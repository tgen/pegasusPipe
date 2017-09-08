#!/usr/bin/env bash
#####################################################################
# Copyright (c) 2011 by The Translational Genomics Research
# Institute. All rights reserved. This License is limited to, and you may
# use the Software solely for, your own internal and non-commercial use
# for academic and research purposes. Without limiting the foregoing, you
# may not use the Software as part of, or in any way in connection with
# the production, marketing, sale or support of any commercial product or
# service or for any governmental purposes. For commercial or governmental
# use, please contact dcraig@tgen.org. By installing this Software you are
# agreeing to the terms of the LICENSE file distributed with this
# software.
#####################################################################

thisStep="pegasus_nextJob_mergeFastqs.txt"
nxtStep1="pegasus_nextJob_rnaAlign.txt"
nxtStep2="pegasus_nextJob_detectFusion.txt"
nxtStep3="pegasus_nextJob_salmon.txt"
nxtStep4="pegasus_nextJob_kallisto.txt"

constants=~/jetstream/constants/constants.txt
constantsDir=~/jetstream/constants/
myName=`basename $0 | cut -d_ -f2`

time=`date +%d-%m-%Y-%H-%M`
echo "Starting $0 at $time"
if [ "$1" == "" ] ; then
	echo "### Please provide runfolder as the only parameter"
	echo "### Exiting!!!"
	exit
fi
runDir=$1
qsubFails=0
projName=`basename $runDir | awk -F'_ps20' '{print $1}'`
configFile=$runDir/$projName.config
if [ ! -e $configFile ] ; then
	echo "### Config file not found at $configFile!!!"
	echo "### Exiting!!!"
	exit
else
	echo "### Config file found."
fi

recipe=`cat $configFile | grep "^RECIPE=" | cut -d= -f2 | head -1 | tr -d [:space:]`
debit=`cat $configFile | grep "^DEBIT=" | cut -d= -f2 | head -1 | tr -d [:space:]`

nCores=`grep @@${myName}_CORES= $constantsDir/$recipe | cut -d= -f2`


echo "### projName: $projName"
echo "### confFile: $configFile"
incFastq=`grep "@@"$recipe"@@" $constants | grep @@INCFASTQ= | cut -d= -f2`
#incFastq=`cat $configFile | grep ^INCFASTQ | cut -d= -f2`
#for fqLine in `cat $configFile | grep ^FQ=`
#do
#	fqFile=`echo $fqLine | cut -d= -f2 | cut -d, -f2`
#	r2File=${fqFile/_R1/_R2}
#done
#if [[ "$incFastq" != "yes" && "$incFastq" != "no" ]] ; then
#	echo "### Valid values for INCFASTQ is either yes or no"
#	echo "### Exiting!!!"
#	exit
#else
#	echo "### Inc fastq is good: $incFastq"
#fi 
skipLines=1
count=0
sampleCount=0
missingFastqTotal=0
d=`echo $runDir | cut -c 2-`

for configLine in `cat $configFile`
do
	if [ "$configLine" == "=START" ] ; then
		skipLines=0
		continue
	fi
	if [ $skipLines -eq 0 ] ; then
		if [[ $configLine == SAMPLE* || $configLine == =END* ]] ; then
			#echo "config line is $configLine"
			arrayCount=${#mergeArray[@]}
			if [ $arrayCount -gt 0 ] ; then
				echo "### Starting with $samName"
				missingFastqSample=0
				accountFastqSample=0
				fastqList1=""
				fastqList2=""
				read1Count=0
				read2Count=0
				((sampleCount++))
				#echo "### Sample with $arrayCount rows found for kit: $kitName, sample: $samName."
				lastItemIndex=`echo ${#mergeArray[@]}-1 | bc`
				for (( i=0; i<=$lastItemIndex; i++ ))
				do
					#echo "array with index $i is::: ${mergeArray[$i]}"
					((accountFastqSample++));((accountFastqSample++))
					thisFq=`echo ${mergeArray[$i]} | cut -d= -f2 | cut -d, -f2`
					r2File=${thisFq/_R1/_R2}
					sourceName="$runDir/$kitName/$samName/$samName"`printf "_%03d" "$i"`"_R1.fastq.gz"
					targetName="$runDir/$kitName/$samName/$samName.proj.R1.fastq.gz"
					sourR2Name="$runDir/$kitName/$samName/$samName"`printf "_%03d" "$i"`"_R2.fastq.gz"
					targR2Name="$runDir/$kitName/$samName/$samName.proj.R2.fastq.gz"
					#echo "target name is $sourceName"
					if [[ ! -e $sourceName.cpFqPass || ! -e $sourceName ]] ; then
						echo "### File missing: $sourceName or its .cpFqPass file"
						((missingFastqTotal++))
						((missingFastqSample++))
					else
						echo "### File found $sourceName "	
						fastqList1="$sourceName $fastqList1"
						((read1Count++))
					fi
					if [[ ! -e $sourR2Name.cpPass && ! -e $sourR2Name ]] ; then
						echo "### File missing: $sourR2Name or its .cpFqPass file"
						((missingFastqTotal++))
						((missingFastqSample++))
					else
						echo "### File found $sourR2Name "	
						fastqList2="$sourR2Name $fastqList2"
						((read2Count++))
					fi
				done
				echo "### Done with $samName with $missingFastqSample missing files out of $accountFastqSample" 
				if [ $missingFastqSample -eq 0 ] ; then
					if [ "$assayID" == "RNA" ] ; then
						echo "### Assay ID is RNA"
						echo "### All fastqs were accounted for $samName"
						echo "### R1 files($read1Count): $fastqList1"
						echo "### R2 files($read2Count): $fastqList2"
						if [[ -e $targetName.mergeFastqPass || -e $targetName.mergeFastqInQueue || -e $targetName.mergeFastqFail ]] ; then
							echo "### Already passed, inQueue, or failed"
						else
							echo "### Ready to submit for read 1..."
							sbatch -n 1 -N 1 --cpus-per-task $nCores -v CNT=$read1Count,RUNDIR=$runDir,NTX1=$nxtStep1,NXT2=$nxtStep2,NXT3=$nxtStep3,NXT4=$nxtStep4,FASTQLIST="$fastqList1",MERGEDFASTQ=$targetName,D=$d $pegasusPbsHome/pegasus_mergeFastqs.pbs
							if [ $? -eq 0 ] ; then
								touch $targetName.mergeFastqInQueue
							else
								echo "### There was a failure in qsub command"
								((qsubFails++))
							fi
						fi
						if [[ -e $targR2Name.mergeFastqPass || -e $targR2Name.mergeFastqInQueue || -e $targR2Name.mergeFastqFail ]] ; then
							echo "### Already passed, inQueue, or failed"
						else
							echo "### Ready to submit for read 2..."
							sbatch -n 1 -N 1 --cpus-per-task $nCores -v CNT=$read1Count,RUNDIR=$runDir,NXT1=$nxtStep1,NXT2=$nxtStep2,NXT3=$nxtStep3,NXT4=$nxtStep4,FASTQLIST="$fastqList2",MERGEDFASTQ=$targR2Name,D=$d $pegasusPbsHome/pegasus_mergeFastqs.pbs
							if [ $? -eq 0 ] ; then
								touch $targR2Name.mergeFastqInQueue
							else
								echo "### There was a failure in qsub command"
								((qsubFails++))
							fi
						fi
					else
						echo "### Assay ID is not RNA. $assayID"
					fi
				else
					echo "### Some files were missing for $samName"
				fi
			fi
			kitName=`echo $configLine | cut -d= -f2 | cut -d, -f1`
			samName=`echo $configLine | cut -d= -f2 | cut -d, -f2`
			assayID=`echo $configLine | cut -d= -f2 | cut -d, -f3`
			unset mergeArray
			count=0
			continue
		else #doesnt start with =, add to mergeArray
			#echo "adding $configLine to mergeArray"
			mergeArray[$count]=$configLine	
			((count++))
		fi
	else
		continue
	fi
done

if [ $qsubFails -eq 0 ] ; then
#all jobs submitted succesffully, remove this dir from messages
	echo "### Removing $thisStep from $runDir."
	rm -f $runDir/$thisStep
else
#qsub failed at some point, this runDir must stay in messages
	echo "### Failure in qsub. Keeping $thisStep"
fi

time=`date +%d-%m-%Y-%H-%M`
echo "Ending $0 at $time"
