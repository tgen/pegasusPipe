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

thisStep="pegasus_nextJob_splitFastqs.txt"
nxtStep1="pegasus_nextJob_dnaAlignParts.txt"

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
echo "### projName: $projName"
echo "### confFile: $configFile"
if [ ! -e $configFile ] ; then
	echo "### Config file not found at $configFile!!!"
	echo "### Exiting!!!"
	exit
else
	echo "### Config file found."
fi
debit=`cat $configFile | grep "^DEBIT=" | cut -d= -f2 | head -1 | tr -d [:space:]`
recipe=`cat $configFile | grep "^RECIPE=" | cut -d= -f2 | head -1 | tr -d [:space:]`

#commented this out until its fixed
#nCores=`grep @@${myName}_CORES= $constantsDir/$recipe | cut -d= -f2`
nCores=8

incFastq=`grep "@@"$recipe"@@" $constants | grep @@INCFASTQ= | cut -d= -f2`

d=`echo $runDir | cut -c 2-`
skipLines=1
count=0
sampleCount=0

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
				((sampleCount++))
				#echo "### Sample with $arrayCount rows found for kit: $kitName, sample: $samName."
				lastItemIndex=`echo ${#mergeArray[@]}-1 | bc`
				for (( i=0; i<=$lastItemIndex; i++ ))
				do
					#echo "array with index $i is::: ${mergeArray[$i]}"
					targetName="$runDir/$kitName/$samName/$samName"`printf "_%03d" "$i"`"_R1.fastq.gz"
					targR2Name="$runDir/$kitName/$samName/$samName"`printf "_%03d" "$i"`"_R2.fastq.gz"
					splitDir=$runDir/$kitName/$samName/splitFastqs
					mkdir -p $splitDir
					if [[ -e $targetName.needsToBeSplit ]] ; then
						echo "### Needs to be split."
						if [[ -e $targetName.fastqSplitInQueue || -e $targetName.fastqSplitPass || -e $targetName.fastqSplitFail ]] ; then
							echo "### This fastq file already in queue, failed, or passed for splitting"
						else
							echo "### Submitting to split: $targetName"		
							prefix="$samName"`printf "_%03d" "$i"`"_R1.split.fastq.gz"
							sbatch -n 1 -N 1 --cpus-per-task $nCores -v PF=$prefix,DIR=$splitDir,D=$d,FQ=$targetName,NXT1=$nxtStep1,RUNDIR=$runDir $pegasusPbsHome/pegasus_splitFastq.pbs
							if [ $? -eq 0 ] ; then
								touch $targetName.fastqSplitInQueue
							else
								((qsubFails++))
							fi
							sleep 1
						fi
						if [[ -e $targR2Name.fastqSplitInQueue || -e $targR2Name.fastqSplitPass || -e $targR2Name.fastqSplitFail ]] ; then
							echo "### This fastq file already in queue, failed, or passed for splitting"
						else
							echo "### Submitting to split: $targR2Name"		
							prefix="$samName"`printf "_%03d" "$i"`"_R2.split.fastq.gz"
							sbatch -n 1 -N 1 --cpus-per-task $nCores -v PF=$prefix,DIR=$splitDir,D=$d,FQ=$targR2Name,NXT1=$nxtStep1,RUNDIR=$runDir $pegasusPbsHome/pegasus_splitFastq.pbs
							if [ $? -eq 0 ] ; then
								touch $targR2Name.fastqSplitInQueue
							else
								((qsubFails++))
							fi
							sleep 1
						fi
					else
						echo "### Does not need to be split."
					fi
				done
			fi
			kitName=`echo $configLine | cut -d= -f2 | cut -d, -f1`
			samName=`echo $configLine | cut -d= -f2 | cut -d, -f2`
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
