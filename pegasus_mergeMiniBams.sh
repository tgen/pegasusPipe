#!/bin/bash
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

thisStep="pegasus_nextJob_mergeMiniBams.txt"
nxtStep1="pegasus_nextJob_mergeBams.txt"
pbsHome="~/pegasus-pipe/jobScripts"
constants="~/central-pipe/constants/constants.txt"
constantsDir="~/central-pipe/constants"
myName=`basename $0 | cut -d_ -f2`

time=`date +%d-%m-%Y-%H-%M`
echo "Starting $0 at $time"
if [ "$1" == "" ] ; then
	echo "### Please provide runfolder as the only parameter"
	echo "### Exiting!!!"
	exit
fi
runDir=$1
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

dnaAligner=`grep "@@"$recipe"@@" $constants | grep @@DNAALIGNER= | cut -d= -f2`
ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
gatkPath=`grep "@@"$recipe"@@" $constants | grep @@GATKPATH= | cut -d= -f2`
picardPath=`grep "@@"$recipe"@@" $constants | grep @@PICARDPATH= | cut -d= -f2`
known=`grep @@"$recipe"@@ $constants | grep @@SNPS= | cut -d= -f2`
bwaPath=`grep "@@"$recipe"@@" $constants | grep @@BWAPATH= | cut -d= -f2`
irRequested=`grep "@@"$recipe"@@" $constants | grep @@INDELREALIGN= | cut -d= -f2`
rcRequested=`grep "@@"$recipe"@@" $constants | grep @@RECALIBRATE= | cut -d= -f2`
faiFile="$ref.fai"
optionSum=0
if [ "$irRequested" == "yes" ] ; then
	optionSum=`echo "$optionSum + 1" | bc `
fi
if [ "$rcRequested" == "yes" ] ; then
	optionSum=`echo "$optionSum + 2" | bc `
fi
echo "### Option sum is $optionSum"
case $optionSum in
0) echo "### No ir, no rc"
	inputNextExt="bam"
	prevStepExt="dnaAlignPass"
	prevStepFile="bam"
;;
1) echo "### Yes ir, no rc"
	inputNextExt="ir.bam"
	prevStepExt="indelRealignPass"
	prevStepFile="bam"
;;
2) echo "### No ir, yes rc"
	inputNextExt="rc.bam"
	prevStepExt="recalibratePass"
	prevStepFile="bam"
;;
3) echo "### Yes ir, yes rc"
	inputNextExt="ir.rc.bam"
	prevStepExt="recalibratePass"
	prevStepFile="ir.bam"
;;
*) echo "### I shouldnt be here"
;;
esac

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### dnaAlign: $dnaAligner"
d=`echo $runDir | cut -c 2-`

skipLines=1
skipGroup=0
qsubFails=0

###
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
				echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID"
				if [ "$assayID" == "RNA" ] ; then
					echo "### Assay ID is $assayID. Will skip."
					skipGroup=1
				fi
				lastItemIndex=`echo ${#mergeArray[@]}-1 | bc`
				if [ $skipGroup -ne 1 ] ; then
					for (( i=0; i<=$lastItemIndex; i++ ))
					do
						#echo "array with index $i is::: ${mergeArray[$i]}"
						bamPre="$runDir/$kitName/$samName/$samName.preMerge."`printf "%03d" "$i"`
						bamName="$runDir/$kitName/$samName/$samName.preMerge."`printf "%03d" "$i"`".bam"
						read1="$runDir/$kitName/$samName/$samName"`printf "_%03d" "$i"`"_R1.fastq.gz"
						read2="$runDir/$kitName/$samName/$samName"`printf "_%03d" "$i"`"_R2.fastq.gz"
						if [ ! -e $bamName.processSplit ] ; then
							echo "### This bam file does NOT need to be processed in pieces."
							echo "### bam: $bamName"
							continue
						fi
						if [[ ! -e $read1.fastqSplitPass || ! -e $read2.fastqSplitPass ]] ; then
							echo "### Perhaps splitting is not done yet"
							continue
						fi
						d=`echo $runDir | cut -c 2-`
						mergedBamName="$runDir/$kitName/$samName/$samName.preMerge."`printf "%03d" "$i"`".$inputNextExt"
						mergedBamShort="$samName.preMerge."`printf "%03d" "$i"`".$inputNextExt"
						missingMiniBams=0
						bamMiniCount=0
						bamMiniCountListed=0
						bamList=""
						for thisMiniFq in `ls $runDir/$kitName/$samName/splitFastqs/$samName*R1*part.fastq.gz`
						do
							((bamMiniCountListed++))
							echo "### Mini fastq is $thisMiniFq"
							thisMiniR2=`echo $thisMiniFq | sed 's/\(.*\)_R1/\1_R2/'`
							bamMiniPre=${thisMiniFq/.fastq.gz/}
							bamMiniName=$bamMiniPre.bam
							echo "### Mini bam name is $bamMiniName"

							if [[ ! -e $bamMiniPre.$prevStepFile.$prevStepExt || ! -e $bamMiniPre.$inputNextExt ]] ; then
								echo "### This bam file is not done with prev. step. or input file is missing: $prevStepFile.$prevStepExtn"
								echo "### Missing: $bamMiniPre.$prevStepFile.$prevStepExt"
								echo "### Missing: $bamMiniPre.$inputNextExt"
								((missingMiniBamTotal++))
								((missingMiniBams++))
							else
								echo "### File found."
								bamList="I=$bamMiniPre.$inputNextExt $bamList"
								((bamMiniCount++))
							fi
						done #looping around mini fastq inside splitFastq folder
						echo "### Done with $mergedBamShort with $missingMiniBams missing files out of $bamMiniCountListed" 
						if [ $missingMiniBams -eq 0 ] ; then
							echo "### All bams were accounted for $mergedBamShort"
							echo "### Bam files($bamMiniCount): $bamList"
							if [[ -e $mergedBamName.mergeBamPass || -e $mergedBamName.mergeBamInQueue || -e $mergedBamName.mergeBamFail ]] ; then
								echo "### Already passed, inQueue, or failed"
							else
								echo "### Ready to submit to create $mergedBamName"
								sbatch -n 1 -N 1 --cpus-per-task $nCores -v PICARDPATH=$picardPath,SAMTOOLSPATH=$samtoolsPath,CNT=$bamMiniCount,RUNDIR=$runDir,NXT1=$nxtStep1,BAMLIST="$bamList",MERGEDBAM=$mergedBamName,D=$d $pbsHome/pegasus_mergeBams.pbs
								if [ $? -eq 0 ] ; then
									touch $mergedBamName.mergeBamInQueue
								else
									echo "### There was a failure in qsub command"
									((qsubFails++))
								fi
							fi
						else
							echo "### Some files were missing for $mergedBamShort"
							((qsubFails++))
						fi
					done
				else #else of skipGroup 
					echo "### Skipping group."
				fi #end of skipGroup 
			fi
			kitName=`echo $configLine | cut -d= -f2 | cut -d, -f1`
			samName=`echo $configLine | cut -d= -f2 | cut -d, -f2`
			assayID=`echo $configLine | cut -d= -f2 | cut -d, -f3`
			libraID=`echo $configLine | cut -d= -f2 | cut -d, -f4`
			unset mergeArray
			count=0
			skipGroup=0
			continue
		else #doesnt start with =, add to mergeArray
			mergeArray[$count]=$configLine	
			((count++))
		fi
	else
		continue
	fi
done

if [ $qsubFails -eq 0 ] ; then
#all jobs submitted succesffully, remove this dir from messages
	echo "### I should remove $thisStep from $runDir."
	rm -f $runDir/$thisStep
else
#qsub failed at some point, this runDir must stay in messages
	echo "### Failure in qsub. Not touching $thisStep"
fi

time=`date +%d-%m-%Y-%H-%M`
echo "Ending $0 at $time"
