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

thisStep="pegasus_nextJob_mergeJointIRBams.txt"
nxtStep1="pegasus_nextJob_seurat.txt"
nxtStep2="pegasus_nextJob_clonalCov.txt"
nxtStep3="pegasus_nextJob_trn.txt"
nxtStep4="pegasus_nextJob_strelka.txt"
nxtStep5="pegasus_nextJob_mutect.txt"
pbsHome="~/pegasus-pipe/jobScripts"
constants="~/central-pipe/constants/constants.txt"
constantsDir="~/central-pipe/constants"
myName=`basename $0 | cut -d_ -f2`

declare -a chrGroups=(1:11:17:21 2:10:16:22 3:9:15:18:MT 4:7:14:Y 5:X:13:20 6:8:12:19)
declare -a chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)

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


ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
indels=`grep @@"$recipe"@@ $constants | grep @@INDELS= | cut -d= -f2`
picardPath=`grep @@"$recipe"@@ $constants | grep @@PICARDPATH= | cut -d= -f2`
samtoolsPath=`grep @@"$recipe"@@ $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###
#for dnaPairLine in `cat $configFile | grep '^DNAPAIR=\|^DNAFAMI='`
for dnaPairLine in `cat $configFile | grep '^JIRSET='`
do
	echo "### DNA pair line is $dnaPairLine"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
	usableName=${sampleNames//,/-}
	sampleCount=0
	missingSampleCount=0
	sampleList=""
	bigBamDetected=0
	jirDir="$runDir/jointIR"
	workDir="$jirDir/$usableName"
	irIntFile="$runDir/jointIR/$usableName/$usableName.intervals"
	trackName="$runDir/jointIR/$usableName/$usableName"

	for eachSample in ${sampleNames//,/ }
	do
		((sampleCount++))
		#echo "eachsample: $eachSample"
		sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
		kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
		samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
		eachSampleBam=$runDir/$kitName/$samName/$samName.proj.md.bam
		eachSamplePass=$runDir/$kitName/$samName/$samName.proj.bam.mdPass
		bamSize=`stat -c%s $eachSampleBam`
		if [ $bamSize -ge 36000000000 ] ; then
			echo "### This bam is greater than 40GB. Setting flag."
			bigBamDetected=1
		fi
	done
	if [ $bigBamDetected -eq 0 ] ; then
		echo "### None of the bams must be bigger than 40GB."
		echo "### Nothing to do."
	else
		echo "### One of the bams must be bigger than 40GB."
		echo "Must check if all 25 chrs in 6 groups have passed"
		missingGrp=0
                for chrGrp in "${chrGroups[@]}"
                do
			grpName=`echo $chrGrp | cut -d: -f1`
			if [ ! -e $trackName.jointIR-group$grpName-Pass ] ; then 
				echo "### Joint indel realign on big bam group $chrGrp is not done!!!"
				echo "### Missing $trackName.jointIR-group$grpName-Pass"
				((missingGrp++))
			fi 
		done
		if [ $missingGrp -eq 0 ] ; then
			echo "### All 6 groups must have passed for joint indel realign."
			for eachSample in ${sampleNames//,/ }
			do
				echo "### Preparing to submit mergin fo $eachSample"	
				((sampleCount++))
				mergedBam=$workDir/$eachSample.proj.md.jr.bam
				if [[ -e $mergedBam.mergeBamPass || -e $mergedBam.mergeBamFail || -e $mergedBam.mergeBamInQueue ]] ; then
					echo "### Merging of joint indel realigned bam pieces is already done, failed, or in queue."
					continue
				fi
				sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
				kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
				samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
				newLoc=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
				fileList=""
				for thisChr in "${chrs[@]}"
				do
					thisBam=$workDir/$eachSample.proj.md.chr$thisChr.jr.bam
					if [ ! -e $thisBam ] ; then
						echo "### I thought we were good? Missing $thisBam"
					fi
					fileList="$fileList I=$thisBam"
				done
				sbatch -n 1 -N 1 --cpus-per-task $nCores -v NXT1=$nxtStep1,NXT2=$nxtStep2,NXT3=$nxtStep3,NXT4=$nxtStep4,NXT5=$nxtStep5,NEWLOC=$newLoc,PICARDPATH=$picardPath,SAMTOOLSPATH=$samtoolsPath,CNT=25,MERGEDBAM=$mergedBam,BAMLIST="$fileList",RUNDIR=$runDir,D=$d $pbsHome/pegasus_mergeBamsForBigJIR.pbs
				if [ $? -eq 0 ] ; then
					touch $mergedBam.mergeBamInQueue
				else
					((qsubFails++))
				fi
				sleep 2

				done
		else
			echo "### Not all chr groups have passed seurat"
			((qsubFails++))
		fi 
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
