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

thisStep="pegasus_nextJob_seurat.txt"
nxtStep1="pegasus_nextJob_checkProjectComplete.txt"
nxtStep2="pegasus_nextJob_snpEff.txt"
nxtStep3="pegasus_nextJob_alleleCount.txt"
pbsHome="/home/tgenjetstream/pegasus-pipe/jobScripts"
constants="/home/tgenjetstream/central-pipe/constants/constants.txt"
constantsDir="/home/tgenjetstream/central-pipe/constants"
myName=`basename $0 | cut -d_ -f2`

#declare -a chrGroups=(1:11:17:21 2:10:16:22 3:9:15:18:MT 4:7:14:Y 5:X:13:20 6:8:12:19)
#echo "array below"
#echo "${chrGroups[@]}"
#echo "array above"

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
seuratPath=`grep @@"$recipe"@@ $constants | grep @@SEURATPATH= | cut -d= -f2`
chrList=`grep @@"$recipe"@@ $constants | grep @@CHRLIST= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
jirRequested=`grep @@"$recipe"@@ $constants | grep @@JIRREALIGN= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

qsubFails=0
###
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
	echo "### DNA pair line is $dnaPairLine"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
	usableName=${sampleNames//,/-}

	pair1=`echo $sampleNames | cut -d, -f1`
	pair2=`echo $sampleNames | cut -d, -f2`

	pair1SampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$pair1"'"'`
	pair2SampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$pair2"'"'`
	pair1KitName=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f1`
	pair2KitName=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f1`
	pair1SamName=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f2`
	pair2SamName=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f2`

	 if [ "${jirRequested}" == "no" ] ; then
                normalBamFile=$runDir/$pair1KitName/$pair1/$pair1.proj.md.bam
                tumorBamFile=$runDir/$pair2KitName/$pair2/$pair2.proj.md.bam

                normalBaiFile=${normalBamFile/.bam/.bai}
                tumorBaiFile=${tumorBamFile/.bam/.bai}

                normalPass=$runDir/$pair1KitName/$pair1/$pair1.proj.bam.mdPass
                tumorPass=$runDir/$pair2KitName/$pair2/$pair2.proj.bam.mdPass
        else
                normalBamFile=$runDir/$pair1KitName/$pair1/$pair1.proj.md.jr.bam
                tumorBamFile=$runDir/$pair2KitName/$pair2/$pair2.proj.md.jr.bam

                normalBaiFile=${normalBamFile/.bam/.bai}
                tumorBaiFile=${tumorBamFile/.bam/.bai}

                normalPass=$runDir/$pair1KitName/$pair1/$pair1.proj.md.jr.bam.jointIRPass
                tumorPass=$runDir/$pair2KitName/$pair2/$pair2.proj.md.jr.bam.jointIRPass
        fi

        echo "### normal BAM: $normalBamFile"
        echo "### tumor  BAM: $tumorBamFile"

	if [[ ! -e $normalBamFile || ! -e $normalBaiFile || ! -e $normalPass ]] ; then
		echo "Normal bam, bai, or jointIRPass does not exist"
		((qsubFails++))
		continue
	fi
	if [[ ! -e $tumorBamFile || ! -e $tumorBaiFile || ! -e $tumorPass ]] ; then
		echo "Tumor bam, bai, or jointIRPass does not exist"
		((qsubFails++))
		continue
	fi

	seuratDir="$runDir/seurat"
	if [ ! -d $seuratDir ] ; then
		mkdir $seuratDir
	fi
	mkdir -p $seuratDir/$usableName
	trackName="$runDir/seurat/$usableName/$usableName"
	STEP=0
	#STEP_COUNT=24
	STEP_COUNT=`ls $chrList/*list | wc -l`
	echo "### Submitting to queue to run seurat on $wd"
	while [ ${STEP} -lt $STEP_COUNT ]
	do
		(( STEP++ ))
		if [[ -e ${trackName}_Step${STEP}.seuratInQueue || -e ${trackName}_Step${STEP}.seuratPass || -e ${trackName}_Step${STEP}.seuratFail ]] ; then
			echo "### Seurat caller is already done, failed, or inqueue for ${trackName}"
			continue
		fi

		echo Starting Seurat caller Step${STEP}
		sbatch -n 1 -N 1 --cpus-per-task $nCores -v STEPCOUNT=$STEP_COUNT,GATKPATH=$gatkPath,SEURATPATH=$seuratPath,NXT1=$nxtStep1,NXT2=$nxtStep2,NXT3=$nxtStep3,TRK=$trackName,CHRLIST=$chrList,REF=$ref,STEP=${STEP},NORMAL=$normalBamFile,TUMOR=$tumorBamFile,RUNDIR=$runDir,D=$d $pbsHome/pegasus_seurat.pbs
		if [ $? -eq 0 ] ; then
			touch ${trackName}_Step${STEP}.seuratInQueue
		else
			((qsubFails++))
		fi
		sleep 2
	done
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
