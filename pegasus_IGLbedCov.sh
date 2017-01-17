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

thisStep="pegasus_nextJob_IGLbedCov.txt"
nxtStep1="pegasus_nextJob_postPIGLbedCov.txt"
pbsHome="/home/tgenjetstream/pegasus-pipe/jobScripts"
constants="/home/tgenjetstream/central-pipe/constants/constants.txt"
constantsDir="/home/tgenjetstream/central-pipe/constants"
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
iglBedCov=`grep @@IGLBEDCOV= $constantsDir/$recipe | cut -d= -f2`


ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
picardPath=`grep @@"$recipe"@@ $constants | grep @@PICARDPATH= | cut -d= -f2`
samtoolsPath=`grep @@"$recipe"@@ $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
bedtoolsPath=`grep @@"$recipe"@@ $constants | grep @@BEDTOOLSPATH= | cut -d= -f2`
iglListBed=`grep @@"$recipe"@@ $constants | grep @@IGLLISTBED= | cut -d= -f2`
chrListBed=`grep @@"$recipe"@@ $constants | grep @@CHRLISTBED= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

qsubFails=0
echo "### rnaAlign: $rnaAligner"
###
for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
	if [ $iglBedCov != "yes" ] ; then
		echo "### iglBedCov not requested for this recipe"
		continue
	fi

	echo "### Sample is $sampleLine"
	kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
	samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
	assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
	libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
	echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID"
	if [ "$assayID" != "RNA" ] ; then
		echo "### Assay ID is $assayID. Skipping."
		continue
	fi
	case $rnaAligner in 
	tophat) echo "### Tophat case"
		rnaDir="$runDir/$kitName/$samName/$samName.topHatDir"
		rnaBam="$rnaDir/$samName.proj.accepted_hits.md.bam"
		rnaPass="$rnaDir/$samName.proj.accepted_hits.bam.rnaMarkDupPass"
		;;
	star) echo "### Star case"
		rnaDir="$runDir/$kitName/$samName/$samName.starDir"
		rnaBam="$rnaDir/$samName.proj.Aligned.out.sorted.md.bam"
		rnaPass="$rnaDir/$samName.proj.Aligned.out.sorted.bam.rnaMarkDupPass"
		;;
	anotherRNAaligner) echo "### Anoter RNA aligner"
		;;
	*) echo "### I should not be here"
		;;
	esac

	echo "### My bam is $rnaBam"
	ownDir=`dirname $rnaBam`
	if [ ! -e $rnaPass ] ; then
		echo "### Looks like rna align is not done yet. $rnaPass doesnt exist yet"
		((qsubFails++))
		continue
	fi
	if [ ! -e $rnaBam ] ; then
		echo "### Weird. Bam itself is missing: $rnaBam"
		((qsubFails++))
		continue
	fi

	if [[ -e $rnaBam.IGLbedCovPass || -e $rnaBam.IGLbedCovFail || -e $rnaBam.IGLbedCovInQueue ]] ; then 
		echo "### IGL bed cov is already done, failed or inQueue"
		continue
	fi 
	if [ ! -d $runDir/stats ] ; then
		mkdir $runDir/stats
	fi
	echo "### Submitting $rnaBam to queue for picard RNA Metrics..."
	qsub -A $debit -l nodes=1:ppn=$nCores -v DIR=$ownDir,SAMPLE=$samName,CHRLISTBED=$chrListBed,REF=$ref,IGLLISTBED=$iglListBed,BEDTOOLSPATH=$bedtoolsPath,SAMTOOLSPATH=$samtoolsPath,BAMFILE=$rnaBam,RUNDIR=$runDir,D=$d $pbsHome/pegasus_IGLbedCov.pbs
	if [ $? -eq 0 ] ; then
		touch $rnaBam.IGLbedCovInQueue
	else
		((qsubFails++))
	fi
	sleep 2
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
