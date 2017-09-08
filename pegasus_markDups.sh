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

thisStep="pegasus_nextJob_markDups.txt"
nxtStep1="pegasus_nextJob_jointIR.txt"
nxtStep2="pegasus_nextJob_haplotypeCaller.txt"
nxtStep3="pegasus_nextJob_picardMultiMetrics.txt"
nxtStep4="pegasus_nextJob_samtoolsStats.txt"
nxtStep5="pegasus_nextJob_picardHSMetrics.txt"
nxtStep6="pegasus_nextJob_picardGcBiasMetrics.txt"
nxtStep7="pegasus_nextJob_snpSniff.txt"
nxtStep8="pegasus_nextJob_freebayes.txt"
nxtStep9="pegasus_nextJob_samtoolsMpileUp.txt"
nxtStep10="pegasus_nextJob_ancestry.txt"

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
picardPath=`grep @@"$recipe"@@ $constants | grep @@PICARDPATH= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### dnaAlign: $dnaAligner"
d=`echo $runDir | cut -c 2-`

skipLines=1
skipGroup=0
qsubFails=0
###
for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
	echo "### Sample is $sampleLine"
	kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
	samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
	assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
	libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
	#echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID"
	bamName="$runDir/$kitName/$samName/$samName.proj.bam"
	mdBamFile="$runDir/$kitName/$samName/$samName.proj.md.bam"
	if [ "$assayID" == "RNA" ] ; then
		echo "### Assay ID is $assayID. Will skip."
		continue
	fi
	if [[ -e $bamName.mdInQueue || -e $bamName.mdPass || -e $bamName.mdFail ]] ; then
		echo "### This bam already in queue, passed, or failed for mark dups: $bamName"
		continue
	fi
	if [[ ! -e $bamName.mergeBamPass || ! -e $bamName ]] ; then
		echo "### Either mergeBamPass or the bam itself is missing"
		echo "### Missing: $bamName.mergeBamPass"
		echo "### Missing: $bamName"
		((qsubFails++))
		continue
	fi
	echo "### Submitting to queue to mark duplicates: $bamName"
	sbatch -n 1 -N 1 --cpus-per-task $nCores --export PICARDPATH=$picardPath,BAMFILE=$bamName,NXT1=$nxtStep1,NXT2=$nxtStep2,NXT3=$nxtStep3,NXT4=$nxtStep4,NXT5=$nxtStep5,NXT6=$nxtStep6,NXT7=$nxtStep7,NXT8=$nxtStep8,NXT9=$nxtStep9,NXT10=$nxtStep10,RUNDIR=$runDir,OUTPUTBAM=$mdBamFile,D=$d $pegasusPbsHome/pegasus_markDups.pbs
	if [ $? -eq 0 ] ; then
		touch $bamName.mdInQueue
	else
		((qsubFails++))
	fi
	sleep 1

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
