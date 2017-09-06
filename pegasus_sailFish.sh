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

thisStep="pegasus_nextJob_sailFish.txt"
nxtStep1="pegasus_nextJob_postSailFish.txt"
#nxtStep2="pegasus_nextJob_cuffDiff.txt"
#nxtStep3="pegasus_nextJob_htSeq.txt"
#nxtStep4="pegasus_nextJob_picardRNAMetrics.txt"
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
sailFishGTF=`grep @@SAILFISHGTF= $constantsDir/$recipe | cut -d= -f2`
ccdsGTF=`grep @@CCDSGTF= $constantsDir/$recipe | cut -d= -f2`
sailFishIndexDir=`grep @@SAILFISHINDEXDIR= $constantsDir/$recipe | cut -d= -f2`
sailFishCCDSIndex=`grep @@SAILFISHCCDSINDEX= $constantsDir/$recipe | cut -d= -f2`
sailFishPath=`grep @@SAILFISHPATH= $constantsDir/$recipe | cut -d= -f2`
sailFish=`grep @@SAILFISH= $constantsDir/$recipe | cut -d= -f2`

ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
indexbase=`grep @@"$recipe"@@ $constants | grep @@BOWTIE2INDEX= | cut -d= -f2`
transIndex=`grep @@"$recipe"@@ $constants | grep @@TRANSINDEX= | cut -d= -f2`
tophat2Path=`grep @@"$recipe"@@ $constants | grep @@TOPHAT2PATH= | cut -d= -f2`
bowtie2Path=`grep @@"$recipe"@@ $constants | grep @@BOWTIE2PATH= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
picardPath=`grep @@"$recipe"@@ $constants | grep @@PICARDPATH= | cut -d= -f2`
starPath=`grep @@"$recipe"@@ $constants | grep @@STARPATH= | cut -d= -f2`
bwaPath=`grep "@@"$recipe"@@" $constants | grep @@BWAPATH= | cut -d= -f2`
refPreTophat=`grep @@"$recipe"@@ $constants | grep @@REFPRETOPHAT= | cut -d= -f2`
faiFile="$ref.fai"
echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###
for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
	if [ $sailFish != "yes" ] ; then
		echo "### Sailfish not requested for this recipe"
		continue
	fi
	echo "### Sample is $sampleLine"
	kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
	samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
	assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
	libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
	rnaStrand=`grep "@@"$kitName"@@" $constants | cut -d= -f2`
	echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID, rnaStrand: $rnaStrand"
	if [ "$assayID" != "RNA" ] ; then
		echo "### Assay ID is $assayID. Skipping."
		continue
	fi
	read1Name=$runDir/$kitName/$samName/$samName.proj.R1.fastq.gz
	read2Name=$runDir/$kitName/$samName/$samName.proj.R2.fastq.gz
	echo "read 1 name: $read1Name"
	echo "read 2 name: $read2Name"
	echo "### Aligner for recipe $recipe is $rnaAligner"
	thisPath=`dirname $runDir`
	cd $thisPath
	ownDir=${read1Name/.proj.R1.fastq.gz/.sailFishDir}
	if [[ ! -e $read1Name || ! -e $read2Name || ! -e $read1Name.mergeFastqPass || ! -e $read2Name.mergeFastqPass ]] ; then
		echo "### one of the fastq files or read pass files dont exist"
		echo "### read1Pass: $read1Name.mergeFastqPass"
		echo "### read2Pass: $read2Name.mergeFastqPass"
		((qsubFails++))
		continue
	fi
	if [[ -e $ownDir.sailFishPass || -e $ownDir.sailFishFail || -e $ownDir.sailFishInQueue ]] ; then 
		echo "### Sail Fish already done, failed or inQueue"
		continue
	fi 
	if [ ! -d $ownDir ] ; then
		mkdir -p $ownDir
	fi

	#creating linked files to the original reads
	r1Name=`basename $read1Name`
	r2Name=`basename $read2Name`
	cd $ownDir
	ln -s $read1Name $r1Name
	read1Name=$ownDir/$r1Name
	ln -s $read2Name $r2Name
	read2Name=$ownDir/$r2Name
	cd -
	#done creating links. vars for reads changed.
	echo "### read 1 name: $read1Name"
	echo "### read 2 name: $read2Name"
	#lineLength=`gunzip -c $read1Name | head -2 | tail -1 | wc -c` 
	#let "readLength=$lineLength-1"
	#echo "### Read length determined to be $readLength for $ownDir"
	#refGrep="STARREF"$readLength
	#starRef=`grep "@@"$recipe"@@" $constants | grep @@"$refGrep"= | cut -d= -f2`
	#echo "### Star reference is $starRef"
	echo "### submitting $ownDir to queue for Sail Fish... "
	echo "### sail fish path is: $sailFishPath"
#	if [[ $rnaStrand == "FIRST" || $rnaStrand == "SECOND" ]] ; then
#        	echo "##running stranded sailFish case"
#		sbatch -n 1 -N 1 --cpus-per-task $nCores -v SAILFISHPATH=$sailFishPath,SAMPLE=$samName,RNASTRAND=$rnaStrand,SAILFISHGTF=$sailFishGTF,CCDSGTF=$ccdsGTF,SAILFISHINDEXDIR=$sailFishIndexDir,SAILFISHCCDSINDEX=$sailFishCCDSIndex,FASTQ1=$read1Name,FASTQ2=$read2Name,DIR=$ownDir,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pbsHome/pegasus_strandedSailFish.pbs
#		if [ $? -eq 0 ] ; then
#			touch $ownDir.sailFishInQueue
#		else
#			((qsubFails++))
#		fi
#		sleep 2
#	else
		echo "##running unstranded Sail Fish case"
		sbatch -n 1 -N 1 --cpus-per-task $nCores -v SAILFISHPATH=$sailFishPath,SAMPLE=$samName,SAILFISHGTF=$sailFishGTF,CCDSGTF=$ccdsGTF,SAILFISHINDEXDIR=$sailFishIndexDir,SAILFISHCCDSINDEX=$sailFishCCDSIndex,FASTQ1=$read1Name,FASTQ2=$read2Name,DIR=$ownDir,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pbsHome/pegasus_sailFish.pbs
		if [ $? -eq 0 ] ; then
			touch $ownDir.sailFishInQueue
		else
			((qsubFails++))
		fi
		sleep 2

#	fi
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
