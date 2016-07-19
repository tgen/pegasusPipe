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

thisStep="pegasus_nextJob_unifiedGenotyper.txt"
#nxtStep1="pegasus_nextJob_vqsr.txt"
nxtStep1="pegasus_nextJob_snpEff.txt"
pbsHome="/home/mrussell/pegasus-pipe/jobScripts"
constants="/home/mrussell/central-pipe/constants/constants.txt"
constantsDir="/home/mrussell/central-pipe/constants"
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


ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
snps=`grep @@"$recipe"@@ $constants | grep @@SNPS= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0

#DNAPAIR samples have multi sample gatkUG calls and seurat calls
#DNAFAMI samples have multi sample gatkUG calls
#samples in neither group have single sample gatkUG calls 

for dnaPairLine in `cat $configFile | grep '^DNAPAIR=\|^DNAFAMI='`
do
	echo "### DNA pair line is $dnaPairLine"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
	usableName=${sampleNames//,/-}
	sampleCount=0
	missingSampleCount=0
	sampleList=""
	for eachSample in ${sampleNames//,/ }
	do
		((sampleCount++))
		#echo "eachsample: $eachSample"
		sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
		kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
		samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
		eachSampleBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
		eachSamplePass=$runDir/$kitName/$samName/$samName.proj.md.jr.bam.jointIRPass
		if [[ ! -e $eachSampleBam || ! -e $eachSamplePass ]] ; then
			echo "### Can't find the jr bam or the jr pass"
			echo "### BAM: $eachSampleBam"
			echo "### PAS: $eachSamplePass"
			((missingSampleCount++))
		else
			sampleList="-I $eachSampleBam $sampleList"		
		fi
	done
	if [ $missingSampleCount -eq 0 ] ; then
		echo "### All sample bams are found"
		sampleList="${sampleList%?}"
		#echo "### Sample list: $sampleList"
	else
		echo "### There were missing things"
		echo "### Samples missing: $missingSampleCount/$sampleCount"
		((qsubFails++))
		continue
	fi
	ugDir="$runDir/ug"
	if [ ! -d $ugDir ] ; then
		mkdir $ugDir
	fi
	mkdir -p $ugDir/$usableName
	workDir="$ugDir/$usableName"
	trackName="$runDir/ug/$usableName/$usableName"
	if [[ -e $trackName.ugInQueue || -e $trackName.ugPass || -e $trackName.ugFail ]] ; then 
		echo "### Joint unified genotyper is already done, failed or inQueue"
		continue
	fi 
	echo "### Submitting $usableName to queue for joint unified genotyper..."
	qsub -A $debit -l nodes=1:ppn=$nCores -v GATKPATH=$gatkPath,TRK=$trackName,KNOWN=$snps,BAMLIST="'$sampleList'",REF=$ref,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pbsHome/pegasus_unifiedGenotyper.pbs
	if [ $? -eq 0 ] ; then
		touch $trackName.ugInQueue
	else
		((qsubFails++))
	fi
	sleep 2
done
for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
	echo "sample is $sampleLine"
	kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
	samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
	assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
	libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
	echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID"
	if [ "$assayID" != "RNA" ] ; then
		echo "### Assay ID is $assayID. Must be genome or exome."
		pcDir=$runDir/$kitName/$samName
		inBam=$runDir/$kitName/$samName/$samName.proj.bam
		mdBam=$runDir/$kitName/$samName/$samName.proj.md.bam
		jrBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
		jrRequested=`cat $configFile | grep '^DNAPAIR=\|^DNAFAMI=' | grep $samName | wc -l`
		if [ $jrRequested -gt 0 ] ; then
			echo "### Unified genotyper will not operate on a single bam because Joint IR is requested for $samName"
		else
			echo "### Unified genotyper will run on single bam because Joint IR is NOT requested for $samName"
			if [[ ! -e $inBam.mdPass || ! -e $mdBam ]] ; then
				echo "### Either mdPass or the bam itself is missing for $mdBam"
				((qsubFails++))
			else
				if [[ -e $mdBam.ugPass || -e $mdBam.ugInQueue || -e $mdBam.ugFail ]] ; then
					echo "### Unified genotyper already passed, in queue, or failed for $mdBam"
				else
					echo "### Submitting for single bam unified genotyper: $mdBam"
					qsub -A $debit -l nodes=1:ppn=$nCores -v GATKPATH=$gatkPath,TRK=$mdBam,KNOWN=$snps,BAMLIST=$mdBam,REF=$ref,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pbsHome/pegasus_unifiedGenotyperSingle.pbs
					if [ $? -eq 0 ] ; then
						touch $mdBam.ugInQueue
					else
						((qsubFails++))
					fi
					sleep 2
				fi
			fi
		fi
	else
		echo "### Assay ID is $assayID. Must be RNA."
		#code for calling AS metrics on tophat bams here
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
