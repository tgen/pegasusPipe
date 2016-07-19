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

thisStep="pegasus_nextJob_picardGcBiasMetrics.txt"
nxtStep1="pegasus_nextJob_postPicGcBiasMetrics.txt"
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
thFusionRef=`grep "@@"$recipe"@@" $constants | grep @@THFUSIONREF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
dnaAligner=`grep "@@"$recipe"@@" $constants | grep @@DNAALIGNER= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
picardPath=`grep "@@"$recipe"@@" $constants | grep @@PICARDPATH= | cut -d= -f2`
jirRequested=`grep @@"$recipe"@@ $constants | grep @@JIRREALIGN= | cut -d= -f2`
echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0

if [ ! -d $runDir/stats ] ; then
	mkdir -p $runDir/stats
fi

###
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
		#if [[ ! -e $inBam.mergeBamPass || ! -e $inBam ]] ; then
		#	echo "### Either mergeBamPass or the bam itself is missing for $inBam"
		#	((qsubFails++))
		#else
		#	if [[ -e $inBam.picGcBiasMetricsPass || -e $inBam.picGcBiasMetricsInQueue || -e $inBam.picGcBiasMetricsFail ]] ; then
		#		echo "### Picard GC Bias summary metric already passed, in queue, or failed for $inBam"
		#	else
		#		echo "### Submitting for picard AS Metrics: $inBam"
		#		qsub -A $debit -l nodes=1:ppn=$nCores -v PICARDPATH=$picardPath,RUNDIR=$runDir,REF=$ref,BAMFILE=$inBam,DIR=$pcDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_picardGcBiasMetrics.pbs
		#		if [ $? -eq 0 ] ; then
		#			touch $inBam.picGcBiasMetricsInQueue
		#		else
		#			((qsubFails++))
		#		fi
		#		sleep 2
#
#			fi
#		fi
		if [[ ! -e $inBam.mdPass || ! -e $mdBam ]] ; then
			echo "### Either mdPass or the bam itself is missing for $mdBam"
			((qsubFails++))
		else
			if [[ -e $mdBam.picGcBiasMetricsPass || -e $mdBam.picGcBiasMetricsInQueue || -e $mdBam.picGcBiasMetricsFail ]] ; then
				echo "### Picard GC Bias summary metric already passed, in queue, or failed for $mdBam"
			else
				echo "### Submitting for picard GC Bias Metrics: $mdBam"
				qsub -A $debit -l nodes=1:ppn=$nCores -v PICARDPATH=$picardPath,RUNDIR=$runDir,REF=$ref,BAMFILE=$mdBam,DIR=$pcDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_picardGcBiasMetrics.pbs
				if [ $? -eq 0 ] ; then
					touch $mdBam.picGcBiasMetricsInQueue
				else
					((qsubFails++))
				fi
				sleep 2
			fi
		fi
		jrRequested=`cat $configFile | grep '^DNAPAIR=\|^DNAFAMI=' | grep $samName | wc -l`
		if [ $jrRequested -gt 0 ] && [ $jirRequested != "no" ] ; then
			echo "### Joint IR requested for $samName"
			if [[ ! -e $jrBam.jointIRPass || ! -e $jrBam ]] ; then
				echo "### Either jointIRPass or the bam itself is missing for $jrBam"
				((qsubFails++))
			else
				if [[ -e $jrBam.picGcBiasMetricsPass || -e $jrBam.picGcBiasMetricsInQueue || -e $jrBam.picGcBiasMetricsFail ]] ; then
					echo "### Picard GC Bias summary metric already passed, in queue, or failed for $jrBam"
				else
					echo "### Submitting for picard GC Bias Metrics: $jrBam"
					qsub -A $debit -l nodes=1:ppn=$nCores -v PICARDPATH=$picardPath,RUNDIR=$runDir,REF=$ref,BAMFILE=$jrBam,DIR=$pcDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_picardGcBiasMetrics.pbs
					if [ $? -eq 0 ] ; then
						touch $jrBam.picGcBiasMetricsInQueue
					else
						((qsubFails++))
					fi
					sleep 2
				fi
			fi
		else
			echo "### Joint IR is NOT requested for $samName"
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
