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

thisStep="pegasus_nextJob_clonalCov.txt"
nxtStep1="pegasus_nextJob_cna15.txt"

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
ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
samtoolsPath=`grep @@"$recipe"@@ $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
clonalCovPath=`grep @@"$recipe"@@ $constants | grep @@CLONALCOVPATH= | cut -d= -f2`
indels=`grep @@"$recipe"@@ $constants | grep @@INDELS= | cut -d= -f2`
doClonalCov=`grep @@"$recipe"@@ $constants | grep @@CLONALCOV= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###
# clonal coverage takes a single bam as input, so its not dependant on pairwise analysis
# but I search for it via DNAPAIR because its never used outside pairwise alignment
# Szabi made the case that individual clonal coverate will be useful, so change the code to 
# call clonal coverage, per sample, rather than member of a "DNAPAIR"
###
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
	if [ $doClonalCov == "no" ] ; then
                echo "clonal coverage not requested for this recipe"
                continue
        fi
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

	normalBamFile=$runDir/$pair1KitName/$pair1/$pair1.proj.md.jr.bam
	tumorBamFile=$runDir/$pair2KitName/$pair2/$pair2.proj.md.jr.bam

	echo "### normal BAM: $normalBamFile"
	echo "### tumor  BAM: $tumorBamFile"

	normalBaiFile=${normalBamFile/.bam/.bai}
	tumorBaiFile=${tumorBamFile/.bam/.bai}
	if [[ -e $normalBamFile.clonalCovInQueue || -e $normalBamFile.clonalCovFail || -e $normalBamFile.clonalCovPass ]] ; then
		echo "### Clonal cov already in queue, passed, or failed"
	else
		if [[ ! -e $normalBamFile || ! -e $normalBaiFile || ! -e $normalBamFile.jointIRPass ]] ; then
			echo "### Normal bam, bai, or jointIRPass does not exist"
			((qsubFails++))
		else
			echo "### Submitting to queue with $normalBamFile"
			sbatch -v BAMFILE=$normalBamFile,OUTFILE=$normalBamFile.clc,CPATH=$clonalCovPath,SAMPATH=$samtoolsPath,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pegasusPbsHome/pegasus_clonalCov.pbs
			if [ $? -eq 0 ] ; then
				touch $normalBamFile.clonalCovInQueue
			else
				((qsubFails++))
			fi
			sleep 2
		fi
	fi
	if [[ -e $tumorBamFile.clonalCovInQueue || -e $tumorBamFile.clonalCovFail || -e $tumorBamFile.clonalCovPass ]] ; then
		echo "### Clonal cov already in queue, passed, or failed"
	else
		if [[ ! -e $tumorBamFile || ! -e $tumorBaiFile || ! -e $tumorBamFile.jointIRPass ]] ; then
			echo "Tumor bam, bai, or jointIRPass does not exist"
			((qsubFails++))
		else
			echo "### Submitting to queue with $tumorBamFile"
			sbatch -v BAMFILE=$tumorBamFile,OUTFILE=$tumorBamFile.clc,CPATH=$clonalCovPath,SAMPATH=$samtoolsPath,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pegasusPbsHome/pegasus_clonalCov.pbs
			if [ $? -eq 0 ] ; then
				touch $tumorBamFile.clonalCovInQueue
			else
				((qsubFails++))
			fi
			sleep 2
		fi

	fi
done
for dnaFamiLine in `cat $configFile | grep '^DNAFAMI='`
do
	if [ $doClonalCov == "no" ] ; then
                echo "clonal coverage not requested for this recipe"
                continue
        fi
	echo "### DNAFAMI line is $dnaFamiLine"
        sampleNames=`echo $dnaFamiLine | cut -d= -f2 | cut -d\; -f1`
        shorterName=`echo $dnaFamiLine | cut -d= -f2 | cut -d\; -f2`
        if [[ $shorterName == *,* ]] ; then
                shorterName=""
                echo "### No nick name detected for these pairs"
                usableName=${sampleNames//,/-}
        else
                echo "### Nick name for these pairs are: $shorterName"
                usableName="$shorterName"
        fi

	#IFS=","
	sampleList=`echo $sampleNames | sed 's/,/ /g'`
	for ccSample in $sampleList
	do
		echo "ccSample: $ccSample"
		sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$ccSample"'"'`
		sampleKitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
		samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
		bamFile=$runDir/$sampleKitName/$samName/$samName.proj.md.jr.bam
		echo "What I have: sampleLine: $sampleLine, sampleKitName: $sampleKitName, samName: $samName"
		echo "BAM FILE: $bamFile"
		baiFile=${bamFile/.bam/.bai}
		echo "BAI FILE: $baiFile"

		if [[ -e $bamFile.clonalCovInQueue || -e $bamFile.clonalCovFail || -e $bamFile.clonalCovPass ]] ; then
			echo "### Clonal cov already in queue, passed, or failed"
		else
			if [[ ! -e $bamFile || ! -e $baiFile || ! -e $bamFile.jointIRPass ]] ; then
				echo "### bam, bai, or jointIRPass does not exist"
				((qsubFails++))
			else
				echo "### Submitting to queue with $bamFile"
				sbatch -v BAMFILE=$bamFile,OUTFILE=$bamFile.clc,CPATH=$clonalCovPath,SAMPATH=$samtoolsPath,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pegasusPbsHome/pegasus_clonalCov.pbs
				if [ $? -eq 0 ] ; then
					touch $normalBamFile.clonalCovInQueue
				else
					((qsubFails++))
				fi
				sleep 2
			fi
		fi

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
