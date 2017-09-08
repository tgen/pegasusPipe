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

thisStep="pegasus_nextJob_snpSniff.txt"
nxtStep1="pegasus_nextJob_checkProjectComplete.txt"

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

snpSnifferPath=`grep @@SNPSNIFFERPATH= $constantsDir/$recipe | cut -d= -f2`
snpSniff=`grep @@SNPSNIFF= $constantsDir/$recipe | cut -d= -f2`


ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
hapmap=`grep @@"$recipe"@@ $constants | grep @@HAPMAP= | cut -d= -f2`
omni=`grep @@"$recipe"@@ $constants | grep @@OMNI= | cut -d= -f2`

snpeffdb=`grep @@"$recipe"@@ $constants | grep @@SNPEFFDB= | cut -d= -f2`
dbsnp=`grep @@"$recipe"@@ $constants | grep @@SNPS= | cut -d= -f2`
snpeffPath=`grep @@"$recipe"@@ $constants | grep @@SNPEFFPATH= | cut -d= -f2`
samtoolsPath=`grep @@"$recipe"@@ $constants | grep @@SAMTOOLSSNPSNIFFPATH= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
	if [ $snpSniff != "yes" ] ; then
		echo "snpSniff not requested for this recipe"
		continue
	fi
	echo "sample is $sampleLine"
	kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
	samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
	assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
	libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
	echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID"
	if [[ "$assayID" == "Exome" || "$assayID" == "Genome" ]] ; then
		echo "### Assay ID is $assayID. Must be genome or exome."
		pcDir=$runDir/$kitName/$samName
		mdBam=$runDir/$kitName/$samName/$samName.proj.md.bam
		jrBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
		jrPas=$runDir/$kitName/$samName/$samName.proj.md.jr.bam.jointIRPass
		mdPas=$runDir/$kitName/$samName/$samName.proj.bam.mdPass
		finalOut=$runDir/$kitName/$samName/$samName.snpSniffer.vcf
		jrRequested=`cat $configFile | grep '^DNAPAIR=\|^DNAFAMI=' | grep $samName | wc -l`
		if [ $jrRequested -gt 0 ] ; then
			echo "### Will call snpSniffer on the jr bam because JR was requested for:  $samName"
			pasFile=$jrPas
			bamFile=$jrBam
		else
			echo "### Will call snpSniffer on the md bam because JR was requested for:  $samName"
			pasFile=$mdPas
			bamFile=$mdBam

		fi
		if [ ! -e $pasFile ] ; then
			echo "### Looks like bam is not ready for snpSniffer. Missing file $pasFile"
			((qsubFails++))
		else
			if [[ -e $bamFile.snpSniffInQueue || -e $bamFile.snpSniffFail || -e $bamFile.snpSniffPass ]] ; then
				echo "### Looks like snp sniff is alread in queue, failed, or passed"
			else
				echo "### Submitting for snpSniffer: $bamFile"
				sbatch -n 1 -N 1 --cpus-per-task $nCores --export SAMTOOLSPATH=$samtoolsPath,OUTVCF=$finalOut,REF=$ref,SNPSNIFFERPATH=$snpSnifferPath,BAM=$bamFile,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pegasusPbsHome/pegasus_snpSniffer.pbs
				if [ $? -eq 0 ] ; then
					touch $bamFile.snpSniffInQueue
				else
					((qsubFails++))
				fi
				sleep 2
			fi
		fi
	elif [ "$assayID" == "RNA" ] ; then
		echo "### Assay ID is $assayID. Must be RNA."
		case $rnaAligner in 
		tophat) echo "tophat case"
			pasFile="$runDir/$kitName/$samName/$samName.topHatDir/$samName.proj.accepted_hits.bam.rnaMarkDupPass"
			bamFile="$runDir/$kitName/$samName/$samName.topHatDir/$samName.proj.accepted_hits.bam"
			;;
		star) echo "star case"
			pasFile="$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sorted.bam.rnaMarkDupPass"
			bamFile="$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sorted.md.bam"
			;;
		anotherRNAaligner) echo "example RNA aligner"
			;;
		*) echo "I should not be here"
			;;
		esac 
		if [ ! -e $pasFile ] ; then 
			echo "### Looks like rna mark dup is not done yet. Missing $pasFile"
			((qsubFails++))
		else
			if [[ -e $bamFile.snpSniffInQueue || -e $bamFile.snpSniffFail || -e $bamFile.snpSniffPass ]] ; then
				echo "### Looks like snp sniff is alread in queue, failed, or passed"
			else
				finalOut=$runDir/$kitName/$samName/$samName.starDir/$samName.snpSniffer.vcf
				echo "### Submitting for snpSniffer: $bamFile"
				sbatch -n 1 -N 1 --cpus-per-task $nCores --export SAMTOOLSPATH=$samtoolsPath,REF=$ref,OUTVCF=$finalOut,SNPSNIFFERPATH=$snpSnifferPath,BAM=$bamFile,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pegasusPbsHome/pegasus_snpSniffer.pbs
				if [ $? -eq 0 ] ; then
					touch $bamFile.snpSniffInQueue
				else
					((qsubFails++))
				fi
			fi
		fi
	else
		echo "### Assay ID is $assayID"
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
