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

thisStep="pegasus_nextJob_vqsr.txt"
nxtStep1="pegasus_nextJob_snpEff.txt"

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


ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
snps=`grep @@"$recipe"@@ $constants | grep @@SNPS= | cut -d= -f2`
hapmap=`grep @@"$recipe"@@ $constants | grep @@HAPMAP= | cut -d= -f2`
omni=`grep @@"$recipe"@@ $constants | grep @@OMNI= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0

for dnaPairLine in `cat $configFile | grep '^DNAPAIR=\|^DNAFAMI='`
do
	echo "### DNA pair line is $dnaPairLine"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
	for eachSample in ${sampleNames//,/ }
	do
		((sampleCount++))
		#echo "eachsample: $eachSample"
		sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
		kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
		samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
		assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
		libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
	done
	usableName=${sampleNames//,/-}
	sampleCount=0
	missingSampleCount=0
	sampleList=""
	ugDir="$runDir/ug"
	workDir="$ugDir/$usableName"
	trackName="$runDir/ug/$usableName/$usableName"
	if [[ ! -e $trackName.ugPass || ! -e $trackName.UG.vcf ]] ; then
		echo "ugPass or UG.vcf itself does not exist yet for: $trackName"
		continue
	fi
	if [[ -e $trackName.UG.vcf.vqsrInQueue || -e $trackName.UG.vcf.vqsrPass || -e $trackName.UG.vcf.vqsrFail ]] ; then 
		echo "### VQSR is already done, failed or inQueue"
		continue
	fi 
	echo "### Submitting $trackName.UG.vcf for VQSR..."
	sbatch -n 1 -N 1 --cpus-per-task $nCores --export HAPMAP=$hapmap,RSCRIPT=$trackName.UG.vcf.plots.R,OMNI=$omni,RECAL=$trackName.UG.vcf.recal,TRANCHES=$trackName.UG.vcf.tranches,ASSAY=$assayID,GATKPATH=$gatkPath,KNOWN=$snps,VCF=$trackName.UG.vcf,REF=$ref,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_vqsr.sh
	if [ $? -eq 0 ] ; then
		touch $trackName.UG.vcf.vqsrInQueue
	else
		((qsubFails++))
	fi
	sleep 2

#	for eachSample in ${sampleNames//,/ }
#	do
#		((sampleCount++))
#		#echo "eachsample: $eachSample"
#		sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
#		kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
#		samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
#		eachSampleVcf=$runDir/$kitName/$samName/$samName.proj.md.jr.bam.UG.vcf
#		eachSamplePass=$runDir/$kitName/$samName/$samName.proj.md.jr.bam.ugPass
#		if [[ ! -e $eachSampleVcf || ! -e $eachSamplePass ]] ; then
#			echo "### Can't find the UG.vcf or md.jr.bam.ugPass"
#			echo "### BAM: $eachSampleVcf"
#			echo "### PAS: $eachSamplePass"
#			((missingSampleCount++))
#		else
#			if [[ -e $eachSampleVcf.vqsrInQueue || -e $eachSampleVcf.vqsrPass || -e $eachSampleVcf.vqsrFail ]] ; then 
#				echo "### VQSR is already done, failed or inQueue"
#				continue
#			fi 
#			echo "### Submitting $eachSampleVcf for VQSR..."
#			sbatch -n 1 -N 1 --cpus-per-task $nCores --export HAPMAP=$hapmap,RSCRIPT=$eachSampleVcf.plots.R,OMNI=$omni,RECAL=$eachSampleVcf.recal,TRANCHES=$eachSampleVcf.tranches,ASSAY=$assayID,GATKPATH=$gatkPath,KNOWN=$snps,VCF=$eachSampleVcf,REF=$ref,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_vqsr.sh
#			if [ $? -eq 0 ] ; then
#				touch $eachSampleVcf.vqsrInQueue
#			else
#				((qsubFails++))
#			fi
#			sleep 2
#		fi
#	done
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
		mdVcf=$runDir/$kitName/$samName/$samName.proj.md.bam.UG.vcf
		jrBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
		jrRequested=`cat $configFile | grep '^DNAPAIR=\|^DNAFAMI=' | grep $samName | wc -l`
		if [ $jrRequested -gt 0 ] ; then
			echo "### Unified genotyper will not operate on a single bam because Joint IR is requested for $samName"
		else
			echo "### Unified genotyper will run on single bam because Joint IR is NOT requested for $samName"
			if [[ ! -e $mdBam.ugPass || ! -e $mdVcf ]] ; then
				echo "### Either mdPass or the bam itself is missing for $mdVcf"
				((qsubFails++))
			else
				if [[ -e $mdVcf.vqsrPass || -e $mdVcf.vqsrInQueue || -e $mdVcf.vqsrFail ]] ; then
					echo "### Unified genotyper already passed, in queue, or failed for $mdVcf"
				else
					echo "### Submitting for single bam unified genotyper: $mdVcf"
					sbatch -n 1 -N 1 --cpus-per-task $nCores --export HAPMAP=$hapmap,RSCRIPT=$mdVcf.plots.R,OMNI=$omni,RECAL=$mdVcf.recal,TRANCHES=$mdVcf.tranches,ASSAY=$assayID,GATKPATH=$gatkPath,KNOWN=$snps,VCF=$mdVcf,REF=$ref,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_vqsr.sh
					if [ $? -eq 0 ] ; then
						touch $mdVcf.vqsrInQueue
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
