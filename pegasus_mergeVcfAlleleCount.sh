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

thisStep="pegasus_nextJob_mergeVcfAlleleCount.txt "
nxtStep1="pegasus_nextJob_postmergeVcfAlleleCount.txt"
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


snpSift=`grep @@"$recipe"@@ $constants | grep @@SNPEFFPATH= | cut -d= -f2`
samTools=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
varScan=`grep @@VARSCANPATH= $constantsDir/$recipe | cut -d= -f2`
ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
refDict=${ref//fa/dict}
cosmicVcf=`grep "@@"$recipe"@@" $constants | grep @@COSMIC_VCF= | cut -d= -f2`
snps=`grep "@@"$recipe"@@" $constants | grep @@SNPS= | cut -d= -f2`
indels=`grep @@"$recipe"@@ $constants | grep @@INDELS= | cut -d= -f2`

rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
hapmap=`grep @@"$recipe"@@ $constants | grep @@HAPMAP= | cut -d= -f2`
omni=`grep @@"$recipe"@@ $constants | grep @@OMNI= | cut -d= -f2`

snpeffdb=`grep @@"$recipe"@@ $constants | grep @@SNPEFFDB= | cut -d= -f2`
dbsnp=`grep @@"$recipe"@@ $constants | grep @@SNPS= | cut -d= -f2`
snpeffPath=`grep @@"$recipe"@@ $constants | grep @@SNPEFFPATH= | cut -d= -f2`

DBNSFP=/home/tgenref/pecan/bin/vcfMerger/dbNSFP2.4.txt.gz
VCFMERGER=/home/tgenref/pecan/bin/vcfMerger_V2_norm/production.version.VCFMERGER_V2.WithNorm.201600819/pecan.merge.3vcfs.main.sh
VCFMERGER_DIR=/home/tgenref/pecan/bin/vcfMerger_V2_norm/production.version.VCFMERGER_V2.WithNorm.201600819
VCFSORTER=/home/tgenref/pecan/bin/vcfMerger/vcfsorter.pl
RNA_VCF_HEADER=/home/tgenref/pecan/bin/vcfMerger/RNA_VCF_HEADER.vcf
POST_MERGE_VENN=/home/tgenref/pecan/bin/vcfMerger_V2_norm/production.version.VCFMERGER_V2.WithNorm.201600819/pecan.Venn_postMMRF_specific_filtering.sh
KG=/home/tgenref/pecan/bin/vcfMerger/1000G_phase1.snps.high_confidence.b37.vcf
COSMIC=/home/tgenref/pecan/bin/vcfMerger/CosmicCodingMuts_v66_20130725_sorted.vcf
NHLBI=/home/tgenref/pecan/bin/vcfMerger/ESP6500SI-V2_snps_indels.vcf
DBSNP=/home/tgenref/pecan/bin/vcfMerger/dbsnp_137.b37.vcf


echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###first check all these vcfs are complete/passed
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
	echo "### DNA pair line is $dnaPairLine for seurat stuff"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
	alleleCount=`cat $configFile | grep '^TRIPLET4ALLELECOUNT=' | grep ${sampleNames} | head -1`
	if [ -z "$alleleCount" ] ; then
		echo "allele count not requested for this DNAPAIR"
	else
		rnaSample=`echo $alleleCount | cut -d, -f3`
		rnaBam=`find $runDir -name ${rnaSample}.proj.Aligned.out.sorted.md.bam | head -1`
		echo "alle count is requested for $sampleNames"
		echo "the matching rnaSample is $rnaSample"

		for eachSample in ${sampleNames//,/ }
		do
			((sampleCount++))
			#echo "eachsample: $eachSample"
			sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
			kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
			samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
			assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
		done
		usableName=${sampleNames//,/-}
		sampleCount=0
		missingSampleCount=0
		sampleList=""
		mergerDir="$runDir/vcfMerger/$usableName"
		
		if [[ -e ${mergerDir}/${usableName}.vcfMergerACPass || -e ${mergerDir}/${usableName}.vcfMergerACInQueue || -e ${mergerDir}/${usableName}.vcfMergerACFail ]] ; then
			
			echo "### This vcf merger allele count pair already passed, failed, or inQueue."
			continue
                fi
		echo "first checking for merged vcf"
		mergedVcf="$mergerDir/$usableName.merge.sort.clean.f2t.ann.vcf"
		if [ ! -e $mergerDir/$usableName.vcfMergerPass ] ; then
			echo "### Vcf Merger Pass does not exist yet $mergerDir/$usableName.vcfMergerPass"
			((qsubFails++))
			continue
		else	
			echo "### Submitting vcf to queue for vcf merger allele count..."
			sbatch --export SNPEFFPATH=$snpeffPath,SNPSIFT=$snpSift,DBNSFP=$DBNSFP,SAMTOOLS=$samTools,VARSCAN=$varScan,REF=$ref,DICT=$refDict,COSMIC=$COSMIC,KG=$KG,NHLBI=$NHLBI,SNPS=$snps,INDELS=$indels,GATK=$gatkPath,VCFMERGER=$VCFMERGER,BASENAME=$usableName,VCFMERGER_DIR=$VCFMERGER_DIR,VCFSORTER=$VCFSORTER,RNA_VCF_HEADER=$RNA_VCF_HEADER,POST_MERGE_VENN=$POST_MERGE_VENN,DBSNP=$dbsnp,DBVERSION=$snpeffdb,MERGERDIR=$mergerDir,RNABAM=$rnaBam,ASSAYID=$assayID,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_mergeVcfAlleleCount.pbs
			if [ $? -eq 0 ] ; then
				touch ${mergerDir}/${usableName}.vcfMergerACInQueue
			else
				((qsubFails++))
			fi
			sleep 2
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
