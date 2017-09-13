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

thisStep="pegasus_nextJob_germVcfMerger.txt"
#nxtStep1="pegasus_nextJob_mergeVcfAlleleCount.txt"

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
vtPath=`grep @@"$recipe"@@ $constants | grep @@VTPATH= | cut -d= -f2`

DBNSFP=/home/tgenref/pecan/bin/vcfMerger/dbNSFP2.4.txt.gz
VCFMERGER=/home/tgenref/pecan/bin/vcfMerger/pecan.merge.3vcfs.main.sh
VCFMERGER_DIR=/home/tgenref/pecan/bin/vcfMerger
VCFSORTER=/home/tgenref/pecan/bin/vcfMerger/vcfsorter.pl
RNA_VCF_HEADER=/home/tgenref/pecan/bin/vcfMerger/RNA_VCF_HEADER.vcf
POST_MERGE_VENN=/home/tgenref/pecan/bin/vcfMerger/pecan.Venn_postMMRF_specific_filtering.sh
#DBSNP=/home/tgenref/pecan/bin/vcfMerger/dbsnp_137.b37.vcf

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###first check all these vcfs are complete/passed
for dnaPairLine in `cat $configFile | grep '^DNAPAIR=\|^DNAFAMI='`
do
    #rnaBam=""
    echo "### DNA pair line is $dnaPairLine for seurat stuff"
    sampleNames=`echo $dnaPairLine | cut -d= -f2`

        shorterName=`echo $dnaPairLine | cut -d= -f2 | cut -d\; -f2`
        if [[ $shorterName == *,* ]] ; then
                shorterName=""
                echo "### No nick name detected for these pairs"
                usableName=${sampleNames//,/-}
        else
                echo "### Nick name for these pairs are: $shorterName"
                usableName="$shorterName"
        fi
        sampleCount=0
        missingSampleCount=0
        sampleList=""
    echo "$kitName"
    if [[ "$kitName" == "TSE61" ]] ; then
        bedFile="/home/tgenref/pipeline_v0.3/annotations/exome_capture/illumina_truseq/TruSeq_exome_targeted_regions_b37_padded.bed"
    elif [[ "$kitName" == *S5U ]] || [[ "$kitName" == *S5X ]] ; then
        bedFile="/home/tgenref/pipeline_v0.3/ensembl70/Ensembl_v70_hs37d5_exonic_coordinates_touched_v5UTR_padded25.bed"
    elif [[ "$kitName" == *STX ]] ; then
        bedFile="/home/tgenref/pipeline_v0.4/annotations/exome_capture/strexome/Strexome_targets_sortedTabs2_padded25.bed"
    elif [[ "$kitName" == *SCR ]] ; then
                bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_clinical_research_exome/Agilent_Clinical_Research_Exome_hs37d5_TargetsPadded25_sortedTabs2_Picard.txt"
    fi
    #bedFileGrep=$kitName"_CNABED"
        #bedFile=`grep "@@"$recipe"@@" $constants | grep @@"$bedFileGrep"= | cut -d= -f2`
        echo "### BED FILE= $bedFile"

    echo "first checking for haplotype caller snpeff vcf"
    hcTrackName="$runDir/hc/$usableName/$usableName"
    if [ ! -e $hcTrackName.HC_All.vcf.snpEffPass ] ; then
        echo "### HC snpEffPass doesnt exist yet: $hcTrackName.HC_All.vcf.snpEffPass"
        ((qsubFails++))
        exit
    fi
    echo "checking for mpileup snpeff vcfs"
    mpTrackName="$runDir/mpileup/$usableName/$usableName"
    if [ ! -e $mpTrackName.mpileup_All.vcf.snpEffPass ] ; then
                echo "### mpileup snpEffPass doesnt exist yet: $mpTrackName.mpileup_All.vcf.snpEffPass"
                ((qsubFails++))
                exit
        fi
    echo "now checking for freebayes vcfs"
    fbTrackName="$runDir/freebayes/$usableName/$usableName"
    if [ ! -e $fbTrackName.freebayes_All.vcf.snpEffPass  ] ; then
        echo "### freebayes snpEff pass doesnt exist yet: $fbTrackName.freebayes_All.vcf.snpEffPass"
        ((qsubFails++))
        continue
    else

        mergerDir="$runDir/germlineVcfMerger/$usableName"
        trackName="$runDir/germlineVcfMerger/$usableName/$usableName"
        mkdir -p $mergerDir
        hcVcf=$hcTrackName.HC_All.vcf
        mpVcf=$mpTrackName.mpileup_All.vcf
        fbVcf=$fbTrackName.freebayes_All.vcf

        if [[ -e ${mergerDir}/${usableName}.vcfMergerPass || -e ${mergerDir}/${usableName}.vcfMergerInQueue || -e ${mergerDir}/${usableName}.vcfMergerFail ]] ; then
            echo "### This vcf merger pair already passed, failed, or inQueue."
            continue
        fi
                echo "### Submitting vcfs to queue for vcf merger..."

                echo "qsub -A $debit -l nodes=1:ppn=8 -v VTPATH=$vtPath,SNPEFFPATH=$snpeffPath,SNPSIFT=$snpSift,DBNSFP=$DBNSFP,DBNSP=$DBNSP,TRACKNAME=$trackName,SAMTOOLS=$samTools,VARSCAN=$varScan,REF=$ref,DICT=$refDict,COSMIC=$COSMIC,KG=$KG,NHLBI=$NHLBI,SNPS=$snps,INDELS=$indels,GATKPATH=$gatkPath,VCFMERGER=$VCFMERGER,VCFMERGER_DIR=$VCFMERGER_DIR,VCFSORTER=$VCFSORTER,RNA_VCF_HEADER=$RNA_VCF_HEADER,POST_MERGE_VENN=$POST_MERGE_VENN,DBSNP=$dbsnp,DBVERSION=$snpeffdb,HCVCF=$hcVcf,MPVCF=$mpVcf,FBVCF=$fbVcf,MERGERDIR=$mergerDir,ASSAYID=$assayID,BEDFILE=$bedFile,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pegasusPbsHome/pegasus_germVcfMerge.sh"
                sbatch --output $runDir/oeFiles/%x-slurm-%j.out --export VTPATH=$vtPath,SNPEFFPATH=$snpeffPath,SNPSIFT=$snpSift,DBNSFP=$DBNSFP,DBNSP=$DBNSP,TRACKNAME=$trackName,SAMTOOLS=$samTools,VARSCAN=$varScan,REF=$ref,DICT=$refDict,COSMIC=$COSMIC,KG=$KG,NHLBI=$NHLBI,SNPS=$snps,INDELS=$indels,GATKPATH=$gatkPath,VCFMERGER=$VCFMERGER,VCFMERGER_DIR=$VCFMERGER_DIR,VCFSORTER=$VCFSORTER,RNA_VCF_HEADER=$RNA_VCF_HEADER,POST_MERGE_VENN=$POST_MERGE_VENN,DBSNP=$dbsnp,DBVERSION=$snpeffdb,HCVCF=$hcVcf,MPVCF=$mpVcf,FBVCF=$fbVcf,MERGERDIR=$mergerDir,ASSAYID=$assayID,BEDFILE=$bedFile,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pegasusPbsHome/pegasus_germVcfMerge.sh
                if [ $? -eq 0 ] ; then
                        touch ${mergerDir}/${usableName}.vcfMergerInQueue
                else
                        ((qsubFails++))
                fi
                sleep 2
        fi
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
               # pcDir=$runDir/$kitName/$samName
                #inBam=$runDir/$kitName/$samName/$samName.proj.bam
                #mdBam=$runDir/$kitName/$samName/$samName.proj.md.bam
                #jrBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
                jrRequested=`cat $configFile | grep '^DNAPAIR=\|^DNAFAMI=' | grep $samName | wc -l`
                if [ $jrRequested -gt 0 ] ; then
                        echo "### Germline vcf merger will not operate on a single sample because Joint IR is requested for $samName"
                else
                        echo "### Germline vcf merger will run on single sample because Joint IR is NOT requested for $samName"
            echo "Kitname: $kitName"
            if [[ "$kitName" == "TSE61" ]] ; then
                bedFile="/home/tgenref/pipeline_v0.3/annotations/exome_capture/illumina_truseq/TruSeq_exome_targeted_regions_b37_padded.bed"
            elif [[ "$kitName" == *S5U ]] || [[ "$kitName" == *S5X ]] ; then
                bedFile="/home/tgenref/pipeline_v0.3/ensembl70/Ensembl_v70_hs37d5_exonic_coordinates_touched_v5UTR_padded25.bed"
            elif [[ "$kitName" == *STX ]] ; then
                bedFile="/home/tgenref/pipeline_v0.4/annotations/exome_capture/strexome/Strexome_targets_sortedTabs2_padded25.bed"
            elif [[ "$kitName" == *SCR ]] ; then
                bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_clinical_research_exome/Agilent_Clinical_Research_Exome_hs37d5_TargetsPadded25_sortedTabs2_Picard.txt"
            fi
            echo "### BED FILE= $bedFile"

            hcVcf=$runDir/$kitName/$samName/$samName.proj.md.bam.HC_All.vcf
            fbVcf=$runDir/freebayes/$samName/$samName.freebayes_All.vcf
            mpVcf=$runDir/mpileup/$samName/$samName.mpileup_All.vcf
            echo "first checking for haplotype caller snpeff vcf"
            if [ ! -e $hcVcf.snpEffPass ] ; then
                echo "### HC snpEffPass doesnt exist yet: $hcVcf.snpEffPass"
                ((qsubFails++))
                exit
            fi
            echo "checking for mpileup snpeff vcfs"
            if [ ! -e $mpVcf.snpEffPass ] ; then
                echo "### mpileup snpEffPass doesnt exist yet: $mpVcf.snpEffPass"
                ((qsubFails++))
                exit
            fi
            echo "now checking for freebayes vcfs"
            if [ ! -e $fbVcf.snpEffPass  ] ; then
                echo "### freebayes snpEff pass doesnt exist yet: $fbVcf.snpEffPass"
                ((qsubFails++))
                continue
            else

                mergerDir="$runDir/germlineVcfMerger/$samName"
                trackName="$runDir/germlineVcfMerger/$samName/$samName"
                mkdir -p $mergerDir

                if [[ -e ${mergerDir}/${samName}.vcfMergerPass || -e ${mergerDir}/${samName}.vcfMergerInQueue || -e ${mergerDir}/${samName}.vcfMergerFail ]] ; then
                    echo "### This vcf merger already passed, failed, or inQueue."
                    continue
                fi
                echo "### Submitting vcfs to queue for vcf merger..."

                echo "qsub -A $debit -l nodes=1:ppn=8 -v VTPATH=$vtPath,SNPEFFPATH=$snpeffPath,SNPSIFT=$snpSift,DBNSFP=$DBNSFP,DBNSP=$DBNSP,TRACKNAME=$trackName,SAMTOOLS=$samTools,VARSCAN=$varScan,REF=$ref,DICT=$refDict,COSMIC=$COSMIC,KG=$KG,NHLBI=$NHLBI,SNPS=$snps,INDELS=$indels,GATKPATH=$gatkPath,VCFMERGER=$VCFMERGER,VCFMERGER_DIR=$VCFMERGER_DIR,VCFSORTER=$VCFSORTER,RNA_VCF_HEADER=$RNA_VCF_HEADER,POST_MERGE_VENN=$POST_MERGE_VENN,DBSNP=$dbsnp,DBVERSION=$snpeffdb,HCVCF=$hcVcf,MPVCF=$mpVcf,FBVCF=$fbVcf,MERGERDIR=$mergerDir,ASSAYID=$assayID,BEDFILE=$bedFile,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pegasusPbsHome/pegasus_germVcfMerge.sh"
                sbatch --output $runDir/oeFiles/%x-slurm-%j.out --export VTPATH=$vtPath,SNPEFFPATH=$snpeffPath,SNPSIFT=$snpSift,DBNSFP=$DBNSFP,DBNSP=$DBNSP,TRACKNAME=$trackName,SAMTOOLS=$samTools,VARSCAN=$varScan,REF=$ref,DICT=$refDict,COSMIC=$COSMIC,KG=$KG,NHLBI=$NHLBI,SNPS=$snps,INDELS=$indels,GATKPATH=$gatkPath,VCFMERGER=$VCFMERGER,VCFMERGER_DIR=$VCFMERGER_DIR,VCFSORTER=$VCFSORTER,RNA_VCF_HEADER=$RNA_VCF_HEADER,POST_MERGE_VENN=$POST_MERGE_VENN,DBSNP=$dbsnp,DBVERSION=$snpeffdb,HCVCF=$hcVcf,MPVCF=$mpVcf,FBVCF=$fbVcf,MERGERDIR=$mergerDir,ASSAYID=$assayID,BEDFILE=$bedFile,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pegasusPbsHome/pegasus_germVcfMerge.sh
                if [ $? -eq 0 ] ; then
                    touch ${mergerDir}/${samName}.vcfMergerInQueue
                else
                    ((qsubFails++))
                fi
                sleep 2
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
