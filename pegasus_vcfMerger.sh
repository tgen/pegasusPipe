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

thisStep="pegasus_nextJob_vcfMerger.txt"
nxtStep1="pegasus_nextJob_mergeVcfAlleleCount.txt"

constants=${JETSTREAM_HOME}/centralPipe/constants/constants.txt
constantsDir=${JETSTREAM_HOME}/centralPipe/constants/
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
matchedNormal=`cat $configFile | grep "^MATCHEDNORMAL=" | cut -d= -f2 | head -1 | tr -d [:space:]`

nCores=`grep @@${myName}_CORES= $constantsDir/$recipe | cut -d= -f2`


snpSift=`grep @@"$recipe"@@ $constants | grep @@SNPEFFPATH= | cut -d= -f2`
samTools=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
varScan=`grep @@VARSCANPATH= $constantsDir/$recipe | cut -d= -f2`
ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
refDict="${ref%.fa}.dict"
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
picardPath=`grep @@"$recipe"@@ $constants | grep @@PICARDPATH= | cut -d= -f2`

DBNSFP="/home/tgenref/homo_sapiens/grch37_hg19/public_databases/dbnsfp/v2.4/dbNSFP2.4.txt.gz"
VCFMERGER_DIR="/home/tgenref/binaries/vcfMerger/mergeVcf/normalization/vcfMerger/"
VCFMERGER="${VCFMERGER_DIR}/pecan.merge.3vcfs.main.sh"
VCFSORTER="/home/tgenref/binaries/vcfMerger/vcfMerger/vcfsorter.pl"
RNA_VCF_HEADER="/home/tgenref/binaries/vcfMerger/vcfMerger/RNA_VCF_HEADER.vcf"
COSMIC="home/tgenref/binaries/vcfMerger/vcfMerger/CosmicCodingMuts_v66_20130725_sorted.vcf"
KG=`grep "@@"$recipe"@@" $constants | grep @@KG= | cut -d= -f2`
NHLBI=`grep "@@"$recipe"@@" $constants | grep @@NHLBI= | cut -d= -f2`
COSMICC=`grep "@@"$recipe"@@" $constants | grep @@COSMICC= | cut -d= -f2`
COSMICNC=`grep "@@"$recipe"@@" $constants | grep @@COSMICNC= | cut -d= -f2`
EXAC=`grep "@@"$recipe"@@" $constants | grep @@EXAC= | cut -d= -f2`

POST_MERGE_VENN=${VCFMERGER_DIR}/pecan.Venn_postMMRF_specific_filtering.sh

##used if there is no matched normal
DBSNP_DIV_BED="/home/tgenref/homo_sapiens/grch37_hg19/hs37d5/tool_specific_resources/tconut/nonmatched_tumor_normal/Merged_Indels_no_COSMIC_Y.bed"
DBSNP_SNV_BED="/home/tgenref/homo_sapiens/grch37_hg19/hs37d5/tool_specific_resources/tconut/nonmatched_tumor_normal/Merged_SNPs_no_COSMIC_Y.bed"


echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###first check all these vcfs are complete/passed
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
    rnaBam=""
    echo "### DNA pair line is $dnaPairLine for seurat stuff"
    sampleNames=`echo $dnaPairLine | cut -d= -f2`
    alleleCount=`cat $configFile | grep '^TRIPLET4ALLELECOUNT=' | grep ${sampleNames} | head -1`
    if [ -z "$alleleCount" ] ; then
        echo "allele count not requested for this DNAPAIR"
    else
        rnaSample=`echo $alleleCount | cut -d, -f3`
        echo "rnaSample: $rnaSample"
        rnaBam=`find $runDir -name "${rnaSample}.proj.Aligned.out.sorted.md.bam" | head -1`
        echo "rnaBam: $rnaBam"
        if [ -z "$rnaBam" ] ; then
            echo "Couldn't find the RNABAM for $rnaSample"
            continue
        fi
        echo "allele count is requested for $sampleNames"
        echo "the matchng rnaSample is $rnaSample"
    fi

    for eachSample in ${sampleNames//,/ }
    do
        ((sampleCount++))
        #echo "eachsample: $eachSample"
        sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
        kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
        samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
        assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
    done
    control=`echo $dnaPairLine | cut -d= -f2 | cut -d, -f1`
    tumor=`echo $dnaPairLine | cut -d, -f2`
    echo "control = $control tumor = $tumor"
    usableName=${sampleNames//,/-}
    sampleCount=0
    missingSampleCount=0
    sampleList=""

    echo "$kitName"
    if [[ "$kitName" == "TSE61" ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/illumina_truseq/TruSeq_exome_targeted_regions_b37_padded.bed"
    elif [[ "$kitName" == *S5U ]] || [[ "$kitName" == *S5X ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/hs37d5/gene_model/ensembl_v70/tool_specific_resources/vcfmerger/Agilent_SureSelect_V5_plusUTR/Ensembl_v70_hs37d5_exonic_coordinates_touched_v5UTR_padded25.bed"
    elif [[ "$kitName" == *STX ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v5_strexome/Strexome_targets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *SCR ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_cre_v1/Agilent_Clinical_Research_Exome_hs37d5_PaddedTargets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *S2X ]] ; then
        bedFile="/home/tgenref/pipeline_v0.4/annotations/exome_capture/Agilent_V2_hs37d5/Agilent_V2_hs37d5_Targets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *STL ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_custom_strexome_lite/Strexome_Lite_Targets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *S1X ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v1_NA/SureSelectV1_hs37d5_PaddedTargets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *S6X ]] ; then
        bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_v6_noUTR/Agilent_V6_noUTR_hs37d5_Targets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *SXP ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v6r2_plusUTR/Agilent_SureSelect_V6R2_plusUTR_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *S4X ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v4_noUTR/Agilent_V4_noUTR_hs37d5_Targets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *E62 ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/illumina_nextera_expanded/NexteraExpandedExome_hs37d5_Targets_PicardPadded100.bed"
    elif [[ "$kitName" == *SC2 ]] ; then
        bedFile="/home/tgenref/canis_familiaris/canfam3/capture_targets/canine_exonV2/agilent_canine_exonV2_targets_padded100.bed"
    elif [[ "$kitName" == *UBC ]] ; then
        bedFile="/home/tgenref/canis_familiaris/canfam3/capture_targets/Uppsala_Broad_EZ_HX1/120705_CF3_Uppsala_Broad_EZ_HX1_padded100.bed"
    elif [[ "$kitName" == *S6A ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v6r2_strexomeAD/StrAD_targets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *S4U ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v4_plusUTR/Agilent_SureSelect_V4_plusUTR_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *ST2 ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v6r2_strexomeV2/Agilent_SureSelect_V6R2_StrexomeV2_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *CCC ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_clearseq_compCancer_NA/Agilent_ClearSeq_Beta_ComprehensiveCancer_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *STP ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v6r2_strexomeProstateV2/Agilent_SureSelect_V6R2_StxProstate_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *S2U ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v2_NA/Agilent_SureSelect_V2_NA_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *V6C ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v6r2_plusCOSMIC/Agilent_SureSelect_V6R2_plusCOSMIC_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *CR2 ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_cre_v2/Agilent_SureSelect_CREv2_cliRes_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
    elif [[ "$kitName" == *R12 ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/illumina_nextera_v1.2/padded_targets.bed"
    elif [[ "$kitName" == *R37 ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/illumina_nextera_v1.2/NexteraRapidCaptureExomeV1.2_Targets.txt"
    elif [[ "$kitName" == *KBM ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/Baylor_MS_Ig_CS_SNP3/Baylor_MS_Ig_CS_SNP3_Regions.preprocessed.filter.bed"
    elif [[ "$kitName" == *S6B ]] ; then
        bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_human_all_exon_v6s_xt_beta/Human_All_Exon_v6S_XT_beta_targets.preprocessed.filter.bed"
    fi

    echo "### BED FILE= $bedFile"

    echo "first checking for seurat snpeff vcf"
    seuratTrackName="$runDir/seurat/$usableName/$usableName"
    if [ ! -e $seuratTrackName.seurat.vcf.snpEffPass ] ; then
        echo "### Seurat snpEffPass doesnt exist yet: $seuratTrackName.seurat.vcf.snpEffPass"
        ((qsubFails++))
        exit 1
    fi
    echo "checking for strelka snpeff vcfs"
    strelkaTrackName="$runDir/strelka/$usableName"
    vcfPre="$runDir/strelka/$usableName/myAnalysis/results/$usableName"
    if [[ ! -e $vcfPre.strelka.all.somatic.snvs.vcf.snpEffPass || ! -e $vcfPre.strelka.passed.somatic.snvs.vcf.snpEffPass || ! -e $vcfPre.strelka.passed.somatic.indels.vcf.snpEffPass || ! -e $vcfPre.strelka.all.somatic.indels.vcf.snpEffPass ]] ; then
        echo "### strelka snpEff doesn't exist for one of the 4 strelka vcfs"
        ((qsubFails++))
        exit 1
    fi
    echo "now checking for mutect vcfs"
    mutectTrackName="$runDir/mutect/$usableName/$usableName"
    vcf="${mutectTrackName}_MuTect_All.vcf"
    if [ ! -e $vcf.snpEffPass  ] ; then
        echo "### mutect snpEff pass doesnt exist yet: $mutectTrackName.mutectPass"
        ((qsubFails++))
        continue
    else

        mergerDir="$runDir/vcfMerger/$usableName"
        mkdir -p $mergerDir
        seuratVcf="$seuratTrackName.seurat.vcf"
        mutectVcf="${mutectTrackName}_MuTect_All.vcf"
        strelkaIndelVcf="$vcfPre.strelka.passed.somatic.indels.vcf"
        strelkaSnvVcf="$vcfPre.strelka.passed.somatic.snvs.vcf"
        seurat_basename=`basename ${seuratVcf} ".seurat.vcf"`

        if [[ -e ${mergerDir}/${seurat_basename}.vcfMergerPass || -e ${mergerDir}/${seurat_basename}.vcfMergerInQueue || -e ${mergerDir}/${seurat_basename}.vcfMergerFail ]] ; then
            echo "### This vcf merger pair already passed, failed, or inQueue."
            continue
        fi
                echo "### Submitting vcfs to queue for vcf merger..."
                echo "sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out --export ALL,RECIPE=$recipe,PICARDPATH=$picardPath,EXAC=$EXAC,COSMICC=$COSMICC,COSMICNC=$COSMICNC,DBSNP_SNV_BED=$DBSNP_SNV_BED,DBSNP_DIV_BED=$DBSNP_DIV_BED,MATCHEDNORMAL=$matchedNormal,SNPEFFPATH=$snpeffPath,CONTROL=$control,TUMOR=$tumor,SNPSIFT=$snpSift,DBNSFP=$DBNSFP,DBNSP=$DBNSP,SAMTOOLS=$samTools,VARSCAN=$varScan,REF=$ref,DICT=$refDict,COSMIC=$COSMIC,KG=$KG,NHLBI=$NHLBI,SNPS=$snps,INDELS=$indels,GATK=$gatkPath,VCFMERGER=$VCFMERGER,VCFMERGER_DIR=$VCFMERGER_DIR,VCFSORTER=$VCFSORTER,RNA_VCF_HEADER=$RNA_VCF_HEADER,POST_MERGE_VENN=$POST_MERGE_VENN,DBSNP=$dbsnp,DBVERSION=$snpeffdb,SEURAT_VCF=$seuratVcf,MUTECT_VCF=$mutectVcf,STRELKA_SNV_VCF=$strelkaSnvVcf,STRELKA_INDEL_VCF=$strelkaIndelVcf,MERGERDIR=$mergerDir,RNABAM=$rnaBam,ASSAYID=$assayID,BEDFILE=$bedFile,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_vcfMerger.sh"
                sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out --export ALL,RECIPE=$recipe,PICARDPATH=$picardPath,EXAC=$EXAC,COSMICC=$COSMICC,COSMICNC=$COSMICNC,DBSNP_SNV_BED=$DBSNP_SNV_BED,DBSNP_DIV_BED=$DBSNP_DIV_BED,MATCHEDNORMAL=$matchedNormal,SNPEFFPATH=$snpeffPath,CONTROL=$control,TUMOR=$tumor,SNPSIFT=$snpSift,DBNSFP=$DBNSFP,DBNSP=$DBNSP,SAMTOOLS=$samTools,VARSCAN=$varScan,REF=$ref,DICT=$refDict,COSMIC=$COSMIC,KG=$KG,NHLBI=$NHLBI,SNPS=$snps,INDELS=$indels,GATK=$gatkPath,VCFMERGER=$VCFMERGER,VCFMERGER_DIR=$VCFMERGER_DIR,VCFSORTER=$VCFSORTER,RNA_VCF_HEADER=$RNA_VCF_HEADER,POST_MERGE_VENN=$POST_MERGE_VENN,DBSNP=$dbsnp,DBVERSION=$snpeffdb,SEURAT_VCF=$seuratVcf,MUTECT_VCF=$mutectVcf,STRELKA_SNV_VCF=$strelkaSnvVcf,STRELKA_INDEL_VCF=$strelkaIndelVcf,MERGERDIR=$mergerDir,RNABAM=$rnaBam,ASSAYID=$assayID,BEDFILE=$bedFile,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_vcfMerger.sh
                if [ $? -eq 0 ] ; then
                        touch ${mergerDir}/${seurat_basename}.vcfMergerInQueue
                else
                        ((qsubFails++))
                fi
                sleep 2
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
