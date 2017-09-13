#!/usr/bin/env bash
#####################################################################
# Copyright (c) 2011 by The Translational Genomics Research
# Institute Ahmet Kurdoglu. All rights reserved. This License is limited 
# to, and you may use the Software solely for, your own internal and i
# non-commercial use for academic and research purposes. Without limiting 
# the foregoing, you may not use the Software as part of, or in any way 
# in connection with the production, marketing, sale or support of any 
# commercial product or service or for any governmental purposes. For 
# commercial or governmental use, please contact dcraig@tgen.org. By 
# installing this Software you are agreeing to the terms of the LICENSE 
# file distributed with this software.
#####################################################################
# source ~/.bashrc # Why?
time=`date +%d-%m-%Y-%H-%M`
echo "Starting $0 at $time"

# scriptsHome="/home/tgenjetstream/pegasus-pipe" # scriptsHome should always be
# relative to the location of this main script.
# There is no reliable way of resolving the location of a script in bash, see
# http://mywiki.wooledge.org/BashFAQ/028
# https://stackoverflow.com/a/246128/3924113
# but, this should work most of the time
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
scriptsHome="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

export pegasusPbsHome=$scriptsHome/jobScripts/

logs=~/jetstream/pegasusPipe/logs/

topProjDir=${1? Projects directory is required}
myhostname=`hostname`

echo "### ~~Running on $myhostname~~"

findCount=`ps -e | awk '$4=="find"' | wc -l`
if [ $findCount -ge 3 ] ; then
    echo "Too many finds on $myhostname ($findCount) already, quitting for $myhostname!!!"
    exit 1
fi

for messageFile in `find $topProjDir -maxdepth 2 -name [p-P]egasus_nextJob_*txt`
do
    projDir=`dirname $messageFile`
    msgName=`basename $messageFile`
    echo "### Message file: $msgName"
    case $msgName in
    Pegasus_nextJob_copyFastqs.txt)    echo "### Will copy fastqs for $projDir"
        nohup $scriptsHome/pegasus_copyFastqs.sh $projDir >> $projDir/logs/pegasus_copyFastqsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_ancestry.txt)    echo "### Will run ancestry for $projDir"
        nohup $scriptsHome/pegasus_ancestry.sh $projDir >> $projDir/logs/pegasus_ancestryLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_runFastQC.txt)    echo "### Will run fastQC for $projDir"
        nohup $scriptsHome/pegasus_runFastQC.sh $projDir >> $projDir/logs/pegasus_runFastQCLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_sexRelCheck.txt)    echo "### Will run gender/relatedness check for $projDir"
        nohup $scriptsHome/pegasus_sexRelCheck.sh $projDir >> $projDir/logs/pegasus_sexRelCheckLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_splitFastqs.txt)    echo "### Will split fastqs for $projDir"
        nohup $scriptsHome/pegasus_splitFastqs.sh $projDir >> $projDir/logs/pegasus_splitFastqsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_dnaAlign.txt)    echo "### Will run dnaAlign for $projDir"
        nohup $scriptsHome/pegasus_dnaAlign.sh $projDir >> $projDir/logs/pegasus_dnaAlignLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_digar.txt)    echo "### Will run digar for $projDir"
        nohup $scriptsHome/pegasus_digar.sh $projDir >> $projDir/logs/pegasus_digarLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_digarPost.txt)    echo "### Will run digarPost for $projDir"
        nohup $scriptsHome/pegasus_digarPost.sh $projDir >> $projDir/logs/pegasus_digarPostLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_dnaAlignParts.txt)    echo "### Will run dnaAlign in parts for $projDir"
        nohup $scriptsHome/pegasus_dnaAlignParts.sh $projDir >> $projDir/logs/pegasus_dnaAlignPartsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_indelRealign.txt)    echo "### Will run indel realign for $projDir"
        nohup $scriptsHome/pegasus_indelRealign.sh $projDir >> $projDir/logs/pegasus_indelRealignLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_indelRealignParts.txt)    echo "### Will run indel realign in parts for $projDir"
        nohup $scriptsHome/pegasus_indelRealignParts.sh $projDir >> $projDir/logs/pegasus_indelRealignPartsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_recalibrate.txt)    echo "### Will run recalibrate for $projDir"
        nohup $scriptsHome/pegasus_recalibrate.sh $projDir >> $projDir/logs/pegasus_recalibrateLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_RNArecalibrate.txt)    echo "### Will run RNArecalibrate for $projDir"
        nohup $scriptsHome/pegasus_RNArecalibrate.sh $projDir >> $projDir/logs/pegasus_RNArecalibrateLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_recalibrateParts.txt)    echo "### Will run recalibrate in parts for $projDir"
        nohup $scriptsHome/pegasus_recalibrateParts.sh $projDir >> $projDir/logs/pegasus_recalibratePartsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_mergeBams.txt)    echo "### Will merge bams $projDir"
        nohup $scriptsHome/pegasus_mergeBams.sh $projDir >> $projDir/logs/pegasus_mergeBamsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_mergeMiniBams.txt)    echo "### Will merge bams $projDir"
        nohup $scriptsHome/pegasus_mergeMiniBams.sh $projDir >> $projDir/logs/pegasus_mergeMiniBamsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_markDups.txt)    echo "### Will run mark dups for $projDir"
        nohup $scriptsHome/pegasus_markDups.sh $projDir >> $projDir/logs/pegasus_markDupsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_jointIR.txt)    echo "### Will run joint IR for $projDir"
        nohup $scriptsHome/pegasus_jointIR.sh $projDir >> $projDir/logs/pegasus_jointIRLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_mergeJointIRBams.txt)    echo "### Will run merging for joint IR bams for $projDir"
        nohup $scriptsHome/pegasus_mergeJointIRBams.sh $projDir >> $projDir/logs/pegasus_mergeJointIRBamsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_seurat.txt)    echo "### Will run seurat $projDir"
        nohup $scriptsHome/pegasus_seurat.sh $projDir >> $projDir/logs/pegasus_seuratLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_REVseurat.txt)    echo "### Will run REVseurat $projDir"
        nohup $scriptsHome/pegasus_REVseurat.sh $projDir >> $projDir/logs/pegasus_REVseuratLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_samtoolsMpileUp.txt)    echo "### Will run samtools mPileUp $projDir"
        nohup $scriptsHome/pegasus_samtoolsMpileUp.sh $projDir >> $projDir/logs/pegasus_samtoolsMpileUpLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_strelka.txt)    echo "### Will run strelka $projDir"
        nohup $scriptsHome/pegasus_strelka.sh $projDir >> $projDir/logs/pegasus_strelkaLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_mutect.txt)    echo "### Will run mutect $projDir"
        nohup $scriptsHome/pegasus_mutect.sh $projDir >> $projDir/logs/pegasus_mutectLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_clonalCov.txt)    echo "### Will run clonal coverage for $projDir"
        nohup $scriptsHome/pegasus_clonalCov.sh $projDir >> $projDir/logs/pegasus_clonalCovLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_cna15.txt)    echo "### Will run copy number analysis for $projDir"
        nohup $scriptsHome/pegasus_cna15.sh $projDir >> $projDir/logs/pegasus_cna15LOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_cna.txt)    echo "### Will run copy number analysis for $projDir"
        nohup $scriptsHome/pegasus_cna.sh $projDir >> $projDir/logs/pegasus_cnaLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_trn.txt)    echo "### Will run translocations analysis for $projDir"
        nohup $scriptsHome/pegasus_trn.sh $projDir >> $projDir/logs/pegasus_trnLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_haplotypeCaller.txt)    echo "### Will run haplotype caller for $projDir"
        nohup $scriptsHome/pegasus_haplotypeCaller.sh $projDir >> $projDir/logs/pegasus_haplotypeCallerLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_snpEff.txt)    echo "### Will run snpEff for $projDir"
        nohup $scriptsHome/pegasus_snpEff.sh $projDir >> $projDir/logs/pegasus_snpEffLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_snpSniff.txt)    echo "### Will run snpSniff for $projDir"
        nohup $scriptsHome/pegasus_snpSniffer.sh $projDir >> $projDir/logs/pegasus_snpSnifferLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_checkProjectComplete.txt)    echo "### Will check if project complete for $projDir"
        nohup $scriptsHome/pegasus_checkProjectComplete.sh $projDir >> $projDir/logs/pegasus_checkProjectCompleteLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_summaryStats.txt)    echo "### Will run summary stats for $projDir"
        nohup $scriptsHome/pegasus_summaryStats.sh $projDir >> $projDir/logs/pegasus_summaryStatsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_finalize.txt)    echo "### Will run finalize for $projDir"
        nohup $scriptsHome/pegasus_finalize.sh $projDir >> $projDir/logs/pegasus_finalizeLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_picardMultiMetrics.txt)    echo "### Will run picard Multi Metrics for $projDir"
        nohup $scriptsHome/pegasus_picardMultiMetrics.sh $projDir >> $projDir/logs/pegasus_picardMultiMetricsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_picardGcBiasMetrics.txt)    echo "### Will run picard GcBias Metrics for $projDir"
        nohup $scriptsHome/pegasus_picardGcBiasMetrics.sh $projDir >> $projDir/logs/pegasus_picardGcBiasMetricsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_picardHSMetrics.txt)    echo "### Will run picard HS Metrics for $projDir"
        nohup $scriptsHome/pegasus_picardHSMetrics.sh $projDir >> $projDir/logs/pegasus_picardHSMetricsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_phaseBT.txt)    echo "### Will run phase BT for $projDir"
        nohup $scriptsHome/pegasus_phaseBT.sh $projDir >> $projDir/logs/pegasus_phaseBTLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_splitNCigarReads.txt)    echo "### Will run splitNCigarReads for $projDir"
        nohup $scriptsHome/pegasus_splitNCigarReads.sh $projDir >> $projDir/logs/pegasus_splitNCigarReadsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_samtoolsStats.txt)    echo "### Will run samtools stats for $projDir"
        nohup $scriptsHome/pegasus_samtoolsStats.sh $projDir >> $projDir/logs/pegasus_samtoolsStatsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_splitNCigarReads.txt)    echo "### Will run splitNCigarReads for $projDir"
        nohup $scriptsHome/pegasus_splitNCigarReads.sh $projDir >> $projDir/logs/pegasus_splitNCigarLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_mergeFastqs.txt)    echo "### Will run mergeFastqs for $projDir"
        nohup $scriptsHome/pegasus_mergeFastqs.sh $projDir >> $projDir/logs/pegasus_mergeFastqsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_detectFusion.txt)    echo "### Will run detect fusion for $projDir"
        nohup $scriptsHome/pegasus_detectFusion.sh $projDir >> $projDir/logs/pegasus_detectFusionLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_sailFish.txt)    echo "### Will run sail fish for $projDir"
        nohup $scriptsHome/pegasus_sailFish.sh $projDir >> $projDir/logs/pegasus_sailFishLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_salmon.txt)    echo "### Will run salmon for $projDir"
        nohup $scriptsHome/pegasus_salmon.sh $projDir >> $projDir/logs/pegasus_salmonLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_tophatFusionPost.txt)    echo "### Will run tophat fusion post for $projDir"
        nohup $scriptsHome/pegasus_tophatFusionPost.sh $projDir >> $projDir/logs/pegasus_tophatFusionPostLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_rnaAlign.txt)    echo "### Will run rnaAlign for $projDir"
        nohup $scriptsHome/pegasus_rnaAlign.sh $projDir >> $projDir/logs/pegasus_rnaAlignLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_rnaMarkDup.txt)    echo "### Will run rna mark dups for $projDir"
        nohup $scriptsHome/pegasus_rnaMarkDup.sh $projDir >> $projDir/logs/pegasus_rnaMarkDupLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_RNAhaplotypeCaller.txt)    echo "### Will run rna haplotype caller for $projDir"
        nohup $scriptsHome/pegasus_RNAhaplotypeCaller.sh $projDir >> $projDir/logs/pegasus_rnaHaplotypeCallerLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_picardRNAMetrics.txt)    echo "### Will run picardRNAMetrics for $projDir"
        nohup $scriptsHome/pegasus_picardRNAMetrics.sh $projDir >> $projDir/logs/pegasus_picardRNAMetricsLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_IGLbedCov.txt)    echo "### Will run IGL bed cov for $projDir"
        nohup $scriptsHome/pegasus_IGLbedCov.sh $projDir >> $projDir/logs/pegasus_IGLbedCovLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_cuffLink.txt)    echo "### Will run cuff links for $projDir"
        nohup $scriptsHome/pegasus_cuffLink.sh $projDir >> $projDir/logs/pegasus_cuffLinkLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_cuffQuant.txt)    echo "### Will run cuff quant for $projDir"
        nohup $scriptsHome/pegasus_cuffQuant.sh $projDir >> $projDir/logs/pegasus_cuffQuantLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_htSeq.txt)    echo "### Will run HT seq for $projDir"
        nohup $scriptsHome/pegasus_htSeq.sh $projDir >> $projDir/logs/pegasus_htSeqLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_cuffDiff.txt)    echo "### Will run cuff diff for $projDir"
        nohup $scriptsHome/pegasus_cuffDiff.sh $projDir >> $projDir/logs/pegasus_cuffDiffLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_DEXseq.txt)    echo "### Will run DEXseq for $projDir"
        nohup $scriptsHome/pegasus_DEXseq.sh $projDir >> $projDir/logs/pegasus_DEXseqLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_DEXseqCount.txt)    echo "### Will run DEXseqCount for $projDir"
        nohup $scriptsHome/pegasus_DEXseqCount.sh $projDir >> $projDir/logs/pegasus_DEXseqCountLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_sleuth.txt)    echo "### Will run sleuth for $projDir"
        nohup $scriptsHome/pegasus_sleuth.sh $projDir >> $projDir/logs/pegasus_sleuthLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_kallisto.txt)    echo "### Will run kallisto for $projDir"
        nohup $scriptsHome/pegasus_kallisto.sh $projDir >> $projDir/logs/pegasus_kallistoLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_deSeq.txt)    echo "### Will run deSeq for $projDir"
        nohup $scriptsHome/pegasus_deSeq.sh $projDir >> $projDir/logs/pegasus_deSeqLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_freebayes.txt)    echo "### Will run freebayes for $projDir"
        nohup $scriptsHome/pegasus_freebayes.sh $projDir >> $projDir/logs/pegasus_freebayesLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_deSeq2.txt)    echo "### Will run deSeq2 for $projDir"
        nohup $scriptsHome/pegasus_deSeq2.sh $projDir >> $projDir/logs/pegasus_deSeq2LOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_vcfMerger.txt)    echo "### Will run vcf Merger for $projDir"
        nohup $scriptsHome/pegasus_vcfMerger.sh $projDir >> $projDir/logs/pegasus_vcfMergerLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_mergeVcfAlleleCount.txt)    echo "### Will run merge vcf allele count for $projDir"
        nohup $scriptsHome/pegasus_mergeVcfAlleleCount.sh $projDir >> $projDir/logs/pegasus_vcfMergerACLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_alleleCount.txt)    echo "### Will run alleleCount for $projDir"
        nohup $scriptsHome/pegasus_alleleCount.sh $projDir >> $projDir/logs/pegasus_alleleCountLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_deNovoGear.txt)    echo "### Will run deNovoGear for $projDir"
        nohup $scriptsHome/pegasus_deNovoGear.sh $projDir >> $projDir/logs/pegasus_deNovoGearLOG.txt 2>&1 &
        sleep 1
        ;;
    pegasus_nextJob_germVcfMerger.txt)    echo "### Will run germline merger for $projDir"
        nohup $scriptsHome/pegasus_germVcfMerge.sh $projDir >> $projDir/logs/pegasus_germVcfMergerLOG.txt 2>&1 &
        sleep 1
        ;;
    *)     echo "### Nothing to process $msgName with on $myhostname. Skipped."
        sleep 1
        ;;
    esac
done
