# Dependencies

## Modules

These modules are loaded by jobscripts:

```
BEDTools/2.14.0
BEDTools/2.26.0
HTSeq/0.6.0
kallisto/0.42.4
MCR/8.3
MCR/9.0
perl/5.24.1
python/2.7.3
python/3.4
R/2.14.1
R/2.15.2
R/3.1.1
R/3.2.1
  Vennerable
  VennDiagram
  gplots
  DESeq
  DESeq2
  DEXSeq
R/3.4.1
sailfish/0.6.3
salmon/0.6.0
```

## Binaries

These binaries are called directly from a path

```
/home/tgenref/pecan/bin/vt/vt
/home/tgenref/pecan/bin/java/jdk1.6.0_45/bin/java
```

## Other

```
/home/tgenref/pipeline_v0.4/cna/pos35.txt
pegasus_cna.sh:            bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_V5_PlusUTR_hs37d5_PaddedTargets_Picard_DGV_1kg_100_interval.bed"
pegasus_cna.sh:            bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_V6_noUTR_hs37d5_PaddedTargets_Picard_DGV_1kg_100_interval.bed"
pegasus_cna.sh:                        bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_SureSelect_V4_plusUTR/Agilent_SureSelect_V4_plusUTR_hs37d5_Padded_DGV_1kg.cna.bed"
pegasus_cna.sh:                        bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_SureSelect_V6R2_StrexomeProstate/Agilent_SureSelect_V6R2_StrexomeProstate_hs37d5_Padded_DGV_1kg.cna.bed"
pegasus_cna.sh:            bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_ClearSeq_Beta_ComprehensiveCancer/Agilent_ClearSeq_Beta_ComprehensiveCancer_hs37d5_Padded_DGV_1kg.cna.bed"
pegasus_cna.sh:            bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_SureSelect_CREv2_cliRes/Agilent_SureSelect_CREv2_cliRes_hs37d5_Padded_DGV_1kg.cna.bed"
pegasus_cna15.sh:                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_V5_PlusUTR_hs37d5_PaddedTargets_Picard_DGV_1kg_100_interval.bed"
pegasus_cna15.sh:                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_V6_noUTR_hs37d5_PaddedTargets_Picard_DGV_1kg_100_interval.bed"
pegasus_cna15.sh:                                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_SureSelect_V4_plusUTR/Agilent_SureSelect_V4_plusUTR_hs37d5_Padded_DGV_1kg.cna.bed"
pegasus_cna15.sh:                                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_SureSelect_V6R2_StxProstate/Agilent_SureSelect_V6R2_StxProstate_hs37d5_Padded_DGV_1kg.cna.bed"
pegasus_cna15.sh:                                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_SureSelect_V6R2_StrexomeV2/Agilent_SureSelect_V6R2_StrexomeV2_hs37d5_Padded_DGV_1kg.cna.bed"
pegasus_cna15.sh:                                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_SureSelect_V2_NA/Agilent_SureSelect_V2_NA_hs37d5_Padded_DGV_1kg.cna.bed"
pegasus_cna15.sh:                                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_ClearSeq_Beta_ComprehensiveCancer/Agilent_ClearSeq_Beta_ComprehensiveCancer_hs37d5_Padded_DGV_1kg.cna.bed"
pegasus_cna15.sh:                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_SureSelect_CREv2_cliRes/Agilent_SureSelect_CREv2_cliRes_hs37d5_Padded_DGV_1kg.cna.bed"
pegasus_cna15.sh:                    bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_V5_PlusUTR_hs37d5_PaddedTargets_Picard_DGV_1kg_100_interval.bed"
pegasus_cna15.sh:                    bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_V6_noUTR_hs37d5_PaddedTargets_Picard_DGV_1kg_100_interval.bed"
pegasus_detectFusion.sh:        spConfig=/home/tgenref/pecan/bin/SOAPfuse-v1.26/config/config.txt
pegasus_germVcfMerge.sh:DBNSFP=/home/tgenref/pecan/bin/vcfMerger/dbNSFP2.4.txt.gz
pegasus_germVcfMerge.sh:VCFMERGER=/home/tgenref/pecan/bin/vcfMerger/pecan.merge.3vcfs.main.sh
pegasus_germVcfMerge.sh:VCFMERGER_DIR=/home/tgenref/pecan/bin/vcfMerger
pegasus_germVcfMerge.sh:VCFSORTER=/home/tgenref/pecan/bin/vcfMerger/vcfsorter.pl
pegasus_germVcfMerge.sh:RNA_VCF_HEADER=/home/tgenref/pecan/bin/vcfMerger/RNA_VCF_HEADER.vcf
pegasus_germVcfMerge.sh:POST_MERGE_VENN=/home/tgenref/pecan/bin/vcfMerger/pecan.Venn_postMMRF_specific_filtering.sh
pegasus_germVcfMerge.sh:#DBSNP=/home/tgenref/pecan/bin/vcfMerger/dbsnp_137.b37.vcf
pegasus_germVcfMerge.sh:        bedFile="/home/tgenref/pipeline_v0.3/annotations/exome_capture/illumina_truseq/TruSeq_exome_targeted_regions_b37_padded.bed"
pegasus_germVcfMerge.sh:        bedFile="/home/tgenref/pipeline_v0.3/ensembl70/Ensembl_v70_hs37d5_exonic_coordinates_touched_v5UTR_padded25.bed"
pegasus_germVcfMerge.sh:        bedFile="/home/tgenref/pipeline_v0.4/annotations/exome_capture/strexome/Strexome_targets_sortedTabs2_padded25.bed"
pegasus_germVcfMerge.sh:                bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_clinical_research_exome/Agilent_Clinical_Research_Exome_hs37d5_TargetsPadded25_sortedTabs2_Picard.txt"
pegasus_germVcfMerge.sh:                bedFile="/home/tgenref/pipeline_v0.3/annotations/exome_capture/illumina_truseq/TruSeq_exome_targeted_regions_b37_padded.bed"
pegasus_germVcfMerge.sh:                bedFile="/home/tgenref/pipeline_v0.3/ensembl70/Ensembl_v70_hs37d5_exonic_coordinates_touched_v5UTR_padded25.bed"
pegasus_germVcfMerge.sh:                bedFile="/home/tgenref/pipeline_v0.4/annotations/exome_capture/strexome/Strexome_targets_sortedTabs2_padded25.bed"
pegasus_germVcfMerge.sh:                bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_clinical_research_exome/Agilent_Clinical_Research_Exome_hs37d5_TargetsPadded25_sortedTabs2_Picard.txt"
pegasus_germlineACmerge.sh:DBNSFP=/home/tgenref/pecan/bin/vcfMerger/dbNSFP2.4.txt.gz
pegasus_germlineACmerge.sh:VCFMERGER=/home/tgenref/pecan/bin/vcfMerger_V2_norm/production.version.VCFMERGER_V2.WithNorm.20150603/pecan.merge.3vcfs.main.sh
pegasus_germlineACmerge.sh:VCFMERGER_DIR=/home/tgenref/pecan/bin/vcfMerger_V2_norm/production.version.VCFMERGER_V2.WithNorm.20150603
pegasus_germlineACmerge.sh:VCFSORTER=/home/tgenref/pecan/bin/vcfMerger/vcfsorter.pl
pegasus_germlineACmerge.sh:RNA_VCF_HEADER=/home/tgenref/pecan/bin/vcfMerger/RNA_VCF_HEADER.vcf
pegasus_germlineACmerge.sh:POST_MERGE_VENN=/home/tgenref/pecan/bin/vcfMerger_V2_norm/production.version.VCFMERGER_V2.WithNorm.20150603/pecan.Venn_postMMRF_specific_filtering.sh
pegasus_germlineACmerge.sh:KG=/home/tgenref/pecan/bin/vcfMerger/1000G_phase1.snps.high_confidence.b37.vcf
pegasus_germlineACmerge.sh:COSMIC=/home/tgenref/pecan/bin/vcfMerger/CosmicCodingMuts_v66_20130725_sorted.vcf
pegasus_germlineACmerge.sh:NHLBI=/home/tgenref/pecan/bin/vcfMerger/ESP6500SI-V2_snps_indels.vcf
pegasus_germlineACmerge.sh:DBSNP=/home/tgenref/pecan/bin/vcfMerger/dbsnp_137.b37.vcf
pegasus_mergeVcfAlleleCount.sh:DBNSFP=/home/tgenref/pecan/bin/vcfMerger/dbNSFP2.4.txt.gz
pegasus_mergeVcfAlleleCount.sh:VCFMERGER=/home/tgenref/pecan/bin/vcfMerger_V2_norm/production.version.VCFMERGER_V2.WithNorm.201600819/pecan.merge.3vcfs.main.sh
pegasus_mergeVcfAlleleCount.sh:VCFMERGER_DIR=/home/tgenref/pecan/bin/vcfMerger_V2_norm/production.version.VCFMERGER_V2.WithNorm.201600819
pegasus_mergeVcfAlleleCount.sh:VCFSORTER=/home/tgenref/pecan/bin/vcfMerger/vcfsorter.pl
pegasus_mergeVcfAlleleCount.sh:RNA_VCF_HEADER=/home/tgenref/pecan/bin/vcfMerger/RNA_VCF_HEADER.vcf
pegasus_mergeVcfAlleleCount.sh:POST_MERGE_VENN=/home/tgenref/pecan/bin/vcfMerger_V2_norm/production.version.VCFMERGER_V2.WithNorm.201600819/pecan.Venn_postMMRF_specific_filtering.sh
pegasus_mergeVcfAlleleCount.sh:KG=/home/tgenref/pecan/bin/vcfMerger/1000G_phase1.snps.high_confidence.b37.vcf
pegasus_mergeVcfAlleleCount.sh:COSMIC=/home/tgenref/pecan/bin/vcfMerger/CosmicCodingMuts_v66_20130725_sorted.vcf
pegasus_mergeVcfAlleleCount.sh:NHLBI=/home/tgenref/pecan/bin/vcfMerger/ESP6500SI-V2_snps_indels.vcf
pegasus_mergeVcfAlleleCount.sh:DBSNP=/home/tgenref/pecan/bin/vcfMerger/dbsnp_137.b37.vcf
pegasus_vcfMerger.sh:DBNSFP=/home/tgenref/pecan/bin/vcfMerger/dbNSFP2.4.txt.gz
pegasus_vcfMerger.sh:VCFMERGER_DIR=/home/tgenref/pecan/bin/mergeVcf/normalization/vcfMerger
pegasus_vcfMerger.sh:VCFSORTER=/home/tgenref/pecan/bin/vcfMerger/vcfsorter.pl
pegasus_vcfMerger.sh:RNA_VCF_HEADER=/home/tgenref/pecan/bin/vcfMerger/RNA_VCF_HEADER.vcf
pegasus_vcfMerger.sh:COSMIC=/home/tgenref/pecan/bin/vcfMerger/CosmicCodingMuts_v66_20130725_sorted.vcf
pegasus_vcfMerger.sh:#DBSNP=/home/tgenref/pecan/bin/vcfMerger/dbsnp_137.b37.vcf
pegasus_vcfMerger.sh:DBSNP_DIV_BED=/home/tgenref/pecan/annotations/dbsnp/Merged_Indels_no_COSMIC_Y.bed
pegasus_vcfMerger.sh:DBSNP_SNV_BED=/home/tgenref/pecan/annotations/dbsnp/Merged_SNPs_no_COSMIC_Y.bed
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pipeline_v0.3/annotations/exome_capture/illumina_truseq/TruSeq_exome_targeted_regions_b37_padded.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pipeline_v0.3/ensembl70/Ensembl_v70_hs37d5_exonic_coordinates_touched_v5UTR_padded25.bed"
pegasus_vcfMerger.sh:        #bedFile="/home/tgenref/pipeline_v0.4/annotations/exome_capture/strexome/Strexome_targets_sortedTabs2_padded150.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pipeline_v0.4/annotations/exome_capture/strexome/Strexome_targets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        #bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_clinical_research_exome/Agilent_Clinical_Research_Exome_hs37d5_TargetsPadded25_sortedTabs2_Picard.txt"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_clinical_research_exome/Agilent_Clinical_Research_Exome_hs37d5_PaddedTargets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        #bedFile="/home/tgenref/pipeline_v0.4/annotations/exome_capture/Agilent_V2_hs37d5/Agilent_V2_hs37d5_TargetsPadded25sorted_Picard.txt"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pipeline_v0.4/annotations/exome_capture/Agilent_V2_hs37d5/Agilent_V2_hs37d5_Targets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        #bedFile="/home/tgenref/pecan/annotations/exome_capture/strexome_lite/temp/Strexome_Lite_Targets_padded25.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/strexome_lite/Strexome_Lite_Targets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        #bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_SureSelectV1/temp/SureSelectV1_hs37d5_TargetsPadded25_Picard.txt"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_SureSelectV1/SureSelectV1_hs37d5_PaddedTargets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        #bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_v6_noUTR/Agilent_V6_noUTR_hs37d5_TargetsPadded25.txt"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_v6_noUTR/Agilent_V6_noUTR_hs37d5_Targets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        #bedFile="/home/tgenref/pecan/annotations/exome_capture/prostateStrexome/prostateStrexome.targetsPadded25.txt"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/Agilent_SureSelect_V6R2_plusUTR/Agilent_SureSelect_V6R2_plusUTR_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_v4_noUTR/Agilent_V4_noUTR_hs37d5_Targets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/annotations/dog/canfam3/vcfMergerBed/agilent_canine_exonV2_targets.padded100.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/StrAD/StrAD_targets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/Agilent_SureSelect_V4_plusUTR/Agilent_SureSelect_V4_plusUTR_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/illumina_nextera_expanded/NexteraExpandedExome_hs37d5_Targets_PicardPadded100.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/Agilent_SureSelect_V6R2_StxProstate/Agilent_SureSelect_V6R2_StxProstate_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/Agilent_SureSelect_V6R2_StrexomeV2/Agilent_SureSelect_V6R2_StrexomeV2_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/Agilent_SureSelect_V2_NA/Agilent_SureSelect_V2_NA_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/Agilent_SureSelect_V6R2_plusCOSMIC/Agilent_SureSelect_V6R2_plusCOSMIC_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/Agilent_ClearSeq_Beta_ComprehensiveCancer/Agilent_ClearSeq_Beta_ComprehensiveCancer_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/Agilent_SureSelect_CREv2_cliRes/Agilent_SureSelect_CREv2_cliRes_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
pegasus_vcfMerger.sh:        bedFile="/home/tgenref/pecan/annotations/exome_capture/illumina_nextera_v1.2/padded_targets.bed"
```
