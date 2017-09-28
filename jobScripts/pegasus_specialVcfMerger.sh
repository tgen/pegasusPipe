#!/usr/bin/env bash
#SBATCH --job-name="pegasus_vcfMerger"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL

beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### SNPSIFT: ${SNPSIFT}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### SEURAT_VCF: ${SEURAT_VCF}"
echo "### MUTECT_VCF: ${MUTECT_VCF}"
echo "### MERGERDIR: ${MERGERDIR}"
echo "### SNPEFFPATH= $SNPEFFPATH"
echo "### SNPSIFT= $SNPSIFT"
echo "### SAMTOOLS= $SAMTOOLS"
echo "### VARSCAN= $VARSCAN"
echo "### MATCHEDNORMAL= $MATCHEDNORMAL"
echo "### DBSNP_SNV_BED= $DBSNP_SNV_BED"
echo "### DBSNP_DIV_BED= $DBSNP_DIV_BED"
echo "### REF= $REF"
echo "### DICT= $DICT"
echo "### COSMIC= $COSMIC"
echo "### COSMICC= $COSMICC"
echo "### COSMICNC= $COSMICNC"
echo "### EXAC= ${EXAC}"
echo "### KG= ${KG}"
echo "### NHLBI= ${NHLBI}"
echo "### SNPS= $SNPS"
echo "### INDELS= $INDELS"
echo "### GATK= $GATK"
echo "### VCFMERGER= $VCFMERGER"
echo "### VCFMERGER_DIR= $VCFMERGER_DIR"
echo "### VCFSORTER= $VCFSORTER"
echo "### RNA_VCF_HEADER= $RNA_VCF_HEADER"
echo "### POST_MERGE_VENN= $POST_MERGE_VENN"
echo "### DBSNP= $DBSNP"
echo "### DBVERSION= $DBVERSION"
echo "### SEURAT_VCF= $SEURAT_VCF"
echo "### MUTECT_VCF= $MUTECT_VCF"
echo "### STRELKA_SNV_VCF= $STRELKA_SNV_VCF"
echo "### STRELKA_INDEL_VCF= $STRELKA_INDEL_VCF"
echo "### MERGERDIR= $MERGERDIR"
echo "### RNABAM= $RNABAM"
echo "### ASSAYID= $ASSAYID"
echo "### BEDFILE= $BEDFILE"
echo "### RUNDIR= $RUNDIR"
echo "### CONTROL= $CONTROL"
echo "### TUMOR= $TUMOR"
echo "### D= $D"

module load BEDTools/2.14.0 
module load R/3.1.1

SEURAT_BASENAME=`basename ${SEURAT_VCF} ".seurat.vcf"`


# Merge vcf from the 3 callers
cd ${MERGERDIR}

# Annotate the merged vcf using GATK
echo "Annotate the merged VCF"
echo "..."

#Cleanup VCF for testing to remove seurat indels
java -Xmx24g -jar ${GATK}/GenomeAnalysisTK.jar -R ${REF} \
    -T VariantAnnotator \
    -nt 4 \
    -o ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.vcf \
    --variant ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.vcf \
    --dbsnp ${DBSNP} \
    --comp:EXAC ${EXAC} \
    --comp:NHLBI ${NHLBI} \
    --comp:1000G ${KG} \
    --comp:COSMIC_NC ${COSMICNC} \
    --comp:COSMIC_C ${COSMICC}

if [ $? -ne 0 ] ; then
    echo "### vcf merger failed at annotate merged vcf stage"
    mv ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerFail
    exit 1
fi

echo "Finished Annotating Merged VCF"


##### Determine if you need to run RNA allele counts or not
if [ -z "${RNABAM}" ] ; then
    echo "allele count was not requested for this pair"

    ## add dbNSFP annotaions
    java -jar ${SNPSIFT}/SnpSift.jar dbnsfp \
        -v ${DBNSFP} \
        -a \
        -f Interpro_domain,Polyphen2_HVAR_pred,GERP++_NR,GERP++_RS,LRT_score,MutationTaster_score,MutationAssessor_score,FATHMM_score,Polyphen2_HVAR_score,SIFT_score,Polyphen2_HDIV_score \
        ${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.vcf > ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.vcf

    if [ $? -ne 0 ] ; then
        echo "### vcf merger failed at annotate with dbNSFP stage"
        mv ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerFail
        exit 1
    fi

    # add snpEFF annotations
    java -Xmx4G -jar ${SNPEFFPATH}/snpEff.jar -canon -c ${SNPEFFPATH}/snpEff.config -v -lof ${DBVERSION} ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.vcf > ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.se74lofcan.vcf
    java -Xmx4G -jar ${SNPEFFPATH}/snpEff.jar -c ${SNPEFFPATH}/snpEff.config -v -lof ${DBVERSION} ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.vcf > ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.se74lof.vcf

    if [ $? -ne 0 ] ; then
        echo "### vcf merger failed at snpEff annotation stage"
        mv ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerFail
        exit 1
    fi

    # Make final call list venn
    ${POST_MERGE_VENN} --vcf ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.se74lofcan.vcf --outprefix  ${MERGERDIR}/${SEURAT_BASENAME}_finalVenn  --maintitle ${SEURAT_BASENAME} --
    if [ $? -ne 0 ] ; then
        echo "### vcf merger failed at venn diagram stage"
        mv ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerFail
        exit 1
    fi

    ##clean up final vcfs to save back
    mv ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.se74lofcan.vcf ${MERGERDIR}/${SEURAT_BASENAME}.merged.canonicalOnly.final.vcf
    mv ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.se74lof.vcf ${MERGERDIR}/${SEURAT_BASENAME}.merged.allTranscripts.final.vcf
    rm ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.vcf
    rm ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.vcf
    rm ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.vcf.idx
    rm ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.vcf
    rm ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.vcf
else
    echo "need to run allele count script and finalize for this DNAPAIR"
    touch ${RUNDIR}/${NXT1}
fi

rm ${MERGERDIR}/${STRELKA_SNV_VCF_BN}
rm ${MERGERDIR}/${STRELKA_INDEL_VCF_BN}

touch ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerPass

rm -rf ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:VCFMERGER:$hours:$mins" > ${MERGERDIR}/${SEURAT_BASENAME}.vcfMerger.totalTime
time=`date +%d-%m-%Y-%H-%M` 
echo "Ending snpEff Annotator."
