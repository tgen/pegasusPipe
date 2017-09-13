#!/usr/bin/env bash
#SBATCH --job-name="pegasus_germVCFmerge"
#SBATCH --time=0-96:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine" >> ${TRACKNAME}.vcfMergerInQueue
echo "### JOBID: ${SLURM_JOB_ID}" >> ${TRACKNAME}.vcfMergerInQueue
echo "### REF: ${REF}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### SAMTOOLSPATH: ${SAMTOOLSPATH}"
echo "### BCFTOOLSPATH: ${BCFTOOLSPATH}"
echo "### BAMFILE: ${BAMFILE}"
echo "### D: ${D}"
echo "### BEDFILE: ${BEDFILE}"
echo "### CHRLIST: ${CHRLIST}"
echo "### GATKPATH: ${GATKPATH}"
echo "### HCVCF: $HCVCF"
echo "### FBVCF: $FBVCF"
echo "### MPVCF: $MPVCF"
echo "### VTPATH: $VTPATH"
echo "### TRACKNAME: ${TRACKNAME}"

echo "### Germline vcf merger started at $time."
${SAMTOOLSPATH}/samtools view -x BD -x BI ${BAM} | head

##Run the actual merging script
java -jar ${GATKPATH}/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R ${REF} \
    --variant:haplotypeCaller ${TRACKNAME}.HC.norm.vcf \
    --variant:freebayes ${TRACKNAME}.freebayes.norm.vcf \
    --variant:mpileup ${TRACKNAME}.mpileup.norm.vcf \
    --genotypemergeoption UNIQUIFY \
    --out ${TRACKNAME}.merged.vcf

if [ $? -eq 0 ] ; then
    echo "the 3 variant callers were merged successfully"
    touch ${TRACKNAME}.vcfMergerPass
else
    touch ${TRACKNAME}.vcfMergerFail
fi

rm -f ${TRACKNAME}.vcfMergerInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:GERMVCFMERGE:$hours:$mins" > ${TRACKNAME}.germVcfMerge.totalTime
time=`date +%d-%m-%Y-%H-%M`
echo "Germline vcf Merger finished at $time."
