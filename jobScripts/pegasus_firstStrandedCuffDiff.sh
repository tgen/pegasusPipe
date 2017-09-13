#!/usr/bin/env bash
#SBATCH --job-name="pegasus_FScuffDiff"
#SBATCH --time=0-240:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL

 
time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### DIRNAME: ${DIRNAME}"
echo "### RUNDIR: ${RUNDIR}"
echo "### CDPATH: ${CUFFDIFFPATH}"
echo "### GTF: ${GTF}"
echo "### MASK: ${MASK}"
echo "### NXT1: ${NXT1}"
echo "### BAM1: ${BAM1}"
echo "### BAM2: ${BAM2}"

echo "TIME:$time starting cuffdiff on ${BAM1} and ${BAM2}"
echo "running first stranded cufflinks options for this RNA"
cd ${DIRNAME}
newName=`basename ${DIRNAME}`
newName=${newName/_cdDir/}
#update-check ${GTF} ${BAM1} ${BAM2}
${CUFFDIFFPATH}/cuffdiff -p 16 -N -M ${MASK} -b ${REF} --library-type fr-secondstrand -L Control,Tumor ${GTF} ${BAM1} ${BAM2} > ${DIRNAME}.cuffDiffOut
if [ $? -eq 0 ] ; then
    mv ${DIRNAME}/tss_groups.fpkm_tracking ${DIRNAME}/${newName}_tss_groups.fpkm_tracking
    mv ${DIRNAME}/isoforms.fpkm_tracking ${DIRNAME}/${newName}_isoforms.fpkm_tracking
    mv ${DIRNAME}/genes.fpkm_tracking ${DIRNAME}/${newName}_genes.fpkm_tracking
    mv ${DIRNAME}/cds.fpkm_tracking ${DIRNAME}/${newName}_cds.fpkm_tracking
    mv ${DIRNAME}/tss_group_exp.diff ${DIRNAME}/${newName}_tss_group_exp.diff
    mv ${DIRNAME}/splicing.diff ${DIRNAME}/${newName}_splicing.diff
    mv ${DIRNAME}/promoters.diff ${DIRNAME}/${newName}_promoters.diff
    mv ${DIRNAME}/isoform_exp.diff ${DIRNAME}/${newName}_isoform_exp.diff
    mv ${DIRNAME}/gene_exp.diff ${DIRNAME}/${newName}_gene_exp.diff
    mv ${DIRNAME}/cds_exp.diff ${DIRNAME}/${newName}_cds_exp.diff
    mv ${DIRNAME}/cds.diff ${DIRNAME}/${newName}_cds.diff

    mv ${DIRNAME}/isoforms.count_tracking ${DIRNAME}/${newName}_isoforms.count_tracking
    mv ${DIRNAME}/tss_groups.count_tracking ${DIRNAME}/${newName}_tss_groups.count_tracking
    mv ${DIRNAME}/cds.count_tracking ${DIRNAME}/${newName}_cds.count_tracking
    mv ${DIRNAME}/genes.count_tracking ${DIRNAME}/${newName}_genes.count_tracking
    mv ${DIRNAME}/isoforms.read_group_tracking ${DIRNAME}/${newName}_isoforms.read_group_tracking
    mv ${DIRNAME}/tss_groups.read_group_tracking ${DIRNAME}/${newName}_tss_groups.read_group_tracking
    mv ${DIRNAME}/cds.read_group_tracking ${DIRNAME}/${newName}_cds.read_group_tracking
    mv ${DIRNAME}/genes.read_group_tracking ${DIRNAME}/${newName}_genes.read_group_tracking
    mv ${DIRNAME}/read_groups.info ${DIRNAME}/${newName}_read_groups.info
    mv ${DIRNAME}/run.info ${DIRNAME}/${newName}_run.info

    echo "starting 3 external scripts"
    ${PROCESSCDLISTPATH}/processCuffDiffLists.sh ${DIRNAME}/${newName}_genes.fpkm_tracking ${DIRNAME}/${newName}_gene_exp.diff
    ${PROCESSCDLISTPATH}/processCuffDiffLists.sh ${DIRNAME}/${newName}_isoforms.fpkm_tracking ${DIRNAME}/${newName}_isoform_exp.diff
    ${CUFFDIFF2VCFPATH}/cuffdiff2vcf.pl ${DIRNAME}/${newName}_gene_exp.diff ${BAM1}
    echo "done with 3 external scripts"

    mv ${DIRNAME}.cuffDiffOut ${DIRNAME}.cuffDiffPass
    #touch ${RUNDIR}/${NXT1}
else
    mv ${DIRNAME}.cuffDiffOut ${DIRNAME}.cuffDiffFail
fi
rm -f ${DIRNAME}.cuffDiffInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:CUFFDIFF:$hours:$mins" > ${DIRNAME}.cuffdiff.totalTime
time=`date +%d-%m-%Y-%H-%M`
echo "TIME:$time finished cuffdiff on ${BAM1} and ${BAM2}"
