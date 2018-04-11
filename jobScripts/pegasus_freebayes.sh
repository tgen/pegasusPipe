#!/usr/bin/env bash
#SBATCH --job-name="pegasus_freebayes"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### KNOWN: ${KNOWN}"
echo "### STEP: ${STEP}"
echo "### FBBAM: ${FBBAM}"
echo "### FREEBAYESPATH: ${FREEBAYESPATH}"
echo "### TRACKNAME: ${TRACKNAME}"
echo "### CHRLIST: ${CHRLIST}"

echo "### freebayes started at $time."
echo "${FREEBAYESPATH}/freebayes -f ${REF} -b ${FBBAM} -t ${CHRLIST}/Step${STEP}.bed --ploidy 2 --min-repeat-entropy 1 > ${TRACKNAME}_Step${STEP}.freebayes.vcf"
${FREEBAYESPATH}/freebayes -f ${REF} -b ${FBBAM} -t ${CHRLIST}/Step${STEP}.bed --ploidy 2 --min-repeat-entropy 1 > ${TRACKNAME}_Step${STEP}.freebayes.vcf

# Freebayes jobs are split over each chromosome, here
# we discover the progress of the entire job by counting
# the number of freebayesPass files. There should be one 
# for each step in STEPCOUNT.
if [ $? -eq 0 ] ; then
    echo Freebayes_${STEP}.Done
    echo "${SLURM_JOB_ID}" > ${TRACKNAME}_Step${STEP}.freebayesPass
    PROGRESS=$(ls ${TRACKNAME}*freebayesPass | wc -l)
else
    touch ${TRACKNAME}_Step${STEP}.freebayesFail
    rm -f ${TRACKNAME}_Step${STEP}.freebayesInQueue
    exit 1
fi

# Here we make a look to create the list of vcfs based on STEPCOUNT
vcfList=""
for i in `seq 1 ${STEPCOUNT}`;
do
    thisVcf="-V ${TRACKNAME}_Step$i.freebayes.vcf "
    vcfList="$vcfList $thisVcf"
done

# IF the progress count equals the step count merge to single vcf
if [ ${PROGRESS} -eq ${STEPCOUNT} ]; then
    # Concatenate VCF with GATK
    java -cp ${GATKPATH}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R ${REF} $vcfList -out ${TRACKNAME}.freebayes_All.vcf -assumeSorted
    if [ $? -eq 0 ] ; then
        touch ${TRACKNAME}.freebayesPass
        touch ${RUNDIR}/${NXT1}
        touch ${RUNDIR}/${NXT2}
    else
        touch ${TRACKKNAME}.freebayesFail
    fi
fi

rm -f ${TRACKNAME}_Step${STEP}.freebayesInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:FREEBAYES:$hours:$mins" > ${TRACKNAME}_Step${STEP}.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "Freebayes finished at $time."

