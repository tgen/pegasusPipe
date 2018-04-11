#!/usr/bin/env bash
#SBATCH --job-name="pegasus_hc"
#SBATCH --time=0-170:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL


time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### CHRLIST: ${CHRLIST}"
echo "### STEP: ${STEP}"
echo "### STEPCOUNT: ${STEPCOUNT}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### GATKPATH: ${GATKPATH}"
echo "### KNOWN: ${KNOWN}"
echo "### BAMLIST: ${BAMLIST}"

echo "### Haplotype caller started for multiple bams at $time."
echo "java -Djava.io.tmpdir=$TMPDIR -jar -Xmx24g ${GATKPATH}/GenomeAnalysisTK.jar \
-l INFO \
-R ${REF} \
-L ${CHRLIST}/Step${STEP}.list \
-nct 8 \
-T HaplotypeCaller \
${BAMLIST} \
-D ${KNOWN} \
-mbq 10 \
-o ${TRK}_Step${STEP}.HC.vcf"

java -Djava.io.tmpdir=$TMPDIR -jar -Xmx24g ${GATKPATH}/GenomeAnalysisTK.jar \
-l INFO \
-R ${REF} \
-L ${CHRLIST}/Step${STEP}.list \
-nct 8 \
-T HaplotypeCaller \
${BAMLIST} \
-D ${KNOWN} \
-mbq 10 \
-o ${TRK}_Step${STEP}.HC.vcf

# HaplotypeCaller jobs are split over each chromosome, here 
# we discover the progress of the entire job by counting
# the number of hcPass files. There should be one for each
# step in STEPCOUNT.
if [ $? -eq 0 ] ; then
    echo "${SLURM_JOB_ID}" > ${TRK}_Step${STEP}.hcPass
    PROGRESS=$(ls ${TRK}*hcPass | wc -l)
else
    echo "${SLURM_JOB_ID}" > ${TRK}_Step${STEP}.hcFail
    rm -f ${TRK}_Step${STEP}.hcInQueue
    exit 1
fi

# Here vcfList is the arguments built for the following merger
# step by appending "-V <vcf>" for each vcf in STEPCOUNT.
vcfList=""
for i in `seq 1 ${STEPCOUNT}`;
do
        thisVcf="-V ${TRK}_Step$i.HC.vcf "
        vcfList="$vcfList $thisVcf"
done

# IF the progress count equals the step count merge to single vcf
if [ ${PROGRESS} -eq ${STEPCOUNT} ]; then
    echo "PROGRESS: ${PROGRESS} equals STEPCOUNT: ${STEPCOUNT}"
    echo HapCaller_${STEP}.Done

    # Concatenate VCF with GATK
     java -cp ${GATKPATH}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R ${REF} $vcfList -out ${TRK}.HC_All.vcf -assumeSorted

    if [ $? -eq 0 ] ; then
        touch ${TRK}.hcPass
        touch ${RUNDIR}/${NXT1}
        touch ${RUNDIR}/${NXT2}
        touch ${RUNDIR}/${NXT3}
    else
        touch ${TRK}.hcFail
    fi

else
    echo HapCaller_${STEP}.Done
fi

rm -f ${TRK}_Step${STEP}.hcInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:GATKHC:$hours:$mins" > ${TRK}_Step${STEP}.hapCal.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "Haplotypecaller finished at $time."

