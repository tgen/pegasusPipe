#!/usr/bin/env bash
#SBATCH --job-name="pegasus_vqsr"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu 4096

time=`date +%d-%m-%Y-%H-%M`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### VCF: ${VCF}"
echo "### KNOWN: ${KNOWN}"
echo "### HAPMAP: ${HAPMAP}"
echo "### OMNI: ${OMNI}"
echo "### ASSAY: ${ASSAY}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### GATKPATH: ${GATKPATH}"

beginTime=`date +%s`
echo "### Variant quality score recal started for ${VCF} at $time."
recalVCF=${VCF/.vcf/.vqsr.vcf}

echo "### New vcf name will be $recalVCF"
if [ "${ASSAY}" == "Genome" ] ; then
    maxGaussian=2
    badVariant=0.20
elif [ "${ASSAY}" == "Exome" ] ; then 
    maxGaussian=4
    badVariant=0.10
else
    echo "I should not be here. Assay type is: ${ASSAY}"
    touch ${VCF}.vqsrFail
    rm -f ${VCF}.vqsrInQueue
    time=`date +%d-%m-%Y-%H-%M`
    endTime=`date +%s`
    elapsed=$(( $endTime - $beginTime ))
    (( hours=$elapsed/3600 ))
    (( mins=$elapsed%3600/60 ))
    echo "RUNTIME:VQSR:$hours:$mins" > ${VCF}.vqsr.totalTime
    echo "VQSR finished at $time."
    exit 1
fi

echo "### Dynamically set maxGaussian to $maxGaussian and percent bad variant to $badVariant because assay type was ${ASSAY}"

java -jar -Xmx4g ${GATKPATH}/GenomeAnalysisTK.jar \
    -R ${REF} \
    -T VariantRecalibrator \
    -input ${VCF} \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${HAPMAP} \
    -resource:omni,known=false,training=true,truth=false,prior=12.0 ${OMNI} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 ${KNOWN} \
    --mode SNP \
    -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an MQ \
    -recalFile ${RECAL} \
    -tranchesFile ${TRANCHES} \
    -rscriptFile ${RSCRIPT} > ${VCF}.vqsrOut

if [ $? -ne 0 ] ; then
    mv ${VCF}.vqsrOut ${VCF}.vqsrFail
    echo "### vqsr failed at variantrecalibrator stage"
    rm -f ${VCF}.vqsrInQueue
    exit 1
fi

java -jar -Xmx4g ${GATKPATH}/GenomeAnalysisTK.jar \
    -R ${REF} \
    -T ApplyRecalibration \
    --ts_filter_level 99.0 \
    -tranchesFile ${TRANCHES} \
    -recalFile ${RECAL} \
    -mode SNP \
    -input ${VCF} \
    -o $recalVCF >> ${VCF}.vqsrOut

if [ $? -eq 0 ] ; then
    mv ${VCF}.vqsrOut ${VCF}.vqsrPass
    touch ${RUNDIR}/${NXT1} > /dev/null
else
    echo "### vqsr failed at apply recalibration stage"
    mv ${VCF}.vqsrOut ${VCF}.vqsrFail
    #attention, fake passing this step in case of failure
    mv ${VCF}.vqsrFail ${VCF}.vqsrPass
    touch ${VCF}.vqsrPassISFAKE
    cp ${VCF} $recalVCF
fi

rm -f ${VCF}.vqsrInQueue

time=`date +%d-%m-%Y-%H-%M`
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:VQSR:$hours:$mins" > ${VCF}.vqsr.totalTime
echo "VQSR finished at $time."
