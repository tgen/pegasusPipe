#!/usr/bin/env bash
#SBATCH --job-name="pegasus_RNAhc"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8

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
java -Djava.io.tmpdir=$TMPDIR -jar -Xmx32g ${GATKPATH}/GenomeAnalysisTK.jar \
    -l INFO \
    -R ${REF} \
    -L ${CHRLIST}/Step${STEP}.list \
    -nct 8 \
    -T HaplotypeCaller \
    -I ${BAMLIST} \
    -D ${KNOWN} \
    -dontUseSoftClippedBases \
    -stand_call_conf 20 \
    -stand_emit_conf 20 \
    -mbq 10 \
    -o ${TRK}_Step${STEP}.rnaHC.vcf > ${TRK}_Step${STEP}.RNAhcOut

if [ $? -eq 0 ] ; then
    echo "${STEP} Completed" >> ${TRK}_RNAhcStatus.txt
    PROGRESS=`wc -l ${TRK}_RNAhcStatus.txt | awk '{print $1}'`
    mv ${TRK}_Step${STEP}.RNAhcOut ${TRK}_Step${STEP}.RNAhcPass
else
    mv ${TRK}_Step${STEP}.RNAhcOut ${TRK}_Step${STEP}.RNAhcFail
    rm -f ${TRK}_Step${STEP}.RNAhcInQueue
    exit 1
fi

vcfList=""
# Here we make a look to create the list of vcfs based on STEPCOUNT
for i in `seq 1 ${STEPCOUNT}`;
do
    thisVcf="-V ${TRK}_Step$i.rnaHC.vcf "
    vcfList="$vcfList $thisVcf"
done

# IF the progress count equals the step count merge to single vcf
if [ ${PROGRESS} -eq ${STEPCOUNT} ]
then
    echo HapCaller_${STEP}.Done
    # Concatenate VCF with GATK
     java -cp ${GATKPATH}/GenomeAnalysisTK.jar \
         org.broadinstitute.gatk.tools.CatVariants \
         -R ${REF} \
         $vcfList -out \
         ${TRK}.rnaHC_All.vcf -assumeSorted
    if [ $? -eq 0 ] ; then
        touch ${TRK}.RNAhcPass
        touch ${RUNDIR}/${NXT1}
    else
        touch ${TRK}.RNAhcFail
    fi

    mv ${TRK}_RNAhcStatus.txt ${TRK}_RNAhcStatus.txt.used
else
    echo
    echo HapCaller_${STEP}.Done
fi

rm -f ${TRK}_Step${STEP}.RNAhcInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:GATKrnaHC:$hours:$mins" > ${TRK}_Step${STEP}.hapCal.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "Haplotypecaller finished at $time."
