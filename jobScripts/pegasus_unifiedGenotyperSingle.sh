#!/usr/bin/env bash
#SBATCH --job-name="pegasus_uniGenoSin"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### GATKPATH: ${GATKPATH}"
echo "### KNOWN: ${KNOWN}"
echo "### BAMLIST: ${BAMLIST}"

echo "### UnifiedGenotyper started for multiple bams at $time."
java -Djava.io.tmpdir=$TMPDIR -jar -Xmx32g ${GATKPATH}/GenomeAnalysisTK.jar \
    -l INFO \
    -R ${REF} \
    -T UnifiedGenotyper \
    -rf BadCigar \
    -nct 3 \
    -nt 8 \
    -glm BOTH \
    -I ${BAMLIST} \
    -D ${KNOWN} \
    -mbq 10 \
    -metrics ${TRK}.metrics.txt \
    -o ${TRK}.UG.vcf > ${TRK}.ugOut

if [ $? -eq 0 ] ; then
    mv ${TRK}.ugOut ${TRK}.ugPass
    touch ${RUNDIR}/${NXT1}
else
    mv ${TRK}.ugOut ${TRK}.ugFail
fi

rm -f ${TRK}.ugInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:GATKUG:$hours:$mins" > ${TRK}.uniGeno.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "UnifiedGenotyper finished at $time."
