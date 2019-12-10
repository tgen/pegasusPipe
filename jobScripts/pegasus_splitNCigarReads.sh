#!/usr/bin/env bash
#SBATCH --job-name="pegasus_splitNCigar"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem-per-cpu 4096
#SBATCH --cpus-per-task 14

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### GATKPATH: ${GATKPATH}"
echo "### RNABAM: ${RNABAM}"
echo "### OUTBAM: ${OUTBAM}"

echo "### GATK splitNCigarReads started at $time."

# All mapping qualities of 255 will be reassigned to 60 in the bam
java -Xmx48G -Djava.io.tmpdir=$TMPDIR -jar ${GATKPATH}/GenomeAnalysisTK.jar \
    -l INFO \
    -R ${REF} \
    -T SplitNCigarReads \
    -I ${RNABAM} \
    -o ${OUTBAM} \
    -rf ReassignOneMappingQuality \
    -RMQF 255 \
    -RMQT 60 \
    -U ALLOW_N_CIGAR_READS > ${RNABAM}.splitNCigarOut

if [ $? -eq 0 ] ; then
    mv ${RNABAM}.splitNCigarOut ${RNABAM}.splitNCigarPass
    touch ${RUNDIR}/${NXT1}
else
    mv ${RNABAM}.splitNCigarOut ${RNABAM}.splitNCigarFail
    rm -f ${RNABAM}.splitNCigarInQueue
    exit 1
fi

rm -f ${RNABAM}.splitNCigarInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:GATKsplitNCigar:$hours:$mins" > ${RNABAM}.splitNCigar.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "splitNCigar finished at $time."
