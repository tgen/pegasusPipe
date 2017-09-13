#!/usr/bin/env bash
#SBATCH --job-name="pegasus_reassignMap"
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
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### GATKPATH: ${GATKPATH}"
echo "### RNABAM: ${RNABAM}"
echo "### OUTBAM: ${OUTBAM}"

echo "### GATK reassign Mapping Quality started at $time."
java -Djava.io.tmpdir=$TMPDIR -jar -Xmx32g ${GATKPATH}/GenomeAnalysisTK.jar \
    -l INFO \
    -R ${REF} \
    -T PrintReads \
    -rf ReassignMappingQuality \
    -I ${RNABAM} \
    -o ${OUTBAM} > ${RNABAM}.reassignMapOut

if [ $? -eq 0 ] ; then
    mv ${RNABAM}.reassignMapOut ${RNABAM}.reassignMapPass
    touch ${RUNDIR}/${NXT1}
else
    mv ${RNABAM}.reassignMapOut ${RNABAM}.reassignMapFail
    rm -f ${RNABAM}.reassignMapInQueue
    exit 1
fi

rm -f ${RNABAM}.reassignMapInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:GATKreassignMap:$hours:$mins" > ${RNABAM}.reassignMap.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "reassignMap finished at $time."
