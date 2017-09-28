#!/usr/bin/env bash
#SBATCH --job-name="pegasus_htSeq"
#SBATCH --time=0-96:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL


module load python/2.7.3

beginTime=`date +%s`
time=`date +%d-%m-%Y-%H-%M`
machine=`hostname`
echo "### NODE: $machine"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### PICARDPATH: ${PICARDPATH}"
echo "### BAM: ${BAM}"
echo "### SAM: ${SAM}"

base=`basename ${SAM}`
anotherName=${base/.proj.Aligned.out.sam}
DIR=$(dirname "${SAM}")

htseq-count -q --format=bam --stranded=no --mode=union ${BAM} ${GTF} > ${DIR}/${anotherName}.htSeqCounts

if [ $? -eq 0 ] ; then
    touch ${SAM}.htSeqPass
else
    mv ${SAM}.htSeqOut ${SAM}.htSeqFail
fi
rm -f ${SAM}.htSeqInQueue

time=`date +%d-%m-%Y-%H-%M`
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:HTSEQ:$hours:$mins" > ${SAM}.htSeq.totalTime
echo "TIME:$time finished htseq on ${BAM}"
