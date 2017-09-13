#!/usr/bin/env bash
#SBATCH --job-name="pegasus_strandhtSeq4Star"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
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

echo "### Now running htseq on SAM"
htseq-count -q --stranded=yes --mode=union ${SAM} ${GTF} > ${DIR}/${anotherName}.htSeqCounts

if [ $? -eq 0 ] ; then
    touch ${SAM}.htSeqPass
    echo "Deleted by htSeq to save on space at $time" > ${SAM}
else
    touch ${SAM}.htSeqFail
fi

rm -f ${SAM}.htSeqInQueue

time=`date +%d-%m-%Y-%H-%M`
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:HTSEQ:$hours:$mins" > ${SAM}.htSeq.totalTime

echo "TIME:$time finished htseq on ${BAM}"
