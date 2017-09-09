#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_DEXseqCount"
#SBATCH --time=0-60:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
##PBS -e /${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err

module load R/3.2.1
module load HTSeq/0.6.0
module load python/2.7.3 

beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### DIR: ${DIR}"
echo "### NXT1: ${NXT1}"
echo "### RUNDIR: ${RUNDIR}"
echo "### DEXSEQCOUNTPATH: ${DEXSEQCOUNTPATH}"
echo "### DEXSEQGFF: ${DEXSEQGFF}"
echo "### RNABAM: ${RNABAM}"
echo "### DEXSEQOUTDIR: ${DEXSEQOUTDIR}"
echo "### DEXSEQOUT: ${DEXSEQOUT}"
echo "### DEXSEQCOUNTOUT: ${DEXSEQCOUNTOUT}"


perf stat python ${DEXSEQCOUNTPATH} -p yes -a 10 -f bam -r pos -s no ${DEXSEQGFF} ${RNABAM} ${DEXSEQOUT}.dexSeqCounts.txt 2> ${DEXSEQOUT}.DEXseqCount.perfOut
if [ $? -eq 0 ] ; then
	touch ${RUNDIR}/${NXT1}
	touch ${DEXSEQOUT}.DEXseqCountPass
else
	touch ${DEXSEQOUT}.DEXseqCountFail
fi

rm ${DEXSEQOUT}.DEXseqCountInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:DEXSEQCOUNT:$hours:$mins" > ${DEXSEQOUT}.DEXseqCount.totalTime
echo "### Ending DEXseqCount."
