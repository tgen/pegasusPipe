#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_DEXseq"
#SBATCH --time=0-60:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
##PBS -e /${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err

module load R/3.2.1

beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### DIR: ${DIR}"
echo "### NXT1: ${NXT1}"
echo "### RUNDIR: ${RUNDIR}"
echo "### DEXSEQPATH: ${DEXSEQPATH}"
echo "### DEXSEQCONFIG: ${DEXSEQCONFIG}"
echo "### DEXSEQOUTDIR: ${DEXSEQOUTDIR}"
echo "### OBJECTDATA: ${OBJECTDATA}"
echo "### DEXSEQOUTFILE: ${DEXSEQOUTFILE}"
echo "### KALLISTOOUT: ${KALLISTOOUT}"

base_Dir=${KALLISTOOUT}
configFile=${DEXSEQCONFIG}
outData=${DEXSEQOUTFILE}
objectData=${OBJECTDATA}

##use DEXSeq R script to perform differential exon usage analysis.

perf stat Rscript --vanilla ${DEXSEQPATH} \
	${DEXSEQCONFIG} \
	${DEXSEQGFF} \
	${DEXSEQOUTFILE} \
	${OBJECTDATA} 2> ${DEXSEQOUTDIR}.DEXseq.perfOut
	if [ $? -eq 0 ] ; then
		touch ${DEXSEQOUTDIR}.DEXseqPass
	else
		touch ${DEXSEQOUTDIR}.DEXseqFail
	fi


rm ${DEXSEQOUTDIR}.DEXseqInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:DEXSEQ:$hours:$mins" > ${DIR}.DEXseq.totalTime
echo "### Ending DEXseq."
