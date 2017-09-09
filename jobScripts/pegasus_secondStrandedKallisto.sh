#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_kallisto"
#SBATCH --time=0-60:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
##PBS -e /${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err

#updated module 12/8/15
module load kallisto/0.42.4
#module load kallisto/0.42.3
#module load kallisto/0.43.0

beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### DIR: ${DIR}"
echo "### NXT1: ${NXT1}"
echo "### RUNDIR: ${RUNDIR}"
echo "### FASTQL: ${FASTQL}"
echo "### KALLISTO_INDEX_CDNA: ${KALLISTO_INDEX_CDNA}"
echo "### KALLISTO_INDEX_GTF: ${KALLISTO_INDEX_GTF}"

FASTQL=`echo "${FASTQL}" | tr ',' ' '`
echo "FASTQL: ${FASTQL}"
echo "### Starting Kallisto quant cDNA"

perf stat kallisto quant \
	--index=${KALLISTO_INDEX_CDNA} \
	--output-dir=${DIR}/cDNA \
	--bias \
	--rf-stranded \
	--bootstrap-samples=100 \
	--threads=14 \
	--seed=42 \
	${FASTQL} 2> ${DIR}.kallistoCDNA.perfOut
	if [ $? -eq 0 ] ; then
		echo "kallisto cDNA Passed"
	else
        	echo "kallisto cDNA Failed"
        	touch ${DIR}.kallistoFail
	fi


perf stat kallisto quant \
	--index=${KALLISTO_INDEX_GTF} \
	--output-dir=${DIR}/GTF \
	--bias \
	--bootstrap-samples=100 \
	--rf-stranded \
	--threads=14 \
	--seed=42 \
	${FASTQL} 2> ${DIR}.kallistoGTF.perfOut
	if [ $? -eq 0 ] ; then
		echo "kallisto GTF Passed"
		touch ${DIR}.kallistoPass
	else
		echo "kallisto GTF Failed"
		touch ${DIR}.kallistoFail
	fi

rm ${DIR}.kallistoInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:KALLISTOquant:$hours:$mins" > ${DIR}.kallistoQuant.totalTime
echo "### Ending Kallisto."
