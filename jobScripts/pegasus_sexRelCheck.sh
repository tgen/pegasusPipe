#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_sexRelCheck"
#SBATCH --time=0-96:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### VCF: ${VCF}"
echo "### OUTVCF: ${OUTVCF}"
echo "### GREGORPATH: ${GREGORPATH}"
echo "### PLINK2PATH: ${PLINK2PATH}"
echo "### OUTTRACKNAME: ${OUTTRACKNAME}"

module load python/3.4

cd ${OUTTRACK}

perf stat ${GREGORPATH}/gregor.py --vcf ${VCF} --plink2 ${PLINK2PATH} 2> ${OUTTRACKNAME}.sexRelCheck.perfOut
if [ $? -eq 0 ] ; then
	touch ${OUTTRACKNAME}.sexRelCheckPass
	mv ${VCF}.gregor2.txt ${OUTTRACKNAME}.gregor2.txt
else	
	touch${OUTTRACKNAME}.sexRelCheckFail
fi

rm -f ${OUTTRACKNAME}.sexRelCheckInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:sexRelCheck:$hours:$mins" > ${OUTTRACKNAME}.sexRelCheck.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "Sex and Relationship Check finished at $time."
