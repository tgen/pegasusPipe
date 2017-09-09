#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_phaseBT"
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
echo "### REF: ${REF}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### GATKPATH: ${GATKPATH}"
echo "### VCF: ${VCF}"
echo "### OUTVCF: ${OUTVCF}"
echo "### PED: ${PED}"
echo "### OUTTRACKNAME: ${OUTTRACKNAME}"

echo "### Phase by transmission started at $time."
perf stat java -Djava.io.tmpdir=/scratch/tgenjetstream/tmp/ -jar -Xmx32g ${GATKPATH}/GenomeAnalysisTK.jar \
-R ${REF} \
-T PhaseByTransmission \
-V ${VCF} \
-ped ${PED} \
-o ${OUTVCF} > ${OUTTRACKNAME}.phaseBTOut 2> ${OUTTRACKNAME}.phaseBT.perfOut
if [ $? -eq 0 ] ; then
	mv ${OUTTRACKNAME}.phaseBTOut ${OUTTRACKNAME}.phaseBTPass
	touch ${RUNDIR}/${NXT1}
else	
	mv ${OUTTRACKNAME}.phaseBTOut ${OUTTRACKNAME}.phaseBTFail
fi
rm -f ${OUTTRACKNAME}.phaseBTInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:GATKPHASEBT:$hours:$mins" > ${OUTTRACKNAME}.phaseBT.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "Phase by transmission finished at $time."
