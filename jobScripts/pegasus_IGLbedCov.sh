#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_IGLbc"
#SBATCH --time=0-32:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"
 
module load R/2.15.2
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### DIR: ${DIR}"
echo "### IGLLISTBED: ${IGLLISTBED}"
echo "### SAMTOOLS: ${SAMTOOLSPATH}"
echo "### BEDTOOLS: ${BEDTOOLSPATH}"
echo "### BAMFILE: ${BAMFILE}"
echo "### PICARDPATH: ${PICARDPATH}"

echo "### Starting igl bed cov"
#perf stat ${SAMTOOLSPATH}/samtools bedcov ${IGLLISTBED} ${BAMFILE} > ${BAMFILE}.IGLbedCovOut 2> ${BAMFILE}.iglBedCov.perfOut
#perf stat ${BEDTOOLSPATH}/samtools bedcov ${IGLLISTBED} ${CHRLISTBED} ${BAMFILE} > ${BAMFILE}.IGLbedCovOut 2> ${BAMFILE}.iglBedCov.perfOut
#${SAMTOOLSPATH}/samtools view -b ${BAMFILE} | ${BEDTOOLSPATH}/bedtools coverage -split -counts -abam stdin -b ${CHRLISTBED} > ${BAMFILE}.bedToolsChrCount.txt
${SAMTOOLSPATH}/samtools view -b ${BAMFILE} | ${BEDTOOLSPATH}/bedtools coverage -split -counts -abam stdin -b ${IGLLISTBED} > ${BAMFILE}.bedToolsIglCount.txt
if [ $? -eq 0 ] ; then
	touch ${BAMFILE}.IGLbedCovPass
else
	touch ${BAMFILE}.IGLbedCovFail
fi
rm -f ${BAMFILE}.IGLbedCovInQueue
#a little organizing
if [ -d ${RUNDIR}/stats/ ] ; then
	echo "### Moving files into stats folder"
	#mv ${BAMFILE}.bedToolsChrCount.txt ${RUNDIR}/stats/
	mv ${BAMFILE}.bedToolsIglCount.txt ${RUNDIR}/stats/
fi
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:IGLBEDCOV:$hours:$mins" > ${BAMFILE}.iglBedCov.totalTime
echo "### Ending igl bed cov"
