#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_picHSMetrics"
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
echo "### BAMFILE: ${BAMFILE}"
echo "### PICARDPATH: ${PICARDPATH}"
echo "### BAITS: ${BAITS}"
echo "### TARGETS: ${TARGETS}"

cd ${DIR}
echo "### Starting picard rna metrics"
perf stat java -Xmx15g -jar ${PICARDPATH}/picard.jar CalculateHsMetrics \
	REFERENCE_SEQUENCE=${REF} \
	BAIT_INTERVALS=${BAITS} \
	TARGET_INTERVALS=${TARGETS} \
	INPUT=${BAMFILE} \
	OUTPUT=${BAMFILE}.picHSMetrics \
	PER_TARGET_COVERAGE=${BAMFILE}.picStats.HsPerTargetCov \
	TMP_DIR=/scratch/tgenjetstream/tmp/ \
	VALIDATION_STRINGENCY=SILENT > ${BAMFILE}.picHSMetricsOut 2> ${BAMFILE}.picardHS.perfOut
if [ $? -eq 0 ] ; then
	mv ${BAMFILE}.picHSMetricsOut ${BAMFILE}.picHSMetricsPass
else
	mv ${BAMFILE}.picHSMetricsOut ${BAMFILE}.picHSMetricsFail
fi
rm -f ${BAMFILE}.picHSMetricsInQueue
#a little organizing
if [ -d ${RUNDIR}/stats/ ] ; then
	echo "moving files into stats folder"
	mv ${BAMFILE}.picHSMetrics ${RUNDIR}/stats/	
	mv ${BAMFILE}.picHSMetrics.pdf ${RUNDIR}/stats/	
	mv ${BAMFILE}.picStats.HsPerTargetCov ${RUNDIR}/stats/
fi
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:PICHSMET:$hours:$mins" > ${BAMFILE}.picHSMet.totalTime
echo "ending picard HS metrics"
