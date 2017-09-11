#!/usr/bin/env bash
#SBATCH --job-name="pegasus_picMltMtrcs"
#SBATCH --time=0-32:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL

module load R/2.15.2
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### BAMFILE: ${BAMFILE}"
echo "### PICARDPATH: ${PICARDPATH}"

cd ${DIR}
echo "### Starting picard multi metrics"
java -Djava.io.tmpdir=/scratch/tgenjetstream/tmp/ -Xmx15g -jar ${PICARDPATH}/picard.jar CollectMultipleMetrics \
	INPUT=${BAMFILE} \
	REFERENCE_SEQUENCE=${REF} \
	PROGRAM=CollectInsertSizeMetrics \
	PROGRAM=CollectAlignmentSummaryMetrics \
	PROGRAM=QualityScoreDistribution \
	PROGRAM=MeanQualityByCycle \
	OUTPUT=${BAMFILE}.picMultiMetrics \
	TMP_DIR=/scratch/tgenjetstream/tmp/ \
	ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=SILENT > ${BAMFILE}.picMultiMetricsOut

if [ $? -eq 0 ] ; then
	mv ${BAMFILE}.picMultiMetricsOut ${BAMFILE}.picMultiMetricsPass
else
	mv ${BAMFILE}.picMultiMetricsOut ${BAMFILE}.picMultiMetricsFail
fi

rm -f ${BAMFILE}.picMultiMetricsInQueue

# A little organizing
if [ -d ${RUNDIR}/stats/ ] ; then
	echo "### Moving files into stats folder"
	mv ${BAMFILE}.picMultiMetrics ${RUNDIR}/stats/	
	mv ${BAMFILE}.picMultiMetrics.insert_size_Histogram.pdf ${BAMFILE}.picMultiMetrics.insert_size_histogram.pdf	
	mv ${BAMFILE}.picMultiMetrics.* ${RUNDIR}/stats/	
fi

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:PICMULTIMET:$hours:$mins" > ${BAMFILE}.picMultiMet.totalTime
echo "### Ending picard MULTI metrics"
