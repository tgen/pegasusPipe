##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_SSpicRNAMetrics"
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
echo "### REFFLAT: ${REFFLAT}"
echo "### RIBINTS: ${RIBINTS}"
echo "### BAMFILE: ${BAMFILE}"
echo "### PICARDPATH: ${PICARDPATH}"

echo "### Starting first stranded picard rna metrics"
perf stat java -Xmx15g -Djava.io.tmpdir=/scratch/tgenjetstream/tmp/ -jar ${PICARDPATH}/picard.jar CollectRnaSeqMetrics \
	REF_FLAT=${REFFLAT} \
	REFERENCE_SEQUENCE=${REF} \
	RIBOSOMAL_INTERVALS=${RIBINTS} \
	STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
	INPUT=${BAMFILE} \
	OUTPUT=${BAMFILE}.picRNAMetrics \
	CHART_OUTPUT=${BAMFILE}.picRNAMetrics.pdf \
	TMP_DIR=/scratch/tgenjetstream/tmp/ \
	VALIDATION_STRINGENCY=SILENT > ${BAMFILE}.picRNAMetricsOut 2> ${BAMFILE}.picardRNA.perfOut
if [ $? -eq 0 ] ; then
	mv ${BAMFILE}.picRNAMetricsOut ${BAMFILE}.picRNAMetricsPass
else
	mv ${BAMFILE}.picRNAMetricsOut ${BAMFILE}.picRNAMetricsFail
fi
rm -f ${BAMFILE}.picRNAMetricsInQueue
#a little organizing
if [ -d ${RUNDIR}/stats/ ] ; then
	echo "moving files into stats folder"
	mv ${BAMFILE}.picRNAMetrics ${RUNDIR}/stats/
	mv ${BAMFILE}.picRNAMetrics.pdf ${RUNDIR}/stats/
fi
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:PICRNAMET:$hours:$mins" > ${BAMFILE}.picRnaMet.totalTime
echo "ending picard rna metrics"
