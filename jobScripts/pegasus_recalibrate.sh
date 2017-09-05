##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_RC"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 14
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"

time=`date +%d-%m-%Y-%H-%M` 
beginTime=`date +%s`
machine=`hostname`

echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### GATK: ${GATKPATH}"
echo "### KNOWN: ${KNOWN}"
echo "gatk base recalibration started on $machine"

perf stat java -Xmx44g -jar ${GATKPATH}/GenomeAnalysisTK.jar \
		-T BaseRecalibrator \
		-nct 8 \
		-l INFO \
		-R ${REF} \
		-knownSites ${KNOWN} \
		-I ${BAMFILE} \
		-cov ReadGroupCovariate \
		-cov QualityScoreCovariate \
		-cov CycleCovariate \
		-cov ContextCovariate \
		--disable_indel_quals \
		-o ${BAMFILE}.recal_data.grp 2> ${BAMFILE}.baseRecal.perfOut > ${BAMFILE}.recalibrateOut
if [ $? -ne 0 ] ; then
	mv ${BAMFILE}.recalibrateOut ${BAMFILE}.recalibrateFail
	echo "recal failed at base recalibrator"
	rm -rf ${BAMFILE}.recalibrateInQueue
	exit
fi
echo "gatk base recalibration print reads stage started"
perf stat java -Xmx44g -jar ${GATKPATH}/GenomeAnalysisTK.jar \
		-l INFO \
		-nct 14 \
		-R ${REF} \
		-I ${BAMFILE} \
		-T PrintReads \
		--out ${RECALBAM} \
		--disable_indel_quals \
		-BQSR ${BAMFILE}.recal_data.grp 2> ${BAMFILE}.recalPrint.perfOut >> ${BAMFILE}.recalibrateOut
if [ $? -eq 0 ] ; then
	mv ${BAMFILE}.recalibrateOut ${BAMFILE}.recalibratePass
	echo "Automatically removed by recalibration step to save on space" > ${BAMFILE}
	touch ${RUNDIR}/${NXT1}
else
	mv ${BAMFILE}.recalibrateOut ${BAMFILE}.recalibrateFail
fi
rm -rf ${BAMFILE}.recalibrateInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:RECAL:$hours:$mins" > ${BAMFILE}.rc.totalTime
time=`date +%d-%m-%Y-%H-%M` 
echo "gatk recalibration ended"
