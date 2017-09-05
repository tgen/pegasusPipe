##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_markDups"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"
 
cd ${DIR}
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### PICARDPATH: ${PICARDPATH}"
echo "### SAMTOOLSPATH: ${SAMTOOLSPATH}"

echo "### Starting picard mark duplicates"
#echo "faking it" > ${BAMFILE}.mdOut
perf stat java -Xmx22g -jar ${PICARDPATH}/picard.jar MarkDuplicates ASSUME_SORTED=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT TMP_DIR=/scratch/tgenjetstream/tmp INPUT=${BAMFILE} OUTPUT=${OUTPUTBAM} METRICS_FILE=${BAMFILE}.picStats.MarkDupMetrics MAX_RECORDS_IN_RAM=18000000 CREATE_INDEX=true 2> ${BAMFILE}.markDups.perfOut > ${BAMFILE}.mdOut
if [ $? -eq 0 ] ; then
	mv ${BAMFILE}.mdOut ${BAMFILE}.mdPass
	echo "Automatically removed by mark duplicates step to save on space" > ${BAMFILE}
	#a little organizing
	if [ ! -d ${RUNDIR}/stats/ ] ; then
		mkdir -p ${RUNDIR}/stats
	fi
	mv ${BAMFILE}.picStats.MarkDupMetrics ${RUNDIR}/stats/	
	touch ${RUNDIR}/${NXT1}
	touch ${RUNDIR}/${NXT2}
	touch ${RUNDIR}/${NXT3}
	touch ${RUNDIR}/${NXT4}
	touch ${RUNDIR}/${NXT5}
	touch ${RUNDIR}/${NXT6}
	touch ${RUNDIR}/${NXT7}
	touch ${RUNDIR}/${NXT8}
	touch ${RUNDIR}/${NXT9}
	touch ${RUNDIR}/${NXT10}
else
	mv ${BAMFILE}.mdOut ${BAMFILE}.mdFail
fi
rm -f ${BAMFILE}.mdInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:MARKDUPS:$hours:$mins" > ${OUTPUTBAM}.md.totalTime
echo "### Ending picard mark duplicates"
