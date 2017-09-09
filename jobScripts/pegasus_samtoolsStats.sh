#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_samStats"
#SBATCH --time=0-32:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"
 
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### RUNDIR: ${RUNDIR}"
echo "### BAMFILE: ${BAMFILE}"
echo "### SAMTOOLSPATH: ${SAMTOOLSPATH}"

cd ${DIR}
echo "### Starting samtools stats"
perf stat ${SAMTOOLSPATH}/samtools idxstats ${BAMFILE} > ${BAMFILE}.idxStats 2> ${BAMFILE}.idxStats.perfOut
perf stat ${SAMTOOLSPATH}/samtools flagstat ${BAMFILE} > ${BAMFILE}.samStats 2> ${BAMFILE}.samStats.perfOut
if [ $? -eq 0 ] ; then
	touch ${BAMFILE}.samtoolsStatsPass
else
	touch ${BAMFILE}.samtoolsStatsFail
fi
rm -f ${BAMFILE}.samtoolsStatsInQueue
#a little organizing
if [ -d ${RUNDIR}/stats/ ] ; then
	echo "moving files into stats folder"
	mv ${BAMFILE}.idxStats ${RUNDIR}/stats/	
	mv ${BAMFILE}.samStats ${RUNDIR}/stats/	
fi
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:SAMSTATS:$hours:$mins" > ${BAMFILE}.samStats.totalTime
echo "ending samtools stats"
