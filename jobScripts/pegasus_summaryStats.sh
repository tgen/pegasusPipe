#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_Summarystats"
#SBATCH --time=0-12:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### RUNDIR: ${RUNDIR}"
echo "### SUMSTATSPATH: ${SUMSTATSPATH}"
echo "### NXT1: ${NXT1}"

echo "Summary stats starting for ${RUNDIR} at $time."
perf stat java -jar ${SUMSTATSPATH}/parseDirforStats_JK.jar ${RUNDIR} ${SUMSTATSPATH}/Stats_v3.txt ${SUMSTATSPATH}/TGENAssayCodes.txt > ${RUNDIR}/summaryStatsOut 2> ${RUNDIR}/summaryStat.perfOut
if [ $? -eq 0 ] ; then
	mv ${RUNDIR}/summaryStatsOut ${RUNDIR}/summaryStatsPass
	touch ${RUNDIR}/${NXT1}

	#cp ${RUNDIR}/Summary_allStats.txt ${TARGETDIR}/Summary_allStats.txt
	#echo "Summary of all stats is attached for project at ${TARGETDIR}." >> ${TARGETDIR}/mini_Summary_allStats.txt
	#grep MEAN_TARGET ${RUNDIR}/Summary_allStats.txt >> ${TARGETDIR}/mini_Summary_allStats.txt
	#grep PCT_TARGET_BASES_20X ${RUNDIR}/Summary_allStats.txt >> ${TARGETDIR}/mini_Summary_allStats.txt

	#grep samStats.reads.mapped ${RUNDIR}/Summary_allStats.txt >> ${TARGETDIR}/mini_Summary_allStats.txt
	#grep PCT_MRNA_BASES ${RUNDIR}/Summary_allStats.txt >> ${TARGETDIR}/mini_Summary_allStats.txt
	#grep MEDIAN_5PRIME_BIAS ${RUNDIR}/Summary_allStats.txt >> ${TARGETDIR}/mini_Summary_allStats.txt
	#grep MEDIAN_3PRIME_BIAS ${RUNDIR}/Summary_allStats.txt >> ${TARGETDIR}/mini_Summary_allStats.txt
	#grep MEDIAN_5PRIME_TO_3PRIME_BIAS ${RUNDIR}/Summary_allStats.txt >> ${TARGETDIR}/mini_Summary_allStats.txt

	#grep samStats.reads.mapped ${RUNDIR}/Summary_allStats.txt >> ${TARGETDIR}/mini_Summary_allStats.txt
	#grep samStats.reads.properly ${RUNDIR}/Summary_allStats.txt >> ${TARGETDIR}/mini_Summary_allStats.txt
	#grep MEAN_INSERT_SIZE ${RUNDIR}/Summary_allStats.txt >> ${TARGETDIR}/mini_Summary_allStats.txt
	#grep STANDARD_DEVIATION ${RUNDIR}/Summary_allStats.txt >> ${TARGETDIR}/mini_Summary_allStats.txt
	#grep MEDIAN_INSERT_SIZE ${RUNDIR}/Summary_allStats.txt >> ${TARGETDIR}/mini_Summary_allStats.txt
	#grep MEDIAN_ABSOLUTE_DEVIATION ${RUNDIR}/Summary_allStats.txt >> ${TARGETDIR}/mini_Summary_allStats.txt
	#this is where we grep things out of file above and make it the body of the email below.
	#cat ${TARGETDIR}/mini_Summary_allStats.txt | mail -s "pnap pipeline: stats for project ${TARGETDIR}" -a ${RUNDIR}/Summary_allStats.txt ${EMAIL}
	#cat ${TARGETDIR}/mini_Summary_allStats.txt | mail -s "pnap pipeline: stats for project ${TARGETDIR}" -a ${RUNDIR}/Summary_allStats.txt tgenjetstream@tgen.org
	#cat ${TARGETDIR}/mini_Summary_allStats.txt | mail -s "pnap pipeline: stats for project ${TARGETDIR}" -a ${RUNDIR}/Summary_allStats.txt snasser@tgen.org
	#cat ${TARGETDIR}/mini_Summary_allStats.txt | mail -s "pnap pipeline: stats for project ${TARGETDIR}" -a ${RUNDIR}/Summary_allStats.txt jaldrich@tgen.org
elif [ $? -eq 1 ] ; then
	mv ${RUNDIR}/summaryStatsOut ${RUNDIR}/summaryStatsWARNING
else
	mv ${RUNDIR}/summaryStatsOut ${RUNDIR}/summaryStatsFail
fi
rm -f ${RUNDIR}/summaryStatsInQueue
time=`date +%d-%m-%Y-%H-%M`
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:SUMMARYSTAT:$hours:$mins" > ${RUNDIR}/summaryStats.totalTime
echo "Summary stats finished at $time."
