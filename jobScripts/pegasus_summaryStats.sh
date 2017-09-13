#!/usr/bin/env bash
#SBATCH --job-name="pegasus_Summarystats"
#SBATCH --time=0-12:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### RUNDIR: ${RUNDIR}"
echo "### SUMSTATSPATH: ${SUMSTATSPATH}"
echo "### NXT1: ${NXT1}"

echo "Summary stats starting for ${RUNDIR} at $time."
java -jar ${SUMSTATSPATH}/parseDirforStats_JK.jar \
    ${RUNDIR} \
    ${SUMSTATSPATH}/Stats_v3.txt \
    ${SUMSTATSPATH}/TGENAssayCodes.txt > ${RUNDIR}/summaryStatsOut

if [ $? -eq 0 ] ; then
    mv ${RUNDIR}/summaryStatsOut ${RUNDIR}/summaryStatsPass
    touch ${RUNDIR}/${NXT1}
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
