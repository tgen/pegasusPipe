#!/usr/bin/env bash
#SBATCH --job-name="pegasus_samStats"
#SBATCH --time=0-32:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL

beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### RUNDIR: ${RUNDIR}"
echo "### BAMFILE: ${BAMFILE}"
echo "### SAMTOOLSPATH: ${SAMTOOLSPATH}"

cd ${DIR}
echo "### Starting samtools stats"
${SAMTOOLSPATH}/samtools idxstats ${BAMFILE} > ${BAMFILE}.idxStats
${SAMTOOLSPATH}/samtools flagstat ${BAMFILE} > ${BAMFILE}.samStats

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
