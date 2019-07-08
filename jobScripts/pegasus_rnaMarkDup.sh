#!/usr/bin/env bash
#SBATCH --job-name="pegasus_rnaMD"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH --mem 64G

cd ${DIR}
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### PICARDPATH: ${PICARDPATH}"
echo "### BAMFILE: ${BAMFILE}"
echo "### SAMTOOLSPATH: ${SAMTOOLSPATH}"

echo "### Starting picard mark duplicates"
java -Xmx22g -jar ${PICARDPATH}/picard.jar MarkDuplicates \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=false \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=/scratch/tgenjetstream/tmp \
    INPUT=${BAMFILE} OUTPUT=${OUTPUTBAM} \
    METRICS_FILE=${BAMFILE}.picStats.MarkDupMetrics \
    MAX_RECORDS_IN_RAM=18000000 \
    CREATE_INDEX=true > ${BAMFILE}.rnaMarkDupOut

if [ $? -eq 0 ] ; then
    mv ${BAMFILE}.rnaMarkDupOut ${BAMFILE}.rnaMarkDupPass
    # A little organizing
    if [ ! -d ${RUNDIR}/stats/ ] ; then
        mkdir -p ${RUNDIR}/stats
    fi
    ${SAMTOOLSPATH}/samtools idxstats ${OUTPUTBAM} > ${OUTPUTBAM}.idxStats
    mv ${BAMFILE}.picStats.MarkDupMetrics ${RUNDIR}/stats/
    mv ${OUTPUTBAM}.idxStats ${RUNDIR}/stats/
    touch ${RUNDIR}/${NXT1}
    touch ${RUNDIR}/${NXT2}
    touch ${RUNDIR}/${NXT3}
    touch ${RUNDIR}/${NXT4}
    touch ${RUNDIR}/${NXT5}
else
    mv ${BAMFILE}.rnaMarkDupOut ${BAMFILE}.rnaMarkDupFail
fi

rm -f ${BAMFILE}.rnaMarkDupInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:MARKDUPS:$hours:$mins" > ${OUTPUTBAM}.rnaMarkDup.totalTime
echo "### Ending picard mark duplicates"
