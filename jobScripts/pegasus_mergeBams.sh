#!/usr/bin/env bash
#SBATCH --job-name="pegasus_MergeBams"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL

cd ${RUNDIR}
time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### PICARDPATH: ${PICARDPATH}"
echo "### SAMTOOLSPATH: ${SAMTOOLSPATH}"
echo "### CNT: ${CNT}"
echo "### RUNDIR: ${RUNDIR}"
echo "### MERGEDBAM: ${MERGEDBAM}"
echo "### BAMLIST: ${BAMLIST}"
echo "### NXT1: ${NXT1}"

echo "### TIME:$time starting picard merge bams to create ${MERGEDBAM}"

onlyBamFile=${BAMLIST/I=/}
onlyBaiFile=${onlyBamFile/.bam/.bai}
mergedBai=${MERGEDBAM/.bam/.bai}

if [ ${CNT} -eq 1 ] ; then
    #nothing really merged, only copied
	echo "### Just copying $onlyBamFile to ${MERGEDBAM}" > ${MERGEDBAM}.mergeBamOut
	cp $onlyBaiFile $mergedBai
	cp $onlyBamFile ${MERGEDBAM}
	if [ $? -ne 0 ] ; then #bad cp
		mv ${MERGEDBAM}.mergeBamOut ${MERGEDBAM}.mergeBamFail
	else #good cp
		mv ${MERGEDBAM}.mergeBamOut ${MERGEDBAM}.mergeBamPass
		echo "Automatically removed by merge bam step to save on space" > $onlyBamFile
		touch ${RUNDIR}/${NXT1}
	fi
else
    #actually merged with picard
	java -Xmx42g -jar ${PICARDPATH}/picard.jar MergeSamFiles \
	    ASSUME_SORTED=true \
	    USE_THREADING=true \
	    VALIDATION_STRINGENCY=SILENT \
	    TMP_DIR=/scratch/tgenjetstream/tmp \
	    OUTPUT=${MERGEDBAM} \
	    ${BAMLIST} > ${MERGEDBAM}.mergeBamOut

	if [ $? -ne 0 ] ; then #bad merge
		mv ${MERGEDBAM}.mergeBamOut ${MERGEDBAM}.mergeBamFail
	else #good merge
		echo "### Starting indexing of bam with samtools now that merge finished OK"
		${SAMTOOLSPATH}/samtools index ${MERGEDBAM}
		mv ${MERGEDBAM}.bai $mergedBai
		echo "### Ended indexing of bam after merging."
		mv ${MERGEDBAM}.mergeBamOut ${MERGEDBAM}.mergeBamPass
		for bam in ${BAMLIST}
		do
			bamPath=`echo $bam | cut -d= -f2`
			echo "Automatically removed by merge bam step to save on space" > $bamPath
		done
		touch ${RUNDIR}/${NXT1}
	fi
fi

rm ${MERGEDBAM}.mergeBamInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:MERGEBAMS:$hours:$mins" > ${MERGEDBAM}.mergeBam.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "TIME:$time picard merge bams finished"
