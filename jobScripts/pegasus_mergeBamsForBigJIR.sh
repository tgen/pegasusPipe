#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_MergeForJIR"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"


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
echo "### NEWLOC: ${NEWLOC}"
echo "### BAMLIST: ${BAMLIST}"
echo "### NXT1: ${NXT1}"
echo "### NXT2: ${NXT2}"
echo "### NXT3: ${NXT3}"
echo "### NXT4: ${NXT4}"
echo "### NXT5: ${NXT5}"

echo "### TIME:$time starting picard merge bams to create ${MERGEDBAM}"

onlyBamFile=${BAMLIST/I=/}
onlyBaiFile=${onlyBamFile/.bam/.bai}
mergedBai=${MERGEDBAM/.bam/.bai}
newLocBai=${NEWLOC/.bam/.bai}

if [ ${CNT} -eq 1 ] ; then #nothing really merged, only copied
	echo "just copying $onlyBamFile to ${MERGEDBAM}" > ${MERGEDBAM}.mergeBamOut
	cp $onlyBaiFile $mergedBai
	perf stat cp $onlyBamFile ${MERGEDBAM} 2> ${MERGEDBAM}.mergeBam.perfOut
	if [ $? -ne 0 ] ; then #bad cp
		mv ${MERGEDBAM}.mergeBamOut ${MERGEDBAM}.mergeBamFail
	else #good cp
		mv ${MERGEDBAM}.mergeBamOut ${MERGEDBAM}.mergeBamPass
		echo "Automatically removed by merge bam step to save on space" > $onlyBamFile
		touch ${RUNDIR}/${NXT1}
	fi
else #actually merged with picard
	perf stat java -Xmx42g -jar ${PICARDPATH}/picard.jar MergeSamFiles ASSUME_SORTED=true USE_THREADING=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/scratch/tgenjetstream/tmp OUTPUT=${MERGEDBAM} ${BAMLIST} 2> ${MERGEDBAM}.mergeBam.perfOut > ${MERGEDBAM}.mergeBamOut
	if [ $? -ne 0 ] ; then #bad merge
		mv ${MERGEDBAM}.mergeBamOut ${MERGEDBAM}.mergeBamFail
	else #good merge
		echo "### Starting indexing of bam with samtools now that merge finished OK"
		perf stat ${SAMTOOLSPATH}/samtools index ${MERGEDBAM} 2> ${MERGEDBAM}.samindex.perfOut
		mv ${MERGEDBAM}.bai $mergedBai
		echo "### Ended indexing of bam after merging."
		mv ${MERGEDBAM}.mergeBamOut ${MERGEDBAM}.mergeBamPass
		echo "### Moving ${MERGEDBAM}"
		if [ -e ${NEWLOC} ] ; then
			echo "### ${NEWLOC} already exists on target, possibly from another joint IR"
		else
			echo "### ${NEWLOC} does not exist on target, copying now..."
			cp ${MERGEDBAM} ${NEWLOC}
			cp $mergedBai $newLocBai
			echo "### Moved out of here to its own dir at ${NEWLOC}" > ${MERGEDBAM}
			echo "### Moved out of here to its own dir at $newLocBai" > $mergedBai
			touch ${NEWLOC}.jointIRPass
		fi
		for bam in ${BAMLIST}
		do
			bamPath=`echo $bam | cut -d= -f2`
			echo "Automatically removed by merge bam step to save on space" > $bamPath
		done

		touch ${RUNDIR}/${NXT1}
		touch ${RUNDIR}/${NXT2}
		touch ${RUNDIR}/${NXT3}
		touch ${RUNDIR}/${NXT4}
		touch ${RUNDIR}/${NXT5}
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
