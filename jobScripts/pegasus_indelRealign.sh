#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_IR"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"
#originalBam=${IRBAMFILE/%.md.ir.bam/.bam}
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### GATK: ${GATKPATH}"

echo "### Starting indel realigning of file ${BAMFILE}"
echo "### Starting step 1, target creator"
perf stat java -Xmx4g -Djava.io.tmpdir=/scratch/tgenjetstream/tmp/ \
	-jar ${GATKPATH}/GenomeAnalysisTK.jar \
	-I ${BAMFILE} \
	-R ${REF} \
	-T RealignerTargetCreator \
	-nt 16 \
	--maxIntervalSize 350 \
	-DBQ 1 \
	-o ${INTS} \
	--disable_auto_index_creation_and_locking_when_reading_rods \
	--known ${INDELS} 2> ${BAMFILE}.irStep1.perfOut > ${BAMFILE}.indelRealignOut
if [ $? -ne 0 ] ; then
	mv ${BAMFILE}.indelRealignOut ${BAMFILE}.indelRealignFail
	rm -f ${BAMFILE}.indelRealignInQueue
	exit
fi
echo "### Starting step 2, indel realignment"
perf stat java -Xmx44g -Djava.io.tmpdir=/scratch/tgenjetstream/tmp/ \
	-jar ${GATKPATH}/GenomeAnalysisTK.jar \
	-I ${BAMFILE} \
	-R ${REF} \
	-T IndelRealigner \
	-DBQ 1 \
	-targetIntervals ${INTS} \
	--maxReadsInMemory 5000000 \
	--maxConsensuses 24 \
	--maxReadsForConsensuses 80 \
	--maxReadsForRealignment 12000 \
	-o ${IRBAMFILE} \
	-model KNOWNS_ONLY \
	--disable_auto_index_creation_and_locking_when_reading_rods \
	-known ${INDELS} 2> ${BAMFILE}.irStep2.perfOut >> ${BAMFILE}.indelRealignOut
if [ $? -eq 0 ] ; then
	mv ${BAMFILE}.indelRealignOut ${BAMFILE}.indelRealignPass
	echo "Automatically removed by indel realignment step to save on space" > ${BAMFILE}
	touch ${RUNDIR}/${NXT1}
else
	mv ${BAMFILE}.indelRealignOut ${BAMFILE}.indelRealignFail
fi

rm -f ${BAMFILE}.indelRealignInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:GATKIR:$hours:$mins" > ${BAMFILE}.ir.totalTime
