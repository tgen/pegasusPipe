#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_FSthFusion"
#SBATCH --time=0-240:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"

module load python/2.7.3
 
echo "### Variables coming in:"
echo "### SAMPLE=${SAMPLE}"
echo "### REF=${REF}"
echo "### FAI=${FAI}"
echo "### FASTQ1=${FASTQ1}"
echo "### FASTQ2=${FASTQ2}"
echo "### DIR=${DIR}"
echo "### INDEXBASE=${INDEXBASE}"
echo "### TOPHAT2PATH: ${TOPHAT2PATH}"
echo "### THFUSIONPATH: ${THFUSIONPATH}"
echo "### SAMTOOLSPATH: ${SAMTOOLSPATH}"
echo "### TRIMFASTQPATH: ${TRIMFASTQPATH}"
echo "### REFPRETOPHAT: ${REFPRETOPHAT}"
echo "### PICARDPATH: ${PICARDPATH}"
echo "### BWAPATH: ${BWAPATH}"
export PATH=${BOWTIE1PATH}:$PATH

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`

echo "TIME:$time starting tophat fusion on ${FASTQ1}"
base=`basename ${FASTQ1}`
anotherName=${base/.R1.fastq.gz}
anotherName2=${base/.proj.R1.fastq.gz}

fastq1tmp=${FASTQ1/.fastq.gz/.TEMP.fastq}
fastq2tmp=${FASTQ2/.fastq.gz/.TEMP.fastq}
tempBamPrefix=${FASTQ1/.R1.fastq.gz}

echo "fastq1tmp = ${fastq1tmp} fastq2tmp = ${fastq2tmp} tempBamPrefix = ${tempBamPrefix}"
#trimmedRead1=${FASTQ1/.R1.fastq.gz/.50mer.R1.fastq.gz}
#trimmedRead2=${FASTQ2/.R2.fastq.gz/.50mer.R2.fastq.gz}

cd ${DIR}
#echo "Starting trim section"
#echo "Trim values: ${PRE1},${POS1},${PRE2},${POS2}"
#r1Trimmed=0
#r2Trimmed=0
#if [[ ${PRE1} != 0 || ${POS1} != 0 ]] ; then
#	r1Trimmed=1
#	echo "### we need trimming first on ${FASTQ1}"
#	/home/tgenjetstream/bin/DPtrimFastq.pl ${FASTQ1} ${PRE1} ${POS1} | gzip > $trimmedRead1 &
#	f1id=$! #throw in background capture pid
#else
#	echo "### no trimming needed on ${FASTQ1}"
#	ln -s ${FASTQ1} $trimmedRead1
#	touch ${FASTQ1}.trimPass
#fi
#if [[ ${PRE2} != 0 || ${POS2} != 0 ]] ; then
#	r2Trimmed=1
#	echo "### we need trimming first on ${FASTQ2}"
#	/home/tgenjetstream/bin/DPtrimFastq.pl ${FASTQ2} ${PRE2} ${POS2} | gzip > $trimmedRead2
#	if [ $? -eq 0 ] ; then
#		echo "### trim finished OK"
		#juggling files
		#mv ${FASTQ2} ${FASTQ2}.untrimmed
		#mv ${FASTQ2}.trimmed ${FASTQ2}
#		touch ${FASTQ2}.trimPass
#	else
#		echo "### trim failed"
#		touch ${FASTQ2}.trimFail
#	fi
#else
#	echo "### no trimming needed on ${FASTQ2}"
#	ln -s ${FASTQ2} $trimmedRead2
#	touch ${FASTQ2}.trimPass
#fi

#if [[ ${PRE1} != 0 || ${POS1} != 0 ]] ; then
#	wait $f1id #better wait for  trimming job on read 1 to end
#	if [ $? -eq 0 ] ; then
#		echo "### trim finished OK"
		#juggling files
		#mv ${FASTQ1} ${FASTQ1}.untrimmed
		#mv ${FASTQ1}.trimmed ${FASTQ1}
#		touch ${FASTQ1}.trimPass
#	else
#		echo "### trim failed"
#		touch ${FASTQ1}.trimFail
#	fi
#fi

#now set names
#if [ $r1Trimmed -eq 1 ] ; then
#	FASTQ1=$trimmedRead1
#fi
#if [ $r2Trimmed -eq 1 ] ; then
#	FASTQ2=$trimmedRead2
#fi
#end of setting names
#
#echo "Starting reducing fastq to 60 million reads"
#gunzip -c ${FASTQ1} | head -240000000 | gzip > ${FASTQ1}.60mil.gz
#gunzip -c ${FASTQ2} | head -240000000 | gzip > ${FASTQ2}.60mil.gz
#echo "End of reducing fastq to 60 million reads"

echo "### Start of getting insert size section"
#
echo "### Getting first 2 million reads"
gunzip -c ${FASTQ1} | head -n 8000000 > $fastq1tmp &
gunzip -c ${FASTQ2} | head -n 8000000 > $fastq2tmp &
wait
f1size=`stat -c%s $fastq1tmp`
f2size=`stat -c%s $fastq2tmp`
echo "### Created temp fastq files with size: $f1size, $f2size"
#
echo "### Aligning with BWA mem"
echo "### Command: ${BWAPATH}/bwa mem -t8 ${REFPRETOPHAT} $fastq1tmp $fastq2tmp | ${SAMTOOLSPATH}/samtools view -S -h -b -t ${FAI} - | ${SAMTOOLSPATH}/samtools sort - $tempBamPrefix"
${BWAPATH}/bwa mem -t8 ${REFPRETOPHAT} $fastq1tmp $fastq2tmp | ${SAMTOOLSPATH}/samtools view -S -h -b -t ${FAI} - | ${SAMTOOLSPATH}/samtools sort - $tempBamPrefix
echo "### BWA mem is done"
bsize=`stat -c%s $tempBamPrefix.bam`
echo "### Bam size: $bsize"
#
echo "### Getting insert size metrics with picard"
module load R/2.14.1
statsOutName=${tempBamPrefix/.proj}
java -jar ${PICARDPATH}/picard.jar CollectInsertSizeMetrics INPUT=$tempBamPrefix.bam OUTPUT=$statsOutName.bwa.transcriptome.picInsertMetrics.txt HISTOGRAM_FILE=$statsOutName.bwa.transcriptome.picInsertMetrics.pdf VALIDATION_STRINGENCY=SILENT TMP_DIR=${DIR} LEVEL=ALL_READS 2>&1
echo "### End of getting insert size section"
#
echo "### Calculating combined read length"
read1Length=`gunzip -c ${FASTQ1} | head -2 | tail -1 | wc -c`
read2Length=`gunzip -c ${FASTQ2} | head -2 | tail -1 | wc -c`
combinedLength=`echo "$read1Length + $read2Length - 2" | bc`
echo "### Combined read length is $combinedLength"
#
##echo "### Calculating allowed mismatches"
#mmAllowed=`echo "$read1Length/25 + 1" | bc`
#echo "calculated mismatches to be: $mmAllowed"

echo "### Reading insert size from picard output"
#INNERDIST=`head -n 8 $tempBamPrefix.bam.picInsertMetrics.txt | tail -n 1 | cut -f5 | awk -v var1="$combinedLength" '{printf "%.0f\n", $1-var1}'`
#STDEV=`head -n 8 $tempBamPrefix.bam.picInsertMetrics.txt | tail -n 1 | cut -f6 | awk '{printf "%.0f\n", $1}'`
INNERDIST=`head -n 8 $statsOutName.bwa.transcriptome.picInsertMetrics.txt | tail -n 1 | cut -f5 | awk -v var1="$combinedLength" '{printf "%.0f\n", $1-var1}'`
STDEV=`head -n 8 $statsOutName.bwa.transcriptome.picInsertMetrics.txt | tail -n 1 | cut -f6 | awk '{printf "%.0f\n", $1}'`

echo "### COMBINED READ LENGTH=$combinedLength (from $read1Length and $read2Length - 2)"
echo "### Got INNERDIST=${INNERDIST} and STDEV=${STDEV}"
#
echo "###End of getting insert size section"

#setting path
PATH=${BOWTIE1PATH}:$PATH
	#--b2-very-fast \
perf stat ${THFUSIONPATH}/tophat2 \
	-p 16 \
	-N 3 \
	--library-type fr-secondstrand \
	--read-edit-dist 3 \
	--fusion-search \
	--bowtie1 \
	--no-coverage-search \
	-r ${INNERDIST} \
	--mate-std-dev ${STDEV} \
	--fusion-min-dist 100000 \
	--max-intron-length 100000 \
	--keep-fasta-order \
	--fusion-anchor-length 20 \
	--fusion-ignore-chromosomes MT \
	-o ${DIR} \
	${INDEXBASE} ${FASTQ1} ${FASTQ2} > ${DIR}.thFusionOut 2> ${DIR}.thFusion.perfOut
if [ $? -eq 0 ] ; then
	echo "success."
	echo "renaming..."
	mv ${DIR}/accepted_hits.bam ${DIR}/$anotherName.accepted_hits.bam
	mv ${DIR}/unmapped.bam ${DIR}/$anotherName.unmapped.bam
	mv ${DIR}/junctions.bed ${DIR}/$anotherName.junctions.bed
	mv ${DIR}/insertions.bed ${DIR}/$anotherName.insertions.bed
	mv ${DIR}/deletions.bed ${DIR}/$anotherName.deletions.bed
	mv ${DIR}/fusions.out ${DIR}/$anotherName.fusions.out
	mv ${DIR}/$tempBamPrefix.bam ${DIR}/$tempBamPrefix.bam.2MilReads
	echo "renaming done"
	echo "Now making bam index and flagstat for ${DIR}/accepted_hits.bam"
	${SAMTOOLSPATH}/samtools index ${DIR}/$anotherName.accepted_hits.bam
	${SAMTOOLSPATH}/samtools flagstat ${DIR}/$anotherName.accepted_hits.bam > ${DIR}/$anotherName.accepted_hits.bam.samStats
	echo "bam indexing and flagstat finished"


	mv ${DIR}.thFusionOut ${DIR}.thFusionPass	
	touch ${RUNDIR}/${NXT1}
else
	mv ${DIR}.thFusionOut ${DIR}.thFusionFail
fi
rm -f ${DIR}.thFusionInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:THFUSION:$hours:$mins" > ${DIR}/$anotherName.accepted_hits.bam.thFusion.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "TIME:$time tophat fusion finished"

