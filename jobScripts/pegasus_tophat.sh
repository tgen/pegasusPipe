##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_tophat"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.out"
#SBATCH --error="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.err"

echo "### Variables coming in:"
echo "### SAMPLE=${SAMPLE}"
echo "### REF=${REF}"
echo "### FASTQ1=${FASTQ1}"
echo "### FASTQ2=${FASTQ2}"
echo "### DIR=${DIR}"
echo "### INDEXBASE=${INDEXBASE}"
echo "### USEGTF=${USEGTF}"
echo "### TRANSINDEX=${TRANSINDEX}"
echo "### RGID: ${RGID}"
echo "### TOPHAT2PATH: ${TOPHAT2PATH}"
echo "### PICARDPATH: ${PICARDPATH}"
echo "### SAMTOOLSPATH: ${SAMTOOLSPATH}"
echo "### BWAPATH: ${BWAPATH}"
echo "### BOWTIE2PATH: ${BOWTIE2PATH}"
echo "### DNAALIGNER: ${DNAALIGNER}"
export PATH=${BOWTIE2PATH}:$PATH

base=`basename ${FASTQ1}`
anotherName=${base/.R1.fastq.gz}
tempBamPrefix=${FASTQ1/.R1.fastq.gz}
trimmedRead1=${FASTQ1/.R1.fastq.gz/.trimmed.R1.fastq.gz}
trimmedRead2=${FASTQ2/.R2.fastq.gz/.trimmed.R2.fastq.gz}
fastq1tmp=${FASTQ1/.fastq.gz/.TEMP.fastq}
fastq2tmp=${FASTQ2/.fastq.gz/.TEMP.fastq}
cd ${DIR}

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`

echo "Start of getting insert size section"

#	echo "calculating combined read length"
#	read1Length=`gunzip -c ${FASTQ1} | head -2 | tail -1 | wc -c`
#	read2Length=`gunzip -c ${FASTQ2} | head -2 | tail -1 | wc -c`
#	combinedLength=`echo "$read1Length + $read2Length - 2" | bc`
#	echo "combined read length is $combinedLength"

#	echo "calculating allowed mismatches"
#	mmAllowed=`echo "$read1Length/25" | bc`
#	echo "calculated mismatches to be: $mmAllowed"

echo "getting first 2 million reads"
gunzip -c ${FASTQ1} | head -n 8000000 > $fastq1tmp &
gunzip -c ${FASTQ2} | head -n 8000000 > $fastq2tmp &
wait
f1size=`stat -c%s $fastq1tmp`
f2size=`stat -c%s $fastq2tmp`
echo "created temp fastq files with size: $f1size, $f2size"

if [ "${DNAALIGNER}" == "bwamem" ] ; then
	echo "aligning with BWA mem"
	${BWAPATH}/bwa mem -M -t8 ${REFPRETOPHAT} $fastq1tmp $fastq2tmp | ${SAMTOOLSPATH}/samtools view -S -h -b -t ${FAI} - | ${SAMTOOLSPATH}/samtools sort - $tempBamPrefix
	echo "BWA mem is done"
else
	echo "aligning with BWA regular"
	${BWAPATH}/bwa aln -t6 ${REF} $fastq1tmp > $fastq1tmp.sai &
	${BWAPATH}/bwa aln -t6 ${REF} $fastq2tmp > $fastq2tmp.sai &
	wait
	s1size=`stat -c%s $fastq1tmp.sai`
	s2size=`stat -c%s $fastq2tmp.sai`
	echo "aligning with BWA is done: created sai files with size: $s1size, $s2size"

	echo "starting bwa sampe"
	${BWAPATH}/bwa sampe -P ${REF} $fastq1tmp.sai $fastq2tmp.sai ${FASTQ1} ${FASTQ2} | ${SAMTOOLSPATH}/samtools view -S -h -b - | ${SAMTOOLSPATH}/samtools sort - $tempBamPrefix
	echo "bwa sampe ended"
fi

bsize=`stat -c%s $tempBamPrefix.bam`
echo "bam size: $bsize"

echo "getting insert size metrics with picard"
module load R/2.14.1
java -jar ${PICARDPATH}/picard.jar CollectInsertSizeMetrics INPUT=$tempBamPrefix.bam OUTPUT=$tempBamPrefix.bam.picInsertMetrics.txt HISTOGRAM_FILE=$tempBamPrefix.bam.picInsertMetrics.pdf VALIDATION_STRINGENCY=SILENT TMP_DIR=${DIR} LEVEL=ALL_READS 2>&1
echo "End of getting insert size section"

echo "calculating combined read length"
read1Length=`gunzip -c ${FASTQ1} | head -2 | tail -1 | wc -c`
read2Length=`gunzip -c ${FASTQ2} | head -2 | tail -1 | wc -c`
combinedLength=`echo "$read1Length + $read2Length - 2" | bc`
echo "combined read length is $combinedLength"

echo "calculating allowed mismatches"
mmAllowed=`echo "$read1Length/25 + 1" | bc`
echo "calculated mismatches to be: $mmAllowed"

echo "reading insert size from picard output"
INNERDIST=`head -n 8 $tempBamPrefix.bam.picInsertMetrics.txt | tail -n 1 | cut -f5 | awk -v var1="$combinedLength" '{printf "%.0f\n", $1-var1}'`
STDEV=`head -n 8 $tempBamPrefix.bam.picInsertMetrics.txt | tail -n 1 | cut -f6 | awk '{printf "%.0f\n", $1}'`
echo "COMBINED READ LENGTH=$combinedLength (from $read1Length and $read2Length - 2)"
echo "got INNERDIST=${INNERDIST} and STDEV=${STDEV}"

echo "End of getting insert size section"

ISODATE=`date --iso-8601`
echo "TIME:$time starting tophat on ${FASTQ1}"
if [ "${USEGTF}" == "no" ] ; then
	echo "tophat with no gtf..."
	#--rg-id ${RGID} --rg-date ${ISODATE} --rg-sample ${SAMPLE} --rg-platform ILLUMINA 
	perf stat ${TOPHAT2PATH}/tophat2 --keep-fasta-order -p 8 -r ${INNERDIST} --mate-std-dev ${STDEV}  -N $mmAllowed --read-edit-dist $mmAllowed --read-gap-length ${mmAllowed} --max-insertion-length 3 --max-deletion-length 3 --b2-very-sensitive -o ${DIR} ${INDEXBASE} ${FASTQ1} ${FASTQ2} > ${DIR}.thOut 2> ${DIR}.tophat.perfOut

else
	echo "tophat with gtf"
	#--rg-id ${RGID} --rg-date ${ISODATE} --rg-sample ${SAMPLE} --rg-platform ILLUMINA 
	perf stat ${TOPHAT2PATH}/tophat2 --keep-fasta-order -p 8 -r ${INNERDIST} --mate-std-dev ${STDEV} -N $mmAllowed --read-edit-dist ${mmAllowed} --read-gap-length ${mmAllowed} --max-insertion-length 3 --max-deletion-length 3 --b2-very-sensitive -o ${DIR} --transcriptome-index=${TRANSINDEX} ${INDEXBASE} ${FASTQ1} ${FASTQ2} > ${DIR}.thOut 2> ${DIR}.tophat.perfOut
fi
if [ $? -eq 0 ] ; then
	mv ${DIR}.thOut ${DIR}.thPass	
	echo "success."
	echo "renaming..."
	mv ${DIR}/accepted_hits.bam ${DIR}/$anotherName.accepted_hits.bam
	mv ${DIR}/unmapped.bam ${DIR}/$anotherName.unmapped.bam
	mv ${DIR}/junctions.bed ${DIR}/$anotherName.junctions.bed
	mv ${DIR}/insertions.bed ${DIR}/$anotherName.insertions.bed
	mv ${DIR}/deletions.bed ${DIR}/$anotherName.deletions.bed
	echo "renaming done"
	echo "Now making bam index and flagstat for ${DIR}/accepted_hits.bam"
	${SAMTOOLSPATH}/samtools index ${DIR}/$anotherName.accepted_hits.bam
	#${SAMTOOLSPATH}/samtools flagstat ${DIR}/$anotherName.accepted_hits.bam > ${DIR}/$anotherName.accepted_hits.bam.samStats
	echo "bam indexing and flagstat finished"

	echo "clean up, removing temp bam"
	rm -f $tempBamPrefix.bam

	touch ${RUNDIR}/${NXT1}
	touch ${RUNDIR}/${NXT2}
	touch ${RUNDIR}/${NXT3}
	touch ${RUNDIR}/${NXT4}
	touch ${RUNDIR}/${NXT6}
	touch ${RUNDIR}/${NXT7}
	touch ${RUNDIR}/${NXT8}
	touch ${RUNDIR}/${NXT9}
	touch ${RUNDIR}/${NXT10}
else
	mv ${DIR}.thOut ${DIR}.thFail
fi
rm -f ${DIR}.thInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))

(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:TOPHAT:$hours:$mins" > ${DIR}/$anotherName.accepted_hits.bam.tophat.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "TIME:$time tophat finished"
