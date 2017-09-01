##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_strandedStar"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.out"
#SBATCH --error="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.err"
#
time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### STARREF=${STARREF}"
echo "### STARGTF=${STARGTF}"
echo "### RGTAGLIST=${RGTAGLIST}"
echo "### FASTQL1=${FASTQL1}"
echo "### FASTQL2=${FASTQL2}"
echo "### DIR=${DIR}"
echo "### SAMTOOLSPATH=${SAMTOOLSPATH}"

base=`basename ${FASTQ1}`
anotherName=${base/.R1.fastq.gz}
anotherName2=${base/.proj.R1.fastq.gz}
tempBamPrefix=${FASTQ1/.R1.fastq.gz}
cd ${DIR}
##--outSAMattributes All \
echo "TIME:$time starting stranded RNA-star on ${FASTQ1}"
perf stat ${STARPATH}/STAR --genomeDir ${STARREF} \
			--runMode alignReads \
			--sjdbGTFfile ${STARGTF} \
			--limitOutSAMoneReadBytes 90000000 \
			--readFilesCommand zcat \
			--readFilesIn ${FASTQL1} ${FASTQL2} \
			--outSAMtype SAM \
			--outFilterType BySJout \
			--outFilterMultimapNmax 10 \
			--outFilterMismatchNmax 10 \
			--outFilterMismatchNoverLmax 0.1 \
			--alignIntronMin 20 \
			--alignIntronMax 1000000 \
			--alignMatesGapMax 1000000 \
			--alignSJoverhangMin 8 \
			--alignSJDBoverhangMin 1 \
			--seedSearchStartLmax 30 \
			--chimSegmentMin 15 \
			--chimJunctionOverhangMin 15 \
			--runThreadN 8 \
			--outReadsUnmapped Fastx \
			--genomeLoad NoSharedMemory \
			--outSAMunmapped Within \
			--outSAMmapqUnique 255 \
			--outSAMattrRGline ${RGTAGLIST} \
			--outSAMmode Full > ${DIR}.starOut 2> ${DIR}/${anotherName}.star.perfOut
if [ $? -eq 0 ] ; then
	echo "### Success. Star finished OK."
	#check output from next four commands
	echo "### Starting sam to bam for Aligned and Chimeric sams"
	perf stat ${SAMTOOLSPATH}/samtools view -bS Aligned.out.sam > Aligned.out.bam 2> ${DIR}/${anotherName}.sam2bam1.perfOut
	perf stat ${SAMTOOLSPATH}/samtools view -bS Chimeric.out.sam > Chimeric.out.bam 2> ${DIR}/${anotherName}.sam2bam2.perfOut

	echo "### Starting sorting for Aligned and Chimeric sams"
	perf stat ${SAMTOOLSPATH}/samtools sort -@4 -m8G Aligned.out.bam Aligned.out.sorted 2> ${DIR}/${anotherName}.sortBam1.perfOut
	perf stat ${SAMTOOLSPATH}/samtools sort -@4 -m8G Chimeric.out.bam Chimeric.out.sorted 2> ${DIR}/${anotherName}.sortBam2.perfOut

	echo "### Starting bam indexing for Aligned and Chimeric bams"
	perf stat ${SAMTOOLSPATH}/samtools index Aligned.out.sorted.bam 2> ${DIR}/${anotherName}.bamIndex1.perfOut
	perf stat ${SAMTOOLSPATH}/samtools index Chimeric.out.sorted.bam 2> ${DIR}/${anotherName}.bamIndex2.perfOut

	echo "### Starting to move the files to their new name"
	mv Aligned.out.sam ${anotherName}.Aligned.out.sam
	mv Aligned.out.bam ${anotherName}.Aligned.out.bam

	mv Aligned.out.sorted.bam ${anotherName}.Aligned.out.sorted.bam
	mv Aligned.out.sorted.bam.bai ${anotherName}.Aligned.out.sorted.bai

	mv Chimeric.out.sam ${anotherName}.Chimeric.out.sam
	mv Chimeric.out.bam ${anotherName}.Chimeric.out.bam

	mv Chimeric.out.sorted.bam ${anotherName}.Chimeric.out.sorted.bam
	mv Chimeric.out.sorted.bam.bai ${anotherName}.Chimeric.out.sorted.bai

	mv Chimeric.out.junction ${anotherName2}.starChimeric.junctions
	mv SJ.out.tab ${anotherName2}.starAligned.junctions

	touch ${RUNDIR}/${NXT1}
	touch ${RUNDIR}/${NXT2}
	touch ${RUNDIR}/${NXT3}
	touch ${RUNDIR}/${NXT4}
	touch ${RUNDIR}/${NXT5}
	touch ${RUNDIR}/${NXT6}
	touch ${RUNDIR}/${NXT7}
	touch ${RUNDIR}/${NXT8}
	touch ${RUNDIR}/${NXT9}
	touch ${RUNDIR}/${NXT10}
	touch ${RUNDIR}/${NXT11}
	mv ${DIR}.starOut ${DIR}.starPass	
else
	echo "### Fail. Star failed."
	mv ${DIR}.starOut ${DIR}.starFail
fi
rm -f ${DIR}.starInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:STAR:$hours:$mins" > ${DIR}/$anotherName.star.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "TIME:$time star finished"
