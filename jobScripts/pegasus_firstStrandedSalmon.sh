#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_salmon"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"

beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### SALMONPATH: ${SALMONPATH}"
echo "### SAMPLE: ${SAMPLE}"
echo "### GTF: ${GTF}"
echo "### SALMON_INDEX_cDNA: ${SALMON_INDEX_cDNA}"
echo "### SALMON_INDEX_GTF: ${SALMON_INDEX_GTF}"
echo "### FASTQ1: ${FASTQ1}"
echo "### FASTQ2: ${FASTQ2}"
echo "### DIR: ${DIR}"
echo "### NXT1: ${NXT1}"
echo "### RUNDIR: ${RUNDIR}"

echo "### Starting salmon..."

module load salmon/0.6.0
#module load salmon/0.5.0
##unstranded salmon

mkdir -p ${DIR}/cDNA
mkdir -p ${DIR}/gtf

perf stat ${SALMONPATH}/salmon quant \
        --index ${SALMON_INDEX_cDNA} \
        --libType ISF \
        --mates1 <(zcat ${FASTQ1}) \
        --mates2 <(zcat ${FASTQ2}) \
        --threads 16 \
        --biasCorrect \
        --geneMap ${GTF} \
        --output ${DIR}/cDNA 2> ${DIR}.salmoncDNA.perfOut
if [ $? -eq 0 ] ; then
        echo "salmon cDNA passed"
	cd ${DIR}/cDNA
	#mv quant_bias_corrected.genes.sf ${SAMPLE}_salmon_bc_cDNA_genes.sf
	#mv quant_bias_corrected.sf ${SAMPLE}_salmon_bc_cDNA_transcripts.sf
	mv quant.genes.sf ${SAMPLE}_salmon_bc_cDNA_genes.sf
	mv quant.sf ${SAMPLE}_salmon_bc_cDNA_transcripts.sf

	perf stat ${SALMONPATH}/salmon quant \
		--index ${SALMON_INDEX_GTF} \
		--libType ISF \
		--mates1 <(zcat ${FASTQ1}) \
		--mates2 <(zcat ${FASTQ2}) \
		--threads 16 \
		--biasCorrect \
		--geneMap ${GTF} \
		--output ${DIR}/gtf 2> ${DIR}.salmonGTF.perfOut
	if [ $? = 0 ] ; then
		echo "Salmon GTF Passed"
		cd ${DIR}/gtf
		#mv quant_bias_corrected.genes.sf ${SAMPLE}_salmon_bc_gtf_genes.sf
		#mv quant_bias_corrected.sf ${SAMPLE}_salmon_bc_gtf_transcripts.sf
		mv quant.genes.sf ${SAMPLE}_salmon_bc_gtf_genes.sf
		mv quant.sf ${SAMPLE}_salmon_bc_gtf_transcripts.sf
		touch ${DIR}.salmonPass
	else
		echo "Salmon GTF Failed"
		touch ${DIR}.salmonFail
	fi
else
	echo "salmon cDNA failed"
	touch ${DIR}.salmonFail
fi
rm ${DIR}.salmonInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:SALMON:$hours:$mins" > ${DIR}.salmon.totalTime
echo "### Ending salmon"
