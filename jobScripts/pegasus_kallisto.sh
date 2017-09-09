#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_kallisto"
#SBATCH --time=0-60:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
##PBS -e /${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err

#updated module 12/8/15
module load kallisto/0.42.4
#module load kallisto/0.42.3
#module load kallisto/0.43.0

beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### DIR: ${DIR}"
echo "### NXT1: ${NXT1}"
echo "### RUNDIR: ${RUNDIR}"
echo "### FASTQL: ${FASTQL}"
echo "### KALLISTO_INDEX_CDNA: ${KALLISTO_INDEX_CDNA}"
echo "### KALLISTO_INDEX_GTF: ${KALLISTO_INDEX_GTF}"

FASTQL=`echo "${FASTQL}" | tr ',' ' '`
echo "FASTQL: ${FASTQL}"
#DIR=/scratch/sszeling/C4RCD/CRDC0154/C4RCD_0073/
#FASTQ_DIR=/scratch/tgenjetstream/chiaPipe/runFolders/150914_7001452_0258_BC7B0WACXX/Project_C4RCD_0073_2_TSMRU/Sample_C4RCD_0073_2_PB_Whole_F0147G2S2M0082P0086_TSMRU_D00932_C7B0WACXX/
#FASTQ1=C4RCD_0073_2_PB_Whole_F0147G2S2M0082P0086_TSMRU_D00932_C7B0WACXX_CTTGTA_L003_R1_001.fastq.gz
#FASTQ2=C4RCD_0073_2_PB_Whole_F0147G2S2M0082P0086_TSMRU_D00932_C7B0WACXX_CTTGTA_L003_R2_001.fastq.gz
#FASTQ3=C4RCD_0073_2_PB_Whole_F0147G2S2M0082P0086_TSMRU_D00932_C7B0WACXX_CTTGTA_L007_R1_001.fastq.gz
#FASTQ4=C4RCD_0073_2_PB_Whole_F0147G2S2M0082P0086_TSMRU_D00932_C7B0WACXX_CTTGTA_L007_R2_001.fastq.gz
#kallisto_index=/home/sszelinger/resources/rnaseq/kallisto_cDNA_index/Homo_sapiens.GRCh37.74.cdna.hs37d5.idx
#KALLISTOPATH=/packages/kallisto/0.42.3/bin/kallisto
#kallisto_cDNA_fasta=/home/tgenref/pecan/salmon_index/salmon42_index_cDNA_E74hs37d5/Homo_sapiens.GRCh37.74.cdna.hs37d5.fasta
#kallisto_cDNA_index=/home/sszelinger/resources/rnaseq/kallisto_cDNA_index/Homo_sapiens.GRCh37.74.cdna.hs37d5
#cd $DIR
#
#	echo -e "\n##Starting kallisto at `date +%d-%m-%Y-%H-%M`...\n"

echo "### Starting Kallisto quant cDNA"

perf stat kallisto quant \
	--index=${KALLISTO_INDEX_CDNA} \
	--output-dir=${DIR}/cDNA \
	--bias \
	--bootstrap-samples=100 \
	--threads=14 \
	--seed=42 \
	${FASTQL} 2> ${DIR}.kallistoCDNA.perfOut
	if [ $? -eq 0 ] ; then
		echo "kallisto cDNA Passed"
	else
                echo "kallisto cDNA Failed"
                touch ${DIR}.kallistoFail
        fi
	
	
perf stat kallisto quant \
	--index=${KALLISTO_INDEX_GTF} \
	--output-dir=${DIR}/GTF \
	--bias \
	--bootstrap-samples=100 \
	--threads=14 \
	--seed=42 \
	${FASTQL} 2> ${DIR}.kallistoGTF.perfOut
	if [ $? -eq 0 ] ; then
		echo "kallisto GTF Passed"
		touch ${DIR}.kallistoPass
		#touch ${RUNDIR}/${NXT1}
	else
		echo "kallisto GTF Failed"
		touch ${DIR}.kallistoFail
	fi


rm ${DIR}.kallistoInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:KALLISTOquant:$hours:$mins" > ${DIR}.kallistoQuant.totalTime
echo "### Ending Kallisto."
