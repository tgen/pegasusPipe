##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_sailFish"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.out"
#SBATCH --error="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.err"

beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### SAILFISHPATH: ${SAILFISHPATH}"
echo "### SAMPLE: ${SAMPLE}"
echo "### SAILFISHGTF: ${SAILFISHGTF}"
echo "### CCDSGTF: ${CCDSGTF}"
echo "### SAILFISHINDEXDIR: ${SAILFISHINDEXDIR}"
echo "### SAILFISHCCDSINDEX: ${SAILFISHCCDSINDEX}"
echo "### FASTQ1: ${FASTQ1}"
echo "### FASTQ2: ${FASTQ2}"
echo "### DIR: ${DIR}"
echo "### NXT1: ${NXT1}"
echo "### RUNDIR: ${RUNDIR}"

echo "### Starting sail fish..."

# Note: SLURM defaults to running jobs in the directory
# where they are submitted, no need for $PBS_O_WORKDIR

module load sailfish/0.6.3

OUTPUT_DIR="${DIR}/full"
OUTPUT_DIR_CCDS="${DIR}/ccds"

perf stat sailfish quant -i ${SAILFISHINDEXDIR} --libtype "T=PE:O=><:S=U" --mates1 <(zcat ${FASTQ1}) --mates2 <(zcat ${FASTQ2}) --out ${OUTPUT_DIR} --threads 12 --polya 2> ${DIR}.sailfish1.perfOut 

cd ${OUTPUT_DIR}

#Rename and process output to generate gene level expression estimates
awk '/ERCC-/{gsub(/+/, "")};{print}' quant_bias_corrected.sf > ${SAMPLE}_sailfish_all_transcripts.expr
${SAILFISHPATH}/TranscriptsToGenes.sh --exp-file ${SAMPLE}_sailfish_all_transcripts.expr --gtf-file ${SAILFISHGTF} --res-file ${SAMPLE}_sailfish_all_genes.expr
mv reads.count_info ${SAMPLE}.reads.count_info

cd ..

#### This is for the CCDS transcripts

#Use Process substitution to run with gzip compressed fastq
#The libtype will need to be dynamic for Single End runs, or stranded runs.  Currently assumes Illumina PE, unstranded library
sailfish quant -i ${SAILFISHCCDSINDEX} --libtype "T=PE:O=><:S=U" --mates1 <(zcat ${FASTQ1}) --mates2 <(zcat ${FASTQ2}) --out ${OUTPUT_DIR_CCDS} --threads 12 --polya

cd ${OUTPUT_DIR_CCDS}

#Rename and process output to generate gene level expression estimates
mv quant_bias_corrected.sf ${SAMPLE}_sailfish_ccds_transcripts.expr
mv reads.count_info ${SAMPLE}.reads.count_info
${SAILFISHPATH}/TranscriptsToGenes.sh --exp-file ${SAMPLE}_sailfish_ccds_transcripts.expr --gtf-file ${CCDSGTF} --res-file ${SAMPLE}_sailfish_ccds_genes.expr 2> ${DIR}.sailfish2.perfOut 
if [ $? -eq 0 ] ; then
	touch ${DIR}.sailFishPass
	#touch ${RUNDIR}/${NXT1}
else
	touch ${DIR}.sailFishFail
fi
rm ${DIR}.sailFishInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:SAILFISH:$hours:$mins" > ${DIR}.sailFish.totalTime
echo "### Ending sail fish"
