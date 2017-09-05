##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_deSeq2"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"
 
time=`date +%d-%m-%Y-%H-%M` 
beginTime=`date +%s`
machine=`hostname`

module add R/3.1.1

echo "### NODE: $machine"
echo "### RUNDIR: ${RUNDIR}"
echo "### GTF: ${GTF}"
echo "### NORMLIST: ${NORMLIST}"
echo "### TUMORLIST: ${TUMORLIST}"
echo "### NXT1: ${NXT1}"
echo "### DIRNAME: ${DIRNAME}"

cd ${DIRNAME}
newName=`basename ${DIRNAME}`
newName=${newName/.dsDir/}

echo "### Starting DESeq3.R"
echo "perf stat /packages/R/3.1.1/bin/Rscript --vanilla ${DESEQ2PATH}/DESeq_v2.R ${NORMLIST} ${TUMORLIST} 2> ${DIRNAME}.deSeq2.perfOut"
perf stat /packages/R/3.1.1/bin/Rscript --vanilla ${DESEQ2PATH}/DESeq_v2.R ${NORMLIST} ${TUMORLIST} 2> ${DIRNAME}.deSeq2.perfOut
if [ $? -eq 0 ] ; then
	touch ${DIRNAME}.deSeq2Pass
	#echo "running this cmd to rename: mv ${DIRNAME}/DESeq_results.txt ${DIRNAME}/$newName.DESeq_results.txt"
	mv ${DIRNAME}/DESeq_results.txt ${DIRNAME}/$newName.DESeq2_results.txt
	#start conversion
	echo "${DESEQ2PATH}/deseq2vcfNewGtf.pl ${GTF} ${DIRNAME}/$newName.DESeq2_results.txt ${NORMLIST}"
	perf stat ${DESEQ2PATH}/deseq2vcfNewGtf.pl ${GTF} ${DIRNAME}/$newName.DESeq2_results.txt ${NORMLIST} 2> ${DIRNAME}.deSeq2vcf.perfOut
	#end conversion
else
	touch ${DIRNAME}.deSeq2Fail
fi
rm -f ${DIRNAME}/masterNor.txt
rm -f ${DIRNAME}/masterTum.txt
rm -f ${DIRNAME}/CountsTable.txt

rm -f ${DIRNAME}.deSeq2InQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:DESEQ2:$hours:$mins" > ${DIRNAME}.deseq2.totalTime
echo "ending deSeq"
