#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_deSeq"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"
 
time=`date +%d-%m-%Y-%H-%M` 
beginTime=`date +%s`
machine=`hostname`
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
echo "perf stat /packages/R/2.15.2/bin/Rscript --vanilla ${DESEQPATH}/DESeq3.R ${NORMLIST} ${TUMORLIST} 2> ${DIRNAME}.deSeq.perfOut"
perf stat /packages/R/2.15.2/bin/Rscript --vanilla ${DESEQPATH}/DESeq3.R ${NORMLIST} ${TUMORLIST} 2> ${DIRNAME}.deSeq.perfOut
if [ $? -eq 0 ] ; then
	touch ${DIRNAME}.deSeqPass
	#echo "running this cmd to rename: mv ${DIRNAME}/DESeq_results.txt ${DIRNAME}/$newName.DESeq_results.txt"
	mv ${DIRNAME}/DESeq_results.txt ${DIRNAME}/$newName.DESeq_results.txt
	#start conversion
	${DESEQPATH}/deseq2vcf.pl ${GTF} ${DIRNAME}/$newName.DESeq_results.txt ${NORMLIST}
	#end conversion
else
	touch ${DIRNAME}.deSeqFail
fi
rm -f ${DIRNAME}/masterNor.txt
rm -f ${DIRNAME}/masterTum.txt
rm -f ${DIRNAME}/CountsTable.txt

rm -f ${DIRNAME}.deSeqInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:DESEQ:$hours:$mins" > ${DIRNAME}.deseq.totalTime
echo "ending deSeq"
