#!/usr/bin/env bash
#SBATCH --job-name="pegasus_bwaMem"
#SBATCH --time=0-72:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL


time=`date +%d-%m-%Y-%H-%M` 
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### FASTQ1: ${FASTQ1}"
echo "### FASTQ2: ${FASTQ2}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### BWAPATH: ${BWAPATH}"

echo "### BWA mem started for ${FASTQ1} at $time"
${BWAPATH}/bwa mem -R ${RGTAG} -M -t8 ${REF} ${FASTQ1} ${FASTQ2} | ${SAMTOOLSPATH}/samtools view -S -h -b -t ${FAI} - | ${SAMTOOLSPATH}/samtools sort - ${BAMPRE}
if [ $? -eq 0 ] ; then
    ${SAMTOOLSPATH}/samtools index ${BAMPRE}.bam
    if [ $? -eq 0 ] ; then
        touch ${RUNDIR}/${NXT1}
        touch ${RUNDIR}/${NXT2}
        touch ${BAMPRE}.bam.dnaAlignPass
    else
        touch ${BAMPRE}.bam.dnaAlignFail
    fi
else
    touch ${BAMPRE}.bam.dnaAlignFail
fi
rm ${BAMPRE}.bam.dnaAlignInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:BWAMEM:$hours:$mins" > ${BAMPRE}.bwaMem.totalTime
time=`date +%d-%m-%Y-%H-%M` 
echo "bwa mem ended at $time"
echo "BWA MEM completed at $time for ${BAMPRE}" >> ${RUNDIR}/ProjectRunSummary.txt
