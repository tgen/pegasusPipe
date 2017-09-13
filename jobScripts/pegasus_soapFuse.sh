#!/usr/bin/env bash
#SBATCH --job-name="pegasus_soapFuse"
#SBATCH --time=0-96:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL

echo "### Variables coming in:"
echo "### SAMPLE=${SAMPLE}"
echo "### SPCONFIG=${SPCONFIG}"
echo "### REF=${REF}"
echo "### FAI=${FAI}"
echo "### FASTQ1=${FASTQ1}"
echo "### FASTQ2=${FASTQ2}"
echo "### DIR=${DIR}"
echo "### SOAPFUSEPATH: ${SOAPFUSEPATH}"
echo "### SAMTOOLSPATH: ${SAMTOOLSPATH}"

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`

echo "TIME:$time starting soap fuse on ${FASTQ1}"
base=`basename ${FASTQ1}`
anotherName=${base/.R1.fastq.gz}

${SOAPFUSEPATH}/SOAPfuse-RUN.pl \
    -c ${SPCONFIG} \
    -fd ${DIR} \
    -l ${SLFILE} \
    -o ${DIR} \
    -fs 1 \
    -es 9 \
    -tp ${SAMPLE} > ${DIR}.soapFuseOut

if [ $? -eq 0 ] ; then
    echo "### Success."
    mv ${DIR}.soapFuseOut ${DIR}.soapFusePass
    touch ${RUNDIR}/${NXT1}
else
    mv ${DIR}.soapFuseOut ${DIR}.soapFuseFail
fi

rm -f ${DIR}.soapFuseInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:SOAPFUSE:$hours:$mins" > ${DIR}/$anotherName.accepted_hits.bam.thFusion.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "TIME:$time soap fuse finished"
