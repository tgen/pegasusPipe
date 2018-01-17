#!/usr/bin/env bash
#SBATCH --job-name="pegasus_splitFq"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL

time=`date +%d-%m-%Y-%H-%M`
machine=`hostname`
echo "### NODE: $machine"
echo "### FQ: ${FQ}"
echo "### DIR: ${DIR}"
echo "### PF: ${PF}"

echo "TIME:$time starting split fastq on ${FQ} on machine $machine"
cd ${DIR}
badExitStatus=0

zcat ${FQ} | split -a 3 -d -l 400000000 - ${PF}
if [ $? -eq 0 ] ; then
    echo "zcatting and splitting is over"
    echo "now have to gzip them again"
    for splitFile in `ls ${PF}*`
    do
        newName=${splitFile/split.fastq.gz/}
        newName=$newName.part.fastq
        echo "now $splitFile is becoming $newName"
        mv $splitFile $newName
        echo "now gzipping $newName"
        gzip $newName
        if [ $? -ne 0 ] ; then
            ((badExitStatus++))
        fi
    done
    if [ $badExitStatus -eq 0 ] ; then
        touch ${FQ}.fastqSplitPass
        touch ${RUNDIR}/${NXT1}
    else
        touch ${FQ}.splitFastqFail
    fi
else
    touch ${FQ}.splitFastqFail
fi

rm -f ${FQ}.fastqSplitInQueue

time=`date +%d-%m-%Y-%H-%M`
echo "TIME:$time finished splitFastq on ${FASTQ}"
