#!/usr/bin/env bash
#SBATCH --job-name="pegasus_mergeFQ"
#SBATCH --time=0-24:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL

time=`date +%d-%m-%Y-%H-%M` 
beginTime=`date +%s`
echo "merging of fastqs started at $time"
machine=`hostname`

echo "### NODE: $machine"
echo "### FAST: ${FASTQLIST}"
echo "### MERG: ${MERGEDFASTQ}"
echo "### CNT : ${CNT}"
echo "### RDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### NXT2: ${NXT2}"
echo "### NXT3: ${NXT3}"

failCount=0
cd ${DIR}
echo "TIME:$time starting fastq merging to create ${MERGEDFASTQ}"
if [ ${CNT} -eq 1 ] ; then
    echo "only one thing to merge, commands are:"
    echo "cp ${FASTQLIST} ${MERGEDFASTQ}"
    cp ${FASTQLIST} ${MERGEDFASTQ}
    if [ $? -ne 0 ] ; then #check if foreground finished OK
        ((failCount++))
    fi
else
    echo "more than one thing to merge, commands are:"
    echo "cat ${FASTQLIST} >> ${MERGEDFASTQ}"
    cat ${FASTQLIST} >> ${MERGEDFASTQ}
    if [ $? -ne 0 ] ; then #check if background finished OK
        ((failCount++))
    fi
fi

if [ $failCount -eq 0 ] ; then
    touch ${MERGEDFASTQ}.mergeFastqPass
    touch ${RUNDIR}/${NXT1}
    touch ${RUNDIR}/${NXT2}
    touch ${RUNDIR}/${NXT3}
    touch ${RUNDIR}/${NXT4}

else
    touch ${MERGEDFASTQ}.mergeFastqFail
fi

rm ${MERGEDFASTQ}.mergeFastqInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:MERGEFASTQ:$hours:$mins" > ${MERGEDFASTQ}.mergeFastq.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "TIME:$time fastq merging finished"
