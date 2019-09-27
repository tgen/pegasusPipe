#!/usr/bin/env bash
#SBATCH --job-name="pegasus_jointIRsplit"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu 4096


beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### TRK: ${TRK}"
echo "### INDELS: ${INDELS}"
echo "### CHRGRP: ${CHRGRP}"
echo "### GRPNAME: ${GRPNAME}"
echo "### INTERVALS: ${INTS}"

cd ${WORKDIR}

colons="${CHRGRP//[^:]}"
colonCount="${#colons}"
declare -a bgJobIDs=()
chrCount=`echo $colonCount + 1 | bc`
for (( i=1; i<=$chrCount; i++ ))
do
    thisChr=`echo ${CHRGRP} | cut -d: -f$i`
    echo "### My i is $i and it is chr $thisChr"
    echo "### This is chr $thisChr, from group ${GRPNAME}, loop index is $i" > ${TRK}.jointIROut.chr$thisChr
    echo "### Starting indel realigning of file ${BAMLIST}"
    intervals=${INTS/.intervals/.chr$thisChr.intervals}
    echo "### Step 1, target creator..."
    java -Xmx2g -Djava.io.tmpdir=$TMPDIR \
        -jar ${GATKPATH}/GenomeAnalysisTK.jar \
        ${BAMLIST} \
        -R ${REF} \
        -T RealignerTargetCreator \
        -L $thisChr \
        -nt 4 \
        --maxIntervalSize 350 \
        -DBQ 1 \
        -o $intervals \
        -known ${INDELS} >> ${TRK}.jointIROut.chr$thisChr &
    bgJobIDs[$i]=$!
done

echo "### All jobs submitted, now time to wait..."
for (( i=1; i<=$chrCount; i++ ))
do
        thisChr=`echo ${CHRGRP} | cut -d: -f$i`
        echo "### Loop for step 1 bg jobs: My i $i and bg job id is ${bgJobIDs[$i]}"
        wait ${bgJobIDs[$i]}
done

echo "### Waiting for bg jobs is over for target creator (step 1) is done."
for (( i=1; i<=$chrCount; i++ ))
do
    thisChr=`echo ${CHRGRP} | cut -d: -f$i`
    echo "### My i is $i and it is chr $thisChr"
    echo "### This is chr $thisChr, from group ${GRPNAME}, loop index is $i" >> ${TRK}.jointIROut.chr$thisChr
    echo "### Starting indel realigning of file ${BAMLIST}"
    echo "### Starting step 2, indel realignment"
    java -Xmx16g -Djava.io.tmpdir=$TMPDIR \
        -jar ${GATKPATH}/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        ${BAMLIST} \
        -R ${REF} \
        -DBQ 1 \
        -L $thisChr \
        -targetIntervals $intervals \
        --maxReadsInMemory 5000000 \
        --maxConsensuses 24 \
        --maxReadsForConsensuses 80 \
        --maxReadsForRealignment 12000 \
        --nWayOut .chr$thisChr.jr.bam \
        -model KNOWNS_ONLY \
        -known ${INDELS} >> ${TRK}.jointIROut.chr$thisChr
    bgJobIDs[$i]=$!
done

echo "### All jobs submitted, now time to wait..."
for (( i=1; i<=$chrCount; i++ ))
do
    thisChr=`echo ${CHRGRP} | cut -d: -f$i`
    echo "### Loop for step 2 bg jobs: My i $i and bg job id is ${bgJobIDs[$i]}"
    wait ${bgJobIDs[$i]}

    if [ $? -eq 0 ] ; then
        mv ${TRK}.jointIROut.chr$thisChr ${TRK}.jointIRPass.chr$thisChr
    else
        mv ${TRK}.jointIROut.chr$thisChr ${TRK}.jointIRFail.chr$thisChr
    fi
done

echo "### Waiting for bg jobs is over for indel realigner is done."

echo "### All jobs finished. Now time to check if all passed"
failedChr=0
for (( i=1; i<=$chrCount; i++ ))
do
    thisChr=`echo ${CHRGRP} | cut -d: -f$i`
    echo "### Checking if pass exists on $i and bg job id is ${bgJobIDs[$i]}"
    if [ ! -e ${TRK}.jointIRPass.chr$thisChr ] ; then
            failedChr=1
    fi
done

if [ $failedChr -eq 0 ] ; then  # This checks if all chromosomes submitted to this node finished OK
    echo "joint IR split calling passed on ${GRPNAME}" >> ${TRK}.jointIR-group${GRPNAME}-InQueue
    mv ${TRK}.jointIR-group${GRPNAME}-InQueue ${TRK}.jointIR-group${GRPNAME}-Pass
    touch ${RUNDIR}/${NXT1}
else
    echo "joint IR split calling failed on ${GRPNAME}" >> ${TRK}.jointIR-group${GRPNAME}-InQueue
    mv ${TRK}.jointIR-group${GRPNAME}-InQueue ${TRK}.jointIR-group${GRPNAME}-Fail
fi
rm -f ${TRK}.jointIR-group${GRPNAME}-InQueue


endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:GATKJIR:$hours:$mins" > ${TRK}.jir.totalTime
