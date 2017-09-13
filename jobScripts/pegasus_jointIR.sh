#!/usr/bin/env bash
#SBATCH --job-name="pegasus_jointIR"
#SBATCH --time=0-96:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL


beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### NXT2: ${NXT2}"
echo "### NXT3: ${NXT3}"
echo "### NXT4: ${NXT4}"
echo "### NXT5: ${NXT5}"
echo "### NXT6: ${NXT6}"
echo "### NXT7: ${NXT7}"
echo "### NXT8: ${NXT8}"
echo "### TRK: ${TRK}"
echo "### INDELS: ${INDELS}"
echo "### GATKPATH: ${GATKPATH}"

cd ${WORKDIR}

echo "### Starting indel realigning of file ${BAMLIST}"
echo "### Step 1, target creator..."

test=0
if [ $test -eq 0 ] ; then
java -Xmx15g -Djava.io.tmpdir=$TMPDIR \
    -jar ${GATKPATH}/GenomeAnalysisTK.jar \
    ${BAMLIST} \
    -R ${REF} \
    -T RealignerTargetCreator \
    -nt 16 \
    --maxIntervalSize 350 \
    -DBQ 1 \
    -o ${INTS} \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -known ${INDELS} > ${TRK}.jointIROut

if [ $? -ne 0 ] ; then
    echo "### JIR failed at RealignerTargetCreator stage"
    mv ${TRK}.jointIROut ${TRK}.jointIRFail
    rm -f ${TRK}.jointIRInQueue
    exit

fi

echo "### Starting step 2, indel realignment"
java -Xmx44g -Djava.io.tmpdir=$TMPDIR \
    -jar ${GATKPATH}/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    ${BAMLIST} \
    -R ${REF} \
    -DBQ 1 \
    -targetIntervals ${INTS} \
    --maxReadsInMemory 5000000 \
    --maxConsensuses 24 \
    --maxReadsForConsensuses 80 \
    --maxReadsForRealignment 12000 \
    --nWayOut .jr.bam \
    -model KNOWNS_ONLY \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -known ${INDELS} >> ${TRK}.jointIROut


if [ $? -eq 0 ] ; then
    mv ${TRK}.jointIROut ${TRK}.jointIRPass
    echo "### Starting jr.bam moving"
    for item in ${BAMLIST}
    do
        if [ "$item" == "-I" ] ; then
            continue
        fi
        itemDir=`dirname $item`
        bamName=`basename $item`
        newName=${bamName/.proj.md.bam/.proj.md.jr.bam}
        newBai=${newName/.md.jr.bam/.md.jr.bai}
        hereName=${WORKDIR}/$newName
        hereBai=${hereName/.md.jr.bam/.md.jr.bai}
        echo "### Moving $hereName"
        if [ -e $itemDir/$newName ] ; then
            echo "### Already exists on target, possibly from another joint IR"
        else
            echo "### Does not exist on target, copying now..."
            mv $hereName $itemDir/$newName
            mv $hereBai $itemDir/$newbai
            echo "### Moved out of here to its own dir at $itemDir" > $hereName
            echo "### Moved out of here to its own dir at $itemDir" > $hereBai
            touch $itemDir/$newName.jointIRPass
        fi
    done
    echo "### Done jr.bam moving"
    touch ${RUNDIR}/${NXT1}
    touch ${RUNDIR}/${NXT2}
    touch ${RUNDIR}/${NXT3}
    touch ${RUNDIR}/${NXT4}
    touch ${RUNDIR}/${NXT5}
    touch ${RUNDIR}/${NXT6}
    touch ${RUNDIR}/${NXT7}
    touch ${RUNDIR}/${NXT8}
else
    mv ${TRK}.jointIROut ${TRK}.jointIRFail
fi

rm -f ${TRK}.jointIRInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:GATKJIR:$hours:$mins" > ${TRK}.jir.totalTime
