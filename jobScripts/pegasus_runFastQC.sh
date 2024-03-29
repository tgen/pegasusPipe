#!/usr/bin/env bash
#SBATCH --job-name="pegasus_runFastQC"
#SBATCH --time=0-24:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL


time=`date +%d-%m-%Y-%H-%M` 
beginTime=`date +%s`
echo "fastQC started at $time"
machine=`hostname`
echo "### NODE: $machine"
echo "### FQ  : ${FQ}"
echo "### RDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### PATH: ${FASTQCPATH}"
qcDir=${FQ/.fastq.gz/_fastqc}
fqName=`basename ${FQ}`
newName=${fqName/.fastq.gz/}

${FASTQCPATH}/fastqc ${FQ}

if [ $? -eq 0 ] ; then
    touch ${FQ}.fastqcPass
    mv $qcDir/Images/per_base_quality.png $qcDir/Images/${newName}_per_base_quality.png
    mv $qcDir/Images/per_base_quality.png $qcDir/Images/${newName}_sequence_length_distribution.png
    mv $qcDir/Images/per_base_quality.png $qcDir/Images/${newName}_per_sequence_quality.png
    mv $qcDir/Images/per_base_quality.png $qcDir/Images/${newName}_per_sequence_gc_content.png
    mv $qcDir/Images/per_base_quality.png $qcDir/Images/${newName}_per_base_sequence_content.png
    mv $qcDir/Images/per_base_quality.png $qcDir/Images/${newName}_per_base_n_content.png
    mv $qcDir/Images/per_base_quality.png $qcDir/Images/${newName}_per_base_gc_content.png
    mv $qcDir/Images/per_base_quality.png $qcDir/Images/${newName}_kmer_profiles.png
    mv $qcDir/Images/per_base_quality.png $qcDir/Images/${newName}_duplication_levels.png
else
    touch ${FQ}.fastqcFail
fi

rm -f ${FQ}.fastqcInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:FASTQC:$hours:$mins" > ${FQ}.fastqc.totalTime
time=`date +%d-%m-%Y-%H-%M` 
echo "fastQC ended at $time"
