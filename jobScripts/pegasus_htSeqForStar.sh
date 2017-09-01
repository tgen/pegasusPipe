##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_htSeq4Star"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.out"
#SBATCH --error="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.err"

module load python/2.7.3
#module load HTSeq/0.5.3p9

beginTime=`date +%s`
time=`date +%d-%m-%Y-%H-%M`
machine=`hostname`
echo "### NODE: $machine"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### PICARDPATH: ${PICARDPATH}"
echo "### BAM: ${BAM}"
echo "### SAM: ${SAM}"

base=`basename ${SAM}`
anotherName=${base/.proj.Aligned.out.sam}
DIR=$(dirname "${SAM}")
echo "### Now running htseq on SAM"
perf stat htseq-count -q --stranded=no --mode=union ${SAM} ${GTF} 2> ${SAM}.htseq.perfOut > ${DIR}/${anotherName}.htSeqCounts
if [ $? -eq 0 ] ; then
	touch ${SAM}.htSeqPass	
	#touch ${RUNDIR}/${NXT1}
	echo "Deleted by htSeq to save on space at $time" > ${SAM}
else
	touch ${SAM}.htSeqFail
fi
rm -f ${SAM}.htSeqInQueue
time=`date +%d-%m-%Y-%H-%M`
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:HTSEQ:$hours:$mins" > ${SAM}.htSeq.totalTime

echo "TIME:$time finished htseq on ${BAM}"
