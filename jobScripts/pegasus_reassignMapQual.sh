##### Author: Megan Russell #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_reassignMap"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.out"
#SBATCH --error="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.err"

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### GATKPATH: ${GATKPATH}"
echo "### RNABAM: ${RNABAM}"
echo "### OUTBAM: ${OUTBAM}"

echo "### GATK reassign Mapping Quality started at $time."
perf stat java -Djava.io.tmpdir=/scratch/tgenjetstream/tmp/ -jar -Xmx32g ${GATKPATH}/GenomeAnalysisTK.jar \
-l INFO \
-R ${REF} \
-T PrintReads \
-rf ReassignMappingQuality \
-I ${RNABAM} \
-o ${OUTBAM} > ${RNABAM}.reassignMapOut 2> ${RNABAM}.reassignMap.perfOut
if [ $? -eq 0 ] ; then
	mv ${RNABAM}.reassignMapOut ${RNABAM}.reassignMapPass
	touch ${RUNDIR}/${NXT1}
else	
	mv ${RNABAM}.reassignMapOut ${RNABAM}.reassignMapFail
	rm -f ${RNABAM}.reassignMapInQueue
	exit
fi
rm -f ${RNABAM}.reassignMapInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:GATKreassignMap:$hours:$mins" > ${RNABAM}.reassignMap.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "reassignMap finished at $time."