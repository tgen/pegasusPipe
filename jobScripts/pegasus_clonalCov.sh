##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_clonalCov"
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
echo "### OUTFILE: ${OUTFILE}"
echo "### BAMFILE: ${BAMFILE}"
echo "### NXT1: ${NXT1}"
echo "### CPATH: ${CPATH}"
echo "### SAMPATH: ${SAMPATH}"

echo "Starting clonal coverage for ${BAMFILE} at $time."
#perf stat ${CPATH}/tgen_CloneCov.v007.pl
perf stat ${CPATH}/tgen_CloneCov.v0092.pl \
	I=${BAMFILE} \
	O=${OUTFILE} \
	M=RG: \
	S=${SAMPATH}/samtools > ${BAMFILE}.clonalCovOut 2> ${BAMFILE}.clonalCov.perfOut
if [ $? -eq 0 ] ; then
	mv ${BAMFILE}.clonalCovOut ${BAMFILE}.clonalCovPass
	touch ${RUNDIR}/${NXT1}
else	
	mv ${BAMFILE}.clonalCovOut ${BAMFILE}.clonalCovFail
fi
rm -f ${BAMFILE}.clonalCovInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:CLNCOV:$hours:$mins" > ${BAMFILE}.clonCov.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "Finished clonal coverage at $time."