#$#### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_digarPost"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"
 
time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### DIGARDIR: ${DIGARDIR}"
echo "### RUNDIR: ${RUNDIR}"
echo "### BAM: ${BAM}"
echo "### ANN: ${ANN}"
echo "### GTF: ${GTF}"
echo "### SAMTOOLSPATH: ${SAMTOOLSPATH}"
echo "### BWAPATH: ${BWAPATH}"
echo "### DIGARPATH: ${DIGARPATH}"


echo "TIME:$time starting digar post on ${DIGARDIR}"
cd ${DIGARDIR}

perf stat ${DIGARPATH}/digarPost.current.pl \
	--gtf ${GTF} \
	--ann ${ANN} \
	--list ${GENEFILE} \
        --samPath ${SAMTOOLSPATH} \
        --path ${DIGARPATH} \
	--dir ${DIGARDIR} 2> ${DIGARDIR}.digarPost.perfOut > ${DIGARDIR}.digarPostOut 2>&1
if [ $? -eq 0 ] ; then
	#finished successfully
	mv ${DIGARDIR}.digarPostOut ${DIGARDIR}.digarPostPass
else
	#failed
	mv ${DIGARDIR}.digarPostOut ${DIGARDIR}.digarPostFail
fi	

rm -f ${DIGARDIR}.digarPostInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:DIGARPOST.$hours:$mins" > ${DIGARDIR}.digar.post.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "TIME:$time finished digar post on ${DIGARDIR}."
