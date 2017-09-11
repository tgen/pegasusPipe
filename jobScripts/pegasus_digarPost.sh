#!/usr/bin/env bash
#SBATCH --job-name="pegasus_digarPost"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8

 
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

${DIGARPATH}/digarPost.current.pl \
	--gtf ${GTF} \
	--ann ${ANN} \
	--list ${GENEFILE} \
        --samPath ${SAMTOOLSPATH} \
        --path ${DIGARPATH} \
	--dir ${DIGARDIR} > ${DIGARDIR}.digarPostOut 2>&1
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
