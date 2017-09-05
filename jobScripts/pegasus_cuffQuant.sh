#$#### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_cuffQuant"
#SBATCH --time=0-240:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"
 
time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### DIRNAME: ${DIRNAME}"
echo "### RUNDIR: ${RUNDIR}"
echo "### BAM: ${BAM}"
echo "### USEGTF: ${USEGTF}"
echo "### USEMASK: ${USEMASK}"
echo "### CQPATH: ${CUFFQUANTPATH}"
echo "### CLGTF: ${CUFFLINKGTF}"
echo "### CLMASK: ${CUFFLINKMASK}"
echo "### NXT1: ${NXT1}"
echo "### PARAMS: ${PARAMS}"

echo "TIME:$time starting cuff quant on ${DIRNAME}"
cd ${DIRNAME}

PARAMS=${PARAMS//\#/ }
echo "### params is $params"
if [ ${USEMASK} == "no" ] ; then
	perf stat ${CUFFQUANTPATH}/cuffquant ${PARAMS} --frag-bias-correct ${REF} ${CUFFLINKGTF} ${BAM} 2> ${DIRNAME}.cuffQuant.perfOut > ${DIRNAME}.cuffQuantOut 2>&1
        if [ $? -eq 0 ] ; then
                newName=`basename ${BAM}`
                newName=${newName/.proj.Aligned.out.sorted.md.bam}
                mv ${DIRNAME}.cuffQuantOut ${DIRNAME}.cuffQuantPass
                mv ${DIRNAME}/abundances.cxb ${DIRNAME}/$newName.cuffQuant.abundances.cxb
        else
                mv ${DIRNAME}.cuffQuantOut ${DIRNAME}.cuffQuantFail
        fi

else
	perf stat ${CUFFQUANTPATH}/cuffquant ${PARAMS} --frag-bias-correct ${REF} --mask-file ${CUFFLINKMASK} ${CUFFLINKGTF} ${BAM} 2> ${DIRNAME}.cuffQuant.perfOut > ${DIRNAME}.cuffQuantOut 2>&1
	if [ $? -eq 0 ] ; then
		newName=`basename ${BAM}`
		newName=${newName/.proj.Aligned.out.sorted.md.bam}
		mv ${DIRNAME}.cuffQuantOut ${DIRNAME}.cuffQuantPass	
		mv ${DIRNAME}/abundances.cxb ${DIRNAME}/$newName.cuffQuant.abundances.cxb
	else
		mv ${DIRNAME}.cuffQuantOut ${DIRNAME}.cuffQuantFail
	fi
fi
rm -f ${DIRNAME}.cuffQuantInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:CUFFQuant:$hours:$mins" > ${DIRNAME}.cuffquant.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "TIME:$time finished cuffQuant on ${DIRNAME}"
