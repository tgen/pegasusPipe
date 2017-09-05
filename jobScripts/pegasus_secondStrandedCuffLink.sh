#$#### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_SScuffLink"
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
echo "### CLPATH: ${CUFFLINKSPATH}"
echo "### CLGTF: ${CUFFLINKGTF}"
echo "### CLMASK: ${CUFFLINKMASK}"
echo "### NXT1: ${NXT1}"
echo "### NXT2: ${NXT2}"
echo "### PARAMS: ${PARAMS}"

echo "TIME:$time starting cufflinks on ${DIRNAME}"
cd ${DIRNAME}
PARAMS=${PARAMS//\@/\$}
echo "### 1. params are: ${PARAMS}"
PARAMS=${PARAMS//\#/ }
echo "### 2. params are: $PARAMS"

echo "Running cuffLinks for second stranded RNA capture"
#perf stat ${CUFFLINKSPATH}/cufflinks ${PARAMS} 2> ${DIRNAME}.cuffLink.perfOut > ${DIRNAME}.cuffLinkOut 2>&1
if [ ${USEMASK} == "no" ] ; then

        perf stat ${CUFFLINKSPATH}/cufflinks ${PARAMS} --library-type fr-firststrand --frag-bias-correct ${REF} --GTF ${CUFFLINKGTF} ${BAM} 2> ${DIRNAME}.cuffLink.perfOut > ${DIRNAME}.cuffLinkOut 2>&1
        if [ $? -eq 0 ] ; then
                newName=`basename ${BAM}`
                newName=${newName/.proj.Aligned.out.sorted.md.bam}
                mv ${DIRNAME}.cuffLinkOut ${DIRNAME}.cuffLinkPass
                mv ${DIRNAME}/transcripts.gtf ${DIRNAME}/$newName.cufflinks.transcripts.gtf
                mv ${DIRNAME}/skipped.gtf ${DIRNAME}/$newName.cufflinks.skipped.gtf
                mv ${DIRNAME}/genes.fpkm_tracking ${DIRNAME}/$newName.cufflinks.genes.fpkm_tracking
                mv ${DIRNAME}/isoforms.fpkm_tracking ${DIRNAME}/$newName.cufflinks.isoforms.fpkm_tracking
        else
                mv ${DIRNAME}.cuffLinkOut ${DIRNAME}.cuffLinkFail
        fi
else

	perf stat ${CUFFLINKSPATH}/cufflinks ${PARAMS} --library-type fr-firststrand --frag-bias-correct ${REF} --GTF ${CUFFLINKGTF} -M ${CUFFLINKMASK} ${BAM} 2> ${DIRNAME}.cuffLink.perfOut > ${DIRNAME}.cuffLinkOut 2>&1
	if [ $? -eq 0 ] ; then
		newName=`basename ${BAM}`
		newName=${newName/.proj.Aligned.out.sorted.md.bam}
		mv ${DIRNAME}.cuffLinkOut ${DIRNAME}.cuffLinkPass
		mv ${DIRNAME}/transcripts.gtf ${DIRNAME}/$newName.cufflinks.transcripts.gtf
		mv ${DIRNAME}/skipped.gtf ${DIRNAME}/$newName.cufflinks.skipped.gtf
		mv ${DIRNAME}/genes.fpkm_tracking ${DIRNAME}/$newName.cufflinks.genes.fpkm_tracking
		mv ${DIRNAME}/isoforms.fpkm_tracking ${DIRNAME}/$newName.cufflinks.isoforms.fpkm_tracking
		#mv ${DIR}/transcripts.expr ${DIR}/${NEWNAME}_transcripts.expr
		#mv ${DIR}/genes.expr ${DIR}/${NEWNAME}_genes.expr
		#mv ${DIR}/isoforms.fpkm_tracking ${DIR}/${NEWNAME}_isoforms.fpkm_tracking
		#touch ${RUNDIR}/${NXT1}
	else
		mv ${DIRNAME}.cuffLinkOut ${DIRNAME}.cuffLinkFail
	fi
fi
rm -f ${DIRNAME}.cuffLinkInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:CUFFLINKS:$hours:$mins" > ${DIRNAME}.cufflinks.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "TIME:$time finished cufflinks on ${DIRNAME}"
