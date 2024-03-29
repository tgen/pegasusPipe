#!/usr/bin/env bash
#SBATCH --job-name="pegasus_FScuffLink"
#SBATCH --time=0-240:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL

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

echo "Running cuffLinks for first stranded RNA capture"
if [ ${USEMASK} == "no" ] ; then

        ${CUFFLINKSPATH}/cufflinks ${PARAMS} --library-type fr-secondstrand --frag-bias-correct ${REF} --GTF ${CUFFLINKGTF} ${BAM} > ${DIRNAME}.cuffLinkOut 2>&1
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

    ${CUFFLINKSPATH}/cufflinks ${PARAMS} --library-type fr-secondstrand --frag-bias-correct ${REF} --GTF ${CUFFLINKGTF} -M ${CUFFLINKMASK} ${BAM} > ${DIRNAME}.cuffLinkOut 2>&1
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
fi
rm -f ${DIRNAME}.cuffLinkInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:CUFFLINKS:$hours:$mins" > ${DIRNAME}.cufflinks.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "TIME:$time finished cufflinks on ${DIRNAME}"
