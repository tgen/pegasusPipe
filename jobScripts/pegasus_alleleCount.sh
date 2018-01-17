#!/usr/bin/env bash
#SBATCH --job-name="pegasus_alleleCount"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 8
 
cd ${RUNDIR}
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### ALCOUNTPATH: ${ALCOUNTPATH}"
echo "### DNABAM: ${DNABAM}"
echo "### RNABAM: ${RNABAM}"
echo "### VCF: ${VCF}"
echo "### OUT: ${OUT}"
echo "### TRACK: ${TRACK}"

echo "### Starting allele count"
echo "${ALCOUNTPATH}/bam_allele_counts_to_vcf.sh -a ${DNABAM} -b ${RNABAM} -v ${VCF} -r ${REF} > ${TRACK}.alleleCountOut"
${ALCOUNTPATH}/bam_allele_counts_to_vcf.sh -a ${DNABAM} -b ${RNABAM} -v ${VCF} -r ${REF} > ${TRACK}.alleleCountOut

if [ $? -eq 0 ] ; then
    mv ${VCF}.allele_counts.vcf ${OUT}
    mv ${TRACK}.alleleCountOut ${TRACK}.alleleCountPass
else
    mv ${TRACK}.alleleCountOut ${TRACK}.alleleCountFail
fi

rm -f ${TRACK}.alleleCountInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:ALLELECOUNT:$hours:$mins" > ${TRACK}.alleleCount.totalTime
echo "### Ending allele count "

