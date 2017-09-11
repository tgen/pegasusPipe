#!/usr/bin/env bash
#SBATCH --job-name="pegasus_IGLbc"
#SBATCH --time=0-32:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL

module load R/2.15.2
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### DIR: ${DIR}"
echo "### IGLLISTBED: ${IGLLISTBED}"
echo "### SAMTOOLS: ${SAMTOOLSPATH}"
echo "### BEDTOOLS: ${BEDTOOLSPATH}"
echo "### BAMFILE: ${BAMFILE}"
echo "### PICARDPATH: ${PICARDPATH}"

echo "### Starting igl bed cov"

${SAMTOOLSPATH}/samtools view -b ${BAMFILE} | ${BEDTOOLSPATH}/bedtools coverage -split -counts -abam stdin -b ${IGLLISTBED} > ${BAMFILE}.bedToolsIglCount.txt
if [ $? -eq 0 ] ; then
	touch ${BAMFILE}.IGLbedCovPass
else
	touch ${BAMFILE}.IGLbedCovFail
fi
rm -f ${BAMFILE}.IGLbedCovInQueue

# A little organizing
if [ -d ${RUNDIR}/stats/ ] ; then
	echo "### Moving files into stats folder"
	mv ${BAMFILE}.bedToolsIglCount.txt ${RUNDIR}/stats/
fi
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:IGLBEDCOV:$hours:$mins" > ${BAMFILE}.iglBedCov.totalTime
echo "### Ending igl bed cov"
