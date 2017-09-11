#!/usr/bin/env bash
#SBATCH --job-name="pegasus_deNovoGear"
#SBATCH --time=0-96:00:00
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
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### GATKPATH: ${GATKPATH}"
echo "### DENOVOPATH: ${DENOVOPATH}"
echo "### VCF: ${VCF}"
echo "### OUTVCF: ${OUTVCF}"
echo "### PED: ${PED}"
echo "### OUTTRACKNAME: ${OUTTRACKNAME}"
echo "### BCFTOOLSPATH: ${BCFTOOLSPATH}"
echo "### BAMLIST: ${BAMLIST}"
echo "### BAMFILE: ${BAMFILE}"

echo "### denovogear started at $time."

cd ${TRACKNAME}
echo "${SAMTOOLSPATH}/samtools mpileup -gDf ${REF} --bam-list ${BAMFILE} | ${DENOVOPATH}/denovogear dnm auto --ped ${PED} --bcf - "
${SAMTOOLSPATH}/samtools mpileup -gDf ${REF} ${BAMFILE} | ${DENOVOPATH}/denovogear dnm auto --ped ${PED} --bcf -

if [ $? -eq 0 ] ; then
	touch ${OUTTRACKNAME}.deNovoGearPass
else
	touch ${OUTTRACKNAME}.deNovoGearFail
fi

rm -f ${OUTTRACKNAME}.deNovoGearInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:DENOVOGEAR:$hours:$mins" > ${OUTTRACKNAME}.deNovoGear.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "DenovoGear finished at $time."
