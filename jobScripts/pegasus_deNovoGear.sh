##### Author: Megan Russell #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_deNovoGear"
#SBATCH --time=0-96:00:00
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
echo "perf stat ${SAMTOOLSPATH}/samtools mpileup -gDf ${REF} --bam-list ${BAMFILE} | ${DENOVOPATH}/denovogear dnm auto --ped ${PED} --bcf - 2> ${OUTTRACKNAME}.denovogear.perfOut"
perf stat ${SAMTOOLSPATH}/samtools mpileup -gDf ${REF} ${BAMFILE} | ${DENOVOPATH}/denovogear dnm auto --ped ${PED} --bcf - 2> ${OUTTRACKNAME}.denovogear.perfOut
#${BCFTOOLSPATH}/bcftools convert ${VCF} -Ou -o ${OUTTRACKNAME}.hc.bcf
#if [ $? -eq 0 ] ; then
	#echo "hc vcf converted to a bcf OK"
	#perf stat ${DENOVOPATH}/denovogear dnm auto --ped ${PED} --bcf ${OUTTRACKNAME}.hc.bcf --output_vcf ${OUTTRACKNAME}.denovogear.vcf 2> ${OUTTRACKNAME}.denovogear.perfOut
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