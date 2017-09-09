#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_strelka"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"
 
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### NXT1: ${NXT1}"
echo "### WD: ${WD}"
echo "### TUMOR: ${TUMOR}"
echo "### NORMAL: ${NORMAL}"
echo "### STRELKAPATH: ${STRELKAPATH}"

justName=`basename ${WD}`
echo "### Step 1, moving to work dir"
cd ${WD}
if [ $? -ne 0 ] ; then
	touch ${WD}.strelkaFail
	rm -f ${WD}.strelkaInQueue
	exit
fi

echo "### Step 2, copying config file here"
cp ${STRELKAPATH}/etc/strelka_config_bwa_default.ini config.ini
if [ $? -ne 0 ] ; then
	touch ${WD}.strelkaFail
	rm -f ${WD}.strelkaInQueue
	exit
fi

if [ ${ASSAY} == "Exome" ] ; then
	echo "### Exome detected. Changing skipdepthfilter to 1"
	sed -i 's/isSkipDepthFilters = 0/isSkipDepthFilters = 1/' config.ini
fi
echo "### Step 3, configure"
${STRELKAPATH}/bin/configureStrelkaWorkflow.pl \
	--normal=${NORMAL} \
	--tumor=${TUMOR} \
	--ref=${REF} \
	--config=config.ini --output-dir=./myAnalysis
if [ $? -ne 0 ] ; then
	touch ${WD}.strelkaFail
	rm -f ${WD}.strelkaInQueue
	exit
fi

echo "### Step4, Run analysis"
cd ./myAnalysis
perf stat make -j8 > ${WD}.strelkaOut 2> ${WD}.strelka.perfOut
if [ $? -eq 0 ] ; then
	mv ${WD}.strelkaOut ${WD}.strelkaPass
	mv ${WD}/myAnalysis/results/passed.somatic.snvs.vcf ${WD}/myAnalysis/results/$justName.strelka.passed.somatic.snvs.vcf
	mv ${WD}/myAnalysis/results/all.somatic.snvs.vcf ${WD}/myAnalysis/results/$justName.strelka.all.somatic.snvs.vcf
	mv ${WD}/myAnalysis/results/passed.somatic.indels.vcf ${WD}/myAnalysis/results/$justName.strelka.passed.somatic.indels.vcf
	mv ${WD}/myAnalysis/results/all.somatic.indels.vcf ${WD}/myAnalysis/results/$justName.strelka.all.somatic.indels.vcf
	touch ${RUNDIR}/${NXT1}
else
	mv ${WD}.strelkaOut ${WD}.strelkaFail
fi
rm -f ${WD}.strelkaInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:STRELKA:$hours:$mins" > ${WD}.strelka.totalTime
