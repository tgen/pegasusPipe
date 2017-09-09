#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_hc"
#SBATCH --time=0-170:00:00
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
echo "### CHRLIST: ${CHRLIST}"
echo "### STEP: ${STEP}"
echo "### STEPCOUNT: ${STEPCOUNT}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### GATKPATH: ${GATKPATH}"
echo "### KNOWN: ${KNOWN}"
echo "### BAMLIST: ${BAMLIST}"

echo "### Haplotype caller started for multiple bams at $time."
perf stat java -Djava.io.tmpdir=/scratch/tgenjetstream/tmp/ -jar -Xmx24g ${GATKPATH}/GenomeAnalysisTK.jar \
-l INFO \
-R ${REF} \
-L ${CHRLIST}/Step${STEP}.list \
-nct 8 \
-T HaplotypeCaller \
${BAMLIST} \
-D ${KNOWN} \
-mbq 10 \
-o ${TRK}_Step${STEP}.HC.vcf > ${TRK}_Step${STEP}.hcOut 2> ${TRK}_Step${STEP}.hapCal.perfOut
if [ $? -eq 0 ] ; then
	echo "${STEP} Completed" >> ${TRK}_hcStatus.txt
	PROGRESS=`wc -l ${TRK}_hcStatus.txt | awk '{print $1}'`
	mv ${TRK}_Step${STEP}.hcOut ${TRK}_Step${STEP}.hcPass
else	
	mv ${TRK}_Step${STEP}.hcOut ${TRK}_Step${STEP}.hcFail
	rm -f ${TRK}_Step${STEP}.hcInQueue
	exit
fi

vcfList=""
#here we make a look to create the list of vcfs based on STEPCOUNT
for i in `seq 1 ${STEPCOUNT}`;
do
        thisVcf="-V ${TRK}_Step$i.HC.vcf "
        vcfList="$vcfList $thisVcf"
done

#IF the progress count equals the step count merge to single vcf
if [ ${PROGRESS} -eq ${STEPCOUNT} ]
then
	echo "PROGRESS: ${PROGRESS} equals STEPCOUNT: ${STEPCOUNT}"
	echo HapCaller_${STEP}.Done
	#Concatenate VCF with GATK
 	java -cp ${GATKPATH}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R ${REF} $vcfList -out ${TRK}.HC_All.vcf -assumeSorted
	#java -cp ${GATKPATH}/GenomeAnalysisTK.jar org.broadinstitute.sting.tools.CatVariants \
	#	-R ${REF} \
	#	-V ${TRK}_Step1.HC.vcf \
	#	-V ${TRK}_Step2.HC.vcf \
	#	-V ${TRK}_Step3.HC.vcf \
	#	-V ${TRK}_Step4.HC.vcf \
	#	-V ${TRK}_Step5.HC.vcf \
	#	-V ${TRK}_Step6.HC.vcf \
	#	-V ${TRK}_Step7.HC.vcf \
	#	-V ${TRK}_Step8.HC.vcf \
	#	-V ${TRK}_Step9.HC.vcf \
	#	-V ${TRK}_Step10.HC.vcf \
	#	-V ${TRK}_Step11.HC.vcf \
	#	-V ${TRK}_Step12.HC.vcf \
	#	-V ${TRK}_Step13.HC.vcf \
	#	-V ${TRK}_Step14.HC.vcf \
	#	-V ${TRK}_Step15.HC.vcf \
	#	-V ${TRK}_Step16.HC.vcf \
	#	-V ${TRK}_Step17.HC.vcf \
	#	-V ${TRK}_Step18.HC.vcf \
	#	-V ${TRK}_Step19.HC.vcf \
	#	-V ${TRK}_Step20.HC.vcf \
	#	-V ${TRK}_Step21.HC.vcf \
	#	-V ${TRK}_Step22.HC.vcf \
	#	-V ${TRK}_Step23.HC.vcf \
	#	-V ${TRK}_Step24.HC.vcf \
	#	-out ${TRK}.HC_All.vcf \
	#	-assumeSorted
		if [ $? -eq 0 ] ; then
			touch ${TRK}.hcPass
			touch ${RUNDIR}/${NXT1}
			touch ${RUNDIR}/${NXT2}
			touch ${RUNDIR}/${NXT3}
		else
			touch ${TRK}.hcFail
		fi
		mv ${TRK}_hcStatus.txt ${TRK}_hcStatus.txt.used
else
	echo
	echo HapCaller_${STEP}.Done
fi

rm -f ${TRK}_Step${STEP}.hcInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:GATKHC:$hours:$mins" > ${TRK}_Step${STEP}.hapCal.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "Haplotypecaller finished at $time."
