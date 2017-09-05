##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_REVseurat"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 4
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.out"
#SBATCH --error="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.err"

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### TRK: ${TRK}"
echo "### CHRLIST: ${CHRLIST}"
echo "### STEP: ${STEP}"
echo "### STEPCOUNT: ${STEPCOUNT}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### NXT2: ${NXT2}"
echo "### NXT3: ${NXT3}"
echo "### GATKPATH: ${GATKPATH}"
echo "### SEURATPATH: ${SEURATPATH}"

echo "### Seurat caller started for bams at $time."
perf stat java -Djava.io.tmpdir=/scratch/tgenjetstream/tmp/ -jar -Xmx8g ${SEURATPATH}/Seurat.jar \
-T Seurat \
-l INFO \
-R ${REF} \
-I:dna_normal ${TUMOR} \
-I:dna_tumor ${NORMAL} \
--both_strands \
-L ${CHRLIST}/Step${STEP}.list \
--metrics \
--indels \
--allele_metrics \
-o ${TRK}_Step${STEP}.Seurat.vcf \
-go ${TRK}_Step${STEP}.perChr.Seurat.txt \
--pileup_info > ${TRK}_Step${STEP}.REVseuratOut 2> ${TRK}_Step${STEP}.REVseurat.perfOut

if [ $? -eq 0 ] ; then
	#Clean-up the produced VCF to exclude lines where the REF and ALT are identical
	grep "#" ${TRK}_Step${STEP}.Seurat.vcf > ${TRK}_Step${STEP}.Seurat.header
	grep -v "#" ${TRK}_Step${STEP}.Seurat.vcf | awk '{if($4 != $5) print $0}' > ${TRK}_Step${STEP}.Seurat.calls

	cat ${TRK}_Step${STEP}.Seurat.header ${TRK}_Step${STEP}.Seurat.calls > ${TRK}_Step${STEP}.Seurat.vcf

	#Clean-up directory to remove temp files
	rm ${TRK}_Step${STEP}.Seurat.header 
	rm ${TRK}_Step${STEP}.Seurat.calls

	echo "${STEP} Completed" >> ${TRK}_REVseuratStatus.txt
	PROGRESS=`wc -l ${TRK}_REVseuratStatus.txt | awk '{print $1}'`
	mv ${TRK}_Step${STEP}.REVseuratOut ${TRK}_Step${STEP}.REVseuratPass
	touch ${RUNDIR}/${NXT1}
else	
	mv ${TRK}_Step${STEP}.REVseuratOut ${TRK}_Step${STEP}.REVseuratFail
	rm -f ${TRK}_Step${STEP}.REVseuratInQueue
	exit
fi
vcfList=""
#here we make a look to create the list of vcfs based on STEPCOUNT
for i in `seq 1 ${STEPCOUNT}`;
do
        thisVcf="-V ${TRK}_Step$i.Seurat.vcf "
        vcfList="$vcfList $thisVcf"

done
#IF the progress count equals the step count merge to single vcf
if [ ${PROGRESS} -eq ${STEPCOUNT} ]
then
	echo SeuratCaller_${STEP}.Done
	#Concatenate VCF with GATK
 	java -cp ${GATKPATH}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R ${REF} $vcfList -out ${TRK}.REVseurat.vcf -assumeSorted
 	#java -cp ${GATKPATH}/GenomeAnalysisTK.jar org.broadinstitute.sting.tools.CatVariants \
	#	-R ${REF} \
	#	-V ${TRK}_Step1.Seurat.vcf \
	#	-V ${TRK}_Step2.Seurat.vcf \
	#	-V ${TRK}_Step3.Seurat.vcf \
	#	-V ${TRK}_Step4.Seurat.vcf \
	#	-V ${TRK}_Step5.Seurat.vcf \
	#	-V ${TRK}_Step6.Seurat.vcf \
	#	-V ${TRK}_Step7.Seurat.vcf \
	#	-V ${TRK}_Step8.Seurat.vcf \
	#	-V ${TRK}_Step9.Seurat.vcf \
	#	-V ${TRK}_Step10.Seurat.vcf \
	#	-V ${TRK}_Step11.Seurat.vcf \
	#	-V ${TRK}_Step12.Seurat.vcf \
	#	-V ${TRK}_Step13.Seurat.vcf \
	#	-V ${TRK}_Step14.Seurat.vcf \
	#	-V ${TRK}_Step15.Seurat.vcf \
	#	-V ${TRK}_Step16.Seurat.vcf \
	#	-V ${TRK}_Step17.Seurat.vcf \
	#	-V ${TRK}_Step18.Seurat.vcf \
	#	-V ${TRK}_Step19.Seurat.vcf \
	#	-V ${TRK}_Step20.Seurat.vcf \
	#	-V ${TRK}_Step21.Seurat.vcf \
	#	-V ${TRK}_Step22.Seurat.vcf \
	#	-V ${TRK}_Step23.Seurat.vcf \
	#	-V ${TRK}_Step24.Seurat.vcf \
	#	-out ${TRK}.seurat.vcf \
	#	-assumeSorted
		if [ $? -eq 0 ] ; then
			touch ${TRK}.REVseuratPass
			touch ${RUNDIR}/${NXT1}
			touch ${RUNDIR}/${NXT2}
		else
			touch ${TRK}.REVseuratFail
		fi
		mv ${TRK}_REVseuratStatus.txt ${TRK}_REVseuratStatus.txt.used
else
	echo
	echo SeuratCaller_${STEP}.Done
fi

rm -f ${TRK}_Step${STEP}.REVseuratInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:Seurat:$hours:$mins" > ${TRK}_Step${STEP}.REVseurat.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "Seuratcaller finished at $time."
