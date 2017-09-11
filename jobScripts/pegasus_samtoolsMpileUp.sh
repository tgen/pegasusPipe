#!/usr/bin/env bash
#SBATCH --job-name="pegasus_samtoolsMpileUp"
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
echo "### NXT2: ${NXT2}"
echo "### SAMTOOLSPATH: ${SAMTOOLSPATH}"
echo "### BCFTOOLSPATH: ${BCFTOOLSPATH}"
echo "### TRACKNAME: ${TRACKNAME}"
echo "### BAMFILE: ${BAMFILE}"
echo "### D: ${D}"
echo "### BEDFILE: ${BEDFILE}"
echo "### CHRLIST: ${CHRLIST}"
echo "### STEP: ${STEP}"
echo "### STEPCOUNT: ${STEPCOUNT}"
echo "### GATKPATH: ${GATKPATH}"


echo "### Samtools mPileUP started at $time."

${SAMTOOLSPATH}/samtools mpileup \
    -DsSOg \
    -C 50 \
    -F 0.01 \
    -l ${CHRLIST}/Step${STEP}.bed \
    -f ${REF} \
    ${BAMFILE} | \
${BCFTOOLSPATH}/bcftools call \
    -vmO v \
    -o ${TRACKNAME}_Step${STEP}.mpileup.vcf

if [ $? -eq 0 ] ; then
	echo "${STEP} Completed" >> ${TRACKNAME}_spStatus.txt
        PROGRESS=`wc -l ${TRACKNAME}_spStatus.txt | awk '{print $1}'`
	touch ${TRACKNAME}_Step${STEP}.samtoolsMpileUpPass
else	
	touch ${TRACKNAME}_Step${STEP}.samtoolsMpileUpFail
	rm -f ${TRACKNAME}_Step${STEP}.samtoolsMpileUpInQueue
    exit 1
fi

# Here we make a look to create the list of vcfs based on STEPCOUNT
vcfList=""
for i in `seq 1 ${STEPCOUNT}`; do
    thisVcf="-V ${TRACKNAME}_Step$i.mpileup.vcf "
    vcfList="$vcfList $thisVcf"
done

# IF the progress count equals the step count merge to single vcf
if [ ${PROGRESS} -eq ${STEPCOUNT} ]
then
    echo mpileup_${STEP}.Done
	#Concatenate VCF with GATK
	echo "java -cp ${GATKPATH}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants --reference ${REF} $vcfList -out ${TRACKNAME}.mpileup_All.vcf -assumeSorted"
	java -cp ${GATKPATH}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants --reference ${REF} $vcfList -out ${TRACKNAME}.mpileup_All.vcf -assumeSorted
    if [ $? -eq 0 ] ; then
        touch ${TRACKNAME}.samtoolsMpileUpPass
        touch ${RUNDIR}/${NXT1}
        touch ${RUNDIR}/${NXT2}
    else
        touch ${TRACKNAME}.samtoolsMpileUpFail
    fi
    mv ${TRACKNAME}_spStatus.txt ${TRACKNAME}_spStatus.txt.used

else
    echo
    echo mpileup_${STEP}.Done
fi

rm -f ${TRACKNAME}_Step${STEP}.samtoolsMpileUpInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:SAMTOOLSMPILEUP:$hours:$mins" > ${TRACKNAME}_Step${STEP}.samtoolsPileUp.totalTime
time=`date +%d-%m-%Y-%H-%M`
echo "SamtoolsMPileUp finished at $time."
