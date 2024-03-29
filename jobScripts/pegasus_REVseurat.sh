#!/usr/bin/env bash
#SBATCH --job-name="pegasus_REVseurat"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 16G

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
echo "/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.91-2.6.2.3.el7.x86_64/jre/bin/java -Djava.io.tmpdir=$TMPDIR -jar -Xmx8g ${SEURATPATH}/Seurat.jar \
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
    --pileup_info > ${TRK}_Step${STEP}.REVseuratOut"

/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.91-2.6.2.3.el7.x86_64/jre/bin/java -Djava.io.tmpdir=$TMPDIR -jar -Xmx8g ${SEURATPATH}/Seurat.jar \
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
    --pileup_info > ${TRK}_Step${STEP}.REVseuratOut

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
    exit 1
fi

# Here we make a look to create the list of vcfs based on STEPCOUNT
vcfList=""
for i in `seq 1 ${STEPCOUNT}`;
do
    thisVcf="-V ${TRK}_Step$i.Seurat.vcf "
    vcfList="$vcfList $thisVcf"
done

# IF the progress count equals the step count merge to single vcf
if [ ${PROGRESS} -eq ${STEPCOUNT} ]
then
    echo SeuratCaller_${STEP}.Done
    #Concatenate VCF with GATK
     java -cp ${GATKPATH}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R ${REF} $vcfList -out ${TRK}.REVseurat.vcf -assumeSorted

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
