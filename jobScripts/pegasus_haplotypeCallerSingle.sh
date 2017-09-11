#!/usr/bin/env bash
#SBATCH --job-name="pegasus_hcSin"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### CHRLIST: ${CHRLIST}"
echo "### STEP: ${STEP}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### GATKPATH: ${GATKPATH}"
echo "### KNOWN: ${KNOWN}"
echo "### BAMLIST: ${BAMLIST}"

echo "### Haplotype caller started for multiple bams at $time."
java -Djava.io.tmpdir=/scratch/tgenjetstream/tmp/ -jar -Xmx24g ${GATKPATH}/GenomeAnalysisTK.jar \
-l INFO \
-R ${REF} \
-T HaplotypeCaller \
-L ${CHRLIST}/Step${STEP}.list \
-nct 8 \
-I ${BAMLIST} \
-D ${KNOWN} \
-mbq 10 \
-o ${TRK}_Step${STEP}.HC.vcf > ${TRK}_Step${STEP}.hcOut

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
	echo HapCaller_${STEP}.Done
	java -cp ${GATKPATH}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R ${REF} $vcfList -out ${TRK}.HC_All.vcf -assumeSorted
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
echo "RUNTIME:GATKHC:$hours:$mins" > ${TRK}.hapCal.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "UnifiedGenotyper finished at $time."
