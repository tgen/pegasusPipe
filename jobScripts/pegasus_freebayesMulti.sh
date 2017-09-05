##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_freebayesMulti"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### KNOWN: ${KNOWN}"
echo "### STEP: ${STEP}"
echo "### FBBAM: ${FBBAM}"
echo "### FREEBAYESPATH: ${FREEBAYESPATH}"
echo "### TRACKNAME: ${TRACKNAME}"
echo "### CHRLIST: ${CHRLIST}"
echo "### BAMLIST: ${BAMLIST}"

##REGION=`cat ${CHRLIST}/Step${STEP}.list`

echo "### freebayes started at $time."
perf stat ${FREEBAYESPATH}/freebayes -f ${REF} ${BAMLIST} -t ${CHRLIST}/Step${STEP}.bed --ploidy 2 --min-repeat-entropy 1 > ${TRACKNAME}_Step${STEP}.freebayes.vcf 2> ${TRACKNAME}_Step${STEP}.freebayes.perfOut
if [ $? -eq 0 ] ; then
        echo "${STEP} Completed" >> ${TRACKNAME}_fbStatus.txt
        PROGRESS=`wc -l ${TRACKNAME}_fbStatus.txt | awk '{print $1}'`
	touch ${TRACKNAME}_Step${STEP}.freebayesPass
else	
	touch ${TRACKNAME}_Step${STEP}.freebayesFail
	rm -f ${TRACKNAME}_Step${STEP}.freebayesInQueue
	exit
fi
vcfList=""
#here we make a look to create the list of vcfs based on STEPCOUNT
for i in `seq 1 ${STEPCOUNT}`;
do
        thisVcf="-V ${TRACKNAME}_Step$i.freebayes.vcf "
        vcfList="$vcfList $thisVcf"
done
#IF the progress count equals the step count merge to single vcf
if [ ${PROGRESS} -eq ${STEPCOUNT} ]
then
        echo Freebayes_${STEP}.Done
	#Concatenate VCF with GATK
	java -cp ${GATKPATH}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R ${REF} $vcfList -out ${TRACKNAME}.freebayes_All.vcf -assumeSorted
		if [ $? -eq 0 ] ; then
                        touch ${TRACKNAME}.freebayesPass
                        touch ${RUNDIR}/${NXT1}
                        touch ${RUNDIR}/${NXT2}
                else
                        touch ${TRACKKNAME}.freebayesFail
                fi
                mv ${TRACKNAME}_fbStatus.txt ${TRACKNAME}_fbStatus.txt.used
else
	echo
        echo HapCaller_${STEP}.Done
fi


rm -f ${TRACKNAME}_Step${STEP}.freebayesInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:FREEBAYES:$hours:$mins" > ${TRACKNAME}_Step${STEP}.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "Freebayes finished at $time."
