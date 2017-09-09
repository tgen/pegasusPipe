#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_trn"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### NORMAL: ${NORMAL}"
echo "### BEDFILE: ${BEDFILE}"
echo "### GTF: ${GTF}"
echo "### ASSAY: ${ASSAY}"
echo "### TUMOR: ${TUMOR}"
echo "### SAMPATH: ${SAMPATH}"
echo "### NXT1: ${NXT1}"
echo "### OUTFILE: ${OUTFILE}"
echo "### RECIPE: ${RECIPE}"

cd ${RUNDIR}
echo "Starting translocations for ${OUTFILE} at $time."
if [ "${RECIPE}" == "trial01" ] ; then
echo "Running new strexome version of TRN 02/06/15"
perf stat ${TRNPATH}/tgen_somaticSV.pl \
        T=${TUMOR} \
        R=${NORMAL} \
        O=${OUTFILE} \
        S=${SAMPATH}/samtools > ${OUTFILE}.trnOut 2> ${OUTFILE}.trn.perfOut
else

	if [ "${ASSAY}" == "Genome" ] ; then
	echo "Running as a Genome trn"
	perf stat ${TRNPATH}/tgen_somaticSV.pl \
		T=${TUMOR} \
		R=${NORMAL} \
		O=${OUTFILE} \
		G=${GTF} \
		S=${SAMPATH}/samtools > ${OUTFILE}.trnOut 2> ${OUTFILE}.trn.perfOut

	elif [ "${ASSAY}" == "Exome" ] ; then 
	echo "Running as Exome trn"
	perf stat ${TRNPATH}/tgen_somaticSV.pl \
		T=${TUMOR} \
		R=${NORMAL} \
		O=${OUTFILE} \
		G=${GTF} \
		F=${BEDFILE} \
		S=${SAMPATH}/samtools > ${OUTFILE}.trnOut 2> ${OUTFILE}.trn.perfOut
	else
		echo "### I should not be here"
		touch "some fail file"
	fi
fi 
if [ $? -eq 0 ] ; then
	mv ${OUTFILE}.trnOut ${OUTFILE}.trnPass
	touch ${RUNDIR}/${NXT1}
else	
	mv ${OUTFILE}.trnOut ${OUTFILE}.trnFail
fi
rm -f ${OUTFILE}.trnInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:TRNS:$hours:$mins" > ${OUTFILE}.trnslocs.totalTime
time=`date +%d-%m-%Y-%H-%M`
echo "Finished trnations at $time."
