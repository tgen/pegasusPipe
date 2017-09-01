##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_circos"
#SBATCH --time=0-96:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#SBATCH --output="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.out"
#SBATCH --error="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.out"

time=`date +%d-%m-%Y-%H-%M` 
beginTime=`date +%s`
machine=`hostname`

echo "starting circos on this machine: $machine for $outDir"
	
echo "###OUTFILE is: ${OUTFILE}"
echo "###OUTDIR is: ${OUTDIR}"
echo "###SEURATVCF is ${SEURATVCF}"
echo "###TRNVCF is ${TRNVCF}"
echo "###CNVTSV is ${CNVTSV}"
echo "###RUNDIR is ${RUNDIR}"
echo "###COSMIC is ${COSMIC}"
echo "###CONF is ${CONF}"

cnvsFile=${CNVTSV}
transFile=${TRNVCF}
snvsFile=${SEURATVCF}
cosmic=1
outDir=${OUTDIR}

module add perl
cd $outDir
mkdir ${OUTDIR}/data
/home/akurdoglu/pnap-pipe/transExtract.sh $transFile $outDir
/home/akurdoglu/pnap-pipe/cnvsExtract.sh $cnvsFile $outDir
/home/akurdoglu/pnap-pipe/snvsExtract.sh $snvsFile $outDir $cosmic

/home/akurdoglu/pnap-pipe/circosPrep.sh $outDir


/packages/circos-0.62-1/bin/circos -conf ${CONF} 
if [ $? -eq 0 ] ; then
        touch ${OUTFILE}.circosPass
        #grep "${RUNDIR}" ${MSG1} > /dev/null
        #if [ $? -eq 0 ] ; then
        #        echo "${RUNDIR} already in ${MSG1}. Not appending."
        #else
        #        echo "${RUNDIR}" >> ${MSG1}
        #fi
else
        touch ${OUTFILE}.circosFail
fi
rm ${OUTFILE}.circosInQueue
cd -
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:CIRCOS:$hours:$mins" > ${OUTFILE}.circos.totalTime
