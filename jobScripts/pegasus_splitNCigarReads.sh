##### Author: Megan Russell #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_splitNCigar"
#SBATCH --time=0-48:00:00
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
echo "### RNABAM: ${RNABAM}"
echo "### OUTBAM: ${OUTBAM}"

#eecho "### GATK splitNCigarReads started at $time."
#perf stat java -Djava.io.tmpdir=/scratch/tgenjetstream/tmp/ -jar -Xmx32g ${GATKPATH}/GenomeAnalysisTK.jar \
#-l INFO \
#-R ${REF} \
#-T SplitNCigarReads \
#--filter_reads_with_N_cigar \
#-I ${RNABAM} \
#-o ${OUTBAM} > ${RNABAM}.splitNCigarOut 2> ${RNABAM}.splitNCigar.perfOut
echo "### GATK splitNCigarReads started at $time."
##All mapping qualities of 255 will be reassigned to 60 in the bam
#perf stat java -Djava.io.tmpdir=/scratch/tgenjetstream/tmp/ -jar -Xmx44g ${GATKPATH}/GenomeAnalysisTK.jar \
perf stat java -Xmx40G -Djava.io.tmpdir=/scratch/tgenjetstream/tmp/ -jar ${GATKPATH}/GenomeAnalysisTK.jar \
-l INFO \
-R ${REF} \
-T SplitNCigarReads \
-I ${RNABAM} \
-o ${OUTBAM} \
-rf ReassignOneMappingQuality \
-RMQF 255 \
-RMQT 60 \
-U ALLOW_N_CIGAR_READS > ${RNABAM}.splitNCigarOut 2> ${RNABAM}.splitNCigar.perfOut
if [ $? -eq 0 ] ; then
	mv ${RNABAM}.splitNCigarOut ${RNABAM}.splitNCigarPass
	touch ${RUNDIR}/${NXT1}
else	
	mv ${RNABAM}.splitNCigarOut ${RNABAM}.splitNCigarFail
	rm -f ${RNABAM}.splitNCigarInQueue
	exit
fi
rm -f ${RNABAM}.splitNCigarInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:GATKsplitNCigar:$hours:$mins" > ${RNABAM}.splitNCigar.totalTime

time=`date +%d-%m-%Y-%H-%M`
echo "splitNCigar finished at $time."
