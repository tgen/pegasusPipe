##### Author: Megan #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_ancestry"
#SBATCH --time=0-60:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.out"
##PBS -e /${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.err

module load python/2.7.3
module load R/3.2.1

beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### DIR: ${DIR}"
echo "### NXT1: ${NXT1}"
echo "### REF: ${REF}"
echo "### RUNDIR: ${RUNDIR}"
echo "### LASERPATH: ${LASERPATH}"
echo "### HGDPPATH: ${HGDPPATH}"
echo "### SAMTOOLSPATH: ${SAMTOOLSPATH}"
echo "### BAMFILE: ${BAMFILE}"
echo "### BEDFILE: ${BEDFILE}"
echo "### TRACKNAME: ${TRACKNAME}"


cd ${ANCESTRYDIR}
echo "perf stat ${SAMTOOLSPATH}/samtools mpileup -q 30 -Q 20 -f ${REF} -l ${HGDPPATH}/HGDP_938.bed ${BAMFILE} > ${TRACKNAME}.pileup 2> ${TRACKNAME}.mpileup.perfOut"
perf stat ${SAMTOOLSPATH}/samtools mpileup -q 30 -Q 20 -f ${REF} -l ${HGDPPATH}/HGDP_938.bed ${BAMFILE} > ${TRACKNAME}.pileup 2> ${TRACKNAME}.mpileup.perfOut
if [ $? -eq 0 ] ; then
	perf stat python ${LASERPATH}/pileup2seq/pileup2seq.py \
	-m ${HGDPPATH}/HGDP_938.site \
	-o ${TRACKNAME} \
	${TRACKNAME}.pileup 2> ${TRACKNAME}.pileup2seq.perfOut
	if [ $? -eq 0 ] ; then
		perf stat ${LASERPATH}/laser -g ${HGDPPATH}/HGDP_938.geno -c ${HGDPPATH}/HGDP_938.RefPC.coord -s ${TRACKNAME}.seq -K 20 -k 4 -o ${TRACKNAME} 2> ${TRACKNAME}.laser.perfOut
		if [ $? -eq 0 ] ; then
			
			perf stat Rscript ${LASERPATH}/plot/plotHGDPMegan.r ${TRACKNAME}.SeqPC.coord ${HGDPPATH}/HGDP_938.RefPC.coord 2> ${TRACKNAME}.plot.perfOut
			mv ${ANCESTRYDIR}/Ancestry_on_HGDP.pdf ${TRACKNAME}.ancestry_on_HGDP.pdf

			touch ${TRACKNAME}.ancestryPass
		else
			touch ${TRACKNAME}.ancestryFail
		fi
	else
		touch ${TRACKNAME}.ancestryFail
	fi
else
	touch ${TRACKNAME}.ancestryFail
fi


rm ${TRACKNAME}.ancestryInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:ANCESTRY:$hours:$mins" > ${TRACKNAME}.ancestry.totalTime
echo "### Ending LASER ancestry script."
