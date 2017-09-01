##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_snpEff"
#SBATCH --time=0-16:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.out"
#SBATCH --error="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.err"
 
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### DBVERSION: ${DBVERSION}"
echo "### VCF: ${VCF}"
echo "### SNPEFFPATH: ${SNPEFFPATH}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"

echo "### Starting snpEff Annotator of vcf file: ${VCF}"
OUT=${VCF/.proj.md.bam}
snpEffOut=${OUT/.vcf/.snpEff.vcf}
snpEffInt=${OUT/.vcf/.snpEffInt.vcf}
snpEffTxt=${OUT/.vcf/.snpEff.txt}
summaryOut=${OUT/.vcf/.snpEff.summary_html}
	##-hgvs \
	##-hgvs \
java -Xmx6g -jar ${SNPEFFPATH}/snpEff.jar eff \
	-v \
	-i vcf \
	-o txt \
	-noLog \
	-s ${summaryOut} \
	-c ${SNPEFFPATH}/snpEff.config \
	${DBVERSION} \
	${VCF} > $snpEffTxt
java -Xmx6g -jar ${SNPEFFPATH}/snpEff.jar eff \
	-v \
	-i vcf \
	-o vcf \
	-noLog \
	-s ${summaryOut} \
	-c ${SNPEFFPATH}/snpEff.config \
	${DBVERSION} \
	${VCF} > $snpEffInt
if [ $? -ne 0 ] ; then
	echo "snpEff first part failed." >> ${VCF}.snpEffOut
	mv ${VCF}.snpEffOut ${VCF}.snpEffFail
else
	echo "snpEff first part complete." >> ${VCF}.snpEffOut
	echo "snpEff second part (snpSift) starting." >> ${VCF}.snpEffOut
	java -Xmx6g -jar ${SNPEFFPATH}/SnpSift.jar annotate \
	${DBSNP} \
	$snpEffInt > $snpEffOut
	if [ $? -eq 0 ] ; then
		echo "snpEff second part (snpSift) complete." >> ${VCF}.snpEffOut
		mv ${VCF}.snpEffOut ${VCF}.snpEffPass
		touch ${RUNDIR}/${NXT1}
	else
		echo "snpEff second part (snpSift) failed." >> ${VCF}.snpEffOut
		mv ${VCF}.snpEffOut ${VCF}.snpEffFail
		rm -f $snpEffOut
	fi
fi
rm -f $snpEffInt
rm -rf ${VCF}.snpEffInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:SNPEFF:$hours:$mins" > ${VCF}.snpeff.totalTime
time=`date +%d-%m-%Y-%H-%M` 
echo "Ending snpEff Annotator."
