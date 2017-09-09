#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_seuratSWAPCheck"
#SBATCH --time=0-16:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"
 
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### DBVERSION: ${DBVERSION}"
echo "### VCF: ${VCF}"
echo "### SNPEFFPATH: ${SNPEFFPATH}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"

echo "### Starting swap check"
#############################################################


#    my $file=shift;
    threshold=20
    #if((index($assay, "A1STX")!=-1) || (index($assay, "A1STL")!=-1)){
     #   $threshold = 20;
     #   chomp ($file);
     #   print "\n\nFinalizing STX/STL VCFs : $file\n";
      #  if ($file=~/(.*).(snpEff.vcf)/) {
if [ $ASSAY == "Exome" ] ; then
	$cleanFile="$1.clean.vcf";
	$snpeffFile="$1.clean.snpeff.vcf";
	$finalFile="$1.snpeff.final.vcf";
	$run=`java -jar $snpEffPath/SnpSift.jar rmInfo $file EFF > $cleanFile`;
	 print "Command Run: java -jar $snpEffPath/SnpSift.jar rmInfo $file EFF to $cleanFile\n";
	$run=`java -Xmx4G -jar $snpEffPath/snpEff.jar -c $snpEffPath/snpEff.config -no-downstream  -no-intergenic -no-intron -no-upstream -canon -hgvs GRCh37.74 -o vcf $cleanFile > $snpeffFile`;
	$run=`cat $snpeffFile | java -jar $snpEffPath/SnpSift.jar filter "(AR2 > 0.05 ) & (AR1 < 0.02) & (DP2 > $threshold) & (DP1 > $threshold) & (QUAL > 15) & (FILTER = 'PASS')" | java -jar $snpEffPath/SnpSift.jar intervals $intervals | java -jar $snpEffPath/SnpSift.jar filter "EFF[*].BIOTYPE = 'protein_coding'" > $finalFile`;
	$run=`rm $cleanFile`;
	$run=`rm $snpeffFile`;
       
elif [ $ASSAY == "Genome" ] ; then 
	$cleanFile="$1.clean.vcf";
	$snpeffFile="$1.clean.snpeff.vcf";
	$finalFile="$1.snpeff.final.vcf";
	$run=`java -jar $snpEffPath/SnpSift.jar rmInfo $file EFF > $cleanFile`;
	$run=`java -Xmx4G -jar $snpEffPath/snpEff.jar -c $snpEffPath/snpEff.config -no-downstream  -no-intergenic -no-intron -no-upstream -canon -hgvs GRCh37.74 -o vcf $cleanFile > $snpeffFile`;
	$run = `mv  $snpeffFile $finalFile`;
	$run=`rm $cleanFile`;
	$run=`rm $snpeffFile`;
else
	echo "I should not be here, not exome or genome"

#############################################################
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
	-s ${summaryOut} \
	-c ${SNPEFFPATH}/snpEff.config \
	${DBVERSION} \
	${VCF} > $snpEffTxt
java -Xmx6g -jar ${SNPEFFPATH}/snpEff.jar eff \
	-v \
	-i vcf \
	-o vcf \
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
