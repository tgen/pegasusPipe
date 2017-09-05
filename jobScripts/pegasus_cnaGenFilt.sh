##### Author: Ahmet Kurdoglu #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_cnaGenFilt"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"

MCRPATH=/packages/MCR/7.14/v714

module load perl
module load MCR/7.14
module load R/3.0.0

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### CNAPATH: ${CNAPATH}"
echo "### GTF: ${GTF}"
echo "### OFILE: ${OFILE}"
echo "### ASSAY: ${ASSAY}"
echo "### NORMALDAT: ${NORMALDAT}"
echo "### TUMORDAT: ${TUMORDAT}"
echo "### RUNDIR: ${RUNDIR}"
echo "### MYPATH: ${MYPATH}"
echo "### MERGEDVCF: ${MERGEDVCF}"

echo "Starting cna at ${MYPATH} at $time."
cd ${MYPATH}

## add error catching for all of the steps
##parse merged VCF
echo "### Start parseMergeVCF.pl"
${CNAPATH}/parseMergeVCF.pl ${MERGEDVCF} ${NORMALSAMPLE} ${TUMORSAMPLE}
echo "### End parseMergeVCF.pl"


##Long-inserts only: aggregate het SNPs to increase confidence
${CNAPATH}/aggregateMergeVCF.pl merged.vcf.txt
mv merged.vcf.txt.ag.txt merged.vcf.txt 



##CNA  - Copy Number Analysis (MATLAB)
HETFILE=merged.vcf.txt

smWin=10             #   <<<< THIS CAN BE ADJUSTED - smoothing window size >>>>
fcThresh=0.21       #   <<<< THIS CAN BE ADJUSTED - fold-change threshold - plot >>>>

res=10               #   <<<< THIS CAN BE ADJUSTED - min resolution >>>>
readDepth=100       #   <<<< THIS CAN BE ADJUSTED - min number of counts >>>>
maxGap=100000000   #   <<<< THIS CAN BE ADJUSTED - max distance between consecutive postions >>>>

nhetDepth=0     #   <<<< THIS CAN BE ADJUSTED - min depth of diploid het position >>>>
thetDepth=0     #   <<<< THIS CAN BE ADJUSTED - min depth of diploid het position >>>>

hetDev=0.075          #   <<<< THIS CAN BE ADJUSTED - allowable deviation from ref allele frequency of 0.5


genomeTarget=/home/tgenref/pipeline_v0.4/cna/pos35.txt
echo "### Start run_ngsCNA.sh"
echo "perf stat ${CNAPATH}/run_ngsCNA.sh ${MCRPATH} ${NORMALDAT} ${TUMORDAT} ${OFILE} ${HETFILE} ${smWin} ${fcThresh} ${ASSAY} ${res} ${readDepth} ${maxGap} ${nhetDepth} ${hetDev} $genomeTarget 2> ${OFILE}.runNgsCna.perfOut "
perf stat ${CNAPATH}/run_ngsCNA.sh ${MCRPATH} ${NORMALDAT} ${TUMORDAT} ${OFILE} ${HETFILE} ${smWin} ${fcThresh} ${ASSAY} ${res} ${readDepth} ${maxGap} ${nhetDepth} ${thetDepth} ${hetDev} $genomeTarget 2> ${OFILE}.runNgsCna.perfOut 
if [ $? -ne 0 ] ; then
	echo "### CNA failed at run_ngsCNA.sh"
	touch ${TRACKNAME}.cnaGenFiltFail
	rm -f ${TRACKNAME}.cnaGenFiltInQueue
	exit
else
	echo "### Renaming cnaStats files"
	nFileName=`basename ${NORMALDAT}`
	tFileName=`basename ${TUMORDAT}`
	mv ${NORMALDAT}.cnaStats ${RUNDIR}/stats/$nFileName.genFilt.cnaStats
	mv ${TUMORDAT}.cnaStats ${RUNDIR}/stats/$tFileName.genFilt.cnaStats
	echo "### Renaming/moving cnaStats file is done."
fi
echo "### End run_ngsCNA.sh"

## CBS segmentation
echo "### Start runDNAcopy.R"
Rscript --vanilla ${CNAPATH}/runDNAcopyV2.R ${OFILE}.cna.tsv ${OFILE}.seg
echo "### End runDNAcopy.R"
##plotting

echo "### Start plotCGH.R"
Rscript --vanilla ${CNAPATH}/plotCGH.R ${OFILE}.cna.tsv ${OFILE}.amp.tsv ${OFILE}.del.tsv ${OFILE}	
if [ $? -ne 0 ] ; then
	echo "### CNA failed at plotCGH.R"
	touch ${TRACKNAME}.cnaGenFiltFail
	rm -f ${TRACKNAME}.cnaGenFiltInQueue
	exit
fi
echo "### End plotCGH.R"

if [ -e ${OFILE}.hets.tsv ] ; then
	echo "### Running plots with hets"
	Rscript --vanilla ${CNAPATH}/plotCGHwithHets.R ${OFILE}.cna.tsv ${OFILE}.amp.tsv ${OFILE}.del.tsv ${OFILE}.hets.tsv ${OFILE}_withhets
	if [ $? -ne 0 ] ; then
		echo "### CNA failed at plotCGHwithHets.R"
		touch ${TRACKNAME}.cnaGenFiltFail
		rm -f ${TRACKNAME}.cnaGenFiltInQueue
		exit
	fi
	echo "### End running plots with hets"
fi

##Annotate and convert SEG file to gVCF 
DUPTHRESH=0.58     #   <<<< THIS CAN BE ADJUSTED - Amplification Threshold - log2 fold-change >>>>
DELTHRESH=-0.99    #   <<<< THIS CAN BE ADJUSTED - Deletion Threshold - log2 fold-change >>>>

echo "### Running annotSeg.pl"
${CNAPATH}/annotSeg.pl ${GTF} ${OFILE}.cna.seg ${DUPTHRESH} ${DELTHRESH}
if [ $? -ne 0 ] ; then
	echo "### CNA failed at annotSeg.pl"
	touch ${TRACKNAME}.cnaGenFiltFail
	rm -f ${TRACKNAME}.cnaGenFiltInQueue
	exit
fi
echo "### End running annotSeg.pl"


touch ${TRACKNAME}.cnaGenFiltPass
rm -f ${TRACKNAME}.cnaGenFiltInQueue
time=`date +%d-%m-%Y-%H-%M`
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:CNAGENFILT:$hours:$mins" > ${OFILE}.cnaGenFilt.totalTime
echo "Finished cna gen filt at $time."
