#!/usr/bin/env bash
#SBATCH --job-name="pegasus_cnaGenUnfi"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org


module load VCFtools/0.1.10
module load perl
module load MCR/8.3
module load R/3.0.0
MCRPATH=/packages/MCR/8.3/v83

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
##${CNAPATH}/aggregateMergeVCF.pl merged.vcf.txt 
##mv merged.vcf.txt.ag.txt merged.vcf.txt

##CNA  - Copy Number Analysis (MATLAB)
HETFILE=merged.vcf.txt

smWin=10             #   <<<< THIS CAN BE ADJUSTED - smoothing window size >>>>
fcThresh=0.21       #   <<<< THIS CAN BE ADJUSTED - fold-change threshold - plot >>>>

res=10               #   <<<< THIS CAN BE ADJUSTED - min resolution >>>>
readDepth=100       #   <<<< THIS CAN BE ADJUSTED - min number of counts >>>>
maxGap=100000000   #   <<<< THIS CAN BE ADJUSTED - max distance between consecutive postions >>>>

nhetDepth=0     #   <<<< THIS CAN BE ADJUSTED - min depth of diploid het position >>>>
thetDepth=0     #   <<<< THIS CAN BE ADJUSTED - min depth of diploid het position >>>>
hetDev=0.05          #   <<<< THIS CAN BE ADJUSTED - allowable deviation from ref allele frequency of 0.5

#genomeTarget=/home/tgenref/pipeline_v0.4/cna/pos35.txt
echo "### Start run_ngsCNA.sh"
echo "${CNAPATH}/pegasusCNA/run_ngsCNA.sh ${MCRPATH} ${NORMALDAT} ${TUMORDAT} ${OFILE} ${HETFILE} ${smWin} ${fcThresh} ${ASSAY} ${res} ${readDepth} ${maxGap} ${nhetDepth} ${thetDepth} ${hetDev}"
${CNAPATH}/pegasusCNA/run_ngsCNA.sh ${MCRPATH} ${NORMALDAT} ${TUMORDAT} ${OFILE} ${HETFILE} ${smWin} ${fcThresh} ${ASSAY} ${res} ${readDepth} ${maxGap} ${nhetDepth} ${thetDepth} ${hetDev}
if [ $? -ne 0 ] ; then
	echo "### CNA failed at run_ngsCNA.sh"
	touch ${TRACKNAME}.cnaGenUnfiFail
	rm -f ${TRACKNAME}.cnaGenUnfiInQueue
	exit
else
	echo "### Renaming cnaStats files"
	nFileName=`basename ${NORMALDAT}`
	tFileName=`basename ${TUMORDAT}`
	mv ${NORMALDAT}.cnaStats ${RUNDIR}/stats/$nFileName.genUnfi.cnaStats
	mv ${TUMORDAT}.cnaStats ${RUNDIR}/stats/$tFileName.genUnfi.cnaStats
	echo "### Renaming/moving cnaStats file is done."
fi

echo "### End run_ngsCNA.sh"

## CBS segmentation
echo "### Start runDNAcopy.R"
Rscript --vanilla ${CNAPATH}/runDNAcopyV2.R ${OFILE}.cna.tsv ${OFILE}.seg
echo "### End runDNAcopy.R"
##plotting

Rscript --vanilla ${CNAPATH}/plotCGH.R ${OFILE}.cna.tsv ${OFILE}.amp.tsv ${OFILE}.del.tsv ${OFILE}	

if [ -f ${OFILE}.hets.tsv ] ; then
	echo "### Running plots with hets"
	Rscript --vanilla ${CNAPATH}/plotCGHwithHets.R ${OFILE}.cna.tsv ${OFILE}.amp.tsv ${OFILE}.del.tsv ${OFILE}.hets.tsv ${OFILE}_withhets
	echo "### End running plots with hets"
fi

echo "### Running plotBAF.R"
Rscript --vanilla ${CNAPATH}/plotBAF.R baf.txt ${OFILE}.baf
if [ $? -ne 0 ] ; then
        echo "### CNA failed at plotBAF.R"
        touch ${TRACKNAME}.cnaGenFiltFail
        rm -f ${TRACKNAME}.cnaGenFiltInQueue
        exit
fi
echo "### End running plotBAF.R"

echo "### Running runDNAcopyBAF.R"
Rscript --vanilla ${CNAPATH}/runDNAcopyBAF.R baf.txt ${OFILE}.baf
if [ $? -ne 0 ] ; then
        echo "### CNA failed at runDNAcopyBAF.R"
        touch ${TRACKNAME}.cnaGenFiltFail
        rm -f ${TRACKNAME}.cnaGenFiltInQueue
        exit
fi
echo "### End running runDNAcopyBAF.R"


##Annotate and convert SEG file to gVCF 
DUPTHRESH=0.58     #   <<<< THIS CAN BE ADJUSTED - Amplification Threshold - log2 fold-change >>>>
DELTHRESH=-0.99    #   <<<< THIS CAN BE ADJUSTED - Deletion Threshold - log2 fold-change >>>>

${CNAPATH}/annotSeg.pl ${GTF} ${OFILE}.cna.seg ${DUPTHRESH} ${DELTHRESH}
if [ $? -ne 0 ] ; then
        echo "### CNA failed at annotSeg.pl"
        touch ${TRACKNAME}.cnaGenFiltFail
        rm -f ${TRACKNAME}.cnaGenFiltInQueue
        exit
fi
echo "### End running annotSeg.pl"

echo "### Running validateCNAVariantsVCF.pl"
${CNAPATH}/validateCNAVariantsVCF.pl ${OFILE}.cna.seg.vcf baf.txt ${CNAPATH}/ztable.txt
if [ $? -ne 0 ] ; then
        echo "### CNA failed at validateCNAVariantsVCF.pl"
        touch ${TRACKNAME}.cnaGenFiltFail
        rm -f ${TRACKNAME}.cnaGenFiltInQueue
        exit
fi
echo "### End running validateCNAVariantsVCF.pl"

touch ${TRACKNAME}.cnaGenUnfiPass
rm -f ${TRACKNAME}.cnaGenUnfiInQueue
time=`date +%d-%m-%Y-%H-%M`
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:CNAGENUNFILT:$hours:$mins" > ${OFILE}.cnaGenUnFilt.totalTime
echo "Finished cna gen filt at $time."
