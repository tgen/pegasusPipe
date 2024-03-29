#!/usr/bin/env bash
#SBATCH --job-name="pegasus_cnaExo15"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=jetstream@tgen.org

module load perl/5.24.1
module load MCR/8.3
module load R/3.4.1 ## NEEDS DNAcopy https://bioconductor.org/packages/release/bioc/html/DNAcopy.html
MCRPATH=/packages/MCR/8.3/

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### NCORES: ${NCORES}"
echo "### CNAPATH: ${CNAPATH}"
echo "### CNAEXOMETARGET: ${CNAEXOMETARGET}"
echo "### GTF: ${GTF}"
echo "### OFILE: ${OFILE}"
echo "### ASSAY: ${ASSAY}"
echo "### NORMALDAT: ${NORMALDAT}"
echo "### TUMORDAT: ${TUMORDAT}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NHETDEPTH: ${NHETDEPTH}"
echo "### THETDEPTH: ${THETDEPTH}"
echo "### MYPATH: ${MYPATH}"
echo "### MERGEDVCF: ${MERGEDVCF}"
echo "### CCDSLIST: ${CCDSLIST}"
echo "### NORMALSAMPLE: ${NORMALSAMPLE}"
echo "### TUMORSAMPLE: ${TUMORSAMPLE}"

##echo "### ZTABLE: ${ZTABLE}"

echo "Starting cna at ${MYPATH} at $time."
cd ${MYPATH}

## add error catching for all of the steps
##parse merged VCF
THETDEPTHINT=${THETDEPTH%.*}
if [ ${THETDEPTHINT} -gt 49 ] ; then
    THETDEPTHCALC=50
else
    #send low coverage email
    echo -e "Coverage for your sample: ${RUNDIR} was lower than 50.  Copy number will use HS metrics calculated target mean for read depth on filtering SNPS" | mail -s "Pegasus Pipeline: Copy Number Alert" jetstream@tgen.org
    #touch file of low cov
    touch ${TRACKNAME}.cnaLOWCOVERAGEALERT.txt
    THETDEPTHCALC=${THETDEPTH}
fi
echo "### THETDEPTHCALC: ${THETDEPTHCALC}"

echo "### Start parseMergeVCF.pl"
echo "${CNAPATH}/parseMergeVCF.pl ${MERGEDVCF} ${NORMALSAMPLE} ${TUMORSAMPLE}"
${CNAPATH}/parseMergeVCF.pl ${MERGEDVCF} ${NORMALSAMPLE} ${TUMORSAMPLE} ${THETDEPTHCALC}
echo "### End parseMergeVCF.pl"
##rename baf.txt
mv ${MYPATH}/baf.txt ${MYPATH}/${OFILE}.baf.txt

##CNA  - Copy Number Analysis (MATLAB)
HETFILE=merged.vcf.txt

smWin=6             #   <<<< THIS CAN BE ADJUSTED - smoothing window size >>>>
fcThresh=0.75       #   <<<< THIS CAN BE ADJUSTED - fold-change threshold - plot >>>>

res=2               #   <<<< THIS CAN BE ADJUSTED - min resolution >>>>
#readDepth=300       #   <<<< THIS CAN BE ADJUSTED - min number of counts >>>>
maxGap=1000   #   <<<< THIS CAN BE ADJUSTED - max distance between consecutive postions >>>>

#hetDepth=100        #   <<<< THIS CAN BE ADJUSTED - min depth of diploid het position >>>>
hetDev=0.025         #   <<<< THIS CAN BE ADJUSTED - allowable deviation from ref allele frequency of 0.5

readDepth=`echo "${NHETDEPTH} * 3" | bc `

echo "### Start run_ngsCNA.sh"
echo "${CNAPATH}/run_ngsCNA.sh ${MCRPATH} ${NORMALDAT} ${TUMORDAT} ${OFILE} ${HETFILE} ${smWin} ${fcThresh} ${ASSAY} ${res} ${readDepth} ${maxGap} ${NHETDEPTH} ${THETDEPTH} ${hetDev} ${CNAEXOMETARGET}"
${CNAPATH}/pegasusCNA/run_ngsCNA.sh ${MCRPATH} ${NORMALDAT} ${TUMORDAT} ${OFILE} ${HETFILE} ${smWin} ${fcThresh} ${ASSAY} ${res} ${readDepth} ${maxGap} ${NHETDEPTH} ${THETDEPTH} ${hetDev} ${CNAEXOMETARGET}
if [ $? -ne 0 ] ; then
    echo "### CNA failed at run_ngsCNA.sh"
    touch ${TRACKNAME}.cnaExomeFail
    rm -f ${TRACKNAME}.cnaExomeInQueue
    exit 1
else
    echo "### Renaming cnaStats files"
    nFileName=`basename ${NORMALDAT}`
    tFileName=`basename ${TUMORDAT}`
    mv ${NORMALDAT}.cnaStats ${RUNDIR}/stats/$nFileName.exo.cnaStats
    mv ${TUMORDAT}.cnaStats ${RUNDIR}/stats/$tFileName.exo.cnaStats
    echo "### Renaming/moving cnaStats file is done."
fi
echo "### End run_ngsCNA.sh"

## CBS segmentation
echo "### Start runDNAcopyExome.R"
Rscript --vanilla ${CNAPATH}/runDNAcopyExomeV2.R ${OFILE}.cna.tsv ${OFILE}.seg
if [ $? -ne 0 ] ; then
    echo "### CNA failed at runDNAcopyExome.R"
    touch ${TRACKNAME}.cnaExomeFail
    rm -f ${TRACKNAME}.cnaExomeInQueue
    exit 1
fi

echo "### End runDNAcopyExome.R"
##plotting

Rscript --vanilla ${CNAPATH}/plotCGH_EXOME.R ${OFILE}.cna.tsv ${OFILE}.amp.tsv ${OFILE}.del.tsv ${OFILE}

if [ -f ${OFILE}.hets.tsv ] ; then
    echo "### Running plots with hets"
    Rscript --vanilla ${CNAPATH}/plotCGHwithHets.R ${OFILE}.cna.tsv ${OFILE}.amp.tsv ${OFILE}.del.tsv ${OFILE}.hets.tsv ${OFILE}_withhets
    if [ $? -ne 0 ] ; then
        echo "### CNA failed at plotCGHwithHets.R"
        touch ${TRACKNAME}.cnaExomeFail
        rm -f ${TRACKNAME}.cnaExomeInQueue
        exit 1
    fi

    echo "### End running plots with hets"
fi

##BAF segmentation and figure
Rscript --vanilla ${CNAPATH}/plotBAF.R ${OFILE}.baf.txt ${OFILE}.baf
Rscript --vanilla ${CNAPATH}/runDNAcopyBAF.R ${OFILE}.baf.txt ${OFILE}.baf
${CNAPATH}/annotBAFseg.pl ${CCDSLIST} ${OFILE}.baf.seg 0.15
if [ $? -ne 0 ] ; then
    echo "### CNA failed at annotBAFSeg.pl"
    touch ${TRACKNAME}.cnaExomeFail
    rm -f ${TRACKNAME}.cnaExomeInQueue
    exit 1
fi

###These may be switched to "optional" if they become a problem
module load MCR/9.0
MCRPATH9=/packages/MCR/9.0/
${CNAPATH}/plotLinearCNA/run_plotLinearCNAandBAF.sh ${MCRPATH9} ${OFILE}.cna.tsv ${OFILE}.baf.txt ${OFILE}.cnaBAF.png
if [ $? -ne 0 ] ; then
    echo "### CNA failed at plotLinearCNAandBAF.sh"
    touch ${TRACKNAME}.cnaExomeFail
    rm -f ${TRACKNAME}.cnaExomeInQueue
    exit 1
fi

${CNAPATH}/plotLinearCNAabsBAF/run_plotLinearCNAandAbsBAF.sh ${MCRPATH9} ${OFILE}.cna.tsv ${OFILE}.baf.txt ${OFILE}.cnaAbsBAF.png
if [ $? -ne 0 ] ; then
    echo "### CNA failed at plotLinearCNAandAbsBAF"
    touch ${TRACKNAME}.cnaExomeFail
    rm -f ${TRACKNAME}.cnaExomeInQueue
    exit 1
fi

##Annotate and convert SEG file to gVCF 
DUPTHRESH=0.58     #   <<<< THIS CAN BE ADJUSTED - Amplification Threshold - log2 fold-change >>>>
DELTHRESH=-0.99    #   <<<< THIS CAN BE ADJUSTED - Deletion Threshold - log2 fold-change >>>>

echo "### Start annotSeg.pl"
${CNAPATH}/annotSeg.pl ${GTF} ${OFILE}.cna.seg ${DUPTHRESH} ${DELTHRESH}
if [ $? -ne 0 ] ; then
    echo "### CNA failed at annotSeg.pl"
    touch ${TRACKNAME}.cnaExomeFail
    rm -f ${TRACKNAME}.cnaExomeInQueue
    exit 1
fi
echo "### End annotSeg.pl"

echo "### START validateCNAVariatsVCF.pl"
${CNAPATH}/validateCNAVariantsVCF.pl ${OFILE}.cna.seg.vcf ${OFILE}.baf.txt ${CNAPATH}/ztable.txt
if [ $? -ne 0 ] ; then
    echo "### CNA failed at validateCNAvariantsVCF.pl"
    touch ${TRACKNAME}.cnaExomeFail
    rm -f ${TRACKNAME}.cnaExomeInQueue
    exit 1
fi
echo "### End validateCNAVariantsVCF.pl"


touch ${TRACKNAME}.cnaExomePass
rm -f ${TRACKNAME}.cnaExomeInQueue
time=`date +%d-%m-%Y-%H-%M`
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:CNAEXO:$hours:$mins" > ${OFILE}.cnaExo.totalTime
echo "Finished cna exo at $time."
