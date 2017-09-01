##### Author: Megan Russell #####
##### Parameterized PBS Script ####
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_vcfMerger"
#SBATCH --time=0-96:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.out"
#SBATCH --error="/${D}/oeFiles/${PBS_JOBNAME}_${PBS_JOBID}.err"
 
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
#echo "### DBVERSION: ${DBVERSION}"
#echo "### VCF: ${VCF}"
echo "### SNPSIFT: ${SNPSIFT}"
echo "### RUNDIR: ${RUNDIR}"
echo "### NXT1: ${NXT1}"
echo "### SEURAT_VCF: ${SEURAT_VCF}"
echo "### MUTECT_VCF: ${MUTECT_VCF}"
echo "### MERGERDIR: ${MERGERDIR}"
echo "### SNPEFFPATH= $SNPEFFPATH"
echo "### SNPSIFT= $SNPSIFT"
echo "### SAMTOOLS= $SAMTOOLS"
echo "### VARSCAN= $VARSCAN"
echo "### MATCHEDNORMAL= $MATCHEDNORMAL"
echo "### DBSNP_SNV_BED= $DBSNP_SNV_BED"
echo "### DBSNP_DIV_BED= $DBSNP_DIV_BED"
echo "### REF= $REF"
echo "### DICT= $DICT"
echo "### COSMIC= $COSMIC"
echo "### COSMICC= $COSMICC"
echo "### COSMICNC= $COSMICNC"
echo "### EXAC= ${EXAC}"
echo "### KG= ${KG}"
echo "### NHLBI= ${NHLBI}"
echo "### SNPS= $SNPS"
echo "### INDELS= $INDELS"
echo "### GATK= $GATK"
echo "### VCFMERGER= $VCFMERGER"
echo "### VCFMERGER_DIR= $VCFMERGER_DIR"
echo "### VCFSORTER= $VCFSORTER"
echo "### RNA_VCF_HEADER= $RNA_VCF_HEADER"
echo "### POST_MERGE_VENN= $POST_MERGE_VENN"
echo "### DBSNP= $DBSNP"
echo "### DBVERSION= $DBVERSION"
echo "### SEURAT_VCF= $SEURAT_VCF"
echo "### MUTECT_VCF= $MUTECT_VCF"
echo "### STRELKA_SNV_VCF= $STRELKA_SNV_VCF"
echo "### STRELKA_INDEL_VCF= $STRELKA_INDEL_VCF"
echo "### MERGERDIR= $MERGERDIR"
echo "### RNABAM= $RNABAM"
echo "### ASSAYID= $ASSAYID"
echo "### BEDFILE= $BEDFILE"
echo "### RUNDIR= $RUNDIR"
echo "### CONTROL= $CONTROL"
echo "### TUMOR= $TUMOR"
echo "### REIPE= $RECIPE"
echo "### D= $D"

module load BEDTools/2.14.0 
module load R/3.1.1

SEURAT_BASENAME=`basename ${SEURAT_VCF} ".seurat.vcf"`

#filter the seurat vcf
if [[ ${MATCHEDNORMAL} == "No" ]] ; then
        echo "Sample does not have a matched normal, will filter with bed file ${DBSNP_SNV_bed} for seurat snps"
	cat ${SEURAT_VCF} | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar filter "( TYPE='somatic_SNV' )" | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar intervals -x ${DBSNP_SNV_BED} > ${MERGERDIR}/${SEURAT_BASENAME}_seurat_snv.vcf
else
	echo "matched normal detected running seurat SNV filter"
	cat ${SEURAT_VCF} | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar filter "( TYPE='somatic_SNV' )" > ${MERGERDIR}/${SEURAT_BASENAME}_seurat_snv.vcf
fi
##filter the seurat INDELS, use bed if no matched normal
if [[ ${MATCHEDNORMAL} == "No" ]] ; then
	echo "Sample does not have a matched normal, will filter with bed file ${DBSNP_DIV_bed} for seurat indels"
        cat ${SEURAT_VCF} | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar filter "(( TYPE='somatic_deletion' ) | ( TYPE='somatic_insertion' ))"  | /home/tgenref/pecan/bin/vt/vt normalize - -r ${REF} | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar intervals -x ${DBSNP_DIV_BED} > ${MERGERDIR}/${SEURAT_BASENAME}_seurat_indel.vcf
else
	echo "Matched normal detected running seurat INDEL filter"
        cat ${SEURAT_VCF} | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar filter "(( TYPE='somatic_deletion' ) | ( TYPE='somatic_insertion' ))" > ${MERGERDIR}/${SEURAT_BASENAME}_seurat_indel.vcf
fi

# Filter the INDELs from SEURAT
cat ${MERGERDIR}/${SEURAT_BASENAME}_seurat_snv.vcf | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar filter "(( QUAL >= 15 ) & ( DP1 >= 10 ) & ( DP2 >= 10 ) & ( AR1 <= 0.02 ) & ( AR2 >= 0.05 ))" > ${MERGERDIR}/${SEURAT_BASENAME}.seurat_snv_filt.vcf
cat ${MERGERDIR}/${SEURAT_BASENAME}_seurat_indel.vcf | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar filter "( QUAL >= 25)" > ${MERGERDIR}/${SEURAT_BASENAME}.seurat_indel_filt.vcf
SEURAT_SNV_PATH=${MERGERDIR}/${SEURAT_BASENAME}.seurat_snv_filt.vcf
SEURAT_INDEL_PATH=${MERGERDIR}/${SEURAT_BASENAME}.seurat_indel_filt.vcf

#filter the mutect vcf
MUTECT_BASENAME=`basename ${MUTECT_VCF} "_MuTect_All.vcf"`
#Check the Mutect sample order, leverage our standardized naming to find if its a C or T (STUDY_PATIENT_VISIT_SOURCE_FRACTIONincrement_ASSAYCODE_LIBRARY)
#Need to test the first charcter of the "FRACTIONincrement" is a C or T
echo "The following line contains the Mutect VCF header:"
grep -v "##" ${MUTECT_VCF} | grep "#"
echo "Extracting the first genotype column to determine if it is tumor or normal"
FIRST_GENOTYPE_COLUMN=`grep -v "##" ${MUTECT_VCF} | grep "#" | cut -f10`
echo "The first genotype column header is: ${FIRST_GENOTYPE_COLUMN}"
#FRACTION_LETTER=`echo ${FIRST_GENOTYPE_COLUMN} | cut -d_ -f6 | cut -c1`
#echo "The fraction letter in the first genotype column is: ${FRACTION_LETTER}"
echo "Filtering Mutect calls for ${MUTECT_BASENAME}"
if [ "${FIRST_GENOTYPE_COLUMN}" == "${TUMOR}" ] 
	then
	# Filter the MUTECT calls
	echo "Found expected genotype order - Proceeding with filtering"
	if [[ ${MATCHEDNORMAL} == "No" ]] ; then
		echo "Sample does not have a matched normal, will filter with bed file ${DBSNP_SNV_bed} for mutect snvs"
                cat ${MUTECT_VCF} | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar filter "(( FILTER = 'PASS') & ( GEN[1].FA <= 0.02 ) & ( GEN[0].FA >= 0.05 ))" | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar intervals -x ${DBSNP_SNV_BED} > ${MERGERDIR}/${MUTECT_BASENAME}.mutect_snv_filt.vcf
        else
		echo "Matched normal detected running mutect SNV filter"
		cat ${MUTECT_VCF} | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar filter "(( FILTER = 'PASS') & ( GEN[1].FA <= 0.02 ) & ( GEN[0].FA >= 0.05 ))" > ${MERGERDIR}/${MUTECT_BASENAME}.mutect_snv_filt.vcf
	fi
elif [ "${FIRST_GENOTYPE_COLUMN}" == "${CONTROL}" ] 
	then
	# Reorder the genotype columns and then filter
	##Filtering with a bed file if nor matched normal
	if [[ ${MATCHEDNORMAL} == "No" ]] ; then
	echo "Sample does not have a matched normal, will filter with bed file ${DBSNP_SNV_bed} for mutect snvs"
		awk -F'\t' '{OFS="\t" ; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$10}' ${MUTECT_VCF} | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar filter "(( FILTER = 'PASS') & ( GEN[1].FA <= 0.02 ) & ( GEN[0].FA >= 0.05 ))" | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar intervals -x ${DBSNP_SNV_BED} > ${MERGERDIR}/${MUTECT_BASENAME}.mutect_snv_filt.vcf
	else
		echo "Matched normal detected running mutect SNV filter"
		echo "Found the wrong genotype order - Reordering genotype columns and then Proceeding with filtering"
		awk '{ FS = "\t" ; OFS = "\t" ; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $11, $10}' ${MUTECT_VCF} | java -jar ${SNPSIFT}/SnpSift.jar filter "(( FILTER = 'PASS') & ( GEN[1].FA <= 0.02 ) & ( GEN[0].FA >= 0.05 ))" > ${MERGERDIR}/${MUTECT_BASENAME}.mutect_snv_filt.vcf
	fi
else
	#This should not happen
	echo "ERROR - The mutect vcf did not contain the tumor or the control listed first.  Something is wrong here."
fi
	
sed -i 's/\t*$//g' ${MERGERDIR}/${MUTECT_BASENAME}.mutect_snv_filt.vcf
MUTECT_SNV_VCF=${MERGERDIR}/${MUTECT_BASENAME}.mutect_snv_filt.vcf
#finding STRELKA VCF and copying over, if no matched normal adding in filtering
STRELKA_BASENAME=`basename ${STRELKA_INDEL_VCF} ".strelka.passed.somatic.indels.vcf"`
if [[ ${MATCHEDNORMAL} == "No" ]] ; then
        #STLKA_SNV_VCF=`ls *.passed.somatic.snvs.vcf`
        #STLKA_INDEL_VCF=`ls *.passed.somatic.indels.vcf`
        echo "Sample does not have a matched normal, will filter with bed file ${DBSNP_DIV_bed} and ${DBSNP_SNP_BED}for STRELKA INDELS/SNPS"
        cat ${STRELKA_SNV_VCF} | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar intervals -x ${DBSNP_SNV_BED} > ${MERGERDIR}/${STRELKA_BASENAME}.strelka.passed.somatic.snvs.filt.vcf
        cat ${STRELKA_INDEL_VCF} | /home/tgenref/pecan/bin/vt/vt normalize - -r ${REF} | java -jar -Xmx20g ${SNPSIFT}/SnpSift.jar intervals -x ${DBSNP_DIV_BED} > ${MERGERDIR}/${STRELKA_BASENAME}.strelka.passed.somatic.indels.filt.vcf

        STRELKA_SNV_VCF="${MERGERDIR}/${STRELKA_BASENAME}.strelka.passed.somatic.snvs.filt.vcf"
        STRELKA_INDEL_VCF="${MERGERDIR}/${STRELKA_BASENAME}.strelka.passed.somatic.indels.filt.vcf"
else
	echo "Matched normal detected will not filter for STRELKA SNVs or INDELS"
	#get the strelka vcfs into the MERGERDIR
	cp ${STRELKA_SNV_VCF} ${MERGERDIR}
	cp ${STRELKA_INDEL_VCF} ${MERGERDIR}
fi

echo "Now merging the filtered vcfs from the 3 callers"

######
#Merge vcf from the 3 callers
######

cd ${MERGERDIR}

SEURAT_SNV_PATH_BN=`basename $SEURAT_SNV_PATH`
SEURAT_INDEL_PATH_BN=`basename $SEURAT_INDEL_PATH`
STRELKA_SNV_VCF_BN=`basename $STRELKA_SNV_VCF`
STRELKA_INDEL_VCF_BN=`basename $STRELKA_INDEL_VCF`
MUTECT_SNV_VCF_BN=`basename $MUTECT_SNV_VCF`

#perf stat /home/jkeats/local/merge3vcfs_v20140803/pecan.merge.3vcfs.main.sh --dirscript ${VCFMERGER_DIR} \
perf stat ${VCFMERGER} --dirscript ${VCFMERGER_DIR} \
        --seusnv ${SEURAT_SNV_PATH_BN} \
        --seuindel ${SEURAT_INDEL_PATH_BN} \
        --slksnv ${STRELKA_SNV_VCF_BN} \
        --slkindel ${STRELKA_INDEL_VCF_BN} \
        --mtcsnv ${MUTECT_SNV_VCF_BN} \
	--refgenfa ${REF} \
	--force \
        --outprefix ${SEURAT_BASENAME} 2> ${MERGERDIR}/${SEURAT_BASENAME}.vcfMerger.perfOut
        if [ $? -ne 0 ] ; then
                echo "### vcf merger failed at merge 3 vcfs stage"
                mv ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerFail
                exit
        fi

echo "Finished Merging successfully"

# Sort the Merged VCF for GATK compatibility
echo "Sorting the Merged VCF"
${VCFSORTER} ${DICT} ${MERGERDIR}/${SEURAT_BASENAME}.merge.norm.vcf > ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.vcf
echo "Finished Sorting"
echo "..."

 # Remove unwanted INFO keys
 echo "Removing unwanted INFO keys"
 echo "..."

perf stat java -jar ${SNPSIFT}/SnpSift.jar rmInfo ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.vcf \
	SEURAT_DNA_ALT_ALLELE_FORWARD_FRACTION \
	SEURAT_DNA_ALT_ALLELE_FORWARD \
	SEURAT_DNA_ALT_ALLELE_REVERSE_FRACTION \
	SEURAT_DNA_ALT_ALLELE_REVERSE \
	SEURAT_DNA_ALT_ALLELE_TOTAL_FRACTION \
	SEURAT_DNA_ALT_ALLELE_TOTAL \
	SEURAT_DNA_REF_ALLELE_FORWARD \
	SEURAT_DNA_REF_ALLELE_REVERSE \
	SEURAT_DNA_REF_ALLELE_TOTAL > ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.vcf 2> ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerRmInfoKeys.perfOut
        if [ $? -ne 0 ] ; then
                echo "### vcf merger failed at remove unwanted info keys stage"
                mv ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerFail
                exit
        fi
echo "Finished Removing Unwanted INFO keys"

# Filter to target regions
echo "Filtering calls to respective target regions"
echo "..."

#if [ "${EXOME_TYPE}" == "TSE61" ]
#then
if [ "${ASSAYID}" == "Exome" ] ; then

	perf stat cat ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.vcf | java -jar ${SNPSIFT}/SnpSift.jar intervals ${BEDFILE} > ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.vcf 2> ${MERGERDIR}/${SEURAT_BASENAME}.vcftargetRegions.perfOut              
	if [ $? -ne 0 ] ; then
                echo "### vcf merger failed filtering to target regions"
        	mv ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerFail
                exit
        fi
else	
	echo "This is not an exome, skipping filtering to target regions"
	cp ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.vcf ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.vcf
fi

#elif [ "${EXOME_TYPE}" == "S5U" ]
#then
#	cat ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.vcf | java -jar ${SNPSIFT} intervals ${S5U_BED} > ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.vcf
#fi

echo "Finished Filtering Calls to Targets"

# Annotate the merged vcf using GATK
echo "Annotate the merged VCF"
echo "..."
#Cleanup VCF for testing to remove seurat indels
#grep -v "SEURAT_TYPE=somatic_deletion" ${SEURAT_BASENAME}.merge.sort.clean.f2t.vcf | grep -v "SEURAT_TYPE=somatic_insertion" > ${SEURAT_BASENAME}.merge.sort.clean.f2t2.vcf

if [ "${RECIPE}" == "cerberus" ] ; then

	perf stat java -Xmx24g -jar ${GATK}/GenomeAnalysisTK.jar -R ${REF} \
		-T VariantAnnotator \
		-nt 4 \
		-o ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.vcf \
		-U ALLOW_SEQ_DICT_INCOMPATIBILITY \
		--variant ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.vcf \
		--dbsnp ${DBSNP} \
		--disable_auto_index_creation_and_locking_when_reading_rods 2> ${MERGERDIR}/${SEURAT_BASENAME}.annotateMerge.perfOut 
		if [ $? -ne 0 ] ; then
			echo "### vcf merger failed at annotate merged vcf stage"
			mv ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerFail
			exit
		fi
else
	perf stat java -Xmx24g -jar ${GATK}/GenomeAnalysisTK.jar -R ${REF} \
		-T VariantAnnotator \
		-nt 4 \
		-o ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.vcf \
		-U ALLOW_SEQ_DICT_INCOMPATIBILITY \
		--variant ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.vcf \
		--dbsnp ${DBSNP} \
		--comp:EXAC ${EXAC} \
		--comp:NHLBI ${NHLBI} \
		--comp:1000G ${KG} \
		--comp:COSMIC_NC ${COSMICNC} \
		--comp:COSMIC_C ${COSMICC} \
		--disable_auto_index_creation_and_locking_when_reading_rods 2> ${MERGERDIR}/${SEURAT_BASENAME}.annotateMerge.perfOut 
		if [ $? -ne 0 ] ; then
			echo "### vcf merger failed at annotate merged vcf stage"
			mv ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerFail
			exit
		fi
fi
echo "Finished Annotating Merged VCF"


##### Determine if you need to run RNA allele counts or not
if [ -z "${RNABAM}" ] ; then
	echo "allele count was not requested for this pair"
#	cat ${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.vcf > ${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.rna.vcf
	## add dbNSFP annotaions
	if [ "${RECIPE}" == "cerberus" ] ; then
	
		mv ${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.vcf ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.vcf
	
	else
		perf stat java -jar ${SNPSIFT}/SnpSift.jar dbnsfp \
			-v ${DBNSFP} \
			-a \
			-f Interpro_domain,Polyphen2_HVAR_pred,GERP++_NR,GERP++_RS,LRT_score,MutationTaster_score,MutationAssessor_score,FATHMM_score,Polyphen2_HVAR_score,SIFT_score,Polyphen2_HDIV_score \
			${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.vcf > ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.vcf 2> ${MERGERDIR}/${SEURAT_BASENAME}.dbnsfpannote.perfOut
		if [ $? -ne 0 ] ; then
			echo "### vcf merger failed at annotate with dbNSFP stage"
			mv ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerFail
			exit
		fi
	fi

	# add snpEFF annotations
	java -Xmx4G -jar ${SNPEFFPATH}/snpEff.jar -canon -c ${SNPEFFPATH}/snpEff.config -v -noLog -lof ${DBVERSION} ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.vcf > ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.se74lofcan.vcf
	java -Xmx4G -jar ${SNPEFFPATH}/snpEff.jar -c ${SNPEFFPATH}/snpEff.config -v -noLog -lof ${DBVERSION} ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.vcf > ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.se74lof.vcf
        if [ $? -ne 0 ] ; then
                echo "### vcf merger failed at snpEff annotation stage"
                mv ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerFail
                exit
        fi
	
	# Make final call list venn
	${POST_MERGE_VENN} --vcf ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.se74lofcan.vcf --outprefix  ${MERGERDIR}/${SEURAT_BASENAME}_finalVenn  --maintitle ${SEURAT_BASENAME} --
	if [ $? -ne 0 ] ; then
                echo "### vcf merger failed at venn diagram stage"
                mv ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerFail
                exit
        fi
	
	##clean up final vcfs to save back
	mv ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.se74lofcan.vcf ${MERGERDIR}/${SEURAT_BASENAME}.merged.canonicalOnly.final.vcf
	mv ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.se74lof.vcf ${MERGERDIR}/${SEURAT_BASENAME}.merged.allTranscripts.final.vcf
	rm ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.dbnsfp.vcf
	rm ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.vcf
	rm ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.f2t.ann.vcf.idx
	rm ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.clean.vcf
	rm ${MERGERDIR}/${SEURAT_BASENAME}.merge.sort.vcf
else
	echo "need to run allele count script and finalize for this DNAPAIR"	
	touch ${RUNDIR}/${NXT1}
fi

rm ${MERGERDIR}/${STRELKA_SNV_VCF_BN}
rm ${MERGERDIR}/${STRELKA_INDEL_VCF_BN}

touch ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerPass
rm -rf ${MERGERDIR}/${SEURAT_BASENAME}.vcfMergerInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:VCFMERGER:$hours:$mins" > ${MERGERDIR}/${SEURAT_BASENAME}.vcfMerger.totalTime
time=`date +%d-%m-%Y-%H-%M` 
echo "Ending snpEff Annotator."