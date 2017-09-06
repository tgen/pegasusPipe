#!/bin/bash
#####################################################################
# Copyright (c) 2011 by The Translational Genomics Research
# Institute. All rights reserved. This License is limited to, and you may
# use the Software solely for, your own internal and non-commercial use
# for academic and research purposes. Without limiting the foregoing, you
# may not use the Software as part of, or in any way in connection with
# the production, marketing, sale or support of any commercial product or
# service or for any governmental purposes. For commercial or governmental
# use, please contact dcraig@tgen.org. By installing this Software you are
# agreeing to the terms of the LICENSE file distributed with this
# software.
#####################################################################

thisStep="pegasus_nextJob_checkProjectComplete.txt"
nxtStep1="pegasus_nextJob_summaryStats.txt"
pbsHome="~/pegasus-pipe/jobScripts"
constants="~/central-pipe/constants/constants.txt"
constantsDir="~/central-pipe/constants"
myName=`basename $0 | cut -d_ -f2`

time=`date +%d-%m-%Y-%H-%M`
echo "Starting $0 at $time"
if [ "$1" == "" ] ; then
	echo "### Please provide runfolder as the only parameter"
	echo "### Exiting!!!"
	exit
fi
runDir=$1
projName=`basename $runDir | awk -F'_ps20' '{print $1}'`
configFile=$runDir/$projName.config
if [ ! -e $configFile ] ; then
	echo "### Config file not found at $configFile!!!"
	echo "### Exiting!!!"
	exit
else
	echo "### Config file found."
fi
pedFile=$runDir/$projName.ped
if [ ! -e $pedFile ] ; then
        echo "### Ped file not found at $pedFile!!!"
        echo "### Must not be a family study"
	familyStudy="NO"
else
        echo "### Ped file found."
	familyStudy="YES"
fi
recipe=`cat $configFile | grep "^RECIPE=" | cut -d= -f2 | head -1 | tr -d [:space:]`
debit=`cat $configFile | grep "^DEBIT=" | cut -d= -f2 | head -1 | tr -d [:space:]`

nCores=`grep @@${myName}_CORES= $constantsDir/$recipe | cut -d= -f2`
cuffQuant=`grep @@CUFFQUANT= $constantsDir/$recipe | cut -d= -f2`
digar=`grep @@DIGAR= $constantsDir/$recipe | cut -d= -f2`
snpSniff=`grep @@SNPSNIFF= $constantsDir/$recipe | cut -d= -f2`
strelka=`grep @@STRELKA= $constantsDir/$recipe | cut -d= -f2`
mutect=`grep @@MUTECT= $constantsDir/$recipe | cut -d= -f2`
salmon=`grep @@SALMON= $constantsDir/$recipe | cut -d= -f2`
iglBedCov=`grep @@IGLBEDCOV= $constantsDir/$recipe | cut -d= -f2`

ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
fusionDetector=`grep "@@"$recipe"@@" $constants | grep @@FUSIONDETECT= | cut -d= -f2`
cnaGenomeFilter=`grep @@"$recipe"@@ $constants | grep @@CNAGENOMEFILTER= | cut -d= -f2`
cnaGenomeUnFilter=`grep @@"$recipe"@@ $constants | grep @@CNAGENOMEUNFILTER= | cut -d= -f2`
cna=`grep @@"$recipe"@@ $constants | grep @@CNA= | cut -d= -f2`
jiRReq=`grep @@"$recipe"@@ $constants | grep @@JIRREALIGN= | cut -d= -f2`
trnReq=`grep @@"$recipe"@@ $constants | grep @@TRNREQ= | cut -d= -f2`
dexSeqReq=`grep @@"$recipe"@@ $constants | grep @@DEXSEQ= | cut -d= -f2`
species=`grep @@"$recipe"@@ $constants | grep @@SPECIES= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
totalMissing=0
totalChecked=0
missingGermLineSnpEff=0
missingSeuratSnpEff=0
missingMutect=0
missingStrelka=0
missingVcfMerger=0
#missingMpileup=0
#missingFreebayes=0
missingGermVcfMerger=0
missingMultiSampleHCSnpEff=0
missingMultiSampleFBSnpEff=0
missingMultiSampleMPSnpEff=0
missingSamStat=0
missingPicMulti=0
missingPicHS=0
missingPicHS2=0
missingTrn=0
missingRNAalign=0
missingRNAmetric=0
missingFusionDetect=0
missingHtSeq=0
missingCuffLink=0
missingCuffQuant=0
missingSnpSniff=0
missingSalmon=0
missingKallisto=0
missingIglBedCov=0
missingAlleleCount=0
missingvcfMergerAC=0
missingCnaGenFilt=0
missingCnaGenUnfi=0
missingCnaExome=0
missingDeSeq=0
missingDEXseq=0
missingsleuth=0
missingDigar=0
missingCuffDiff=0
missingPhaseBT=0
missingdeNovoGear=0
missingRNAHCSnpEff=0

##This checks for outputs from family studies that had a ped file
if [ "$familyStudy" = "YES" ] ;then
	for dnaPairLine in `cat $configFile | grep '^DNAFAMI='`
	do
		#echo "### DNA pair line is $dnaPairLine"
		sampleNames=`echo $dnaPairLine | cut -d= -f2`
		shorterName=`echo $dnaPairLine | cut -d= -f2 | cut -d\; -f2`
		if [[ $shorterName == *,* ]] ; then
			shorterName=""
			#echo "### No nick name detected for these pairs"
			usableName=${sampleNames//,/-}
		else
			#echo "### Nick name for these pairs are: $shorterName"
			usableName="$shorterName"
		fi
		phaseBTPass="$runDir/phaseByTransmission/$usableName/$usableName.phaseBTPass"
		deNovoGearPass="$runDir/denovoGear/$usableName/$usableName.deNovoGearPass"
		((totalChecked++))
		if [ ! -e $phaseBTPass ] ; then
			echo "### phaseBTPass does not exist"
			echo "### Missing file: $phaseBTPass"
			##Still let projects save back without phaseBT
			#((missingPhaseBT++))
		fi
		#((totalChecked++))
                #if [ ! -e $deNovoGearPass ] ; then
                #        echo "### deNovoGearPass does not exist"
                #        echo "### Missing file: $deNovoGearPass"
                #        ((missingdeNovoGear++))
                #fi 

	done
else
	echo "This isn't a family study, will skip these checks"
fi
## This is a generic loop for each sample 
for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
	kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
	samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
	assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
	libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
	#echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID, insSize: $insSize"
	if [ "$assayID" != "RNA" ] ; then
		#echo "### Assay ID is $assayID. Must be genome or exome."
		pcDir=$runDir/$kitName/$samName
		inBam=$runDir/$kitName/$samName/$samName.proj.bam
		#mdVcf=$runDir/$kitName/$samName/$samName.proj.md.bam.UG.vcf
		#vqVcf=$runDir/$kitName/$samName/$samName.proj.md.bam.UG.vqsr.vcf
		hcVcf=$runDir/$kitName/$samName/$samName.proj.md.bam.HC.vcf
		jrBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
		samStat1=$runDir/stats/$samName.proj.bam.samStats
		samStat2=$runDir/stats/$samName.proj.md.bam.samStats
		picardHSMetric=$runDir/$kitName/$samName/$samName.proj.md.bam.picHSMetricsPass
		picardMulti2=$runDir/$kitName/$samName/$samName.proj.md.bam.picMultiMetricsPass
		((totalChecked++))
		#if [ ! -e $samStat1 ] ; then
		#	echo "### samStats does not exist for proj.bam"
		#	echo "### Missing file: $samStat1"
		#	((missingSamStat++))	
		#fi
		((totalChecked++))
		if [ ! -e $samStat2 ] ; then
			echo "### samStats does not exist for proj.md.bam"
			echo "### Missing file: $samStat2"
			((missingSamStat++))	
		fi
		((totalChecked++))
		if [ ! -e $picardMulti2 ] ; then
			echo "### picardMultiMetricsPass does not exist for proj.md.bam"
			echo "### Missing file: $picardMulti2"
			((missingPicMulti++))
		fi
		if [ "$assayID" == "Exome" ] ; then
			((totalChecked++))
			if [ ! -e $picardHSMetric ] ; then
				echo "### picard HS metric does not exist for proj.bam"
				echo "### Missing file: $picardHSMetric"
				((missingPicHS++))
			fi

		fi
	else
		#echo "### Assay ID is $assayID. Must be RNA."
		#code for calling Multi metrics on tophat or star bams here
		case $rnaAligner in 
		tophat) #echo "### Tophat case"
			rnaAlignPass=$runDir/$kitName/$samName/$samName.topHatDir.thPass
			htSeqPass=$runDir/$kitName/$samName/$samName.topHatDir/$samName.proj.accepted_hits.sam.htSeqPass
			cuffLinkPass=$runDir/$kitName/$samName/$samName.topHatDir.cuffLinkPass
			rnaMtricPass=$runDir/$kitName/$samName/$samName.topHatDir/$samName.proj.accepted_hits.bam.picRNAMetricsPass
			rnaSamStat=$runDir/stats/$samName.proj.accepted_hits.bam.samStats
			((totalChecked++))
			((totalChecked++))
			((totalChecked++))
			((totalChecked++))
			((totalChecked++))
			;;
		star) #echo "### Star case"
			rnaAlignPass=$runDir/$kitName/$samName/$samName.starDir.starPass
			htSeqPass=$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sam.htSeqPass
			cuffLinkPass=$runDir/$kitName/$samName/$samName.starDir.cuffLinkPass
			rnaMtricPass=$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sorted.md.bam.picRNAMetricsPass
			rnaSamStat=$runDir/stats/$samName.proj.Aligned.out.sorted.md.bam.samStats
			((totalChecked++))
			((totalChecked++))
			((totalChecked++))
			((totalChecked++))
			((totalChecked++))
			;;
		anotherRNAaligner) echo "#### A new rna aligner case"
			rnaAlignPass=$runDir/$kitName/$samName/$samName.newAlignerDir.newAlignerPass
			((totalChecked++))
			;;
		*) echo "### I should not be here."
			;;
		esac
		if [ ! -e $rnaAlignPass ] ; then
			echo "### Missing rnaAlign pass file for $rnaAligner"
			echo "### Missing file: $rnaAlignPass"
			((missingRNAalign++))
		fi
		if [ ! -e $htSeqPass ] ; then
			echo "### Missing htSeq pass file for $rnaAligner"
			echo "### Missing file: $htSeqPass"
			((missingHtSeq++))
		fi
		if [ ! -e $cuffLinkPass ] ; then
			echo "### Missing cuffLinkPass file for $rnaAligner"
			echo "### Missing file: $cuffLinkPass"
			((missingCuffLink++))
		fi
		if [ ! -e $rnaMtricPass ] ; then
			echo "### Missing rnaMetric pass file for $rnaAligner"
			echo "### Missing file: $rnaMtricPass"
			((missingRNAmetric++))
		fi
		if [ ! -e $rnaSamStat ] ; then
			echo "### samStats does not exist for rna bam"
			echo "### Missing file: $rnaSamStat"
			((missingSamStat++))	
		fi
		case $fusionDetector in 
		tophatfusion) #echo "### Tophat case"
			detectFusionPass=$runDir/$kitName/$samName/$samName.topHatFusionDir.thFPostPass
			((totalChecked++))
			;;
		soapFuse) #echo "### Star case"
			detectFusionPass=$runDir/$kitName/$samName/$samName.soapFuseDir.soapFusePass
			((totalChecked++))
			;;
		none) #echo "### Star case"
			detectFusionPass=$configFile
			((totalChecked++))
			;;
		anotherFusionDetector) echo "#### A new fusion detector case"
			detectFusionPass=$runDir/$kitName/$samName/$samName.newAlignerDir.newAlignerPass
			((totalChecked++))
			;;
		*) echo "### I should not be here."
			;;
		esac
		((totalChecked++))
		if [ ! -e $detectFusionPass ] ; then
			echo "### Missing fusion pass file for $samName"
			echo "### Missing file: $detectFusionPass"
			((missingFusionDetect++))
		fi

	fi
done
## This is a generic loop for each sample for non default options such as cuff quant
for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
	kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
	samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
	assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
	libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`

	mdBam=$runDir/$kitName/$samName/$samName.proj.md.bam
	jrBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
	jrPas=$runDir/$kitName/$samName/$samName.proj.md.jr.bam.jointIRPass
	mdPas=$runDir/$kitName/$samName/$samName.proj.bam.mdPass

	#echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID, insSize: $insSize"
	if [ "$assayID" != "RNA" ] ; then
		#echo "### Assay ID is $assayID. Must be genome or exome."
		pcDir=$runDir/$kitName/$samName
		jrRequested=`cat $configFile | grep '^DNAPAIR=\|^DNAFAMI=\|^JIRSET=' | grep $samName | wc -l`
		if [ $jrRequested -gt 0 ] ; then
			pasFile=$jrPass
			bamFile=$jrBam
		else
			pasFile=$mdPas
			bamFile=$mdBam
		fi
		if [ "$snpSniff" == "yes" ] ; then
			((totalChecked++))
			if [  ! -e $bamFile.snpSniffPass ] ; then
				echo "### Missing snpSniff file for $bamFile"	
				echo "### Missing file: $bamFile.snpSniffPass"
				((missingSnpSniff++))
			fi
		fi

		#if [ "$digar" == "yes" ] ; then
		#	((totalChecked++))
		#	if [  ! -e $bamFile.snpSniffPass ] ; then
		#	echo "### Missing digar post file for $bamFile"   
		#	echo "### Missing file: $bamFile.snpSniffPass"
		#	((missingSnpSniff++))
		#	fi
                #fi
	else
		#echo "### Assay ID is $assayID. Must be RNA."
		#code for calling Multi metrics on tophat or star bams here
		case $rnaAligner in 
		tophat) #echo "### Tophat case"
			cuffQuantPass=$runDir/$kitName/$samName/$samName.topHatDir.cuffQuantPass
			snpSniffPass=$runDir/$kitName/$samName/$samName.topHatDir/$samName.proj.accpeted_hits.md.bam.snpSniffPass
			salmonPass=$runDir/$kitName/$samName/$samName.topHatDir.salmonPass
			iglBedCovPass=$runDir/$kitName/$samName/$samName.topHatDir/$samName.proj.accpeted_hits.md.bam.IGLbedCovPass
			;;
		star) #echo "### Star case"
			cuffQuantPass=$runDir/$kitName/$samName/$samName.starDir.cuffQuantPass
			snpSniffPass=$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sorted.md.bam.snpSniffPass
			salmonPass=$runDir/$kitName/$samName/$samName.salmonDir.salmonPass
			iglBedCovPass=$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sorted.md.bam.IGLbedCovPass
			kallistoPass=$runDir/kallisto/$samName.kallistoPass
			;;
		anotherRNAaligner) echo "#### A new rna aligner case"
			rnaAlignPass=$runDir/$kitName/$samName/$samName.newAlignerDir.newAlignerPass
			;;
		*) echo "### I should not be here."
			;;
		esac
		if [ "$cuffQuant" == "yes" ] ; then
			((totalChecked++))
			if [ ! -e $cuffQuantPass ] ; then
				echo "### Missing cuffQuantPass file for $rnaAligner"
				echo "### Missing file: $cuffQuantPass"
				((missingCuffQuant++))
			fi
		fi
		if [ "$snpSniff" == "yes" ] ; then
			((totalChecked++))
			if [ ! -e $snpSniffPass ] ; then
				echo "### Missing snpSniffPass file for $rnaAligner"
				echo "### Missing file: $snpSniffPass"
				((missingSnpSniff++))
			fi
		fi
		if [ "$salmon" == "yes" ] ; then
			((totalChecked++))
			if [ ! -e $salmonPass ] ; then
				echo "### Missing salmonPass file for $rnaAligner"
				echo "### Missing file: $salmonPass"
				((missingSalmon++))
			fi
		fi
		if [ "$kallisto" == "yes" ] ; then
                        ((totalChecked++))
                        if [ ! -e $kallistoPass ] ; then
                                echo "### Missing kallistoPass file for $rnaAligner"
                                echo "### Missing file: $kallistoPass"
                                ((missingKallisto++))
                        fi
                fi

		if [ "$iglBedCov" == "yes" ] ; then
			((totalChecked++))
			if [ ! -e $iglBedCovPass ] ; then
				echo "### Missing igl bed cov file for $rnaAligner"
				echo "### Missing file: $iglBedCovPass"
				((missingIglBedCov++))
			fi
		fi
	fi
done

####this loop is for finding allele count output
for acLine in `cat $configFile | grep '^TRIPLET4ALLELECOUNT='`
do
	echo "### Checking allele count output $acLine"
        sampleNames=`echo $acLine | cut -d= -f2`
        usableName=${sampleNames//,/-}
	trackName="$runDir/alCount/$usableName/$usableName"
	((totalChecked++))
	if [ ! -e $trackName.alleleCountPass ] ; then 
		echo "### alleleCountPass does not exist on allele output!"
		echo "### Missing file: $trackName.alleleCountPass"
		((missingAlleleCount++))
	fi 
	usableName2=`echo $usableName | cut -d '-' -f1-2`
	vcfACmergeTrack=$runDir/vcfMerger/$usableName2/$usableName2
        if [ ! -e $vcfACmergeTrack.vcfMergerACPass ] ; then
                echo "### vcfMerger alleleCountPass does not exist on allele output!"
                echo "### Missing file: $vcfACmergeTrack.vcfMergerACPass"
                ((missingvcfMergerAC++))
        fi
done

###this for loop is for finding translocations 
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
        if [ $trnReq == "no" ] ; then
                echo "### TRN is not requested for this recipe"
                continue
        fi
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
        pair1=`echo $sampleNames | cut -d, -f1`
        pair2=`echo $sampleNames | cut -d, -f2`
	outName="$pair1-$pair2"
        pair1SampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$pair1"'"'`
        pair2SampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$pair2"'"'`
        pair1KitName=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f1`
        pair2KitName=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f1`
        pair1AssayNo=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f3`
        pair2AssayNo=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f3`
	case $pair1AssayNo$pair2AssayNo in 
	GenomeGenome) echo "### Both normal and tumor are Genome"
		echo "### Checking translocations output $dnaPairLine"
		sampleNames=`echo $dnaPairLine | cut -d= -f2`
		usableName=${sampleNames//,/-}
		trackName="$runDir/trn/$usableName/$usableName"
		((totalChecked++))
		if [ ! -e $trackName.trnPass ] ; then 
			echo "### trnPass does not exist!!"
			echo "### Missing file: $trackName.trnPass"
			((missingTrn++))
		fi 
	;;
	ExomeExome) echo "### Both normal and tumor are Exome"
		echo " Checking translocations on output $dnaPairLine"
		sampleNames=`echo $dnaPairLine | cut -d= -f2`
		usableName=${sampleNames//,/-}
		trackName="$runDir/trn/$usableName/$usableName"
		((totalChecked++))
		if [ ! -e $trackName.trnPass ] ; then 
			echo "### trnPass does not exist!!"
			echo "### Missing file: $trackName.trnPass"
			((missingTrn++))
		fi 
	;;
	*) echo "### Something else"
	;;
	esac
done
###this for loop is for finding cna outputs
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
	if [ $cna != "yes" ] ; then
		echo "### CNA is not requested for this recipe"
		continue
	fi

	echo "### Checking CNA output for $dnaPairLine"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
        pair1=`echo $sampleNames | cut -d, -f1`
        pair2=`echo $sampleNames | cut -d, -f2`
	outName="$pair1-$pair2"
        pair1SampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$pair1"'"'`
        pair2SampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$pair2"'"'`
        pair1KitName=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f1`
        pair2KitName=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f1`
        pair1AssayNo=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f3`
        pair2AssayNo=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f3`
	case $pair1AssayNo$pair2AssayNo in 
	GenomeGenome) echo "### Both normal and tumor are Genome"
		if [ "$cnaGenomeFilter" = "yes" ] ; then
			((totalChecked++))
			outName1=$outName"_filt"
			trackName=$runDir/cna/$outName1/$outName1
			if [ ! -e $trackName.cnaGenFiltPass ] ; then
				echo "### cnaGenFiltPass does not exist on cna output!"
				echo "### Missing file: $trackName.cnaGenFiltPass" 
				((missingCnaGenFilt++))
			fi
		fi
		if [ "$cnaGenomeUnFilter" = "yes" ] ; then
			((totalChecked++))
			outName2=$outName"_unfi"
			trackName=$runDir/cna/$outName2/$outName2
			if [ ! -e $trackName.cnaGenUnfiPass ] ; then
				echo "### cnaGenUnfiPass does not exist on cna output!"
				echo "### Missing file: $trackName.cnaGenUnfiPass" 
				((missingCnaGenUnfi++))
			fi
		fi
	;;
	ExomeExome) echo "### Both normal and tumor are Exome"
		((totalChecked++))
		outName3=$outName"_exo"
		trackName=$runDir/cna/$outName3/$outName3
		if [ ! -e $trackName.cnaExomePass ] ; then
			echo "### cnaExomePass does not exist on cna output!"
			echo "### Missing file: $trackName.cnaExomePass" 
			((missingCnaExome++))
		fi
	;;
	*) echo "### Something else"
	;;
	esac
done

###this loop is for finding seurat vcfs
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
	echo "### Checking Seurat output for $dnaPairLine"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
	usableName=${sampleNames//,/-}
	trackName="$runDir/seurat/$usableName/$usableName"
	((totalChecked++))
	if [ ! -e $trackName.seurat.vcf.snpEffPass ] ; then 
		echo "### snpEffPass does not exist on Seurat output!"
		echo "### Missing file: $trackName.seurat.vcf.snpEffPass"
		((missingSeuratSnpEff++))
	fi 
done
###this loop is for finding mutect vcfs and somatic vcf merger vcfs
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
	if [ "$mutect" == "yes" ] ; then
		echo "### Checking mutect output for $dnaPairLine"
		sampleNames=`echo $dnaPairLine | cut -d= -f2`
		usableName=${sampleNames//,/-}
		trackName="$runDir/mutect/$usableName/$usableName"
		((totalChecked++))
		if [ ! -e ${trackName}_MuTect_All.vcf.snpEffPass ] ; then 
			echo "### snpEff pass does not exist for final mutect output "
			echo "### Missing file: ${trackName}_MuTect_All.vcf.snpEffPass"
			((missingMutect++))
		fi 
		echo "Checking vcf merger output for $dnaPairLine"
		vcfMergerTrack="$runDir/vcfMerger/$usableName/$usableName"
		if [ ! -e ${vcfMergerTrack}.vcfMergerPass ] ; then
			echo "### vcfMerger pass does not exist for the 3 somatic callers"
			echo "### missing file: ${vcfMergerTrack}.vcfMergerPass"
			((missingVcfMerger++))
		fi
	fi
done
###this loop is for finding strelka vcfs
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
	if [ "$strelka" == "yes" ] ; then
		echo "### Checking strelka output for $dnaPairLine"
		sampleNames=`echo $dnaPairLine | cut -d= -f2`
		usableName=${sampleNames//,/-}
		trackName="$runDir/strelka/$usableName"
		trackName="$runDir/strelka/$usableName"
		vcfPre="$runDir/strelka/$usableName/myAnalysis/results/$usableName"
		((totalChecked++))
		((totalChecked++))
		((totalChecked++))
		((totalChecked++))
		#if [ ! -e $trackName.strelkaPass ] ; then 
		#	echo "### strelkaPass does not exist "
		#	echo "### Missing file: $trackName.strelkaPass"
		#	((missingStrelka++))
		#fi 
		if [ ! -e $vcfPre.strelka.passed.somatic.snvs.vcf.snpEffPass ] ; then 
			echo "### snpEff is already done, failed or inQueue for passed somatic snvs"
			echo "### Missing file: $vcfPre.strelka.passed.somatic.snvs.vcf.snpEffPass"
			((missingStrelka++))
		fi 
		if [ ! -e $vcfPre.strelka.all.somatic.snvs.vcf.snpEffPass ] ; then 
			echo "### snpEff is already done, failed or inQueue for all somatic snvs"
			echo "### Missing file: $vcfPre.strelka.all.somatic.snvs.vcf.snpEffPass"
			((missingStrelka++))
		fi 
		if [ ! -e $vcfPre.strelka.passed.somatic.indels.vcf.snpEffPass ] ; then 
			echo "### snpEff is already done, failed or inQueue for passed somatic indels"
			echo "### Missing file: $vcfPre.strelka.passed.somatic.indels.vcf.snpEffPass"
			((missingStrelka++))
		fi 
		if [ ! -e $vcfPre.strelka.all.somatic.indels.vcf.snpEffPass ] ; then 
			echo "### snpEff is already done, failed or inQueue for all somatic indels"
			echo "### Missing file: $vcfPre.strelka.all.somatic.indels.vcf.snpEffPass"
			((missingStrelka++))
		fi 
	fi
done
###this loop is for finding picard HS metrics for jr.bam's
for dnaPairLine in `cat $configFile | grep '^DNAPAIR=\|^DNAFAMI='`
do
	echo "### Checking for HS metrics output for samples in $dnaPairLine"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
	for eachSample in ${sampleNames//,/ }
	do
                sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
                kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
                samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
		assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
		if [ "$assayID" == "Exome" ] ; then
			((totalChecked++))
			eachSamplePass=$runDir/$kitName/$samName/$samName.proj.md.jr.bam.picHSMetricsPass
			if [ ! -e $eachSamplePass ] && [ $jiRReq != "no" ] ; then 
				echo "### HS Metrics on jr bam does not exist "
				echo "### Missing file: $eachSamplePass"
				((missingPicHS2++))
			fi 

		fi
	done	
done

## This loop is for finding multi variant called vcf files from GATK
for dnaPairLine in `cat $configFile | grep '^DNAPAIR=\|^DNAFAMI='`
do
	echo "### Checking multi variant calling for $dnaPairLine"
	#sampleNames=`echo $dnaPairLine | cut -d= -f2`
	sampleNames=`echo $dnaPairLine | cut -d= -f2 | cut -d\; -f1`
	shorterName=`echo $dnaPairLine | cut -d= -f2 | cut -d\; -f2`
	if [[ $shorterName == *,* ]] ; then
		shorterName=""
		echo "### No nick name detected for these pairs"
		usableName=${sampleNames//,/-}
	else
		#echo "### Nick name for these pairs are: $shorterName"
		usableName="$shorterName"
	fi
	#usableName=${sampleNames//,/-}
	#trackName="$runDir/ug/$usableName/$usableName"
	trackName="$runDir/hc/$usableName/$usableName"
	((totalChecked++))
	if [ ! -e $trackName.HC_All.vcf.snpEffPass ] ; then 
		echo "### snpEffPass does not exist on multi sample HC output!"
		echo "### Missing file: $trackName.HC_All.vcf.snpEffPass"
		((missingMultiSampleHCSnpEff++))
	fi
	fbtrackName="$runDir/freebayes/$usableName/$usableName"
        ((totalChecked++))
        if [ ! -e $fbtrackName.freebayes_All.vcf.snpEffPass ] ; then
                echo "### snpEffPass does not exist on multi sample freebayes output!"
                echo "### Missing file: $fbtrackName.freebayes_All.vcf.snpEffPass"
                ((missingMultiSampleFBSnpEff++))
        fi
	mptrackName="$runDir/mpileup/$usableName/$usableName"
        ((totalChecked++))
        if [ ! -e $mptrackName.mpileup_All.vcf.snpEffPass ] ; then
                echo "### snpEffPass does not exist on multi sample mpileup output!"
                echo "### Missing file: $mptrackName.mpileup_All.vcf.snpEffPass"
                ((missingMultiSampleMPSnpEff++))
        fi
	((totalChecked++))
	echo "### Checking germline vcf merger output for $dnaPairLine"
	germVcfMergerTrack="$runDir/germlineVcfMerger/$usableName/$usableName"
	if [ ! -e ${germVcfMergerTrack}.vcfMergerPass ] ; then
		echo "### germline vcfMerger pass does not exist for the 3 germline callers"
		echo "### missing file: ${germVcfMergerTrack}.vcfMergerPass"
		((missingGermVcfMerger++))
	fi 
done
## This loop for finding deSeq pairs
for rnaPairLine in `cat $configFile | grep '^RNAPAIR'`
do
	((totalChecked++))
	echo "### Checking deSeq and cuffDiff for $rnaPairLine"
	sampleNames=`echo $rnaPairLine | cut -d= -f2`
	normls=`echo $sampleNames | cut -d, -f1`
	tumors=`echo $sampleNames | cut -d, -f2`
	cdName=`echo $sampleNames | cut -d, -f3`
	dsName=`echo $sampleNames | cut -d, -f3`
	if [ "$cdName" == "" ] ; then
		normlsWithDash=${normls//;/-}
		tumorsWithDash=${tumors//;/-}
		cdName="$normlsWithDash-VS-$tumorsWithDash"
	fi 
	if [ "$dsName" == "" ] ; then
		normlsWithDash=${normls//;/-}
		tumorsWithDash=${tumors//;/-}
		dsName="$normlsWithDash-VS-$tumorsWithDash"
	fi 
	normlList2=""
	tumorList2=""
	for eachNorml in ${normls//;/ }
	do
		normlList2="$eachNorml-$normlList2"		
	done
	for eachTumor in ${tumors//;/ }
	do
		tumorList2="$eachTumor-$tumorList2"		
	done
	normlList2="${normlList2%?}"
	tumorList2="${tumorList2%?}"
	deSeqDir=$runDir/deSeq/$dsName
	if [[ ! -e $deSeqDir.deSeqPass ]] ; then 
		echo "### deSeqPass does not exist for RNAPAIR"
		echo "### Missing file: $deSeqDir.deSeqPass"
		((missingDeSeq++))
	fi
        DEXseqDir=$runDir/DEXseq/$dsName
	if [[ "$sampleNames" != *";"* ]] ; then
                echo "There is only 2 RNA samples, at least 3 are needed to run DEXseq at this time"
                echo "will skip DEXseq check for $sampleNames"
        else
                echo "There are at least 3 RNA samples, will check for DEXseq pass"
        	if [[ ! -e $DEXseqDir.DEXseqPass ]] ; then
			if [ $dexSeqReq == "no" ] ; then
				echo "### DEXSEQ is not requested for this recipe"
			else
				echo "### DEXseqPass does not exist for RNAPAIR"
				echo "### Missing file: $DEXseqDir.DEXseqPass"
				((missingDEXseq++))
			fi
        	fi

	
	fi
	sleuthDir=$runDir/sleuth/$dsName
        #if [[ ! -e $sleuthDir.sleuthPass ]] ; then
        #        echo "### sleuthPass does not exist for RNAPAIR"
        #        echo "### Missing file: $sleuthDir.sleuthPass"
        #        ((missingsleuth++))
        #fi
	cuffDiffDir=$runDir/cuffDiff/$cdName
	if [[ ! -e $cuffDiffDir.cuffDiffPass ]] ; then
		echo "### cuffDiffPass does not exist for RNAPAIR"
		echo "### Missing file: $cuffDiffDir.cuffDiffPass"
		((missingCuffDiff++))
	fi 
done
## This loop is for finding single bam called vcf files from GATK, freebayes, and samtools mpileup
for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
	kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
	samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
	assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
	libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
	if [ "$assayID" != "RNA" ] ; then
		#echo "### Assay ID is $assayID. Must be genome or exome."
		pcDir=$runDir/$kitName/$samName
		inBam=$runDir/$kitName/$samName/$samName.proj.bam
		hcVcf=$runDir/$kitName/$samName/$samName.proj.md.bam.HC_All.vcf
		fbVcf=$runDir/freebayes/$samName/$samName.freebayes_All.vcf
		mpVcf=$runDir/mpileup/$samName/$samName.mpileup_All.vcf
		jrBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
		jrRequested=`cat $configFile | grep '^DNAPAIR=\|^DNAFAMI\|^JIRSET=' | grep $samName | wc -l`
		if [ $jrRequested -gt 0 ] ; then
			echo "### File check will not operate on a single vcf because Joint IR is requested for $samName"
		else
			#echo "### File check will run on single vcf because Joint IR is NOT requested for $samName"
			((totalChecked++))
			if [ ! -e $hcVcf.snpEffPass ] ; then
				echo "### snpEffPass does not exist for germline!"
				echo "### Missing file: $hcVcf.snpEffPass"
				((missingGermLineSnpEff++))
			fi
			((totalChecked++))
                        if [ ! -e $fbVcf.snpEffPass ] ; then
                                echo "### snpEffPass does not exist for germline!"
                                echo "### Missing file: $fbVcf.snpEffPass"
                                ((missingGermLineSnpEff++))
                        fi
                        ((totalChecked++))
                        if [ ! -e $mpVcf.snpEffPass ] ; then
                                echo "### snpEffPass does not exist for germline!"
                                echo "### Missing file: $mpVcf.snpEffPass"
                                ((missingGermLineSnpEff++))
                        fi
			echo "### Checking germline vcf merger output for $dnaPairLine"

			germVcfMergerTrack="$runDir/germlineVcfMerger/$samName/$samName"
			if [ ! -e ${germVcfMergerTrack}.vcfMergerPass ] ; then
				echo "### germline vcfMerger pass does not exist for the 3 germline callers"
				echo "### missing file: ${germVcfMergerTrack}.vcfMergerPass"
				((missingGermVcfMerger++))
			fi

		fi
	##This checks for single sample merged germline vcf

	else
		echo "### Assay ID is $assayID. Must be RNA." > /dev/null
                hcPas=$runDir/rnaHC/$samName/$samName.rnaHC_All.vcf.snpEffPass
                jrRequested=`cat $configFile | grep '^DNAPAIR=\|^DNAFAMI=' | grep $samName | wc -l`
                if [ $jrRequested -gt 0 ] ; then
                        echo "### File check will not operate on a single bam because Joint IR is requested for $samName"
                else
			((totalChecked++))
                        if [ ! -e $hcPas ] ; then
				echo "### snpEffPass does not exist for germline!"
                                echo "### Missing file: $hcPas"
                                ((missingRNAHCSnpEff++))
                        fi
                fi
	fi
done

totalMissing=$(( missingKallisto + missingsleuth + missingDEXseq + missingCnaExome + missingRNAHCSnpEff + missingPhaseBT + missingGermVcfMerger + missingdeNovoGear + missingCnaGenUnfi + missingCnaGenFilt + missingHtSeq + missingRNAalign + missingRNAmetric + missingGermLineSnpEff + missingSeuratSnpEff + missingMultiSampleHCSnpEff + missingMultiSampleFBSnpEff + missingMultiSampleMPSnpEff + missingSamStat + missingPicMulti + missingPicHS + missingPicHS2 + missingTrn + missingAlleleCount + missingvcfMergerAC + missingDeSeq + missingCuffDiff + missingCuffLink + missingCuffQuant + missingSnpSniff + missingSalmon + missingIglBedCov + missingFusionDetect + missingMutect + missingStrelka + missingVcfMerger ))
if [ $totalMissing -eq 0 ] ; then
	echo "### Checked for $totalChecked things and found them all."
	touch $runDir/project.finished	
	rm $runDir/project.incomplete
	echo "### Marking summary stats to run..."
	touch $runDir/$nxtStep1
	rm -f $runDir/$thisStep
else
	echo "### There were $totalMissing missing things out of $totalChecked looked for."
	echo "### There were $totalMissing missing things out of $totalChecked looked for." > $runDir/project.incomplete
fi

time=`date +%d-%m-%Y-%H-%M`
echo "Ending $0 at $time"
