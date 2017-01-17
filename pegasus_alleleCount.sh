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

thisStep="pegasus_nextJob_alleleCount.txt"
nxtStep1="pegasus_nextJob_postAlleleCount.txt"
pbsHome="/home/tgenjetstream/pegasus-pipe/jobScripts"
constants="/home/tgenjetstream/central-pipe/constants/constants.txt"
constantsDir="/home/tgenjetstream/central-pipe/constants"
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
recipe=`cat $configFile | grep "^RECIPE=" | cut -d= -f2 | head -1 | tr -d [:space:]`
debit=`cat $configFile | grep "^DEBIT=" | cut -d= -f2 | head -1 | tr -d [:space:]`

nCores=`grep @@${myName}_CORES= $constantsDir/$recipe | cut -d= -f2`


ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
alPath=`grep "@@"$recipe"@@" $constants | grep @@ALPATH= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
alleleCount=`grep "@@"$recipe"@@" $constants | grep @@ALLELECOUNT= | cut -d= -f2`
jirRequested=`grep @@"$recipe"@@ $constants | grep @@JIRREALIGN= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###
for acLine in `cat $configFile | grep ^TRIPLET4ALLELECOUNT=`
do
	if [ $alleleCount != "yes" ] ; then
		echo "### Allele count is not requested for this recipe"
		continue
	fi
	echo "### Allele count line is $acLine"
	sampleNames=`echo $acLine | cut -d= -f2`
	usableName=${sampleNames//,/-}
	norml=`echo $sampleNames | cut -d, -f1`
	tumor=`echo $sampleNames | cut -d, -f2`
	rnaIn=`echo $sampleNames | cut -d, -f3`

	normlLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$norml"'"'`
	tumorLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$tumor"'"'`
	rnaInLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$rnaIn"'"'`
	normlKit=`echo $normlLine | cut -d= -f2 | cut -d, -f1`
	tumorKit=`echo $tumorLine | cut -d= -f2 | cut -d, -f1`
	rnaInKit=`echo $rnaInLine | cut -d= -f2 | cut -d, -f1`

	echo "### Normals: $normlKit, $norml"
	echo "### Tumors : $tumorKit, $tumor"
	echo "### RNA in : $rnaInKit, $rnaIn"

	seuratVcf=$runDir/seurat/$norml-$tumor/$norml-$tumor.seurat.snpEff.vcf
	snpEffPass=$runDir/seurat/$norml-$tumor/$norml-$tumor.seurat.vcf.snpEffPass
	if [[ ! -e $seuratVcf || ! -e $snpEffPass ]] ; then
		echo "### Seurat vcf or its snpEffPass does not exist yet"
		((qsubFails++))
		continue
	fi
	
	if [ "${jirRequested}" == "no" ] ; then
                tumorBam=$runDir/$tumorKit/$tumor/$tumor.proj.md.bam
		tumorPass=$runDir/$tumorKit/$tumor/$tumor.proj.bam.mdPass
	else
		tumorBam=$runDir/$tumorKit/$tumor/$tumor.proj.md.jr.bam
		tumorPass=$runDir/$tumorKit/$tumor/$tumor.proj.md.jr.bam.jointIRPass
	fi
	if [[ ! -e $tumorBam || ! -e $tumorPass ]] ; then
		echo "### Tumor bam or its Pass does not exist yet"
		((qsubFails++))
		continue
	fi
	#RNA is from the tumor	

	echo "### Aligner for recipe $recipe is $rnaAligner"
	case $rnaAligner in 
	tophat) echo "tophat case"
		rnaInBam=$runDir/$rnaInKit/$rnaIn/$rnaIn.topHatDir/$rnaIn.proj.accepted_hits.bam
		topHPass=$runDir/$rnaInKit/$rnaIn/$rnaIn.topHatDir.thPass
		if [[ ! -e $rnaInBam || ! -e $topHPass ]] ; then
			echo "### RNA bam or its tophat pass does not exist yet"
			echo "### $rnaInBam"
			echo "### $topHPass"
			((qsubFails++))
			continue
		fi
		;;
	star) echo "star case"
		rnaInBam=$runDir/$rnaInKit/$rnaIn/$rnaIn.starDir/$rnaIn.proj.Aligned.out.sorted.bam
		starPass=$runDir/$rnaInKit/$rnaIn/$rnaIn.starDir.starPass
		if [[ ! -e $rnaInBam || ! -e $starPass ]] ; then
			echo "### RNA bam or its star pass does not exist yet"
			echo "### $rnaInBam"
			echo "### $starPass"
			((qsubFails++))
			continue
		fi
		;;
	anotherRNAaligner) echo "example RNA aligner"
		;;
	*) echo "I should not be here"
		;;
	esac 
	alCountDir="$runDir/alCount"
	trackName=$alCountDir/$usableName
	if [ ! -d $trackName ] ; then
		echo "### Making directory $trackName"
		mkdir -p $trackName
	fi
	outFile=$trackName/$usableName.alleleCount.vcf
	if [[ -e $trackName/$usableName.alleleCountInQueue ||  -e $trackName/$usableName.alleleCountPass || -e $trackName/$usableName.alleleCountFail ]] ; then
		echo "### Alelle count already passed, failed, or inqueue... "
		continue
	fi
	echo "### Submitting to queue for allele count..."
	qsub -v OUT=$outFile,ALCOUNTPATH=$alPath,REF=$ref,TRACK=$trackName/$usableName,VCF=$seuratVcf,RNABAM=$rnaInBam,DNABAM=$tumorBam,RUNDIR=$runDir,D=$d $pbsHome/pegasus_alleleCount.pbs
	if [ $? -eq 0 ] ; then
		touch $trackName/$usableName.alleleCountInQueue
	else
		((qsubFails++))
	fi
	sleep 1

done

if [ $qsubFails -eq 0 ] ; then
#all jobs submitted succesffully, remove this dir from messages
	echo "### I should remove $thisStep from $runDir."
	rm -f $runDir/$thisStep
else
#qsub failed at some point, this runDir must stay in messages
	echo "### Failure in qsub. Not touching $thisStep"
fi

time=`date +%d-%m-%Y-%H-%M`
echo "Ending $0 at $time"
