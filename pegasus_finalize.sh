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

thisStep="pegasus_nextJob_finalize.txt"
nxtStep1="central_nextJob_saveToIsilon.txt"
nxtStep2="central_nextJob_saveStats.txt"
nxtStep3="central_nextJob_saveToReport.txt"
pbsHome="/home/mrussell/pegasus-pipe/jobScripts"
constants="/home/mrussell/central-pipe/constants/constants.txt"
constantsDir="/home/mrussell/central-pipe/constants"
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


email=`cat $configFile | grep "^EMAIL=" | cut -d= -f2 | head -1 | tr -d [:space:]`

ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

qsubFails=0
missingStuff=0
if [ ! -e $runDir/project.finished ] ; then
	echo "Can't find $runDir/project.finished"
	((missingStuff++))
fi
if [ $missingStuff -ne 0 ] ; then
	echo "Exiting!!!"	
	exit
fi

alreadyDone=0
if [ -e $runDir/finalizePass ] ; then
	echo "Summary stats already passed!"
	touch $runDir/$nxtStep1
	touch $runDir/$nxtStep2
	touch $runDir/$nxtStep3
	alreadyDone=1
fi
if [ -e $runDir/finalizeFail ] ; then
	echo "Summary stats already failed!"
fi
if [ -e $runDir/finalizeInQueue ] ; then
	echo "Summary stats still in queue!"
	alreadyDone=1
	continue
fi
if [ $alreadyDone -eq 0 ] ; then
	touch $runDir/finalizeInQueue
	echo "Create final BAMS and VCFS"
	for sampleLine in `cat $configFile | grep ^SAMPLE=`
	do
        	echo "sample is $sampleLine"
		echo "choosing final bams and creating md5sums for them"
		kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
		samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
		assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
		libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
		insSize=`echo $sampleLine | cut -d= -f2 | cut -d, -f6`
		echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID"
		if [ "$assayID" != "RNA" ] ; then
			echo "### Assay ID is $assayID. Must be genome or exome."
			jirBam="$runDir/$kitName/$samName/${samName}.proj.md.jr.bam"
			mdBam="$runDir/$kitName/$samName/${samName}.proj.md.bam"
			if [ -e $jirBam ] ; then
				echo "JIRBAM found, will make this the final bam $jirBam"
				mv $runDir/$kitName/$samName/${samName}.proj.md.jr.bam $runDir/$kitName/$samName/${samName}.bwa.final.bam
				mv $runDir/$kitName/$samName/${samName}.proj.md.jr.bai $runDir/$kitName/$samName/${samName}.bwa.final.bai
				md5sum $runDir/$kitName/$samName/${samName}.bwa.final.bam > $runDir/$kitName/$samName/${samName}.bwa.final.bam.md5
				
			elif [ -e $mdBam ] ; then
				echo "JIRBAM not found, will make md BAM the final bam $mdBam"
				mv $runDir/$kitName/$samName/${samName}.proj.md.bam $runDir/$kitName/$samName/${samName}.bwa.final.bam
				mv $runDir/$kitName/$samName/${samName}.proj.md.bai $runDir/$kitName/$samName/${samName}.bwa.final.bai
				md5sum $runDir/$kitName/$samName/${samName}.bwa.final.bam $runDir/$kitName/$samName/${samName}.bwa.final.bam.md5
			else
				echo "###ERROR: THERE IS NOT JIRBAM OR MD BAM FOR THIS SAMPLE: $samName"
				touch $runDir/finalizeFail
				exit
			fi
	
		elif [ "$assayID" == "RNA" ] ; then
			echo "### Assay ID is $assayID. Must be RNA sample"
			starFinalBam="$runDir/$kitName/$samName/${samName}.starDir/${samName}.proj.Aligned.out.sorted.md.bam"
			starChimericBam="$runDir/$kitName/$samName/${samName}.starDir/${samName}.proj.Chimeric.out.sorted.bam"
			thFusionBam="$runDir/$kitName/$samName/${samName}.topHatFusionDir/${samName}.proj.accepted_hits.bam"
		
			if [ -e $starFinalBam ] ; then
				echo "STAR final bam found.  will finalize: $starFinalBam"
				mv $runDir/$kitName/$samName/${samName}.starDir/${samName}.proj.Aligned.out.sorted.md.bam $runDir/$kitName/$samName/${samName}.starDir/${samName}.starAligned.final.bam
				mv $runDir/$kitName/$samName/${samName}.starDir/${samName}.proj.Aligned.out.sorted.md.bai $runDir/$kitName/$samName/${samName}.starDir/${samName}.starAligned.final.bai
				md5sum $runDir/$kitName/$samName/${samName}.starDir/${samName}.starAligned.final.bam > $runDir/$kitName/$samName/${samName}.starDir/${samName}.starAligned.final.bam.md5
			else
				echo "###ERROR: THE FINAL STAR BAM IS MISSING: $starFinalBam"
			fi
			if [ -e $starChimericBam ] ; then
                                echo "STAR final chimeric bam found.  will finalize: $starChimericBam"
				mv $runDir/$kitName/$samName/${samName}.starDir/${samName}.proj.Chimeric.out.sorted.bam $runDir/$kitName/$samName/${samName}.starDir/${samName}.starChimeric.final.bam
                        	mv $runDir/$kitName/$samName/${samName}.starDir/${samName}.proj.Chimeric.out.sorted.bai $runDir/$kitName/$samName/${samName}.starDir/${samName}.starChimeric.final.bai
				md5sum $runDir/$kitName/$samName/${samName}.starDir/${samName}.starChimeric.final.bam > $runDir/$kitName/$samName/${samName}.starDir/${samName}.starChimeric.final.bam.md5
			else
                                echo "###ERROR: THE FINAL STAR CHIMERIC BAM IS MISSING: $starChimericBam"
                        fi
                        if [ -e $thFusionBam ] ; then
                                echo "TH FUSION final bam found.  will finalize: $thFusionBam"
				mv $runDir/$kitName/$samName/${samName}.topHatFusionDir/${samName}.proj.accepted_hits.bam $runDir/$kitName/$samName/${samName}.topHatFusionDir/${samName}.thFusion.final.bam
                        	mv $runDir/$kitName/$samName/${samName}.topHatFusionDir/${samName}.proj.accepted_hits.bam.bai $runDir/$kitName/$samName/${samName}.topHatFusionDir/${samName}.thFusion.final.bai
				md5sum $runDir/$kitName/$samName/${samName}.topHatFusionDir/${samName}.thFusion.final.bam > $runDir/$kitName/$samName/${samName}.topHatFusionDir/${samName}.thFusion.final.bam.md5
			else
                                echo "###ERROR: THE FINAL TH FUSION BAM IS MISSING: $thFusionBam"
                        fi

		else
			echo "Assay ID didn't fall into RNA or not RNA. This should never happen"
		fi

	done

	echo "### Create links for BAM files"
	if [ ! -d $runDir/bams ] ; then
		mkdir -p $runDir/bams
	fi
	cd $runDir/bams
	find ../ -name *final.bam -not -path *jointIR* -exec ln -s {} \;	
	#find ../ -name *md.bam -not -path *jointIR* -exec ln -s {} \;	
	#find ../ -name *proj.accepted_hits.bam -not -path *jointIR* -exec ln -s {} \;	
	#find ../ -name *proj.Chimeric.out.sorted.bam -not -path *jointIR* -exec ln -s {} \;	
	#find ../ -name *proj.Aligned.out.sorted.md.bam -not -path *jointIR* -exec ln -s {} \;	

	find ../ -name *final.bai -not -path *jointIR* -exec ln -s {} \;	
	#find ../ -name *md.bai -not -path *jointIR* -exec ln -s {} \;	
	#find ../ -name *proj.accepted_hits.bam.bai -not -path *jointIR* -exec ln -s {} \;	
	#find ../ -name *proj.Chimeric.out.sorted.bai -not -path *jointIR* -exec ln -s {} \;	
	#find ../ -name *proj.Aligned.out.sorted.md.bai -not -path *jointIR* -exec ln -s {} \;	

	find ../ -name *md5 -exec ln -s {} \;

	echo "### Convert pdfs to pngs"
	cd $runDir
	for pdfFile in `find $runDir -name \*pdf`
	do
		pngName="$pdfFile.png"
		convert $pdfFile $pngName
	done
	rm $runDir/finalizeInQueue
	touch $runDir/finalizePass
	touch $runDir/$nxtStep1
	touch $runDir/$nxtStep2
	touch $runDir/$nxtStep3
fi

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
