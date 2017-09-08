#!/usr/bin/env bash
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

thisStep="pegasus_nextJob_cuffDiff.txt"
nxtStep1="pegasus_nextJob_postCuffDiff.txt"

constants=~/jetstream/constants/constants.txt
constantsDir=~/jetstream/constants/
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


usegtf=`cat $configFile | grep "^CUFFLINKUSEGTF=" | cut -d= -f2 | head -1 | tr -d [:space:]`
usemask=`cat $configFile | grep "^CUFFLINKUSEMASK=" | cut -d= -f2 | head -1 | tr -d [:space:]`

ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
cuffdiffPath=`grep @@"$recipe"@@ $constants | grep @@CUFFDIFFPATH= | cut -d= -f2`
gtfmask=`grep @@"$recipe"@@ $constants | grep @@GTFMASK= | cut -d= -f2`
gtf=`grep @@"$recipe"@@ $constants | grep @@GTF= | cut -d= -f2`
cuffdiff2vcfPath=`grep @@"$recipe"@@ $constants | grep @@CUFFDIFF2VCFPATH= | cut -d= -f2`
processcdlistPath=`grep @@"$recipe"@@ $constants | grep @@PROCESSCDLISTPATH= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###
for rnaPairLine in `cat $configFile | grep ^RNAPAIR=`
do
	echo "### RNA pair line is $rnaPairLine"
	sampleNames=`echo $rnaPairLine | cut -d= -f2`
	normls=`echo $sampleNames | cut -d, -f1`
	tumors=`echo $sampleNames | cut -d, -f2`
	cdName=`echo $sampleNames | cut -d, -f3`
	if [ "$cdName" == "" ] ; then
		normlsWithDash=${normls//;/-}
		tumorsWithDash=${tumors//;/-}
		cdName="$normlsWithDash-VS-$tumorsWithDash"
	fi 
	echo "### Normals: $normls"
	echo "### Tumors : $tumors"
	tumorCount=0
	normlCount=0
	missingTumorCount=0
	missingNormlCount=0
	normlList=""
	tumorList=""
	case $rnaAligner in 
	tophat) echo "tophat case"
		for eachNorml in ${normls//;/ }
		do
			((normlCount++))
			echo "eachnorml: $eachNorml"
			sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachNorml"'"'`
			kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
			samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
			eachNormlBam=$runDir/$kitName/$samName/$samName.topHatDir/$samName.proj.accepted_hits.bam
			eachNormlPas=$runDir/$kitName/$samName/$samName.topHatDir.thPass
			if [[ ! -e $eachNormlBam || ! -e $eachNormlPas ]] ; then
				echo "### Can't find accp hits bam or tophat pass"
				((missingNormlCount++))
			else
				normlList="$eachNormlBam,$normlList"		
			fi
		done
		for eachTumor in ${tumors//;/ }
		do
			((tumorCount++))
			echo "eachtumor: $eachTumor"
			sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachTumor"'"'`
			kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
			samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
			eachTumorBam=$runDir/$kitName/$samName/$samName.topHatDir/$samName.proj.accepted_hits.bam
			eachTumorPas=$runDir/$kitName/$samName/$samName.topHatDir.thPass
			if [[ ! -e $eachTumorBam || ! -e $eachTumorPas ]] ; then
				echo "### Can't find accp hits bam or tophat pass"
				((missingTumorCount++))
			else
				tumorList="$eachTumorBam,$tumorList"		
			fi
		done
		if [[ $missingNormlCount -eq 0 && $missingTumorCount -eq 0 ]] ; then
			echo "### All of normal and tumor bams are found"
			echo "### Norml list: $normlList"
			echo "### Tumor list: $tumorList"
		else
			echo "### There were missing things"
			echo "### Normls missing: $missingNormlCount/$normlCount"
			echo "### Tumors missing: $missingTumorCount/$tumorCount"
			((qsubFails++))
			continue
		fi
		cdDir=$runDir/cuffDiff/$cdName
		if [ ! -d $cdDir ] ; then
			mkdir -p $cdDir
		fi
		if [[ -e $cdDir.cuffDiffPass || -e $cdDir.cuffDiffFail || -e $cdDir.cuffDiffInQueue ]] ; then 
			echo "### Cuff diff is already done, failed or inQueue"
			continue
		fi 

		echo "### Submitting $cdName to queue for cuff diff..."
		sbatch -n 1 -N 1 --cpus-per-task $nCores --export CUFFDIFF2VCFPATH=$cuffdiff2vcfPath,PROCESSCDLISTPATH=$processcdlistPath,CUFFDIFFPATH=$cuffdiffPath,RUNDIR=$runDir,DIRNAME=$cdDir,GTF=$gtf,BAM1="'"$normlList"'",BAM2="'"$tumorList"'",REF=$ref,MASK=$gtfmask,REF=$ref,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_cuffDiff.pbs
		if [ $? -eq 0 ] ; then
			touch $cdDir.cuffDiffInQueue
		else
			((qsubFails++))
		fi
		sleep 2

		;;
	star) echo "star case"
		for eachNorml in ${normls//;/ }
		do
			((normlCount++))
			echo "eachnorml: $eachNorml"
			sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachNorml"'"'`
			kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
			samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
			rnaNormStrand=`grep @@"$kitName" $constants | cut -d= -f2`
			eachNormlBam=$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sorted.md.bam
			eachNormlPas=$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sorted.bam.rnaMarkDupPass
			if [[ ! -e $eachNormlBam || ! -e $eachNormlPas ]] ; then
				echo "### Can't find aligned out sorted bam or star pass"
				((missingNormlCount++))
			else
				normlList="$eachNormlBam,$normlList"		
			fi
		done
		for eachTumor in ${tumors//;/ }
		do
			((tumorCount++))
			echo "eachtumor: $eachTumor"
			sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachTumor"'"'`
			kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
			samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
			rnaTumorStrand=`grep "@@"$kitName"@@" $constants | cut -d= -f2`
			eachTumorBam=$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sorted.md.bam
			eachTumorPas=$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sorted.bam.rnaMarkDupPass
			if [[ ! -e $eachTumorBam || ! -e $eachTumorPas ]] ; then
				echo "### Can't find aligned out sorted bam or star pass"
				((missingTumorCount++))
			else
				tumorList="$eachTumorBam,$tumorList"		
			fi
		done
		if [[ $missingNormlCount -eq 0 && $missingTumorCount -eq 0 ]] ; then
			echo "### All of normal and tumor bams are found"
			echo "### Norml list: $normlList"
			echo "### Tumor list: $tumorList"
		else
			echo "### There were missing things"
			echo "### Normls missing: $missingNormlCount/$normlCount"
			echo "### Tumors missing: $missingTumorCount/$tumorCount"
			((qsubFails++))
			continue
		fi
		cdDir=$runDir/cuffDiff/$cdName
		if [ ! -d $cdDir ] ; then
			mkdir -p $cdDir
		fi
		if [[ -e $cdDir.cuffDiffPass || -e $cdDir.cuffDiffFail || -e $cdDir.cuffDiffInQueue ]] ; then 
			echo "### Cuff diff is already done, failed or inQueue"
			continue
		fi 
		echo "### Submitting $cdName to queue for cuff diff..."
		if [ $rnaTumorStrand == "FIRST" ] ; then
                        echo "##running stranded cuffDiff case"
			sbatch -n 1 -N 1 --cpus-per-task $nCores --export CUFFDIFF2VCFPATH=$cuffdiff2vcfPath,PROCESSCDLISTPATH=$processcdlistPath,CUFFDIFFPATH=$cuffdiffPath,RUNDIR=$runDir,DIRNAME=$cdDir,GTF=$gtf,BAM1="'"$normlList"'",BAM2="'"$tumorList"'",REF=$ref,MASK=$gtfmask,REF=$ref,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_firstStrandedCuffDiff.pbs
			if [ $? -eq 0 ] ; then
				touch $cdDir.cuffDiffInQueue
			else
				((qsubFails++))
			fi
			sleep 2
                elif [ $rnaTumorStrand == "SECOND" ] ; then
                        echo "##running second stranded cuffDiff case"
                        sbatch -n 1 -N 1 --cpus-per-task $nCores --export CUFFDIFF2VCFPATH=$cuffdiff2vcfPath,PROCESSCDLISTPATH=$processcdlistPath,CUFFDIFFPATH=$cuffdiffPath,RUNDIR=$runDir,DIRNAME=$cdDir,GTF=$gtf,BAM1="'"$normlList"'",BAM2="'"$tumorList"'",REF=$ref,MASK=$gtfmask,REF=$ref,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_secondStrandedCuffDiff.pbs
                        if [ $? -eq 0 ] ; then
                                touch $cdDir.cuffDiffInQueue
                        else
                                ((qsubFails++))
                        fi
                        sleep 2
		else
			echo "running unstranded cuffDiff case"
			sbatch -n 1 -N 1 --cpus-per-task $nCores --export CUFFDIFF2VCFPATH=$cuffdiff2vcfPath,PROCESSCDLISTPATH=$processcdlistPath,CUFFDIFFPATH=$cuffdiffPath,RUNDIR=$runDir,DIRNAME=$cdDir,GTF=$gtf,BAM1="'"$normlList"'",BAM2="'"$tumorList"'",REF=$ref,MASK=$gtfmask,REF=$ref,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_cuffDiff.pbs
			if [ $? -eq 0 ] ; then
				touch $cdDir.cuffDiffInQueue
			else
				((qsubFails++))
			fi
			sleep 2
		fi
		;;
	anotherRNAaligner) echo "example RNA aligner"
		;;
	*) echo "I should not be here"
		;;
	esac 
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
