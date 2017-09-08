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

thisStep="pegasus_nextJob_mergeSeuratVcfs.txt"
nxtStep1="pegasus_nextJob_checkProjectComplete.txt"
nxtStep2="pegasus_nextJob_snpEff.txt"
nxtStep3="pegasus_nextJob_alleleCount.txt"

constants=~/jetstream/constants/constants.txt
constantsDir=~/jetstream/constants/
myName=`basename $0 | cut -d_ -f2`

declare -a chrGroups=(1:11:17:21 2:10:16:22 3:9:15:18:MT 4:7:14:Y 5:X:13:20 6:8:12:19)
declare -a chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
echo "array below"
echo "${chrGroups[@]}"
echo "array above"

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
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
seuratPath=`grep @@"$recipe"@@ $constants | grep @@SEURATPATH= | cut -d= -f2`
indels=`grep @@"$recipe"@@ $constants | grep @@INDELS= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
	echo "### DNA pair line is $dnaPairLine"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
	usableName=${sampleNames//,/-}

	pair1=`echo $sampleNames | cut -d, -f1`
	pair2=`echo $sampleNames | cut -d, -f2`

	trackName="$runDir/seurat/$usableName/$usableName"
	mergedVCF=$trackName.seurat.vcf
	if [[ -e $mergedVCF.mergeVcfPass || -e $mergedVCF.mergeVcfFail || -e $mergedVCF.mergeVcfInQueue ]] ; then
		echo "Merge of seurat output already passed, failed or in queue"
		continue
	fi
	missingGrp=0
	for chrGrp in "${chrGroups[@]}"
	do
		grpName=`echo $chrGrp | cut -d: -f1`
		outVCFname=$trackName-group$grpName.seuratOnJIR.vcf
		if [ ! -e $outVCFname.seuratOnJIR-group$grpName-Pass ] ; then 
			echo "missing this: $outVCFname.seuratOnJIR-group$grpName-Pass"
			((missingGrp++))
		fi 
	done
	if [ $missingGrp -eq 0 ] ; then
		echo "All 6 groups must have passed"
		fileList=""
		for thisChr in "${chrs[@]}"
		do
			thisVCF=`ls $trackName-group*.seuratOnJIR.vcf-chr$thisChr`
			#echo "thisVCF is $thisVCF"
			fileList="$fileList I=$thisVCF"
		done
		sbatch -n 1 -N 1 --cpus-per-task $nCores --export FILELIST="$fileList",NXT1=$nxtStep1,NXT2=$nxtStep2,NXT3=$nxtStep3,MERGEDVCF=$mergedVCF,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_mergeVCFs.pbs
		if [ $? -eq 0 ] ; then
			touch $mergedVCF.mergeVcfInQueue
		else
			((qsubFails++))
		fi
		sleep 2

	else
		echo "Not all chr groups have passed seurat"
		((qsubFails++))
	fi 
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
