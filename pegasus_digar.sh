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

thisStep="pegasus_nextJob_digar.txt"
nxtStep1="pegasus_nextJob_digarPost.txt"
pbsHome="~/pegasus-pipe/jobScripts"
constants="~/central-pipe/constants/constants.txt"
constantsDir="~/central-pipe/constants"
myName=`basename $0 |  cut -d_ -f2`

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
#params=`grep @@${myName}_PARAMS= $constantsDir/$recipe | cut -d= -f2`


#usegtf=`cat $configFile | grep "^CUFFLINKUSEGTF=" | cut -d= -f2 | head -1 | tr -d [:space:]`
#usemask=`cat $configFile | grep "^CUFFLINKUSEMASK=" | cut -d= -f2 | head -1 | tr -d [:space:]`

ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
#cufflinksPath=`grep @@"$recipe"@@ $constants | grep @@CUFFLINKSPATH= | cut -d= -f2`
#gtfmask=`grep @@"$recipe"@@ $constants | grep @@GTFMASK= | cut -d= -f2`
gtf=`grep @@"$recipe"@@ $constants | grep @@GTF= | cut -d= -f2`
bwaPath=`grep "@@"$recipe"@@" $constants | grep @@BWAPATH= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
digarPath=`grep "@@"$recipe"@@" $constants | grep @@DIGARPATH= | cut -d= -f2`
digarAnn=`grep "@@"$recipe"@@" $constants | grep @@DIGARANN= | cut -d= -f2`
listOfGenes=`grep "@@"$recipe"@@" $constants | grep @@DIGARGENES= | cut -d= -f2`
trinityPath=`grep "@@"$recipe"@@" $constants | grep @@TRINITYPATH= | cut -d= -f2`
digar=`grep @@DIGAR= $constantsDir/$recipe | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
echo "### list of genes is: $listOfGenes"
echo "### digarPath is: $digarPath"
echo "### digarAnn is: $digarAnn"
echo "### digar is $digar"

d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###
for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
       if [[ "$digar" != "yes" ]] ; then
                echo "digar not requested for this recipe"
		echo "### I should remove $thisStep from $runDir."
        	rm -f $runDir/$thisStep
		echo "### Exiting!!!"
		exit
       else
		echo "### Sample is $sampleLine"
		kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
		samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
		assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
		libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
		echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID"
		if [ "$assayID" != "RNA" ] ; then
			echo "### Assay ID is $assayID. Skipping."
			continue
		fi
		read1Name=$runDir/$kitName/$samName/$samName.proj.R1.fastq.gz
		read2Name=$runDir/$kitName/$samName/$samName.proj.R2.fastq.gz
		echo "read 1 name: $read1Name"
		echo "read 2 name: $read2Name"
		
		case $rnaAligner in 
		tophat) echo "tophat case"
			echo "digar not meant for tophat bams"
			continue
			;;
		star) echo "star case"
			starDir="$runDir/$kitName/$samName/$samName.starDir"
			digarDir="$runDir/$kitName/$samName/$samName.digarDir"
			starBam="$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sorted.md.bam"
			rnaDir="$runDir/$kitName/$samName"
			if [ ! -e $starDir/$samName.proj.Aligned.out.sorted.bam.rnaMarkDupPass ] ; then
				echo "### Looks like rna mark dup is not done yet. $starDir/$samName.proj.Aligned.out.sorted.bam.rnaMarkDupPass"
				((qsubFails++))
				continue
			fi
			if [ ! -d $digarDir ] ; then
				mkdir -p $digarDir
			fi

			#cp $read1Name $digarDir
			#cp $read2Name $digarDir

			#creating linked files to the original reads
			#ir1Name=`basename $read1Name`
			#r2Name=`basename $read2Name`
			#cd $digarDir
			#ln -s $read1Name $r1Name
			#read1Name=$digarDir/$r1Name
			#ln -s $read2Name $r2Name
			#read2Name=$digarDir/$r2Name
			#cd -
			#done creating links. vars for reads changed.
			echo "### read 1 name: $read1Name"
			echo "### read 2 name: $read2Name"
			
			while read line	
			do
				geneName=$line
				echo "The gene is now $geneName"
				if [[ -e ${digarDir}/$samName.$geneName.digarPass || -e ${digarDir}/$samName.$geneName.digarFail || -e ${digarDir}/$samName.$geneName.digarInQueue ]] ; then 
					echo "### Digar is already done, failed or inQueue"
					continue
				fi 

				echo "### Submitting $digarDir to queue for digar..."
				echo "sbatch -n 1 -N 1 --cpus-per-task $nCores -v FASTQ1=$read1Name,FASTQ2=$read2Name,SAMNAME=$samName,TRINITYPATH=$trinityPath,DIGARPATH=$digarPath,ANN=$digarAnn,SAMTOOLSPATH=$samtoolsPath,BWAPATH=$bwaPath,GENENAME=$geneName,DIGARDIR=$digarDir,REF=$ref,BAM=$starBam,GTF=$gtf,RUNDIR=$runDir,D=$d $pbsHome/pegasus_digar.pbs"
				sbatch -n 1 -N 1 --cpus-per-task $nCores -v FASTQ1=$read1Name,FASTQ2=$read2Name,LISTOFGENES=$listOfGenes,SAMNAME=$samName,TRINITYPATH=$trinityPath,DIGARPATH=$digarPath,ANN=$digarAnn,SAMTOOLSPATH=$samtoolsPath,BWAPATH=$bwaPath,GENENAME=$geneName,DIGARDIR=$digarDir,REF=$ref,NXT1=$nxtStep1,BAM=$starBam,GTF=$gtf,RUNDIR=$runDir,D=$d $pbsHome/pegasus_digar.pbs
				if [ $? -eq 0 ] ; then
					touch ${digarDir}/$samName.$geneName.digarInQueue
				else
					((qsubFails++))
				fi
				sleep 2

			done < $listOfGenes
			
		;;
		anotherRNAaligner) echo "example RNA aligner"
			;;
		*) echo "I should not be here"
			;;
		esac 
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
