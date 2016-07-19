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

thisStep="pegasus_nextJob_seuratSwapCheck.txt"
nxtStep1="pegasus_nextJob_checkProjectComplete.txt"
#nxtStep2="pegasus_nextJob_seuratSwapCheck.txt"
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
ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
seuratPath=`grep @@"$recipe"@@ $constants | grep @@SEURATPATH= | cut -d= -f2`
chrList=`grep @@"$recipe"@@ $constants | grep @@CHRLIST= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

qsubFails=0
###
##first for loop is for finding REVERSE seurat vcfs
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
        echo "### DNA pair line is $dnaPairLine for reverse seurat stuff"
        sampleNames=`echo $dnaPairLine | cut -d= -f2`
        for eachSample in ${sampleNames//,/ }
        do
                ((sampleCount++))
                #echo "eachsample: $eachSample"
                sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
                kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
                samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
                assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
        done
        usableName=${sampleNames//,/-}
        usableName=${usableName}_REVERSE
        sampleCount=0
        missingSampleCount=0
        sampleList=""
        trackName="$runDir/seurat/$usableName/$usableName"
        if [ ! -e $trackName.REVseuratPass ] ; then
                echo "### Seurat Pass doesnt exist yet: $trackName.REVseuratPass"
                ((qsubFails++))
                continue
        fi
        if [[ -e $trackName.REVseurat.vcf.seuratSwapCheckInQueue || -e $trackName.REVseurat.vcf.seuratSwapCheckPass || -e $trackName.REVseurat.vcf.seuratSwapCheckFail ]] ; then
                echo "### seuratSwapCheck is already done, failed or inQueue"
                continue
        fi
        echo "### Submitting $trackName.REVseurat.vcf to queue for seuratSwapCheck..."
        qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$trackName.REVseurat.vcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_seuratSwapCheck.pbs
        if [ $? -eq 0 ] ; then
                touch $trackName.REVseurat.vcf.seuratSwapCheckInQueue
        else
                ((qsubFails++))
        fi
        sleep 2
done
####first for loop is for finding seurat vcfs
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
        echo "### DNA pair line is $dnaPairLine for seurat stuff"
        sampleNames=`echo $dnaPairLine | cut -d= -f2`
        for eachSample in ${sampleNames//,/ }
        do
                ((sampleCount++))
                #echo "eachsample: $eachSample"
                sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
                kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
                samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
                assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
        done
        usableName=${sampleNames//,/-}
        sampleCount=0
        missingSampleCount=0
        sampleList=""
        trackName="$runDir/seurat/$usableName/$usableName"
        if [ ! -e $trackName.seuratPass ] ; then
                echo "### Seurat Pass doesnt exist yet: $trackName.seuratPass"
                ((qsubFails++))
                continue
        fi
        if [[ -e $trackName.seurat.vcf.seuratSwapCheckInQueue || -e $trackName.seurat.vcf.seuratSwapCheckPass || -e $trackName.seurat.vcf.seuratSwapCheckFail ]] ; then
                echo "### seuratSwapCheck is already done, failed or inQueue"
                continue
        fi
        echo "### Submitting $trackName.seurat.vcf to queue for seuratSwapCheck..."
        qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$trackName.seurat.vcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_seuratSwapCheck.pbs
        if [ $? -eq 0 ] ; then
                touch $trackName.seurat.vcf.seuratSwapCheckInQueue
        else
                ((qsubFails++))
        fi
        sleep 2
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
