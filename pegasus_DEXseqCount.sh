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

thisStep="pegasus_nextJob_DEXseqCount.txt"
nxtStep1="pegasus_nextJob_DEXseq.txt"

constantsDir=${JETSTREAM_HOME}/centralPipe/constants/
constants=${JETSTREAM_HOME}/centralPipe/constants/constants.txt

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
DEXseqPath=`grep "@@"$recipe"@@" $constants | grep @@DEXSEQPATH= | cut -d= -f2`
DEXseqGff=`grep "@@"$recipe"@@" $constants | grep @@DEXSEQGFF= | cut -d= -f2`
DEXseqCountPath=`grep "@@"$recipe"@@" $constants | grep @@DEXSEQCOUNTPATH= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
picardPath=`grep @@"$recipe"@@ $constants | grep @@PICARDPATH= | cut -d= -f2`
dexSeq=`grep @@"$recipe"@@ $constants | grep @@DEXSEQ= | cut -d= -f2`
echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
d=`echo $runDir | cut -c 2-`

qsubFails=0

for rnaPairLine in `cat $configFile | grep ^RNAPAIR=`
do
    if [[ $dexSeq == no ]] ; then
                echo "### DEXSEQ not requested for this recipe"
                continue
        fi
    echo "### RNA pair line is $rnaPairLine"
        sampleNames=`echo $rnaPairLine | cut -d= -f2`
        normls=`echo $sampleNames | cut -d, -f1`
        tumors=`echo $sampleNames | cut -d, -f2`
        altNam=`echo $sampleNames | cut -d, -f3`
        DEXseqCountName=`echo $sampleNames | cut -d, -f3`
        if [ "$DEXseqCountName" == "" ] ; then
                normlsWithDash=${normls//;/-}
                tumorsWithDash=${tumors//;/-}
                DEXseqCountName="$normlsWithDash-VS-$tumorsWithDash"
        fi

        echo "### Normals: $normls"
        echo "### Tumors : $tumors"
    tumorCount=0
        normlCount=0
        missingTumorCount=0
        missingNormlCount=0
        normlList=""
        tumorList=""
        normlList2=""
        tumorList2=""

    DEXseqCountDir=$runDir/DEXseq/$DEXseqCountName
    if [ ! -d $DEXseqCountDir ] ; then
        mkdir -p $DEXseqCountDir
    fi
    countsPath=${DEXseqCountDir}/${DEXseqCountName}.can.dexSeqCounts.txt

    for eachNorml in ${normls//;/ }
    do
        ((normlCount++))
        echo "eachnorml: $eachNorml"
        sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachNorml"'"'`
        kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
        samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
        rnaBam=$runDir/$kitName/$samName/${samName}.starDir/${samName}.proj.Aligned.out.sorted.md.bam
        rnaMarkDupPass=$runDir/$kitName/$samName/${samName}.starDir/${samName}.proj.Aligned.out.sorted.bam.rnaMarkDupPass
        if [ ! -e $rnaMarkDupPass ] ; then
            echo "### Can't find rna mark dup pass"
            ((missingNormlCount++))
        else
            DEXseqOut=$DEXseqCountDir/$samName
            DEXseqCountOut=$DEXseqCountDir/$samName.dexSeqCounts.txt

            if [[ -e $DEXseqOut.DEXseqCountPass || -e $DEXseqOut.DEXseqCountFail || -e $DEXseqOut.DEXseqCountInQueue ]] ; then
                echo "### DEXseqCount is already done, failed or inQueue"
            else
                echo "### Submitting $samName to queue for DEXseqCount..."
                sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --mem-per-cpu 4096 --export ALL,DEXSEQOUT=$DEXseqOut,DEXSEQCOUNTPATH=$DEXseqCountPath,DEXSEQGFF=$DEXseqGff,RNABAM=$rnaBam,RUNDIR=$runDir,DEXSEQOUTDIR=$DEXseqCountDir,DEXSEQCOUNTOUT=$DEXseqCountOut,NXT1=$nxtStep1,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_DEXseqCount.sh
                if [ $? -eq 0 ] ; then
                    touch $DEXseqOut.DEXseqCountInQueue
                else
                    ((qsubFails++))
                fi
            fi
        fi
    done
    for eachTumor in ${tumors//;/ }
    do
        ((tumorCount++))
        echo "eachtumor: $eachTumor"
        sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachTumor"'"'`
        kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
        samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
        rnaBam=$runDir/$kitName/$samName/${samName}.starDir/${samName}.proj.Aligned.out.sorted.md.bam
        rnaMarkDupPass=$runDir/$kitName/$samName/${samName}.starDir/${samName}.proj.Aligned.out.sorted.bam.rnaMarkDupPass
        if [ ! -e $rnaMarkDupPass ] ; then
            echo "### Can't find tumor rna mark dup pass"
            ((missingTumorCount++))
        else
            DEXseqOut=$DEXseqCountDir/$samName
            DEXseqCountOut=$DEXseqCountDir/$samName.dexSeqCounts.txt

            if [[ -e $DEXseqOut.DEXseqCountPass || -e $DEXseqOut.DEXseqCountFail || -e $DEXseqOut.DEXseqCountInQueue ]] ; then
                echo "### DEXseqCount is already done, failed or inQueue"
            else
                echo "### Submitting $samName to queue for DEXseqCount..."
                sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ALL,DEXSEQOUT=$DEXseqOut,DEXSEQCOUNTPATH=$DEXseqCountPath,DEXSEQGFF=$DEXseqGff,RUNDIR=$runDir,RNABAM=$rnaBam,DEXSEQOUTDIR=$DEXseqCountDir,DEXSEQCOUNTOUT=$DEXseqCountOut,NXT1=$nxtStep1,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_DEXseqCount.sh

                if [ $? -eq 0 ] ; then
                    touch $DEXseqOut.DEXseqCountInQueue
                else
                    ((qsubFails++))
                fi
            fi
        fi
    done
    if [[ $missingNormlCount -eq 0 && $missingTumorCount -eq 0 ]] ; then
        echo "### All of normal and tumor bams are found"
    else
        echo "### There were missing things"
        ((qsubFails++))
        continue
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
