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

thisStep="pegasus_nextJob_DEXseq.txt"
nxtStep1="pegasus_nextJob_postDEXseq.txt"

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

ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
DEXseqPath=`grep "@@"$recipe"@@" $constants | grep @@DEXSEQPATH= | cut -d= -f2`
DEXseqGff=`grep "@@"$recipe"@@" $constants | grep @@DEXSEQGFF= | cut -d= -f2`
DEXseqCountPath=`grep "@@"$recipe"@@" $constants | grep @@DEXSEQCOUNTPATH= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
picardPath=`grep @@"$recipe"@@ $constants | grep @@PICARDPATH= | cut -d= -f2`
echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
d=`echo $runDir | cut -c 2-`

qsubFails=0

for rnaPairLine in `cat $configFile | grep ^RNAPAIR=`
do

    echo "### RNA pair line is $rnaPairLine"
        sampleNames=`echo $rnaPairLine | cut -d= -f2`
        normls=`echo $sampleNames | cut -d, -f1`
        tumors=`echo $sampleNames | cut -d, -f2`
        altNam=`echo $sampleNames | cut -d, -f3`
        DEXseqName=`echo $sampleNames | cut -d, -f3`
        
    if [[ "$sampleNames" != *";"* ]] ; then
        echo "There is only 2 RNA samples, at least 3 are needed to run DEXseq at this time"
        echo "will skip DEXseq for $sampleNames"
        continue
    else
        echo "There are at least 3 RNA samples, will continue to run DEXseq"
    fi

    if [ "$DEXseqName" == "" ] ; then
                normlsWithDash=${normls//;/-}
                tumorsWithDash=${tumors//;/-}
                DEXseqName="$normlsWithDash-VS-$tumorsWithDash"
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

    DEXseqDir=$runDir/DEXseq/$DEXseqName
    if [ ! -d $DEXseqDir ] ; then
        mkdir -p $DEXseqDir
    fi
    DEXseqConfig=$DEXseqDir/${DEXseqName}_DEXseq.config
    echo "sampleID    condition    libType    exonCountModel    countFilePath" > $DEXseqConfig

    for eachNorml in ${normls//;/ }
                do
                        ((normlCount++))
                        echo "eachnorml: $eachNorml"
            countsPath=${DEXseqDir}/${eachNorml}.dexSeqCounts.txt
            countsPass=${DEXseqDir}/${eachNorml}.DEXseqCountPass
                        if [ ! -e $countsPass ] ; then
                                echo "### Can't find normal counts pass"
                                ((missingNormlCount++))
                continue
                        else
                echo "$eachNorml    unaffected    paired-end    noAggregate    $countsPath" >> $DEXseqConfig
                        fi
                done
                for eachTumor in ${tumors//;/ }
                do
                        ((tumorCount++))
                        echo "eachtumor: $eachTumor"
            countsPath=${DEXseqDir}/${eachTumor}.dexSeqCounts.txt
                        countsPass=${DEXseqDir}/${eachTumor}.DEXseqCountPass
                        if [ ! -e $countsPass ] ; then
                                echo "### Can't find tumor counts pass"
                                ((missingTumorCount++))
                continue
                        else
                echo "$eachTumor    affected    paired-end    noAggregate    $countsPath" >> $DEXseqConfig
                        fi
                done
                if [[ $missingNormlCount -eq 0 && $missingTumorCount -eq 0 ]] ; then
                        echo "### All of normal and tumor bams are found"
                        echo "### Norml list: $normlList"
                        echo "### Tumor list: $tumorList"
                        normlList="${normlList%?}"
                        tumorList="${tumorList%?}"
                        normlList2="${normlList2%?}"
                        tumorList2="${tumorList2%?}"
                else
                        echo "### There were missing things"
                        echo "### Normls missing: $missingNormlCount/$normlCount"
                        echo "### Tumors missing: $missingTumorCount/$tumorCount"
                        ((qsubFails++))
                        continue
                fi

        DEXseqOut=$runDir/DEXseq/$DEXseqName/${DEXseqName}_DEXseq.txt
        objectData=$runDir/DEXseq/$DEXseqName/${DEXseqName}_DEXseq.rda

                if [[ -e $DEXseqDir.DEXseqPass || -e $DEXseqDir.DEXseqFail || -e $DEXseqDir.DEXseqInQueue ]] ; then
                        echo "### DEXseq is already done, failed or inQueue"
                        continue
                fi
                echo "### Submitting $normlList2-VS-$tumorList2 to queue for DEXseq..."
                sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export DEXSEQPATH=$DEXseqPath,DEXSEQGFF=$DEXseqGff,DEXSEQCONFIG=$DEXseqConfig,RUNDIR=$runDir,DEXSEQOUTDIR=$DEXseqDir,OBJECTDATA=$objectData,DEXSEQOUTFILE=$DEXseqOut,KALLISTOOUT=$kallistoDir,NXT1=$nxtStep1,D=$d $pegasusPbsHome/pegasus_DEXseq.sh
                if [ $? -eq 0 ] ; then
                       touch $DEXseqDir.DEXseqInQueue
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
