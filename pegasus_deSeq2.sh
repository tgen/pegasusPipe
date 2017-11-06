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

thisStep="pegasus_nextJob_deSeq2.txt"
nxtStep1="pegasus_nextJob_postDeSeq2.txt"

constants=${JETSTREAM_HOME}/centralPipe/constants/constants.txt
constantsDir=${JETSTREAM_HOME}/centralPipe/constants/
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
deseq2Path=`grep @@"$recipe"@@ $constants | grep @@DESEQ2PATH= | cut -d= -f2`
gtfmask=`grep @@"$recipe"@@ $constants | grep @@GTFMASK= | cut -d= -f2`
gtf=`grep @@"$recipe"@@ $constants | grep @@GTF= | cut -d= -f2`
cuffdiff2vcfPath=`grep @@"$recipe"@@ $constants | grep @@CUFFDIFF2VCFPATH= | cut -d= -f2`
processcdlistPath=`grep @@"$recipe"@@ $constants | grep @@PROCESSCDLISTPATH= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
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
    altNam=`echo $sampleNames | cut -d, -f3`
    deSeqName=`echo $sampleNames | cut -d, -f3`
    if [ "$deSeqName" == "" ] ; then
        normlsWithDash=${normls//;/-}
        tumorsWithDash=${tumors//;/-}
        deSeqName="$normlsWithDash-VS-$tumorsWithDash"
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
    case $rnaAligner in
    tophat) echo "tophat case"
        for eachNorml in ${normls//;/ }
        do
            ((normlCount++))
            echo "eachnorml: $eachNorml"
            sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachNorml"'"'`
            kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
            samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
            eachNormlHTseqOut=$runDir/$kitName/$samName/$samName.topHatDir/$samName.proj.accepted_hits.sam.htSeqCounts
            eachNormlHTseqPas=$runDir/$kitName/$samName/$samName.topHatDir/$samName.proj.accepted_hits.sam.htSeqPass
            if [ ! -e $eachNormlHTseqPas ] ; then
                echo "### Can't find normal ht seq pass"
                ((missingNormlCount++))
            else
                normlList="$eachNormlHTseqOut;$normlList"
                normlList2="$eachNorml-$normlList2"
            fi
        done
        for eachTumor in ${tumors//;/ }
        do
            ((tumorCount++))
            echo "eachtumor: $eachTumor"
            sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachTumor"'"'`
            kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
            samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
            eachTumorHTseqOut=$runDir/$kitName/$samName/$samName.topHatDir/$samName.proj.accepted_hits.sam.htSeqCounts
            if [ ! -e $eachTumorHTseqOut ] ; then
                echo "### Can't find tumor ht seq pass"
                ((missingTumorCount++))
            else
                tumorList="$eachTumorHTseqOut;$tumorList"
                tumorList2="$eachTumor-$tumorList2"
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
        #deSeqDir=$runDir/deSeq/$normlList2-VS-$tumorList2.dsDir
        #if [ $altNam != "" ] ; then
        #    deSeqDir=$runDir/deSeq/$altNam
        #fi
        deSeqDir=$runDir/deSeq2/$deSeqName
        if [ ! -d $deSeqDir ] ; then
            mkdir -p $deSeqDir
        fi
        if [[ -e $deSeqDir.deSeq2Pass || -e $deSeqDir.deSeq2Fail || -e $deSeqDir.deSeq2InQueue ]] ; then
            echo "### DeSeq is already done, failed or inQueue"
            continue
        fi

        echo "### Submitting $normlList2-VS-$tumorList2 to queue for deSeq..."
        sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ALL,DESEQ2PATH=$deseq2Path,RUNDIR=$runDir,DIRNAME=$deSeqDir,GTF=$gtf,NORMLIST="$normlList",TUMORLIST="$tumorList",REF=$ref,NXT1=$nxtStep1,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_deSeq2.sh
        if [ $? -eq 0 ] ; then
            touch $deSeqDir.deSeq2InQueue
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
            eachNormlHTseqOut=$runDir/$kitName/$samName/$samName.starDir/$samName.htSeqCounts
            eachNormlHTseqPas=$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sam.htSeqPass
            if [ ! -e $eachNormlHTseqPas ] ; then
                echo "### Can't find normal ht seq pass"
                ((missingNormlCount++))
            else
                normlList="$eachNormlHTseqOut;$normlList"
                normlList2="$eachNorml-$normlList2"
            fi
        done
        for eachTumor in ${tumors//;/ }
        do
            ((tumorCount++))
            echo "eachtumor: $eachTumor"
            sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachTumor"'"'`
            kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
            samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
            eachTumorHTseqOut=$runDir/$kitName/$samName/$samName.starDir/$samName.htSeqCounts
            eachTumorHTseqPas=$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sam.htSeqPass
            if [ ! -e $eachTumorHTseqPas ] ; then
                echo "### Can't find tumor ht seq pass"
                ((missingTumorCount++))
            else
                tumorList="$eachTumorHTseqOut;$tumorList"
                tumorList2="$eachTumor-$tumorList2"
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

        deSeqDir=$runDir/deSeq2/$deSeqName
        if [ ! -d $deSeqDir ] ; then
            mkdir -p $deSeqDir
        fi
        if [[ -e $deSeqDir.deSeq2Pass || -e $deSeqDir.deSeq2Fail || -e $deSeqDir.deSeq2InQueue ]] ; then
            echo "### DeSeq2 is already done, failed or inQueue"
            continue
        fi
        if [ ! -d $deSeqDir ] ; then
            mkdir -p $deSeqDir
        fi
        echo "### Submitting $normlList2-VS-$tumorList2 to queue for deSeq..."
        sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ALL,DESEQ2PATH=$deseq2Path,RUNDIR=$runDir,DIRNAME=$deSeqDir,GTF=$gtf,NORMLIST="$normlList",TUMORLIST="$tumorList",REF=$ref,NXT1=$nxtStep1,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_deSeq2.sh
        if [ $? -eq 0 ] ; then
            touch $deSeqDir.deSeq2InQueue
        else
            ((qsubFails++))
        fi
        sleep 2

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
