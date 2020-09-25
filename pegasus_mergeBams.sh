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

thisStep="pegasus_nextJob_mergeBams.txt"
nxtStep1="pegasus_nextJob_markDups.txt"

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
qsubFails=0
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

picardPath=`grep @@"$recipe"@@ $constants | grep @@PICARDPATH= | cut -d= -f2`
samtoolsPath=`grep @@"$recipe"@@ $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
irRequested=`grep "@@"$recipe"@@" $constants | grep @@INDELREALIGN= | cut -d= -f2`
rcRequested=`grep "@@"$recipe"@@" $constants | grep @@RECALIBRATE= | cut -d= -f2`

#optionSum=0
#if [ "$irRequested" == "yes" ] ; then
#    optionSum=`echo "$optionSum + 1" | bc `
#fi
#if [ "$rcRequested" == "yes" ] ; then
#    optionSum=`echo "$optionSum + 2" | bc `
#fi
#echo "### Option sum is $optionSum"
#case $optionSum in
#0) echo "### No ir, no rc"
#    inputNextExt="bam"
#    prevStepExt="dnaAlignPass"
#    prevStepFile="bam"
#;;
#1) echo "### Yes ir, no rc"
#    inputNextExt="ir.bam"
#    prevStepExt="indelRealignPass"
#    prevStepFile="bam"
#;;
#2) echo "### No ir, yes rc"
#    inputNextExt="rc.bam"
#    prevStepExt="recalibratePass"
#    prevStepFile="bam"
#;;
#3) echo "### Yes ir, yes rc"
#    inputNextExt="ir.rc.bam"
#    prevStepExt="recalibratePass"
#    prevStepFile="ir.bam"
#;;
#*) echo "### I shouldnt be here"
#;;
#esac
echo "### projName: $projName"
echo "### confFile: $configFile"

skipLines=1
count=0
sampleCount=0
missingBamTotal=0
d=`echo $runDir | cut -c 2-`

for configLine in `cat $configFile`
do
    if [ "$configLine" == "=START" ] ; then
        skipLines=0
        continue
    fi
    if [ $skipLines -eq 0 ] ; then
        if [[ $configLine == SAMPLE* || $configLine == =END* ]] ; then

            optionSum=0
            if [ "$irRequested" == "yes" ] ; then
                optionSum=`echo "$optionSum + 1" | bc `
            fi
            if [ "$rcRequested" == "yes" ] ; then
                optionSum=`echo "$optionSum + 2" | bc `
            fi
            echo "### Option sum is $optionSum"
            case $optionSum in
            0) echo "### No ir, no rc"
                inputNextExt="bam"
                prevStepExt="dnaAlignPass"
                prevStepFile="bam"
            ;;
            1) echo "### Yes ir, no rc"
                inputNextExt="ir.bam"
                prevStepExt="indelRealignPass"
                prevStepFile="bam"
            ;;
            2) echo "### No ir, yes rc"
                inputNextExt="rc.bam"
                prevStepExt="recalibratePass"
                prevStepFile="bam"
            ;;
            3) echo "### Yes ir, yes rc"
                inputNextExt="ir.rc.bam"
                prevStepExt="recalibratePass"
                prevStepFile="ir.bam"
            ;;
            *) echo "### I shouldnt be here"
            ;;
            esac

            #echo "config line is $configLine"
            arrayCount=${#mergeArray[@]}
            if [ $arrayCount -gt 0 ] ; then
                echo "### Starting with $samName"
                missingBamSample=0
                accountBamSample=0
                bamList=""
                bamCount=0
                ((sampleCount++))
                #echo "### Sample with $arrayCount rows found for kit: $kitName, sample: $samName."
                lastItemIndex=`echo ${#mergeArray[@]}-1 | bc`
                mergedBamName="$runDir/$kitName/$samName/$samName.proj.bam"
                checkIfSplit="$runDir/$kitName/$samName/$samName.preMerge.000.bam.processSplit"
                inAPair=`cat $configFile | grep 'DNAPAIR' | grep $samName`
                inAFam=`cat $configFile | grep 'DNAFAMI' | grep $samName`
                if [ -z "$inAPair" ] && [ -z "$inAFam" ] ; then
                    inputNextExt="ir.rc.bam"
                    prevStepExt="recalibratePass"
                    prevStepFile="ir.bam"
                fi
                if [ -e $checkIfSplit ] ; then
                    echo "### Setting previous step extension to mergeBamPass because this bam was processed in pieces"
                    prevStepExt="mergeBamPass"
                    #prevStepFile="$prevStepFile"
                    prevStepFile="$inputNextExt"
                fi
                for (( i=0; i<=$lastItemIndex; i++ ))
                do
                    #echo "array with index $i is::: ${mergeArray[$i]}"
                    ((accountBamSample++))
                    bamPre="$runDir/$kitName/$samName/$samName.preMerge."`printf "%03d" "$i"`
                    #bamName="$runDir/$kitName/$samName/$samName.preMerge."`printf "%03d" "$i"`".bam"
                    if [[ ! -e $bamPre.$prevStepFile.$prevStepExt || ! -e $bamPre.$inputNextExt ]] ; then
                        echo "### File missing: $bamPre.$prevStepFile.$prevStepExt OR $bamPre.$inputNextExt "
                        ((missingBamTotal++))
                        ((missingBamSample++))
                    else
                        echo "### File found $bamPre.$inputNextExt "
                        bamList="I=$bamPre.$inputNextExt $bamList"
                        ((bamCount++))
                    fi
                done
                echo "### Done with $samName with $missingBamSample missing files out of $accountBamSample"
                if [ $missingBamSample -eq 0 ] ; then
                    echo "### All bams were accounted for $samName"
                    echo "### Bam files($bamCount): $bamList"
                    if [[ -e $mergedBamName.mergeBamPass || -e $mergedBamName.mergeBamInQueue || -e $mergedBamName.mergeBamFail ]] ; then
                        echo "### Already passed, inQueue, or failed"
                    else
                        echo "### Ready to submit to create $mergedBamName"
                        sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ALL,PICARDPATH=$picardPath,SAMTOOLSPATH=$samtoolsPath,CNT=$bamCount,RUNDIR=$runDir,NXT1=$nxtStep1,BAMLIST="$bamList",MERGEDBAM=$mergedBamName,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_mergeBams.sh
                        if [ $? -eq 0 ] ; then
                            touch $mergedBamName.mergeBamInQueue
                        else
                            echo "### There was a failure in qsub command"
                            ((qsubFails++))
                        fi
                    fi
                else
                    echo "### Some files were missing for $samName"
                fi
            fi
            kitName=`echo $configLine | cut -d= -f2 | cut -d, -f1`
            samName=`echo $configLine | cut -d= -f2 | cut -d, -f2`
            unset mergeArray
            count=0
            continue
        else #doesnt start with =, add to mergeArray
            #echo "adding $configLine to mergeArray"
            mergeArray[$count]=$configLine
            ((count++))
        fi
    else
        continue
    fi
done

if [ $qsubFails -eq 0 ] ; then
#all jobs submitted succesffully, remove this dir from messages
    echo "### Removing $thisStep from $runDir."
    rm -f $runDir/$thisStep
else
#qsub failed at some point, this runDir must stay in messages
    echo "### Failure in qsub. Keeping $thisStep"
fi

time=`date +%d-%m-%Y-%H-%M`
echo "Ending $0 at $time"
