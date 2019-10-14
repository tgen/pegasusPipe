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

thisStep="pegasus_nextJob_indelRealign.txt"
nxtStep1="pegasus_nextJob_recalibrate.txt"

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

dnaAligner=`grep "@@"$recipe"@@" $constants | grep @@DNAALIGNER= | cut -d= -f2`
ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
gatkPath=`grep "@@"$recipe"@@" $constants | grep @@GATKPATH= | cut -d= -f2`
indels=`grep @@"$recipe"@@ $constants | grep @@INDELS= | cut -d= -f2`
bwaPath=`grep "@@"$recipe"@@" $constants | grep @@BWAPATH= | cut -d= -f2`
irRequested=`grep "@@"$recipe"@@" $constants | grep @@INDELREALIGN= | cut -d= -f2`
faiFile="$ref.fai"
indel=`grep @@INDELREALIGN= $constantsDir/$recipe | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### dnaAlign: $dnaAligner"
d=`echo $runDir | cut -c 2-`

skipLines=1
skipGroup=0
qsubFails=0
###
for configLine in `cat $configFile`
do
    if [ "$configLine" == "=START" ] ; then
        skipLines=0
        continue
    fi
    if [ $skipLines -eq 0 ] ; then
        if [[ $configLine == SAMPLE* || $configLine == =END* ]] ; then
            #echo "config line is $configLine"
            arrayCount=${#mergeArray[@]}
            if [ $arrayCount -gt 0 ] ; then
                ((sampleCount++))
                echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID, insSize: $insSize"
                if [ "$assayID" == "RNA" ] ; then
                    echo "### Assay ID is $assayID. Will skip."
                    skipGroup=1
                fi
                if [ "$irRequested" == "no" ] ; then
                    echo "### Indel realign is $irRequested."
                    inAPair=`cat $configFile | grep 'DNAPAIR' | grep $samName`
                    inAFam=`cat $configFile | grep 'DNAFAMI' | grep $samName`
                    if [ -z "$inAPair" ] && [ -z "$inAFam" ] ; then
                        echo "This sample isn't paired so will do a solo indel realign"
                    else
                        skipGroup=1
                    fi
                fi
                lastItemIndex=`echo ${#mergeArray[@]}-1 | bc`
                if [ $skipGroup -ne 1 ] ; then
                    for (( i=0; i<=$lastItemIndex; i++ ))
                    do
                        #echo "array with index $i is::: ${mergeArray[$i]}"
                        bamPre="$runDir/$kitName/$samName/$samName.preMerge."`printf "%03d" "$i"`
                        bamName="$runDir/$kitName/$samName/$samName.preMerge."`printf "%03d" "$i"`".bam"
                        irBamFile=${bamName/.bam/.ir.bam}
                        irIntFile=${bamName/.bam/.bam.intervals}
                        #echo "### BAMPRE: $bamPre"
                        if [ ! -e $bamName.dnaAlignPass ] ; then
                            echo "### This bam file is not done with alignment yet"
                            continue
                        fi
                        if [ -e $bamName.processSplit ] ; then
                            echo "### This bam file needs to be processed in pieces."
                            continue
                        fi
                        if [[ -e $bamName.indelRealignPass || -e $bamName.indelRealignInQueue || -e $bamName.indelRealignFail ]] ; then
                            echo "### This bam file already passed, failed, or inQueue"
                            continue
                        fi
                        d=`echo $runDir | cut -c 2-`
			if [ $indel != "yes" ] ; then
        		    echo "Indel realignment not requested for this recipe, skipping job submission..."
			    echo "Renaming bam file to bypass naming scheme validation in later steps..."
			    echo "mv $bamName $irBamFile"
			    mv $bamName $irBamFile
			    echo "Touching ${bamName}.indelRealignPass"
			    touch ${bamName}.indelRealignPass
			    echo "Touching ${runDir}/${nxtStep1}"
    			    touch ${runDir}/${nxtStep1}
        		    continue
    			fi
                        echo "### Submitting to indel realign to create $bamName"
			echo "sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --mem 128G --export ALL,GATKPATH=$gatkPath,INTS=$irIntFile,IRBAMFILE=$irBamFile,D=$d,INDELS=$indels,REF=$ref,BAMFILE=$bamName,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_indelRealign.sh"
                        sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --mem 128G --export ALL,GATKPATH=$gatkPath,INTS=$irIntFile,IRBAMFILE=$irBamFile,D=$d,INDELS=$indels,REF=$ref,BAMFILE=$bamName,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_indelRealign.sh
                        if [ $? -eq 0 ] ; then
                            touch $bamName.indelRealignInQueue
                        else
                            ((qsubFails++))
                        fi
                        sleep 1
                    done
                else #else of skipGroup
                    echo "### Skipping group."
                fi #end of skipGroup
            fi
            kitName=`echo $configLine | cut -d= -f2 | cut -d, -f1`
            samName=`echo $configLine | cut -d= -f2 | cut -d, -f2`
            assayID=`echo $configLine | cut -d= -f2 | cut -d, -f3`
            libraID=`echo $configLine | cut -d= -f2 | cut -d, -f4`
            unset mergeArray
            count=0
            skipGroup=0
            continue
        else #doesnt start with =, add to mergeArray
            mergeArray[$count]=$configLine
            ((count++))
        fi
    else
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
