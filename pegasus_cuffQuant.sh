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

thisStep="pegasus_nextJob_cuffQuant.txt"
nxtStep1="pegasus_nextJob_postCuffQuant.txt"

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
params=`grep @@${myName}_PARAMS= $constantsDir/$recipe | cut -d= -f2`
cuffQuant=`grep @@CUFFQUANT= $constantsDir/$recipe | cut -d= -f2`
usegtf=`cat $configFile | grep "^CUFFLINKUSEGTF=" | cut -d= -f2 | head -1 | tr -d [:space:]`
usemask=`grep @@"$recipe"@@ $constants | grep @@CUFFLINKUSEMASK= | cut -d= -f2`
ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
cuffQuantPath=`grep @@"$recipe"@@ $constants | grep @@CUFFQUANTPATH= | cut -d= -f2`
gtfmask=`grep @@"$recipe"@@ $constants | grep @@GTFMASK= | cut -d= -f2`
gtf=`grep @@"$recipe"@@ $constants | grep @@GTF= | cut -d= -f2`

echo "gtfmask is now $gtfmask"

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###
for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
    if [ $cuffQuant != "yes" ] ; then
        echo "### CuffQuant not requested for this recipe"
        continue
    fi
    echo "### Sample is $sampleLine"
    kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
    samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
    assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
    libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
    rnaStrand=`grep "@@"$kitName"@@" $constants | cut -d= -f2`

    echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID, rnaStrand: $rnaStrand"
    if [ "$assayID" != "RNA" ] ; then
        echo "### Assay ID is $assayID. Skipping."
        continue
    fi

    #need to add option to run on star output
    case $rnaAligner in
    tophat) echo "tophat case"
        topHatDir="$runDir/$kitName/$samName/$samName.topHatDir"
        accHitsBam="$topHatDir/$samName.proj.accepted_hits.bam"
        echo "### My bam is $accHitsBam"
        if [ ! -e $topHatDir.thPass ] ; then
            echo "### Looks like tophat is not done yet. $topHatDir.thPass doesnt exist yet"
            ((qsubFails++))
            continue
        fi
        if [[ -e $topHatDir.cuffQuantPass || -e $topHatDir.cuffQuantFail || -e $topHatDir.cuffQuantInQueue ]] ; then
            echo "### Cuffquant is already done, failed or inQueue"
            continue
        fi

        echo "### Submitting $topHatDir to queue for cuff quant..."
        sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export PARAMS=${params},DIRNAME=$topHatDir,CUFFQUANTPATH=$cuffQuantPath,REF=$ref,BAM=$accHitsBam,USEGTF=$usegtf,USEMASK=$usemask,CUFFLINKGTF=$gtf,CUFFLINKMASK=$gtfmask,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_cuffQuant.sh
        if [ $? -eq 0 ] ; then
            touch $topHatDir.cuffQuantInQueue
        else
            ((qsubFails++))
        fi
        sleep 2
        ;;
    star) echo "star case"
        starDir="$runDir/$kitName/$samName/$samName.starDir"
        starBam="$runDir/$kitName/$samName/$samName.starDir/$samName.proj.Aligned.out.sorted.md.bam"
        if [ ! -e $starDir/$samName.proj.Aligned.out.sorted.bam.rnaMarkDupPass ] ; then
            echo "### Looks like rna mark dup is not done yet. $starDir/$samName.proj.Aligned.out.sorted.bam.rnaMarkDupPass"
            ((qsubFails++))
            continue
        fi
        if [[ -e $starDir.cuffQuantPass || -e $starDir.cuffQuantFail || -e $starDir.cuffQuantInQueue ]] ; then
            echo "### Cuffquant is already done, failed or inQueue"
            continue
        fi

        echo "### Submitting $starDir to queue for cuff quant..."
        if [ $rnaStrand == "FIRST" ] ; then
                        echo "##running stranded cuffQuant case"
            sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export PARAMS=$params,DIRNAME=$starDir,CUFFQUANTPATH=$cuffQuantPath,REF=$ref,BAM=$starBam,USEGTF=$usegtf,USEMASK=$usemask,CUFFLINKGTF=$gtf,CUFFLINKMASK=$gtfmask,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_firstStrandedCuffQuant.sh
            if [ $? -eq 0 ] ; then
                touch $starDir.cuffQuantInQueue
            else
                ((qsubFails++))
            fi
            sleep 2
                elif [ $rnaStrand == "SECOND" ] ; then
                        echo "##running second stranded cuffQuant case"
                        sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export PARAMS=$params,DIRNAME=$starDir,CUFFQUANTPATH=$cuffQuantPath,REF=$ref,BAM=$starBam,USEGTF=$usegtf,USEMASK=$usemask,CUFFLINKGTF=$gtf,CUFFLINKMASK=$gtfmask,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_secondStrandedCuffQuant.sh
                        if [ $? -eq 0 ] ; then
                                touch $starDir.cuffQuantInQueue
                        else
                                ((qsubFails++))
                        fi
                        sleep 2
        else
            echo "##running unstranded cuffQuant case"
            sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export PARAMS=$params,DIRNAME=$starDir,CUFFQUANTPATH=$cuffQuantPath,REF=$ref,BAM=$starBam,USEGTF=$usegtf,USEMASK=$usemask,CUFFLINKGTF=$gtf,CUFFLINKMASK=$gtfmask,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_cuffQuant.sh
            if [ $? -eq 0 ] ; then
                touch $starDir.cuffQuantInQueue
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
