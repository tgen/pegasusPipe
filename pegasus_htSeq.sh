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

thisStep="pegasus_nextJob_htSeq.txt"
nxtStep1="pegasus_nextJob_postHtSeq.txt"

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


ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
picardPath=`grep @@"$recipe"@@ $constants | grep @@PICARDPATH= | cut -d= -f2`
samtoolsPath=`grep @@"$recipe"@@ $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
gtf=`grep @@"$recipe"@@ $constants | grep @@GTF= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###
for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
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
    #enter case statement here
    case $rnaAligner in
    tophat) echo "### Tophat case"
            topHatDir="$runDir/$kitName/$samName/$samName.topHatDir"
            accHitsBam="$topHatDir/$samName.proj.accepted_hits.bam"
            accHitsSam="$topHatDir/$samName.proj.accepted_hits.sam"
            echo "### My bam is $accHitsBam"
            if [ ! -e $topHatDir.thPass ] ; then
                echo "### Looks like tophat is not done yet. $topHatDir.thPass doesnt exist yet"
                ((qsubFails++))
                continue
            fi
            if [ ! -e $accHitsBam ] ; then
                echo "### Weird. Bam itself is missing: $accHitsBam"
                ((qsubFails++))
                continue
            fi
            if [[ -e $accHitsSam.htSeqPass || -e $accHitsSam.htSeqFail || -e $accHitsSam.htSeqInQueue ]] ; then
                echo "### HT Seq is already done, failed or inQueue"
                continue
            fi
            echo "### Submitting $accHitsSam to queue for HT Seq..."
            sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ALL,SAMTOOLSPATH=$samtoolsPath,PICARDPATH=$picardPath,SAM=$accHitsSam,BAM=$accHitsBam,GTF=$gtf,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_htSeq.sh
            if [ $? -eq 0 ] ; then
                touch $accHitsSam.htSeqInQueue
            else
                ((qsubFails++))
            fi
        sleep 2
        ;;
    star) echo "### Star case"
            starDir="$runDir/$kitName/$samName/$samName.starDir"
            rmDupPass="$starDir/$samName.proj.Aligned.out.sorted.bam.rnaMarkDupPass"
            alignedBam="$starDir/$samName.proj.Aligned.out.sorted.md.bam"
            alignedSam="$starDir/$samName.proj.Aligned.out.sam"
            echo "### My bam is $alignedBam"
            if [ ! -e $starDir.starPass ] ; then
                echo "### Looks like star alignment is not done yet. $starPass doesnt exist yet"
                ((qsubFails++))
                continue
            fi
            if [ ! -e $alignedSam ] ; then
                echo "### Weird. Sam itself is missing: $alignedSam"
                ((qsubFails++))
                continue
            fi
            if [[ -e $alignedSam.htSeqPass || -e $alignedSam.htSeqFail || -e $alignedSam.htSeqInQueue ]] ; then
                echo "### HT Seq is already done, failed or inQueue"
                continue
            fi
            echo "### Submitting $alignedSam to queue for HT Seq..."
             if [[ $rnaStrand == "FIRST" ]] ; then
                echo "##running stranded STAR case"
                sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ALL,SAMTOOLSPATH=$samtoolsPath,BAM=$alignedBam,SAM=$alignedSam,GTF=$gtf,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_strandedHtSeqForStar.sh
                if [ $? -eq 0 ] ; then
                    touch $alignedSam.htSeqInQueue
                else
                    ((qsubFails++))
                fi
                sleep 2
            elif  [[ $rnaStrand == "SECOND"  ]] ; then
                                echo "##running stranded STAR case"
                                sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ALL,SAMTOOLSPATH=$samtoolsPath,BAM=$alignedBam,SAM=$alignedSam,GTF=$gtf,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_revStrandedHtSeqForStar.sh
                                if [ $? -eq 0 ] ; then
                                        touch $alignedSam.htSeqInQueue
                                else
                                        ((qsubFails++))
                                fi
                                sleep 2

            else
                echo "##running unstranded STAR case"
                sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ALL,SAMTOOLSPATH=$samtoolsPath,BAM=$alignedBam,SAM=$alignedSam,GTF=$gtf,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_htSeqForStar.sh
                if [ $? -eq 0 ] ; then
                    touch $alignedSam.htSeqInQueue
                else
                    ((qsubFails++))
                fi
                sleep 2
            fi
        ;;
    anotherRNAaligner) echo "### Anoter RNA aligner"
        sleep 2
        ;;
    *) echo "### I should not be here"
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
