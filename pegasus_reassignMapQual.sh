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

thisStep="pegasus_nextJob_reassignMapQual.sh"
#nxtStep1="pegasus_nextJob_snpEff.txt"

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

rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
snps=`grep @@"$recipe"@@ $constants | grep @@SNPS= | cut -d= -f2`
chrList=`grep @@"$recipe"@@ $constants | grep @@CHRLIST= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
echo "### reference: $ref"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0

for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
    echo "sample is $sampleLine"
    kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
    samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
    assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
    libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
    echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID"
    if [ "$assayID" == "RNA" ] ; then
        echo "### Assay ID is $assayID."
        rnaDir="$runDir/$kitName/$samName/$samName.starDir"
                rnaBam="$rnaDir/$samName.proj.Aligned.out.sorted.md.splitN.bam"
                rnaPass="$rnaDir/$samName.proj.Aligned.out.sorted.md.bam.splitNCigarPass"
        outBam="$rnaDir/$samName.proj.Aligned.out.sorted.md.splitN.reassignMap.bam"

        echo "### My bam is $rnaBam"
        if [ ! -e $rnaPass ] ; then
            echo "### Looks like splitNCigar is not done yet. $rnaPass doesnt exist yet"
            ((qsubFails++))
            continue
        fi
        if [ ! -e $rnaBam ] ; then
            echo "### Weird. Bam itself is missing: $rnaBam"
            ((qsubFails++))
            continue
        else
            if [[ -e $rnaBam.reassignMapPass || -e $rnaBam.reassignMapInQueue || -e $rnaBam.reassignMapFail ]] ; then
                            echo "### This bam file already passed, failed, or inQueue"
                                continue
                        fi
            echo "Starting reassignMappingQuality"
            sbatch --output $runDir/oeFiles/%x-slurm-%j.out --export GATKPATH=$gatkPath,OUTBAM=$outBam,RNABAM=$rnaBam,REF=$ref,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_reassignMapQual.sh
            if [ $? -eq 0 ] ; then
                touch $rnaBam.reassignMapInQueue
            else
                ((qsubFails++))
            fi

            sleep 2
        fi
    else
        echo "### Assay ID is $assayID. Must not be RNA reassignMappingQuality will only operate on RNA"
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
