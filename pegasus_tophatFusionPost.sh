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

thisStep="pegasus_nextJob_tophatFusionPost.txt"
nxtStep1="pegasus_nextJob_checkProjectComplete.txt"

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
tophat2Path=`grep @@"$recipe"@@ $constants | grep @@TOPHAT2PATH= | cut -d= -f2`
thFusion2vcfPath=`grep @@"$recipe"@@ $constants | grep @@THFUSION2VCFPATH= | cut -d= -f2`
thFusionRef=`grep "@@"$recipe"@@" $constants | grep @@THFUSIONREF= | cut -d= -f2`
gtfmask=`grep @@"$recipe"@@ $constants | grep @@GTFMASK= | cut -d= -f2`
gtf=`grep @@"$recipe"@@ $constants | grep @@GTF= | cut -d= -f2`
thfDir=`grep @@"$recipe"@@ $constants | grep @@THFUSIONPOST= | cut -d= -f2`

indexbase=${thFusionRef/.fa}

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
    topHatFDir="$runDir/$kitName/$samName/$samName.topHatFusionDir"
    if [ ! -e $topHatFDir.thFusionPass ] ; then
        echo "tophat fusion not complete yet"
        continue
    fi
    echo $topHatFDir
    if [[ -e $topHatFDir.thFPostInQueue || -e $topHatFDir.thFPostPass || -e $topHatFDir.thFPostFail ]] ; then
            echo "$topHatFDir already in queue, passed or failed"
            continue
    fi
    #making that directory for tophat.. madness
    thatDir=`basename $topHatFDir`
    thatDir=${thatDir/.topHatFusionDir}
    thatName=$thatDir
    thatDir="tophat_$thatDir"
    thatDir=$topHatFDir/$thatDir
    mkdir $thatDir
    cd $thatDir
    ln -s ../$thatName.proj.accepted_hits.bam accepted_hits.bam
    ln -s ../$thatName.proj.fusions.out fusions.out
    cd -
    echo "that directory is $thatDir"
    #end making that directory tophat
    #linking stuff
    if [ ! -e $topHatFDir/refGene.txt ] ; then
        ln -s $thfDir/refGene.txt $topHatFDir/refGene.txt
    fi
    if [ ! -e $topHatFDir/ensGene.txt ] ; then
        ln -s $thfDir/ensGene.txt $topHatFDir/ensGene.txt
    fi
    if [ ! -e $topHatFDir/mcl ] ; then
        ln -s $thfDir/mcl $topHatFDir/mcl
    fi
    if [ ! -d $topHatFDir/blast ] ; then
        ln -s $thfDir/blast $topHatFDir/blast
    fi
    #end linking stuff

    sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ALL,THFUSION2VCFPATH=$thFusion2vcfPath,TOPHAT2PATH=$tophat2Path,RUNDIR=$runDir,DIR=$topHatFDir,INDEXBASE=$indexbase,NXT1=$nxtStep1,REF=$thFusionRef,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_tophatFusionPost.sh
    if [ $? -eq 0 ] ; then
        touch $topHatFDir.thFPostInQueue
    else
        ((qsubFails++))
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
