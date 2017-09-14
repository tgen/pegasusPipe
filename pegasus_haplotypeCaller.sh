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

thisStep="pegasus_nextJob_haplotypeCaller.txt"
nxtStep1="pegasus_nextJob_snpEff.txt"
nxtStep2="pegasus_nextJob_phaseBT.txt"
nxtStep3="pegasus_nextJob_sexRelCheck.txt"

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

#ref=`cat $configFile | grep "^REF=" | cut -d= -f2 | head -1 | tr -d [:space:]`
#if [[ -z "$ref"  ]] ; then
ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
#fi

rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
snps=`grep @@"$recipe"@@ $constants | grep @@SNPS= | cut -d= -f2`
chrList=`grep @@"$recipe"@@ $constants | grep @@CHRLIST= | cut -d= -f2`
jirRequested=`grep @@"$recipe"@@ $constants | grep @@JIRREALIGN= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
echo "### reference: $ref"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0

#DNAPAIR samples have multi sample gatkHC calls and seurat calls
#DNAFAMI samples have multi sample gatkHC calls
#samples in neither group have single sample gatkHC calls 

for dnaPairLine in `cat $configFile | grep '^DNAPAIR=\|^DNAFAMI='`
do
    echo "### DNA pair line is $dnaPairLine"
    sampleNames=`echo $dnaPairLine | cut -d= -f2 | cut -d\; -f1`
    shorterName=`echo $dnaPairLine | cut -d= -f2 | cut -d\; -f2`
    if [[ $shorterName == *,* ]] ; then
        shorterName=""
        echo "### No nick name detected for these pairs"
        usableName=${sampleNames//,/-}
    else
        echo "### Nick name for these pairs are: $shorterName"
        usableName="$shorterName"
    fi
    sampleCount=0
    missingSampleCount=0
    sampleList=""
    for eachSample in ${sampleNames//,/ }
    do
        ((sampleCount++))
        #echo "eachsample: $eachSample"
        sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
        kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
        samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
        if [ "${jirRequested}" == "no" ] ; then

                        eachSampleBam=$runDir/$kitName/$samName/$samName.proj.md.bam
                        eachSamplePass=$runDir/$kitName/$samName/$samName.proj.bam.mdPass
                else
            eachSampleBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
            eachSamplePass=$runDir/$kitName/$samName/$samName.proj.md.jr.bam.jointIRPass
        fi
        echo "sampleLine: $sampleLine, kitName: $kitName, samName: $samName"
        if [[ ! -e $eachSampleBam || ! -e $eachSamplePass ]] ; then
            echo "### Can't find the jr bam or the jr pass"
            echo "### BAM: $eachSampleBam"
            echo "### PAS: $eachSamplePass"
            ((missingSampleCount++))
        else
            sampleList="-I $eachSampleBam $sampleList"
        fi
    done
    if [ $missingSampleCount -eq 0 ] ; then
        echo "### All sample bams are found"
        sampleList="${sampleList%?}"
        #echo "### Sample list: $sampleList"
    else
        echo "### There were missing things"
        echo "### Samples missing: $missingSampleCount/$sampleCount"
        ((qsubFails++))
        continue
    fi
    hcDir="$runDir/hc"
    if [ ! -d $hcDir ] ; then
        mkdir $hcDir
    fi
    mkdir -p $hcDir/$usableName
    workDir="$hcDir/$usableName"
    trackName="$runDir/hc/$usableName/$usableName"
    STEP=0
    #STEP_COUNT=24
    STEP_COUNT=`ls $chrList/*list | wc -l`
    echo "### Submitting to queue to run hc on $wd"
    while [ ${STEP} -lt $STEP_COUNT ]
    do
        (( STEP++ ))
        if [[ -e ${trackName}_Step${STEP}.hcInQueue || -e ${trackName}_Step${STEP}.hcPass || -e ${trackName}_Step${STEP}.hcFail ]] ; then
            echo "### Haplotype caller is already done, failed, or inqueue for ${trackName}"
            continue
        fi

        echo Starting Haplotype caller Step${STEP}
        sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export GATKPATH=$gatkPath,STEPCOUNT=$STEP_COUNT,TRK=$trackName,KNOWN=$snps,BAMLIST="$sampleList",TRK=$trackName,CHRLIST=$chrList,REF=$ref,STEP=${STEP},NXT1=$nxtStep1,NXT2=$nxtStep2,NXT3=$nxtStep3,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_haplotypeCaller.sh
        if [ $? -eq 0 ] ; then
            touch ${trackName}_Step${STEP}.hcInQueue
        else
            ((qsubFails++))
        fi

        sleep 2
    done
done
for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
    echo "sample is $sampleLine"
    kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
    samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
    assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
    libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
    echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID"
    if [ "$assayID" != "RNA" ] ; then
        echo "### Assay ID is $assayID. Must be genome or exome."
        pcDir=$runDir/$kitName/$samName
        inBam=$runDir/$kitName/$samName/$samName.proj.bam
        mdBam=$runDir/$kitName/$samName/$samName.proj.md.bam
        jrBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
        jrRequested=`cat $configFile | grep '^DNAPAIR=\|^DNAFAMI=' | grep $samName | wc -l`
        if [ $jrRequested -gt 0 ] ; then
            echo "### Haplotype caller will not operate on a single bam because Joint IR is requested for $samName"
        else
            echo "### Haplotype caller will run on single bam because Joint IR is NOT requested for $samName"
            if [[ ! -e $inBam.mdPass || ! -e $mdBam ]] ; then
                echo "### Either mdPass or the bam itself is missing for $mdBam"
                ((qsubFails++))
            else
                STEP=0
                STEP_COUNT=`ls $chrList/*list | wc -l`
                echo "### Submitting to queue to run hc on $wd"
                while [ ${STEP} -lt $STEP_COUNT ]
                do
                    (( STEP++ ))
                    if [[ -e ${mdBam}_Step${STEP}.hcInQueue || -e ${mdBam}_Step${STEP}.hcPass || -e ${mdBam}_Step${STEP}.hcFail ]] ; then
                        echo "### Haplotype caller is already done, failed, or inqueue for ${mdBam}"
                        continue
                    fi

                    echo Starting Haplotype caller for Step${STEP}
                    sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export GATKPATH=$gatkPath,STEPCOUNT=$STEP_COUNT,TRK=$mdBam,KNOWN=$snps,BAMLIST=$mdBam,CHRLIST=$chrList,REF=$ref,STEP=${STEP},NXT1=$nxtStep1,NXT2=$nxtStep2,NXT3=$nxtStep3,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_haplotypeCallerSingle.sh
                    if [ $? -eq 0 ] ; then
                        touch ${mdBam}_Step${STEP}.hcInQueue
                    else
                        ((qsubFails++))
                    fi

                    sleep 2
                done
            fi
        fi
    else
        echo "### Assay ID is $assayID. Must be RNA."
        #code for calling AS metrics on tophat bams here
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
