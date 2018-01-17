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

thisStep="pegasus_nextJob_deNovoGear.txt"
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
pedFile=$runDir/$projName.ped
if [ ! -e $pedFile ] ; then
        echo "### Ped file not found at $pedFile!!!"
        echo "### Exiting!!!"
    echo "### Must not be a family study"
        echo "### I should remove $thisStep from $runDir."
        rm -f $runDir/$thisStep

        exit
else
        echo "### Ped file found."
fi
recipe=`cat $configFile | grep "^RECIPE=" | cut -d= -f2 | head -1 | tr -d [:space:]`
debit=`cat $configFile | grep "^DEBIT=" | cut -d= -f2 | head -1 | tr -d [:space:]`

nCores=`grep @@${myName}_CORES= $constantsDir/$recipe | cut -d= -f2`
ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
hapmap=`grep @@"$recipe"@@ $constants | grep @@HAPMAP= | cut -d= -f2`
omni=`grep @@"$recipe"@@ $constants | grep @@OMNI= | cut -d= -f2`
deNovoGearPath=`grep @@"$recipe"@@ $constants | grep @@DENOVOGEARPATH= | cut -d= -f2`
bcfToolsPath=`grep @@"$recipe"@@ $constants | grep @@BCFTOOLSPATH= | cut -d= -f2`

snpeffdb=`grep @@"$recipe"@@ $constants | grep @@SNPEFFDB= | cut -d= -f2`
dbsnp=`grep @@"$recipe"@@ $constants | grep @@SNPS= | cut -d= -f2`
snpeffPath=`grep @@"$recipe"@@ $constants | grep @@SNPEFFPATH= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0

##only for projects with PED file
## This loop if for finding multi variant called vcf files from GATK
for dnaPairLine in `cat $configFile | grep '^DNAFAMI='`
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
    sampleList=""
        missingSampleCount=0
    echo "$sampleNames"
    echo "${sampleNames//,/} "
        trackName="$runDir/denovoGear/$usableName"
        hcTrackName="$runDir/hc/$usableName/$usableName"
        outTrackName="$runDir/denovoGear/$usableName/$usableName"
        outTrack="$runDir/denovoGear"
        outVcf="$runDir/denovoGear/$usableName/$usableName"
        if [ ! -d $outTrack ] ; then
                mkdir $outTrack
        fi
        mkdir -p $outTrack/$usableName
    bamText="$outTrack/${usableName}_bams.txt"
    for eachSample in ${sampleNames//,/ }
    do
        ((sampleCount++))
        echo "eachsample: $eachSample"
        sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
        kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
        samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
        assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
        echo "sampleLine: $sampleLine, kitName: $kitName, samName: $samName"
        eachSampleBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
                eachSamplePass=$runDir/$kitName/$samName/$samName.proj.md.jr.bam.jointIRPass
                if [[ ! -e $eachSampleBam || ! -e $eachSamplePass ]] ; then
                        echo "### Can't find the jr bam or the jr pass"
                        echo "### BAM: $eachSampleBam"
                        echo "### PAS: $eachSamplePass"
                        ((missingSampleCount++))
                else
                        #sampleList="$eachSampleBam $sampleList"
                    echo "$eachSampleBam" >> $bamText
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


    if [ ! -e $hcTrackName.HC_All.vcf.snpEffPass ] ; then
        echo "### HC snpEff pass is not done yet. Missing: $hcTrackName.HC_All.vcf.snpEffPass"
        ((qsubFails++))
        continue
    fi

    if [[ -e $outTrackName.deNovoGearInQueue || -e $outTrackName.deNovoGearPass || -e $outTrackName.deNovoGearFail ]] ; then
        echo "### deNovoGear is already done, failed or inQueue"
        continue
    fi
    echo "### Submitting $trackName.HC_All.vcf to queue for deNovoGear..."
    sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out --export ALL,DENOVOPATH=$deNovoGearPath,TRACKNAME=$trackname,SAMTOOLSPATH=$samtoolsPath,BAMLIST="$sampleList",BAMFILE=$bamText,BCFTOOLSPATH=$bcfToolsPath,GATKPATH=$gatkPath,PED=$pedFile,OUTVCF=$outVcf,REF=$ref,OUTTRACKNAME=$outTrackName,SNPEFFPATH=$snpeffPath,VCF=${hcTrackName}.HC_All.vcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_deNovoGear.sh
    if [ $? -eq 0 ] ; then
        touch $outTrackName.deNovoGearInQueue
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
