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

thisStep="pegasus_nextJob_strelka.txt"
nxtStep1="pegasus_nextJob_snpEff.txt"
#nxtStep2="pegasus_nextJob_cuffDiff.txt"
#nxtStep3="pegasus_nextJob_htSeq.txt"
#nxtStep4="pegasus_nextJob_picardRNAMetrics.txt"
#nxtStep5="pegasus_nextJob_deSeq.txt"
#nxtStep6="pegasus_nextJob_checkProjectComplete.txt"
#nxtStep7="pegasus_nextJob_samtoolsStats.txt"

constants=~/jetstream/constants/constants.txt
constantsDir=~/jetstream/constants/
myName=`basename $0 | cut -d_ -f2`

declare -a chrGroups=(1:11:17:21 2:10:16:22 3:9:15:18:MT 4:7:14:Y 5:X:13:20 6:8:12:19)

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
strelka=`grep @@STRELKA= $constantsDir/$recipe | cut -d= -f2`

ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
starGTF=`grep "@@"$recipe"@@" $constants | grep @@STARGTF= | cut -d= -f2`
usegtf=`grep "@@"$recipe"@@" $constants | grep @@TOPHATGTF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
dnaAligner=`grep "@@"$recipe"@@" $constants | grep @@DNAALIGNER= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
picardPath=`grep @@"$recipe"@@ $constants | grep @@PICARDPATH= | cut -d= -f2`
starPath=`grep @@"$recipe"@@ $constants | grep @@STARPATH= | cut -d= -f2`
bwaPath=`grep "@@"$recipe"@@" $constants | grep @@BWAPATH= | cut -d= -f2`
seuratPath=`grep @@"$recipe"@@ $constants | grep @@SEURATPATH= | cut -d= -f2`
strelkaPath=`grep @@"$recipe"@@ $constants | grep @@STRELKAPATH= | cut -d= -f2`
jirRequested=`grep @@"$recipe"@@ $constants | grep @@JIRREALIGN= | cut -d= -f2`
faiFile="$ref.fai"
echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### somaticC: $somaticCaller"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
    if [ $strelka != "yes" ] ; then
                echo "strelka not requested for this recipe"
                continue
        fi

    echo "### DNA pair line is $dnaPairLine"
    sampleNames=`echo $dnaPairLine | cut -d= -f2`
    usableName=${sampleNames//,/-}

    pair1=`echo $sampleNames | cut -d, -f1`
    pair2=`echo $sampleNames | cut -d, -f2`

    pair1SampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$pair1"'"'`
    pair2SampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$pair2"'"'`
    pair1KitName=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f1`
    pair2KitName=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f1`
    pair1SamName=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f2`
    pair2SamName=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f2`
    pair1AssayNo=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f3`
    pair2AssayNo=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f3`

    case $pair1AssayNo$pair2AssayNo in
    GenomeGenome) echo "### Both normal and tumor are Genome"
            assay="Genome"
            ;;
    ExomeExome) echo "### Both normal and tumor are Exome"
            assay="Exome"
            ;;
    *)    echo "### Something else"
        echo "### Do not mix and match exomes and genomes!!!"
            assay="Mixed"
            ;;
    esac

    if [ "${jirRequested}" == "no" ] ; then
        normalBamFile=$runDir/$pair1KitName/$pair1/$pair1.proj.md.bam
        tumorBamFile=$runDir/$pair2KitName/$pair2/$pair2.proj.md.bam

        normalBaiFile=${normalBamFile/.bam/.bai}
        tumorBaiFile=${tumorBamFile/.bam/.bai}

        normalPass=$runDir/$pair1KitName/$pair1/$pair1.proj.bam.mdPass
        tumorPass=$runDir/$pair2KitName/$pair2/$pair2.proj.bam.mdPass
    else
        normalBamFile=$runDir/$pair1KitName/$pair1/$pair1.proj.md.jr.bam
        tumorBamFile=$runDir/$pair2KitName/$pair2/$pair2.proj.md.jr.bam

        normalBaiFile=${normalBamFile/.bam/.bai}
        tumorBaiFile=${tumorBamFile/.bam/.bai}

        normalPass=$runDir/$pair1KitName/$pair1/$pair1.proj.md.jr.bam.jointIRPass
        tumorPass=$runDir/$pair2KitName/$pair2/$pair2.proj.md.jr.bam.jointIRPass
    fi

    echo "### normal BAM: $normalBamFile"
    echo "### tumor  BAM: $tumorBamFile"

    if [[ ! -e $normalBamFile || ! -e $normalBaiFile || ! -e $normalPass ]] ; then
        echo "Normal bam, bai, or Pass does not exist"
        ((qsubFails++))
        continue
    fi
    if [[ ! -e $tumorBamFile || ! -e $tumorBaiFile || ! -e $tumorPass ]] ; then
        echo "Tumor bam, bai, or Pass does not exist"
        ((qsubFails++))
        continue
    fi
    strelkaDir="$runDir/strelka"
    if [ ! -d $strelkaDir ] ; then
        mkdir $strelkaDir
    fi
    mkdir -p $strelkaDir/$usableName

    wd="$runDir/strelka/$usableName"
    if [[ -e $wd.strelkaInQueue || -e $wd.strelkaPass || -e $wd.strelkaFail ]] ; then
        echo "### Strelka is already done, failed, or inqueue for $wd"
        continue
    fi
    echo "### Submitting to queue to run strelka on $wd"
    sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ASSAY=$assay,STRELKAPATH=$strelkaPath,WD=$wd,REF=$ref,NXT1=$nxtStep1,NORMAL=$normalBamFile,TUMOR=$tumorBamFile,D=$d $pegasusPbsHome/pegasus_strelka.sh
    if [ $? -eq 0 ] ; then
        touch $wd.strelkaInQueue
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
