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

thisStep="pegasus_nextJob_trn.txt"
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


ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
gtf=`grep @@"$recipe"@@ $constants | grep @@GTF= | cut -d= -f2`
samtoolsPath=`grep @@"$recipe"@@ $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
trnPath=`grep @@"$recipe"@@ $constants | grep @@TRNPATH= | cut -d= -f2`
indels=`grep @@"$recipe"@@ $constants | grep @@INDELS= | cut -d= -f2`
trn=`grep @@"$recipe"@@ $constants | grep @@TRNREQ= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
        if [[ $trn == "no" ]] ; then
                echo "### TRN is not requested for this recipe"
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
            #continue
            ;;
    *)    echo "### Something else"
        echo "### Do not mix and match exomes and genomes!!!"
            assay="Mixed"
            continue
            ;;
    esac

    normalBamFile=$runDir/$pair1KitName/$pair1/$pair1.proj.md.jr.bam
    tumorBamFile=$runDir/$pair2KitName/$pair2/$pair2.proj.md.jr.bam

    echo "normal BAM: $normalBamFile"
    echo "tumor  BAM: $tumorBamFile"

    normalBaiFile=${normalBamFile/.bam/.bai}
    tumorBaiFile=${tumorBamFile/.bam/.bai}

    if [[ ! -e $normalBamFile || ! -e $normalBaiFile || ! -e $normalBamFile.jointIRPass ]] ; then
        echo "Normal bam, bai, or jointIRPass does not exist"
        ((qsubFails++))
        continue
    fi
    if [[ ! -e $tumorBamFile || ! -e $tumorBaiFile || ! -e $tumorBamFile.jointIRPass ]] ; then
        echo "Tumor bam, bai, or jointIRPass does not exist"
        ((qsubFails++))
        continue
    fi

    bedFileGrep=$pair1KitName"_CNABED"
    bedFile=`grep "@@"$recipe"@@" $constants | grep @@"$bedFileGrep"= | cut -d= -f2`
    echo "### BED FILE= $bedFile"


    trnDir="$runDir/trn"
    if [ ! -d $trnDir ] ; then
        mkdir $trnDir
    fi
    mkdir -p $trnDir/$usableName
    trackName="$runDir/trn/$usableName/$usableName"
    trnDatFile="$runDir/trn/$usableName/$usableName"
    if [[ -e $trnDatFile.trnInQueue || -e $trnDatFile.trnFail || -e $trnDatFile.trnPass ]] ; then
        echo "Translocations are already passed, failed, or in queue for $trnDatFile"
        continue
    fi
    echo "### Submitting to queue with $normalBamFile"
    sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --mem-per-cpu 6000 --export ALL,ASSAY=$assay,BEDFILE=$bedFile,GTF=$gtf,NORMAL=$normalBamFile,TUMOR=$tumorBamFile,OUTFILE=$trnDatFile,TRNPATH=$trnPath,SAMPATH=$samtoolsPath,RUNDIR=$runDir,RECIPE=$recipe,NXT1=$nxtStep1,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_trn.sh
    if [ $? -eq 0 ] ; then
        touch $trnDatFile.trnInQueue
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
