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

thisStep="pegasus_nextJob_detectFusion.txt"
nxtStep1="pegasus_nextJob_tophatFusionPost.txt"

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
thfPath=`grep @@THFUSIONPATH= $constantsDir/$recipe | cut -d= -f2`
soapFusePath=`grep @@SOAPFUSEPATH= $constantsDir/$recipe | cut -d= -f2`


ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
thFusionRef=`grep "@@"$recipe"@@" $constants | grep @@THFUSIONREF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
fusionDetect=`grep "@@"$recipe"@@" $constants | grep @@FUSIONDETECT= | cut -d= -f2`
dnaAligner=`grep "@@"$recipe"@@" $constants | grep @@DNAALIGNER= | cut -d= -f2`
#indexbase=`grep @@"$recipe"@@ $constants | grep @@BOWTIE2INDEX= | cut -d= -f2`
tophat2Path=`grep @@"$recipe"@@ $constants | grep @@TOPHAT2PATH= | cut -d= -f2`
#bowtie2Path=`grep @@"$recipe"@@ $constants | grep @@BOWTIE2PATH= | cut -d= -f2`
bowtie1Path=`grep @@"$recipe"@@ $constants | grep @@BOWTIE1PATH= | cut -d= -f2`
bwaPath=`grep @@"$recipe"@@ $constants | grep @@BWAPATH= | cut -d= -f2`
picardPath=`grep @@"$recipe"@@ $constants | grep @@PICARDPATH= | cut -d= -f2`
trimFastqPath=`grep @@"$recipe"@@ $constants | grep @@TRIMFASTQPATH= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
faiFile="$ref.fai"
indexbase=${thFusionRef/.fa}
echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###
for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
    echo "sample is $sampleLine"
    kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
    samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
    assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
    libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
    rnaStrand=`grep "@@"$kitName"@@" $constants | cut -d= -f2`
    echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID, insSize: $insSize, rnaStrand: $rnaStrand"
    if [ "$assayID" != "RNA" ] ; then
        echo "### Assay ID is $assayID. Skipping."
        continue
    fi
    read1Name=$runDir/$kitName/$samName/$samName.proj.R1.fastq.gz
    read2Name=$runDir/$kitName/$samName/$samName.proj.R2.fastq.gz
    echo "read 1 name: $read1Name"
    echo "read 2 name: $read2Name"
    echo "### Fusion detector for recipe $recipe is $fusionDetect"
    case $fusionDetect in
    tophatfusion) echo "### Tophatfusion fusion detector"
        RG_CN=TGen
        RG_PL=ILLUMINA
        RG_ID_fromConfig="$libraID"
        #RG_ID_fromConfig=`echo ${mergeArray[$i]} | cut -d= -f2 | cut -d, -f1`
        FCID=`zcat $read1Name | head -n1 | cut -d: -f3`
        LANE=`zcat $read1Name | head -n1 | cut -d: -f4`
        INDEX=`zcat $read1Name | head -n1 | cut -d: -f10`
        RG_LB=`echo $RG_ID_fromConfig | cut -d_ -f3`
        RG_PU="${FCID}_${LANE}"
        RG_ID="${FCID}_${LANE}_${RG_LB}"
        rgTag="ID:${RG_ID}\tSM:$samName\tPL:${RG_PL}\tCN:${RG_CN}\tPU:${RG_PU}\tLB:${RG_LB}\tKS:${INDEX}"
        #rgTag="@RG\tID:$rgID\tSM:$samName\tPL:ILLUMINA\tLB:$libraID\tPI:$insertSize"
        thisPath=`dirname $runDir`
        cd $thisPath
        ownDir=${read1Name/.proj.R1.fastq.gz/.topHatFusionDir}
        if [[ ! -e $read1Name || ! -e $read2Name || ! -e $read1Name.mergeFastqPass || ! -e $read2Name.mergeFastqPass ]] ; then
            echo "one of the fastq files or read pass files dont exist"
            echo "read1Pass: $read1Name.mergeFastqPass"
            echo "read2Pass: $read2Name.mergeFastqPass"
            ((qsubFails++))
            continue
        fi
         if [[ -e $ownDir.thFusionPass || -e $ownDir.thFusionFail || -e $ownDir.thFusionInQueue ]] ; then
            echo "tophat fusion already done, failed or inQueue"
            continue
        fi
        if [ ! -d $ownDir ] ; then
            mkdir -p $ownDir
        fi

        #creating linked files to reads in tophatDir
        r1Name=`basename $read1Name`
        r2Name=`basename $read2Name`
        cd $ownDir
        ln -s $read1Name $r1Name
        read1Name=$ownDir/$r1Name
        ln -s $read2Name $r2Name
        read2Name=$ownDir/$r2Name
        cd -
        #done creating links. vars for reads changed.
        echo "### Read 1 name: $read1Name"
        echo "### Read 2 name: $read2Name"

        #manually set trim values to make 50mers??
        read1Length=`gunzip -c $read1Name | head -2 | tail -1 | wc -c`
        read2Length=`gunzip -c $read2Name | head -2 | tail -1 | wc -c`
        ((read1Length--));((read2Length--))
        pre1=0
        pre2=0
        pos1=`echo "$read1Length - 50" | bc`
        pos2=`echo "$read2Length - 50" | bc`
        echo "### Read1Length: ${read1Length}, Read2Length: $read2Length"
        echo "### Trim values are: $pre1, $pos1, $pre2, $pos2"

        echo "### Submitting $ownDir to queue for tophat fusion..."
        if [ $rnaStrand == "FIRST" ] ; then
                        echo "##running stranded tophatfusion case"
            sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ALL,FAI=$faiFile,PICARDPATH=$picardPath,REFPRETOPHAT=$ref,BWAPATH=$bwaPath,BOWTIE1PATH=$bowtie1Path,INDEXBASE=$indexBase,TOPHAT2PATH=$tophat2Path,THFUSIONPATH=$thfPath,SAMTOOLSPATH=$samtoolsPath,REF=$ref,FASTQ1=$read1Name,FASTQ2=$read2Name,DIR=$ownDir,INDEXBASE=$indexbase,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d,PRE1=$pre1,POS1=$pos1,PRE2=$pre2,POS2=$pos2 ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_firstStrandedtophatFusion.sh
            if [ $? -eq 0 ] ; then
                touch $ownDir.thFusionInQueue
            else
                ((qsubFails++))
            fi
            sleep 2

        elif [ $rnaStrand == "SECOND" ] ; then
                        echo "##running second  stranded tophatfusion case"
                        sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ALL,FAI=$faiFile,PICARDPATH=$picardPath,REFPRETOPHAT=$ref,BWAPATH=$bwaPath,BOWTIE1PATH=$bowtie1Path,INDEXBASE=$indexBase,TOPHAT2PATH=$tophat2Path,THFUSIONPATH=$thfPath,SAMTOOLSPATH=$samtoolsPath,REF=$ref,FASTQ1=$read1Name,FASTQ2=$read2Name,DIR=$ownDir,INDEXBASE=$indexbase,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d,PRE1=$pre1,POS1=$pos1,PRE2=$pre2,POS2=$pos2 ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_secondStrandedtophatFusion.sh
                        if [ $? -eq 0 ] ; then
                                touch $ownDir.thFusionInQueue
                        else
                                ((qsubFails++))
                        fi
                        sleep 2
        else
            echo "###running unstranded tophatfusion case"
            sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ALL,FAI=$faiFile,PICARDPATH=$picardPath,REFPRETOPHAT=$ref,BWAPATH=$bwaPath,BOWTIE1PATH=$bowtie1Path,INDEXBASE=$indexBase,TOPHAT2PATH=$tophat2Path,THFUSIONPATH=$thfPath,SAMTOOLSPATH=$samtoolsPath,REF=$ref,FASTQ1=$read1Name,FASTQ2=$read2Name,DIR=$ownDir,INDEXBASE=$indexbase,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d,PRE1=$pre1,POS1=$pos1,PRE2=$pre2,POS2=$pos2 ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_tophatFusion.sh
            if [ $? -eq 0 ] ; then
                touch $ownDir.thFusionInQueue
            else
                ((qsubFails++))
            fi
            sleep 2
        fi
        ;;
    soapFuse) echo "### Soap Fuse fusion detector"
        RG_CN=TGen
                RG_PL=ILLUMINA
                RG_ID_fromConfig="$libraID"
                #RG_ID_fromConfig=`echo ${mergeArray[$i]} | cut -d= -f2 | cut -d, -f1`
                FCID=`zcat $read1Name | head -n1 | cut -d: -f3`
                LANE=`zcat $read1Name | head -n1 | cut -d: -f4`
                INDEX=`zcat $read1Name | head -n1 | cut -d: -f10`
                RG_LB=`echo $RG_ID_fromConfig | cut -d_ -f3`
                RG_PU="${FCID}_${LANE}"
                RG_ID="${FCID}_${LANE}_${RG_LB}"
                rgTag="ID:${RG_ID}\tSM:$samName\tPL:${RG_PL}\tCN:${RG_CN}\tPU:${RG_PU}\tLB:${RG_LB}\tKS:${INDEX}"

        #rgTag="@RG\tID:$rgID\tSM:$samName\tPL:ILLUMINA\tLB:$libraID\tPI:$insertSize"
        thisPath=`dirname $runDir`
        cd $thisPath
        ownDir=${read1Name/.proj.R1.fastq.gz/.soapFuseDir}
        if [[ ! -e $read1Name || ! -e $read2Name || ! -e $read1Name.mergeFastqPass || ! -e $read2Name.mergeFastqPass ]] ; then
            echo "### one of the fastq files or read pass files dont exist"
            echo "### read1Pass: $read1Name.mergeFastqPass"
            echo "### read2Pass: $read2Name.mergeFastqPass"
            ((qsubFails++))
            continue
        fi
         if [[ -e $ownDir.soapFusePass || -e $ownDir.soapFuseFail || -e $ownDir.soapFuseInQueue ]] ; then
            echo "### soap fuse already done, failed or inQueue"
            continue
        fi
        if [ ! -d $ownDir ] ; then
            mkdir -p $ownDir
        fi
        #mkdir -p $ownDir/WHOLE_SEQ-DATA_DIR/$samName/$libraID/
        mkdir -p $ownDir/$samName/$libraID/

        read1Length=`gunzip -c $read1Name | head -2 | tail -1 | wc -c`
        ((read1Length--))

        spConfig=/home/tgenref/pecan/bin/SOAPfuse-v1.26/config/config.txt
        echo "$samName $libraID $samName $read1Length" > $ownDir/soapFuse.sampleList

        #creating linked files to reads in tophatDir
        r1Name=`basename $read1Name`
        r1Name=${r1Name/.proj.R1.fastq.gz/_1.fastq.gz}
        r2Name=`basename $read2Name`
        r2Name=${r2Name/.proj.R2.fastq.gz/_2.fastq.gz}
        #cd $ownDir/WHOLE_SEQ-DATA_DIR/$samName/$libraID/
        cd $ownDir/$samName/$libraID/
        ln -s $read1Name $r1Name
        read1Name=$ownDir/$r1Name
        ln -s $read2Name $r2Name
        read2Name=$ownDir/$r2Name
        cd -
        #done creating links. vars for reads changed.
        echo "### Read 1 name: $read1Name"
        echo "### Read 2 name: $read2Name"

        echo "### Submitting $ownDir to queue for soap fuse..."
        sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export ALL,SAMPLE=$samName,SLFILE=$ownDir/soapFuse.sampleList,SPCONFIG=$spConfig,SOAPFUSEPATH=$soapFusePath,SAMTOOLSPATH=$samtoolsPath,REF=$ref,FASTQ1=$read1Name,FASTQ2=$read2Name,DIR=$ownDir,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pecan_soapFuse.sh
        if [ $? -eq 0 ] ; then
            touch $ownDir.soapFuseInQueue
        else
            ((qsubFails++))
        fi
        sleep 2
        ;;
    anotherFusionDetector) echo "example fusion detector"
        sleep 2
        ;;
    none) echo "No fusion detection requested for this recipe."
        sleep 2
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
