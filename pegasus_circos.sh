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

thisStep="pegasus_nextJob_circos.txt"
nxtStep1="pegasus_nextJob_postCircos.txt"

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
ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
samtoolsPath=`grep @@"$recipe"@@ $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
#circosPath=`grep @@"$recipe"@@ $constants | grep @@TRNPATH= | cut -d= -f2`
indels=`grep @@"$recipe"@@ $constants | grep @@INDELS= | cut -d= -f2`
cosmicVcf=`grep "@@"$recipe"@@" $constants | grep @@COSMIC_VCF= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
###
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
    echo "### DNA pair line is $dnaPairLine"
    ## 1. runs on whole genome based on DNAPAIR
    ## 2. requires seurat output from an exome pair (that means exome pair names required)
    ## 3. requires CGH output, namely *cgh.tsv file from this pair
    ## 4. requires translocations output from this pair

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
    assayID=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f3`

    normalBamFile=$runDir/$pair1KitName/$pair1/$pair1.proj.md.jr.bam
    tumorBamFile=$runDir/$pair2KitName/$pair2/$pair2.proj.md.jr.bam
    normalBaiFile=${normalBamFile/.bam/.bai}
    tumorBaiFile=${tumorBamFile/.bam/.bai}

    echo "first checking for seurat snpeff vcf"
    seuratTrackName="$runDir/seurat/$usableName/$usableName"
    if [ ! -e $seuratTrackName.seurat.vcf.snpEffPass ] ; then
        echo "### Seurat snpEffPass doesnt exist yet: $seuratTrackName.seurat.vcf.snpEffPass"
        ((qsubFails++))
        exit
    fi

    if [ $assayID == "Exome"  ] ; then
        cnvsTSV="$runDir/cna/${usableName}_exo/${usableName}_exo.cna.tsv"
        cnvsPass="$runDir/cna/${usableName}_exo/${usableName}_exo.cnaExomePass"
    elif [ $assayID == "Genome" ] ; then
        cnvsTSV="$runDir/cna/${usableName}_filt/${usableName}_filt.cna.tsv"
        cnvsPass="$runDir/cna/${usableName}_filt/${usableName}_filt.cnaGenFiltPass"
    else
        echo "I should not be here, assay ID was not Exome or Genome"
    fi

    seuratVcf="$seuratTrackName.seurat.snpEff.txt"
    trnVcf="$runDir/trn/$usableName/$usableName.trn.vcf"

    if [[ ! -e $seuratTrackName.seurat.vcf.snpEffPass || ! -e $cnvsPass || ! -e $runDir/trn/$usableName/$usableName.trnPass  ]] ; then
        echo "One of these does not exist:"
        echo "$seuratTrackName.seurat.vcf.snpEffPass"
        echo "$cnvsPass"
        echo "$runDir/trn/$usableName/$usableName.trnPass"
        exit
    fi

    echo "seurat vcf: $seuratVcf"
    echo "cnvsTSV: $cnvsTSV"
    echo "trnVcf: $trnVcf"

    circosDir="$runDir/circos"
    if [ ! -d $circosDir ] ; then
        mkdir $circosDir
    fi

    mkdir -p $circosDir/$usableName
    #trackName="$runDir/circos/$usableName/$usableName"
    outDir="$circosDir/$usableName"
    circosSamples="$runDir/circos/$usableName/$usableName"
    if [[ -e $circosSamples.circosInQueue || -e $circosSamples.circosFail || -e $circosSamples.circosPass ]] ; then
        echo "Circos are already passed, failed, or in queue for $circosSamples"
        continue
    fi

    conf=$outDir/template_circos.conf
    cp -r /home/tizatt/circosTemplateFolder/* $outDir/
    cat $outDir/template_circos.part1.conf > $outDir/template_circos.conf
    echo "#########pipeline insertion start************" >> $outDir/template_circos.conf
    echo "dir = $outDir" >> $outDir/template_circos.conf
    echo "file = $projName.circos1.png" >> $outDir/template_circos.conf
    echo "#########pipeline insertion end**************" >> $outDir/template_circos.conf
    cat $outDir/template_circos.part2.conf >> $outDir/template_circos.conf

    echo "### Submitting to queue with $normalBamFile"
    sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out --export ALL,OUTFILE=$circosSamples,CONF=$conf,OUTDIR=$outDir,SEURATVCF=$seuratVcf,COSMIC=$cosmicVcf,TRNVCF=$trnVcf,CNVTSV=$cnvsTSV,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_circos.sh
    if [ $? -eq 0 ] ; then
        touch $circosSamples.circosInQueue
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
