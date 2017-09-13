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

thisStep="pegasus_nextJob_cna15.txt"
nxtStep1="pegasus_nextJob_postCna15.txt"

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
matchedNormal=`cat $configFile | grep "^MATCHEDNORMAL=" | cut -d= -f2 | head -1 | tr -d [:space:]`
nCores=`grep @@${myName}_CORES= $constantsDir/$recipe | cut -d= -f2`

ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
cnaPath=`grep @@"$recipe"@@ $constants | grep @@CNA15PATH= | cut -d= -f2`
ccdsPath=`grep @@"$recipe"@@ $constants | grep @@CNA15CCDS= | cut -d= -f2`
seuratPath=`grep @@"$recipe"@@ $constants | grep @@SEURATPATH= | cut -d= -f2`
indels=`grep @@"$recipe"@@ $constants | grep @@INDELS= | cut -d= -f2`
cnaGenomeFilter=`grep @@"$recipe"@@ $constants | grep @@CNAGENOMEFILTER= | cut -d= -f2`
cnaGenomeUnFilter=`grep @@"$recipe"@@ $constants | grep @@CNAGENOMEUNFILTER= | cut -d= -f2`
gtf=`grep @@"$recipe"@@ $constants | grep @@GTF= | cut -d= -f2`
cna=`grep @@"$recipe"@@ $constants | grep @@CNA= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
# find the vcf from UG separately, just the UG.vcf ran on jr.bam
# NO
# find the merged VCF created by UG
###
for dnaPairLine in `cat $configFile | grep '^DNAPAIR=\|^DNAFAMI='`
do
    if [[ $cna != "yes" ]] ; then
        echo "### CNA is not requested for this recipe"
        continue
    fi
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

#    sampleNames=`echo $dnaPairLine | cut -d= -f2`
#    usableName=${sampleNames//,/-}

    famORpair=`echo $dnaPairLine | cut -d= -f1`
    if [ $famORpair == "DNAPAIR" ] ; then
        echo "RUNNING WITH DNAPAIR LINE"
        pair1=`echo $sampleNames | cut -d, -f1`
        pair2=`echo $sampleNames | cut -d, -f2`

        #usableName="$pair1-$pair2"
        outName=$usableName
        #ugDir="$runDir/ug/$usableName"
        hcDir="$runDir/hc/$usableName"
        #ugVcf=$usableName.UG.vcf
        #ugPas=$usableName.ugPass
        hcVcf=$usableName.HC_All.snpEff.vcf
        hcPas=$usableName.HC_All.vcf.snpEffPass
        vcfFile=$hcDir/$hcVcf
        vcfPass=$hcDir/$hcPas


        pair1SampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$pair1"'"'`
        pair2SampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$pair2"'"'`
        pair1KitName=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f1`
        pair2KitName=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f1`
        pair1SamName=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f2`
        pair2SamName=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f2`
        pair1AssayNo=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f3`
        pair2AssayNo=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f3`

        normalDatFile=$runDir/$pair1KitName/$pair1/$pair1.proj.md.jr.bam.clc.cln.dat
        tumorDatFile=$runDir/$pair2KitName/$pair2/$pair2.proj.md.jr.bam.clc.cln.dat
        normalDatPass=$runDir/$pair1KitName/$pair1/$pair1.proj.md.jr.bam.clonalCovPass
        tumorDatPass=$runDir/$pair2KitName/$pair2/$pair2.proj.md.jr.bam.clonalCovPass

        nHsMetricPass=$runDir/$pair1KitName/$pair1/$pair1.proj.md.jr.bam.picHSMetricsPass
        tHsMetricPass=$runDir/$pair2KitName/$pair2/$pair2.proj.md.jr.bam.picHSMetricsPass

        nHsMetric=$runDir/stats/$pair1.proj.md.jr.bam.picHSMetrics
        tHsMetric=$runDir/stats/$pair2.proj.md.jr.bam.picHSMetrics

        echo "### normal DAT: $normalDatFile, $pair1AssayNo"
        echo "### tumor  DAT: $tumorDatFile, $pair2AssayNo"

        bedFileGrep=$pair1KitName"_CNABED"
        bedFile=`grep "@@"$recipe"@@" $constants | grep @@"$bedFileGrep"= | cut -d= -f2`
        echo "matchedNormal = ${matchedNormal}"
        if [[ ${matchedNormal} == "No" ]] ; then
            echo "BED file will be swapped out since it has no matched normal"
            if [[ "$pair1KitName" == *S5U ]]  ; then
                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_V5_PlusUTR_hs37d5_PaddedTargets_Picard_DGV_1kg_100_interval.bed"
            elif [[ "$pair1KitName" == *S6X ]]  ; then
                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_V6_noUTR_hs37d5_PaddedTargets_Picard_DGV_1kg_100_interval.bed"
            elif [[ "$pair1KitName" == *S4U ]]  ; then
                                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_SureSelect_V4_plusUTR/Agilent_SureSelect_V4_plusUTR_hs37d5_Padded_DGV_1kg.cna.bed"
                elif [[ "$pair1KitName" == *S4U ]]  ; then
                                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_SureSelect_V6R2_StxProstate/Agilent_SureSelect_V6R2_StxProstate_hs37d5_Padded_DGV_1kg.cna.bed"
            elif [[ "$pair1KitName" == *ST2 ]]  ; then
                                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_SureSelect_V6R2_StrexomeV2/Agilent_SureSelect_V6R2_StrexomeV2_hs37d5_Padded_DGV_1kg.cna.bed"
            elif [[ "$pair1KitName" == *S2U ]]  ; then
                                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_SureSelect_V2_NA/Agilent_SureSelect_V2_NA_hs37d5_Padded_DGV_1kg.cna.bed"
            elif [[ "$pair1KitName" == *CCC ]]  ; then
                                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_ClearSeq_Beta_ComprehensiveCancer/Agilent_ClearSeq_Beta_ComprehensiveCancer_hs37d5_Padded_DGV_1kg.cna.bed"
            elif [[ "$pair1KitName" == *CR2 ]]  ; then
                bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_SureSelect_CREv2_cliRes/Agilent_SureSelect_CREv2_cliRes_hs37d5_Padded_DGV_1kg.cna.bed"
            else
                echo "We don't have a bed filter for this kit for no matched normal at this time"
            fi
        fi
        echo "### BED FILE= $bedFile"

        if [[ ! -e $normalDatFile || ! -e $tumorDatFile || ! -e $normalDatPass || ! -e $tumorDatPass ]] ; then
            echo "### Normal, tumor dat file or pass does not exist"
            ((qsubFails++))
            continue
        fi
        if [[ ! -e $vcfFile || ! -e $vcfPass ]] ; then
            echo "### HC vcf file called on DNAPAIR not found or not passed yet: $vcfFile"
            ((qsubFails++))
            continue
        fi
        cnaDir="$runDir/cna"
        if [ ! -d $cnaDir ] ; then
            mkdir -p $cnaDir
        fi
        #mkdir -p $cnaDir/$usableName
        #trackName="$runDir/cna/$usableName/$usableName"
        case $pair1AssayNo$pair2AssayNo in
        GenomeGenome) echo "### Both normal and tumor are Genome"
                assay="Genome"
                if [ "$cnaGenomeFilter" = "yes" ] ; then
                    outName1=$outName"_filt"
                    mkdir -p $cnaDir/$outName1
                    trackName="$runDir/cna/$outName1/$outName1"
                    if [[ -e $trackName.cnaGenFiltInQueue || -e $trackName.cnaGenFiltPass || -e $trackName.cnaGenFiltFail ]] ; then
                        echo "CNA gen filt plot already in queue, passed, or failed"
                    else
                        echo "### Submitting to queue with $trackName"
                        sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export TRACKNAME=$trackName,MYPATH=$cnaDir/$outName1,CCDSLIST=$ccdsPath,CNAPATH=$cnaPath,MERGEDVCF=$vcfFile,NORMALSAMPLE=$pair1SamName,TUMORSAMPLE=$pair2SamName,NORMALDAT=$normalDatFile,TUMORDAT=$tumorDatFile,OFILE=$outName1,ASSAY=$assay,GTF=$gtf,NXT1=$nxtStep1,REF=$ref,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_cnaGenFilt15.sh
                        if [ $? -eq 0 ] ; then
                            touch $trackName.cnaGenFiltInQueue
                            sleep 1
                        else
                            ((qsubFails++))
                        fi
                        sleep 2
                    fi
                fi
                if [ "$cnaGenomeUnFilter" = "yes" ] ; then
                    outName2=$outName"_unfi"
                    mkdir -p $cnaDir/$outName2
                    trackName="$runDir/cna/$outName2/$outName2"
                    if [[ -e $trackName.cnaGenUnfiInQueue || -e $trackName.cnaGenUnfiPass || -e $trackName.cnaGenUnfiFail ]] ; then
                        echo "CNA gen filt plot already in queue, passed, or failed"
                    else
                        echo "### Submitting to queue with $trackName"
                        sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export TRACKNAME=$trackName,MYPATH=$cnaDir/$outName2,CCDSLIST=$ccdsPath,CNAPATH=$cnaPath,MERGEDVCF=$vcfFile,NORMALSAMPLE=$pair1SamName,TUMORSAMPLE=$pair2SamName,NORMALDAT=$normalDatFile,TUMORDAT=$tumorDatFile,OFILE=$outName2,ASSAY=$assay,GTF=$gtf,NXT1=$nxtStep1,REF=$ref,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_cnaGenUnfi15.sh
                        if [ $? -eq 0 ] ; then
                            touch $trackName.cnaGenUnfiInQueue
                            sleep 1
                        else
                            ((qsubFails++))
                        fi
                        sleep 2
                    fi
                fi

                ;;
        ExomeExome) echo "### Both normal and tumor are Exome"
                assay="Exome"
                outName3=$outName"_exo"

                if [[ ! -e $nHsMetricPass || ! -e $tHsMetricPass || ! -e $tHsMetric || ! -e $nHsMetric ]] ; then
                    echo "### Normal, tumor Hs Metrics Pass or the hs metric file itself does not exist"
                    echo "$nHsMetricPass"
                    echo "$tHsMetricPass"
                    ((qsubFails++))
                    continue
                fi
                nHetDepth=`cat $nHsMetric | grep -A 1 BAIT_SET | tail -1 | cut -f22`
                tHetDepth=`cat $tHsMetric | grep -A 1 BAIT_SET | tail -1 | cut -f22`

                mkdir -p $cnaDir/$outName3
                trackName="$runDir/cna/$outName3/$outName3"
                if [[ -e $trackName.cnaExomeInQueue || -e $trackName.cnaExomePass || -e $trackName.cnaExomeFail ]] ; then
                    echo "CNA gen filt plot already in queue, passed, or failed"
                else
                    echo "### Submitting to queue with $trackName"
                    sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export NHETDEPTH=$nHetDepth,THETDEPTH=$tHetDepth,TRACKNAME=$trackName,CCDSLIST=$ccdsPath,MYPATH=$cnaDir/$outName3,CNAPATH=$cnaPath,MERGEDVCF=$vcfFile,NORMALSAMPLE=$pair1SamName,TUMORSAMPLE=$pair2SamName,CNAEXOMETARGET=$bedFile,NORMALDAT=$normalDatFile,TUMORDAT=$tumorDatFile,OFILE=$outName3,ASSAY=$assay,GTF=$gtf,NXT1=$nxtStep1,REF=$ref,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_cnaExo15.sh
                    if [ $? -eq 0 ] ; then
                        touch $trackName.cnaExomeInQueue
                        sleep 1
                    else
                        ((qsubFails++))
                    fi
                    sleep 2
                fi
                ;;
        *)    echo "### Something else"
            echo "### Do not mix and match exomes and genomes!!!"
                assay="Mixed"
                ;;
        esac
    elif [ $famORpair == "DNAFAMI" ] ; then
        pedFile="$runDir/$projName.ped"
        if [[ -z $pedFile ]] ; then
            echo "### FAMILIES NEED PED FILE TO RUN COPY NUMBER"
            echo "### EXITING"
            exit
        fi
        echo "RUNNING WITH DNAFAMI LINE"
        echo "$dnaPairLine"
        #goes through ped file line by line
        while read sampleLine;
        do
            pair1=`echo $sampleLine | cut -d ' ' -f2`
            pair2=`echo $sampleLine | cut -d ' ' -f3`
            pair3=`echo $sampleLine | cut -d ' ' -f4`
            echo "sampleLine= $sampleLine"
            echo "pair1= $pair1"
            echo "pair2= $pair2"
            echo "pair3= $pair3"

            if [[ "$pair2" == 0 && "$pair3" == 0 ]] ; then
                echo "Sample: $pair1 doesn't have a mom or dad listed, will continue"
                continue
            fi
            #make sure each person is in the current DNAPAIR line and thus hc vcf
            pair1Check=`echo $dnaPairLine | grep ${pair1} | head -1`
            pair2Check=`echo $dnaPairLine | grep ${pair2} | head -1`
            pair3Check=`echo $dnaPairLine | grep ${pair3} | head -1`
            echo "pair1Check: $pair1Check"
            echo "pair2Check: $pair2Check"
            echo "pair3Check: $pair3Check"

            if [[ -z $pair1Check || -z $pair2Check || -z $pair3Check ]] ; then
                echo "One of these samples wasn't in the DNAFAMI line  pair1: $pair1 pair2: $pair2 pair3: $pair3 "
                echo "Will continue"
                continue
            fi
            #usableName="$pair1-$pair2"
            #ugDir="$runDir/ug/$usableName"
            hcDir="$runDir/hc/$usableName"
            #ugVcf=$usableName.UG.vcf
            #ugPas=$usableName.ugPass
            hcVcf=$usableName.HC_All.snpEff.vcf
            hcPas=$usableName.HC_All.vcf.snpEffPass
            vcfFile=$hcDir/$hcVcf
            vcfPass=$hcDir/$hcPas


            pair1SampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$pair1"'"'`
            pair2SampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$pair2"'"'`
            pair3SampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$pair3"'"'`
            pair1KitName=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f1`
            pair2KitName=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f1`
            pair3KitName=`echo $pair3SampleLine | cut -d= -f2 | cut -d, -f1`
            pair1SamName=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f2`
            pair2SamName=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f2`
            pair3SamName=`echo $pair3SampleLine | cut -d= -f2 | cut -d, -f2`
            pair1AssayNo=`echo $pair1SampleLine | cut -d= -f2 | cut -d, -f3`
            pair2AssayNo=`echo $pair2SampleLine | cut -d= -f2 | cut -d, -f3`
            pair3AssayNo=`echo $pair3SampleLine | cut -d= -f2 | cut -d, -f3`

            normalDatFile1=$runDir/$pair2KitName/$pair2/$pair2.proj.md.jr.bam.clc.cln.dat
            normalDatFile2=$runDir/$pair3KitName/$pair3/$pair3.proj.md.jr.bam.clc.cln.dat
            tumorDatFile=$runDir/$pair1KitName/$pair1/$pair1.proj.md.jr.bam.clc.cln.dat
            normalDatPass1=$runDir/$pair2KitName/$pair2/$pair2.proj.md.jr.bam.clonalCovPass
            normalDatPass2=$runDir/$pair3KitName/$pair3/$pair3.proj.md.jr.bam.clonalCovPass
            tumorDatPass=$runDir/$pair1KitName/$pair1/$pair1.proj.md.jr.bam.clonalCovPass

            nHsMetricPass1=$runDir/$pair2KitName/$pair2/$pair2.proj.md.jr.bam.picHSMetricsPass
            nHsMetricPass2=$runDir/$pair3KitName/$pair3/$pair3.proj.md.jr.bam.picHSMetricsPass
            tHsMetricPass=$runDir/$pair1KitName/$pair1/$pair1.proj.md.jr.bam.picHSMetricsPass

            nHsMetric1=$runDir/stats/$pair2.proj.md.jr.bam.picHSMetrics
            nHsMetric2=$runDir/stats/$pair3.proj.md.jr.bam.picHSMetrics
            tHsMetric=$runDir/stats/$pair1.proj.md.jr.bam.picHSMetrics

            echo "### normal DAT: DAD: $normalDatFile1, $pair2AssayNo Mom: $normalDatFile2, $pair3AssayNo"
            echo "### tumor  DAT: $tumorDatFile, $pair1AssayNo"

            bedFileGrep=$pair1KitName"_CNABED"
            bedFile=`grep "@@"$recipe"@@" $constants | grep @@"$bedFileGrep"= | cut -d= -f2`
            echo "matchedNormal = ${matchedNormal}"
            if [[ ${matchedNormal} == "No" ]] ; then
                echo "BED file will be swapped out since it has no matched normal"
                if [[ "$pair1KitName" == *S5U ]]  ; then
                    bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_V5_PlusUTR_hs37d5_PaddedTargets_Picard_DGV_1kg_100_interval.bed"

                elif [[ "$pair1KitName" == *S6X ]]  ; then
                    bedFile="/home/tgenref/pecan/annotations/CN_filter/Agilent_V6_noUTR_hs37d5_PaddedTargets_Picard_DGV_1kg_100_interval.bed"
                else
                    echo "We don't have a bed filter for this kit for no matched normal at this time"
                fi
            fi
            echo "### BED FILE= $bedFile"

            if [[ ! -e $normalDatFile1 || ! -e $tumorDatFile || ! -e $normalDatPass1 || ! -e $tumorDatPass || ! -e $normalDatPass2 || ! -e $normalDatFile2 ]] ; then
                echo "### Normal, tumor dat file or pass does not exist"
                ((qsubFails++))
                #continue
            fi
            if [[ ! -e $vcfFile || ! -e $vcfPass ]] ; then
                echo "### HC vcf file called on DNAPAIR not found or not passed yet: $vcfFile"
                ((qsubFails++))
                continue
            fi
            cnaDir="$runDir/cna"
            if [ ! -d $cnaDir ] ; then
                mkdir -p $cnaDir
            fi
            #mkdir -p $cnaDir/$usableName
            #trackName="$runDir/cna/$usableName/$usableName"
            case $pair1AssayNo$pair2AssayNo in
            GenomeGenome) echo "### Both normal and tumor are Genome"
                    assay="Genome"
                    if [ "$cnaGenomeFilter" = "yes" ] ; then
                        outName1=$outName"_filt"
                        mkdir -p $cnaDir/$outName1
                        trackName="$runDir/cna/$outName1/$outName1"
                        if [[ -e $trackName.cnaGenFiltInQueue || -e $trackName.cnaGenFiltPass || -e $trackName.cnaGenFiltFail ]] ; then
                            echo "CNA gen filt plot already in queue, passed, or failed"
                        else
                            echo "### Submitting to queue with $trackName"
                            sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export TRACKNAME=$trackName,MYPATH=$cnaDir/$outName1,CCDSLIST=$ccdsPath,CNAPATH=$cnaPath,MERGEDVCF=$vcfFile,NORMALSAMPLE=$pair1SamName,TUMORSAMPLE=$pair2SamName,NORMALDAT=$normalDatFile,TUMORDAT=$tumorDatFile,OFILE=$outName1,ASSAY=$assay,GTF=$gtf,NXT1=$nxtStep1,REF=$ref,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_cnaGenFilt15.sh
                            if [ $? -eq 0 ] ; then
                                touch $trackName.cnaGenFiltInQueue
                                sleep 1
                            else
                                ((qsubFails++))
                            fi
                            sleep 2
                        fi
                        if [[ -e $trackName.cnaGenFiltInQueue || -e $trackName.cnaGenFiltPass || -e $trackName.cnaGenFiltFail ]] ; then
                                                        echo "CNA gen filt plot already in queue, passed, or failed"
                                                else
                                                        echo "### Submitting to queue with $trackName"
                                                        sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export TRACKNAME=$trackName,MYPATH=$cnaDir/$outName1,CCDSLIST=$ccdsPath,CNAPATH=$cnaPath,MERGEDVCF=$vcfFile,NORMALSAMPLE=$pair1SamName,TUMORSAMPLE=$pair2SamName,NORMALDAT=$normalDatFile,TUMORDAT=$tumorDatFile,OFILE=$outName1,ASSAY=$assay,GTF=$gtf,NXT1=$nxtStep1,REF=$ref,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_cnaGenFilt15.sh
                                                        if [ $? -eq 0 ] ; then
                                                                touch $trackName.cnaGenFiltInQueue
                                                                sleep 1
                                                        else
                                                                ((qsubFails++))
                                                        fi
                                                        sleep 2
                                                fi

                    fi
                    if [ "$cnaGenomeUnFilter" = "yes" ] ; then
                        outName2=$outName"_unfi"
                        mkdir -p $cnaDir/$outName2
                        trackName="$runDir/cna/$outName2/$outName2"
                        if [[ -e $trackName.cnaGenUnfiInQueue || -e $trackName.cnaGenUnfiPass || -e $trackName.cnaGenUnfiFail ]] ; then
                            echo "CNA gen filt plot already in queue, passed, or failed"
                        else
                            echo "### Submitting to queue with $trackName"
                            sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export TRACKNAME=$trackName,MYPATH=$cnaDir/$outName2,CCDSLIST=$ccdsPath,CNAPATH=$cnaPath,MERGEDVCF=$vcfFile,NORMALSAMPLE=$pair1SamName,TUMORSAMPLE=$pair2SamName,NORMALDAT=$normalDatFile,TUMORDAT=$tumorDatFile,OFILE=$outName2,ASSAY=$assay,GTF=$gtf,NXT1=$nxtStep1,REF=$ref,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_cnaGenUnfi15.sh
                            if [ $? -eq 0 ] ; then
                                touch $trackName.cnaGenUnfiInQueue
                                sleep 1
                            else
                                ((qsubFails++))
                            fi
                            sleep 2
                        fi
                        if [[ -e $trackName.cnaGenUnfiInQueue || -e $trackName.cnaGenUnfiPass || -e $trackName.cnaGenUnfiFail ]] ; then
                                                        echo "CNA gen filt plot already in queue, passed, or failed"
                                                else
                                                        echo "### Submitting to queue with $trackName"
                                                        sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export TRACKNAME=$trackName,MYPATH=$cnaDir/$outName2,CCDSLIST=$ccdsPath,CNAPATH=$cnaPath,MERGEDVCF=$vcfFile,NORMALSAMPLE=$pair1SamName,TUMORSAMPLE=$pair2SamName,NORMALDAT=$normalDatFile,TUMORDAT=$tumorDatFile,OFILE=$outName2,ASSAY=$assay,GTF=$gtf,NXT1=$nxtStep1,REF=$ref,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_cnaGenUnfi15.sh
                                                        if [ $? -eq 0 ] ; then
                                                                touch $trackName.cnaGenUnfiInQueue
                                                                sleep 1
                                                        else
                                                                ((qsubFails++))
                                                        fi
                                                        sleep 2
                                                fi
                    fi

                    ;;
            ExomeExome) echo "### Both normal and tumor are Exome"
                    assay="Exome"
                    outNameDad=$pair2SamName"-"$pair1SamName"_exo"
                    outNameMom=$pair3SamName"-"$pair1SamName"_exo"
                    outName3=$outName"_exo"

                    if [[ ! -e $nHsMetricPass1 || ! -e $tHsMetricPass || ! -e $tHsMetric || ! -e $nHsMetric1 ]] ; then
                        echo "### Normal, tumor Hs Metrics Pass or the hs metric file itself does not exist"
                        ((qsubFails++))
                        continue
                    fi
                    nHetDepth1=`cat $nHsMetric1 | grep -A 1 BAIT_SET | tail -1 | cut -f22`
                    tHetDepth=`cat $tHsMetric | grep -A 1 BAIT_SET | tail -1 | cut -f22`

                    mkdir -p $cnaDir/$outNameDad
                    trackName="$runDir/cna/$outNameDad/$outNameDad"
                    if [[ -e $trackName.cnaExomeInQueue || -e $trackName.cnaExomePass || -e $trackName.cnaExomeFail ]] ; then
                        echo "CNA gen filt plot already in queue, passed, or failed"
                    else
                        echo "### Submitting to queue with $trackName"
                        sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export NHETDEPTH=$nHetDepth1,THETDEPTH=$tHetDepth,TRACKNAME=$trackName,CCDSLIST=$ccdsPath,MYPATH=$cnaDir/$outNameDad,CNAPATH=$cnaPath,MERGEDVCF=$vcfFile,NORMALSAMPLE=$pair2SamName,TUMORSAMPLE=$pair1SamName,CNAEXOMETARGET=$bedFile,NORMALDAT=$normalDatFile1,TUMORDAT=$tumorDatFile,OFILE=$outNameDad,ASSAY=$assay,GTF=$gtf,NXT1=$nxtStep1,REF=$ref,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_cnaExo15.sh
                        if [ $? -eq 0 ] ; then
                            touch $trackName.cnaExomeInQueue
                            sleep 1
                        else
                            ((qsubFails++))
                        fi
                        sleep 2
                    fi

                    if [[ ! -e $nHsMetricPass2 || ! -e $tHsMetricPass || ! -e $tHsMetric || ! -e $nHsMetric2 ]] ; then
                        echo "### Normal, tumor Hs Metrics Pass or the hs metric file itself does not exist"
                        ((qsubFails++))
                        continue
                    fi
                    nHetDepth2=`cat $nHsMetric2 | grep -A 1 BAIT_SET | tail -1 | cut -f22`
                    tHetDepth=`cat $tHsMetric | grep -A 1 BAIT_SET | tail -1 | cut -f22`

                    mkdir -p $cnaDir/$outNameMom
                    trackName="$runDir/cna/$outNameMom/$outNameMom"
                    if [[ -e $trackName.cnaExomeInQueue || -e $trackName.cnaExomePass || -e $trackName.cnaExomeFail ]] ; then
                        echo "CNA gen filt plot already in queue, passed, or failed"
                    else
                        echo "### Submitting to queue with $trackName"
                        sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export NHETDEPTH=$nHetDepth2,THETDEPTH=$tHetDepth,TRACKNAME=$trackName,CCDSLIST=$ccdsPath,MYPATH=$cnaDir/$outNameMom,CNAPATH=$cnaPath,MERGEDVCF=$vcfFile,NORMALSAMPLE=$pair3SamName,TUMORSAMPLE=$pair1SamName,CNAEXOMETARGET=$bedFile,NORMALDAT=$normalDatFile2,TUMORDAT=$tumorDatFile,OFILE=$outNameMom,ASSAY=$assay,GTF=$gtf,NXT1=$nxtStep1,REF=$ref,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_cnaExo15.sh
                        if [ $? -eq 0 ] ; then
                            touch $trackName.cnaExomeInQueue
                            sleep 1
                        else
                            ((qsubFails++))
                        fi
                        sleep 2
                    fi
                    ;;
            *)    echo "### Something else"
                echo "### Do not mix and match exomes and genomes!!!"
                    assay="Mixed"
                    ;;
            esac
        done<$pedFile
        else
            echo "###Something is wrong. Didn't find a DNAFAMI or DNAPAIR line"
            exit
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
