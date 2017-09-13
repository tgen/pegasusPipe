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

thisStep="pegasus_nextJob_kallisto.txt"
nxtStep1="pegasus_nextJob_postKallisto.txt"
##nxtStep1="pegasus_nextJob_sleuth.txt"

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
kallistoIndexCDNA=`grep @@KALLISTO_INDEX_cDNA= $constantsDir/$recipe | cut -d= -f2` 
kallistoIndexGTF=`grep @@KALLISTO_INDEX_GTF= $constantsDir/$recipe | cut -d= -f2` 
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
picardPath=`grep @@"$recipe"@@ $constants | grep @@PICARDPATH= | cut -d= -f2`
kallisto=`grep @@KALLISTO= $constantsDir/$recipe | cut -d= -f2`
echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
skipGroup=0

for configLine in `cat $configFile`
do
    if [[ $kallisto == no ]] ; then
                echo "### Salmon not requested for this recipe"
                continue
        fi
        if [ "$configLine" == "=START" ] ; then
                skipLines=0
                continue
        fi
        if [ $skipLines -eq 0 ] ; then
                if [[ $configLine == SAMPLE* || $configLine == =END* ]] ; then
                        echo "config line is $configLine"
                        arrayCount=${#mergeArray[@]}
            echo "arrayCount is: $arrayCount"
                        if [ $arrayCount -gt 0 ] ; then
                                echo "### Starting with $samName"
                                missingFastqSample=0
                                accountFastqSample=0
                                fastqList=""
                                #fastqList2=""
                rgTagList=""
                                read1Count=0
                                read2Count=0
                                ((sampleCount++))
                                echo "### Sample with $arrayCount rows found for kit: $kitName, sample: $samName."
                echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID"
                if [[ "$assayID" == "Exome"  ||  "$assayID" == "Genome" || "$assayID" == "exome" || "$assayID" == "genome" ]] ; then
                                        echo "### Assay ID is $assayID. Will skip."
                                        skipGroup=1
                                fi
                                lastItemIndex=`echo ${#mergeArray[@]}-1 | bc`
                if [ $skipGroup -ne 1 ] ; then
                    echo "Assay is RNA will continue"
                    for (( i=0; i<=$lastItemIndex; i++ ))
                    do
                        #echo "array with index $i is::: ${mergeArray[$i]}"
                        ((accountFastqSample++));((accountFastqSample++))
                        thisFq=`echo ${mergeArray[$i]} | cut -d= -f2 | cut -d, -f2`
                        r2File=${thisFq/_R1/_R2}
                        sourceName="$runDir/$kitName/$samName/$samName"`printf "_%03d" "$i"`"_R1.fastq.gz"
                        #targetName="$runDir/$kitName/$samName/$samName.proj.R1.fastq.gz"
                        sourR2Name="$runDir/$kitName/$samName/$samName"`printf "_%03d" "$i"`"_R2.fastq.gz"
                        #targR2Name="$runDir/$kitName/$samName/$samName.proj.R2.fastq.gz"
                        #RG_CN=TGen
                        #RG_PL=ILLUMINA
                        #RG_ID_fromConfig=`echo ${mergeArray[$i]} | cut -d= -f2 | cut -d, -f1`
                        #FCID=`zcat $sourceName | head -n1 | cut -d: -f3`
                        #LANE=`zcat $sourceName | head -n1 | cut -d: -f4`
                        #INDEX=`zcat $sourceName | head -n1 | cut -d: -f10`
                        #RG_LB=`echo $RG_ID_fromConfig | cut -d_ -f3`
                        #RG_PU="${FCID}_${LANE}"
                        #RG_ID="${FCID}_${LANE}_${RG_LB}"
                        #rgTag="@RG\tID:${RG_ID}\tSM:$samName\tPL:${RG_PL}\tCN:${RG_CN}\tPU:${RG_PU}\tLB:${RG_LB}\tKS:${INDEX}\t"
                        #rgTag="ID:${RG_ID}    SM:$samName    PL:${RG_PL}    CN:${RG_CN}    PU:${RG_PU}    LB:${RG_LB}    KS:${INDEX}"
                        #echo "RG TAG is: $rgTag"
                        #rgTagList="${rgTag} , ${rgTagList}"

                        #echo "target name is $sourceName"
                        if [[ ! -e $sourceName.cpFqPass || ! -e $sourceName ]] ; then
                            echo "### File missing: $sourceName or its .cpFqPass file"
                            ((missingFastqTotal++))
                            ((missingFastqSample++))
                        else
                            echo "### File found $sourceName "
                            fastqList="${fastqList},${sourceName}"
                            ((read1Count++))
                        fi
                        if [[ ! -e $sourR2Name.cpPass && ! -e $sourR2Name ]] ; then
                            echo "### File missing: $sourR2Name or its .cpFqPass file"
                            ((missingFastqTotal++))
                            ((missingFastqSample++))
                        else
                            echo "### File found $sourR2Name "
                            fastqList="${fastqList},${sourR2Name}"
                            ((read2Count++))
                        fi

                        thisPath=`dirname $runDir`
                        cd $thisPath
                        #ownDir=$runDir/$kitName/$samName/${samName}.kallistoDir
                        ownDir=$runDir/kallisto/$samName
                        if [ ! -d $ownDir ] ; then
                            mkdir -p $ownDir
                        fi
                        #creating linked files to the original reads
                        r1Name=`basename $sourceName`
                        r2Name=`basename $sourR2Name`
                        cd $ownDir
                        ln -s $sourceName $r1Name
                        read1Name=$ownDir/$r1Name
                        ln -s $sourR2Name $r2Name
                        read2Name=$ownDir/$r2Name
                        cd -

                    done
                    echo "### Done with $samName with $missingFastqSample missing files out of $accountFastqSample"
                    if [ $missingFastqSample -eq 0 ] ; then
                        #if [ "$assayID" == "RNA" ] ; then
                            echo "### Assay ID is $assayID"
                            echo "### All fastqs were accounted for $samName"
                            echo "### FASTQ files($read1Count): $fastqList"
                            #echo "### R2 files($read2Count): $fastqList"
                            #repeats the rg tag so both R1 and R1 have an RG tag
                            #gets rid of the trailing comma on the rg tag list and fastqList2
                            #rgTagList="${rgTagList}${rgTagList}"
                            #rgTagList=`echo "${rgTagList%?}"`
                            #fastqList=`echo "${fastqList%?}"`
                            #fastqList2=`echo "${fastqList2%?}"`
                            #echo "### RG TAG LIST: $rgTagList"

                            echo "fastqList: ${fastqList}"
                            if [[ -e $ownDir.kallistoPass || -e $ownDir.kallistoFail || -e $ownDir.kallistoInQueue ]] ; then
                                echo "### STAR already done, failed or inQueue"
                                continue
                            fi
                            echo "### submitting $ownDir to queue for Kallisto... "

                             if [[ $rnaStrand == "FIRST" ]] ; then
                                        echo "##running first stranded kallisto case"
                                sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export KALLISTO_INDEX_CDNA=$kallistoIndexCDNA,KALLISTO_INDEX_GTF=$kallistoIndexGTF,SAMNAME=$samName,FASTQL="'"$fastqList"'",DIR=$ownDir,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_firstStrandedKallisto.sh
                                if [ $? -eq 0 ] ; then
                                    touch $ownDir.kallistoInQueue
                                else
                                    ((qsubFails++))
                                    sleep 2
                                fi
                            elif [[ $rnaStrand == "SECOND" ]] ; then
                                        echo "##running second stranded kallisto case"
                                sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export KALLISTO_INDEX_CDNA=$kallistoIndexCDNA,KALLISTO_INDEX_GTF=$kallistoIndexGTF,SAMNAME=$samName,FASTQL="'"$fastqList"'",DIR=$ownDir,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_secondStrandedKallisto.sh
                                                                if [ $? -eq 0 ] ; then
                                                                        touch $ownDir.kallistoInQueue
                                                                else
                                                                        ((qsubFails++))
                                                                        sleep 2
                                                                fi
                            else
                                        echo "##running unstranded kallisto case"
                                sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export KALLISTO_INDEX_CDNA=$kallistoIndexCDNA,KALLISTO_INDEX_GTF=$kallistoIndexGTF,SAMNAME=$samName,FASTQL="'"$fastqList"'",DIR=$ownDir,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_kallisto.sh
                                                                if [ $? -eq 0 ] ; then
                                                                        touch $ownDir.kallistoInQueue
                                                                else
                                                                        ((qsubFails++))
                                                                        sleep 2
                                                                fi
                            fi
                    else
                        echo "### Some files were missing for $samName"
                    fi
                else #else of skipGroup
                    echo "##Skipping Group"
                fi #end of skipGroup
                        fi
                        kitName=`echo $configLine | cut -d= -f2 | cut -d, -f1`
                        samName=`echo $configLine | cut -d= -f2 | cut -d, -f2`
                        assayID=`echo $configLine | cut -d= -f2 | cut -d, -f3`
            rnaStrand=`grep "@@"$kitName"@@" $constants | cut -d= -f2`

                        unset mergeArray
                        count=0
            skipGroup=0
                        continue
                else #doesnt start with =, add to mergeArray
                        #echo "adding $configLine to mergeArray"
                        mergeArray[$count]=$configLine
            #echo "mergeArray[$count]=$configLine"
                        ((count++))
                fi
        else
                continue
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
