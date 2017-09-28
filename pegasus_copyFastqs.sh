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

thisStep="Pegasus_nextJob_copyFastqs.txt"
#nxtStep1="pegasus_nextJob_runFastQC.txt"
nxtStep2="pegasus_nextJob_dnaAlign.txt"
nxtStep3="pegasus_nextJob_mergeFastqs.txt"
nxtStep4="pegasus_nextJob_splitFastqs.txt"


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
qsubFails=0
projName=`basename $runDir | awk -F'_ps20' '{print $1}'`
configFile=$runDir/$projName.config
echo "### projName: $projName"
echo "### confFile: $configFile"
if [ ! -e $configFile ] ; then
    echo "### Config file not found at $configFile!!!"
    echo "### Exiting!!!"
    exit
else
    echo "### Config file found."
fi
recipe=`cat $configFile | grep "^RECIPE=" | cut -d= -f2 | head -1 | tr -d [:space:]`
echo "### Recipe is $recipe"
#for fqLine in `cat $configFile | grep ^FQ=`
#do
#    fqFile=`echo $fqLine | cut -d= -f2 | cut -d, -f2`
#    r2File=${fqFile/_R1/_R2}
#done
email=`cat $configFile | grep "^EMAIL=" | cut -d= -f2 | tr -d [:space:]`
incFastq=`grep "@@"$recipe"@@" $constants | grep @@INCFASTQ= | cut -d= -f2`
if [[ "$incFastq" != "yes" && "$incFastq" != "no" ]] ; then
    echo "### Valid values for INCFASTQ is either yes or no. You have: $incFastq"
    echo "### Exiting!!!"
    exit
else
    echo "### Inc fastq is good: $incFastq"
fi 
skipLines=1
count=0
sampleCount=0

#Creating a file for a summary of the pipeline run
touch $runDir/ProjectRunSummary.txt

for configLine in `cat $configFile`
do
    if [ "$configLine" == "=START" ] ; then
        skipLines=0
        continue
    fi
    if [ $skipLines -eq 0 ] ; then
        if [[ $configLine == SAMPLE* || $configLine == =END* ]] ; then
            #echo "config line is $configLine"
            arrayCount=${#mergeArray[@]}
            if [ $arrayCount -gt 0 ] ; then
                ((sampleCount++))
                #echo "### Sample with $arrayCount rows found for kit: $kitName, sample: $samName."
                lastItemIndex=`echo ${#mergeArray[@]}-1 | bc`
                for (( i=0; i<=$lastItemIndex; i++ ))
                do
                    ##3/28/16 Megan added check to make sure fastqs have the right extension and R2 exists
                    #echo "array with index $i is::: ${mergeArray[$i]}"
                    thisFq=`echo ${mergeArray[$i]} | cut -d= -f2 | cut -d, -f2`
                    fqExt="_R1_001.fastq.gz"
                    if [[ "$thisFq" != *"$fqExt" ]] ; then
                        echo "###This fastq doesn't have the proper extension: $thisFq"
                        echo "###Exiting"
                        if [ ! -e $runDir/checkConfig.wrongFQextension ] ; then
                            echo "###This fastq doesn't have the proper extension: $thisFq, the proper extension is _R1_001.fastq.gz" >> ~/mailtmp-$$.txt
                            cat ~/mailtmp-$$.txt | mail -s "Pegasus pipeline ERROR: Your project $projName has an ERROR" "jetstream@tgen.org"
                            cat ~/mailtmp-$$.txt | mail -s "Pegasus pipeline ERROR: Your project $projName has an ERROR" "${email}"
                            mv ~/mailtmp-$$.txt $runDir/checkConfig.wrongFQextension
                        fi
                        exit
                    fi
                    #r2File=${thisFq/_R1/_R2}
                    r2File=`echo $thisFq | sed 's/\(.*\)_R1_/\1_R2_/'`
                    if [[ ! -f $r2File ]] ; then
                        echo "FILE: $r2File not found, EXITING"
                        if [ ! -e $runDir/checkConfig.copyFastqFail ] ; then
                            echo "###FILE: $r2File not found" >> ~/mailtmp-$$.txt
                            cat ~/mailtmp-$$.txt | mail -s "Pegasus pipeline ERROR: Your project $projName has an ERROR" "jetstream@tgen.org"
                            cat ~/mailtmp-$$.txt | mail -s "Pegasus pipeline ERROR: Your project $projName has an ERROR" "${email}"
                            mv ~/mailtmp-$$.txt $runDir/checkConfig.copyFastqFail
                        fi
                        exit
                    fi
                    targetName="$runDir/$kitName/$samName/$samName"`printf "_%03d" "$i"`"_R1.fastq.gz"
                    targR2Name="$runDir/$kitName/$samName/$samName"`printf "_%03d" "$i"`"_R2.fastq.gz"

                    fastqSize=`stat -Lc%s $thisFq`
                    if [ $fastqSize -ge 36000000000 ] ; then
                        echo "This fastq is quite large, needs to be split, then aligned: $fastqSize"
                        touch $targetName.needsToBeSplit
                    fi
                    #echo "target name is $targetName"
                    if [[ ! -e $targetName.cpFqInQueue && ! -e $targetName.cpFqPass ]] ; then
                        if [ "$incFastq" == "yes" ] ; then
                            echo "### Copying to $targetName"
                            touch $targetName.cpFqInQueue
                            cp $thisFq $targetName
                            if [ $? -eq 0 ] ; then
                                touch $targetName.cpFqPass
                            else
                                touch $targetName.cpFqFail
                                ((qsubFails++))
                            fi
                            rm -f $targetName.cpFqInQueue
                        else #must be no, not copying, just linking
                            echo "### Linking for $thisFq to $targetName"
                            ln -s $thisFq $targetName
                            if [ $? -eq 0 ] ; then
                                touch $targetName.cpFqPass
                            else
                                touch $targetName.cpFqFail
                                ((qsubFails++))
                            fi
                        fi
                    else
                        echo "### Copy already in queue or passed $targetName "
                    fi
                    if [[ ! -e $targR2Name.cpFqInQueue && ! -e $targR2Name.cpFqPass ]] ; then
                        if [ "$incFastq" == "yes" ] ; then
                            echo "### Copying to $targR2Name"
                            touch $targR2Name.cpFqInQueue
                            cp $r2File $targR2Name
                            if [ $? -eq 0 ] ; then
                                touch $targR2Name.cpFqPass
                            else
                                touch $targR2Name.cpFqFail
                                ((qsubFails++))
                            fi
                            rm -f $targR2Name.cpFqInQueue
                        else #must be no, not copying, just linking
                            echo "### Linking for $r2File to $targR2Name"
                            ln -s $r2File $targR2Name
                            if [ $? -eq 0 ] ; then
                                touch $targR2Name.cpFqPass
                            else
                                touch $targR2Name.cpFqFail
                                ((qsubFails++))
                            fi
                        fi
                    else
                        echo "### Copy already in queue or passed $targR2Name "
                    fi

                done

            fi
            kitName=`echo $configLine | cut -d= -f2 | cut -d, -f1`
            samName=`echo $configLine | cut -d= -f2 | cut -d, -f2`
            unset mergeArray
            count=0
            continue
        else #doesnt start with =, add to mergeArray
            #echo "adding $configLine to mergeArray"
            mergeArray[$count]=$configLine
            ((count++))
        fi
    else
        continue
    fi
done

if [ $qsubFails -eq 0 ] ; then
#all jobs submitted succesffully, remove this dir from messages
    echo "### Removing $thisStep from $runDir."
    rm -f $runDir/$thisStep
    #echo "### Touching $nxtStep1 under $runDir"
    #touch $runDir/$nxtStep1
    echo "### Touching $nxtStep2 under $runDir"
    touch $runDir/$nxtStep2
    echo "### Touching $nxtStep3 under $runDir"
    touch $runDir/$nxtStep3
    echo "### Touching $nxtStep4 under $runDir"
    touch $runDir/$nxtStep4

else
#qsub failed at some point, this runDir must stay in messages
    echo "### Failure in qsub. Keeping $thisStep"
fi

time=`date +%d-%m-%Y-%H-%M`
echo "Ending $0 at $time"
