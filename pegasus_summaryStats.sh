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

thisStep="pegasus_nextJob_summaryStats.txt"
#nxtStep1="pegasus_nextJob_saveToIsilon.txt"
nxtStep1="pegasus_nextJob_finalize.txt"

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

sumStatsPath=`grep @@"$recipe"@@ $constants | grep @@SUMSTATSPATH= | cut -d= -f2`

nCores=`grep @@${myName}_CORES= $constantsDir/$recipe | cut -d= -f2`


email=`cat $configFile | grep "^EMAIL=" | cut -d= -f2 | head -1 | tr -d [:space:]`

ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

qsubFails=0
missingStuff=0
if [ ! -e $runDir/project.finished ] ; then
    echo "Can't find $runDir/project.finished"
    ((missingStuff++))
fi
if [ $missingStuff -ne 0 ] ; then
    echo "Exiting!!!"
    exit
fi

alreadyDone=0
if [ -e $runDir/summaryStatsPass ] ; then
    echo "Summary stats already passed!"
    touch $runDir/pegasus_nextJob_finalize.txt
    alreadyDone=1
fi
if [ -e $runDir/summaryStatsFail ] ; then
    echo "Summary stats already failed!"
fi
if [ -e $runDir/summaryStatsInQueue ] ; then
    echo "Summary stats still in queue!"
    alreadyDone=1
fi
if [ $alreadyDone -eq 0 ] ; then
    echo "submitting $runDir to queue for Summary stats"
    sbatch --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores --export SUMSTATSPATH=$sumStatsPath,RUNDIR=$runDir,EMAIL=$email,NXT1=$nxtStep1,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_summaryStats.sh
    if [ $? -eq 0 ] ; then
        touch $runDir/summaryStatsInQueue
    else
        ((qsubFails++))
    fi
    sleep 2
fi

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
