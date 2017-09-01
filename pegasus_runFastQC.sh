#!/bin/bash
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

thisStep="pegasus_nextJob_runFastQC.txt"
nxtStep1="pegasus_nextJob_postFastQC.txt"
pbsHome="/home/tgenjetstream/pegasus-pipe/jobScripts"
constants="/home/tgenjetstream/central-pipe/constants/constants.txt"
constantsDir="/home/tgenjetstream/central-pipe/constants"
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
echo "$nCores"
echo "$myName"

fastqcPath=`grep "@@"$recipe"@@" $constants | grep @@FASTQCPATH= | cut -d= -f2`
echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### fastqcpt: $fastqcPath"

qsubFails=0
for cpFqPassFile in `find $runDir -name *cpFqPass`
do
	thisFq=${cpFqPassFile/.cpFqPass/}
	if [ ! -e $thisFq ] ; then
		echo "### Fastq itself doesnt exist. Weird: $thisFq"
		((qsubFails++))
		continue
	fi
	if [[ -e $thisFq.fastqcPass || -e $thisFq.fastqcInQueue || -e $thisFq.fastqcFail ]] ; then
		echo "### This fastq file already passed, failed, or inQueue"
		continue
	fi
	echo "### Submitting to fastQC for $thisFq"
	d=`echo $runDir | cut -c 2-`
	sbatch -n 1 -N 1 --cpus-per-task $nCores -v FQ=$thisFq,FASTQCPATH=$fastqcPath,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_runFastQC.pbs
	if [ $? -eq 0 ] ; then
		touch $thisFq.fastqcInQueue
	else
		((qsubFails++))
	fi
	sleep 1
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
