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

thisStep="pegasus_nextJob_sleuth.txt"
nxtStep1="pegasus_nextJob_postSleuth.txt"

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
sleuthPath=`grep "@@"$recipe"@@" $constants | grep @@SLEUTHPATH= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
picardPath=`grep @@"$recipe"@@ $constants | grep @@PICARDPATH= | cut -d= -f2`
echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
d=`echo $runDir | cut -c 2-`

qsubFails=0

for rnaPairLine in `cat $configFile | grep ^RNAPAIR=`
do
	
	echo "### RNA pair line is $rnaPairLine"
        sampleNames=`echo $rnaPairLine | cut -d= -f2`
        normls=`echo $sampleNames | cut -d, -f1`
        tumors=`echo $sampleNames | cut -d, -f2`
        altNam=`echo $sampleNames | cut -d, -f3`
        sleuthName=`echo $sampleNames | cut -d, -f3`
        if [ "$sleuthName" == "" ] ; then
                normlsWithDash=${normls//;/-}
                tumorsWithDash=${tumors//;/-}
                sleuthName="$normlsWithDash-VS-$tumorsWithDash"
        fi

        echo "### Normals: $normls"
        echo "### Tumors : $tumors"
	tumorCount=0
        normlCount=0
        missingTumorCount=0
        missingNormlCount=0
        normlList=""
        tumorList=""
        normlList2=""
        tumorList2=""
	
	sleuthDir=$runDir/sleuth/$sleuthName
	if [ ! -d $sleuthDir ] ; then
		mkdir -p $sleuthDir
	fi
	
	sleuthConfig=$sleuthDir/$sleuthName.config
	echo "sample	condition" > $sleuthConfig
	
	for eachNorml in ${normls//;/ }
                do
                        ((normlCount++))
                        echo "eachnorml: $eachNorml"
			echo "$eachNorml	unaffected" >> $sleuthConfig
                        sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachNorml"'"'`
                        kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
                        samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
                        kalllistoOut=$runDir/kallisto/$samName
                        kallistoPass=$runDir/kallisto/$samName.kallistoPass
                        if [ ! -e $kallistoPass ] ; then
                                echo "### Can't find normal kallisto pass"
                                ((missingNormlCount++))
                        else
                                normlList="$kallistoOut;$normlList"
                                normlList2="$eachNorml-$normlList2"
                        fi
                done
                for eachTumor in ${tumors//;/ }
                do
                        ((tumorCount++))
                        echo "eachtumor: $eachTumor"
			echo "$eachTumor	affected" >> $sleuthConfig
                        sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachTumor"'"'`
                        kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
                        samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
			kalllistoOut=$runDir/kallisto/$samName
                        kallistoPass=$runDir/kallisto/$samName.kallistoPass
                        if [ ! -e $kallistoPass ] ; then
                                echo "### Can't find tumor kallisto pass"
                                ((missingTumorCount++))
                        else
                                tumorList="$kallistoOut;$tumorList"
                                tumorList2="$eachTumor-$tumorList2"
                        fi
                done
                if [[ $missingNormlCount -eq 0 && $missingTumorCount -eq 0 ]] ; then
                        echo "### All of normal and tumor bams are found"
                        echo "### Norml list: $normlList"
                        echo "### Tumor list: $tumorList"
                        normlList="${normlList%?}"
                        tumorList="${tumorList%?}"
                        normlList2="${normlList2%?}"
                        tumorList2="${tumorList2%?}"
                else
                        echo "### There were missing things"
                        echo "### Normls missing: $missingNormlCount/$normlCount"
                        echo "### Tumors missing: $missingTumorCount/$tumorCount"
                        ((qsubFails++))
                        continue
                fi

		kallistoDir=$runDir/kallisto
		sleuthOut=$runDir/sleuth/$sleuthName/${sleuthName}_sleuth.txt
		objectData=$runDir/sleuth/$sleuthName/${sleuthName}_sleuth.rda

                if [[ -e $sleuthDir.sleuthPass || -e $sleuthDir.sleuthFail || -e $sleuthDir.sleuthInQueue ]] ; then
                        echo "### Sleuth is already done, failed or inQueue"
                        continue
                fi
                echo "### Submitting $normlList2-VS-$tumorList2 to queue for sleuth..."
                sbatch -n 1 -N 1 --cpus-per-task $nCores --export SLEUTHPATH=$sleuthPath,SLEUTHCONFIG=$sleuthConfig,RUNDIR=$runDir,SLEUTHOUTDIR=$sleuthDir,OBJECTDATA=$objectData,SLEUTHOUTFILE=$sleuthOut,KALLISTOOUT=$kallistoDir,NXT1=$nxtStep1,D=$d $pegasusPbsHome/pegasus_sleuth.pbs
                if [ $? -eq 0 ] ; then
                       touch $sleuthDir.sleuthInQueue
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
