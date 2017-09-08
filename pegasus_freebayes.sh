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

thisStep="pegasus_nextJob_freebayes.txt"
nxtStep1="pegasus_nextJob_snpEff.txt"
nxtStep2="pegasus_nextJob_germVcfMerger.txt"

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
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
freebayesPath=`grep @@"$recipe"@@ $constants | grep @@FREEBAYESPATH= | cut -d= -f2`
snps=`grep @@"$recipe"@@ $constants | grep @@SNPS= | cut -d= -f2`
chrListBed=`grep @@"$recipe"@@ $constants | grep @@CHRBED= | cut -d= -f2`
jirRequested=`grep @@"$recipe"@@ $constants | grep @@JIRREALIGN= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
echo "### reference: $ref"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
for dnaPairLine in `cat $configFile | grep '^DNAPAIR=\|^DNAFAMI='`
do
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
        sampleCount=0
        missingSampleCount=0
        sampleList=""
        for eachSample in ${sampleNames//,/ }
        do
                ((sampleCount++))
                #echo "eachsample: $eachSample"
                sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
                kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
                samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
		if [ "${jirRequested}" == "no" ] ; then

                        eachSampleBam=$runDir/$kitName/$samName/$samName.proj.md.bam
                        eachSamplePass=$runDir/$kitName/$samName/$samName.proj.bam.mdPass
                else
                	eachSampleBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
                	eachSamplePass=$runDir/$kitName/$samName/$samName.proj.md.jr.bam.jointIRPass
		fi

                if [[ ! -e $eachSampleBam || ! -e $eachSamplePass ]] ; then
                        echo "### Can't find the jr bam or the jr pass"
                        echo "### BAM: $eachSampleBam"
                        echo "### PAS: $eachSamplePass"
                        ((missingSampleCount++))
                else
                        sampleList="-b $eachSampleBam $sampleList"
                fi
        done
        if [ $missingSampleCount -eq 0 ] ; then
                echo "### All sample bams are found"
                sampleList="${sampleList%?}"
                echo "### Sample list: $sampleList"
        else
                echo "### There were missing things"
                echo "### Samples missing: $missingSampleCount/$sampleCount"
                ((qsubFails++))
                continue
        fi
 	fbDir="$runDir/freebayes"
	if [ ! -d $fbDir ] ; then
		mkdir $fbDir
	fi
	mkdir -p $fbDir/$usableName
	workDir="$fbDir/$usableName"
	trackName="$runDir/freebayes/$usableName/$usableName"

        STEP=0
        STEP_COUNT=`ls $chrListBed/*bed | wc -l`
        echo "### Submitting to queue to run freebayes on $wd"
        while [ ${STEP} -lt $STEP_COUNT ]
        do
                (( STEP++ ))
                if [[ -e ${trackName}_Step${STEP}.freebayesInQueue || -e ${trackName}_Step${STEP}.freebayesPass || -e ${trackName}_Step${STEP}.freebayesFail ]] ; then
                        echo "### Freebayes is already done, failed, or inqueue for ${trackName}"
                        continue
                fi
	        echo Starting freebayes caller Step${STEP}
		sbatch --export FREEBAYESPATH=$freebayesPath,GATKPATH=$gatkPath,BAMLIST="'$sampleList'",TRACKNAME=$trackName,KNOWN=$snps,STEP=${STEP},STEPCOUNT=$STEP_COUNT,CHRLIST=$chrListBed,FBBAM=$fbBam,REF=$ref,NXT1=$nxtStep1,NXT2=$nxtStep2,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_freebayesMulti.pbs
		if [ $? -eq 0 ] ; then
			touch ${trackName}_Step${STEP}.freebayesInQueue
		else
			((qsubFails++))
		fi
                sleep 2
        done
done

for sampleLine in `cat $configFile | grep ^SAMPLE=`
do
	echo "sample is $sampleLine"
	kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
	samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
	assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
	libraID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
	echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID"
	if [ "$assayID" != "RNA" ] ; then
		echo "### Assay ID is $assayID. Must be genome or exome."
		pcDir=$runDir/$kitName/$samName
		inBam=$runDir/$kitName/$samName/$samName.proj.bam
		mdBam=$runDir/$kitName/$samName/$samName.proj.md.bam
		jrBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
		jrRequested=`cat $configFile | grep '^DNAPAIR=\|^DNAFAMI=' | grep $samName | wc -l`
		if [ $jrRequested -gt 0 ] ; then
			echo "### freebayes will not operate on a single bam because Joint IR is requested for $samName"
			echo "setting to JIRbam"
			fbBam=${jrBam}
			fbPass=${jrBam}.jointIRPass
		else
			
			echo "### freebayes will run on single bam because Joint IR is NOT requested for $samName"
			echo "setting to mdBAM"
			fbBam=${mdBam}
			fbPass=$inBam.mdPass
		#fi
			if [[ ! -e $fbPass || ! -e $fbBam ]] ; then
				echo "### Either mdPass or the bam itself is missing for $fbBam"
				((qsubFails++))
			else
				fbDir="$runDir/freebayes"
				if [ ! -d $fbDir ] ; then
					mkdir $fbDir
				fi
				mkdir -p $fbDir/$samName
				workDir="$fbDir/$samName"
				trackName="$runDir/freebayes/$samName/$samName"
				STEP=0
				STEP_COUNT=`ls $chrListBed/*bed | wc -l`
				##echo "### Submitting to queue to run freebayes on $wd"
				while [ ${STEP} -lt $STEP_COUNT ]
				do
					(( STEP++ ))
					
					echo "### Submitting to queue to run freebayes on STEP: ${STEP}  $wd"
					if [[ -e ${trackName}_Step${STEP}.freebayesInQueue || -e ${trackName}_Step${STEP}.freebayesPass || -e ${trackName}_Step${STEP}.freebayesFail ]] ; then
						echo "### freebayes is already done, failed, or inqueue for ${fbBam}"
						continue
					fi

					echo Starting freebayes for ${fbBam}
					sbatch --export FREEBAYESPATH=$freebayesPath,GATKPATH=$gatkPath,TRACKNAME=$trackName,KNOWN=$snps,STEP=${STEP},STEPCOUNT=$STEP_COUNT,CHRLIST=$chrListBed,FBBAM=$fbBam,REF=$ref,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pegasusPbsHome/pegasus_freebayes.pbs
					if [ $? -eq 0 ] ; then
						touch ${trackName}_Step${STEP}.freebayesInQueue
					else
						((qsubFails++))
					fi

					sleep 2
				done
			fi
		fi
	else
		echo "### Assay ID is $assayID. Must be RNA."
		#code for calling AS metrics on tophat bams here
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
