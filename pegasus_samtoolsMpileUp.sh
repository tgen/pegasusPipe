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

thisStep="pegasus_nextJob_samtoolsMpileUp.txt"
nxtStep1="pegasus_nextJob_snpEff.txt"
nxtStep2="pegasus_nextJob_germVcfMerger.txt"
pbsHome="/home/mrussell/pegasus-pipe/jobScripts"
constants="/home/mrussell/central-pipe/constants/constants.txt"
constantsDir="/home/mrussell/central-pipe/constants"
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
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
bcftoolsPath=`grep "@@"$recipe"@@" $constants | grep @@BCFTOOLSPATH= | cut -d= -f2`
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
        mpDir="$runDir/mpileup"
        if [ ! -d $mpDir ] ; then
                mkdir $mpDir
        fi
	bamText="$mpDir/${usableName}_bams.txt"
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
                        echo "### Can't find the bam or the pass"
                        echo "### BAM: $eachSampleBam"
                        echo "### PAS: $eachSamplePass"
                        ((missingSampleCount++))
                else
                        echo "$eachSampleBam" >> $bamText
			#sampleList="$eachSampleBam $sampleList"
                fi
        done
        if [ $missingSampleCount -eq 0 ] ; then
                echo "### All sample bams are found"
                sampleList="${sampleList%?}"
                #echo "### Sample list: $sampleList"
        else
                echo "### There were missing things"
                echo "### Samples missing: $missingSampleCount/$sampleCount"
                ((qsubFails++))
                continue
        fi
        mkdir -p $mpDir/$usableName
        workDir="$mpDir/$usableName"
        trackName="$runDir/mpileup/$usableName/$usableName"
        STEP=0
        #STEP_COUNT=24
        STEP_COUNT=`ls $chrListBed/*bed | wc -l`
        echo "### Submitting to queue to run mpileup on $wd"
        while [ ${STEP} -lt $STEP_COUNT ]
        do
                (( STEP++ ))
                if [[ -e ${trackName}_Step${STEP}.samtoolsMpileUpInQueue || -e ${trackName}_Step${STEP}.samtoolsMpileUpPass || -e ${trackName}_Step${STEP}.samtoolsMpileUpFail ]] ; then
                        echo "### mpileup is already done, failed, or inqueue for ${trackName}"
                        continue
                fi
                        echo Starting mpileup Step${STEP}
                        ##qsub -A $debit -l nodes=1:ppn=$nCores -v GATKPATH=$gatkPath,STEPCOUNT=$STEP_COUNT,TRK=$trackName,KNOWN=$snps,BAMLIST="'$sampleList'",TRK=$trackName,CHRLIST=$chrList,REF=$ref,STEP=${STEP},NXT1=$nxtStep1,NXT2=$nxtStep2,RUNDIR=$runDir,D=$d $pbsHome/pegasus_haplotypeCaller.pbs
			qsub -A $debit -l nodes=1:ppn=8 -v BEDFILE=$targets,GATKPATH=$gatkPath,SAMTOOLSPATH=$samtoolsPath,BCFTOOLSPATH=$bcftoolsPath,CHRLIST=$chrListBed,TRACKNAME=$trackName,STEP=${STEP},STEPCOUNT=${STEP_COUNT},KNOWN=$snps,BAMFILE=$bamText,REF=$ref,NXT1=$nxtStep1,NXT2=$nxtStep2,RUNDIR=$runDir,D=$d $pbsHome/pegasus_samtoolsMpileUpMulti.pbs
                        if [ $? -eq 0 ] ; then
                                touch ${trackName}_Step${STEP}.samtoolsMpileUpInQueue
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
	trgtGrep=$kitName"_T"
	targets=`grep "@@"$recipe"@@" $constants | grep @@"$trgtGrep"= | cut -d= -f2`
	echo "### TARGETS: $targets"

	echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID"
	if [ "$assayID" != "RNA" ] ; then
		echo "### Assay ID is $assayID. Must be genome or exome."
		pcDir=$runDir/$kitName/$samName
		inBam=$runDir/$kitName/$samName/$samName.proj.bam
		mdBam=$runDir/$kitName/$samName/$samName.proj.md.bam
		jrBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
		jrRequested=`cat $configFile | grep '^DNAPAIR=\|^DNAFAMI=' | grep $samName | wc -l`
		if [ $jrRequested -gt 0 ] ; then
			echo "### samtoolsMpileUp will not operate on a single bam because Joint IR is requested for $samName"
			echo "setting to JIRbam"
			bamFile=${jrBam}
			bamFilePass=${jrBam}.jointIRPass
		else
			
			echo "### samtoolsMpileUp will run on single bam because Joint IR is NOT requested for $samName"
			echo "setting to mdBAM"
			bamFile=${mdBam}
			bamFilePass=$inBam.mdPass
#		fi
			if [[ ! -e $bamFilePass || ! -e $bamFile ]] ; then
				echo "### Either mdPass or the bam itself is missing for $bamFile"
				((qsubFails++))
			else
				STEP=0
				STEP_COUNT=`ls $chrListBed/*bed | wc -l`
				##echo "### Submitting to queue to run freebayes on $wd"
				while [ ${STEP} -lt $STEP_COUNT ]
				do
					(( STEP++ ))

					echo "### Submitting to queue to run samtoolsMpileUp on $wd"
					#trackName=$runDir/$kitName/$samName/$samName
					mpDir="$runDir/mpileup"
					if [ ! -d $mpDir ] ; then
						mkdir $mpDir
					fi
					mkdir -p $mpDir/$samName
					workDir="$mpDir/$samName"
					trackName="$runDir/mpileup/$samName/$samName"

					if [[ -e ${trackName}_Step${STEP}.samtoolsMpileUpInQueue || -e ${trackName}_Step${STEP}.samtoolsMpileUpPass || -e ${trackName}_Step${STEP}.samtoolsMpileUpFail ]] ; then
						echo "### samtoolsMpileUp is already done, failed, or inqueue for ${bamFile}"
						continue
					fi

					echo Starting samtoolsMpileUp for ${bamFile}
					qsub -A $debit -l nodes=1:ppn=8 -v BEDFILE=$targets,GATKPATH=$gatkPath,SAMTOOLSPATH=$samtoolsPath,BCFTOOLSPATH=$bcftoolsPath,CHRLIST=$chrListBed,TRACKNAME=$trackName,STEP=${STEP},STEPCOUNT=${STEP_COUNT},KNOWN=$snps,BAMFILE=$bamFile,REF=$ref,NXT1=$nxtStep1,NXT2=$nxtStep2,RUNDIR=$runDir,D=$d $pbsHome/pegasus_samtoolsMpileUp.pbs
					if [ $? -eq 0 ] ; then
						touch ${trackName}_Step${STEP}.samtoolsMpileUpInQueue
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
