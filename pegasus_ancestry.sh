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

thisStep="pegasus_nextJob_ancestry.txt"
nxtStep1="pegasus_nextJob_postAncestry.txt"
pbsHome="~/pegasus-pipe/jobScripts"
constants="~/central-pipe/constants/constants.txt"
constantsDir="~/central-pipe/constants"
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
laserPath=`grep "@@"$recipe"@@" $constants | grep @@LASERPATH= | cut -d= -f2`
hgdpPath=`grep "@@"$recipe"@@" $constants | grep @@HGDPPATH= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
bcftoolsPath=`grep "@@"$recipe"@@" $constants | grep @@BCFTOOLSPATH= | cut -d= -f2`
snps=`grep @@"$recipe"@@ $constants | grep @@SNPS= | cut -d= -f2`
chrListBed=`grep @@"$recipe"@@ $constants | grep @@CHRBED= | cut -d= -f2`
species=`grep @@SPECIES= $constantsDir/$recipe | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### rnaAlign: $rnaAligner"
echo "### reference: $ref"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0

if [ "$species" != "HUMAN" ] ; then

	echo "The ancestry script only works on HUMAN samples"
	echo "This species is $species"
	echo "Will exit this script"
	exit
fi

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
		trackName="$runDir/ancestry/$samName/$samName"
		bamFile="$runDir/$kitName/$samName/$samName.proj.md.bam"
		bamFilePass="$runDir/$kitName/$samName/$samName.proj.bam.mdPass"
		if [[ ! -e $bamFilePass || ! -e $bamFile ]] ; then
			echo "### Either mdPass or the bam itself is missing for $bamFile"
			((qsubFails++))
		else

			if [[ -e ${trackName}.ancestryInQueue || -e ${trackName}.ancestryPass || -e ${trackName}.ancestryFail ]] ; then
				echo "### ancestry is already done, failed, or inqueue for ${bamFile}"
				continue
			fi
			ancestryDir="$runDir/ancestry/$samName"
			if [ ! -d $ancestryDir ] ; then
               			mkdir -p $ancestryDir
        		fi
			echo "Starting ancestry for ${bamFile}"
			sbatch --export ANCESTRYDIR=$ancestryDir,BEDFILE=$targets,GATKPATH=$gatkPath,SAMTOOLSPATH=$samtoolsPath,LASERPATH=$laserPath,HGDPPATH=$hgdpPath,TRACKNAME=$trackName,KNOWN=$snps,BAMFILE=$bamFile,REF=$ref,NXT1=$nxtStep1,RUNDIR=$runDir,D=$d $pbsHome/pegasus_ancestry.sh
			if [ $? -eq 0 ] ; then
				touch ${trackName}.ancestryInQueue
			else
				((qsubFails++))
				sleep 2
			fi
		fi
	else
		echo "### Assay ID is $assayID. This is for DNA only."
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
