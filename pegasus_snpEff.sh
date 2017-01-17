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

thisStep="pegasus_nextJob_snpEff.txt"
nxtStep1="pegasus_nextJob_checkProjectComplete.txt"
##nxtStep2="pegasus_nextJob_germVcfMerger.txt"
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


ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
rnaAligner=`grep "@@"$recipe"@@" $constants | grep @@RNAALIGNER= | cut -d= -f2`
gatkPath=`grep @@"$recipe"@@ $constants | grep @@GATKPATH= | cut -d= -f2`
hapmap=`grep @@"$recipe"@@ $constants | grep @@HAPMAP= | cut -d= -f2`
omni=`grep @@"$recipe"@@ $constants | grep @@OMNI= | cut -d= -f2`

snpeffdb=`grep @@"$recipe"@@ $constants | grep @@SNPEFFDB= | cut -d= -f2`
dbsnp=`grep @@"$recipe"@@ $constants | grep @@SNPS= | cut -d= -f2`
snpeffPath=`grep @@"$recipe"@@ $constants | grep @@SNPEFFPATH= | cut -d= -f2`

echo "### projName: $projName"
echo "### confFile: $configFile"
d=`echo $runDir | cut -c 2-`

skipLines=1
qsubFails=0
####first for loop is for finding seurat vcfs
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
	echo "### DNA pair line is $dnaPairLine for seurat stuff"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
	for eachSample in ${sampleNames//,/ }
	do
		((sampleCount++))
		#echo "eachsample: $eachSample"
		sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
		kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
		samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
		assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
	done
	usableName=${sampleNames//,/-}
	sampleCount=0
	missingSampleCount=0
	sampleList=""
	trackName="$runDir/seurat/$usableName/$usableName"
	if [ ! -e $trackName.seuratPass ] ; then
		echo "### Seurat Pass doesnt exist yet: $trackName.seuratPass"
		((qsubFails++))
		continue
	fi
	if [[ -e $trackName.seurat.vcf.snpEffInQueue || -e $trackName.seurat.vcf.snpEffPass || -e $trackName.seurat.vcf.snpEffFail ]] ; then 
		echo "### snpEff is already done, failed or inQueue"
		continue
	fi 
	echo "### Submitting $trackName.seurat.vcf to queue for snpEff..."
	qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$trackName.seurat.vcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
	if [ $? -eq 0 ] ; then
		touch $trackName.seurat.vcf.snpEffInQueue
	else
		((qsubFails++))
	fi
	sleep 2
done
##first for loop is for finding REVERSE seurat vcfs
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
	echo "### DNA pair line is $dnaPairLine for reverse seurat stuff"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
	for eachSample in ${sampleNames//,/ }
	do
		((sampleCount++))
		#echo "eachsample: $eachSample"
		sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
		kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
		samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
		assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f3`
	done
	usableName=${sampleNames//,/-}
	usableName=${usableName}_REVERSE
	sampleCount=0
	missingSampleCount=0
	sampleList=""
	trackName="$runDir/seurat/$usableName/$usableName"
	if [ ! -e $trackName.REVseuratPass ] ; then
		echo "### Seurat Pass doesnt exist yet: $trackName.REVseuratPass"
		((qsubFails++))
		continue
	fi
	if [[ -e $trackName.REVseurat.vcf.snpEffInQueue || -e $trackName.REVseurat.vcf.snpEffPass || -e $trackName.REVseurat.vcf.snpEffFail ]] ; then 
		echo "### snpEff is already done, failed or inQueue"
		continue
	fi 
	echo "### Submitting $trackName.REVseurat.vcf to queue for snpEff..."
	qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$trackName.REVseurat.vcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
	if [ $? -eq 0 ] ; then
		touch $trackName.REVseurat.vcf.snpEffInQueue
	else
		((qsubFails++))
	fi
	sleep 2
done
###second for loop is for finding strelka vcfs
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
	echo "### DNA pair line is $dnaPairLine for strelka stuff"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
	for eachSample in ${sampleNames//,/ }
	do
		((sampleCount++))
		#echo "eachsample: $eachSample"
		sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
		kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
		samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
		assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
	done
	usableName=${sampleNames//,/-}
	sampleCount=0
	missingSampleCount=0
	sampleList=""
	trackName="$runDir/strelka/$usableName"
	vcfPre="$runDir/strelka/$usableName/myAnalysis/results/$usableName"
	if [ ! -e $trackName.strelkaPass ] ; then
		echo "### Strelka pass doesnt exist yet: $trackName.strelkaPass"
		((qsubFails++))
		continue
	fi
	if [[ -e $vcfPre.strelka.passed.somatic.snvs.vcf.snpEffInQueue || -e $vcfPre.strelka.passed.somatic.snvs.vcf.snpEffPass || -e $vcfPre.strelka.passed.somatic.snvs.vcf.snpEffFail ]] ; then 
		echo "### snpEff is already done, failed or inQueue for passed somatic snvs"
	else
		echo "### Submitting $vcfPre.passed.somatic.snvs.vcf to queue for snpEff..."
		qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$vcfPre.strelka.passed.somatic.snvs.vcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
		if [ $? -eq 0 ] ; then
			touch $vcfPre.strelka.passed.somatic.snvs.vcf.snpEffInQueue 
		else
			((qsubFails++))
		fi
		sleep 2
	fi 
	if [[ -e $vcfPre.strelka.all.somatic.snvs.vcf.snpEffInQueue || -e $vcfPre.strelka.all.somatic.snvs.vcf.snpEffPass || -e $vcfPre.strelka.all.somatic.snvs.vcf.snpEffFail ]] ; then 
		echo "### snpEff is already done, failed or inQueue for all somatic snvs"
	else
		echo "### Submitting $vcfPre.all.somatic.snvs.vcf to queue for snpEff..."
		qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$vcfPre.strelka.all.somatic.snvs.vcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
		if [ $? -eq 0 ] ; then
			touch $vcfPre.strelka.all.somatic.snvs.vcf.snpEffInQueue 
		else
			((qsubFails++))
		fi
		sleep 2
	fi 
	if [[ -e $vcfPre.strelka.passed.somatic.indels.vcf.snpEffInQueue || -e $vcfPre.strelka.passed.somatic.indels.vcf.snpEffPass || -e $vcfPre.strelka.passed.somatic.indels.vcf.snpEffFail ]] ; then 
		echo "### snpEff is already done, failed or inQueue for passed somatic indels"
	else
		echo "### Submitting $vcfPre.passed.somatic.indels.vcf to queue for snpEff..."
		qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$vcfPre.strelka.passed.somatic.indels.vcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
		if [ $? -eq 0 ] ; then
			touch $vcfPre.strelka.passed.somatic.indels.vcf.snpEffInQueue 
		else
			((qsubFails++))
		fi
		sleep 2
	fi 
	if [[ -e $vcfPre.strelka.all.somatic.indels.vcf.snpEffInQueue || -e $vcfPre.strelka.all.somatic.indels.vcf.snpEffPass || -e $vcfPre.strelka.all.somatic.indels.vcf.snpEffFail ]] ; then 
		echo "### snpEff is already done, failed or inQueue for all somatic indels"
	else
		echo "### Submitting $vcfPre.all.somatic.indels.vcf to queue for snpEff..."
		qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$vcfPre.strelka.all.somatic.indels.vcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
		if [ $? -eq 0 ] ; then
			touch $vcfPre.strelka.all.somatic.indels.vcf.snpEffInQueue 
		else
			((qsubFails++))
		fi
		sleep 2
	fi 
done
###third loop is for finding mutect vcfs
for dnaPairLine in `cat $configFile | grep '^DNAPAIR='`
do
	echo "### DNA pair line is $dnaPairLine for mutect stuff"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
	for eachSample in ${sampleNames//,/ }
	do
		((sampleCount++))
		#echo "eachsample: $eachSample"
		sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
		kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
		samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
		assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
	done
	usableName=${sampleNames//,/-}
	sampleCount=0
	missingSampleCount=0
	sampleList=""
	trackName="$runDir/mutect/$usableName/$usableName"
	vcf="${trackName}_MuTect_All.vcf"
	if [ ! -e $trackName.mutectPass ] ; then
		echo "### Strelka pass doesnt exist yet: $trackName.mutectPass"
		((qsubFails++))
		continue
	fi
	if [[ -e $vcf.snpEffInQueue || -e $vcf.snpEffPass || -e $vcf.snpEffFail ]] ; then 
		echo "### snpEff is already done, failed or inQueue for mutect vcf"
		continue
	fi 
	echo "### Submitting $vcf to queue for snpEff..."
	qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$vcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
	if [ $? -eq 0 ] ; then
		touch $vcf.snpEffInQueue 
	else
		((qsubFails++))
	fi
	sleep 2
done

## This loop if for finding multi variant called vcf files from GATK
for dnaPairLine in `cat $configFile | grep '^DNAPAIR=\|^DNAFAMI='`
do
	echo "### DNA pair line is $dnaPairLine"
	sampleNames=`echo $dnaPairLine | cut -d= -f2`
	shorterName=`echo $dnaPairLine | cut -d= -f2 | cut -d\; -f2`
        if [[ $shorterName == *,* ]] ; then
                shorterName=""
                echo "### No nick name detected for these pairs"
                usableName=${sampleNames//,/-}
        else
                echo "### Nick name for these pairs are: $shorterName"
                usableName="$shorterName"
        fi

	for eachSample in ${sampleNames//,/ }
	do
		((sampleCount++))
		#echo "eachsample: $eachSample"
		sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
		kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
		samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
		assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
	done
	sampleCount=0
	missingSampleCount=0
	sampleList=""
	trackName="$runDir/hc/$usableName/$usableName"
	if [ ! -e $trackName.hcPass ] ; then
		echo "### HC pass is not done yet. Missing: $trackName.hcPass"
		((qsubFails++))
		continue
	fi

	if [[ -e $trackName.HC_All.vcf.snpEffInQueue || -e $trackName.HC_All.vcf.snpEffPass || -e $trackName.HC_All.vcf.snpEffFail ]] ; then 
		echo "### snpEff is already done, failed or inQueue"
		continue
	fi 
	echo "### Submitting $trackName.HC_All.vcf to queue for snpEff..."
	qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$trackName.HC_All.vcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
	if [ $? -eq 0 ] ; then
		touch $trackName.HC_All.vcf.snpEffInQueue
	else
		((qsubFails++))
	fi
	sleep 2
done
## This loop if for finding multi variant called vcf files from samtools mpileup
for dnaPairLine in `cat $configFile | grep '^DNAPAIR=\|^DNAFAMI='`
do
        echo "### DNA pair line is $dnaPairLine"
        sampleNames=`echo $dnaPairLine | cut -d= -f2`
        shorterName=`echo $dnaPairLine | cut -d= -f2 | cut -d\; -f2`
        if [[ $shorterName == *,* ]] ; then
                shorterName=""
                echo "### No nick name detected for these pairs"
                usableName=${sampleNames//,/-}
        else
                echo "### Nick name for these pairs are: $shorterName"
                usableName="$shorterName"
        fi

        for eachSample in ${sampleNames//,/ }
        do
                ((sampleCount++))
                #echo "eachsample: $eachSample"
                sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
                kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
                samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
                assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
        done
        sampleCount=0
        missingSampleCount=0
        sampleList=""
        trackName="$runDir/mpileup/$usableName/$usableName"
        if [ ! -e $trackName.samtoolsMpileUpPass ] ; then
                echo "### mpileup pass is not done yet. Missing: $trackName.samtoolsMpileUpPass"
                ((qsubFails++))
                continue
        fi

        if [[ -e $trackName.mpileup_All.vcf.snpEffInQueue || -e $trackName.mpileup_All.vcf.snpEffPass || -e $trackName.mpileup_All.vcf.snpEffFail ]] ; then
                echo "### snpEff is already done, failed or inQueue"
                continue
        fi
        echo "### Submitting $trackName.mpileup_All.vcf to queue for snpEff..."
        qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$trackName.mpileup_All.vcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
        if [ $? -eq 0 ] ; then
                touch $trackName.mpileup_All.vcf.snpEffInQueue
        else
                ((qsubFails++))
        fi
        sleep 2
done
## This loop if for finding multi variant called vcf files from freebayes
for dnaPairLine in `cat $configFile | grep '^DNAPAIR=\|^DNAFAMI='`
do
        echo "### DNA pair line is $dnaPairLine"
        sampleNames=`echo $dnaPairLine | cut -d= -f2`
        shorterName=`echo $dnaPairLine | cut -d= -f2 | cut -d\; -f2`
        if [[ $shorterName == *,* ]] ; then
                shorterName=""
                echo "### No nick name detected for these pairs"
                usableName=${sampleNames//,/-}
        else
                echo "### Nick name for these pairs are: $shorterName"
                usableName="$shorterName"
        fi

        for eachSample in ${sampleNames//,/ }
        do
                ((sampleCount++))
                #echo "eachsample: $eachSample"
                sampleLine=`cat $configFile | awk '/^SAMPLE=/' | awk 'BEGIN{FS=","} $2=="'"$eachSample"'"'`
                kitName=`echo $sampleLine | cut -d= -f2 | cut -d, -f1`
                samName=`echo $sampleLine | cut -d= -f2 | cut -d, -f2`
                assayID=`echo $sampleLine | cut -d= -f2 | cut -d, -f4`
        done
        sampleCount=0
        missingSampleCount=0
        sampleList=""
        trackName="$runDir/freebayes/$usableName/$usableName"
        if [ ! -e $trackName.freebayesPass ] ; then
                echo "### freebayes pass is not done yet. Missing: $trackName.freebayesPass"
                ((qsubFails++))
                continue
        fi

        if [[ -e $trackName.freebayes_All.vcf.snpEffInQueue || -e $trackName.freebayes_All.vcf.snpEffPass || -e $trackName.freebayes_All.vcf.snpEffFail ]] ; then
                echo "### snpEff is already done, failed or inQueue"
                continue
        fi
        echo "### Submitting $trackName.freebayes_All.vcf to queue for snpEff..."
        qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$trackName.freebayes_All.vcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
        if [ $? -eq 0 ] ; then
                touch $trackName.freebayes_All.vcf.snpEffInQueue
        else
                ((qsubFails++))
        fi
        sleep 2
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
		hcVcf=$runDir/$kitName/$samName/$samName.proj.md.bam.HC_All.vcf
		hcPas=$runDir/$kitName/$samName/$samName.proj.md.bam.hcPass
		mpVcf=$runDir/mpileup/$samName/$samName.mpileup_All.vcf
		mpPas=$runDir/mpileup/$samName/$samName.samtoolsMpileUpPass
		fbVcf=$runDir/freebayes/$samName/$samName.freebayes_All.vcf
		fbPas=$runDir/freebayes/$samName/$samName.freebayesPass
		jrBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
		jrRequested=`cat $configFile | grep '^DNAPAIR=\|^DNAFAMI=' | grep $samName | wc -l`
		if [ $jrRequested -gt 0 ] ; then
			echo "### Single sample snpEff will not operate on a single bam because Joint IR is requested for $samName"
		else
			echo "### Single sample snpeff will run on single bam because Joint IR is NOT requested for $samName"
			if [[ ! -e $hcPas || ! -e $hcVcf ]] ; then
				echo "### Either hcPass or the vcf itself is missing for $hcVcf"
				((qsubFails++))
			else
				if [[ -e $hcVcf.snpEffPass || -e $hcVcf.snpEffInQueue || -e $hcVcf.snpEffFail ]] ; then
					echo "### snpEff already passed, in queue, or failed for $hcVcf"
				else
					echo "### Submitting for hc vcf for snpEff: $hcVcf"
					qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$hcVcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
					if [ $? -eq 0 ] ; then
						touch $hcVcf.snpEffInQueue
					else
						((qsubFails++))
					fi
					sleep 2
				fi
			fi
			if [[ ! -e $mpPas || ! -e $mpVcf ]] ; then
                                echo "### Either mpPass or the vcf itself is missing for $mpVcf"
                                ((qsubFails++))
                        else
                                if [[ -e $mpVcf.snpEffPass || -e $mpVcf.snpEffInQueue || -e $mpVcf.snpEffFail ]] ; then
                                        echo "### snpEff already passed, in queue, or failed for $mpVcf"
                                else
                                        echo "### Submitting for mp vcf for snpEff: $mpVcf"
                                        qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$mpVcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
                                        if [ $? -eq 0 ] ; then
                                                touch $mpVcf.snpEffInQueue
                                        else
                                                ((qsubFails++))
                                        fi
                                        sleep 2
                                fi
                        fi
			if [[ ! -e $fbPas || ! -e $fbVcf ]] ; then
                                echo "### Either fbPass or the vcf itself is missing for $fbVcf"
                                ((qsubFails++))
                        else
                                if [[ -e $fbVcf.snpEffPass || -e $fbVcf.snpEffInQueue || -e $fbVcf.snpEffFail ]] ; then
                                        echo "### snpEff already passed, in queue, or failed for $fbVcf"
                                else
                                        echo "### Submitting for fb vcf for snpEff: $fbVcf"
                                        qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$fbVcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
                                        if [ $? -eq 0 ] ; then
                                                touch $fbVcf.snpEffInQueue
                                        else
                                                ((qsubFails++))
                                        fi
                                        sleep 2
                                fi
                        fi

		fi
	else
		echo "### Assay ID is $assayID. Must be RNA."
		#pcDir=$runDir/$kitName/$samName
		#inBam=$runDir/$kitName/$samName/$samName.proj.bam
		hcVcf=$runDir/rnaHC/$samName/$samName.rnaHC_All.vcf
		hcPas=$runDir/rnaHC/$samName/$samName.RNAhcPass
		jrBam=$runDir/$kitName/$samName/$samName.proj.md.jr.bam
		jrRequested=`cat $configFile | grep '^DNAPAIR=\|^DNAFAMI=' | grep $samName | wc -l`
		if [ $jrRequested -gt 0 ] ; then
			echo "### HC will not operate on a single bam because Joint IR is requested for $samName"
		else
			echo "### HC will run on single bam because Joint IR is NOT requested for $samName"
			if [[ ! -e $hcPas || ! -e $hcVcf ]] ; then
				echo "### Either hcPass or the vcf itself is missing for $hcVcf"
				((qsubFails++))
			else
				if [[ -e $hcVcf.snpEffPass || -e $hcVcf.snpEffInQueue || -e $hcVcf.snpEffFail ]] ; then
					echo "### snpEff already passed, in queue, or failed for $hcVcf"
				else
					echo "### Submitting for hc vcf for snpEff: $hcVcf"
					qsub -A $debit -l nodes=1:ppn=$nCores -v SNPEFFPATH=$snpeffPath,VCF=$hcVcf,DBSNP=$dbsnp,DBVERSION=$snpeffdb,RUNDIR=$runDir,NXT1=$nxtStep1,D=$d $pbsHome/pegasus_snpEff.pbs
					if [ $? -eq 0 ] ; then
						touch $hcVcf.snpEffInQueue
					else
						((qsubFails++))
					fi
					sleep 2
				fi
			fi
		fi

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
