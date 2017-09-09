#!/usr/bin/env bash
#PBS -S /bin/bash
#SBATCH --job-name="pegasus_thFPost"
#SBATCH --time=0-96:00:00
#SBATCH --mail-user=tgenjetstream@tgen.org
#SBATCH --mail-type=FAIL
#PBS -j oe
#SBATCH --output="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
#SBATCH --error="/${D}/oeFiles/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"

cd ${DIR}

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### DIR: ${DIR}"
echo "### INDEXBASE: ${INDEXBASE}"
echo "### 2VCFPATH: ${THFUSION2VCFPATH}"

newName=`basename ${DIR}`
newName=${newName/.topHatFusionDir}
echo "TIME:$time starting tophat fusion post on ${DIR} with indexbase of ${INDEXBASE}"
#tophat-fusion-post -p 8 --num-fusion-reads 1 --num-fusion-pairs 2 --num-fusion-both 5 ~/references/bowtie/Homo_sapiens.GRCh37.6*
#tophat-fusion-post -p 8 --num-fusion-reads 1 --num-fusion-pairs 2 --num-fusion-both 5 ${BOWTIE1_INDEX} 
perf stat ${TOPHAT2PATH}/tophat-fusion-post -p 16 --num-fusion-reads 3 --num-fusion-pairs 2 --num-fusion-both 5 --skip-read-dist --fusion-read-mismatches 3 ${INDEXBASE} > ${DIR}.thFPostOut 2> ${DIR}.thFPost.perfOut
if [ $? -eq 0 ] ; then
	echo "success."
	echo "renaming..."
	mv ${DIR}/tophatfusion_out/fusion_seq.fa ${DIR}/tophatfusion_out/$newName.fusion_seq.fa
	mv ${DIR}/tophatfusion_out/fusion_seq.bwtout ${DIR}/tophatfusion_out/$newName.thFusion.fusion_seq.bwtout
	mv ${DIR}/tophatfusion_out/fusion_seq.map ${DIR}/tophatfusion_out/$newName.thFusion.fusion_seq.map
	mv ${DIR}/tophatfusion_out/potential_fusion.txt ${DIR}/tophatfusion_out/$newName.thFusion.potential_fusion.txt
	mv ${DIR}/tophatfusion_out/result.txt ${DIR}/tophatfusion_out/$newName.thFusion.result.txt
	mv ${DIR}/tophatfusion_out/result.html ${DIR}/tophatfusion_out/$newName.thFusion.result.html
	echo "renaming done."
	#start Legendre's script here
	cd ${DIR}/tophatfusion_out/
	${THFUSION2VCFPATH}/tophatFusion2vcf.sh ${DIR}/tophatfusion_out/$newName.thFusion.result.txt $newName ${REF}
	cd -
	#end Legendre's script here
	mv ${DIR}.thFPostOut ${DIR}.thFPostPass	
else
	mv ${DIR}.thFPostOut ${DIR}.thFPostFail
fi
rm -f ${DIR}.thFPostInQueue
endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:THFPOST:$hours:$mins" > ${DIR}.thFPost.totalTime
time=`date +%d-%m-%Y-%H-%M`
echo "TIME:$time tophat fusion post finished"
