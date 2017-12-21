#!/usr/bin/env bash
#SBATCH --job-name="pegasus_thFPost"
#SBATCH --time=0-96:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL

cd ${DIR}

time=`date +%d-%m-%Y-%H-%M`
beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### REF: ${REF}"
echo "### DIR: ${DIR}"
echo "### INDEXBASE: ${INDEXBASE}"
echo "### 2VCFPATH: ${THFUSION2VCFPATH}"

module load bowtie/0.12.9

PATH=/home/tgenref/binaries/blast/ncbi-blast-2.2.29+/bin:$PATH
export PATH

newName=`basename ${DIR}`
newName=${newName/.topHatFusionDir}
echo "TIME:$time starting tophat fusion post on ${DIR} with indexbase of ${INDEXBASE}"
echo "${TOPHAT2PATH}/tophat-fusion-post \
    -p 4 \
    --num-fusion-reads 3 \
    --num-fusion-pairs 2 \
    --num-fusion-both 5 \
    --skip-read-dist \
    --fusion-read-mismatches 3 \
    ${INDEXBASE} > ${DIR}.thFPostOut
"

${TOPHAT2PATH}/tophat-fusion-post \
    -p 4 \
    --num-fusion-reads 3 \
    --num-fusion-pairs 2 \
    --num-fusion-both 5 \
    --skip-read-dist \
    --fusion-read-mismatches 3 \
    ${INDEXBASE} > ${DIR}.thFPostOut

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

    cd ${DIR}/tophatfusion_out/
    ${THFUSION2VCFPATH}/tophatFusion2vcf.sh ${DIR}/tophatfusion_out/$newName.thFusion.result.txt $newName ${REF}
    cd -

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
