#!/usr/bin/env bash
#SBATCH --job-name="pegasus_sleuth"
#SBATCH --time=0-60:00:00
#SBATCH --mail-user=jetstream@tgen.org
#SBATCH --mail-type=FAIL

module load R/3.2.1

beginTime=`date +%s`
machine=`hostname`
echo "### NODE: $machine"
echo "### NXT1: ${NXT1}"
echo "### RUNDIR: ${RUNDIR}"
echo "### SLEUTHPATH: ${SLEUTHPATH}"
echo "### SLEUTHCONFIG: ${SLEUTHCONFIG}"
echo "### SLEUTHOUTDIR: ${SLEUTHOUTDIR}"
echo "### OBJECTDATA: ${OBJECTDATA}"
echo "### SLEUTHOUTFILE: ${SLEUTHOUTFILE}"
echo "### KALLISTOOUT: ${KALLISTOOUT}"

base_Dir=${KALLISTOOUT}
configFile=${SLEUTHCONFIG}
outData=${SLEUTHOUTFILE}
objectData=${OBJECTDATA}

Rscript --vanilla ${SLEUTHPATH} \
    $base_Dir \
    $configFile \
    $outData \
    $objectData

if [ $? -eq 0 ] ; then
    touch ${SLEUTHOUTDIR}.sleuthPass
else
    touch ${SLEUTHOUTDIR}.sleuthFail
fi

rm ${SLEUTHOUTDIR}.sleuthInQueue

endTime=`date +%s`
elapsed=$(( $endTime - $beginTime ))
(( hours=$elapsed/3600 ))
(( mins=$elapsed%3600/60 ))
echo "RUNTIME:SLEUTH:$hours:$mins" > ${SLEUTHOUTDIR}.sleuth.totalTime
echo "### Ending Sleuth."
