#!/bin/bash
#SBATCH -p atlas_all
# #SBATCH --distribution=cyclic
###SBATCH --exclude=c-12-19,c-12-23,c-12-35,c-12-15
##SBATCH --exclude=c-12-35,c-12-39,c-12-15,c-12-23,c-12-27,c-12-31
#SBATCH -N 1 -n 2
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=8:00:00

#echo 'S_ANA_NAME     : '${S_ANA_NAME}
#echo 'S_MODE         : '${S_MODE}
#echo 'S_IN_DIRECTORY : '${S_IN_DIRECTORY}  
#echo 'S_WORKDIR      : '${S_WORKDIR}  
#echo 'S_STARTDIR     : '${S_STARTDIR}  
#echo 'S_SYSTEMATIC   : '${S_SYSTEMATIC}  
#echo 'S_OUTPUT_DIR   : '${S_OUTPUT_DIR}

echo 'S_ANA_NAME     : '${S_ANA_NAME}
echo 'S_STARTDIR     : '${S_STARTDIR}
echo 'S_WORKDIR      : '${S_WORKDIR}
echo 'S_OUTPUT_DIR   : '${S_OUTPUT_DIR}
echo 'S_IN_FILELIST  : '${S_IN_FILELIST}
echo 'S_PROXY        : '${S_PROXY}
echo 'CERT_DIR       : '${X509_CERT_DIR}
echo 'USER_PROX      : '${X509_USER_PROXY}
echo 'MY_USER_PROX   : '${S_PROXY}
echo 'VOMS-DIR       : '${X509_VOMS_DIR}
echo 'VOMSES         : '${X509_VOMSES}

echo ""
echo ""
echo ""


echo ""
echo "Starting on `hostname`, `date`"
echo ""

RANNUM=$RANDOM
SCRATCH=/scratch/${USER}/${SLURM_JOB_ID}_$RANNUM
sleep 1
mkdir -p ${SCRATCH}
if [ ! -d "${SCRATCH}" ]; then
    SCRATCH=/scratch/${USER}/${SLURM_JOB_ID}_$RANNUM$RANNUM
    sleep 1
fi
cd       ${SCRATCH}
echo "Working from ${PWD}"

pwd
echo ""

#export PATH=$PATH:/gdata/atlas/dantrim/SusyAna/Super/Superflow/bin/ 
echo "setting up atlas local root base"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh


echo "setting up fax"
localSetupFAX --skipConfirm
#
echo "calling voms-proxy-init -voms atlas"
voms-proxy-init -voms atlas --noregen


echo "calling: ${S_ANA_NAME} -i ${S_IN_FILELIST} -n 20"
${S_ANA_NAME} -i ${S_IN_FILELIST} -n 20
#${S_ANA_NAME} -i ${S_IN_DIRECTORY} -n 20 # -${S_MODE}  

echo ""

pwd
pwd
pwd

echo ""
echo ""

sleep 1
SUFFIX=$RANDOM$RANDOM
echo "${PWD} contents"
ls -ltrh
#cd ${SCRATCH}
#rm *entrylist*  #remove entry list (to make the next step easier)
#
## check whether or not the destination file(name) exists
## if so, rename with a new random number. TODO: make this a WHILE LOOP
#for filename in *.root; do
#    if [ -f "${S_OUTPUT_DIR}${SUFFIX}_$filename" ]; then
#        sleep 1
#        SUFFIX2=$RANDOM$RANDOM
#        echo Output file exists. Renaming current file.
#        echo From ${SUFFIX}_$filename to ${SUFFIX2}_$filename
#        mv "$filename" "${SUFFIX2}_$filename"
#    else
#        mv "$filename" "${SUFFIX}_$filename"
#    fi;
#done
#for filename in *.root; do mv "$filename" "${SUFFIX}_$filename"; done;

echo "Contents of"
echo "  ${PWD}"
ls -ltrh
cp -p ${SCRATCH}/*.root ${S_OUTPUT_DIR}

echo ${PWD}
echo "Done, `date`"
echo "Cleaning up ${SCRATCH}"
echo ${PWD}
date
date
date
rm -rf $SCRATCH || exit $?
