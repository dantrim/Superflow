#!/bin/bash

R_ANA_NAME='toyNtAna'
R_START_DIR='/gdata/atlas/dantrim/SusyAna/mt2Tail/'
R_OUTPUT_DIR='/gdata/atlas/dantrim/SusyAna/mt2Tail/'
R_LOG_DIR='/gdata/atlas/dantrim/SusyAna/mt2Tail/'

R_WORK_BASE='/scratch/dantrim/'

R_RUN='/gdata/atlas/dantrim/SusyAna/mt2Tail/Superflow/run/'
R_GEN=${R_RUN}
R_LIST_DIR='/gdata/atlas/dantrim/SusyAna/mt2Tail/filelists/'
R_LIST_POST='.txt'

R_PROXY='/data7/atlas/dantrim/SusyAna/mt2Tail/Superflow/run/proxy/x509up_u674'
R_SAMPLES_LIST=("ttbar")
#R_SAMPLES_LIST=("ttbar wjets ww wz zz singletop zee zmm ztt")
#R_SAMPLES_LIST=("singletop" "ttbar" "wjets_sherpa" "ww_powheg" "wz_powheg" "zz_powheg" "zjets_powheg")
#R_SAMPLES_LIST=("ww_powheg" "ttbar")
#R_SAMPLES_LIST=("wjets_sherpa")
#R_SAMPLES_LIST=("bwn")
#R_SAMPLES_LIST=("datalist")

for file_ in ${R_SAMPLES_LIST[@]}; do
    for TXT_FILE in ${R_LIST_DIR}${file_}/*txt; do
        export S_ANA_NAME=${R_ANA_NAME}
        export S_STARTDIR=${R_START_DIR}
        export S_PROXY=${R_PROXY}
        export X509_USER_PROXY=${R_PROXY}
        sleep 0.1
		RANNUM=$RANDOM$RANDOM$RANDOM
        export S_WORKDIR=${R_WORK_BASE}${RANNUM}
        export S_IN_FILELIST=${TXT_FILE}

        export S_OUTPUT_DIR=${R_OUTPUT_DIR}
        RED_line=${TXT_FILE%?}
        lfilelistname=$(basename $RED_line)
        strip_one=${lfilelistname#user.*.}
        strip_two=${strip_one%.SusyNt*}
        echo ${strip_two}
        sbatch -J 'test ttbar' -o ${R_LOG_DIR}${lfilelistname}_slurm-%j.log ${R_RUN}SuperflowFAX.sh
    done
done        



#for file_ in ${R_SAMPLES_LIST[@]}; do
#	FILE_NAME=${R_LIST_DIR}${R_LIST_PREFIX}${file_}${R_LIST_POST}
#	FILE_LINES=`cat $FILE_NAME`
#	
#	for line in $FILE_LINES ; do
#        export S_ANA_NAME=${R_ANA_NAME}
#        export S_MODE='a'
#		export S_STARTDIR=${R_START_DIR}
#		
#		sleep 0.1
#		RANNUM=$RANDOM$RANDOM$RANDOM
#		
#		export S_WORKDIR=${R_WORK_BASE}${RANNUM}
#		export S_IN_DIRECTORY=${line%?}'/'
#	#	export S_IN_DIRECTORY=${line%?}'t'
#		export S_SYSTEMATIC='NONE'
#		export S_OUTPUT_DIR=${R_OUTPUT_DIR}
#		
#		RED_line=${line%?}
#		lFileName=$(basename $RED_line)
#		strip_one=${lFileName#user.*.}
#		strip_two=${strip_one%.SusyNt*}
#		
#		echo $strip_two
#		
#		sbatch -J 'wjets '${strip_two} -o ${R_LOG_DIR}${lFileName}_slurm-%j.log ${R_RUN}Superflow.sh 
#	done
#done
