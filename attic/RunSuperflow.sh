#!/bin/bash

R_ANA_NAME='wwlikeAna'
R_START_DIR='/gdata/atlas/dantrim/SusyAna/n0211val/'
R_OUTPUT_DIR='/gdata/atlas/dantrim/SusyAna/histoAna/run2early/n0211/mc/Raw/'
R_LOG_DIR='/gdata/atlas/dantrim/SusyAna/histoAna/run2early/n0211/mc/logs/'
#R_OUTPUT_DIR='/gdata/atlas/dantrim/SusyAna/histoAna/run2early/n0211/data/Raw/'
#R_LOG_DIR='/gdata/atlas/dantrim/SusyAna/histoAna/run2early/n0211/data/logs/'

R_WORK_BASE='/scratch/dantrim/'

R_RUN='/gdata/atlas/dantrim/SusyAna/n0211val/Superflow/run/'
R_GEN=${R_RUN}
R_LIST_DIR=${R_GEN}"filelists/"
R_LIST_POST='_n0211.txt'

#R_SAMPLES_LIST=("ttbar wjets ww wz zz singletop zee zmm ztt")
#R_SAMPLES_LIST=("singletop" "ttbar" "wjets_sherpa" "ww_powheg" "wz_powheg" "zz_powheg" "zjets_powheg")
#R_SAMPLES_LIST=("ww_powheg" "ttbar")
R_SAMPLES_LIST=("wjets_sherpa")
#R_SAMPLES_LIST=("bwn")
#R_SAMPLES_LIST=("datalist")


for file_ in ${R_SAMPLES_LIST[@]}; do
	FILE_NAME=${R_LIST_DIR}${R_LIST_PREFIX}${file_}${R_LIST_POST}
	FILE_LINES=`cat $FILE_NAME`
	
	for line in $FILE_LINES ; do
        export S_ANA_NAME=${R_ANA_NAME}
        export S_MODE='a'
		export S_STARTDIR=${R_START_DIR}
		
		sleep 0.1
		RANNUM=$RANDOM$RANDOM$RANDOM
		
		export S_WORKDIR=${R_WORK_BASE}${RANNUM}
		export S_IN_DIRECTORY=${line%?}'/'
	#	export S_IN_DIRECTORY=${line%?}'t'
		export S_SYSTEMATIC='NONE'
		export S_OUTPUT_DIR=${R_OUTPUT_DIR}
		
		RED_line=${line%?}
		lFileName=$(basename $RED_line)
		strip_one=${lFileName#user.*.}
		strip_two=${strip_one%.SusyNt*}
		
		echo $strip_two
		
		sbatch -J 'wjets '${strip_two} -o ${R_LOG_DIR}${lFileName}_slurm-%j.log ${R_RUN}Superflow.sh 
	done
done
