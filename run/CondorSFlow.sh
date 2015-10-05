#!/bin/bash

cd ${3} # startidr
echo "Starting at ${PWD}"
echo "Starting at ${PWD}"
echo "Starting at ${PWD}"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
#echo "source susynt-read/bash/setup_root.sh"
#source susynt-read/bash/setup_root.sh
lsetup root
echo "source RootCore/scripts/setup.sh"
source RootCore/scripts/setup.sh
echo "rc find_packages"
rc find_packages
##echo "rc compile"
##rc compile
echo "source RestFrames/setup_RestFrames.sh"
source RestFrames/setup_RestFrames.sh

echo "cd ${4}"
cd ${4} # outdir

echo "Working from ${PWD}"
echo "Working from ${PWD}"
echo "Working from ${PWD}"
echo "input ${2}"
echo "input ${2}"
echo "input ${2}"
echo "executable ${1}"
echo "executable ${1}"
echo "executable ${1}"

echo "${1} -i ${2} -c"
${1} -i ${2} -c
