#!/bin/bash

cd ${3} # startidr

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
source susynt-read/bash/setup_root.sh
source RootCore/scripts/setup.sh
rc find_packages

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

${1} -i ${2} -c -n 500
