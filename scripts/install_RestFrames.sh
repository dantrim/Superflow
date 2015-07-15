#!/bin bash

# This script will checkout the RestFrames libraries
# and install them to your ${HOME} directory.

cd ${ROOTCOREDIR}
cd ..
startDir=${PWD}
echo "Starting from ${startDir}"
# checkout the master branch
echo " > Checking out RestFrames from GitHub"
git clone -b master git@github.com:crogan/RestFrames.git

echo "cd RestFrames"
cd RestFrames
installDir=${PWD}
echo "Installing RestFrames libraries here: ${installDir}"
echo "./configure --prefix ${installDir}"
./configure --prefix ${installDir}

make &&

make install

cd ${startDir}
source RestFrames/setup_RestFrames.sh
#cd Superflow/cmt/
#
#cxxFlags=$(restframes-config --cxxflags)
#libFlags=$(restframes-config --libs)
#echo ${libFlags}
#inText="PACKAGE"
#outText="PACKAGE\_LIBFLAGS = ${libFlags}"
#sed -i "s%PACKAGE\_LIBFLAGS =%$outText%g" Makefile.RootCore
#outText="PACKAGE\_CXXFLAGS = ${cxxFlags}"
#sed -i "s%PACKAGE\_CXXFLAGS =%$outText%g" Makefile.RootCore
#
#cd ${startDir}

echo ""
tput setaf 2
echo "Note:                                                                                        "
tput setaf 1
echo " > Heed these words:"
tput setaf 2
echo " > The script RestFrames/setup_RestFrames.sh(csh) sets up the necessary environment variables"
echo " > in order to properly run using the RestFrames libraries.                                  "
echo " >"
echo " > You must add the following line to your .bashrc file:                                     "
echo " >       . \$(restframes-config --prefix)/libexec/setup_RestFrames.sh                        "
echo " > Or have something equivalent to call inside of your work area.                             "
echo " >"
echo " > In cmt/Makefile.RootCore for your Superflow package (or package that depends on Superflow) "
echo " > you must add the following to the library and cxx flags lines:                            "
echo " >        PACKAGE_CXXFLAGS =  \$(restframes-config --cxxflags)                               "
echo " > and                                                                                       "
echo " >        PACKAGE_LIBFLAGS = \$(restframes-config --libs)                                    "
echo " > Where you should put the explict, expanded forms of these environment variables, e.g.  copy"
echo " > the output from 'echo \$(restframes-config --cxxflags)' to the PACKAGE_CXXFLAGS line      "
echo " > and the output from 'echo \$(restframes-config --libs)' to the PACKAGE_LIBFLAGS line.     "
tput setaf 1
echo " > Do this NOW!                                                                              "
tput sgr0
echo ""

