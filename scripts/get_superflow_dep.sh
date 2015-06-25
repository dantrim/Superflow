#!/bin/bash

#
# This will check out several of the packages that Superflow
# depends on out of the box.
#

SVNPHYS="svn+ssh://svn.cern.ch/reps/atlasphys/"

qflipURL="$SVNPHYS/Physics/SUSY/Analyses/WeakProduction/ChargeFlip/tags/ChargeFlip-00-00-19-01"

echo ""
echo "Checking out ChargeFlip package from SVN."
echo " > Before use, make sure that the q-flip probability map for your"
echo " > analysis is present (you may need another tag or to make the "
echo " > new map)."
echo ""
svn co $qflipURL ChargeFlip || return || exit

echo ""
echo "Checking out DileptonMatrixMethod from github."
echo " > Before use, make sure that the FakeMatrix (in DileptonMatrixMethod/data)"
echo " > for your analysis is present."
echo ""
git clone git@github.com:gerbaudo/DileptonMatrixMethod.git

echo "Done checking out dependencies of Superflow. Now compile."
