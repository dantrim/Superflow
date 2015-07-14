#include "Superflow/SuperTools.h"

#include <iostream>
#include <iomanip>
#include <cstddef>


namespace SuperTools {

std::fstream outMatrix;
std::string debugMatrixName = "debugMatrixMethod.txt";

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// from [1] https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/SUSYPhys/SUSYTools/trunk/Root/SUSYObjDef.cxx#L2833
//
double ttbar_powheg_differentialxsec(double ttbarpt)
{
    // Note: ttbarpt = (top + antitop).Pt() and expects Pt to be in MeV
    double weight = 1.0;
    
    if (ttbarpt/1000. < 40.)
        weight = (1./1.011850 + 1./0.994193)/2.;
    else if (ttbarpt/1000. < 170.)
         weight = (1./1.095920 + 1./1.034480)/2.;
    else if (ttbarpt/1000. < 340)
      weight = (1./1.407280 + 1./1.319110)/2.;
    else
      weight = (1./1.799380 + 1./1.710780)/2.;
 
    return weight;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

} // close SuperTools namespace
