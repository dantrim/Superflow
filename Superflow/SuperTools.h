#pragma once

#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include "DileptonMatrixMethod/DileptonMatrixMethod.h"
#include "SusyNtuple/SusyNtTools.h"
#include "Superflow/Superlink.h"
#include "Superflow/Superflow.h"



namespace SuperTools {


double ttbar_powheg_differentialxsec(double ttbarpt);
double getFakeWeight(sflow::Superlink* sl, std::string fakeRegion, susy::fake::Systematic::Value sys, bool doSyst);
//double getFakeWeight(sflow::Superlink* sl, const std::string fakeRegion, int sys);


}
