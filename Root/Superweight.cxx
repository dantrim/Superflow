#include "Superflow/Superweight.h"

namespace sflow {

//////////////////////////////////////////////////////////////////////////////
Superweight::Superweight()
{
    reset();
}
//////////////////////////////////////////////////////////////////////////////
void Superweight::reset()
{
    susynt = 1.0;
    pileup = 1.0;
    lepSf = 1.0;
    trigSf = 1.0;
    btagSf = 1.0;
    jvtSf = 1.0;
}
//////////////////////////////////////////////////////////////////////////////
double Superweight::product()
{
    // dantrim July 2017 -- the user should be smart about using the product method,
    // I personally don't recommend it since the SF are analysis/selection dependent
    #warning need to commit to git adding back the btagSf or fix how it is handled in the Superflow weight syst
    return susynt * pileup * lepSf * jvtSf * btagSf;
}

}; // namespace
