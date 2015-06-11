#include "Superflow/EventFlags.h"

#include <sstream>      // std::ostringstream

using sflow::EventFlags;

//-----------------------------------------
std::string EventFlags::str() const
{
    std::ostringstream oss;
    oss << " grl: " << grl
        << " larErr: " << larErr
        << " tileErr: " << tileErr
        << " ttcVeto: " << ttcVeto
        << " goodVtx: " << goodVtx
        << " tileTrip: " << tileTrip
        << " lAr: " << lAr
        << " badJet: " << badJet
        << " deadRegions: " << deadRegions
        << " badMuon: " << badMuon
        << " cosmicMuon: " << cosmicMuon
        << " hfor: " << hfor
        << " ge2blep: " << ge2blep
        << " eq2blep: " << eq2blep
        << " mllMin: " << mllMin;
    return oss.str();
}
//-----------------------------------------
