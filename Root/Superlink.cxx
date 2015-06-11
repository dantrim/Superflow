// Superlink.cxx
//

using namespace std;

#include "Superflow/Superlink.h"

namespace sflow {

    Superlink::Superlink()
    {
        isMC = false;
        isData = false;
        doFake = false;

        tools = nullptr;

        anaType = Ana_2Lep;

        nt = nullptr;
        weights = nullptr;
        nt_sys = NtSys_NOM;

        preElectrons = nullptr;
        preMuons = nullptr;
        preJets = nullptr;
        baseElectrons = nullptr;
        baseMuons = nullptr;
        baseLeptons = nullptr;
        baseTaus = nullptr;
        baseJets = nullptr;

        electrons = nullptr;
        muons = nullptr;
        leptons = nullptr;
        taus = nullptr;
        jets = nullptr;
        jets2Lep = nullptr;

        mediumTaus = nullptr;
        tightTaus = nullptr;

        met = nullptr;

        dileptonTrigger = nullptr;

        jvfTool = nullptr;
    }

    Superlink::~Superlink()
    {
        isMC = false;
        isData = false;
        doFake = false;

        tools = nullptr;

        anaType = Ana_2Lep;

        nt = nullptr;
        weights = nullptr;
        nt_sys = NtSys_NOM;

        preElectrons = nullptr;
        preMuons = nullptr;
        preJets = nullptr;
        baseElectrons = nullptr;
        baseMuons = nullptr;
        baseLeptons = nullptr;
        baseTaus = nullptr;
        baseJets = nullptr;

        electrons = nullptr;
        muons = nullptr;
        leptons = nullptr;
        taus = nullptr;
        jets = nullptr;
        jets2Lep = nullptr;

        mediumTaus = nullptr;
        tightTaus = nullptr;

        met = nullptr;

        dileptonTrigger = nullptr;

        jvfTool = nullptr;
    }
}
