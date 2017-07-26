// Superlink.cxx
//

using namespace std;

#include "Superflow/Superlink.h"


namespace sflow {

    Superlink::Superlink()
    {
        isMC = false;
        isData = false;

        tools = nullptr;

        anaType = Susy::AnalysisType::Ana_2Lep;

        nt = nullptr;
        weights = nullptr;
        nt_sys = Susy::NtSys::NOM;

        preElectrons = nullptr;
        preMuons = nullptr;
        preTaus = nullptr;
        preJets = nullptr;
        baseElectrons = nullptr;
        baseMuons = nullptr;
        baseLeptons = nullptr;
        baseTaus = nullptr;
        baseJets = nullptr;

        electrons = nullptr;
        muons = nullptr;
        leptons = nullptr;
        signalLeptons = nullptr;
        taus = nullptr;
        jets = nullptr;

        met = nullptr;
        trackMet = nullptr;

    }

    Superlink::~Superlink()
    {
        isMC = false;
        isData = false;

        tools = nullptr;

        anaType = Susy::AnalysisType::Ana_2Lep;

        nt = nullptr;
        weights = nullptr;
        nt_sys = Susy::NtSys::NOM;

        preElectrons = nullptr;
        preMuons = nullptr;
        preTaus = nullptr;
        preJets = nullptr;
        baseElectrons = nullptr;
        baseMuons = nullptr;
        baseLeptons = nullptr;
        baseTaus = nullptr;
        baseJets = nullptr;

        electrons = nullptr;
        muons = nullptr;
        leptons = nullptr;
        signalLeptons = nullptr;
        taus = nullptr;
        jets = nullptr;

        met = nullptr;
        trackMet = nullptr;
    }
}
