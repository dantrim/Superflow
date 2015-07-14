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

        ntTrig = nullptr;

        anaType = Susy::AnalysisType::Ana_2Lep;

        nt = nullptr;
        weights = nullptr;
        nt_sys = Susy::NtSys::NOM;

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


      //  jvfTool = nullptr;
    }

    Superlink::~Superlink()
    {
        isMC = false;
        isData = false;

        tools = nullptr;

        ntTrig = nullptr;

        anaType = Susy::AnalysisType::Ana_2Lep;

        nt = nullptr;
        weights = nullptr;
        nt_sys = Susy::NtSys::NOM;

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


     //   jvfTool = nullptr;
    }
}
