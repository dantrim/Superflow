#pragma once

#include <functional>
#include <vector>

#include "SusyNtuple/SusyNt.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/SusyNtObject.h"
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/DilTrigLogic.h"

#include "Superflow/Supersys.h"
#include "Superflow/DataDefinitions.h"
#include "JVFUncertaintyTool/JVFUncertaintyTool.h"

#include "DileptonMatrixMethod/DileptonMatrixMethod.h"

using namespace DataDefinitions;

namespace sflow {

    class Superlink {

    public:
        Superlink();
        ~Superlink();

        bool isMC;
        bool isData;
        bool doFake;

        SusyNtTools* tools;

        AnalysisType anaType;

        Susy::SusyNtObject* nt;
        Superweight* weights;
        SusyNtSys nt_sys;

        ElectronVector* preElectrons; ///< selected electrons before OR
        MuonVector* preMuons; ///< selected muons before OR
        JetVector* preJets; ///< selected jets before OR

        ElectronVector* baseElectrons; ///< baseline electrons
        MuonVector* baseMuons; ///< baseline muons
        LeptonVector* baseLeptons; ///< baseline leptons
        TauVector* baseTaus; ///< baseline taus
        JetVector* baseJets; ///< baseline jets

        ElectronVector* electrons; ///< signal electrons
        MuonVector* muons; ///< signal muons
        LeptonVector* leptons; ///< signal leptons
        TauVector* taus; ///< signal taus
        JetVector* jets; ///< signal jets
        JetVector* jets2Lep; ///< signal jets for 2 Lep

        // New organization of tau selections
        TauVector* mediumTaus; ///< taus with medium ID
        TauVector* tightTaus; ///< taus with tight ID

        const Susy::Met* met; ///< Met

        DilTrigLogic* dileptonTrigger;

        JVFUncertaintyTool* jvfTool;

        susy::fake::DileptonMatrixMethod* fakeMatrix;
    };
};
