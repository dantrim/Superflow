#pragma once

#include <functional>
#include <vector>

#include "SusyNtuple/SusyNt.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/SusyNtObject.h"
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/AnalysisType.h"
#include "SusyNtuple/Trigger.h"

#include "Superflow/Supersys.h"
#include "Superflow/DataDefinitions.h"
#include "JVFUncertaintyTool/JVFUncertaintyTool.h"


using namespace DataDefinitions;

namespace sflow {

    class Superlink {

    public:
        Superlink();
        ~Superlink();

        bool isMC;
        bool isData;
        
        SusyNtTools* tools;

        Trigger* ntTrig;

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

        JVFUncertaintyTool* jvfTool;

    };
};
