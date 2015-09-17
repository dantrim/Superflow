#pragma once

#include <functional>
#include <vector>

#include "SusyNtuple/SusyNt.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/SusyNtObject.h"
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/AnalysisType.h"
#include "SusyNtuple/TriggerTools.h"

#include "Superflow/Supersys.h"
#include "Superflow/DataDefinitions.h"


using namespace DataDefinitions;

namespace sflow {

    class Superlink {

    public:
        Superlink();
        ~Superlink();

        bool isMC;
        bool isData;
        
        SusyNtTools* tools;

        AnalysisType anaType;

        Susy::SusyNtObject* nt;
        Superweight*        weights;
        SusyNtSys           nt_sys;

        ElectronVector* preElectrons;   ///< pre electrons prior to baseline selection and OR 
        MuonVector*     preMuons;       ///< pre muons prior to baseline selection and OR
        TauVector*      preTaus;        ///< pre taus prior to baseline selection and OR
        JetVector*      preJets;        ///< pre jets prior to baseline selection and OR

        ElectronVector* baseElectrons;  ///< baseline electrons with OR
        MuonVector*     baseMuons;      ///< baseline muons with OR
        LeptonVector*   baseLeptons;    ///< baseline leptons with OR
        TauVector*      baseTaus;       ///< baseline taus with OR
        JetVector*      baseJets;       ///< baseline jets with OR

        ElectronVector* electrons;      ///< signal electrons
        MuonVector*     muons;          ///< signal muons
        LeptonVector*   leptons;        ///< signal leptons
        TauVector*      taus;           ///< signal taus
        JetVector*      jets;           ///< signal jets

        const Susy::Met* met;           ///< Met
        const Susy::TrackMet* trackMet; ///< TrackMet


    };
};
