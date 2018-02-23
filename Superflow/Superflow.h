#pragma once

// std
#include <vector>

//Superflow
#include "Superflow/SuperflowBase.h"

namespace sflow {

    class Superflow : public SuperflowBase {

    public:
        Superflow();
        ~Superflow();

        // TSelector Methods
        void Begin(TTree *tree);
        void Init(TTree *tree);
        Bool_t Notify();
        Bool_t Process(Long64_t entry);
        void Terminate();
        
        void setLumi(const float lumi = 1000. /*pb-1*/);
        void setCountWeights(bool value);
        void setRunMode(SuperflowRunMode run_mode_);
        void setSingleEventSyst(SusyNtSys nt_syst_);
        void setChain(TChain* input_chain_);
        void setFileSuffix(string suffix);

    protected:

        void attach_superlink(Superlink* sl_);
        
        MCWeighter* m_mcWeighter;
        bool initialize_mc_weighter(TTree* tree);

        bool computeWeights(Susy::SusyNtObject &ntobj, MCWeighter &weighter, const LeptonVector& leptons,
            const JetVector& jets, Supersys* super_sys, Superweight* weightComponents);

        double computeLeptonEfficiencySf(const Susy::Lepton &lep, const SupersysWeight sys);

        float m_luminosity;

        Superweight* m_weights;


    };

}
