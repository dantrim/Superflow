#pragma once

#include <vector>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TEntryList.h"
#include "TGraphAsymmErrors.h"

#include "TVector2.h"

#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/SusyNtTools.h"

#include "SusyNtuple/MCWeighter.h"
#include "SusyNtuple/DilTrigLogic.h"

#include "Superflow/DataDefinitions.h"

#include "Superflow/Cut.h"
#include "Superflow/Superlink.h"
#include "Superflow/Supervar.h"
#include "Superflow/Supersys.h"

#include "DileptonMatrixMethod/DileptonMatrixMethod.h"
#include "ChargeFlip/chargeFlip.h"

using namespace DataDefinitions;

namespace sflow {

    enum class SuperflowRunMode {
        nominal,
        nominal_and_weight_syst,
        single_event_syst,
        all_syst,
        data,
        fakes,
        null
    };

    enum class ATLAS_period { A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z, null };
    enum class ATLAS_stream { Main, Egamma, Muons, null };
    //enum class ATLAS_stream { Main, null };

    class Superflow : public SusyNtAna {

    public:
        Superflow();
        ~Superflow();

        // Cut Operators
        Superflow& operator<<(CutName cut_);
        Superflow& operator<<(std::function<bool(Superlink*)> cut_);

        // Var Operators
        Superflow& operator<<(NewVar new_var_name);
        Superflow& operator<<(HFTname hft_name);
        Superflow& operator<<(SaveVar save_var);

        Superflow& operator<<(std::function<double(Superlink*, var_float*)> var_);
        Superflow& operator<<(std::function<double(Superlink*, var_double*)> var_);
        Superflow& operator<<(std::function<int(Superlink*, var_int*)> var_);
        Superflow& operator<<(std::function<bool(Superlink*, var_bool*)> var_);
        Superflow& operator<<(std::function<void(Superlink*, var_void*)> var_);

        // Systematics Operators
        Superflow& operator<<(NewSystematic new_sys);
        Superflow& operator<<(TreeName tree_name);
        Superflow& operator<<(EventSystematic obj_);
        Superflow& operator<<(WeightSystematic obj_);
        Superflow& operator<<(SaveSystematic save_var);

        // TSelector Methods
        void Begin(TTree *tree);        ///< called before looping on entries
        void Init(TTree *tree);         ///< called when the TChain is attached
        Bool_t Notify();                ///< called at each event
        Bool_t Process(Long64_t entry); ///< called at each event
        void Terminate();               ///< called after looping is finished
        
        void setAnaName(string name) { app_name = name + "    "; }
        void setLumi(const float lumi = LUMI_A_A3); ///< Set the MC normalization lumi
        void setCountWeights(bool value); ///< Toggle the display of the weighted cuts. (default off)
        void setRunMode(SuperflowRunMode run_mode_);
        void setSingleEventSyst(SusyNtSys nt_syst_);
        void setChain(TChain* input_chain_);
        void setFakeRegion(string fk_reg);
        void setFake2dParam(bool use2d);
        void setQFlip();
        float computeChargeFlipWeight(const LeptonVector &leptons, const SupersysWeight sys);
        bool isGenuineSS(const LeptonVector& leptons);
        bool hasQFlip(const LeptonVector& leptons); 
        bool isMM(const LeptonVector& leptons);

    protected:
        string app_name;
        void attach_superlink(Superlink* sl_);

        DilTrigLogic* m_trigObj; ///< trigger logic class
        MCWeighter* m_mcWeighter; ///< tool to determine the normalization

        bool m_do_qflip;
        chargeFlip*     m_chargeFlip; ///< chargeFlip tool

        // Matrix Method tool and configurables
        susy::fake::DileptonMatrixMethod* m_matrix;
        std::string m_matrixFilename;
        bool m_use2dparametrization;
        bool m_allconfigured;
        bool initMatrixTool();
        bool selectBaseLineLeptons;

        bool computeWeights(
            Susy::SusyNtObject &ntobj,
            MCWeighter &weighter,
            const LeptonVector& leptons,
            const JetVector& jets,
            Supersys* super_sys,
            Superweight* weightComponents);

        //double computeDileptonTriggerWeight(const LeptonVector &leptons, const SusyNtSys sys);

        double computeBtagWeight(const JetVector& jets, const Susy::Event* evt, SupersysWeight sys);

        double computeLeptonEfficiencySf(const Susy::Lepton &lep, const SupersysWeight sys);

        double get_isr_weights(sflow::Superlink* sl, const SupersysWeight sys);

        vector<double> m_RawCounter;
        vector<double> m_WeightCounter; // indexed by cut #

        bool m_countWeights;
        float m_luminosity;

        Superweight* m_weights;

        string m_outputFileName;
        string m_entry_list_FileName;
        string m_tree_name_auto;
        TFile* m_outputFile;
        TFile* m_entryListFile;
        TTree* m_HFT;

        TFile** m_output_array;
        TTree** m_HFT_array;

        TEntryList* m_entry_list_total;
        TEntryList* m_entry_list_single_tree;

        ATLAS_period m_period;
        ATLAS_stream m_stream;

        vector<std::function<bool(Superlink*)>> m_CutStore;
        vector<string> m_CutStoreNames;
        int m_CutStoreUntitled;
        bool m_CutStore_Name_Exists;

        std::function<double(Superlink*, var_float*)> m_nullExprFloat;
        std::function<double(Superlink*, var_double*)> m_nullExprDouble;
        std::function<int(Superlink*, var_int*)> m_nullExprInt;
        std::function<bool(Superlink*, var_bool*)> m_nullExprBool;
        std::function<void(Superlink*, var_void*)> m_nullExprVoid;

        vector<std::function<double(Superlink*, var_float*)>> m_varExprFloat;
        vector<std::function<double(Superlink*, var_double*)>> m_varExprDouble;
        vector<std::function<int(Superlink*, var_int*)>> m_varExprInt;
        vector<std::function<bool(Superlink*, var_bool*)>> m_varExprBool;
        vector<std::function<void(Superlink*, var_void*)>> m_varExprVoid;

        Float_t* m_varFloat;
        Double_t* m_varDouble;
        Int_t* m_varInt;
        Bool_t* m_varBool;

        Float_t** m_varFloat_array;
        Double_t** m_varDouble_array;
        Int_t** m_varInt_array;
        Bool_t** m_varBool_array; // the arrays are not initialized unless needed

        SupervarState m_varState;

        vector<SupervarType> m_varType;
        vector<string> m_varNiceName;
        vector<string> m_varHFTName;
        bool m_superVar_hasFunction;
        bool m_superVar_hasNiceName;
        bool m_superVar_hasHFTName;
        int m_superVar_Untitled;

        SuperflowRunMode m_runMode;
        SupersysState m_sysState;
        string m_fake_region;

        Supersys m_sysTemplate;
        vector<Supersys> m_sysStore;
        bool m_sys_hasNiceName;
        bool m_sys_hasTreeName;
        bool m_sys_hasType;
        bool m_sys_hasSystematic;

        SusyNtSys m_singleEventSyst;
        Supersys* m_RunSyst;

        vector<int> index_weight_sys;
        vector<int> index_event_sys;

        int m_tree_leafs_size;
        int m_weight_leaf_offset;

    private:
        bool initMcWeighter(TTree *tree);

        map<ATLAS_stream, map<ATLAS_period, string>> m_data_periods;
        map<ATLAS_stream, string> m_data_stream2string;
        map<SusyNtSys, string> m_NtSys_to_string;

        TChain* m_input_chain;

        const double m_epsilon = 1e-12;
    };

};
