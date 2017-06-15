#ifndef SUPERFLOWBASE_H
#define SUPERFLOWBASE_H

//SusyNtuple
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/TauId.h"
#include "SusyNtuple/MCWeighter.h"
#include "SusyNtuple/TriggerTools.h"

//Superflow
#include "Superflow/Cut.h"
#include "Superflow/Superlink.h"
#include "Superflow/Supervar.h"
#include "Superflow/Supersys.h"

//ROOT
#include "TSelector.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TEntryList.h"
#include "TGraphAsymmErrors.h"
#include "TVector2.h"


//std/stl
#include <functional>
#include <map>
#include <string>



namespace sflow {

    /////////////////////////////////////
    // Superflow Run Modes
    /////////////////////////////////////
    enum class SuperflowRunMode {
        nominal,
        nominal_and_weight_syst,
        single_event_syst,
        all_syst,
        data,
        null
    };

    /////////////////////////////////////
    // ATLAS stream enum
    /////////////////////////////////////
    enum class ATLAS_stream {
        Main,
        null
    };

    class SuperflowBase : public SusyNtAna {

        public :
            SuperflowBase();
            ~SuperflowBase(){};

            virtual void setAnaName(string name) { app_name = name + "    "; }
            virtual void setDebug(bool dbg) { m_dbg = dbg; }

            // TSelector methods
            virtual void Begin(TTree *tree);
            virtual void Init(TTree *tree);
            virtual Bool_t Notify();
            virtual void Terminate();

            void initialize_output_files(TString input = "");
            void determine_output_filename(TString input = "");
            void write_output_files();

            const double m_epsilon = 1e-12;

            /////////////////////////////////////////////
            // Operators
            /////////////////////////////////////////////

            // Cut operators
            SuperflowBase& operator<<(CutName cut_);
            SuperflowBase& operator<<(std::function<bool(Superlink*)> cut_);

            // Variable operators
            SuperflowBase& operator<<(NewVar new_var_name);
            SuperflowBase& operator<<(HFTname hft_name);
            SuperflowBase& operator<<(SaveVar save_var);

            SuperflowBase& operator<<(std::function<double(Superlink*, var_float*)> var_);
            SuperflowBase& operator<<(std::function<double(Superlink*, var_double*)> var_);
            SuperflowBase& operator<<(std::function<vector<double>(Superlink*, var_float_array*)> var_);
            SuperflowBase& operator<<(std::function<vector<bool>(Superlink*, var_bool_array*)> var_);
            SuperflowBase& operator<<(std::function<int(Superlink*, var_int*)> var_);
            SuperflowBase& operator<<(std::function<bool(Superlink*, var_bool*)> var_);
            SuperflowBase& operator<<(std::function<void(Superlink*, var_void*)> var_);

            // Systematic operators
            SuperflowBase& operator<<(NewSystematic new_sys);
            SuperflowBase& operator<<(TreeName tree_name);
            SuperflowBase& operator<<(EventSystematic obj_);
            SuperflowBase& operator<<(WeightSystematic obj_);
            SuperflowBase& operator<<(SaveSystematic save_var);


        protected :

            std::string app_name;

            bool m_dbg;

            void initialize_output_arrays();
            void initialize_sys_state();
            void initialize_data_periods();

            SusyNtSys m_singleEventSyst;
            Supersys* m_RunSyst;

            /////////////////////////////////////////
            // input TChain that we loop over
            /////////////////////////////////////////
            TChain* m_input_chain;

            /////////////////////////////////////////
            // output
            /////////////////////////////////////////
            string m_outputFileName;
            string m_outputFileNameSuffix;
            string m_entry_list_FileName;
            string m_tree_name_auto;
            TFile* m_outputFile;
            TFile* m_entryListFile;
            TTree* m_HFT;

            TFile** m_output_array;
            TTree** m_HFT_array;
            
            /////////////////////////////////////////
            // entry lists
            /////////////////////////////////////////
            TEntryList* m_entry_list_total;
            TEntryList* m_entry_list_single_tree;

            /////////////////////////////////////////
            // counters
            /////////////////////////////////////////
            vector<double> m_RawCounter;
            vector<double> m_WeightCounter;
            bool m_countWeights;

            ATLAS_stream m_stream;
            map<ATLAS_stream, string> m_data_stream2string;

            /////////////////////////////////////////
            // flags
            /////////////////////////////////////////
            int m_CutStoreUntitled;
            bool m_CutStore_Name_Exists;
            SuperflowRunMode m_runMode;
            SupersysState m_sysState;
            SupervarState m_varState;

            /////////////////////////////////////////
            // holders for loaded stuff
            /////////////////////////////////////////

            vector<std::function<bool(Superlink*)>> m_CutStore;
            vector<string> m_CutStoreNames;

            Supersys m_sysTemplate;
            vector<Supersys> m_sysStore;
            bool m_sys_hasNiceName;
            bool m_sys_hasTreeName;
            bool m_sys_hasType;
            bool m_sys_hasSystematic;

            vector<SupervarType> m_varType;
            vector<string> m_varNiceName;
            vector<string> m_varHFTName;
            bool m_superVar_hasFunction;
            bool m_superVar_hasNiceName;
            bool m_superVar_hasHFTName;
            int m_superVar_Untitled;

            /////////////////////////////////////////
            // function placeholders
            /////////////////////////////////////////
            std::function<double(Superlink*, var_float*)>               m_nullExprFloat;
            std::function<vector<double>(Superlink*, var_float_array*)> m_nullExprFloatArray;
            std::function<vector<bool>(Superlink*, var_bool_array*)>    m_nullExprBoolArray;
            std::function<double(Superlink*, var_double*)>              m_nullExprDouble;
            std::function<int(Superlink*, var_int*)>                    m_nullExprInt;
            std::function<bool(Superlink*, var_bool*)>                  m_nullExprBool;
            std::function<void(Superlink*, var_void*)>                  m_nullExprVoid;

            /////////////////////////////////////////
            // arrays of filled functions
            /////////////////////////////////////////
            vector<std::function<double(Superlink*, var_float*)>>       m_varExprFloat;
            vector<std::function<vector<double>(Superlink*, var_float_array*)>> m_varExprFloatArray;
            vector<std::function<vector<bool>(Superlink*, var_bool_array*)>> m_varExprBoolArray;
            vector<std::function<double(Superlink*, var_double*)>> m_varExprDouble;
            vector<std::function<int(Superlink*, var_int*)>> m_varExprInt;
            vector<std::function<bool(Superlink*, var_bool*)>> m_varExprBool;
            vector<std::function<void(Superlink*, var_void*)>> m_varExprVoid;
            

            /////////////////////////////////////////
            // arrays of branch values
            /////////////////////////////////////////
            Float_t* m_varFloat;
            vector<vector<double>> m_varFloatArray;
            vector<vector<bool>> m_varBoolArray;
            Double_t* m_varDouble;
            Int_t* m_varInt;
            Bool_t* m_varBool;

            vector<vector<vector<double>>> m_varFloatArray_array;
            vector<vector<vector<bool>>> m_varBoolArray_array;
            Float_t** m_varFloat_array;
            Double_t** m_varDouble_array;
            Int_t** m_varInt_array;
            Bool_t** m_varBool_array; // the arrays are not initialized unless needed
            
            vector<int> index_weight_sys;
            vector<int> index_event_sys;

            int m_tree_leafs_size;
            int m_weight_leaf_offset;


    }; // class SuperflowBase


}; // namespace sflow




#endif
