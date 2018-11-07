#include "Superflow/SuperflowBase.h"


namespace sflow {

///////////////////////////////////////////////////////////////////////////////
SuperflowBase::SuperflowBase()
{
    m_dbg = false;

    app_name = "Superflow    ";

    m_runMode = SuperflowRunMode::null;

    m_CutStore_Name_Exists = false;
    m_CutStoreUntitled     = 1;

    m_countWeights = true;

    m_outputFileName        = "";
    m_outputFileNameSuffix  = "";
    m_tree_name_auto        = "";
    m_outputFile            = nullptr;
    m_HFT                   = nullptr;

    m_HFT_array             = nullptr;

    m_output_array          = nullptr;

    //m_entry_list_total       = nullptr;
    m_entry_list_single_tree = nullptr;

    m_stream = ATLAS_stream::null;

    m_varState = SupervarState::closed;
    m_superVar_hasFunction = false;
    m_superVar_hasNiceName = false;
    m_superVar_hasHFTName  = false;
    m_superVar_Untitled    = 1;

    initialize_output_arrays();

    initialize_sys_state();

    m_tree_leafs_size    = 0;
    m_weight_leaf_offset = 0;

    m_input_chain = nullptr;

    initialize_data_periods();

}

///////////////////////////////////////////////////////////////////////////////
void SuperflowBase::initialize_output_arrays()
{
    m_nullExprFloat  = [](Superlink* sl, var_float*) -> double { return 0.0; };
    m_nullExprFloatArray = [](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> null;
            for(int i = 0; i < 25; i++) {
                null.push_back(0.0);
            }
            return null;
    };
    m_nullExprDouble = [](Superlink* sl, var_double*) -> double { return 0.0; };
    m_nullExprInt    = [](Superlink* sl, var_int*) -> int { return 0; };
    m_nullExprBool   = [](Superlink* sl, var_bool*) -> bool { return false; };
    m_nullExprVoid   = [](Superlink* sl, var_void*) {};

    m_varFloat  = nullptr;
    m_varDouble = nullptr;
    m_varInt    = nullptr;
    m_varBool   = nullptr;

    m_varFloatArray_array.clear();
    m_varBoolArray_array.clear();
    m_varFloat_array  = nullptr;
    m_varDouble_array = nullptr;
    m_varInt_array    = nullptr;
    m_varBool_array   = nullptr;
}
///////////////////////////////////////////////////////////////////////////////
void SuperflowBase::initialize_sys_state()
{
    m_sysState = SupersysState::closed;

    m_sysTemplate.reset();
    m_sys_hasNiceName   = false;
    m_sys_hasTreeName   = false;
    m_sys_hasType       = false;
    m_sys_hasSystematic = false;

    m_singleEventSyst = NtSys::NOM;
    m_RunSyst         = new Supersys(SupersysType::central);

}
///////////////////////////////////////////////////////////////////////////////
void SuperflowBase::initialize_data_periods()
{
    m_data_stream2string[ATLAS_stream::Main] = "physics_Main";

}
///////////////////////////////////////////////////////////////////////////////
void SuperflowBase::Begin(TTree* /*tree*/)
{
    // determine number of weight systematics
    if (m_runMode == SuperflowRunMode::nominal_and_weight_syst || m_runMode == SuperflowRunMode::all_syst) {
        for (int i = 0; i < m_sysStore.size(); i++) {
            if (m_sysStore[i].type == SupersysType::weight) {
                index_weight_sys.push_back(i);
            }
        }
    }

    // determine number of event systematics
    if (m_runMode == SuperflowRunMode::all_syst) {
        for (int i = 0; i < m_sysStore.size(); i++) {
            if (m_sysStore[i].type == SupersysType::event) {
                index_event_sys.push_back(i);
            }
        }
    }
    
    // if run mode single_event_syst, set which syst to run on
    if (m_runMode == SuperflowRunMode::single_event_syst) {
        for (int i = 0; i < m_sysStore.size(); i++) {
            if (m_sysStore[i].type == SupersysType::event && m_sysStore[i].event_syst == m_singleEventSyst) {
                // initially set to nominal
                delete m_RunSyst;

                m_RunSyst = &m_sysStore[i]; // don't delete!!
            }
        }
    } // end if run mode single_event_syst
    

}
///////////////////////////////////////////////////////////////////////////////
void SuperflowBase::Init(TTree* tree)
{

}
///////////////////////////////////////////////////////////////////////////////
//Bool_t SuperflowBase::Notify()
//{
//    cout << app_name << "Superflow::Notify" << endl;
//    static int tree_counter;
//
//    //if (m_entry_list_single_tree != nullptr) m_entry_list_total->Add(m_entry_list_single_tree);
//    //delete m_entry_list_single_tree;
//
//    string new_list_name = m_tree_name_auto + "_" + to_string(tree_counter);
//    tree_counter++;
//
//    m_entry_list_single_tree = new TEntryList();
//    m_entry_list_single_tree->SetTree(m_input_chain->GetTree());
//
//    return kTRUE;
//
//}
///////////////////////////////////////////////////////////////////////////////
void SuperflowBase::determine_output_filename(TString input_sample)
{
    // output file name
    stringstream sfile_name_;

    // determine output file name
    if (m_runMode == SuperflowRunMode::data) {
        m_countWeights = true;

        // output file name
        stringstream sfile_name_;
        sfile_name_ << "CENTRAL_";

        bool data15 = input_sample.Contains("data15_13TeV");
        bool data16 = input_sample.Contains("data16_13TeV");
        bool data17 = input_sample.Contains("data17_13TeV");
        bool is_data = (data15 || data16 || data17);

        if(is_data) {
            TString data_run; data_run.Form("%d",nt.evt()->run);
            string year = (data15 ? "2015" : data16 ? "2016" : data17 ? "2017" : "");
            if(year=="") {
                cout << app_name << "SuperflowBase::determine_output_filename    ERROR Could not determine data year from sample input container: " << endl; 
                cout << app_name << "SuperflowBase::determine_output_filename    \t" << input_sample << endl;
                cout << app_name << "SuperflowBase::determine_output_filename    (seen as: 2015? " << data15 << ", 2016? " << data16 << ", 2017? " << data17 << ")" << endl;
                exit(1);
            }
            if(input_sample.Contains("physics_Main")) {
                cout << app_name << "SuperflowBase::determine_output_filename    ====== Determined input data specifications ====== " << endl;
                cout << app_name << "SuperflowBase::determine_output_filename     stream: Main " << endl;
                cout << app_name << "SuperflowBase::determine_output_filename     run   : " << data_run << endl;
                cout << app_name << "SuperflowBase::determine_output_filename    ==================================================" << endl;
                cout << endl;
                m_stream = ATLAS_stream::Main;
            } // if main
            else {
                cout << app_name << "SuperflowBase::determine_output_filename    ERROR    Could not determine data stream for sample from the input container: " << endl;
                cout << app_name << "SuperflowBase::determine_output_filename    ERROR    \t" << input_sample << endl;
                cout << app_name << "SuperflowBase::determine_output_filename    ERROR    The only supported stream is 'Main' ('physics_Main')." << endl;
                cout << app_name << "SuperflowBase::determine_output_filename    ERROR    >>> Exiting." << endl;
                exit(1);
            } 

            stringstream suffix;
            if(m_outputFileNameSuffix!="")
                suffix << "_" << m_outputFileNameSuffix;
            else suffix << "";

            sfile_name_ << m_data_stream2string[m_stream] << "_" << year << "_" << data_run << suffix.str() << ".root";
            cout << app_name << "SuperflowBase::determine_output_filename    Setting output file name to: " << sfile_name_.str() << endl;
            
            m_outputFileName = sfile_name_.str();

        }
        else {
            cout << app_name << "SuperflowBase::determine_output_filename    ERROR    The input container name does not appear to be a data sample." << endl;
            cout << app_name << "SuperflowBase::determine_output_filename    ERROR    It does not contain either the 'data15_13TeV', 'data16_13TeV', or 'data17_13TeV' grouping and the run mode is" << endl;
            cout << app_name << "SuperflowBase::determine_output_filename    ERROR    set for data (SuperflowRunMode::data). " << endl;
            cout << app_name << "SuperflowBase::determine_output_filename    ERROR    >>> Exiting." << endl;
            exit(1);
        } 

    }
    else if (m_runMode == SuperflowRunMode::single_event_syst) {
        if (Susy::NtSys::SusyNtSysNames.count(m_singleEventSyst) != 0) {

            stringstream suffix;
            if(m_outputFileNameSuffix!="")
                suffix << "_" << m_outputFileNameSuffix;
            else suffix << "";

            stringstream sfile_name_; // output file name
            sfile_name_ << Susy::NtSys::SusyNtSysNames.at(m_singleEventSyst) << "_" << nt.evt()->mcChannel
                        << suffix.str() << ".root";
            cout << app_name << "SuperflowBase::determine_output_filename    Run mode: SuperflowRunMode::single_event_syst" << endl;
            cout << app_name << "SuperflowBase::determine_output_filename    Setting output file name to: " << sfile_name_.str() << endl;
            m_outputFileName = sfile_name_.str();
        }
        else {
            cout << app_name << "SuperflowBase::determine_output_filename    ERROR (Fatal)    Unknown event systematic! Code: " << static_cast<int>(m_singleEventSyst) << endl;
            exit(1);
        }
    }
    else if (nt.evt()->isMC) {
        stringstream suffix;
        if(m_outputFileNameSuffix!="")
            suffix << "_" << m_outputFileNameSuffix;
        else suffix << "";

        stringstream sfile_name_; // output file name
        sfile_name_ << "CENTRAL_" << nt.evt()->mcChannel << suffix.str() << ".root";
        if (m_runMode == SuperflowRunMode::nominal_and_weight_syst) {
            cout << app_name << "SuperflowBase::determine_output_filename    Run mode: SuperflowRunMode::nominal_and_weight_syst" << endl;
        }
        else {
            cout << app_name << "SuperflowBase::determine_output_filename    Run mode: SuperflowRunMode::nominal (weighted)" << endl;
        }
        cout << app_name << "SuperflowBase::determine_output_filename    Setting output file name to: " << sfile_name_.str() << endl;
        m_outputFileName = sfile_name_.str();
    }
    else {
        cout << app_name << "SuperflowBase::determine_output_filename    ERROR (Fatal)    Inconsistent setup." << endl;
        exit(1);
    }

    m_outputFile = new TFile(m_outputFileName.data(), "RECREATE");


}
///////////////////////////////////////////////////////////////////////////////
void SuperflowBase::initialize_output_files(TString input)
{

    determine_output_filename(input);

    stringstream tree_name;
    tree_name << "superNt";
    cout << app_name << "SuperflowBase::determine_output_filename    Setting output tree name to: " << tree_name.str() << endl; 
    m_tree_name_auto = tree_name.str();

    // initialize output tree
    m_HFT = new TTree(tree_name.str().data(), tree_name.str().data());
    m_HFT->SetDirectory(m_outputFile);
    m_HFT->SetAutoFlush(-16777216L);

    // define number of trees
    m_tree_leafs_size = m_varType.size() + 2 * index_weight_sys.size();// 2nd term may be zero
    m_weight_leaf_offset = m_varType.size();
   

    m_varFloat = new Float_t[m_tree_leafs_size]; // this one is larger to hold the syst_WEIGHT
    //m_varFloatArray = new Double_t*[m_varType.size()];
    m_varDouble = new Double_t[m_varType.size()];
    m_varInt = new Int_t[m_varType.size()];
    m_varBool = new Bool_t[m_varType.size()];

    // initialize HFT
    for (int i = 0; i < m_tree_leafs_size; i++) m_varFloat[i] = 1.0;
    for (int i = 0; i < m_varType.size(); i++) m_varDouble[i] = 1.0;
    for (int i = 0; i < m_varType.size(); i++) m_varInt[i] = 0;
    for (int i = 0; i < m_varType.size(); i++) m_varBool[i] = false;
    for (int i = 0; i < m_varType.size(); i++) { //m_varFloatArray[i] = 1.0; }
        vector<double> null(25, -999);
        vector<bool> f(25, false);
        m_varFloatArray.push_back(null);
        m_varBoolArray.push_back(f);
    }

    for (int i = 0; i < m_varType.size(); i++) {
        switch (m_varType[i]) {
            case SupervarType::sv_void: break;
            case SupervarType::sv_float: {
                string leaflist_ = m_varHFTName[i] + "/F";
                m_HFT->Branch(m_varHFTName[i].data(), m_varFloat + i, leaflist_.data(), 65536);
                break;
            }
            case SupervarType::sv_double: {
                string leaflist_ = m_varHFTName[i] + "/D";
                m_HFT->Branch(m_varHFTName[i].data(), m_varDouble + i, leaflist_.data(), 65536);
                break;
            }
            case SupervarType::sv_int: {
                string leaflist_ = m_varHFTName[i] + "/I";
                m_HFT->Branch(m_varHFTName[i].data(), m_varInt + i, leaflist_.data(), 65536);
                break;
            }
            case SupervarType::sv_bool: {
                string leaflist_ = m_varHFTName[i] + "/O";
                m_HFT->Branch(m_varHFTName[i].data(), m_varBool + i, leaflist_.data(), 65536);
                break;
            }
            case SupervarType::sv_float_array: {
                //string leaflist_ = m_varHFTName[i] + "[25]/D";
                //m_HFT->Branch(m_varHFTName[i].data(), m_varFloatArray.at(i), leaflist_.data(), 65536);
                m_HFT->Branch(m_varHFTName[i].data(), &m_varFloatArray[i]);
                break;
            }
            case SupervarType::sv_bool_array: {
                m_HFT->Branch(m_varHFTName[i].data(), &m_varBoolArray[i]);
                break;
            }
        }
    }

    for (int i = 0; i < index_weight_sys.size(); i++) {
        string syst_var_name_up = weight_prefix + m_sysStore[index_weight_sys[i]].tree_name + weight_suffix_up;
        string syst_var_name_down = weight_prefix + m_sysStore[index_weight_sys[i]].tree_name + weight_suffix_down;

        string leaflist_up = syst_var_name_up + "/F";
        string leaflist_down = syst_var_name_down + "/F";

        cout << app_name << "Weight var trees: " << syst_var_name_up << ", " << syst_var_name_down << endl;

        m_HFT->Branch(syst_var_name_up.data(), m_varFloat + m_weight_leaf_offset + 2 * i, leaflist_up.data(), 65536);
        m_HFT->Branch(syst_var_name_down.data(), m_varFloat + m_weight_leaf_offset + 2 * i + 1, leaflist_down.data(), 65536);
    }

    // Make an output file for each event systematic.
    if (m_runMode == SuperflowRunMode::all_syst) {

        m_output_array = new TFile*[index_event_sys.size()];
        for (int i = 0; i < index_event_sys.size(); i++) {

            stringstream suffix;
            if(m_outputFileNameSuffix!="")
                suffix << "_" << m_outputFileNameSuffix;
            else suffix << "";

            stringstream sfile_name_; // output file name
            sfile_name_ << Susy::NtSys::SusyNtSysNames.at(m_sysStore[index_event_sys[i]].event_syst) << "_" << nt.evt()->mcChannel << suffix.str() << ".root";
            m_output_array[i] = new TFile(sfile_name_.str().data(), "RECREATE");
        }

        m_HFT_array = new TTree*[index_event_sys.size()];
        for (int i = 0; i < index_event_sys.size(); i++) {
            m_HFT_array[i] = new TTree(tree_name.str().data(), tree_name.str().data());
            m_HFT_array[i]->SetDirectory(m_output_array[i]);
            m_HFT_array[i]->SetAutoFlush(-16777216L);
        }

        m_varFloat_array = new Float_t*[index_event_sys.size()];
        m_varDouble_array = new Double_t*[index_event_sys.size()];
        m_varInt_array = new Int_t*[index_event_sys.size()];
        m_varBool_array = new Bool_t*[index_event_sys.size()];

        m_varFloatArray_array.clear();
        m_varBoolArray_array.clear();
        vector<vector<vector<double>>> f(index_event_sys.size(), vector<vector<double>>(m_varType.size(),
                        vector<double>(25,1.0)));
        m_varFloatArray_array = f;

        vector<vector<vector<bool>>> fb(index_event_sys.size(), vector<vector<bool>>(m_varType.size(),
                        vector<bool>(25,false)));
        m_varBoolArray_array = fb;


        for (int i = 0; i < index_event_sys.size(); i++) {
            m_varFloat_array[i] = new Float_t[m_varType.size()];
            m_varDouble_array[i] = new Double_t[m_varType.size()];
            m_varInt_array[i] = new Int_t[m_varType.size()];
            m_varBool_array[i] = new Bool_t[m_varType.size()];

            for (int j = 0; j < m_varType.size(); j++) m_varFloat_array[i][j] = 1.0;
            for (int j = 0; j < m_varType.size(); j++) m_varDouble_array[i][j] = 1.0;
            for (int j = 0; j < m_varType.size(); j++) m_varInt_array[i][j] = 0;
            for (int j = 0; j < m_varType.size(); j++) m_varBool_array[i][j] = false;
        }

        for (int i = 0; i < index_event_sys.size(); i++) {
            for (int j = 0; j < m_varType.size(); j++) {
                switch (m_varType[j]) {
                    case SupervarType::sv_void: break;
                    case SupervarType::sv_float: {
                        string leaflist_ = m_varHFTName[j] + "/F";
                        m_HFT_array[i]->Branch(m_varHFTName[j].data(), m_varFloat_array[i] + j, leaflist_.data(), 65536);
                        break;
                    }
                    case SupervarType::sv_float_array: { 
                        m_HFT_array[i]->Branch(m_varHFTName[j].data(), &m_varFloatArray_array[i][j]); 
                        break;
                    }
                    case SupervarType::sv_bool_array: {
                        m_HFT_array[i]->Branch(m_varHFTName[j].data(), &m_varBoolArray_array[i][j]);
                        break;
                    }
                    case SupervarType::sv_double: {
                        string leaflist_ = m_varHFTName[j] + "/D";
                        m_HFT_array[i]->Branch(m_varHFTName[j].data(), m_varDouble_array[i] + j, leaflist_.data(), 65536);
                        break;
                    }
                    case SupervarType::sv_int: {
                        string leaflist_ = m_varHFTName[j] + "/I";
                        m_HFT_array[i]->Branch(m_varHFTName[j].data(), m_varInt_array[i] + j, leaflist_.data(), 65536);
                        break;
                    }
                    case SupervarType::sv_bool: {
                        string leaflist_ = m_varHFTName[j] + "/O";
                        m_HFT_array[i]->Branch(m_varHFTName[j].data(), m_varBool_array[i] + j, leaflist_.data(), 65536);
                        break;
                    }
                }
            }
        }
    }

}
///////////////////////////////////////////////////////////////////////////////
void SuperflowBase::write_output_files()
{
    // write nominal tree
    bool ok = true;
    if(m_outputFile->Write() == 0) ok = false;
    // write tree for each of the event systematics
    for (int i = 0; i < index_event_sys.size(); i++) { if( m_output_array[i]->Write() == 0) ok = false; }

    // close nominal file
    m_outputFile->Close();
    // close systematics files
    for (int i = 0; i < index_event_sys.size(); i++) m_output_array[i]->Close();

    if(!ok) {
        cout << app_name << "SuperflowBase::write_output_files    WARNING Problem with output files" << endl;
    }
    else {
        cout << app_name << "SuperflowBase::write_output_files    Files OK." << endl;
    }
}
///////////////////////////////////////////////////////////////////////////////
void SuperflowBase::Terminate()
{
    write_output_files();

    delete m_outputFile;
    for (int i = 0; i < index_event_sys.size(); i++) delete m_output_array[i];
    delete[] m_output_array;

    delete[] m_varFloat;
    delete[] m_varDouble;
    delete[] m_varInt;
    delete[] m_varBool;

    if (m_runMode == SuperflowRunMode::all_syst) {
        for (int i = 0; i < index_event_sys.size(); i++) delete[] m_varFloat_array[i];
        delete m_varFloat_array;
        for (int i = 0; i < index_event_sys.size(); i++) delete[] m_varDouble_array[i];
        delete m_varDouble_array;
        for (int i = 0; i < index_event_sys.size(); i++) delete[] m_varInt_array[i];
        delete m_varInt_array;
        for (int i = 0; i < index_event_sys.size(); i++) delete[] m_varBool_array[i];
        delete m_varBool_array;
    }
    
    if (m_runMode != SuperflowRunMode::single_event_syst) delete m_RunSyst;

    cout << app_name << "SuperflowBase::Terminate    Done." << endl;

}
///////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(CutName cut_)
{
    if (m_sysState == SupersysState::closed && m_varState == SupervarState::closed) {
        m_CutStore_Name_Exists = true;
        m_CutStoreNames.push_back(cut_.name);

        cout << app_name << "New cut: " << cut_.name << endl;
    }
    else {
        cout << app_name << "ERROR (Fatal): Cutflow operations are incorrectly ordered.";
        exit(1);
    }
    return *this;
}
///////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(std::function<bool(Superlink*)> cut_)
{
    if (m_sysState == SupersysState::closed && m_varState == SupervarState::closed) {
        m_CutStore.push_back(cut_);

        if (m_CutStore_Name_Exists) {
            m_CutStore_Name_Exists = false;
        }
        else {
            m_CutStoreNames.push_back("Untitled-" + to_string(m_CutStoreUntitled));
            cout << app_name << "New cut: " << "Untitled-" << m_CutStoreUntitled << endl;
            m_CutStoreUntitled++;
        }

        m_RawCounter.push_back(0.0);
        m_WeightCounter.push_back(0.0);
    }
    else {
        cout << app_name << "ERROR (Fatal): Cutflow operations are incorrectly ordered.";
        exit(1);
    }
    return *this;
}
///////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(NewVar new_var_name) // this is the NiceName
{
    if (m_varState == SupervarState::closed && m_sysState == SupersysState::closed) {
        m_varNiceName.push_back(new_var_name.name);
        m_superVar_hasNiceName = true;
        m_varState = SupervarState::open;
    }
    else {
        cout << app_name << "ERROR (Fatal): Close the Var using SaveVar().";
        exit(1);
    }
    return *this;
}
///////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(HFTname hft_name)
{
    if (m_varState == SupervarState::open && m_sysState == SupersysState::closed && !m_superVar_hasHFTName) {
        m_varHFTName.push_back(hft_name.name);
        m_superVar_hasHFTName = true;
    }
    else {
        cout << app_name << "ERROR (Fatal): Open a new Var using NewVar().";
        exit(1);
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(std::function<double(Superlink*, var_float*)> var_)
{
    if (m_varState == SupervarState::open && m_sysState == SupersysState::closed && !m_superVar_hasFunction) {
        m_varExprFloat.push_back(var_); // fill
        m_varExprFloatArray.push_back(m_nullExprFloatArray);
        m_varExprBoolArray.push_back(m_nullExprBoolArray);
        m_varExprDouble.push_back(m_nullExprDouble);
        m_varExprInt.push_back(m_nullExprInt);
        m_varExprBool.push_back(m_nullExprBool);
        m_varExprVoid.push_back(m_nullExprVoid);

        m_varType.push_back(SupervarType::sv_float);
        m_superVar_hasFunction = true;
    }
    else {
        cout << app_name << "ERROR (Fatal): First open a new Var using NewVar().";
        exit(1);
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(std::function<vector<double>(Superlink*, var_float_array*)> var_)
{
    if (m_varState == SupervarState::open && m_sysState == SupersysState::closed && !m_superVar_hasFunction) {
        m_varExprFloat.push_back(m_nullExprFloat);
        m_varExprDouble.push_back(m_nullExprDouble);
        m_varExprBoolArray.push_back(m_nullExprBoolArray);
        m_varExprFloatArray.push_back(var_);
        m_varExprInt.push_back(m_nullExprInt);
        m_varExprBool.push_back(m_nullExprBool);
        m_varExprVoid.push_back(m_nullExprVoid);

        m_varType.push_back(SupervarType::sv_float_array);
        m_superVar_hasFunction = true;
    }
    else {
        cout << app_name << "ERROR (Fatal): First open a new Var using NewVar().";
        exit(1);
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(std::function<vector<bool>(Superlink*, var_bool_array*)> var_)
{
    if (m_varState == SupervarState::open && m_sysState == SupersysState::closed && !m_superVar_hasFunction) {
        m_varExprFloat.push_back(m_nullExprFloat);
        m_varExprDouble.push_back(m_nullExprDouble);
        m_varExprBoolArray.push_back(var_);
        m_varExprFloatArray.push_back(m_nullExprFloatArray);
        m_varExprInt.push_back(m_nullExprInt);
        m_varExprBool.push_back(m_nullExprBool);
        m_varExprVoid.push_back(m_nullExprVoid);

        m_varType.push_back(SupervarType::sv_bool_array);
        m_superVar_hasFunction = true;
    }
    else {
        cout << app_name << "ERROR (Fatal): First open a new Var using NewVar().";
        exit(1);
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(std::function<double(Superlink*, var_double*)> var_)
{
    if (m_varState == SupervarState::open && m_sysState == SupersysState::closed && !m_superVar_hasFunction) {
        m_varExprFloat.push_back(m_nullExprFloat);
        m_varExprFloatArray.push_back(m_nullExprFloatArray);
        m_varExprBoolArray.push_back(m_nullExprBoolArray);
        m_varExprDouble.push_back(var_); // fill
        m_varExprInt.push_back(m_nullExprInt);
        m_varExprBool.push_back(m_nullExprBool);
        m_varExprVoid.push_back(m_nullExprVoid);

        m_varType.push_back(SupervarType::sv_double);
        m_superVar_hasFunction = true;
    }
    else {
        cout << app_name << "ERROR (Fatal): First open a new Var using NewVar().";
        exit(1);
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(std::function<int(Superlink*, var_int*)> var_)
{
    if (m_varState == SupervarState::open && m_sysState == SupersysState::closed && !m_superVar_hasFunction) {
        m_varExprFloat.push_back(m_nullExprFloat);
        m_varExprDouble.push_back(m_nullExprDouble);
        m_varExprFloatArray.push_back(m_nullExprFloatArray);
        m_varExprBoolArray.push_back(m_nullExprBoolArray);
        m_varExprInt.push_back(var_); // fill
        m_varExprBool.push_back(m_nullExprBool);
        m_varExprVoid.push_back(m_nullExprVoid);

        m_varType.push_back(SupervarType::sv_int);
        m_superVar_hasFunction = true;
    }
    else {
        cout << app_name << "ERROR (Fatal): First open a new Var using NewVar().";
        exit(1);
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(std::function<bool(Superlink*, var_bool*)> var_)
{
    if (m_varState == SupervarState::open && m_sysState == SupersysState::closed && !m_superVar_hasFunction) {
        m_varExprFloat.push_back(m_nullExprFloat);
        m_varExprDouble.push_back(m_nullExprDouble);
        m_varExprFloatArray.push_back(m_nullExprFloatArray);
        m_varExprBoolArray.push_back(m_nullExprBoolArray);
        m_varExprInt.push_back(m_nullExprInt);
        m_varExprBool.push_back(var_); // fill
        m_varExprVoid.push_back(m_nullExprVoid);

        m_varType.push_back(SupervarType::sv_bool);
        m_superVar_hasFunction = true;
    }
    else {
        cout << app_name << "ERROR (Fatal): First open a new Var using NewVar().";
        exit(1);
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(std::function<void(Superlink*, var_void*)> var_)
{
    m_varExprFloat.push_back(m_nullExprFloat);
    m_varExprDouble.push_back(m_nullExprDouble);
    m_varExprFloatArray.push_back(m_nullExprFloatArray);
    m_varExprBoolArray.push_back(m_nullExprBoolArray);
    m_varExprInt.push_back(m_nullExprInt);
    m_varExprBool.push_back(m_nullExprBool);
    m_varExprVoid.push_back(var_);// fill

    m_varNiceName.push_back("void");
    m_varHFTName.push_back("void");
    m_varType.push_back(SupervarType::sv_void);
    return *this;
}
//////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(SaveVar save_var)
{
    if (m_varState == SupervarState::open && m_superVar_hasFunction && m_superVar_hasNiceName) {
        if (!m_superVar_hasHFTName) {
            m_varHFTName.push_back("Untitled_" + to_string(m_superVar_Untitled));
            m_superVar_Untitled++;
            m_superVar_hasHFTName = true;
        }

        m_varState = SupervarState::closed;
        m_superVar_hasFunction = false;
        m_superVar_hasNiceName = false;
        m_superVar_hasHFTName = false;

        cout << app_name << "New var: " << m_varNiceName.back() << endl;
        cout << app_name << "    HFT: " << m_varHFTName.back() << endl;
    }
    else {
        if (m_varState != SupervarState::open) {
            cout << app_name << "ERROR (Fatal): First open a new Var using NewVar().";
        }
        else {
            cout << app_name << "ERROR (Fatal): A lambda-expression is required.";
        }
        exit(1);
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(NewSystematic new_sys) // this is the NiceName
{
    if (m_sysState == SupersysState::closed && m_varState == SupervarState::closed) {
        m_sysTemplate.name = new_sys.name;
        m_sys_hasNiceName = true;
        m_sysState = SupersysState::open;
    }
    else {
        cout << app_name << "ERROR (Fatal): Close using SaveSystematic().";
        exit(1);
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(TreeName tree_name)
{
    if (m_sysState == SupersysState::open && m_varState == SupervarState::closed && !m_sys_hasTreeName) {
        m_sysTemplate.tree_name = tree_name.name;
        m_sys_hasTreeName = true;
    }
    else {
        cout << app_name << "ERROR (Fatal): Open a NewSystematic().";
        exit(1);
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(EventSystematic obj_)
{
    if (m_sysState == SupersysState::open && m_varState == SupervarState::closed && !m_sys_hasSystematic) {
        m_sysTemplate.event_syst = obj_.event_syst_;
        m_sysTemplate.weight_syst = SupersysWeight::null;
        m_sys_hasSystematic = true;

        m_sysTemplate.type = SupersysType::event;
        m_sys_hasType = true;
    }
    else {
        cout << app_name << "ERROR (Fatal): Open a NewSystematic().";
        exit(1);
    }
    return *this;
}

SuperflowBase& SuperflowBase::operator<<(WeightSystematic obj_)
{
    if (m_sysState == SupersysState::open && m_varState == SupervarState::closed && !m_sys_hasSystematic) {
        m_sysTemplate.event_syst = Susy::NtSys::NOM;
        m_sysTemplate.weight_syst_up = obj_.weight_syst_up;
        m_sysTemplate.weight_syst_down = obj_.weight_syst_down;
        m_sys_hasSystematic = true;

        m_sysTemplate.type = SupersysType::weight;
        m_sys_hasType = true;
    }
    else {
        cout << app_name << "ERROR (Fatal): Open a NewSystematic().";
        exit(1);
    }
    return *this;
}
//////////////////////////////////////////////////////////////////////////////
SuperflowBase& SuperflowBase::operator<<(SaveSystematic save_var)
{
    if (m_sysState == SupersysState::open
        && m_varState == SupervarState::closed
        && m_sys_hasNiceName
        && m_sys_hasTreeName
        && m_sys_hasType
        && m_sys_hasSystematic
        ) {
        if (m_sysTemplate.name == "") m_sysTemplate.name = m_sysTemplate.tree_name;

        m_sysState = SupersysState::closed;
        m_sys_hasNiceName = false;
        m_sys_hasTreeName = false;
        m_sys_hasSystematic = false;
        m_sys_hasType = false;

        m_sysStore.push_back(m_sysTemplate);
        m_sysTemplate.reset();

        cout << app_name << "New systematic: " << m_sysStore.back().name << endl;
        if (m_sysStore.back().type == SupersysType::event) {
            cout << app_name << "    event systematic: " << m_sysStore.back().tree_name << endl;
        }
        else if (m_sysStore.back().type == SupersysType::weight) {
            cout << app_name << "    weight systematics: " << m_sysStore.back().tree_name << weight_suffix_up
                << "/" << m_sysStore.back().tree_name << weight_suffix_down << endl;
        }
        else {
            cout << app_name << "ERROR (Fatal): Impossible SupersysType.";
            exit(1);
        }

    }
    else {
        if (m_sysState != SupersysState::open) {
            cout << app_name << "ERROR (Fatal): First open a NewSystematic().";
        }
        else {
            cout << app_name << "ERROR (Fatal): Can't save incomplete systematic object.";
        }
        exit(1);
    }
    return *this;
}

}; //namespace sflow
