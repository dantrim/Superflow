// std
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <random>
#include <fstream>

// Superflow
#include "Superflow/Superflow.h"
#include "Superflow/StringTools.h"

using namespace std;

namespace sflow {

///////////////////////////////////////////////////////////////////////////////
Superflow::Superflow() :
    m_mcWeighter(nullptr)
{
}
///////////////////////////////////////////////////////////////////////////////
Superflow::~Superflow()
{}
///////////////////////////////////////////////////////////////////////////////
void Superflow::attach_superlink(Superlink* sl_)
{
    sl_->tools = &m_nttools;

    sl_->anaType = m_nttools.getAnaType();

    sl_->nt = &nt;
    sl_->weights = m_weights;
    sl_->nt_sys = m_RunSyst->event_syst;

    sl_->preElectrons = &m_preElectrons;
    sl_->preMuons = &m_preMuons;
    sl_->preTaus = &m_preTaus;
    sl_->preJets = &m_preJets;

    sl_->baseLeptons = &m_baseLeptons;
    sl_->baseElectrons = &m_baseElectrons;
    sl_->baseMuons = &m_baseMuons;
    sl_->baseTaus = &m_baseTaus;
    sl_->baseJets = &m_baseJets;

    sl_->leptons = &m_signalLeptons;
    sl_->electrons = &m_signalElectrons;
    sl_->muons = &m_signalMuons;
    sl_->taus = &m_signalTaus;
    sl_->jets = &m_signalJets;

    sl_->met = m_met;
    sl_->trackMet = m_trackMet;

    if (nt.evt()->isMC) {
        sl_->isMC = true;
    }
    else {
        sl_->isData = true;
    }
}
///////////////////////////////////////////////////////////////////////////////
void Superflow::Begin(TTree* /*tree*/)
{
    cout << app_name << "Superflow::Begin" << endl;
    SusyNtAna::Begin(0);

    cout << app_name << "Superflow::Begin    Input sample: " << m_sample << endl;

    if (m_runMode == SuperflowRunMode::null) {
        cout << app_name << "Superflow::Begin ERROR    SuperflowRunMode is \"null\"! -- Missing call to Superflow::setRunMode(), exiting" << endl;
        exit(1);
    }

    if (m_varState != SupervarState::closed || m_sysState != SupersysState::closed) {
        cout << app_name << "Superflow::Begin ERROR    Invalid setup of variables in executable -- Close the Var using SaveVar(), exiting" << endl;;
        exit(1);
    }

    if (m_runMode == SuperflowRunMode::single_event_syst && m_singleEventSyst == Susy::NtSys::NOM) {
        cout << app_name << "Superflow::Begin ERROR    SuperflowRunMode::single_event_syst: Call setSingleEventSyst(SusyNtSys nt_syst_)." << endl;
        exit(1);
    }

    SuperflowBase::Begin(0);

}
///////////////////////////////////////////////////////////////////////////////
void Superflow::Init(TTree* tree)
{
    cout << app_name << "Superflow::Init" << endl;
    SusyNtAna::Init(tree);
    
    TString input_sample = m_input_chain->GetFile()->Get("inputContainerName")->GetTitle();
    TString output_sample = m_input_chain->GetFile()->Get("outputContainerName")->GetTitle();
    TString nt_tag = m_input_chain->GetFile()->Get("productionTag")->GetTitle();
    TString prod_command = m_input_chain->GetFile()->Get("productionCommand")->GetTitle();

    cout << endl;
    cout << app_name << "Superflow::Init    ================ Sample Information ==============" << endl;
    cout << app_name << "Superflow::Init     Input container name  : " << input_sample  << endl;
    cout << app_name << "Superflow::Init     Output container name : " << output_sample << endl;
    cout << app_name << "Superflow::Init     SusyNt production tag : " << nt_tag << endl;
    cout << app_name << "Superflow::Init     SusyNt production cmd : " << prod_command << endl; 
    cout << app_name << "Superflow::Init    ==================================================" << endl;

    if (!nt.evt()->isMC) {
        cout << app_name << "Superflow::Init    Switching run mode to SuperflowRunMode::data" << endl;
        m_runMode = SuperflowRunMode::data;
    }
    else if (nt.evt()->isMC) {
        initialize_mc_weighter(tree);
    }

    initialize_output_files(input_sample);
}
///////////////////////////////////////////////////////////////////////////////
//Bool_t Superflow::Notify()
//{
//    SuperflowBase::Notify();
//    return kTRUE;
//}
///////////////////////////////////////////////////////////////////////////////
Bool_t Superflow::Process(Long64_t entry)
{
    GetEntry(entry);

    m_chainEntry++; // SusyNtAna counter

    //////////////////////////////////////////////
    // process counter
    //////////////////////////////////////////////
    if (m_chainEntry % 5000 == 0) {
        cout << app_name << "Superflow::Process    **** Processing entry " << setw(6) << m_chainEntry
            << " run "   << setw(6) << nt.evt()->run
            << " event " << setw(7) << nt.evt()->eventNumber << " ****" << endl;
    }


    /////////////////////////////////////////////////////////////////////////
    // select the objects for the selected run modes
    ////////////////////////////////////////////////////////////////////////

    // these are used to select the signature of the lambdas
    var_float* vf_          = nullptr;
    var_double* vd_         = nullptr;
    var_float_array* vfa_   = nullptr;
    var_int_array* via_     = nullptr;
    var_bool_array* vba_    = nullptr;
    var_int* vi_            = nullptr;
    var_bool* vb_           = nullptr;
    var_void* vv_           = nullptr;

    switch (m_runMode) {
        ////////////////////////////////////////////////////////////////////
        // data
        ////////////////////////////////////////////////////////////////////
        case SuperflowRunMode::data: {
            SusyNtAna::clearObjects();
            SusyNtAna::selectObjects(m_RunSyst->event_syst);

            m_weights = new Superweight();
            Superlink* sl_ = new Superlink;
            attach_superlink(sl_);

            // loop over the loaded cuts
            bool pass_cuts = true;
            if (m_CutStore.size() > 0) {
                for (int i = 0; i < (int)m_CutStore.size(); i++) {
                    pass_cuts = m_CutStore[i](sl_);
                    if (pass_cuts) {
                        m_RawCounter[i]++;
                    }
                    else {
                        break;
                    }
                }
            }

            // we have passed the cuts, now fill the trees
            if (pass_cuts) {
                for (int v_ = 0; v_ < (int)m_varType.size(); v_++) {
                    switch (m_varType[v_]) {
                        case SupervarType::sv_float: {
                            m_varFloat[v_] = m_varExprFloat[v_](sl_, vf_);
                            break;
                        }
                        case SupervarType::sv_double: {
                            m_varDouble[v_] = m_varExprDouble[v_](sl_, vd_);
                            break;
                        }
                        case SupervarType::sv_float_array: { 
                            m_varFloatArray[v_] = m_varExprFloatArray[v_](sl_,vfa_);
                            break;
                        }
                        case SupervarType::sv_int_array: { 
                            m_varIntArray[v_] = m_varExprIntArray[v_](sl_,via_);
                            break;
                        }
                        case SupervarType::sv_bool_array: {
                            m_varBoolArray[v_] = m_varExprBoolArray[v_](sl_,vba_);
                            break;
                        }
                        case SupervarType::sv_int: {
                            m_varInt[v_] = m_varExprInt[v_](sl_, vi_);
                            break;
                        }
                        case SupervarType::sv_bool: {
                            m_varBool[v_] = m_varExprBool[v_](sl_, vb_);
                            break;
                        }
                        case SupervarType::sv_void: {
                            m_varExprVoid[v_](sl_, vv_);
                            break;
                        }
                    }
                }
                m_HFT->Fill();
            }
            delete sl_;
            delete m_weights;

        } break;
        
        ////////////////////////////////////////////////////////////////////
        // nominal, single_event_syst
        ////////////////////////////////////////////////////////////////////
        case SuperflowRunMode::nominal:
        case SuperflowRunMode::single_event_syst: {

            SusyNtAna::clearObjects();
            SusyNtAna::selectObjects(m_RunSyst->event_syst);

            m_weights = new Superweight();
            Superlink* sl_ = new Superlink;
            attach_superlink(sl_);

            bool pass_cuts = true;

            if (m_countWeights) {
                computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, m_weights);
            }

            if (m_CutStore.size() > 0) {
                for (int i = 0; i < (int)m_CutStore.size(); i++) {
                    pass_cuts = m_CutStore[i](sl_); // run the cut function
                    if (pass_cuts) {
                        m_RawCounter[i]++;
                        if (m_countWeights) m_WeightCounter[i] += m_weights->product();
                    }
                    else {
                        break;
                    }
                }
            }

            if (pass_cuts) {
                if (!m_countWeights) {
                    computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, m_weights);
                    m_WeightCounter[m_CutStore.size() - 1] += m_weights->product();
                }

                // FILL_HFTs
                for (int v_ = 0; v_ < (int)m_varType.size(); v_++) {
                    switch (m_varType[v_]) {
                        case SupervarType::sv_float: {
                            m_varFloat[v_] = m_varExprFloat[v_](sl_, vf_);
                            break;
                        }
                        case SupervarType::sv_float_array: { 
                            m_varFloatArray[v_] = m_varExprFloatArray[v_](sl_,vfa_);
                            break;
                        }
                        case SupervarType::sv_int_array: { 
                            m_varIntArray[v_] = m_varExprIntArray[v_](sl_,via_);
                            break;
                        }
                        case SupervarType::sv_bool_array: {
                            m_varBoolArray[v_] = m_varExprBoolArray[v_](sl_,vba_);
                            break;
                        }
                        case SupervarType::sv_double: {
                            m_varDouble[v_] = m_varExprDouble[v_](sl_, vd_);
                            break;
                        }
                        case SupervarType::sv_int: {
                            m_varInt[v_] = m_varExprInt[v_](sl_, vi_);
                            break;
                        }
                        case SupervarType::sv_bool: {
                            m_varBool[v_] = m_varExprBool[v_](sl_, vb_);
                            break;
                        }
                        case SupervarType::sv_void: {
                            m_varExprVoid[v_](sl_, vv_);
                            break;
                        }
                    }
                }
                m_HFT->Fill();
            }
            delete sl_;
            delete m_weights;

        } break;
        ////////////////////////////////////////////////////////////////////
        // processing all systematics
        ////////////////////////////////////////////////////////////////////
        case SuperflowRunMode::all_syst: {
            for (int i = 0; i < (int)index_event_sys.size(); i++) { // loop over event systematics
                SusyNtAna::clearObjects();
                delete m_RunSyst;

                m_RunSyst = &m_sysStore[index_event_sys[i]]; // don't delete!!
                SusyNtAna::selectObjects(m_RunSyst->event_syst);

                m_weights = new Superweight();
                Superlink* sl_ = new Superlink;
                attach_superlink(sl_);

                bool pass_cuts = true;

                if (m_CutStore.size() > 0) {
                    for (int k = 0; k < (int)m_CutStore.size(); k++) {
                        pass_cuts = m_CutStore[k](sl_);
                        if (!pass_cuts) break;
                    }
                }

                if (pass_cuts) {
                    computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, m_weights);
                    // FILL_HFTs
                    for (int v_ = 0; v_ < (int)m_varType.size(); v_++) {
                        switch (m_varType[v_]) {
                            case SupervarType::sv_float: {
                                m_varFloat_array[i][v_] = m_varExprFloat[v_](sl_, vf_);
                                break;
                            }
                            case SupervarType::sv_double: {
                                m_varDouble_array[i][v_] = m_varExprDouble[v_](sl_, vd_);
                                break;
                            }
                            case SupervarType::sv_float_array: { 
                                m_varFloatArray_array[i][v_] = m_varExprFloatArray[v_](sl_,vfa_);
                                break;
                            }
                            case SupervarType::sv_int_array: { 
                                m_varIntArray_array[i][v_] = m_varExprIntArray[v_](sl_,via_);
                                break;
                            }
                            case SupervarType::sv_bool_array: {
                                m_varBoolArray_array[i][v_] = m_varExprBoolArray[v_](sl_,vba_);
                                break;
                            }
                            case SupervarType::sv_int: {
                                m_varInt_array[i][v_] = m_varExprInt[v_](sl_, vi_);
                                break;
                            }
                            case SupervarType::sv_bool: {
                                m_varBool_array[i][v_] = m_varExprBool[v_](sl_, vb_);
                                break;
                            }
                            case SupervarType::sv_void: {
                                m_varExprVoid[v_](sl_, vv_);
                                break;
                            }
                        }
                    }
                    m_HFT_array[i]->Fill();
                }
                delete sl_;
                delete m_weights;

                m_RunSyst = nullptr;
            }
        }
        ////////////////////////////////////////////////////////////////////
        // processing nominal and weight systematics
        ////////////////////////////////////////////////////////////////////
        case SuperflowRunMode::nominal_and_weight_syst: {
            delete m_RunSyst;
            m_RunSyst = new Supersys(SupersysType::central);

            SusyNtAna::clearObjects();
            SusyNtAna::selectObjects(m_RunSyst->event_syst);

            m_weights = new Superweight();
            Superlink* sl_ = new Superlink;
            attach_superlink(sl_);

            if (m_countWeights) {
                computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, m_weights);
            }

            bool pass_cuts = true;

            if (m_CutStore.size() > 0) {
                for (int i = 0; i < (int)m_CutStore.size(); i++) {
                    pass_cuts = m_CutStore[i](sl_); // run the cut function

                    if (pass_cuts) {
                        m_RawCounter[i]++;
                        if (m_countWeights) m_WeightCounter[i] += m_weights->product();
                    }
                    else {
                        break;
                    }
                }
            }

            if (pass_cuts) {
                computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, m_weights);

                double nom_eventweight = m_weights->product();
                if(!m_countWeights)
                    m_WeightCounter[m_CutStore.size() - 1] += m_weights->product();

                // FILL HFTs
                for (int v_ = 0; v_ < (int)m_varType.size(); v_++) {
                    switch (m_varType[v_]) {
                        case SupervarType::sv_float: {
                            m_varFloat[v_] = m_varExprFloat[v_](sl_, vf_);
                            break;
                        }
                        case SupervarType::sv_double: {
                            m_varDouble[v_] = m_varExprDouble[v_](sl_, vd_);
                            break;
                        }
                        case SupervarType::sv_float_array: { 
                            m_varFloatArray[v_] = m_varExprFloatArray[v_](sl_,vfa_);
                            break;
                        }
                        case SupervarType::sv_int_array: { 
                            m_varIntArray[v_] = m_varExprIntArray[v_](sl_,via_);
                            break;
                        }
                        case SupervarType::sv_bool_array: {
                            m_varBoolArray[v_] = m_varExprBoolArray[v_](sl_,vba_);
                            break;
                        }
                        case SupervarType::sv_int: {
                            m_varInt[v_] = m_varExprInt[v_](sl_, vi_);
                            break;
                        }
                        case SupervarType::sv_bool: {
                            m_varBool[v_] = m_varExprBool[v_](sl_, vb_);
                            break;
                        }
                        case SupervarType::sv_void: {
                            m_varExprVoid[v_](sl_, vv_);
                            break;
                        }
                    }
                }

                /////////////////////////////////////////
                // fill the weight variations
                /////////////////////////////////////////
                for (int w_ = 0; w_ < (int)index_weight_sys.size(); w_++) {
                    Superweight* weightComponents_copy = new Superweight(*m_weights);

                    // Up variation
                    m_RunSyst->weight_syst = m_sysStore[index_weight_sys[w_]].weight_syst_up; // do up variation
                    computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, weightComponents_copy);
                    double up_weight = weightComponents_copy->product();
                    delete weightComponents_copy;

                    weightComponents_copy = new Superweight(*m_weights);
                    // Down variation
                    m_RunSyst->weight_syst = m_sysStore[index_weight_sys[w_]].weight_syst_down; // do down variation
                    computeWeights(nt, *m_mcWeighter, m_signalLeptons, m_baseJets, m_RunSyst, weightComponents_copy);
                    double down_weight = weightComponents_copy->product();

                    delete weightComponents_copy;
                    weightComponents_copy = nullptr;

                    if (up_weight < nom_eventweight && down_weight > nom_eventweight) {
                        double temp = up_weight; // reason to swap these
                        up_weight = down_weight;
                        down_weight = temp;
                    }
                    if (nom_eventweight > m_epsilon) {
                        up_weight = up_weight / nom_eventweight;
                        down_weight = down_weight / nom_eventweight;
                    }
                    else {
                        up_weight = 1.0;
                        down_weight = 1.0;
                    }

                    *(m_varFloat + m_weight_leaf_offset + 2 * w_) = up_weight; // put in TBranch
                    *(m_varFloat + m_weight_leaf_offset + 2 * w_ + 1) = down_weight;
                }
                m_RunSyst->weight_syst = SupersysWeight::null; // must reset this value!
                m_HFT->Fill();
            }
            delete sl_;
            delete m_weights;
        } break;
        default: break;
    }
    return kTRUE;
}
///////////////////////////////////////////////////////////////////////////////
void Superflow::Terminate()
{
    cout << app_name << "Superflow::Terminate" << endl;

    cout << app_name << "------------------------------ ------------------------------" << endl;
    cout << app_name << " Raw cutflow " << endl;
    cout << app_name << "----------------------------- -------------------------------" << endl;
    cout << std::fixed;
    cout << std::setprecision(0);
    for (int i = 0; i < (int)m_CutStore.size(); i++) {
        cout << app_name << "Cut " << pad_width(to_string(i), 2) << ": " << pad_width(m_CutStoreNames[i], 32) << ": " << m_RawCounter[i] << endl;
    }
    cout << app_name << endl << app_name << endl;

    cout << app_name << "------------------------------ ------------------------------" << endl;
    cout << app_name << " Weighted cutflow " << endl;
    cout << app_name << "------------------------------- -----------------------------" << endl;
    cout << std::resetiosflags(std::ios::floatfield);
    cout << std::resetiosflags(std::ios::adjustfield);
    cout << std::setprecision(6);
    for (int i = 0; i < (int)m_CutStore.size(); i++) {
        cout << app_name << "Cut " << pad_width(to_string(i), 2) << ": " << pad_width(m_CutStoreNames[i], 32) << ": " << m_WeightCounter[i] << endl;
    }
    cout << app_name << endl << app_name << endl;

    SuperflowBase::Terminate();
    SusyNtAna::Terminate();

    //dantrim -- use MCWeighter from SusyNtAna
    //if (m_mcWeighter) delete m_mcWeighter;
    //if (m_trigObj) delete m_trigObj;

}
///////////////////////////////////////////////////////////////////////////////
bool Superflow::initialize_mc_weighter(TTree *tree)
{
    cout << app_name << "Superflow::initialize_mc_weighter    Initializing MCWeighter" << endl;
    bool success = false;
    if (tree) {
        m_mcWeighter = &mcWeighter(); // use MCWeighter instance from SusyNtAna
        m_mcWeighter->setVerbose(m_dbg);
    
        //m_mcWeighter->printSumwMap();
        m_mcWeighter->printSumwMapPeriod();
        if (m_dbg) {
            cout << app_name << "Superflow::initialize_mc_weighter    MCWeighter has been initialized." << endl;
        }
    }
    else {
        cout << app_name << "Superflow::initialize_mc_weighter ERROR    Invalid input tree, cannot initialize MCWeighter." << endl;
        cout << app_name << "Superflow::initialize_mc_weighter ERROR    >>> Exiting." << endl;
        exit(1);
    }
    return success;
}
///////////////////////////////////////////////////////////////////////////////
bool Superflow::computeWeights( Susy::SusyNtObject &ntobj, MCWeighter &weighter,
            const LeptonVector& leptons, const JetVector& jets,
            Supersys* super_sys, Superweight* weights_)
{

    // MCWeighter's susynt-weight calculation
    if (ntobj.evt()->isMC) {
        //MCWeighter::WeightSys wSys = MCWeighter::Sys_NOM;
        NtSys::SusyNtSys wSys = NtSys::NOM;


        bool do_susynt_w = false;

        switch (super_sys->weight_syst) {
            case SupersysWeight::PILEUP_UP : {
                wSys = NtSys::PILEUP_UP;
                do_susynt_w = true;
                break;
            }
            case SupersysWeight::PILEUP_DN : {
                wSys = NtSys::PILEUP_DN;
                do_susynt_w = true;
                break;
            }
            case SupersysWeight::null : {
                do_susynt_w = true;
            }
            default: break;
        } // switch
        if (do_susynt_w) {
            weights_->susynt = weighter.getMCWeight(ntobj.evt(), m_luminosity, wSys, /*multi-period*/false, /*include PRW*/false);
            weights_->susynt_multi = weighter.getMCWeight(ntobj.evt(), m_luminosity, wSys, /*multi-period*/true, /*include PRW*/false);
        }

        // JVT eff
        bool do_jvt = false;
        switch (super_sys->weight_syst) {
            case SupersysWeight::JET_JVTEff_UP :
            case SupersysWeight::JET_JVTEff_DN : {
                do_jvt = true;
                break;
            } default: break;
        } // switch
        if(do_jvt && jets.size()>0) {
            double outsf = 1.0;
            for(int i = 0; i < (int)jets.size(); i++) {
                double sf = jets[i]->jvtEff;
                double delta = 0.0;
                if(super_sys->weight_syst == SupersysWeight::JET_JVTEff_UP) {
                    delta = jets[i]->jvtEff_up;
                }
                else if(super_sys->weight_syst == SupersysWeight::JET_JVTEff_DN) {
                    delta = jets[i]->jvtEff_dn;
                }
                outsf *= (sf + delta);
            } // i
            weights_->jvtSf = outsf;
        } // do jvt 

        if(leptons.size()>=1) {
            bool do_lepSf_ = false;
            switch(super_sys->weight_syst) {
                case SupersysWeight::EL_EFF_ChargeIDSel_DN :
                case SupersysWeight::EL_EFF_ChargeIDSel_UP :
                case SupersysWeight::EL_EFF_ID_TOTAL_Uncorr_DN :
                case SupersysWeight::EL_EFF_ID_TOTAL_Uncorr_UP :
                case SupersysWeight::EL_EFF_Iso_TOTAL_Uncorr_DN :
                case SupersysWeight::EL_EFF_Iso_TOTAL_Uncorr_UP :
                case SupersysWeight::EL_EFF_Reco_TOTAL_Uncorr_DN :
                case SupersysWeight::EL_EFF_Reco_TOTAL_Uncorr_UP :
                case SupersysWeight::EL_EFF_TriggerEff_TOTAL_DN :
                case SupersysWeight::EL_EFF_TriggerEff_TOTAL_UP :
                case SupersysWeight::EL_EFF_Trigger_TOTAL_DN :
                case SupersysWeight::EL_EFF_Trigger_TOTAL_UP :
                case SupersysWeight::MUON_EFF_BADMUON_STAT_DN :
                case SupersysWeight::MUON_EFF_BADMUON_STAT_UP :
                case SupersysWeight::MUON_EFF_BADMUON_SYS_DN :
                case SupersysWeight::MUON_EFF_BADMUON_SYS_UP :
                case SupersysWeight::MUON_EFF_ISO_STAT_DN :
                case SupersysWeight::MUON_EFF_ISO_STAT_UP :
                case SupersysWeight::MUON_EFF_ISO_SYS_DN :
                case SupersysWeight::MUON_EFF_ISO_SYS_UP :
                case SupersysWeight::MUON_EFF_RECO_STAT_DN :
                case SupersysWeight::MUON_EFF_RECO_STAT_UP :
                case SupersysWeight::MUON_EFF_RECO_SYS_DN :
                case SupersysWeight::MUON_EFF_RECO_SYS_UP :
                case SupersysWeight::MUON_EFF_RECO_STAT_LOWPT_DN :
                case SupersysWeight::MUON_EFF_RECO_STAT_LOWPT_UP :
                case SupersysWeight::MUON_EFF_RECO_SYS_LOWPT_DN :
                case SupersysWeight::MUON_EFF_RECO_SYS_LOWPT_UP :
                case SupersysWeight::MUON_EFF_TTVA_STAT_DN :
                case SupersysWeight::MUON_EFF_TTVA_STAT_UP :
                case SupersysWeight::MUON_EFF_TTVA_SYS_DN :
                case SupersysWeight::MUON_EFF_TTVA_SYS_UP :
                case SupersysWeight::MUON_EFF_TrigStat_DN :
                case SupersysWeight::MUON_EFF_TrigStat_UP :
                case SupersysWeight::MUON_EFF_TrigSys_DN :
                case SupersysWeight::MUON_EFF_TrigSys_UP :
                case SupersysWeight::null : {
                    do_lepSf_ = true;
                }; break;
                default : break;
            } // switch
            if(do_lepSf_) {
                float outSF = 1.;
                for(unsigned int iL = 0; iL < leptons.size(); iL++) {
                    const Lepton &lep = *(leptons[iL]);
                    outSF *= computeLeptonEfficiencySf(lep, super_sys->weight_syst);
                }
                weights_->lepSf = outSF;
            } // do_lepSf
        } // lepts >= 1

        //flavor tagging systematics

        bool do_btag_sys = false;
        bool do_btag_nom = false;

        switch (super_sys->weight_syst) {
            case SupersysWeight::FT_EFF_B_UP :
            case SupersysWeight::FT_EFF_B_DN :
            case SupersysWeight::FT_EFF_C_UP :
            case SupersysWeight::FT_EFF_C_DN :
            case SupersysWeight::FT_EFF_LT_UP :
            case SupersysWeight::FT_EFF_LT_DN :
            case SupersysWeight::FT_EFF_EXTRAP_UP :
            case SupersysWeight::FT_EFF_EXTRAP_DN :
            case SupersysWeight::FT_EFF_EXTRAPC_UP :
            case SupersysWeight::FT_EFF_EXTRAPC_DN : {
                do_btag_sys = true;
                do_btag_nom = false;
                break;
            };
            case SupersysWeight::null : {
                do_btag_sys = false;
                do_btag_nom = true;
                break;
            };
            default: break;
        } //switch
        if (do_btag_nom) {
            double btagSF = 1.0;
            if(jets.size()>0) { btagSF = m_nttools.bTagSF(jets); } 
            weights_->btagSf = btagSF;
        }
        if (do_btag_sys) {
            double btagSF = 1.0;
            NtSys::SusyNtSys sys = NtSys::SusyNtSys::NOM;
            {
                switch(super_sys->weight_syst) {
                    case SupersysWeight::FT_EFF_B_UP :
                        sys = NtSys::SusyNtSys::FT_EFF_B_systematics_UP;
                        break;
                    case SupersysWeight::FT_EFF_B_DN : 
                        sys = NtSys::SusyNtSys::FT_EFF_B_systematics_DN;
                        break;
                    case SupersysWeight::FT_EFF_C_UP :
                        sys = NtSys::SusyNtSys::FT_EFF_C_systematics_UP;
                        break;
                    case SupersysWeight::FT_EFF_C_DN :
                        sys = NtSys::SusyNtSys::FT_EFF_C_systematics_DN;
                        break;
                    case SupersysWeight::FT_EFF_LT_UP :
                        sys = NtSys::SusyNtSys::FT_EFF_Light_systematics_UP;
                        break;
                    case SupersysWeight::FT_EFF_LT_DN :
                        sys = NtSys::SusyNtSys::FT_EFF_Light_systematics_DN;
                        break;
                    case SupersysWeight::FT_EFF_EXTRAP_UP :
                        sys = NtSys::SusyNtSys::FT_EFF_extrapolation_UP;
                        break; 
                    case SupersysWeight::FT_EFF_EXTRAP_DN :
                        sys = NtSys::SusyNtSys::FT_EFF_extrapolation_DN;
                        break; 
                    case SupersysWeight::FT_EFF_EXTRAPC_UP :
                        sys = NtSys::SusyNtSys::FT_EFF_extrapolation_charm_UP;
                        break;
                    case SupersysWeight::FT_EFF_EXTRAPC_DN :
                        sys = NtSys::SusyNtSys::FT_EFF_extrapolation_charm_DN;
                        break;
                    default: break;
                } // switch
            }
            if(jets.size()>0) { btagSF = m_nttools.bTagSFError(jets, sys); }
            weights_->btagSf = btagSF;
        }

    } //isMC
    return true;
}
///////////////////////////////////////////////////////////////////////////////
double Superflow::computeLeptonEfficiencySf(const Susy::Lepton &lep, const SupersysWeight sys)
{
    double out_SF = 1.0;
    double sf = 1.0;
    double delta = 0.0;

    if (lep.isEle()) {
        const Electron* el = dynamic_cast<const Electron*> (&lep);

        // get the ElectronId for the analysis type to select the correct errSF
        ElectronId id = m_nttools.electronSelector().signalId();

        // select the nominal SF
        sf = el->eleEffSF[id];

        // not storing charge id selection SF
        //if ( sys == SupersysWeight::EL_EFF_ChargeIDSel_DN ) {
        //}
        //else if( sys == SupersysWeight::EL_EFF_ChargeIDSel_UP ) {
        //}
        if( sys == SupersysWeight::EL_EFF_ID_TOTAL_Uncorr_DN ) {
            delta = el->errEffSF_id_dn[id];
        }
        else if( sys == SupersysWeight::EL_EFF_ID_TOTAL_Uncorr_UP ) {
            delta = el->errEffSF_id_up[id];
        }
        else if( sys == SupersysWeight::EL_EFF_Iso_TOTAL_Uncorr_DN ) {
            delta = el->errEffSF_iso_dn[id];
        }
        else if( sys == SupersysWeight::EL_EFF_Iso_TOTAL_Uncorr_UP ) {
            delta = el->errEffSF_iso_up[id];
        }
        else if( sys == SupersysWeight::EL_EFF_Reco_TOTAL_Uncorr_DN ) {
            delta = el->errEffSF_reco_dn[id];
        }
        else if( sys == SupersysWeight::EL_EFF_Reco_TOTAL_Uncorr_UP ) {
            delta = el->errEffSF_reco_up[id];
        }
        //else if( sys == SupersysWeight::EL_EFF_TriggerEff_TOTAL_DN ) {
        //}
        //else if( sys == SupersysWeight::EL_EFF_TriggerEff_TOTAL_UP ) {
        //}
        //else if( sys == SupersysWeight::EL_EFF_Trigger_TOTAL_DN ) {
        //}
        //else if( sys == SupersysWeight::EL_EFF_Trigger_TOTAL_UP ) {
        //}
    } // isEle
    else if (lep.isMu()) {
        const Muon* mu = dynamic_cast<const Muon*> (&lep);

        // get the MuonId for the analysis type to select the correct errSF
        MuonId id = m_nttools.muonSelector().signalId();

        // select the nominal SF
        sf = mu->muoEffSF[id];

        if ( sys == SupersysWeight::MUON_EFF_BADMUON_STAT_DN ) {
            delta = mu->errEffSF_badmu_stat_dn[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_BADMUON_STAT_UP ) {
            delta = mu->errEffSF_badmu_stat_up[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_BADMUON_SYS_DN ) {
            delta = mu->errEffSF_badmu_syst_dn[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_BADMUON_SYS_UP ) {
            delta = mu->errEffSF_badmu_syst_up[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_ISO_STAT_DN ) {
            delta = mu->errEffSF_iso_stat_dn[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_ISO_STAT_UP ) {
            delta = mu->errEffSF_iso_stat_up[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_ISO_SYS_DN ) {
            delta = mu->errEffSF_iso_syst_dn[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_ISO_SYS_UP ) {
            delta = mu->errEffSF_iso_syst_up[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_RECO_STAT_DN ) {
            delta = mu->errEffSF_reco_stat_dn[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_RECO_STAT_UP ) {
            delta = mu->errEffSF_reco_stat_up[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_RECO_SYS_DN ) {
            delta = mu->errEffSF_reco_syst_dn[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_RECO_SYS_UP ) {
            delta = mu->errEffSF_reco_syst_up[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_RECO_STAT_LOWPT_DN ) {
            delta = mu->errEffSF_reco_lowpt_stat_dn[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_RECO_STAT_LOWPT_UP ) {
            delta = mu->errEffSF_reco_lowpt_stat_up[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_RECO_SYS_LOWPT_DN ) {
            delta = mu->errEffSF_reco_lowpt_syst_dn[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_RECO_SYS_LOWPT_UP ) {
            delta = mu->errEffSF_reco_lowpt_syst_up[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_TTVA_STAT_DN ) {
            delta = mu->errEffSF_ttva_stat_dn[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_TTVA_STAT_UP ) {
            delta = mu->errEffSF_ttva_stat_up[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_TTVA_SYS_DN ) {
            delta = mu->errEffSF_ttva_syst_dn[id];
        }
        else if( sys == SupersysWeight::MUON_EFF_TTVA_SYS_UP ) {
            delta = mu->errEffSF_ttva_syst_up[id];
        }
        //else if( sys == SupersysWeight::MUON_EFF_TrigStat_DN ) {
        //    delta = mu->errEffSF_trig_stat_dn[id];
        //}
        //else if( sys == SupersysWeight::MUON_EFF_TrigStat_UP ) {
        //    delta = mu->errEffSF_trig_stat_up[id];
        //}
        //else if( sys == SupersysWeight::MUON_EFF_TrigSys_DN ) {
        //    delta = mu->errEffSF_trig_syst_dn[id];
        //}
        //else if( sys == SupersysWeight::MUON_EFF_TrigSys_UP ) {
        //    delta = mu->errEffSF_trig_syst_up[id];
        //}
    } // isMu

    out_SF = (sf + delta);
    return out_SF;
}
///////////////////////////////////////////////////////////////////////////////
void Superflow::setLumi(const float lumi) {
    m_luminosity = lumi;
    cout << app_name << "------------------------------ ------------------------------" << endl;
    cout << app_name << " Setting MC normalization (luminosity) to " << m_luminosity << " pb^-1." << endl;
    cout << app_name << "----------------------------- -------------------------------" << endl;
}
///////////////////////////////////////////////////////////////////////////////
void Superflow::setCountWeights(bool value) ///> public function, if set true it prints the weighted cutflow
{
    m_countWeights = value;
}
///////////////////////////////////////////////////////////////////////////////
void Superflow::setRunMode(SuperflowRunMode run_mode_) ///> public function
{
    m_runMode = run_mode_;
}
///////////////////////////////////////////////////////////////////////////////
void Superflow::setSingleEventSyst(SusyNtSys nt_syst_)
{
    m_singleEventSyst = nt_syst_;
}
///////////////////////////////////////////////////////////////////////////////
void Superflow::setChain(TChain* input_chain_)
{
    m_input_chain = input_chain_;
}
///////////////////////////////////////////////////////////////////////////////
void Superflow::setFileSuffix(string suffix)
{
    m_outputFileNameSuffix = suffix;
}

} // namespace sflow
