// SuperflowAna.cxx
//

// std
#include <cstdlib>
#include <cmath>
#include <fstream> 
#include <iostream>
#include <string>
#include <getopt.h>

// ROOT
#include "TChain.h"
#include "TVectorD.h"

// SusyNtuple
#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"
#include "SusyNtuple/SusyNtSys.h"
#include "SusyNtuple/TriggerTools.h"
#include "SusyNtuple/KinematicTools.h"

// Superflow
#include "Superflow/Superflow.h"
#include "Superflow/Superlink.h"
#include "Superflow/Cut.h"
#include "Superflow/StringTools.h"
#include "Superflow/PhysicsTools.h"
#include "Superflow/LeptonTruthDefinitions.h"

// RestFrames
#include "RestFrames/RestFrames.hh"


using namespace std;
using namespace sflow;
using namespace RestFrames;

// function prototypes
void print_usage(const char *exeName);
void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_, SuperflowRunMode& run_mode_, SusyNtSys& nt_sysnt_);

const string analysis_name = "wwlikeAna";

int main(int argc, char* argv[])
{
    // START read-in
    int n_skip_ = 0;
    int num_events_ = -1;
    string sample_;
    SuperflowRunMode run_mode = SuperflowRunMode::nominal;
    SusyNtSys nt_sys_ = NtSys::NOM;

    TChain* chain = new TChain("susyNt");
    chain->SetDirectory(0);
    
    ////////////////////////////////////////////////////////////
    // Read in the command-line options (input file, num events, etc...)
    ////////////////////////////////////////////////////////////
    read_options(argc, argv, chain, n_skip_, num_events_, sample_, run_mode, nt_sys_); // defined below

    ////////////////////////////////////////////////////////////
    // Initialize & configure the analysis
    //  > Superflow inherits from SusyNtAna : TSelector
    ////////////////////////////////////////////////////////////
    Superflow* cutflow = new Superflow(); // initialize the cutflow
    cutflow->setAnaName(analysis_name);
    cutflow->setAnaType(AnalysisType::Ana_Stop2L); 
    cutflow->setLumi(523.3); // set the MC normalized to lumi periods A1-A3
    cutflow->setSampleName(sample_);
    cutflow->setRunMode(run_mode);
    cutflow->setCountWeights(true); // print the weighted cutflows
    cutflow->setChain(chain);

    cout << analysis_name <<"    Total Entries   : " << chain->GetEntries() << endl;
    if(num_events_ > 0) {
    cout << analysis_name <<"    Process Entries : " << num_events_ << endl;
    } else {
    cout << analysis_name <<"    Process Entries : " << chain->GetEntries() << endl;
    }

    if (run_mode == SuperflowRunMode::single_event_syst) cutflow->setSingleEventSyst(nt_sys_);

  //  TriggerTools* triggers = new TriggerTools(chain, true);

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    //
    //  Superflow methods [BEGIN]
    //
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    
    *cutflow << CutName("read in") << [](Superlink* sl) -> bool { return true; };

    ////////////////////////////////////////////////////////////
    //  Cleaning Cuts
    ////////////////////////////////////////////////////////////
    
    int cutflags = 0;
    
    *cutflow << CutName("Pass GRL") << [&](Superlink* sl) -> bool {
        cutflags = sl->nt->evt()->cutFlags[sl->nt_sys];
        return (sl->tools->passGRL(cutflags));
    };
    
    *cutflow << CutName("LAr error") << [&](Superlink* sl) -> bool {
        return (sl->tools->passLarErr(cutflags));
    };
    
    *cutflow << CutName("Tile error") << [&](Superlink* sl) -> bool {
        return (sl->tools->passTileErr(cutflags));
    };
    
    *cutflow << CutName("TTC veto") << [&](Superlink* sl) -> bool {
        return (sl->tools->passTTC(cutflags));
    };

    *cutflow << CutName("pass Good Vertex") << [&](Superlink* sl) -> bool {
        return (sl->tools->passGoodVtx(cutflags));
    };
    
    *cutflow << CutName("pass bad muon veto") << [&](Superlink* sl) -> bool {
        return (sl->tools->passBadMuon(sl->preMuons));
    };
    
    *cutflow << CutName("pass cosmic muon veto") << [&](Superlink* sl) -> bool {
        return (sl->tools->passCosmicMuon(sl->baseMuons));
    };
    
    *cutflow << CutName("pass jet cleaning") << [&](Superlink* sl) -> bool {
        return (sl->tools->passJetCleaning(sl->baseJets));
    };

    ////////////////////////////////////////////////////////////
    //  Analysis Cuts
    ////////////////////////////////////////////////////////////
    
    *cutflow << CutName("tau veto") << [](Superlink* sl) -> bool {
        return sl->taus->size() == 0;
    };
    *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
        return sl->leptons->size() == 2;
    };

    ////////////////////////////////////////////////////////////
    //  Output Ntuple Setup
    //      > Ntuple variables
    ////////////////////////////////////////////////////////////

    *cutflow << NewVar("event weight"); {
        *cutflow << HFTname("eventweight");
        *cutflow << [](Superlink* sl, var_double*) -> double { 
            return sl->weights->product();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mcChannel (dsid)"); {
        *cutflow << HFTname("dsid");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            if (sl->isMC) {                        
                return sl->nt->evt()->mcChannel;
            }
            else{
                return 0.0;
            }
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Monte-Carlo generator event weight"); {
        *cutflow << HFTname("w");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->w; 
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Pile-up weight"); {
        *cutflow << HFTname("pupw");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->wPileup;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Event run number"); {
        *cutflow << HFTname("runNumber");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->run; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Event number"); {
        *cutflow << HFTname("eventNumber");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->eventNumber; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is Monte Carlo"); {
        *cutflow << HFTname("isMC");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->nt->evt()->isMC ? true : false; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of primary vertices"); {
        *cutflow << HFTname("nVtx");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->nVtx; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("average interactions per b.c."); {
        *cutflow << HFTname("avgMu");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->nt->evt()->avgMu; };
        *cutflow << SaveVar();
    }


    /* ------------------------------------------------------------- */
    /*  Leptons     leptons                                          */
    /*  Leptons     leptons                                          */
    /*  Leptons     leptons                                          */
    /* ------------------------------------------------------------- */

    LeptonVector leptons;
    ElectronVector electrons;
    MuonVector muons;
    *cutflow << [&](Superlink* sl, var_void*) { leptons = *sl->leptons; };
    *cutflow << [&](Superlink* sl, var_void*) { electrons = *sl->electrons; };
    *cutflow << [&](Superlink* sl, var_void*) { muons = *sl->muons; };

    *cutflow << NewVar("number of leptons"); {
        *cutflow << HFTname("nLeptons");
        *cutflow << [&](Superlink* sl, var_int*) -> int { return leptons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of electrons"); {
        *cutflow << HFTname("nElectrons");
        *cutflow << [&](Superlink* sl, var_int*) -> int {return electrons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of muons"); {
        *cutflow << HFTname("nMuons");
        *cutflow <<[&](Superlink* sl, var_int*) -> int {return muons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton pt"); {
        *cutflow << HFTname("l_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> lep_pt;
            for(int i = 0; i< leptons.size(); i++) {
                lep_pt.push_back(leptons.at(i)->Pt());
            }
            return lep_pt;
            };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton eta"); {
        *cutflow << HFTname("l_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> lep_eta;
            for(int i = 0; i < leptons.size(); i++) {
                lep_eta.push_back(leptons.at(i)->Eta());
            }
            return lep_eta;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton phi"); {
        *cutflow << HFTname("l_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> lep_phi;
            for(int i = 0; i < leptons.size(); i++) {
                lep_phi.push_back(leptons.at(i)->Phi());
            }
            return lep_phi;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton charge"); {
        *cutflow << HFTname("l_q");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->q);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton energy"); {
        *cutflow << HFTname("l_e");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->E());
            }
            return out;
            };
        *cutflow << SaveVar();
    }

    // muons
    *cutflow << NewVar("mll leptons"); {
        *cutflow << HFTname("mll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double mll = -999.0;
            if(leptons.size()==2) {
                Lepton* l0 = leptons.at(0);
                Lepton* l1 = leptons.at(1);
                mll = (*l0 + *l1).M();
            }
            return mll;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("dilepton transverse momentum"); {
        *cutflow << HFTname("pTll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double pTll = 0.0;
            if(leptons.size()==2) {
                Lepton* l0 = leptons.at(0);
                Lepton* l1 = leptons.at(1);
                pTll = (*l0 + *l1).Pt();
            } 
            return pTll;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta phi between two leptons"); {
        *cutflow << HFTname("dphi_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double dphi = -999.0;
            if(leptons.size()==2) {
                Lepton l0 = *leptons.at(0);
                Lepton l1 = *leptons.at(1);
                dphi = l0.DeltaPhi(l1);
            }
            return dphi;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta eta between two leptons"); {
        *cutflow << HFTname("deta_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double deta = -999.0;
            if(leptons.size()==2) {
                Lepton l0 = *leptons.at(0);
                Lepton l1 = *leptons.at(1);
                deta = l0.Eta() - l1.Eta();
            }
            return deta;
        };
        *cutflow << SaveVar();
    }

    // JETS
    // JETS
    // JETS
    JetVector jets;
    JetVector bjets;
    JetVector sjets;

    *cutflow << [&](Superlink* sl, var_void*) { jets = *sl->jets; };

    *cutflow << [&](Superlink* sl, var_void*) {
        for(int i = 0; i < jets.size(); i++) {
            Jet* j = jets[i];
            if( sl->tools->jetSelector().isB(j)) bjets.push_back(j);
            if(!sl->tools->jetSelector().isB(j)) sjets.push_back(j);
        } // i
    };
    *cutflow << NewVar("number of jets"); {
        *cutflow << HFTname("nJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return jets.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of sjets"); {
        *cutflow << HFTname("nSJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return sjets.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of bjets"); {
        *cutflow << HFTname("nBJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return bjets.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("jet pt"); {
        *cutflow << HFTname("j_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->Pt());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet pt"); {
        *cutflow << HFTname("sj_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->Pt());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet pt"); {
        *cutflow << HFTname("bj_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->Pt());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("jet eta"); {
        *cutflow << HFTname("j_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->Eta());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet eta"); {
        *cutflow << HFTname("sj_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->Eta());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet eta"); {
        *cutflow << HFTname("bj_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->Eta());
            }
            return out;
            };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet phi"); {
        *cutflow << HFTname("j_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->Phi());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet phi"); {
        *cutflow << HFTname("sj_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->Phi());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet phi"); {
        *cutflow << HFTname("bj_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->Phi());
            }
            return out;
            };
        *cutflow << SaveVar();
    }

    // MET 
    // MET 
    // MET 
    Met met;
    *cutflow << [&](Superlink* sl, var_void*) { met = *sl->met; };
    
    *cutflow << NewVar("transverse missing energy (Etmiss)"); {
        *cutflow << HFTname("met");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return met.lv().Pt(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("phi coord. of Etmiss"); {
        *cutflow << HFTname("metPhi");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return met.lv().Phi(); };
        *cutflow << SaveVar();
    }
    

    ////////////////////////////////////////
    // WW-like analysis variables
    ////////////////////////////////////////

    *cutflow << NewVar("cosine(Theta_b)"); {
        *cutflow << HFTname("cosThetaB");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double cosThetaB = -999.0;
            if(leptons.size()==2){
                TLorentzVector lp;
                TLorentzVector lm;
                for(unsigned int iL = 0; iL < leptons.size(); iL++){
                    if(leptons.at(iL)->q < 0) lm = *leptons.at(iL);
                    else if (leptons.at(iL)->q > 0) lp = *leptons.at(iL);
                }
                TLorentzVector ll = lp + lm;
                TVector3 boost = ll.BoostVector();
                lp.Boost(-boost);
                lm.Boost(-boost);
                cosThetaB = tanh((lp.Eta()-lm.Eta())/2.);
            }
            return cosThetaB;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("stransverse mass (mT2)"); {
        *cutflow << HFTname("mt2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double mt2 = -999.0;
            if(leptons.size()==2) {

                mt2 = kin::getMT2(*sl->leptons, *sl->met);


                //Mt2_bisect::mt2 mt2_event;
                //Double *pa, *pb, *pmiss;
                //Pa = new double[3]; pa[0] = sl->leptons->at(0)->M(); pa[1] = sl->leptons->at(0)->Px(), pa[2] = sl->leptons->at(0)->Py();
                //Pb = new double[3]; pb[0] = sl->leptons->at(1)->M(); pb[1] = sl->leptons->at(1)->Px(), pb[2] = sl->leptons->at(1)->Py();
                //Pmiss = new double[3]; pmiss[0] = 0.0; pmiss[1] = sl->met->Et * cos(sl->met->phi); pmiss[2] = sl->met->Et * sin(sl->met->phi);
 
                //Mt2_event.set_momenta(pa, pb, pmiss);
                //Mt2_event.set_mn(0.0); // LSP mass = 0 is Generic
 
                //Mt2 = mt2_event.get_mt2();
 
                //Delete[] pa;
                //Delete[] pb;
                //Delete[] pmiss;
            }
            return mt2;
        };
        *cutflow << SaveVar();
    }
    double meff;
    *cutflow << NewVar("meff : scalar sum pt of all jets, leptons, and met"); {
        *cutflow << HFTname("meff");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            meff = 0.0;
            // met
            meff += met.lv().Pt();
            // jets
            for(unsigned int ij=0; ij<jets.size(); ij++){
                meff += jets.at(ij)->Pt();
            }
            // leptons
            for(unsigned int il=0; il<leptons.size(); il++){
                meff += leptons.at(il)->Pt();
            }
            return meff;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("R1: met/meff"); {
        *cutflow << HFTname("R1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double R1 = -99.0;
            if (meff>0.0) {
                R1 = met.lv().Pt() / meff * 1.0;
            }
            return R1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("R2: met / (met + l0pt + l1pt)"); {
        *cutflow << HFTname("R2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double R2 = -99.0;
            if(leptons.size()==2){
                double denom = met.lv().Pt() + leptons.at(0)->Pt() + leptons.at(1)->Pt();
                R2 = met.lv().Pt() / denom * 1.0;
            }
            return R2;
        };
        *cutflow << SaveVar();
    }

    // dilepton super-razor variables
    double MDR, shatr, cosThetaRp1, DPB, dphi_l1_l2, gamma_r;
    double dphi_vBeta_R_vBeta_T;
    TVector3 vBeta_z, pT_CM, vBeta_T_CMtoR, vBeta_r;
    *cutflow << NewVar("Super-razor variables -- shatr"); {
        *cutflow << HFTname("shatr");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            MDR = shatr = cosThetaRp1 = DPB = dphi_l1_l2 = gamma_r = -999.0;
            dphi_vBeta_R_vBeta_T = -999.0;
            if(leptons.size()==2){
                kin::superRazor(leptons, &met, vBeta_z, pT_CM,
                                  vBeta_T_CMtoR, vBeta_r, shatr, DPB,
                                  dphi_l1_l2, gamma_r, dphi_vBeta_R_vBeta_T,
                                  MDR, cosThetaRp1);
            }
            return shatr;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("Super-razor variables -- dphi_l1_l2"); {
        *cutflow << HFTname("dphi_l1_l2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return dphi_l1_l2;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("Super-razor variables -- MDR"); {
        *cutflow << HFTname("MDR");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("Super-razor variables -- CosThetaR+1"); {
        *cutflow << HFTname("cosThetaRp1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return cosThetaRp1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("Super-razor variables -- DPB"); {
        *cutflow << HFTname("DPB");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPB;
        };
        *cutflow << SaveVar();
    } 
    ///////////////////////////////////////////////////////////////////////////////////////////
    //  RestFrames - 0 b (no requirements on b) -- SUPERRAZOR CASE
    ///////////////////////////////////////////////////////////////////////////////////////////

    double shat_0;
    double pTT_t;
    double pTT_z;
    double RPT_0;
    double RPZ_0;
    double gamma_rp1_0;
    double MDR_i1_t1_0;
    double MDR_i1_t2_0;
    double MDR_i2_t1_0;
    double MDR_i2_t2_0;
    double MDR_v1_t1_0;
    double MDR_v1_t2_0;
    double MDR_v2_t2_0;
    double MDR_v2_t1_0;
    double COS_tt_0;
    double COS_t1_0;
    double COS_t2_0;
    double DPD_tt_t1_0;
    double DPD_tt_t2_0;
    double DPD_t1_tt_0;
    double DPD_t2_tt_0;
    double DPD_t1_t2_0;
    double DPD_lab_tt_0;
    double DPB_vTT_0;
    double v_i1_0;
    

    *cutflow << [&](Superlink* sl, var_void*) {
        LabRecoFrame lab("lab", "lab");
        DecayRecoFrame tt("tt", "#tilde{t} #tilde{#bar{t}}");
        DecayRecoFrame t1("t1", "#tilde{t}_{1}");
        DecayRecoFrame t2("t2", "#tilde{t}_{2}");
        VisibleRecoFrame v1("v1", "vis_{1}");
        VisibleRecoFrame v2("v2", "vis_{2}");
        InvisibleRecoFrame i1("i1", "invis_{1}");
        InvisibleRecoFrame i2("i2", "invis_{2}");


        //
        /// connect the frames
        //
        lab.SetChildFrame(tt);
        tt.AddChildFrame(t1);
        tt.AddChildFrame(t2);
        t1.AddChildFrame(v1);
        t1.AddChildFrame(i1);
        t2.AddChildFrame(v2);
        t2.AddChildFrame(i2);

        // check that the decay tree is connected correctly
        if(!lab.InitializeTree()) {
            cout << "RestFrames::InitializeTree ERROR   Unable to initialize tree from lab frame. Exitting." << endl;
            exit(1);
        }

        //
        /// Define groups
        //

        // invisible group
        InvisibleGroup inv("inv", "wimp jigsaws");
        inv.AddFrame(i1);
        inv.AddFrame(i2);

        // combinatoric group is list of visible objects that end up in our hemispheres
        CombinatoricGroup vis("vis", "visible object jigsaws");
        vis.AddFrame(v1);
        vis.SetNElementsForFrame(v1,1,false);
        vis.AddFrame(v2);
        vis.SetNElementsForFrame(v2,1,false);


        SetMassInvJigsaw MinMassJigsaw("MinMass", "Invisible system mass jigsaw");
        inv.AddJigsaw(MinMassJigsaw);

        SetRapidityInvJigsaw RapidityJigsaw("RapidityJigsaw", "Invisible system rapidity jigsaw");
        inv.AddJigsaw(RapidityJigsaw);
        RapidityJigsaw.AddVisibleFrames(lab.GetListVisibleFrames());

        ContraBoostInvJigsaw ContraBoostJigsaw("ContraBoostJigsaw", "ContraBoost Invariant Jigsaw");
        inv.AddJigsaw(ContraBoostJigsaw);
        ContraBoostJigsaw.AddVisibleFrames((t1.GetListVisibleFrames()), 0);
        ContraBoostJigsaw.AddVisibleFrames((t2.GetListVisibleFrames()), 1);
        ContraBoostJigsaw.AddInvisibleFrame(i1, 0);
        ContraBoostJigsaw.AddInvisibleFrame(i2, 1);

        MinMassesCombJigsaw HemiJigsaw("hemi_jigsaw", "Minimize m_{vis_{1,2}} jigsaw");
        vis.AddJigsaw(HemiJigsaw);
        HemiJigsaw.AddFrame(v1,0);
        HemiJigsaw.AddFrame(v2,1);
        

        // check that the jigsaws are in place
        if(!lab.InitializeAnalysis()) {
            cout << "RestFrames::InitializeAnalysis ERROR   Unable to initialize analysis from labframe. Exitting." << endl;
            exit(1);
        }

        lab.ClearEvent();

        // set the met
        TVector3 met3vector(sl->met->lv().Px(), sl->met->lv().Py(), sl->met->lv().Pz());
        inv.SetLabFrameThreeVector(met3vector);

        // add the leptons to the vis group
        vis.AddLabFrameFourVector(*leptons.at(0));
        vis.AddLabFrameFourVector(*leptons.at(1));

        // analyze the event
        lab.AnalyzeEvent();

        ///////// MTT
        shat_0 = tt.GetMass();

        //////////////////////////////////
        // RATIO OF CM pT to CM mass
        TVector3 vPTT_0 = tt.GetFourVector(lab).Vect();
        RPT_0 = vPTT_0.Pt() / (vPTT_0.Pt() + shat_0/4.);
        // RATIO OF CM pZ to CM mass
        RPZ_0 = vPTT_0.Pz() / (vPTT_0.Pz() + shat_0/4.);

        pTT_t = vPTT_0.Pt();
        pTT_z = vPTT_0.Pz();

        ///////// shapes
        gamma_rp1_0 = tt.GetVisibleShape();

        ///////// energies
        MDR_i1_t1_0 = i1.GetEnergy(t1);
        MDR_i1_t2_0 = i1.GetEnergy(t2);
        MDR_i2_t1_0 = i2.GetEnergy(t1);
        MDR_i2_t2_0 = i2.GetEnergy(t2);

        MDR_v1_t1_0 = 2.0 * v1.GetEnergy(t1);
        MDR_v1_t2_0 = 2.0 * v1.GetEnergy(t2);
        MDR_v2_t2_0 = 2.0 * v2.GetEnergy(t2);
        MDR_v2_t1_0 = 2.0 * v2.GetEnergy(t1);

        ////////// angles
        // decay angle of tt (~cosThetaB)
        COS_tt_0 = tt.GetCosDecayAngle();
        COS_t1_0 = t1.GetCosDecayAngle();
        COS_t2_0 = t2.GetCosDecayAngle();
        

        // dphi between tt and t1 decay planes
        DPD_tt_t1_0 = tt.GetDeltaPhiDecayPlanes(t1);
        DPD_tt_t2_0 = tt.GetDeltaPhiDecayPlanes(t2);
        // dphi between lab and tt decay planes
        DPD_lab_tt_0 = lab.GetDeltaPhiDecayPlanes(tt);
        
        //////////////////////////////////////
        // BOOST ANGLES
        // delta phi between tt visible decay products and tt momentum
        DPB_vTT_0 = tt.GetDeltaPhiBoostVisible();

        //////////////////////////////////////
        // DECAY PLANES
        // delta phi between t1 decay frame and tt decay frame
        DPD_t1_tt_0 = t1.GetDeltaPhiDecayPlanes(tt);
        // delta phi between t2 decay frame and tt decay frame
        DPD_t2_tt_0 = t2.GetDeltaPhiDecayPlanes(tt);
        // delta phi between t1 and t2 decay frames in lab frame
        DPD_t1_t2_0 = t1.GetDeltaPhiDecayPlanes(t2);

        // "velocity" of invisible
        const RestFrame& i1_prdframe = i1.GetProductionFrame();
        v_i1_0 = i1.GetMomentum(i1_prdframe) / t1.GetMass();

    }; // end var_void

    *cutflow << NewVar("RFSR -- shat_0"); {
        *cutflow << HFTname("shat_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return shat_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- pTT_t"); {
        *cutflow << HFTname("pTT_t_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return pTT_t;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- pTT_z"); {
        *cutflow << HFTname("pTT_z_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return pTT_z;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- RPT_0"); {
        *cutflow << HFTname("RPT_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPT_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- RPZ_0"); {
        *cutflow << HFTname("RPZ_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPZ_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- gamma_rp1_0"); {
        *cutflow << HFTname("gamma_rp1_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return gamma_rp1_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- MDR_i1_t1_0"); {
        *cutflow << HFTname("MDR_i1_t1_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR_i1_t1_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- MDR_i1_t2_0"); {
        *cutflow << HFTname("MDR_i1_t2_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR_i1_t2_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- MDR_i2_t1_0"); {
        *cutflow << HFTname("MDR_i2_t1_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR_i2_t1_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- MDR_i2_t2_0"); {
        *cutflow << HFTname("MDR_i2_t2_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR_i2_t2_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- MDR_v1_t1_0"); {
        *cutflow << HFTname("MDR_v1_t1_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR_v1_t1_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- MDR_v1_t2_0"); {
        *cutflow << HFTname("MDR_v1_t2_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR_v1_t2_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- MDR_v2_t2_0"); {
        *cutflow << HFTname("MDR_v2_t2_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR_v2_t2_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- MDR_v2_t1_0"); {
        *cutflow << HFTname("MDR_v2_t1_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR_v2_t1_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- COS_tt_0"); {
        *cutflow << HFTname("COS_tt_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return COS_tt_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- COS_t1_0"); {
        *cutflow << HFTname("COS_t1_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return COS_t1_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- COS_t2_0"); {
        *cutflow << HFTname("COS_t2_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return COS_t2_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- DPD_tt_t1_0"); {
        *cutflow << HFTname("DPD_tt_t1_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPD_tt_t1_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- DPD_tt_t2_0"); {
        *cutflow << HFTname("DPD_tt_t2_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPD_tt_t2_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- DPD_t1_tt_0"); {
        *cutflow << HFTname("DPD_t1_tt_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPD_t1_tt_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- DPD_t2_tt_0"); {
        *cutflow << HFTname("DPD_t2_tt_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPD_t2_tt_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- DPD_t1_t2_0"); {
        *cutflow << HFTname("DPD_t1_t2_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPD_t1_t2_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- DPD_lab_tt_0"); {
        *cutflow << HFTname("DPD_lab_tt_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPD_lab_tt_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- DPB_vTT_0"); {
        *cutflow << HFTname("DPB_vTT_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPB_vTT_0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RFSR -- v_i1_0"); {
        *cutflow << HFTname("v_i1_0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return v_i1_0;
        };
        *cutflow << SaveVar();
    }




    ///////////////////////////////////////////////////////////////////////////////////////////
    //  RestFrames - 1 b 
    ///////////////////////////////////////////////////////////////////////////////////////////

    // rf variables
    double shat;
    double RPT;
    double RPZ;
    double dphi_b_w1_t1;
    double dphi_b_w2_t1;
    double dphi_b_w1_t2;
    double dphi_b_w2_t2;
    double DPB_vTT;
    double DPB_vT1;
    double DPB_vT2;
    double DPD_t1_tt;
    double DPD_t2_tt;
    double DPD_t1_t2;
    double COS_t1_lab;
    double COS_t2_lab;
    double COS_t1_tt;
    double COS_t2_tt;
    double COS_tt_lab;
    double MDR_l1_t1;
    double MDR_l2_t2;
    double MDR_b1_t1;
    double MDR_l1_w1;
    double MDR_l2_w2;
    double MTT;

    *cutflow << [&](Superlink* sl, var_void*) {
        if(bjets.size()>0) {
        LabRecoFrame lab("lab", "lab");
        DecayRecoFrame tt("tt", "#tilde{t} #tilde{#bar{t}}");
        DecayRecoFrame t1("t1", "#tilde{t}_{1}");
        DecayRecoFrame t2("t2", "#tilde{t}_{2}");
        DecayRecoFrame w1("w1", "W_{1}");
        DecayRecoFrame w2("w2", "W_{2}");
        VisibleRecoFrame b1("b1", "b_{1}");
        //VisibleFrame b2("b2", "b_{2}");
        InvisibleRecoFrame c1("c1", "(#nu + #tilde{#chi})_{1}");
        InvisibleRecoFrame c2("c2", "(#nu + #tilde{#chi})_{2}");
        VisibleRecoFrame l1("l1", "#it{l}_{1}");
        VisibleRecoFrame l2("l2", "#it{l}_{2}");


        //
        /// connect the frames
        //
        lab.SetChildFrame(tt);
        tt.AddChildFrame(t1);
        tt.AddChildFrame(t2);
        t1.AddChildFrame(b1);
        t1.AddChildFrame(w1);
        //t2.AddChildFrame(b2);
        t2.AddChildFrame(w2);
        w1.AddChildFrame(c1);
        w1.AddChildFrame(l1);
        w2.AddChildFrame(c2);
        w2.AddChildFrame(l2);

        // check that the decay tree is connected correctly
        if(!lab.InitializeTree()) {
            cout << "RestFrames::InitializeTree ERROR   Unable to initialize tree from lab frame. Exitting." << endl;
            exit(1);
        }

        //
        /// Define groups
        //

        // invisible group
        InvisibleGroup inv("inv", "wimp jigsaws");
        inv.AddFrame(c1);
        inv.AddFrame(c2);

        // b-jet combinatoric group
        //CombinatoricGroup BJETS("bjets", "b-jet jigsaws");
        //BJETS.AddFrame(b1);
        //BJETS.SetNElementsForFrame(b1, 1, true);

        SetMassInvJigsaw MinMassJigsaw("MinMass", "Invisible system mass jigsaw");
        inv.AddJigsaw(MinMassJigsaw);

        SetRapidityInvJigsaw RapidityJigsaw("RapidityJigsaw", "Invisible system rapidity jigsaw");
        inv.AddJigsaw(RapidityJigsaw);
        RapidityJigsaw.AddVisibleFrames(lab.GetListVisibleFrames());

        ContraBoostInvJigsaw ContraBoostJigsaw("ContraBoostJigsaw", "ContraBoost Invariant Jigsaw");
        inv.AddJigsaw(ContraBoostJigsaw);
        ContraBoostJigsaw.AddVisibleFrames((t1.GetListVisibleFrames()), 0);
        ContraBoostJigsaw.AddVisibleFrames((t2.GetListVisibleFrames()), 1);
        ContraBoostJigsaw.AddInvisibleFrame(c1, 0);
        ContraBoostJigsaw.AddInvisibleFrame(c2, 1);

        //MinMassesCombJigsaw HemiJigsaw("HemiJigsaw", "minimize m(b #it{l}) jigsaw");
        //BJETS.AddJigsaw(HemiJigsaw);
        //HemiJigsaw.AddFrames((t1.GetListVisibleFrames()), 0);
        //HemiJigsaw.AddFrames((t2.GetListVisibleFrames()), 1);

        // check that the jigsaws are in place
        if(!lab.InitializeAnalysis()) {
            cout << "RestFrames::InitializeAnalysis ERROR   Unable to initialize analysis from labframe. Exitting." << endl;
            exit(1);
        }

        lab.ClearEvent();

        // set the met
        TVector3 met3vector(sl->met->lv().Px(), sl->met->lv().Py(), sl->met->lv().Pz());
        inv.SetLabFrameThreeVector(met3vector);

        // set the leptons
        l1.SetLabFrameFourVector(*leptons.at(0));
        l2.SetLabFrameFourVector(*leptons.at(1));

        // set the bjet
        b1.SetLabFrameFourVector(*bjets.at(0));


        // analyze the event
        lab.AnalyzeEvent();

        //////////////////////////////////
        // MASS IN CM FRAME
        shat = tt.GetMass();

        //////////////////////////////////
        // RATIO OF CM pT to CM mass
        TVector3 vPTT = tt.GetFourVector(lab).Vect();
        RPT = vPTT.Pt() / (vPTT.Pt() + shat/4.);
        // RATIO OF CM pZ to CM mass
        RPZ = vPTT.Pz() / (vPTT.Pz() + shat/4.);

        //////////////////////////////////
        // BJET ANGLES
        // angle between b-jet and W's in the t frame
        // -> t1
        TLorentzVector btlv_t1 = b1.GetFourVector(t1);
        TLorentzVector w1tlv_t1 = w1.GetFourVector(t1);
        TLorentzVector w2tlv_t1 = w2.GetFourVector(t1);
        dphi_b_w1_t1 = btlv_t1.DeltaPhi(w1tlv_t1);
        dphi_b_w2_t1 = btlv_t1.DeltaPhi(w2tlv_t1);
        // -> t2
        TLorentzVector btlv_t2 = b1.GetFourVector(t2);
        TLorentzVector w1tlv_t2 = w1.GetFourVector(t2);
        TLorentzVector w2tlv_t2 = w2.GetFourVector(t2);
        dphi_b_w1_t2 = btlv_t2.DeltaPhi(w1tlv_t2);
        dphi_b_w2_t2 = btlv_t2.DeltaPhi(w2tlv_t2);

  //      cout << "dphi_b_w1_t1 : " << dphi_b_w1_t1 << endl;
  //      cout << "dphi_b_w2_t1 : " << dphi_b_w2_t1 << endl;
  //      cout << "dphi_b_w1_t2 : " << dphi_b_w1_t2 << endl;
  //      cout << "dphi_b_w2_t2 : " << dphi_b_w2_t2 << endl;


        //////////////////////////////////////
        // BOOST ANGLES
        // delta phi between tt visible decay products and tt momentum
        DPB_vTT = tt.GetDeltaPhiBoostVisible();
        // delta phi between t1 visible decay products and t1 momentum
        DPB_vT1 = t1.GetDeltaPhiBoostVisible();
        // delta phi between t2 visible decay products and t2 momentum
        DPB_vT2 = t2.GetDeltaPhiBoostVisible();

   //     cout << "DPB_vTT : " << DPB_vTT << endl;
   //     cout << "DPB_vT1 : " << DPB_vT1 << endl;
   //     cout << "DPB_vT2 : " << DPB_vT2 << endl;


        //////////////////////////////////////
        // DECAY PLANES
        // delta phi between t1 decay frame and tt decay frame
        DPD_t1_tt = t1.GetDeltaPhiDecayPlanes(tt);
        // delta phi between t2 decay frame and tt decay frame
        DPD_t2_tt = t2.GetDeltaPhiDecayPlanes(tt);
        // delta phi between t1 and t2 decay frames in lab frame
        DPD_t1_t2 = t1.GetDeltaPhiDecayPlanes(t2);
        // cosine of t1 decay plane as seen in lab frame
        COS_t1_lab = t1.GetCosDecayAngle();
        COS_t2_lab = t2.GetCosDecayAngle();
        // cosine of t1 decay plane as seen in tt frame
        COS_t1_tt  = t1.GetCosDecayAngle(tt);
        COS_t2_tt  = t2.GetCosDecayAngle(tt);
        // cosine of tt decay plane as seen in lab frame
        COS_tt_lab = tt.GetCosDecayAngle();

    //    cout << "DPD_t1_tt : " << DPD_t1_tt << endl;
    //    cout << "DPD_t2_tt : " << DPD_t2_tt << endl;
    //    cout << "DPD_t1_t2 : " << DPD_t1_t2 << endl;
    //    cout << "COS_t1_lab : " << COS_t1_lab << endl;
    //    cout << "COS_t2_lab : " << COS_t2_lab << endl;
    //    cout << "COS_t1_tt : " << COS_t1_tt << endl;
    //    cout << "COS_t2_tt : " << COS_t2_tt << endl;
    //    cout << "COS_tt_lab : " << COS_tt_lab << endl;


        /////////////////////////////////////////
        // MASS VARIABLES
        // energy of l1 in t1 frame
        MDR_l1_t1 = 2.0 * l1.GetEnergy(t1);
        MDR_l2_t2 = 2.0 * l2.GetEnergy(t2);
        MDR_b1_t1 = 2.0 * b1.GetEnergy(t1);
        MDR_l1_w1 = 2.0 * l1.GetEnergy(w1);
        MDR_l2_w2 = 2.0 * l2.GetEnergy(w2);

        // mtt
        double PT = t1.GetMomentum(tt);
        double PB = b1.GetMomentum(t1);
        double EB = b1.GetEnergy(t1);
        double MW1 = 2. * l1.GetEnergy(w1);
        double MW2 = 2. * l2.GetEnergy(w2);
        double MT1 = EB + sqrt(PB*PB + MW1*MW1); 
        double MT2 = MW2;
        MTT = sqrt(MT1*MT1 + PT*PT) + sqrt(MT2*MT2 + PT*PT);

     //   cout << "MDR_l1_t1 : " << MDR_l1_t1 << endl;
     //   cout << "MDR_l2_t2 : " << MDR_l2_t2 << endl;
     //   cout << "MDR_b1_t1 : " << MDR_b1_t1 << endl;
     //   cout << "MTT       : " << MTT << endl;
     //   cout << endl;
        



    } // bjets size > 0
    };  /// end RestFrames
    *cutflow << NewVar("Super-razor variables -- dphi_l1_l2"); {
        *cutflow << HFTname("dphi_l1_l2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return dphi_l1_l2;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("RF -- shat"); {
        *cutflow << HFTname("shat");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return shat;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- RPT"); {
        *cutflow << HFTname("RPT");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPT;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- RPZ"); {
        *cutflow << HFTname("RPZ");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return RPZ;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- dphi_b_w1_t1"); {
        *cutflow << HFTname("dphi_b_w1_t1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return dphi_b_w1_t1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- dphi_b_w2_t1"); {
        *cutflow << HFTname("dphi_b_w2_t1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return dphi_b_w2_t1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- dphi_b_w1_t2"); {
        *cutflow << HFTname("dphi_b_w1_t2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return dphi_b_w2_t2;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- dphi_b_w2_t2"); {
        *cutflow << HFTname("dphi_b_w2_t2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return dphi_b_w2_t2;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- DPB_vTT"); {
        *cutflow << HFTname("DPB_vTT");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPB_vTT;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- DPB_vT1"); {
        *cutflow << HFTname("DPB_vT1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPB_vT1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- DPB_vT2"); {
        *cutflow << HFTname("DPB_vT2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPB_vT2;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- DPD_t1_tt"); {
        *cutflow << HFTname("DPD_t1_tt");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPD_t1_tt;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- DPD_t2_tt"); {
        *cutflow << HFTname("DPD_t2_tt");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPD_t2_tt;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- DPD_t1_t2"); {
        *cutflow << HFTname("DPD_t1_t2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return DPD_t1_t2;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- COS_t1_lab"); {
        *cutflow << HFTname("COS_t1_lab");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return COS_t1_lab;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- COS_t2_lab"); {
        *cutflow << HFTname("COS_t2_lab");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return COS_t2_lab;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- COS_t1_tt"); {
        *cutflow << HFTname("COS_t1_tt");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return COS_t1_tt;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- COS_t2_tt"); {
        *cutflow << HFTname("COS_t2_tt");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return COS_t2_tt;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- COS_tt_lab"); {
        *cutflow << HFTname("COS_tt_lab");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return COS_tt_lab;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- MDR_l1_t1"); {
        *cutflow << HFTname("MDR_l1_t1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR_l1_t1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- MDR_l2_t2"); {
        *cutflow << HFTname("MDR_l2_t2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR_l2_t2;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- MDR_b1_t1"); {
        *cutflow << HFTname("MDR_b1_t1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR_b1_t1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- MDR_l1_w1"); {
        *cutflow << HFTname("MDR_l1_w1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR_l1_w1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- MDR_l2_w2"); {
        *cutflow << HFTname("MDR_l2_w2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MDR_l2_w2;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("RF -- MTT"); {
        *cutflow << HFTname("MTT");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return MTT;
        };
        *cutflow << SaveVar();
    }









    *cutflow << [&](Superlink* sl, var_void*) { leptons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { electrons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { muons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { jets.clear(); }; 
    *cutflow << [&](Superlink* sl, var_void*) { bjets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { sjets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { met.clear(); };
    
    ////////////////////////////////////////////////////////////
    //  Output Ntuple Setup
    //      > Setup the output systematic ntuples
    ////////////////////////////////////////////////////////////

    //
    //// Weight variation systematics
    //

    *cutflow << NewSystematic("shift in data mu by +/- 10 %"); {
        *cutflow << WeightSystematic(SupersysWeight::PILEUP_UP, SupersysWeight::PILEUP_DOWN);
        *cutflow << TreeName("PILEUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in electron ID -SF"); {
        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_ID_UP, SupersysWeight::EL_EFF_ID_DOWN);
        *cutflow << TreeName("ELEFFID");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in electron reco - SF"); {
        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_RECO_UP, SupersysWeight::EL_EFF_RECO_DOWN);
        *cutflow << TreeName("ELEFFRECO");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in muon eff - stat"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_STAT_UP, SupersysWeight::MUON_EFF_STAT_DOWN);
        *cutflow << TreeName("MUEFFSTAT");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in muon eff - syst"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_SYST_UP, SupersysWeight::MUON_EFF_SYST_DOWN);
        *cutflow << TreeName("MUEFFSYST");
        *cutflow << SaveSystematic();
    }


/*
    //
    // Weight variation systematics
    //  > stored in nominal output tree
    //
    *cutflow << NewSystematic("positve shift due to background estimation method"); {
        *cutflow << WeightSystematic(SupersysWeight::BKGMETHODUP, SupersysWeight::BKGMETHODDOWN);
        *cutflow << TreeName("BKGMETHOD");    
        *cutflow << SaveSystematic();
    }


    *cutflow << NewSystematic("shift in electron trigger weights"); {
        *cutflow << WeightSystematic(SupersysWeight::ETRIGREWUP, SupersysWeight::ETRIGREWDOWN);
        *cutflow << TreeName("ETRIGREW");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in muon trigger weights"); {
        *cutflow << WeightSystematic(SupersysWeight::MTRIGREWUP, SupersysWeight::MTRIGREWDOWN);
        *cutflow << TreeName("MTRIGREW");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in b-tag scale factor"); {
        *cutflow << WeightSystematic(SupersysWeight::BJETUP, SupersysWeight::BJETDOWN);
        *cutflow << TreeName("BJET");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in c-tag scale factor"); {
        *cutflow << WeightSystematic(SupersysWeight::CJETUP, SupersysWeight::CJETDOWN);
        *cutflow << TreeName("CJET");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in light-tag scale factor"); {
        *cutflow << WeightSystematic(SupersysWeight::BMISTAGUP, SupersysWeight::BMISTAGDOWN);
        *cutflow << TreeName("BMISTAG");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in electron efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::ESFUP, SupersysWeight::ESFDOWN);
        *cutflow << TreeName("ESF");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in muon efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::MEFFUP, SupersysWeight::MEFFDOWN);
        *cutflow << TreeName("MEFF");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in pileup"); {
        *cutflow << WeightSystematic(SupersysWeight::PILEUPUP, SupersysWeight::PILEUPDOWN);
        *cutflow << TreeName("PILEUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in cross section"); {
        *cutflow << WeightSystematic(SupersysWeight::XSUP, SupersysWeight::XSDOWN);
        *cutflow << TreeName("XS");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in ISR uncertainty (MG5 scale variation)"); {
        *cutflow << WeightSystematic(SupersysWeight::ISRUP, SupersysWeight::ISRDOWN);
        *cutflow << TreeName("ISR");
        *cutflow << SaveSystematic();
    }


    //
    // Event/object systematics
    //  > stored in separate output ntuples
    //

    *cutflow << NewSystematic("Electron Scale Zsys + sigma"); {
        *cutflow << EventSystematic(NtSys::EES_Z_UP);
        *cutflow << TreeName("EESZUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Zsys - sigma"); {
        *cutflow << EventSystematic(NtSys::EES_Z_DN);
        *cutflow << TreeName("EESZDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Material + sigma"); {
        *cutflow << EventSystematic(NtSys::EES_MAT_UP);
        *cutflow << TreeName("EESMATUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Material - sigma"); {
        *cutflow << EventSystematic(NtSys::EES_MAT_DN);
        *cutflow << TreeName("EESMATDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Presampler + sigma"); {
        *cutflow << EventSystematic(NtSys::EES_PS_UP);
        *cutflow << TreeName("EESPSUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Presampler - sigma"); {
        *cutflow << EventSystematic(NtSys::EES_PS_DN);
        *cutflow << TreeName("EESPSDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Low Pt + sigma"); {
        *cutflow << EventSystematic(NtSys::EES_LOW_UP);
        *cutflow << TreeName("EESLOWUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Low Pt - sigma"); {
        *cutflow << EventSystematic(NtSys::EES_LOW_DN);
        *cutflow << TreeName("EESLOWDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Resolution + sigma"); {
        *cutflow << EventSystematic(NtSys::EER_UP);
        *cutflow << TreeName("EERUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Resolution - sigma"); {
        *cutflow << EventSystematic(NtSys::EER_DN);
        *cutflow << TreeName("EERDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Muon MS track + sigma"); {
        *cutflow << EventSystematic(NtSys::MS_UP);
        *cutflow << TreeName("MSUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Muon MS track - sigma"); {
        *cutflow << EventSystematic(NtSys::MS_DN);
        *cutflow << TreeName("MSDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Muon ID track + sigma"); {
        *cutflow << EventSystematic(NtSys::ID_UP);
        *cutflow << TreeName("IDUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Muon ID track - sigma"); {
        *cutflow << EventSystematic(NtSys::ID_DN);
        *cutflow << TreeName("IDDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Jet Energy Scale + sigma"); {
        *cutflow << EventSystematic(NtSys::JES_UP);
        *cutflow << TreeName("JESUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Jet Energy Scale - sigma"); {
        *cutflow << EventSystematic(NtSys::JES_DN);
        *cutflow << TreeName("JESDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Jet Energy Resolution (gaussian)"); {
        *cutflow << EventSystematic(NtSys::JER);
        *cutflow << TreeName("JER");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Met scale soft term + sigma"); {
        *cutflow << EventSystematic(NtSys::SCALEST_UP);
        *cutflow << TreeName("SCALESTUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Met scale soft term - sigma"); {
        *cutflow << EventSystematic(NtSys::SCALEST_DN);
        *cutflow << TreeName("SCALESTDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Met resolution soft term + sigma"); {
        *cutflow << EventSystematic(NtSys::RESOST);
        *cutflow << TreeName("RESOST");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Tau energy scale + sigma"); {
        *cutflow << EventSystematic(NtSys::TES_UP);
        *cutflow << TreeName("TESUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Tau energy scale - sigma"); {
        *cutflow << EventSystematic(NtSys::TES_DN);
        *cutflow << TreeName("TESDOWN");
        *cutflow << SaveSystematic();
    }

     *cutflow << NewSystematic("Jet JVF cut + sigma"); { // THIS systematic is erroring out.
         *cutflow << EventSystematic(NtSys::JVF_UP);
         *cutflow << TreeName("JVFUP");
         *cutflow << SaveSystematic();
     }
     
     *cutflow << NewSystematic("Jet JVF cut - sigma"); { // THIS systematic is erroring out.
         *cutflow << EventSystematic(NtSys::JVF_DN);
         *cutflow << TreeName("JVFDOWN");
         *cutflow << SaveSystematic();
     }
*/
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    //
    //  Superflow methods [END]
    //
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    // END Setup systematics
    // END Setup systematics
    // END Setup systematics

    // Initialize the cutflow and start the event loop.
    chain->Process(cutflow, sample_.c_str(), num_events_, n_skip_);

    delete cutflow;
    delete chain;

    cout << "Done." << endl;
    exit(0);
}

void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_,
    SuperflowRunMode& run_mode_, SusyNtSys& nt_sys)
{
    bool nominal_ = false;
    bool nominal_and_weight_syst_ = false;
    bool all_syst_ = false;
    bool single_event_syst_ = false;

    string systematic_ = "undefined";

    string input;

    /** Read inputs to program */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0)
            num_events_ = atoi(argv[++i]);
        else if (strcmp(argv[i], "-c") == 0)
            nominal_ = true;
        else if (strcmp(argv[i], "-w") == 0)
            nominal_and_weight_syst_ = true;
        else if (strcmp(argv[i], "-e") == 0) {
            single_event_syst_ = true;
        }
        else if (strcmp(argv[i], "-a") == 0) {
            all_syst_ = true;
        }
        else if (strcmp(argv[i], "-i") == 0) {
            input = argv[++i];
        }
        else if (strcmp(argv[i], "-s") == 0) {
            systematic_ = argv[++i];
        }
        else {
            cout << analysis_name <<"    Error (fatal): Bad arguments." << endl;
            exit(1);
        }
    }

    bool inputIsFile = Susy::utils::endswith(input, ".root");
    bool inputIsList = Susy::utils::endswith(input, ".txt");
    bool inputIsDir = Susy::utils::endswith(input, "/");
    bool validInput(inputIsFile || inputIsList || inputIsDir);
    if (!validInput) {
        cout << analysis_name <<"    invalid input '" << input << "'" << endl;
        exit(1);
    }
    if (inputIsFile) {
        ChainHelper::addFile(chain, input);
        cout << analysis_name <<"    file: " << input << endl;
        cout << analysis_name <<"    file: " << input << endl;
        cout << analysis_name <<"    file: " << input << endl;
        sample_ = input;
    }
    if (inputIsList) {
        ChainHelper::addFileList(chain, input);
        cout << analysis_name << "    list: " << input << endl;
        cout << analysis_name << "    list: " << input << endl;
        cout << analysis_name << "    list: " << input << endl;
        ifstream infile(input.c_str());
        if (infile.good()) {
            string sLine;
            getline(infile, sLine);
            sample_ = sLine;
        }
        else {
            sample_ = input;
        }
        infile.close();
    }
    if (inputIsDir) {
        ChainHelper::addFileDir(chain, input);
        cout << analysis_name << "    dir: " << input << endl;
        cout << analysis_name << "    dir: " << input << endl;
        cout << analysis_name << "    dir: " << input << endl;
        sample_ = input;
    }
    Long64_t tot_num_events = chain->GetEntries();
    num_events_ = (num_events_ < 0 ? tot_num_events : num_events_);
    // if (debug) chain->ls();

    if (nominal_) {
        run_mode_ = SuperflowRunMode::nominal;
        cout << analysis_name <<"    run mode: SuperflowRunMode::nominal" << endl;
    }
    if (nominal_and_weight_syst_) {
        run_mode_ = SuperflowRunMode::nominal_and_weight_syst;
        cout << analysis_name <<"    run mode: SuperflowRunMode::nominal_and_weight_syst" << endl;
    }
    if (single_event_syst_) {
        run_mode_ = SuperflowRunMode::single_event_syst;
        cout << analysis_name <<"    run mode: SuperflowRunMode::single_event_syst" << endl;
    }

    if (all_syst_) {
        run_mode_ = SuperflowRunMode::all_syst;
        cout << analysis_name <<"    run mode: SuperflowRunMode::all_syst" << endl;
    }

    map <string, SusyNtSys> event_syst_map;
/*    event_syst_map["EESZUP"] = NtSys::EES_Z_UP;
    event_syst_map["EESZDOWN"] = NtSys::EES_Z_DN;
    event_syst_map["EESMATUP"] = NtSys::EES_MAT_UP;
    event_syst_map["EESMATDOWN"] = NtSys::EES_MAT_DN;
    event_syst_map["EESPSUP"] = NtSys::EES_PS_UP;
    event_syst_map["EESPSDOWN"] = NtSys::EES_PS_DN;
    event_syst_map["EESLOWUP"] = NtSys::EES_LOW_UP;
    event_syst_map["EESLOWDOWN"] = NtSys::EES_LOW_DN;
    event_syst_map["EERUP"] = NtSys::EER_UP;
    event_syst_map["EERDOWN"] = NtSys::EER_DN;
    event_syst_map["MSUP"] = NtSys::MS_UP;
    event_syst_map["MSDOWN"] = NtSys::MS_DN;
    event_syst_map["IDUP"] = NtSys::ID_UP;
    event_syst_map["IDDOWN"] = NtSys::ID_DN;
    event_syst_map["JESUP"] = NtSys::JES_UP;
    event_syst_map["JESDOWN"] = NtSys::JES_DN;
    event_syst_map["JER"] = NtSys::JER;
    event_syst_map["SCALESTUP"] = NtSys::SCALEST_UP;
    event_syst_map["SCALESTDOWN"] = NtSys::SCALEST_DN;
    event_syst_map["RESOST"] = NtSys::RESOST;
    event_syst_map["TRIGSFELUP"] = NtSys::TRIGSF_EL_UP;
    event_syst_map["TRIGSFELDN"] = NtSys::TRIGSF_EL_DN;
    event_syst_map["TRIGSFMUUP"] = NtSys::TRIGSF_MU_UP;
    event_syst_map["TRIGSFMUDN"] = NtSys::TRIGSF_MU_DN;
    event_syst_map["TESUP"] = NtSys::TES_UP;
    event_syst_map["TESDOWN"] = NtSys::TES_DN;
    event_syst_map["JVFUP"] = NtSys::JVF_UP;
    event_syst_map["JVFDOWN"] = NtSys::JVF_DN;
*/
    if (single_event_syst_) {
        if (event_syst_map.count(systematic_) == 1) {
            nt_sys = event_syst_map[systematic_];
        }
        else {
            cout << analysis_name << " ERROR (fatal)    Event systematic option /s " << systematic_ << " -> not found." << endl;
            exit(1);
        }
    }
}
