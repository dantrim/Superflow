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


using namespace std;
using namespace sflow;

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
    cutflow->setLumi(78.3); // set the MC normalized to lumi periods A1-A3
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
    *cutflow << CutName("at least 1 signal lepton") << [](Superlink* sl) -> bool {
        return sl->leptons->size() >= 1;
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

    *cutflow << NewVar("mu passes HLT_mu24_iloose_L1MU15"); {
        *cutflow << HFTname("muPass_mu24_iloose_L1MU15");
        *cutflow << [&](Superlink* sl, var_bool_array*) -> vector<bool> {
            vector<bool> out;
            for(uint imu=0; imu<muons.size(); imu++){
                if(sl->tools->triggerTool().passTrigger(muons[imu]->trigBits, "HLT_mu24_iloose_L1MU15")) { out.push_back(true); }
                else { out.push_back(false); }
            }
            return out;
            };
        *cutflow << SaveVar();
    }

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
    *cutflow << NewVar("lepton flavor (0: e, 1: m)"); {
        *cutflow << HFTname("l_flav");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> lep_flav;
            for(int i = 0; i < leptons.size(); i++) {
                lep_flav.push_back(leptons.at(i)->isEle() ? 0 : 1);
            }
            return lep_flav;
            };
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
    *cutflow << NewVar("lepton d0"); {
        *cutflow << HFTname("l_d0");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->d0);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton errD0"); {
        *cutflow << HFTname("l_errD0");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->errD0);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton d0sig"); {
        *cutflow << HFTname("l_d0sig");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> d0sig;
            for(int i = 0; i < leptons.size(); i++) {
                d0sig.push_back(leptons.at(i)->d0Sig());
            }
            return d0sig;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton d0sig (BSCorr)"); {
        *cutflow << HFTname("l_d0sigBSCorr");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> d0sigBSCorr;
            for(int i = 0; i < leptons.size(); i++) {
                d0sigBSCorr.push_back(leptons.at(i)->d0sigBSCorr);
            }
            return d0sigBSCorr;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton z0"); {
        *cutflow << HFTname("l_z0");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->z0);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton errZ0"); {
        *cutflow << HFTname("l_errZ0");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->errZ0);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton z0sinTheta"); {
        *cutflow << HFTname("l_z0sinTheta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> z0;
            for(int i = 0; i < leptons.size(); i++) {
                z0.push_back(leptons.at(i)->z0SinTheta());
            }
            return z0;
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
    *cutflow << NewVar("lepton ptvarcone20"); {
        *cutflow << HFTname("l_ptvarcone20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->ptvarcone20);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton ptvarcone30"); {
        *cutflow << HFTname("l_ptvarcone30");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->ptvarcone30);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton etconetopo20"); {
        *cutflow << HFTname("l_etconetopo20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->etconetopo20);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton etconetopo30"); {
        *cutflow << HFTname("l_etconetopo30");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->etconetopo30);
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
/*
    TriggerTools* trigtool = new TriggerTools(chain, true);
    *cutflow << NewVar("m_mumu - invariant mass of di-muon events with muons matched to HLT_mu14"); {
        *cutflow << HFTname("m_mumu_mu14");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            MuonVector matchedMuons;
            for(int i = 0; i < muons.size(); i++){
                Muon* mu = muons.at(i);
                if(sl->ntTrig->passTriggerTools(mu->trigBits, "HLT_mu14")) { matchedMuons.push_back(mu); }
            }
            double mll = -999.0;
            if(matchedMuons.size()==2) {
                Muon* mu0 = matchedMuons.at(0);
                Muon* mu1 = matchedMuons.at(1);
                mll = (*mu0 + *mu1).M();
            }
            return mll;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("m_mumu - invariant mass of di-muon events with muons matched to HLT_mu26_imedium"); {
        *cutflow << HFTname("m_mumu_mu26");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            MuonVector matchedMuons;
            for(int i = 0; i < muons.size(); i++){
                Muon* mu = muons.at(i);
                if(sl->ntTrig->passTriggerTools(mu->trigBits, "HLT_mu26_imedium")) { matchedMuons.push_back(mu); }
            }
            double mll = -999.0;
            if(matchedMuons.size()==2) {
                Muon* mu0 = matchedMuons.at(0);
                Muon* mu1 = matchedMuons.at(1);
                mll = (*mu0 + *mu1).M();
            }
            return mll;
        };
        *cutflow << SaveVar();
    }
    delete trigtool;
*/
/*
*/

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
 
    *cutflow << NewVar("dphi between lead jet and dilepton system (Z-balance)"); {
        *cutflow << HFTname("phi_balance");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double balance = -99.0;
            if(leptons.size()==2 && jets.size()>0) {
                Lepton* l0 = leptons.at(0);
                Lepton* l1 = leptons.at(1);
                TLorentzVector ll = *l0 + *l1;
                Jet* j0 = jets.at(0);

                balance = j0->DeltaPhi(ll);
            }
            return balance;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("pT balance (pT-leadJet / pTll)"); {
        *cutflow << HFTname("pT_balance");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double balance = -99.0;
            if(leptons.size()==2 && jets.size()>0) {
                Lepton* l0 = leptons.at(0);
                Lepton* l1 = leptons.at(1);
                TLorentzVector ll = *l0 + *l1;
                Jet* j0 = jets.at(0);

                balance = (j0->Pt() / ll.Pt()) * 1.0;
            }
            return balance;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("pT balance of sub-lead (pT-sublead /pTll)"); {
        *cutflow << HFTname("pT_sublead_balance");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double balance = -99.0;
            if(leptons.size()==2 && jets.size()>1) {
                Lepton* l0 = leptons.at(0);
                Lepton* l1 = leptons.at(1);
                TLorentzVector ll = *l0 + *l1;
                Jet* j1 = jets.at(1);

                balance = (j1->Pt() / ll.Pt()) * 1.0;
            }
            return balance;
        };
        *cutflow << SaveVar();
    }
    
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

    *cutflow << NewVar("jet jvf"); {
        *cutflow << HFTname("j_jvf");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->jvf);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet jvf"); {
        *cutflow << HFTname("sj_jvf");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->jvf);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet jvf"); {
        *cutflow << HFTname("bj_jvf");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->jvf);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("jet jvt"); {
        *cutflow << HFTname("j_jvt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->jvt);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet jvt"); {
        *cutflow << HFTname("sj_jvt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->jvt);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet jvt"); {
        *cutflow << HFTname("bj_jvt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->jvt);
            }
            return out;
            };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet nTracks"); {
        *cutflow << HFTname("j_nTracks");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->nTracks);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet nTracks"); {
        *cutflow << HFTname("sj_nTracks");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->nTracks);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet nTracks"); {
        *cutflow << HFTname("bj_nTracks");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->nTracks);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("jet mv2c20"); {
        *cutflow << HFTname("j_mv2c20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->mv2c20);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet mv2c20"); {
        *cutflow << HFTname("sj_mv2c20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->mv2c20);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet mv2c20"); {
        *cutflow << HFTname("bj_mv2c20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->mv2c20);
            }
            return out;
            };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet isBJet"); {
        *cutflow << HFTname("j_isBJet");
        *cutflow << [&](Superlink* sl, var_bool_array*) -> vector<bool> {
            vector<bool> out;
            for(int i = 0; i < jets.size(); i++) {
                if(jets.at(i)->bjet) { out.push_back(true); }
                else { out.push_back(false); }
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
    
    *cutflow << NewVar("x-coord. of Etmiss"); {
        *cutflow << HFTname("metX");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return (met.lv().Pt() * cos(met.lv().Phi())); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("y-coord. of Etmiss"); {
        *cutflow << HFTname("metY");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return (met.lv().Pt() * sin(met.lv().Phi())); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("met rel"); {
        *cutflow << HFTname("Etmiss rel");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return kin::getMetRel(&met, *sl->leptons, *sl->jets);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ref electron term et"); {
        *cutflow << HFTname("refEle_et");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refEle_et;
        };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("Ref gamma term et"); {
        *cutflow << HFTname("refGamma_et");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refGamma_et;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ref tau term et"); {
        *cutflow << HFTname("refTau_et");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refTau_et;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ref jet term et"); {
        *cutflow << HFTname("refJet_et");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refJet_et;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Met soft term et"); {
        *cutflow << HFTname("softTerm_et");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.softTerm_et;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ref muon term et"); {
        *cutflow << HFTname("refMuo_et");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refMuo_et;
        };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("sumet"); {
        *cutflow << HFTname("sumet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.sumet;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ref electron term sumet"); {
        *cutflow << HFTname("refEle_sumet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refEle_sumet;
        };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("Ref gamma term sumet"); {
        *cutflow << HFTname("refGamma_sumet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refGamma_sumet;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ref tau term sumet"); {
        *cutflow << HFTname("refTau_sumet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refTau_sumet;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ref jet term sumet"); {
        *cutflow << HFTname("refJet_sumet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refJet_sumet;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Met soft term sumet"); {
        *cutflow << HFTname("softTerm_sumet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.softTerm_sumet;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ref muon term sumet"); {
        *cutflow << HFTname("refMuo_sumet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refMuo_sumet;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ref electron term phi"); {
        *cutflow << HFTname("refEle_phi");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refEle_phi;
        };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("Ref gamma term phi"); {
        *cutflow << HFTname("refGamma_phi");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refGamma_phi;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ref tau term phi"); {
        *cutflow << HFTname("refTau_phi");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refTau_phi;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ref jet term phi"); {
        *cutflow << HFTname("refJet_phi");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refJet_phi;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Met soft term phi"); {
        *cutflow << HFTname("softTerm_phi");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.softTerm_phi;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ref muon term phi"); {
        *cutflow << HFTname("refMuo_phi");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.refMuo_phi;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("dphi met ref ele"); {
        *cutflow << HFTname("dphi_met_eleTerm");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector refTLV;
            float et  = met.refEle_et;
            float phi = met.refEle_phi;
            refTLV.SetPtEtaPhiE(et, 0, phi, et);
            return met.lv().DeltaPhi(refTLV);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("dphi met ref gamma"); {
        *cutflow << HFTname("dphi_met_gammaTerm");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector refTLV;
            float et  = met.refGamma_et;
            float phi = met.refGamma_phi;
            refTLV.SetPtEtaPhiE(et, 0, phi, et);
            return met.lv().DeltaPhi(refTLV);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar( "dphi met ref tau"); {
        *cutflow << HFTname("dphi_met_tauTerm");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector refTLV;
            float et  = met.refTau_et;
            float phi = met.refTau_phi;
            refTLV.SetPtEtaPhiE(et, 0, phi, et);
            return met.lv().DeltaPhi(refTLV);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar( "dphi met ref jet"); {
        *cutflow << HFTname("dphi_met_jetTerm");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector refTLV;
            float et  = met.refJet_et;
            float phi = met.refJet_phi;
            refTLV.SetPtEtaPhiE(et, 0, phi, et);
            return met.lv().DeltaPhi(refTLV);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar( "dphi met ref soft"); {
        *cutflow << HFTname("dphi_met_softTerm");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector refTLV;
            float et  = met.softTerm_et;
            float phi = met.softTerm_phi;
            refTLV.SetPtEtaPhiE(et, 0, phi, et);
            return met.lv().DeltaPhi(refTLV);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar( "dphi met ref muon"); {
        *cutflow << HFTname("dphi_met_muonTerm");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector refTLV;
            float et  = met.refMuo_et;
            float phi = met.refMuo_phi;
            refTLV.SetPtEtaPhiE(et, 0, phi, et);
            return met.lv().DeltaPhi(refTLV);
        };
        *cutflow << SaveVar();
    }

    TrackMet trackMet;
    *cutflow << [&](Superlink* sl, var_void*) { trackMet = *sl->trackMet; };
    *cutflow << NewVar("transverse missing energy, Track (Track MET)"); {
        *cutflow << HFTname("trackMet");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return trackMet.lv().Pt(); };
        *cutflow << SaveVar();
    }
    

    *cutflow << NewVar("phi coord. of Etmiss, Track (Track MET phi)"); {
        *cutflow << HFTname("trackMet_phi");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return trackMet.lv().Phi(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("sumet of Etmiss, Track (Track MET sumet)"); {
        *cutflow << HFTname("trackMet_sumet");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return trackMet.sumet; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("x-coord. of Etmiss Track"); {
        *cutflow << HFTname("trackMetX");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return (trackMet.lv().Pt() * cos(trackMet.lv().Phi())); };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("y-coord. of Etmiss Track"); {
        *cutflow << HFTname("trackMetY");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return (trackMet.lv().Pt() * sin(trackMet.lv().Phi())); };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("transverse mass"); {
        *cutflow << HFTname("mT0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double mT = -999;
            if(leptons.size()>0) {
                mT = sqrt(2 * leptons.at(0)->Pt() * met.Et * (1 - cos(leptons.at(0)->DeltaPhi((met.lv())))));
            }
            return mT;
        };
        *cutflow << SaveVar();
    }


    ////////////////////////////////////////
    // WW-like analysis variables
    ////////////////////////////////////////
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

    *cutflow << NewVar("deltaX: pz0 + pz1 / sqrt(s)"); {
        *cutflow << HFTname("deltaX");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double deltaX = -999.0;
            if(leptons.size()==2){
                double sqrtS= 13000.0;
                double num = leptons.at(0)->Pz() + leptons.at(1)->Pz();
                deltaX = num / sqrtS * 1.0;
            }
            return deltaX;
        };
        *cutflow << SaveVar();
    }

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
    *cutflow << NewVar("cosine(Theta)_ll"); {
        *cutflow << HFTname("cosThetaLL");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double cosThetaLL = -999.0;
            if(leptons.size()==2){
                TLorentzVector lp;
                TLorentzVector lm;
                for(unsigned int iL = 0; iL < leptons.size(); iL++){
                    if(leptons.at(iL)->q < 0) lm = *leptons.at(iL);
                    else if (leptons.at(iL)->q > 0) lp = *leptons.at(iL);
                }
                cosThetaLL = tanh((lp.Eta()-lm.Eta())/2.);
            }
            return cosThetaLL;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("stransverse mass (mT2)"); {
        *cutflow << HFTname("mt2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double mt2 = -999.0;
            if(leptons.size()==2) {

                mt2 = kin::getMT2(*sl->leptons, *sl->met);

                //mt2_bisect::mt2 mt2_event;
                //double *pa, *pb, *pmiss;
                //pa = new double[3]; pa[0] = sl->leptons->at(0)->M(); pa[1] = sl->leptons->at(0)->Px(), pa[2] = sl->leptons->at(0)->Py();
                //pb = new double[3]; pb[0] = sl->leptons->at(1)->M(); pb[1] = sl->leptons->at(1)->Px(), pb[2] = sl->leptons->at(1)->Py();
                //pmiss = new double[3]; pmiss[0] = 0.0; pmiss[1] = sl->met->Et * cos(sl->met->phi); pmiss[2] = sl->met->Et * sin(sl->met->phi);
 
                //mt2_event.set_momenta(pa, pb, pmiss);
                //mt2_event.set_mn(0.0); // LSP mass = 0 is Generic
 
                //mt2 = mt2_event.get_mt2();
 
                //delete[] pa;
                //delete[] pb;
                //delete[] pmiss;
            }
            return mt2;
        };
        *cutflow << SaveVar();
    }

    TLorentzVector pb_ll;
    *cutflow << NewVar("pb_ll"); {
        *cutflow << HFTname("pBll");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            double pBll_pt = -999.0;
            if(leptons.size()==2) {
                pb_ll = (met.lv() + *leptons.at(0) + *leptons.at(1));
                pBll_pt = pb_ll.Pt();
            }
            return pBll_pt;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("dphi pBll-mET"); {
        *cutflow << HFTname("dphi_met_pbll");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            TLorentzVector met_tlv = met.lv();
            double dphi = -999.0;
            if(leptons.size()==2) {
                dphi = fabs(met_tlv.DeltaPhi(pb_ll));
            }
            return dphi;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta phi between met and closest jet"); {
        *cutflow << HFTname("dphi_met_j");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            double dphi = -999.0;
            TLorentzVector met_tlv = met.lv();
            double dR = 999;
            for(unsigned int ij = 0; ij < jets.size(); ij++){
                Jet j = *jets.at(ij);
                if(fabs(met_tlv.DeltaR(j)) < dR) {
                    dR = fabs(met_tlv.DeltaR(j));
                    dphi = fabs(met_tlv.DeltaPhi(j));
                }
            }
            return dphi;
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


    *cutflow << [&](Superlink* sl, var_void*) { leptons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { electrons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { muons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { jets.clear(); }; 
    *cutflow << [&](Superlink* sl, var_void*) { bjets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { sjets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { met.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { trackMet.clear(); };
    
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
    bool verbose = true;
    ChainHelper::addInput(chain, input, verbose);
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
