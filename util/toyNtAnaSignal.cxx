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
    cutflow->setAnaName("SuperflowAna");
    cutflow->setAnaType(AnalysisType::Ana_2Lep); 
    cutflow->setLumi(6.6); // set the MC normalized to lumi periods A1-A3
    cutflow->setSampleName(sample_);
    cutflow->setRunMode(run_mode);
    cutflow->setCountWeights(true); // print the weighted cutflows
    cutflow->setChain(chain);

    cout << "Analysis    Total Entries: " << chain->GetEntries() << endl;

    if (run_mode == SuperflowRunMode::single_event_syst) cutflow->setSingleEventSyst(nt_sys_);

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
        return (cutflags & ECut_GRL);
    };
    
    *cutflow << CutName("LAr error") << [&](Superlink* sl) -> bool {
        return (cutflags & ECut_LarErr);
    };
    
    *cutflow << CutName("Tile error") << [&](Superlink* sl) -> bool {
        return (cutflags & ECut_TileErr);
    };
    
    *cutflow << CutName("TTC veto") << [&](Superlink* sl) -> bool {
        return (cutflags & ECut_TTC);
    };

    *cutflow << CutName("bad muon veto") << [&](Superlink* sl) -> bool {
        return (cutflags & ECut_BadMuon);
    };
    
    *cutflow << CutName("jet cleaning") << [&](Superlink* sl) -> bool {
        return (cutflags & ECut_BadJet);
    };
    
    *cutflow << CutName("pass good vertex") << [&](Superlink* sl) -> bool {
        return (cutflags & ECut_GoodVtx);
    };
    
    *cutflow << CutName("pass cosmic veto") << [&](Superlink* sl) -> bool {
        return (cutflags & ECut_Cosmic);
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
    
    // base electrons
    *cutflow << NewVar("ele is signal (SUSYTools 'signal' flag)"); {
        *cutflow << HFTname("el_STsignal");
        *cutflow << [&](Superlink* sl, var_bool_array*) -> vector<bool> {
            vector<bool> out;
            for(int i = 0; i < electrons.size(); i++) {
                out.push_back(electrons.at(i)->isSignal);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ele is looseLH"); {
        *cutflow << HFTname("el_isLooseLH");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> isLooseLH;
            for(int i = 0; i < electrons.size(); i++) {
                isLooseLH.push_back(electrons.at(i)->looseLH);
            }
            return isLooseLH;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ele is mediumLH"); {
        *cutflow << HFTname("el_isMediumLH");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> isMediumLH;
            for(int i = 0; i < electrons.size(); i++) {
                    isMediumLH.push_back(electrons.at(i)->mediumLH);
            }
            return isMediumLH;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ele is tightLH"); {
        *cutflow << HFTname("el_isTightLH");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> istightLH;
            for(int i = 0; i < electrons.size(); i++) {
                istightLH.push_back(electrons.at(i)->tightLH);
            }
            return istightLH;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ele is looseLH (no d0 cut)"); {
        *cutflow << HFTname("el_isLooseLH_noD0");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> istightLH;
            for(int i = 0; i < electrons.size(); i++) {
                istightLH.push_back(electrons.at(i)->looseLH_nod0);
            }
            return istightLH;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ele is mediumLH (no d0 cut)"); {
        *cutflow << HFTname("el_isMediumLH_noD0");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> istightLH;
            for(int i = 0; i < electrons.size(); i++) {
                istightLH.push_back(electrons.at(i)->mediumLH_nod0);
            }
            return istightLH;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ele is tightLH (no d0 cut)"); {
        *cutflow << HFTname("el_isTightLH_noD0");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> istightLH;
            for(int i = 0; i < electrons.size(); i++) {
                istightLH.push_back(electrons.at(i)->tightLH_nod0);
            }
            return istightLH;
            };
        *cutflow << SaveVar();
    }
/*
    *cutflow << NewVar("ele passes isoGradientLoose"); {
        *cutflow << HFTname("el_isoGradientLoose");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> isGL;
            for(int i = 0; i < electrons.size(); i++) {
                isGL.push_back(electrons.at(i)->isoGradientLoose);
            }
            return isGL;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ele passes isoGradient"); {
        *cutflow << HFTname("el_isoGradient");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> isGL;
            for(int i = 0; i < electrons.size(); i++) {
                isGL.push_back(electrons.at(i)->isoGradient);
            }
            return isGL;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ele passes isoVeryLoose"); {
        *cutflow << HFTname("el_isoVeryLoose");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> isGL;
            for(int i = 0; i < electrons.size(); i++) {
                isGL.push_back(electrons.at(i)->isoVeryLoose);
            }
            return isGL;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ele passes isoLoose"); {
        *cutflow << HFTname("el_isoLoose");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> isGL;
            for(int i = 0; i < electrons.size(); i++) {
                isGL.push_back(electrons.at(i)->isoLoose);
            }
            return isGL;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ele passes isoTight"); {
        *cutflow << HFTname("el_isoTight");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> isGL;
            for(int i = 0; i < electrons.size(); i++) {
                isGL.push_back(electrons.at(i)->isoTight);
            }
            return isGL;
            };
        *cutflow << SaveVar();
    }
*/
    // muons
    *cutflow << NewVar("mu is signal (SUSYTools 'signal' flag)"); {
        *cutflow << HFTname("mu_STsignal");
        *cutflow << [&](Superlink* sl, var_bool_array*) -> vector<bool> {
            vector<bool> out;
            for(int i = 0; i < muons.size(); i++) {
                out.push_back(muons.at(i)->isSignal);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("mu passes veryLoose"); {
        *cutflow << HFTname("mu_isVeryLoose");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < muons.size(); i++) {
                out.push_back(muons.at(i)->veryLoose);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("mu passes loose quality"); {
        *cutflow << HFTname("mu_isLoose");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < muons.size(); i++) {
                out.push_back(muons.at(i)->loose);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("mu passes medium quality"); {
        *cutflow << HFTname("mu_isMedium");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < muons.size(); i++) {
                out.push_back(muons.at(i)->medium);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("mu passes tight quality"); {
        *cutflow << HFTname("mu_isTight");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < muons.size(); i++) {
                out.push_back(muons.at(i)->tight);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
/*
    *cutflow << NewVar("mu passes isoGradientLoose"); {
        *cutflow << HFTname("mu_isoGradientLoose");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < muons.size(); i++) {
                out.push_back(muons.at(i)->isoGradientLoose);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("mu passes isoGradient"); {
        *cutflow << HFTname("mu_isoGradient");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < muons.size(); i++) {
                out.push_back(muons.at(i)->isoGradient);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("mu passes isoVeryLoose"); {
        *cutflow << HFTname("mu_isoVeryLoose");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < muons.size(); i++) {
                out.push_back(muons.at(i)->isoVeryLoose);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("mu passes isoLoose"); {
        *cutflow << HFTname("mu_isoLoose");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < muons.size(); i++) {
                out.push_back(muons.at(i)->isoLoose);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("mu passes isoTight"); {
        *cutflow << HFTname("mu_isoTight");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < muons.size(); i++) {
                out.push_back(muons.at(i)->isoTight);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
*/
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
    *cutflow << [&](Superlink* sl, var_void*) {
        JetVector susyJets = *sl->baseJets;
        for(int i = 0; i < susyJets.size(); i++) {
            Jet* j = susyJets.at(i);
            bool passPt = j->Pt()>25 ? true : false;
            bool passEta = j->Eta()<2.8 ? true : false;
            if(passPt && passEta) { jets.push_back(j); }
        }
    };

    *cutflow << [&](Superlink* sl, var_void*) {
        JetVector susyJets = *sl->baseJets;
        for(int i = 0; i < susyJets.size(); i++) {
            Jet* j = susyJets.at(i);
            bool passPt = j->Pt()> 25 ? true: false;
            bool passEta = j->Eta() < 2.8 ? true : false;
            bool isBjet = j->bjet ? true : false;
            if(passPt && passEta && isBjet) { bjets.push_back(j); }
            if(passPt && passEta && !isBjet) { sjets.push_back(j); }
        }
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
            cout << "Analysis    Error (fatal): Bad arguments." << endl;
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
        cout << "Analysis    run mode: SuperflowRunMode::nominal" << endl;
    }
    if (nominal_and_weight_syst_) {
        run_mode_ = SuperflowRunMode::nominal_and_weight_syst;
        cout << "Analysis    run mode: SuperflowRunMode::nominal_and_weight_syst" << endl;
    }
    if (single_event_syst_) {
        run_mode_ = SuperflowRunMode::single_event_syst;
        cout << "Analysis    run mode: SuperflowRunMode::single_event_syst" << endl;
    }

    if (all_syst_) {
        run_mode_ = SuperflowRunMode::all_syst;
        cout << "Analysis    run mode: SuperflowRunMode::all_syst" << endl;
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
            cout << "Analysis" << "    ERROR (fatal): Event systematic option /s " << systematic_ << " -> not found." << endl;
            exit(1);
        }
    }
}
