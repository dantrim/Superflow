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
#include "SusyNtuple/KinematicTools.h"

// Superflow
#include "Superflow/Superflow.h"
#include "Superflow/Superlink.h"
#include "Superflow/Cut.h"
#include "Superflow/StringTools.h"
#include "Superflow/PhysicsTools.h"
#include "Superflow/LeptonTruthDefinitions.h"

// Mt2
#include "Mt2/mt2_bisect.h"


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
    cutflow->setLumi(LUMI_A_A4); // set the MC normalized to lumi periods A1-A3
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
    
    *cutflow << CutName("exactly two base leptons") << [](Superlink* sl) -> bool {
        return sl->baseLeptons->size() == 2;
    };

    *cutflow << CutName("m_ll > 20 GeV") << [](Superlink* sl) -> bool {
        return (*sl->baseLeptons->at(0) + *sl->baseLeptons->at(1)).M() > 20.0;
    };

    *cutflow << CutName("tau veto") << [](Superlink* sl) -> bool {
        return sl->taus->size() == 0;
    };

    *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
        return (sl->leptons->size() == 2);
    };

    *cutflow << CutName("opposite sign") << [](Superlink* sl) -> bool {
        return (sl->leptons->at(0)->q * sl->leptons->at(1)->q < 0);
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

    // LEPTONS
    // LEPTONS
    // LEPTONS

    *cutflow << NewVar("is genuine same-sign"); {
        *cutflow << HFTname("isGenSS");
        *cutflow << [](Superlink *sl, var_bool*) -> bool { return PhysicsTools::isGenuineSS(sl->leptons); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is e + e"); {
        *cutflow << HFTname("isElEl");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is e + mu"); {
        *cutflow << HFTname("isElMu");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->leptons->at(0)->isEle() ^ sl->leptons->at(1)->isEle(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is mu + mu"); {
        *cutflow << HFTname("isMuMu");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isMu(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is opposite-sign"); {
        *cutflow << HFTname("isOS");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->leptons->at(0)->q * sl->leptons->at(1)->q < 0; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is Mu (lead) + E (sub)"); {
        *cutflow << HFTname("isME");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isEle(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is E (lead) + Mu (sub)"); {
        *cutflow << HFTname("isEM");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isMu(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 Pt"); {
        *cutflow << HFTname("lept1Pt");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(0)->Pt(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 Eta"); {
        *cutflow << HFTname("lept1Eta");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(0)->Eta(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 Phi"); {
        *cutflow << HFTname("lept1Phi");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(0)->Phi(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 Energy"); {
        *cutflow << HFTname("lept1E");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(0)->E(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 charge"); {
        *cutflow << HFTname("lept1q");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->leptons->at(0)->q; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 flavor"); { // 0=el, 1=mu, see HistFitterTree.h
        *cutflow << HFTname("lept1Flav");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->leptons->at(0)->isEle() ? 0 : 1; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 Pt"); {
        *cutflow << HFTname("lept2Pt");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(1)->Pt(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 Eta"); {
        *cutflow << HFTname("lept2Eta");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(1)->Eta(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 Phi"); {
        *cutflow << HFTname("lept2Phi");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(1)->Phi(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 Energy"); {
        *cutflow << HFTname("lept2E");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(1)->E(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 charge"); {
        *cutflow << HFTname("lept2q");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->leptons->at(1)->q; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 flavor"); {
        *cutflow << HFTname("lept2Flav");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->leptons->at(1)->isEle() ? 0 : 1; };
        *cutflow << SaveVar();
    }

    // JETS
    // JETS
    // JETS

    JetVector central_light_jets; // local variable!

    *cutflow << NewVar("number of central light jets"); {
        *cutflow << HFTname("nCentralLightJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            for (int i = 0; i < sl->jets->size(); i++) {
                if( sl->tools->jetSelector().isCentralLightJet(sl->jets->at(i))){
                    central_light_jets.push_back(sl->jets->at(i));
                }
            }
            return central_light_jets.size();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of central b jets"); {
        *cutflow << HFTname("nCentralBJets");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfCBJets(*sl->jets); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of forward jets"); {
        *cutflow << HFTname("nForwardJets");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfFJets(*sl->jets); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet-1 Pt"); {
        *cutflow << HFTname("jet1Pt");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return central_light_jets.size() >= 1 ? central_light_jets[0]->Pt() : 0.0;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet-1 Eta"); {
        *cutflow << HFTname("jet1Eta");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return central_light_jets.size() >= 1 ? central_light_jets[0]->Eta() : 0.0;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet-1 Phi"); {
        *cutflow << HFTname("jet1Phi");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return central_light_jets.size() >= 1 ? central_light_jets[0]->Phi() : 0.0;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet-2 Pt"); {
        *cutflow << HFTname("jet2Pt");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return central_light_jets.size() >= 2 ? central_light_jets[1]->Pt() : 0.0;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet-2 Eta"); {
        *cutflow << HFTname("jet2Eta");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return central_light_jets.size() >= 2 ? central_light_jets[1]->Eta() : 0.0;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet-2 Phi"); {
        *cutflow << HFTname("jet2Phi");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return central_light_jets.size() >= 2 ? central_light_jets[1]->Phi() : 0.0;
        };
        *cutflow << SaveVar();
    }

    *cutflow << [&](Superlink* sl, var_void*) { central_light_jets.clear(); };

    // MET
    // MET
    // MET

    *cutflow << NewVar("transverse missing energy (Et)"); {
        *cutflow << HFTname("met");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->met->Et; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("transverse missing energy (Phi)"); {
        *cutflow << HFTname("metPhi");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->met->phi; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Etmiss Rel"); {
        *cutflow << HFTname("metrel");
        *cutflow << [](Superlink* sl, var_float*) -> double { return kin::getMetRel(sl->met, *sl->leptons, *sl->jets); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta Phi of leading lepton and met"); {
        *cutflow << HFTname("deltaPhi_met_l1");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->leptons->at(0)->Phi() - sl->met->phi); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta Phi of subleading lepton and met"); {
        *cutflow << HFTname("deltaPhi_met_l2");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->leptons->at(1)->Phi() - sl->met->phi); };
        *cutflow << SaveVar();
    }

    // VARS
    // VARS
    // VARS

    TLorentzVector local_ll; // local variable!

    *cutflow << NewVar("mass of di-lepton system, M_ll"); {
        *cutflow << HFTname("mll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            local_ll = (*sl->leptons->at(0) + *sl->leptons->at(1));
            return local_ll.M();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Pt of di-lepton system, Pt_ll"); {
        *cutflow << HFTname("pTll");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return local_ll.Pt(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta Eta of di-lepton system"); {
        *cutflow << HFTname("deta_ll");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->leptons->at(0)->Eta() - sl->leptons->at(1)->Eta()); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta Phi of di-lepton system"); {
        *cutflow << HFTname("dphi_ll");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->leptons->at(0)->Phi() - sl->leptons->at(1)->Phi()); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ht (m_Eff: lep + met + jet)"); {
        *cutflow << HFTname("ht");
        *cutflow << [](Superlink* sl, var_float*) -> double {
            double ht = 0.0;

            ht += sl->leptons->at(0)->Pt() + sl->leptons->at(1)->Pt();
            ht += sl->met->Et;
            for (int i = 0; i < sl->jets->size(); i++) {
                if (sl->jets->at(i)->Pt() > 20.0) {
                    ht += sl->jets->at(i)->Pt();
                }
            }
            return ht;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("stransverse mass"); {
        *cutflow << HFTname("MT2");
        *cutflow << [](Superlink* sl, var_float*) -> double {

            mt2_bisect::mt2 mt2_event;
            double *pa, *pb, *pmiss;
            pa = new double[3]; pa[0] = sl->leptons->at(0)->M(); pa[1] = sl->leptons->at(0)->Px(), pa[2] = sl->leptons->at(0)->Py();
            pb = new double[3]; pb[0] = sl->leptons->at(1)->M(); pb[1] = sl->leptons->at(1)->Px(), pb[2] = sl->leptons->at(1)->Py();
            pmiss = new double[3]; pmiss[0] = 0.0; pmiss[1] = sl->met->Et * cos(sl->met->phi); pmiss[2] = sl->met->Et * sin(sl->met->phi);

            mt2_event.set_momenta(pa, pb, pmiss);
            mt2_event.set_mn(0.0); // LSP mass = 0 is Generic

            double mt2_ = mt2_event.get_mt2();
            // SUPRESSED messages "Deltasq_high not found at event 0" (~1 per 10000 events)

            delete[] pa;
            delete[] pb;
            delete[] pmiss;

            return mt2_;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Super-Razor Var: shatr"); {
        *cutflow << HFTname("shatr");
        *cutflow << [](Superlink* sl, var_float*) -> double {
            double   mDeltaR              = 0.0;
            double   SHATR                = 0.0;
            double   cosThetaRp1          = 0.0; 
            double   dphi_LL_vBETA_T      = 0.0;
            double   dphi_L1_L2           = 0.0;
            double   gamma_R              = 0.0;
            double   dphi_vBETA_R_vBETA_T = 0.0;
            TVector3 vBETA_z;
            TVector3 pT_CM;
            TVector3 vBETA_T_CMtoR;
            TVector3 vBETA_R;
            kin::superRazor(*sl->leptons, sl->met, vBETA_z, pT_CM,
                                    vBETA_T_CMtoR, vBETA_R, SHATR, dphi_LL_vBETA_T,
                                    dphi_L1_L2, gamma_R, dphi_vBETA_R_vBETA_T,
                                    mDeltaR, cosThetaRp1);
            return SHATR;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Super-Razor Var: mDeltaR"); {
        *cutflow << HFTname("mDeltaR");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            double   mDeltaR              = 0.0;
            double   SHATR                = 0.0;
            double   cosThetaRp1          = 0.0; 
            double   dphi_LL_vBETA_T      = 0.0;
            double   dphi_L1_L2           = 0.0;
            double   gamma_R              = 0.0;
            double   dphi_vBETA_R_vBETA_T = 0.0;
            TVector3 vBETA_z;
            TVector3 pT_CM;
            TVector3 vBETA_T_CMtoR;
            TVector3 vBETA_R;
            kin::superRazor(*sl->leptons, sl->met, vBETA_z, pT_CM,
                                    vBETA_T_CMtoR, vBETA_R, SHATR, dphi_LL_vBETA_T,
                                    dphi_L1_L2, gamma_R, dphi_vBETA_R_vBETA_T,
                                    mDeltaR, cosThetaRp1);
            return mDeltaR;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Super-Razor Var: cos(Theta_{R+1})"); {
        *cutflow << HFTname("cosThetaRp1");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            double   mDeltaR              = 0.0;
            double   SHATR                = 0.0;
            double   cosThetaRp1          = 0.0; 
            double   dphi_LL_vBETA_T      = 0.0;
            double   dphi_L1_L2           = 0.0;
            double   gamma_R              = 0.0;
            double   dphi_vBETA_R_vBETA_T = 0.0;
            TVector3 vBETA_z;
            TVector3 pT_CM;
            TVector3 vBETA_T_CMtoR;
            TVector3 vBETA_R;
            kin::superRazor(*sl->leptons, sl->met, vBETA_z, pT_CM,
                                    vBETA_T_CMtoR, vBETA_R, SHATR, dphi_LL_vBETA_T,
                                    dphi_L1_L2, gamma_R, dphi_vBETA_R_vBETA_T,
                                    mDeltaR, cosThetaRp1);
            return cosThetaRp1;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Super-Razor Var: Delta Phi_{beta}^{R}"); {
        *cutflow << HFTname("dphi_ll_vBetaT");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            double   mDeltaR              = 0.0;
            double   SHATR                = 0.0;
            double   cosThetaRp1          = 0.0; 
            double   dphi_LL_vBETA_T      = 0.0;
            double   dphi_L1_L2           = 0.0;
            double   gamma_R              = 0.0;
            double   dphi_vBETA_R_vBETA_T = 0.0;
            TVector3 vBETA_z;
            TVector3 pT_CM;
            TVector3 vBETA_T_CMtoR;
            TVector3 vBETA_R;
            kin::superRazor(*sl->leptons, sl->met, vBETA_z, pT_CM,
                                    vBETA_T_CMtoR, vBETA_R, SHATR, dphi_LL_vBETA_T,
                                    dphi_L1_L2, gamma_R, dphi_vBETA_R_vBETA_T,
                                    mDeltaR, cosThetaRp1);
            return dphi_LL_vBETA_T;
        };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("R1"); {
        *cutflow << HFTname("R1");
        *cutflow << [](Superlink* sl, var_float*) -> double {
            double mEff = 0.0;

            mEff += sl->leptons->at(0)->Pt() + sl->leptons->at(1)->Pt();
            mEff += sl->met->Et;
            for (int i = 0; i < sl->jets->size(); i++) {
                if (sl->jets->at(i)->Pt() > 20.0) {
                    mEff += sl->jets->at(i)->Pt();
                }
            }
            return sl->met->Et/mEff;
        };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("R2"); {
        *cutflow << HFTname("R2");
        *cutflow << [](Superlink* sl, var_float*) -> double {
            return sl->met->lv().Pt()/(sl->met->lv().Pt() + sl->leptons->at(0)->Pt() + sl->leptons->at(1)->Pt());
        };
        *cutflow << SaveVar();
    } 

    *cutflow << NewVar("DeltaX, ~parton momentum fraction difference"); {
        *cutflow << HFTname("deltaX");
        *cutflow << [](Superlink* sl, var_float*) -> double {
            return abs(2*(sl->leptons->at(0)->Pz() + sl->leptons->at(1)->Pz())/(8000));
        };
        *cutflow << SaveVar();
    }

    ////////////////////////////////////////////////////////////
    //  Output Ntuple Setup
    //      > Setup the output systematic ntuples
    ////////////////////////////////////////////////////////////
/*
    //
    // Weight variation systematics
    //  > stored in nominal output tree
    //

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

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    //
    //  Superflow methods [END]
    //
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
*/
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

    bool inputIsFile = Susy::utils::endswith(input, ".root");
    bool inputIsList = Susy::utils::endswith(input, ".txt");
    bool inputIsDir = Susy::utils::endswith(input, "/");
    bool validInput(inputIsFile || inputIsList || inputIsDir);
    if (!validInput) {
        cout << "Analysis    invalid input '" << input << "'" << endl;
        exit(1);
    }
    if (inputIsFile) {
        ChainHelper::addFile(chain, input);
        cout << "Analysis    file: " << input << endl;
        cout << "Analysis    file: " << input << endl;
        cout << "Analysis    file: " << input << endl;
        sample_ = input;
    }
    if (inputIsList) {
        ChainHelper::addFileList(chain, input);
        cout << "Analysis    list: " << input << endl;
        cout << "Analysis    list: " << input << endl;
        cout << "Analysis    list: " << input << endl;
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
        cout << "Analysis    dir: " << input << endl;
        cout << "Analysis    dir: " << input << endl;
        cout << "Analysis    dir: " << input << endl;
        sample_ = input;
    }
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
