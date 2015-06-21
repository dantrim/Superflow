// SuperflowAna.cxx
//

#include <cstdlib>
#include <cmath>
#include <fstream> 
#include <iostream>
#include <string>
#include <getopt.h>

#include "TChain.h"
#include "TVectorD.h"
#include "TVector3.h"

#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"
#include "SusyNtuple/MCWeighter.h"
#include "SusyNtuple/SusyNtSys.h"

#include "Superflow/Superflow.h"
#include "Superflow/Superlink.h"
#include "Superflow/Cut.h"
#include "Superflow/StringTools.h"
#include "Superflow/PhysicsTools.h"
#include "Superflow/LeptonTruthDefinitions.h"

// RestFrames
#include "RestFrames/RestFrame.hh"
#include "RestFrames/RFrame.hh"
#include "RestFrames/RLabFrame.hh"
#include "RestFrames/RDecayFrame.hh"
#include "RestFrames/RVisibleFrame.hh"
#include "RestFrames/RInvisibleFrame.hh"
#include "RestFrames/RSelfAssemblingFrame.hh"
#include "RestFrames/InvisibleMassJigsaw.hh"
#include "RestFrames/InvisibleRapidityJigsaw.hh"
#include "RestFrames/ContraBoostInvariantJigsaw.hh"
#include "RestFrames/MinimizeMassesCombinatoricJigsaw.hh"
#include "RestFrames/InvisibleGroup.hh"
#include "RestFrames/CombinatoricGroup.hh"
#include "RestFrames/FramePlot.hh"

#include "TROOT.h"
#include "TSystem.h"



#include "Mt2/mt2_bisect.h"


using namespace std;
using namespace sflow;
using namespace RestFrames;

// constants
const double GeV_to_MeV = 1000.0;
//const double GeV_to_MeV = 1.0;  // rather than comment out all of the places where we thought MeV was reasonable

// function prototypes
void print_usage(const char *exeName);
void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_, SuperflowRunMode& run_mode_, SusyNtSys& nt_sysnt_);


int main(int argc, char* argv[])
{
  //  gSystem->Load("/export/home/dantrim/RestFrames/lib/libRestFrames.so");
    // START read-in
    int n_skip_ = 0;
    int num_events_ = -1;
    string sample_;
    SuperflowRunMode run_mode = SuperflowRunMode::nominal;
    SusyNtSys nt_sys_ = NtSys::NOM;

    TChain* chain = new TChain("susyNt");
    chain->SetDirectory(0);

    read_options(argc, argv, chain, n_skip_, num_events_, sample_, run_mode, nt_sys_); // defined below
    // END read-in

    Superflow* cutflow = new Superflow(); // initialize the cutflow
    cutflow->setAnaType(AnalysisType::Ana_2Lep); // Ana_2Lep Ana_2LepWH 
    cutflow->setSampleName(sample_);
    cutflow->setRunMode(run_mode);
    cutflow->setCountWeights(true); // print the weighted cutflows
    cutflow->setChain(chain);

    string xsecDir = gSystem->ExpandPathName("$ROOTCOREBIN/data/SUSYTools/mc12_8TeV/");
    MCWeighter* weighter = new MCWeighter(chain, xsecDir);

    cout << "Analysis    Total Entries: " << chain->GetEntries() << endl;

    if (run_mode == SuperflowRunMode::single_event_syst) cutflow->setSingleEventSyst(nt_sys_);

    // START Setup cuts
    // START Setup cuts
    // START Setup cuts
    
    *cutflow << CutName("read in") << [](Superlink* sl) -> bool { return true; };

    int cutFlags = 0;
    *cutflow << CutName("pass GRL") << [&](Superlink* sl) -> bool { 
        cutFlags = sl->nt->evt()->cutFlags[sl->nt_sys];
        bool pass_grl(cutFlags & ECut_GRL);
        return pass_grl;
    };
    
    *cutflow << CutName("error flags") << [&](Superlink* sl) -> bool {
        bool passLarErr(cutFlags & ECut_LarErr);
        bool pass_tile(cutFlags & ECut_TileErr);
        bool pass_TTC(cutFlags & ECut_TTC);
        bool passFlags(passLarErr && pass_tile && pass_TTC);
        return passFlags;
    };

    *cutflow << CutName("bad muon") << [&](Superlink* sl) -> bool {
        bool pass_bad_moun(cutFlags & ECut_BadMuon);
        return pass_bad_moun;
    };

    *cutflow << CutName("jet cleaning") << [&](Superlink *sl) -> bool {
        bool passJetCleaning(cutFlags & ECut_BadJet);
        return passJetCleaning;
    };

    *cutflow << CutName("pass good vtx") << [&](Superlink *sl) -> bool {
        bool pass_goodpv(cutFlags & ECut_GoodVtx);
        return pass_goodpv;
    }; 
    
    *cutflow << CutName("cosmics") << [&](Superlink *sl) -> bool {
        bool pass_cosmic(cutFlags & ECut_Cosmic);
        return pass_cosmic;
    };

    *cutflow << CutName("exactly two base leptons") << [](Superlink* sl) -> bool {
        return sl->baseLeptons->size() == 2;
    };

 //   *cutflow << CutName("m_ll > 20 GeV") << [](Superlink* sl) -> bool {
 //       return (*sl->baseLeptons->at(0) + *sl->baseLeptons->at(1)).M() > 20.0;
 //   };

 //   *cutflow << CutName("muon eta < 2.4") << [](Superlink* sl) -> bool {
 //       return (!sl->leptons->at(0)->isMu() || abs(sl->leptons->at(0)->Eta()) < 2.4)
 //           && (!sl->leptons->at(1)->isMu() || abs(sl->leptons->at(1)->Eta()) < 2.4);
 //   };

 //   *cutflow << CutName("tau veto") << [](Superlink* sl) -> bool {
 //       return sl->taus->size() == 0;
 //   };

//    *cutflow << CutName("pass dilepton trigger") << [](Superlink* sl) -> bool {
//        return sl->dileptonTrigger->passDilTrig(*sl->leptons, sl->met->lv().Pt(), sl->nt->evt());
//    };

//    *cutflow << CutName("pass MET trigger") << [](Superlink* sl) -> bool {
//        return sl->nt->evt()->trigFlags & TRIG_xe80T_tclcw_loose;
//    };

    *cutflow << CutName("exactly 2 signal leptons") << [](Superlink* sl) -> bool {
        return (sl->leptons->size()==2);
    };


 //   *cutflow << CutName("prompt leptons") << [&](Superlink* sl) -> bool {       // DANTRIM (capture)
 //       bool pass_ = true;
 //       if (sl->isMC) {                        
 //           for (int l_ = 0; l_ < sl->leptons->size(); l_++) {
 //               bool isReal = sl->leptons->at(l_)->truthType == LeptonTruthType::PROMPT;

 //               bool isChargeFlip = sl->leptons->at(l_)->isEle()
 //                   && static_cast<Electron*>(sl->leptons->at(l_))->isChargeFlip;

 //               if (!isReal || isChargeFlip) pass_ = false;
 //           }
 //       }
 //       return pass_;
 //   };
    
//    *cutflow << CutName("same-sign") << [](Superlink *sl) -> bool {
//        return sl->leptons->at(0)->q * sl->leptons->at(1)->q > 0;
//    };


    *cutflow << CutName("exactly 2 signal leptons") << [](Superlink* sl) -> bool {
        return (sl->leptons->size()==2);
    };

    *cutflow << CutName("opposite sign") << [](Superlink* sl) -> bool {
        return (sl->leptons->at(0)->q * sl->leptons->at(1)->q < 0);
    };


    // END Setup cuts
    // END Setup cuts
    // END Setup cuts

    // GAP //
    // GAP //
    // GAP //

    // START Setup output trees
    // START Setup output trees
    // START Setup output trees




    *cutflow << NewVar("event weight"); {
        *cutflow << HFTname("eventweight");
        *cutflow << [](Superlink* sl, var_double*) -> double { 
             return sl->weights->product();
        };
        *cutflow << SaveVar();
    }

    // VARIABLES FOR DEBUGGING "NEW" AND "OLD" ALPGEN+PYTHIA Z+JETS [BEGIN]
    // VARIABLES FOR DEBUGGING "NEW" AND "OLD" ALPGEN+PYTHIA Z+JETS [BEGIN]
    // VARIABLES FOR DEBUGGING "NEW" AND "OLD" ALPGEN+PYTHIA Z+JETS [BEGIN]

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

    *cutflow << NewVar("sumw"); {
        *cutflow << HFTname("sumw");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return weighter->getSumw(sl->nt->evt());
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("monte carlo generator event weight"); {
        *cutflow << HFTname("w");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->w; 
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pile-up weight"); {
        *cutflow << HFTname("pupw");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->wPileup;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("xsec times eff"); {
        *cutflow << HFTname("xsec");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return weighter->getXsecTimesEff(sl->nt->evt(), MCWeighter::Sys_NOM);
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("run number"); {
        *cutflow << HFTname("runNumber");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->run; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("event number"); {
        *cutflow << HFTname("eventNumber");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->event; };
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
                if( sl->tools->m_jetSelector.isCentralLightJet(sl->jets->at(i))){
                //if (sl->tools->isCentralLightJet(sl->jets->at(i), sl->jvfTool, sl->nt_sys, sl->anaType)) {
                    central_light_jets.push_back(sl->jets->at(i));
                }
            }
            return central_light_jets.size();
        };
        *cutflow << SaveVar();
    }

 //   *cutflow << NewVar("number of central b jets"); {
 //       *cutflow << HFTname("nCentralBJets");
 //       *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfCBJets(*sl->jets); };
 //       *cutflow << SaveVar();
 //   }

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
            return central_light_jets.size() >= 2 ? central_light_jets[1]->Pt()  : 0.0;
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
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->tools->getMetRel(sl->met, *sl->leptons, *sl->jets); };
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

    // Super-Razor variables
    double mDeltaR, shatr, cosThetaRp1, dphi_ll_vBetaT, dphi_l1_l2, gamma_r = 0.0;
    double dphi_vBeta_R_vBeta_T = 0.0;
    TVector3 vBeta_z, pT_CM, vBeta_T_CMtoR, vBeta_r;
    *cutflow << NewVar("Super-Razor variables -- shatr"); {
        *cutflow << HFTname("shatr");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            sl->tools->superRazor(*sl->leptons, sl->met, vBeta_z, pT_CM,
                                    vBeta_T_CMtoR, vBeta_r, shatr, dphi_ll_vBetaT,
                                    dphi_l1_l2, gamma_r, dphi_vBeta_R_vBeta_T,
                                    mDeltaR, cosThetaRp1);
            return shatr;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("Super-Razor variables -- mDeltaR"); {
        *cutflow << HFTname("mDeltaR");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return mDeltaR;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("Super-Razor variables -- cos(Theta_{R+1})"); {
        *cutflow << HFTname("cosThetaRp1");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return cosThetaRp1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("Super-Razor variables -- Delta Phi_{beta}^{R}"); {
        *cutflow << HFTname("dphi_ll_vBetaT");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return dphi_ll_vBetaT;
        };
        *cutflow << SaveVar();
    } 

/*    ////////////////////////////////////////////////
    // Setup RestFrames 
    ////////////////////////////////////////////////

    // frames of interest
    RLabFrame LAB("LAB", "lab");
    RDecayFrame SS("SS", "SS");
    RSelfAssemblingFrame S1("S1", "#tilde{S}_{a}");
    RSelfAssemblingFrame S2("S2", "#tilde{S}_{b}");
    RVisibleFrame V1("V1", "V_{a}");
    RVisibleFrame V2("V2", "V_{b}");
    RInvisibleFrame I1("I1", "I_{a}");
    RInvisibleFrame I2("I2", "I_{b}");
    
    /////////////////////////
    // groups of particles
    /////////////////////////

    // The invisible group is all of the weakly interacting particles
    InvisibleGroup INV("INV", "Invisible State Jigsaws");
    INV.AddFrame(I1);
    INV.AddFrame(I2);
    
    // the combinatric group is the list of visible objects
    // that go into ourhemispheres
    CombinatoricGroup VIS("VIS", "Visible Object Jigsaws");
    VIS.AddFrame(V1);
    VIS.SetNElementsForFrame(V1,1,false);
    VIS.AddFrame(V2);
    VIS.SetNElementsForFrame(V2,1,false);

    //////////////////////////
    // Define how all the frames connect
    // in the decay tree
    //////////////////////////
    LAB.SetChildFrame(SS);
    
    SS.AddChildFrame(S1);
    SS.AddChildFrame(S2);
    
    S1.AddChildFrame(V1);
    S1.AddChildFrame(I1);
    S2.AddChildFrame(V2);
    S2.AddChildFrame(I2);

    ///////////////////////////
    // Make sure the tree is logical
    ///////////////////////////
    LAB.InitializeTree();

    ////////////////////////////////////////////
    // now define the 'jigsaw rules' that tell the
    // tree how to define the objects in our groups
    ////////////////////////////////////////////
    InvisibleMassJigsaw MinMassJigsaw("MINMASS_JIGSAW", "Invisible system mass Jigsaw");
    INV.AddJigsaw(MinMassJigsaw);
    
    InvisibleRapidityJigsaw RapidityJigsaw("RAPIDITY_JIGSAW", "Invisible system rapidity Jigsaw");
    INV.AddJigsaw(RapidityJigsaw);
    RapidityJigsaw.AddVisibleFrame((LAB.GetListVisibleFrames()));

    ContraBoostInvariantJigsaw ContraBoostJigsaw("CB_JIGSAW", "Contraboost invariant Jigsaw");
    INV.AddJigsaw(ContraBoostJigsaw);
    ContraBoostJigsaw.AddVisibleFrame((S1.GetListVisibleFrames()), 0);
    ContraBoostJigsaw.AddVisibleFrame((S2.GetListVisibleFrames()), 1);
    ContraBoostJigsaw.AddInvisibleFrame((S1.GetListInvisibleFrames()), 0);
    ContraBoostJigsaw.AddInvisibleFrame((S2.GetListInvisibleFrames()), 1);

    MinimizeMassesCombinatoricJigsaw HemiJigsaw("HEM_JIGSAW", "Minimize m_{V_{a,b}} Jigsaw");
    VIS.AddJigsaw(HemiJigsaw);
    HemiJigsaw.AddFrame(V1,0);
    HemiJigsaw.AddFrame(V2,1);
*/    
    ////////////////////////////////
    // check that jigsaws are logical
    ////////////////////////////////
    vector<GroupElementID> jetID;

    double shatr_rr;
    double inv_gammaRp1_rr;
    double inv_E1_in_S2;
    int n_jets_in_Hem1;
    int n_jets_in_Hem2;
    int jet1_in_Hem1; 
    double cos_theta_SS;
    int wimp_depth_S1;
    int wimp_depth_S2;
    double mDeltaR_rr;
    double dphi_ll_vBetaT_rr;


    *cutflow << [&](Superlink* sl, var_void*) {

        ////////////////////////////////////////////////
        // Setup RestFrames 
        ////////////////////////////////////////////////

        // frames of interest
        RLabFrame LAB("LAB", "lab");
        RDecayFrame SS("SS", "SS");
        RSelfAssemblingFrame S1("S1", "#tilde{S}_{a}");
        RSelfAssemblingFrame S2("S2", "#tilde{S}_{b}");
        RVisibleFrame V1("V1", "V_{a}");
        RVisibleFrame V2("V2", "V_{b}");
        RInvisibleFrame I1("I1", "I_{a}");
        RInvisibleFrame I2("I2", "I_{b}");
        
        /////////////////////////
        // groups of particles
        /////////////////////////

        // The invisible group is all of the weakly interacting particles
        InvisibleGroup INV("INV", "Invisible State Jigsaws");
        INV.AddFrame(I1);
        INV.AddFrame(I2);
        
        // the combinatric group is the list of visible objects
        // that go into ourhemispheres
        CombinatoricGroup VIS("VIS", "Visible Object Jigsaws");
        VIS.AddFrame(V1);
        VIS.SetNElementsForFrame(V1,1,false);
        VIS.AddFrame(V2);
        VIS.SetNElementsForFrame(V2,1,false);

        //////////////////////////
        // Define how all the frames connect
        // in the decay tree
        //////////////////////////
        LAB.SetChildFrame(SS);
        
        SS.AddChildFrame(S1);
        SS.AddChildFrame(S2);
        
        S1.AddChildFrame(V1);
        S1.AddChildFrame(I1);
        S2.AddChildFrame(V2);
        S2.AddChildFrame(I2);

        ///////////////////////////
        // Make sure the tree is logical
        ///////////////////////////
        LAB.InitializeTree();

        ////////////////////////////////////////////
        // now define the 'jigsaw rules' that tell the
        // tree how to define the objects in our groups
        ////////////////////////////////////////////
        InvisibleMassJigsaw MinMassJigsaw("MINMASS_JIGSAW", "Invisible system mass Jigsaw");
        INV.AddJigsaw(MinMassJigsaw);
        
        InvisibleRapidityJigsaw RapidityJigsaw("RAPIDITY_JIGSAW", "Invisible system rapidity Jigsaw");
        INV.AddJigsaw(RapidityJigsaw);
        RapidityJigsaw.AddVisibleFrame((LAB.GetListVisibleFrames()));

        ContraBoostInvariantJigsaw ContraBoostJigsaw("CB_JIGSAW", "Contraboost invariant Jigsaw");
        INV.AddJigsaw(ContraBoostJigsaw);
        ContraBoostJigsaw.AddVisibleFrame((S1.GetListVisibleFrames()), 0);
        ContraBoostJigsaw.AddVisibleFrame((S2.GetListVisibleFrames()), 1);
        ContraBoostJigsaw.AddInvisibleFrame((S1.GetListInvisibleFrames()), 0);
        ContraBoostJigsaw.AddInvisibleFrame((S2.GetListInvisibleFrames()), 1);

        MinimizeMassesCombinatoricJigsaw HemiJigsaw("HEM_JIGSAW", "Minimize m_{V_{a,b}} Jigsaw");
        VIS.AddJigsaw(HemiJigsaw);
        HemiJigsaw.AddFrame(V1,0);
        HemiJigsaw.AddFrame(V2,1);

        LAB.InitializeAnalysis();

        // Clear the event (reset all frames)
        LAB.ClearEvent();
     
        // add the leptons to the "VIS" group
        for(int i = 0; i < int(sl->leptons->size()) ; i++)
            jetID.push_back(VIS.AddLabFrameFourVector(*sl->leptons->at(i)));
        // Set the "INV" group momentum in the lab frame
        //   (i.e. the MET)
        TVector3 met3vector(sl->met->lv().Px(), sl->met->lv().Py(), sl->met->lv().Pz());
        INV.SetLabFrameThreeVector(met3vector);
        LAB.AnalyzeEvent();
      //  if(!LAB.AnalyzeEvent()) return;
    
        // Can get our two visible hemispheres' four vectors in
        // any references frame in our decay tree
        TLorentzVector HEM1_inLab = V1.GetFourVector(LAB);
        TLorentzVector HEM2_inLab = V2.GetFourVector(LAB);
        TLorentzVector HEM1_inS2 = V1.GetFourVector(S2);

        // can get our friendly Super-Razor variables
        // shatr
        shatr_rr = SS.GetMass();
        // 1/ gamma_R+1
        inv_gammaRp1_rr = SS.GetVisibleShape();
        // invsible energy INV1 in S2
        inv_E1_in_S2 = I1.GetEnergy(S2);

        // mDeltaR
        mDeltaR_rr = 2.0 * V1.GetEnergy(S1);

        // dphi_ll_vBetaT_rr
        dphi_ll_vBetaT_rr = SS.GetDeltaPhiBoostVisible();
        
        // where did the jets end up?
        n_jets_in_Hem1 = VIS.GetNElementsInFrame(V1);
        n_jets_in_Hem2 = VIS.GetNElementsInFrame(V2);

        // can find out where a specific jet ended up
        jet1_in_Hem1 = V1.IsSame(VIS.GetFrame(jetID[0])) ? 1 : 0;
    
        // can get decay angles of different frames
        cos_theta_SS = SS.GetCosDecayAngle();
        
        // can find out 'depth' of the WIMPS in ecah decay hemisphere
        wimp_depth_S1 = S1.GetFrameDepth(I1);
        wimp_depth_S2 = S2.GetFrameDepth(I2); 
    };

   

    *cutflow << NewVar("shatr_rr"); {
        *cutflow << HFTname("shatr_rr");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return shatr_rr;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("1/gamma_R+1"); {
        *cutflow << HFTname("inv_gammaRp1_rr");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return inv_gammaRp1_rr;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("invisible energy INV1 in S2"); {
        *cutflow << HFTname("inv_E1_in_S2");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return inv_E1_in_S2;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mDeltaR_rr"); {
        *cutflow << HFTname("mDeltaR_rr");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return mDeltaR_rr;
        };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("dphi_ll_vBetaT_rr"); {
        *cutflow << HFTname("dphi_ll_vBetaT_rr");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return dphi_ll_vBetaT_rr;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("nJets in Hem1"); {
        *cutflow << HFTname("n_jets_in_Hem1");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return n_jets_in_Hem1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("nJets in Hem2"); {
        *cutflow << HFTname("n_jets_in_Hem2");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return n_jets_in_Hem2;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("decay angle of SS"); {
        *cutflow << HFTname("CosDecayAngle_SS");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return cos_theta_SS;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("depth of wimp 1 in S1"); {
        *cutflow << HFTname("wimp_depth_S1");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return wimp_depth_S1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("depth of wimp 2 in S2"); {
        *cutflow << HFTname("wimp_depth_S2");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return wimp_depth_S2;
        };
        *cutflow << SaveVar();
    }

 
    *cutflow << [&](Superlink* sl, var_void*) { central_light_jets.clear(); };
    // END Setup output trees
    // END Setup output trees
    // END Setup output trees

    // GAP //
    // GAP //
    // GAP //

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
    bool inputIsDir  = Susy::utils::endswith(input, "/");
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
    event_syst_map["EESZUP"] = NtSys::EES_Z_UP;
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
