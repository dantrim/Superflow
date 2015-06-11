// SAME.cxx
//

#include <cstdlib>
#include <cmath>
#include <fstream> 
#include <iostream>
#include <iomanip>
#include <string>

#include "TChain.h"
#include "TVectorD.h"
//#include "Cintex/Cintex.h"

#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"

#include "Superflow/Superflow.h"
#include "Superflow/Superlink.h"
#include "Superflow/Cut.h"
#include "Superflow/StringTools.h"
#include "Superflow/PhysicsTools.h"
#include "Superflow/LeptonTruthDefinitions.h"

#include "Mt2/mt2_bisect.h"

using namespace std;
using namespace sflow;

// constants
const double GeV_to_MeV = 1000.0;
bool use2dParam_ = false;  // option "/t" toggles whether to set fake 2d parametrization (pt, eta) vs. pt-only 

// function prototypes
void print_usage(const char *exeName);
void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_, SuperflowRunMode& run_mode_, SusyNtSys& nt_sysnt_, string& fake_region_);

int main(int argc, char* argv[])
{
    // START read-in
    int n_skip_ = 0;
    int num_events_ = -1;
    string sample_;
    string fake_region_;
    SuperflowRunMode run_mode = SuperflowRunMode::nominal;
    SusyNtSys nt_sys_ = NtSys_NOM;

    TChain* chain = new TChain("susyNt");
    chain->SetDirectory(0);

    read_options(argc, argv, chain, n_skip_, num_events_, sample_, run_mode, nt_sys_, fake_region_); // defined below
    // END read-in

    Superflow* cutflow = new Superflow(); // initialize the cutflow
    cutflow->setAnaType(Ana_2Lep); // Ana_2Lep Ana_2LepWH 
    cutflow->setSampleName(sample_);
    cutflow->setRunMode(run_mode);
    cutflow->setChain(chain);
    cutflow->setCountWeights(true);
    if (run_mode == SuperflowRunMode::fakes) {
        cutflow->setFakeRegion(fake_region_);
    }
    cutflow->setFake2dParam(use2dParam_);

    cout << "Analysis    Total Entries: " << chain->GetEntries() << endl;

    if (run_mode == SuperflowRunMode::single_event_syst) cutflow->setSingleEventSyst(nt_sys_);

    // START Setup cuts
    // START Setup cuts
    // START Setup cuts

    *cutflow << CutName("read in") << [](Superlink* sl) -> bool { return true; };

    *cutflow << CutName("HFOR") << [](Superlink *sl) -> bool {
        bool pass_ = true;
        if(sl->nt->evt()->hfor==4){
            pass_ = false;
        }
        return pass_;
    };

    *cutflow << CutName("Mll Cut") << [](Superlink *sl) -> bool {
        bool pass_ = true;

        if(sl->nt->evt()->isMC){
        if((sl->nt->evt()->mcChannel>=147770 && sl->nt->evt()->mcChannel<=147772) ||
           (sl->nt->evt()->mcChannel>=128975 && sl->nt->evt()->mcChannel<=128977) ||
           (sl->nt->evt()->mcChannel>=146820 && sl->nt->evt()->mcChannel<=146822)
          ){
           if(sl->nt->evt()->mllMcTruth > 60){ pass_ = false; }
           }
        if((sl->nt->evt()->mcChannel>=110805 && sl->nt->evt()->mcChannel<=110828) ||
           (sl->nt->evt()->mcChannel>=117650 && sl->nt->evt()->mcChannel<=117675)
          ){
           if(sl->nt->evt()->mllMcTruth < 60){ pass_ = false; }
           }
        } 
        
        return pass_;
    };
        


    *cutflow << CutName("at least two signal leptons") << [](Superlink* sl) -> bool {
        return !(sl->leptons->size() < 2);
    };

    *cutflow << CutName("remove higgsino events") << [](Superlink* sl) -> bool {
        return !sl->nt->evt()->eventWithSusyProp;
    };

    int cutFlags = 0;

    *cutflow << CutName("GRL, tile trip, and LAr error") << [&](Superlink* sl) -> bool {
        cutFlags = sl->tools->cleaningCutFlags(sl->nt->evt()->cutFlags[sl->nt_sys], *sl->preMuons, *sl->baseMuons, *sl->preJets, *sl->baseJets);
        return sl->tools->passGRL(cutFlags)
            && sl->tools->passTileTripCut(cutFlags)
            && sl->tools->passLarErr(cutFlags);
    };

    *cutflow << CutName("bad jets") << [&](Superlink* sl) -> bool {
        JetVector jets = sl->tools->getPreJets(sl->nt, sl->nt_sys);
        sl->tools->e_j_overlap(*sl->baseElectrons, jets, J_E_DR, true);
        sl->tools->t_j_overlap(*sl->taus, jets, J_T_DR, true);
        return !sl->tools->hasBadJet(jets);
    };

    *cutflow << CutName("dead regions") << [&](Superlink* sl) -> bool {
        return sl->tools->passDeadRegions(*sl->preJets, sl->met, sl->nt->evt()->run, sl->nt->evt()->isMC);
    };

    *cutflow << CutName("bad muons") << [&](Superlink* sl) -> bool {
        return !sl->tools->hasBadMuon(*sl->preMuons);
    };

    *cutflow << CutName("cosmic muons") << [&](Superlink* sl) -> bool {
        return !sl->tools->hasCosmicMuon(*sl->baseMuons);
    };

    *cutflow << CutName("hotspot jets") << [&](Superlink* sl) -> bool {
        return !sl->tools->hasHotSpotJet(*sl->preJets);
    };

    *cutflow << CutName("TTC Veto and Good Vertex") << [&](Superlink* sl) -> bool {
        return sl->tools->passTTCVeto(cutFlags) && sl->tools->passGoodVtx(cutFlags);
    };

    *cutflow << CutName("exactly two base leptons") << [](Superlink* sl) -> bool {
        return sl->baseLeptons->size() == 2;
    };

    *cutflow << CutName("m_ll > 20 GeV") << [](Superlink* sl) -> bool {
        return (*sl->baseLeptons->at(0) + *sl->baseLeptons->at(1)).M() > 20.0;
    };

    *cutflow << CutName("muon eta < 2.4") << [](Superlink* sl) -> bool {
        return (!sl->leptons->at(0)->isMu() || abs(sl->leptons->at(0)->Eta()) < 2.4)
            && (!sl->leptons->at(1)->isMu() || abs(sl->leptons->at(1)->Eta()) < 2.4);
    };

    *cutflow << CutName("tau veto") << [](Superlink* sl) -> bool {
        return sl->taus->size() == 0;
    };

    *cutflow << CutName("pass dilepton trigger") << [](Superlink* sl) -> bool {
        return sl->dileptonTrigger->passDilTrig(*sl->leptons, sl->met->lv().Pt(), sl->nt->evt());
    };

    *cutflow << CutName("prompt leptons") << [&](Superlink* sl) -> bool {       // DANTRIM (capture)
        bool pass_ = true;

        if (sl->isMC) {                      
            for (int l_ = 0; l_ < sl->leptons->size(); l_++) {
                bool isReal = sl->leptons->at(l_)->truthType == LeptonTruthType::PROMPT;

                bool isChargeFlip = sl->leptons->at(l_)->isEle()
                    && static_cast<Electron*>(sl->leptons->at(l_))->isChargeFlip;

                if (!isReal || isChargeFlip) pass_ = false;
            }
        }
        return pass_;
    };

    *cutflow << CutName("same sign") << [](Superlink* sl) -> bool {
        return sl->leptons->at(0)->q * sl->leptons->at(1)->q > 0;
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
        *cutflow << [](Superlink* sl, var_double*) -> double { please return sl->weights->product(); };
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
    
    *cutflow << NewVar("is tight-loose"); {
        *cutflow << HFTname("isTL");
        *cutflow << [](Superlink *sl, var_bool*) -> bool {
                bool isMC               = sl->nt->evt()->isMC;
                unsigned int nVtx       = sl->nt->evt()->nVtx;
                const Lepton &l0        = *sl->leptons->at(0);
                const Lepton &l1        = *sl->leptons->at(1);
                bool l0IsSig(sl->tools->isSignalLepton(&l0, *sl->baseElectrons, *sl->baseMuons, nVtx, isMC));
                bool l1IsSig(sl->tools->isSignalLepton(&l1, *sl->baseElectrons, *sl->baseMuons, nVtx, isMC));
                
                bool isTL(l0IsSig && !l1IsSig);
                return isTL;
            };
        *cutflow << SaveVar();
        }

    *cutflow << NewVar("is loose-tight"); {
        *cutflow << HFTname("isLT");
        *cutflow << [](Superlink *sl, var_bool*) -> bool {
                bool isMC               = sl->nt->evt()->isMC;
                unsigned int nVtx       = sl->nt->evt()->nVtx;
                const Lepton &l0        = *sl->leptons->at(0);
                const Lepton &l1        = *sl->leptons->at(1);
                bool l0IsSig(sl->tools->isSignalLepton(&l0, *sl->baseElectrons, *sl->baseMuons, nVtx, isMC));
                bool l1IsSig(sl->tools->isSignalLepton(&l1, *sl->baseElectrons, *sl->baseMuons, nVtx, isMC));
                
                bool isLT(!l0IsSig && l1IsSig);
                return isLT;
            };
        *cutflow << SaveVar();
        }
                
    *cutflow << NewVar("is tight-tight"); {
        *cutflow << HFTname("isTT");
        *cutflow << [](Superlink *sl, var_bool*) -> bool {
                bool isMC               = sl->nt->evt()->isMC;
                unsigned int nVtx       = sl->nt->evt()->nVtx;
                const Lepton &l0        = *sl->leptons->at(0);
                const Lepton &l1        = *sl->leptons->at(1);
                bool l0IsSig(sl->tools->isSignalLepton(&l0, *sl->baseElectrons, *sl->baseMuons, nVtx, isMC));
                bool l1IsSig(sl->tools->isSignalLepton(&l1, *sl->baseElectrons, *sl->baseMuons, nVtx, isMC));
                
                bool isTT(l0IsSig && l1IsSig);
                return isTT;
            };
        *cutflow << SaveVar();
        }
                
    *cutflow << NewVar("is loose-loose"); {
        *cutflow << HFTname("isLL");
        *cutflow << [](Superlink *sl, var_bool*) -> bool {
                bool isMC               = sl->nt->evt()->isMC;
                unsigned int nVtx       = sl->nt->evt()->nVtx;
                const Lepton &l0        = *sl->leptons->at(0);
                const Lepton &l1        = *sl->leptons->at(1);
                bool l0IsSig(sl->tools->isSignalLepton(&l0, *sl->baseElectrons, *sl->baseMuons, nVtx, isMC));
                bool l1IsSig(sl->tools->isSignalLepton(&l1, *sl->baseElectrons, *sl->baseMuons, nVtx, isMC));
                
                bool isLL(!l0IsSig && !l1IsSig);
                return isLL;
            };
        *cutflow << SaveVar();
        }
                
            

    // LEPTONS
    // LEPTONS
    // LEPTONS

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

    *cutflow << NewVar("is genuine same-same"); {
        *cutflow << HFTname("isGenSS");
        *cutflow << [](Superlink *sl, var_bool*) -> bool { return PhysicsTools::isGenuineSS(sl->leptons); };
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
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(0)->Pt() * GeV_to_MeV; };
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
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(0)->E() * GeV_to_MeV; };
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
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(1)->Pt() * GeV_to_MeV; };
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
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(1)->E() * GeV_to_MeV; };
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
                if (sl->tools->isCentralLightJet(sl->jets->at(i), sl->jvfTool, sl->nt_sys, sl->anaType)) {
                    central_light_jets.push_back(sl->jets->at(i));
                }
            }
            please return central_light_jets.size();
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
            return central_light_jets.size() >= 1 ? central_light_jets[0]->Pt() * GeV_to_MeV : 0.0;
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
            return central_light_jets.size() >= 2 ? central_light_jets[1]->Pt() * GeV_to_MeV  : 0.0;
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
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->met->Et * GeV_to_MeV ; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("transverse missing energy (Phi)"); {
        *cutflow << HFTname("metPhi");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->met->phi; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Etmiss Rel"); {
        *cutflow << HFTname("metrel");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->tools->getMetRel(sl->met, *sl->leptons, *sl->jets) * GeV_to_MeV ; };
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
            return local_ll.M() * GeV_to_MeV ;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Pt of di-lepton system, Pt_ll"); {
        *cutflow << HFTname("pTll");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return local_ll.Pt() * GeV_to_MeV ; };
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
            return ht * GeV_to_MeV;
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

            return mt2_ * GeV_to_MeV;
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
            double   cosTheta_b           = 0.0;
            TVector3 vBETA_z;
            TVector3 pT_CM;
            TVector3 vBETA_T_CMtoR;
            TVector3 vBETA_R;
            sl->tools->superRazor(*sl->leptons, sl->met, vBETA_z, pT_CM,
                                    vBETA_T_CMtoR, vBETA_R, SHATR, dphi_LL_vBETA_T,
                                    dphi_L1_L2, gamma_R, dphi_vBETA_R_vBETA_T,
                                    mDeltaR, cosThetaRp1, cosTheta_b);
            return SHATR * GeV_to_MeV;
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
            double   cosTheta_b           = 0.0;
            TVector3 vBETA_z;
            TVector3 pT_CM;
            TVector3 vBETA_T_CMtoR;
            TVector3 vBETA_R;
            sl->tools->superRazor(*sl->leptons, sl->met, vBETA_z, pT_CM,
                                    vBETA_T_CMtoR, vBETA_R, SHATR, dphi_LL_vBETA_T,
                                    dphi_L1_L2, gamma_R, dphi_vBETA_R_vBETA_T,
                                    mDeltaR, cosThetaRp1, cosTheta_b);
            return mDeltaR * GeV_to_MeV;
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
            double   cosTheta_b           = 0.0;
            TVector3 vBETA_z;
            TVector3 pT_CM;
            TVector3 vBETA_T_CMtoR;
            TVector3 vBETA_R;
            sl->tools->superRazor(*sl->leptons, sl->met, vBETA_z, pT_CM,
                                    vBETA_T_CMtoR, vBETA_R, SHATR, dphi_LL_vBETA_T,
                                    dphi_L1_L2, gamma_R, dphi_vBETA_R_vBETA_T,
                                    mDeltaR, cosThetaRp1, cosTheta_b);
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
            double   cosTheta_b           = 0.0;
            TVector3 vBETA_z;
            TVector3 pT_CM;
            TVector3 vBETA_T_CMtoR;
            TVector3 vBETA_R;
            sl->tools->superRazor(*sl->leptons, sl->met, vBETA_z, pT_CM,
                                    vBETA_T_CMtoR, vBETA_R, SHATR, dphi_LL_vBETA_T,
                                    dphi_L1_L2, gamma_R, dphi_vBETA_R_vBETA_T,
                                    mDeltaR, cosThetaRp1, cosTheta_b);
            return dphi_LL_vBETA_T;
        };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("cosine(Theta_{b})"); {
        *cutflow << HFTname("cosTheta_b");
        *cutflow << [](Superlink* sl, var_float*) -> double {
            double   mDeltaR              = 0.0;
            double   SHATR                = 0.0;
            double   cosThetaRp1          = 0.0; 
            double   dphi_LL_vBETA_T      = 0.0;
            double   dphi_L1_L2           = 0.0;
            double   gamma_R              = 0.0;
            double   dphi_vBETA_R_vBETA_T = 0.0;
            double   cosTheta_b           = 0.0;
            TVector3 vBETA_z;
            TVector3 pT_CM;
            TVector3 vBETA_T_CMtoR;
            TVector3 vBETA_R;
            sl->tools->superRazor(*sl->leptons, sl->met, vBETA_z, pT_CM,
                                    vBETA_T_CMtoR, vBETA_R, SHATR, dphi_LL_vBETA_T,
                                    dphi_L1_L2, gamma_R, dphi_vBETA_R_vBETA_T,
                                    mDeltaR, cosThetaRp1, cosTheta_b);
            return cosTheta_b;
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
            return sl->met->Et/mEff * GeV_to_MeV;
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

    // END Setup output trees
    // END Setup output trees
    // END Setup output trees

    // GAP //
    // GAP //
    // GAP //
    
    // START Setup systematics
    // START Setup systematics
    // START Setup systematics

    // weight variation systematics
    *cutflow << NewSystematic("shift in electron fake rate"); {
        *cutflow << WeightSystematic(SupersysWeight::ELFRUP, SupersysWeight::ELFRDOWN);
        *cutflow << TreeName("ELFR");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in electron real efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::ELREUP, SupersysWeight::ELREDOWN);
        *cutflow << TreeName("ELRE");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in muon fake rate"); {
        *cutflow << WeightSystematic(SupersysWeight::MUFRUP, SupersysWeight::MUFRDOWN);
        *cutflow << TreeName("MUFR");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in muon real efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::MUREUP, SupersysWeight::MUREDOWN);
        *cutflow << TreeName("MURE");
        *cutflow << SaveSystematic();
    }


    // Initialize the cutflow and start the event loop.
    chain->Process(cutflow, sample_.c_str(), num_events_, n_skip_);

    delete cutflow;
    delete chain;

    cout << "Done." << endl;
    exit(0);
}

void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_,
                  SuperflowRunMode& run_mode_, SusyNtSys& nt_sys, string& fake_region_)
{
    bool nominal_ = false;
    bool nominal_and_weight_syst_ = false;
    bool all_syst_ = false;
    bool single_event_syst_ = false;

    bool do_fakes_ = false;

    string systematic_ = "undefined";

    string input = "";

    /** Read inputs to program */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "/n") == 0)
            num_events_ = atoi(argv[++i]);
        else if (strcmp(argv[i], "/c") == 0)
            nominal_ = true;
        else if (strcmp(argv[i], "/w") == 0)
            nominal_and_weight_syst_ = true;
        else if (strcmp(argv[i], "/e") == 0) {
            single_event_syst_ = true;
        }
        else if (strcmp(argv[i], "/a") == 0) {
            all_syst_ = true;
        }
        else if (strcmp(argv[i], "/f") == 0) {
            do_fakes_ = true;
            fake_region_ = argv[++i];
        }
        else if (strcmp(argv[i], "/i") == 0) {
            input = argv[++i];
        }
        else if (strcmp(argv[i], "/s") == 0) {
            systematic_ = argv[++i];
        }
        else if (strcmp(argv[i], "/t") == 0) {
            use2dParam_ = true;
        }
        else {
            cout << "Analysis    Error (fatal): Bad arguments." << endl;
            exit(1);
        }
    }

    bool inputIsFile = susy::utils::endswith(input, ".root");
    bool inputIsList = susy::utils::endswith(input, ".txt");
    bool inputIsDir = susy::utils::endswith(input, "/");
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

    if (do_fakes_) {
        run_mode_ = SuperflowRunMode::fakes;
        cout << "Analysis    run mode: SuperflowRunMode::fakes, region: " << fake_region_ << endl;
    }

    map <string, SusyNtSys> event_syst_map;
    event_syst_map["EESZUP"] = NtSys_EES_Z_UP;
    event_syst_map["EESZDOWN"] = NtSys_EES_Z_DN;
    event_syst_map["EESMATUP"] = NtSys_EES_MAT_UP;
    event_syst_map["EESMATDOWN"] = NtSys_EES_MAT_DN;
    event_syst_map["EESPSUP"] = NtSys_EES_PS_UP;
    event_syst_map["EESPSDOWN"] = NtSys_EES_PS_DN;
    event_syst_map["EESLOWUP"] = NtSys_EES_LOW_UP;
    event_syst_map["EESLOWDOWN"] = NtSys_EES_LOW_DN;
    event_syst_map["EERUP"] = NtSys_EER_UP;
    event_syst_map["EERDOWN"] = NtSys_EER_DN;
    event_syst_map["MSUP"] = NtSys_MS_UP;
    event_syst_map["MSDOWN"] = NtSys_MS_DN;
    event_syst_map["IDUP"] = NtSys_ID_UP;
    event_syst_map["IDDOWN"] = NtSys_ID_DN;
    event_syst_map["JESUP"] = NtSys_JES_UP;
    event_syst_map["JESDOWN"] = NtSys_JES_DN;
    event_syst_map["JER"] = NtSys_JER;
    event_syst_map["SCALESTUP"] = NtSys_SCALEST_UP;
    event_syst_map["SCALESTDOWN"] = NtSys_SCALEST_DN;
    event_syst_map["RESOST"] = NtSys_RESOST;
    event_syst_map["TRIGSFELUP"] = NtSys_TRIGSF_EL_UP;
    event_syst_map["TRIGSFELDN"] = NtSys_TRIGSF_EL_DN;
    event_syst_map["TRIGSFMUUP"] = NtSys_TRIGSF_MU_UP;
    event_syst_map["TRIGSFMUDN"] = NtSys_TRIGSF_MU_DN;
    event_syst_map["TESUP"] = NtSys_TES_UP;
    event_syst_map["TESDOWN"] = NtSys_TES_DN;
    event_syst_map["JVFUP"] = NtSys_JVF_UP;
    event_syst_map["JVFDOWN"] = NtSys_JVF_DN;

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


// List of event systematics (for scripting)
// 
// "EESZUP",
// "EESZDOWN",
// "EESMATUP",
// "EESMATDOWN",
// "EESPSUP",
// "EESPSDOWN",
// "EESLOWUP",
// "EESLOWDOWN",
// "EERUP",
// "EERDOWN",
// "MSUP",
// "MSDOWN",
// "IDUP",
// "IDDOWN",
// "JESUP",
// "JESDOWN",
// "JER",
// "SCALESTUP",
// "SCALESTDOWN",
// "RESOST",
// "TRIGSFELUP",
// "TRIGSFELDN",
// "TRIGSFMUUP",
// "TRIGSFMUDN",
// "TESUP",
// "TESDOWN",
// "JVFUP",
// "JVFDOWN",
// 


// A SELECTION OF CUTS and Variables
// A SELECTION OF CUTS and Variables
// A SELECTION OF CUTS and Variables
//
// *cutflow << CutName("exactly two base leptons") << [](Superlink* sl) -> bool {
//     return sl->baseLeptons->size() == 2;
// };
//
// *cutflow << CutName("m_ll > 20 GeV") << [](Superlink* sl) -> bool {
//     double m_ll = (*sl->baseLeptons->at(0) + *sl->baseLeptons->at(1)).M();
//     return m_ll > 20.0;
// };
//
// *cutflow << CutName("is MM") << [](Superlink* sl) -> bool {
//     return sl->baseLeptons->at(0)->isMu() && sl->baseLeptons->at(1)->isMu();
// }; // debug only !!!
//
// *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
//     return sl->leptons->size() == 2;
// };
//
// *cutflow << NewVar("mCT"); {
//     *cutflow << HFTname("mct");
//     *cutflow << [&](Superlink* sl, var_float*) -> double { please return PhysicsTools::mCT(*sl->leptons->at(0), *sl->leptons->at(1)) * GeV_to_MeV; };
//     *cutflow << SaveVar();
// }
// 
// *cutflow << NewVar("mCT perpendicular"); {
//     *cutflow << HFTname("mctPerp");
//     *cutflow << [&](Superlink* sl, var_float*) -> double { return PhysicsTools::mCTperp(*sl->leptons->at(0), *sl->leptons->at(1), sl->met->lv()) * GeV_to_MeV; };
//     *cutflow << SaveVar();
// }
//
// *cutflow << NewSystematic("Trigger Scale factor + error for el"); {
//     *cutflow << EventSystematic(NtSys_TRIGSF_EL_UP);
//     *cutflow << TreeName("TRIGSFELUP");
//     *cutflow << SaveSystematic();
// } // These are correctly included as event variation systematics
// 
// *cutflow << NewSystematic("Trigger Scale factor - error for el"); {
//     *cutflow << EventSystematic(NtSys_TRIGSF_EL_DN);
//     *cutflow << TreeName("TRIGSFELDN");
//     *cutflow << SaveSystematic();
// } // These are correctly included as event variation systematics
// 
// *cutflow << NewSystematic("Trigger Scale factor + error for mu"); {
//     *cutflow << EventSystematic(NtSys_TRIGSF_MU_UP);
//     *cutflow << TreeName("TRIGSFMUUP");
//     *cutflow << SaveSystematic();
// } // These are correctly included as event variation systematics
// 
// *cutflow << NewSystematic("Trigger Scale factor - error for mu"); {
//     *cutflow << EventSystematic(NtSys_TRIGSF_MU_DN);
//     *cutflow << TreeName("TRIGSFMUDN");
//     *cutflow << SaveSystematic();
// } // These are correctly included as event variation systematics
//
// A SELECTION OF CUTS and Variables
// A SELECTION OF CUTS and Variables
// A SELECTION OF CUTS and Variables


//rem// Cut* is_e_mu = new IsEMu(); // initialize cut
//rem// Cut* is_two_lepton = new Is2Lepton();
//rem// Cut* is_e_e_and_os = new IsEEOS();
//rem// *cutflow << (new IsEEOS()); // push cut to cutflow
//rem// class Is2Lepton : public Cut {
//rem// public:
//rem//     Is2Lepton()
//rem//     {
//rem//         name = "base is 2-lepton";
//rem//     }
//rem//     bool operator() (Superlink* sl) // return true to pass the cut
//rem//     {
//rem//         return sl->baseLeptons->size() == 2; // exactly two base leptons
//rem//     }
//rem// };
//rem// 
//rem// class IsEMu : public Cut {
//rem// public:
//rem//     IsEMu()
//rem//     {
//rem//         name = "base is e + mu";
//rem//     }
//rem//     bool operator() (Superlink* sl) // return true to pass the cut
//rem//     {
//rem//         if (sl->baseLeptons->size() == 2) { // exactly two base leptons
//rem//             return sl->baseLeptons->at(0)->isEle() ^ sl->baseLeptons->at(1)->isEle(); // e + mu
//rem//         }
//rem//         else {
//rem//             return false;
//rem//         }
//rem//     }
//rem// };
//rem// 
//rem// class IsEEOS : public Cut {
//rem// public:
//rem//     IsEEOS()
//rem//     {
//rem//         name = "base is e + e + OS";
//rem//     }
//rem//     bool operator() (Superlink* sl) // return true to pass the cut
//rem//     {
//rem//         if (sl->baseLeptons->size() == 2) { // exactly two base leptons
//rem//             return (sl->baseLeptons->at(0)->q * sl->baseLeptons->at(1)->q < 0)
//rem//                 && (sl->baseLeptons->at(0)->isEle() && sl->baseLeptons->at(1)->isEle()); // e + mu
//rem//         }
//rem//         else {
//rem//             return false;
//rem//         }
//rem//     }
//rem// };


/*
*cutflow << CutName("all cleaning") << [](Superlink* sl) -> bool {
bool pass_cleaning = true;
int cutFlags = sl->tools->cleaningCutFlags(sl->nt->evt()->cutFlags[sl->nt_sys], *sl->preMuons, *sl->baseMuons, *sl->preJets, *sl->baseJets);

pass_cleaning &=
!sl->nt->evt()->eventWithSusyProp
&& sl->tools->passGRL(cutFlags)
&& sl->tools->passTileTripCut(cutFlags)
&& sl->tools->passLarErr(cutFlags)
;

if (pass_cleaning) {
JetVector jets = sl->tools->getPreJets(sl->nt, sl->nt_sys);
sl->tools->e_j_overlap(*sl->baseElectrons, jets, J_E_DR, true);
sl->tools->t_j_overlap(*sl->taus, jets, J_T_DR, true);

pass_cleaning &=
!sl->tools->hasBadJet(jets)
&& sl->tools->passDeadRegions(*sl->preJets, sl->met, sl->nt->evt()->run, sl->nt->evt()->isMC)
&& !sl->tools->hasBadMuon(*sl->preMuons)
&& !sl->tools->hasCosmicMuon(*sl->baseMuons)
&& !sl->tools->hasHotSpotJet(*sl->preJets)
&& sl->tools->passTTCVeto(cutFlags)
&& sl->tools->passGoodVtx(cutFlags)
&& sl->baseLeptons->size() == 2
&& (*sl->baseLeptons->at(0) + *sl->baseLeptons->at(1)).M() > 20.0
;
}

return pass_cleaning;
};
*/
