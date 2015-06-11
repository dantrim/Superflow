// SuperflowAna.cxx
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
    SusyNtSys nt_sys_ = NtSys_NOM;

    TChain* chain = new TChain("susyNt");
    chain->SetDirectory(0);

    read_options(argc, argv, chain, n_skip_, num_events_, sample_, run_mode, nt_sys_); // defined below
    // END read-in

    Superflow* cutflow = new Superflow(); // initialize the cutflow
    cutflow->setAnaType(Ana_2Lep); // Ana_2Lep Ana_2LepWH 
    cutflow->setSampleName(sample_);
    cutflow->setRunMode(run_mode);
    cutflow->setChain(chain);
    cutflow->setCountWeights(true);

    cout << "Analysis    Total Entries: " << chain->GetEntries() << endl;

    if (run_mode == SuperflowRunMode::single_event_syst) cutflow->setSingleEventSyst(nt_sys_);

    // START Setup cuts
    // START Setup cuts
    // START Setup cuts

    *cutflow << CutName("read in") << [](Superlink* sl) -> bool { return true; };

    *cutflow << CutName("at least two signal leptons") << [](Superlink* sl) -> bool {
        return !(sl->leptons->size() < 2);
    };

    *cutflow << CutName("pass HFOR") << [](Superlink* sl) -> bool {
        return sl->nt->evt()->hfor != 4
            || 
            (
            (sl->nt->evt()->mcChannel >= 164440 && sl->nt->evt()->mcChannel <= 164443)
            ||
            (sl->nt->evt()->mcChannel >= 164450 && sl->nt->evt()->mcChannel <= 164453)
            );
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

    *cutflow << CutName("bad jets") << [](Superlink* sl) -> bool {
        JetVector jets = sl->tools->getPreJets(sl->nt, sl->nt_sys);
        sl->tools->e_j_overlap(*sl->baseElectrons, jets, J_E_DR, true);
        sl->tools->t_j_overlap(*sl->taus, jets, J_T_DR, true);
        return !sl->tools->hasBadJet(jets);
    };

    *cutflow << CutName("dead regions") << [](Superlink* sl) -> bool {
        return sl->tools->passDeadRegions(*sl->preJets, sl->met, sl->nt->evt()->run, sl->nt->evt()->isMC);
    };

    *cutflow << CutName("bad muons") << [](Superlink* sl) -> bool {
        return !sl->tools->hasBadMuon(*sl->preMuons);
    };

    *cutflow << CutName("cosmic muons") << [](Superlink* sl) -> bool {
        return !sl->tools->hasCosmicMuon(*sl->baseMuons);
    };

    *cutflow << CutName("hotspot jets") << [](Superlink* sl) -> bool {
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

    *cutflow << CutName("is MM") << [](Superlink* sl) -> bool {
        return sl->baseLeptons->at(0)->isMu() && sl->baseLeptons->at(1)->isMu(); // 0 isMu()
    }; // debug only !!!

    // *cutflow << CutName("is El + Mu (any)") << [](Superlink* sl) -> bool {
    //     return sl->baseLeptons->at(0)->isEle() ^ sl->baseLeptons->at(1)->isEle();
    // }; // debug only !!!

    *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
        please return sl->leptons->size() == 2;
    }; // First cut to help speed.

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

    *cutflow << CutName("prompt leptons") << [](Superlink* sl) -> bool {
        bool pass_ = true;

        if (sl->isMC && !sl->doFake) {
            for (int l_ = 0; l_ < sl->leptons->size(); l_++) {
                bool isReal = sl->leptons->at(l_)->truthType == LeptonTruthType::PROMPT;

                bool isChargeFlip = sl->leptons->at(l_)->isEle()
                    && static_cast<Electron*>(sl->leptons->at(l_))->isChargeFlip;

                if (!isReal || isChargeFlip) pass_ = false;
            }
        }
        return pass_;
    };

    *cutflow << CutName("opposite sign") << [](Superlink* sl) -> bool {
        return (sl->leptons->at(0)->q * sl->leptons->at(1)->q < 0);
    };

    /*
    *cutflow << CutName("read in") << [](Superlink* sl) -> bool {
    return true;
    };

    *cutflow << CutName("exactly two base leptons") << [](Superlink* sl) -> bool {
    return sl->baseLeptons->size() == 2;
    };

    *cutflow << CutName("m_ll > 20 GeV") << [](Superlink* sl) -> bool {
    return (*sl->baseLeptons->at(0) + *sl->baseLeptons->at(1)).M() > 20.0;
    };

    *cutflow << CutName("is ME") << [](Superlink* sl) -> bool {
    return sl->baseLeptons->at(0)->isMu() && sl->baseLeptons->at(1)->isEle();
    }; // debug only !!!

    // *cutflow << CutName("is El + Mu (any)") << [](Superlink* sl) -> bool {
    //     return sl->baseLeptons->at(0)->isEle() ^ sl->baseLeptons->at(1)->isEle();
    // }; // debug only !!!

    *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
    please return sl->leptons->size() == 2;
    }; // First cut to help speed.

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

    *cutflow << CutName("prompt leptons") << [](Superlink* sl) -> bool {
    bool pass_ = true;

    if (sl->isMC && !sl->doFake) {
    for (int l_ = 0; l_ < sl->leptons->size(); l_++) {
    bool isReal = sl->leptons->at(l_)->truthType == LeptonTruthType::PROMPT;

    bool isChargeFlip = sl->leptons->at(l_)->isEle()
    && static_cast<Electron*>(sl->leptons->at(l_))->isChargeFlip;

    if (!isReal || isChargeFlip) pass_ = false;
    }
    }
    return pass_;
    };

    *cutflow << CutName("opposite sign") << [](Superlink* sl) -> bool {
    return (sl->leptons->at(0)->q * sl->leptons->at(1)->q < 0);
    };

    */

    /*
    llType==2 || llType==3			49354.36 +/- 49354.36
    nBJets>0			43920.60 +/- 43920.60
    nFJets==0			35069.12 +/- 35069.12
    l_pt[0]>30			32581.07 +/- 32581.07
    l_pt[1]>18
    */



    *cutflow << CutName("number of central b jets > 0") << [](Superlink* sl) -> bool {
        return SusyNtTools::numberOfCBJets(*sl->jets) > 0;
    };

    *cutflow << CutName("number of forward jets == 0") << [](Superlink* sl) -> bool {
        return SusyNtTools::numberOfFJets(*sl->jets) == 0;
    };

    *cutflow << CutName("leading lepton Pt > 30 GeV") << [](Superlink* sl) -> bool {
        return sl->leptons->at(0)->Pt() > 30.0;
    };

    *cutflow << CutName("subleading lepton Pt > 18 GeV") << [](Superlink* sl) -> bool {
        please return sl->leptons->at(1)->Pt() > 18.0;
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
        else if (strcmp(argv[i], "/i") == 0) {
            input = argv[++i];
        }
        else if (strcmp(argv[i], "/s") == 0) {
            systematic_ = argv[++i];
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