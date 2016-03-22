// rjigsawAna.cxx


// std
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>

// ROOT
#include "TChain.h"
#include "TVectorD.h"
#include "TRandom.h"

// SusyNtuple
#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"
#include "SusyNtuple/TriggerTools.h"
#include "SusyNtuple/KinematicTools.h"

// Superflow
#include "Superflow/Superflow.h"
#include "Superflow/Superlink.h"
#include "Superflow/Cut.h"


using namespace std;
using namespace sflow;

const string analysis_name = "newSysAna";

////////////////////////////////////////
// Function Prototypes
////////////////////////////////////////
void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_, SuperflowRunMode& run_mode_, SusyNtSys& nt_sys_, bool dbg);


////////////////////////////////////////
// MAIN
////////////////////////////////////////
int main(int argc, char* argv[])
{
    ///////////////////////
    // READ IN
    ///////////////////////
    int n_skip_ = 0;
    int num_events_ = -1;
    bool m_dbg = false;
    string sample_;
    SuperflowRunMode run_mode_ = SuperflowRunMode::nominal;
    SusyNtSys nt_sys_ = NtSys::NOM;

    TChain* chain = new TChain("susyNt");
    chain->SetDirectory(0);

    read_options(argc, argv, chain, n_skip_, num_events_, sample_, run_mode_, nt_sys_, m_dbg);

    ////////////////////////////////////////////////////
    // Construct and configure the Superflow object
    ////////////////////////////////////////////////////
    Superflow* cutflow = new Superflow();
    cutflow->setAnaName(analysis_name);
    cutflow->setAnaType(AnalysisType::Ana_Stop2L);
    cutflow->setLumi(3209); // 3.34/fb
    cutflow->setSampleName(sample_);
    cutflow->setRunMode(run_mode_);
    cutflow->setCountWeights(true);
    cutflow->setChain(chain);

    // print some useful
    cout << analysis_name << "    Total Entries    : " << chain->GetEntries() << endl;
    if(num_events_ > 0) {
    cout << analysis_name << "    Process Entries  : " << num_events_ << endl;
    } else {
    cout << analysis_name << "    Process Entries  : " << chain->GetEntries() << endl;
    }

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    //
    // Superflow methods [BEGIN]
    //
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    *cutflow << CutName("read in ") << [](Superlink* sl) -> bool { return true; };

    ////////////////////////////////////////////////////
    // Cleaning cuts
    ////////////////////////////////////////////////////
    int cutflags = 0;
    *cutflow << CutName("Pass GRL") << [&](Superlink* sl) -> bool {
        cutflags = sl->nt->evt()->cutFlags[NtSys::NOM];
        //cutflags = sl->nt->evt()->cutFlags[sl->nt_sys];
        return (sl->tools->passGRL(cutflags));
    };
    *cutflow << CutName("LAr error") << [&](Superlink* sl) -> bool {
        return (sl->tools->passLarErr(cutflags));
    };
    *cutflow << CutName("Tile Error") << [&](Superlink* sl) -> bool {
        return (sl->tools->passTileErr(cutflags));
    };
    *cutflow << CutName("TTC veto") << [&](Superlink* sl) -> bool {
        return (sl->tools->passTTC(cutflags));
    };
    *cutflow << CutName("pass Good Vertex") << [&](Superlink * sl) -> bool {
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
    ///////////////////////////////////////////////////
    // Analysis Cuts
    ///////////////////////////////////////////////////
    *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
        return sl->leptons->size() == 2;
    };

    *cutflow << CutName("opposite sign") << [](Superlink* sl) -> bool {
        return ((sl->leptons->at(0)->q * sl->leptons->at(1)->q) < 0);
    };
    

    ///////////////////////////////////////////////////
    // Ntuple Setup
    ///////////////////////////////////////////////////

    // standard variables
    *cutflow << NewVar("event weight"); {
        *cutflow << HFTname("eventweight");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product();
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

    *cutflow << NewVar("is MC"); {
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

    // fill containers
    // fill containers
    // fill containers
    //
    Met met;
    *cutflow << [&](Superlink* sl, var_void*) { met = *sl->met; };

    LeptonVector leptons;
    ElectronVector electrons;
    MuonVector muons;
    *cutflow << [&](Superlink* sl, var_void*) { leptons = *sl->leptons; };
    *cutflow << [&](Superlink* sl, var_void*) { electrons = *sl->electrons; };
    *cutflow << [&](Superlink* sl, var_void*) { muons = *sl->muons; };

    JetVector jets;
    JetVector bjets;
    JetVector sjets;
    *cutflow << [&](Superlink* sl, var_void*) { jets = *sl->jets; };
    *cutflow << [&](Superlink* sl, var_void*) {
        for(int i = 0; i < jets.size(); i++) {
            Jet* j = jets[i];
            if(sl->tools->jetSelector().isB(j)) bjets.push_back(j);
            if(!sl->tools->jetSelector().isB(j)) sjets.push_back(j);
        } //i
    };

    // lepton variables
    // lepton variables
    // lepton variables
    

    *cutflow << NewVar("number of leptons"); {
        *cutflow << HFTname("nLeptons");
        *cutflow << [&](Superlink* sl, var_int*) -> int { return leptons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of electrons"); {
        *cutflow << HFTname("nElectrons");
        *cutflow << [&](Superlink* sl, var_int*) -> int { return electrons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of muons"); {
        *cutflow << HFTname("nMuons");
        *cutflow << [&](Superlink* sl, var_int*) -> int { return muons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead lep is E"); {
        *cutflow << HFTname("leadIsE");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return (leptons.at(0)->isEle() ? 1 : 0);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lep1q"); {
        *cutflow << HFTname("l0q");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return leptons.at(0)->q;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lep2q"); {
        *cutflow << HFTname("l1q");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return leptons.at(1)->q;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lept1pt"); {
        *cutflow << HFTname("l0pt");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return leptons.at(0)->Pt();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lept2pt"); {
        *cutflow << HFTname("l1pt");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return leptons.at(1)->Pt();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lept1eta"); {
        *cutflow << HFTname("l0eta");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return leptons.at(0)->Eta();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lept2eta"); {
        *cutflow << HFTname("l1eta");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return leptons.at(1)->Eta();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("nJets"); {
        *cutflow << HFTname("nJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return (int)jets.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("nBJets"); {
        *cutflow << HFTname("nBJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return (int)bjets.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("j0pt"); {
        *cutflow << HFTname("j0pt");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return (jets.size()>0 ? jets.at(0)->Pt() : -999.);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("j0eta"); {
        *cutflow << HFTname("j0eta");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return (jets.size()>0 ? jets.at(0)->Eta() : -999.);
        };
        *cutflow << SaveVar();
    }

    // clear the wectors
    *cutflow << [&](Superlink* sl, var_void*) { leptons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { electrons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { muons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { jets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { bjets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { sjets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { met.clear(); };

    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //
    // Sysystematics [BEGIN]
    //
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////
    // weight systematics
    ////////////////////////////////////

    // electron eff
    *cutflow << NewSystematic("shift in electron ID efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_ID_UP, SupersysWeight::EL_EFF_ID_DN);
        *cutflow << TreeName("EL_EFF_ID");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in electron ISO efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_ISO_UP, SupersysWeight::EL_EFF_ISO_DN);
        *cutflow << TreeName("EL_EFF_Iso");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in electron RECO efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_RECO_UP, SupersysWeight::EL_EFF_RECO_DN);
        *cutflow << TreeName("EL_EFF_Reco");
        *cutflow << SaveSystematic();
    }

    // muon eff
    *cutflow << NewSystematic("muon eff stat uncertainty"); {
        *cutflow << WeightSystematic(SupersysWeight::MU_EFF_STAT_UP, SupersysWeight::MU_EFF_STAT_DN);
        *cutflow << TreeName("MUON_EFF_STAT");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff stat uncertainty (low pt)"); {
        *cutflow << WeightSystematic(SupersysWeight::MU_EFF_STAT_LOWPT_UP, SupersysWeight::MU_EFF_STAT_LOWPT_DN);
        *cutflow << TreeName("MUON_EFF_STAT_LOWPT");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff syst uncertainty"); {
        *cutflow << WeightSystematic(SupersysWeight::MU_EFF_SYS_UP, SupersysWeight::MU_EFF_SYS_DN);
        *cutflow << TreeName("MUON_EFF_SYS");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff syst uncertainty (low pt"); {
        *cutflow << WeightSystematic(SupersysWeight::MU_EFF_SYS_LOWPT_UP, SupersysWeight::MU_EFF_SYS_LOWPT_DN);
        *cutflow << TreeName("MUON_EFF_SYS_LOWPT");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff iso stat uncertainty"); {
        *cutflow << WeightSystematic(SupersysWeight::MU_EFF_ISO_STAT_UP, SupersysWeight::MU_EFF_ISO_STAT_DN);
        *cutflow << TreeName("MUON_ISO_STAT");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff iso syst uncertainty"); {
        *cutflow << WeightSystematic(SupersysWeight::MU_EFF_ISO_SYS_UP, SupersysWeight::MU_EFF_ISO_SYS_DN);
        *cutflow << TreeName("MUON_ISO_SYS");
        *cutflow << SaveSystematic();
    }

    // flavor tagging eff
    *cutflow << NewSystematic("shift in b-tag efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_B_UP, SupersysWeight::FT_EFF_B_DN);
        *cutflow << TreeName("FT_EFF_B");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in c-tag efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_C_UP, SupersysWeight::FT_EFF_C_DN);
        *cutflow << TreeName("FT_EFF_C");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in light tag (i.e. mis-tag) efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_LT_UP, SupersysWeight::FT_EFF_LT_DN);
        *cutflow << TreeName("FT_EFF_Light");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift flavor tagging extrapolation ?"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_EXTRAP_UP, SupersysWeight::FT_EFF_EXTRAP_DN);
        *cutflow << TreeName("FT_EFF_extrapolation");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift flavor tagging extrapolation - charm ?"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_EXTRAPC_UP, SupersysWeight::FT_EFF_EXTRAPC_DN);
        *cutflow << TreeName("FT_EFF_extrapolation_charm");
        *cutflow << SaveSystematic();
    }

    // jvt eff
    *cutflow << NewSystematic("shift in JVT efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::JVT_EFF_UP, SupersysWeight::JVT_EFF_DN);
        *cutflow << TreeName("JET_JVTEff");
        *cutflow << SaveSystematic();
    }

    // pileup
    *cutflow << NewSystematic("shift in data mu (pile-up)"); {
        *cutflow << WeightSystematic(SupersysWeight::PILEUP_UP, SupersysWeight::PILEUP_DN);
        *cutflow << TreeName("PILEUP");
        *cutflow << SaveSystematic();
    }


    ////////////////////////////////////
    // shape systematics
    ////////////////////////////////////

    // egamma
    *cutflow << NewSystematic("shift in e-gamma resolution (UP)"); {
        *cutflow << EventSystematic(NtSys::EG_RESOLUTION_ALL_UP);
        *cutflow << TreeName("EG_RESOLUTION_ALL_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in e-gamma resolution (DOWN)"); {
        *cutflow << EventSystematic(NtSys::EG_RESOLUTION_ALL_DN);
        *cutflow << TreeName("EG_RESOLUTION_ALL_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in e-gamma scale (UP)"); {
        *cutflow << EventSystematic(NtSys::EG_SCALE_ALL_UP);
        *cutflow << TreeName("EG_SCALE_ALL_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in e-gamma scale (DOWN)"); {
        *cutflow << EventSystematic(NtSys::EG_SCALE_ALL_DN);
        *cutflow << TreeName("EG_SCALE_ALL_DN");
        *cutflow << SaveSystematic();
    }
    // muon
    *cutflow << NewSystematic("muon ID (UP)"); {
        *cutflow << EventSystematic(NtSys::MUONS_ID_UP);
        *cutflow << TreeName("MUONS_ID_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon ID (DOWN)"); {
        *cutflow << EventSystematic(NtSys::MUONS_ID_DN);
        *cutflow << TreeName("MUONS_ID_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon MS (UP)"); {
        *cutflow << EventSystematic(NtSys::MUONS_MS_UP);
        *cutflow << TreeName("MUONS_MS_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon MS (DOWN)"); {
        *cutflow << EventSystematic(NtSys::MUONS_MS_DN);
        *cutflow << TreeName("MUONS_MS_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon scale shift (UP)"); {
        *cutflow << EventSystematic(NtSys::MUONS_SCALE_UP);
        *cutflow << TreeName("MUONS_SCALE_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon scale shift (DN)"); {
        *cutflow << EventSystematic(NtSys::MUONS_SCALE_DN);
        *cutflow << TreeName("MUONS_SCALE_DN");
        *cutflow << SaveSystematic();
    }

    // jet
    *cutflow << NewSystematic("JER"); {
        *cutflow << EventSystematic(NtSys::JER);
        *cutflow << TreeName("JER");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 1 (up)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_1_UP);
        *cutflow << TreeName("JET_GroupedNP_1_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 1 (down)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_1_DN);
        *cutflow << TreeName("JET_GroupedNP_1_DN");
        *cutflow << SaveSystematic();
    }

    // met
    *cutflow << NewSystematic("MET TST Soft-Term resolution (parallel)"); {
        *cutflow << EventSystematic(NtSys::MET_SoftTrk_ResoPara);
        *cutflow << TreeName("MET_SoftTrk_ResoPara");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("MET TST Soft-Term resolution (perpendicular)"); {
        *cutflow << EventSystematic(NtSys::MET_SoftTrk_ResoPerp);
        *cutflow << TreeName("MET_SoftTrk_ResoPerp");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("MET TST Soft-Term shift in scale (UP)"); {
        *cutflow << EventSystematic(NtSys::MET_SoftTrk_ScaleUp);
        *cutflow << TreeName("MET_SoftTrk_ScaleUp");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("MET TST Soft-Term shift in scale (DOWN)"); {
        *cutflow << EventSystematic(NtSys::MET_SoftTrk_ScaleDown);
        *cutflow << TreeName("MET_SoftTrk_ScaleDown");
        *cutflow << SaveSystematic();
    }

    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //
    // Superflow methods [END]
    //
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    // initialize the cutflow and start the event loop
    chain->Process(cutflow, sample_.c_str(), num_events_, n_skip_);
    delete cutflow;
    delete chain;
    cout << "La Fin." << endl;
    exit(0);

}

void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_,
    SuperflowRunMode& run_mode_, SusyNtSys& nt_sys, bool dbg)
{
    bool nominal_ = false;
    bool nominal_and_weight_sys_ = false;
    bool all_sys_ = false;
    string input;

    for(int i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-n") == 0)
            num_events_ = atoi(argv[++i]);
        else if(strcmp(argv[i], "-i") == 0) {
            input = argv[++i];
        }
        else if(strcmp(argv[i], "-c") == 0) {
            nominal_ = true;
        }
        else if(strcmp(argv[i], "-w") == 0) {
            nominal_and_weight_sys_ = true;
        }
        else if(strcmp(argv[i], "-a") == 0) {
            all_sys_ = true;
        }
        else if(strcmp(argv[i], "-d") == 0) {
            dbg = true;
        }
        else {
            cout << analysis_name << "    Error (fatal) : Bad arguments." << endl;
            exit(1);
        }
    } // i

    cout << analysis_name << "    input: " << input << endl;
    cout << analysis_name << "    input: " << input << endl;
    cout << analysis_name << "    input: " << input << endl;
    ChainHelper::addInput(chain, input, dbg); 

    Long64_t tot_num_events = chain->GetEntries();
    num_events_ = (num_events_ < 0 ? tot_num_events : num_events_);
    if(nominal_) {
        run_mode_ = SuperflowRunMode::nominal;
        cout << analysis_name << "    run mode: SuperflowRunMode::nominal" << endl;
    }
    if(nominal_and_weight_sys_) {
        run_mode_ = SuperflowRunMode::nominal_and_weight_syst;
        cout << analysis_name << "    run mode: SuperflowRunMode::nominal_and_weight_syst" << endl;
    }
    if(all_sys_) {
        run_mode_ = SuperflowRunMode::all_syst;
        cout << analysis_name << "    run mode: SuperflowRunMode::all_syst" << endl;
    }
}
