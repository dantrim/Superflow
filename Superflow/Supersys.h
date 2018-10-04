#pragma once

#include <string>

#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/SusyNtSys.h"

using namespace std;

namespace sflow {

    ////////////////////////////////////////////////////////
    // variables for defining output tree names
    ////////////////////////////////////////////////////////
    const string weight_prefix = "syst_";
    const string weight_suffix_up = "UP";
    const string weight_suffix_down = "DOWN";

    //////////////////////////////////////////////////////
    enum class SupersysType {
        central,
        event,
        weight,
        null
    };

    //////////////////////////////////////////////////////
    enum class SupersysState {
        closed,
        open
    };

    //////////////////////////////////////////////////////
    enum class SupersysWeight {
        // electron
        EL_EFF_ChargeIDSel_DN,
        EL_EFF_ChargeIDSel_UP,
        EL_EFF_ID_TOTAL_Uncorr_DN,
        EL_EFF_ID_TOTAL_Uncorr_UP,
        EL_EFF_Iso_TOTAL_Uncorr_DN,
        EL_EFF_Iso_TOTAL_Uncorr_UP,
        EL_EFF_Reco_TOTAL_Uncorr_DN,
        EL_EFF_Reco_TOTAL_Uncorr_UP,
        EL_EFF_TriggerEff_TOTAL_DN,
        EL_EFF_TriggerEff_TOTAL_UP,
        EL_EFF_Trigger_TOTAL_DN,
        EL_EFF_Trigger_TOTAL_UP,
        // flavor taggin
        FT_EFF_B_UP,
        FT_EFF_B_DN,
        FT_EFF_C_UP,
        FT_EFF_C_DN,
        FT_EFF_LT_UP,
        FT_EFF_LT_DN,
        FT_EFF_EXTRAP_UP,
        FT_EFF_EXTRAP_DN,
        FT_EFF_EXTRAPC_UP,
        FT_EFF_EXTRAPC_DN,
        // jvt
        JET_JVTEff_UP,
        JET_JVTEff_DN,
        // muon
        MUON_EFF_BADMUON_STAT_DN,
        MUON_EFF_BADMUON_STAT_UP,
        MUON_EFF_BADMUON_SYS_DN,
        MUON_EFF_BADMUON_SYS_UP,
        MUON_EFF_ISO_STAT_DN,
        MUON_EFF_ISO_STAT_UP,
        MUON_EFF_ISO_SYS_DN,
        MUON_EFF_ISO_SYS_UP,
        MUON_EFF_RECO_STAT_DN,
        MUON_EFF_RECO_STAT_UP,
        MUON_EFF_RECO_SYS_DN,
        MUON_EFF_RECO_SYS_UP,
        MUON_EFF_RECO_STAT_LOWPT_DN,
        MUON_EFF_RECO_STAT_LOWPT_UP,
        MUON_EFF_RECO_SYS_LOWPT_DN,
        MUON_EFF_RECO_SYS_LOWPT_UP,
        MUON_EFF_TTVA_STAT_DN,
        MUON_EFF_TTVA_STAT_UP,
        MUON_EFF_TTVA_SYS_DN,
        MUON_EFF_TTVA_SYS_UP,
        MUON_EFF_TrigStat_DN,
        MUON_EFF_TrigStat_UP,
        MUON_EFF_TrigSys_DN,
        MUON_EFF_TrigSys_UP,
        // pileup reweighting
        PILEUP_UP,
        PILEUP_DN,
        null,
        block
    };
    

    //////////////////////////////////////////////////////
    class Supersys {
        public:
                Supersys();
                Supersys(SupersysType type_);

                void reset();

                std::string name;
                std::string tree_name;

                double weight_up;
                double weight_down;

                SupersysType type;

                Susy::NtSys::SusyNtSys event_syst; /// The event systematics are from SusyDefs
                SupersysWeight weight_syst;
                SupersysWeight weight_syst_up; // The weight systematics are given above.
                SupersysWeight weight_syst_down; // The weight systematics are given above.
    };

    //////////////////////////////////////////////////////
    class NewSystematic {
        public:
            NewSystematic(std::string sys_name);

            std::string name;
    };


    //////////////////////////////////////////////////////
    class TreeName {
        public:
            TreeName(std::string var_name);

            std::string name;
    };

    //////////////////////////////////////////////////////
    class EventSystematic {
        public:
            EventSystematic(Susy::NtSys::SusyNtSys syst_val);

            Susy::NtSys::SusyNtSys event_syst_;
    };

    //////////////////////////////////////////////////////
    class WeightSystematic {
        public:
            WeightSystematic(SupersysWeight syst_val_up, SupersysWeight syst_val_down);

            SupersysWeight weight_syst_up;
            SupersysWeight weight_syst_down;
    };

    //////////////////////////////////////////////////////
    class SaveSystematic {
        public:
            SaveSystematic() {};
    };

    //////////////////////////////////////////////////////
    BTagSys supersys_to_btagsys(SupersysWeight sys);

};
