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
        EL_EFF_ID_UP,           //EL_EFF_ID_TotalCorrUncertainty_UP
        EL_EFF_ID_DN,           //EL_EFF_ID_TotalCorrUncertainty_DN
        EL_EFF_ISO_UP,          //EL_EFF_Iso_TotalCorrUncertainty_UP
        EL_EFF_ISO_DN,          //EL_EFF_Iso_TotalCorrUncertainty_DN
        EL_EFF_RECO_UP,         //EL_EFF_Reco_TotalCorrUncertainty_UP
        EL_EFF_RECO_DN,         //EL_EFF_Reco_TotalCorrUncertainty_DN
        EL_EFF_TRIG_UP,         //EL_EFF_Trigger_TOTAL_UP
        EL_EFF_TRIG_DN,         //EL_EFF_Trigger_TOTAL_DN
        // flavor tagging
        FT_EFF_B_UP,            //FT_EFF_B_systematics_UP
        FT_EFF_B_DN,            //FT_EFF_B_systematics_DN
        FT_EFF_C_UP,            //FT_EFF_C_systematics_UP
        FT_EFF_C_DN,            //FT_EFF_C_systematics_DN
        FT_EFF_LT_UP,           //FT_EFF_Light_systematics_UP
        FT_EFF_LT_DN,           //FT_EFF_Light_systematics_DN
        FT_EFF_EXTRAP_UP,       //FT_EFF_extrapolation_UP
        FT_EFF_EXTRAP_DN,       //FT_EFF_extrapolation_DN
        FT_EFF_EXTRAPC_UP,      //FT_EFF_extrapolation_charm_UP
        FT_EFF_EXTRAPC_DN,      //FT_EFF_extrapolation_charm_DN
        // jvt
        JVT_EFF_UP,             //JET_JVTEff_UP
        JVT_EFF_DN,             //JET_JVTEff_DN
        // muon
        MUON_EFF_STAT_UP,         //MUON_EFF_STAT_UP
        MUON_EFF_STAT_DN,         //MUON_EFF_STAT_DN
        MUON_EFF_STAT_LOWPT_UP,   //MUON_EFF_STAT_LOWPT_UP
        MUON_EFF_STAT_LOWPT_DN,   //MUON_EFF_STAT_LOWPT_DN
        MUON_EFF_SYS_UP,          //MUON_EFF_SYS_UP
        MUON_EFF_SYS_DN,          //MUON_EFF_SYS_DN
        MUON_EFF_SYS_LOWPT_UP,    //MUON_EFF_SYS_LOWPT_UP
        MUON_EFF_SYS_LOWPT_DN,    //MUON_EFF_SYS_LOWPT_DN
        MUON_EFF_ISO_STAT_UP,     //MUON_ISO_STAT_UP
        MUON_EFF_ISO_STAT_DN,     //MUON_ISO_STAT_DN
        MUON_EFF_ISO_SYS_UP,      //MUON_ISO_SYS_UP
        MUON_EFF_ISO_SYS_DN,      //MUON_ISO_SYS_DN
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

}
