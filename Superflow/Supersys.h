#pragma once

#include <string>

#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/SusyNtSys.h"

using namespace std;

namespace sflow {

    const string weight_prefix = "syst_";
    const string weight_suffix_up = "UP";
    const string weight_suffix_down = "DOWN";

    enum class SupersysType { central, event, weight, null };

    enum class SupersysState {
        closed,
        open
    };

    enum class SupersysWeight {
        EL_EFF_ID_UP,
        EL_EFF_ID_DOWN,
        EL_EFF_RECO_UP,
        EL_EFF_RECO_DOWN,
        MUON_EFF_STAT_UP,
        MUON_EFF_STAT_DOWN,
        MUON_EFF_SYST_UP,
        MUON_EFF_SYST_DOWN,
        PILEUP_UP,
        PILEUP_DOWN,
        null,
        block
    };
    

/*
 // run 1 sys
    enum class SupersysWeight {
        ETRIGREWUP,         ///< Positive shift in electron trigger weights
        ETRIGREWDOWN,       ///< Negative shift in electron trigger weights
        MTRIGREWUP,         ///< Positive shift in muon trigger weights
        MTRIGREWDOWN,       ///< Negative shift in muon trigger weights
        BJETUP,             ///< Positive shift in btag scale factor
        BJETDOWN,           ///< Negative shift in btag scale factor
        CJETUP,             ///< Positive shift in ctag scale factor
        CJETDOWN,           ///< Negative shift in ctag scale factor
        BMISTAGUP,          ///< Positive shift in ltag scale factor
        BMISTAGDOWN,        ///< Negative shift in ltag scale factor
        ESFUP,              ///< Positive shift in electron efficiency
        ESFDOWN,            ///< Negative shift in electron efficiency
        MEFFUP,             ///< Positive shift in muon efficiency
        MEFFDOWN,           ///< Negative shift in muon efficiency
        PILEUPUP,           ///< Positive shift in pileup
        PILEUPDOWN,         ///< Negative shift in pileup
        XSUP,               ///< Positive shift in xsec
        XSDOWN,             ///< Negative shift in xsec
        null,               ///< Nominal value
        block               ///< Prevent (re)calculation of the weight
    };
*/
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

    class Superweight {
    public:
        Superweight();

        double product();
        void reset();

        double susynt; ///< from SusyNtTools::getEventWeight: includes gen, pileup, xs, lumi, sumw
        double gen;
        double pileup;
        double norm; ///< breakdown of the above; norm is xs*lumi/sumw
        double lepSf;
        double btag;
        double trigger;
        double qflip;
        double fake;
        double isr;
    };

    class NewSystematic {
    public:
        NewSystematic(std::string sys_name);

        std::string name;
    };


    class TreeName {
    public:
        TreeName(std::string var_name);

        std::string name;
    };

    class EventSystematic {
    public:
        EventSystematic(Susy::NtSys::SusyNtSys syst_val);

        Susy::NtSys::SusyNtSys event_syst_;
    };

    class WeightSystematic {
    public:
        WeightSystematic(SupersysWeight syst_val_up, SupersysWeight syst_val_down);

        SupersysWeight weight_syst_up;
        SupersysWeight weight_syst_down;
    };

    class SaveSystematic {
    public:
        SaveSystematic() {};
    };

    BTagSys supersys_to_btagsys(SupersysWeight sys);

};
