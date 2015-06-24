#pragma once

#include <string>

#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/SusyNtSys.h"

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
        BKGMETHODUP,        ///< Positive shift due to background estimation method (chargeflip, jetveto, fakeweights,...)
        BKGMETHODDOWN,      ///< Negative shift due to background estimation method (chargeflip, jetveto, fakeweights,...)
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
        ELFRUP,             ///< Positive shift in electron fake rate
        ELFRDOWN,           ///< Negative shift in electron fake rate
        ELREUP,             ///< Positive shift in electron real efficiency
        ELREDOWN,           ///< Negative shift in electron real efficiency
        MUFRUP,             ///< Positive shift in muon fake rate
        MUFRDOWN,           ///< Negative shift in muon fake rate
        MUREUP,             ///< Positive shift in muon real efficiency
        MUREDOWN,           ///< Negative shift in muon real efficiency
        ELFRACUP,           ///< Positive shift in electron fake fraction (see Davide about what this systematic is)
        ELFRACDOWN,         ///< Negative shift in electron fake fraction
        MUFRACUP,           ///< Positive shift in muon fake fraction
        MUFRACDOWN,         ///< Negative shift in muon fake fraction
        ISRUP,              ///< Positive shift in ISR uncertainties (MG5 scale variations up)
        ISRDOWN,            ///< Negative shift in ISR uncertainties (MG5 scale variations down)
        null,               ///< Nominal value
        block               ///< Prevent (re)calculation of the weight
    };

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
