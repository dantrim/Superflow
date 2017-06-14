// Supersys.cxx
//
using namespace std;

#include "Superflow/Supersys.h"

namespace sflow {

    //////////////////////////////////////////////////////////////////////////
    Supersys::Supersys()
    {
        reset();
    }

    Supersys::Supersys(SupersysType type_)
    {
        if (type_ == SupersysType::central) {
            name = "CENTRAL";
            tree_name = "NOM";

            type = SupersysType::central;
            event_syst = Susy::NtSys::NOM;
            weight_syst = SupersysWeight::null;
        }
        else {
            reset();
        }
    }

    void Supersys::reset()
    {
        name = "";
        tree_name = "";

        weight_up = 0.0;
        weight_down = 0.0;

        type = SupersysType::null;
        event_syst = Susy::NtSys::NOM;
        weight_syst = SupersysWeight::null;
        weight_syst_up = SupersysWeight::null;
        weight_syst_down = SupersysWeight::null;
    }

    //////////////////////////////////////////////////////////////////////////
    Superweight::Superweight()
    {
        reset();
    }

    double Superweight::product()
    {
        return susynt * lepSf * btag * jvt;
    }

    void Superweight::reset()
    {
        susynt = 1.0;
        pileup = 1.0;
        lepSf = 1.0;
        btag = 1.0;
        jvt = 1.0;
    }

    //////////////////////////////////////////////////////////////////////////
    NewSystematic::NewSystematic(std::string sys_name)
    {
        name = sys_name;
    }

    //////////////////////////////////////////////////////////////////////////
    TreeName::TreeName(std::string tree_name)
    {
        name = tree_name;
    }

    //////////////////////////////////////////////////////////////////////////
    EventSystematic::EventSystematic(Susy::NtSys::SusyNtSys syst_val)
    {
        event_syst_ = syst_val;
    }

    //////////////////////////////////////////////////////////////////////////
    WeightSystematic::WeightSystematic(SupersysWeight syst_val_up, SupersysWeight syst_val_down)
    {
        weight_syst_up = syst_val_up;
        weight_syst_down = syst_val_down;
    }

    //////////////////////////////////////////////////////////////////////////
    BTagSys supersys_to_btagsys(SupersysWeight sys) // BTagSys is from SusyDefs
    {
        BTagSys bsys = BTag_NOM;

        switch (sys) {
          //  case SupersysWeight::BJETUP: bsys = BTag_BJet_UP; break;
          //  case SupersysWeight::BJETDOWN: bsys = BTag_BJet_DN; break;
          //  case SupersysWeight::CJETUP: bsys = BTag_CJet_UP; break;
          //  case SupersysWeight::CJETDOWN: bsys = BTag_CJet_DN; break;
          //  case SupersysWeight::BMISTAGUP: bsys = BTag_LJet_UP; break;
          //  case SupersysWeight::BMISTAGDOWN: bsys = BTag_LJet_DN; break;
            default: break;
        }
        return bsys;
    }
}
