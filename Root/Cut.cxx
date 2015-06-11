// Cut.cxx
//
#include "Superflow/Cut.h"

namespace sflow {

    CutName::CutName(std::string cut_name)
    {
        name = cut_name;
    }

    void CutName::operator<<(std::function<bool(Superlink*)> cut_)
    {
        stored_cut = cut_;
    }
}