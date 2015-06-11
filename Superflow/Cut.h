#pragma once

#include <string>

#include "Superflow/Superlink.h"

namespace sflow {

    class Cut {
    public:
        Cut() {};
        virtual ~Cut() {};

        std::string name;
        virtual bool operator() (Superlink* sl) = 0;// return true to pass the cut
    };

    class CutName {
    public:
        std::string name;
        std::function<bool(Superlink*)> stored_cut;

        CutName(std::string cut_name);
        ~CutName() {};

        void operator<<(std::function<bool(Superlink*)> cut_);
    };

};
