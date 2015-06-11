#pragma once

#include <string>

namespace sflow {

    class Supervar {
    public:
        Supervar() {}; 
        ~Supervar() {};
        // this class is used for type identification 
        // in operator overloading
    };

    class var_float {
        // this class is used for type identification 
        // in operator overloading
    };

    class var_double {
        // this class is used for type identification 
        // in operator overloading
    };

    class var_int {
        // this class is used for type identification 
        // in operator overloading
    };

    class var_bool {
        // this class is used for type identification 
        // in operator overloading
    };

    class var_void {
        // this class is used for type identification 
        // in operator overloading
    };

    enum class SupervarState {
        closed,
        open
    };

    enum class SupervarType {
        sv_float,
        sv_double,
        sv_int,
        sv_bool,
        sv_void
    };

    class NewVar {
    public:
        std::string name;

        NewVar(std::string var_name);
        ~NewVar() {};
    };

    class HFTname {
    public:
        std::string name;

        HFTname(std::string hft_name);
    };

    class SaveVar {
    public:
        SaveVar() {};
    };

};

