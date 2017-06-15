#ifndef INPUT_OPTIONS_H
#define INPUT_OPTIONS_H

#include <string>

#include "Superflow/Superflow.h"


// struct to hold the options
struct SFOptions {

    SFOptions(int in_argc, char** in_argv) : argc(in_argc), argv(in_argv),  ana_name(""), input(""), n_events_to_process(-1), suffix_name(""), sumw_file_name(""), run_mode(sflow::SuperflowRunMode::nominal), dbg(false) {}

    int argc;
    char** argv;

    std::string ana_name;
    std::string input;
    int n_events_to_process;

    std::string suffix_name;
    std::string sumw_file_name;

    sflow::SuperflowRunMode run_mode;

    bool dbg;


};

void help(std::string ana);
bool read_options(SFOptions& options);


#endif
