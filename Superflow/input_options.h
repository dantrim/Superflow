#ifndef INPUT_OPTIONS_H
#define INPUT_OPTIONS_H

#include <string>

#include "Superflow/Superflow.h"
//enum class sflow::SuperflowRunMode;

//read_options(std::std::string, int, char**, std::std::string&, int&, std::std::string&, sflow::SuperflowRunMode&, std::std::string&, bool)
void help(std::string ana);
bool read_options(std::string ana, int argc, char* argv[], std::string& input, int& num_events_, std::string& suffix_name_,
            sflow::SuperflowRunMode& run_mode_, std::string& sumw_file, bool dbg);



#endif
