#include "Superflow/input_options.h"

#include <iostream>
#include <fstream>
using namespace std;



void help(string ana)
{
    cout << "---------------------------------------------------------------" << endl;
    cout << " " << ana << endl;
    cout << endl;
    cout << " Options: " << endl;
    cout << "    -i             : input (ROOT file, *.txt file, or dir) [REQUIRED]" << endl;
    cout << "    -n             : number of events to process [default: -1, all]" << endl;
    cout << "    -d             : set debug true [default: false]" << endl;
    cout << "    --suffix       : provide a string to be added as a suffix to all output files [default: \"\"]" << endl;
    cout << "    --sumw         : provide a SUMW file to provide to SusyNtuple::MCWeighter for eventweight calculation" << endl;
    cout << "    -h|--help      : print this help message" << endl;
    cout << " - - - - - - - " << endl;
    cout << " Run Modes (one is REQUIRED) " << endl;
    cout << "    -c             : CENTRAL (nominal) running only" << endl; 
    cout << "    -w             : CENTRAL + WEIGHT SYS running" << endl;
    cout << "    -a             : CENTRAL + WEIGHT + KINEMATIC SYS running" << endl;
    cout << endl;
    cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
    cout << "  Usage: " << endl;
    cout << endl;
    cout << "    " << ana << " -i susyNt.root -c --sumw ./sumw_file.txt" << endl;
    cout << "---------------------------------------------------------------" << endl;


}

bool read_options(string ana_name, int argc, char* argv[], string& input, int& num_events, string& suffix_name,
            sflow::SuperflowRunMode& run_mode, string& sumw_file, bool dbg)
{
    bool nominal = false;
    bool nominal_and_weight_sys = false;
    bool all_sys = false;

    for(int i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-n") == 0) {
            num_events = atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-i") == 0) {
            input = argv[++i];
        }
        else if(strcmp(argv[i], "-c") == 0) {
            nominal = true;
        }
        else if(strcmp(argv[i], "-w") == 0) {
            nominal_and_weight_sys = true;
        }
        else if(strcmp(argv[i], "-a") == 0) {
            all_sys = true;
        }
        else if(strcmp(argv[i], "-d") == 0) {
            dbg = true;
        }
        else if(strcmp(argv[i], "--suffix") == 0) {
            suffix_name = argv[++i];
        }
        else if(strcmp(argv[i], "--sumw") == 0) {
            sumw_file = argv[++i];
        }
        else if( (strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0) ) {
            help(ana_name);
            return false;
        }
        else {
            cout << ana_name << "::input_options    Unknown command line argument '" << argv[i] << "', exiting" << endl;
            help(ana_name);
            return false;
        }
    }


    if(nominal) {
        run_mode = sflow::SuperflowRunMode::nominal;
        cout << ana_name << "    Run mode: nominal" << endl;
    }
    else if(nominal_and_weight_sys) {
        run_mode = sflow::SuperflowRunMode::nominal_and_weight_syst;
        cout << ana_name << "    Run mode: nominal_and_weight_syst" << endl;
    }
    else if(all_sys) {
        run_mode = sflow::SuperflowRunMode::all_syst;
        cout << ana_name << "    Run mode: all_syst" << endl;
    }

    // check input sumw file if it exists
    bool sumw_file_exists = true;
    if(sumw_file != "") {
        sumw_file_exists = std::ifstream(sumw_file).good();
    }
    if(!sumw_file_exists) {
        cout << ana_name << "    ERROR Provided SUMW file '" << sumw_file << "' cannot be found!" << endl;
        return false;
    }

    // check that the user provided an input
    if(input=="") {
        cout << ana_name << "    ERROR You must provide an input file, dir, or filelist!" << endl;
        help(ana_name);
        return false;
    }

    return true;
}
