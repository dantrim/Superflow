// Standard
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <dirent.h>

// Root
#include "TROOT.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"

using namespace std;

// prototypes
void list_files(vector<string>& list, string dir);
map<string, string> get_listnames();

// MC_TEXT_DIR : directory containing the filelists
const string MC_TEXT_DIR     = "/gdata/atlas/dantrim/SusyAna/userSusyNt/Superflow/run/filelists/n0206pup/";
// RAW_SAMPLES_DIR : directory where the "raw" ntuples are located
const string RAW_SAMPLES_DIR = "/gdata/atlas/dantrim/SusyAna/histoAna/run2early/n0206pup/Jul7_signal/mc/Raw/";
// NEW_SAMPLES_DIR : output location for the merged ntuple
const string NEW_SAMPLES_DIR = "/gdata/atlas/dantrim/SusyAna/histoAna/run2early/n0206pup/Jul7_signal/mc/Processed/";
// OUT_FILENAME : name of output, merged ntuple
const string OUT_FILENAME    = "mc15_13TeV.root";

struct hft_process {
    string in_file_string;
    string output_tree_name;
};
vector<hft_process> defineTrees(bool do_data) {
    vector<hft_process> out_process;
    hft_process process;
   
    if(!do_data) { 
        process.in_file_string   = "powhegpythia_W_n0206pup.txt";
        process.output_tree_name = "W";
        out_process.push_back(process);
        
        process.in_file_string   = "powhegpythia_Z_n0206pup.txt";
        process.output_tree_name = "Z";
        out_process.push_back(process);
        
        process.in_file_string   = "powheg_ttbar_n0206pup.txt";
        process.output_tree_name = "ttbar";
        out_process.push_back(process);
    }
    else if(do_data) {
        process.in_file_string   = "dummy_n0206pup.txt";
        process.output_tree_name = "Data";
        out_process.push_back(process);
    }

    return out_process;
}
    
enum hft_sys_type {
    k_sys_object,
    k_sys_one_sided_object,
    k_sys_weight,
    k_sys_adhoc,
    k_N_sys_type
};
struct hft_systematic {
    hft_sys_type sys_type;
    string basename;
    string up_name;
    string down_name;
    string preface;
};

vector<hft_systematic> defineSystematics()
{
    vector<hft_systematic> hft_syst;
    hft_systematic syst;

    /////////////////////
    // Nominal tree
    /////////////////////
    syst.sys_type = k_sys_one_sided_object;
    syst.basename = "CENTRAL";
    syst.up_name = "";
    syst.down_name = "";
    syst.preface = "";
    hft_syst.push_back(syst);


    return hft_syst;
}


string trim(string s)
{
    string val = s.erase(s.find_last_not_of(" \n\r\t") + 1);
    val.erase(0, s.find_first_not_of(" \n\r\t"));
    return val;
}

int main(int argc, char** argv)
{
    bool do_data = false;
    for(int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--data") == 0) do_data = true;
    }

    vector<string> sample_dir_files;
    list_files(sample_dir_files, RAW_SAMPLES_DIR); 

    vector<hft_process> hft_trees = defineTrees(do_data);

    vector<hft_systematic> hft_syst = defineSystematics();
    vector<string> syst_strings;

    string output_filename = string(NEW_SAMPLES_DIR) + OUT_FILENAME;
    TFile* output_file = new TFile(output_filename.data(), "RECREATE");
    output_file->Close();
    delete output_file;
    
    for (int i = 0; i < hft_syst.size(); i++) {
        if( hft_syst[i].sys_type == k_sys_object ) {
            syst_strings.push_back(hft_syst[i].basename + hft_syst[i].up_name);
            syst_strings.push_back(hft_syst[i].basename + hft_syst[i].down_name);
        }
        else if (hft_syst[i].sys_type == k_sys_one_sided_object) {
            syst_strings.push_back(hft_syst[i].basename);
        }
    }

    for (int s = 0; s < syst_strings.size(); s++) {
        for (int g = 0; g < hft_trees.size(); g++) {
            if(syst_strings[s].compare("CENTRAL") != 0 && hft_trees[g].output_tree_name == 0) continue;
            
            vector<string> sample_list;
            stringstream in_filename;
            in_filename << MC_TEXT_DIR << hft_trees[g].in_file_string;
            ifstream mc_in_(in_filename.str().data());

            if (!mc_in_.is_open()) {
                cout << "File (" << in_filename.str() << ") not open." << endl;
                cout << " > Exitting." << endl;
                return 0;
            }
            else {
                cout << endl;
                cout << "Opened " << in_filename.str() << endl;
            }

            while (!mc_in_.eof()) {
                string line_ = "";
                getline(mc_in_, line_);
                line_ = trim(line_);
                
                string find_this = "";
                int jump1 = 0;
                int jump2 = 0;
                if(do_data) { find_this = "data15_13TeV"; jump1 = 15; jump2 = 6; }
                else { find_this = "mc15_13TeV"; jump1 = 11; jump2 = 6; }
                
                size_t found = line_.find(find_this);
                if(found != string::npos && line_ != "") {
                    string sample = line_.substr(found + jump1, jump2);
                    sample_list.push_back(sample);
                    cout << "sample : " << sample << endl;
                }
            }
            
            string tree_name = hft_trees[g].output_tree_name + "_" + syst_strings[s];
            cout << "output tree : " << tree_name << endl;
            
            output_file = new TFile(output_filename.data(), "UPDATE");
            output_file->cd();

            TChain* merge_chain = new TChain(tree_name.data());
            TTree::SetMaxTreeSize(137438953472LL);
            
            int sum_trees = 0;
            
            for (int i = 0; i < sample_list.size(); i++) {
                string sample_treename = "";
                if(do_data) { sample_treename = "id_physics_Main"; }
                else {
                    sample_treename = string("id_") + sample_list[i];
                }
                string sample_filename = "";
                if(do_data) { sample_filename = syst_strings[s] + "_physics_Main_" + sample_list[i] + ".root"; }
                else {
                    sample_filename = syst_strings[s] + "_" + sample_list[i] + ".root";
                }
                bool flag_found_root = false;

                auto file_nm = find(sample_dir_files.begin(), sample_dir_files.end(), sample_filename);
                if(file_nm != sample_dir_files.end()) {
                    flag_found_root = true;
                }
                else {
                    cout << sample_filename << "  not found      !#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#" << endl;
                }
                if(!flag_found_root) continue;

                string root_file_name = string(RAW_SAMPLES_DIR) + sample_filename;

                TFile* in_file = new TFile(root_file_name.data());
                TTree* in_tree = static_cast<TTree*>(in_file->Get(sample_treename.data()));

                int n_entries = 0;
                if (in_tree != nullptr) {
                    n_entries = in_tree->GetEntries();
                    cout << sample_filename << "   " << to_string(n_entries) << endl;
                    sum_trees += n_entries;
                }
                if(n_entries == 0) continue;
                merge_chain->AddFile(root_file_name.data(), 0, sample_treename.data());
                delete in_tree;
                delete in_file;
            }

            output_file->cd();
            merge_chain->Merge(output_file, 0, "fast");

            delete merge_chain;
            cout << "\nsum : " << sum_trees << endl;
            mc_in_.close();


        }
    }
    return 0;
}

void list_files(vector<string>& list, string dir)
{
    DIR *pDIR;
    struct dirent *entry;
    if((pDIR = opendir(dir.data()))) {
        while (( entry = readdir(pDIR) )) {
            if(strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0)
                list.push_back(entry->d_name);
        }
        closedir(pDIR);
    }
}
                
                
