#include <cstdlib>
#include <string>
#include <cassert>
#include <cmath> // isnan
#include <cfloat> // FLT_MAX, FLT_MIN
#include <iomanip> // setw, setprecision
#include <fstream> 
#include <sstream>  // std::ostringstream
#include <dirent.h> // UNIX

#include "TChain.h"
#include "TCanvas.h"
#include "Cintex/Cintex.h"

using namespace std;

void listFiles(vector<string>& list_, string dir_);

#define MC_TEXT_DIR "/gdata/atlas/suneetu/Documents/LFV_Higgs2014/generation/filelists/"
#define MC_SAMPLES_DIR "/gdata/atlas/ucintprod/SusyNt/mc12_n0154b/"
#define OUT_PREFIX "file_list_"

const int n_files = 5;
string files[] = { "Higgs", "top_MCNLO", "WW_Sherpa", "WZ_ZZ_Sherpa", "Zjets_AlpgenPythia" };
string ext = ".txt";

string ren_files[] = { "Higgs", "ttbar", "WW", "ZV", "ZPlusJets" };

int main(int argc, char** argv)
{
    ROOT::Cintex::Cintex::Enable();

    vector<string> sample_numbers;
    vector<string> files_in_mc;

    listFiles(files_in_mc, MC_SAMPLES_DIR);

    for (int f_ = 0; f_ < n_files; f_++) {
        vector<string> new_list;

        stringstream in_filename;
        in_filename << MC_TEXT_DIR << files[f_] << ext;

        ifstream mc_in_(in_filename.str().data());

        if (!mc_in_.is_open()) {
            cout << "File not open. Exiting." << endl;
            return 0;
        }
        else {
            cout << endl << "Opened " << in_filename.str() << endl;
        }

        while (!mc_in_.eof()) {
            string line_ = "";
            getline(mc_in_, line_);

            size_t found_ = line_.find_first_of('.');

            if (found_ != string::npos && line_ != "") {
                int mc_found = 0;
                string sample_ = line_.substr(found_ + 1, 6);

                vector<int> near_index;
                vector<string> near_list;

                for (uint i = 0; i < files_in_mc.size(); i++) {
                    string line_ = files_in_mc[i];
                    size_t found = line_.find(sample_);

                    if (found != string::npos) {
                        mc_found++;

                        cout << sample_ << " -> " << line_;
                        if (mc_found > 1) cout << " x" << mc_found << " #.#.#.#";
                        cout << endl;

                        near_list.push_back(line_);
                        size_t find_e_ = line_.find("SusyNt.e");// 8 chars
                        if (find_e_ != string::npos) {
                            near_index.push_back(atoi(line_.substr(find_e_ + 8, 4).data()));
                        }
                        else cout << "ERROR: can't find e-number" << endl;
                    }
                }

                if (near_index.size() != 0) {
                    int large = 0;
                    int one_index = 0;

                    for (uint j = 0; j < near_index.size(); j++) {
                        if (near_index[j] > large) {
                            large = near_index[j];
                            one_index = j;
                        }
                    }

                    new_list.push_back(near_list[one_index]);
                }
                else {
                    cout << "ERROR: can't find process " << sample_ << endl;
                }
            }
        }
        mc_in_.close();

        stringstream out_filename;
        out_filename << MC_TEXT_DIR << OUT_PREFIX << ren_files[f_] << ext;

        ofstream out_file;
        out_file.open(out_filename.str().data(), ofstream::trunc);

        if (!out_file.is_open()) {
            cout << "Did not create output file!" << endl;
            return 0;
        }

        for (uint i = 0; i < new_list.size(); i++) {
            out_file << MC_SAMPLES_DIR << new_list[i] << "/" << "\r\n";
            //out_file << new_list[i] << "/" << "\r\n";
        }
        out_file.close();
    }

    return 0;
}

void listFiles(vector<string>& list_, string dir_)
{
    DIR *pDIR;
    struct dirent *entry;

    if ((pDIR = opendir(dir_.data()))) {
        while ((entry = readdir(pDIR))) {
            if (strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0)
                list_.push_back(entry->d_name);
        }
        closedir(pDIR);
    }
}