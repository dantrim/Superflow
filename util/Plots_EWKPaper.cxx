// Copyright 2014, The ATLAS Collaboration
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <vector>
#include <map>

#include "TROOT.h"
//#include "TDirectory.h"
#include "TMath.h"

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCut.h"
#include "TEventList.h"

#include "THStack.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TColor.h"
#include "TStyle.h"

#include "Superflow/StringTools.h"

#include "SusyNtuple/TGuiUtils.h"

//#include "Cintex/Cintex.h"

using namespace std;
using namespace sflow;

#define DEBUG_NTC true


//
//MC Backgrounds + DATA
//
#define DIR_DATA "/scratch/dantrim/n0209/"
#define DATA_FILE "data15_13TeV.root"

//#define DIR_MC "/gdata/atlas/dantrim/SusyAna/histoAna/run2early/n0206pup/Jun30_signal/mc/Processed/"
#define DIR_MC "/scratch/dantrim/n0209/"
#define MC_FILE "mc15_13TeV.root"

#define TT_FILE "ttbar_wjets_mc15_13TeV.root"

#define DIR_1 "/gdata/atlas/dantrim/SusyAna/n0209val/plots/"


const string m_central = "CENTRAL";
//const string m_central = "";

const int n_hists = 1;
string hists_entries[] = {};

const double Pi = 3.141592653589793238462;

struct plot {
    string root_get_member;
    string root_get_member_name;
    string root_get_factor;
    string signal_region;
    string x_label;
    string y_label;
    string request_cuts;
    string reweight;
    string blind_data;
    bool is_Log;
    bool leg_is_left;  // for putting legend on left hand side (MOVE ANNOYING LABELS TO RHS!!)
    double x_bin_width;
    double x_range_max;
    double x_range_min;
    double y_range_max;
    double y_range_min;
    string long_name;
    int run_group;
    double x_red_arrow;
};

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
    string adhoc_tcut;
    string adhoc_tcut2;
    int binding_group;// inside a binding group, *trees* add linearly
};

enum processed_tree_type {
    k_tree_mc,
    k_tree_signal,
    k_tree_fakes,
    k_tree_data,
    k_N_tree_type
};

struct processed_tree {
    processed_tree_type tree_type;
    string treename;
    string displayname;
    string blinding_cuts;
    string special_cuts;
    string file;
    int color;
    vector<int> systematic_groups;
    TEventList central_list;
    map<string, TEventList> all_lists;
    bool all_lists_exist;
};

void convertErrorsToPoisson(TH1* histo, TGraphAsymmErrors* graph)
{
    // Needed variables
    double value = 0;
    double error_poisson_up = 0;
    double error_poisson_down = 0;
    double alpha = 0.158655, beta = 0.158655; // 68%

    // loop over bins and overwrite values
    for (int i = 1; i < histo->GetNbinsX() + 1; i++) {
        value = histo->GetBinContent(i);
        if (value != 0) {
            error_poisson_up = 0.5*TMath::ChisquareQuantile(1 - beta, 2 * (value + 1)) - value;
            error_poisson_down = value - 0.5*TMath::ChisquareQuantile(alpha, 2 * value);
            graph->SetPoint(i - 1, histo->GetBinCenter(i), value);
            graph->SetPointError(i - 1, 0., 0., error_poisson_down, error_poisson_up);
        }
        else {
            graph->SetPoint(i - 1, histo->GetBinCenter(i), -10);
            graph->SetPointError(i - 1, 0., 0., 0., 0.);
        }
    }
}

void myText(Double_t x, Double_t y, Color_t color, const char *text)
{
  //  Double_t tsize = 0.055;
    Double_t tsize = 0.055;
    TLatex l; //l.SetTextAlign(12); 
    l.SetTextSize(tsize);
    l.SetTextFont(42);
    l.SetNDC();
    l.SetTextColor(color);
    l.DrawLatex(x, y, text);
}

vector<processed_tree> defineTrees(string region)
{
    cout << "defineTrees()" << endl;
    vector<processed_tree> hft_trees;

    processed_tree tree_;

    vector<int> grouplist;
    grouplist.clear();

/*    size_t is_vr = region.find("vr");
    size_t is_sr0  = region.find("sr0");
    size_t is_sr1a = region.find("sr1a");
    size_t is_sr1b = region.find("sr1b");
    size_t is_sr1c  = region.find("sr1c");

    if(is_sr1a!=string::npos || is_sr1c!=string::npos){
        tree_.file = SIG_FILE;
        tree_.tree_type = k_tree_signal;
        tree_.blinding_cuts = "";
        tree_.special_cuts  = "";
        tree_.treename = "SMCwslep8TeV_100.0_80.0";
        tree_.displayname = "(m_{#tilde{#chi}_{1}^{ #pm}}, m_{#tilde{#chi}_{1}^{0}}) = (100, 80) GeV";
        tree_.color = kBlack;
        grouplist.clear();
        tree_.systematic_groups = grouplist;
        hft_trees.push_back(tree_);
        
        tree_.treename = "SMCwslep8TeV_100.0_35.0";
        tree_.displayname = "(m_{#tilde{#chi}_{1}^{ #pm}}, m_{#tilde{#chi}_{1}^{0}}) = (100, 35) GeV";
        tree_.color = kGreen-4;
        grouplist.clear();
        tree_.systematic_groups = grouplist;
        hft_trees.push_back(tree_);
    }
    else if(is_sr1b!=string::npos){
        tree_.file = SIG_FILE;
        tree_.tree_type = k_tree_signal;
        tree_.blinding_cuts = "";
        tree_.special_cuts  = "";
        tree_.treename = "SMCwslep8TeV_130.0_30.0";
        tree_.displayname = "(130, 30)";
        tree_.color = kBlack;
        grouplist.clear();
        tree_.systematic_groups = grouplist;
        hft_trees.push_back(tree_);
        
        tree_.treename = "SMCwslep8TeV_112.5_12.5";
        tree_.displayname = "(112.5, 12.5)";
        tree_.color = kGreen-4;
        grouplist.clear();
        tree_.systematic_groups = grouplist;
        hft_trees.push_back(tree_);
    }
*/


    TEventList blank;                 // COMMENTED THESE OUT WHEN ADDED NEW TOP
    tree_.central_list = blank;
    tree_.all_lists_exist = false;

    tree_.file = MC_FILE;

    tree_.tree_type = k_tree_mc;
    tree_.treename = "Z";
    tree_.displayname = "Z";
    tree_.color = TColor::GetColor("#82DE68");
    tree_.special_cuts = "0.94"; // 0.32
    grouplist.clear();
    grouplist.push_back(0);
    tree_.systematic_groups = grouplist;
    hft_trees.push_back(tree_);

    // Top
    tree_.treename = "ttbar";
    tree_.displayname = "ttbar";
    tree_.color = TColor::GetColor("#E67067");
    tree_.special_cuts = "0.94";
    grouplist.clear();
    grouplist.push_back(0);
    tree_.systematic_groups = grouplist;
    hft_trees.push_back(tree_);

    // WW
    tree_.treename = "W";
    tree_.displayname = "W";
    tree_.color = TColor::GetColor("#5E9AD6");
    tree_.special_cuts = "0.94";
    grouplist.clear();
    grouplist.push_back(0);
//    grouplist.push_back(4);// XS
    tree_.systematic_groups = grouplist;
    hft_trees.push_back(tree_);

    // DATA trees
    tree_.file = DATA_FILE;
    tree_.treename = "Data";
    tree_.displayname = "Data";
    tree_.tree_type = k_tree_data;
    tree_.blinding_cuts = "";//   (mcollCorr > 150 || mcollCorr < 100) 
    tree_.special_cuts = "";
    tree_.color = kBlack;
    grouplist.clear();// no systematics on data
    tree_.systematic_groups = grouplist;
    hft_trees.push_back(tree_);


    return hft_trees;
}

vector<hft_systematic> defineSystematics(string region)
{
    cout << "defineSystematics" << endl;
    vector<hft_systematic> hft_syst;
/*
    hft_systematic syst;

    syst.adhoc_tcut2 = "";

    size_t is_sr1a = region.find("sr1a");
    size_t is_sr1b = region.find("sr1b");
    size_t is_sr1c = region.find("sr1c");
    size_t is_crWWb = region.find("crWWb");

   // THEORY/GEN UNCERTAINTIES FROM SARAH
    // (SEE NOTE)
    syst.sys_type = k_sys_adhoc;
    syst.preface  = "syst_";
   
    if(is_sr1a!=string::npos || is_sr1b!=string::npos || is_sr1c!=string::npos || is_crWWb!=string::npos){ 
        // WW
        syst.binding_group = 2;
        if     (is_sr1a!=string::npos) { syst.basename = "SR1a_WW_Theory"; syst.adhoc_tcut = "1.07"; syst.adhoc_tcut2 = "0.93"; }
        else if(is_sr1b!=string::npos) { syst.basename = "SR1b_WW_Theory"; syst.adhoc_tcut = "1.11"; syst.adhoc_tcut2 = "0.89"; }
        else if(is_sr1c!=string::npos) { syst.basename = "SR1c_WW_Theory"; syst.adhoc_tcut = "1.08"; syst.adhoc_tcut2 = "0.92"; }
       // syst.adhoc_tcut  = "1.07";
       // syst.adhoc_tcut2 = "0.93";
        hft_syst.push_back(syst);

        // ZV
        syst.binding_group = 3;
        if     (is_sr1a!=string::npos) { syst.basename = "SR1a_ZV_Theory"; syst.adhoc_tcut = "1.28"; syst.adhoc_tcut2 = "0.72"; }
        else if(is_sr1b!=string::npos) { syst.basename = "SR1b_ZV_Theory"; syst.adhoc_tcut = "1.45"; syst.adhoc_tcut2 = "0.55"; }
        else if(is_sr1c!=string::npos) { syst.basename = "SR1c_ZV_Theory"; syst.adhoc_tcut = "1.46"; syst.adhoc_tcut2 = "0.54"; }
       // syst.basename = "SR1a_ZV_Theory";
       // syst.adhoc_tcut  = "1.281";
       // syst.adhoc_tcut2 = "0.719";
        hft_syst.push_back(syst);

        // Top
        syst.binding_group = 4;
        if     (is_sr1a!=string::npos)  { syst.basename = "SR1a_Top_PS_Theory";  syst.adhoc_tcut = "1.23"; syst.adhoc_tcut2 = "0.77"; }
        else if(is_sr1b!=string::npos)  { syst.basename = "SR1b_Top_PS_Theory";  syst.adhoc_tcut = "1.25"; syst.adhoc_tcut2 = "0.75"; }
        else if(is_sr1c!=string::npos)  { syst.basename = "SR1c_Top_PS_Theory";  syst.adhoc_tcut = "1.27"; syst.adhoc_tcut2 = "0.73"; }
        else if(is_crWWb!=string::npos) { syst.basename = "crWWb_Top_PS_Theory"; syst.adhoc_tcut = "1.29"; syst.adhoc_tcut2 = "0.71"; }
        //syst.basename = "SR1a_Top_Theory";
        //syst.adhoc_tcut  = "1.259";
        //syst.adhoc_tcut2 = "0.741";
        hft_syst.push_back(syst);

        // Top
        syst.binding_group = 5;
        if     (is_sr1a!=string::npos)  { syst.basename = "SR1a_Top_OTHER_Theory";  syst.adhoc_tcut = "1.12"; syst.adhoc_tcut2 = "0.88"; }
        else if(is_sr1b!=string::npos)  { syst.basename = "SR1b_Top_OTHER_Theory";  syst.adhoc_tcut = "1.16"; syst.adhoc_tcut2 = "0.84"; }
        else if(is_sr1c!=string::npos)  { syst.basename = "SR1c_Top_OTHER_Theory";  syst.adhoc_tcut = "1.12"; syst.adhoc_tcut2 = "0.88"; }
        else if(is_crWWb!=string::npos) { syst.basename = "crWWb_Top_OTHER_Theory"; syst.adhoc_tcut = "1.13"; syst.adhoc_tcut2 = "0.87"; }
        hft_syst.push_back(syst);
    }


     // OBJECT SYSTEMATICS
     // these apply to all
     syst.sys_type = k_sys_object;
     syst.up_name = "UP";
     syst.down_name = "DOWN";
     syst.preface = "";
     syst.adhoc_tcut = "";
     syst.binding_group = 0;
     
     // EESZ
     syst.basename = "EESZ";
     hft_syst.push_back(syst);
     
     // EER
     syst.basename = "EER";
     hft_syst.push_back(syst);
     
     // EESLOW
     syst.basename = "EESLOW";
     hft_syst.push_back(syst);
     
     // EESMAT
     syst.basename = "EESMAT";
     hft_syst.push_back(syst);
     
     // EESPS
     syst.basename = "EESPS";
     hft_syst.push_back(syst);
     
     // ID
     syst.basename = "ID";
     hft_syst.push_back(syst);
     
     // JES
     syst.basename = "JES";
     hft_syst.push_back(syst);
     
     // MS
     syst.basename = "MS";
     hft_syst.push_back(syst);
     
     // SCALEST
     syst.basename = "SCALEST";
     hft_syst.push_back(syst);
     
     // TES
     syst.basename = "TES";
     hft_syst.push_back(syst);
     
     // ONE-SIDED SYSTEMATICS
     // these apply to all
     syst.sys_type = k_sys_one_sided_object;
     syst.up_name = "";
     syst.down_name = "";
     syst.preface = "";
     syst.adhoc_tcut = "";
     syst.binding_group = 0;
     
     // JER
     syst.basename = "JER";
     hft_syst.push_back(syst);
     
     // RESOST
     syst.basename = "RESOST";
     hft_syst.push_back(syst);
     
     // WEIGHT SYSTEMATICS
     // these apply to all
     syst.sys_type = k_sys_weight;
     syst.up_name = "UP";
     syst.down_name = "DOWN";
     syst.preface = "syst_";
     syst.adhoc_tcut = "";
     syst.binding_group = 0;
     
     // BJET
     syst.basename = "BJET";
     hft_syst.push_back(syst);
     
     // CJET
     syst.basename = "CJET";
     hft_syst.push_back(syst);
     
     // BMISTAG
     syst.basename = "BMISTAG";
     hft_syst.push_back(syst);
     
     // ETRIGREW
     syst.basename = "ETRIGREW";
     hft_syst.push_back(syst);
     
     // MTRIGREW
     syst.basename = "MTRIGREW";
     hft_syst.push_back(syst);
     
     // ESF
     syst.basename = "ESF";
     hft_syst.push_back(syst);
     
     // MEFF
     syst.basename = "MEFF";
     hft_syst.push_back(syst);
     
     // PILEUP
     syst.basename = "PILEUP";
     hft_syst.push_back(syst);
     
     // MEFF
     syst.basename = "MEFF";
     hft_syst.push_back(syst);

     syst.binding_group = 1;    // move to fakes
     // ELFR
     syst.basename = "ELFR";
     hft_syst.push_back(syst);

     // ELRE
     syst.basename = "ELRE";
     hft_syst.push_back(syst);

     // MUFR
     syst.basename = "MUFR";
     hft_syst.push_back(syst);

     // MURE
     syst.basename = "MURE";
     hft_syst.push_back(syst);

    // MUFRAC
      syst.basename = "ELFRAC";
      hft_syst.push_back(syst);

    // ELFRAC
      syst.basename = "MUFRAC";
      hft_syst.push_back(syst);
*/     
/*     // SPLIT XS
     // these apply to all XS
     syst.sys_type = k_sys_weight;
     syst.up_name = "UP";
     syst.down_name = "DOWN";
     syst.preface = "syst_";
     syst.adhoc_tcut = "";
     
     // XS
     syst.basename = "XS";
     syst.binding_group = 3;
     hft_syst.push_back(syst);
     
     syst.basename = "XS";
     syst.binding_group = 4;
     hft_syst.push_back(syst);
     
     syst.basename = "XS";
     syst.binding_group = 5;
     hft_syst.push_back(syst);
     
     // no higgs!
     
     syst.basename = "XS";
     syst.binding_group = 7;
     hft_syst.push_back(syst);
     
     // XS on Z+jets
     syst.sys_type = k_sys_adhoc;
     syst.basename = "XS";
     syst.binding_group = 6;
     syst.preface = "syst_";
     syst.adhoc_tcut = "1.0014";
     syst.adhoc_tcut2 = "0.9986"; // note
     hft_syst.push_back(syst);

    
    // SPLIT GEN
    // these apply to all GEN
    syst.sys_type = k_sys_weight;
    syst.up_name = "UP";
    syst.down_name = "DOWN";
    syst.preface = "syst_";
    syst.adhoc_tcut = "";

    //GEN
    syst.basename = "GEN";
    syst.binding_group = 3;
    hft_syst.push_back(syst);

    syst.basename = "GEN";
    syst.binding_group = 4;
    hft_syst.push_back(syst);

    syst.basename = "GEN";
    syst.binding_group = 5;
    hft_syst.push_back(syst);

    syst.basename = "GEN";
    syst.binding_group = 6;
    hft_syst.push_back(syst);

    syst.basename = "GEN";
    syst.binding_group = 7;
    hft_syst.push_back(syst);

    //FAKE SYSTEMATICS (object)
    //these apply to all
    syst.sys_type = k_sys_object;
    syst.up_name = "UP";
    syst.down_name = "DOWN";
    syst.preface = "";
    syst.adhoc_tcut = "";
    syst.binding_group = 1;

    // ELFR
    syst.basename = "ELFR";
    hft_syst.push_back(syst);

    // MUFR
    syst.basename = "MUFR";
    hft_syst.push_back(syst);

    // ELRE
    syst.basename = "ELRE";
    hft_syst.push_back(syst);

    // MURE
    syst.basename = "MURE";
    hft_syst.push_back(syst);

    // ELFRAC
    syst.basename = "ELFRAC";
    hft_syst.push_back(syst);

    // MUFRAC
    syst.basename = "MUFRAC";
    hft_syst.push_back(syst);

    */

    // DONE //
    return hft_syst;
}

vector<plot> definePlots()
{
    cout << "definePlots()" << endl;
    vector<plot> vect_plot;
    plot plot_;

    plot_.run_group = 0;
    
    plot_.request_cuts = "";
    plot_.blind_data = "";
    plot_.is_Log = false;
        /////////////
        // z_like
        /////////////
        plot_.root_get_member   = "mll";
        plot_.root_get_factor   = "";
        plot_.x_label           = "m_{ee} [GeV]";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_mll";
            plot_.y_label       = "Events / 2 GeV";
            plot_.x_bin_width   = 2;
            plot_.x_range_min   = 70;
            plot_.x_range_max   = 110;
            plot_.y_range_min   = 0.0;
            plot_.y_range_max   = 1000;
            vect_plot.push_back(plot_);

        plot_.is_Log = true;
        plot_.root_get_member   = "nJets";
        plot_.root_get_factor   = "";
        plot_.x_label           = "Number of jets";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_nJets";
            plot_.y_label       = "Events";
            plot_.x_bin_width   = 1;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 8;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);


        plot_.root_get_member   = "met";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} [GeV]";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_met";
            plot_.y_label       = "Events / 5 GeV";
            plot_.x_bin_width   = 5;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 180;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "metPhi";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} #phi";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_metPhi";
            plot_.y_label       = "Events";
            plot_.x_bin_width   = 0.2;
            plot_.x_range_min   = -3.2;
            plot_.x_range_max   = 3.2;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);


        plot_.root_get_member   = "l_pt[0]";
        plot_.root_get_factor   = "";
        plot_.x_label           = "Lead lepton p_{T} [GeV]";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_lpt0";
            plot_.y_label       = "Events / 5 GeV";
            plot_.x_bin_width   = 5;
            plot_.x_range_min   = 25;
            plot_.x_range_max   = 220;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "l_pt[1]";
        plot_.root_get_factor   = "";
        plot_.x_label           = "Sub-lead lepton p_{T} [GeV]";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_lpt1";
            plot_.y_label       = "Events / 5 GeV";
            plot_.x_bin_width   = 5;
            plot_.x_range_min   = 25;
            plot_.x_range_max   = 180;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);



        plot_.root_get_member   = "refEle_et";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} electron term [GeV]";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_refEle_et";
            plot_.y_label       = "Events / 10 GeV";
            plot_.x_bin_width   = 10;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 300;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "refEle_phi";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} electron term #phi";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_refEle_phi";
            plot_.y_label       = "Events";
            plot_.x_bin_width   = 0.2;
            plot_.x_range_min   = -3.2;
            plot_.x_range_max   = 3.2;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);



        plot_.root_get_member   = "refGamma_et";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} photon term [GeV]";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_refGamma_et";
            plot_.y_label       = "Events / 2 GeV";
            plot_.x_bin_width   = 2;
            plot_.x_range_min   = 10;
            plot_.x_range_max   = 55;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e1;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "refGamma_phi";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} photon term #phi";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_refGamma_phi";
            plot_.y_label       = "Events";
            plot_.x_bin_width   = 0.2;
            plot_.x_range_min   = -3.2;
            plot_.x_range_max   = 3.2;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);



        plot_.root_get_member   = "refTau_et";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} tau term [GeV]";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_refTau_et";
            plot_.y_label       = "Events / 2 GeV";
            plot_.x_bin_width   = 2;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 55;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "refJet_et";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} jet term [GeV]";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_refJet_et";
            plot_.y_label       = "Events / 10 GeV";
            plot_.x_bin_width   = 10;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 400;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "refJet_phi";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} jet term #phi";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_refJet_phi";
            plot_.y_label       = "Events";
            plot_.x_bin_width   = 0.2;
            plot_.x_range_min   = -3.2;
            plot_.x_range_max   = 3.2;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "refMuo_et";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} muon term [GeV]";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_refMuo_et";
            plot_.y_label       = "Events / 2 GeV";
            plot_.x_bin_width   = 2;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 55;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);


        plot_.root_get_member   = "softTerm_et";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} soft-term [GeV]";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_softTerm";
            plot_.y_label       = "Events / 2 GeV";
            plot_.x_bin_width   = 2;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 55;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);


        plot_.root_get_member   = "softTerm_phi";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} soft-term #phi";
            plot_.signal_region = "z_like_e";
            plot_.long_name     = "z_like_e_softTerm_phi";
            plot_.y_label       = "Events";
            plot_.x_bin_width   = 0.2;
            plot_.x_range_min   = -3.2;
            plot_.x_range_max   = 3.2;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

//        plot_.root_get_member   = "j_jvf[0]";
//        plot_.root_get_factor   = "";
//        plot_.x_label           = "JVF of leading jet";
//            plot_.signal_region = "z_like_e";
//            plot_.long_name     = "z_like_e_jvf0";
//            plot_.y_label       = "Events";
//            plot_.x_bin_width   = 0.05;
//            plot_.x_range_min   = 0;
//            plot_.x_range_max   = 1;
//            plot_.y_range_min   = 1e-1;
//            plot_.y_range_max   = 1e4;
//            vect_plot.push_back(plot_);
//
//        plot_.root_get_member   = "j_jvf[1]";
//        plot_.root_get_factor   = "";
//        plot_.x_label           = "JVF of second-leading jet";
//            plot_.signal_region = "z_like_e";
//            plot_.long_name     = "z_like_e_jvf1";
//            plot_.y_label       = "Events";
//            plot_.x_bin_width   = 0.05;
//            plot_.x_range_min   = 0;
//            plot_.x_range_max   = 1;
//            plot_.y_range_min   = 1e-1;
//            plot_.y_range_max   = 1e4;
//            vect_plot.push_back(plot_);
//
//        plot_.root_get_member   = "j_jvf[2]";
//        plot_.root_get_factor   = "";
//        plot_.x_label           = "JVF of third-leading jet";
//            plot_.signal_region = "z_like_e";
//            plot_.long_name     = "z_like_e_jvf2";
//            plot_.y_label       = "Events";
//            plot_.x_bin_width   = 0.05;
//            plot_.x_range_min   = 0;
//            plot_.x_range_max   = 1;
//            plot_.y_range_min   = 1e-1;
//            plot_.y_range_max   = 1e4;
//            vect_plot.push_back(plot_);


        plot_.is_Log = false;
        plot_.root_get_member   = "mll";
        plot_.root_get_factor   = "";
        plot_.x_label           = "m_{#mu#mu} [GeV]";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_mll";
            plot_.y_label       = "Events / 2 GeV";
            plot_.x_bin_width   = 2;
            plot_.x_range_min   = 70;
            plot_.x_range_max   = 110;
            plot_.y_range_min   = 0.0;
            plot_.y_range_max   = 1400;
            vect_plot.push_back(plot_);

        plot_.is_Log = true;
        plot_.root_get_member   = "nJets";
        plot_.root_get_factor   = "";
        plot_.x_label           = "Number of jets";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_nJets";
            plot_.y_label       = "Events";
            plot_.x_bin_width   = 1;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 8;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e5;
            vect_plot.push_back(plot_);


        plot_.root_get_member   = "met";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} [GeV]";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_met";
            plot_.y_label       = "Events / 5 GeV";
            plot_.x_bin_width   = 5;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 180;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "metPhi";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} #phi";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_metPhi";
            plot_.y_label       = "Events";
            plot_.x_bin_width   = 0.2;
            plot_.x_range_min   = -3.2;
            plot_.x_range_max   = 3.2;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "l_pt[0]";
        plot_.root_get_factor   = "";
        plot_.x_label           = "Lead lepton p_{T} [GeV]";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_lpt0";
            plot_.y_label       = "Events / 5 GeV";
            plot_.x_bin_width   = 5;
            plot_.x_range_min   = 25;
            plot_.x_range_max   = 220;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "l_pt[1]";
        plot_.root_get_factor   = "";
        plot_.x_label           = "Sub-lead lepton p_{T} [GeV]";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_lpt1";
            plot_.y_label       = "Events / 5 GeV";
            plot_.x_bin_width   = 5;
            plot_.x_range_min   = 25;
            plot_.x_range_max   = 180;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "refEle_et";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} electron term [GeV]";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_refEle_et";
            plot_.y_label       = "Events / 2 GeV";
            plot_.x_bin_width   = 2;
            plot_.x_range_min   = -3.2;
            plot_.x_range_max   = 300;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);



        plot_.root_get_member   = "refGamma_et";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} photon term [GeV]";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_refGamma_et";
            plot_.y_label       = "Events / 2 GeV";
            plot_.x_bin_width   = 2;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 55;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "refTau_et";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} tau term [GeV]";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_refTau_et";
            plot_.y_label       = "Events / 2 GeV";
            plot_.x_bin_width   = 2;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 55;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "refJet_et";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} jet term [GeV]";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_refJet_et";
            plot_.y_label       = "Events / 10 GeV";
            plot_.x_bin_width   = 10;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 400;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "refJet_phi";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} jet term #phi";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_refJet_phi";
            plot_.y_label       = "Events";
            plot_.x_bin_width   = 0.2;
            plot_.x_range_min   = -3.2;
            plot_.x_range_max   = 3.2;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "refMuo_et";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} muon term [GeV]";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_refMuo_et";
            plot_.y_label       = "Events / 10 GeV";
            plot_.x_bin_width   = 10;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 300;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

  //      plot_.root_get_member   = "refMuo_phi";
  //      plot_.root_get_factor   = "";
  //      plot_.x_label           = "TST #slash{E}_{T} muon term #phi";
  //          plot_.signal_region = "z_like_m";
  //          plot_.long_name     = "z_like_m_refMuo_phi";
  //          plot_.y_label       = "Events";
  //          plot_.x_bin_width   = 0.1;
  //          plot_.x_range_min   = 0;
  //          plot_.x_range_max   = 3.2;
  //          plot_.y_range_min   = 1e-1;
  //          plot_.y_range_max   = 1e4;
  //          vect_plot.push_back(plot_);


      //  plot_.is_Log = false;
        plot_.root_get_member   = "softTerm_et";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} soft-term [GeV]";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_softTerm";
            plot_.y_label       = "Events / 2 GeV";
            plot_.x_bin_width   = 2;
            plot_.x_range_min   = 0;
            plot_.x_range_max   = 55;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);

        plot_.root_get_member   = "softTerm_phi";
        plot_.root_get_factor   = "";
        plot_.x_label           = "TST #slash{E}_{T} soft-term #phi";
            plot_.signal_region = "z_like_m";
            plot_.long_name     = "z_like_m_softTerm_phi";
            plot_.y_label       = "Events";
            plot_.x_bin_width   = 0.2;
            plot_.x_range_min   = -3.2;
            plot_.x_range_max   = 3.2;
            plot_.y_range_min   = 1e-1;
            plot_.y_range_max   = 1e4;
            vect_plot.push_back(plot_);


//        plot_.root_get_member   = "j_jvf[0]";
//        plot_.root_get_factor   = "";
//        plot_.x_label           = "JVF of leading jet";
//            plot_.signal_region = "z_like_m";
//            plot_.long_name     = "z_like_m_jvf0";
//            plot_.y_label       = "Events";
//            plot_.x_bin_width   = 0.05;
//            plot_.x_range_min   = 0;
//            plot_.x_range_max   = 1;
//            plot_.y_range_min   = 1e-1;
//            plot_.y_range_max   = 1e4;
//            vect_plot.push_back(plot_);
//
//        plot_.root_get_member   = "j_jvf[1]";
//        plot_.root_get_factor   = "";
//        plot_.x_label           = "JVF of second-leading jet";
//            plot_.signal_region = "z_like_m";
//            plot_.long_name     = "z_like_m_jvf1";
//            plot_.y_label       = "Events";
//            plot_.x_bin_width   = 0.05;
//            plot_.x_range_min   = 0;
//            plot_.x_range_max   = 1;
//            plot_.y_range_min   = 1e-1;
//            plot_.y_range_max   = 1e4;
//            vect_plot.push_back(plot_);
//
//        plot_.root_get_member   = "j_jvf[2]";
//        plot_.root_get_factor   = "";
//        plot_.x_label           = "JVF of third-leading jet";
//            plot_.signal_region = "z_like_m";
//            plot_.long_name     = "z_like_m_jvf2";
//            plot_.y_label       = "Events";
//            plot_.x_bin_width   = 0.05;
//            plot_.x_range_min   = 0;
//            plot_.x_range_max   = 1;
//            plot_.y_range_min   = 1e-1;
//            plot_.y_range_max   = 1e4;
//            vect_plot.push_back(plot_);




















































    return vect_plot;
}

int main(int argc, char* argv[])
{
    string request_plot = "";
    bool plot_requested = false;

    int request_group = 0;
    bool group_requested = false;

    bool no_systematics = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "/p") == 0) {
            request_plot = argv[++i];
            if (request_plot != "") plot_requested = true;
        }
        else if (strcmp(argv[i], "/g") == 0) {
            request_group = atoi(argv[++i]);
            if (request_group > 0) group_requested = true;
        }
        else if (strcmp(argv[i], "/nosys") == 0) {
            no_systematics = true;
        }
        else {
            cout << "Super-Razor    Error (fatal): Bad arguments." << endl;
            exit(1);
        }
    }

    vector<plot> request_plots = definePlots();
    vector<processed_tree> hft_trees = defineTrees(request_plot);
    vector<hft_systematic> hft_syst = defineSystematics(request_plot);


    // SIGNAL REGION DEFS
    // SELECTIONS

    // define the all-important TCut strings
    map <string, string> tcuts;
    stringstream tcuts_join;

    // ---------------------------------------------- //
    //  signal region definitions                       
    // ---------------------------------------------- //
    tcuts_join << "isOS && isElEl && lept1Pt>25 && lept2Pt>20 && abs(mll-91.2)<10.";
    tcuts["val_ZEE"] = tcuts_join.str(); tcuts_join.str("");

    tcuts_join << "isOS && isMuMu && lept1Pt>20 && lept2Pt>20 && abs(mll-91.2)<15.";
    tcuts["val_ZMM"] = tcuts_join.str(); tcuts_join.str("");

    tcuts_join << "isOS && isElMu && nCentralBJets>=1";
    tcuts["ttlike"] = tcuts_join.str(); tcuts_join.str("");

    tcuts_join << "nLeptons==2 && nMuons==2 && abs(m_mumu_mu14-91.2)<20";
    tcuts["z_mumu_14"] = tcuts_join.str(); tcuts_join.str("");

    tcuts_join << "nLeptons==2 && nMuons==2 && abs(m_mumu_mu26-91.2)<20";
    tcuts["z_mumu_26"] = tcuts_join.str(); tcuts_join.str("");

    tcuts_join << "nLeptons==2 && nElectrons==2 && abs(mll-91.2)<20";
    tcuts["z_ee"] = tcuts_join.str(); tcuts_join.str("");

    tcuts_join << "nLeptons==2 && nElectrons==2 && nJets>=1";
    tcuts["z_like_e"] = tcuts_join.str(); tcuts_join.str("");

    tcuts_join << "nLeptons==2 && nMuons==2 && nJets>=0";
    tcuts["z_like_m"] = tcuts_join.str(); tcuts_join.str("");

    tcuts_join << "nLeptons==1 && nElectrons==1 && met>50";
    tcuts["w_like_e"] = tcuts_join.str(); tcuts_join.str("");

    tcuts_join << "nLeptons==1 && nMuons==1 && met>50";
    tcuts["w_like_m"] = tcuts_join.str(); tcuts_join.str("");


    tcuts["none"] = "1";

    TCanvas* ttemp = new TCanvas("canvas_1", "", 768, 768);

    TGuiUtils* _utils = new TGuiUtils();

    vector<vector<vector<double>>> table_columns;
    vector<string> sr;

    string last_sr = "";

    for (int r_ = 0; r_ < request_plots.size(); r_++) {

        plot plot_ = request_plots[r_];
        sr.push_back(plot_.signal_region);

        if (plot_requested && plot_.long_name.compare(request_plot) != 0) { continue; };
        if (group_requested && plot_.run_group != request_group) { continue; };

        if (plot_.signal_region.compare(last_sr) != 0) {
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                hft_trees[t_].all_lists_exist = false;
            }
            last_sr = plot_.signal_region;
        }

        vector<vector<vector<double>>> table_col_;// column is per process

        cout << setprecision(4);
        cout << std::fixed;
        cout << "+++++++++++++ " << plot_.signal_region << " +++++++++++++" << endl;

        // Complete/finalize TCut strings.
        stringstream signal_region_stream;
        signal_region_stream << "( ";
        signal_region_stream << tcuts[plot_.signal_region];
        if (plot_.request_cuts.compare("") != 0) signal_region_stream << " && " << plot_.request_cuts;
        signal_region_stream << " )";// Don't forget the close parentheses

        string signal_tcut = signal_region_stream.str();

        TH1F** hist_central = new TH1F*[hft_trees.size()];
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            stringstream hist_;
            if (plot_.root_get_member_name.compare("") == 0) {
                hist_ << "h_" << plot_.signal_region << "_" << plot_.root_get_member << "_" << hft_trees[t_].treename;
            }
            else {
                hist_ << "h_" << plot_.signal_region << "_" << plot_.root_get_member_name << "_" << hft_trees[t_].treename;
            }
            // calculate the number of bins
            int n_bins = floor((plot_.x_range_max - plot_.x_range_min) / (plot_.x_bin_width) + 0.5);
            hist_central[t_] = new TH1F(hist_.str().data(), hist_.str().data(), n_bins, plot_.x_range_min, plot_.x_range_max);
            hist_central[t_]->Sumw2();
        }

        int num_bins = hist_central[0]->GetNbinsX();

        double*** systematic_up = new double**[num_bins];
        for (int b_ = 0; b_ < num_bins; b_++) {
            systematic_up[b_] = new double*[hft_trees.size()];
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                systematic_up[b_][t_] = new double[hft_syst.size()];
                for (int s_ = 0; s_ < hft_syst.size(); s_++) {
                    systematic_up[b_][t_][s_] = 0;
                }
            }
        }

        double*** systematic_down = new double**[num_bins];
        for (int b_ = 0; b_ < num_bins; b_++) {
            systematic_down[b_] = new double*[hft_trees.size()];
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                systematic_down[b_][t_] = new double[hft_syst.size()];
                for (int s_ = 0; s_ < hft_syst.size(); s_++) {
                    systematic_down[b_][t_][s_] = 0;
                }
            }
        }

        vector<TH1F> total_ststematics;
        vector<vector<double>> tree_total;
        vector<vector<double>> tree_stat_error;

        vector<double> single_tree_total;

        for (int b_ = 0; b_ < num_bins; b_++) {
            vector<double> entry;
            tree_total.push_back(entry);
            tree_stat_error.push_back(entry);
        }

        // define central yields table
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            stringstream filename_stream;
            if(hft_trees[t_].file == MC_FILE && hft_trees[t_].tree_type == k_tree_mc) {
                filename_stream << DIR_MC << MC_FILE;
            }
            else{
                filename_stream << DIR_DATA << hft_trees[t_].file;
            }








            string filename_ = filename_stream.str();

            stringstream chain_name;// select CENTRAL
            chain_name << hft_trees[t_].treename << "_" << m_central;

            TChain* chain_;
            chain_ = new TChain(chain_name.str().data());
            chain_->Add(filename_.data());

            if (chain_->IsZombie()) {
                cout << chain_name.str() << "not found." << endl;
                continue;
            }

            TCut sel = TCut(signal_tcut.data());// key TCut

            TCut weight = TCut("eventweight");// selects eventweight leaf
            TCut blind = TCut(hft_trees[t_].blinding_cuts.data());

            // START scale factor (SF) re-weighting
            TCut scale_factor = TCut("1.0");
            size_t is_sr1 = plot_.signal_region.find("sr1");
            size_t is_sr0 = plot_.signal_region.find("sr0");
            if(is_sr1!=string::npos){
                if     (hft_trees[t_].treename == "ttbar") scale_factor = TCut("1.54");
                else if(hft_trees[t_].treename == "W")  scale_factor = TCut("1.54");
                else if(hft_trees[t_].treename == "Z")  scale_factor = TCut("1.54");
            }
            if (hft_trees[t_].tree_type == k_tree_mc && hft_trees[t_].special_cuts.compare("") != 0) {
                scale_factor = TCut(hft_trees[t_].special_cuts.data());
            }
            // END scale factor re-weighting

            // START reweighting segment
            TCut tc_reweight = TCut("1.0");
            if (plot_.reweight.compare("") != 0 && hft_trees[t_].tree_type == k_tree_mc) {
                tc_reweight = TCut(plot_.reweight.data());
            }
            // END reweighting segment

            // START selection and event-list
            if (!hft_trees[t_].all_lists_exist) {
                string list_name = string("list_") + hist_central[t_]->GetName();

                stringstream draw_list;
                draw_list << ">> " << list_name; // hist_central[t_]->GetName();
                stringstream draw_hist;
                draw_hist << plot_.root_get_member << plot_.root_get_factor << " >>+ " << hist_central[t_]->GetName();


                chain_->Draw(draw_list.str().data(), (sel + blind)); // drawing list

                TEventList* eventList = (TEventList*)gDirectory->Get(list_name.data());
                chain_->SetEventList(eventList);

                hft_trees[t_].central_list = *eventList;

                chain_->Draw(draw_hist.str().data(), weight * tc_reweight * scale_factor); // drawing histogram

                // Draw overflow
                _utils->moveUnderOverFlow(hist_central[t_], 2);
            }
            else {
                stringstream draw_hist;
                draw_hist << plot_.root_get_member << plot_.root_get_factor << " >>+ " << hist_central[t_]->GetName();

                chain_->SetEventList(&hft_trees[t_].central_list);
                chain_->Draw(draw_hist.str().data(), weight * tc_reweight * scale_factor); // drawing histogram

                // Draw overflow
                _utils->moveUnderOverFlow(hist_central[t_], 2);
            }
            // END selection and event-list

            double error_ = 0;
            double integral_ = hist_central[t_]->IntegralAndError(0, hist_central[t_]->GetNbinsX() + 1, error_);

            cout << hist_central[t_]->GetName() << ":" << integral_ << endl;

            for (int b_ = 1; b_ <= num_bins; b_++) {
                tree_total[b_ - 1].push_back(hist_central[t_]->GetBinContent(b_));
                tree_stat_error[b_ - 1].push_back(hist_central[t_]->GetBinError(b_));
            }

            single_tree_total.push_back(integral_);

            delete chain_;
            // delete eventList; in gDirectory, not needed.draw_hist
        }

        cout << ">";

        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            stringstream filename_stream;
            #warning HAVE NOT SET UP SYSTMEATIC TREE/file LOOP
            filename_stream << DIR_DATA << hft_trees[t_].file;




            //filename_stream << DIR_0 << hft_trees[t_].file;
            string filename_ = filename_stream.str();

            if (!DEBUG_NTC) {
                cout << "<< " << hft_trees[t_].treename;
            }
            int lcount = 0;

            for (int s_ = 0; s_ < hft_syst.size(); s_++) {

                if (no_systematics) continue;

                if (hft_trees[t_].systematic_groups.size() == 0) continue;

                bool add_this_systematic = false;
                for (int i = 0; i < hft_trees[t_].systematic_groups.size(); i++) {
                    if (hft_trees[t_].systematic_groups[i] == hft_syst[s_].binding_group) add_this_systematic = true;
                }

                if (!add_this_systematic) continue;

                stringstream chain_name_up;
                stringstream chain_name_down;

                // classify systematics
                if (hft_syst[s_].sys_type == k_sys_object) {
                    chain_name_up << hft_trees[t_].treename << "_" << hft_syst[s_].basename << hft_syst[s_].up_name;
                    chain_name_down << hft_trees[t_].treename << "_" << hft_syst[s_].basename << hft_syst[s_].down_name;
                }
                else if (hft_syst[s_].sys_type == k_sys_one_sided_object) {
                    chain_name_up << hft_trees[t_].treename << "_" << hft_syst[s_].basename;
                    chain_name_down << hft_trees[t_].treename << "_" << hft_syst[s_].basename;
                }
                else if (hft_syst[s_].sys_type == k_sys_weight) {
                    chain_name_up << hft_trees[t_].treename << "_" << m_central;
                    chain_name_down << hft_trees[t_].treename << "_" << m_central;
                    //reweight_flag = true;
                }
                else if (hft_syst[s_].sys_type == k_sys_adhoc) {
                    chain_name_up << hft_trees[t_].treename << "_" << m_central;
                    chain_name_down << hft_trees[t_].treename << "_" << m_central;
                    //reweight_flag = true;
                }
                else {
                    continue;
                }

                TChain* chain_up;
                chain_up = new TChain(chain_name_up.str().data());
                chain_up->Add(filename_.data());

                TChain* chain_down;
                chain_down = new TChain(chain_name_down.str().data());
                chain_down->Add(filename_.data());

                if (chain_up->IsZombie()) {
                    cout << chain_name_up.str() << "not found." << endl;
                    continue;
                }
                if (chain_down->IsZombie()) {
                    cout << chain_name_down.str() << "not found." << endl;
                    continue;
                }

                TCut sel = TCut(signal_tcut.data());// key TCut
                TCut weight = TCut("eventweight");// selects eventweight leaf
                TCut blind = TCut(hft_trees[t_].blinding_cuts.data());

                // apply direct weighting
                TCut syst_var_up;
                TCut syst_var_down;
                if (hft_syst[s_].sys_type == k_sys_weight) {
                    stringstream weight_name_up;
                    weight_name_up << hft_syst[s_].preface << hft_syst[s_].basename << hft_syst[s_].up_name;
                    stringstream weight_name_down;
                    weight_name_down << hft_syst[s_].preface << hft_syst[s_].basename << hft_syst[s_].down_name;

                    syst_var_up = TCut(weight_name_up.str().data());
                    syst_var_down = TCut(weight_name_down.str().data());
                }
                else if (hft_syst[s_].sys_type == k_sys_adhoc) {
                    syst_var_up = TCut(hft_syst[s_].adhoc_tcut.data());
                    syst_var_down = TCut(hft_syst[s_].adhoc_tcut2.data());
                }
                else {
                    syst_var_up = TCut("1.0");
                    syst_var_down = TCut("1.0");
                }
            
                // START scale factor (SF) re-weighting (systematic side)
                TCut scale_factor = TCut("1.0");
                size_t is_sr1 = plot_.signal_region.find("sr1");
                size_t is_sr0 = plot_.signal_region.find("sr0");
                if(is_sr1!=string::npos){
                    if     (hft_trees[t_].treename == "Top") scale_factor = TCut("1.06");
                    else if(hft_trees[t_].treename == "WW")  scale_factor = TCut("1.04");
                    else if(hft_trees[t_].treename == "ZV")  scale_factor = TCut("1.19");
                }
             //   if (hft_trees[t_].tree_type == k_tree_mc && hft_trees[t_].special_cuts.compare("") != 0) {
             //       scale_factor = TCut(hft_trees[t_].special_cuts.data());
             //   }
                // END scale factor re-weighting

                // START reweighting segment
                TCut tc_reweight = TCut("1.0");
                if (plot_.reweight.compare("") != 0 && hft_trees[t_].tree_type == k_tree_mc) {
                    tc_reweight = TCut(plot_.reweight.data());
                }
                // END reweighting segment

                int n_bins = floor((plot_.x_range_max - plot_.x_range_min) / (plot_.x_bin_width) + 0.5);

                // temp histogram to get the variation
                TH1F* temp_hist_up = new TH1F("temp_hist_up", "temp_hist_up", n_bins, plot_.x_range_min, plot_.x_range_max);
                TH1F* temp_hist_down = new TH1F("temp_hist_down", "temp_hist_down", n_bins, plot_.x_range_min, plot_.x_range_max);

                stringstream draw_up;
                draw_up << plot_.root_get_member << plot_.root_get_factor << " >> +" << "temp_hist_up";
                stringstream draw_down;
                draw_down << plot_.root_get_member << plot_.root_get_factor << " >> +" << "temp_hist_down";

                if (hft_syst[s_].sys_type == k_sys_weight || hft_syst[s_].sys_type == k_sys_adhoc) {
                    chain_up->SetEventList(&hft_trees[t_].central_list);
                    chain_down->SetEventList(&hft_trees[t_].central_list);
                    chain_up->Draw(draw_up.str().data(), weight * syst_var_up * tc_reweight * scale_factor);
                    chain_down->Draw(draw_down.str().data(), weight * syst_var_down * tc_reweight * scale_factor);

                    // Draw overflow
                    _utils->moveUnderOverFlow(temp_hist_up, 2);
                    _utils->moveUnderOverFlow(temp_hist_down, 2);
                    
                }
                // else if (hft_syst[s_].sys_type == k_sys_adhoc) {
                //     chain_up->Draw(draw_up.str().data(), sel * weight * syst_var_up * tc_reweight);
                //     chain_down->Draw(draw_down.str().data(), sel * weight * syst_var_down * tc_reweight);
                // } // needless duplication
                else if (hft_trees[t_].all_lists_exist) {
                    string list_name_up = string("eventlist_") + hft_syst[s_].basename + hft_syst[s_].up_name;
                    string list_name_down = string("eventlist_") + hft_syst[s_].basename + hft_syst[s_].down_name;

                    chain_up->SetEventList(&hft_trees[t_].all_lists[list_name_up]);
                    chain_down->SetEventList(&hft_trees[t_].all_lists[list_name_down]);

                    chain_up->Draw(draw_up.str().data(), weight * tc_reweight * scale_factor);
                    chain_down->Draw(draw_down.str().data(), weight * tc_reweight * scale_factor);
                    
                    // Draw overflow
                    _utils->moveUnderOverFlow(temp_hist_up, 2);
                    _utils->moveUnderOverFlow(temp_hist_down, 2);
                }
                else {
                    string list_name_up = string("eventlist_") + hft_syst[s_].basename + hft_syst[s_].up_name;
                    string list_name_down = string("eventlist_") + hft_syst[s_].basename + hft_syst[s_].down_name;

                    string draw_list_up = string(">> ") + list_name_up;
                    string draw_list_down = string(">> ") + list_name_down;

                    chain_up->Draw(draw_list_up.data(), sel); // drawing list
                    chain_down->Draw(draw_list_down.data(), sel); // drawing list

                    TEventList* eventList_up = (TEventList*)gDirectory->Get(list_name_up.data());
                    TEventList* eventList_down = (TEventList*)gDirectory->Get(list_name_down.data());

                    chain_up->SetEventList(eventList_up);
                    chain_down->SetEventList(eventList_down);

                    hft_trees[t_].all_lists[list_name_up] = *eventList_up;
                    hft_trees[t_].all_lists[list_name_down] = *eventList_down;

                    chain_up->Draw(draw_up.str().data(), weight * tc_reweight * scale_factor);
                    chain_down->Draw(draw_down.str().data(), weight * tc_reweight * scale_factor);

                    // Draw overflow
                    _utils->moveUnderOverFlow(temp_hist_up, 2);
                    _utils->moveUnderOverFlow(temp_hist_down, 2);
                }

                double int_up = temp_hist_up->Integral(0, -1) - single_tree_total[t_];
                double int_down = temp_hist_down->Integral(0, -1) - single_tree_total[t_];

                bool use_up = false;
                bool use_down = false;
                // sign and sidedness adjustments
                if (hft_syst[s_].sys_type == k_sys_one_sided_object) {
                    if (int_up > 0) {
                        int_down = 0.0;
                        use_up = true;
                    }
                    else {
                        int_up = 0.0;
                        use_down = true;
                    }
                }
                else if (hft_syst[s_].sys_type == k_sys_adhoc) {
                    int_up = abs(int_up);
                    int_down = -abs(int_down);
                }

                // systematic_up[t_][s_] = int_up;
                // systematic_down[t_][s_] = int_down;

                for (int b_ = 1; b_ <= num_bins; b_++) {
                    systematic_up[b_ - 1][t_][s_] = temp_hist_up->GetBinContent(b_) - tree_total[b_ - 1][t_];
                    systematic_down[b_ - 1][t_][s_] = temp_hist_down->GetBinContent(b_) - tree_total[b_ - 1][t_];

                    if (use_up)  systematic_down[b_ - 1][t_][s_] = 0;
                    if (use_down) systematic_up[b_ - 1][t_][s_] = 0;
                }

                if (DEBUG_NTC) {
                    // cout << hft_trees[t_].treename << " > " << hft_syst[s_].basename << " ";
                    // cout << tree_total[t_] << ",+" << int_up / tree_total[t_] << ",-" << int_down / tree_total[t_] << endl;
                    //cout << hft_trees[t_].treename << " > " << hft_syst[s_].basename << " ";
                    //cout << single_tree_total[t_] << ",+" << int_up << ",-" << int_down << endl;

                    cout << pad_width(hft_trees[t_].treename + " > " + hft_syst[s_].basename, 20);
                    cout << pad_width(to_string(single_tree_total[t_]) + ",", 16);
                    cout << pad_width("+" + to_string(int_up), 12) << pad_width("-" + to_string(int_down), 12) << endl;

                }

                delete temp_hist_up;
                temp_hist_up = 0;
                delete temp_hist_down;
                temp_hist_down = 0;
                delete chain_up;
                delete chain_down;

                lcount++;

                if (!(lcount % 10 == 0 || lcount == 0)) cout << "+";
                else cout << "/";
                cout << flush;
            }

            hft_trees[t_].all_lists_exist = true; // set when you've run through the tree once

            if (!DEBUG_NTC) {
                cout << endl;
            }
        }

        for (int b_ = 1; b_ <= num_bins; b_++) {
            vector<vector<double> > entry;
            // for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            // 
            //     vector<double> single_entry(4, 0.0);
            //     entry.push_back(single_entry);
            // 
            // }
            table_col_.push_back(entry);
        }

        // per grouping table
        for (int b_ = 0; b_ < num_bins; b_++) {
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                vector<double> single_entry(4, 0.0);

                double syst_up = 0;
                double syst_down = 0;

                for (int s_ = 0; s_ < hft_syst.size(); s_++) {
                    double up_var = systematic_up[b_][t_][s_];
                    double down_var = systematic_down[b_][t_][s_];

                    if (up_var > 0) syst_up += up_var * up_var;
                    else syst_down += up_var * up_var;

                    if (down_var < 0) syst_down += down_var * down_var;
                    else syst_up += down_var * down_var;
                }

                syst_up = sqrt(syst_up);
                syst_down = sqrt(syst_down);

                single_entry[0] = tree_total[b_][t_];
                single_entry[1] = tree_stat_error[b_][t_];
                single_entry[2] = syst_up;
                single_entry[3] = syst_down;

                table_col_[b_].push_back(single_entry);
            }
        }

        // gather systematic groups
        vector<int> syst_groups;
        for (int s_ = 0; s_ < hft_syst.size(); s_++) {
            bool add_group = true;
            for (int g_ = 0; g_ < syst_groups.size(); g_++) {
                if (syst_groups[g_] == hft_syst[s_].binding_group) add_group = false;
            }

            if (add_group) syst_groups.push_back(hft_syst[s_].binding_group);
        }

        vector<double> total_stat_up;
        vector<double> total_stat_down;

        for (int b_ = 0; b_ < num_bins; b_++) {
            total_stat_up.push_back(0);
            total_stat_down.push_back(0);
        }

        // totals
        for (int b_ = 0; b_ < num_bins; b_++) {
            for (int g_ = 0; g_ < syst_groups.size(); g_++) {
                for (int s_ = 0; s_ < hft_syst.size(); s_++) {
                    if (hft_syst[s_].binding_group != syst_groups[g_]) continue;

                    double syst_up = 0;
                    double syst_down = 0;

                    for (int t_ = 0; t_ < hft_trees.size(); t_++) {

                        bool add_this_systematic = false;
                        for (int i = 0; i < hft_trees[t_].systematic_groups.size(); i++) {
                            if (hft_trees[t_].systematic_groups[i] == hft_syst[s_].binding_group) add_this_systematic = true;
                        }

                        if (!add_this_systematic) continue;

                        double up_var = systematic_up[b_][t_][s_];
                        double down_var = systematic_down[b_][t_][s_];

                        syst_up += up_var;
                        syst_down += down_var;
                    }

                    //total_stat_up += syst_up * syst_up;
                    //total_stat_down += syst_down * syst_down;

                    if (syst_up > 0) {
                        total_stat_up[b_] += syst_up * syst_up;
                    }
                    else {
                        total_stat_down[b_] += syst_up * syst_up;
                    }

                    if (syst_down < 0) {
                        total_stat_down[b_] += syst_down * syst_down;
                    }
                    else {
                        total_stat_up[b_] += syst_down * syst_down;
                    }
                }
            }
        }

        for (int b_ = 0; b_ < num_bins; b_++) {
            total_stat_up[b_] = sqrt(total_stat_up[b_]);
            total_stat_down[b_] = sqrt(total_stat_down[b_]);
        }

        vector<double> realtotal;
        vector<double> staterror;
        for (int b_ = 0; b_ < num_bins; b_++) {
            realtotal.push_back(0);
            staterror.push_back(0);
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
//                if (hft_trees[t_].tree_type == k_tree_data) continue;
                if ((hft_trees[t_].tree_type == k_tree_data) || (hft_trees[t_].tree_type == k_tree_signal)) continue;
                realtotal[b_] += tree_total[b_][t_];
                staterror[b_] += tree_stat_error[b_][t_] * tree_stat_error[b_][t_];
            }
            staterror[b_] = sqrt(staterror[b_]);
        }

        //cout << "Total" << " > " << realtotal << " \xB1" << staterror;
        //cout << " + " << total_stat_up << " - " << total_stat_down << endl;

        vector<double> total_entry(4, 0.0);

        // realtotal[b_];        // these are the ultimate results.
        // staterror[b_];        // these are the ultimate results.
        // total_stat_up[b_];    // these are the ultimate results.
        // total_stat_down[b_];  // these are the ultimate results.

        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///

        double* nominal_yield = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) nominal_yield[j] = 0.0;
        // calculate and store, per bin nominal_yield;
        TH1F* hist_prediction = new TH1F();
        delete hist_prediction;
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if ((hft_trees[t_].tree_type == k_tree_mc) || (hft_trees[t_].tree_type == k_tree_fakes)) {
                hist_prediction = (TH1F*)hist_central[t_]->Clone();
                break;
            }
        }
        // ---------------------------------
        //  dealing with negative bin yields
        // ---------------------------------
        for (int t_=0; t_ < hft_trees.size(); t_++) {
            if ((hft_trees[t_].tree_type == k_tree_fakes)) {
                for(int j=1; j <= hist_central[t_]->GetNbinsX(); j++) {
                    double bin_yield = 0.0;
                    bin_yield = hist_central[t_]->GetBinContent(j);
                    if (bin_yield < 0) {
                        hist_central[t_]->SetBinContent(j,0.0);
                    }
                } // over bins
            } // if fakes
        } // loop over trees

        bool skip_init = true;
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if ((hft_trees[t_].tree_type == k_tree_mc) || (hft_trees[t_].tree_type == k_tree_fakes)) {
                if (skip_init) {
                    skip_init = false;
                }
                else {
                    hist_prediction->Add(hist_central[t_], 1.0);
                }
            }
        }
        for (int j = 1; j <= hist_central[0]->GetNbinsX(); j++) {// note bin numbering convention
            nominal_yield[j - 1] = hist_prediction->GetBinContent(j);
        }
        delete hist_prediction;

        // X BIN CENTER
        double* x_bin_center = new double[hist_central[0]->GetNbinsX()];
        for (int j = 1; j <= hist_central[0]->GetNbinsX(); j++) x_bin_center[j - 1] = hist_central[0]->GetBinCenter(j);
        // X BIN DOWN
        double* x_bin_down = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) x_bin_down[j] = hist_central[0]->GetBinWidth(j) / 2;
        // X BIN UP
        double* x_bin_up = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) x_bin_up[j] = hist_central[0]->GetBinWidth(j) / 2;
        // Y BIN DOWN
        double* y_bin_down = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) y_bin_down[j] = sqrt(total_stat_down[j] * total_stat_down[j] + staterror[j] * staterror[j]);
        // Y BIN UP
        double* y_bin_up = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) y_bin_up[j] = sqrt(total_stat_up[j] * total_stat_up[j] + staterror[j] * staterror[j]);

        // order histograms smallest to largest
        vector<int> tree_ordering;
        {
            map<double, int> tree_map;
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                tree_map[single_tree_total[t_]] = t_;
            }
            for (auto t_ = tree_map.begin(); t_ != tree_map.end(); t_++) {
                tree_ordering.push_back(t_->second);
            }
        }

        cout << "-------------" << endl;

        // Output to pdf or eps

        gROOT->SetStyle("Plain");
        gStyle->SetOptStat(false);
        gStyle->SetTitleBorderSize(0);
        gStyle->SetLineWidth(1);
        gStyle->SetTextFont(42);

        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);

        gStyle->SetErrorX(0.0);

        stringstream out_name;
        out_name << plot_.signal_region << "_" << plot_.root_get_member;
        stringstream out_file;
        out_file << DIR_1 << plot_.long_name << ".eps";


 //       stringstream out_file_image;
 //       out_file_image << DIR_1 << plot_.long_name << ".png";

        TCanvas* tc = new TCanvas("canvas", "", 500, 500);
   //     TLegend* tl = new TLegend(0.68, 0.49, 0.88, 0.87);
        //TLegend* tl_1 = new TLegend(0.54,0.58,0.72,0.88);
        //TLegend* tl_2 = new TLegend(0.72,0.58,0.91,0.88);
        TLegend* tl_1 = new TLegend(0.55,0.74,0.70,0.91);
        //TLegend* tl_2 = new TLegend(0.76,0.64,0.91,0.91);
        TLegend* tl_2 = new TLegend(0.74,0.74,0.89,0.91);

        tc->Divide(1, 2);

        TPad* canvas_up = (TPad*)tc->GetListOfPrimitives()->FindObject("canvas_1");
        TPad* canvas_dw = (TPad*)tc->GetListOfPrimitives()->FindObject("canvas_2");

        // define the size
        double up_height = 0.78;
        double font_size_dw = 0.73;
        double dw_height = 0.30;

        // set pad size
        canvas_up->SetPad(0., 1 - up_height, 1., 1.);
        canvas_dw->SetPad(0., 0., 1., dw_height);
//        canvas_dw->SetPad(0., 0., 1., dw_height-0.2);   // DANTRIM
        canvas_up->SetFrameFillColor(0);
        canvas_up->SetFillColor(0);
        canvas_dw->SetFrameFillColor(0);
        canvas_dw->SetFillColor(0);

        //canvas_up->SetLeftMargin  (0.1  );
        //canvas_up->SetRightMargin (0.075);
        //canvas_up->SetBottomMargin(0.15 );
        //canvas_dw->SetLeftMargin  (0.1  );
        //canvas_dw->SetRightMargin (0.075);
        //canvas_dw->SetBottomMargin(0.4  );
        //canvas_dw->SetTopMargin   (0.05 );
        canvas_up->SetLeftMargin  (0.14  );
        canvas_up->SetRightMargin (0.05  );
        canvas_up->SetBottomMargin(0.12 );
        canvas_up->SetTopMargin(0.5 * canvas_up->GetTopMargin());
        canvas_dw->SetLeftMargin  (0.14  );
        canvas_dw->SetRightMargin (0.05  );
        canvas_dw->SetBottomMargin(0.4  );
    //    canvas_dw->SetTopMargin   (0.05 );

        // draw top figure
        canvas_up->cd();

        if (plot_.is_Log) canvas_up->SetLogy(true);

        tl_1->SetBorderSize(0);
        tl_1->SetFillColor(0);
        tl_1->SetTextFont(42);
        tl_1->SetTextSize(0.04);
        tl_1->SetLineWidth(0);
        tl_2->SetBorderSize(0);
        tl_2->SetFillColor(0);
        tl_2->SetTextFont(42);
        tl_2->SetTextSize(0.04);
        tl_2->SetLineWidth(0);

        // declare early for legend
        TGraphAsymmErrors* Data = new TGraphAsymmErrors();
        // data legend FIRST
        tl_1->AddEntry(Data, "Data", "p");

        // fill legend 2 for bkg
        TH1F* mcError = new TH1F("mcError", "mcError", 2,0,2);
        mcError->SetFillStyle(3354);
        mcError->SetFillColor(kBlack);
        mcError->SetLineColor(kRed);
        gStyle->SetHatchesSpacing(0.9);
        gStyle->SetHatchesLineWidth(2);
     //   asym_errors->SetFillStyle(3354); // 3004
        tl_2->AddEntry(mcError, "Total SM", "fl");

        THStack* hist_stack = new THStack(out_name.str().data(), "");

        //for (int t_ = 0; t_ < hft_trees.size(); t_++) {
        for (int i = 0; i < tree_ordering.size(); i++) {
            int t_ = tree_ordering[i];

            if ((hft_trees[t_].tree_type == k_tree_mc) || (hft_trees[t_].tree_type == k_tree_fakes)){
                hist_central[t_]->SetMarkerStyle(1);
                hist_central[t_]->SetLineColor(hft_trees[t_].color);
                hist_central[t_]->SetLineWidth(0);
                hist_central[t_]->SetFillColor(hft_trees[t_].color);
                hist_stack->Add(hist_central[t_]);
            }
        }
        //}
        
        // draw Total SM red line
        TH1F*  hist_sm = (TH1F*)hist_stack->GetStack()->Last()->Clone("hist_sm");
        hist_sm->SetLineColor(kRed);
        hist_sm->SetLineWidth(2);
        hist_sm->SetLineStyle(1);
        hist_sm->SetFillStyle(0);
    
        // Signal
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if (hft_trees[t_].tree_type == k_tree_signal) {
                hist_central[t_]->SetMarkerStyle(1);
                hist_central[t_]->SetLineColor(hft_trees[t_].color);
                hist_central[t_]->SetLineWidth(3);
                hist_central[t_]->SetLineStyle(2);
            }
        }

        // MC Legend
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if(t_%2==0 && (hft_trees[t_].tree_type == k_tree_mc || hft_trees[t_].tree_type == k_tree_fakes)) tl_2->AddEntry(hist_central[t_], hft_trees[t_].displayname.data(), "f");
            else if ( hft_trees[t_].tree_type == k_tree_mc || hft_trees[t_].tree_type == k_tree_fakes) tl_1->AddEntry(hist_central[t_], hft_trees[t_].displayname.data(), "f");
        }

        // Signal legend
        // signal trees are filled first
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if (hft_trees[t_].tree_type == k_tree_signal) {
                if(t_==0) { tl_1->AddEntry(hist_central[t_], hft_trees[t_].displayname.data(), "l"); }
                if(t_==1) { tl_1->AddEntry(hist_central[t_], hft_trees[t_].displayname.data(), "l"); }
            }

        }

        hist_stack->Draw("HIST");
        hist_stack->GetXaxis()->SetTitle(plot_.x_label.data());
        hist_stack->GetYaxis()->SetTitle(plot_.y_label.data());
        hist_stack->GetXaxis()->SetTitleFont(42);
        hist_stack->GetYaxis()->SetTitleFont(42);
        hist_stack->GetXaxis()->SetLabelFont(42);
        hist_stack->GetYaxis()->SetLabelFont(42);
        hist_stack->GetXaxis()->SetLabelOffset(999);
        hist_stack->GetXaxis()->SetTitleOffset(999);
        hist_stack->GetYaxis()->SetTitleOffset(0.85 * 1.28);
        hist_stack->GetYaxis()->SetLabelOffset(0.013);
        hist_stack->SetMinimum(plot_.y_range_min);
        hist_stack->SetMaximum(plot_.y_range_max);
        hist_stack->GetXaxis()->SetLabelSize(0.046);
        hist_stack->GetYaxis()->SetLabelSize(0.05);
        hist_stack->GetXaxis()->SetTitleSize(0.13);
        hist_stack->GetXaxis()->SetTitleSize(0.048);
        hist_stack->GetYaxis()->SetTitleSize(0.055);

        // draw Total SM
        hist_sm->Draw("hist same");

        TGraphAsymmErrors* asym_errors = new TGraphAsymmErrors(hist_central[0]->GetNbinsX(), x_bin_center, nominal_yield, x_bin_down, x_bin_up, y_bin_down, y_bin_up);
        gStyle->SetHatchesSpacing(0.9);
        gStyle->SetHatchesLineWidth(1);
        asym_errors->SetFillStyle(3354); // 3004
        asym_errors->SetFillColor(kGray + 3);
//        tl_2->AddEntry(asym_errors, "Bkg. Uncert.", "f");
        asym_errors->Draw("option same 02");

        // Decoration
        char annoyingLabel1[100] = "#bf{#it{ATLAS}}  Internal";
        char annoyingLabel2[100] = "#sqrt{s} = 13 TeV, 6.6 pb^{-1}";
      //  myText(0.14, 0.82, kBlack, annoyingLabel1);
      //  myText(0.14, 0.74, kBlack, annoyingLabel2);
        myText(0.18, 0.88, kBlack, annoyingLabel1);
        myText(0.16, 0.82, kBlack, annoyingLabel2);

    //    char annoyingLabel3;
        plot plt = request_plots[0];
        string channel = plot_.signal_region;
        size_t is_eeSR1a = channel.find("sr1a_EE");
        size_t is_mmSR1a = channel.find("sr1a_MM");
        size_t is_emSR1a = channel.find("sr1a_DF");
        size_t is_sfSR1a = channel.find("sr1a_SF");
         
        size_t is_eeSR1b = channel.find("sr1b_EE");
        size_t is_mmSR1b = channel.find("sr1b_MM");
        size_t is_emSR1b = channel.find("sr1b_DF");
        size_t is_sfSR1b = channel.find("sr1b_SF");

        size_t is_eeSR1c = channel.find("sr1c_EE");
        size_t is_mmSR1c = channel.find("sr1c_MM");
        size_t is_emSR1c = channel.find("sr1c_DF");
        size_t is_sfSR1c = channel.find("sr1c_SF");

        size_t is_crTopa = channel.find("crTopa");
        size_t is_crWWa  = channel.find("crWWa");
        size_t is_crZVa  = channel.find("crZVa");

        size_t is_crTopb = channel.find("crTopb");
        size_t is_crWWb  = channel.find("crWWb");
        size_t is_crZVb  = channel.find("crZVb");

        if     (is_eeSR1a!=string::npos) { char annoyingLabel3[100] = "SR2l-1a ee channel"; myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_mmSR1a!=string::npos) { char annoyingLabel3[100] = "SR2l-1a mm channel"; myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_emSR1a!=string::npos) { char annoyingLabel3[100] = "SR2l-1a DF channel"; myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_sfSR1a!=string::npos) { char annoyingLabel3[100] = "SR2l-1a SF channel"; myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_eeSR1b!=string::npos) { char annoyingLabel3[100] = "SR2l-1b ee channel"; myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_mmSR1b!=string::npos) { char annoyingLabel3[100] = "SR2l-1b mm channel"; myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_emSR1b!=string::npos) { char annoyingLabel3[100] = "SR2l-1b DF channel"; myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_sfSR1b!=string::npos) { char annoyingLabel3[100] = "SR2l-1b SF channel"; myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_eeSR1c!=string::npos) { char annoyingLabel3[100] = "SR2l-1b ee channel"; myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_mmSR1c!=string::npos) { char annoyingLabel3[100] = "SR2l-1b mm channel"; myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_emSR1c!=string::npos) { char annoyingLabel3[100] = "SR2l-1b DF channel"; myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_sfSR1c!=string::npos) { char annoyingLabel3[100] = "SR2l-1b SF channel"; myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_crTopa!=string::npos) { char annoyingLabel3[100] = "CRTop_{a}"; myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_crTopb!=string::npos) { char annoyingLabel3[100] =  "CR2l-Top"; myText(0.24, 0.73, kBlack, annoyingLabel3); }
        else if(is_crWWa !=string::npos) { char annoyingLabel3[100] = "CRWW_{a}";  myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_crWWb !=string::npos) { char annoyingLabel3[100] =  "CR2l-WW";  myText(0.24, 0.73, kBlack, annoyingLabel3); }
        else if(is_crZVa !=string::npos) { char annoyingLabel3[100] = "CRZV_{a}";  myText(0.18, 0.73, kBlack, annoyingLabel3); }
        else if(is_crZVb !=string::npos) { char annoyingLabel3[100] =  "CR2l-ZV";  myText(0.24, 0.73, kBlack, annoyingLabel3); }

  //      myText(0.20, 0.79, kBlack, annoyingLabel3);
         
        string var_name = plot_.x_label;
        size_t find_g = var_name.find("[GeV]");
        if (find_g != string::npos) {
            var_name = var_name.substr(0, find_g);
        }
        string title_ = plot_.signal_region + ": " + var_name;

        // this draws the label above the plots
        //myText(0.14, 0.940, kBlack, title_.c_str());

        tl_1->Draw();
        tl_2->Draw();


        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if (hft_trees[t_].tree_type == k_tree_data) {
                convertErrorsToPoisson(hist_central[t_], Data);
                break;
            }
        }

        // Signal
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if (hft_trees[t_].tree_type == k_tree_signal) {
                hist_central[t_]->Draw("same HIST");
            }
        }

        Data->SetLineWidth(2);
        Data->SetMarkerStyle(20);
        Data->SetMarkerSize(1.1);
        Data->SetLineColor(1);
        Data->Draw("option same pz");

        // draw bottom figure
        // draw bottom figure
        // draw bottom figure

        canvas_dw->cd();

        TLine line_ = TLine(plot_.x_range_min, 1.0, plot_.x_range_max, 1.0);
        TLine line_up = TLine(plot_.x_range_min, 1.5, plot_.x_range_max, 1.5);
        TLine line_down = TLine(plot_.x_range_min, 0.5, plot_.x_range_max, 0.5);

        line_up.SetLineStyle(3);
        line_down.SetLineStyle(3);

        line_.SetLineColor(kRed);
        line_.SetLineStyle(2);

        //
        double* nominal_one = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) nominal_one[j] = 1.0;

        double* value_zero = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) value_zero[j] = 0.0;

        double* y_bin_down_ratio = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) {
            if (nominal_yield[j] != 0) {
                y_bin_down_ratio[j] = y_bin_down[j] / nominal_yield[j];
            }
            else {
                y_bin_down_ratio[j] = 0;
            }
        }

        double* y_bin_up_ratio = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) {
            if (nominal_yield[j] != 0) {
                y_bin_up_ratio[j] = y_bin_up[j] / nominal_yield[j];
            }
            else {
                y_bin_up_ratio[j] = 0;
            }
        }

        double* data_ratio = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) {
            if (nominal_yield[j] > 0 && hist_central[hft_trees.size() - 1]->GetBinContent(j + 1) > 0) {
                data_ratio[j] = hist_central[hft_trees.size() - 1]->GetBinContent(j + 1) / nominal_yield[j];
            }
            else {
                data_ratio[j] = 100.0;
            }
        }

        double* data_err_down_ratio = new double[hist_central[0]->GetNbinsX()];
        double* data_err_up_ratio = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) {
            if (nominal_yield[j] > 0) {

                data_err_down_ratio[j] = 1.0 - ((nominal_yield[j] - Data->GetErrorYlow(j)) / nominal_yield[j]);

                data_err_up_ratio[j] = ((nominal_yield[j] + Data->GetErrorYhigh(j)) / nominal_yield[j]) - 1.0;
            }
            else {
                data_err_down_ratio[j] = 0.0;
                data_err_up_ratio[j] = 0.0;
            }
        }

        TGraphAsymmErrors* asym_errors2 = new TGraphAsymmErrors(hist_central[0]->GetNbinsX(), x_bin_center, nominal_one, x_bin_down, x_bin_up, y_bin_down_ratio, y_bin_up_ratio);

        asym_errors2->SetFillStyle(3354); // 3004
        asym_errors2->SetFillColor(kGray + 3);

        TH1F* hist_ratio = (TH1F*)hist_central[0]->Clone();
        TGraphAsymmErrors* Data2 = new TGraphAsymmErrors(hist_central[0]->GetNbinsX(), x_bin_center, data_ratio, value_zero, value_zero, data_err_down_ratio, data_err_up_ratio);

        hist_ratio->SetName("ratio_");
        hist_ratio->SetMaximum(2.0);
        hist_ratio->SetMinimum(0.0);
        hist_ratio->GetXaxis()->SetTitleFont(42);
        hist_ratio->GetXaxis()->SetTitleFont(42);
        hist_ratio->GetYaxis()->SetTitleFont(42);
        hist_ratio->GetXaxis()->SetLabelFont(42);
        hist_ratio->GetYaxis()->SetLabelFont(42);
        hist_ratio->GetXaxis()->SetTitle(plot_.x_label.data());
        hist_ratio->GetYaxis()->SetTitle("Data/SM");
        hist_ratio->GetXaxis()->SetTickLength(0.06);
        hist_ratio->GetXaxis()->SetLabelSize(0.13);
        hist_ratio->GetYaxis()->SetLabelSize(0.13);
        hist_ratio->GetXaxis()->SetTitleSize(0.14);
        hist_ratio->GetYaxis()->SetTitleSize(0.14);
        hist_ratio->GetYaxis()->SetTitleOffset(0.45);
        hist_ratio->GetYaxis()->SetLabelOffset(0.98 * 0.013);
        hist_ratio->GetXaxis()->SetLabelOffset(0.02);
        //hist_ratio->GetYaxis()->SetNdivisions(4, false);
        hist_ratio->GetYaxis()->SetNdivisions(5);
        hist_ratio->Draw("AXIS");

        asym_errors2->Draw("option same 02");

        Data2->SetTitle("");
        Data2->SetLineWidth(2);
        Data2->SetMarkerStyle(20);
        Data2->SetMarkerSize(1.1);
        Data2->SetLineColor(1);
        Data2->Draw("option same 0pz");

        line_.Draw("same");
        line_up.Draw("same");
        line_down.Draw("same");


       // stringstream out_file_root;
       // out_file_root << DIR_1 << plot_.long_name << ".root";

        tc->Print(out_file.str().data());
//        tc->Print(out_file_image.str().data());
//        tc->SaveAs(out_file_root.str().data());
        tc->Close();

        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///

        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///


        delete asym_errors;
        delete asym_errors2;
        delete Data;
        delete Data2;
        delete hist_stack;
        delete hist_sm;

        delete[] data_ratio;
        delete[] data_err_up_ratio;
        delete[] data_err_down_ratio;

        delete[] nominal_yield;
        delete[] x_bin_center;
        delete[] x_bin_down;
        delete[] x_bin_up;
        delete[] y_bin_down;
        delete[] y_bin_up;

        delete[] nominal_one;
        delete[] value_zero;
        delete[] y_bin_down_ratio;
        delete[] y_bin_up_ratio;

        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            delete hist_central[t_];
        }
        delete[] hist_central;

        for (int b_ = 0; b_ < num_bins; b_++) {
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                delete[] systematic_up[b_][t_];
            }
            delete[] systematic_up[b_];
        }
        delete[] systematic_up;

        for (int b_ = 0; b_ < num_bins; b_++) {
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                delete[] systematic_down[b_][t_];
            }
            delete[] systematic_down[b_];
        }
        delete[] systematic_down;


        delete tl_1;
        delete tl_2;
        delete tc;

        ttemp->Clear();
    }

    delete ttemp;

    cout << "Done." << endl;
    return 0;
}






// // // VR_Etmiss_Zmm
// // // Plot: MT2
// // plot_.root_get_member = "MT2";
// // plot_.root_get_factor = "";
// // plot_.signal_region = "VR_Etmiss_Zmm";
// // plot_.request_cuts = "";
// // plot_.blind_data = "";
// // plot_.x_label = "m_{T2} [GeV]";
// // plot_.y_label = "Events / 10 GeV";
// // plot_.x_bin_width = 10.0;
// // plot_.x_range_min = 0.0;
// // plot_.x_range_max = 300.0;
// // plot_.y_range_min = 1.0e-1;
// // plot_.y_range_max = 1.0e7;
// // plot_.long_name = "VR_Etmiss_Zmm_MT2";
// // vect_plot.push_back(plot_);
