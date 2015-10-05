#include <TLorentzVector.h>
#include <TVector3.h>
#include <TCanvas.h>
#include <iostream>
#include "RestFrames/RestFrames.hh"


using namespace std;
using namespace RestFrames;

////////////////////////////////////////////////
// This is a script setting up the RestFrames 
// for the Super-Razor set-up where you group
// together two legs of invisible decay objects
// and only include the two leptons for each of
// the visible legs' mega-jets.
//
// daniel.joseph.antrim@cern.ch
// September 2015
////////////////////////////////////////////////


void testRestFrames() {

    /////////////////////////////
    // Set up the RestFrames
    /////////////////////////////
    LabRecoFrame LAB("lab", "LAB");
    DecayRecoFrame tt("SS", "SS");
    DecayRecoFrame t1("S1", "S_{1}");
    DecayRecoFrame t2("S2", "S_{2}");
    VisibleRecoFrame v1("v1", "vis_{1}");
    VisibleRecoFrame v2("v2", "vis_{2}");
    InvisibleRecoFrame i1("i1", "invis_{1}");
    InvisibleRecoFrame i2("i2", "invis_{2}");


    ////////////////////////////
    // Connect the frames
    ////////////////////////////
    LAB.SetChildFrame(tt);
    tt.AddChildFrame(t1);
    tt.AddChildFrame(t2);
    t1.AddChildFrame(v1);
    t1.AddChildFrame(i1);
    t2.AddChildFrame(v2);
    t2.AddChildFrame(i2);


    // check that the decay tree is connected correctly
    if(!LAB.InitializeTree()) {
        cout << "testRestFrames::InitializeTree ERROR    Unable to initialize tree from lab frame. Exitting." << endl;
        exit(1);
    }

    ////////////////////////////
    // Define the groups
    ////////////////////////////

    // invisible group for weakly interacting particles
    InvisibleGroup INV("inv", "WIMP Jigsaws");
    INV.AddFrame(i1);
    INV.AddFrame(i2);

    // combinatoric group (list of visible objects that end up in our hemispheres)
    CombinatoricGroup VIS("vis", "Visible object Jigsaws");
    VIS.AddFrame(v1);
    VIS.SetNElementsForFrame(v1, 1, false);
    VIS.AddFrame(v2);
    VIS.SetNElementsForFrame(v2, 1, false);

    ////////////////////////////
    // Setup the jigsaw rules
    ////////////////////////////
    SetMassInvJigsaw MinMassJigsaw("MinMass", "Invisible system mass Jigsaw");
    INV.AddJigsaw(MinMassJigsaw);

    SetRapidityInvJigsaw RapidityJigsaw("RapidityJigsaw", "Invisible system rapidity Jigsaw");
    INV.AddJigsaw(RapidityJigsaw);
    RapidityJigsaw.AddVisibleFrames(LAB.GetListVisibleFrames());

    ContraBoostInvJigsaw ContraBoostJigsaw("ContraBoostJigsaw", "ContraBoost Invariant Jigsaw");
    INV.AddJigsaw(ContraBoostJigsaw);
    ContraBoostJigsaw.AddVisibleFrames((t1.GetListVisibleFrames()), 0);
    ContraBoostJigsaw.AddVisibleFrames((t2.GetListVisibleFrames()), 1);
    ContraBoostJigsaw.AddInvisibleFrame(i1, 0);
    ContraBoostJigsaw.AddInvisibleFrame(i2, 1);

    MinMassesCombJigsaw HemiJigsaw("HemiJigsaw", "Minimize m_{Vis_{1,2}} Jigsaw");
    VIS.AddJigsaw(HemiJigsaw);
    HemiJigsaw.AddFrame(v1, 0);
    HemiJigsaw.AddFrame(v2, 1);


    // check that the Jigsaws are in place
    if(!LAB.InitializeAnalysis()) {
        cout << "testRestFrames::InitializeAnalysis ERROR    Unable to initialize analysis from lab-frame. Exitting." << endl;
        exit(1);
    }

    LAB.ClearEvent();


    //////////////////////////////
    // Draw some pictures of our decay tree
    //////////////////////////////
    TreePlot* treePlot = new TreePlot("tree", "Decay Tree");

    treePlot->SetFrameTree(LAB);
    treePlot->Draw("Decay Tree", "Decay Tree");

    treePlot->SetGroupTree(INV);
    treePlot->Draw("InvTree", "Invisible Jigsaws");

    treePlot->SetGroupTree(VIS);
    treePlot->Draw("VisTree", "Visible Jigsaws");
    

};

    
