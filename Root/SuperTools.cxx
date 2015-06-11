#include "Superflow/SuperTools.h"
#include "DileptonMatrixMethod/DileptonMatrixMethod.h"

#include <iostream>
#include <iomanip>
#include <cstddef>


namespace SuperTools {

std::fstream outMatrix;
std::string debugMatrixName = "debugMatrixMethod.txt";

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// from [1] https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/SUSYPhys/SUSYTools/trunk/Root/SUSYObjDef.cxx#L2833
//
double ttbar_powheg_differentialxsec(double ttbarpt)
{
    // Note: ttbarpt = (top + antitop).Pt() and expects Pt to be in MeV
    double weight = 1.0;
    
    if (ttbarpt/1000. < 40.)
        weight = (1./1.011850 + 1./0.994193)/2.;
    else if (ttbarpt/1000. < 170.)
         weight = (1./1.095920 + 1./1.034480)/2.;
    else if (ttbarpt/1000. < 340)
      weight = (1./1.407280 + 1./1.319110)/2.;
    else
      weight = (1./1.799380 + 1./1.710780)/2.;
 
    return weight;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// see https://github.com/gerbaudo/SusyntHlfv/blob/master/Root/MatrixPrediction.cxx

double getFakeWeight(sflow::Superlink* sl, std::string fakeRegion, susy::fake::Systematic::Value sys, bool computeSyst)
{
 //   std::string test = susy::fake::Systematic::str(sys);
 //   std::cout << "\t " << test << std::endl;
    
    if(!computeSyst){
        outMatrix.open((fakeRegion + "_" + debugMatrixName).c_str(), ios::app | ios::out);
    }
    double weight=1.0;

    unsigned int nVtx = sl->nt->evt()->nVtx;
    bool isMC         = sl->nt->evt()->isMC;
    const Lepton &l0  = *sl->leptons->at(0);
    const Lepton &l1  = *sl->leptons->at(1);
    float metRel      = 0.0; // dummy value 
    bool l0IsSig(sl->tools->isSignalLepton(&l0, *sl->baseElectrons, *sl->baseMuons, nVtx, isMC));
    bool l1IsSig(sl->tools->isSignalLepton(&l1, *sl->baseElectrons, *sl->baseMuons, nVtx, isMC));
    std::size_t iRegion = sl->fakeMatrix->getIndexRegion(fakeRegion);
    susy::fake::Lepton fl0(l0IsSig, l0.isEle(), l0.Pt(), l0.Eta());
    susy::fake::Lepton fl1(l1IsSig, l1.isEle(), l1.Pt(), l1.Eta());
    if(!computeSyst){        
        weight = sl->fakeMatrix->getTotalFake(fl0, fl1, iRegion, metRel, susy::fake::Systematic::SYS_NOM);
    }
    else{
    //    double nominal = sl->fakeMatrix->getTotalFake(fl0, fl1, iRegion, metRel, susy::fake::Systematic::SYS_NOM);
    //    if(nominal!=0){
    //        double in = 1.0/nominal;
    //        weight = (sl->fakeMatrix->getTotalFake(fl0, fl1, iRegion, metRel, sys))*in;
        weight = sl->fakeMatrix->getTotalFake(fl0, fl1, iRegion, metRel, sys);
        }
   // }//if(nominal)
 
    if(!computeSyst){
            int isElectron1 = l0.isEle() ? 1 : 0;
            int isElectron2 = l1.isEle() ? 1 : 0;
            int isTight1 = l0IsSig ? 1 : 0;
            int isTight2 = l1IsSig ? 1 : 0;
            double pt1 = l0.Pt();
            double pt2 = l1.Pt();
            double eta1 = l0.Eta();
            double eta2 = l1.Eta();
            outMatrix << isElectron1 << "\t"
                            << isElectron2 << "\t"
                            << isTight1    << "\t"
                            << isTight2    << "\t"
                            << pt1         << "\t"
                            << pt2         << "\t"
                            << eta1        << "\t"
                            << eta2        << "\n";
            outMatrix.close();
    }
    return weight;
}


} // close SuperTools namespace
