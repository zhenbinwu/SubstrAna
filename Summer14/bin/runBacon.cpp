#include <iostream>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "fastjet/contrib/GenericSubtractor.hh"

#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/ClusterSequence.hh"
//#include "fastjet/contrib/JetCleanser.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/SafeSubtractor.hh"
#include "fastjet/contrib/SoftKiller.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TClonesArray.h"
//#include "puppiContainer.hh"

#include "BaconAna/DataFormats/interface/TPFPart.hh"

using namespace std;
using namespace fastjet;
using namespace baconhep;

int main(){

  std::cout << "hello world" << std::endl;

  TFile* fIn = new TFile("/eos/uscms/store/user/ntran/PUPPI/bacon/rsgww1000_62x_PU40BX50/ntuple_1_1_VQC.root");
  TTree* tree = (TTree*) fIn->Get("Events");
  TClonesArray *fPFPart = new TClonesArray("baconhep::TPFPart");
  tree->SetBranchAddress("PFPart",       &fPFPart);

  for(int i0 = 0; i0 < tree->GetEntriesFast(); i0++) { 
        
      if (i0 > 10) break;
      tree->GetEntry(i0);
        
      std::vector<fastjet::PseudoJet> particles;
      particles.clear();
        
      std::cout << "i0 = " << i0 << ", and N PF candidates = " << fPFPart->GetEntriesFast() << std::endl;
        
      for( int i1 = 0; i1 < fPFPart->GetEntriesFast(); i1++){//9,entries loop,fill the vector particles with PF particles

        baconhep::TPFPart *pPartTmp = (baconhep::TPFPart*)((*fPFPart)[i1]);

        double Px = pPartTmp->pt*cos(pPartTmp->phi);
        double Py = pPartTmp->pt*sin(pPartTmp->phi);
        double theta = 2*atan(exp(-pPartTmp->eta)); //eta = -ln(tan(theta/2))
        double Pz = pPartTmp->pt/tan(theta);
        double E  = pPartTmp->e;
        //double pdgId = pPartTmp->pfType;
        //int charge = pPartTmp->q;
        fastjet::PseudoJet tmp_psjet(Px, Py, Pz, E);
        particles.push_back(tmp_psjet);
  
    }//9,entries loop ,fill the vector particles with PFparticles
  }
}




