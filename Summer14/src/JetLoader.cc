#include "../include/JetLoader.hh"
#include <iostream>
#include "TMath.h"

using namespace baconhep;

JetLoader::JetLoader(TTree *iTree) { 
  fJets  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress("Jet04",       &fJets);
  fJetBr  = iTree->GetBranch("Jet04");
}
JetLoader::~JetLoader() { 
  delete fJets;
  delete fJetBr;
}
void JetLoader::load(int iEvent) { 
  fJets   ->Clear();
  fJetBr ->GetEntry(iEvent);
}
std::vector<fastjet::PseudoJet> JetLoader::genjets() { 
  std::vector<fastjet::PseudoJet> lVec;
  for  (int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    TJet *pJet = (TJet*)((*fJets)[i0]);
    if(pJet->genpt < 5) continue;
    TLorentzVector pLVec; pLVec.SetPtEtaPhiM(pJet->genpt,pJet->geneta,pJet->genphi,pJet->genm);
    fastjet::PseudoJet pVec(pLVec.Px(),pLVec.Py(),pLVec.Pz(),pLVec.E());
    lVec.push_back(pVec);
  }
  return lVec;
}
