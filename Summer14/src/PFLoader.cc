#include "../include/PFLoader.hh"
#include "../include/RecoObj.hh"
#include "../include/puppiContainer.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include <fastjet/tools/GridMedianBackgroundEstimator.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
//#include :
#include "TMath.h"

using namespace baconhep;

PFLoader::PFLoader(TTree *iTree) { 
  fPFCands  = new TClonesArray("baconhep::TPFPart");
  iTree->SetBranchAddress("PFPart",       &fPFCands);
  fPFCandBr  = iTree->GetBranch("PFPart");
  fPuppi.clear(); 
}
PFLoader::~PFLoader() { 
  delete fPFCands;
  delete fPFCandBr;
}
void PFLoader::reset() { 
  fMet     = 0; 
  fMetPhi  = 0; 
  fPuppi.clear();
}
void PFLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("met"    ,&fMet   ,"fMet/F");
  fTree->Branch("metphi" ,&fMetPhi,"fMetPhi/F");
}
void PFLoader::load(int iEvent) { 
  fPFCands  ->Clear();
  fPFCandBr ->GetEntry(iEvent);
}
TLorentzVector PFLoader::met() { 
  TLorentzVector lVec(0,0,0,0);
  for(int i0 = 0; i0 < fPFCands->GetEntriesFast(); i0++) { 
    TPFPart *pPart = (TPFPart*)((*fPFCands)[i0]);    
    TLorentzVector pVec(0,0,0,0);
    pVec.SetPtEtaPhiM(pPart->pt,0.,pPart->phi,0.);
    lVec += pVec;
  }
  return lVec;
}
std::vector<fastjet::PseudoJet> PFLoader::puppiFetch() {
  bool lIsCh = false;
  bool lIsPV = false;
  std::vector<RecoObj> lPuppi;
  for(int i0 = 0; i0 < fPFCands->GetEntriesFast(); i0++) { 
    TPFPart *pPart = (TPFPart*)((*fPFCands)[i0]);    
    lIsCh   = (pPart->pfType == 1 || pPart->pfType == 2 || pPart->pfType == 3) && (pPart->vtxId > -1 || fabs(pPart->dz) < 0.2) ;
    lIsPV   = (pPart->vtxId  == 0  || (fabs(pPart->dz) < 0.2 && lIsCh));
    int lID = -1;
    if (!lIsCh) lID = 1;
    if (lIsCh &&  lIsPV) lID = 2;
    if (lIsCh && !lIsPV) lID = 3;
    RecoObj pJet;
    pJet.pt  = pPart->pt;
    pJet.eta = pPart->eta;
    pJet.phi = pPart->phi;
    pJet.m   = pPart->m;
    pJet.id  = lID;
    pJet.vtxId = pPart->vtxId;
    pJet.trkChi2 = pPart->trkChi2;
    pJet.vtxChi2 = pPart->vtxChi2;
    pJet.pfType  = pPart->pfType;
    pJet.depth   = pPart->depth;
    pJet.time    = pPart->time;
    pJet.d0        = pPart->d0;
    pJet.dZ        = pPart->dz;
    lPuppi.push_back(pJet); 
  }
  puppiContainer curEvent(lPuppi);
  fPuppi = curEvent.puppiFetch(7,0.5);
  return fPuppi;
}
std::vector<fastjet::PseudoJet> PFLoader::puppiJets(std::vector<TLorentzVector> iVetoes) { 
  std::vector < fastjet::PseudoJet > lJets;
  getJets(fPuppi,lJets,iVetoes);
  return lJets;
}
std::vector<fastjet::PseudoJet> PFLoader::pfJets(std::vector<TLorentzVector> iVetoes) { 
  std::vector < fastjet::PseudoJet > lConstits;
  std::vector < fastjet::PseudoJet > lJets;
  for(int i0 = 0; i0 < fPFCands->GetEntriesFast(); i0++) { 
    TPFPart *pPart = (TPFPart*)((*fPFCands)[i0]);    
    fastjet::PseudoJet pJet;
    pJet.reset_PtYPhiM(pPart->pt,pPart->eta,pPart->phi,pPart->m);
    lConstits.push_back(pJet);
  }
  getJets(lConstits,lJets,iVetoes);
  return lJets;
}
std::vector<fastjet::PseudoJet> PFLoader::chsJets(std::vector<TLorentzVector> iVetoes) { 
  std::vector < fastjet::PseudoJet > lConstits;
  std::vector < fastjet::PseudoJet > lJets;
  for(int i0 = 0; i0 < fPFCands->GetEntriesFast(); i0++) { 
    TPFPart *pPart = (TPFPart*)((*fPFCands)[i0]);    
    if(pPart->vtxId > 0) continue;
    fastjet::PseudoJet pJet;
    pJet.reset_PtYPhiM(pPart->pt,pPart->eta,pPart->phi,pPart->m);
    lConstits.push_back(pJet);
  }
  getJets(lConstits,lJets,iVetoes);
  return lJets;
}
void PFLoader::getJets(std::vector < fastjet::PseudoJet > &constits,std::vector < fastjet::PseudoJet > &jets,std::vector<TLorentzVector> iVetoes) { 
  double rParam = 0.7;
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);    
  int activeAreaRepeats = 1;
  double ghostArea = 0.01;
  double ghostEtaMax = 7.0;
  fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
  fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
  fastjet::ClusterSequenceArea* thisClustering_ = new fastjet::ClusterSequenceArea(constits, jetDef, fjAreaDefinition);
  std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering_->inclusive_jets(5.0));
  for(unsigned int i0 = 0; i0 < out_jets.size(); i0++) {
    bool pMatch = false;
    for(unsigned int i1 = 0; i1 < iVetoes.size(); i1++) {
      double pDEta = out_jets[i0].eta() - iVetoes[i1].Eta();
      double pDPhi = fabs(out_jets[i0].phi() - iVetoes[i1].Phi());
      if(fabs(pDPhi) > 2.*TMath::Pi()-fabs(pDPhi)) pDPhi =  2.*TMath::Pi()-fabs(pDPhi);
      if(sqrt(pDPhi*pDPhi+pDEta*pDEta) < 0.5) continue;
      pMatch = true;
    } 
    if(pMatch) continue;
    jets.push_back(out_jets[i0]); 
  }
}
