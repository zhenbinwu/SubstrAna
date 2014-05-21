#include "../include/MuonLoader.hh"
#include "TMath.h"

using namespace baconhep;

MuonLoader::MuonLoader(TTree *iTree) { 
  fMuons  = new TClonesArray("baconhep::TMuon");
  iTree->SetBranchAddress("Muon",       &fMuons);
  fMuonBr  = iTree->GetBranch("Muon");
}
MuonLoader::~MuonLoader() { 
  delete fMuons;
  delete fMuonBr;
}
void MuonLoader::reset() { 
  fPt   = 0; 
  fEta  = 0; 
  fPhi  = 0; 
}
void MuonLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("pt_1"  ,&fPt ,"fPt/F");
  fTree->Branch("eta_1" ,&fEta,"fEta/F");
  fTree->Branch("phi_1" ,&fPhi,"fPhi/F");
}
void MuonLoader::load(int iEvent) { 
  fMuons   ->Clear();
  fMuonBr ->GetEntry(iEvent);
}
bool MuonLoader::selectSingleMu() {
  reset(); 
  TMuon *lMuon = 0; 
  int lCount = 0; 
  for  (int i0 = 0; i0 < fMuons->GetEntriesFast(); i0++) { 
    TMuon *pMuon = (TMuon*)((*fMuons)[i0]);
    if(!passLoose(pMuon)) continue;
    lMuon = pMuon;
    lCount++;
    if(lCount > 1) return false;
  }
  if(lMuon == 0) return false;
  fPt  = lMuon->pt;
  fEta = lMuon->eta;
  fPhi = lMuon->phi;
  return true;
}
bool MuonLoader::vetoMu() {
  for  (int i0 = 0; i0 < fMuons->GetEntriesFast(); i0++) { 
    TMuon *pMuon = (TMuon*)((*fMuons)[i0]);
    if(passLoose(pMuon)) return true;
  }
  return false;
}
//H=>ZZ Mu Id
bool MuonLoader::passLoose(TMuon *muon) { 
  if(!(muon->typeBits     & kGlobal || muon->typeBits & kTracker))  return false;
  if(!(muon->selectorBits & kAllArbitrated))                        return false;
  if(!(muon->typeBits     & kPFMuon))                               return false;
  if(fabs(muon->dz)> 1.0)                                           return false;
  double chargedIso = muon->chHadIso04;
  double neutralIso = TMath::Max(muon->gammaIso04 + muon->neuHadIso04 - 0.5 * muon->puIso04, 0.0);
  double totalIso   = chargedIso+neutralIso;
  if(totalIso/muon->pt > 0.4) return false;
  return true;
}
bool MuonLoader::passTight(TMuon *iMuon) { 
  if(!(iMuon->typeBits & kGlobal))  return false;
  if(fabs(iMuon->dz)  > 0.2)        return false;
  if(fabs(iMuon->d0)  > 0.045)      return false;
  if(iMuon->muNchi2        > 10)    return false;
  if(iMuon->nValidHits     < 1)     return false;
  if(iMuon->nMatchStn      < 2)     return false;
  if(iMuon->nPixHits       < 1)     return false;
  if(iMuon->nTkLayers      < 6)     return false;
  if(!(iMuon->typeBits & kPFMuon))  return false;


  double chargedIso = iMuon->chHadIso04;
  double neutralIso = TMath::Max(iMuon->gammaIso04 + iMuon->neuHadIso04 - 0.5 * iMuon->puIso04, 0.0);
  double totalIso   = chargedIso+neutralIso;
  if(totalIso/iMuon->pt > 0.15) return false;
  return true;
}
TLorentzVector MuonLoader::muon() { 
  TLorentzVector lMuon; 
  lMuon.SetPtEtaPhiM(fPt,fEta,fPhi,0.105);
  return lMuon;
}
