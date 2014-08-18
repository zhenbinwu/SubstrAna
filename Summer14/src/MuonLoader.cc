#include "../include/MuonLoader.hh"
#include <iostream>
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
  fPt1   = 0; 
  fEta1  = 0; 
  fPhi1  = 0; 
  fPt2   = 0; 
  fEta2  = 0; 
  fPhi2  = 0; 
  fZPt   = 0; 
  fZY    = 0; 
  fZM    = 0; 
}
void MuonLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("pt_z"  ,&fZPt ,"fZPt/F");
  fTree->Branch("y_z" ,&fZY,"fZY/F");
  fTree->Branch("m_z" ,&fZM,"fZM/F");
  fTree->Branch("pt_1"  ,&fPt1 ,"fPt1/F");
  fTree->Branch("eta_1" ,&fEta1,"fEta1/F");
  fTree->Branch("phi_1" ,&fPhi1,"fPhi1/F");
  fTree->Branch("pt_2"  ,&fPt2 ,"fPt2/F");
  fTree->Branch("eta_2" ,&fEta2,"fEta2/F");
  fTree->Branch("phi_2" ,&fPhi2,"fPhi2/F");
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
  fPt1  = lMuon->pt;
  fEta1 = lMuon->eta;
  fPhi1 = lMuon->phi;
  return true;
}
bool MuonLoader::selectZ(std::vector<TLorentzVector> &iVetoes) {
  reset(); 
  int lCount = 0; 
  for  (int i0 = 0; i0 < fMuons->GetEntriesFast(); i0++) { 
    TMuon *pMuon = (TMuon*)((*fMuons)[i0]);
    if(!passLoose(pMuon)) continue;
    lCount++;
    bool lPassZ = false;
    for  (int i1 = i0+1; i1 < fMuons->GetEntriesFast(); i1++) { 
      TMuon *pMuon1 = (TMuon*)((*fMuons)[i1]);
      if(!passLoose(pMuon1)) continue;
      TLorentzVector lVec ; lVec .SetPtEtaPhiM(pMuon ->pt,pMuon ->eta,pMuon ->phi,0.105);
      TLorentzVector lVec1; lVec1.SetPtEtaPhiM(pMuon1->pt,pMuon1->eta,pMuon1->phi,0.105);
      if((lVec+lVec1).M() < 60 || (lVec+lVec1).M() > 120) continue; 
      iVetoes.push_back(lVec);
      iVetoes.push_back(lVec1);
      lPassZ = true;
      break;
    }
    if(lPassZ) break;
  }
  if(iVetoes.size() < 2) return false;
  TLorentzVector lZ; lZ = iVetoes[0] + iVetoes[1];
  fPt1  = iVetoes[0].Pt();
  fEta1 = iVetoes[0].Eta();
  fPhi1 = iVetoes[0].Phi();
  fPt2  = iVetoes[1].Pt();
  fEta2 = iVetoes[1].Eta();
  fPhi2 = iVetoes[1].Phi();
  fZPt  = lZ.Pt();
  fZY   = lZ.Rapidity();
  fZM   = lZ.M();
  fZEta = lZ.Eta();
  fZPhi = lZ.Phi();
  
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
  //if(!(muon->typeBits     & kGlobal || muon->typeBits & kTracker))  return false;
  //if(!(muon->selectorBits & kAllArbitrated))                        return false;
  if(!(muon->typeBits     & kPFMuon))                               return false;
  //if(fabs(muon->dz)> 1.0)                                           return false;
  double chargedIso = muon->chHadIso04;
  //double neutralIso = TMath::Max(muon->gammaIso04 + muon->neuHadIso04 - 0.5 * muon->puIso04, 0.0);
  double totalIso   = chargedIso;//+neutralIso;
  if(totalIso/muon->pt > 0.25) return false;
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
  lMuon.SetPtEtaPhiM(fPt1,fEta1,fPhi1,0.105);
  return lMuon;
}
TLorentzVector MuonLoader::boson() { 
  TLorentzVector lMuon; 
  lMuon.SetPtEtaPhiM(fZPt,fZEta,fZPhi,fZM);
  return lMuon;
}
