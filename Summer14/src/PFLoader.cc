#include "../include/PFLoader.hh"
#include "Dummy/Puppi/interface/RecoObj.hh"
#include "../include/puppiContainer.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include <fastjet/tools/GridMedianBackgroundEstimator.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "TMath.h"
#include "boost/shared_ptr.hpp"

using namespace baconhep;

PFLoader::PFLoader(TTree *iTree,std::string iName) { 
  fPFCands  = new TClonesArray("baconhep::TPFPart");
  fPFCandBr = 0; 
  iTree->SetBranchAddress("PFPart",       &fPFCands, &fPFCandBr);
  //fPFCandBr  = iTree->GetBranch("PFPart");
  boost::shared_ptr<edm::ParameterSet> lConfig = edm::readPSetsFrom(iName.c_str());
  //const edm::ParameterSet& lConfig = edm::readPSetsFrom(iName.c_str())->getParameter<edm::ParameterSet>("puppi"); 
  //edm::ParameterSet lConfig1 = *lConfig;
  const edm::ParameterSet& lConfig1 = lConfig->getParameter<edm::ParameterSet>("puppi"); 
  fPuppiContainer = new PuppiContainer(lConfig1);
  fPuppi.clear(); 
}
PFLoader::~PFLoader() { 
  delete fPFCands;
  delete fPFCandBr;
}
void PFLoader::reset() { 
  fPt     = 0; 
  fEta    = 0; 
  fPhi    = 0; 
  fEcalE  = 0; 
  fHcalE  = 0; 
  fPFType = 0; 
  fPuppi.clear();

  fSumEt     = 0; 
  fMet       = 0; 
  fMetPhi    = 0; 
  fU1        = 0; 
  fU2        = 0; 
  fPupSumEt  = 0; 
  fPupMet    = 0; 
  fPupMetPhi = 0; 
  fPupU1     = 0; 
  fPupU2     = 0;   
  fCHSSumEt  = 0; 
  fCHSMet    = 0; 
  fCHSMetPhi = 0; 
  fCHSU1     = 0; 
  fCHSU2     = 0;   
}
void PFLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("sumet"        ,&fSumEt     ,"fSumEt/F");
  fTree->Branch("met"          ,&fMet       ,"fMet/F");
  fTree->Branch("metphi"       ,&fMetPhi    ,"fMetPhi/F");
  fTree->Branch("u1"           ,&fU1        ,"fU1/F");
  fTree->Branch("u2"           ,&fU2        ,"fU2/F");
  fTree->Branch("pupisumet"    ,&fPupSumEt  ,"fPupSumEt/F");
  fTree->Branch("puppimet"     ,&fPupMet    ,"fPupMet/F");
  fTree->Branch("puppimetphi"  ,&fPupMetPhi ,"fPupMetPhi/F");
  fTree->Branch("puppiu1"      ,&fPupU1     ,"fPupU1/F");
  fTree->Branch("puppiu2"      ,&fPupU2     ,"fPupU2/F");
  fTree->Branch("chssumet"     ,&fCHSSumEt  ,"fCHSSumEt/F");
  fTree->Branch("chsmet"       ,&fCHSMet    ,"fCHSMet/F");
  fTree->Branch("chsmetphi"    ,&fCHSMetPhi ,"fCHSMetPhi/F");
  fTree->Branch("chsu1"        ,&fCHSU1     ,"fCHSU1/F");
  fTree->Branch("chsu2"        ,&fCHSU2     ,"fCHSU2/F");
  //fTree->Branch("pt"     ,&fPt    ,"fPt/F");
  //fTree->Branch("eta"    ,&fEta   ,"fEta/F");
  //fTree->Branch("phi"    ,&fPhi   ,"fPhi/F");
  //fTree->Branch("ecalE"  ,&fEcalE ,"fEcalE/F");
  //fTree->Branch("hcalE"  ,&fHcalE ,"fHcalE/F");
  //fTree->Branch("pfType" ,&fPFType,"fPFType/F");
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
void PFLoader::load(int iEvent,TLorentzVector &iVec) { 
  fPFCands  ->Clear();
  fPFCandBr ->GetEntry(iEvent);
  fetch(iVec);
}
void PFLoader::fetch(TLorentzVector &iVec) { 
  fAllParticles  .resize(0);
  fPFParticles   .resize(0);
  fPFCHSParticles.resize(0);
  std::vector<RecoObj> lPuppi;

  TLorentzVector lVec   (0,0,0,0);
  TLorentzVector lCHSVec(0,0,0,0);
  fSumEt    = 0;
  fCHSSumEt = 0; 
  for(int i0 = 0; i0 < fPFCands->GetEntriesFast(); i0++) { 
    TPFPart  *pPart = (TPFPart*)((*fPFCands)[i0]);    
    //if(pPart->eta > 0 ) cout << "---> " << pPart->pt << " -- " << pPart->eta << " -- " << pPart->phi << " -- " << pPart->pfType << " -- " << pPart->hcalE << " -- " << pPart->ecalE   << endl;
    //if(pPart->pfType == 5) continue;
    RecoObj   pObj  = convert(pPart);
    PseudoJet pPar  = convert(&pObj);
    fAllParticles  .push_back(pObj); 
    fPFParticles   .push_back(pPar);
    if(pPar.user_index() != 3) fPFCHSParticles.push_back(pPar);
    TLorentzVector pVec(0,0,0,0);
    pVec.SetPtEtaPhiM(pPart->pt,0.,pPart->phi,0.);
    lVec    -= pVec;
    if(pPar.user_index() != 3) lCHSVec -= pVec;
    fSumEt    += pPart->pt;
    if(pPar.user_index() != 3) fCHSSumEt += pPart->pt;
    fPt     = pPart->pt;
    fEta    = pPart->eta;
    fPhi    = pPart->phi;
    fEcalE  = pPart->ecalE;
    fHcalE  = pPart->hcalE;
    fPFType = float(pPart->pfType);
    //fTree->Fill();
  }
  fMet       = lVec.Pt();
  fMetPhi    = lVec.Phi();
  lVec += iVec;
  lVec.RotateZ(-iVec.Phi());
  fU1        = lVec.Px();
  fU2        = lVec.Py();

  fCHSMet    = lCHSVec.Pt();
  fCHSMetPhi = lCHSVec.Phi();
  lCHSVec += iVec;
  lCHSVec.RotateZ(-iVec.Phi());
  fCHSU1     = lCHSVec.Px();
  fCHSU2     = lCHSVec.Py();
}
std::vector<fastjet::PseudoJet> PFLoader::puppiFetch(TLorentzVector &iVec) {
  //puppiContainer curEvent(fAllParticlesPuppi);
  //fPuppi = curEvent.puppiFetch(7,0.5);
  fPuppiContainer->initialize(fAllParticles);
  fPuppiContainer->puppiWeights();
  fPuppi = fPuppiContainer->puppiParticles();
  TLorentzVector lVec   (0,0,0,0);
  fPupSumEt    = 0;
  for(unsigned int i0 = 0; i0 < fPuppi.size(); i0++) { 
    TLorentzVector pVec;
    pVec.SetPtEtaPhiM(fPuppi[i0].pt(),0.,fPuppi[i0].phi(),0.);
    lVec -= pVec;
    fPupSumEt  += fPuppi[i0].pt();
  }
  fPupMet    = lVec.Pt();
  fPupMetPhi = lVec.Phi();
  lVec      += iVec;
  lVec.RotateZ(-iVec.Phi());
  fPupU1     = lVec.Px();
  fPupU2     = lVec.Py();
  return fPuppi;
}
std::vector<fastjet::PseudoJet> PFLoader::pfchsFetch(double iPtCut) { 
  if(iPtCut < 0) return fPFCHSParticles;
  std::vector<PseudoJet> lParts;
  for(unsigned int i0 = 0; i0 < fPFCHSParticles.size(); i0++) { 
    int charge_tmp = fPFCHSParticles[i0].user_index() > 1;
    if(charge_tmp && fPFCHSParticles[i0].pt() < iPtCut) continue;
    lParts.push_back(fPFCHSParticles[i0]);
  }
  return lParts;
}
RecoObj PFLoader::convert(TPFPart *iPart) { 
    bool lIsCh   = (iPart->pfType == 1 || iPart->pfType == 2 || iPart->pfType == 3) && (iPart->vtxId > -1 || fabs(iPart->dz) < 0.3) ;
    bool lIsPV   = (iPart->vtxId  == 0 || (fabs(iPart->dz) < 0.3 && lIsCh));
    int lID = -1;
    // if(fabs(iPart->eta) > 2.5) lIsCh = false;
    if (!lIsCh) lID = 1;
    if (lIsCh &&  lIsPV) lID = 2;
    if (lIsCh && !lIsPV) lID = 3;
    RecoObj pJet;
    pJet.pt      = iPart->pt;
    pJet.eta     = iPart->eta;
    pJet.phi     = iPart->phi;
    pJet.m       = iPart->m;
    pJet.id      = lID;
    pJet.vtxId   = iPart->vtxId;
    pJet.trkChi2 = iPart->trkChi2;
    pJet.vtxChi2 = iPart->vtxChi2;
    pJet.pfType  = iPart->pfType;
    pJet.depth   = iPart->depth;
    pJet.time    = iPart->time;
    pJet.d0      = iPart->d0;
    pJet.dZ      = iPart->dz;
    //if(!lIsCh && pJet.pfType == 1) pJet.pfType  = 5;//iPart->d0;
    // if(!lIsCh) pJet.d0      = 0;//iPart->d0;
    //if(!lIsCh) pJet.dZ      = 0;//iPart->dz;
    return pJet;
}
PseudoJet PFLoader::convert(RecoObj *iObj) { 
  double Px    = iObj->pt*cos(iObj->phi);
  double Py    = iObj->pt*sin(iObj->phi);
  double theta = 2*atan(exp( -iObj->eta)); //eta = -ln(tan(theta/2))
  double Pz    = iObj->pt/tan(theta);
  double E     = sqrt(Px*Px+Py*Py+Pz*Pz+iObj->m*iObj->m);
  fastjet::PseudoJet tmp_psjet(Px, Py, Pz, E);
  tmp_psjet.set_user_index(iObj->id);
  return tmp_psjet;
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
