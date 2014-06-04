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
void PFLoader::load(int iEvent) { 
  fPFCands  ->Clear();
  fPFCandBr ->GetEntry(iEvent);
  fetch();
}
void PFLoader::fetch() { 
  fAllParticles  .resize(0);
  fPFParticles   .resize(0);
  fPFCHSParticles.resize(0);
  std::vector<RecoObj> lPuppi;
  for(int i0 = 0; i0 < fPFCands->GetEntriesFast(); i0++) { 
    TPFPart  *pPart = (TPFPart*)((*fPFCands)[i0]);    
    RecoObj   pObj  = convert(pPart);
    PseudoJet pPar  = convert(&pObj);
    fAllParticles  .push_back(pObj); 
    fPFParticles   .push_back(pPar);
    if(pPar.user_index() != 3) fPFCHSParticles.push_back(pPar);
  }
}
std::vector<fastjet::PseudoJet> PFLoader::puppiFetch() {
  //puppiContainer curEvent(fAllParticlesPuppi);
  //fPuppi = curEvent.puppiFetch(7,0.5);
  fPuppiContainer->initialize(fAllParticles);
  fPuppiContainer->puppiWeights();
  fPuppi = fPuppiContainer->puppiParticles();
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
    bool lIsCh   = (iPart->pfType == 1 || iPart->pfType == 2 || iPart->pfType == 3) && (iPart->vtxId > -1 || fabs(iPart->dz) < 0.2) ;
    bool lIsPV   = (iPart->vtxId  == 0 || (fabs(iPart->dz) < 0.2 && lIsCh));
    int lID = -1;
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
