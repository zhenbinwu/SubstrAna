#include "../include/GenLoader.hh"
#include "../include/MuonLoader.hh"
#include "../include/PFLoader.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/JetCleanser.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/SafeSubtractor.hh"
#include "fastjet/contrib/SoftKiller.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/Selector.hh"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>

using namespace std;
using namespace fastjet;
using namespace contrib;

//Object Processors
GenLoader       *fGen      = 0; 
MuonLoader      *fMuon     = 0; 
PFLoader        *fPFCand   = 0; 

TClonesArray *fJet;
TBranch      *fJetBr;

TTree* load(std::string iName) { 
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  return lTree;
}


struct JetInfo {
  int npu ;

  vector<float> pt;
  vector<float> ptcorr;
  vector<float> ptraw;
  vector<float> ptclean;
  vector<float> pttrim;
  vector<float> pttrimsafe;
  vector<float> ptconst;
  vector<float> ptunc;
  vector<float> eta;
  vector<float> phi;
  vector<float> m;
  vector<float> mraw;
  vector<float> mclean;
  vector<float> mtrim;
  vector<float> mtrimsafe;
  vector<float> mconst;
  vector<int>   nparticles;
  vector<int>   nneutrals;
  vector<int>   ncharged;

  // gen level info
  vector<float> ptgen;
  vector<float> etagen;
  vector<float> phigen;
  vector<float> mgen;
  vector<float> mrawgen;
  vector<float> mtrimgen;
  vector<float> mtrimsafegen;
  vector<float> mcleangen; //needed?
  vector<float> mconstgen;//needed?
  vector<int>   imatch;
};


void getConstitsForCleansing(vector<PseudoJet> inputs, vector<PseudoJet> &oNeutrals, vector<PseudoJet> &oChargedLV, vector<PseudoJet> &oChargedPU){
  for (unsigned int i = 0; i < inputs.size(); i++){
    if (inputs[i].user_index() <= 1) oNeutrals.push_back(inputs[i]);
    if (inputs[i].user_index() == 2) oChargedLV.push_back(inputs[i]);
    if (inputs[i].user_index() == 3) oChargedPU.push_back(inputs[i]);
  }
}


class SW_IsPupCharged : public SelectorWorker {
public:
  SW_IsPupCharged(){}
  virtual bool pass(const PseudoJet & jet) const {
    return (jet.user_index() > 1);
  }
};

Selector SelectorIsPupCharged(){
  return Selector(new SW_IsPupCharged());
}


class SW_IsPupVertex : public SelectorWorker {
public:
  SW_IsPupVertex(){}
  virtual bool pass(const PseudoJet & jet) const {
    return (jet.user_index() == 2);
  }
};


Selector SelectorIsPupVertex(){
  return Selector(new SW_IsPupVertex());
}


double correction( PseudoJet &iJet,FactorizedJetCorrector *iJetCorr,double iRho) { 
  iJetCorr->setJetPt (iJet.pt());
  iJetCorr->setJetEta(iJet.eta());
  iJetCorr->setJetPhi(iJet.phi());
  iJetCorr->setJetE  (iJet.e());
  iJetCorr->setJetA  (iJet.area());
  iJetCorr->setRho(iRho);
  iJetCorr->setJetEMF(-99.0);
  double jetcorr= iJetCorr->getCorrection();
  return jetcorr;
}


double unc( PseudoJet &iJet,JetCorrectionUncertainty *iJetUnc) { 
  if(fabs(iJet.eta()) > 5. || fabs(iJet.pt()) < 10.) return 1.;
  iJetUnc->setJetPt ( iJet.pt()  );
  iJetUnc->setJetEta( iJet.eta() );
  double jetunc = iJetUnc->getUncertainty(true);
  return jetunc;
}


int matchingIndex(PseudoJet jet, vector<PseudoJet> genjets) {
  float rmin = 9999.;  
  int imatch = -1;
  for(unsigned int i = 0; i < genjets.size(); i++) {
    float rtemp = jet.delta_R(genjets[i]);
    if ( rtemp > 0.3 ) continue;
    if ( rtemp < rmin ){
      rmin =  rtemp;
      imatch = i;
    }
  }
  return (imatch);  
}



int matchingIndexFromJetInfo(PseudoJet jet, JetInfo jetInfo) {
  float rmin = 9999.;  
  int imatch = -1;
  for(unsigned int i = 0; i < (jetInfo.eta).size(); i++) {
    double dEta = fabs(jet.eta() - (jetInfo.eta)[i]);
    double dPhi = fabs(jet.phi() - (jetInfo.phi)[i]);
    if(dPhi > 2.*TMath::Pi()-dPhi) dPhi =  2.*TMath::Pi()-dPhi;
    float rtemp = sqrt(dEta*dEta+dPhi*dPhi);
    if ( rtemp > 0.3 ) continue;
    if ( rtemp < rmin ){
      rmin =  rtemp;
      imatch = i;
    }
  }
  return (imatch);  
}




void setupTree(TTree *iTree, JetInfo &iJet, std::string iName) {
  iTree->Branch((iName+"npu"       ).c_str(),&iJet.npu       );

  iTree->Branch((iName+"pt"        ).c_str(),&iJet.pt        );
  iTree->Branch((iName+"ptcorr"    ).c_str(),&iJet.ptcorr    );
  iTree->Branch((iName+"ptraw"     ).c_str(),&iJet.ptraw     );
  iTree->Branch((iName+"ptclean"   ).c_str(),&iJet.ptclean   );
  iTree->Branch((iName+"pttrim"    ).c_str(),&iJet.pttrim    );
  iTree->Branch((iName+"pttrimsafe").c_str(),&iJet.pttrimsafe);
  iTree->Branch((iName+"ptconst"   ).c_str(),&iJet.ptconst   );
  iTree->Branch((iName+"ptunc"     ).c_str(),&iJet.ptunc     );
  iTree->Branch((iName+"eta"       ).c_str(),&iJet.eta       );
  iTree->Branch((iName+"phi"       ).c_str(),&iJet.phi       );
  iTree->Branch((iName+"m"         ).c_str(),&iJet.m         );
  iTree->Branch((iName+"mraw"      ).c_str(),&iJet.mraw      );
  iTree->Branch((iName+"mtrim"     ).c_str(),&iJet.mtrim     );
  iTree->Branch((iName+"mtrimsafe" ).c_str(),&iJet.mtrimsafe );
  iTree->Branch((iName+"mclean"    ).c_str(),&iJet.mclean    );
  iTree->Branch((iName+"mconst"    ).c_str(),&iJet.mconst    );
  iTree->Branch((iName+"nparticles").c_str(),&iJet.nparticles);
  iTree->Branch((iName+"nneutrals" ).c_str(),&iJet.nneutrals);
  iTree->Branch((iName+"ncharged"  ).c_str(),&iJet.ncharged);
  // gen info
  iTree->Branch((iName+"ptgen"       ).c_str(),&iJet.ptgen       );
  iTree->Branch((iName+"etagen"      ).c_str(),&iJet.etagen      );
  iTree->Branch((iName+"phigen"      ).c_str(),&iJet.phigen      );
  iTree->Branch((iName+"mgen"        ).c_str(),&iJet.mgen        );
  iTree->Branch((iName+"mrawgen"     ).c_str(),&iJet.mrawgen     );//needed?
  iTree->Branch((iName+"mtrimgen"    ).c_str(),&iJet.mtrimgen    );//needed?
  iTree->Branch((iName+"mtrimsafegen").c_str(),&iJet.mtrimsafegen);//needed?
  iTree->Branch((iName+"mcleangen"   ).c_str(),&iJet.mcleangen   );//needed?
  iTree->Branch((iName+"mconstgen"   ).c_str(),&iJet.mconstgen   );//needed?

  iTree->Branch((iName+"imatch"      ).c_str(),&iJet.imatch      );
}


void clear(JetInfo &iJet) {
  iJet.npu  = -1;

  iJet.pt         .clear();
  iJet.ptcorr     .clear();
  iJet.ptraw      .clear();
  iJet.ptclean    .clear();
  iJet.pttrim     .clear();
  iJet.pttrimsafe .clear();
  iJet.eta        .clear();
  iJet.phi        .clear();
  iJet.m          .clear();
  iJet.mraw       .clear();
  iJet.mtrim      .clear();
  iJet.mtrimsafe  .clear();
  iJet.mclean     .clear();
  iJet.mconst     .clear();
  iJet.nparticles .clear();
  iJet.nneutrals  .clear();
  iJet.ncharged   .clear();

  iJet.ptgen       .clear();
  iJet.etagen      .clear();
  iJet.phigen      .clear();
  iJet.mgen        .clear();
  iJet.mrawgen     .clear();
  iJet.mtrimgen    .clear();
  iJet.mtrimsafegen.clear();
  iJet.mcleangen   .clear();
  iJet.mconstgen   .clear();
  iJet.imatch      .clear();
}



void setJet(PseudoJet &iJet, JetInfo &iJetI,JetMedianBackgroundEstimator bge_rho, JetMedianBackgroundEstimator bge_rhom, JetMedianBackgroundEstimator bge_rhoC, 
	    bool isGEN, bool isCHS, FactorizedJetCorrector *iJetCorr, JetCorrectionUncertainty *iJetUnc, JetCleanser &gsn_cleanser, 
	    bool doGenMatching, vector<PseudoJet> genJets) {

  // -- area-median subtractor  ( safe area subtractor )
  contrib::SafeAreaSubtractor *area_subtractor = 0;
  if(!isCHS) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom);
  if( isCHS) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom,SelectorIsPupCharged(),SelectorIsPupVertex());
  PseudoJet lCorr =  (*area_subtractor)(iJet);
  
  // -- constituent subtractor
  contrib::ConstituentSubtractor *const_subtractor = 0;
  const_subtractor = new contrib::ConstituentSubtractor(&bge_rhoC);
  (*const_subtractor).use_common_bge_for_rho_and_rhom(true);
  PseudoJet lConstit = (*const_subtractor)(iJet);

  // -- cleansing 
  vector<PseudoJet> neutrals,chargedLV,chargedPU;
  getConstitsForCleansing(iJet.constituents(),neutrals,chargedLV,chargedPU);
  PseudoJet     lClean = gsn_cleanser(neutrals,chargedLV,chargedPU);
  
  // -- trimming
  fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05)));
  PseudoJet lTrim     = (trimmer)(iJet);
  trimmer.set_subtractor(area_subtractor);
  PseudoJet lTrimSafe = (trimmer)(iJet);
    
  // -- apply the JEC
  double lJEC = 1.;
  double lUnc = 0 ;
  if (!isGEN){
    lJEC = correction(iJet,iJetCorr,bge_rho.rho());  
    lUnc = unc       (iJet,iJetUnc);
  }

  // -- find the gen jet matched to this reco jet
  int imatch = -1;
  if (doGenMatching) imatch = matchingIndex(iJet,genJets);
  
  // -- fill jet info
  (iJetI.pt        ).push_back(lCorr     .pt());
  (iJetI.ptcorr    ).push_back(iJet      .pt()*lJEC);
  (iJetI.ptraw     ).push_back(iJet      .pt());
  (iJetI.ptclean   ).push_back(lClean    .pt());
  (iJetI.pttrim    ).push_back(lTrim     .pt());
  (iJetI.pttrimsafe).push_back(lTrimSafe .pt());
  (iJetI.ptconst   ).push_back(lConstit  .pt());
  (iJetI.ptunc     ).push_back(lUnc);
  (iJetI.eta       ).push_back(iJet      .eta());
  (iJetI.phi       ).push_back(iJet      .phi());
  (iJetI.mraw      ).push_back(iJet      .m());
  (iJetI.m         ).push_back(lCorr     .m());
  (iJetI.mclean    ).push_back(lClean    .m());
  (iJetI.mtrim     ).push_back(lTrim     .m());
  (iJetI.mtrimsafe ).push_back(lTrimSafe .m());
  (iJetI.mconst    ).push_back(lConstit  .m());
  (iJetI.nparticles).push_back((iJet.constituents()).size());
  (iJetI.nneutrals ).push_back(neutrals.size());
  (iJetI.ncharged  ).push_back(chargedLV.size()+chargedPU.size());
  
  if (imatch > -1){
    (iJetI.imatch  ).push_back(imatch);
    (iJetI.ptgen    ).push_back(genJets[imatch].pt());
    (iJetI.etagen   ).push_back(genJets[imatch].eta());
    (iJetI.phigen   ).push_back(genJets[imatch].phi());
    (iJetI.mgen     ).push_back(genJets[imatch].m());
  }
  else {
    (iJetI.imatch   ).push_back(imatch);
    (iJetI.ptgen    ).push_back(-999.);
    (iJetI.etagen   ).push_back(-999.);
    (iJetI.phigen   ).push_back(-999.);
    (iJetI.mgen     ).push_back(-999.);
  }
  
}


void setRecoJet(PseudoJet &iJet, JetInfo &iJetI, JetInfo iGenJetI, JetMedianBackgroundEstimator bge_rho, JetMedianBackgroundEstimator bge_rhom, JetMedianBackgroundEstimator bge_rhoC, 
		bool isCHS, FactorizedJetCorrector *iJetCorr, JetCorrectionUncertainty *iJetUnc, JetCleanser &gsn_cleanser) {

  // -- area-median subtractor  ( safe area subtractor )
  contrib::SafeAreaSubtractor *area_subtractor = 0;
  if(!isCHS) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom);
  if( isCHS) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom,SelectorIsPupCharged(),SelectorIsPupVertex());
  PseudoJet lCorr =  (*area_subtractor)(iJet);
  
  // -- constituent subtractor
  contrib::ConstituentSubtractor *const_subtractor = 0;
  const_subtractor = new contrib::ConstituentSubtractor(&bge_rhoC);
  (*const_subtractor).use_common_bge_for_rho_and_rhom(true);
  PseudoJet lConstit = (*const_subtractor)(iJet);

  // -- cleansing 
  vector<PseudoJet> neutrals,chargedLV,chargedPU;
  getConstitsForCleansing(iJet.constituents(),neutrals,chargedLV,chargedPU);
  PseudoJet     lClean = gsn_cleanser(neutrals,chargedLV,chargedPU);
  
  // -- trimming
  fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05)));
  PseudoJet lTrim     = (trimmer)(iJet);
  trimmer.set_subtractor(area_subtractor);
  PseudoJet lTrimSafe = (trimmer)(iJet);
    
  // -- apply the JEC
  double lJEC = correction(iJet,iJetCorr,bge_rho.rho());  
  double lUnc = unc       (iJet,iJetUnc);
  

  // -- find the gen jet matched to this reco jet
  //int imatch = matchingIndex(iJet,genJets);
  int imatch = matchingIndexFromJetInfo(iJet, iGenJetI);
  
  // -- fill jet info
  (iJetI.pt        ).push_back(lCorr     .pt());
  (iJetI.ptcorr    ).push_back(iJet      .pt()*lJEC);
  (iJetI.ptraw     ).push_back(iJet      .pt());
  (iJetI.ptclean   ).push_back(lClean    .pt());
  (iJetI.pttrim    ).push_back(lTrim     .pt());
  (iJetI.pttrimsafe).push_back(lTrimSafe .pt());
  (iJetI.ptconst   ).push_back(lConstit  .pt());
  (iJetI.ptunc     ).push_back(lUnc);
  (iJetI.eta       ).push_back(iJet      .eta());
  (iJetI.phi       ).push_back(iJet      .phi());
  (iJetI.mraw      ).push_back(iJet      .m());
  (iJetI.m         ).push_back(lCorr     .m());
  (iJetI.mclean    ).push_back(lClean    .m());
  (iJetI.mtrim     ).push_back(lTrim     .m());
  (iJetI.mtrimsafe ).push_back(lTrimSafe .m());
  (iJetI.mconst    ).push_back(lConstit  .m());
  (iJetI.nparticles).push_back((iJet.constituents()).size());
  (iJetI.nneutrals ).push_back(neutrals.size());
  (iJetI.ncharged  ).push_back(chargedLV.size()+chargedPU.size());
  
  if (imatch > -1){
    (iJetI.imatch).push_back(imatch);
    (iJetI.ptgen    ).push_back((iGenJetI.pt)[imatch]);
    (iJetI.etagen   ).push_back((iGenJetI.eta)[imatch]);
    (iJetI.phigen   ).push_back((iGenJetI.phi)[imatch]);
    
    (iJetI.mgen     ).push_back((iGenJetI.m)[imatch]);
    (iJetI.mrawgen     ).push_back((iGenJetI.mraw)[imatch]);
    (iJetI.mtrimgen    ).push_back((iGenJetI.mtrim)[imatch]);
    (iJetI.mtrimsafegen).push_back((iGenJetI.mtrimsafe)[imatch]);
    (iJetI.mcleangen   ).push_back((iGenJetI.mclean)[imatch]);
    (iJetI.mconstgen   ).push_back((iGenJetI.mconst)[imatch]);
  }
  else { //dummy values ?
    (iJetI.imatch).push_back(imatch);
    (iJetI.ptgen    ).push_back(-999.);
    (iJetI.etagen   ).push_back(-999.);
    (iJetI.phigen   ).push_back(-999.);
    (iJetI.mgen     ).push_back(-999.);
    (iJetI.mrawgen     ).push_back(-999.);
    (iJetI.mtrimgen    ).push_back(-999.);
    (iJetI.mtrimsafegen).push_back(-999.);
    (iJetI.mcleangen   ).push_back(-999.);
    (iJetI.mconstgen   ).push_back(-999.);
  }
  
}


void setGenJet(PseudoJet &iJet, JetInfo &iJetI,JetMedianBackgroundEstimator bge_rho, JetMedianBackgroundEstimator bge_rhom, JetMedianBackgroundEstimator bge_rhoC, JetCleanser &gsn_cleanser) {

  // -- area-median subtractor  ( safe area subtractor )
  contrib::SafeAreaSubtractor *area_subtractor = 0;
  area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom);
  PseudoJet lCorr =  (*area_subtractor)(iJet);
  
  // -- constituent subtractor
  contrib::ConstituentSubtractor *const_subtractor = 0;
  const_subtractor = new contrib::ConstituentSubtractor(&bge_rhoC);
  (*const_subtractor).use_common_bge_for_rho_and_rhom(true);
  PseudoJet lConstit = (*const_subtractor)(iJet);

  // -- cleansing 
  vector<PseudoJet> neutrals,chargedLV,chargedPU;
  getConstitsForCleansing(iJet.constituents(),neutrals,chargedLV,chargedPU);
  PseudoJet     lClean = gsn_cleanser(neutrals,chargedLV,chargedPU);
  
  // -- trimming
  fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05)));
  PseudoJet lTrim     = (trimmer)(iJet);
  trimmer.set_subtractor(area_subtractor);
  PseudoJet lTrimSafe = (trimmer)(iJet);
    
  // -- fill jet info
  (iJetI.pt        ).push_back(lCorr     .pt());
  (iJetI.ptcorr    ).push_back(iJet      .pt());
  (iJetI.ptraw     ).push_back(iJet      .pt());
  (iJetI.ptclean   ).push_back(lClean    .pt());
  (iJetI.pttrim    ).push_back(lTrim     .pt());
  (iJetI.pttrimsafe).push_back(lTrimSafe .pt());
  (iJetI.ptconst   ).push_back(lConstit  .pt());
  (iJetI.ptunc     ).push_back(0.);
  (iJetI.eta       ).push_back(iJet      .eta());
  (iJetI.phi       ).push_back(iJet      .phi());
  (iJetI.mraw      ).push_back(iJet      .m());
  (iJetI.m         ).push_back(lCorr     .m());
  (iJetI.mclean    ).push_back(lClean    .m());
  (iJetI.mtrim     ).push_back(lTrim     .m());
  (iJetI.mtrimsafe ).push_back(lTrimSafe .m());
  (iJetI.mconst    ).push_back(lConstit  .m());
  (iJetI.nparticles).push_back((iJet.constituents()).size());
  (iJetI.nneutrals ).push_back(neutrals.size());
  (iJetI.ncharged  ).push_back(chargedLV.size()+chargedPU.size());
    
  // some duplicate infos.... should be fixed.
  (iJetI.imatch      ).push_back((iJetI.pt).size()-1);
  (iJetI.ptgen       ).push_back(lCorr     .pt());
  (iJetI.etagen      ).push_back(iJet      .eta());
  (iJetI.phigen      ).push_back(iJet      .phi());
  (iJetI.mgen        ).push_back(lCorr     .m());
  (iJetI.mrawgen     ).push_back(iJet      .m());
  (iJetI.mtrimgen    ).push_back(lTrim     .m());
  (iJetI.mtrimsafegen).push_back(lTrimSafe .m());
  (iJetI.mcleangen   ).push_back(lClean    .m());
  (iJetI.mconstgen   ).push_back(lConstit  .m());
}


// ------------------------------------------------------------------------------------------
void fillGenJetsInfo(vector<PseudoJet> &iJets, vector<PseudoJet> &iParticles, JetInfo &iJetInfo, JetCleanser &gsn_cleanser, int nPU ){

  // -- Compute rho, rho_m for SafeAreaSubtraction
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  JetDefinition jet_def_for_rho(kt_algorithm, 0.4);
  Selector rho_range =  SelectorAbsRapMax(5.0);
  ClusterSequenceArea clust_seq_rho(iParticles, jet_def_for_rho, area_def);
  JetMedianBackgroundEstimator bge_rho(rho_range, clust_seq_rho);
  JetMedianBackgroundEstimator bge_rhom(rho_range, clust_seq_rho);
  BackgroundJetPtMDensity m_density;
  bge_rhom.set_jet_density_class(&m_density);
  
  //// use GridMedianBackgroundEstimator (faster), however doesn't work with SafeAreaSubtraction
  //GridMedianBackgroundEstimator bge_rho(5.0,0.8);
  //bge_rho.set_particles(iParticles);
  //GridMedianBackgroundEstimator bge_rhom(5.0,0.8);
  //bge_rhom.set_particles(iParticles);
  
  // -- Background estimator for constituents subtractor
  JetMedianBackgroundEstimator bge_rhoC(rho_range,jet_def_for_rho, area_def);
  BackgroundJetScalarPtDensity *scalarPtDensity = new BackgroundJetScalarPtDensity();
  bge_rhoC.set_jet_density_class(scalarPtDensity);
  bge_rhoC.set_particles(iParticles);
    
  // -- Clear jet info for each event                                                                                                                                                                                                  
  clear(iJetInfo);
  
  iJetInfo.npu = nPU;

  // -- Loop over jets in the event and set jets variables                                                                                                                                                                      
  for (unsigned int j = 0; j < iJets.size(); j++){
    setGenJet( iJets[j], iJetInfo, bge_rho, bge_rhom, bge_rhoC, gsn_cleanser);
    //cout << iTree.GetName() << "  " << (iJetInfo.pt)[j] << "  "<< (iJetInfo.ptcorr)[j] <<endl;                                                                                                                                      
  }

}
// ------------------------------------------------------------------------------------------


void fillRecoJetsInfo(vector<PseudoJet> &iJets, vector<PseudoJet> &iParticles, JetInfo &iJetInfo, JetInfo iGenJetInfo, bool isCHS, FactorizedJetCorrector *jetCorr, JetCorrectionUncertainty *ijetUnc,
		      JetCleanser &gsn_cleanser, int nPU ){
  
  // -- Compute rho, rho_m for SafeAreaSubtraction
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  JetDefinition jet_def_for_rho(kt_algorithm, 0.4);
  Selector rho_range =  SelectorAbsRapMax(5.0);
  ClusterSequenceArea clust_seq_rho(iParticles, jet_def_for_rho, area_def);
  JetMedianBackgroundEstimator bge_rho(rho_range, clust_seq_rho);
  JetMedianBackgroundEstimator bge_rhom(rho_range, clust_seq_rho);
  BackgroundJetPtMDensity m_density;
  bge_rhom.set_jet_density_class(&m_density);
  
  //// use GridMedianBackgroundEstimator (faster), however doesn't work with SafeAreaSubtraction
  //GridMedianBackgroundEstimator bge_rho(5.0,0.8);
  //bge_rho.set_particles(iParticles);
  //GridMedianBackgroundEstimator bge_rhom(5.0,0.8);
  //bge_rhom.set_particles(iParticles);
  
  // -- Background estimator for constituents subtractor
  JetMedianBackgroundEstimator bge_rhoC(rho_range,jet_def_for_rho, area_def);
  BackgroundJetScalarPtDensity *scalarPtDensity = new BackgroundJetScalarPtDensity();
  bge_rhoC.set_jet_density_class(scalarPtDensity);
  bge_rhoC.set_particles(iParticles);
    
  // -- Clear jet info for each event                                                                                                                                                                                                  
  clear(iJetInfo);
  
  iJetInfo.npu = nPU;

  // -- Loop over jets in the event and set jets variables                                                                                                                                                                      
  for (unsigned int j = 0; j < iJets.size(); j++){
    setRecoJet( iJets[j], iJetInfo, iGenJetInfo, bge_rho, bge_rhom, bge_rhoC, isCHS, jetCorr, ijetUnc, gsn_cleanser);
    //cout << iTree.GetName() << "  " << (iJetInfo.pt)[j] << "  "<< (iJetInfo.ptcorr)[j] <<endl;                                                                                                                                      
  }

}

// ------------------------------------------------------------------------------------------
void setupCMSSWJetReadOut(TTree *iTree, float R ) {
  
  cout << "Setting up to read jet collection : " << Form("Jet0%d",int(R*10)) << endl;
  fJet  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress(Form("Jet0%d",int(R*10)), &fJet);
  fJetBr  = iTree->GetBranch(Form("Jet0%d",int(R*10)));
}
// ------------------------------------------------------------------------------------------


// ------------------------------------------------------------------------------------------
void readCMSSWJet(int entry, TTree *iTree, TTree &oTree,  std::vector<fastjet::PseudoJet> genJets, JetInfo &iJetI) {


  // -- Clear jet info for each event
  clear(iJetI);

  // -- Read event and fill jet info
  iTree->GetEntry(entry);

  for (int i = 0; i < fJet->GetEntriesFast(); i++){
    TJet *pJet = (TJet*)((*fJet)[i]);
    
    // -- fill jet info                                                                                                                                                                                                       
    (iJetI.pt        ).push_back(pJet->pt);
    (iJetI.ptcorr    ).push_back(pJet->pt);
    (iJetI.ptraw     ).push_back(pJet->ptRaw);
    (iJetI.eta       ).push_back(pJet->eta);
    (iJetI.phi       ).push_back(pJet->phi);
    (iJetI.m         ).push_back(pJet->mass);
    (iJetI.nparticles).push_back(pJet->nParticles);
    (iJetI.nneutrals ).push_back(pJet->nNeutrals);
    (iJetI.ncharged  ).push_back(pJet->nCharged);
    

    // for now fill this branches with dummy value
    (iJetI.ptclean   ).push_back(-999.);
    (iJetI.pttrim    ).push_back(-999.);
    (iJetI.pttrimsafe).push_back(-999.);
    (iJetI.ptconst   ).push_back(-999.);
    (iJetI.ptunc     ).push_back(-999.);
    (iJetI.eta       ).push_back(-999.);
    (iJetI.phi       ).push_back(-999.);
    (iJetI.mraw      ).push_back(-999.);
    (iJetI.mclean    ).push_back(-999.);
    (iJetI.mtrim     ).push_back(-999.);
    (iJetI.mtrimsafe ).push_back(-999.);
    (iJetI.mconst    ).push_back(-999.);
  
    //-- gen matching
    int imatch = -1;
    double mindr = 0.3;
    TLorentzVector *recojet = new TLorentzVector();
    recojet->SetPtEtaPhiM(pJet->pt, pJet->eta, pJet->phi, pJet->mass);
    for (unsigned int ig = 0; ig < genJets.size(); ig++){
      TLorentzVector *genjet = new TLorentzVector();
      genjet->SetPtEtaPhiM(genJets[ig].pt(), genJets[ig].eta(), genJets[ig].phi(), genJets[ig].m());
      double dr = recojet->DeltaR(*genjet);
      if (dr < mindr){
	mindr = dr;
	imatch = ig;
      }
      delete genjet;
    }
    
    delete recojet;

    if (imatch > -1){
      (iJetI.ptgen    ).push_back(genJets[imatch].pt());
      (iJetI.etagen   ).push_back(genJets[imatch].eta());
      (iJetI.phigen   ).push_back(genJets[imatch].phi());
      (iJetI.mgen     ).push_back(genJets[imatch].m());
      (iJetI.imatch   ).push_back(imatch);
    }
    else {
      (iJetI.ptgen    ).push_back(-999.);
      (iJetI.etagen   ).push_back(-999.);
      (iJetI.phigen   ).push_back(-999.);
      (iJetI.mgen     ).push_back(-999.);
      (iJetI.imatch   ).push_back(imatch);
    }
  }

  // --- fill tree 
  oTree.Fill();
  
}
// ------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------
bool FillChain(TChain& chain, const std::string& inputFileList)
{
  std::ifstream inFile(inputFileList.c_str());
  std::string buffer;

  if(!inFile.is_open())
    {
      std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
      return false;
    }
  
  while(1)
    {
      inFile >> buffer;
      if(!inFile.good()) break;
      chain.Add(buffer.c_str());
    }

  return true;
}
// ------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------------
//--- MAIN PROGRAM
//---------------------------------------------------------------------------------------------------------------
int main (int argc, char ** argv) {

  // args 
  std::string inputFilesList = argv[1];        // input file name
  int maxEvents              = atoi(argv[2]);  // max events
  std::string fOut           = argv[3];        // output name
  float jetR                 = atof(argv[4]);  // jet cone size      
  bool doCMSSWJets           = atoi(argv[5]);  // if want to analyze PF Jets from CMSSW

  // --- Read list of files to be analyzed and fill TChain
  TChain* lTree = new TChain("Events");
  FillChain(*lTree, inputFilesList);
  if (lTree->GetEntries() < maxEvents || maxEvents == -1) maxEvents = lTree->GetEntries(); 

  cout << "This analysis will run on "<< maxEvents << " events" <<endl; 

  // --- Load branches
  fPFCand = new PFLoader (lTree,"Puppi_cff.py");
  fGen    = new GenLoader(lTree);
  if (doCMSSWJets) setupCMSSWJetReadOut(lTree, jetR);

  TEventInfo *eventInfo = new TEventInfo();
  lTree->SetBranchAddress("Info",&eventInfo);

  // --- Setup JEC on the fly
  std::string cmsenv = "/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_7_patch2/src/";
  std::vector<JetCorrectorParameters> corrParams;
  corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L1FastJet_AK5PF.txt"));
  corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L2Relative_AK5PF.txt"));
  corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L3Absolute_AK5PF.txt"));
  //corrParams.push_back(JetCorrectorParameter(cmsenv+'BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt'));
  JetCorrectorParameters     param(cmsenv+"BaconProd/Utils/data/Summer13_V1_DATA_Uncertainty_AK5PF.txt");
  FactorizedJetCorrector   *jetCorr = new FactorizedJetCorrector(corrParams);
  JetCorrectionUncertainty *jetUnc  = new JetCorrectionUncertainty(param);

  // --- Setup JetAlgos
  JetDefinition jet_def(antikt_algorithm,jetR);         // the jet definition....
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  
  // --- Setup cleansing
  JetDefinition subjet_def(kt_algorithm,0.2);
  JetCleanser gsn_cleanser(subjet_def,JetCleanser::gaussian_cleansing,JetCleanser::input_nc_separate);
  gsn_cleanser.SetGaussianParameters(0.617,0.62,0.15,0.22);

  // --- Setup output trees
  TFile *fout = new TFile(fOut.c_str(),"RECREATE");

  TTree *genTree   = new TTree("gen"  , "gen"  );
  TTree *pfTree    = new TTree("pf"   , "pf"   );
  TTree *chsTree   = new TTree("chs"  , "chs"  );
  TTree *puppiTree = new TTree("puppi", "puppi");
  TTree *cmsswTree = new TTree("cmsswpf", "cmsswpf");

  JetInfo JGenInfo, JPFInfo, JCHSInfo, JPuppiInfo, JCMSSWPFInfo;    
  setupTree(genTree,   JGenInfo    , "" );
  setupTree(pfTree,    JPFInfo     , "" );
  setupTree(chsTree,   JCHSInfo    , "" );
  setupTree(puppiTree, JPuppiInfo  , "" );
  if (doCMSSWJets) setupTree(cmsswTree, JCMSSWPFInfo, "" );

  
  // --- start loop over events
  for(int ientry = 0; ientry < maxEvents; ientry++) { 

    if(ientry % 50 == 0) 
      std::cout << "===> Processed " << ientry << " - Done : " << (float(ientry)/float(maxEvents))*100 << "%" << std::endl;
    
    // -- For each event build collections of particles (gen, puppi, etc..) to cluster
    fPFCand->load(ientry);
    fGen   ->load(ientry); 
    vector<PseudoJet> gen_event       = fGen   ->genFetch();
    vector<PseudoJet> pf_event        = fPFCand->pfFetch();
    vector<PseudoJet> chs_event       = fPFCand->pfchsFetch(-1);
    vector<PseudoJet> puppi_event     = fPFCand->puppiFetch();

    // -- Cluster jets
    ClusterSequenceArea pGen    (gen_event    , jet_def, area_def);
    ClusterSequenceArea pPup    (puppi_event  , jet_def, area_def);
    ClusterSequenceArea pPF     (pf_event     , jet_def, area_def);
    ClusterSequenceArea pCHS    (chs_event    , jet_def, area_def);

    vector<PseudoJet> genJets     = sorted_by_pt(pGen    .inclusive_jets(25.));
    vector<PseudoJet> puppiJets   = sorted_by_pt(pPup    .inclusive_jets(25.));
    vector<PseudoJet> pfJets      = sorted_by_pt(pPF     .inclusive_jets(25.));
    vector<PseudoJet> chsJets     = sorted_by_pt(pCHS    .inclusive_jets(25.));

    lTree->GetEntry(ientry);
    int nPU = eventInfo->nPU;

    // save jet info in a tree
    fillGenJetsInfo(genJets, gen_event, JGenInfo, gsn_cleanser, nPU);
    fillRecoJetsInfo(puppiJets, puppi_event, JPuppiInfo, JGenInfo, false, jetCorr, jetUnc, gsn_cleanser,nPU);
    fillRecoJetsInfo(pfJets   , pf_event   , JPFInfo   , JGenInfo, false, jetCorr, jetUnc, gsn_cleanser,nPU);
    fillRecoJetsInfo(chsJets  , chs_event  , JCHSInfo  , JGenInfo, true , jetCorr, jetUnc, gsn_cleanser,nPU);

    genTree->Fill();
    puppiTree->Fill();
    pfTree->Fill();
    chsTree->Fill();

    if (doCMSSWJets)
      readCMSSWJet(ientry, lTree, *cmsswTree, genJets, JCMSSWPFInfo);
    
  }


  // --- Write trees 
  fout->cd();

  genTree  ->Write();
  pfTree   ->Write();
  chsTree  ->Write();
  puppiTree->Write();

  if (doCMSSWJets)  cmsswTree->Write();
}  
