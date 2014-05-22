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
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "../include/puppiContainer.hh"

#include "BaconAna/DataFormats/interface/TPFPart.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

using namespace std;
using namespace fastjet;
using namespace baconhep;

PseudoJet convert(TPFPart *iPart) { 
  double Px    = iPart->pt*cos(iPart->phi);
  double Py    = iPart->pt*sin(iPart->phi);
  double theta = 2*atan(exp(-iPart->eta)); //eta = -ln(tan(theta/2))
  double Pz    = iPart->pt/tan(theta);
  double E     = iPart->e;
  bool lIsCh   = (iPart->pfType == 1 || iPart->pfType == 2 || iPart->pfType == 3) && (iPart->vtxId > -1 || fabs(iPart->dz) < 0.2) ;
  bool lIsPV   = (iPart->vtxId  == 0 || (fabs(iPart->dz) < 0.2 && lIsCh));
  int lID = -1;
  if (!lIsCh) lID = 1;
  if (lIsCh &&  lIsPV) lID = 2;
  if (lIsCh && !lIsPV) lID = 3;
  fastjet::PseudoJet tmp_psjet(Px, Py, Pz, E);
  tmp_psjet.set_user_index(lID);
  return tmp_psjet;
}
PseudoJet convert(TGenParticle *iPart) { 
  double Px    = iPart->pt*cos(iPart->phi);
  double Py    = iPart->pt*sin(iPart->phi);
  double theta = 2*atan(exp(-iPart->eta)); //eta = -ln(tan(theta/2))
  double Pz    = iPart->pt/tan(theta);
  double E     = sqrt(Px*Px+Py*Py+Pz*Pz+iPart->mass*iPart->mass);
  fastjet::PseudoJet tmp_psjet(Px, Py, Pz, E);
  return tmp_psjet;
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
  iJetUnc->setJetPt ( iJet.pt()  );
  iJetUnc->setJetEta( iJet.eta() );
  double jetunc = iJetUnc->getUncertainty(true);
  return jetunc;
}
int main(){
  TFile* fIn =  TFile::Open("root://eoscms.cern.ch//store/group/phys_jetmet/ntran/PUPPI/miniSamples/62x/rsgww1000_62x_PU40BX50/ntuple_1_1_VQC.root");
  TTree* tree = (TTree*) fIn->Get("Events");
  TClonesArray *fPFPart  = new TClonesArray("baconhep::TPFPart");
  TClonesArray *fGenPart = new TClonesArray("baconhep::TGenParticle");
  tree->SetBranchAddress("GenParticle",  &fGenPart);
  tree->SetBranchAddress("PFPart",       &fPFPart);
  //Setup the JetEnergy Corrections
  std::string cmsenv = "/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_7_patch2/src/";
  std::vector<JetCorrectorParameters> corrParams;
  corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L1FastJet_AK5PF.txt"));
  corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L2Relative_AK5PF.txt"));
  corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L3Absolute_AK5PF.txt"));
  //corrParams.push_back(JetCorrectorParameter(cmsenv+'BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt'));
  JetCorrectorParameters     param(cmsenv+"BaconProd/Utils/data/Summer13_V1_DATA_Uncertainty_AK5PF.txt");
  FactorizedJetCorrector   *jetCorr = new FactorizedJetCorrector(corrParams);
  JetCorrectionUncertainty *jetUnc  = new JetCorrectionUncertainty(param);
  
  //Setup the Collections
  std::vector<PseudoJet> particles;
  std::vector<PseudoJet> genparticles;  
  for(int i0 = 0; i0 < tree->GetEntriesFast(); i0++) { 
    if (i0 > 10) break;
      tree->GetEntry(i0);
      particles   .clear();
      genparticles.clear();
        
      std::cout << "i0 = " << i0 << ", and N PF candidates = " << fPFPart->GetEntriesFast() << std::endl;
      //First Lets get the PF Candidates
      for( int i1 = 0; i1 < fPFPart->GetEntriesFast(); i1++){//9,entries loop,fill the vector particles with PF particles
        baconhep::TPFPart *pPartTmp = (baconhep::TPFPart*)((*fPFPart)[i1]);
	//Convert from Particle Flow to PseudoJet with id
	fastjet::PseudoJet pFastJet = convert(pPartTmp);
	//Build the collection 
        particles.push_back(pFastJet);
      }
      //Second Lets get the Gen Particles
      for( int i1 = 0; i1 < fGenPart->GetEntriesFast(); i1++){//9,entries loop,fill the vector particles with PF particles
        baconhep::TGenParticle *pPartTmp = (baconhep::TGenParticle*)((*fGenPart)[i1]);
	if(pPartTmp->status != 1) continue;
	//Convert gen particle to PseudoJet
	fastjet::PseudoJet pFastJet = convert(pPartTmp);
	//Build the collection 
        genparticles.push_back(pFastJet);
      }
      //Lets set the particles collections in the Puppi container
      puppiContainer curEvent(particles); 
      curEvent.setGen(genparticles); 
      
      //Lets fetch the particle collections 
      vector<PseudoJet> gen_event       = curEvent.genFetch();
      vector<PseudoJet> puppi_event     = curEvent.puppiFetch(40.);
      vector<PseudoJet> pf_event        = curEvent.pfFetch();
      vector<PseudoJet> chs_event       = curEvent.pfchsFetch(-1);
      //Now lets define a jet collection we would like to cluster
      JetDefinition jet_def(kt_algorithm, 0.4);
      AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
      //With this set of candidates we can cluster
     
      ClusterSequenceArea pGen    (gen_event    ,jet_def,area_def);
      ClusterSequenceArea pPup    (puppi_event  ,jet_def,area_def);
      ClusterSequenceArea pPF     (pf_event     ,jet_def,area_def);
      ClusterSequenceArea pCHS    (chs_event    ,jet_def,area_def);
      //Now lets define a selector to get the leading jet
      Selector selector = SelectorNHardest(1); 
      //And lets select the leading jet from the above collections
      vector<PseudoJet> genJets     = selector(sorted_by_pt(pGen    .inclusive_jets()));
      vector<PseudoJet> pupJets     = selector(sorted_by_pt(pPup    .inclusive_jets()));
      vector<PseudoJet> pfJets      = selector(sorted_by_pt(pPF     .inclusive_jets()));
      vector<PseudoJet> chsJets     = selector(sorted_by_pt(pCHS    .inclusive_jets()));
      //Now we've got jets so lets calculate rho for the two that need it (CHS/PF)
      AreaDefinition area_def_for_rho(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
      JetDefinition jet_def_for_rho(kt_algorithm, 0.6);
      Selector rho_range    =  SelectorAbsRapMax(5.0);
      //Selector rho_rangeCHS =  SelectorAbsRapMax(2.5);
      ClusterSequenceArea clust_seq_rho   (pf_event , jet_def_for_rho, area_def_for_rho);
      ClusterSequenceArea clust_seq_rhoCHS(chs_event, jet_def_for_rho, area_def_for_rho);
      // the two background estimators
      JetMedianBackgroundEstimator rho   (rho_range, clust_seq_rho);
      JetMedianBackgroundEstimator rhoCHS(rho_range, clust_seq_rhoCHS);
      //Finally we can run the JEC on the first jet
      double genPt   = genJets[0].pt() * 1.;
      double pupPt   = pupJets[0].pt() * correction(pupJets[0],jetCorr,0);         //Assume rho is zero for this
      double pfPt    = pfJets [0].pt() * correction(pfJets [0],jetCorr,rho.rho()); 
      double chsPt   = chsJets[0].pt() * correction(chsJets[0],jetCorr,rhoCHS.rho()); 
      double pfPtUnc = unc(pfJets[0],jetUnc);
      //Lets cout the jet pts 
      cout << "===> Gen Pt : " << genPt <<  "  -- Puppi Pt : " << pupPt << " -- pf Pt : " << pfPt << " +/- " << pfPtUnc << " -- chs Pt : " << chsPt << endl;
  }
}




