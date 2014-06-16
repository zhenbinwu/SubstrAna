///////
//
// runHats.cpp
// code to run on JMEtuples (Bacon) for HATS tutorial and Boost2014 studies
// usage: runHats  <Cone size>  <clustering algorithm (kt, ca, ak)>  <first event> <number of events>
// example: runHats 0.8 ca 0 1000
//
///////

// standard
#include <iostream>
#include <fstream> 

// Fastjet
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/GenericSubtractor.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

// Groomers
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/ModifiedMassDropTagger.hh"

// Pileup reduction
#include "fastjet/contrib/SafeSubtractor.hh"
#include "fastjet/contrib/SoftKiller.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"
//#include "fastjet/contrib/JetCleanser.hh"

// Taggers 
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "CMSTopTagger.hh"
#include <fastjet/tools/JHTopTagger.hh>
#include "QjetsPlugin.h"
#include "Qjets.h"
#include "fastjet/contrib/EnergyCorrelator.hh"
//#include "HEPTopTagger.hh"
//#include "HEPTopTaggerWrapper.hh"

// Jet Corrections
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//ROOT
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TClonesArray.h"

// Bacon (JME tuples)
#include "BaconAna/DataFormats/interface/TPFPart.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

// SoftKiller
#include "fastjet/contrib/SoftKiller.hh"

// PUPPI - make sure to include last
#include "../include/puppiContainer.hh"



using namespace std;
using namespace fastjet;
using namespace baconhep;
using namespace contrib;
using namespace fastjet::contrib;


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
// default CMS jet corrections
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
// CMS top tagger minmass calculation
double calculate_minmass(vector<PseudoJet> subjets)
{
  double mmin = -1;

  //minimum pariwise mass is only defined for jets with at least 3 subjets. 
  if ( subjets.size() >=3 )
  {
    // Select the leading 3 subjets
    Selector subjet_selector = SelectorNHardest(3); 
    vector<PseudoJet> leading3subjets      = subjet_selector(sorted_by_pt(subjets));
    
    PseudoJet sum01 = leading3subjets[0]+leading3subjets[1];
    PseudoJet sum02 = leading3subjets[0]+leading3subjets[2];
    PseudoJet sum12 = leading3subjets[1]+leading3subjets[2];
 
    vector<double> v_subjet_pair_masses;
    v_subjet_pair_masses.push_back(sum01.m());
    v_subjet_pair_masses.push_back(sum02.m());
    v_subjet_pair_masses.push_back(sum12.m());

    vector<double>::const_iterator smallestMassPair_iter = min_element(v_subjet_pair_masses.begin(),v_subjet_pair_masses.end());
    mmin  = *smallestMassPair_iter;
  }
  return mmin;
}

double deltaR(PseudoJet &iJet,PseudoJet &jJet) {
   
        double pEta = fabs(iJet.eta()-jJet.eta());
        double pPhi = fabs(iJet.phi()-jJet.phi());
        if(pPhi > 2.*TMath::Pi()-pPhi) pPhi =  2.*TMath::Pi()-pPhi;
        double deltaR = sqrt(pEta*pEta+pPhi*pPhi);
        return deltaR;
}
PseudoJet match(PseudoJet &iJet,vector<PseudoJet> &iJets) {
    for(unsigned int i0 = 0; i0 < iJets.size(); i0++) { 
        double pEta = fabs(iJet.eta()-iJets[i0].eta()); 
        double pPhi = fabs(iJet.phi() - iJets[i0].phi());
        if(pPhi > 2.*TMath::Pi()-pPhi) pPhi =  2.*TMath::Pi()-pPhi;
        if(sqrt(pEta*pEta+pPhi*pPhi) > 0.3) continue;
        return iJets[i0];
    }
    return PseudoJet();
}
// returns a vector containing: [0] uncorr mmdt mass [1] corr mmdt mass [2] uncorr mmdt filtered mass [3] corr mmdt filtered mass
vector<double> mmdt_corr_and_filter(PseudoJet &iJet, FactorizedJetCorrector *iJetCorr,double iRho)
{
  vector<double> vec;
  if (iJet!=0)
    {
      // if (verbose) cout<<"mmdt jet pt "<<mmdtjet.perp()<<" mass "<<mmdtjet.m()<<endl;
      // if (verbose) cout << "  delta_R between subjets: " << mmdtjet.structure_of<MMDT>().delta_R() << endl;
      // if (verbose) cout << "  symmetry measure(z):     " << mmdtjet.structure_of<MMDT>().symmetry() << endl;
      // if (verbose) cout << "  mass drop(mu):           " << mmdtjet.structure_of<MMDT>().mu() << endl;
      
      double mass_mmdt_corr    = iJet.m() * correction(iJet,iJetCorr,iRho); 

      // filter jet dynamically based on deltaR between subjets (arXiv:0802.2470)
      double dyn_Rfilt = min(0.3, iJet.structure_of<contrib::ModifiedMassDropTagger>().delta_R()*0.5);
      int    dyn_nfilt = 3;
      Filter filter(dyn_Rfilt, SelectorNHardest(dyn_nfilt));
      PseudoJet filtered_mmdtjet = filter(iJet);
      double mass_mmdt_filtered_corr    = filtered_mmdtjet.m() * correction(filtered_mmdtjet,iJetCorr,iRho); 
      

      vec.push_back( iJet.m() );
      vec.push_back( mass_mmdt_corr );
      vec.push_back( filtered_mmdtjet.m() );
      vec.push_back( mass_mmdt_filtered_corr );
    }
    // else {
    //   cout<<" mmdt_corr_and_filter failed to get jet"<<endl;
    // }
    return vec;
}

//Q jets stuff
// double getQjetVolatility(std::vector < fastjet::PseudoJet > constits, int QJetsN, int seed){
    
//     std::vector< float > qjetmasses;
    
//     double zcut(0.1), dcut_fctr(0.5), exp_min(0.), exp_max(0.), rigidity(0.1), truncationFactor(0.01);
    
//     for(unsigned int ii = 0 ; ii < (unsigned int) QJetsN ; ii++){
//         QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity, truncationFactor);
//         qjet_plugin.SetRandSeed(seed+ii); // new feature in Qjets to set the random seed
//         fastjet::JetDefinition qjet_def(&qjet_plugin);
//         fastjet::ClusterSequence qjet_seq(constits, qjet_def);
//         vector<fastjet::PseudoJet> inclusive_jets2 = sorted_by_pt(qjet_seq.inclusive_jets(5.0));
        
//         if (inclusive_jets2.size()>0) { qjetmasses.push_back( inclusive_jets2[0].m() ); }
        
//     }
    
//     // find RMS of a vector
//     float qjetsRMS = FindRMS( qjetmasses );
//     // find mean of a vector
//     float qjetsMean = FindMean( qjetmasses );
//     float qjetsVolatility = qjetsRMS/qjetsMean;
//     return qjetsVolatility;
// }
// float FindRMS( std::vector< float > qjetmasses ){
    
//     float total = 0.;
//     float ctr = 0.;
//     for (unsigned int i = 0; i < qjetmasses.size(); i++){
//         total = total + qjetmasses[i];
//         ctr++;
//     }
//     float mean = total/ctr;
    
//     float totalsquared = 0.;
//     for (unsigned int i = 0; i < qjetmasses.size(); i++){
//         totalsquared += (qjetmasses[i] - mean)*(qjetmasses[i] - mean) ;
//     }
//     float RMS = sqrt( totalsquared/ctr );
//     return RMS;
// }

// float FindMean( std::vector< float > qjetmasses ){
//     float total = 0.;
//     float ctr = 0.;
//     for (unsigned int i = 0; i < qjetmasses.size(); i++){
//         total = total + qjetmasses[i];
//         ctr++;
//     }
//     return total/ctr;
// }


int main( int argc, char *argv[] ){

  /////////////////////////////////////////////////////////////////////
  //  Get input parameters
  /////////////////////////////////////////////////////////////////////


  if (argc <= 5)
  {
    cout << "Usage: " << argv[0] << " <Cone size>  <clustering algorithm (kt, ca, ak)>  <first event> <number of events> <1 for ttbar 0 for qcd>"<< endl; 
    exit(1);
  }

  char *ConeSize = argv[1];
  std::string algo(argv[2]);
  char *FirstEventChar = argv[3];
  char *NeventsChar = argv[4];
  char *SampleTypeChar = argv[5];

  double R = atof (ConeSize);
  int FirstEvent = atoi (FirstEventChar);
  int Nevents = atoi (NeventsChar);
  int Sample = atoi (SampleTypeChar);  //top 1 qcd 0

  bool verbose = false;
  cout<<"Cluster jets with algorithm = "<<algo<<", R = "<<R<<endl;

  //int jetalgo = -1;
  fastjet::JetAlgorithm jetalgo;

  if (algo=="kt")         jetalgo = fastjet::kt_algorithm        ;         
  else if (algo=="ca")    jetalgo = fastjet::cambridge_algorithm ;
  else if (algo=="ak")    jetalgo = fastjet::antikt_algorithm    ;
  else if (algo=="KT")    jetalgo = fastjet::kt_algorithm        ;
  else if (algo=="CA")    jetalgo = fastjet::cambridge_algorithm ;
  else if (algo=="AK")    jetalgo = fastjet::antikt_algorithm    ;
  else if (algo=="kt_algorithm")        jetalgo = fastjet::kt_algorithm        ;
  else if (algo=="cambridge_algorithm") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo=="antikt_algorithm")    jetalgo = fastjet::antikt_algorithm    ;
  else if (algo=="0") jetalgo = fastjet::kt_algorithm        ;
  else if (algo=="1") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo=="2") jetalgo = fastjet::antikt_algorithm    ;
  else 
  {
    cout << "The entered jet algorithm name is not kt, ca or ak" << endl; 
    exit(1);
  }
  /////////////////////////////////////////////////////////////////////
  //  Outfile and ttree
  /////////////////////////////////////////////////////////////////////


  Float_t PFjetMassUncorr                    ;        
  Float_t CHSjetMassUncorr                   ;        
  Float_t PUPjetMassUncorr                   ;        
  Float_t GENjetMassUncorr                   ;        
  Float_t PFjetMass                          ;        
  Float_t CHSjetMass                         ;        
  Float_t PUPjetMass                         ;        
  Float_t GENjetMass                         ;        
  Float_t PFjetPt                            ;        
  Float_t CHSjetPt                           ;        
  Float_t PUPjetPt                           ;        
  Float_t GENjetPt                           ;        
  Float_t PFjetPtUnc                         ;        
  Float_t CHSjetPtUnc                        ;        
  Float_t PUPjetPtUnc                        ;        
  Float_t GENjetPtUnc                        ;        
  Float_t PFjetMassTrimmedUncorr             ;               
  Float_t CHSjetMassTrimmedUncorr            ;               
  Float_t PUPjetMassTrimmedUncorr            ;               
  Float_t GENjetMassTrimmedUncorr            ;               
  Float_t PFjetMassTrimmed                   ;               
  Float_t CHSjetMassTrimmed                  ;               
  Float_t PUPjetMassTrimmed                  ;               
  Float_t GENjetMassTrimmed                  ;               
  Float_t PFjetMassPrunedUncorr              ;              
  Float_t CHSjetMassPrunedUncorr             ;              
  Float_t PUPjetMassPrunedUncorr             ;              
  Float_t GENjetMassPrunedUncorr             ;              
  Float_t PFjetMassPruned                    ;              
  Float_t CHSjetMassPruned                   ;              
  Float_t PUPjetMassPruned                   ;              
  Float_t GENjetMassPruned                   ;              
  Float_t PFjetMassFilteredUncorr            ;                 
  Float_t CHSjetMassFilteredUncorr           ;                 
  Float_t PUPjetMassFilteredUncorr           ;                 
  Float_t GENjetMassFilteredUncorr           ;                 
  Float_t PFjetMassFiltered                  ;      
  Float_t CHSjetMassFiltered                 ;      
  Float_t PUPjetMassFiltered                 ;      
  Float_t GENjetMassFiltered                 ;      
  Float_t PFjetMassMMDTUncorr                ;        
  Float_t PFjetMassMMDT                      ;        
  Float_t PFjetMassMMDTFilteredUncorr        ;        
  Float_t PFjetMassMMDTFiltered              ;        
  Float_t CHSjetMassMMDTUncorr               ;         
  Float_t CHSjetMassMMDT                     ;         
  Float_t CHSjetMassMMDTFilteredUncorr       ;         
  Float_t CHSjetMassMMDTFiltered             ;         
  Float_t PUPjetMassMMDTUncorr               ;         
  Float_t PUPjetMassMMDT                     ;         
  Float_t PUPjetMassMMDTFilteredUncorr       ;         
  Float_t PUPjetMassMMDTFiltered             ;         
  Float_t GENjetMassMMDTUncorr               ;         
  Float_t GENjetMassMMDT                     ;         
  Float_t GENjetMassMMDTFilteredUncorr       ; 
  Float_t GENjetMassMMDTFiltered             ;  
  Float_t PFjetMassSoftDropUncorr            ;           
  Float_t PFjetMassSoftDrop                  ;           
  Float_t PFjetSoftDropSymmetry              ;           
  Float_t PFjetSoftDropDR                    ;           
  Float_t PFjetSoftDropMassDrop              ;           
  Float_t PFjetSoftDropEnergyLoss            ;           
  Float_t CHSjetMassSoftDropUncorr           ;           
  Float_t CHSjetMassSoftDrop                 ;           
  Float_t CHSjetSoftDropSymmetry             ;           
  Float_t CHSjetSoftDropDR                   ;           
  Float_t CHSjetSoftDropMassDrop             ;           
  Float_t CHSjetSoftDropEnergyLoss           ;           
  Float_t PUPjetMassSoftDropUncorr           ;            
  Float_t PUPjetMassSoftDrop                 ;            
  Float_t PUPjetSoftDropSymmetry             ;            
  Float_t PUPjetSoftDropDR                   ;            
  Float_t PUPjetSoftDropMassDrop             ;            
  Float_t PUPjetSoftDropEnergyLoss           ;            
  Float_t GENjetMassSoftDropUncorr           ;           
  Float_t GENjetMassSoftDrop                 ;           
  Float_t GENjetSoftDropSymmetry             ;           
  Float_t GENjetSoftDropDR                   ;           
  Float_t GENjetSoftDropMassDrop             ;           
  Float_t GENjetSoftDropEnergyLoss           ;   



  Float_t PFjetArea                          ;        
  Float_t CHSjetArea                         ;        
  Float_t PUPjetArea                         ;        
  Float_t GENjetArea                         ;                              
  Float_t PFjetAreaTrimmed                   ;               
  Float_t CHSjetAreaTrimmed                  ;               
  Float_t PUPjetAreaTrimmed                  ;               
  Float_t GENjetAreaTrimmed                  ;                           
  Float_t PFjetAreaPruned                    ;              
  Float_t CHSjetAreaPruned                   ;              
  Float_t PUPjetAreaPruned                   ;              
  Float_t GENjetAreaPruned                   ;                          
  Float_t PFjetAreaFiltered                  ;      
  Float_t CHSjetAreaFiltered                 ;      
  Float_t PUPjetAreaFiltered                 ;      
  Float_t GENjetAreaFiltered                 ;      
  Float_t PFjetAreaMMDT                      ;   
  Float_t CHSjetAreaMMDT                     ;         
  Float_t PUPjetAreaMMDT                     ;         
  Float_t GENjetAreaMMDT                     ;         
  Float_t PFjetAreaMMDTFiltered              ;        
  Float_t CHSjetAreaMMDTFiltered             ;         
  Float_t PUPjetAreaMMDTFiltered             ;         
  Float_t GENjetAreaMMDTFiltered             ;  
  Float_t PFjetAreaSoftDrop                  ;              
  Float_t CHSjetAreaSoftDrop                 ;                   
  Float_t PUPjetAreaSoftDrop                 ;                 
  Float_t GENjetAreaSoftDrop                 ;     

      
  Float_t PFjetNconst                          ;        
  Float_t CHSjetNconst                         ;        
  Float_t PUPjetNconst                         ;        
  Float_t GENjetNconst                         ;                              
  Float_t PFjetNconstTrimmed                   ;               
  Float_t CHSjetNconstTrimmed                  ;               
  Float_t PUPjetNconstTrimmed                  ;               
  Float_t GENjetNconstTrimmed                  ;                           
  Float_t PFjetNconstPruned                    ;              
  Float_t CHSjetNconstPruned                   ;              
  Float_t PUPjetNconstPruned                   ;              
  Float_t GENjetNconstPruned                   ;                          
  Float_t PFjetNconstFiltered                  ;      
  Float_t CHSjetNconstFiltered                 ;      
  Float_t PUPjetNconstFiltered                 ;      
  Float_t GENjetNconstFiltered                 ;      
  Float_t PFjetNconstMMDT                      ;   
  Float_t CHSjetNconstMMDT                     ;         
  Float_t PUPjetNconstMMDT                     ;         
  Float_t GENjetNconstMMDT                     ;         
  Float_t PFjetNconstMMDTFiltered              ;        
  Float_t CHSjetNconstMMDTFiltered             ;         
  Float_t PUPjetNconstMMDTFiltered             ;         
  Float_t GENjetNconstMMDTFiltered             ;  
  Float_t PFjetNconstSoftDrop                  ;              
  Float_t CHSjetNconstSoftDrop                 ;                   
  Float_t PUPjetNconstSoftDrop                 ;                 
  Float_t GENjetNconstSoftDrop                 ;           
        


  Float_t PFtau1                             ;   
  Float_t PFtau2                             ;   
  Float_t PFtau3                             ;   
  Float_t CHStau1                            ;    
  Float_t CHStau2                            ;    
  Float_t CHStau3                            ;    
  Float_t PUPtau1                            ;    
  Float_t PUPtau2                            ;    
  Float_t PUPtau3                            ;    
  Float_t GENtau1                            ;    
  Float_t GENtau2                            ;    
  Float_t GENtau3                            ;    
  Float_t PFcmsttJetMass                     ;            
  Float_t PFcmsttWMass                       ;            
  Float_t PFcmsttHelicity                    ;            
  Float_t PFcmsttNsubjets                    ;            
  Float_t PFcmsttMinMass                     ;            
  Float_t CHScmsttJetMass                    ;            
  Float_t CHScmsttWMass                      ;            
  Float_t CHScmsttHelicity                   ;            
  Float_t CHScmsttNsubjets                   ;            
  Float_t CHScmsttMinMass                    ;            
  Float_t PUPcmsttJetMass                    ;             
  Float_t PUPcmsttWMass                      ;             
  Float_t PUPcmsttHelicity                   ;             
  Float_t PUPcmsttNsubjets                   ;             
  Float_t PUPcmsttMinMass                    ;             
  Float_t GENcmsttJetMass                    ;            
  Float_t GENcmsttWMass                      ;            
  Float_t GENcmsttHelicity                   ;            
  Float_t GENcmsttNsubjets                   ;            
  Float_t GENcmsttMinMass                    ;            

  stringstream tempj;
  tempj << R*10;
  string Rstring = tempj.str();

  stringstream tempk;
  tempk << FirstEvent;
  string FirstEventString = tempk.str();

  stringstream templ;
  templ << FirstEvent+Nevents-1;
  string LastEventString = templ.str();

  string outname;
  if (Sample==2) outname = "out_0615_WW_"+algo+"_R"+Rstring+"_"+FirstEventString+"_"+LastEventString+".root";
  if (Sample==1) outname = "out_0615_RS3000T_"+algo+"_R"+Rstring+"_"+FirstEventString+"_"+LastEventString+".root";
  if (Sample==0) outname = "out_0615_QCD_"+algo+"_R"+Rstring+"_"+FirstEventString+"_"+LastEventString+".root";
  TFile *outfile = new TFile(outname.c_str(), "RECREATE");
  TTree * JetTree = new TTree("JetTree","Tree for saving jet info");
  JetTree->Branch("PFjetPt"                         ,&PFjetPt                      ,"PFjetPt"                      );
  JetTree->Branch("CHSjetPt"                        ,&CHSjetPt                     ,"CHSjetPt"                     );
  JetTree->Branch("PUPjetPt"                        ,&PUPjetPt                     ,"PUPjetPt"                     );
  JetTree->Branch("GENjetPt"                        ,&GENjetPt                     ,"GENjetPt"                     );
  JetTree->Branch("PFjetMassUncorr"                 ,&PFjetMassUncorr              ,"PFjetMassUncorr"              );
  JetTree->Branch("CHSjetMassUncorr"                ,&CHSjetMassUncorr             ,"CHSjetMassUncorr"             );
  JetTree->Branch("PUPjetMassUncorr"                ,&PUPjetMassUncorr             ,"PUPjetMassUncorr"             );
  JetTree->Branch("GENjetMassUncorr"                ,&GENjetMassUncorr             ,"GENjetMassUncorr"             );
  JetTree->Branch("PFjetPtUnc"                      ,&PFjetPtUnc                   ,"PFjetPtUnc"                   );
  JetTree->Branch("CHSjetPtUnc"                     ,&CHSjetPtUnc                  ,"CHSjetPtUnc"                  );
  JetTree->Branch("PUPjetPtUnc"                     ,&PUPjetPtUnc                  ,"PUPjetPtUnc"                  );
  JetTree->Branch("GENjetPtUnc"                     ,&GENjetPtUnc                  ,"GENjetPtUnc"                  );
  JetTree->Branch("PFjetMassTrimmedUncorr"          ,&PFjetMassTrimmedUncorr       ,"PFjetMassTrimmedUncorr"       );
  JetTree->Branch("CHSjetMassTrimmedUncorr"         ,&CHSjetMassTrimmedUncorr      ,"CHSjetMassTrimmedUncorr"      );
  JetTree->Branch("PUPjetMassTrimmedUncorr"         ,&PUPjetMassTrimmedUncorr      ,"PUPjetMassTrimmedUncorr"      );
  JetTree->Branch("GENjetMassTrimmedUncorr"         ,&GENjetMassTrimmedUncorr      ,"GENjetMassTrimmedUncorr"      );
  JetTree->Branch("PFjetMassPrunedUncorr"           ,&PFjetMassPrunedUncorr        ,"PFjetMassPrunedUncorr"        );
  JetTree->Branch("CHSjetMassPrunedUncorr"          ,&CHSjetMassPrunedUncorr       ,"CHSjetMassPrunedUncorr"       );
  JetTree->Branch("PUPjetMassPrunedUncorr"          ,&PUPjetMassPrunedUncorr       ,"PUPjetMassPrunedUncorr"       );
  JetTree->Branch("GENjetMassPrunedUncorr"          ,&GENjetMassPrunedUncorr       ,"GENjetMassPrunedUncorr"       );
  JetTree->Branch("PFjetMassFilteredUncorr"         ,&PFjetMassFilteredUncorr      ,"PFjetMassFilteredUncorr"      );
  JetTree->Branch("CHSjetMassFilteredUncorr"        ,&CHSjetMassFilteredUncorr     ,"CHSjetMassFilteredUncorr"     );
  JetTree->Branch("PUPjetMassFilteredUncorr"        ,&PUPjetMassFilteredUncorr     ,"PUPjetMassFilteredUncorr"     );
  JetTree->Branch("GENjetMassFilteredUncorr"        ,&GENjetMassFilteredUncorr     ,"GENjetMassFilteredUncorr"     );
  JetTree->Branch("PFjetMassMMDTUncorr"             ,&PFjetMassMMDTUncorr          ,"PFjetMassMMDTUncorr"          );
  JetTree->Branch("CHSjetMassMMDTUncorr"            ,&CHSjetMassMMDTUncorr         ,"CHSjetMassMMDTUncorr"         );
  JetTree->Branch("PUPjetMassMMDTUncorr"            ,&PUPjetMassMMDTUncorr         ,"PUPjetMassMMDTUncorr"         );
  JetTree->Branch("GENjetMassMMDTUncorr"            ,&GENjetMassMMDTUncorr         ,"GENjetMassMMDTUncorr"         );
  JetTree->Branch("PFjetMassMMDTFilteredUncorr"     ,&PFjetMassMMDTFilteredUncorr  ,"PFjetMassMMDTFilteredUncorr"  );
  JetTree->Branch("CHSjetMassMMDTFilteredUncorr"    ,&CHSjetMassMMDTFilteredUncorr ,"CHSjetMassMMDTFilteredUncorr" );
  JetTree->Branch("PUPjetMassMMDTFilteredUncorr"    ,&PUPjetMassMMDTFilteredUncorr ,"PUPjetMassMMDTFilteredUncorr" );
  JetTree->Branch("GENjetMassMMDTFilteredUncorr"    ,&GENjetMassMMDTFilteredUncorr ,"GENjetMassMMDTFilteredUncorr" );
  JetTree->Branch("PFjetMassSoftDropUncorr"         ,&PFjetMassSoftDropUncorr      ,"PFjetMassSoftDropUncorr"      );
  JetTree->Branch("CHSjetMassSoftDropUncorr"        ,&CHSjetMassSoftDropUncorr     ,"CHSjetMassSoftDropUncorr"     );
  JetTree->Branch("PUPjetMassSoftDropUncorr"        ,&PUPjetMassSoftDropUncorr     ,"PUPjetMassSoftDropUncorr"     );
  JetTree->Branch("GENjetMassSoftDropUncorr"        ,&GENjetMassSoftDropUncorr     ,"GENjetMassSoftDropUncorr"     );
  JetTree->Branch("PFjetMass"                       ,&PFjetMass                    ,"PFjetMass"                    );
  JetTree->Branch("CHSjetMass"                      ,&CHSjetMass                   ,"CHSjetMass"                   );
  JetTree->Branch("PUPjetMass"                      ,&PUPjetMass                   ,"PUPjetMass"                   );
  JetTree->Branch("GENjetMass"                      ,&GENjetMass                   ,"GENjetMass"                   );
  JetTree->Branch("PFjetMassTrimmed"                ,&PFjetMassTrimmed             ,"PFjetMassTrimmed"             );
  JetTree->Branch("CHSjetMassTrimmed"               ,&CHSjetMassTrimmed            ,"CHSjetMassTrimmed"            );
  JetTree->Branch("PUPjetMassTrimmed"               ,&PUPjetMassTrimmed            ,"PUPjetMassTrimmed"            );
  JetTree->Branch("GENjetMassTrimmed"               ,&GENjetMassTrimmed            ,"GENjetMassTrimmed"            );
  JetTree->Branch("PFjetMassPruned"                 ,&PFjetMassPruned              ,"PFjetMassPruned"              );
  JetTree->Branch("CHSjetMassPruned"                ,&CHSjetMassPruned             ,"CHSjetMassPruned"             );
  JetTree->Branch("PUPjetMassPruned"                ,&PUPjetMassPruned             ,"PUPjetMassPruned"             );
  JetTree->Branch("GENjetMassPruned"                ,&GENjetMassPruned             ,"GENjetMassPruned"             );
  JetTree->Branch("PFjetMassFiltered"               ,&PFjetMassFiltered            ,"PFjetMassFiltered"            );
  JetTree->Branch("CHSjetMassFiltered"              ,&CHSjetMassFiltered           ,"CHSjetMassFiltered"           );
  JetTree->Branch("PUPjetMassFiltered"              ,&PUPjetMassFiltered           ,"PUPjetMassFiltered"           );
  JetTree->Branch("GENjetMassFiltered"              ,&GENjetMassFiltered           ,"GENjetMassFiltered"           );
  JetTree->Branch("PFjetMassMMDT"                   ,&PFjetMassMMDT                ,"PFjetMassMMDT"                );
  JetTree->Branch("CHSjetMassMMDT"                  ,&CHSjetMassMMDT               ,"CHSjetMassMMDT"               );
  JetTree->Branch("PUPjetMassMMDT"                  ,&PUPjetMassMMDT               ,"PUPjetMassMMDT"               );
  JetTree->Branch("GENjetMassMMDT"                  ,&GENjetMassMMDT               ,"GENjetMassMMDT"               );
  JetTree->Branch("PFjetMassMMDTFiltered"           ,&PFjetMassMMDTFiltered        ,"PFjetMassMMDTFiltered"        );
  JetTree->Branch("CHSjetMassMMDTFiltered"          ,&CHSjetMassMMDTFiltered       ,"CHSjetMassMMDTFiltered"       );
  JetTree->Branch("PUPjetMassMMDTFiltered"          ,&PUPjetMassMMDTFiltered       ,"PUPjetMassMMDTFiltered"       );
  JetTree->Branch("GENjetMassMMDTFiltered"          ,&GENjetMassMMDTFiltered       ,"GENjetMassMMDTFiltered"       );
  JetTree->Branch("PFjetMassSoftDrop"               ,&PFjetMassSoftDrop            ,"PFjetMassSoftDrop"            );
  JetTree->Branch("CHSjetMassSoftDrop"              ,&CHSjetMassSoftDrop           ,"CHSjetMassSoftDrop"           );
  JetTree->Branch("PUPjetMassSoftDrop"              ,&PUPjetMassSoftDrop           ,"PUPjetMassSoftDrop"           );
  JetTree->Branch("GENjetMassSoftDrop"              ,&GENjetMassSoftDrop           ,"GENjetMassSoftDrop"           );
  JetTree->Branch("PFjetSoftDropSymmetry"           ,&PFjetSoftDropSymmetry        ,"PFjetSoftDropSymmetry"        );
  JetTree->Branch("PFjetSoftDropDR"                 ,&PFjetSoftDropDR              ,"PFjetSoftDropDR"              );
  JetTree->Branch("PFjetSoftDropMassDrop"           ,&PFjetSoftDropMassDrop        ,"PFjetSoftDropMassDrop"        );
  JetTree->Branch("PFjetSoftDropEnergyLoss"         ,&PFjetSoftDropEnergyLoss      ,"PFjetSoftDropEnergyLoss"      );
  JetTree->Branch("CHSjetSoftDropSymmetry"          ,&CHSjetSoftDropSymmetry       ,"CHSjetSoftDropSymmetry"       );
  JetTree->Branch("CHSjetSoftDropDR"                ,&CHSjetSoftDropDR             ,"CHSjetSoftDropDR"             );
  JetTree->Branch("CHSjetSoftDropMassDrop"          ,&CHSjetSoftDropMassDrop       ,"CHSjetSoftDropMassDrop"       );
  JetTree->Branch("CHSjetSoftDropEnergyLoss"        ,&CHSjetSoftDropEnergyLoss     ,"CHSjetSoftDropEnergyLoss"     );
  JetTree->Branch("PUPjetSoftDropSymmetry"          ,&PUPjetSoftDropSymmetry       ,"PUPjetSoftDropSymmetry"       );
  JetTree->Branch("PUPjetSoftDropDR"                ,&PUPjetSoftDropDR             ,"PUPjetSoftDropDR"             );
  JetTree->Branch("PUPjetSoftDropMassDrop"          ,&PUPjetSoftDropMassDrop       ,"PUPjetSoftDropMassDrop"       );
  JetTree->Branch("PUPjetSoftDropEnergyLoss"        ,&PUPjetSoftDropEnergyLoss     ,"PUPjetSoftDropEnergyLoss"     );
  JetTree->Branch("GENjetSoftDropSymmetry"          ,&GENjetSoftDropSymmetry       ,"GENjetSoftDropSymmetry"       );
  JetTree->Branch("GENjetSoftDropDR"                ,&GENjetSoftDropDR             ,"GENjetSoftDropDR"             );
  JetTree->Branch("GENjetSoftDropMassDrop"          ,&GENjetSoftDropMassDrop       ,"GENjetSoftDropMassDrop"       );
  JetTree->Branch("GENjetSoftDropEnergyLoss"        ,&GENjetSoftDropEnergyLoss     ,"GENjetSoftDropEnergyLoss"     );
  JetTree->Branch("PFjetArea"                       ,&PFjetArea                    ,"PFjetArea"                    );
  JetTree->Branch("CHSjetArea"                      ,&CHSjetArea                   ,"CHSjetArea"                   );
  JetTree->Branch("PUPjetArea"                      ,&PUPjetArea                   ,"PUPjetArea"                   );
  JetTree->Branch("GENjetArea"                      ,&GENjetArea                   ,"GENjetArea"                   );
  JetTree->Branch("PFjetAreaTrimmed"                ,&PFjetAreaTrimmed             ,"PFjetAreaTrimmed"             );
  JetTree->Branch("CHSjetAreaTrimmed"               ,&CHSjetAreaTrimmed            ,"CHSjetAreaTrimmed"            );
  JetTree->Branch("PUPjetAreaTrimmed"               ,&PUPjetAreaTrimmed            ,"PUPjetAreaTrimmed"            );
  JetTree->Branch("GENjetAreaTrimmed"               ,&GENjetAreaTrimmed            ,"GENjetAreaTrimmed"            );
  JetTree->Branch("PFjetAreaPruned"                 ,&PFjetAreaPruned              ,"PFjetAreaPruned"              );
  JetTree->Branch("CHSjetAreaPruned"                ,&CHSjetAreaPruned             ,"CHSjetAreaPruned"             );
  JetTree->Branch("PUPjetAreaPruned"                ,&PUPjetAreaPruned             ,"PUPjetAreaPruned"             );
  JetTree->Branch("GENjetAreaPruned"                ,&GENjetAreaPruned             ,"GENjetAreaPruned"             );
  JetTree->Branch("PFjetAreaFiltered"               ,&PFjetAreaFiltered            ,"PFjetAreaFiltered"            );
  JetTree->Branch("CHSjetAreaFiltered"              ,&CHSjetAreaFiltered           ,"CHSjetAreaFiltered"           );
  JetTree->Branch("PUPjetAreaFiltered"              ,&PUPjetAreaFiltered           ,"PUPjetAreaFiltered"           );
  JetTree->Branch("GENjetAreaFiltered"              ,&GENjetAreaFiltered           ,"GENjetAreaFiltered"           );
  JetTree->Branch("PFjetAreaMMDT"                   ,&PFjetAreaMMDT                ,"PFjetAreaMMDT"                );
  JetTree->Branch("CHSjetAreaMMDT"                  ,&CHSjetAreaMMDT               ,"CHSjetAreaMMDT"               );
  JetTree->Branch("PUPjetAreaMMDT"                  ,&PUPjetAreaMMDT               ,"PUPjetAreaMMDT"               );
  JetTree->Branch("GENjetAreaMMDT"                  ,&GENjetAreaMMDT               ,"GENjetAreaMMDT"               );
  JetTree->Branch("PFjetAreaMMDTFiltered"           ,&PFjetAreaMMDTFiltered        ,"PFjetAreaMMDTFiltered"        );
  JetTree->Branch("CHSjetAreaMMDTFiltered"          ,&CHSjetAreaMMDTFiltered       ,"CHSjetAreaMMDTFiltered"       );
  JetTree->Branch("PUPjetAreaMMDTFiltered"          ,&PUPjetAreaMMDTFiltered       ,"PUPjetAreaMMDTFiltered"       );
  JetTree->Branch("GENjetAreaMMDTFiltered"          ,&GENjetAreaMMDTFiltered       ,"GENjetAreaMMDTFiltered"       );
  JetTree->Branch("PFjetAreaSoftDrop"               ,&PFjetAreaSoftDrop            ,"PFjetAreaSoftDrop"            );
  JetTree->Branch("CHSjetAreaSoftDrop"              ,&CHSjetAreaSoftDrop           ,"CHSjetAreaSoftDrop"           );
  JetTree->Branch("PUPjetAreaSoftDrop"              ,&PUPjetAreaSoftDrop           ,"PUPjetAreaSoftDrop"           );
  JetTree->Branch("GENjetAreaSoftDrop"              ,&GENjetAreaSoftDrop           ,"GENjetAreaSoftDrop"           );
  JetTree->Branch("PFjetNconst"                     ,&PFjetNconst                  ,"PFjetNconst"                  );
  JetTree->Branch("CHSjetNconst"                    ,&CHSjetNconst                 ,"CHSjetNconst"                 );
  JetTree->Branch("PUPjetNconst"                    ,&PUPjetNconst                 ,"PUPjetNconst"                 );
  JetTree->Branch("GENjetNconst"                    ,&GENjetNconst                 ,"GENjetNconst"                 );
  JetTree->Branch("PFjetNconstTrimmed"              ,&PFjetNconstTrimmed           ,"PFjetNconstTrimmed"           );
  JetTree->Branch("CHSjetNconstTrimmed"             ,&CHSjetNconstTrimmed          ,"CHSjetNconstTrimmed"          );
  JetTree->Branch("PUPjetNconstTrimmed"             ,&PUPjetNconstTrimmed          ,"PUPjetNconstTrimmed"          );
  JetTree->Branch("GENjetNconstTrimmed"             ,&GENjetNconstTrimmed          ,"GENjetNconstTrimmed"          );
  JetTree->Branch("PFjetNconstPruned"               ,&PFjetNconstPruned            ,"PFjetNconstPruned"            );
  JetTree->Branch("CHSjetNconstPruned"              ,&CHSjetNconstPruned           ,"CHSjetNconstPruned"           );
  JetTree->Branch("PUPjetNconstPruned"              ,&PUPjetNconstPruned           ,"PUPjetNconstPruned"           );
  JetTree->Branch("GENjetNconstPruned"              ,&GENjetNconstPruned           ,"GENjetNconstPruned"           );
  JetTree->Branch("PFjetNconstFiltered"             ,&PFjetNconstFiltered          ,"PFjetNconstFiltered"          );
  JetTree->Branch("CHSjetNconstFiltered"            ,&CHSjetNconstFiltered         ,"CHSjetNconstFiltered"         );
  JetTree->Branch("PUPjetNconstFiltered"            ,&PUPjetNconstFiltered         ,"PUPjetNconstFiltered"         );
  JetTree->Branch("GENjetNconstFiltered"            ,&GENjetNconstFiltered         ,"GENjetNconstFiltered"         );
  JetTree->Branch("PFjetNconstMMDT"                 ,&PFjetNconstMMDT              ,"PFjetNconstMMDT"              );
  JetTree->Branch("CHSjetNconstMMDT"                ,&CHSjetNconstMMDT             ,"CHSjetNconstMMDT"             );
  JetTree->Branch("PUPjetNconstMMDT"                ,&PUPjetNconstMMDT             ,"PUPjetNconstMMDT"             );
  JetTree->Branch("GENjetNconstMMDT"                ,&GENjetNconstMMDT             ,"GENjetNconstMMDT"             );
  JetTree->Branch("PFjetNconstMMDTFiltered"         ,&PFjetNconstMMDTFiltered      ,"PFjetNconstMMDTFiltered"      );
  JetTree->Branch("CHSjetNconstMMDTFiltered"        ,&CHSjetNconstMMDTFiltered     ,"CHSjetNconstMMDTFiltered"     );
  JetTree->Branch("PUPjetNconstMMDTFiltered"        ,&PUPjetNconstMMDTFiltered     ,"PUPjetNconstMMDTFiltered"     );
  JetTree->Branch("GENjetNconstMMDTFiltered"        ,&GENjetNconstMMDTFiltered     ,"GENjetNconstMMDTFiltered"     );
  JetTree->Branch("PFjetNconstSoftDrop"             ,&PFjetNconstSoftDrop          ,"PFjetNconstSoftDrop"          );
  JetTree->Branch("CHSjetNconstSoftDrop"            ,&CHSjetNconstSoftDrop         ,"CHSjetNconstSoftDrop"         );
  JetTree->Branch("PUPjetNconstSoftDrop"            ,&PUPjetNconstSoftDrop         ,"PUPjetNconstSoftDrop"         );
  JetTree->Branch("GENjetNconstSoftDrop"            ,&GENjetNconstSoftDrop         ,"GENjetNconstSoftDrop"         );
  JetTree->Branch("PFtau1"                          ,&PFtau1                       ,"PFtau1"                       );
  JetTree->Branch("PFtau2"                          ,&PFtau2                       ,"PFtau2"                       );
  JetTree->Branch("PFtau3"                          ,&PFtau3                       ,"PFtau3"                       );
  JetTree->Branch("CHStau1"                         ,&CHStau1                      ,"CHStau1"                      );
  JetTree->Branch("CHStau2"                         ,&CHStau2                      ,"CHStau2"                      );
  JetTree->Branch("CHStau3"                         ,&CHStau3                      ,"CHStau3"                      );
  JetTree->Branch("PUPtau1"                         ,&PUPtau1                      ,"PUPtau1"                      );
  JetTree->Branch("PUPtau2"                         ,&PUPtau2                      ,"PUPtau2"                      );
  JetTree->Branch("PUPtau3"                         ,&PUPtau3                      ,"PUPtau3"                      );
  JetTree->Branch("GENtau1"                         ,&GENtau1                      ,"GENtau1"                      );
  JetTree->Branch("GENtau2"                         ,&GENtau2                      ,"GENtau2"                      );
  JetTree->Branch("GENtau3"                         ,&GENtau3                      ,"GENtau3"                      );
  JetTree->Branch("PFcmsttJetMass"                  ,&PFcmsttJetMass               ,"PFcmsttJetMass"               );
  JetTree->Branch("PFcmsttWMass"                    ,&PFcmsttWMass                 ,"PFcmsttWMass"                 );
  JetTree->Branch("PFcmsttHelicity"                 ,&PFcmsttHelicity              ,"PFcmsttHelicity"              );
  JetTree->Branch("PFcmsttNsubjets"                 ,&PFcmsttNsubjets              ,"PFcmsttNsubjets"              );
  JetTree->Branch("PFcmsttMinMass"                  ,&PFcmsttMinMass               ,"PFcmsttMinMass"               );
  JetTree->Branch("CHScmsttJetMass"                 ,&CHScmsttJetMass              ,"CHScmsttJetMass"              );
  JetTree->Branch("CHScmsttWMass"                   ,&CHScmsttWMass                ,"CHScmsttWMass"                );
  JetTree->Branch("CHScmsttHelicity"                ,&CHScmsttHelicity             ,"CHScmsttHelicity"             );
  JetTree->Branch("CHScmsttNsubjets"                ,&CHScmsttNsubjets             ,"CHScmsttNsubjets"             );
  JetTree->Branch("CHScmsttMinMass"                 ,&CHScmsttMinMass              ,"CHScmsttMinMass"              );
  JetTree->Branch("PUPcmsttJetMass"                 ,&PUPcmsttJetMass              ,"PUPcmsttJetMass"              );
  JetTree->Branch("PUPcmsttWMass"                   ,&PUPcmsttWMass                ,"PUPcmsttWMass"                );
  JetTree->Branch("PUPcmsttHelicity"                ,&PUPcmsttHelicity             ,"PUPcmsttHelicity"             );
  JetTree->Branch("PUPcmsttNsubjets"                ,&PUPcmsttNsubjets             ,"PUPcmsttNsubjets"             );
  JetTree->Branch("PUPcmsttMinMass"                 ,&PUPcmsttMinMass              ,"PUPcmsttMinMass"              );
  JetTree->Branch("GENcmsttJetMass"                 ,&GENcmsttJetMass              ,"GENcmsttJetMass"              );
  JetTree->Branch("GENcmsttWMass"                   ,&GENcmsttWMass                ,"GENcmsttWMass"                );
  JetTree->Branch("GENcmsttHelicity"                ,&GENcmsttHelicity             ,"GENcmsttHelicity"             );
  JetTree->Branch("GENcmsttNsubjets"                ,&GENcmsttNsubjets             ,"GENcmsttNsubjets"             );
  JetTree->Branch("GENcmsttMinMass"                 ,&GENcmsttMinMass              ,"GENcmsttMinMass"              );

  
  /////////////////////////////////////////////////////////////////////
  //  Get PF particles and Gen particles from the input files
  /////////////////////////////////////////////////////////////////////

  TChain fIn("Events");  
  std::ifstream input_file_list;
  if (Sample==2) input_file_list.open("fileListRSGWW1000.txt");
  if (Sample==1) input_file_list.open("fileListRS3000.txt");
  if (Sample==0) input_file_list.open("fileListQCD.txt");
  std::string tmp;

  if(!input_file_list.is_open())
  {
    cout << "Did not find input file list" << endl;
    return false;
  }
  while(1)
  {
    input_file_list >> tmp;
    if(!input_file_list.good()) break;
    fIn.Add(tmp.c_str());
  }

  TClonesArray *fPFPart  = new TClonesArray("baconhep::TPFPart");
  TClonesArray *fGenPart = new TClonesArray("baconhep::TGenParticle");
  fIn.SetBranchAddress("GenParticle",  &fGenPart);
  fIn.SetBranchAddress("PFPart",       &fPFPart);
  
  /////////////////////////////////////////////////////////////////////
  //  Setup Jet Energy Corrections
  /////////////////////////////////////////////////////////////////////

  std::string cmsenv = "/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_7_patch2/src/";
  std::vector<JetCorrectorParameters> corrParams;
  corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L1FastJet_AK5PF.txt"));
  corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L2Relative_AK5PF.txt"));
  corrParams.push_back(JetCorrectorParameters(cmsenv+"BaconProd/Utils/data/Summer13_V1_MC_L3Absolute_AK5PF.txt"));
  //corrParams.push_back(JetCorrectorParameter(cmsenv+'BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt'));
  JetCorrectorParameters     param(cmsenv+"BaconProd/Utils/data/Summer13_V1_DATA_Uncertainty_AK5PF.txt");
  FactorizedJetCorrector   *jetCorr = new FactorizedJetCorrector(corrParams);
  JetCorrectionUncertainty *jetUnc  = new JetCorrectionUncertainty(param);
  

  /////////////////////////////////////////////////////////////////////
  //  Setup N-subjettiness
  /////////////////////////////////////////////////////////////////////

  double beta = 1.0; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
  //double R0 = 0.8; // Characteristic jet radius for normalization
  double Rcut = 10000.0; // maximum R particles can be from axis to be included in jet (large value for no cutoff)   
  NsubParameters paraNsub(beta, R, Rcut);

  Njettiness nSubKT(Njettiness::kt_axes,paraNsub);
  Njettiness nSubMin(Njettiness::min_axes,paraNsub);
  Njettiness nSubOnePass(Njettiness::onepass_kt_axes,paraNsub);

  /////////////////////////////////////////////////////////////////////
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  //  START EVENT LOOP
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  /////////////////////////////////////////////////////////////////////

  std::vector<PseudoJet> particles;
  std::vector<PseudoJet> genparticles;  
  std::vector<PseudoJet> gentops;  
  std::vector<PseudoJet> genWs;  

  int nEventsInChain = fIn.GetEntries();
  cout<<" nEventsInChain "<<nEventsInChain<<" FirstEvent "<<FirstEvent<<" Nevents "<<Nevents<<endl;

  if (Nevents ==-1 ) Nevents = nEventsInChain;
  if (Nevents+FirstEvent > nEventsInChain ) {cout << "Not enough events in chain" << endl; exit(1);}

  for(int i0 = FirstEvent; i0 < Nevents+FirstEvent; i0++) { 
    fIn.GetEntry(i0);
    particles   .clear();
    genparticles.clear();
    gentops.clear();
    genWs.clear();

    if (i0%100==0) std::cout <<"Event: "<< i0 <<"  (Running from "<< FirstEvent << " -> " <<Nevents+FirstEvent-1 <<") - (total # of events in chain = "<<nEventsInChain<<")"<< std::endl;
    if (verbose) std::cout <<"Event: "<< i0 <<"  (Running from "<< FirstEvent << " -> " <<Nevents+FirstEvent-1 <<") - (total # of events in chain = "<<nEventsInChain<<") -  N PF candidates = " << fPFPart->GetEntriesFast() << std::endl;

    //Get the PF Candidates
    for( int i1 = 0; i1 < fPFPart->GetEntriesFast(); i1++){//9,entries loop,fill the vector particles with PF particles
      baconhep::TPFPart *pPartTmp = (baconhep::TPFPart*)((*fPFPart)[i1]);
      //Convert from Particle Flow to PseudoJet with id
      fastjet::PseudoJet pFastJet = convert(pPartTmp);
      //Build the collection 
      particles.push_back(pFastJet);
    }

    //Get the Gen Particles
    for( int i1 = 0; i1 < fGenPart->GetEntriesFast(); i1++){
      baconhep::TGenParticle *pPartTmp = (baconhep::TGenParticle*)((*fGenPart)[i1]);
      if ( fabs(pPartTmp->pdgId)  ==6 && pPartTmp->status >20 && pPartTmp->status<30 )
      {
        fastjet::PseudoJet topholder = convert(pPartTmp);
        gentops.push_back(topholder);
      }  
      if ( fabs(pPartTmp->pdgId)  ==24 && pPartTmp->status >20 && pPartTmp->status<30 )
      {
        fastjet::PseudoJet Wholder = convert(pPartTmp);
        genWs.push_back(Wholder);
      }  
      if(pPartTmp->status != 1) continue;
      //Convert gen particle to PseudoJet
      fastjet::PseudoJet pFastJet = convert(pPartTmp);
      //Build the collection 
      genparticles.push_back(pFastJet);
    }
    if (verbose) cout<<"found "<<genparticles.size()<<" final state gen particles and "<<gentops.size() <<" initial tops "<<endl;

    if (Sample==1 && gentops.size()<1 )
    {
      if (verbose) cout << "************************ Did not find tops - Event "<<i0<<" ************************" << endl;
      for( int i1 = 0; i1 < fGenPart->GetEntriesFast(); i1++){
        baconhep::TGenParticle *pPartTmp = (baconhep::TGenParticle*)((*fGenPart)[i1]);
        if (pPartTmp->pt>1000 && verbose) cout<<"     Particle "<<i1<<" pdgId "<<pPartTmp->pdgId<<" status "<<pPartTmp->status<<" parent "<<pPartTmp->parent<<" pt "<<pPartTmp->pt<<" eta "<<pPartTmp->eta<<" phi "<<pPartTmp->phi<<endl;
      }
      continue;
    } 
    if (Sample==2 && genWs.size()<1 )
    {
      if (verbose) cout << "************************ Did not find Ws - Event "<<i0<<" ************************" << endl;
      for( int i1 = 0; i1 < fGenPart->GetEntriesFast(); i1++){
        baconhep::TGenParticle *pPartTmp = (baconhep::TGenParticle*)((*fGenPart)[i1]);
        if (pPartTmp->pt>1000 && verbose) cout<<"     Particle "<<i1<<" pdgId "<<pPartTmp->pdgId<<" status "<<pPartTmp->status<<" parent "<<pPartTmp->parent<<" pt "<<pPartTmp->pt<<" eta "<<pPartTmp->eta<<" phi "<<pPartTmp->phi<<endl;
      }
      continue;
    } 

    /////////////////////////////////////////////////////////////////////
    //  Event PF particles
    //  - gen particles
    //  - pf particles
    //  - charge hadron subtracted particles
    //  - puppi weighted particles 
    /////////////////////////////////////////////////////////////////////

    //Lets set the particles collections in the Puppi container
    puppiContainer curEvent(particles); 
    curEvent.setGen(genparticles); 


    // Setup soft killer
    SoftKiller soft_killer   (0.4,0.4);
    //SoftKiller soft_killerCHS(4.0,0.5, !SelectorIsPupCharged());

  
    //Lets fetch the particle collections 
    vector<PseudoJet> gen_event       = curEvent.genFetch();
    vector<PseudoJet> puppi_event     = curEvent.puppiFetch(40.);
    vector<PseudoJet> pf_event        = curEvent.pfFetch();
    vector<PseudoJet> chs_event       = curEvent.pfchsFetch(-1);
    vector<PseudoJet> chs_event2GeV   = curEvent.pfchsFetch( 2.);
    vector<PseudoJet> soft_event      = soft_killer   (pf_event);
    //vector<PseudoJet> softCHS_event   = soft_killerCHS(chs_event);
    if (verbose ) cout<<"got particle collections"<<endl;
    // for(unsigned int i0 = 0; i0 < genJets.size(); i0++) {
    //   lIndex = i0;
    //   PseudoJet puppiJet   = match(genJets[i0],puppiJets);
    //   PseudoJet pfJet      = match(genJets[i0],pfJets   );
    //   PseudoJet chsJet     = match(genJets[i0],chsJets  );
    //   PseudoJet chs2GeVJet = match(genJets[i0],chs2GeVJets);
    //   PseudoJet softJet    = match(genJets[i0],softJets);
    //   // PseudoJet softCHSJet = match(genJets[i0],softCHSJets);
    //   setJet(genJets[i0],JGen    ,gen_event   ,false,jetCorr,jetUnc,gsn_cleanser);
    //   if(pfJet.pt()      != 0) setJet(pfJet ,     JPF     ,pf_event     ,false,jetCorr,jetUnc,gsn_cleanser);
    //   if(chsJet.pt()     != 0) setJet(chsJet,     JCHS    ,chs_event    ,true ,jetCorr,jetUnc,gsn_cleanser);
    //   if(chs2GeVJet.pt() != 0) setJet(chs2GeVJet, JCHS2GeV,chs_event2GeV,true ,jetCorr,jetUnc,gsn_cleanser);
    //   if(puppiJet.pt()   != 0) setJet(puppiJet  , JPup    ,puppi_event  ,true ,jetCorr,jetUnc,gsn_cleanser);
    //   if(softJet.pt()    != 0) setJet(softJet   , JSoft   ,soft_event   ,false,jetCorr,jetUnc,gsn_cleanser);
    //   // if(softCHSJet.pt() != 0) setJet(softCHSJet, JSoftCHS,softCHS_event,true ,jetCorr,jetUnc,gsn_cleanser);
    //   lOut->Fill();
    // }

    /////////////////////////////////////////////////////////////////////
    //  Cluster jets from different input particle collections
    /////////////////////////////////////////////////////////////////////

    AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));

    //Now lets define a jet collection we would like to cluster
    JetDefinition jet_def (jetalgo, R);

    // Do the actual clustering 
    //    - > get CHS collection first and apply pT cuts (just to speed things up)
    ClusterSequenceArea clus_seq_Gen     (gen_event    , jet_def  , area_def);

    Selector selectorquick = SelectorNHardest(1); 
    vector<PseudoJet> quickjets      = selectorquick(sorted_by_pt(clus_seq_Gen    .inclusive_jets()));
    
    if (quickjets[0].pt()<400){  
      if (verbose) cout<<"Event failed minimum Pt cut. genjet pt = "<< quickjets[0].pt() << endl;
      continue;
    }
    //    - > now cluster the rest
    ClusterSequenceArea clus_seq_PF      (pf_event     , jet_def  , area_def);
    ClusterSequenceArea clus_seq_Pup     (puppi_event  , jet_def  , area_def);
    //ClusterSequenceArea clus_seq_Gen     (gen_event    , jet_def  , area_def);
    ClusterSequenceArea clus_seq_CHS     (chs_event    , jet_def  , area_def);
    ClusterSequenceArea clus_seq_CHS2GeV (chs_event2GeV, jet_def  , area_def);
    ClusterSequenceArea clus_seq_Soft    (soft_event   , jet_def  , area_def);
    //lusterSequenceArea pSoftCHS (softCHS_event, jet_def  , area_def);
    

    //Now lets define a selector to get the leading jet
    Selector selector = SelectorNHardest(4); 

    //And lets select the leading jet from the above collections
    vector<PseudoJet> genJets       = selector(sorted_by_pt(clus_seq_Gen   .inclusive_jets()));
    vector<PseudoJet>  pfJets       = selector(sorted_by_pt(clus_seq_PF    .inclusive_jets()));
    vector<PseudoJet> chsJets       = selector(sorted_by_pt(clus_seq_CHS   .inclusive_jets()));
    vector<PseudoJet> pupJets       = selector(sorted_by_pt(clus_seq_Pup   .inclusive_jets()));
    vector<PseudoJet> sofJets       = selector(sorted_by_pt(clus_seq_Soft  .inclusive_jets()));
    // vector<PseudoJet> genJets       = sorted_by_pt(clus_seq_Gen   .inclusive_jets());
    // vector<PseudoJet>  pfJets       = sorted_by_pt(clus_seq_PF    .inclusive_jets());
    // vector<PseudoJet> chsJets       = sorted_by_pt(clus_seq_CHS   .inclusive_jets());
    // vector<PseudoJet> pupJets       = sorted_by_pt(clus_seq_Pup   .inclusive_jets());

    if (verbose ) cout<<"clustered jets"<<endl;

    // find gen jet matched top top quark
    PseudoJet genJet0 ;
    genWs       = sorted_by_pt(genWs);
    gentops       = sorted_by_pt(gentops);
    if (Sample==2) genJet0 = match(genWs[0],genJets);  
    if (Sample==1) genJet0 = match(gentops[0],genJets);  
    else genJet0 = genJets[0];



    if(genJet0.pt() < 1) 
    {
      cout<<"************************  did not find genjet matched to top ************************"<<endl;
      cout<<"   genWs.size() "<<genWs.size()<<" pt0 "<<genWs[0].perp()<<" pt1 "<<genWs[1].perp()<<endl;
      cout<<"   gentops.size() "<<gentops.size()<<" pt0 "<<gentops[0].perp()<<" pt1 "<<gentops[1].perp()<<endl;
      cout<<"   genJets.size() "<<genJets.size()<<" pt0 "<<genJets[0].perp()<<" pt1 "<<genJets[1].perp()<<" dr1 "<< deltaR(gentops[0],genJets[0])<<" dr2 "<<deltaR(gentops[0],genJets[1])<<endl;
      continue;
    }


    // find reco jets matched to gen jets
    PseudoJet chsJet0 = match(genJet0,chsJets);
    PseudoJet pfJet0  = match(genJet0,pfJets);
    PseudoJet pupJet0 = match(genJet0,pupJets);

    if( pfJet0.pt() < 1  || chsJet0.pt() < 1  ||  pupJet0.pt() < 1 ) 
    {
      if (verbose) cout<<"************************  did not find genjet matched to recojets ************************"<<endl;
      continue;
    }

    if(verbose ) {
      cout<<"genJet0 mass  "<<genJet0 .m()<<endl;
      cout<<"chsJet0  mass "<<chsJet0 .m()<<endl;
      cout<<"pfJet0   mass "<<pfJet0  .m()<<endl;
      cout<<"pupJet0  mass "<<pupJet0 .m()<<endl;
    }
   
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
    double pfPt    =  pfJet0.pt() * correction( pfJet0 ,jetCorr,rho.rho()); 
    double chsPt   = chsJet0.pt() * correction(chsJet0,jetCorr,rhoCHS.rho()); 
    double pupPt   = pupJet0.pt() * correction(pupJet0,jetCorr,0);         //Assume rho is zero for this
    double genPt   = genJet0.pt() * 1.;

    double pfPtUnc  = unc( pfJet0,jetUnc);
    double chsPtUnc = unc(chsJet0,jetUnc);
    double pupPtUnc = unc(chsJet0,jetUnc);
    double genPtUnc = unc(chsJet0,jetUnc);

    double pfMass    =  pfJet0.m() * correction( pfJet0, jetCorr,rho.rho()); 
    double chsMass   = chsJet0.m() * correction(chsJet0,jetCorr,rhoCHS.rho()); 
    double pupMass   = pupJet0.m() * correction(pupJet0,jetCorr,0);         //Assume rho is zero for this
    double genMass   = genJet0.m() * 1.;

    vector<PseudoJet> pfConstituents   = pfJet0.constituents();
    vector<PseudoJet> chsConstituents = chsJet0.constituents();
    vector<PseudoJet> pupConstituents = pupJet0.constituents();
    vector<PseudoJet> genConstituents = genJet0.constituents();

    // Fill tree and output some info
    PFjetMassUncorr  =  pfJet0.m() ;
    CHSjetMassUncorr = chsJet0.m();
    PUPjetMassUncorr = pupJet0.m();
    GENjetMassUncorr = genJet0.m();

    PFjetMass        = pfMass  ;
    CHSjetMass       = chsMass ;
    PUPjetMass       = pupMass ;
    GENjetMass       = genMass ;

    PFjetPt          = pfPt   ;
    CHSjetPt         = chsPt  ;
    PUPjetPt         = pupPt  ;
    GENjetPt         = genPt  ;

    PFjetPtUnc       = pfPtUnc  ;
    CHSjetPtUnc      = chsPtUnc ;
    PUPjetPtUnc      = pupPtUnc ;
    GENjetPtUnc      = genPtUnc ;

     PFjetArea       =  pfJet0.area() ;
    CHSjetArea       = chsJet0.area() ;
    PUPjetArea       = pupJet0.area() ;
    GENjetArea       = genJet0.area() ;

     PFjetNconst       = pfConstituents .size() ;
    CHSjetNconst       = chsConstituents.size() ;
    PUPjetNconst       = pupConstituents.size() ;
    GENjetNconst       = genConstituents.size() ;
    
    if (verbose ) cout<<" pfPt     "<<pfPt    <<endl;
    if (verbose ) cout<<" chsPt    "<<chsPt   <<endl;
    if (verbose ) cout<<" pupPt    "<<pupPt   <<endl;
    if (verbose ) cout<<" genPt    "<<genPt   <<endl;
    if (verbose ) cout<<" pfPtUnc  "<<pfPtUnc <<endl;
    if (verbose ) cout<<" chsPtUnc "<<chsPtUnc<<endl;
    if (verbose ) cout<<" pupPtUnc "<<pupPtUnc<<endl;
    if (verbose ) cout<<" genPtUnc "<<genPtUnc<<endl;
    if (verbose ) cout<<" pfMass   "<<pfMass  <<endl;
    if (verbose ) cout<<" chsMass  "<<chsMass <<endl;
    if (verbose ) cout<<" pupMass  "<<pupMass <<endl;
    if (verbose ) cout<<" genMass  "<<genMass <<endl;

    /////////////////////////////////////////////////////////////////////
    //  Jet Grooming 
    /////////////////////////////////////////////////////////////////////
  
    // HATS groups can experiment with adding new groomers with different grooming parameters
    // -  Compare grooming with the cambridge_algorithm, kt_algorithm, and antikt_algorithm 
    // -  Change the grooming parameters (example Rtrim, ptfrac etc.)
    // -  Compare jet mass, jet area, measure jet mass response for each
    // -  Apply jet corrections to groomed jets
    // -  Experiment with using all jet constituents, CHS constituents, PUPPI constituents

    // -- Trimming ------------------------------------------
    double Rtrim = 0.2;
    double ptfrac = 0.05;
    Filter trimmer1(JetDefinition(cambridge_algorithm, Rtrim), SelectorPtFractionMin(ptfrac) );
    PseudoJet trimmed_pfJet  = trimmer1( pfJet0 );
    PseudoJet trimmed_chsJet = trimmer1(chsJet0 );
    PseudoJet trimmed_pupJet = trimmer1(pupJet0 );
    PseudoJet trimmed_genJet = trimmer1(genJet0 );

    double mass_corr_trimmed_pfJet     = trimmed_pfJet .m() * correction(trimmed_pfJet ,jetCorr,rho.rho()); 
    double mass_corr_trimmed_chsJet    = trimmed_chsJet.m() * correction(trimmed_chsJet,jetCorr,rho.rho()); 
    double mass_corr_trimmed_pupJet    = trimmed_pupJet.m() * correction(trimmed_pupJet,jetCorr,rho.rho()); 
    double mass_corr_trimmed_genJet    = trimmed_genJet.m() * correction(trimmed_genJet,jetCorr,rho.rho()); 

    PFjetMassTrimmedUncorr  = trimmed_pfJet .m();
    CHSjetMassTrimmedUncorr = trimmed_chsJet.m();
    PUPjetMassTrimmedUncorr = trimmed_pupJet.m();
    GENjetMassTrimmedUncorr = trimmed_genJet.m();

    PFjetMassTrimmed        = mass_corr_trimmed_pfJet  ;
    CHSjetMassTrimmed       = mass_corr_trimmed_chsJet ;
    PUPjetMassTrimmed       = mass_corr_trimmed_pupJet ;
    GENjetMassTrimmed       = mass_corr_trimmed_genJet ;

     PFjetAreaTrimmed       = trimmed_pfJet .area() ;
    CHSjetAreaTrimmed       = trimmed_chsJet.area() ;
    PUPjetAreaTrimmed       = trimmed_pupJet.area() ;
    GENjetAreaTrimmed       = trimmed_genJet.area() ;

     PFjetNconstTrimmed       = trimmed_pfJet .constituents().size() ;
    CHSjetNconstTrimmed       = trimmed_chsJet.constituents().size() ;
    PUPjetNconstTrimmed       = trimmed_pupJet.constituents().size() ;
    GENjetNconstTrimmed       = trimmed_genJet.constituents().size() ;

    if (verbose) cout<<"done trimming"<<endl;
    // -- Pruning ------------------------------------------
    double  prune_zcut=0.1;
    double Dcut_pfJets  =  pfJet0 .m() /  pfJet0 .perp();
    double Dcut_chsJets = chsJet0 .m() / chsJet0.perp();
    double Dcut_pupJets = pupJet0 .m() / pupJet0.perp();
    double Dcut_genJets = genJet0 .m() / genJet0.perp();

    Pruner prune_Dcut_pfJets (cambridge_algorithm, prune_zcut, Dcut_pfJets);
    Pruner prune_Dcut_chsJets(cambridge_algorithm, prune_zcut, Dcut_chsJets);
    Pruner prune_Dcut_pupJets(cambridge_algorithm, prune_zcut, Dcut_pupJets);
    Pruner prune_Dcut_genJets(cambridge_algorithm, prune_zcut, Dcut_genJets);

    PseudoJet pruned_pfJet  = prune_Dcut_pfJets ( pfJet0);
    PseudoJet pruned_chsJet = prune_Dcut_chsJets(chsJet0);
    PseudoJet pruned_pupJet = prune_Dcut_pupJets(pupJet0);
    PseudoJet pruned_genJet = prune_Dcut_genJets(genJet0);

    double mass_corr_pruned_pfJet     = pruned_pfJet .m() * correction(pruned_pfJet ,jetCorr,rho.rho()); 
    double mass_corr_pruned_chsJet    = pruned_chsJet.m() * correction(pruned_chsJet,jetCorr,rho.rho()); 
    double mass_corr_pruned_pupJet    = pruned_pupJet.m() * correction(pruned_pupJet,jetCorr,rho.rho()); 
    double mass_corr_pruned_genJet    = pruned_genJet.m() * correction(pruned_genJet,jetCorr,rho.rho()); 

    PFjetMassPrunedUncorr  = pruned_pfJet .m();
    CHSjetMassPrunedUncorr = pruned_chsJet.m();
    PUPjetMassPrunedUncorr = pruned_pupJet.m();
    GENjetMassPrunedUncorr = pruned_genJet.m();

     PFjetMassPruned        = mass_corr_pruned_pfJet  ;
    CHSjetMassPruned       = mass_corr_pruned_chsJet ;
    PUPjetMassPruned       = mass_corr_pruned_pupJet ;
    GENjetMassPruned       = mass_corr_pruned_genJet ;

     PFjetAreaPruned       = pruned_pfJet .area() ;
    CHSjetAreaPruned       = pruned_chsJet.area() ;
    PUPjetAreaPruned       = pruned_pupJet.area() ;
    GENjetAreaPruned       = pruned_genJet.area() ;

     PFjetNconstPruned       = pruned_pfJet .constituents().size() ;
    CHSjetNconstPruned       = pruned_chsJet.constituents().size() ;
    PUPjetNconstPruned       = pruned_pupJet.constituents().size() ;
    GENjetNconstPruned       = pruned_genJet.constituents().size() ;


    // -- Filtering ------------------------------------------
    double Rfilt = 0.3;
    unsigned int nfilt = 3;
    Filter filter_CA3(JetDefinition(cambridge_algorithm, Rfilt), SelectorNHardest(nfilt));
    
    PseudoJet filtered_pfJet  = filter_CA3( pfJet0);
    PseudoJet filtered_chsJet = filter_CA3(chsJet0);
    PseudoJet filtered_pupJet = filter_CA3(pupJet0);
    PseudoJet filtered_genJet = filter_CA3(genJet0);

    double mass_corr_filtered_pfJet     = filtered_pfJet .m() * correction(filtered_pfJet ,jetCorr,rho.rho()); 
    double mass_corr_filtered_chsJet    = filtered_chsJet.m() * correction(filtered_chsJet,jetCorr,rho.rho()); 
    double mass_corr_filtered_pupJet    = filtered_pupJet.m() * correction(filtered_pupJet,jetCorr,rho.rho()); 
    double mass_corr_filtered_genJet    = filtered_genJet.m() * correction(filtered_genJet,jetCorr,rho.rho()); 

    PFjetMassFilteredUncorr  = filtered_pfJet .m();
    CHSjetMassFilteredUncorr = filtered_chsJet.m();
    PUPjetMassFilteredUncorr = filtered_pupJet.m();
    GENjetMassFilteredUncorr = filtered_genJet.m();

    PFjetMassFiltered        = mass_corr_filtered_pfJet  ;
    CHSjetMassFiltered       = mass_corr_filtered_chsJet ;
    PUPjetMassFiltered       = mass_corr_filtered_pupJet ;
    GENjetMassFiltered       = mass_corr_filtered_genJet ;

     PFjetAreaFiltered       = filtered_pfJet .area() ;
    CHSjetAreaFiltered       = filtered_chsJet.area() ;
    PUPjetAreaFiltered       = filtered_pupJet.area() ;
    GENjetAreaFiltered       = filtered_genJet.area() ;

     PFjetNconstFiltered       = filtered_pfJet .constituents().size() ;
    CHSjetNconstFiltered       = filtered_chsJet.constituents().size() ;
    PUPjetNconstFiltered       = filtered_pupJet.constituents().size() ;
    GENjetNconstFiltered       = filtered_genJet.constituents().size() ;



    // -- Modified mass drop ------------------------------------------
    //typedef contrib::ModifiedMassDropTagger MMDT;
    
    // use just a symmetry cut, with no mass-drop requirement
    double z_cut = 0.10;
    contrib::ModifiedMassDropTagger mmdt(z_cut);

    PseudoJet mmdt_pfJet  = mmdt( pfJet0 );
    PseudoJet mmdt_chsJet = mmdt(chsJet0 );
    PseudoJet mmdt_pupJet = mmdt(pupJet0 );
    PseudoJet mmdt_genJet = mmdt(genJet0 );
        if (verbose) cout<<"somemmdt"<<endl;

    //returns [0] uncorr mmdt mass [1] corr mmdt mass [2] uncorr mmdt filtered mass [3] corr mmdt filtered mass
    vector<double> masses_mmdt_pfJet  = mmdt_corr_and_filter(mmdt_pfJet, jetCorr,rho.rho());
    vector<double> masses_mmdt_chsJet = mmdt_corr_and_filter(mmdt_chsJet, jetCorr,rho.rho());
    vector<double> masses_mmdt_pupJet = mmdt_corr_and_filter(mmdt_pupJet, jetCorr,rho.rho());
    vector<double> masses_mmdt_genJet = mmdt_corr_and_filter(mmdt_genJet, jetCorr,rho.rho());
    
    if (verbose) cout<<"mmdt func"<<endl;
    if (verbose) cout<<" masses_mmdt_chsJet.size() "<< masses_mmdt_chsJet.size() <<endl;
    if (verbose) cout<<" masses_mmdt_pupJet.size() "<< masses_mmdt_pupJet.size() <<endl;
    if (verbose) cout<<" masses_mmdt_genJet.size() "<< masses_mmdt_genJet.size() <<endl;
    if (verbose) cout<<" masses_mmdt_pfJet .size() "<< masses_mmdt_pfJet .size() <<endl;

    if (masses_mmdt_pfJet.size() >0 && masses_mmdt_chsJet.size() >0  && masses_mmdt_pupJet.size() >0 && masses_mmdt_genJet.size() >0 )
   { // save to tree
    PFjetMassMMDTUncorr           = masses_mmdt_pfJet[0];
    PFjetMassMMDT                 = masses_mmdt_pfJet[1];
    PFjetMassMMDTFilteredUncorr   = masses_mmdt_pfJet[2];
    PFjetMassMMDTFiltered         = masses_mmdt_pfJet[3];

    CHSjetMassMMDTUncorr          = masses_mmdt_chsJet[0];
    CHSjetMassMMDT                = masses_mmdt_chsJet[1];
    CHSjetMassMMDTFilteredUncorr  = masses_mmdt_chsJet[2];
    CHSjetMassMMDTFiltered        = masses_mmdt_chsJet[3];

    PUPjetMassMMDTUncorr          = masses_mmdt_pupJet[0];
    PUPjetMassMMDT                = masses_mmdt_pupJet[1];
    PUPjetMassMMDTFilteredUncorr  = masses_mmdt_pupJet[2];
    PUPjetMassMMDTFiltered        = masses_mmdt_pupJet[3];

    GENjetMassMMDT                = masses_mmdt_genJet[0];
    GENjetMassMMDTFiltered        = masses_mmdt_genJet[2];

    if (verbose) cout<<"mmdt  "<<masses_mmdt_pfJet[0]<<" "<<masses_mmdt_pfJet[1]<<" "<<masses_mmdt_pfJet[2]<<" "<<masses_mmdt_pfJet[3]<<" "<<endl;
    }

    // -- Soft drop ------------------------------------------
    double beta = 1.0;
    double zcut = 0.1;
    double mu   = 1.0;
    SoftDrop soft_drop(beta, zcut, mu);

    PseudoJet sd_pfJet  = soft_drop(pfJets[0]);
    PseudoJet sd_chsJet = soft_drop(chsJets[0]);
    PseudoJet sd_pupJet = soft_drop(pupJets[0]);
    PseudoJet sd_genJet = soft_drop(genJets[0]);
    
   
    if (sd_pfJet!=0)
    {
      double mass_sd_corr     = sd_pfJet.m() * correction(sd_pfJet,jetCorr,rho.rho()); 
      PFjetMassSoftDropUncorr = sd_pfJet.m();
      PFjetMassSoftDrop       = mass_sd_corr;
      PFjetSoftDropSymmetry   = sd_pfJet.structure_of<contrib::SoftDrop>().symmetry();
      PFjetSoftDropDR         = sd_pfJet.structure_of<contrib::SoftDrop>().delta_R();
      PFjetSoftDropMassDrop   = sd_pfJet.structure_of<contrib::SoftDrop>().mu();
      PFjetSoftDropEnergyLoss = 1-sd_pfJet.pt()/pfJets[0].pt();
      PFjetAreaSoftDrop       = sd_pfJet .area() ;
      PFjetNconstSoftDrop     = sd_pfJet .constituents().size() ;

    }
    if (sd_chsJet!=0)
    {
      double mass_sd_corr      = sd_chsJet.m() * correction(sd_chsJet,jetCorr,rho.rho()); 
      CHSjetMassSoftDropUncorr = sd_chsJet.m();
      CHSjetMassSoftDrop       = mass_sd_corr;
      CHSjetSoftDropSymmetry   = sd_chsJet.structure_of<contrib::SoftDrop>().symmetry();
      CHSjetSoftDropDR         = sd_chsJet.structure_of<contrib::SoftDrop>().delta_R();
      CHSjetSoftDropMassDrop   = sd_chsJet.structure_of<contrib::SoftDrop>().mu();
      CHSjetSoftDropEnergyLoss = 1-sd_chsJet.pt()/chsJets[0].pt();
      CHSjetAreaSoftDrop       = sd_chsJet .area() ;
      CHSjetNconstSoftDrop     = sd_chsJet .constituents().size() ;
    } 
    if (sd_pupJet!=0)
    {
      double mass_sd_corr      = sd_pupJet.m() * correction(sd_pupJet,jetCorr,rho.rho()); 
      PUPjetMassSoftDropUncorr = sd_pupJet.m();
      PUPjetMassSoftDrop       = mass_sd_corr;
      PUPjetSoftDropSymmetry   = sd_pupJet.structure_of<contrib::SoftDrop>().symmetry();
      PUPjetSoftDropDR         = sd_pupJet.structure_of<contrib::SoftDrop>().delta_R();
      PUPjetSoftDropMassDrop   = sd_pupJet.structure_of<contrib::SoftDrop>().mu();
      PUPjetSoftDropEnergyLoss = 1-sd_pupJet.pt()/pupJets[0].pt();
      PUPjetAreaSoftDrop       = sd_pupJet .area() ;
      PUPjetNconstSoftDrop     = sd_pupJet .constituents().size() ;
    }
    if (sd_genJet!=0)
    {
      double mass_sd_corr      = sd_genJet.m() * correction(sd_genJet,jetCorr,rho.rho()); 
      GENjetMassSoftDropUncorr = sd_genJet.m();
      GENjetMassSoftDrop       = mass_sd_corr;
      GENjetSoftDropSymmetry   = sd_genJet.structure_of<contrib::SoftDrop>().symmetry();
      GENjetSoftDropDR         = sd_genJet.structure_of<contrib::SoftDrop>().delta_R();
      GENjetSoftDropMassDrop   = sd_genJet.structure_of<contrib::SoftDrop>().mu();
      GENjetSoftDropEnergyLoss = 1-sd_genJet.pt()/genJets[0].pt();
      GENjetAreaSoftDrop       = sd_genJet .area() ;
      GENjetNconstSoftDrop     = sd_genJet .constituents().size() ;
    }   
    if (verbose) cout<<" done soft drop"<<endl;

    /////////////////////////////////////////////////////////////////////
    //  Jet tagging 
    /////////////////////////////////////////////////////////////////////

    // -- N-subjettiness ------------------------------------------
    // ----- HATS groups can experiment with the following:
    // ----- a. different methods to find the axis (kt,onepass)
    // ----- b. modify beta (see NsubParameters above)
    // ----- Questions:
    // ----- 1. How are N-subjettiness and jet mass correlated? Plot tau3/tau2 before and after a top mass selection.



    vector<PseudoJet> jet_constituents_pf  = pfJets[0].constituents();      
    vector<PseudoJet> jet_constituents_chs = chsJets[0].constituents();      
    vector<PseudoJet> jet_constituents_pup = pupJets[0].constituents();      
    vector<PseudoJet> jet_constituents_gen = genJets[0].constituents(); 
    if (verbose) cout<<" N-subjett"<<endl;

    double tau1_pf  = nSubOnePass.getTau(1,jet_constituents_pf );
    double tau1_chs = nSubOnePass.getTau(1,jet_constituents_chs);
    double tau1_pup = nSubOnePass.getTau(1,jet_constituents_pup);
    double tau1_gen = nSubOnePass.getTau(1,jet_constituents_gen);
    double tau2_pf  = nSubOnePass.getTau(2,jet_constituents_pf );
    double tau2_chs = nSubOnePass.getTau(2,jet_constituents_chs);
    double tau2_pup = nSubOnePass.getTau(2,jet_constituents_pup);
    double tau2_gen = nSubOnePass.getTau(2,jet_constituents_gen);
    double tau3_pf  = nSubOnePass.getTau(3,jet_constituents_pf );
    double tau3_chs = nSubOnePass.getTau(3,jet_constituents_chs);
    double tau3_pup = nSubOnePass.getTau(3,jet_constituents_pup);
    double tau3_gen = nSubOnePass.getTau(3,jet_constituents_gen);

    PFtau1 = tau1_pf;
    PFtau2 = tau2_pf;
    PFtau3 = tau3_pf;
    CHStau1 = tau1_chs;
    CHStau2 = tau2_chs;
    CHStau3 = tau3_chs;
    PUPtau1 = tau1_pup;
    PUPtau2 = tau2_pup;
    PUPtau3 = tau3_pup;
    GENtau1 = tau1_gen;
    GENtau2 = tau2_gen;
    GENtau3 = tau3_gen;

    // Get N-subjettiness subjets
    // vector<fastjet::PseudoJet> onepass1jets = nSubOnePass.getJets(jet_constituents);
    // vector<fastjet::PseudoJet> onepass2jets = nSubOnePass.getJets(jet_constituents);
    // vector<fastjet::PseudoJet> onepass3jets = nSubOnePass.getJets(jet_constituents);

    // -- Pruned N-subjettiness ------------------------------------------
    // ------ HATS Groups: calculated N-subjettiness for pruned jets

    // -- Trimmed N-subjettiness ------------------------------------------
    // ------ HATS Groups: calculated N-subjettiness for trimmed jets

    // -- PUPPI N-subjettiness ------------------------------------------
    // ------ HATS Groups: calculated N-subjettiness for PUPPI jets


    // // -- JHU Top Tagger ------------------------------------------

    // double delta_p = 0.05;
    // double delta_r=0.19;
    // double cos_theta_W_max = 0.7; // the maximal allowed value of the W helicity angle

    // JHTopTagger jhu_top_tagger(delta_p, delta_r, cos_theta_W_max);
    // jhu_top_tagger.set_top_selector(SelectorMassRange(0,10000));  //140,250));
    // jhu_top_tagger.set_W_selector  (SelectorMassRange( 0,10000));  // 50,110));

    // PseudoJet jhu_top_candidate = jhu_top_tagger(chsJets_ca10[0]);

    // //JHU_Pt_Denom -> Fill( jets[0].perp() );

    // if (jhu_top_candidate != 0)
    // {
    //   //JHU_Pt_Numer ->Fill( tagged.perp() );
    //   double jhu_top_mass = jhu_top_candidate.m();
    //   double jhu_W_mass   = jhu_top_candidate.structure_of<JHTopTagger>().W().m();
    //   double jhu_cosTheta   = jhu_top_candidate.structure_of<JHTopTagger>().cos_theta_W();
    //   PseudoJet subjet1 = jhu_top_candidate.structure_of<JHTopTagger>().W1();
    //   PseudoJet subjet2 = jhu_top_candidate.structure_of<JHTopTagger>().W1();
    //   PseudoJet subjet3 = jhu_top_candidate.structure_of<JHTopTagger>().non_W();

    //   vector<PseudoJet> kept_subjets0 = jhu_top_candidate.structure_of<JHTopTagger>().W().pieces();
    //   vector<PseudoJet> kept_subjets1 = jhu_top_candidate.structure_of<JHTopTagger>().non_W().pieces();
    //   if (verbose) cout<<" JHU Nsubjets0 "<<kept_subjets0.size()<<endl;
    //   if (verbose) cout<<" JHU Nsubjets1 "<<kept_subjets1.size()<<endl;
    //   if (verbose) cout<<" JHU mass "<<jhu_top_mass <<" jhu_W_mass "<<jhu_W_mass<<" jhu_cosTheta "<< jhu_cosTheta<<endl;
    // } 

    if (verbose) cout<<" start CMS tt"<<endl;

    // // -- CMS Top Tagger ------------------------------------------
    double cms_delta_p = 0.05;
    double cms_delta_r=0.4;
    double A=0.0004;

    CMSTopTagger cms_top_tagger(cms_delta_p, cms_delta_r, A);
    // cms_top_tagger.set_top_selector(SelectorMassRange(140,250));
    // cms_top_tagger.set_W_selector  (SelectorMassRange( 50,1000));

    PseudoJet cms_top_candidate_pf  = cms_top_tagger(pfJets[0]);
    PseudoJet cms_top_candidate_chs = cms_top_tagger(chsJets[0]);
    PseudoJet cms_top_candidate_pup = cms_top_tagger(pupJets[0]);
    PseudoJet cms_top_candidate_gen = cms_top_tagger(genJets[0]);

    if (cms_top_candidate_pf != 0)
    {
      vector<PseudoJet> kept_subjets0 = cms_top_candidate_pf.structure_of<CMSTopTagger>().W().pieces();
      vector<PseudoJet> kept_subjets1 = cms_top_candidate_pf.structure_of<CMSTopTagger>().non_W().pieces();
      vector<PseudoJet> all_subjets = kept_subjets0;
      all_subjets.insert( all_subjets.end(), kept_subjets1.begin(), kept_subjets1.end() );

      PFcmsttJetMass  = cms_top_candidate_pf.m();
      PFcmsttWMass    = cms_top_candidate_pf.structure_of<CMSTopTagger>().W().m();
      PFcmsttHelicity = cms_top_candidate_pf.structure_of<CMSTopTagger>().cos_theta_W();
      PFcmsttNsubjets = all_subjets.size();
      PFcmsttMinMass  = calculate_minmass(all_subjets);

    } 

    if (cms_top_candidate_chs != 0)
    {

      vector<PseudoJet> kept_subjets0 = cms_top_candidate_chs.structure_of<CMSTopTagger>().W().pieces();
      vector<PseudoJet> kept_subjets1 = cms_top_candidate_chs.structure_of<CMSTopTagger>().non_W().pieces();
      vector<PseudoJet> all_subjets = kept_subjets0;
      all_subjets.insert( all_subjets.end(), kept_subjets1.begin(), kept_subjets1.end() );

      CHScmsttJetMass  = cms_top_candidate_chs.m();
      CHScmsttWMass    = cms_top_candidate_chs.structure_of<CMSTopTagger>().W().m();
      CHScmsttHelicity = cms_top_candidate_chs.structure_of<CMSTopTagger>().cos_theta_W();
      CHScmsttNsubjets = all_subjets.size();
      CHScmsttMinMass  = calculate_minmass(all_subjets);

    } 

    if (cms_top_candidate_pup != 0)
    {

      vector<PseudoJet> kept_subjets0 = cms_top_candidate_pup.structure_of<CMSTopTagger>().W().pieces();
      vector<PseudoJet> kept_subjets1 = cms_top_candidate_pup.structure_of<CMSTopTagger>().non_W().pieces();
      vector<PseudoJet> all_subjets = kept_subjets0;
      all_subjets.insert( all_subjets.end(), kept_subjets1.begin(), kept_subjets1.end() );

      PUPcmsttJetMass  = cms_top_candidate_pup.m();
      PUPcmsttWMass    = cms_top_candidate_pup.structure_of<CMSTopTagger>().W().m();
      PUPcmsttHelicity = cms_top_candidate_pup.structure_of<CMSTopTagger>().cos_theta_W();
      PUPcmsttNsubjets = all_subjets.size();
      PUPcmsttMinMass  = calculate_minmass(all_subjets);

    } 

    if (cms_top_candidate_gen != 0)
    {

      vector<PseudoJet> kept_subjets0 = cms_top_candidate_gen.structure_of<CMSTopTagger>().W().pieces();
      vector<PseudoJet> kept_subjets1 = cms_top_candidate_gen.structure_of<CMSTopTagger>().non_W().pieces();
      vector<PseudoJet> all_subjets = kept_subjets0;
      all_subjets.insert( all_subjets.end(), kept_subjets1.begin(), kept_subjets1.end() );

      GENcmsttJetMass  = cms_top_candidate_gen.m();
      GENcmsttWMass    = cms_top_candidate_gen.structure_of<CMSTopTagger>().W().m();
      GENcmsttHelicity = cms_top_candidate_gen.structure_of<CMSTopTagger>().cos_theta_W();
      GENcmsttNsubjets = all_subjets.size();
      GENcmsttMinMass  = calculate_minmass(all_subjets);

    } 


    // -- HEP Top Tagger ------------------------------------------

    //     double topmass=172.3;
    //     double wmass=80.4;
    //     HEPTopTagger::HEPTopTagger cm_toptag(clus_seq_PF,pfJets[0],topmass,wmass);
    //     cm_toptag.set_top_range(150.,200.);
    //     cout<< "========= Top Tagger ============" << endl;
    //     cm_toptag.run_tagger();
    //     cout<< "-------- setting  --------" << endl;
    //     cm_toptag.get_setting();
    //     cout<< "-------- resutls  --------" << endl;
    //     cm_toptag.get_info();

    //       if(cm_toptag.is_masscut_passed()){
    //   cout << "### masscut_passed ###" << endl;
    //   PseudoJet top=cm_toptag.top_candidate();
    //   PseudoJet b=cm_toptag.top_subjets().at(0);
    //   PseudoJet W1=cm_toptag.top_subjets().at(1);
    //   PseudoJet W2=cm_toptag.top_subjets().at(2);
    //   cout << "top mass: " << top.m() << endl;
    //   cout << "bottom mass: "<< b.m() << endl;
    //   cout << "W mass: "<< (W1+W2).m() << endl;
    // }
    // -- Qjets ------------------------------------------

    // int QJetsPreclustering = 999;
    // std::vector<fastjet::PseudoJet> constits;
    // unsigned int nqjetconstits = chsJets_ca10[0].constituents().size();
    // if (nqjetconstits < (unsigned int) QJetsPreclustering) constits = chsJets_ca10[0].constituents();
    // else constits = chsJets_ca10[0].associated_cluster_sequence()->exclusive_subjets_up_to(chsJets_ca10[0],QJetsPreclustering);
    // double qjet_vol = getQjetVolatility(constits, 25, i0*25) );
    // constits.clear();

     // double zcut =0.1; 
     // double dcut_fctr = 0.5;
     // double exp_min = 0.0;
     // double exp_max = 0.0;
     // double rigidity = 0.1;
     // double truncationFactor = 0.01;
    
    //  QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity, truncationFactor);

    // -- Shower deconstruction ------------------------------------------

    // -- Jet Charge ------------------------------------------

    // -- pruned mass drop ------------------------------------------

    // -- ECF ------------------------------------------

    //------------------------------------
    // EnergyCorrelatorRatio
    //------------------------------------

    // EnergyCorrelatorRatio C2beta0 (2,0. ,fastjet::EnergyCorrelator::pt_R);
    // EnergyCorrelatorRatio C2beta02(2,0.2,fastjet::EnergyCorrelator::pt_R);
    // EnergyCorrelatorRatio C2beta05(2,0.5,fastjet::EnergyCorrelator::pt_R);
    // EnergyCorrelatorRatio C2beta10(2,1.0,fastjet::EnergyCorrelator::pt_R);
    // EnergyCorrelatorRatio C2beta20(2,2.0,fastjet::EnergyCorrelator::pt_R);
    // double ca8_C2beta0  = C2beta0 ( chsJets_ca10[0] );
    // double ca8_C2beta02 = C2beta02( chsJets_ca10[0] );
    // double ca8_C2beta05 = C2beta05( chsJets_ca10[0] );
    // double ca8_C2beta10 = C2beta10( chsJets_ca10[0] );
    // double ca8_C2beta20 = C2beta20( chsJets_ca10[0] );

    // -- Planar Flow ------------------------------------------


    JetTree->Fill();


  }//end event loop

  outfile->cd();
  JetTree ->Write();
  outfile->Write();
  outfile->Close();
}