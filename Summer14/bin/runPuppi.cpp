//#include "../include/puppiContainer.hh"
#include "../include/GenLoader.hh"
#include "../include/MuonLoader.hh"
#include "../include/JetLoader.hh"
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
#include "boost/tokenizer.hpp"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TChain.h"
#include "TGraph.h"

#include <cstdlib>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace fastjet;
using namespace contrib;
const bool verbose = false;

//Object Processors
GenLoader       *fGen      = 0; 
JetLoader       *fJet      = 0; 
MuonLoader      *fMuon     = 0; 
PFLoader        *fPFCand   = 0; 
std::vector<TGraph*> iPuppiCorr;

TChain* load(std::string iName) { 
  TChain *chain = new TChain("Events");
  typedef boost::tokenizer<boost::char_separator<char> > passwdTokenizer;

  if (iName.find("list")!= std::string::npos)
  {
    std::cout <<  "input file " << iName << std::endl;
    std::fstream input(iName);
    for(std::string line; getline(input, line);)
    {
      if (line[0] == '#') continue;
      std::cout << "Add File: " << line << std::endl;
      chain->Add(line.c_str());
    }
  } else {
    boost::char_separator<char> tokenSep(",");
    passwdTokenizer tok(iName, tokenSep);
    for(passwdTokenizer::iterator curTok=tok.begin(); curTok!=tok.end(); ++curTok)
    {
      std::cout << "Add File: " << *curTok<< std::endl;
      chain->Add(curTok->c_str());
    }
  }

  return chain;
}

struct JetInfo {
    float pt;
    float ptcorr;
    float ptraw;
    //float ptclean;
    //float pttrim;
    //float pttrimsafe;
    //float ptconst;
    float ptunc;
    float eta;
    float phi;
    float m;
    float mraw;
    //float mclean;
    //float mtrim;
    //float mtrimsafe;
    //float mconst;

    float Beta;
    float BetaStar;     
    float MeanSqDeltaR; 
    float NCharged;     
    float NNeutrals;    

};

struct MetInfo {
  float fSumEt;
  float fMet;
  float fMetPhi;
  float fU1;
  float fU2;
};


double correctPhil(double iPt, double iEta);
void loadPhilJEC(const std::string &globalTag, std::vector<TGraph*> &iCorr);
bool DefaultJet(PseudoJet jet);
std::vector<TLorentzVector> GetCorJets(bool IsPuppi, std::vector<PseudoJet> &iJets, std::vector<TLorentzVector> &lVetoes, ClusterSequenceArea &clust_seq_rho, FactorizedJetCorrector *iJetCorr);
bool setMET(std::vector<PseudoJet> &iJets, MetInfo& iMET, std::vector<TLorentzVector> &lVetoes);
bool setMET(bool IsPuppi, std::vector<PseudoJet> &iJets, MetInfo& iMET, std::vector<TLorentzVector> &lVetoes, ClusterSequenceArea &clust_seq_rho, FactorizedJetCorrector *iJetCorr);
bool setupMETTree(TTree *iTree,MetInfo &iMet,std::string iName);
bool RemoveMuonPFCand(std::vector<PseudoJet> &iJets, std::vector<TLorentzVector> &lVetoes);

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
  if (iJetCorr == NULL) return 1;
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

// ===  FUNCTION  ============================================================
//         Name:  GetPUJetID
//  Description:  
// ===========================================================================
std::map<std::string, float> GetPUJetID(PseudoJet &iJet, std::vector<PseudoJet> &neutrals, std::vector<PseudoJet> &chargedLV, std::vector<PseudoJet> &chargedPU)
{
  std::map<std::string, float> outmap;

  float sumpt = 0.;
  float sumptch = 0.;
  float sumptchpv = 0.;
  float sumptchpu = 0.;
  float sumdrsqptsq = 0.;
  float sumptsq = 0.;
  float nc = 0;
  float nn = 0;

  // For Neutrals
  for (unsigned int i = 0; i < neutrals.size(); ++i)
  {
    double pEta = fabs(iJet.eta()-neutrals[i].eta());
    double pPhi = fabs(iJet.phi()-neutrals[i].phi());
    if(pPhi > 2.*TMath::Pi()-pPhi) pPhi =  2.*TMath::Pi()-pPhi;
    double dr = sqrt(pEta*pEta+pPhi*pPhi);
    double pt =  neutrals[i].pt();

    sumpt += pt;
    sumdrsqptsq += dr*dr*pt*pt;
    sumptsq += pt*pt;
    nn++;
  }

  // For charged hadron LV
  for (unsigned int i = 0; i < chargedLV.size(); ++i)
  {
    double pEta = fabs(iJet.eta()-chargedLV[i].eta());
    double pPhi = fabs(iJet.phi()-chargedLV[i].phi());
    if(pPhi > 2.*TMath::Pi()-pPhi) pPhi =  2.*TMath::Pi()-pPhi;
    double dr = sqrt(pEta*pEta+pPhi*pPhi);
    double pt =chargedLV[i].pt();

    sumpt += pt;
    sumptch += pt;
    sumptchpv += pt;
    sumdrsqptsq += dr*dr*pt*pt;
    sumptsq += pt*pt;
    nc++;
  }


  // For charged hadron PU
  for (unsigned int i = 0; i < chargedPU.size(); ++i)
  {
    double pEta = fabs(iJet.eta()-chargedPU[i].eta());
    double pPhi = fabs(iJet.phi()-chargedPU[i].phi());
    if(pPhi > 2.*TMath::Pi()-pPhi) pPhi =  2.*TMath::Pi()-pPhi;
    double dr = sqrt(pEta*pEta+pPhi*pPhi);
    double pt =chargedPU[i].pt();

    sumpt += pt;
    sumptch += pt;
    sumptchpu += pt;
    sumdrsqptsq += dr*dr*pt*pt;
    sumptsq += pt*pt;
    nc++;
  }



  if (sumptch > 0.) {
    outmap["Beta"] = sumptchpv/sumptch;
    outmap["BetaStar"] = sumptchpu/sumptch;
  } else {
    outmap["Beta"] = -999.;
    outmap["BetaStar"] = -999.;
  }

  if (sumptsq > 0.) {
    outmap["MeanSqDeltaR"] =  sumdrsqptsq/sumptsq;
  } else {
    outmap["MeanSqDeltaR"] = -999.;
  }

  outmap["NCharged"] = nc;
  outmap["NNeutrals"] = nn;

  return outmap;
}       // -----  end of function GetPUJetID  -----

//void setJet(PseudoJet &iJet,JetInfo &iJetI, ClusterSequenceArea &clust_seq_rho, bool iCHS,FactorizedJetCorrector *iJetCorr,JetCorrectionUncertainty *iJetUnc,JetCleanser &gsn_cleanser) {
void setJet(bool IsPuppi, PseudoJet &iJet,JetInfo &iJetI, ClusterSequenceArea &clust_seq_rho, bool iCHS,FactorizedJetCorrector *iJetCorr,JetCorrectionUncertainty *iJetUnc) {
    vector<PseudoJet> neutrals,chargedLV,chargedPU;
    getConstitsForCleansing(iJet.constituents(),neutrals,chargedLV,chargedPU);
    std::map<std::string, float> PUJetID = GetPUJetID(iJet, neutrals, chargedLV, chargedPU);

    //PseudoJet     lClean = gsn_cleanser(neutrals,chargedLV,chargedPU);
//----------------------------------------------------------------------------
//  PU Jet ID
//----------------------------------------------------------------------------
    
    // define safeAreaSub (PF)
    //AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
    //JetDefinition jet_def_for_rho(kt_algorithm, 0.4);
    Selector rho_range =  SelectorAbsRapMax(5.0);
    //ClusterSequenceArea clust_seq_rho(iParticles, jet_def_for_rho, area_def);
    // the two background estimators
    JetMedianBackgroundEstimator bge_rho (rho_range, clust_seq_rho);
    JetMedianBackgroundEstimator bge_rhom(rho_range, clust_seq_rho);
    BackgroundJetPtMDensity m_density;
    bge_rhom.set_jet_density_class(&m_density);
    
    // declare an area-median subtractor from this
    contrib::SafeAreaSubtractor *area_subtractor = 0;
    if(!iCHS) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom);
    if( iCHS) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom,SelectorIsPupCharged(),SelectorIsPupVertex());
    //iGMBE->set_particles(iParticles);
    PseudoJet lCorr =  (*area_subtractor)(iJet);

    //fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05)));
    //PseudoJet lTrim     = (trimmer)(iJet);
    //trimmer.set_subtractor(area_subtractor);
    //PseudoJet lTrimSafe = (trimmer)(iJet);
    //JetMedianBackgroundEstimator bge_rhoC(rho_range,jet_def_for_rho, area_def);
    //BackgroundJetScalarPtDensity *scalarPtDensity = new BackgroundJetScalarPtDensity();
    //bge_rhoC.set_jet_density_class(scalarPtDensity);
    //bge_rhoC.set_particles(iParticles);
    //contrib::ConstituentSubtractor subtractor(&bge_rhoC);
    //subtractor.use_common_bge_for_rho_and_rhom(true);
    //PseudoJet lConstit = subtractor(iJet);

    //Finally apply the JEC
    double lJEC = correction(iJet,iJetCorr,bge_rho.rho());  
    if (IsPuppi) 
    {
      lJEC = correction(iJet,iJetCorr,1.);  
      lJEC *= correctPhil(iJet.pt()*lJEC, iJet.eta()) ;
    }

    double lUnc = -999.;
    if (iJetUnc != NULL) lUnc = unc       (iJet,iJetUnc);
    iJetI.pt          = lCorr     .pt();
    iJetI.ptcorr      = iJet      .pt()*lJEC;
    iJetI.ptraw       = iJet      .pt();
    //iJetI.ptclean     = lClean    .pt();
    //iJetI.pttrim      = lTrim     .pt();
    //iJetI.pttrimsafe  = lTrimSafe .pt();
    //iJetI.ptconst     = lConstit  .pt();
    iJetI.ptunc       = lUnc;
    iJetI.eta         = iJet      .eta();
    iJetI.phi         = iJet      .phi();
    iJetI.mraw        = iJet      .m();
    iJetI.m           = lCorr     .m();
    //iJetI.mclean      = lClean    .m();
    //iJetI.mtrim       = lTrim     .m();
    //iJetI.mtrimsafe   = lTrimSafe .m();
    //iJetI.mconst      = lConstit  .m();

    iJetI.Beta         = PUJetID["Beta"];
    iJetI.BetaStar     = PUJetID["BetaStar"];
    iJetI.MeanSqDeltaR = PUJetID["MeanSqDeltaR"];
    iJetI.NCharged     = PUJetID["NCharged"];
    iJetI.NNeutrals    = PUJetID["NNeutrals"];

  
    delete area_subtractor;
}

void setupTree(TTree *iTree,JetInfo &iJet,std::string iName) {
    iTree->Branch((iName+"pt"           ) .c_str() ,&iJet.pt           ,(iName+"pt/F"           ) .c_str( ) ) ;
    iTree->Branch((iName+"ptcorr"       ) .c_str() ,&iJet.ptcorr       ,(iName+"ptcorr/F"       ) .c_str( ) ) ;
    iTree->Branch((iName+"ptraw"        ) .c_str() ,&iJet.ptraw        ,(iName+"ptraw/F"        ) .c_str( ) ) ;
    iTree->Branch((iName+"ptunc"        ) .c_str() ,&iJet.ptunc        ,(iName+"ptunc/F"        ) .c_str( ) ) ;
    iTree->Branch((iName+"eta"          ) .c_str() ,&iJet.eta          ,(iName+"eta/F"          ) .c_str( ) ) ;
    iTree->Branch((iName+"phi"          ) .c_str() ,&iJet.phi          ,(iName+"phi/F"          ) .c_str( ) ) ;
    iTree->Branch((iName+"m"            ) .c_str() ,&iJet.m            ,(iName+"m/F"            ) .c_str( ) ) ;
    iTree->Branch((iName+"mraw"         ) .c_str() ,&iJet.mraw         ,(iName+"mraw/F"         ) .c_str( ) ) ;
    //                                  PUJet ID       Variable
    iTree->Branch((iName+"Beta"         ) .c_str() ,&iJet.Beta         ,(iName+"Beta/F"         ) .c_str( ) ) ;
    iTree->Branch((iName+"BetaStar"     ) .c_str() ,&iJet.BetaStar     ,(iName+"BetaStar/F"     ) .c_str( ) ) ;
    iTree->Branch((iName+"MeanSqDeltaR" ) .c_str() ,&iJet.MeanSqDeltaR ,(iName+"MeanSqDeltaR/F" ) .c_str( ) ) ;
    iTree->Branch((iName+"NCharged"     ) .c_str() ,&iJet.NCharged     ,(iName+"NCharged/F"     ) .c_str( ) ) ;
    iTree->Branch((iName+"NNeutrals"    ) .c_str() ,&iJet.NNeutrals    ,(iName+"NNeutrals/F"    ) .c_str( ) ) ;


    //iTree->Branch((iName+"mtrim"     ).c_str(),&iJet.mtrim     ,(iName+"mtrim/F"     ).c_str());
    //iTree->Branch((iName+"mtrimsafe" ).c_str(),&iJet.mtrimsafe ,(iName+"mtrimsafe/F" ).c_str());
    //iTree->Branch((iName+"mclean"    ).c_str(),&iJet.mclean    ,(iName+"mclean/F"    ).c_str());
    //iTree->Branch((iName+"mconst"    ).c_str(),&iJet.mconst    ,(iName+"mconst/F"    ).c_str());
    //iTree->Branch((iName+"ptclean"   ).c_str(),&iJet.ptclean,   (iName+"ptclean/F"   ).c_str());
    //iTree->Branch((iName+"pttrim"    ).c_str(),&iJet.pttrim    ,(iName+"pttrim/F"    ).c_str());
    //iTree->Branch((iName+"pttrimsafe").c_str(),&iJet.pttrimsafe,(iName+"pttrimsafe/F").c_str());
    //iTree->Branch((iName+"ptconst"   ).c_str(),&iJet.ptconst   ,(iName+"ptconst/F"   ).c_str());
}

PseudoJet match(PseudoJet &iJet,vector<PseudoJet> &iJets) {

    std::list< std::pair<double, int> > lJetEta; //Sorting jet energy 
    for(unsigned int i0 = 0; i0 < iJets.size(); i0++)
    {
      // Ignore those soft pt jets
      if (iJets[i0].pt() < 0.5) continue; 

      double pEta = fabs(iJet.eta()-iJets[i0].eta());
      double pPhi = fabs(iJet.phi() - iJets[i0].phi());
      if(pPhi > 2.*TMath::Pi()-pPhi) pPhi =  2.*TMath::Pi()-pPhi;
      double deltaR = sqrt(pEta*pEta+pPhi*pPhi);
      if( deltaR > 0.3) continue;
      else
        lJetEta.push_back(std::make_pair(deltaR, i0));
    }

    lJetEta.sort();
    if (lJetEta.size() != 0)
    {
      int MatchedIdx = lJetEta.front().second;
      PseudoJet retJet = iJets[MatchedIdx];
      iJets.erase(iJets.begin()+MatchedIdx);
      return retJet;
    }
    else 
      return PseudoJet();
}

void clear(JetInfo &iJet) {
    iJet.pt         = -999;
    iJet.ptraw      = -999;
    //iJet.ptclean    = -999;
    //iJet.pttrim     = -999;
    //iJet.pttrimsafe = -999;
    iJet.eta        = -999;
    iJet.phi        = -999;
    iJet.m          = -999;
    iJet.mraw       = -999;
    //iJet.mtrim      = -999;
    //iJet.mtrimsafe  = -999;
    //iJet.mclean     = -999;
    //iJet.mconst     = -999;
    iJet.Beta       = -999;
    iJet.BetaStar   = -999;
    iJet.MeanSqDeltaR   = -999;
    iJet.NNeutrals = -999;
    iJet.NCharged = -999;
}

void clear(MetInfo &iMet) {
  iMet.fSumEt  = -999.;
  iMet.fMet    = -999.;
  iMet.fMetPhi = -999.;
  iMet.fU1     = -999.;
  iMet.fU2     = -999.;
}


int main (int argc, char ** argv) {

  //gROOT->ProcessLine("#include <vector>");          
  int maxEvents               = atoi(argv[1]);
  const std::string globalTag = argv[2];
  bool        lGen            = atoi(argv[3]);
  bool        lUseZ           = atoi(argv[4]);
  bool        lMet            = atoi(argv[5]);
  std::string lName           = argv[6];
//----------------------------------------------------------------------------
//  Setup JEC on the fly
//----------------------------------------------------------------------------
  std::string cmsenv = getenv("CMSSW_BASE");
  const std::string JetType = "AK4PF";

  //AK5PF
  std::vector<JetCorrectorParameters> PFcorrParams;
  PFcorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L1FastJet_"+JetType+".txt"));
  PFcorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L2Relative_"+JetType+".txt"));
  PFcorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L3Absolute_"+JetType+".txt"));
  //corrParams.push_back(JetCorrectorParameter(cmsenv+'BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_"+JetType+".txt'));
  FactorizedJetCorrector   *PFjetCorr = new FactorizedJetCorrector(PFcorrParams);
  JetCorrectionUncertainty *PFjetUnc  = NULL;
  //JetCorrectorParameters     PFparam(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_Uncertainty_"+JetType+".txt");
  //JetCorrectionUncertainty *PFjetUnc  = new JetCorrectionUncertainty(PFparam);

  //AK5PFCHS
  std::vector<JetCorrectorParameters> CHScorrParams;
  CHScorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L1FastJet_"+JetType+"chs.txt"));
  CHScorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L2Relative_"+JetType+"chs.txt"));
  CHScorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L3Absolute_"+JetType+"chs.txt"));
  //corrParams.push_back(JetCorrectorParameter(cmsenv+'BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_"+JetType+".txt'));
  FactorizedJetCorrector   *CHSjetCorr = new FactorizedJetCorrector(CHScorrParams);
  JetCorrectionUncertainty *CHSjetUnc  = NULL;
  //JetCorrectorParameters    CHSparam(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_Uncertainty_"+JetType+"chs.txt");
  //JetCorrectionUncertainty *CHSjetUnc  = new JetCorrectionUncertainty(CHSparam);

  //Puppi (L2L3)
  std::vector<JetCorrectorParameters> PuppicorrParams;

  PuppicorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/POSTLS170_V6_L1FastJet_AK7PF.txt"));
  PuppicorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/POSTLS170_V6_L2Relative_AK7PF.txt"));
  PuppicorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/POSTLS170_V6_L3Absolute_AK7PF.txt"));
  //corrParams.push_back(JetCorrectorParameter(cmsenv+'BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_"+JetType+".txt'));
  FactorizedJetCorrector   *PuppijetCorr = new FactorizedJetCorrector(PuppicorrParams);
  JetCorrectionUncertainty *PuppijetUnc  = NULL;
  loadPhilJEC(globalTag, iPuppiCorr);

  //JetCorrectorParameters    Puppiparam(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_Uncertainty_"+JetType+"chs.txt");
  //JetCorrectionUncertainty *PuppijetUnc  = new JetCorrectionUncertainty(Puppiparam);

  //Setup JetAlgos
  double R = 0.4;
  JetDefinition jet_def(antikt_algorithm,R);         // the jet definition....
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  Selector selector = SelectorNHardest(3);   // definition of a selector for the three hardest jets

  //Now setup cleansing
  //JetDefinition subjet_def(kt_algorithm,0.2);
  //JetCleanser gsn_cleanser(subjet_def,JetCleanser::gaussian_cleansing,JetCleanser::input_nc_separate);
  //gsn_cleanser.SetGaussianParameters(0.617,0.62,0.15,0.22);

  //Now read a file
  //root://eoscms.cern.ch//store/group/phys_jetmet/ntran/PUPPI/miniSamples/62x/rsgww1000_62x_PU40BX50/ntuple_1_1_VQC.root";
  TChain *lTree = load(lName); 
  if(lTree->GetEntries() < maxEvents || maxEvents == -1) 
    maxEvents = lTree->GetEntries(); 

  //Declare Readers
  //fEvt      = new EvtLoader     (lTree);
  fMuon     = new MuonLoader    (lTree);
  if (globalTag == "DES19_V1_MC" || globalTag == "AGE1K_V1_MC")
    fPFCand   = new PFLoader      (lTree, cmsenv+"/src/SubstrAna/Summer14/data/Puppi_PhaseI_cff.py");
  else if (globalTag == "PH2_1K_FB")
    fPFCand   = new PFLoader      (lTree, cmsenv+"/src/SubstrAna/Summer14/data/Puppi_PhaseII_cff.py");
  fGen      = new GenLoader     (lTree);

  //fJet      = new JetLoader     (lTree);

//----------------------------------------------------------------------------
//  Setup Output
//----------------------------------------------------------------------------
  std::size_t pos = lName.find_last_of("/");
  std::string outfile = lName.substr(pos+1, lName.size() - pos);
  TString outname(outfile);
  outname.ReplaceAll(".list", ".root");
  std::stringstream ss;
  ss <<"Puppi" << lGen << lUseZ << lMet << "_" << outname;
  std::cout << " Outfile name " << ss.str() << std::endl;
  TFile *lFile = new TFile(ss.str().c_str(),"RECREATE");
  TTree *lOut  = new TTree("Tree","Tree");

  //Setup Tree
  //fEvt    ->setupTree      (lOut);
  fMuon   ->setupTree      (lOut);
  fPFCand ->setupTree      (lOut);
  if(lGen) 
    fGen ->setupTree      (lOut);

  int lIndex = 0; 
  lOut->Branch("index",&lIndex,"lIndex/F");

  JetInfo JGen; setupTree(lOut,    JGen ,"Gen"     );
  JetInfo JPF;  setupTree(lOut,    JPF  ,"PF"      );
  JetInfo JPup; setupTree(lOut,    JPup ,"Puppi");
  JetInfo JCHS; setupTree(lOut,    JCHS ,"CHS"     );
  MetInfo MPF;  setupMETTree(lOut, MPF  ,"PF"      );
  MetInfo MPup; setupMETTree(lOut, MPup ,"Puppi");
  MetInfo MCHS; setupMETTree(lOut, MCHS ,"CHS"     );
  MetInfo MRawPF;  setupMETTree(lOut, MRawPF  ,"PFRaw"      );
  MetInfo MRawPup; setupMETTree(lOut, MRawPup ,"PuppiRaw");
  MetInfo MRawCHS; setupMETTree(lOut, MRawCHS ,"CHSRaw"     );
  //JetInfo JCHS2GeV; setupTree(lOut,JCHS2GeV,"CHS2GeV");
  //JetInfo JSoft;    setupTree(lOut,JSoft   ,"SK"  );
  //JetInfo JSoftCHS; setupTree(lOut,JSoftCHS,"SKCHS");
  std::vector<TLorentzVector> lVetoes;
  for(int i0 = 0; i0 < maxEvents; i0++)
  { 
    lVetoes.clear();
    if(i0 % 500 == 0) 
      std::cout << "===> Processed " << i0 << " - Done : " << (float(i0)/float(maxEvents)) << " with " << lTree->GetCurrentFile()->GetName() << std::endl;

    clear(JGen);
    clear(JPF);
    clear(JPup);
    clear(JCHS);
    //clear(JCHS2GeV);
    clear(MPF);
    clear(MPup);
    clear(MCHS);
    clear(MRawPF);
    clear(MRawPup);
    clear(MRawCHS);

    lTree->GetEntry(i0);
    //fJet    ->load(i0); 
    fMuon   ->load(i0); 
    fGen    ->load(i0);
    bool lFindZ = fMuon->selectZ(lVetoes);
    if(!lFindZ && lUseZ) continue;

    vector<PseudoJet> gen_event;
    gen_event  = fGen   ->genFetch();
    //bool lFind20Jet = (genJets[2].pt() > 20);
    //if(!lFind20Jet) continue;
    //////////////////////////////////////////////////////
    TLorentzVector lZ = fMuon->boson();
    fPFCand->load(i0,lVetoes);
    vector<PseudoJet> puppi_event     = fPFCand->puppiFetch(lVetoes);
    vector<PseudoJet> pf_event        = fPFCand->pfFetch();
    vector<PseudoJet> chs_event       = fPFCand->pfchsFetch(-1);

    if (verbose)
    {
      int EventCount[3][3] = {};
      for(unsigned int i=0; i < pf_event.size(); ++i)
      {
        PseudoJet jet = pf_event.at(i);
        if (jet.user_index() <= 1) EventCount[0][0]++;
        if (jet.user_index() == 2) EventCount[0][1]++;
        if (jet.user_index() == 3) EventCount[0][2]++;
        assert(jet.user_index() < 4);
      }
      for(unsigned int i=0; i < chs_event.size(); ++i)
      {
        PseudoJet jet = chs_event.at(i);
        if (jet.user_index() <= 1) EventCount[1][0]++;
        if (jet.user_index() == 2) EventCount[1][1]++;
        if (jet.user_index() == 3) EventCount[1][2]++;
        assert(jet.user_index() < 4);
      }
      for(unsigned int i=0; i < puppi_event.size(); ++i)
      {
        PseudoJet jet = puppi_event.at(i);
        if (jet.user_index() <= 1) EventCount[2][0]++;
        if (jet.user_index() == 2) EventCount[2][1]++;
        if (jet.user_index() == 3) EventCount[2][2]++;
        assert(jet.user_index() < 4);
      }
      std::cout << "PF size " << pf_event.size() <<" CHS size " << chs_event.size() << " Puppisize " << puppi_event.size() << std::endl;
      std::cout << "PF size : " << EventCount[0][0] <<" " << EventCount[0][1] <<" "<< EventCount[0][2] <<" "<< std::endl;
      std::cout << "CHS size : " << EventCount[1][0] <<" " << EventCount[1][1] <<" "<< EventCount[1][2] <<" "<< std::endl;
      std::cout << "Puppi size : " << EventCount[2][0] <<" " << EventCount[2][1] <<" "<< EventCount[2][2] <<" "<< std::endl;
    }



    //----------------------------------------------------------------------------
    //  For MET recoil
    //----------------------------------------------------------------------------
    if (lFindZ && lUseZ && lMet)
    {
      if(!RemoveMuonPFCand(pf_event, lVetoes)) continue;
      if(!RemoveMuonPFCand(chs_event, lVetoes)) continue;
      if(!RemoveMuonPFCand(puppi_event, lVetoes)) continue;
    }


    //vector<PseudoJet> chs_event2GeV   = fPFCand->pfchsFetch( 2.);
    ClusterSequenceArea pGen    (gen_event    ,jet_def,area_def);
    ClusterSequenceArea pPup    (puppi_event  ,jet_def,area_def);
    ClusterSequenceArea pPF     (pf_event     ,jet_def,area_def);
    ClusterSequenceArea pCHS    (chs_event    ,jet_def,area_def);
    //ClusterSequenceArea pCHS2GeV(chs_event2GeV,jet_def,area_def);
    //

    vector<PseudoJet> genJets     = sorted_by_pt(pGen    .inclusive_jets());
    vector<PseudoJet> puppiJets   = sorted_by_pt(pPup    .inclusive_jets());
    vector<PseudoJet> pfJets      = sorted_by_pt(pPF     .inclusive_jets());
    vector<PseudoJet> chsJets     = sorted_by_pt(pCHS    .inclusive_jets());
    //vector<PseudoJet> chs2GeVJets = sorted_by_pt(pCHS2GeV.inclusive_jets());

    if (verbose)
    {
      std::cout << "GenJet "<< count_if(genJets.begin(), genJets.end(), DefaultJet) 
        << " PFJet "<< count_if(pfJets.begin(), pfJets.end(), DefaultJet) 
        << " CHSJet "<< count_if(chsJets.begin(), chsJets.end(), DefaultJet) << std::endl; 
    }


    //////////////////////////////////////////////////////
    //SoftKiller soft_killer   (0.4,0.4);
    //SoftKiller soft_killerCHS(4.0,0.5, !SelectorIsPupCharged());
    //vector<PseudoJet> soft_event    = soft_killer   (pf_event);
    //vector<PseudoJet> softCHS_event = soft_killerCHS(chs_event);
    //ClusterSequenceArea pSoft   (soft_event   ,jet_def,area_def);
    //ClusterSequenceArea pSoftCHS(softCHS_event,jet_def,area_def);
    //vector<PseudoJet> softJets    = sorted_by_pt(pSoft   .inclusive_jets());
    //vector<PseudoJet> softCHSJets = sorted_by_pt(pSoftCHS.inclusive_jets());
    
    // Rho
    JetDefinition jet_def_for_rho(kt_algorithm, 0.4);
    ClusterSequenceArea pGen_Rho    (gen_event  , jet_def_for_rho, area_def);
    ClusterSequenceArea pPup_Rho    (puppi_event, jet_def_for_rho, area_def);
    ClusterSequenceArea pPF_Rho     (pf_event   , jet_def_for_rho, area_def);
    ClusterSequenceArea pCHS_Rho    (chs_event  , jet_def_for_rho, area_def);


    if(lMet)
    {
      setMET(false, pfJets,    MPF,  lVetoes, pPF_Rho,  PFjetCorr);
      setMET(false, chsJets,   MCHS, lVetoes, pCHS_Rho, CHSjetCorr);
      setMET(true, puppiJets, MPup, lVetoes, pPup_Rho, NULL);
      setMET(pfJets,    MRawPF,  lVetoes);
      setMET(chsJets,   MRawCHS, lVetoes);
      setMET(puppiJets, MRawPup, lVetoes);
      lOut->Fill();
      continue;
    }

    vector<PseudoJet> loopJets = pfJets;
    if(lGen) loopJets = genJets;

    for(unsigned int i0 = 0; i0 < loopJets.size(); i0++) 
    {
      lIndex = i0;
      if(loopJets[i0].pt() < 1) continue;

      PseudoJet genJet     = PseudoJet();
      if (lGen) genJet = loopJets[i0];
      else genJet = match(loopJets[i0],genJets);
    
      PseudoJet pfJet      = match(loopJets[i0],pfJets   );
      PseudoJet chsJet     = match(loopJets[i0],chsJets  );
      PseudoJet puppiJet   = match(loopJets[i0],puppiJets);

      //PseudoJet chs2GeVJet = match(loopJets[i0],chs2GeVJets);
      //std::cout << "==> " << genJet.pt() << " -- " << genJets[i0].pt() << " -- " << loopJets[i0].pt() << std::endl;
      if(genJet.pt()     != 0) setJet(false, genJet  , JGen, pGen_Rho, false, PFjetCorr   , PFjetUnc   );
      if(pfJet.pt()      != 0) setJet(false, pfJet   , JPF , pPF_Rho , false, PFjetCorr   , PFjetUnc   );
      if(chsJet.pt()     != 0) setJet(false, chsJet  , JCHS, pCHS_Rho, true , CHSjetCorr  , CHSjetUnc  );
      if(puppiJet.pt()   != 0) setJet(false, puppiJet, JPup, pPup_Rho, true , PuppijetCorr, PuppijetUnc);
      //if(chs2GeVJet.pt() != 0) setJet(chs2GeVJet, JCHS2GeV, chs_event2GeV, true , CHSjetCorr  , CHSjetUnc  , gsn_cleanser);

      //////////////////////////////////////////////////////
      //PseudoJet softJet    = match(loopJets[i0],softJets);
      //PseudoJet softCHSJet = match(loopJets[i0],softCHSJets);
      //if(softJet.pt()    != 0) setJet(softJet   , JSoft   , soft_event   , false, PFjetCorr   , PFjetUnc   , gsn_cleanser);
      //if(softCHSJet.pt() != 0) setJet(softCHSJet, JSoftCHS, softCHS_event, true , PFjetCorr   , PFjetUnc   , gsn_cleanser);

      lOut->Fill();
    }
  }
  lFile->cd();
  lOut->Write();
  lFile->Close();

  //delete PFjetCorr;
  //delete PFjetUnc;
  //delete CHSjetCorr;
  //delete CHSjetUnc;
  //delete PuppijetCorr;
  //delete PuppijetUnc;
  //delete fMuon;
  //delete fPFCand;
  //delete fGen;
  //delete lFile;
  //delete lOut;

  //delete FactorizedJetCorrector;
  //delete JetCorrectionUncertainty;
} 

// ===  FUNCTION  ============================================================
//         Name:  DefaultJet
//  Description:  
// ===========================================================================
bool DefaultJet(PseudoJet jet)
{
  return jet.pt() > 3;
}       // -----  end of function DefaultJet  -----

// ===  FUNCTION  ============================================================
//         Name:  GetCorJets
//  Description:  
// ===========================================================================
std::vector<TLorentzVector> GetCorJets(bool IsPuppi, std::vector<PseudoJet> &iJets, std::vector<TLorentzVector> &lVetoes, ClusterSequenceArea &clust_seq_rho, FactorizedJetCorrector *iJetCorr)
{
  std::vector<TLorentzVector> CorJets;
  Selector rho_range =  SelectorAbsRapMax(5.0);
  JetMedianBackgroundEstimator bge_rho (rho_range, clust_seq_rho);

  for(unsigned int i=0; i < iJets.size(); ++i)
  {
    PseudoJet ijet = iJets.at(i);
    double lJEC = correction(ijet,iJetCorr,bge_rho.rho());  
    if (IsPuppi) 
    {
      lJEC = correction(ijet,iJetCorr,1.);  
      lJEC *= correctPhil(ijet.pt()*lJEC, ijet.eta()) ;
    }
    TLorentzVector pVec(0,0,0,0);
    pVec.SetPtEtaPhiM(ijet.pt()*lJEC, ijet.eta(), ijet.phi(), ijet.m());
    CorJets.push_back(pVec);
  }
  return CorJets;
}       // -----  end of function GetCorJets  -----

// ===  FUNCTION  ============================================================
//         Name:  setupMETTree
//  Description:  
// ===========================================================================
bool setupMETTree(TTree *iTree,MetInfo &iMet,std::string iName)
{
  iTree->Branch((iName+"sumet"  ) .c_str() ,&iMet.fSumEt  ,(iName+"sumet/F"  ) .c_str( ) );
  iTree->Branch((iName+"met"    ) .c_str() ,&iMet.fMet    ,(iName+"met/F"    ) .c_str( ) );
  iTree->Branch((iName+"metphi" ) .c_str() ,&iMet.fMetPhi ,(iName+"metphi/F" ) .c_str( ) );
  iTree->Branch((iName+"u1"     ) .c_str() ,&iMet.fU1     ,(iName+"u1/F"     ) .c_str( ) );
  iTree->Branch((iName+"u2"     ) .c_str() ,&iMet.fU2     ,(iName+"u2/F"     ) .c_str( ) );
  return true;
}       // -----  end of function setupMETTree  -----

// ===  FUNCTION  ============================================================
//         Name:  setMET
//  Description:  
// ===========================================================================
bool setMET(bool IsPuppi, std::vector<PseudoJet> &iJets, MetInfo& iMET, std::vector<TLorentzVector> &lVetoes, ClusterSequenceArea &clust_seq_rho, FactorizedJetCorrector *iJetCorr)
{
  std::vector<TLorentzVector> iCorJets = GetCorJets(IsPuppi, iJets, lVetoes, clust_seq_rho, iJetCorr);

  double SumEt = 0.0;
  TLorentzVector lVec(0,0,0,0);
  TLorentzVector lUT(0,0,0,0);
  TLorentzVector lQT(0,0,0,0);

  for(std::vector<TLorentzVector>::const_iterator it=iCorJets.begin();
    it!=iCorJets.end(); ++it)
  {
    lVec -= *it;
    lUT += *it;
    SumEt += it->Pt();
  }

  for(std::vector<TLorentzVector>::const_iterator it=lVetoes.begin();
    it!=lVetoes.end(); ++it)
  {
    lVec -= *it;
    lQT += *it;
    SumEt += it->Pt();
  }

  iMET.fSumEt  = SumEt;
  iMET.fMet    = lVec.Pt();
  iMET.fMetPhi = lVec.Phi();
  double Dphi =  lUT.DeltaPhi(lQT);
  iMET.fU2 = lUT.Pt() * std::sin(Dphi);
  iMET.fU1 = lUT.Pt() * std::cos(Dphi);

  return true;
}       // -----  end of function setMET  -----

// ===  FUNCTION  ============================================================
//         Name:  setMET
//  Description:  
// ===========================================================================
bool setMET(std::vector<PseudoJet> &iJets, MetInfo& iMET, std::vector<TLorentzVector> &lVetoes)
{
  std::vector<TLorentzVector> iCorJets;
  for(unsigned int i=0; i < iJets.size(); ++i)
  {
    PseudoJet ijet = iJets.at(i);
    TLorentzVector pVec(0,0,0,0);
    pVec.SetPtEtaPhiM(ijet.pt(), ijet.eta(), ijet.phi(), ijet.m());
    iCorJets.push_back(pVec);
  }

  double SumEt = 0.0;
  TLorentzVector lVec(0,0,0,0);
  TLorentzVector lUT(0,0,0,0);
  TLorentzVector lQT(0,0,0,0);

  for(std::vector<TLorentzVector>::const_iterator it=iCorJets.begin();
    it!=iCorJets.end(); ++it)
  {
    lVec -= *it;
    lUT += *it;
    SumEt += it->Pt();
  }

  for(std::vector<TLorentzVector>::const_iterator it=lVetoes.begin();
    it!=lVetoes.end(); ++it)
  {
    lQT += *it;
    lVec -= *it;
    SumEt += it->Pt();
  }

  iMET.fSumEt  = SumEt;
  iMET.fMet    = lVec.Pt();
  iMET.fMetPhi = lVec.Phi();
  double Dphi =  lUT.DeltaPhi(lQT);
  iMET.fU2 = lUT.Pt() * std::sin(Dphi);
  iMET.fU1 = lUT.Pt() * std::cos(Dphi);

  return true;
}       // -----  end of function setMET  -----



// ===  FUNCTION  ============================================================
//         Name:  RemoveMuonPFCand
//  Description:  
// ===========================================================================
bool RemoveMuonPFCand(std::vector<PseudoJet> &iJets, std::vector<TLorentzVector> &lVetoes)
{
  int icoutn = 0;
  assert(iJets.size() != 0);
  for(std::vector<PseudoJet>::iterator it=iJets.begin(); it!=iJets.end(); )
  {
    PseudoJet pPar = *it;
    // recoil
    bool lMatch = false;
    for(unsigned int i1 = 0; i1 < lVetoes.size(); i1++) 
    { 
      double pDPhi = fabs(lVetoes[i1].Phi()-pPar.phi());
      if(pDPhi > TMath::Pi()*2.-pDPhi) pDPhi = TMath::Pi()*2.-pDPhi;
      if (pDPhi > 0.01) continue;
      double pDEta = fabs(lVetoes[i1].Eta()-pPar.eta());
      if (pDEta > 0.01 ) continue;
      if ( fabs(lVetoes[i1].Pt()-pPar.pt()) > 0.1) continue;
      if(sqrt(pDPhi*pDPhi+pDEta*pDEta)  < 0.05) lMatch = true;
    }

    if (lMatch)
    {
      it = iJets.erase(it);
      icoutn++;
    } else it++;
  }
  return icoutn == 2;
}       // -----  end of function RemoveMuonPFCand  -----

void loadPhilJEC(const std::string &globalTag, std::vector<TGraph*> &iCorr) {
  std::string cmsenv = getenv("CMSSW_BASE");
  std::string JECname = cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_Puppi.root";
  TFile *lFile = TFile::Open(JECname.c_str());

  for(int i0 = 0; i0 < 20; i0++) {
    std::stringstream pSS0;
    pSS0 << "puppi" << i0;
    TGraph* lF0 = (TGraph*) lFile->FindObjectAny(pSS0.str().c_str());
    iCorr.push_back(lF0);
  }
  return;
}

double correctPhil(const double iPt, const double iEta) {
  double lPt = iPt;
  if(lPt > 1000) return 1.;
  if(iPt < 30) lPt = 30;
  if (fabs(iEta) > 4.99) iEta = iEta/fabs(iEta) * 4.99;
  if(iPt < 10) lPt = 10;

  int iId  = 0.0;
  int lEta = int((iEta + 5)/0.5); 
  iId += lEta;
  double pCorr = iPuppiCorr[iId]->Eval(lPt);
  pCorr/=lPt;
  return pCorr;
}
