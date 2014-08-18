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
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include <cstdlib>

using namespace std;
using namespace fastjet;
using namespace contrib;


typedef std::map<double, std::map<double, std::vector<int> > > JetMAP;
typedef std::map<double, std::map<double, std::vector<int> > >::iterator JetMAPIt;

//Object Processors
GenLoader       *fGen      = 0; 
JetLoader       *fJet      = 0; 
MuonLoader      *fMuon     = 0; 
PFLoader        *fPFCand   = 0; 

TTree* load(std::string iName) { 
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  return lTree;
}

struct JetInfo {
    float pt;
    float ptcorr;
    float ptraw;
    float ptclean;
    float pttrim;
    float pttrimsafe;
    float ptconst;
    float ptunc;
    float eta;
    float phi;
    float m;
    float mraw;
    float mclean;
    float mtrim;
    float mtrimsafe;
    float mconst;
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

//// ===  FUNCTION  ============================================================
////         Name:  GetPUJetID
////  Description:  
//// ===========================================================================
//std::vector<double> GetPUJetID(PseudoJet &iJet, std::vector<PseudoJet> &neutrals, std::vector<PseudoJet> &chargedLV, std::vector<PseudoJet> &chargedPU)
//{
     //std::vector<double> output;
  
    //float sumpt = 0.;
    //float sumptchpv = 0.;
    //float sumptchpu = 0.;
    //float sumdrsqptsq = 0.;
    //float sumptsq = 0.;

    //// For Neutrals
    //for (int i = 0; i < neutrals.size(); ++i)
    //{
      //sumpt += neutrals.pt();
      //sumdrsqptsq += dr*dr*pt*pt;
      //sumptsq += pt*pt;
      
    //}

  //return true;
//}       // -----  end of function GetPUJetID  -----

//void setJet(PseudoJet &iJet,JetInfo &iJetI, ClusterSequenceArea &clust_seq_rho, bool iCHS,FactorizedJetCorrector *iJetCorr,JetCorrectionUncertainty *iJetUnc,JetCleanser &gsn_cleanser) {
void setJet(PseudoJet &iJet,JetInfo &iJetI, ClusterSequenceArea &clust_seq_rho, bool iCHS,FactorizedJetCorrector *iJetCorr,JetCorrectionUncertainty *iJetUnc) {
    vector<PseudoJet> neutrals,chargedLV,chargedPU;
    getConstitsForCleansing(iJet.constituents(),neutrals,chargedLV,chargedPU);

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
    double lUnc = unc       (iJet,iJetUnc);
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
    delete area_subtractor;
}

void setupTree(TTree *iTree,JetInfo &iJet,std::string iName) {
    iTree->Branch((iName+"pt"        ).c_str(),&iJet.pt        ,(iName+"pt/F"        ).c_str());
    iTree->Branch((iName+"ptcorr"    ).c_str(),&iJet.ptcorr    ,(iName+"ptcorr/F"    ).c_str());
    iTree->Branch((iName+"ptraw"     ).c_str(),&iJet.ptraw     ,(iName+"ptraw/F"     ).c_str());
    iTree->Branch((iName+"ptclean"   ).c_str(),&iJet.ptclean,   (iName+"ptclean/F"   ).c_str());
    iTree->Branch((iName+"pttrim"    ).c_str(),&iJet.pttrim    ,(iName+"pttrim/F"    ).c_str());
    iTree->Branch((iName+"pttrimsafe").c_str(),&iJet.pttrimsafe,(iName+"pttrimsafe/F").c_str());
    iTree->Branch((iName+"ptconst"   ).c_str(),&iJet.ptconst   ,(iName+"ptconst/F"   ).c_str());
    iTree->Branch((iName+"ptunc"     ).c_str(),&iJet.ptunc     ,(iName+"ptunc/F"     ).c_str());
    iTree->Branch((iName+"eta"       ).c_str(),&iJet.eta       ,(iName+"eta/F"       ).c_str());
    iTree->Branch((iName+"phi"       ).c_str(),&iJet.phi       ,(iName+"phi/F"       ).c_str());
    iTree->Branch((iName+"m"         ).c_str(),&iJet.m         ,(iName+"m/F"         ).c_str());
    iTree->Branch((iName+"mraw"      ).c_str(),&iJet.mraw      ,(iName+"mraw/F"      ).c_str());
    iTree->Branch((iName+"mtrim"     ).c_str(),&iJet.mtrim     ,(iName+"mtrim/F"     ).c_str());
    iTree->Branch((iName+"mtrimsafe" ).c_str(),&iJet.mtrimsafe ,(iName+"mtrimsafe/F" ).c_str());
    iTree->Branch((iName+"mclean"    ).c_str(),&iJet.mclean    ,(iName+"mclean/F"    ).c_str());
    iTree->Branch((iName+"mconst"    ).c_str(),&iJet.mconst    ,(iName+"mconst/F"    ).c_str());
}
//vector<PseudoJet> threeHardest(vector<PseudoJet> &iParts, JetDefinition &iJetDef, Selector &iSelector,std::vector<ClusterSequence> &iCSs) {
// cluster full event (hard + pileup)
//  vector<PseudoJet> threehardest = iSelector(sorted_by_pt(cs.inclusive_jets()));
//  iCSs.push_back(cs);
//  return threehardest;
//}
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
    iJet.pt         = -1;
    iJet.ptraw      = -1;
    iJet.ptclean    = -1;
    iJet.pttrim     = -1;
    iJet.pttrimsafe = -1;
    iJet.eta        = -1;
    iJet.phi        = -1;
    iJet.m          = -1;
    iJet.mraw       = -1;
    iJet.mtrim      = -1;
    iJet.mtrimsafe  = -1;
    iJet.mclean     = -1;
    iJet.mconst     = -1;
}

int main (int argc, char ** argv) {

  //gROOT->ProcessLine("#include <vector>");          
  int maxEvents     = atoi(argv[1]);
  std::string lName = argv[2];
  bool        lGen  = atoi(argv[3]);
  bool        lUseZ = atoi(argv[4]);
  bool        lMet  = atoi(argv[5]);

//----------------------------------------------------------------------------
//  Setup JEC on the fly
//----------------------------------------------------------------------------
  std::string cmsenv = getenv("CMSSW_BASE");
  const std::string globalTag = "DES19_62_V8";

  //AK5PF
  std::vector<JetCorrectorParameters> PFcorrParams;
  PFcorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L1FastJet_AK5PF.txt"));
  PFcorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L2Relative_AK5PF.txt"));
  PFcorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L3Absolute_AK5PF.txt"));
  //corrParams.push_back(JetCorrectorParameter(cmsenv+'BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt'));
  JetCorrectorParameters     PFparam(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_Uncertainty_AK5PF.txt");
  FactorizedJetCorrector   *PFjetCorr = new FactorizedJetCorrector(PFcorrParams);
  JetCorrectionUncertainty *PFjetUnc  = new JetCorrectionUncertainty(PFparam);

  //AK5PFCHS
  std::vector<JetCorrectorParameters> CHScorrParams;
  CHScorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L1FastJet_AK5PFchs.txt"));
  CHScorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L2Relative_AK5PFchs.txt"));
  CHScorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L3Absolute_AK5PFchs.txt"));
  //corrParams.push_back(JetCorrectorParameter(cmsenv+'BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt'));
  JetCorrectorParameters    CHSparam(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_Uncertainty_AK5PFchs.txt");
  FactorizedJetCorrector   *CHSjetCorr = new FactorizedJetCorrector(CHScorrParams);
  JetCorrectionUncertainty *CHSjetUnc  = new JetCorrectionUncertainty(CHSparam);

  //Puppi (L2L3)
  std::vector<JetCorrectorParameters> PuppicorrParams;
  PuppicorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L2Relative_AK5PFchs.txt"));
  PuppicorrParams.push_back(JetCorrectorParameters(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_L3Absolute_AK5PFchs.txt"));
  //corrParams.push_back(JetCorrectorParameter(cmsenv+'BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt'));
  JetCorrectorParameters    Puppiparam(cmsenv+"/src/SubstrAna/Summer14/data/"+globalTag+"_Uncertainty_AK5PFchs.txt");
  FactorizedJetCorrector   *PuppijetCorr = new FactorizedJetCorrector(PuppicorrParams);
  JetCorrectionUncertainty *PuppijetUnc  = new JetCorrectionUncertainty(Puppiparam);

  //Setup JetAlgos
  double R = 0.4;
  JetDefinition jet_def(antikt_algorithm,R);         // the jet definition....
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  Selector selector = SelectorNHardest(3);   // definition of a selector for the three hardest jets

  //Now setup cleansing
  JetDefinition subjet_def(kt_algorithm,0.2);
  //JetCleanser gsn_cleanser(subjet_def,JetCleanser::gaussian_cleansing,JetCleanser::input_nc_separate);
  //gsn_cleanser.SetGaussianParameters(0.617,0.62,0.15,0.22);

  //Now read a file
  //root://eoscms.cern.ch//store/group/phys_jetmet/ntran/PUPPI/miniSamples/62x/rsgww1000_62x_PU40BX50/ntuple_1_1_VQC.root";
  TTree *lTree = load(lName); 
  if(lTree->GetEntries() < maxEvents || maxEvents == -1) 
    maxEvents = lTree->GetEntries(); 

  //Declare Readers
  //fEvt      = new EvtLoader     (lTree);
  fMuon     = new MuonLoader    (lTree);
  fPFCand   = new PFLoader      (lTree,"Puppi_cff.py");
  fGen      = new GenLoader     (lTree);
  //fJet      = new JetLoader     (lTree);

//----------------------------------------------------------------------------
//  Setup Output
//----------------------------------------------------------------------------
  std::size_t pos = lName.find_last_of("/");
  std::string outfile = lName.substr(pos+1, lName.size() - pos);
  outfile = "Puppi_" +outfile;
  TFile *lFile = new TFile(outfile.c_str(),"RECREATE");
  TTree *lOut  = new TTree("Tree","Tree");

  //Setup Tree
  //fEvt    ->setupTree      (lOut);
  fMuon   ->setupTree      (lOut);
  fPFCand ->setupTree      (lOut);
  if(lGen) 
    fGen ->setupTree      (lOut);

  int lIndex = 0; 
  lOut->Branch("index",&lIndex,"lIndex/F");

  JetInfo JGen;     setupTree(lOut,JGen    ,"Gen"  );
  JetInfo JPF;      setupTree(lOut,JPF     ,"PF"   );
  JetInfo JPup;     setupTree(lOut,JPup    ,"Puppi");
  JetInfo JCHS;     setupTree(lOut,JCHS    ,"CHS"  );
  //JetInfo JCHS2GeV; setupTree(lOut,JCHS2GeV,"CHS2GeV");
  //JetInfo JSoft;    setupTree(lOut,JSoft   ,"SK"  );
  //JetInfo JSoftCHS; setupTree(lOut,JSoftCHS,"SKCHS");

  std::vector<TLorentzVector> lVetoes;
  for(int i0 = 0; i0 < maxEvents; i0++)
  { 
    lVetoes.clear();
    if(i0 % 500 == 0) 
      std::cout << "===> Processed " << i0 << " - Done : " << (float(i0)/float(maxEvents)) << std::endl;

    clear(JGen);
    clear(JPF);
    clear(JPup);
    clear(JCHS);
    //clear(JCHS2GeV);

    //fJet    ->load(i0); 
    fMuon   ->load(i0); 
    fGen    ->load(i0);
    bool lFindZ = fMuon->selectZ(lVetoes);
    if(!lFindZ && lUseZ) continue;
    if(lUseZ) continue;

    vector<PseudoJet> gen_event;
    gen_event  = fGen   ->genFetch();
    //bool lFind20Jet = (genJets[2].pt() > 20);
    //if(!lFind20Jet) continue;
    //////////////////////////////////////////////////////
    TLorentzVector lZ = fMuon->boson();
    //TLorentzVector lZ(0, 0, 0, 0);
    fPFCand->load(i0,lZ);
    vector<PseudoJet> puppi_event     = fPFCand->puppiFetch(lZ);
    if(lMet) lOut->Fill();
    if(lMet) continue;
    vector<PseudoJet> pf_event        = fPFCand->pfFetch();
    vector<PseudoJet> chs_event       = fPFCand->pfchsFetch(-1);
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
    ClusterSequenceArea pGen_Rho (gen_event    , jet_def_for_rho, area_def);
    ClusterSequenceArea pPup_Rho   (puppi_event, jet_def_for_rho, area_def);
    ClusterSequenceArea pPF_Rho     (pf_event   , jet_def_for_rho, area_def);
    ClusterSequenceArea pCHS_Rho    (chs_event  , jet_def_for_rho, area_def);

    vector<PseudoJet> loopJets = pfJets;
    if(lGen) loopJets = genJets;

    for(unsigned int i0 = 0; i0 < loopJets.size(); i0++) 
    {
      lIndex = i0;
      bool lMatch = false;
      for(unsigned int i1 = 0; i1 < lVetoes.size(); i1++) 
      { 
        double pDPhi = fabs(lVetoes[i1].Phi()-loopJets[i0].phi());
        if(pDPhi > TMath::Pi()*2.-pDPhi) pDPhi = TMath::Pi()*2.-pDPhi;
        double pDEta = fabs(lVetoes[i1].Eta()-loopJets[i0].eta());
        if(sqrt(pDPhi*pDPhi+pDEta*pDEta)  < 0.5) lMatch = true;
      }
      if(lMatch) continue;
      if(loopJets[i0].pt() < 1) continue;


      PseudoJet genJet     = PseudoJet();
      if (lGen) genJet = loopJets[i0];
      else genJet = match(loopJets[i0],genJets);
    
      PseudoJet pfJet      = match(loopJets[i0],pfJets   );
      PseudoJet chsJet     = match(loopJets[i0],chsJets  );
      PseudoJet puppiJet   = match(loopJets[i0],puppiJets);

      //PseudoJet chs2GeVJet = match(loopJets[i0],chs2GeVJets);
      //std::cout << "==> " << genJet.pt() << " -- " << genJets[i0].pt() << " -- " << loopJets[i0].pt() << std::endl;
      if(genJet.pt()     != 0) setJet(genJet  , JGen, pGen_Rho, false, PFjetCorr   , PFjetUnc   );
      if(pfJet.pt()      != 0) setJet(pfJet   , JPF , pPF_Rho , false, PFjetCorr   , PFjetUnc   );
      if(chsJet.pt()     != 0) setJet(chsJet  , JCHS, pCHS_Rho, true , CHSjetCorr  , CHSjetUnc  );
      if(puppiJet.pt()   != 0) setJet(puppiJet, JPup, pPup_Rho, true , PuppijetCorr, PuppijetUnc);
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
