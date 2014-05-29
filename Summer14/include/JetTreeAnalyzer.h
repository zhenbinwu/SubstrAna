#ifndef JetTreeAnalyzer_h
#define JetTreeAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TStyle.h"

#include "TMath.h"
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

using namespace std;

class JetTreeAnalyzer{

 public :

  JetTreeAnalyzer(TTree *tree);
  virtual ~JetTreeAnalyzer();

  virtual void Init(TTree *tree);
  virtual int  GetEntry(Long64_t entry);

  virtual void bookHistograms(std::string suffix="", float maxpt=200.);
  virtual void fillHistograms(int maxEntries, float minPt);
  virtual void saveHistograms(TFile *file, std::string dir);

  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain                                                                                                                                          
  Int_t           fCurrent; //!current Tree number in a TChain                                                                                                                                                  

  // Declaration of leaf types                                                                                                                                                                                  
  int npu;
  std::vector<float>   *pt;
  std::vector<float>   *ptcorr;
  std::vector<float>   *ptraw;
  std::vector<float>   *ptclean;
  std::vector<float>   *pttrim;
  std::vector<float>   *pttrimsafe;
  std::vector<float>   *ptconst;
  std::vector<float>   *ptunc;
  std::vector<float>   *eta;
  std::vector<float>   *phi;
  std::vector<float>   *m;
  std::vector<float>   *mraw;
  std::vector<float>   *mtrim;
  std::vector<float>   *mtrimsafe;
  std::vector<float>   *mclean;
  std::vector<float>   *mconst;
  std::vector<int>     *nparticles;
  std::vector<int>     *nneutrals;
  std::vector<int>     *ncharged;
  std::vector<float>   *ptgen;
  std::vector<float>   *etagen;
  std::vector<float>   *phigen;
  std::vector<float>   *mgen;

  std::vector<float>   *mrawgen;
  std::vector<float>   *mtrimgen;
  std::vector<float>   *mtrimsafegen;
  std::vector<float>   *mcleangen;
  std::vector<float>   *mconstgen;

  std::vector<int>     *imatch;


  // List of branches
  TBranch        *b_npu;   //!                                                                                                                                                                                                     
  TBranch        *b_pt;   //!                                                                                                                                                                                                      
  TBranch        *b_ptcorr;   //!                                                                                                                                                                                                  
  TBranch        *b_ptraw;   //!                                                                                                                                                                                                   
  TBranch        *b_ptclean;   //!                                                                                                                                                                                                 
  TBranch        *b_pttrim;   //!                                                                                                                                                                                                  
  TBranch        *b_pttrimsafe;   //!                                                                                                                                                                                              
  TBranch        *b_ptconst;   //!                                                                                                                                                                                                 
  TBranch        *b_ptunc;   //!                                                                                                                                                                                                   
  TBranch        *b_eta;   //!                                                                                                                                                                                                     
  TBranch        *b_phi;   //!                                                                                                                                                                                                     
  TBranch        *b_m;   //!                                                                                                                                                                                                       
  TBranch        *b_mraw;   //!                                                                                                                                                                                                    
  TBranch        *b_mtrim;   //!                                                                                                                                                                                                   
  TBranch        *b_mtrimsafe;   //!                                                                                                                                                                                               
  TBranch        *b_mclean;   //!                                                                                                                                                                                                  
  TBranch        *b_mconst;   //!                                                                                                                                                                                                  
  TBranch        *b_nparticles;   //!                                                                                                                                                                                              
  TBranch        *b_nneutrals;   //!                                                                                                                                                                                               
  TBranch        *b_ncharged;   //!                                                                                                                                                                                                
  TBranch        *b_ptgen;   //!                                                                                                                                                                                                   
  TBranch        *b_etagen;   //!                                                                                                                                                                   
  TBranch        *b_phigen;   //!                                                                                                                                                                                                  
  TBranch        *b_mgen;   //!
  TBranch        *b_mrawgen;   //!
  TBranch        *b_mtrimgen;   //!
  TBranch        *b_mtrimsafegen;   //!
  TBranch        *b_mcleangen;   //!
  TBranch        *b_mconstgen;   //!
  TBranch        *b_imatch;   //!              

  // histograms declaration
  TH1F *hnjets;

  TH1F* hptgen;
  TH1F* hptgen_pu;
  TH1F* hptgen_good;

  TH1F* hptraw;
  TH1F* hptraw_pu;
  TH1F* hptraw_good;
  TH1F* hptraw_response;
  
  TH1F* hpt;
  TH1F* hpt_pu;
  TH1F* hpt_good;
  TH1F* hpt_response;

  TH1F* hptcorr;
  TH1F* hptcorr_pu;
  TH1F* hptcorr_good;
  TH1F* hptcorr_response;

  TH1F* heta;
  TH1F* heta_pu;
  TH1F* heta_good;

  TH1F* hnpu;
  TH1F* hnpu_pu;
  TH1F* hnpu_good;

  TH1F* hmraw;
  TH1F* hmraw_response;

  TH1F* hm;
  TH1F* hm_response;  

  TH1F* hmtrim;
  TH1F* hmtrim_response;

  TH1F* hmtrimsafe;
  TH1F* hmtrimsafe_response;

  TH1F* hmclean;
  TH1F* hmclean_response;

  TH1F* hmconst;
  TH1F* hmconst_response;

  TH1F* hnparticles;
  TH1F* hnneutrals;
  TH1F* hncharged;

  TH2F* hpt_response_vs_pt;
  TH2F* hptraw_response_vs_pt;
  TH2F* hmraw_response_vs_pt;
  TH2F* hm_response_vs_pt;
  TH2F* hmtrim_response_vs_pt;
  TH2F* hmtrimsafe_response_vs_pt;
  TH2F* hmclean_response_vs_pt;
  TH2F* hmconst_response_vs_pt;

  TH2F* hpt_response_vs_eta;
  TH2F* hptraw_response_vs_eta;
  TH2F* hmraw_response_vs_eta;
  TH2F* hm_response_vs_eta;
  TH2F* hmtrim_response_vs_eta;
  TH2F* hmtrimsafe_response_vs_eta;
  TH2F* hmclean_response_vs_eta;
  TH2F* hmconst_response_vs_eta;

  TH2F* hpt_response_vs_npu;
  TH2F* hptraw_response_vs_npu;
  TH2F* hmraw_response_vs_npu;
  TH2F* hm_response_vs_npu;
  TH2F* hmtrim_response_vs_npu;
  TH2F* hmtrimsafe_response_vs_npu;
  TH2F* hmclean_response_vs_npu;
  TH2F* hmconst_response_vs_npu;

  // leading jet
  TH1F* hptraw_leadjet;
  TH1F* hptraw_pu_leadjet;
  TH1F* hptraw_good_leadjet;
  TH1F* hptraw_response_leadjet;

  TH1F* hpt_leadjet;
  TH1F* hpt_pu_leadjet;
  TH1F* hpt_good_leadjet;
  TH1F* hpt_response_leadjet;
 
  TH1F* heta_leadjet;
  TH1F* heta_pu_leadjet;
  TH1F* heta_good_leadjet;

  TH1F* hmraw_leadjet;
  TH1F* hmraw_response_leadjet;

  TH1F* hm_leadjet;
  TH1F* hm_response_leadjet;

  TH1F* hmtrim_leadjet;
  TH1F* hmtrim_response_leadjet;

  TH1F* hmtrimsafe_leadjet;
  TH1F* hmtrimsafe_response_leadjet;

  TH1F* hmclean_leadjet;
  TH1F* hmclean_response_leadjet;

  TH1F* hmconst_leadjet;
  TH1F* hmconst_response_leadjet;

  TH1F* hnparticles_leadjet;
  TH1F* hnneutrals_leadjet;
  TH1F* hncharged_leadjet;



 private:

};
#endif
