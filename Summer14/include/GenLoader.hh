#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "fastjet/PseudoJet.hh"
#include <vector>
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

using namespace baconhep;

class GenLoader { 
public:
  GenLoader(TTree *iTree);
  ~GenLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load (int iEvent);
  void fillGenEvent();
  //Fill specific Gen Info
  void selectBoson();
  //Helpters
  TGenParticle* getStatus1(int iId);
  bool isNeutrino(TGenParticle *iPart);
  //Get gen particles
  std::vector<fastjet::PseudoJet>  genFetch();
  fastjet::PseudoJet convert(TGenParticle *iPart);
protected: 
  TClonesArray  *fGens;
  TBranch       *fGenBr;
  TGenEventInfo *fGenInfo;
  TBranch       *fGenInfoBr;
  TTree         *fTree;
  float fQ;
  float fPId1;
  float fX1;
  float fPdf1;
  float fPId2;
  float fX2;
  float fPdf2;

  float fVPt;
  float fVEta;
  float fVPhi;
  float fVM;
  int   fVId;
  
  float fPt1;
  float fEta1;
  float fPhi1;
  float fM1;
  int   fId1;

  float fPt2;
  float fEta2;
  float fPhi2;
  float fM2;
  int   fId2;
};
