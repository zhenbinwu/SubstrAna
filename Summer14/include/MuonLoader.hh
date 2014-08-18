#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "BaconAna/DataFormats/interface/TMuon.hh"

using namespace baconhep;

class MuonLoader { 
public:
  MuonLoader(TTree *iTree);
  ~MuonLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load (int iEvent);
  bool selectSingleMu();
  bool selectZ(std::vector<TLorentzVector> &iVetoes);
  bool vetoMu();
  bool passLoose(TMuon *iMuon);
  bool passTight(TMuon *iMuon);
  TLorentzVector muon();
  TLorentzVector boson();

protected: 
  TClonesArray *fMuons;
  TBranch      *fMuonBr;
  TTree        *fTree;
  float fZPt;
  float fZY;
  float fZEta;
  float fZM;
  float fZPhi;
  float fPt1;
  float fEta1;
  float fPhi1;
  float fPt2;
  float fEta2;
  float fPhi2;
};
