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
  bool vetoMu();
  bool passLoose(TMuon *iMuon);
  bool passTight(TMuon *iMuon);
  TLorentzVector muon();

protected: 
  TClonesArray *fMuons;
  TBranch      *fMuonBr;
  TTree        *fTree;
  float fPt;
  float fEta;
  float fPhi;
};
