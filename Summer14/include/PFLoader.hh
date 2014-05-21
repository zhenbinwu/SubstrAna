#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "BaconAna/DataFormats/interface/TPFPart.hh"

using namespace baconhep;

class PFLoader { 
public:
  PFLoader(TTree *iTree);
  ~PFLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load (int iEvent);
  TLorentzVector met();
  TLorentzVector puppiFetch();
  std::vector<PseudoJet> puppiJets(std::vector<TLorentzVector> lVetoes);
  std::vector<PseudoJet> pfJets   (std::vector<TLorentzVector> lVetoes);
  std::vector<PseudoJet> chsJets  (std::vector<TLorentzVector> lVetoes);
  void getJets(std::vector < fastjet::PseudoJet > &constits,std::vector < fastjet::PseudoJet > &jets);
 
protected: 
  TClonesArray *fPFCands;
  TBranch      *fPFPartBr;
  TTree        *fTree;
  float fMet;
  float fMetPhi;
};
