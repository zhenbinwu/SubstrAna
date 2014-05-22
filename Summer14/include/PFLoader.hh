#include "fastjet/PseudoJet.hh"
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
  std::vector<fastjet::PseudoJet> puppiFetch();
  std::vector<fastjet::PseudoJet> puppiJets(std::vector<TLorentzVector> lVetoes);
  std::vector<fastjet::PseudoJet> pfJets   (std::vector<TLorentzVector> lVetoes);
  std::vector<fastjet::PseudoJet> chsJets  (std::vector<TLorentzVector> lVetoes);
  void getJets(std::vector < fastjet::PseudoJet > &constits,std::vector < fastjet::PseudoJet > &jets,std::vector<TLorentzVector> lVetoes);
 
protected: 
  TClonesArray *fPFCands;
  TBranch      *fPFCandBr;
  TTree        *fTree;
  float fMet;
  float fMetPhi;
  std::vector<fastjet::PseudoJet> fPuppi;
};
