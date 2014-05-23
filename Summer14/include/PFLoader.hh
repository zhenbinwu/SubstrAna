#include "fastjet/PseudoJet.hh"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "Dummy/Puppi/interface/RecoObj.hh"
#include "Dummy/Puppi/interface/PuppiContainer.h"
#include "BaconAna/DataFormats/interface/TPFPart.hh"

using namespace baconhep;

class PFLoader { 
public:
  PFLoader(TTree *iTree,std::string iName);
  ~PFLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load (int iEvent);
  TLorentzVector met();
  std::vector<fastjet::PseudoJet> puppiFetch();
  std::vector<fastjet::PseudoJet> pfFetch   (){ return fPFParticles; }
  std::vector<fastjet::PseudoJet> pfchsFetch(double iPt);
  void                            fetch();
  RecoObj            convert(TPFPart *iPart);
  fastjet::PseudoJet convert(RecoObj *iObj);  
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
  std::vector<RecoObj>            fAllParticles;
  std::vector<fastjet::PseudoJet> fPuppi;
  std::vector<fastjet::PseudoJet> fPFParticles;
  std::vector<fastjet::PseudoJet> fPFCHSParticles;
  PuppiContainer *fPuppiContainer;
};
