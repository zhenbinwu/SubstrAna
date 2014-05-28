#ifndef  GROOMINGTREEPRODUCER
#define GROOMINGTREEPRODUCER
#include <iostream>
#include <dirent.h>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/GenericSubtractor.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/SafeSubtractor.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/SoftKiller.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "../include/puppiContainer.hh"

#include "BaconAna/DataFormats/interface/TPFPart.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"


using namespace std;
using namespace fastjet;
using namespace baconhep;


class GroomingTreeProducer
{       
    std::vector <string> inputFileNames;
    string inputFileName;
    
	  
  public:
  string prefixSamples, outFileName;
  JetDefinition jet_def, jet_def_filtering;
  GroomingTreeProducer();
  void SetInputFiles(std::vector <string> inputFileNames_);
  void SetInputFiles(string inputFileName_);
  void ProduceTree();
	
};

#endif



