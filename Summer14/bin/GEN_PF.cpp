#include <iostream>
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

#include "SubstrAna/Summer14/src/GroomingTreeProducer.cpp"


int main()
{
              
	std::vector <string> inputFileNames;
	
	inputFileNames.push_back("ntuple_1_1_VQC.root");
	inputFileNames.push_back("ntuple_2_1_WxV.root");
	inputFileNames.push_back("ntuple_3_1_12y.root");
	inputFileNames.push_back("ntuple_4_1_Dks.root");
	inputFileNames.push_back("ntuple_5_1_oHy.root");
	inputFileNames.push_back("ntuple_6_1_Moz.root");
	inputFileNames.push_back("ntuple_7_1_sFX.root");
	
	GroomingTreeProducer treeProducer;
	treeProducer.SetInputFiles(inputFileNames);
	treeProducer.prefixSamples = "/storage/a/ishvetso/JetSubstructureStudies/RSGraviton_JMET_samples/";
	JetDefinition jet_def(antikt_algorithm, 0.8);
	JetDefinition jet_def_filtering(antikt_algorithm, 0.3);
	treeProducer.jet_def = jet_def;
	treeProducer.jet_def_filtering = jet_def_filtering;
	treeProducer.outFileName = "AK8.root";
	treeProducer.ProduceTree();
	
	
}




