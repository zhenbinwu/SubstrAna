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


int main(){
       
	
	
	TFile OutputFile ("output1.root", "RECREATE");
	
	cout << "this is running" << endl;
	
 	Float_t mass_PFjet, mass_Filtered_PFjet, mass_Pruned_PFjet, mass_Trimmed_PFjet, pt_PFjet, mass_GENjet, mass_Filtered_GENjet, mass_Pruned_GENjet, mass_Trimmed_GENjet, pt_GENjet, mass_subtracted;
	Int_t nPU;
	
	
	TTree *PF_tree = new TTree("PF", "PF");
	TTree *GEN_tree = new TTree("GEN", "GEN");
	
	PF_tree -> Branch("mass_PFjet",&mass_PFjet,"mass_PFjet/F");
	PF_tree -> Branch("pt_PFjet",&pt_PFjet,"pt_PFjet/F");
	PF_tree -> Branch("mass_Pruned_PFjet",&mass_Pruned_PFjet,"mass_Pruned_PFjet/F");
	PF_tree -> Branch("mass_Trimmed_PFjet",&mass_Trimmed_PFjet,"mass_Trimmed_PFjet/F");
	PF_tree -> Branch("mass_Filtered_PFjet",&mass_Filtered_PFjet,"mass_Filtered_PFjet/F");
	PF_tree -> Branch("mass_subtracted",&mass_subtracted,"mass_subtracted/F");
	PF_tree -> Branch("nPU",&nPU,"nPU/I");
	
	GEN_tree -> Branch("mass_GENjet",&mass_GENjet,"mass_GENjet/F");
	GEN_tree -> Branch("pt_GENjet",&pt_GENjet,"pt_GENjet/F");
	GEN_tree -> Branch("mass_Pruned_GENjet",&mass_Pruned_GENjet,"mass_Pruned_GENjet/F");
	GEN_tree -> Branch("mass_Trimmed_GENjet",&mass_Trimmed_GENjet,"mass_Trimmed_GENjet/F");
	GEN_tree -> Branch("mass_Filtered_GENjet",&mass_Filtered_GENjet,"mass_Filtered_GENjet/F");
	GEN_tree -> Branch("nPU",&nPU,"nPU/I");

        
	std::vector <string> inputFileNames;
	
	inputFileNames.push_back("ntuple_1_1_VQC.root");
	inputFileNames.push_back("ntuple_2_1_WxV.root");
	inputFileNames.push_back("ntuple_3_1_12y.root");
	inputFileNames.push_back("ntuple_4_1_Dks.root");
	inputFileNames.push_back("ntuple_5_1_oHy.root");
	inputFileNames.push_back("ntuple_6_1_Moz.root");
	inputFileNames.push_back("ntuple_7_1_sFX.root");
	
	
	for (unsigned int jj = 0; jj < inputFileNames.size(); jj++)
	{
	  
	  TFile* fIn = TFile::Open(("/storage/a/ishvetso/JetSubstructureStudies/RSGraviton_JMET_samples/" + inputFileNames.at(jj)).c_str());
	  TTree* tree = (TTree*) fIn->Get("Events");
	  TClonesArray *PF = new TClonesArray("baconhep::TPFPart");
	  TClonesArray *Gen = new TClonesArray("baconhep::TGenParticle");
	  TEventInfo *ev_Info = new TEventInfo();
	  	  
	  tree -> SetBranchAddress("PFPart", &PF);
	  tree -> SetBranchAddress("GenParticle", &Gen);
	  tree -> SetBranchAddress("Info",&ev_Info);
	  
	  for(int i0 =0; i0 <tree -> GetEntriesFast(); i0++) 
	  {
	    //1, event loop
	
	     tree -> GetEntry(i0);
	     bool leptonic = false;
	    
	    std::vector<fastjet::PseudoJet> PFparticles;
	    std::vector<fastjet::PseudoJet> GENparticles;
	    
	    nPU =  ev_Info -> nPU;
	    	    
	    for( int i1 = 0; i1 < PF -> GetEntriesFast(); i1++)
	    {
	      //2,entries loop,fill the vector particles with PF particles
	      
	      baconhep::TPFPart *PFTmp = (baconhep::TPFPart*)((*PF)[i1]);
	    
	      double Px=PFTmp->pt*cos(PFTmp->phi);
	      double Py= PFTmp->pt*sin(PFTmp->phi);
	      double theta = 2*atan(exp(-PFTmp->eta)); //eta = -ln(tan(theta/2))
	      double Pz = PFTmp->pt/tan(theta);
	      double E = PFTmp->e;
	      
	      fastjet::PseudoJet tmp_psjet(Px, Py, Pz, E);
	      PFparticles.push_back(tmp_psjet);
	    
	    }
	  
	   
	   //loop over genparticles
	    for( int i3 = 0; i3 < Gen -> GetEntriesFast(); i3++)
	    {
	      
	      baconhep::TGenParticle *GenParticleTmp = (baconhep::TGenParticle*)((*Gen)[i3]);
	      double Px=GenParticleTmp->pt*cos(GenParticleTmp->phi);
	      double Py= GenParticleTmp->pt*sin(GenParticleTmp->phi);
	      double theta = 2*atan(exp(-GenParticleTmp->eta)); //eta = -ln(tan(theta/2))
	      double Pz = GenParticleTmp->pt/tan(theta);
	      double Genmass = GenParticleTmp -> mass;
	      double E = sqrt(Px*Px + Py*Py + Pz*Pz + Genmass*Genmass);
	      fastjet::PseudoJet tmp_genpsjet(Px, Py, Pz, E);
	      
	    //  cout << "Number  "<< i3 <<   "|status " << GenParticleTmp -> status << "| pdgId " << GenParticleTmp -> pdgId << "| parent " << GenParticleTmp -> parent << endl;
	      
	      //check if W decays leptonically
	      if (fabs(GenParticleTmp -> pdgId) == 11 || fabs(GenParticleTmp -> pdgId) == 13 || fabs(GenParticleTmp -> pdgId) == 15) 
	      {
	      	  if ( fabs(((baconhep::TGenParticle*)((*Gen)[GenParticleTmp -> parent])) -> pdgId) == 24) leptonic = true;
	      }
	      
	      
	      if ((GenParticleTmp -> status) == 1 && fabs(GenParticleTmp -> pdgId)!=12 && fabs(GenParticleTmp -> pdgId)!=14 && fabs(GenParticleTmp -> pdgId)!=16 ) GENparticles.push_back(tmp_genpsjet);//taking particles only with status == 1 and filtering neutrinos
	  
		  
	    }
	    
	    if (leptonic) continue;//use only leptonically decays of W
	    fastjet::AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(4.0)));
	    JetDefinition jet_def(cambridge_algorithm, 0.8);
	    
	     
	    fastjet::ClusterSequenceArea clust_seq_PF(PFparticles, jet_def, area_def);
	    fastjet::ClusterSequenceArea clust_seq_GEN(GENparticles, jet_def, area_def);
	    
	    std::vector<fastjet::PseudoJet> PF_jets_basic = sorted_by_pt(clust_seq_PF.inclusive_jets(300.0));
	    std::vector<fastjet::PseudoJet> GEN_jets_basic = sorted_by_pt(clust_seq_GEN.inclusive_jets(300.0));
	    
	    fastjet::Selector rho_range =  SelectorAbsRapMax(4.0);
	    
	    JetMedianBackgroundEstimator bge_rho (rho_range, clust_seq_PF);
	    JetMedianBackgroundEstimator bge_rhom(rho_range, clust_seq_PF);
	    BackgroundJetPtMDensity m_density;
	    bge_rhom.set_jet_density_class(&m_density);
	    
	    contrib::SafeAreaSubtractor *area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom);
	     
	    // fastjet::JetMedianBackgroundEstimator bge_rho_PF(rho_range, clust_seq_PF);
	   //fastjet::JetMedianBackgroundEstimator bge_rhom_PF(rho_range, clust_seq_PF);
	    
	   //fastjet::JetMedianBackgroundEstimator bge_rho_GEN(rho_range, clust_seq_GEN);
	   //fastjet::JetMedianBackgroundEstimator bge_rhom_GEN(rho_range, clust_seq_GEN);
	   
	   // contrib::SafeAreaSubtractor *area_subtractor;
	   
	   double RCut= 0.5;
	   Pruner pruner(jet_def, 0.1,RCut);
	   Filter trimmer(JetDefinition(cambridge_algorithm, 0.3) ,SelectorPtFractionMin(0.05) );
	   Filter filter(JetDefinition(cambridge_algorithm, 0.3) , SelectorNHardest(3) );
	//   trimmer.set_subtractor(area_subtractor);
	  // filter.set_subtractor(area_subtractor);
	   
	   PseudoJet pruned_jet, trimmed_jet, filtered_jet, subtracted_jet;
	   
	   for (unsigned j =0; j < PF_jets_basic.size();++j)
	   {
	   
	    
	    trimmed_jet = (PseudoJet) trimmer.result(PF_jets_basic.at(j));
	    filtered_jet = (PseudoJet) filter.result(PF_jets_basic.at(j));
	    subtracted_jet = (PseudoJet) (*area_subtractor)(PF_jets_basic.at(j));
	    pruned_jet = pruner(subtracted_jet);
	        
	    mass_PFjet = PF_jets_basic.at(j).m();
	    mass_subtracted = subtracted_jet.m();
	    pt_PFjet = PF_jets_basic.at(j).pt();
	    mass_Pruned_PFjet = pruned_jet.m();
	    mass_Filtered_PFjet = filtered_jet.m();
	    mass_Trimmed_PFjet = trimmed_jet.m();
	    
	    PF_tree -> Fill();
	    
	  }

	   for (unsigned j =0; j < GEN_jets_basic.size();++j)
	   {
	   
	    pruned_jet = pruner(GEN_jets_basic.at(j));
	    trimmed_jet = (PseudoJet) trimmer.result(GEN_jets_basic.at(j));
	    filtered_jet = (PseudoJet) filter.result(GEN_jets_basic.at(j));
	        
	    mass_GENjet = GEN_jets_basic.at(j).m();
	    pt_GENjet = GEN_jets_basic.at(j).pt();
	    mass_Pruned_GENjet = pruned_jet.m();
	    mass_Filtered_GENjet = filtered_jet.m();
	    mass_Trimmed_GENjet = trimmed_jet.m();
	    
	    GEN_tree -> Fill();
	    
	  }
	
	  PF_jets_basic.clear();
	  GEN_jets_basic.clear();
	}//end of event loop
    }
    
   OutputFile.cd();
   GEN_tree -> Write();
   PF_tree -> Write();
}




