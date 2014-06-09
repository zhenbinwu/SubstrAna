#include "../src/JetTreeAnalyzer.cc"


//-------------------------------------------------------                                                                                                              
// MAIN                                                                                                                                                                       
//-------------------------------------------------------                                                                                                                                  
  
int main( int argc, char **argv ) {

  gROOT->ProcessLine("#include <vector>");

  TFile *inputFile = TFile::Open(argv[1]);
  if (inputFile==0){
    std::cout<<"Error: cannot open " << inputFile->GetName() << std::endl;
    exit(0);
  }

  std::string outname = argv[2];

<<<<<<< HEAD
  int maxEntries = -1;
  float minpt = atof(argv[3]);
  float maxpt  = atof(argv[4]);
  bool doCMSSWJets = atoi(argv[5]);

  // -- gen
  TTree *tree_gen   = (TTree *)inputFile->Get("gen");
  JetTreeAnalyzer *genAnalyzer = new JetTreeAnalyzer(tree_gen);
  genAnalyzer->bookHistograms("_gen", maxpt);
  genAnalyzer->fillHistograms(maxEntries,minpt);
  delete tree_gen;

  // -- pf
  TTree *tree_pf    = (TTree *)inputFile->Get("pf");
  JetTreeAnalyzer *pfAnalyzer = new JetTreeAnalyzer(tree_pf);
  pfAnalyzer->bookHistograms("_pf", maxpt);
  pfAnalyzer->fillHistograms(maxEntries,minpt);
  delete tree_pf;

  // -- pfchs
  TTree *tree_pfchs = (TTree *)inputFile->Get("chs");
=======
  TTree *tree_gen   = (TTree *)inputFile->Get("gen");
  TTree *tree_pf    = (TTree *)inputFile->Get("pf");
  TTree *tree_pfchs = (TTree *)inputFile->Get("chs");
  TTree *tree_puppi = (TTree *)inputFile->Get("puppi");
  //TTree *tree_pfcmssw = (TTree *)inputFile->Get("cmsswpf");

  int maxEntries = -1;
  float minpt = atof(argv[3]);
  float maxpt  = atof(argv[4]);

  JetTreeAnalyzer *genAnalyzer = new JetTreeAnalyzer(tree_gen);
  genAnalyzer->bookHistograms("_gen", maxpt);
  genAnalyzer->fillHistograms(maxEntries,minpt);

  JetTreeAnalyzer *pfAnalyzer = new JetTreeAnalyzer(tree_pf);
  pfAnalyzer->bookHistograms("_pf", maxpt);
  pfAnalyzer->fillHistograms(maxEntries,minpt);
  
>>>>>>> 6bfe8398e4508d4b6164356e713476c2c4272a93
  JetTreeAnalyzer *pfchsAnalyzer = new JetTreeAnalyzer(tree_pfchs);
  pfchsAnalyzer->bookHistograms("_pfchs", maxpt);
  pfchsAnalyzer->fillHistograms(maxEntries,minpt);

<<<<<<< HEAD
  // -- puppi
  TTree *tree_puppi = (TTree *)inputFile->Get("puppi");
  JetTreeAnalyzer *puppiAnalyzer = new JetTreeAnalyzer(tree_puppi);
  puppiAnalyzer->bookHistograms("_puppi", maxpt);
  puppiAnalyzer->fillHistograms(maxEntries,minpt);
  delete tree_puppi;

  // -- pf cmssw
  TTree *tree_pfcmssw;
  JetTreeAnalyzer *pfcmsswAnalyzer = 0;
  if (doCMSSWJets){
    tree_pfcmssw = (TTree *)inputFile->Get("cmsswpf");
    pfcmsswAnalyzer = new JetTreeAnalyzer(tree_pfcmssw);
    pfcmsswAnalyzer->bookHistograms("_pfcmssw", maxpt);
    pfcmsswAnalyzer->fillHistograms(maxEntries, minpt);
    delete tree_pfcmssw;
  }

=======
  JetTreeAnalyzer *puppiAnalyzer = new JetTreeAnalyzer(tree_puppi);
  puppiAnalyzer->bookHistograms("_puppi", maxpt);
  puppiAnalyzer->fillHistograms(maxEntries,minpt);

  //JetTreeAnalyzer *pfcmsswAnalyzer = new JetTreeAnalyzer(tree_pfcmssw);
  //pfcmsswAnalyzer->bookHistograms("_pfcmssw", maxpt);
  //pfcmsswAnalyzer->fillHistograms(maxEntries, minpt);
  
>>>>>>> 6bfe8398e4508d4b6164356e713476c2c4272a93
  // save results in file
  TFile *outfile = new TFile(outname.c_str(),"RECREATE");
  genAnalyzer->saveHistograms(outfile,"gen");
  pfAnalyzer->saveHistograms(outfile,"pf");
  pfchsAnalyzer->saveHistograms(outfile,"pfchs");
  puppiAnalyzer->saveHistograms(outfile,"puppi");
<<<<<<< HEAD
  if (doCMSSWJets) pfcmsswAnalyzer->saveHistograms(outfile,"pfcmssw");
=======
  //pfcmsswAnalyzer->saveHistograms(outfile,"pfcmssw");
>>>>>>> 6bfe8398e4508d4b6164356e713476c2c4272a93
  

}
