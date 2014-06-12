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

  int maxEntries = -1;
  float minpt = atof(argv[3]);
  float maxpt  = atof(argv[4]);
  bool doCMSSWJets = atoi(argv[5]);

  // -- gen
  TTree *tree_gen   = (TTree *)inputFile->Get("gen");
  JetTreeAnalyzer *genAnalyzer = new JetTreeAnalyzer(tree_gen);
  genAnalyzer->bookHistograms("_gen");
  genAnalyzer->fillHistograms(maxEntries,minpt, maxpt);
  delete tree_gen;

  // -- pf
  TTree *tree_pf    = (TTree *)inputFile->Get("pf");
  JetTreeAnalyzer *pfAnalyzer = new JetTreeAnalyzer(tree_pf);
  pfAnalyzer->bookHistograms("_pf");
  pfAnalyzer->fillHistograms(maxEntries,minpt,maxpt);
  delete tree_pf;

  // -- pfchs
  TTree *tree_pfchs = (TTree *)inputFile->Get("chs");
  JetTreeAnalyzer *pfchsAnalyzer = new JetTreeAnalyzer(tree_pfchs);
  pfchsAnalyzer->bookHistograms("_pfchs");
  pfchsAnalyzer->fillHistograms(maxEntries,minpt,maxpt);

  // -- puppi
  TTree *tree_puppi = (TTree *)inputFile->Get("puppi");
  JetTreeAnalyzer *puppiAnalyzer = new JetTreeAnalyzer(tree_puppi);
  puppiAnalyzer->bookHistograms("_puppi");
  puppiAnalyzer->fillHistograms(maxEntries,minpt,maxpt);
  delete tree_puppi;

  // -- pf cmssw
  TTree *tree_pfcmssw;
  JetTreeAnalyzer *pfcmsswAnalyzer = 0;
  if (doCMSSWJets){
    tree_pfcmssw = (TTree *)inputFile->Get("cmsswpf");
    pfcmsswAnalyzer = new JetTreeAnalyzer(tree_pfcmssw);
    pfcmsswAnalyzer->bookHistograms("_pfcmssw");
    pfcmsswAnalyzer->fillHistograms(maxEntries, minpt,maxpt);
    delete tree_pfcmssw;
  }

  // save results in file
  TFile *outfile = new TFile(outname.c_str(),"RECREATE");
  genAnalyzer->saveHistograms(outfile,"gen");
  pfAnalyzer->saveHistograms(outfile,"pf");
  pfchsAnalyzer->saveHistograms(outfile,"pfchs");
  puppiAnalyzer->saveHistograms(outfile,"puppi");
  if (doCMSSWJets) pfcmsswAnalyzer->saveHistograms(outfile,"pfcmssw");
  

}
