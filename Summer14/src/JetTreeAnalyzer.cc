#include "../include/JetTreeAnalyzer.h"

// --- constructor ------------------------------------------------------
JetTreeAnalyzer::JetTreeAnalyzer(TTree *tree){
  Init(tree);
}

// ----------------------------------------------------------------------
JetTreeAnalyzer::~JetTreeAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}


// --- Init tree --------------------------------------------------------
void JetTreeAnalyzer::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize                                                                                                   
  // a new tree or chain. Typically here the branch addresses and branch                                                    
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated     
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF 
  // (once per file to be processed).                                                                                                                                            
  
  pt = 0;
  ptcorr = 0;
  ptraw = 0;
  ptclean = 0;
  pttrim = 0;
  pttrimsafe = 0;
  ptconst = 0;
  ptunc = 0;
  eta = 0;
  phi = 0;
  m = 0;
  mraw = 0;
  mtrim = 0;
  mtrimsafe = 0;
  mclean = 0;
  mconst = 0;
  nparticles = 0;
  nneutrals = 0;
  ncharged = 0;
  ptgen = 0;
  etagen = 0;
  phigen = 0;
  mgen = 0;
  mrawgen = 0;
  mtrimgen = 0;
  mtrimsafegen = 0;
  mcleangen = 0;
  mconstgen = 0;
  imatch = 0;

  // Set branch addresses and branch pointers                                                                                                                                                                   
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("npu", &npu, &b_npu);
  fChain->SetBranchAddress("pt", &pt, &b_pt);
  fChain->SetBranchAddress("ptcorr", &ptcorr, &b_ptcorr);
  fChain->SetBranchAddress("ptraw", &ptraw, &b_ptraw);
  fChain->SetBranchAddress("ptclean", &ptclean, &b_ptclean);
  fChain->SetBranchAddress("pttrim", &pttrim, &b_pttrim);
  fChain->SetBranchAddress("pttrimsafe", &pttrimsafe, &b_pttrimsafe);
  fChain->SetBranchAddress("ptconst", &ptconst, &b_ptconst);
  fChain->SetBranchAddress("ptunc", &ptunc, &b_ptunc);
  fChain->SetBranchAddress("eta", &eta, &b_eta);
  fChain->SetBranchAddress("phi", &phi, &b_phi);
  fChain->SetBranchAddress("m", &m, &b_m);
  fChain->SetBranchAddress("mraw", &mraw, &b_mraw);
  fChain->SetBranchAddress("mtrim", &mtrim, &b_mtrim);
  fChain->SetBranchAddress("mtrimsafe", &mtrimsafe, &b_mtrimsafe);
  fChain->SetBranchAddress("mclean", &mclean, &b_mclean);
  fChain->SetBranchAddress("mconst", &mconst, &b_mconst);
  fChain->SetBranchAddress("nparticles", &nparticles, &b_nparticles);
  fChain->SetBranchAddress("nneutrals", &nneutrals, &b_nneutrals);
  fChain->SetBranchAddress("ncharged", &ncharged, &b_ncharged);
  fChain->SetBranchAddress("ptgen", &ptgen, &b_ptgen);
  fChain->SetBranchAddress("etagen", &etagen, &b_etagen);
  fChain->SetBranchAddress("phigen", &phigen, &b_phigen);
  fChain->SetBranchAddress("mgen", &mgen, &b_mgen);
  fChain->SetBranchAddress("mrawgen", &mrawgen, &b_mrawgen);
  fChain->SetBranchAddress("mtrimgen", &mtrimgen, &b_mtrimgen);
  fChain->SetBranchAddress("mtrimsafegen", &mtrimsafegen, &b_mtrimsafegen);
  fChain->SetBranchAddress("mcleangen", &mcleangen, &b_mcleangen);
  fChain->SetBranchAddress("mconstgen", &mconstgen, &b_mconstgen);
  fChain->SetBranchAddress("imatch", &imatch, &b_imatch);

}

// --- get Tree entry ----------------------------------------------------------------
Int_t JetTreeAnalyzer::GetEntry(Long64_t entry)
{
  // Read contents of entry.                                                                                                                                                                           
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}


// --- Book histograms ---------------------------------------------------------------
void JetTreeAnalyzer::bookHistograms(std::string suffix, float maxpt){

  std::cout << "Booking histograms for " << suffix.c_str() << " tree" << std::endl;

  hnjets = new TH1F(("hnjets"+suffix).c_str(), ("hnjets"+suffix).c_str(), 50, 0, 50 );


  // all jets
  hptgen = new TH1F(("hptgen"+suffix).c_str(), ("hptgen"+suffix).c_str(), maxpt, 0, maxpt );
  hptgen_pu = new TH1F(("hptgen_pu"+suffix).c_str(), ("hptgen_pu"+suffix).c_str(), maxpt, 0, maxpt );
  hptgen_good = new TH1F(("hptgen_good"+suffix).c_str(), ("hptgen_good"+suffix).c_str(), maxpt, 0, maxpt );

  hptraw = new TH1F(("hptraw"+suffix).c_str(), ("hptraw"+suffix).c_str(), maxpt, 0, maxpt );
  hptraw_pu = new TH1F(("hptraw_pu"+suffix).c_str(), ("hptraw_pu"+suffix).c_str(), maxpt, 0, maxpt );
  hptraw_good = new TH1F(("hptraw_good"+suffix).c_str(), ("hptraw_good"+suffix).c_str(), maxpt, 0, maxpt );
  hptraw_response = new TH1F(("hptraw_response"+suffix).c_str(), ("hptraw_response"+suffix).c_str(), 200, -100, 100 );

  hpt = new TH1F(("hpt"+suffix).c_str(), ("hpt"+suffix).c_str(), maxpt, 0, maxpt );
  hpt_pu = new TH1F(("hpt_pu"+suffix).c_str(), ("hpt_pu"+suffix).c_str(), maxpt, 0, maxpt );
  hpt_good = new TH1F(("hpt_good"+suffix).c_str(), ("hpt_good"+suffix).c_str(), maxpt, 0, maxpt );
  hpt_response = new TH1F(("hpt_response"+suffix).c_str(), ("hpt_response"+suffix).c_str(), 200, -100, 100 );

  hptcorr = new TH1F(("hptcorr"+suffix).c_str(), ("hptcorr"+suffix).c_str(), maxpt, 0, maxpt );
  hptcorr_pu = new TH1F(("hptcorr_pu"+suffix).c_str(), ("hptcorr_pu"+suffix).c_str(), maxpt, 0, maxpt );
  hptcorr_good = new TH1F(("hptcorr_good"+suffix).c_str(), ("hptcorr_good"+suffix).c_str(), maxpt, 0, maxpt );
  hptcorr_response = new TH1F(("hptcorr_response"+suffix).c_str(), ("hptcorr_response"+suffix).c_str(), 200, -100, 100 );

  heta = new TH1F(("heta"+suffix).c_str(), ("heta"+suffix).c_str(), 100, -5, 5 );
  heta_pu = new TH1F(("heta_pu"+suffix).c_str(), ("heta_pu"+suffix).c_str(), 100, -5, 5 );
  heta_good = new TH1F(("heta_good"+suffix).c_str(), ("heta_good"+suffix).c_str(), 100, -5, 5 );

  hnpu = new TH1F(("hnpu"+suffix).c_str(), ("hnpu"+suffix).c_str(), 100, 0, 100 );
  hnpu_pu = new TH1F(("hnpu_pu"+suffix).c_str(), ("hnpu_pu"+suffix).c_str(), 100, 0, 100 );
  hnpu_good = new TH1F(("hnpu_good"+suffix).c_str(), ("hnpu_good"+suffix).c_str(), 100, 0, 100 );

  hm = new TH1F(("hm"+suffix).c_str(), ("hm"+suffix).c_str(), 200, 0, 200 );
  hm_response = new TH1F(("hm_response"+suffix).c_str(), ("hm_response"+suffix).c_str(), 200, -100, 100 );

  hmraw = new TH1F(("hmraw"+suffix).c_str(), ("hmraw"+suffix).c_str(), 200, 0, 200 );
  hmraw_response = new TH1F(("hmraw_response"+suffix).c_str(), ("hmraw_response"+suffix).c_str(), 200, -100, 100 );
  
  hmtrim = new TH1F(("hmtrim"+suffix).c_str(), ("hmtrim"+suffix).c_str(), 200, 0, 200 );
  hmtrim_response = new TH1F(("hmtrim_response"+suffix).c_str(), ("hmtrim_response"+suffix).c_str(), 200, -100, 100 );

  hmtrimsafe = new TH1F(("hmtrimsafe"+suffix).c_str(), ("hmtrimsafe"+suffix).c_str(), 200, 0, 200 );
  hmtrimsafe_response = new TH1F(("hmtrimsafe_response"+suffix).c_str(), ("hmtrimsafe_response"+suffix).c_str(), 200, -100, 100 );

  hmclean = new TH1F(("hmclean"+suffix).c_str(), ("hmclean"+suffix).c_str(), 200, 0, 200 );
  hmclean_response = new TH1F(("hmclean_response"+suffix).c_str(), ("hmclean_response"+suffix).c_str(), 200, -100, 100 );

  hmconst = new TH1F(("hmconst"+suffix).c_str(), ("hmconst"+suffix).c_str(), 200, 0, 200 );
  hmconst_response = new TH1F(("hmconst_response"+suffix).c_str(), ("hmconst_response"+suffix).c_str(), 200, -100, 100 );

  hnparticles = new TH1F(("hnparticles"+suffix).c_str(), ("hnparticles"+suffix).c_str(), 1000, 0, 1000 );
  hnneutrals = new TH1F(("hnneutrals"+suffix).c_str(), ("hnneutrals"+suffix).c_str(), 1000, 0, 1000 );
  hncharged = new TH1F(("hncharged"+suffix).c_str(), ("hncharged"+suffix).c_str(), 1000, 0, 1000 );

  // leading jet 
  hptraw_leadjet = new TH1F(("hptraw_leadjet"+suffix).c_str(), ("hptraw_leadjet"+suffix).c_str(), maxpt, 0, maxpt );
  hptraw_pu_leadjet = new TH1F(("hptraw_pu_leadjet"+suffix).c_str(), ("hptraw_pu_leadjet"+suffix).c_str(), maxpt, 0, maxpt );
  hptraw_good_leadjet = new TH1F(("hptraw_good_leadjet"+suffix).c_str(), ("hptraw_good_leadjet"+suffix).c_str(), maxpt, 0, maxpt );
  hptraw_response_leadjet = new TH1F(("hptraw_response_leadjet"+suffix).c_str(), ("hptraw_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hpt_leadjet = new TH1F(("hpt_leadjet"+suffix).c_str(), ("hpt_leadjet"+suffix).c_str(), maxpt, 0, maxpt );
  hpt_pu_leadjet = new TH1F(("hpt_pu_leadjet"+suffix).c_str(), ("hpt_pu_leadjet"+suffix).c_str(), maxpt, 0, maxpt );
  hpt_good_leadjet = new TH1F(("hpt_good_leadjet"+suffix).c_str(), ("hpt_good_leadjet"+suffix).c_str(), maxpt, 0, maxpt );
  hpt_response_leadjet = new TH1F(("hpt_response_leadjet"+suffix).c_str(), ("hpt_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hptcorr_leadjet = new TH1F(("hptcorr_leadjet"+suffix).c_str(), ("hptcorr_leadjet"+suffix).c_str(), maxpt, 0, maxpt );
  hptcorr_pu_leadjet = new TH1F(("hptcorr_pu_leadjet"+suffix).c_str(), ("hptcorr_pu_leadjet"+suffix).c_str(), maxpt, 0, maxpt );
  hptcorr_good_leadjet = new TH1F(("hptcorr_good_leadjet"+suffix).c_str(), ("hptcorr_good_leadjet"+suffix).c_str(), maxpt, 0, maxpt );
  hptcorr_response_leadjet = new TH1F(("hptcorr_response_leadjet"+suffix).c_str(), ("hptcorr_response_leadjet"+suffix).c_str(), 200, -100, 100 );
  
  heta_leadjet = new TH1F(("heta_leadjet"+suffix).c_str(), ("heta_leadjet"+suffix).c_str(), 100, -5, 5 );
  heta_pu_leadjet = new TH1F(("heta_pu_leadjet"+suffix).c_str(), ("heta_pu_leadjet"+suffix).c_str(), 100, -5, 5 );
  heta_good_leadjet = new TH1F(("heta_good_leadjet"+suffix).c_str(), ("heta_good_leadjet"+suffix).c_str(), 100, -5, 5 );

  hmraw_leadjet = new TH1F(("hmraw_leadjet"+suffix).c_str(), ("hmraw_leadjet"+suffix).c_str(), 200, 0, 200 );
  hmraw_response_leadjet = new TH1F(("hmraw_response_leadjet"+suffix).c_str(), ("hmraw_response_leadjet"+suffix).c_str(), 200, -100, 100 );
  
  hm_leadjet = new TH1F(("hm_leadjet"+suffix).c_str(), ("hm_leadjet"+suffix).c_str(), 200, 0, 200 );
  hm_response_leadjet = new TH1F(("hm_response_leadjet"+suffix).c_str(), ("hm_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hmtrim_leadjet = new TH1F(("hmtrim_leadjet"+suffix).c_str(), ("hmtrim_leadjet"+suffix).c_str(), 200, 0, 200 );
  hmtrim_response_leadjet = new TH1F(("hmtrim_response_leadjet"+suffix).c_str(), ("hmtrim_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hmtrimsafe_leadjet = new TH1F(("hmtrimsafe_leadjet"+suffix).c_str(), ("hmtrimsafe_leadjet"+suffix).c_str(), 200, 0, 200 );
  hmtrimsafe_response_leadjet = new TH1F(("hmtrimsafe_response_leadjet"+suffix).c_str(), ("hmtrimsafe_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hmclean_leadjet = new TH1F(("hmclean_leadjet"+suffix).c_str(), ("hmclean_leadjet"+suffix).c_str(), 200, 0, 200 );
  hmclean_response_leadjet = new TH1F(("hmclean_response_leadjet"+suffix).c_str(), ("hmclean_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hmconst_leadjet = new TH1F(("hmconst_leadjet"+suffix).c_str(), ("hmconst_leadjet"+suffix).c_str(), 200, 0, 200 );
  hmconst_response_leadjet = new TH1F(("hmconst_response_leadjet"+suffix).c_str(), ("hmconst_response_leadjet"+suffix).c_str(), 200, -100, 100 );

  hnparticles_leadjet = new TH1F(("hnparticles_leadjet "+suffix).c_str(), ("hnparticles_leadjet "+suffix).c_str(), 100, 0, 100 );
  hnneutrals_leadjet  = new TH1F(("hnneutrals_leadjet "+suffix).c_str(), ("hnneutrals_leadjet "+suffix).c_str(), 100, 0, 100 );
  hncharged_leadjet   = new TH1F(("hncharged_leadjet "+suffix).c_str(), ("hncharged_leadjet "+suffix).c_str(), 100, 0, 100 );

  // 2d histograms

  hptraw_response_vs_pt     = new TH2F(("hptraw_response_vs_pt"+suffix).c_str(), ("hptraw_response_vs_pt"+suffix).c_str(), maxpt, 0, maxpt, 200, -100, 100 );
  hpt_response_vs_pt        = new TH2F(("hpt_response_vs_pt"+suffix).c_str(), ("hpt_response_vs_pt"+suffix).c_str(), maxpt, 0, maxpt, 200, -100, 100 );
  hmraw_response_vs_pt      = new TH2F(("hmraw_response_vs_pt"+suffix).c_str(), ("hmraw_response_vs_pt"+suffix).c_str(), maxpt, 0, maxpt, 200, -100, 100 );
  hm_response_vs_pt         = new TH2F(("hm_response_vs_pt"+suffix).c_str(), ("hm_response_vs_pt"+suffix).c_str(), maxpt, 0, maxpt, 200, -100, 100 );
  hmtrim_response_vs_pt     = new TH2F(("hmtrim_response_vs_pt"+suffix).c_str(), ("hmtrim_response_vs_pt"+suffix).c_str(), maxpt, 0, maxpt, 200, -100, 100 );
  hmtrimsafe_response_vs_pt = new TH2F(("hmtrimsafe_response_vs_pt"+suffix).c_str(), ("hmtrimsafe_response_vs_pt"+suffix).c_str(), maxpt, 0, maxpt, 200, -100, 100 );
  hmclean_response_vs_pt    = new TH2F(("hmclean_response_vs_pt"+suffix).c_str(), ("hmclean_response_vs_pt"+suffix).c_str(), maxpt, 0, maxpt, 200, -100, 100 );
  hmconst_response_vs_pt    = new TH2F(("hmconst_response_vs_pt"+suffix).c_str(), ("hmconst_response_vs_pt"+suffix).c_str(), maxpt, 0, maxpt, 200, -100, 100 );

  hptraw_response_vs_eta = new TH2F(("hptraw_response_vs_eta"+suffix).c_str(), ("hptraw_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hpt_response_vs_eta = new TH2F(("hpt_response_vs_eta"+suffix).c_str(), ("hpt_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmraw_response_vs_eta = new TH2F(("hmraw_response_vs_eta"+suffix).c_str(), ("hmraw_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hm_response_vs_eta = new TH2F(("hm_response_vs_eta"+suffix).c_str(), ("hm_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmtrim_response_vs_eta = new TH2F(("hmtrim_response_vs_eta"+suffix).c_str(), ("hmtrim_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmtrimsafe_response_vs_eta = new TH2F(("hmtrimsafe_response_vs_eta"+suffix).c_str(), ("hmtrimsafe_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmclean_response_vs_eta = new TH2F(("hmclean_response_vs_eta"+suffix).c_str(), ("hmclean_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );
  hmconst_response_vs_eta = new TH2F(("hmconst_response_vs_eta"+suffix).c_str(), ("hmconst_response_vs_eta"+suffix).c_str(), 100, -5, 5, 200, -100, 100 );

  hptraw_response_vs_npu = new TH2F(("hptraw_response_vs_npu"+suffix).c_str(), ("hptraw_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hpt_response_vs_npu = new TH2F(("hpt_response_vs_npu"+suffix).c_str(), ("hpt_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmraw_response_vs_npu = new TH2F(("hmraw_response_vs_npu"+suffix).c_str(), ("hmraw_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hm_response_vs_npu = new TH2F(("hm_response_vs_npu"+suffix).c_str(), ("hm_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmtrim_response_vs_npu = new TH2F(("hmtrim_response_vs_npu"+suffix).c_str(), ("hmtrim_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmtrimsafe_response_vs_npu = new TH2F(("hmtrimsafe_response_vs_npu"+suffix).c_str(), ("hmtrimsafe_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmclean_response_vs_npu = new TH2F(("hmclean_response_vs_npu"+suffix).c_str(), ("hmclean_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );
  hmconst_response_vs_npu = new TH2F(("hmconst_response_vs_npu"+suffix).c_str(), ("hmconst_response_vs_npu"+suffix).c_str(), 100, 0, 100, 200, -100, 100 );

}


// --- Fill histograms ---------------------------------------------------------------
void JetTreeAnalyzer::fillHistograms(int maxEntries, float minPt){

  std::cout << "Filling histograms..." << std::endl;
  
  if (fChain==0){
    std::cout<<"Error: cannot open " << fChain->GetName() << std::endl;
    exit(0);
  }

  if (maxEntries == -1)
    maxEntries = fChain->GetEntries();

  
  for (int entry = 0; entry < maxEntries; entry++){
    fChain->GetEntry(entry);

    if (entry%100==0) std::cout << "Analyzing entry : " << entry << "\r" << std::flush;
        
    // --- Loop over jets in this event                                                                                                                                      
    int nj = 0;

    for (unsigned int j = 0; j < ptraw->size(); j++){

      float thispt = pt->at(j); // use pt 
      
      if (thispt < minPt)  continue;

      nj++;
      int matchInd = imatch->at(j);

      
      hptgen    -> Fill(ptgen->at(j));
      hptraw    -> Fill(ptraw->at(j));
      hpt       -> Fill(pt->at(j));
      hptcorr   -> Fill(ptcorr->at(j));
      heta      -> Fill(eta->at(j));
      hnpu      -> Fill(npu);
      hmraw     -> Fill(mraw->at(j));
      hm        -> Fill(m->at(j));
      hmtrim    -> Fill(mtrim->at(j));
      hmtrimsafe-> Fill(mtrimsafe->at(j));
      hmclean   -> Fill(mclean->at(j));
      hmconst   -> Fill(mconst->at(j));

      hnparticles -> Fill(nparticles->at(j));
      hnneutrals  -> Fill(nneutrals->at(j));
      hncharged   -> Fill(ncharged->at(j));

      
      if (j == 0){
	hptraw_leadjet    -> Fill(ptraw->at(j));
	hpt_leadjet       -> Fill(pt->at(j));
	hptcorr_leadjet   -> Fill(ptcorr->at(j));
	heta_leadjet      -> Fill(eta->at(j));
	hmraw_leadjet     -> Fill(mraw->at(j));
	hm_leadjet        -> Fill(m->at(j));
	hmtrim_leadjet    -> Fill(mtrim->at(j));
	hmtrimsafe_leadjet-> Fill(mtrimsafe->at(j));
	hmclean_leadjet   -> Fill(mclean->at(j));
	hmconst_leadjet   -> Fill(mconst->at(j));
      }

      if (matchInd == -1 ) {
	hptgen_pu -> Fill(ptgen->at(j));
	hptraw_pu -> Fill(ptraw->at(j));
	hpt_pu    -> Fill(pt->at(j));
	hptcorr_pu-> Fill(ptcorr->at(j));
	heta_pu   -> Fill(eta->at(j));
	hnpu_pu   -> Fill(npu);
	if (j == 0){
	  hptraw_pu_leadjet -> Fill(ptraw->at(j));
	  hpt_pu_leadjet    -> Fill(pt->at(j));
	  hptcorr_pu_leadjet-> Fill(ptcorr->at(j));
	  heta_pu_leadjet   -> Fill(eta->at(j));
	}
      }
      else {
	hptgen_good -> Fill(ptgen->at(j));
	hptraw_good -> Fill(ptraw->at(j));
	hpt_good    -> Fill(pt->at(j));
	hptcorr_good-> Fill(ptcorr->at(j));
	heta_good   -> Fill(eta->at(j));
	hnpu_good   -> Fill(npu);
	if (j == 0){
	  hptraw_good_leadjet -> Fill(ptraw->at(j));
	  hpt_good_leadjet    -> Fill(pt->at(j));
	  hptcorr_good_leadjet-> Fill(ptcorr->at(j));
	  heta_good_leadjet   -> Fill(eta->at(j));
	}
      }


      // -- response plots
      if (matchInd >- 1){
	hptraw_response     -> Fill(ptraw->at(j)-ptgen->at(j));
	hpt_response        -> Fill(pt->at(j)-ptgen->at(j));
	hptcorr_response    -> Fill(ptcorr->at(j)-ptgen->at(j));
	hmraw_response      -> Fill(mraw->at(j)-mrawgen->at(j));
	hm_response         -> Fill(m->at(j)-mgen->at(j));
	hmtrim_response     -> Fill(mtrim->at(j)-mtrimgen->at(j));
	hmtrimsafe_response -> Fill(mtrimsafe->at(j)-mtrimsafegen->at(j));
	hmclean_response    -> Fill(mclean->at(j)-mcleangen->at(j));
	hmconst_response    -> Fill(mconst->at(j)-mconstgen->at(j));
	
	hptraw_response_vs_pt     -> Fill(pt->at(j),ptraw->at(j)-ptgen->at(j));
	hpt_response_vs_pt        -> Fill(pt->at(j),pt->at(j)-ptgen->at(j));
	hmraw_response_vs_pt      -> Fill(pt->at(j),mraw->at(j)-mrawgen->at(j));
	hm_response_vs_pt         -> Fill(pt->at(j),m->at(j)-mgen->at(j));
	hmtrim_response_vs_pt     -> Fill(pt->at(j),mtrim->at(j)-mtrimgen->at(j));
	hmtrimsafe_response_vs_pt -> Fill(pt->at(j),mtrimsafe->at(j)-mtrimsafegen->at(j));
	hmclean_response_vs_pt    -> Fill(pt->at(j),mclean->at(j)-mcleangen->at(j));
	hmconst_response_vs_pt    -> Fill(pt->at(j),mconst->at(j)-mconstgen->at(j));

	hptraw_response_vs_eta     -> Fill(eta->at(j),ptraw->at(j)-ptgen->at(j));
	hpt_response_vs_eta        -> Fill(eta->at(j),pt->at(j)-ptgen->at(j));
	hmraw_response_vs_eta      -> Fill(eta->at(j),mraw->at(j)-mrawgen->at(j));
	hm_response_vs_eta         -> Fill(eta->at(j),m->at(j)-mgen->at(j));
	hmtrim_response_vs_eta     -> Fill(eta->at(j),mtrim->at(j)-mtrimgen->at(j));
	hmtrimsafe_response_vs_eta -> Fill(eta->at(j),mtrimsafe->at(j)-mtrimsafegen->at(j));
	hmclean_response_vs_eta    -> Fill(eta->at(j),mclean->at(j)-mcleangen->at(j));
	hmconst_response_vs_eta    -> Fill(eta->at(j),mconst->at(j)-mconstgen->at(j));

	hptraw_response_vs_npu     -> Fill(npu,ptraw->at(j)-ptgen->at(j));
	hpt_response_vs_npu        -> Fill(npu,pt->at(j)-ptgen->at(j));
	hmraw_response_vs_npu      -> Fill(npu,mraw->at(j)-mrawgen->at(j));
	hm_response_vs_npu         -> Fill(npu,m->at(j)-mgen->at(j));
	hmtrim_response_vs_npu     -> Fill(npu,mtrim->at(j)-mtrimgen->at(j));
	hmtrimsafe_response_vs_npu -> Fill(npu,mtrimsafe->at(j)-mtrimsafegen->at(j));
	hmclean_response_vs_npu    -> Fill(npu,mclean->at(j)-mcleangen->at(j));
	hmconst_response_vs_npu    -> Fill(npu,mconst->at(j)-mconstgen->at(j));


	if (j == 0){
	  hptraw_response_leadjet     -> Fill(ptraw->at(j)-ptgen->at(j));
	  hpt_response_leadjet        -> Fill(pt->at(j)-ptgen->at(j));
	  hptcorr_response_leadjet    -> Fill(ptcorr->at(j)-ptgen->at(j));
	  hmraw_response_leadjet      -> Fill(mraw->at(j)-mrawgen->at(j));
	  hm_response_leadjet         -> Fill(m->at(j)-mgen->at(j));
	  hmtrim_response_leadjet     -> Fill(mtrim->at(j)-mtrimgen->at(j));
	  hmtrimsafe_response_leadjet -> Fill(mtrimsafe->at(j)-mtrimsafegen->at(j));
	  hmclean_response_leadjet    -> Fill(mclean->at(j)-mcleangen->at(j));
	  hmconst_response_leadjet    -> Fill(mconst->at(j)-mconstgen->at(j));
	}
      }
 
    }// end loop over jets 

    hnjets->Fill(nj);

  }// end loop over entries
}


// --- Save histograms ---------------------------------------------------------------
void JetTreeAnalyzer::saveHistograms(TFile *outfile, std::string dir){

  std::cout << "Saving histograms ... " << std::endl;
  
  outfile->cd();
  TDirectory *thisdir = outfile->mkdir(dir.c_str());
  thisdir->cd();    // make the "thisdir" directory
  
  hnjets->Write();

  hptgen->Write();
  hptgen_pu->Write();
  hptgen_good->Write();

  hptraw->Write();
  hptraw_pu->Write();
  hptraw_good->Write();
  hptraw_response->Write();

  hpt->Write();
  hpt_pu->Write();
  hpt_good->Write();
  hpt_response->Write();

  hptcorr->Write();
  hptcorr_pu->Write();
  hptcorr_good->Write();
  hptcorr_response->Write();

  heta->Write();
  heta_pu->Write();
  heta_good->Write();

  hnpu->Write();
  hnpu_pu->Write();
  hnpu_good->Write();

  hnparticles->Write();
  hnneutrals->Write();
  hncharged->Write();

  hmraw->Write();
  hmraw_response->Write();
  hm->Write();
  hm_response->Write();
  hmtrim->Write();
  hmtrim_response->Write();
  hmtrimsafe->Write();
  hmtrimsafe_response->Write();
  hmclean->Write();
  hmclean_response->Write();
  hmconst->Write();
  hmconst_response->Write();

  // leading jet
  hptraw_leadjet->Write();
  hptraw_pu_leadjet->Write();
  hptraw_good_leadjet->Write();
  hptraw_response_leadjet->Write();

  hpt_leadjet->Write();
  hpt_pu_leadjet->Write();
  hpt_good_leadjet->Write();
  hpt_response_leadjet->Write();
  hptcorr_leadjet->Write();
  hptcorr_pu_leadjet->Write();
  hptcorr_good_leadjet->Write();
  hptcorr_response_leadjet->Write();

  heta_leadjet->Write();
  heta_pu_leadjet->Write();
  heta_good_leadjet->Write();

  hnparticles_leadjet->Write();
  hnneutrals_leadjet->Write();
  hncharged_leadjet->Write();

  hmraw_leadjet->Write();
  hmraw_response_leadjet->Write();
  hm_leadjet->Write();
  hm_response_leadjet->Write();
  hmtrim_leadjet->Write();
  hmtrim_response_leadjet->Write();
  hmtrimsafe_leadjet->Write();
  hmtrimsafe_response_leadjet->Write();
  hmclean_leadjet->Write();
  hmclean_response_leadjet->Write();
  hmconst_leadjet->Write();
  hmconst_response_leadjet->Write();

  // 2d
  hptraw_response_vs_pt->Write();
  hpt_response_vs_pt->Write();
  hmraw_response_vs_pt->Write();
  hm_response_vs_pt->Write();
  hmtrim_response_vs_pt->Write();
  hmtrimsafe_response_vs_pt->Write();
  hmclean_response_vs_pt->Write();
  hmconst_response_vs_pt->Write();

  hptraw_response_vs_eta->Write();
  hpt_response_vs_eta->Write();
  hm_response_vs_eta->Write();
  hmraw_response_vs_eta->Write();
  hmtrim_response_vs_eta->Write();
  hmtrimsafe_response_vs_eta->Write();
  hmclean_response_vs_eta->Write();
  hmconst_response_vs_eta->Write();

  hptraw_response_vs_npu->Write();
  hpt_response_vs_npu->Write();
  hmraw_response_vs_npu->Write();
  hm_response_vs_npu->Write();
  hmtrim_response_vs_npu->Write();
  hmtrimsafe_response_vs_npu->Write();
  hmclean_response_vs_npu->Write();
  hmconst_response_vs_npu->Write();

}
