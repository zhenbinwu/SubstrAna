{
  gROOT->Reset();
  gROOT->ProcessLine(".L compareQuantities.C");

  TCanvas *c1236= new TCanvas("c1236","",200,10,800,700);
  
  // Jet Area
  compareQuantities("samples/QCD_CA8_combined.root", "samples/QCD_AK8_combined.root", "CHSjetArea", "CHSjetArea", 100, 0, 10);
  
  // Jet Mass
  compareQuantities("samples/ttbar_CA8_combined.root", "", "CHSjetMass", "CHSjetMassTrimmed", 100, 0, 500);
}