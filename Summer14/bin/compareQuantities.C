void compareQuantities(string filename, string var1, string var2, int nbins, float min, float max){


	TFile *inputFile = new TFile(filename.c_str());

	TTree *jetTree = (TTree *) inputFile->Get("JetTree");

	TH1F *var1H = new TH1F("var1H", "var1H", nbins, min, max);
	TH1F *var2H = new TH1F("var2H", "var2H", nbins, min, max);
	
	cout << "The tree has " << jetTree->GetEntries() << " entries." << endl;

	cout << var1 << "  " << var2 << endl;	

	jetTree->Draw(Form("%s>>var1H",var1.c_str()), "", "goff");
	jetTree->Draw(Form("%s>>var2H",var2.c_str()), "", "goff");

	cout << var1H->GetEntries() << " " << var2H->GetEntries() << endl;

	var1H->SetLineWidth(2);
	var2H->SetLineWidth(2);

	var1H->SetLineColor(kRed);
	var2H->SetLineColor(kBlue);


	//var1H->SetMaximum( TMath::Max(var1H->GetMaximum(), var2H->GetMaximum()) * 1.25 );

	


	var1H->DrawNormalized();
	var2H->DrawNormalized("same");
	

	TLegend *leg = new TLegend(0.1,0.9,0.9,0.99);
	leg->SetFillColor(0);
	leg->SetNColumns(2);
	leg->AddEntry(var1H, var1.c_str(), "l");
	leg->AddEntry(var2H, var2.c_str(), "l");
	leg->Draw("same");

}

