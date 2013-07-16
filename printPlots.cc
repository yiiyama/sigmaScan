void printPlots(TString const& point = "600_1020")
{
  // EXAMPLE ANALYSIS MACRO

  TFile* source = TFile::Open("/afs/cern.ch/user/y/yiiyama/public/sigmaScan/recoMet_ph40_el30_mt100.root");

  TTree* sigmaTree = (TTree*)source->Get("sigmaTree");

  TString* gridParams = new TString;
  float significance;
  float cut_recoMet;

  sigmaTree->SetBranchAddress("gridParams", &gridParams);
  sigmaTree->SetBranchAddress("significance", &significance);
  sigmaTree->SetBranchAddress("cut_recoMet", &cut_recoMet);

  TH2F* h2_cut = new TH2F("h2_cut", "#slash{E}_{T} cut value;Squark mass (GeV);Gluino mass (GeV);Cut value (GeV)", 17, 350., 2050., 17, 350., 2050.);
  TH2F* h2_significance = new TH2F("h2_significance", "S/#sqrt{S+B} (20/fb);Squark mass (GeV);Gluino mass (GeV)", 17, 350., 2050., 17, 350., 2050.);

  TPRegexp gridParamPat("^([0-9]+)_([0-9]+)$");

  long iPoint(-1);

  long iEntry = 0;
  while(sigmaTree->GetEntry(iEntry++) != 0){
    TObjArray* matches = gridParamPat.MatchS(*gridParams);

    if(matches->GetEntries() == 0){
      delete matches;
      continue;
    }

    TString msqStr = matches->At(1)->GetName();
    double msq = msqStr.Atof();

    TString mglStr = matches->At(2)->GetName();
    double mgl = mglStr.Atof();

    delete matches;

    if(*gridParams == point) iPoint = iEntry - 1;

    int bin = h2_cut->FindBin(msq, mgl);

    h2_cut->SetBinContent(bin, cut_recoMet);
    h2_significance->SetBinContent(bin, significance);
  }

  TCanvas* cgrid = new TCanvas("cgrid", "cgrid", 600, 800);
  cgrid->Divide(1, 2);

  cgrid->cd(1);
  h2_cut->Draw("colz");

  cgrid->cd(2);
  h2_significance->Draw("colz");

  if(iPoint != -1){
    TCanvas* cdist = new TCanvas("cdist", "cdist", 600, 400);
    sigmaTree->SetBranchAddress("dist_recoMet", &cdist);
    sigmaTree->GetEntry(iPoint);
    cdist->SetLogy();

    cdist->Draw();
  }
}
