MatrixMethod_Results_Efficiencies(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile* file_MM_MC_E = new TFile("MatrixMethod_Efficiency_MC_E.root");
  file_MM_MC_E->cd();
  hSignalEfficiency_MC_E_NJets = (TH1F*)gDirectory->Get("SignalEfficiencyNJets");
  hSignalEfficiency_MC_E_Eta = (TH1F*)gDirectory->Get("SignalEfficiencyEta");
  hSignalEfficiency_MC_E_PT = (TH1F*)gDirectory->Get("SignalEfficiencyPT");
  hFakeRate_MC_E_NJets = (TH1F*)gDirectory->Get("FakeRateNJets");
  hFakeRate_MC_E_Eta = (TH1F*)gDirectory->Get("FakeRateEta");
  hFakeRate_MC_E_PT = (TH1F*)gDirectory->Get("FakeRatePT");

  TFile* file_MM_DATA_E = new TFile("MatrixMethod_Efficiency_DATA_E.root");
  file_MM_DATA_E->cd();
  hSignalEfficiency_DATA_E_NJets = (TH1F*)gDirectory->Get("SignalEfficiencyNJets");
  hSignalEfficiency_DATA_E_Eta = (TH1F*)gDirectory->Get("SignalEfficiencyEta");
  hSignalEfficiency_DATA_E_PT = (TH1F*)gDirectory->Get("SignalEfficiencyPT");
  hFakeRate_DATA_E_NJets = (TH1F*)gDirectory->Get("FakeRateNJets");
  hFakeRate_DATA_E_NJets_10 = (TH1F*)gDirectory->Get("FakeRateNJets_10");
  hFakeRate_DATA_E_NJets_20 = (TH1F*)gDirectory->Get("FakeRateNJets_20");
  hFakeRate_DATA_E_NJets_30 = (TH1F*)gDirectory->Get("FakeRateNJets_30");
  hFakeRate_DATA_E_NJets_40 = (TH1F*)gDirectory->Get("FakeRateNJets_40");
  hFakeRate_DATA_E_NJets_50 = (TH1F*)gDirectory->Get("FakeRateNJets_50");
  hFakeRate_DATA_E_NJets_60 = (TH1F*)gDirectory->Get("FakeRateNJets_60");
  hFakeRate_DATA_E_NJets_70 = (TH1F*)gDirectory->Get("FakeRateNJets_70");
  hFakeRate_DATA_E_NJets_80 = (TH1F*)gDirectory->Get("FakeRateNJets_80");
  hFakeRate_DATA_E_NJets_90 = (TH1F*)gDirectory->Get("FakeRateNJets_90");
  hFakeRate_DATA_E_Eta = (TH1F*)gDirectory->Get("FakeRateEta");
  hFakeRate_DATA_E_PT = (TH1F*)gDirectory->Get("FakeRatePT");
  hFakeRate_DATA_E_NVTX = (TH1F*)gDirectory->Get("FakeRateNVTX");




  TFile* file_MM_MC_Mu = new TFile("MatrixMethod_Efficiency_MC_Mu.root");
  file_MM_MC_Mu->cd();
  hSignalEfficiency_MC_Mu_NJets = (TH1F*)gDirectory->Get("SignalEfficiencyNJets");
  hSignalEfficiency_MC_Mu_Eta = (TH1F*)gDirectory->Get("SignalEfficiencyEta");
  hSignalEfficiency_MC_Mu_PT = (TH1F*)gDirectory->Get("SignalEfficiencyPT");
  hFakeRate_MC_Mu_NJets = (TH1F*)gDirectory->Get("FakeRateNJets");
  hFakeRate_MC_Mu_Eta = (TH1F*)gDirectory->Get("FakeRateEta");
  hFakeRate_MC_Mu_PT = (TH1F*)gDirectory->Get("FakeRatePT");

  TFile* file_MM_DATA_Mu_ = new TFile("MatrixMethod_Efficiency_DATA_Mu.root");
  file_MM_DATA_Mu_->cd();
  hSignalEfficiency_DATA_Mu_NJets = (TH1F*)gDirectory->Get("SignalEfficiencyNJets");
  hSignalEfficiency_DATA_Mu_Eta = (TH1F*)gDirectory->Get("SignalEfficiencyEta");
  hSignalEfficiency_DATA_Mu_PT = (TH1F*)gDirectory->Get("SignalEfficiencyPT");
  hFakeRate_DATA_Mu_NJets = (TH1F*)gDirectory->Get("FakeRateNJets");
  hFakeRate_DATA_Mu_NJets_10 = (TH1F*)gDirectory->Get("FakeRateNJets_10");
  hFakeRate_DATA_Mu_NJets_20 = (TH1F*)gDirectory->Get("FakeRateNJets_20");
  hFakeRate_DATA_Mu_NJets_30 = (TH1F*)gDirectory->Get("FakeRateNJets_30");
  hFakeRate_DATA_Mu_NJets_40 = (TH1F*)gDirectory->Get("FakeRateNJets_40");
  hFakeRate_DATA_Mu_NJets_50 = (TH1F*)gDirectory->Get("FakeRateNJets_50");
  hFakeRate_DATA_Mu_NJets_60 = (TH1F*)gDirectory->Get("FakeRateNJets_60");
  hFakeRate_DATA_Mu_NJets_70 = (TH1F*)gDirectory->Get("FakeRateNJets_70");
  hFakeRate_DATA_Mu_NJets_80 = (TH1F*)gDirectory->Get("FakeRateNJets_80");
  hFakeRate_DATA_Mu_NJets_90 = (TH1F*)gDirectory->Get("FakeRateNJets_90");
  hFakeRate_DATA_Mu_Eta = (TH1F*)gDirectory->Get("FakeRateEta");
  hFakeRate_DATA_Mu_PT = (TH1F*)gDirectory->Get("FakeRatePT");
  hFakeRate_DATA_Mu_NVTX = (TH1F*)gDirectory->Get("FakeRateNVTX");



  title = new TLatex(10.,10.,"CMS preliminary");
  title->SetNDC();
  title->SetTextAlign(12);
  title->SetX(0.35);
  title->SetY(0.70);
  title->SetTextFont(42);
  title->SetTextSize(0.04);
  title->SetTextSizePixels(24);


  
  TLegend* qw = 0;
  qw = new TLegend(0.65,0.99,0.93,0.91);
  qw->SetFillColor(kWhite);

  hSignalEfficiency_MC_E_NJets->GetXaxis()->SetTitle("Jet multiplicity");
  hSignalEfficiency_MC_E_NJets->GetYaxis()->SetTitle("Efficiency");
  hSignalEfficiency_MC_E_NJets->SetTitle("#varepsilon_{s}^{e} (tight-to-loose ratio)");
  hSignalEfficiency_MC_E_NJets->GetYaxis()->SetLabelSize(0.025);
  hSignalEfficiency_MC_E_NJets->SetLineColor(kRed);
  hSignalEfficiency_MC_E_NJets->GetYaxis()->SetRangeUser(0.,1.);
  hSignalEfficiency_MC_E_NJets->SetMarkerStyle(8);
  hSignalEfficiency_MC_E_NJets->SetMarkerColor(kRed);
  hSignalEfficiency_DATA_E_NJets->SetLineColor(kBlue);
  hSignalEfficiency_DATA_E_NJets->SetMarkerStyle(8);
  hSignalEfficiency_DATA_E_NJets->SetMarkerColor(kBlue);

  TCanvas *SENJetsCanvas_E = new TCanvas("SENJetsCanvas_E","SENJetsCanvas_E",200,10,600,800);
  SENJetsCanvas_E->cd();
  hSignalEfficiency_MC_E_NJets->Draw();
  hSignalEfficiency_DATA_E_NJets->Draw("same");
  qw->AddEntry(hSignalEfficiency_DATA_E_NJets,"DATA","pl");
  qw->AddEntry(hSignalEfficiency_MC_E_NJets,"MC","pl");
  qw->Draw("same");
  SENJetsCanvas_E->Update();
  title->Draw("same");
  SENJetsCanvas_E->SaveAs("SignalEfficiencyNJets_E.pdf");


 

  hSignalEfficiency_MC_E_Eta->GetXaxis()->SetTitle("#eta");
  hSignalEfficiency_MC_E_Eta->GetYaxis()->SetTitle("Efficiency");
  hSignalEfficiency_MC_E_Eta->SetTitle("#varepsilon_{s}^{e} (tight-to-loose ratio)");
  hSignalEfficiency_MC_E_Eta->GetYaxis()->SetLabelSize(0.025);
  hSignalEfficiency_MC_E_Eta->SetLineColor(kRed);
  hSignalEfficiency_MC_E_Eta->GetYaxis()->SetRangeUser(0.,1.);
  hSignalEfficiency_MC_E_Eta->SetMarkerStyle(8);
  hSignalEfficiency_MC_E_Eta->SetMarkerColor(kRed);
  hSignalEfficiency_DATA_E_Eta->SetLineColor(kBlue);
  hSignalEfficiency_DATA_E_Eta->SetMarkerStyle(8);
  hSignalEfficiency_DATA_E_Eta->SetMarkerColor(kBlue);

  TCanvas *SEEtaCanvas_E = new TCanvas("SEEtaCanvas_E","SEEtaCanvas_E",200,10,600,800);
  SEEtaCanvas_E->cd();
  hSignalEfficiency_MC_E_Eta->Draw();
  hSignalEfficiency_DATA_E_Eta->Draw("same");
  qw->Draw("same");
  SEEtaCanvas_E->Update();
  title->Draw("same");
  SEEtaCanvas_E->SaveAs("SignalEfficiencyEta_E.pdf");

  hSignalEfficiency_MC_E_PT->GetXaxis()->SetTitle("p_{T}");
  hSignalEfficiency_MC_E_PT->GetYaxis()->SetTitle("Efficiency");
  hSignalEfficiency_MC_E_PT->SetTitle("#varepsilon_{s}^{e} (tight-to-loose ratio)");
  hSignalEfficiency_MC_E_PT->GetYaxis()->SetLabelSize(0.025);
  hSignalEfficiency_MC_E_PT->SetLineColor(kRed);
  hSignalEfficiency_MC_E_PT->GetYaxis()->SetRangeUser(0.,1.);
  hSignalEfficiency_MC_E_PT->GetXaxis()->SetRangeUser(0.,150.);
  hSignalEfficiency_MC_E_PT->SetMarkerStyle(8);
  hSignalEfficiency_MC_E_PT->SetMarkerColor(kRed);
  hSignalEfficiency_DATA_E_PT->SetLineColor(kBlue);
  hSignalEfficiency_DATA_E_PT->SetMarkerStyle(8);
  hSignalEfficiency_DATA_E_PT->SetMarkerColor(kBlue);

  TCanvas *SEPTCanvas_E = new TCanvas("SEPTCanvas_E","SEPTCanvas_E",200,10,600,800);
  SEPTCanvas_E->cd();
  hSignalEfficiency_MC_E_PT->Draw();
  hSignalEfficiency_DATA_E_PT->Draw("same");
  qw->Draw("same");
  SEPTCanvas_E->Update();
  title->Draw("same");
  SEPTCanvas_E->SaveAs("SignalEfficiencyPT_E.pdf");



  hFakeRate_MC_E_NJets->GetXaxis()->SetTitle("Jet multiplicity");
  hFakeRate_MC_E_NJets->GetYaxis()->SetTitle("Fake Rate");
  hFakeRate_MC_E_NJets->SetTitle("#varepsilon_{f}^{e} (tight-to-loose ratio)");
  hFakeRate_MC_E_NJets->GetYaxis()->SetLabelSize(0.025);
  hFakeRate_MC_E_NJets->SetLineColor(kRed);
  hFakeRate_MC_E_NJets->GetYaxis()->SetRangeUser(0.,1.);
  hFakeRate_MC_E_NJets->SetMarkerStyle(8);
  hFakeRate_MC_E_NJets->SetMarkerColor(kRed);
  hFakeRate_DATA_E_NJets->SetLineColor(kBlue);
  hFakeRate_DATA_E_NJets->SetMarkerStyle(8);
  hFakeRate_DATA_E_NJets->SetMarkerColor(kBlue);

  TCanvas *FRNJetsCanvas_E = new TCanvas("FRNJetsCanvas_E","FRNJetsCanvas_E",200,10,600,800);
  FRNJetsCanvas_E->cd();
  hFakeRate_MC_E_NJets->Draw();
  hFakeRate_DATA_E_NJets->Draw("same");
  qw->Draw("same");
  FRNJetsCanvas_E->Update();
  title->Draw("same");
  FRNJetsCanvas_E->SaveAs("FakeRateNJets_E.pdf");

  hFakeRate_MC_E_Eta->GetXaxis()->SetTitle("#eta");
  hFakeRate_MC_E_Eta->GetYaxis()->SetTitle("Fake Rate");
  hFakeRate_MC_E_Eta->SetTitle("#varepsilon_{f}^{e} (tight-to-loose ratio)");
  hFakeRate_MC_E_Eta->GetYaxis()->SetLabelSize(0.025);
  hFakeRate_MC_E_Eta->SetLineColor(kRed);
  hFakeRate_MC_E_Eta->GetYaxis()->SetRangeUser(0.,1.);
  hFakeRate_MC_E_Eta->SetMarkerStyle(8);
  hFakeRate_MC_E_Eta->SetMarkerColor(kRed);
  hFakeRate_DATA_E_Eta->SetLineColor(kBlue);
  hFakeRate_DATA_E_Eta->SetMarkerStyle(8);
  hFakeRate_DATA_E_Eta->SetMarkerColor(kBlue);

  TCanvas *FREtaCanvas_E = new TCanvas("FREtaCanvas_E","FREtaCanvas_E",200,10,600,800);
  FREtaCanvas_E->cd();
  hFakeRate_MC_E_Eta->Draw();
  hFakeRate_DATA_E_Eta->Draw("same");
  qw->Draw("same");
  FREtaCanvas_E->Update();
  title->Draw("same");
  FREtaCanvas_E->SaveAs("FakeRateEta_E.pdf");

  hFakeRate_MC_E_PT->GetXaxis()->SetTitle("p_{T}");
  hFakeRate_MC_E_PT->GetYaxis()->SetTitle("Fake Rate");
  hFakeRate_MC_E_PT->SetTitle("#varepsilon_{f}^{e} (tight-to-loose ratio)");
  hFakeRate_MC_E_PT->GetYaxis()->SetLabelSize(0.025);
  hFakeRate_MC_E_PT->SetLineColor(kRed);
  hFakeRate_MC_E_PT->GetYaxis()->SetRangeUser(0.,1.);
  hFakeRate_MC_E_PT->GetXaxis()->SetRangeUser(0.,150.);
  hFakeRate_MC_E_PT->SetMarkerStyle(8);
  hFakeRate_MC_E_PT->SetMarkerColor(kRed);
  hFakeRate_DATA_E_PT->SetLineColor(kBlue);
  hFakeRate_DATA_E_PT->SetMarkerStyle(8);
  hFakeRate_DATA_E_PT->SetMarkerColor(kBlue);

  TCanvas *FRPTCanvas_E = new TCanvas("FRPTCanvas_E","FRPTCanvas_E",200,10,600,800);
  FRPTCanvas_E->cd();
  hFakeRate_MC_E_PT->Draw();
  hFakeRate_DATA_E_PT->Draw("same");
  qw->Draw("same");
  FRPTCanvas_E->Update();
  title->Draw("same");
  FRPTCanvas_E->SaveAs("FakeRatePT_E.pdf");











  hSignalEfficiency_MC_Mu_NJets->GetXaxis()->SetTitle("Jet multiplicity");
  hSignalEfficiency_MC_Mu_NJets->GetYaxis()->SetTitle("Efficiency");
  hSignalEfficiency_MC_Mu_NJets->SetTitle("#varepsilon_{s}^{#mu} (tight-to-loose ratio)");
  hSignalEfficiency_MC_Mu_NJets->GetYaxis()->SetLabelSize(0.025);
  hSignalEfficiency_MC_Mu_NJets->SetLineColor(kRed);
  hSignalEfficiency_MC_Mu_NJets->GetYaxis()->SetRangeUser(0.,1.);
  hSignalEfficiency_MC_Mu_NJets->SetMarkerStyle(8);
  hSignalEfficiency_MC_Mu_NJets->SetMarkerColor(kRed);
  hSignalEfficiency_DATA_Mu_NJets->SetLineColor(kBlue);
  hSignalEfficiency_DATA_Mu_NJets->SetMarkerStyle(8);
  hSignalEfficiency_DATA_Mu_NJets->SetMarkerColor(kBlue);

  TCanvas *SENJetsCanvas_Mu = new TCanvas("SENJetsCanvas_Mu","SENJetsCanvas_Mu",200,10,600,800);
  SENJetsCanvas_Mu->cd();
  hSignalEfficiency_MC_Mu_NJets->Draw();
  hSignalEfficiency_DATA_Mu_NJets->Draw("same");
  qw->Draw("same");
  SENJetsCanvas_Mu->Update();
  title->Draw("same");
  SENJetsCanvas_Mu->SaveAs("SignalEfficiencyNJets_Mu.pdf");

  hSignalEfficiency_MC_Mu_Eta->GetXaxis()->SetTitle("#eta");
  hSignalEfficiency_MC_Mu_Eta->GetYaxis()->SetTitle("Efficiency");
  hSignalEfficiency_MC_Mu_Eta->SetTitle("#varepsilon_{s}^{#mu} (tight-to-loose ratio)");
  hSignalEfficiency_MC_Mu_Eta->GetYaxis()->SetLabelSize(0.025);
  hSignalEfficiency_MC_Mu_Eta->SetLineColor(kRed);
  hSignalEfficiency_MC_Mu_Eta->GetYaxis()->SetRangeUser(0.,1.);
  hSignalEfficiency_MC_Mu_Eta->SetMarkerStyle(8);
  hSignalEfficiency_MC_Mu_Eta->SetMarkerColor(kRed);
  hSignalEfficiency_DATA_Mu_Eta->SetLineColor(kBlue);
  hSignalEfficiency_DATA_Mu_Eta->SetMarkerStyle(8);
  hSignalEfficiency_DATA_Mu_Eta->SetMarkerColor(kBlue);

  TCanvas *SEEtaCanvas_Mu = new TCanvas("SEEtaCanvas_Mu","SEEtaCanvas_Mu",200,10,600,800);
  SEEtaCanvas_Mu->cd();
  hSignalEfficiency_MC_Mu_Eta->Draw();
  hSignalEfficiency_DATA_Mu_Eta->Draw("same");
  qw->Draw("same");
  SEEtaCanvas_Mu->Update();
  title->Draw("same");
  SEEtaCanvas_Mu->SaveAs("SignalEfficiencyEta_Mu.pdf");

  hSignalEfficiency_MC_Mu_PT->GetXaxis()->SetTitle("p_{T}");
  hSignalEfficiency_MC_Mu_PT->GetYaxis()->SetTitle("Efficiency");
  hSignalEfficiency_MC_Mu_PT->SetTitle("#varepsilon_{s}^{#mu} (tight-to-loose ratio)");
  hSignalEfficiency_MC_Mu_PT->GetYaxis()->SetLabelSize(0.025);
  hSignalEfficiency_MC_Mu_PT->SetLineColor(kRed);
  hSignalEfficiency_MC_Mu_PT->GetYaxis()->SetRangeUser(0.,1.);
  hSignalEfficiency_MC_Mu_PT->GetXaxis()->SetRangeUser(0.,150.);
  hSignalEfficiency_MC_Mu_PT->SetMarkerStyle(8);
  hSignalEfficiency_MC_Mu_PT->SetMarkerColor(kRed);
  hSignalEfficiency_DATA_Mu_PT->SetLineColor(kBlue);
  hSignalEfficiency_DATA_Mu_PT->SetMarkerStyle(8);
  hSignalEfficiency_DATA_Mu_PT->SetMarkerColor(kBlue);

  TCanvas *SEPTCanvas_Mu = new TCanvas("SEPTCanvas_Mu","SEPTCanvas_Mu",200,10,600,800);
  SEPTCanvas_Mu->cd();
  hSignalEfficiency_MC_Mu_PT->Draw();
  hSignalEfficiency_DATA_Mu_PT->Draw("same");
  qw->Draw("same");
  SEPTCanvas_Mu->Update();
  title->Draw("same");
  SEPTCanvas_Mu->SaveAs("SignalEfficiencyPT_Mu.pdf");



  hFakeRate_MC_Mu_NJets->GetXaxis()->SetTitle("Jet multiplicity");
  hFakeRate_MC_Mu_NJets->GetYaxis()->SetTitle("Fake Rate");
  hFakeRate_MC_Mu_NJets->SetTitle("#varepsilon_{f}^{#mu} (tight-to-loose ratio)");
  hFakeRate_MC_Mu_NJets->GetYaxis()->SetLabelSize(0.025);
  hFakeRate_MC_Mu_NJets->SetLineColor(kRed);
  hFakeRate_MC_Mu_NJets->GetYaxis()->SetRangeUser(0.,1.);
  hFakeRate_MC_Mu_NJets->SetMarkerStyle(8);
  hFakeRate_MC_Mu_NJets->SetMarkerColor(kRed);
  hFakeRate_DATA_Mu_NJets->SetLineColor(kBlue);
  hFakeRate_DATA_Mu_NJets->SetMarkerStyle(8);
  hFakeRate_DATA_Mu_NJets->SetMarkerColor(kBlue);

  TCanvas *FRNJetsCanvas_Mu = new TCanvas("FRNJetsCanvas_Mu","FRNJetsCanvas_Mu",200,10,600,800);
  FRNJetsCanvas_Mu->cd();
  hFakeRate_MC_Mu_NJets->Draw();
  hFakeRate_DATA_Mu_NJets->Draw("same");
  qw->Draw("same");
  FRNJetsCanvas_Mu->Update();
  title->Draw("same");
  FRNJetsCanvas_Mu->SaveAs("FakeRateNJets_Mu.pdf");

  hFakeRate_MC_Mu_Eta->GetXaxis()->SetTitle("#eta");
  hFakeRate_MC_Mu_Eta->GetYaxis()->SetTitle("Fake Rate");
  hFakeRate_MC_Mu_Eta->SetTitle("#varepsilon_{f}^{#mu} (tight-to-loose ratio)");
  hFakeRate_MC_Mu_Eta->GetYaxis()->SetLabelSize(0.025);
  hFakeRate_MC_Mu_Eta->SetLineColor(kRed);
  hFakeRate_MC_Mu_Eta->GetYaxis()->SetRangeUser(0.,1.);
  hFakeRate_MC_Mu_Eta->SetMarkerStyle(8);
  hFakeRate_MC_Mu_Eta->SetMarkerColor(kRed);
  hFakeRate_DATA_Mu_Eta->SetLineColor(kBlue);
  hFakeRate_DATA_Mu_Eta->SetMarkerStyle(8);
  hFakeRate_DATA_Mu_Eta->SetMarkerColor(kBlue);

  TCanvas *FREtaCanvas_Mu = new TCanvas("FREtaCanvas_Mu","FREtaCanvas_Mu",200,10,600,800);
  FREtaCanvas_Mu->cd();
  hFakeRate_MC_Mu_Eta->Draw();
  hFakeRate_DATA_Mu_Eta->Draw("same");
  qw->Draw("same");
  FREtaCanvas_Mu->Update();
  title->Draw("same");
  FREtaCanvas_Mu->SaveAs("FakeRateEta_Mu.pdf");

  hFakeRate_MC_Mu_PT->GetXaxis()->SetTitle("p_{T}");
  hFakeRate_MC_Mu_PT->GetYaxis()->SetTitle("Fake Rate");
  hFakeRate_MC_Mu_PT->SetTitle("#varepsilon_{f}^{#mu} (tight-to-loose ratio)");
  hFakeRate_MC_Mu_PT->GetYaxis()->SetLabelSize(0.025);
  hFakeRate_MC_Mu_PT->SetLineColor(kRed);
  hFakeRate_MC_Mu_PT->GetYaxis()->SetRangeUser(0.,1.);
  hFakeRate_MC_Mu_PT->GetXaxis()->SetRangeUser(0.,150.);
  hFakeRate_MC_Mu_PT->SetMarkerStyle(8);
  hFakeRate_MC_Mu_PT->SetMarkerColor(kRed);
  hFakeRate_DATA_Mu_PT->SetLineColor(kBlue);
  hFakeRate_DATA_Mu_PT->SetMarkerStyle(8);
  hFakeRate_DATA_Mu_PT->SetMarkerColor(kBlue);

  TCanvas *FRPTCanvas_Mu = new TCanvas("FRPTCanvas_Mu","FRPTCanvas_Mu",200,10,600,800);
  FRPTCanvas_Mu->cd();
  hFakeRate_MC_Mu_PT->Draw();
  hFakeRate_DATA_Mu_PT->Draw("same");
  qw->Draw("same");
  FRPTCanvas_Mu->Update();
  title->Draw("same");
  FRPTCanvas_Mu->SaveAs("FakeRatePT_Mu.pdf");


  hFakeRate_DATA_Mu_NVTX->GetXaxis()->SetTitle("Number of vertices");
  hFakeRate_DATA_Mu_NVTX->GetYaxis()->SetTitle("Fake Rate");
  hFakeRate_DATA_Mu_NVTX->SetTitle("#varepsilon_{f}^{#mu} (tight-to-loose ratio)");
  hFakeRate_DATA_Mu_NVTX->GetYaxis()->SetLabelSize(0.025);
  hFakeRate_DATA_Mu_NVTX->SetLineColor(kRed);
  hFakeRate_DATA_Mu_NVTX->GetYaxis()->SetRangeUser(0.,1.);
  hFakeRate_DATA_Mu_NVTX->GetXaxis()->SetRangeUser(0,15);
  hFakeRate_DATA_Mu_NVTX->SetMarkerStyle(8);
  hFakeRate_DATA_Mu_NVTX->SetMarkerColor(kBlack);

  TCanvas *FRNVTXCanvas_Mu = new TCanvas("FRNVTXCanvas_Mu","FRNVTXCanvas_Mu",200,10,600,800);
  FRNVTXCanvas_Mu->cd();
  hFakeRate_DATA_Mu_NVTX->Draw();
  //hFakeRate_DATA_Mu_PT->Draw("same");
  //  qw->Draw("same");
  FRNVTXCanvas_Mu->Update();
  title->Draw("same");
  FRNVTXCanvas_Mu->SaveAs("FakeRateNVTX_Mu.pdf");



  hFakeRate_DATA_E_NVTX->GetXaxis()->SetTitle("Number of vertices");
  hFakeRate_DATA_E_NVTX->GetYaxis()->SetTitle("Fake Rate");
  hFakeRate_DATA_E_NVTX->SetTitle("#varepsilon_{f}^{e} (tight-to-loose ratio)");
  hFakeRate_DATA_E_NVTX->GetYaxis()->SetLabelSize(0.025);
  hFakeRate_DATA_E_NVTX->SetLineColor(kRed);
  hFakeRate_DATA_E_NVTX->GetYaxis()->SetRangeUser(0.,1.);
  hFakeRate_DATA_E_NVTX->GetXaxis()->SetRangeUser(0,15);
  hFakeRate_DATA_E_NVTX->SetMarkerStyle(8);
  hFakeRate_DATA_E_NVTX->SetMarkerColor(kBlack);

  TCanvas *FRNVTXCanvas_E = new TCanvas("FRNVTXCanvas_E","FRNVTXCanvas_E",200,10,600,800);
  FRNVTXCanvas_E->cd();
  hFakeRate_DATA_E_NVTX->Draw();
  //hFakeRate_DATA_E_PT->Draw("same");
  //  qw->Draw("same");
  FRNVTXCanvas_E->Update();
  title->Draw("same");
  FRNVTXCanvas_E->SaveAs("FakeRateNVTX_E.pdf");




  TLegend* qw2 = 0;
  qw2 = new TLegend(0.65,0.99,0.93,0.81);
  qw2->SetFillColor(kWhite);

  hFakeRate_DATA_E_NJets_20->GetXaxis()->SetTitle("Jet multiplicity");
  hFakeRate_DATA_E_NJets_20->GetYaxis()->SetTitle("Fake Rate");
  hFakeRate_DATA_E_NJets_20->SetTitle("#varepsilon_{f}^{e} (tight-to-loose ratio)");
  hFakeRate_DATA_E_NJets_20->GetYaxis()->SetLabelSize(0.025);
  hFakeRate_DATA_E_NJets_20->GetYaxis()->SetRangeUser(0.,1.);

  hFakeRate_DATA_E_NJets_20->SetMarkerStyle(8);
  hFakeRate_DATA_E_NJets_20->SetMarkerColor(kRed);
  hFakeRate_DATA_E_NJets_20->SetLineColor(kRed);

  hFakeRate_DATA_E_NJets_60->SetMarkerStyle(8);
  hFakeRate_DATA_E_NJets_60->SetMarkerColor(kGreen);
  hFakeRate_DATA_E_NJets_60->SetLineColor(kGreen);

  hFakeRate_DATA_E_NJets_90->SetMarkerStyle(8);
  hFakeRate_DATA_E_NJets_90->SetMarkerColor(kBlue);
  hFakeRate_DATA_E_NJets_90->SetLineColor(kBlue);

  TCanvas *FRNJetsCanvas_SuperImposed_E = new TCanvas("FRNJetsCanvas_SuperImposed_E","FRNJetsCanvas_SuperImposed_E",200,10,600,800);
  FRNJetsCanvas_SuperImposed_E->cd();
  hFakeRate_DATA_E_NJets_20->Draw();
  hFakeRate_DATA_E_NJets_60->Draw("same");
  hFakeRate_DATA_E_NJets_90->Draw("same");
  qw2->AddEntry(hFakeRate_DATA_E_NJets_20,"Jet 20","pl");
  qw2->AddEntry(hFakeRate_DATA_E_NJets_60,"Jet 60","pl");
  qw2->AddEntry(hFakeRate_DATA_E_NJets_90,"Jet 90","pl");
  qw2->Draw("same");
  FRNJetsCanvas_SuperImposed_E->Update();
  title->Draw("same");
  FRNJetsCanvas_SuperImposed_E->SaveAs("FakeRateNJets_SuperImposed_E.pdf");



  hFakeRate_DATA_Mu_NJets_20->GetXaxis()->SetTitle("Jet multiplicity");
  hFakeRate_DATA_Mu_NJets_20->GetYaxis()->SetTitle("Fake Rate");
  hFakeRate_DATA_Mu_NJets_20->SetTitle("#varepsilon_{f}^{#mu} (tight-to-loose ratio)");
  hFakeRate_DATA_Mu_NJets_20->GetYaxis()->SetLabelSize(0.025);
  hFakeRate_DATA_Mu_NJets_20->GetYaxis()->SetRangeUser(0.,1.);

  hFakeRate_DATA_Mu_NJets_20->SetMarkerStyle(8);
  hFakeRate_DATA_Mu_NJets_20->SetMarkerColor(kRed);
  hFakeRate_DATA_Mu_NJets_20->SetLineColor(kRed);

  hFakeRate_DATA_Mu_NJets_60->SetMarkerStyle(8);
  hFakeRate_DATA_Mu_NJets_60->SetMarkerColor(kGreen);
  hFakeRate_DATA_Mu_NJets_60->SetLineColor(kGreen);

  hFakeRate_DATA_Mu_NJets_90->SetMarkerStyle(8);
  hFakeRate_DATA_Mu_NJets_90->SetMarkerColor(kBlue);
  hFakeRate_DATA_Mu_NJets_90->SetLineColor(kBlue);

  TCanvas *FRNJetsCanvas_SuperImposed_Mu = new TCanvas("FRNJetsCanvas_SuperImposed_Mu","FRNJetsCanvas_SuperImposed_Mu",200,10,600,800);
  FRNJetsCanvas_SuperImposed_Mu->cd();
  hFakeRate_DATA_Mu_NJets_20->Draw();
  hFakeRate_DATA_Mu_NJets_60->Draw("same");
  hFakeRate_DATA_Mu_NJets_90->Draw("same");
  qw2->Draw("same");
  FRNJetsCanvas_SuperImposed_Mu->Update();
  title->Draw("same");
  FRNJetsCanvas_SuperImposed_Mu->SaveAs("FakeRateNJets_SuperImposed_Mu.pdf");













}
