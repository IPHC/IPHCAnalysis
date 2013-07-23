
void DrawSelection(TCanvas* canv, TH1F* estimate, THStack* expectation, TH1F** mc, unsigned int* mc_nr,
   const int Nd, TLegend* qw, string* dataset, int NjetsMin, int NjetsMax)
{
  
  estimate->SetMarkerStyle(20);
  estimate->GetXaxis()->SetTitle("Jet multiplicity");
  estimate->GetYaxis()->SetTitle("Number of events");
  estimate->GetYaxis()->SetLabelSize(0.025);
  estimate->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  estimate->Draw();
  for(unsigned int id=0; id<Nd; id++)
  if(mc[mc_nr[id]]) expectation->Add(mc[mc_nr[id]]);
  expectation->Draw("hesame");
  estimate->Draw("same");
  TH1F * tmp  = new TH1F("tmp",  "tmp",  11, -0.5,  10.5);
  tmp->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  for(unsigned int id=0; id<Nd; id++)
  if(mc[mc_nr[id]])
  {
    tmp->Add(mc[mc_nr[id]]);
    qw->AddEntry(mc[mc_nr[id]], dataset[mc_nr[id]].c_str(),"f");
  }

  TGraphErrors *thegraph = new TGraphErrors(tmp);
  thegraph->SetFillStyle(3005);
  thegraph->SetFillColor(1);
  thegraph->Draw("e2same");
  qw->AddEntry(expectation,    " Matrix Method  ","p");
  qw->Draw("same");

}
    

void MatrixMethod_Results_Estimations(string selection, string eventType, string ntupleType, string fast, string other,
 string channel="MuMu"){


  float NjetsMin = -0.5;
  float NjetsMax = 10.5;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile * file_Expected   = new TFile(("MatrixMethod_OutPut_"+channel+"_MC"+other+".root").c_str());
  TFile * file_Estimated  = new TFile(("MatrixMethod_OutPut_"+channel+ntupleType+fast+other+".root").c_str());
  
  string Sel[] = {"Tight", "Medium", "Loose"};
  string Type[] = {"Signal", "Wlike", "QCDlike"};

  const int Nd = 24;
  string dataset[Nd] = {"TTbarSignal", "TTbarSemileptonic", "TTbar", "Zjets", 
                       "DYToMuMu_M-20", "DYToMuMu_M-10To20", "DYToTauTau_M-20","DYToTauTau_M-10To20",
		       "Wjets", "Tt", "Ts", "TtW", "Tbart", "Tbars", "TbartW", "WW", "WZ", "ZZ", 
		       "BBbar_Pt25", "BBbar_Pt15", 
		       "QCD_Pt-50to80_Mu", "QCD_Pt-80to120_Mu", "QCD_Pt-120to150_Mu", "QCD_Pt-150_Mu"};
  int color[Nd] = {kRed+1, kRed-7, kRed-7, kAzure-2, kAzure-2, kAzure-2, kAzure-2, kAzure-2, 
                   kGreen-3, kMagenta-6, kMagenta-6, kMagenta-6, kMagenta-6, kMagenta-6, kMagenta-6, kGray+2, kGray+2, kGray+2, 
		   kYellow+3, kYellow+3, kYellow+3, kYellow+3, kYellow+3, kYellow+3};


  // Get histos of expected # of events
  file_Expected->cd();

  TH1F* mc[3][Nd]; // Tight, Medium, Loose
  for(int isel=0; isel<3; isel++)
  for(unsigned int id=0; id<Nd; id++)
  { 
    string histoname = "MMExpected_"+Sel[isel]+channel+"_"+dataset[id];
    mc[isel][id] = (TH1F*)gROOT->FindObject(histoname.c_str());
    
    if(mc[isel][id])
    {
      mc[isel][id]->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
      mc[isel][id]->SetFillColor(color[id]);
    }
  }

  
  // Get histos of estimated # of events
  file_Estimated->cd();
  TH1F* estimate[3][3]; // Tight, Medium, Loose // Signal, Wlike, QCDlike
  
  for(int isel=0; isel<3; isel++)
  {
    estimate[isel][0]  = (TH1F*)gROOT->FindObject(("MMEstimated_"+Sel[isel]+channel+"_Signal").c_str() );
    estimate[isel][1]  = (TH1F*)gROOT->FindObject(("MMEstimated_"+Sel[isel]+channel+"_WJets").c_str() );
    estimate[isel][2]  = (TH1F*)gROOT->FindObject(("MMEstimated_"+Sel[isel]+channel+"_QCD").c_str() );
  }

  
  const int NdS = 3;
  unsigned int signallike[] = {0, 16, 17}; 
/*  tightSignalLike->Add(tight_TTDilept);
  tightSignalLike->Add(tight_TtW);
  tightSignalLike->Add(tight_TbartW);
  tightSignalLike->Add(tight_ZJets);
  tightSignalLike->Add(tight_DYHigh);
  tightSignalLike->Add(tight_DYLow);
  tightSignalLike->Add(tight_DYTauTauHigh);
  tightSignalLike->Add(tight_DYTauTauLow);
  tightSignalLike->Add(tight_WW);
  tightSignalLike->Add(tight_WZ);
  tightSignalLike->Add(tight_ZZ);*/

  const int NdW = 14;
  unsigned int wlike[] = {1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  //tightWLike->Add(tight_TTSemi);
  //tightWLike->Add(tight_WJets);
  
  const int NdQCD = 7;
  unsigned int qcdlike[] = {2, 18, 19, 20, 21, 22, 23};


  THStack* hstack[3][3];
  unsigned int* datasetlist[3] = {signallike, wlike, qcdlike};
  const int Ndlist[] = {NdS, NdW, NdQCD};
  
  TLegend* qw = 0;
  qw = new TLegend(0.70,0.75,0.95,0.95);
  qw->SetFillColor(kWhite);

  title = new TLatex(10.,10.,"CMS preliminary");
  title->SetNDC();
  title->SetTextAlign(12);
  title->SetX(0.60);
  title->SetY(0.70);
  title->SetTextFont(42);
  title->SetTextSize(0.04);
  title->SetTextSizePixels(24);

  // Draw graphs for each selection and type
  for(int isel=0; isel<3; isel++)
  if(selection == Sel[isel])
  {
    for(int ityp=0; ityp<3; ityp++)
    if(eventType == Type[ityp])
    { 
      string hstackname = Sel[isel]+Type[ityp];
      string hstacktitle = Sel[isel]+" "+Type[ityp];
      hstack[isel][ityp] = new THStack(hstackname.c_str(), hstacktitle.c_str());
      TCanvas * canv = new TCanvas(("canv"+channel+selection+eventType+ntupleType+fast+other).c_str(),
       ("canv"+channel+selection+eventType+ntupleType+fast+other).c_str(), 300, 400);
      string estimatetitle;
      if(channel=="MuMu") estimatetitle = Type[ityp]+" #mu#mu ("+Sel[isel]+")";
      if(channel=="EE") estimatetitle = Type[ityp]+" ee ("+Sel[isel]+")";
      estimate[isel][ityp]->SetTitle(estimatetitle.c_str());
      
      DrawSelection(canv, estimate[isel][ityp], hstack[isel][ityp],
        mc[isel], datasetlist[ityp], Ndlist[ityp],
        qw, dataset, NjetsMin, NjetsMax);
      title->Draw("same");
      canv->SaveAs((Sel[isel]+Type[ityp]+channel+ntupleType+fast+other+".pdf").c_str());
    
    }
  }
 
}
