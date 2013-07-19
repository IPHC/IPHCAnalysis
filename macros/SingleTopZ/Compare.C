
void CompareHisto(char* histo, char* file1, char* file2, char* leg1=" ", char* leg2=" ", bool log=false, char* filename="", bool print=false)
{
 TFile f1(file1,"read");
 TFile f2(file2,"read");
 
 TH1F* h1 = (TH1F*) f1.Get(histo);
 TH1F* h2 = (TH1F*) f2.Get(histo);

 if(!h1 || !h2) return;

 Compare(h1, h2, leg1, leg2, log, filename, print);
}

void Draw2Histos(char* histo1, char* histo2, char* file1, char* leg1=" ", char* leg2=" ", bool log=false, char* filename="", bool print=false, bool ratio=false)
{
 TFile f1(file1,"read");
 TH1F* h1 = (TH1F*) f1.Get(histo1);
 TH1F* h2 = (TH1F*) f1.Get(histo2);
 if(!h1 || !h2) { cout<<"No such histo"<<endl; return;}
 Compare(h1, h2, leg1, leg2, log, filename, print, ratio);
}

void Draw3Histos(char* histo1, char* histo2, char* histo3, 
  char* file1, char* leg1=" ", char* leg2=" ", char* leg3=" ", bool log=false, char* filename="", bool print=false)
{
 TFile f1(file1,"read");
 TH1F* h1 = (TH1F*) f1.Get(histo1);
 TH1F* h2 = (TH1F*) f1.Get(histo2);
 TH1F* h3 = (TH1F*) f1.Get(histo3);
 if(!h1 || !h2|| !h3) { cout<<"No such histo"<<endl; return;}
 Compare(h1, h2, h3, leg1, leg2, leg3, log, filename, print);
}

void Compare(TH1F* h1, TH1F* h2, char* leg1=" ", char* leg2=" ", bool log=false, char* filename="", bool print=false, bool ratio=false)
{

// h1->Sumw2();
// h2->Sumw2();
 
 cout<<"h1 : N "<<h1->Integral()<<" mean "<<h1->GetMean()<<"+/-"<<h1->GetMeanError()<<" rms "<<h1->GetRMS()<<"+/-"<<h1->GetRMSError()<<endl;
 cout<<"h2 : N "<<h2->Integral()<<" mean "<<h2->GetMean()<<"+/-"<<h2->GetMeanError()<<" rms "<<h2->GetRMS()<<"+/-"<<h2->GetRMSError()<<endl;

 h1->Rebin(4);
 h2->Rebin(4);

 h1->Scale(1/h1->Integral());
 h2->Scale(1/h2->Integral());

 TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 700, 500);
 if(h1->GetMaximum()>h2->GetMaximum()) h1->Draw("h");
 else h2->Draw("h");
 h2->SetLineColor(4);
 h2->SetLineWidth(2);
 h2->Draw("h:same");
 h1->SetLineColor(2);
 h1->SetLineWidth(2);
 h1->Draw("h:same");
 
 TLegend leg(0.6, 0.65, 0.85, 0.8);
 leg.AddEntry(h1, leg1, "l");
 leg.AddEntry(h2, leg2, "l");
 if(strcmp(leg1,"") || strcmp(leg2,"")) leg.Draw();

 //h2->Fit("gaus","+","", -0.05, 0.04);
 //h3->Fit("gaus","+","", -0.05, 0.04);
 h1->SetStats(0);
 h2->SetStats(0);

 if(log) c1->SetLogy();

 TCanvas *c2;
 if(ratio)
 {
   c2 = new TCanvas("c2", "c2", 700, 0, 700, 500);
   TH1F* hr = h2->Clone();
   hr->Divide(h1);
   hr->Rebin(2);
   hr->Scale(0.5);
   hr->SetMinimum(0);
   hr->SetMaximum(5);
   hr->Draw("histo:e");
   float xmin = h1->GetMean()-h1->GetRMS();
   float xmax = h1->GetMean()+h1->GetRMS();
   //hr->Fit("pol1", "", "", xmin, xmax);
   //cout<<" fit range : "<<xmin<<" "<<xmax<<endl;  
   TLine l(hr->GetXaxis()->GetXmin(), 1, hr->GetXaxis()->GetXmax(), 1);
   l.Draw();
   c2->Modified();
   c2->Update();
 }

 c1->Modified();
 c1->Update();
 getchar();
 if(strcmp(filename,"") && print) c1->Print(filename);
 c2->Delete();

}



void Compare(TH1F* h1, TH1F* h2, TH1F* h3, char* leg1=" ", char* leg2=" ", char* leg3=" ", bool log=false, char* filename="", bool print=false)
{

 cout<<"h1 : N "<<h1->Integral()<<" mean "<<h1->GetMean()<<"+/-"<<h1->GetMeanError()<<" rms "<<h1->GetRMS()<<"+/-"<<h1->GetRMSError()<<endl;
 cout<<"h2 : N "<<h2->Integral()<<" mean "<<h2->GetMean()<<"+/-"<<h2->GetMeanError()<<" rms "<<h2->GetRMS()<<"+/-"<<h2->GetRMSError()<<endl;
 cout<<"h3 : N "<<h3->Integral()<<" mean "<<h3->GetMean()<<"+/-"<<h3->GetMeanError()<<" rms "<<h3->GetRMS()<<"+/-"<<h3->GetRMSError()<<endl;

 if(h1->GetMaximum()>h2->GetMaximum() && h1->GetMaximum()>h3->GetMaximum()) h1->Draw("h");
 else
 if(h2->GetMaximum()>h1->GetMaximum() && h2->GetMaximum()>h3->GetMaximum()) h2->Draw("h");
 else h3->Draw("h");
 h2->SetLineColor(4);
 h2->SetLineWidth(2);
 h2->Draw("h:same");
 h1->SetLineColor(2);
 h1->SetLineWidth(2);
 h1->Draw("h:same");
 h3->SetLineColor(8);
 h3->SetLineWidth(2);
 h3->Draw("h:same");
 
 TLegend leg(0.6, 0.65, 0.85, 0.8);
 leg.AddEntry(h1, leg1, "l");
 leg.AddEntry(h2, leg2, "l");
 leg.AddEntry(h3, leg3, "l");
 if(strcmp(leg1,"") || strcmp(leg2,"") || strcmp(leg3,"")) leg.Draw();

 h1->SetStats(0);
 h2->SetStats(0);
 h3->SetStats(0);
 h1->SetTitle(0);
 h2->SetTitle(0);
 h3->SetTitle(0);

 if(log) c1->SetLogy();
 c1->Modified();
 c1->Update();
 getchar();
 if(strcmp(filename,"") && print) c1->Print(filename);

}

