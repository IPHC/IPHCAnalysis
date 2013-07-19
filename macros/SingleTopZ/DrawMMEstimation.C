#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TKey.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include <iostream>

using namespace std;

void SetHistograms(TH1F *&hS, TH1F *&hZ, TH1F *&hother, TH1F *&htot, TH1F* h, TString name)
{
  //cout<<hS<<endl;
  // Z + 1 lepton
  if(name.Contains("WZ") || name.Contains("ZZ"))
  {
    if(!hS) 
    hS = (TH1F*) h->Clone();
    else hS->Add(h);
  }
  // Z + 1 fake
  else if(name.Contains("DYToLL") || name.Contains("Zjets"))
      //|| name.Contains("TTbar"))
     // || name.Contains("TtW") || name.Contains("TbartW") || name.Contains("WW"))
  {
    if(!hZ) hZ = (TH1F*) h->Clone();
    else hZ->Add(h);
  }
  // other
  else
  {
    //cout<<name<<endl;
    if(!hother) hother = (TH1F*) h->Clone();
    else hother->Add(h);
  }
  // all
  if(!htot) htot = (TH1F*) h->Clone();
  else htot->Add(h);
  //cout<<hS<<"-"<<endl;

}

//------------------------------------------------------------------------------------------------
// root -l DrawMMEstimation.C+ en raison de crash avec l'interpreteur
void DrawMMEstimation(bool isData=false)
{
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
 
  TFile* f = new TFile("MMethod.root", "read");

  string dataset="";
  string channels[] = {"mumumu", "mumue", "eemu", "eee"};

  TCanvas c1;
  c1.Divide(2,2);

  // Loop over channels
  for(int iChannel=0; iChannel<4; iChannel++)
  {
    cout<<"-------------"<<endl; 
    cout<<"Channel "<<channels[iChannel]<<endl; 
    cout<<"-------------"<<endl; 
  
    if(isData)
    {
      if(iChannel==0) dataset="DataMu";
      if(iChannel==1) dataset="DataMuEG";
      if(iChannel==2) dataset="DataMuEG";
      if(iChannel==3) dataset="DataEG"; 
    }
    else dataset="MC";
       
    // Get estimations histos
    TH1F* hS_est = (TH1F*) f->Get(string("Estimation_"+channels[iChannel]+"_S_"+dataset).c_str());
    TH1F* hZ_est = (TH1F*) f->Get(string("Estimation_"+channels[iChannel]+"_Z_"+dataset).c_str());
    
    if(!hS_est || !hZ_est) {cout<<"No estimation histo found"<<endl; return;}
    TH1F* htot_est = (TH1F*) hS_est->Clone();
    htot_est->Add(hZ_est);
    
    // Get MC truth histos
    TH1F* hS=0;
    TH1F* hZ=0;
    TH1F* hother=0;
    TH1F* htot=0;
    
    TH1F* hS_loose=0;
    TH1F* hZ_loose=0;
    TH1F* hother_loose=0;
    TH1F* htot_loose=0;

    TH1F* hS_err=0;
    TH1F* hZ_err=0;
    TH1F* hother_err=0;
    TH1F* htot_err=0;
    
    TH1F* hS_loose_err=0;
    TH1F* hZ_loose_err=0;
    TH1F* hother_loose_err=0;
    TH1F* htot_loose_err=0;

    // Loop over histos in file to get MC histos
    TIter nextkey( gDirectory->GetListOfKeys() );
    TKey *key;
    while ( (key = (TKey*)nextkey()))
    {
      TObject *obj = key->ReadObj();
      if(! obj->IsA()->InheritsFrom( TH1F::Class() )) continue;

      TString name(obj->GetName());
      if(!name.Contains("CutFlow")) continue;
      if(name.Contains("Data") || name.Contains("FCNC")) continue;

      TH1F* h = (TH1F*) obj;

      if(!name.Contains(channels[iChannel])) continue;
      //cout<<name<<endl;

      if(name.Contains("CutFlow") && !name.Contains("ErrCutFlow"))
      {
	if(name.Contains("tight")) SetHistograms(hS, hZ, hother, htot, h, name);
	if(name.Contains("loose")) SetHistograms(hS_loose, hZ_loose, hother_loose, htot_loose, h, name);
      }
      if(name.Contains("ErrCutFlow"))
      {
	if(name.Contains("tight")) SetHistograms(hS_err, hZ_err, hother_err, htot_err, h, name);
	if(name.Contains("loose")) SetHistograms(hS_loose_err, hZ_loose_err, hother_loose_err, htot_loose_err, h, name);
      }
    }

    
    if(!hS || !hZ) {cout<<"No MC truth histo found"<<endl; continue;}
    if(!hS_loose || !hZ_loose) {cout<<"No MC truth loose histo found"<<endl; continue;}
   
   
   // Draw MC histos and estimations
   //-------------------------------
   
   // Draw histos for 3 leptons datasets
   TVirtualPad* p1=c1.cd(1);
   p1->SetLogy();
      
   hS->SetFillColor(8);
   hS->GetXaxis()->SetRangeUser(0.5,5.5);
   hS->SetTitle(Form("Cut flow, %s channel", channels[iChannel].c_str()));
   hS->GetXaxis()->SetTitle("Cut");
   hS->Draw("histo");
   TH1F* hSd_err = (TH1F*) hS->Clone();
   // Set errors
   if(hS_err)
   for(int ibin=0; ibin<hSd_err->GetNbinsX(); ibin++)
   {
     hSd_err->SetBinError(ibin, sqrt(hS_err->GetBinContent(ibin))); // square of errors in ErrCutFlow 
     hS->SetBinError(ibin, sqrt(hS_err->GetBinContent(ibin)));
   }
   hSd_err->SetFillStyle(3004);
   hSd_err->SetFillColor(12);
   hSd_err->Draw("e2:same");
   hS_est->SetMarkerStyle(20);
   hS_est->Draw("pe:same");
   
   TVirtualPad* p2=c1.cd(2);
   p2->SetLogy(0);
   
   TH1F* hS_est_ratio = (TH1F*) hS_est->Clone();
   hS_est_ratio->SetTitle(Form("Ratio, %s channel", channels[iChannel].c_str()));
   hS_est_ratio->GetXaxis()->SetTitle("Cut");
   hS_est_ratio->GetYaxis()->SetTitle("Ratio");
   hS_est_ratio->Divide(hS);
   hS_est_ratio->GetXaxis()->SetRangeUser(0.5,5.5);
   hS_est_ratio->Draw("pe");
   TLine l(0.5, 1, 6.5, 1);
   l.Draw();
   
   
   // Draw histos for 2 leptons + 1 fake datasets
   TVirtualPad* p3=c1.cd(3);
   p3->SetLogy();
   
   TH1F* hfaketot = (TH1F*) hother->Clone();
   hfaketot->Add(hZ);
   hfaketot->SetTitle(Form("Cut flow, %s channel", channels[iChannel].c_str()));
   hfaketot->GetXaxis()->SetTitle("Cut");
   hfaketot->SetFillColor(45);//16 //kOrange+8
   hfaketot->GetXaxis()->SetRangeUser(0.5,5.5);
   hfaketot->Draw("histo");

   hZ->SetFillColor(46);
   hZ->Draw("same");
      
   TH1F* hfaketot_err = (TH1F*) hfaketot->Clone();
   // Set errors
   if(hZ_err && hother_err)
   for(int ibin=0; ibin<hfaketot_err->GetNbinsX(); ibin++)
   {
     float err=sqrt(hZ_err->GetBinContent(ibin)+hother_err->GetBinContent(ibin)); 
     // sous-estime erreur quand 0 evt dans 1 dataset MC
     hfaketot_err->SetBinError(ibin, err);
     hfaketot->SetBinError(ibin, err);
   }
   hfaketot_err->SetFillStyle(3004);
   hfaketot_err->SetFillColor(12);
   hfaketot_err->Draw("e2:same");
   
   hZ_est->SetMarkerStyle(20);
   hZ_est->Draw("pe:same");
   hZ_est->Print("all");
   
   TVirtualPad* p4=c1.cd(4);
   p4->SetLogy(0);
   
   TH1F* hZ_est_ratio = (TH1F*) hZ_est->Clone();
   hZ_est_ratio->SetTitle(Form("Ratio, %s channel", channels[iChannel].c_str()));
   hZ_est_ratio->GetXaxis()->SetTitle("Cut");
   hZ_est_ratio->GetYaxis()->SetTitle("Ratio");
   hZ_est_ratio->Divide(hfaketot);
   hZ_est_ratio->GetXaxis()->SetRangeUser(0.5,5.5);
   hZ_est_ratio->Draw("pe");
   l.Draw();

   c1.cd();
   c1.Modified();
   c1.Update();
   c1.Print(Form("MMEstimation_%s_%s.pdf", channels[iChannel].c_str(), dataset.c_str()));
   
   // for MC deduce MC truth Eff and Fake rate
   TCanvas c2;
   TCanvas c3;
   if(!isData)
   {
     TGraphErrors *geff = new TGraphErrors();
     TGraphErrors *gfake = new TGraphErrors();
   
     float Eff, Eff_err, Fake_Z, Fake_Z_err, Fake_other, Fake_other_err, Fake_all, Fake_all_err;
     for(int iCut=2; iCut<7; iCut++)
     {
       Eff=0;
       Eff_err=0;
       Fake_Z=0;
       Fake_Z_err=0;
       Fake_other=0;
       Fake_other_err=0;
       Fake_all=0;
       Fake_all_err=0;
       // with binomial error
       if(hS->GetBinContent(iCut) && hS_loose->GetBinContent(iCut)){
        Eff=hS->GetBinContent(iCut)/hS_loose->GetBinContent(iCut);
	Eff_err=(hS_loose->GetBinContent(iCut)-2*hS->GetBinContent(iCut))/hS_loose->GetBinContent(iCut)
	        *(hS_err->GetBinContent(iCut)/hS->GetBinContent(iCut)/hS->GetBinContent(iCut))
		+(hS_loose_err->GetBinContent(iCut)/hS_loose->GetBinContent(iCut)/hS_loose->GetBinContent(iCut));
	Eff_err=sqrt(Eff_err)*Eff;
	geff->SetPoint(iCut-2, iCut-1, Eff);
	geff->SetPointError(iCut-2, 0, Eff_err);
       }
       cout.precision(3);
       cout<<iCut-1<<"  Eff "<<Eff<<"+/-"<<Eff_err;
       if(hZ->GetBinContent(iCut) && hZ_loose->GetBinContent(iCut)){
        Fake_Z=hZ->GetBinContent(iCut)/hZ_loose->GetBinContent(iCut);
	Fake_Z_err=(hZ_loose->GetBinContent(iCut)-2*hZ->GetBinContent(iCut))/hZ_loose->GetBinContent(iCut)
	        *(hZ_err->GetBinContent(iCut)/hZ->GetBinContent(iCut)/hZ->GetBinContent(iCut))
		+(hZ_loose_err->GetBinContent(iCut)/hZ_loose->GetBinContent(iCut)/hZ_loose->GetBinContent(iCut));
	Fake_Z_err=sqrt(Fake_Z_err)*Fake_Z;
       }
       cout<<"\t  Fake Zjets "<<Fake_Z<<"+/-"<<Fake_Z_err;
       if(hother->GetBinContent(iCut) && hother_loose->GetBinContent(iCut)){
        Fake_other=hother->GetBinContent(iCut)/hother_loose->GetBinContent(iCut);
	Fake_other_err=(hother_loose->GetBinContent(iCut)-2*hother->GetBinContent(iCut))/hother_loose->GetBinContent(iCut)
	        *(hother_err->GetBinContent(iCut)/hother->GetBinContent(iCut)/hother->GetBinContent(iCut))
		+(hother_loose_err->GetBinContent(iCut)/hother_loose->GetBinContent(iCut)/hother_loose->GetBinContent(iCut));
	Fake_other_err=sqrt(Fake_other_err)*Fake_other;
       }
       cout<<"\t ttbar "<<Fake_other<<"+/-"<<Fake_other_err;
       if(hZ->GetBinContent(iCut)+hother->GetBinContent(iCut)>0 && hZ_loose->GetBinContent(iCut)+hother_loose->GetBinContent(iCut)){
        float a=hZ->GetBinContent(iCut)+hother->GetBinContent(iCut);
	float b=hZ_loose->GetBinContent(iCut)+hother_loose->GetBinContent(iCut);
        Fake_all=a/b;
	Fake_all_err=(b-2*a)/b
	        *(hZ_err->GetBinContent(iCut)+hother_err->GetBinContent(iCut))/a/a
		+(hZ_loose_err->GetBinContent(iCut)+hother_loose_err->GetBinContent(iCut))/b/b;
	Fake_all_err=sqrt(Fake_all_err)*Fake_all;
	gfake->SetPoint(iCut-2, iCut-1, Fake_all);
	gfake->SetPointError(iCut-2, 0, Fake_all_err);
       }
       cout<<"\t all "<<Fake_all<<"+/-"<<Fake_all_err<<endl;
     }
     
     c2.cd();
     geff->SetTitle(Form("Efficiencies from MC truth,  %s channel", channels[iChannel].c_str()));
     TH1F* hframe_eff = geff->GetHistogram();
     hframe_eff->GetYaxis()->SetRangeUser(0.75,1);
     geff->GetXaxis()->SetTitle("Cut");
     geff->SetMarkerStyle(25);
     geff->Draw("AP");
     c2.Modified();
     c2.Update();
     c2.Print(Form("MCtruth_Eff_%s_%s.pdf", channels[iChannel].c_str(), dataset.c_str()));
     
     c3.cd();
     gfake->SetTitle(Form("Tight to loose ratio from MC truth,  %s channel", channels[iChannel].c_str()));
     TH1F* hframe_fake = gfake->GetHistogram();
     if(iChannel==0 || iChannel==2) hframe_fake->GetYaxis()->SetRangeUser(0,0.2);
     else hframe_fake->GetYaxis()->SetRangeUser(0,0.5);
     gfake->GetXaxis()->SetTitle("Cut");
     gfake->SetMarkerStyle(25);
     gfake->Draw("AP"); 
     c3.Modified();
     c3.Update();
     c3.Print(Form("MCtruth_FakeRate_%s_%s.pdf", channels[iChannel].c_str(), dataset.c_str()));
   }
      
   getchar();

  }

}
