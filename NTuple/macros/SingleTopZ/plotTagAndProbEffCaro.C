#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"


double calculErrorSF(double a, double delta_a, double b, double delta_b){


     double error = pow(
	pow(delta_a/a, 2.) +
	pow(delta_b*a/(b*b), 2.) ,
	0.5);
	return error;
}



double calculSF(double a,  double b){


     double sfact = a/b;
     //cout << "a " << a  << "  b " << b << " SF is " << sfact << endl;
     return sfact;
}





TGraphAsymmErrors* Produce_Tgraph_caro(TH1F* tight_histo_Zjets, TH1F* looseID_histo_Zjets)
{
        // solution temporaire!
        float eff_id[100];
        float bin_center[100];
        float err_eff_id_m[100];
        float err_eff_id_p[100];
        float err_bin_m[100];
        float err_bin_p[100];
        for (int ixx=0; ixx<100; ixx++) {
          eff_id[ixx]=0;
          bin_center[ixx]=0;
          err_eff_id_m[ixx]=0;
          err_eff_id_p[ixx]=0;
          err_bin_m[ixx]=0;
          err_bin_p[ixx]=0;
        }

        int nbinxx=tight_histo_Zjets->GetNbinsX();
        if (nbinxx>100) {
           cout << " probleme de taille !!! tight_histo_Zjets a trop de bins! " << endl;
           nbinxx=100;
        }
        for (int ixx=0; ixx<nbinxx; ixx++) {
          float pass= tight_histo_Zjets->GetBinContent(ixx+1);
          float total1= looseID_histo_Zjets->GetBinContent(ixx+1);
          eff_id[ixx]=pass/total1;
//          float err1=sqrt(total1*eff_id[ixx]*(1-eff_id[ixx]))/total1;  //ceci c'est juste s'il n'y a pas de poids et si eff pas trop proche de 1!
//          malheureusement je n'ai pas la distribution non weightee mais si on suppose que le poids est le meme pour tous les events (style juste
//          une normalisation liee a la section efficace et la lumi --- on oublie un peu le pileup), on peut dire que 
//          N_sans_poids_bin_i = N_avec_poids_bin_i * (N_entries_histo/ N_integral_histo)
            float tt1=total1*(looseID_histo_Zjets->GetEntries()/looseID_histo_Zjets->Integral());
            float err1=sqrt(tt1*eff_id[ixx]*(1-eff_id[ixx]))/tt1;
           err_eff_id_m[ixx]=err1;
           err_eff_id_p[ixx]=err1;
           bin_center[ixx]=tight_histo_Zjets->GetBinCenter(ixx+1);
           err_bin_m[ixx]=tight_histo_Zjets->GetBinCenter(ixx+1)-tight_histo_Zjets->GetBinLowEdge(ixx+1);
           err_bin_p[ixx]=tight_histo_Zjets->GetBinLowEdge(ixx+1)-tight_histo_Zjets->GetBinCenter(ixx+1)+tight_histo_Zjets->GetBinWidth(ixx+1);
        }
        TGraphAsymmErrors* graph_caro= new TGraphAsymmErrors(nbinxx,bin_center,eff_id,err_bin_m,err_bin_p,err_eff_id_m,err_eff_id_p);
        return graph_caro;
}





std::vector<pair<double, double> > plotTagAndProbEff(TString lepton, TString q, TString isoid, TString etaCut, ofstream *ofile, bool computeSF, int systType){
  

  bool normal = false;  // utilise une solution temporaire pour le probleme de TGraghAsymmErrors sur du MC weightes

  std::vector<pair<double, double> > output;
  
  Int_t eventType = 1; // 0 MC, 1 Data
  
  Int_t stati=0;
  Bool_t  fit=1;
  Bool_t logy=0;
  
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  gStyle->SetOptDate(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleColor(1);
  //gStyle->SetTitleW(0.4);
  //gStyle->SetTitleH(0.07);
  //gStyle->SetOptStat(stati);
  gStyle->SetOptStat(0);
  // gStyle->SetFrameFillColor(18);
  //gStyle->SetPalette(51,0);
  gStyle->SetPalette(1);
  
  
  TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
   
  c1->SetFillColor(10);
  c1->SetFillStyle(4000);
  c1->SetBorderSize(2);
  
  //c1.Divide(2,2);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(10);
  
  
  
  TFile *f_Data_El ; 
  TFile *f_Data_Mu ; 
  TFile *f_MC_Zjets; 
  
  if(systType == 0){
    f_Data_El = new TFile("RootFileTagAndProb_good/2011AB/tagAndProbeHisto_DataEG.root"); 
    f_Data_Mu = new TFile("RootFileTagAndProb_good/2011AB/tagAndProbeHisto_DataMu.root"); 
    f_MC_Zjets= new TFile("RootFileTagAndProb_good/2011AB/tagAndProbeHisto_Zjets.root"); 
  }
  if(systType == -1){
    f_Data_El = new TFile("RootFileTagAndProb_good/2011AB_ZMm5/tagAndProbeHisto_DataEG.root"); 
    f_Data_Mu = new TFile("RootFileTagAndProb_good/2011AB_ZMm5/tagAndProbeHisto_DataMu.root"); 
    f_MC_Zjets= new TFile("RootFileTagAndProb_good/2011AB_ZMm5/tagAndProbeHisto_Zjets.root"); 
  }
  if(systType == 1){
    cout << "set syst " << endl;
    f_Data_El = new TFile("RootFileTagAndProb_good/2011AB_ZMp5/tagAndProbeHisto_DataEG.root"); 
    f_Data_Mu = new TFile("RootFileTagAndProb_good/2011AB_ZMp5/tagAndProbeHisto_DataMu.root"); 
    f_MC_Zjets= new TFile("RootFileTagAndProb_good/2011AB_ZMp5/tagAndProbeHisto_Zjets.root"); 
  }
  
  
  
  
  TH1F* tight_histo;
  TH1F* looseID_histo;
  TH1F* looseIso_histo;
  
  TH1F* tight_histo_LS;
  TH1F* looseID_histo_LS;
  TH1F* looseIso_histo_LS;
  
  TH2F* tight_histo_CP;
  TH2F* looseID_histo_CP;
  TH2F* looseIso_histo_CP;
  
  TH1F* tight_histo_TTbar    ;
  TH1F* looseID_histo_TTbar  ;
  TH1F* looseIso_histo_TTbar ;
   
   
  TH1F* tight_histo_Zjets    ;
  TH1F* looseID_histo_Zjets  ;
  TH1F* looseIso_histo_Zjets ;
      
   
  TLegend* qw = 0;
  qw = new TLegend(0.50,0.15,0.85,0.50);
 
  TGraphAsymmErrors * grID;
  TGraphAsymmErrors * grIso;
  TGraphAsymmErrors * grID_Zjets;
  TGraphAsymmErrors * grIso_Zjets;
  
  if(lepton == "Mu"){
    
    if(q == "pt"){
      f_Data_Mu->cd();
      tight_histo    = (TH1F*)gROOT->FindObject( ("DataMumuons_Tight_pt"+etaCut).Data() );
      looseID_histo  = (TH1F*)gROOT->FindObject( ("DataMumuons_looseID_pt"+etaCut).Data() );
      looseIso_histo = (TH1F*)gROOT->FindObject( ("DataMumuons_looseIso_pt"+etaCut).Data() );
     
      f_MC_Zjets->cd();
      tight_histo_Zjets    = (TH1F*)gROOT->FindObject( ("Zjetsmuons_Tight_pt"+etaCut).Data() );
      looseID_histo_Zjets  = (TH1F*)gROOT->FindObject( ("Zjetsmuons_looseID_pt"+etaCut).Data() );
      looseIso_histo_Zjets = (TH1F*)gROOT->FindObject( ("Zjetsmuons_looseIso_pt"+etaCut).Data() );
      
      
      grID  = new TGraphAsymmErrors();
      grID->BayesDivide(tight_histo,looseID_histo);
      grIso = new TGraphAsymmErrors();
      grIso->BayesDivide(tight_histo,looseIso_histo);
      
     
     
// probleme dans le cas de poids!

      if (normal==true) {
       grID_Zjets  = new TGraphAsymmErrors();
       grID_Zjets->BayesDivide(tight_histo_Zjets,looseID_histo_Zjets);
       grIso_Zjets = new TGraphAsymmErrors();
       grIso_Zjets->BayesDivide(tight_histo_Zjets,looseIso_histo_Zjets);
      } 
      else {
        // solution temporaire!
        grID_Zjets= Produce_Tgraph_caro(tight_histo_Zjets, looseID_histo_Zjets);
        grIso_Zjets=Produce_Tgraph_caro(tight_histo_Zjets, looseIso_histo_Zjets);
      }
      
      
      
      looseID_histo->GetXaxis()->SetTitle("muon p_{T} [GeV/c]");
      looseIso_histo->GetXaxis()->SetTitle("muon p_{T} [GeV/c]");
      
      //c1->SaveAs("muon_eff_pt.gif");
      if(computeSF){
        /*tight_histo->Sumw2();
        looseID_histo->Sumw2();
        tight_histo_Zjets->Sumw2();
        looseID_histo_Zjets->Sumw2();
	if(isoid == "id"){
	  cout << "in id " << endl;
          tight_histo->Divide(tight_histo,looseID_histo, 1, 1, "b");
          tight_histo_Zjets->Divide(tight_histo_Zjets,looseID_histo_Zjets, 1, 1, "b");
          tight_histo->Divide(tight_histo,tight_histo_Zjets, 1, 1, "b");
	}
	if(isoid == "iso"){
	  cout << "in iso " << endl;
          tight_histo->Divide(tight_histo,looseIso_histo, 1, 1, "b");
          tight_histo_Zjets->Divide(tight_histo_Zjets,looseIso_histo_Zjets, 1, 1, "b");
          tight_histo->Divide(tight_histo,tight_histo_Zjets, 1, 1, "b");
	}*/
      }
    } 
    
    
    if(q == "eta"){
      
      f_Data_Mu->cd();
      tight_histo    = (TH1F*)gROOT->FindObject(("DataMumuons_Tight_eta"+etaCut).Data());
      looseID_histo  = (TH1F*)gROOT->FindObject(("DataMumuons_looseID_eta"+etaCut).Data());
      looseIso_histo = (TH1F*)gROOT->FindObject(("DataMumuons_looseIso_eta"+etaCut).Data());
      
      
      f_MC_Zjets->cd();
      tight_histo_Zjets    = (TH1F*)gROOT->FindObject(("Zjetsmuons_Tight_eta"+etaCut).Data());
      looseID_histo_Zjets  = (TH1F*)gROOT->FindObject(("Zjetsmuons_looseID_eta"+etaCut).Data());
      looseIso_histo_Zjets = (TH1F*)gROOT->FindObject(("Zjetsmuons_looseIso_eta"+etaCut).Data());
     
     
      
      grID  = new TGraphAsymmErrors();
      grID->BayesDivide(tight_histo,looseID_histo);
      grIso = new TGraphAsymmErrors();
      grIso->BayesDivide(tight_histo,looseIso_histo);
      
      
      if (normal==true) {
       grID_Zjets  = new TGraphAsymmErrors();
       grID_Zjets->BayesDivide(tight_histo_Zjets,looseID_histo_Zjets);
       grIso_Zjets = new TGraphAsymmErrors();
       grIso_Zjets->BayesDivide(tight_histo_Zjets,looseIso_histo_Zjets);
      } 
      else {
        // solution temporaire!
        grID_Zjets= Produce_Tgraph_caro(tight_histo_Zjets, looseID_histo_Zjets);
        grIso_Zjets=Produce_Tgraph_caro(tight_histo_Zjets, looseIso_histo_Zjets);
      }
      
      looseID_histo->GetXaxis()->SetTitle("muon #eta");
      looseIso_histo->GetXaxis()->SetTitle("muon #eta");
      //c1->SaveAs("muon_eff_eta.gif");
    } 
    if(q == "phi"){
      
      f_Data_Mu->cd();
      tight_histo    = (TH1F*)gROOT->FindObject(("DataMumuons_Tight_phi"+etaCut).Data());
      looseID_histo  = (TH1F*)gROOT->FindObject(("DataMumuons_looseID_phi"+etaCut).Data());
      looseIso_histo = (TH1F*)gROOT->FindObject(("DataMumuons_looseIso_phi"+etaCut).Data());
          
      f_MC_Zjets->cd();
      tight_histo_Zjets    = (TH1F*)gROOT->FindObject(("Zjetsmuons_Tight_phi"+etaCut).Data());
      looseID_histo_Zjets  = (TH1F*)gROOT->FindObject(("Zjetsmuons_looseID_phi"+etaCut).Data());
      looseIso_histo_Zjets = (TH1F*)gROOT->FindObject(("Zjetsmuons_looseIso_phi"+etaCut).Data());
     
      
      grID  = new TGraphAsymmErrors();
      grID->BayesDivide(tight_histo,looseID_histo);
      grIso = new TGraphAsymmErrors();
      grIso->BayesDivide(tight_histo,looseIso_histo);
      
      
      if (normal==true) {
       grID_Zjets  = new TGraphAsymmErrors();
       grID_Zjets->BayesDivide(tight_histo_Zjets,looseID_histo_Zjets);
       grIso_Zjets = new TGraphAsymmErrors();
       grIso_Zjets->BayesDivide(tight_histo_Zjets,looseIso_histo_Zjets);
      } 
      else {
        // solution temporaire!
        grID_Zjets= Produce_Tgraph_caro(tight_histo_Zjets, looseID_histo_Zjets);
        grIso_Zjets=Produce_Tgraph_caro(tight_histo_Zjets, looseIso_histo_Zjets);
      }
      
      
      looseID_histo->GetXaxis()->SetTitle("muon #phi");
      looseIso_histo->GetXaxis()->SetTitle("muon #phi");
    } 
    if(q == "njet"){
      
      f_Data_Mu->cd();
      tight_histo    = (TH1F*)gROOT->FindObject(("DataMumuons_Tight_njet"+etaCut).Data());
      looseID_histo  = (TH1F*)gROOT->FindObject(("DataMumuons_looseID_njet"+etaCut).Data());
      looseIso_histo = (TH1F*)gROOT->FindObject(("DataMumuons_looseIso_njet"+etaCut).Data());
      
      
      f_MC_Zjets->cd();
      tight_histo_Zjets    = (TH1F*)gROOT->FindObject(("Zjetsmuons_Tight_njet"+etaCut).Data());
      looseID_histo_Zjets  = (TH1F*)gROOT->FindObject(("Zjetsmuons_looseID_njet"+etaCut).Data());
      looseIso_histo_Zjets = (TH1F*)gROOT->FindObject(("Zjetsmuons_looseIso_njet"+etaCut).Data());
     
      
      grID  = new TGraphAsymmErrors();
      grID->BayesDivide(tight_histo,looseID_histo);
      grIso = new TGraphAsymmErrors();
      grIso->BayesDivide(tight_histo,looseIso_histo);
      
      
      if (normal==true) {
       grID_Zjets  = new TGraphAsymmErrors();
       grID_Zjets->BayesDivide(tight_histo_Zjets,looseID_histo_Zjets);
       grIso_Zjets = new TGraphAsymmErrors();
       grIso_Zjets->BayesDivide(tight_histo_Zjets,looseIso_histo_Zjets);
      } 
      else {
        // solution temporaire!
        grID_Zjets= Produce_Tgraph_caro(tight_histo_Zjets, looseID_histo_Zjets);
        grIso_Zjets=Produce_Tgraph_caro(tight_histo_Zjets, looseIso_histo_Zjets);
      }
      
      
      looseID_histo->GetXaxis()->SetTitle("NJet");
      looseIso_histo->GetXaxis()->SetTitle("NJet");
    } 
    if(isoid == "id" ) qw->AddEntry(looseID_histo,  "muon ID eff. data","P");
    if(isoid == "iso") qw->AddEntry(looseIso_histo, "muon Iso eff. data","P");
    if(isoid == "id" ) qw->AddEntry(looseID_histo_Zjets,  "muon ID eff.  Z jets","P");
    if(isoid == "iso") qw->AddEntry(looseIso_histo_Zjets, "muon iso eff. Z jets","P");
  }
  
  
  
  if(lepton == "El"){
    
    if(q == "pt"){
      
      f_Data_El->cd();
      tight_histo    = (TH1F*)gROOT->FindObject(("DataEGelectrons_Tight_pt"+etaCut).Data());
      looseID_histo  = (TH1F*)gROOT->FindObject(("DataEGelectrons_looseID_pt"+etaCut).Data());
      looseIso_histo = (TH1F*)gROOT->FindObject(("DataEGelectrons_looseIso_pt"+etaCut).Data());
      
      tight_histo_LS    = (TH1F*)gROOT->FindObject("DataEGelectrons_Tight_pt_LS");
      looseID_histo_LS  = (TH1F*)gROOT->FindObject("DataEGelectrons_looseID_pt_LS");
      looseIso_histo_LS = (TH1F*)gROOT->FindObject("DataEGelectrons_looseIso_pt_LS");
      
      
      tight_histo_CP    = (TH2F*)gROOT->FindObject("DataEGelectrons_Tight_pt_CP");
      looseID_histo_CP  = (TH2F*)gROOT->FindObject("DataEGelectrons_looseID_pt_CP");
      looseIso_histo_CP = (TH2F*)gROOT->FindObject("DataEGelectrons_looseIso_pt_CP");
      
      
      /*f_MC_Zjets->cd();
      
      
      tight_histo    = (TH1F*)gROOT->FindObject("Zjetselectrons_Tight_pt");
      looseID_histo  = (TH1F*)gROOT->FindObject("Zjetselectrons_looseID_pt");
      looseIso_histo = (TH1F*)gROOT->FindObject("Zjetselectrons_looseIso_pt");
      
      tight_histo_LS    = (TH1F*)gROOT->FindObject("Zjetselectrons_Tight_pt_LS");
      looseID_histo_LS  = (TH1F*)gROOT->FindObject("Zjetselectrons_looseID_pt_LS");
      looseIso_histo_LS = (TH1F*)gROOT->FindObject("Zjetselectrons_looseIso_pt_LS");
      
      
      tight_histo_CP    = (TH2F*)gROOT->FindObject("Zjetselectrons_Tight_pt_CP");
      looseID_histo_CP  = (TH2F*)gROOT->FindObject("Zjetselectrons_looseID_pt_CP");
      looseIso_histo_CP = (TH2F*)gROOT->FindObject("Zjetselectrons_looseIso_pt_CP");
       */
      /*for(unsigned int i=1; i<=14; i++){
      
	cout << "*******************  "  << endl;
        double qqmean_looseIso = 0;
        qqmean_looseIso        = -1*looseIso_histo_CP->GetBinContent(i,1) + looseIso_histo_CP->GetBinContent(i,3);
	cout << " looseIso_histo_CP->GetBinContent(i,1) " <<  looseIso_histo_CP->GetBinContent(i,1) << endl;
	cout << " looseIso_histo_CP->GetBinContent(i,3) " <<  looseIso_histo_CP->GetBinContent(i,3)  << endl;
	cout << " qqmean_looseIso 1 " <<  qqmean_looseIso << endl;
	qqmean_looseIso        = qqmean_looseIso/(looseIso_histo->GetBinContent(i)+looseIso_histo_LS->GetBinContent(i));
	cout << " qqmean_looseIso 2 " <<  qqmean_looseIso << endl;
	
	double qmiss          = (1- pow(fabs(qqmean_looseIso), 0.5) ) / 2;
	cout << "qmiss            " <<  qmiss << endl;
	cout << "pow(1-qmiss, 2)  " <<  pow(1-qmiss, 2) << endl;
	double the_N_S        = (looseIso_histo->GetBinContent(i) - looseIso_histo_LS->GetBinContent(i))/( pow(1-qmiss, 2)); 
	cout << "looseIso_histo->GetBinContent(i) " << looseIso_histo->GetBinContent(i) << endl;
	looseIso_histo->SetBinContent(i,the_N_S);
	cout <<" the_N_S " << the_N_S << endl;
        
      } */
      
      
      /*f_MC_TTbar->cd();
      tight_histo_TTbar    = (TH1F*)gROOT->FindObject("genEl_Tight_pt_TTbar");
      looseID_histo_TTbar  = (TH1F*)gROOT->FindObject("genEl_LooseID_pt_TTbar");
      looseIso_histo_TTbar = (TH1F*)gROOT->FindObject("genEl_LooseIso_pt_TTbar");
      */
      f_MC_Zjets->cd();
      tight_histo_Zjets    = (TH1F*)gROOT->FindObject(("Zjetselectrons_Tight_pt"+etaCut).Data());
      looseID_histo_Zjets  = (TH1F*)gROOT->FindObject(("Zjetselectrons_looseID_pt"+etaCut).Data());
      looseIso_histo_Zjets = (TH1F*)gROOT->FindObject(("Zjetselectrons_looseIso_pt"+etaCut).Data());
      
      /*
      tight_histo->Rebin(4);
      looseID_histo->Rebin(4);
      looseIso_histo->Rebin(4);
      
      tight_histo_Zjets->Rebin(4);
      looseID_histo_Zjets->Rebin(4);
      looseIso_histo_Zjets->Rebin(4);*/
      grID  = new TGraphAsymmErrors();
      grID->BayesDivide(tight_histo,looseID_histo);
      grIso = new TGraphAsymmErrors();
      grIso->BayesDivide(tight_histo,looseIso_histo);
      
      
      
      if (normal==true) {
       grID_Zjets  = new TGraphAsymmErrors();
       grID_Zjets->BayesDivide(tight_histo_Zjets,looseID_histo_Zjets);
       grIso_Zjets = new TGraphAsymmErrors();
       grIso_Zjets->BayesDivide(tight_histo_Zjets,looseIso_histo_Zjets);
      } 
      else {
        // solution temporaire!
        grID_Zjets= Produce_Tgraph_caro(tight_histo_Zjets, looseID_histo_Zjets);
        grIso_Zjets=Produce_Tgraph_caro(tight_histo_Zjets, looseIso_histo_Zjets);
      }
      
      
      
      looseID_histo->GetXaxis()->SetTitle("electron p_{T} [GeV/c]");
      looseIso_histo->GetXaxis()->SetTitle("electron p_{T} [GeV/c]");
      
      tight_histo->GetXaxis()->SetTitle("electron p_{T} [GeV/c]");
      //c1->SaveAs("electron_eff_pt.gif");
      //compute scale factors
      
      
      
      
      if(computeSF){
        /*tight_histo->Sumw2();
        looseID_histo->Sumw2();
        tight_histo_Zjets->Sumw2();
        looseID_histo_Zjets->Sumw2();
	if(isoid == "id"){
	  
          tight_histo->Divide(tight_histo,looseID_histo, 1, 1, "b");
          tight_histo_Zjets->Divide(tight_histo_Zjets,looseID_histo_Zjets, 1, 1, "b");
          tight_histo->Divide(tight_histo,tight_histo_Zjets, 1, 1, "b");
	}
	if(isoid == "iso"){
          tight_histo->Divide(tight_histo,looseIso_histo, 1, 1, "b");
          tight_histo_Zjets->Divide(tight_histo_Zjets,looseIso_histo_Zjets, 1, 1, "b");
          tight_histo->Divide(tight_histo,tight_histo_Zjets, 1, 1, "b");
	}*/
      }
      
    } 
    if(q == "eta"){
      
      f_Data_El->cd();
      tight_histo    = (TH1F*)gROOT->FindObject(("DataEGelectrons_Tight_eta"+etaCut).Data());
      looseID_histo  = (TH1F*)gROOT->FindObject(("DataEGelectrons_looseID_eta"+etaCut).Data());
      looseIso_histo = (TH1F*)gROOT->FindObject(("DataEGelectrons_looseIso_eta"+etaCut).Data());
      
      tight_histo_LS    = (TH1F*)gROOT->FindObject("DataEGelectrons_Tight_eta_LS");
      looseID_histo_LS  = (TH1F*)gROOT->FindObject("DataEGelectrons_looseID_eta_LS");
      looseIso_histo_LS = (TH1F*)gROOT->FindObject("DataEGelectrons_looseIso_eta_LS");
      
      f_MC_Zjets->cd();
      tight_histo_Zjets    = (TH1F*)gROOT->FindObject(("Zjetselectrons_Tight_eta"+etaCut).Data());
      looseID_histo_Zjets  = (TH1F*)gROOT->FindObject(("Zjetselectrons_looseID_eta"+etaCut).Data());
      looseIso_histo_Zjets = (TH1F*)gROOT->FindObject(("Zjetselectrons_looseIso_eta"+etaCut).Data());
      
      
      
      /*tight_histo->Rebin(4);
      looseID_histo->Rebin(4);
      looseIso_histo->Rebin(4);
      
      tight_histo_Zjets->Rebin(4);
      looseID_histo_Zjets->Rebin(4);
      looseIso_histo_Zjets->Rebin(4);*/
      
      grID  = new TGraphAsymmErrors();
      grID->BayesDivide(tight_histo,looseID_histo);
      grIso = new TGraphAsymmErrors();
      grIso->BayesDivide(tight_histo,looseIso_histo);
      
      if (normal==true) {
       grID_Zjets  = new TGraphAsymmErrors();
       grID_Zjets->BayesDivide(tight_histo_Zjets,looseID_histo_Zjets);
       grIso_Zjets = new TGraphAsymmErrors();
       grIso_Zjets->BayesDivide(tight_histo_Zjets,looseIso_histo_Zjets);
      } 
      else {
        // solution temporaire!
        grID_Zjets= Produce_Tgraph_caro(tight_histo_Zjets, looseID_histo_Zjets);
        grIso_Zjets=Produce_Tgraph_caro(tight_histo_Zjets, looseIso_histo_Zjets);
      }
      
      
      looseID_histo->GetXaxis()->SetTitle("electron #eta");
      looseIso_histo->GetXaxis()->SetTitle("electron #eta");
      tight_histo->GetXaxis()->SetTitle("electron #eta");
      //c1->SaveAs("electron_eff_eta.gif"); 
      /*if(computeSF){
        tight_histo->Sumw2();
        looseID_histo->Sumw2();
        tight_histo_Zjets->Sumw2();
        looseID_histo_Zjets->Sumw2();
        tight_histo->Divide(tight_histo,looseID_histo, 1, 1, "b");
        tight_histo_Zjets->Divide(tight_histo_Zjets,looseID_histo_Zjets, 1, 1, "b");
        tight_histo->Divide(tight_histo,tight_histo_Zjets, 1, 1, "b");
      }*/
    } 
    if(q == "phi"){
      
      f_Data_El->cd();
      tight_histo    = (TH1F*)gROOT->FindObject(("DataEGelectrons_Tight_phi"+etaCut).Data());
      looseID_histo  = (TH1F*)gROOT->FindObject(("DataEGelectrons_looseID_phi"+etaCut).Data());
      looseIso_histo = (TH1F*)gROOT->FindObject(("DataEGelectrons_looseIso_phi"+etaCut).Data());
      
      tight_histo_LS    = (TH1F*)gROOT->FindObject("DataEGelectrons_Tight_phi_LS");
      looseID_histo_LS  = (TH1F*)gROOT->FindObject("DataEGelectrons_looseID_phi_LS");
      looseIso_histo_LS = (TH1F*)gROOT->FindObject("DataEGelectrons_looseIso_phi_LS");
        
      f_MC_Zjets->cd();
      tight_histo_Zjets    = (TH1F*)gROOT->FindObject(("Zjetselectrons_Tight_phi"+etaCut).Data());
      looseID_histo_Zjets  = (TH1F*)gROOT->FindObject(("Zjetselectrons_looseID_phi"+etaCut).Data());
      looseIso_histo_Zjets = (TH1F*)gROOT->FindObject(("Zjetselectrons_looseIso_phi"+etaCut).Data());
      
      /*tight_histo->Rebin(4);
      looseID_histo->Rebin(4);
      looseIso_histo->Rebin(4);
      
      tight_histo_Zjets->Rebin(4);
      looseID_histo_Zjets->Rebin(4);
      looseIso_histo_Zjets->Rebin(4);*/
      
      grID  = new TGraphAsymmErrors();
      grID->BayesDivide(tight_histo,looseID_histo);
      grIso = new TGraphAsymmErrors();
      grIso->BayesDivide(tight_histo,looseIso_histo);
      
      
      if (normal==true) {
       grID_Zjets  = new TGraphAsymmErrors();
       grID_Zjets->BayesDivide(tight_histo_Zjets,looseID_histo_Zjets);
       grIso_Zjets = new TGraphAsymmErrors();
       grIso_Zjets->BayesDivide(tight_histo_Zjets,looseIso_histo_Zjets);
      } 
      else {
        // solution temporaire!
        grID_Zjets= Produce_Tgraph_caro(tight_histo_Zjets, looseID_histo_Zjets);
        grIso_Zjets=Produce_Tgraph_caro(tight_histo_Zjets, looseIso_histo_Zjets);
      }
      
      
      looseID_histo->GetXaxis()->SetTitle("electron #phi");
      looseIso_histo->GetXaxis()->SetTitle("electron #phi");
      tight_histo->GetXaxis()->SetTitle("electron #phi");
      //c1->SaveAs("electron_eff_phi.gif");
      /*if(computeSF){
        
	
	
        tight_histo->Sumw2();
        looseID_histo->Sumw2();
        tight_histo_Zjets->Sumw2();
        looseID_histo_Zjets->Sumw2();
        tight_histo->Divide(tight_histo,looseID_histo, 1, 1, "b");
        tight_histo_Zjets->Divide(tight_histo_Zjets,looseID_histo_Zjets, 1, 1, "b");
        tight_histo->Divide(tight_histo,tight_histo_Zjets, 1, 1, "b");
      }*/
    } 
    if(q == "njet"){
      
      f_Data_El->cd();
      tight_histo    = (TH1F*)gROOT->FindObject(("DataEGelectrons_Tight_njet"+etaCut).Data());
      looseID_histo  = (TH1F*)gROOT->FindObject(("DataEGelectrons_looseID_njet"+etaCut).Data());
      looseIso_histo = (TH1F*)gROOT->FindObject(("DataEGelectrons_looseIso_njet"+etaCut).Data());
      
      
      f_MC_Zjets->cd();
      tight_histo_Zjets    = (TH1F*)gROOT->FindObject(("Zjetselectrons_Tight_njet"+etaCut).Data());
      looseID_histo_Zjets  = (TH1F*)gROOT->FindObject(("Zjetselectrons_looseID_njet"+etaCut).Data());
      looseIso_histo_Zjets = (TH1F*)gROOT->FindObject(("Zjetselectrons_looseIso_njet"+etaCut).Data());
      
      
      grID  = new TGraphAsymmErrors();
      grID->BayesDivide(tight_histo,looseID_histo);
      grIso = new TGraphAsymmErrors();
      grIso->BayesDivide(tight_histo,looseIso_histo);
      
      
      if (normal==true) {
       grID_Zjets  = new TGraphAsymmErrors();
       grID_Zjets->BayesDivide(tight_histo_Zjets,looseID_histo_Zjets);
       grIso_Zjets = new TGraphAsymmErrors();
       grIso_Zjets->BayesDivide(tight_histo_Zjets,looseIso_histo_Zjets);
      } 
      else {
        // solution temporaire!
        grID_Zjets= Produce_Tgraph_caro(tight_histo_Zjets, looseID_histo_Zjets);
        grIso_Zjets=Produce_Tgraph_caro(tight_histo_Zjets, looseIso_histo_Zjets);
      }
      
      
      
      
      
      looseID_histo->GetXaxis()->SetTitle("NJet");
      looseIso_histo->GetXaxis()->SetTitle("NJet");
      //c1->SaveAs("electron_eff_njet.gif");
    } 
    if(isoid == "id" ) qw->AddEntry(looseID_histo,  "electron ID eff. data","P");
    if(isoid == "iso") qw->AddEntry(looseIso_histo, "electron iso eff.data","P");
    //if(isoid == "id" ) qw->AddEntry(looseID_histo_TTbar,  "electron ID eff. in TT MC","P");
    //if(isoid == "iso") qw->AddEntry(looseIso_histo_TTbar, "electron iso eff. in TT MC","P");
    if(isoid == "id" ) qw->AddEntry(looseID_histo_Zjets,  "electron ID eff.  Z jets","P");
    if(isoid == "iso") qw->AddEntry(looseIso_histo_Zjets, "electron iso eff.  Z jets","P");
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  tight_histo->Sumw2();
  looseID_histo->Sumw2();
  looseIso_histo->Sumw2();
  
  looseID_histo->Divide( tight_histo, looseID_histo, 1, 1, "b");
  looseIso_histo->Divide(tight_histo, looseIso_histo, 1, 1, "b");
  
  looseID_histo->SetMarkerStyle(20);
  looseIso_histo->SetMarkerStyle(20);
  looseID_histo->SetMarkerColor(2);
  looseIso_histo->SetMarkerColor(4);
  
  //looseID_histo_TTbar->SetMarkerStyle(24);
  //looseIso_histo_TTbar->SetMarkerStyle(24);
  //looseID_histo_TTbar->SetMarkerColor(2);
  //looseIso_histo_TTbar->SetMarkerColor(4);
  
  looseID_histo_Zjets->SetMarkerStyle(26);
  looseIso_histo_Zjets->SetMarkerStyle(26);
  looseID_histo_Zjets->SetMarkerColor(2);
  looseIso_histo_Zjets->SetMarkerColor(4);
    
  looseID_histo->SetMaximum(1.1);
  looseID_histo->SetMinimum(0.);
  looseIso_histo->SetMaximum(1.1);
  looseIso_histo->SetMinimum(0.);
  if(isoid == "id" ) looseID_histo->Draw("p");
  if(isoid == "iso") looseIso_histo->Draw("P");
  
  grID->SetMarkerStyle(20);
  grIso->SetMarkerStyle(20);
  grID->SetMarkerColor(2);
  grIso->SetMarkerColor(4);
  
  
  //grID_TTbar->SetMarkerStyle(24);
  //grIso_TTbar->SetMarkerStyle(24);
  //grID_TTbar->SetMarkerColor(2);
  //grIso_TTbar->SetMarkerColor(4);
  
  
  grID_Zjets->SetMarkerStyle(26);
  grIso_Zjets->SetMarkerStyle(26);
  grID_Zjets->SetMarkerColor(2);
  grIso_Zjets->SetMarkerColor(4);
  
  
  if(isoid == "iso")   grIso->Draw("psame");
  if(isoid == "id" )   grID->Draw("psame");
  if(isoid == "id" )   grID_Zjets->Draw("psame");
  if(isoid == "iso")   grIso_Zjets->Draw("psame");
   
  
  qw->SetFillColor(0);
  qw->Draw();
  TString outputfilename =  lepton+"_"+ q+"_"+isoid+".gif";
  
  c1->SaveAs(outputfilename.Data());
  
 
  if(computeSF){
     //double ptRange[] = {20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 100};
     double ptRange[] = {20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 100};
     TH1F * thehistSF ;
     if (q=="pt") {
         thehistSF = new TH1F("thehistSF", "thehistSF", 10, ptRange);
     }
     else  {
         int nnbbiinn=tight_histo->GetNbinsX();
         float xmin=tight_histo->GetBinLowEdge(1);
         float xmax=tight_histo->GetBinLowEdge(nnbbiinn)+tight_histo->GetBinWidth(nnbbiinn);
         cout << q << " nbin " << nnbbiinn << " xmin " << xmin  << " xmax" << xmax <<endl;
         thehistSF = new TH1F("thehistSF", "thehistSF", nnbbiinn, xmin, xmax);
     }
  
     for(unsigned int i=0; i<= looseIso_histo->GetNbinsX()-1; i++){
	  std::pair<double, double> thepair;
       if(isoid == "iso" ){
         double eff_data = 0;
         double x_data   = 0;
	 grIso->GetPoint(i, x_data, eff_data);
	 
	 //cout << "iso x_data " << i << " " << x_data << " iso eff_data " << eff_data << endl;
         double eff_Z  = 0;
         double x_Z    = 0;
	 grIso_Zjets->GetPoint(i, x_Z,eff_Z );
	 //cout << "iso x_Z    " << i << " " << x_Z    << " iso eff_Z    " << eff_Z    << endl;
	 //cout << "calculSF(eff_data, eff_Z)  " << calculSF(eff_data, eff_Z) << endl;
         //*ofile <<   calculSF(eff_data, eff_Z)
         // << " $\pm$ " <<  calculErrorSF(eff_data, grIso->GetErrorY(i), eff_Z, grIso_Zjets->GetErrorY(i)) <<  " & " ;
	  thepair.first  = calculSF(eff_data, eff_Z);
	  thepair.second = calculErrorSF(eff_data, grIso->GetErrorY(i), eff_Z, grIso_Zjets->GetErrorY(i));
	  output.push_back(thepair);
	  thehistSF->SetBinContent(i+1, thepair.first);
	  thehistSF->SetBinError(i+1, thepair.second);
	}
       
       if(isoid == "id" ){
         double eff_data = 0;
         double x_data   = 0;
	 grID->GetPoint(i, x_data, eff_data);
	 //cout << "id x_data " << i << " "  << x_data << " id eff_data " << eff_data << endl;
         double eff_Z  = 0;
         double x_Z    = 0;
	 grID_Zjets->GetPoint(i, x_Z,eff_Z );
	 //cout << "id x_Z    " << i << " "  << x_Z    << " id eff_Z    " << eff_Z    << endl;
	 
         //*ofile <<   calculSF(eff_data, eff_Z)
          //<< " $\pm$ " <<  calculErrorSF(eff_data, grID->GetErrorY(i), eff_Z, grID_Zjets->GetErrorY(i)) <<  " & " ;
	  thepair.first  = calculSF(eff_data, eff_Z);
	  thepair.second = calculErrorSF(eff_data, grID->GetErrorY(i), eff_Z, grID_Zjets->GetErrorY(i));
	  output.push_back(thepair);
	  thehistSF->SetBinContent(i+1, thepair.first);
	  thehistSF->SetBinError(i+1, thepair.second);
	}
       
     } 
   thehistSF->SetMaximum(1.1);
   TString outputfilename =  lepton+"_SF_"+ q+"_"+isoid+".gif";
  
   if(isoid == "id" ){ thehistSF->GetYaxis()->SetTitle("SF_{ID}");   }
   if(isoid == "iso" ){ thehistSF->GetYaxis()->SetTitle("SF_{Iso}"); }
   thehistSF->SetMarkerStyle(20);
 
   thehistSF->SetMaximum(1.1);
   thehistSF->SetMinimum(0.85);
   if(lepton == "El") {
       if (q=="pt") thehistSF->GetXaxis()->SetTitle("electron p_{T} [GeV/c]");
       else if (q=="eta") thehistSF->GetXaxis()->SetTitle("electron #eta");
       else if (q=="phi") thehistSF->GetXaxis()->SetTitle("electron #phi");
   }
   else if(lepton == "Mu") {
       if (q=="pt") thehistSF->GetXaxis()->SetTitle("muon p_{T} [GeV/c]");
       else if (q=="eta") thehistSF->GetXaxis()->SetTitle("muon #eta");
       else if (q=="phi") thehistSF->GetXaxis()->SetTitle("muon #phi");
   }

   thehistSF->Draw("ep");
   c1->SaveAs(outputfilename.Data());
  
  }
     
  
  
  if(!computeSF){
  
     //double ptRange[] = {20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 100};
     for(unsigned int i=0; i<= looseIso_histo->GetNbinsX()-1; i++){
       
       if(isoid == "iso" ){
        //*ofile <<   looseIso_histo->GetBinContent(i) << " $\pm$ " <<  looseIso_histo->GetBinError(i) <<  " & " ;
         double eff_data = 0;
         double x_data   = 0;
	 grIso->GetPoint(i, x_data, eff_data);
	 //*ofile <<  eff_data  << " $\pm$ " << grIso->GetErrorY(i)  <<  " & " ;
	 std::pair<double, double> thepair;
	 thepair.first  = eff_data;
	 thepair.second = grIso->GetErrorY(i) ;
	 output.push_back(thepair);
       }
       if(isoid == "id" ){
         //*ofile <<   looseID_histo->GetBinContent(i)  << " $\pm$ " <<  looseID_histo->GetBinError(i) <<  " & " ;
         double eff_data = 0;
         double x_data   = 0;
	 grID->GetPoint(i, x_data, eff_data);
	 //*ofile <<  eff_data  << " $\pm$ " << grID->GetErrorY(i)  <<  " & " ;
	 std::pair<double, double> thepair;
	 thepair.first  = eff_data;
	 thepair.second = grID->GetErrorY(i) ;
	 output.push_back(thepair);
       }
       
     }
     
     
  }
  
  
  
  
/*
  delete tight_histo;
  delete looseID_histo;
  delete looseIso_histo;
  delete tight_histo_LS;
  delete looseID_histo_LS;
  delete looseIso_histo_LS;
  delete tight_histo_CP;
  delete looseID_histo_CP;
  delete looseIso_histo_CP;
  delete tight_histo_Zjets    ;
  delete looseID_histo_Zjets  ;
  delete looseIso_histo_Zjets ;
  delete grID;
  delete grIso;
  delete grID_Zjets;
  delete grIso_Zjets;
  delete c1;
  delete f_Data_El ; 
  delete f_Data_Mu ; 
  delete f_MC_Zjets; 
*/
  tight_histo->Delete();
  looseID_histo->Delete();
  looseIso_histo->Delete();
/*
  tight_histo_LS->Delete();
  looseID_histo_LS->Delete();
  looseIso_histo_LS->Delete();
  tight_histo_CP->Delete();
  looseID_histo_CP->Delete();
  looseIso_histo_CP->Delete();
*/
  tight_histo_Zjets->Delete();
  looseID_histo_Zjets->Delete();
  looseIso_histo_Zjets->Delete();
  grID->Delete();
  grIso->Delete();
  grIso_Zjets->Delete();
  grID_Zjets->Delete();
//  c1->Delete();
  f_Data_El->Close();
  f_Data_Mu->Close();
  f_MC_Zjets->Close();

  
  return output;
  
}


void dumpTagAndProbeValue(ofstream *ofile, std::vector<pair<double, double > > norm, 
std::vector<pair<double, double > > sup, std::vector<pair<double, double > > sdown , bool part1){
  
  if(part1){
   for(unsigned int ibin = 0; ibin < 5; ibin++){
    
    double systUp   = fabs(norm[ibin].first - sup[ibin].first);
    double systDown = fabs(norm[ibin].first - sdown[ibin].first);
    double syst = systUp;
    if(systUp<systDown) syst = systDown;
    *ofile << norm[ibin].first   << " $\\pm$ " << norm[ibin].second << "$\\pm$"  <<  syst ;
    if(ibin < 4 ) *ofile<< " & " ;
    
    
   }
  }else{
  
   for(unsigned int ibin = 5; ibin < norm.size(); ibin++){
    
    double systUp   = fabs(norm[ibin].first - sup[ibin].first);
    double systDown = fabs(norm[ibin].first - sdown[ibin].first);
    double syst = systUp;
    if(systUp<systDown) syst = systDown;
    *ofile << norm[ibin].first   << " $\\pm$ " << norm[ibin].second << "$\\pm$"  <<  syst ;
    if(ibin < (norm.size()-1) ) *ofile<< " & " ;
    
    
   }
  
  }
  
    *ofile << " \\\\ " << endl;


}
void dumpTagAndProbeValueRoot(
std::vector<pair<double, double > > norm_inf15,  std::vector<pair<double, double > > sup_inf15, std::vector<pair<double, double > > sdown_inf15, 
std::vector<pair<double, double > > norm_sup15,  std::vector<pair<double, double > > sup_sup15, std::vector<pair<double, double > > sdown_sup15, 
bool isMuon ){



  double ptRange[]  = {20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 100};
  double etaRange[] = {0., 1.5, 2.5};
  TH2F * thehistSF = new TH2F("thehistSF", "thehistSF", 10, ptRange, 2, etaRange);
  

  for(unsigned int ibin = 0; ibin < norm_inf15.size(); ibin++){
    
    double systUp   = fabs(norm_inf15[ibin].first - sup_inf15[ibin].first);
    double systDown = fabs(norm_inf15[ibin].first - sdown_inf15[ibin].first);
    double syst = systUp;
    if(systUp<systDown) syst = systDown;
    thehistSF->SetBinContent(ibin+1, 1, norm_inf15[ibin].first);
    thehistSF->SetBinError(ibin+1, 1,syst);
    //cout << "norm_inf15[ibin].first " << norm_inf15[ibin].first << endl;
  }
  
  
  for(unsigned int ibin = 0; ibin < norm_sup15.size(); ibin++){
    
    double systUp   = fabs(norm_sup15[ibin].first - sup_sup15[ibin].first);
    double systDown = fabs(norm_sup15[ibin].first - sdown_sup15[ibin].first);
    double syst = systUp;
    if(systUp<systDown) syst = systDown;
    thehistSF->SetBinContent(ibin+1, 2, norm_sup15[ibin].first);
    thehistSF->SetBinError(ibin+1, 2,syst);
    //cout << "norm_sup15[ibin].first " << norm_sup15[ibin].first << endl;
  }
  
  
  TFile * outputSF;
  if(isMuon) outputSF = new TFile("muonSF.root", "recreate");
  else       outputSF = new TFile("electronSF.root", "recreate");
  outputSF->cd();;
  thehistSF->Write();




}

std::vector<pair<double, double > > multiplySF(std::vector<pair<double, double > > vect1, std::vector<pair<double, double > > vect2){
  
  std::vector<pair<double, double > > globalSF;
  for(unsigned int i=0; i< vect1.size(); i++){
    
    
	 std::pair<double, double> thepair;
	 thepair.first  = vect1[i].first*vect2[i].first;
	 thepair.second = pow(pow(vect1[i].second,2.)+ pow(vect2[i].second,2.), 0.5);
         globalSF.push_back(thepair);
  }
  return globalSF;


}

void ProducePlots(){
  
  string ofilenametex = string("Plots")+string(".tex");
  ofstream ofile(ofilenametex.c_str());
 
//  std::vector<pair<double, double > > elecIsoEff_pt = plotTagAndProbEff("El", "pt", "iso","" , &ofile, 0, 0);
//  std::vector<pair<double, double > > elecIDEff_pt  = plotTagAndProbEff("El", "pt", "id",""  , &ofile, 0, 0);
  std::vector<pair<double, double > > elecIsoSF_pt  = plotTagAndProbEff("El", "pt", "iso","" , &ofile, 1, 0);
  std::vector<pair<double, double > > elecIDSF_pt   = plotTagAndProbEff("El", "pt", "id",""  , &ofile, 1, 0);
  
//  std::vector<pair<double, double > > muonIsoEff_pt = plotTagAndProbEff("Mu", "pt", "iso","" , &ofile, 0, 0);
//  std::vector<pair<double, double > > muonIDEff_pt  = plotTagAndProbEff("Mu", "pt", "id",""  , &ofile, 0, 0);
  std::vector<pair<double, double > > muonIsoSF_pt  = plotTagAndProbEff("Mu", "pt", "iso","" , &ofile, 1, 0);
  std::vector<pair<double, double > > muonIDSF_pt   = plotTagAndProbEff("Mu", "pt", "id",""  , &ofile, 1, 0);
  
  
 
//  std::vector<pair<double, double > > elecIsoEff_eta = plotTagAndProbEff("El", "eta", "iso","" , &ofile, 0, 0);
//  std::vector<pair<double, double > > elecIDEff_eta  = plotTagAndProbEff("El", "eta", "id",""  , &ofile, 0, 0);
  std::vector<pair<double, double > > elecIsoSF_eta  = plotTagAndProbEff("El", "eta", "iso","" , &ofile, 1, 0);
  std::vector<pair<double, double > > elecIDSF_eta   = plotTagAndProbEff("El", "eta", "id",""  , &ofile, 1, 0);
  
//  std::vector<pair<double, double > > muonIsoEff_eta = plotTagAndProbEff("Mu", "eta", "iso","" , &ofile, 0, 0);
//  std::vector<pair<double, double > > muonIDEff_eta  = plotTagAndProbEff("Mu", "eta", "id",""  , &ofile, 0, 0);
  std::vector<pair<double, double > > muonIsoSF_eta  = plotTagAndProbEff("Mu", "eta", "iso","" , &ofile, 1, 0);
  std::vector<pair<double, double > > muonIDSF_eta   = plotTagAndProbEff("Mu", "eta", "id",""  , &ofile, 1, 0);

  //JLA
  std::vector<pair<double, double > >  elecSF_eta = multiplySF(elecIsoSF_eta, elecIDSF_eta);
  std::cout<<"Electron SF vs eta"<<std::endl;
  for(unsigned int iv=0; iv<elecSF_eta.size(); iv++)
  {
    std::cout<<iv<<" "<<elecSF_eta[iv].first<<" +/- "<<elecSF_eta[iv].second<<std::endl;
  }
  
  //JLA
  std::vector<pair<double, double > >  muonSF_eta = multiplySF(muonIsoSF_eta, muonIDSF_eta);
  std::cout<<"Muon SF vs eta"<<std::endl;
  for(unsigned int iv=0; iv<muonSF_eta.size(); iv++)
  {
    std::cout<<iv<<" "<<muonSF_eta[iv].first<<" +/- "<<muonSF_eta[iv].second<<std::endl;
  }
  
 
//  std::vector<pair<double, double > > elecIsoEff_phi = plotTagAndProbEff("El", "phi", "iso","" , &ofile, 0, 0);
//  std::vector<pair<double, double > > elecIDEff_phi  = plotTagAndProbEff("El", "phi", "id",""  , &ofile, 0, 0);
  std::vector<pair<double, double > > elecIsoSF_phi  = plotTagAndProbEff("El", "phi", "iso","" , &ofile, 1, 0);
  std::vector<pair<double, double > > elecIDSF_phi   = plotTagAndProbEff("El", "phi", "id",""  , &ofile, 1, 0);
  
//  std::vector<pair<double, double > > muonIsoEff_phi = plotTagAndProbEff("Mu", "phi", "iso","" , &ofile, 0, 0);
//  std::vector<pair<double, double > > muonIDEff_phi  = plotTagAndProbEff("Mu", "phi", "id",""  , &ofile, 0, 0);
  std::vector<pair<double, double > > muonIsoSF_phi  = plotTagAndProbEff("Mu", "phi", "iso","" , &ofile, 1, 0);
  std::vector<pair<double, double > > muonIDSF_phi   = plotTagAndProbEff("Mu", "phi", "id",""  , &ofile, 1, 0);

}


//**************************
// get global SF info
//**************************

vector<float> get_sfinfo(string chan, string type_id){
  TFile *f_Data_El ; 
  TFile *f_Data_Mu ; 
  TFile *f_MC_Zjets; 
  
  f_Data_El = new TFile("RootFileTagAndProb_good/2011AB/tagAndProbeHisto_DataEG.root"); 
  f_Data_Mu = new TFile("RootFileTagAndProb_good/2011AB/tagAndProbeHisto_DataMu.root"); 
  f_MC_Zjets= new TFile("RootFileTagAndProb_good/2011AB/tagAndProbeHisto_Zjets.root"); 

 /* 
  if(systType == -1){
    f_Data_El = new TFile("RootFileTagAndProb/2011A/tagAndProbeHisto_DataEG.root"); 
    f_Data_Mu = new TFile("RootFileTagAndProb/tagAndProbeHisto_DataMu.root"); 
    f_MC_Zjets= new TFile("RootFileTagAndProb/tagAndProbeHisto_Zjets.root"); 
  }
  if(systType == 1){
    cout << "set syst " << endl;
    f_Data_El = new TFile("RootFileTagAndProb/tagAndProbeHisto_DataEG.root"); 
    f_Data_Mu = new TFile("RootFileTagAndProb/tagAndProbeHisto_DataMu.root"); 
    f_MC_Zjets= new TFile("RootFileTagAndProb/tagAndProbeHisto_Zjets.root"); 
  }*/

  
  
  TH1F* tight_histo;
  TH1F* looseID_histo;
  TH1F* looseIso_histo;
  
  TH1F* tight_histo_Zjets    ;
  TH1F* looseID_histo_Zjets  ;
  TH1F* looseIso_histo_Zjets ;
      
   
  if(chan == "mu"){
      f_Data_Mu->cd();
      tight_histo    = (TH1F*)gROOT->FindObject( "DataMumuons_Tight_pt" );
      looseID_histo  = (TH1F*)gROOT->FindObject( "DataMumuons_looseID_pt" );
      looseIso_histo = (TH1F*)gROOT->FindObject( "DataMumuons_looseIso_pt" );
     
      f_MC_Zjets->cd();
      tight_histo_Zjets    = (TH1F*)gROOT->FindObject( "Zjetsmuons_Tight_pt" );
      looseID_histo_Zjets  = (TH1F*)gROOT->FindObject( "Zjetsmuons_looseID_pt" );
      looseIso_histo_Zjets = (TH1F*)gROOT->FindObject( "Zjetsmuons_looseIso_pt" );
  }
  if (chan == "el") {
      f_Data_El->cd();
      tight_histo    = (TH1F*)gROOT->FindObject( "DataEGelectrons_Tight_pt" );
      looseID_histo  = (TH1F*)gROOT->FindObject( "DataEGelectrons_looseID_pt" );
      looseIso_histo = (TH1F*)gROOT->FindObject( "DataEGelectrons_looseIso_pt" );
      
      f_MC_Zjets->cd();
      tight_histo_Zjets    = (TH1F*)gROOT->FindObject( "Zjetselectrons_Tight_pt" );
      looseID_histo_Zjets  = (TH1F*)gROOT->FindObject( "Zjetselectrons_looseID_pt" );
      looseIso_histo_Zjets = (TH1F*)gROOT->FindObject( "Zjetselectrons_looseIso_pt" );
  }
      
    float pass_data=tight_histo->GetEntries();
    float total_data=1.;
    if (type_id=="id") {total_data=looseID_histo->GetEntries(); }
    else if (type_id=="iso") {total_data=looseIso_histo->GetEntries(); }
    float eff_data = pass_data/total_data;
    float err_data = sqrt(total_data*eff_data*(1.-eff_data))/total_data;

    float pass_Zjets=tight_histo_Zjets->Integral();
    float total_Zjets=1.;
    if (type_id=="id") {total_Zjets=looseID_histo_Zjets->Integral(); }
    else if (type_id=="iso") {total_Zjets=looseIso_histo_Zjets->Integral(); }
    float eff_Zjets = pass_Zjets/total_Zjets;
    float tt_Zjets=tight_histo_Zjets->GetEntries(); // sans poids!
    float err_Zjets = sqrt(tt_Zjets*eff_Zjets*(1.-eff_Zjets))/tt_Zjets;
    
    float SF_data_Zjets= eff_data/eff_Zjets;
    float err_SF = pow( pow(err_data/eff_data, 2.) + pow(err_Zjets*eff_data/(eff_Zjets*eff_Zjets),2) , 0.5 );

  
   cout << " eff data : " << eff_data << " eff Zjets : " << eff_Zjets << " SF : " << SF_data_Zjets << endl;    
    vector<float> info_to_be_returned;
    info_to_be_returned.push_back(SF_data_Zjets);
    info_to_be_returned.push_back(err_SF);
    return info_to_be_returned;
}

//****************************
// Write the global SF to be put in the note
//****************************

void PrintSummary(){

// Caro

  vector<float>  sf_mu_id = get_sfinfo("mu","id");
  vector<float>  sf_e_id = get_sfinfo("el","id");

  vector<float>  sf_mu_iso = get_sfinfo("mu","iso");
  vector<float>  sf_e_iso = get_sfinfo("el","iso");

  cout << " SINGLE LEPTON SCALE FACTOR " << endl;
  cout << " ID " << endl;
  cout << " mu " << sf_mu_id[0] << " +/- " << sf_mu_id[1] << endl;
  cout << " e  " << sf_e_id[0] << " +/- " << sf_e_id[1] << endl;
  cout << endl;
  cout << " ISO " << endl;
  cout << " mu " << sf_mu_iso[0] << " +/- " << sf_mu_iso[1] << endl;
  cout << " e  " << sf_e_iso[0] << " +/- " << sf_e_iso[1] << endl;
  cout << " --------------- " << endl;


  cout << " DOUBLE LEPTON SCALE FACTOR " << endl;
  float sf_mumu_id[2];
  float sf_ee_id[2];
  float sf_emu_id[2];

  sf_mumu_id[0]=sf_mu_id[0]*sf_mu_id[0];
  sf_ee_id[0]  =sf_e_id[0]*sf_e_id[0];
  sf_emu_id[0] =sf_mu_id[0]*sf_e_id[0];
  sf_mumu_id[1]=pow(pow(sf_mu_id[1],2.) + pow(sf_mu_id[1],2.), 0.5) ;
  sf_ee_id[1]  =pow(pow(sf_e_id[1],2.) + pow(sf_e_id[1],2.), 0.5) ;
  sf_emu_id[1] =pow(pow(sf_mu_id[1],2.) + pow(sf_e_id[1],2.), 0.5) ;

  float sf_mumu_iso[2];
  float sf_ee_iso[2];
  float sf_emu_iso[2];

  sf_mumu_iso[0]=sf_mu_iso[0]*sf_mu_iso[0];
  sf_ee_iso[0]  =sf_e_iso[0]*sf_e_iso[0];
  sf_emu_iso[0] =sf_mu_iso[0]*sf_e_iso[0];
  sf_mumu_iso[1]=pow(pow(sf_mu_iso[1],2.) + pow(sf_mu_iso[1],2.), 0.5) ;
  sf_ee_iso[1]  =pow(pow(sf_e_iso[1],2.) + pow(sf_e_iso[1],2.), 0.5) ;
  sf_emu_iso[1] =pow(pow(sf_mu_iso[1],2.) + pow(sf_e_iso[1],2.), 0.5) ;

  float sf_mumu[2];
  float sf_ee[2];
  float sf_emu[2];

  sf_mumu[0]=sf_mumu_id[0]*sf_mumu_iso[0];
  sf_ee[0]  =sf_ee_id[0]*sf_ee_iso[0];
  sf_emu[0] =sf_emu_id[0]*sf_emu_iso[0];
  sf_mumu[1]=pow(pow(sf_mumu_id[1],2.) + pow(sf_mumu_iso[1],2.), 0.5) ;
  sf_ee[1]  =pow(pow(sf_ee_id[1],2.) + pow(sf_ee_iso[1],2.), 0.5) ;
  sf_emu[1] =pow(pow(sf_emu_id[1],2.) + pow(sf_emu_iso[1],2.), 0.5) ;

  cout << " ID " << endl;
  cout << " mumu " << sf_mumu_id[0] << " +/- " << sf_mumu_id[1] << endl;
  cout << " ee   " << sf_ee_id[0] << " +/- " << sf_ee_id[1] << endl;
  cout << " emu  " << sf_emu_id[0] << " +/- " << sf_emu_id[1] << endl;
  cout << endl;
  cout << " ISO " << endl;
  cout << " mumu " << sf_mumu_iso[0] << " +/- " << sf_mumu_iso[1] << endl;
  cout << " ee   " << sf_ee_iso[0] << " +/- " << sf_ee_iso[1] << endl;
  cout << " emu  " << sf_emu_iso[0] << " +/- " << sf_emu_iso[1] << endl;
  cout << endl;
  cout << " ID*IDO " << endl;
  cout << " mumu " << sf_mumu[0] << " +/- " << sf_mumu[1] << endl;
  cout << " ee   " << sf_ee[0] << " +/- " << sf_ee[1] << endl;
  cout << " emu  " << sf_emu[0] << " +/- " << sf_emu[1] << endl;
  cout << endl;
}

//*************************************
//Dump tex table with Tag&Probe results
//eff vs pt for 2 eta range,
//isol eff ad ID eff.
//*************************************

void dumpTagAndProbeTable(){
//plotTagAndProbEff(TString lepton, TString q, TString isoid, TString etaCut, ofstream ofile)
  string ofilenametex = string("TagAndProbe")+string(".tex");
  ofstream ofile(ofilenametex.c_str());
  ofile<<"\\documentclass[amsmath,amssymb]{revtex4}"<<endl;
  ofile<<"\\begin{document}"<<endl;
  
  ofile.setf(ios::fixed);
  ofile.precision(3);
  
  
  
  
   cout << "765  " << endl;
  //*********************************************************************************
  //nominal values
  std::vector<pair<double, double > > elecIsoEff_infEta15 = plotTagAndProbEff("El", "pt", "iso","_infEta15" , &ofile, 0, 0);
  std::vector<pair<double, double > > elecIsoEff_supEta15 = plotTagAndProbEff("El", "pt", "iso","_supEta15" , &ofile, 0, 0);
  std::vector<pair<double, double > > elecIDEff_infEta15  = plotTagAndProbEff("El", "pt", "id","_infEta15"  , &ofile, 0, 0);
  std::vector<pair<double, double > > elecIDEff_supEta15  = plotTagAndProbEff("El", "pt", "id","_supEta15"  , &ofile, 0, 0);
  std::vector<pair<double, double > > elecIsoSF_infEta15  = plotTagAndProbEff("El", "pt", "iso","_infEta15" , &ofile, 1, 0);
  std::vector<pair<double, double > > elecIsoSF_supEta15  = plotTagAndProbEff("El", "pt", "iso","_supEta15" , &ofile, 1, 0);
  std::vector<pair<double, double > > elecIDSF_infEta15   = plotTagAndProbEff("El", "pt", "id","_infEta15"  , &ofile, 1, 0);
  std::vector<pair<double, double > > elecIDSF_supEta15   = plotTagAndProbEff("El", "pt", "id","_supEta15"  , &ofile, 1, 0);
  
     cout << "765  " << endl;
  std::vector<pair<double, double > > muonIsoEff_infEta15 = plotTagAndProbEff("Mu", "pt", "iso","_infEta15" , &ofile, 0, 0);
  std::vector<pair<double, double > > muonIsoEff_supEta15 = plotTagAndProbEff("Mu", "pt", "iso","_supEta15" , &ofile, 0, 0);
  std::vector<pair<double, double > > muonIDEff_infEta15  = plotTagAndProbEff("Mu", "pt", "id","_infEta15"  , &ofile, 0, 0);
  std::vector<pair<double, double > > muonIDEff_supEta15  = plotTagAndProbEff("Mu", "pt", "id","_supEta15"  , &ofile, 0, 0);
  std::vector<pair<double, double > > muonIsoSF_infEta15  = plotTagAndProbEff("Mu", "pt", "iso","_infEta15" , &ofile, 1, 0);
  std::vector<pair<double, double > > muonIsoSF_supEta15  = plotTagAndProbEff("Mu", "pt", "iso","_supEta15" , &ofile, 1, 0);
  std::vector<pair<double, double > > muonIDSF_infEta15   = plotTagAndProbEff("Mu", "pt", "id","_infEta15"  , &ofile, 1, 0);
  std::vector<pair<double, double > > muonIDSF_supEta15   = plotTagAndProbEff("Mu", "pt", "id","_supEta15"  , &ofile, 1, 0);
  
  
     cout << "788  " << endl;
  
  //*********************************************************************************
  //syst up
  std::vector<pair<double, double > > elecIsoEff_infEta15_sUp = plotTagAndProbEff("El", "pt", "iso","_infEta15" , &ofile, 0, 1);
  std::vector<pair<double, double > > elecIsoEff_supEta15_sUp = plotTagAndProbEff("El", "pt", "iso","_supEta15" , &ofile, 0, 1);
  std::vector<pair<double, double > > elecIDEff_infEta15_sUp  = plotTagAndProbEff("El", "pt", "id","_infEta15"  , &ofile, 0, 1);
  std::vector<pair<double, double > > elecIDEff_supEta15_sUp  = plotTagAndProbEff("El", "pt", "id","_supEta15"  , &ofile, 0, 1);
  std::vector<pair<double, double > > elecIsoSF_infEta15_sUp  = plotTagAndProbEff("El", "pt", "iso","_infEta15" , &ofile, 1, 1);
  std::vector<pair<double, double > > elecIsoSF_supEta15_sUp  = plotTagAndProbEff("El", "pt", "iso","_supEta15" , &ofile, 1, 1);
  std::vector<pair<double, double > > elecIDSF_infEta15_sUp   = plotTagAndProbEff("El", "pt", "id","_infEta15"  , &ofile, 1, 1);
  std::vector<pair<double, double > > elecIDSF_supEta15_sUp   = plotTagAndProbEff("El", "pt", "id","_supEta15"  , &ofile, 1, 1);
     cout << "800  " << endl;
  std::vector<pair<double, double > > muonIsoEff_infEta15_sUp = plotTagAndProbEff("Mu", "pt", "iso","_infEta15" , &ofile, 0, 1);
  std::vector<pair<double, double > > muonIsoEff_supEta15_sUp = plotTagAndProbEff("Mu", "pt", "iso","_supEta15" , &ofile, 0, 0);
  std::vector<pair<double, double > > muonIDEff_infEta15_sUp  = plotTagAndProbEff("Mu", "pt", "id","_infEta15"  , &ofile, 0, 1);
  std::vector<pair<double, double > > muonIDEff_supEta15_sUp  = plotTagAndProbEff("Mu", "pt", "id","_supEta15"  , &ofile, 0, 1);
  std::vector<pair<double, double > > muonIsoSF_infEta15_sUp  = plotTagAndProbEff("Mu", "pt", "iso","_infEta15" , &ofile, 1, 1);
  std::vector<pair<double, double > > muonIsoSF_supEta15_sUp  = plotTagAndProbEff("Mu", "pt", "iso","_supEta15" , &ofile, 1, 1);
  std::vector<pair<double, double > > muonIDSF_infEta15_sUp   = plotTagAndProbEff("Mu", "pt", "id","_infEta15"  , &ofile, 1, 1);
  std::vector<pair<double, double > > muonIDSF_supEta15_sUp   = plotTagAndProbEff("Mu", "pt", "id","_supEta15"  , &ofile, 1, 1);
     cout << "809 " << endl;
  
  //*********************************************************************************
  //syst down
  std::vector<pair<double, double > > elecIsoEff_infEta15_sDown = plotTagAndProbEff("El", "pt", "iso","_infEta15" , &ofile, 0, -1);
  std::vector<pair<double, double > > elecIsoEff_supEta15_sDown = plotTagAndProbEff("El", "pt", "iso","_supEta15" , &ofile, 0, -1);
  std::vector<pair<double, double > > elecIDEff_infEta15_sDown  = plotTagAndProbEff("El", "pt", "id","_infEta15"  , &ofile, 0, -1);
  std::vector<pair<double, double > > elecIDEff_supEta15_sDown  = plotTagAndProbEff("El", "pt", "id","_supEta15"  , &ofile, 0, -1);
  std::vector<pair<double, double > > elecIsoSF_infEta15_sDown  = plotTagAndProbEff("El", "pt", "iso","_infEta15" , &ofile, 1, -1);
  std::vector<pair<double, double > > elecIsoSF_supEta15_sDown  = plotTagAndProbEff("El", "pt", "iso","_supEta15" , &ofile, 1, -1);
  std::vector<pair<double, double > > elecIDSF_infEta15_sDown   = plotTagAndProbEff("El", "pt", "id","_infEta15"  , &ofile, 1, -1);
  std::vector<pair<double, double > > elecIDSF_supEta15_sDown   = plotTagAndProbEff("El", "pt", "id","_supEta15"  , &ofile, 1, -1);
       cout << "765  " << endl;
  std::vector<pair<double, double > > muonIsoEff_infEta15_sDown = plotTagAndProbEff("Mu", "pt", "iso","_infEta15" , &ofile, 0, -1);
  std::vector<pair<double, double > > muonIsoEff_supEta15_sDown = plotTagAndProbEff("Mu", "pt", "iso","_supEta15" , &ofile, 0, -1);
  std::vector<pair<double, double > > muonIDEff_infEta15_sDown  = plotTagAndProbEff("Mu", "pt", "id","_infEta15"  , &ofile, 0, -1);
  std::vector<pair<double, double > > muonIDEff_supEta15_sDown  = plotTagAndProbEff("Mu", "pt", "id","_supEta15"  , &ofile, 0, -1);
  std::vector<pair<double, double > > muonIsoSF_infEta15_sDown  = plotTagAndProbEff("Mu", "pt", "iso","_infEta15" , &ofile, 1, -1);
  std::vector<pair<double, double > > muonIsoSF_supEta15_sDown  = plotTagAndProbEff("Mu", "pt", "iso","_supEta15" , &ofile, 1, -1);
  std::vector<pair<double, double > > muonIDSF_infEta15_sDown   = plotTagAndProbEff("Mu", "pt", "id","_infEta15"  , &ofile, 1, -1);
  std::vector<pair<double, double > > muonIDSF_supEta15_sDown   = plotTagAndProbEff("Mu", "pt", "id","_supEta15"  , &ofile, 1, -1);
  
  
  
     double ptRange[] = {20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 100};
  
  //***************
  //for electrons 
  //***************
  // int systType
  ofile << "\\begin{table}[p]" << endl;
  ofile << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 20-22 & 22-24 & 24-26 & 26-28 & 28-30  \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIsoEff_infEta15, elecIsoEff_infEta15_sUp, elecIsoEff_infEta15_sDown, 1);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIsoEff_supEta15, elecIsoEff_supEta15_sUp, elecIsoEff_supEta15_sDown, 1);
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 30-32 & 32-34 & 34-36 & 36-40 & 40-100 \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIsoEff_infEta15, elecIsoEff_infEta15_sUp, elecIsoEff_infEta15_sDown, 0);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIsoEff_supEta15, elecIsoEff_supEta15_sUp, elecIsoEff_supEta15_sDown, 0);
  ofile << "\\hline" << endl;
  ofile << "\\end{tabular}" << endl;
  ofile << "\\caption{Electron isolation efficiency as function of electrons $p_T$ [GeV/c] and $\\eta$.}" << endl;
  ofile << "\\end{table}" << endl;

  ofile << endl;
  ofile << endl;
  ofile << endl;
  
  ofile << "\\begin{table}[p]" << endl;
  ofile << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 20-22 & 22-24 & 24-26 & 26-28 & 28-30  \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIDEff_infEta15, elecIDEff_infEta15_sUp, elecIDEff_infEta15_sDown, 1);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIDEff_supEta15, elecIDEff_supEta15_sUp, elecIDEff_supEta15_sDown, 1);
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 30-32 & 32-34 & 34-36 & 36-40 & 40-100 \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIDEff_infEta15, elecIDEff_infEta15_sUp, elecIDEff_infEta15_sDown, 0);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIDEff_supEta15, elecIDEff_supEta15_sUp, elecIDEff_supEta15_sDown, 0);
  ofile << "\\hline" << endl;
  ofile << "\\end{tabular}" << endl;
  ofile << "\\caption{Electron identification efficiency as function of electrons $p_T$ [GeV/c] and $\\eta$.}" << endl;
  ofile << "\\end{table}" << endl;
 



ofile << "\\begin{table}[p]" << endl;
  ofile << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 20-22 & 22-24 & 24-26 & 26-28 & 28-30  \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIsoEff_infEta15, muonIsoEff_infEta15_sUp, muonIsoEff_infEta15_sDown, 1);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIsoEff_supEta15, muonIsoEff_supEta15_sUp, muonIsoEff_supEta15_sDown, 1);
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 30-32 & 32-34 & 34-36 & 36-40 & 40-100 \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIsoEff_infEta15, muonIsoEff_infEta15_sUp, muonIsoEff_infEta15_sDown, 0);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIsoEff_supEta15, muonIsoEff_supEta15_sUp, muonIsoEff_supEta15_sDown, 0);
  ofile << "\\hline" << endl;
  ofile << "\\end{tabular}" << endl;
  ofile << "\\caption{muon isolation efficiency as function of muons $p_T$ [GeV/c] and $\\eta$.}" << endl;
  ofile << "\\end{table}" << endl;

  ofile << endl;
  ofile << endl;
  ofile << endl;
  
  ofile << "\\begin{table}[p]" << endl;
  ofile << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 20-22 & 22-24 & 24-26 & 26-28 & 28-30  \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIDEff_infEta15, muonIDEff_infEta15_sUp, muonIDEff_infEta15_sDown, 1);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIDEff_supEta15, muonIDEff_supEta15_sUp, muonIDEff_supEta15_sDown, 1);
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 30-32 & 32-34 & 34-36 & 36-40 & 40-100 \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIDEff_infEta15, muonIDEff_infEta15_sUp, muonIDEff_infEta15_sDown, 0);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIDEff_supEta15, muonIDEff_supEta15_sUp, muonIDEff_supEta15_sDown, 0);
  ofile << "\\hline" << endl;
  ofile << "\\end{tabular}" << endl;
  ofile << "\\caption{muon identification efficiency as function of muons $p_T$ [GeV/c] and $\\eta$.}" << endl;
  ofile << "\\end{table}" << endl;
 

  ofile << endl;
  ofile << endl;
  ofile << endl;
   ofile << "\\begin{table}[p]" << endl;
  ofile << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 20-22 & 22-24 & 24-26 & 26-28 & 28-30  \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIsoSF_infEta15, elecIsoSF_infEta15_sUp, elecIsoSF_infEta15_sDown, 1);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIsoSF_supEta15, elecIsoSF_supEta15_sUp, elecIsoSF_supEta15_sDown, 1);
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 30-32 & 32-34 & 34-36 & 36-40 & 40-100 \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIsoSF_infEta15, elecIsoSF_infEta15_sUp, elecIsoSF_infEta15_sDown, 0);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIsoSF_supEta15, elecIsoSF_supEta15_sUp, elecIsoSF_supEta15_sDown, 0);
  ofile << "\\hline" << endl;
  ofile << "\\end{tabular}" << endl;
  ofile << "\\caption{Electron isolation SF as function of electrons $p_T$ [GeV/c] and $\\eta$.}" << endl;
  ofile << "\\end{table}" << endl;

  ofile << endl;
  ofile << endl;
  ofile << endl;
  
  ofile << "\\begin{table}[p]" << endl;
  ofile << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 20-22 & 22-24 & 24-26 & 26-28 & 28-30  \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIDSF_infEta15, elecIDSF_infEta15_sUp, elecIDSF_infEta15_sDown, 1);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIDSF_supEta15, elecIDSF_supEta15_sUp, elecIDSF_supEta15_sDown, 1);
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 30-32 & 32-34 & 34-36 & 36-40 & 40-100 \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIDSF_infEta15, elecIDSF_infEta15_sUp, elecIDSF_infEta15_sDown, 0);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, elecIDSF_supEta15, elecIDSF_supEta15_sUp, elecIDSF_supEta15_sDown, 0);
  ofile << "\\hline" << endl;
  ofile << "\\end{tabular}" << endl;
  ofile << "\\caption{Electron identification SF as function of electrons $p_T$ [GeV/c] and $\\eta$.}" << endl;
  ofile << "\\end{table}" << endl;


  ofile << endl;
  ofile << endl;
  ofile << endl;


  ofile << endl;
  ofile << endl;
  ofile << endl;
   ofile << "\\begin{table}[p]" << endl;
  ofile << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 20-22 & 22-24 & 24-26 & 26-28 & 28-30  \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIsoSF_infEta15, muonIsoSF_infEta15_sUp, muonIsoSF_infEta15_sDown, 1);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIsoSF_supEta15, muonIsoSF_supEta15_sUp, muonIsoSF_supEta15_sDown, 1);
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 30-32 & 32-34 & 34-36 & 36-40 & 40-100 \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIsoSF_infEta15, muonIsoSF_infEta15_sUp, muonIsoSF_infEta15_sDown, 0);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIsoSF_supEta15, muonIsoSF_supEta15_sUp, muonIsoSF_supEta15_sDown, 0);
  ofile << "\\hline" << endl;
  ofile << "\\end{tabular}" << endl;
  ofile << "\\caption{muon isolation SF as function of muons $p_T$ [GeV/c] and $\\eta$.}" << endl;
  ofile << "\\end{table}" << endl;

  ofile << endl;
  ofile << endl;
  ofile << endl;
  
  ofile << "\\begin{table}[p]" << endl;
  ofile << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 20-22 & 22-24 & 24-26 & 26-28 & 28-30  \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIDSF_infEta15, muonIDSF_infEta15_sUp, muonIDSF_infEta15_sDown, 1);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIDSF_supEta15, muonIDSF_supEta15_sUp, muonIDSF_supEta15_sDown, 1);
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 30-32 & 32-34 & 34-36 & 36-40 & 40-100 \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIDSF_infEta15, muonIDSF_infEta15_sUp, muonIDSF_infEta15_sDown, 0);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, muonIDSF_supEta15, muonIDSF_supEta15_sUp, muonIDSF_supEta15_sDown, 0);
  ofile << "\\hline" << endl;
  ofile << "\\end{tabular}" << endl;
  ofile << "\\caption{muon identification SF as function of muons $p_T$ [GeV/c] and $\\eta$.}" << endl;
  ofile << "\\end{table}" << endl;


  cout << "compute global scale factor" << endl;
  //***************
  //Global SF
  //***************
  std::vector<pair<double, double > > mergeSFMuons_infEta15          = multiplySF(muonIDSF_infEta15,       muonIsoSF_infEta15);
  cout << "1102 " << mergeSFMuons_infEta15.size() << endl;
  std::vector<pair<double, double > > mergeSFMuonsSystUp_infEta15    = multiplySF(muonIDSF_infEta15_sUp,   muonIsoSF_infEta15_sUp);
  cout << "1104 " << mergeSFMuonsSystUp_infEta15.size() << endl;
  std::vector<pair<double, double > > mergeSFMuonsSystDown_infEta15  = multiplySF(muonIDSF_infEta15_sDown, muonIsoSF_infEta15_sDown);
  
  cout << "1107 " << mergeSFMuonsSystDown_infEta15.size()  << endl;
  std::vector<pair<double, double > > mergeSFMuons_supEta15          = multiplySF(muonIDSF_supEta15,       muonIsoSF_supEta15);
  cout << "1109 " << mergeSFMuons_supEta15.size()  << endl;
  std::vector<pair<double, double > > mergeSFMuonsSystUp_supEta15    = multiplySF(muonIDSF_supEta15_sUp,   muonIsoSF_supEta15_sUp);
  cout << "1111 " << mergeSFMuonsSystUp_supEta15.size()  << endl;
  std::vector<pair<double, double > > mergeSFMuonsSystDown_supEta15  = multiplySF(muonIDSF_supEta15_sDown, muonIsoSF_supEta15_sDown);
  cout << "1113 " << mergeSFMuonsSystDown_supEta15.size()  << endl;

  ofile << "\\begin{table}[p]" << endl;
  ofile << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 20-22 & 22-24 & 24-26 & 26-28 & 28-30  \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, mergeSFMuons_infEta15, mergeSFMuonsSystUp_infEta15, mergeSFMuonsSystDown_infEta15, 1);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, mergeSFMuons_supEta15, mergeSFMuonsSystUp_supEta15, mergeSFMuonsSystDown_supEta15, 1);
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 30-32 & 32-34 & 34-36 & 36-40 & 40-100 \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, mergeSFMuons_infEta15, mergeSFMuonsSystUp_infEta15, mergeSFMuonsSystDown_infEta15, 0);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, mergeSFMuons_supEta15, mergeSFMuonsSystUp_supEta15, mergeSFMuonsSystDown_supEta15, 0);
  ofile << "\\hline" << endl;
  ofile << "\\end{tabular}" << endl;
  ofile << "\\caption{muon overall SF as function of muon $p_T$ [GeV/c] and $\\eta$.}" << endl;
  ofile << "\\end{table}" << endl;

  std::vector<pair<double, double > > mergeSFelectrons_infEta15          = multiplySF(elecIDSF_infEta15,       elecIsoSF_infEta15);
  std::vector<pair<double, double > > mergeSFelectronsSystUp_infEta15    = multiplySF(elecIDSF_infEta15_sUp,   elecIsoSF_infEta15_sUp);
  std::vector<pair<double, double > > mergeSFelectronsSystDown_infEta15  = multiplySF(elecIDSF_infEta15_sDown, elecIsoSF_infEta15_sDown);
  
  std::vector<pair<double, double > > mergeSFelectrons_supEta15          = multiplySF(elecIDSF_supEta15,       elecIsoSF_supEta15);
  std::vector<pair<double, double > > mergeSFelectronsSystUp_supEta15    = multiplySF(elecIDSF_supEta15_sUp,   elecIsoSF_supEta15_sUp);
  std::vector<pair<double, double > > mergeSFelectronsSystDown_supEta15  = multiplySF(elecIDSF_supEta15_sDown, elecIsoSF_supEta15_sDown);

  ofile << "\\begin{table}[p]" << endl;
  ofile << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 20-22 & 22-24 & 24-26 & 26-28 & 28-30  \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, mergeSFelectrons_infEta15, mergeSFelectronsSystUp_infEta15, mergeSFelectronsSystDown_infEta15, 1);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, mergeSFelectrons_supEta15, mergeSFelectronsSystUp_supEta15, mergeSFelectronsSystDown_supEta15, 1);
  ofile << "\\hline" << endl;
  ofile << " pt ranges & 30-32 & 32-34 & 34-36 & 36-40 & 40-100 \\\\" << endl;
  ofile << "\\hline" << endl;
  ofile << "$|\\eta| <$ 1.5 &" ; dumpTagAndProbeValue(&ofile, mergeSFelectrons_infEta15, mergeSFelectronsSystUp_infEta15, mergeSFelectronsSystDown_infEta15, 0);
  ofile << "$|\\eta| >$ 1.5 &" ; dumpTagAndProbeValue(&ofile, mergeSFelectrons_supEta15, mergeSFelectronsSystUp_supEta15, mergeSFelectronsSystDown_supEta15, 0);
  ofile << "\\hline" << endl;
  ofile << "\\end{tabular}" << endl;
  ofile << "\\caption{electron overall SF as function of electron $p_T$ [GeV/c] and $\\eta$.}" << endl;
  ofile << "\\end{table}" << endl;

  dumpTagAndProbeValueRoot(
  	mergeSFelectrons_infEta15, mergeSFelectronsSystUp_infEta15, mergeSFelectronsSystDown_infEta15,
  	mergeSFelectrons_supEta15, mergeSFelectronsSystUp_supEta15, mergeSFelectronsSystDown_supEta15, 0);
  
  
  dumpTagAndProbeValueRoot(
  	mergeSFMuons_infEta15, mergeSFMuonsSystUp_infEta15, mergeSFMuonsSystDown_infEta15,
  	mergeSFMuons_supEta15, mergeSFMuonsSystUp_supEta15, mergeSFMuonsSystDown_supEta15, 1);
  
  
  
  ofile<<"\\end{document}"<<endl;

}






//************************************************
//compare isolation efficiencies in MC ttbar vs Z.
//extract a correction factor out of it
//************************************************


void compareZ_TTbaTagAndProbef(TString lepton, float shift){
   
  Int_t stati=0;
  Bool_t  fit=1;
  Bool_t logy=0;
  
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  gStyle->SetOptDate(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleColor(1);
  //gStyle->SetTitleW(0.4);
  //gStyle->SetTitleH(0.07);
  //gStyle->SetOptStat(stati);
  gStyle->SetOptStat(0);
  // gStyle->SetFrameFillColor(18);
  //gStyle->SetPalette(51,0);
  gStyle->SetPalette(1);
  
  
  TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
   
  c1->SetFillColor(10);
  c1->SetFillStyle(4000);
  c1->SetBorderSize(2);
  
  //c1.Divide(2,2);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(10);
  
   
   
   //TFile *f_MC_Zjets = new TFile("RootFileTagAndProbeta15_rebin/TriggerPlots.root"); 
   //TFile *f_MC_TTbar = new TFile("RootFileTagAndProbeta15_rebin/TriggerPlots.root");  
 
   //TFile *f_MC_Zjets = new TFile("TriggerPlots_eff2Jets.root"); 
   //TFile *f_MC_TTbar = new TFile("TriggerPlots_eff2Jets.root");  
   
   TFile *f_MC_Zjets = new TFile("TriggerPlots_Zjets.root"); 
   TFile *f_MC_TTbar = new TFile("TriggerPlots_eff2Jets.root");  
   
   TH1F* tight_histo_TTbar    ;
   TH1F* looseID_histo_TTbar  ;
   TH1F* looseIso_histo_TTbar ;
   
   
   TH1F* tight_histo_Zjets    ;
   TH1F* looseID_histo_Zjets  ;
   TH1F* looseIso_histo_Zjets ;
   
   if(lepton == "Mu"){
   
           
      f_MC_Zjets->cd();
      tight_histo_Zjets    = (TH1F*)gROOT->FindObject("genMu_Tight_pt_Zjets");
      looseID_histo_Zjets  = (TH1F*)gROOT->FindObject("genMu_LooseID_pt_Zjets");
      looseIso_histo_Zjets = (TH1F*)gROOT->FindObject("genMu_LooseIso_pt_Zjets");
      
      f_MC_TTbar->cd();
      tight_histo_TTbar    = (TH1F*)gROOT->FindObject("genMu_Tight_pt_TTbarSig");
      looseID_histo_TTbar  = (TH1F*)gROOT->FindObject("genMu_LooseID_pt_TTbarSig");
      looseIso_histo_TTbar = (TH1F*)gROOT->FindObject("genMu_LooseIso_pt_TTbarSig");
   
   }
   
   if(lepton == "El"){
   
       f_MC_Zjets->cd();
      tight_histo_Zjets    = (TH1F*)gROOT->FindObject("genEl_Tight_pt_Zjets");
      looseID_histo_Zjets  = (TH1F*)gROOT->FindObject("genEl_LooseID_pt_Zjets");
      looseIso_histo_Zjets = (TH1F*)gROOT->FindObject("genEl_LooseIso_pt_Zjets");
      
      f_MC_TTbar->cd();
      tight_histo_TTbar    = (TH1F*)gROOT->FindObject("genEl_Tight_pt_TTbarSig");
      looseID_histo_TTbar  = (TH1F*)gROOT->FindObject("genEl_LooseID_pt_TTbarSig");
      looseIso_histo_TTbar = (TH1F*)gROOT->FindObject("genEl_LooseIso_pt_TTbarSig");
   
   }
   
   
    
   
   TGraphAsymmErrors * grID_Zjets  = new TGraphAsymmErrors();
   grID_Zjets->BayesDivide(tight_histo_Zjets,looseID_histo_Zjets);
   TGraphAsymmErrors * grIso_Zjets = new TGraphAsymmErrors();
   grIso_Zjets->BayesDivide(tight_histo_Zjets,looseIso_histo_Zjets);
     
   
     
   TGraphAsymmErrors * grID_TTbar  = new TGraphAsymmErrors();
   grID_TTbar->BayesDivide(tight_histo_TTbar,looseID_histo_TTbar);
   TGraphAsymmErrors * grIso_TTbar = new TGraphAsymmErrors();
   grIso_TTbar->BayesDivide(tight_histo_TTbar,looseIso_histo_TTbar);
     
   //float ptRange[] = {20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 65, 80, 100};
  
   float ptRange[] = {20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 100};
   TH1F * thehistSF = new TH1F("thehistSF", "thehistSF", 10, ptRange);
  
   thehistSF->SetMaximum(1.1);
   thehistSF->SetMinimum(0.);
   if(lepton == "El") thehistSF->GetXaxis()->SetTitle("electron p_{T} [GeV/c]");
   if(lepton == "Mu") thehistSF->GetXaxis()->SetTitle("muon p_{T} [GeV/c]");
   
   grID_Zjets->SetMarkerColor(2);
   grID_TTbar->SetMarkerColor(4);
   grID_Zjets->SetMarkerStyle(20);
   grID_TTbar->SetMarkerStyle(20);
   
   thehistSF->Draw();
   grID_Zjets->Draw("epsame");
   grID_TTbar->Draw("epsame");
   
   c1->SaveAs(("comp_tt_Z_ID_"+lepton+"_.gif").Data());

   
   
   grIso_Zjets->SetMarkerColor(2);
   grIso_TTbar->SetMarkerColor(4);
   grIso_Zjets->SetMarkerStyle(20);
   grIso_TTbar->SetMarkerStyle(20);
   
   thehistSF->Draw();
   grIso_Zjets->Draw("epsame");
   grIso_TTbar->Draw("epsame");
   
   TFile * calibLeptModelling = new TFile("calibLeptModel.root", "recreate");
   TH1F *histocalib = new TH1F("histocalib", "histocalib", 10, ptRange);
   
   for(int i=0; i<14; i++){
        double x_data = 0;
        double eff_ID_Zjets= 0;
        double eff_ID_TTbar= 0;
	grID_Zjets->GetPoint(i, x_data, eff_ID_Zjets);
	grID_TTbar->GetPoint(i, x_data, eff_ID_TTbar);
	cout << "ID diff  " << eff_ID_Zjets - eff_ID_TTbar << endl;
   }
   
   
   double maxvalDiffIso = 0;
   double eff_Iso_Zjets[14];
   double eff_Iso_TTbar[14];
   for(int i=0; i<14; i++){
        double x_data = 0;
        double eff_Iso_Zjets= 0;
        double eff_Iso_TTbar= 0;
	grIso_Zjets->GetPoint(i, x_data, eff_Iso_Zjets);
	grIso_TTbar->GetPoint(i, x_data, eff_Iso_TTbar);
	cout << "Iso diff  " << eff_Iso_Zjets - eff_Iso_TTbar << endl;
	if(fabs(eff_Iso_Zjets - eff_Iso_TTbar) > maxvalDiffIso) maxvalDiffIso = fabs(eff_Iso_Zjets - eff_Iso_TTbar);
	histocalib->SetBinContent(  shift*eff_Iso_TTbar/eff_Iso_Zjets, i+1);
   }
   if(lepton == "El") cout << "maxIso El " << maxvalDiffIso << endl;
   if(lepton == "Mu") cout << "maxIso Mu " << maxvalDiffIso << endl;
   c1->SaveAs(("comp_tt_Z_Iso_"+lepton+"_.gif").Data());
   
   

}



//********************************
//Draw correlation plots for
//trigger efficiency estimate
//********************************


void drawCorrelationPlot(TString channel){
Int_t stati=0;
  Bool_t  fit=1;
  Bool_t logy=0;
  
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  gStyle->SetOptDate(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleColor(1);
  //gStyle->SetTitleW(0.4);
  //gStyle->SetTitleH(0.07);
  //gStyle->SetOptStat(stati);
  gStyle->SetOptStat(0);
  // gStyle->SetFrameFillColor(18);
  //gStyle->SetPalette(51,0);
  gStyle->SetPalette(1);
  
  
  TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
   
  c1->SetFillColor(10);
  c1->SetFillStyle(4000);
  c1->SetBorderSize(2);
  //c1->SetLogz(1);
  
  //c1.Divide(2,2);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(10);
  
  TFile *f_trigger = new TFile("TagAndProbe/TriggerPlots.root"); 
 
  //TH2D * histo_ee = (TH2D*)gROOT->FindObject("correlationMET_ee_TTbarSig");
  //TH2D * histo_em = (TH2D*)gROOT->FindObject("correlationMET_emu_TTbarSig");
  //TH2D * histo_mm = (TH2D*)gROOT->FindObject("correlationMET_mumu_TTbarSig");
  
  TH2D * histo_ee = (TH2D*)gROOT->FindObject("correlationMET_ee_Zjets");
  TH2D * histo_em = (TH2D*)gROOT->FindObject("correlationMET_emu_Zjets");
  TH2D * histo_mm = (TH2D*)gROOT->FindObject("correlationMET_mumu_Zjets");
  
  
  
  if(channel == "ee"){
    histo_ee->GetXaxis()->SetBinLabel(1, "fail MET");
    histo_ee->GetXaxis()->SetBinLabel(2, "pass MET");
    histo_ee->GetYaxis()->SetBinLabel(1, "fail Dilept");
    histo_ee->GetYaxis()->SetBinLabel(2, "pass Dilept");
    
    histo_ee->GetXaxis()->SetTitle( "");
    histo_ee->GetYaxis()->SetTitle( "");
    histo_ee->SetTitle("correlations ee");
    
    histo_ee->Draw("colz");
    
    double n_total  = histo_ee->Integral();
    double N_dilept = histo_ee->GetBinContent(1, 2) + histo_ee->GetBinContent(2, 2);
    double N_met    = histo_ee->GetBinContent(2, 1) + histo_ee->GetBinContent(2, 2);
    double N_dilept_met   = histo_ee->GetBinContent(2, 2) ;
    
    double eff_dilept = N_dilept/n_total;
    cout << "eff dilept " << eff_dilept << endl;
    double eff_met = N_met/n_total;
    cout << "eff met " << eff_met << endl;
    double eff_dilept_met = N_dilept_met/n_total;
    cout << "eff dilept meet " <<  eff_dilept_met<< endl;
    cout << "eff_dilept*eff_met " << eff_met*eff_dilept << endl;
    cout << "corr factor is " << eff_met*eff_dilept/eff_dilept_met<< endl;
    
    
    c1->SaveAs("corree.gif");
  }
  
  
  
  if(channel == "emu"){
    histo_em->GetXaxis()->SetBinLabel(1, "fail MET");
    histo_em->GetXaxis()->SetBinLabel(2, "pass MET");
    histo_em->GetYaxis()->SetBinLabel(1, "fail Dilept");
    histo_em->GetYaxis()->SetBinLabel(2, "pass Dilept");
    
    histo_em->GetXaxis()->SetTitle( "");
    histo_em->GetYaxis()->SetTitle( "");
    histo_em->SetTitle("correlations emu");
    
    histo_em->Draw("colz");
    
    double n_total  = histo_em->Integral();
    double N_dilept = histo_em->GetBinContent(1, 2) + histo_em->GetBinContent(2, 2);
    double N_met    = histo_em->GetBinContent(2, 1) + histo_em->GetBinContent(2, 2);
    double N_dilept_met   = histo_em->GetBinContent(2, 2) ;
    
    double eff_dilept = N_dilept/n_total;
    cout << "eff dilept " << eff_dilept << endl;
    double eff_met = N_met/n_total;
    cout << "eff met " << eff_met << endl;
    double eff_dilept_met = N_dilept_met/n_total;
    cout << "eff dilept meet " <<  eff_dilept_met<< endl;
    cout << "eff_dilept*eff_met " << eff_met*eff_dilept << endl;
    cout << "corr factor is " << eff_met*eff_dilept/eff_dilept_met<< endl;
    c1->SaveAs("corremu.gif");
  
  }
  
  
  
  if(channel == "mumu"){
    histo_mm->GetXaxis()->SetBinLabel(1, "fail MET");
    histo_mm->GetXaxis()->SetBinLabel(2, "pass MET");
    histo_mm->GetYaxis()->SetBinLabel(1, "fail Dilept");
    histo_mm->GetYaxis()->SetBinLabel(2, "pass Dilept");
    
    histo_mm->GetXaxis()->SetTitle( "");
    histo_mm->GetYaxis()->SetTitle( "");
    histo_mm->SetTitle("correlations mumu");
    
    histo_mm->Draw("colz");  
    double n_total  = histo_mm->Integral();
    double N_dilept = histo_mm->GetBinContent(1, 2) + histo_mm->GetBinContent(2, 2);
    double N_met    = histo_mm->GetBinContent(2, 1) + histo_mm->GetBinContent(2, 2);
    double N_dilept_met   = histo_mm->GetBinContent(2, 2) ;
    
    double eff_dilept = N_dilept/n_total;
    cout << "eff dilept " << eff_dilept << endl;
    double eff_met = N_met/n_total;
    cout << "eff met " << eff_met << endl;
    double eff_dilept_met = N_dilept_met/n_total;
    cout << "eff dilept meet " <<  eff_dilept_met<< endl;
    cout << "eff_dilept*eff_met " << eff_met*eff_dilept << endl;
    cout << "corr factor is " << eff_met*eff_dilept/eff_dilept_met<< endl;
   
    c1->SaveAs("corrmumu.gif");
  
  }
     

}


void plotDataTagAndProbEff(TString lepton, TString q, TString isoid, TString etaCut){

  
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  gStyle->SetOptDate(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleColor(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  
  TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
   
  c1->SetFillColor(10);
  c1->SetFillStyle(4000);
  c1->SetBorderSize(2);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(10);
  
  
  TFile *f_Data_1 ; 
  TFile *f_Data_2 ; 
  TFile *f_Data_3 ; 
  TFile *f_Data_4 ; 
  
  string datasetname;
  if(lepton == "Mu") datasetname="DataMu";
  else datasetname="DataEG";
  TString filenam1="TagAndProbe_par_periode/tagAndProbeHisto_"+datasetname+"_may.root";
  TString filenam2="TagAndProbe_par_periode/tagAndProbeHisto_"+datasetname+"_prompv4.root";
  TString filenam3="TagAndProbe_par_periode/tagAndProbeHisto_"+datasetname+"_aug.root";
  TString filenam4="TagAndProbe_par_periode/tagAndProbeHisto_"+datasetname+"_prompv6.root";
  f_Data_1 = new TFile(filenam1);
  f_Data_2 = new TFile(filenam2);
  f_Data_3 = new TFile(filenam3);
  f_Data_4 = new TFile(filenam4);
  
  
  TH1F* tight_histo1;
  TH1F* loose_histo1;
  TH1F* tight_histo2;
  TH1F* loose_histo2;
  TH1F* tight_histo3;
  TH1F* loose_histo3;
  TH1F* tight_histo4;
  TH1F* loose_histo4;
  
   
  TLegend* qw = 0;
  qw = new TLegend(0.70,0.70,0.99,0.99);
 
  TGraphAsymmErrors * gTagAndProbef1;
  TGraphAsymmErrors * gTagAndProbef2;
  TGraphAsymmErrors * gTagAndProbef3;
  TGraphAsymmErrors * gTagAndProbef4;

  string histoname1;
  if (lepton=="Mu") histoname1="muon";
  else histoname1="electron";
  string histoname2;
  if (isoid=="id") histoname2="ID";
  else if (isoid=="iso") histoname2="Iso";
    
  
  f_Data_1->cd();
  tight_histo1    = (TH1F*)gROOT->FindObject( (datasetname+histoname1+"s_Tight_"+q+etaCut).Data() );
  loose_histo1  = (TH1F*)gROOT->FindObject( (datasetname+histoname1+"s_loose"+histoname2+"_"+q+etaCut).Data() );
  gTagAndProbef1  = new TGraphAsymmErrors();
  gTagAndProbef1->BayesDivide(tight_histo1,loose_histo1);
      
  f_Data_2->cd();
  tight_histo2    = (TH1F*)gROOT->FindObject( (datasetname+histoname1+"s_Tight_"+q+etaCut).Data() );
  loose_histo2  = (TH1F*)gROOT->FindObject( (datasetname+histoname1+"s_loose"+histoname2+"_"+q+etaCut).Data() );
  gTagAndProbef2  = new TGraphAsymmErrors();
  gTagAndProbef2->BayesDivide(tight_histo2,loose_histo2);

  f_Data_3->cd();
  tight_histo3    = (TH1F*)gROOT->FindObject( (datasetname+histoname1+"s_Tight_"+q+etaCut).Data() );
  loose_histo3  = (TH1F*)gROOT->FindObject( (datasetname+histoname1+"s_loose"+histoname2+"_"+q+etaCut).Data() );
  gTagAndProbef3  = new TGraphAsymmErrors();
  gTagAndProbef3->BayesDivide(tight_histo3,loose_histo3);

  f_Data_4->cd();
  tight_histo4    = (TH1F*)gROOT->FindObject( (datasetname+histoname1+"s_Tight_"+q+etaCut).Data() );
  loose_histo4  = (TH1F*)gROOT->FindObject( (datasetname+histoname1+"s_loose"+histoname2+"_"+q+etaCut).Data() );
  gTagAndProbef4  = new TGraphAsymmErrors();
  gTagAndProbef4->BayesDivide(tight_histo4,loose_histo4);
     
  TString histo_labelx;
  if (q=="pt") { histo_labelx= histoname1+" p_{T} [GeV/c]" ; }  
  else if (q=="eta") { histo_labelx= histoname1+" #eta" ; } 
  else if(q == "phi"){ histo_labelx= histoname1+" #phi" ; }
  else if (q=="njet") { histo_labelx="NJet" ; }
  loose_histo1->GetXaxis()->SetTitle(histo_labelx);
      
  TString gr_label1=histoname1+" "+histoname2+" eff. (May ReReco)";
  TString gr_label2=histoname1+" "+histoname2+" eff. (Prompt v4)";
  TString gr_label3=histoname1+" "+histoname2+" eff. (Aug ReReco)";
  TString gr_label4=histoname1+" "+histoname2+" eff. (Prompt v6)";
  qw->AddEntry(gTagAndProbef1,  gr_label1,"P");
  qw->AddEntry(gTagAndProbef2,  gr_label2,"P");
  qw->AddEntry(gTagAndProbef3,  gr_label3,"P");
  qw->AddEntry(gTagAndProbef4,  gr_label4,"P");
     
  tight_histo1->Sumw2();
  loose_histo1->Sumw2();
  loose_histo1->Divide( tight_histo1, loose_histo1, 1, 1, "b");
  loose_histo1->SetMaximum(1.1);
  float min=loose_histo1->GetMinimum();
  loose_histo1->SetMinimum(min-0.1);

  loose_histo1->SetMarkerStyle(20);
  loose_histo1->SetMarkerColor(2);
  loose_histo1->Draw("p");
  
  gTagAndProbef1->SetMarkerStyle(20);
  gTagAndProbef1->SetMarkerColor(2);
  gTagAndProbef2->SetMarkerStyle(21);
  gTagAndProbef2->SetMarkerColor(4);
  gTagAndProbef3->SetMarkerStyle(26);
  gTagAndProbef3->SetMarkerColor(8);
  gTagAndProbef4->SetMarkerStyle(27);
  gTagAndProbef4->SetMarkerColor(1);
  
  
  gTagAndProbef1->Draw("psame");
  gTagAndProbef2->Draw("psame");
  gTagAndProbef3->Draw("psame");
  gTagAndProbef4->Draw("psame");
   
  
  qw->SetFillColor(0);
  qw->Draw();
  TString outputfilename =  histoname1+"_"+ q+"_"+isoid+"_4ranges.gif";
  
  c1->SaveAs(outputfilename.Data());
  
 
  tight_histo1->Delete();
  loose_histo1->Delete();
  tight_histo2->Delete();
  loose_histo2->Delete();
  tight_histo3->Delete();
  loose_histo3->Delete();
  tight_histo4->Delete();
  loose_histo4->Delete();
  gTagAndProbef1->Delete();
  gTagAndProbef2->Delete();
  gTagAndProbef3->Delete();
  gTagAndProbef4->Delete();
  f_Data_1->Close();
  f_Data_2->Close();
  f_Data_3->Close();
  f_Data_4->Close();

  
  
}

void ProducePlots_par_periode(){
  
  plotDataTagAndProbEff("El", "pt", "iso","");
  plotDataTagAndProbEff("El", "eta", "iso","");
  plotDataTagAndProbEff("El", "phi", "iso","");
  plotDataTagAndProbEff("El", "pt", "id","");
  plotDataTagAndProbEff("El", "eta", "id","");
  plotDataTagAndProbEff("El", "phi", "id","");

  plotDataTagAndProbEff("Mu", "pt", "iso","");
  plotDataTagAndProbEff("Mu", "eta", "iso","");
  plotDataTagAndProbEff("Mu", "phi", "iso","");
  plotDataTagAndProbEff("Mu", "pt", "id","");
  plotDataTagAndProbEff("Mu", "eta", "id","");
  plotDataTagAndProbEff("Mu", "phi", "id","");

}

