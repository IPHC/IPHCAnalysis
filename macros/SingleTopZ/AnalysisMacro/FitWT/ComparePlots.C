
#include <TH1F.h> 
#include <TFile.h> 
#include <TF1.h> 
#include <TLegend.h> 
#include <TMath.h> 
#include <THStack.h> 
#include <TCanvas.h> 
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include "TLegend.h"



void ComparePlots(TString channel, TString quantity, TString cut, TString dataset, bool isData){

  
  //***********************
  //Get histograms
  //***********************
   
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadLeftMargin(0.13); 
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  gStyle->SetStatColor(0); // kWhite

  
  TCanvas *c1 = new TCanvas("c1", "plots",400,400,800,600);
  c1->SetFillColor(10);
  c1->SetFillStyle(4000);
  c1->SetBorderSize(2); 
  c1->SetLogy(0);

  TFile *f_input2 = new TFile("../RootFiles/backup_outputProof22-12-12_11-12-22_nom_allSF_applied/proof_merged_AllSF.root");
  //TFile *f_input1 = new TFile("BackupRoot/proof_ZjetsTemplate.root");
  TFile *f_input1 ;
  
  if(isData) f_input1 = new TFile("../RootFiles/backup_outputProof20-12-12_10-03-23_noSFapplied_Zenriched/proof_merged_noSF_Zenriched.root");
  else       f_input1 = new TFile("../RootFiles/backup_outputProof22-12-12_11-12-22_nom_allSF_applied/proof_merged_AllSF.root");
  
  TString plotname1 = quantity + "_" + channel + "_"  +cut +  "_" + dataset;
  if(cut       == "") plotname1 = quantity + "_" + channel + "_" + dataset;
  if(quantity =="Inv")  plotname1 = quantity  + channel + "MassPair_" + cut +  "_" +dataset;
  if(quantity =="Inv" && cut == "" )  plotname1 = quantity  + channel + "MassPair_" +dataset;
  
  cout << "plotname 1 " << plotname1 << endl;
  
  f_input1->cd();
  
  TH1F * histo1 = (TH1F*)gROOT->FindObject(plotname1.Data());
  
  
  
  
  
  
  TString plotname2 = quantity + "_" + channel + "_"  +cut +  "_" + "WZ";
  if(cut == "") plotname2 = quantity + "_" + channel + "_" + "WZ";
  if(quantity =="Inv")  plotname2 = quantity  + channel + "MassPair_" + cut +  "_" +"WZ";
  if(quantity =="Inv" && cut == "" )  plotname2 = quantity  + channel + "MassPair_" +"WZ";
  
  f_input2->cd();
  
  TH1F * histo2 = (TH1F*)gROOT->FindObject(plotname2.Data());
  
  
  histo1->Rebin(5);
  histo2->Rebin(5);
  
  histo2->SetFillColor(5);
  histo1->SetMarkerSize(1.2);
  histo1->SetMarkerStyle(20);
  
  histo1->Sumw2();
  histo2->Sumw2();
  
  TLegend* qw = 0;
  if(quantity == "Njets" || quantity == "NBjets" || quantity == "DeltaPhiLLept") qw = new TLegend(.15,.70,.35,.85);
  else qw = new TLegend(.60,.70,.85,.85);
  
  qw->SetShadowColor(0);
  qw->SetFillColor(0);
  qw->SetLineColor(0);
  //qw->AddEntry(histo2,        "ttbar MG5" ,    "ep");
  if(isData) qw->AddEntry(histo1,        "data Z enriched","ep" );
  else       qw->AddEntry(histo1,        "signal FCNC gtu (FastSim)","ep" );
  qw->AddEntry(histo2,        "WZ MC" ,    "f");
 
  
  
  histo1->SetMinimum(0.0001);
  //histo1->SetMaximum(0.0001);
  
  if(isData){
  histo1->DrawNormalized("ep");
  histo2->DrawNormalized("hsame");
  }else{
  
  histo2->DrawNormalized("h");
  histo1->DrawNormalized("epsame");
  }
  
  qw->Draw();
  plotname1 = "plots/"+plotname1+".pdf";
  c1->SaveAs(plotname1.Data());
  

} 

void ComparePlots(){


 /*ComparePlots("mumumu", "mWT", "afterleptsel", "DataMu"   , 1);
 ComparePlots("mumue" , "mWT", "afterleptsel", "DataMuEG" , 1);
 ComparePlots("eemu"  , "mWT", "afterleptsel", "DataMuEG" , 1);
 ComparePlots("eee"   , "mWT", "afterleptsel", "DataEG"   , 1);
 
 */
 ComparePlots("mumumu", "mWT", "afterleptsel", "FCNCzut" , 0);
 ComparePlots("mumue" , "mWT", "afterleptsel", "FCNCzut" , 0);
 ComparePlots("eemu"  , "mWT", "afterleptsel", "FCNCzut" , 0);
 ComparePlots("eee"   , "mWT", "afterleptsel", "FCNCzut" , 0);
 
 //mWT_mumumu_afterleptsel_DataMu
 //mWT_mumumu_afterleptcut_DataMu

}



