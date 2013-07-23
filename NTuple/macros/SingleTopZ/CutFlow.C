//******************************************************
// Compute CutFlowTables from Proof output 
//******************************************************



#include <TH1F.h> 
#include <TH2F.h> 
#include <TFile.h> 
#include <TF1.h> 
#include <TLegend.h> 
#include <TMath.h> 
#include <THStack.h> 
#include <TCanvas.h> 
#include <TROOT.h>
#include <TStyle.h>
#include <fstream>
#include <vector>
#include <TLegend.h>
#include <TGraphErrors.h>



#include <iostream>

void CutFlow(){
  
  
  
  
  double SF_DY_mm = 1.;
  double SF_DY_ee = 1.;
  double SF_DY_em = 1.;
  
  double SF_Wjets_ee = 1.;
  double SF_Wjets_mm = 1.;
  double SF_Wjets_em = 1.;
  
  double SF_QCD_ee = 1.;
  double SF_QCD_mm = 1.;
  double SF_QCD_em = 1.;
  
  double TriggError[4]      = {0.029, 0.031, 0.038, 0.038};
  double SF_Lepton_error[4] = {0.001,  0.001, 0.001, 0.001};
  double SF_MET_error[4]    = { 0.0, 0.,  0., 0.};
  
  
  
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  
  //TCanvas *c1 = new TCanvas("c1", "plots",400,400,800,600);
  //c1->SetFillColor(10);
  //c1->SetFillStyle(4000);
  //c1->SetBorderSize(2); 
  //c1->SetLogy(0);
  
  //*********************************
  //get histograms 
  //*********************************
  
  
  TFile *f_data  = new TFile("backup_outputProof21-05-13_17-21/proof_merged.root");
  
  f_data->cd();
  
  
  //define tables [channel][process][cut step]
  double TabFlow1[5][103][101];
  double TabFlow2[5][103][101];
  
  
  for(int k0=0; k0<5; ++k0) {
    for(int k1=0; k1<103; k1++) {
      for(int k2=0; k2<101; ++k2) {
        TabFlow1[k0][k1][k2]=0.;
        TabFlow2[k0][k1][k2]=0.;
      }
    } 
  } 
  
  // mumumu  
  TH1F *  CutFlow_mumumu_DataMu     = (TH1F*)gROOT->FindObject("CutFlow_mumumu_DataMu");
  if ( CutFlow_mumumu_DataMu==NULL ) std::cout<<"WARNING "<<"CutFlow_mumumu_DataMu empty"<<std::endl;
  
  TH1F *  CutFlow_mumumu_FCNCzutFullSim   = (TH1F*)gROOT->FindObject("CutFlow_mumumu_FCNCzutFullSim");
  if ( CutFlow_mumumu_FCNCzutFullSim==NULL) std::cout<<"WARNING "<<"CutFlow_mumumu_FCNCzutFullSim empty"<<std::endl;
  
  TH1F *  CutFlow_mumumu_TTbarBkg   = (TH1F*)gROOT->FindObject("CutFlow_mumumu_TTbarBkg");
  if ( CutFlow_mumumu_TTbarBkg==NULL) std::cout<<"WARNING CutFlow_mumumu_TTbarBkg"<<" empty"<<std::endl;
  
  TH1F *  CutFlow_mumumu_Zjets      = (TH1F*)gROOT->FindObject("CutFlow_mumumu_Zjets");
  if (CutFlow_mumumu_Zjets ==NULL) std::cout<<"WARNING "<<"CutFlow_mumumu_Zjets empty"<<std::endl;
  
  TH1F *  CutFlow_mumumu_Wjets      = (TH1F*)gROOT->FindObject("CutFlow_mumumu_Wjets");
  if ( CutFlow_mumumu_Wjets==NULL) std::cout<<"WARNING "<<"CutFlow_mumumu_Wjets empty"<<std::endl;
  
  TH1F *  CutFlow_mumumu_TtW        = (TH1F*)gROOT->FindObject("CutFlow_mumumu_TtW");
  if ( CutFlow_mumumu_TtW==NULL) std::cout<<"WARNING "<<"CutFlow_mumumu_TtW empty"<<std::endl;
  
  TH1F *  CutFlow_mumumu_TbartW     = (TH1F*)gROOT->FindObject("CutFlow_mumumu_TbartW");
  if ( CutFlow_mumumu_TbartW==NULL) std::cout<<"WARNING "<<"CutFlow_mumumu_TbartW empty"<<std::endl;
  
  TH1F *  CutFlow_mumumu_TZq     = (TH1F*)gROOT->FindObject("CutFlow_mumumu_TZq");
  if ( CutFlow_mumumu_TZq==NULL) std::cout<<"WARNING "<<"CutFlow_mumumu_TZq empty"<<std::endl;
  
  TH1F *  CutFlow_mumumu_TtChan        = (TH1F*)gROOT->FindObject("CutFlow_mumumu_TtChan");
  if ( CutFlow_mumumu_TtChan==NULL) std::cout<<"WARNING "<<"CutFlow_mumumu_TtChan empty"<<std::endl;
  TH1F *  CutFlow_mumumu_TbartChan     = (TH1F*)gROOT->FindObject("CutFlow_mumumu_TbartChan");
  if ( CutFlow_mumumu_TbartChan==NULL) std::cout<<"WARNING "<<"CutFlow_mumumu_TbartChan empty"<<std::endl;
  
  TH1F *  CutFlow_mumumu_TsChan        = (TH1F*)gROOT->FindObject("CutFlow_mumumu_TsChan");
  if ( CutFlow_mumumu_TsChan==NULL) std::cout<<"WARNING "<<"CutFlow_mumumu_TsChan empty"<<std::endl;
  TH1F *  CutFlow_mumumu_TbarsChan     = (TH1F*)gROOT->FindObject("CutFlow_mumumu_TbarsChan");
  if ( CutFlow_mumumu_TbarsChan==NULL) std::cout<<"WARNING "<<"CutFlow_mumumu_TbarsChan empty"<<std::endl;
  
  
  TH1F *  CutFlow_mumumu_WW         = (TH1F*)gROOT->FindObject("CutFlow_mumumu_WW");
  if ( CutFlow_mumumu_WW==NULL) std::cout<<"WARNING "<<"CutFlow_mumumu_WW empty"<<std::endl;
  TH1F *  CutFlow_mumumu_WZ         = (TH1F*)gROOT->FindObject("CutFlow_mumumu_WZ");
  if ( CutFlow_mumumu_WZ==NULL) std::cout<<"WARNING "<<"CutFlow_mumumu_WZ empty"<<std::endl;
  TH1F *  CutFlow_mumumu_ZZ         = (TH1F*)gROOT->FindObject("CutFlow_mumumu_ZZ");
  if ( CutFlow_mumumu_ZZ==NULL) std::cout<<"WARNING "<<"CutFlow_mumumu_ZZ empty"<<std::endl;
  TH1F *  CutFlow_mumumu_DYToLL_M_10To50      = (TH1F*)gROOT->FindObject("CutFlow_mumumu_DYToLL_M10-50");
  if ( CutFlow_mumumu_DYToLL_M_10To50==NULL) std::cout<<"WARNING CutFlow_mumumu_DYToLL_M10-50"<<" empty"<<std::endl;
  
  
  TH1F *  ErrCutFlow_mumumu_DataMu     = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_DataMu");
  if ( ErrCutFlow_mumumu_DataMu==NULL ) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_DataMu empty"<<std::endl;
  TH1F *  ErrCutFlow_mumumu_FCNCzutFullSim   = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_FCNCzutFullSim");
  if ( ErrCutFlow_mumumu_FCNCzutFullSim==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_FCNCzutFullSim empty"<<std::endl;
  TH1F *  ErrCutFlow_mumumu_TTbarBkg   = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_TTbarBkg");
  if ( ErrCutFlow_mumumu_TTbarBkg==NULL) std::cout<<"WARNING ErrCutFlow_mumumu_TTbarBkg"<<" empty"<<std::endl;
  TH1F *  ErrCutFlow_mumumu_Zjets      = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_Zjets");
  if (ErrCutFlow_mumumu_Zjets ==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_Zjets empty"<<std::endl;
  TH1F *  ErrCutFlow_mumumu_Wjets      = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_Wjets");
  if ( ErrCutFlow_mumumu_Wjets==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_Wjets empty"<<std::endl;


  TH1F *  ErrCutFlow_mumumu_TZq        = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_TZq");
  if ( ErrCutFlow_mumumu_TZq==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_TZq empty"<<std::endl;
  
  TH1F *  ErrCutFlow_mumumu_TtW        = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_TtW");
  if ( ErrCutFlow_mumumu_TtW==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_TtW empty"<<std::endl;
  TH1F *  ErrCutFlow_mumumu_TbartW     = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_TbartW");
  if ( ErrCutFlow_mumumu_TbartW==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_TbartW empty"<<std::endl;

  TH1F *  ErrCutFlow_mumumu_TtChan        = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_TtChan");
  if ( ErrCutFlow_mumumu_TtChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_TtChan empty"<<std::endl;
  TH1F *  ErrCutFlow_mumumu_TbartChan     = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_TbartChan");
  if ( ErrCutFlow_mumumu_TbartChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_TbartChan empty"<<std::endl;

  TH1F *  ErrCutFlow_mumumu_TsChan        = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_TsChan");
  if ( ErrCutFlow_mumumu_TsChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_TsChan empty"<<std::endl;
  TH1F *  ErrCutFlow_mumumu_TbarsChan     = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_TbarsChan");
  if ( ErrCutFlow_mumumu_TbarsChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_TbarsChan empty"<<std::endl;


  TH1F *  ErrCutFlow_mumumu_WW         = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_WW");
  if ( ErrCutFlow_mumumu_WW==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_WW empty"<<std::endl;
  TH1F *  ErrCutFlow_mumumu_WZ         = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_WZ");
  if ( ErrCutFlow_mumumu_WZ==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_WZ empty"<<std::endl;
  TH1F *  ErrCutFlow_mumumu_ZZ         = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_ZZ");
  if ( ErrCutFlow_mumumu_ZZ==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumumu_ZZ empty"<<std::endl;
  TH1F *  ErrCutFlow_mumumu_DYToLL_M_10To50      = (TH1F*)gROOT->FindObject("ErrCutFlow_mumumu_DYToLL_M10-50");
  if ( ErrCutFlow_mumumu_DYToLL_M_10To50==NULL) std::cout<<"WARNING ErrCutFlow_mumumu_DYToLL_M10-50"<<" empty"<<std::endl;
  
  
  // mumue  
  TH1F *  CutFlow_mumue_DataMu     = (TH1F*)gROOT->FindObject("CutFlow_mumue_DataMuEG");
  if ( CutFlow_mumue_DataMu==NULL ) std::cout<<"WARNING "<<"CutFlow_mumue_DataMuEG empty"<<std::endl;
  TH1F *  CutFlow_mumue_FCNCzutFullSim   = (TH1F*)gROOT->FindObject("CutFlow_mumue_FCNCzutFullSim");
  if ( CutFlow_mumue_FCNCzutFullSim==NULL) std::cout<<"WARNING "<<"CutFlow_mumue_FCNCzutFullSim empty"<<std::endl;
  TH1F *  CutFlow_mumue_TTbarBkg   = (TH1F*)gROOT->FindObject("CutFlow_mumue_TTbarBkg");
  if ( CutFlow_mumue_TTbarBkg==NULL) std::cout<<"WARNING CutFlow_mumue_TTbarBkg"<<" empty"<<std::endl;
  TH1F *  CutFlow_mumue_Zjets      = (TH1F*)gROOT->FindObject("CutFlow_mumue_Zjets");
  if (CutFlow_mumue_Zjets ==NULL) std::cout<<"WARNING "<<"CutFlow_mumue_Zjets empty"<<std::endl;
  TH1F *  CutFlow_mumue_Wjets      = (TH1F*)gROOT->FindObject("CutFlow_mumue_Wjets");
  if ( CutFlow_mumue_Wjets==NULL) std::cout<<"WARNING "<<"CutFlow_mumue_Wjets empty"<<std::endl;
  
  TH1F *  CutFlow_mumue_TZq        = (TH1F*)gROOT->FindObject("CutFlow_mumue_TZq");
  if ( CutFlow_mumue_TZq==NULL) std::cout<<"WARNING "<<"CutFlow_mumue_TZq empty"<<std::endl;
  
  TH1F *  CutFlow_mumue_TtW        = (TH1F*)gROOT->FindObject("CutFlow_mumue_TtW");
  if ( CutFlow_mumue_TtW==NULL) std::cout<<"WARNING "<<"CutFlow_mumue_TtW empty"<<std::endl;
  TH1F *  CutFlow_mumue_TbartW     = (TH1F*)gROOT->FindObject("CutFlow_mumue_TbartW");
  
  TH1F *  CutFlow_mumue_TtChan        = (TH1F*)gROOT->FindObject("CutFlow_mumue_TtChan");
  if ( CutFlow_mumue_TtChan==NULL) std::cout<<"WARNING "<<"CutFlow_mumue_TtChan empty"<<std::endl;
  TH1F *  CutFlow_mumue_TbartChan     = (TH1F*)gROOT->FindObject("CutFlow_mumue_TbartChan");
  
  TH1F *  CutFlow_mumue_TsChan        = (TH1F*)gROOT->FindObject("CutFlow_mumue_TsChan");
  if ( CutFlow_mumue_TsChan==NULL) std::cout<<"WARNING "<<"CutFlow_mumue_TsChan empty"<<std::endl;
  TH1F *  CutFlow_mumue_TbarsChan     = (TH1F*)gROOT->FindObject("CutFlow_mumue_TbarsChan");
  
  
  
  if ( CutFlow_mumue_TbartW==NULL) std::cout<<"WARNING "<<"CutFlow_mumue_TbartW empty"<<std::endl;
  TH1F *  CutFlow_mumue_WW         = (TH1F*)gROOT->FindObject("CutFlow_mumue_WW");
  if ( CutFlow_mumue_WW==NULL) std::cout<<"WARNING "<<"CutFlow_mumue_WW empty"<<std::endl;
  TH1F *  CutFlow_mumue_WZ         = (TH1F*)gROOT->FindObject("CutFlow_mumue_WZ");
  if ( CutFlow_mumue_WZ==NULL) std::cout<<"WARNING "<<"CutFlow_mumue_WZ empty"<<std::endl;
  TH1F *  CutFlow_mumue_ZZ         = (TH1F*)gROOT->FindObject("CutFlow_mumue_ZZ");
  if ( CutFlow_mumue_ZZ==NULL) std::cout<<"WARNING "<<"CutFlow_mumue_ZZ empty"<<std::endl;
  TH1F *  CutFlow_mumue_DYToLL_M_10To50      = (TH1F*)gROOT->FindObject("CutFlow_mumue_DYToLL_M10-50");
  if ( CutFlow_mumue_DYToLL_M_10To50==NULL) std::cout<<"WARNING CutFlow_mumue_DYToLL_M10-50"<<" empty"<<std::endl;
  
  
  TH1F *  ErrCutFlow_mumue_DataMu     = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_DataMuEG");
  if ( ErrCutFlow_mumue_DataMu==NULL ) std::cout<<"WARNING "<<"ErrCutFlow_mumue_DataMuEG empty"<<std::endl;
  TH1F *  ErrCutFlow_mumue_FCNCzutFullSim   = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_FCNCzutFullSim");
  if ( ErrCutFlow_mumue_FCNCzutFullSim==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumue_FCNCzutFullSim empty"<<std::endl;
  TH1F *  ErrCutFlow_mumue_TTbarBkg   = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_TTbarBkg");
  if ( ErrCutFlow_mumue_TTbarBkg==NULL) std::cout<<"WARNING ErrCutFlow_mumue_TTbarBkg"<<" empty"<<std::endl;
  TH1F *  ErrCutFlow_mumue_Zjets      = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_Zjets");
  if (ErrCutFlow_mumue_Zjets ==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumue_Zjets empty"<<std::endl;
  TH1F *  ErrCutFlow_mumue_Wjets      = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_Wjets");
  if ( ErrCutFlow_mumue_Wjets==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumue_Wjets empty"<<std::endl;
  
  TH1F *  ErrCutFlow_mumue_TZq        = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_TZq");
  if ( ErrCutFlow_mumue_TZq==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumue_TZq empty"<<std::endl;
  
  TH1F *  ErrCutFlow_mumue_TtW        = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_TtW");
  if ( ErrCutFlow_mumue_TtW==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumue_TtW empty"<<std::endl;
  TH1F *  ErrCutFlow_mumue_TbartW     = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_TbartW");
  if ( ErrCutFlow_mumue_TbartW==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumue_TbartW empty"<<std::endl;
  
  
  TH1F *  ErrCutFlow_mumue_TtChan        = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_TtChan");
  if ( ErrCutFlow_mumue_TtChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumue_TtChan empty"<<std::endl;
  TH1F *  ErrCutFlow_mumue_TbartChan     = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_TbartChan");
  if ( ErrCutFlow_mumue_TbartChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumue_TbartChan empty"<<std::endl;
  
  TH1F *  ErrCutFlow_mumue_TsChan        = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_TsChan");
  if ( ErrCutFlow_mumue_TsChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumue_TsChan empty"<<std::endl;
  TH1F *  ErrCutFlow_mumue_TbarsChan     = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_TbarsChan");
  if ( ErrCutFlow_mumue_TbarsChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumue_TbarsChan empty"<<std::endl;
  
  TH1F *  ErrCutFlow_mumue_WW         = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_WW");
  if ( ErrCutFlow_mumue_WW==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumue_WW empty"<<std::endl;
  TH1F *  ErrCutFlow_mumue_WZ         = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_WZ");
  if ( ErrCutFlow_mumue_WZ==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumue_WZ empty"<<std::endl;
  TH1F *  ErrCutFlow_mumue_ZZ         = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_ZZ");
  if ( ErrCutFlow_mumue_ZZ==NULL) std::cout<<"WARNING "<<"ErrCutFlow_mumue_ZZ empty"<<std::endl;
  TH1F *  ErrCutFlow_mumue_DYToLL_M_10To50      = (TH1F*)gROOT->FindObject("ErrCutFlow_mumue_DYToLL_M10-50");
  if ( ErrCutFlow_mumue_DYToLL_M_10To50==NULL) std::cout<<"WARNING ErrCutFlow_mumue_DYToLL_M10-50"<<" empty"<<std::endl;
  
  
  //eemu
  TH1F *  CutFlow_eemu_DataEG     = (TH1F*)gROOT->FindObject("CutFlow_eemu_DataMuEG");
  if ( CutFlow_eemu_DataEG==NULL ) std::cout<<"WARNING "<<"CutFlow_eemu_DataMuEG empty"<<std::endl;
  TH1F *  CutFlow_eemu_FCNCzutFullSim   = (TH1F*)gROOT->FindObject("CutFlow_eemu_FCNCzutFullSim");
  if ( CutFlow_eemu_FCNCzutFullSim==NULL) std::cout<<"WARNING "<<"CutFlow_eemu_FCNCzutFullSim empty"<<std::endl;
  TH1F *  CutFlow_eemu_TTbarBkg   = (TH1F*)gROOT->FindObject("CutFlow_eemu_TTbarBkg");
  if ( CutFlow_eemu_TTbarBkg==NULL) std::cout<<"WARNING CutFlow_eemu_TTbarBkg"<<" empty"<<std::endl;
  TH1F *  CutFlow_eemu_Zjets      = (TH1F*)gROOT->FindObject("CutFlow_eemu_Zjets");
  if (CutFlow_eemu_Zjets ==NULL) std::cout<<"WARNING "<<"CutFlow_eemu_Zjets empty"<<std::endl;
  TH1F *  CutFlow_eemu_Wjets      = (TH1F*)gROOT->FindObject("CutFlow_eemu_Wjets");
  if ( CutFlow_eemu_Wjets==NULL) std::cout<<"WARNING "<<"CutFlow_eemu_Wjets empty"<<std::endl;
  
  
  TH1F *  CutFlow_eemu_TZq        = (TH1F*)gROOT->FindObject("CutFlow_eemu_TZq");
  if ( CutFlow_eemu_TZq==NULL) std::cout<<"WARNING "<<"CutFlow_eemu_TZq empty"<<std::endl;
  
  TH1F *  CutFlow_eemu_TtW        = (TH1F*)gROOT->FindObject("CutFlow_eemu_TtW");
  if ( CutFlow_eemu_TtW==NULL) std::cout<<"WARNING "<<"CutFlow_eemu_TtW empty"<<std::endl;
  TH1F *  CutFlow_eemu_TbartW     = (TH1F*)gROOT->FindObject("CutFlow_eemu_TbartW");  
  if ( CutFlow_eemu_TbartW==NULL) std::cout<<"WARNING "<<"CutFlow_eemu_TbartW empty"<<std::endl;
  
  TH1F *  CutFlow_eemu_TtChan        = (TH1F*)gROOT->FindObject("CutFlow_eemu_TtChan");
  if ( CutFlow_eemu_TtChan==NULL) std::cout<<"WARNING "<<"CutFlow_eemu_TtChan empty"<<std::endl;
  TH1F *  CutFlow_eemu_TbartChan     = (TH1F*)gROOT->FindObject("CutFlow_eemu_TbartChan");  
  if ( CutFlow_eemu_TbartChan==NULL) std::cout<<"WARNING "<<"CutFlow_eemu_TbartChan empty"<<std::endl;
  
  TH1F *  CutFlow_eemu_TsChan        = (TH1F*)gROOT->FindObject("CutFlow_eemu_TsChan");
  if ( CutFlow_eemu_TsChan==NULL) std::cout<<"WARNING "<<"CutFlow_eemu_TsChan empty"<<std::endl;
  TH1F *  CutFlow_eemu_TbarsChan     = (TH1F*)gROOT->FindObject("CutFlow_eemu_TbarsChan");  
  if ( CutFlow_eemu_TbarsChan==NULL) std::cout<<"WARNING "<<"CutFlow_eemu_TbarsChan empty"<<std::endl;
  
  
  
  TH1F *  CutFlow_eemu_WW         = (TH1F*)gROOT->FindObject("CutFlow_eemu_WW");
  if ( CutFlow_eemu_WW==NULL) std::cout<<"WARNING "<<"CutFlow_eemu_WW empty"<<std::endl;
  TH1F *  CutFlow_eemu_WZ         = (TH1F*)gROOT->FindObject("CutFlow_eemu_WZ");
  if ( CutFlow_eemu_WZ==NULL) std::cout<<"WARNING "<<"CutFlow_eemu_WZ empty"<<std::endl;
  TH1F *  CutFlow_eemu_ZZ         = (TH1F*)gROOT->FindObject("CutFlow_eemu_ZZ");
  if ( CutFlow_eemu_ZZ==NULL) std::cout<<"WARNING "<<"CutFlow_eemu_ZZ empty"<<std::endl;
  TH1F *  CutFlow_eemu_DYToLL_M_10To50      = (TH1F*)gROOT->FindObject("CutFlow_eemu_DYToLL_M10-50");
  if ( CutFlow_eemu_DYToLL_M_10To50==NULL) std::cout<<"WARNING CutFlow_eemu_DYToLL_M10-50"<<" empty"<<std::endl;
  
  TH1F *  ErrCutFlow_eemu_DataEG     = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_DataMuEG");
  if ( ErrCutFlow_eemu_DataEG==NULL ) std::cout<<"WARNING "<<"ErrCutFlow_eemu_DataMuEG empty"<<std::endl;
  TH1F *  ErrCutFlow_eemu_FCNCzutFullSim   = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_FCNCzutFullSim");
  if ( ErrCutFlow_eemu_FCNCzutFullSim==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eemu_FCNCzutFullSim empty"<<std::endl;
  TH1F *  ErrCutFlow_eemu_TTbarBkg   = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_TTbarBkg");
  if ( ErrCutFlow_eemu_TTbarBkg==NULL) std::cout<<"WARNING ErrCutFlow_eemu_TTbarBkg"<<" empty"<<std::endl;
  TH1F *  ErrCutFlow_eemu_Zjets      = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_Zjets");
  if (ErrCutFlow_eemu_Zjets ==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eemu_Zjets empty"<<std::endl;
  TH1F *  ErrCutFlow_eemu_Wjets      = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_Wjets");
  if ( ErrCutFlow_eemu_Wjets==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eemu_Wjets empty"<<std::endl;
  
  
  TH1F *  ErrCutFlow_eemu_TZq        = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_TZq");
  if ( ErrCutFlow_eemu_TZq==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eemu_TZq empty"<<std::endl;
  
  
  
  TH1F *  ErrCutFlow_eemu_TtW        = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_TtW");
  if ( ErrCutFlow_eemu_TtW==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eemu_TtW empty"<<std::endl;
  
  
  
  TH1F *  ErrCutFlow_eemu_TbartW     = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_TbartW");
  if ( ErrCutFlow_eemu_TbartW==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eemu_TbartW empty"<<std::endl;
  
  TH1F *  ErrCutFlow_eemu_TtChan        = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_TtChan");
  if ( ErrCutFlow_eemu_TtChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eemu_TtChan empty"<<std::endl;
  TH1F *  ErrCutFlow_eemu_TbartChan     = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_TbartChan");
  if ( ErrCutFlow_eemu_TbartChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eemu_TbartChan empty"<<std::endl;
  
  TH1F *  ErrCutFlow_eemu_TsChan        = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_TsChan");
  if ( ErrCutFlow_eemu_TsChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eemu_TsChan empty"<<std::endl;
  TH1F *  ErrCutFlow_eemu_TbarsChan     = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_TbarsChan");
  if ( ErrCutFlow_eemu_TbarsChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eemu_TbarsChan empty"<<std::endl;
  
  
  TH1F *  ErrCutFlow_eemu_WW         = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_WW");
  if ( ErrCutFlow_eemu_WW==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eemu_WW empty"<<std::endl;
  TH1F *  ErrCutFlow_eemu_WZ         = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_WZ");
  if ( ErrCutFlow_eemu_WZ==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eemu_WZ empty"<<std::endl;
  TH1F *  ErrCutFlow_eemu_ZZ         = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_ZZ");
  if ( ErrCutFlow_eemu_ZZ==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eemu_ZZ empty"<<std::endl;
  TH1F *  ErrCutFlow_eemu_DYToLL_M_10To50      = (TH1F*)gROOT->FindObject("ErrCutFlow_eemu_DYToLL_M10-50");
  if ( ErrCutFlow_eemu_DYToLL_M_10To50==NULL) std::cout<<"WARNING ErrCutFlow_eemu_DYToLL_M10-50"<<" empty"<<std::endl;
  
  // eee   
  TH1F *  CutFlow_eee_DataEG     = (TH1F*)gROOT->FindObject("CutFlow_eee_DataEG");
  if ( CutFlow_eee_DataEG==NULL ) std::cout<<"WARNING "<<"CutFlow_eee_DataEG empty"<<std::endl;
  TH1F *  CutFlow_eee_FCNCzutFullSim   = (TH1F*)gROOT->FindObject("CutFlow_eee_FCNCzutFullSim");
  if ( CutFlow_eee_FCNCzutFullSim==NULL) std::cout<<"WARNING "<<"CutFlow_eee_FCNCzutFullSim empty"<<std::endl;
  TH1F *  CutFlow_eee_TTbarBkg   = (TH1F*)gROOT->FindObject("CutFlow_eee_TTbarBkg");
  if ( CutFlow_eee_TTbarBkg==NULL) std::cout<<"WARNING CutFlow_eee_TTbarBkg"<<" empty"<<std::endl;
  TH1F *  CutFlow_eee_Zjets      = (TH1F*)gROOT->FindObject("CutFlow_eee_Zjets");
  if (CutFlow_eee_Zjets ==NULL) std::cout<<"WARNING "<<"CutFlow_eee_Zjets empty"<<std::endl;
  TH1F *  CutFlow_eee_Wjets      = (TH1F*)gROOT->FindObject("CutFlow_eee_Wjets");
  if ( CutFlow_eee_Wjets==NULL) std::cout<<"WARNING "<<"CutFlow_eee_Wjets empty"<<std::endl;
  
  TH1F *  CutFlow_eee_TZq        = (TH1F*)gROOT->FindObject("CutFlow_eee_TZq");
  if ( CutFlow_eee_TZq==NULL) std::cout<<"WARNING "<<"CutFlow_eee_TZq empty"<<std::endl;
  
  TH1F *  CutFlow_eee_TtW        = (TH1F*)gROOT->FindObject("CutFlow_eee_TtW");
  if ( CutFlow_eee_TtW==NULL) std::cout<<"WARNING "<<"CutFlow_eee_TtW empty"<<std::endl;
  TH1F *  CutFlow_eee_TbartW     = (TH1F*)gROOT->FindObject("CutFlow_eee_TbartW");
  if ( CutFlow_eee_TbartW==NULL) std::cout<<"WARNING "<<"CutFlow_eee_TbartW empty"<<std::endl;
  
  
  TH1F *  CutFlow_eee_TtChan        = (TH1F*)gROOT->FindObject("CutFlow_eee_TtChan");
  if ( CutFlow_eee_TtChan==NULL) std::cout<<"WARNING "<<"CutFlow_eee_TtChan empty"<<std::endl;
  TH1F *  CutFlow_eee_TbartChan     = (TH1F*)gROOT->FindObject("CutFlow_eee_TbartChan");
  if ( CutFlow_eee_TbartChan==NULL) std::cout<<"WARNING "<<"CutFlow_eee_TbartChan empty"<<std::endl;
  
  
  
  TH1F *  CutFlow_eee_TsChan        = (TH1F*)gROOT->FindObject("CutFlow_eee_TsChan");
  if ( CutFlow_eee_TsChan==NULL) std::cout<<"WARNING "<<"CutFlow_eee_TsChan empty"<<std::endl;
  TH1F *  CutFlow_eee_TbarsChan     = (TH1F*)gROOT->FindObject("CutFlow_eee_TbarsChan");
  if ( CutFlow_eee_TbarsChan==NULL) std::cout<<"WARNING "<<"CutFlow_eee_TbarsChan empty"<<std::endl;
  
  TH1F *  CutFlow_eee_WW         = (TH1F*)gROOT->FindObject("CutFlow_eee_WW");
  if ( CutFlow_eee_WW==NULL) std::cout<<"WARNING "<<"CutFlow_eee_WW empty"<<std::endl;
  TH1F *  CutFlow_eee_WZ         = (TH1F*)gROOT->FindObject("CutFlow_eee_WZ");
  if ( CutFlow_eee_WZ==NULL) std::cout<<"WARNING "<<"CutFlow_eee_WZ empty"<<std::endl;
  TH1F *  CutFlow_eee_ZZ         = (TH1F*)gROOT->FindObject("CutFlow_eee_ZZ");
  if ( CutFlow_eee_ZZ==NULL) std::cout<<"WARNING "<<"CutFlow_eee_ZZ empty"<<std::endl;
  TH1F *  CutFlow_eee_DYToLL_M_10To50      = (TH1F*)gROOT->FindObject("CutFlow_eee_DYToLL_M10-50");
  if ( CutFlow_eee_DYToLL_M_10To50==NULL) std::cout<<"WARNING CutFlow_eee_DYToLL_M10-50"<<" empty"<<std::endl;
  
  TH1F *  ErrCutFlow_eee_DataEG     = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_DataEG");
  if ( ErrCutFlow_eee_DataEG==NULL ) std::cout<<"WARNING "<<"ErrCutFlow_eee_DataEG empty"<<std::endl;
  TH1F *  ErrCutFlow_eee_FCNCzutFullSim   = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_FCNCzutFullSim");
  if ( ErrCutFlow_eee_FCNCzutFullSim==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eee_FCNCzutFullSim empty"<<std::endl;
  TH1F *  ErrCutFlow_eee_TTbarBkg   = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_TTbarBkg");
  if ( ErrCutFlow_eee_TTbarBkg==NULL) std::cout<<"WARNING ErrCutFlow_eee_TTbarBkg"<<" empty"<<std::endl;
  TH1F *  ErrCutFlow_eee_Zjets      = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_Zjets");
  if (ErrCutFlow_eee_Zjets ==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eee_Zjets empty"<<std::endl;
  TH1F *  ErrCutFlow_eee_Wjets      = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_Wjets");
  if ( ErrCutFlow_eee_Wjets==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eee_Wjets empty"<<std::endl;
  
  
  TH1F *  ErrCutFlow_eee_TZq        = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_TZq");
  if ( ErrCutFlow_eee_TZq==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eee_TZq empty"<<std::endl;
  
  
  TH1F *  ErrCutFlow_eee_TtW        = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_TtW");
  if ( ErrCutFlow_eee_TtW==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eee_TtW empty"<<std::endl;
  
  TH1F *  ErrCutFlow_eee_TbartW     = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_TbartW");
  if ( ErrCutFlow_eee_TbartW==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eee_TbartW empty"<<std::endl;
  
  TH1F *  ErrCutFlow_eee_TtChan        = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_TtChan");
  if ( ErrCutFlow_eee_TtChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eee_TtChan empty"<<std::endl;
  TH1F *  ErrCutFlow_eee_TbartChan     = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_TbartChan");
  if ( ErrCutFlow_eee_TbartChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eee_TbartChan empty"<<std::endl;
  
  TH1F *  ErrCutFlow_eee_TsChan        = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_TsChan");
  if ( ErrCutFlow_eee_TsChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eee_TsChan empty"<<std::endl;
  TH1F *  ErrCutFlow_eee_TbarsChan     = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_TbarsChan");
  if ( ErrCutFlow_eee_TbarsChan==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eee_TbarsChan empty"<<std::endl;
  
  
  TH1F *  ErrCutFlow_eee_WW         = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_WW");
  if ( ErrCutFlow_eee_WW==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eee_WW empty"<<std::endl;
  TH1F *  ErrCutFlow_eee_WZ         = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_WZ");
  if ( ErrCutFlow_eee_WZ==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eee_WZ empty"<<std::endl;
  TH1F *  ErrCutFlow_eee_ZZ         = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_ZZ");
  if ( ErrCutFlow_eee_ZZ==NULL) std::cout<<"WARNING "<<"ErrCutFlow_eee_ZZ empty"<<std::endl;
  TH1F *  ErrCutFlow_eee_DYToLL_M_10To50      = (TH1F*)gROOT->FindObject("ErrCutFlow_eee_DYToLL_M10-50");
  if ( ErrCutFlow_eee_DYToLL_M_10To50==NULL) std::cout<<"WARNING ErrCutFlow_eee_DYToLL_M10-50"<<" empty"<<std::endl;
  
  // SFtrigger  
  TH1F *  SFtrigger     = (TH1F*)gROOT->FindObject("SFtrigger");
  if ( SFtrigger==NULL ) std::cout<<"WARNING "<<"SFtrigger empty"<<std::endl;
  
  std::cout<<"================================="<<std::endl;
  std::cout<<"SF_DY_mm = "<<SF_DY_mm<<std::endl;
  std::cout<<"SF_DY_ee = "<<SF_DY_ee<<std::endl;
  std::cout<<"SF_DY_em = "<<SF_DY_em<<std::endl;
  std::cout<<"================================="<<std::endl;
  std::string rep = "";
  std::cout << "Inputs : SF DY corrects ? [y/n] "<<std::endl;
  cin >> rep ;
  
  
  
  
      CutFlow_mumumu_FCNCzutFullSim->Scale(0.1);
      CutFlow_mumue_FCNCzutFullSim->Scale(0.1);
      CutFlow_eemu_FCNCzutFullSim->Scale(0.1);
      CutFlow_eee_FCNCzutFullSim->Scale(0.1);
      ErrCutFlow_mumumu_FCNCzutFullSim->Scale(0.1);
      ErrCutFlow_mumue_FCNCzutFullSim->Scale(0.1);
      ErrCutFlow_eemu_FCNCzutFullSim->Scale(0.1);
      ErrCutFlow_eee_FCNCzutFullSim->Scale(0.1);
      
      
      CutFlow_mumumu_TZq->Scale(0.27);
      CutFlow_mumue_TZq->Scale(0.27);
      CutFlow_eemu_TZq->Scale(0.27);
      CutFlow_eee_TZq->Scale(0.27);
      ErrCutFlow_mumumu_TZq->Scale(0.27);
      ErrCutFlow_mumue_TZq->Scale(0.27);
      ErrCutFlow_eemu_TZq->Scale(0.27);
      ErrCutFlow_eee_TZq->Scale(0.27);
      
      
  if ( rep=="y" ){ 
    
    
    //*********************************
    // 
    //*********************************
    
    for(int i=0; i<12; ++i) {
      
      
      
      // mumumu channel  
      if ( CutFlow_mumumu_FCNCzutFullSim!=NULL ) { 
        TabFlow1[0][0][i] = CutFlow_mumumu_FCNCzutFullSim->GetBinContent(i+1);
	TabFlow2[0][0][i] = ErrCutFlow_mumumu_FCNCzutFullSim->GetBinContent(i+1);
        TabFlow1[4][0][i] += CutFlow_mumumu_FCNCzutFullSim->GetBinContent(i+1);
	TabFlow2[4][0][i] += ErrCutFlow_mumumu_FCNCzutFullSim->GetBinContent(i+1);
      }	
      if ( CutFlow_mumumu_TTbarBkg!=NULL ) { 
        TabFlow1[0][1][i] = CutFlow_mumumu_TTbarBkg->GetBinContent(i+1);
        TabFlow2[0][1][i] = ErrCutFlow_mumumu_TTbarBkg->GetBinContent(i+1);
        TabFlow1[4][1][i] += CutFlow_mumumu_TTbarBkg->GetBinContent(i+1);
        TabFlow2[4][1][i] += ErrCutFlow_mumumu_TTbarBkg->GetBinContent(i+1);
      }	
      
      if ( CutFlow_mumumu_Zjets!=NULL ) {
        TabFlow1[0][2][i] = CutFlow_mumumu_Zjets->GetBinContent(i+1);
	TabFlow2[0][2][i] = ErrCutFlow_mumumu_Zjets->GetBinContent(i+1);
        TabFlow1[4][2][i] += CutFlow_mumumu_Zjets->GetBinContent(i+1);
	TabFlow2[4][2][i] += ErrCutFlow_mumumu_Zjets->GetBinContent(i+1);
      }		 
      if ( CutFlow_mumumu_DYToLL_M_10To50!=NULL ) {
	TabFlow1[0][2][i] += CutFlow_mumumu_DYToLL_M_10To50->GetBinContent(i+1) ;
	TabFlow2[0][2][i] += ErrCutFlow_mumumu_DYToLL_M_10To50->GetBinContent(i+1) ; 
	TabFlow1[4][2][i] += CutFlow_mumumu_DYToLL_M_10To50->GetBinContent(i+1) ;
	TabFlow2[4][2][i] += ErrCutFlow_mumumu_DYToLL_M_10To50->GetBinContent(i+1) ; 
      }	 
      
      
      
      if ( CutFlow_mumumu_Wjets!=NULL ) { 
        TabFlow1[0][3][i] = CutFlow_mumumu_Wjets->GetBinContent(i+1);
	TabFlow2[0][3][i] = ErrCutFlow_mumumu_Wjets->GetBinContent(i+1);
        TabFlow1[4][3][i] += CutFlow_mumumu_Wjets->GetBinContent(i+1);
	TabFlow2[4][3][i] += ErrCutFlow_mumumu_Wjets->GetBinContent(i+1);
      }	
      
      
      if ( CutFlow_mumumu_TtW!=NULL ) { 
        TabFlow1[0][4][i] = CutFlow_mumumu_TtW->GetBinContent(i+1);
	TabFlow2[0][4][i] = ErrCutFlow_mumumu_TtW->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_mumumu_TtW->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_mumumu_TtW->GetBinContent(i+1);
      }	
      if ( CutFlow_mumumu_TbartW!=NULL ) {
        TabFlow1[0][4][i] += CutFlow_mumumu_TbartW->GetBinContent(i+1);
	TabFlow2[0][4][i] += ErrCutFlow_mumumu_TbartW->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_mumumu_TbartW->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_mumumu_TbartW->GetBinContent(i+1);
      }
      
      
      if ( CutFlow_mumumu_TtChan!=NULL ) { 
        TabFlow1[0][4][i] = CutFlow_mumumu_TtChan->GetBinContent(i+1);
	TabFlow2[0][4][i] = ErrCutFlow_mumumu_TtChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_mumumu_TtChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_mumumu_TtChan->GetBinContent(i+1);
      }	
      if ( CutFlow_mumumu_TbartChan!=NULL ) {
        TabFlow1[0][4][i] += CutFlow_mumumu_TbartChan->GetBinContent(i+1);
	TabFlow2[0][4][i] += ErrCutFlow_mumumu_TbartChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_mumumu_TbartChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_mumumu_TbartChan->GetBinContent(i+1);
      }
      
      if ( CutFlow_mumumu_TsChan!=NULL ) { 
        TabFlow1[0][4][i] = CutFlow_mumumu_TsChan->GetBinContent(i+1);
	TabFlow2[0][4][i] = ErrCutFlow_mumumu_TsChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_mumumu_TsChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_mumumu_TsChan->GetBinContent(i+1);
      }	
      if ( CutFlow_mumumu_TbarsChan!=NULL ) {
        TabFlow1[0][4][i] += CutFlow_mumumu_TbarsChan->GetBinContent(i+1);
	TabFlow2[0][4][i] += ErrCutFlow_mumumu_TbarsChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_mumumu_TbarsChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_mumumu_TbarsChan->GetBinContent(i+1);
      }
      
      
      
      	
      if ( CutFlow_mumumu_WZ!=NULL ) {
        TabFlow1[0][5][i] = CutFlow_mumumu_WZ->GetBinContent(i+1) ;
	TabFlow2[0][5][i] = ErrCutFlow_mumumu_WZ->GetBinContent(i+1) ;
        TabFlow1[4][5][i] += CutFlow_mumumu_WZ->GetBinContent(i+1) ;
	TabFlow2[4][5][i] += ErrCutFlow_mumumu_WZ->GetBinContent(i+1) ;
      }	
      if ( CutFlow_mumumu_ZZ!=NULL ) { 
        TabFlow1[0][5][i] += CutFlow_mumumu_ZZ->GetBinContent(i+1) ;
	TabFlow2[0][5][i] += ErrCutFlow_mumumu_ZZ->GetBinContent(i+1) ;
        TabFlow1[4][5][i] += CutFlow_mumumu_ZZ->GetBinContent(i+1) ;
	TabFlow2[4][5][i] += ErrCutFlow_mumumu_ZZ->GetBinContent(i+1) ;
      }	
      if ( CutFlow_mumumu_WW!=NULL ) { 
        TabFlow1[0][5][i] += CutFlow_mumumu_WW->GetBinContent(i+1) ;
	TabFlow2[0][5][i] += ErrCutFlow_mumumu_WW->GetBinContent(i+1) ;
        TabFlow1[4][5][i] += CutFlow_mumumu_WW->GetBinContent(i+1) ;
	TabFlow2[4][5][i] += ErrCutFlow_mumumu_WW->GetBinContent(i+1) ;
      }	
      if ( CutFlow_mumumu_TZq!=NULL ) { 
        TabFlow1[0][6][i] = CutFlow_mumumu_TZq->GetBinContent(i+1);
	TabFlow2[0][6][i] = ErrCutFlow_mumumu_TZq->GetBinContent(i+1);
        TabFlow1[4][6][i] += CutFlow_mumumu_TZq->GetBinContent(i+1);
	TabFlow2[4][6][i] += ErrCutFlow_mumumu_TZq->GetBinContent(i+1);
      }	
      
      if ( CutFlow_mumumu_DataMu!=NULL ) {
        TabFlow1[0][100][i] = CutFlow_mumumu_DataMu->GetBinContent(i+1);
	TabFlow2[0][100][i] = ErrCutFlow_mumumu_DataMu->GetBinContent(i+1);
        TabFlow1[4][100][i] += CutFlow_mumumu_DataMu->GetBinContent(i+1);
	TabFlow2[4][100][i] += ErrCutFlow_mumumu_DataMu->GetBinContent(i+1);
      }	
      
      
      // mumue channel  
      if ( CutFlow_mumue_FCNCzutFullSim!=NULL ) { 
        TabFlow1[1][0][i] = CutFlow_mumue_FCNCzutFullSim->GetBinContent(i+1);
	TabFlow2[1][0][i] = ErrCutFlow_mumue_FCNCzutFullSim->GetBinContent(i+1);
        TabFlow1[4][0][i] += CutFlow_mumue_FCNCzutFullSim->GetBinContent(i+1);
	TabFlow2[4][0][i] += ErrCutFlow_mumue_FCNCzutFullSim->GetBinContent(i+1);
      }	
      if ( CutFlow_mumue_TTbarBkg!=NULL ) { 
        TabFlow1[1][1][i] = CutFlow_mumue_TTbarBkg->GetBinContent(i+1);
        TabFlow2[1][1][i] = ErrCutFlow_mumue_TTbarBkg->GetBinContent(i+1);
        TabFlow1[4][1][i] += CutFlow_mumue_TTbarBkg->GetBinContent(i+1);
        TabFlow2[4][1][i] += ErrCutFlow_mumue_TTbarBkg->GetBinContent(i+1);
      }	
      
      if ( CutFlow_mumue_Zjets!=NULL ) {
        TabFlow1[1][2][i] = CutFlow_mumue_Zjets->GetBinContent(i+1);
	TabFlow2[1][2][i] = ErrCutFlow_mumue_Zjets->GetBinContent(i+1);
        TabFlow1[4][2][i] += CutFlow_mumue_Zjets->GetBinContent(i+1);
	TabFlow2[4][2][i] += ErrCutFlow_mumue_Zjets->GetBinContent(i+1);
      }		 
      if ( CutFlow_mumue_DYToLL_M_10To50!=NULL ) {
	TabFlow1[1][2][i] += CutFlow_mumue_DYToLL_M_10To50->GetBinContent(i+1) ;
	TabFlow2[1][2][i] += ErrCutFlow_mumue_DYToLL_M_10To50->GetBinContent(i+1) ; 
	TabFlow1[4][2][i] += CutFlow_mumue_DYToLL_M_10To50->GetBinContent(i+1) ;
	TabFlow2[4][2][i] += ErrCutFlow_mumue_DYToLL_M_10To50->GetBinContent(i+1) ; 
      }	 
      
      
      
      if ( CutFlow_mumue_Wjets!=NULL ) { 
        TabFlow1[1][3][i] = CutFlow_mumue_Wjets->GetBinContent(i+1);
	TabFlow2[1][3][i] = ErrCutFlow_mumue_Wjets->GetBinContent(i+1);
        TabFlow1[4][3][i] += CutFlow_mumue_Wjets->GetBinContent(i+1);
	TabFlow2[4][3][i] += ErrCutFlow_mumue_Wjets->GetBinContent(i+1);
      }	
      
      
      if ( CutFlow_mumue_TtW!=NULL ) { 
        TabFlow1[1][4][i] = CutFlow_mumue_TtW->GetBinContent(i+1);
	TabFlow2[1][4][i] = ErrCutFlow_mumue_TtW->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_mumue_TtW->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_mumue_TtW->GetBinContent(i+1);
      }	
      if ( CutFlow_mumue_TbartW!=NULL ) {
        TabFlow1[1][4][i] += CutFlow_mumue_TbartW->GetBinContent(i+1);
	TabFlow2[1][4][i] += ErrCutFlow_mumue_TbartW->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_mumue_TbartW->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_mumue_TbartW->GetBinContent(i+1);
      }	
      
      
      if ( CutFlow_mumue_TtChan!=NULL ) { 
        TabFlow1[1][4][i] = CutFlow_mumue_TtChan->GetBinContent(i+1);
	TabFlow2[1][4][i] = ErrCutFlow_mumue_TtChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_mumue_TtChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_mumue_TtChan->GetBinContent(i+1);
      }	
      if ( CutFlow_mumue_TbartChan!=NULL ) {
        TabFlow1[1][4][i] += CutFlow_mumue_TbartChan->GetBinContent(i+1);
	TabFlow2[1][4][i] += ErrCutFlow_mumue_TbartChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_mumue_TbartChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_mumue_TbartChan->GetBinContent(i+1);
      }	
      
      if ( CutFlow_mumue_TsChan!=NULL ) { 
        TabFlow1[1][4][i] = CutFlow_mumue_TsChan->GetBinContent(i+1);
	TabFlow2[1][4][i] = ErrCutFlow_mumue_TsChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_mumue_TsChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_mumue_TsChan->GetBinContent(i+1);
      }	
      if ( CutFlow_mumue_TbarsChan!=NULL ) {
        TabFlow1[1][4][i] += CutFlow_mumue_TbarsChan->GetBinContent(i+1);
	TabFlow2[1][4][i] += ErrCutFlow_mumue_TbarsChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_mumue_TbarsChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_mumue_TbarsChan->GetBinContent(i+1);
      }	
      
      
      if ( CutFlow_mumue_WZ!=NULL ) {
        TabFlow1[1][5][i] = CutFlow_mumue_WZ->GetBinContent(i+1) ;
	TabFlow2[1][5][i] = ErrCutFlow_mumue_WZ->GetBinContent(i+1) ;
        TabFlow1[4][5][i] += CutFlow_mumue_WZ->GetBinContent(i+1) ;
	TabFlow2[4][5][i] += ErrCutFlow_mumue_WZ->GetBinContent(i+1) ;
      }	
      if ( CutFlow_mumue_ZZ!=NULL ) { 
        TabFlow1[1][5][i] += CutFlow_mumue_ZZ->GetBinContent(i+1) ;
	TabFlow2[1][5][i] += ErrCutFlow_mumue_ZZ->GetBinContent(i+1) ;
        TabFlow1[4][5][i] += CutFlow_mumue_ZZ->GetBinContent(i+1) ;
	TabFlow2[4][5][i] += ErrCutFlow_mumue_ZZ->GetBinContent(i+1) ;
      }	
      if ( CutFlow_mumue_WW!=NULL ) { 
        TabFlow1[1][5][i] += CutFlow_mumue_WW->GetBinContent(i+1) ;
	TabFlow2[1][5][i] += ErrCutFlow_mumue_WW->GetBinContent(i+1) ;
        TabFlow1[4][5][i] += CutFlow_mumue_WW->GetBinContent(i+1) ;
	TabFlow2[4][5][i] += ErrCutFlow_mumue_WW->GetBinContent(i+1) ;
      }	
      if ( CutFlow_mumue_TZq!=NULL ) { 
        TabFlow1[1][6][i] = CutFlow_mumue_TZq->GetBinContent(i+1);
	TabFlow2[1][6][i] = ErrCutFlow_mumue_TZq->GetBinContent(i+1);
        TabFlow1[4][6][i] += CutFlow_mumue_TZq->GetBinContent(i+1);
	TabFlow2[4][6][i] += ErrCutFlow_mumue_TZq->GetBinContent(i+1);
      }	
      
      if ( CutFlow_mumue_DataMu!=NULL ) {
        TabFlow1[1][100][i] = CutFlow_mumue_DataMu->GetBinContent(i+1);
	TabFlow2[1][100][i] = ErrCutFlow_mumue_DataMu->GetBinContent(i+1);
        TabFlow1[4][100][i] += CutFlow_mumue_DataMu->GetBinContent(i+1);
	TabFlow2[4][100][i] += ErrCutFlow_mumue_DataMu->GetBinContent(i+1);
      }	
      
      
      // eemuchannel
      if ( CutFlow_eemu_FCNCzutFullSim!=NULL ) { 
        TabFlow1[2][0][i] = CutFlow_eemu_FCNCzutFullSim->GetBinContent(i+1);
	TabFlow2[2][0][i] = ErrCutFlow_eemu_FCNCzutFullSim->GetBinContent(i+1);
        TabFlow1[4][0][i] += CutFlow_eemu_FCNCzutFullSim->GetBinContent(i+1);
	TabFlow2[4][0][i] += ErrCutFlow_eemu_FCNCzutFullSim->GetBinContent(i+1);
      }	
      if ( CutFlow_eemu_TTbarBkg!=NULL ) { 
        TabFlow1[2][1][i] = CutFlow_eemu_TTbarBkg->GetBinContent(i+1);
        TabFlow2[2][1][i] = ErrCutFlow_eemu_TTbarBkg->GetBinContent(i+1);
        TabFlow1[4][1][i] += CutFlow_eemu_TTbarBkg->GetBinContent(i+1);
        TabFlow2[4][1][i] += ErrCutFlow_eemu_TTbarBkg->GetBinContent(i+1);
      }	
      
      if ( CutFlow_eemu_Zjets!=NULL ) {
        TabFlow1[2][2][i] = CutFlow_eemu_Zjets->GetBinContent(i+1);
	TabFlow2[2][2][i] = ErrCutFlow_eemu_Zjets->GetBinContent(i+1);
        TabFlow1[4][2][i] += CutFlow_eemu_Zjets->GetBinContent(i+1);
	TabFlow2[4][2][i] += ErrCutFlow_eemu_Zjets->GetBinContent(i+1);
      }	
      if ( CutFlow_eemu_DYToLL_M_10To50!=NULL ) {
	TabFlow1[2][2][i] += CutFlow_eemu_DYToLL_M_10To50->GetBinContent(i+1) ;
	TabFlow2[2][2][i] += ErrCutFlow_eemu_DYToLL_M_10To50->GetBinContent(i+1) ; 
	TabFlow1[4][2][i] += CutFlow_eemu_DYToLL_M_10To50->GetBinContent(i+1) ;
	TabFlow2[4][2][i] += ErrCutFlow_eemu_DYToLL_M_10To50->GetBinContent(i+1) ; 
      } 
      
      if ( CutFlow_eemu_Wjets!=NULL ) { 
        TabFlow1[2][3][i] = CutFlow_eemu_Wjets->GetBinContent(i+1);
	TabFlow2[2][3][i] = ErrCutFlow_eemu_Wjets->GetBinContent(i+1);
        TabFlow1[4][3][i] += CutFlow_eemu_Wjets->GetBinContent(i+1);
	TabFlow2[4][3][i] += ErrCutFlow_eemu_Wjets->GetBinContent(i+1);
      }	
      
      
      if ( CutFlow_eemu_TtW!=NULL ) { 
        TabFlow1[2][4][i] = CutFlow_eemu_TtW->GetBinContent(i+1);
	TabFlow2[2][4][i] = ErrCutFlow_eemu_TtW->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_eemu_TtW->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_eemu_TtW->GetBinContent(i+1);
      }	
      if ( CutFlow_eemu_TbartW!=NULL ) {
        TabFlow1[2][4][i] += CutFlow_eemu_TbartW->GetBinContent(i+1);
	TabFlow2[2][4][i] += ErrCutFlow_eemu_TbartW->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_eemu_TbartW->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_eemu_TbartW->GetBinContent(i+1);
      }	
      
      if ( CutFlow_eemu_TtChan!=NULL ) { 
        TabFlow1[2][4][i] = CutFlow_eemu_TtChan->GetBinContent(i+1);
	TabFlow2[2][4][i] = ErrCutFlow_eemu_TtChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_eemu_TtChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_eemu_TtChan->GetBinContent(i+1);
      }	
      if ( CutFlow_eemu_TbartChan!=NULL ) {
        TabFlow1[2][4][i] += CutFlow_eemu_TbartChan->GetBinContent(i+1);
	TabFlow2[2][4][i] += ErrCutFlow_eemu_TbartChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_eemu_TbartChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_eemu_TbartChan->GetBinContent(i+1);
      }	
      
      if ( CutFlow_eemu_TsChan!=NULL ) { 
        TabFlow1[2][4][i] = CutFlow_eemu_TsChan->GetBinContent(i+1);
	TabFlow2[2][4][i] = ErrCutFlow_eemu_TsChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_eemu_TsChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_eemu_TsChan->GetBinContent(i+1);
      }	
      if ( CutFlow_eemu_TbarsChan!=NULL ) {
        TabFlow1[2][4][i] += CutFlow_eemu_TbarsChan->GetBinContent(i+1);
	TabFlow2[2][4][i] += ErrCutFlow_eemu_TbarsChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_eemu_TbarsChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_eemu_TbarsChan->GetBinContent(i+1);
      }	
      
      
      
      
      if ( CutFlow_eemu_WZ!=NULL ) {
        TabFlow1[2][5][i] = CutFlow_eemu_WZ->GetBinContent(i+1) ;
	TabFlow2[2][5][i] = ErrCutFlow_eemu_WZ->GetBinContent(i+1) ;
        TabFlow1[4][5][i] += CutFlow_eemu_WZ->GetBinContent(i+1) ;
	TabFlow2[4][5][i] += ErrCutFlow_eemu_WZ->GetBinContent(i+1) ;
      }	
      if ( CutFlow_eemu_ZZ!=NULL ) { 
        TabFlow1[2][5][i] += CutFlow_eemu_ZZ->GetBinContent(i+1) ;
	TabFlow2[2][5][i] += ErrCutFlow_eemu_ZZ->GetBinContent(i+1) ;
        TabFlow1[4][5][i] += CutFlow_eemu_ZZ->GetBinContent(i+1) ;
	TabFlow2[4][5][i] += ErrCutFlow_eemu_ZZ->GetBinContent(i+1) ;
      }	
      if ( CutFlow_eemu_WW!=NULL ) { 
        TabFlow1[2][5][i] += CutFlow_eemu_WW->GetBinContent(i+1) ;
	TabFlow2[2][5][i] += ErrCutFlow_eemu_WW->GetBinContent(i+1) ;
        TabFlow1[4][5][i] += CutFlow_eemu_WW->GetBinContent(i+1) ;
	TabFlow2[4][5][i] += ErrCutFlow_eemu_WW->GetBinContent(i+1) ;
      }	
      
      if ( CutFlow_eemu_TZq!=NULL ) { 
        TabFlow1[2][6][i] = CutFlow_eemu_TZq->GetBinContent(i+1);
	TabFlow2[2][6][i] = ErrCutFlow_eemu_TZq->GetBinContent(i+1);
        TabFlow1[4][6][i] += CutFlow_eemu_TZq->GetBinContent(i+1);
	TabFlow2[4][6][i] += ErrCutFlow_eemu_TZq->GetBinContent(i+1);
      }	
      if ( CutFlow_eemu_DataEG!=NULL ) {
        TabFlow1[2][100][i] = CutFlow_eemu_DataEG->GetBinContent(i+1);
	TabFlow2[2][100][i] = ErrCutFlow_eemu_DataEG->GetBinContent(i+1);
        TabFlow1[4][100][i] += CutFlow_eemu_DataEG->GetBinContent(i+1);
	TabFlow2[4][100][i] += ErrCutFlow_eemu_DataEG->GetBinContent(i+1);
      }	
      
      // eee channel  
      if ( CutFlow_eee_FCNCzutFullSim!=NULL ) { 
        TabFlow1[3][0][i] = CutFlow_eee_FCNCzutFullSim->GetBinContent(i+1);
	TabFlow2[3][0][i] = ErrCutFlow_eee_FCNCzutFullSim->GetBinContent(i+1);
        TabFlow1[4][0][i] += CutFlow_eee_FCNCzutFullSim->GetBinContent(i+1);
	TabFlow2[4][0][i] += ErrCutFlow_eee_FCNCzutFullSim->GetBinContent(i+1);
      }	
      if ( CutFlow_eee_TTbarBkg!=NULL ) { 
        TabFlow1[3][1][i] = CutFlow_eee_TTbarBkg->GetBinContent(i+1);
        TabFlow2[3][1][i] = ErrCutFlow_eee_TTbarBkg->GetBinContent(i+1);
        TabFlow1[4][1][i] += CutFlow_eee_TTbarBkg->GetBinContent(i+1);
        TabFlow2[4][1][i] += ErrCutFlow_eee_TTbarBkg->GetBinContent(i+1);
      }	
      
      if ( CutFlow_eee_Zjets!=NULL ) {
        TabFlow1[3][2][i] = CutFlow_eee_Zjets->GetBinContent(i+1);
	TabFlow2[3][2][i] = ErrCutFlow_eee_Zjets->GetBinContent(i+1);
        TabFlow1[4][2][i] += CutFlow_eee_Zjets->GetBinContent(i+1);
	TabFlow2[4][2][i] += ErrCutFlow_eee_Zjets->GetBinContent(i+1);
      }	
      if ( CutFlow_eee_DYToLL_M_10To50!=NULL ) {
	TabFlow1[3][2][i] += CutFlow_eee_DYToLL_M_10To50->GetBinContent(i+1) ;
	TabFlow2[3][2][i] += ErrCutFlow_eee_DYToLL_M_10To50->GetBinContent(i+1) ; 
	TabFlow1[4][2][i] += CutFlow_eee_DYToLL_M_10To50->GetBinContent(i+1) ;
	TabFlow2[4][2][i] += ErrCutFlow_eee_DYToLL_M_10To50->GetBinContent(i+1) ; 
      }	
      if ( CutFlow_eee_Wjets!=NULL ) { 
        TabFlow1[3][3][i] = CutFlow_eee_Wjets->GetBinContent(i+1);
	TabFlow2[3][3][i] = ErrCutFlow_eee_Wjets->GetBinContent(i+1);
        TabFlow1[4][3][i] += CutFlow_eee_Wjets->GetBinContent(i+1);
	TabFlow2[4][3][i] += ErrCutFlow_eee_Wjets->GetBinContent(i+1);
      }	
      
      
      if ( CutFlow_eee_TtW!=NULL ) { 
        TabFlow1[3][4][i] = CutFlow_eee_TtW->GetBinContent(i+1);
	TabFlow2[3][4][i] = ErrCutFlow_eee_TtW->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_eee_TtW->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_eee_TtW->GetBinContent(i+1);
      }	
      if ( CutFlow_eee_TbartW!=NULL ) {
        TabFlow1[3][4][i] += CutFlow_eee_TbartW->GetBinContent(i+1);
	TabFlow2[3][4][i] += ErrCutFlow_eee_TbartW->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_eee_TbartW->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_eee_TbartW->GetBinContent(i+1);
      }	
      
      if ( CutFlow_eee_TtChan!=NULL ) { 
        TabFlow1[3][4][i] = CutFlow_eee_TtChan->GetBinContent(i+1);
	TabFlow2[3][4][i] = ErrCutFlow_eee_TtChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_eee_TtChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_eee_TtChan->GetBinContent(i+1);
      }	
      if ( CutFlow_eee_TbartChan!=NULL ) {
        TabFlow1[3][4][i] += CutFlow_eee_TbartChan->GetBinContent(i+1);
	TabFlow2[3][4][i] += ErrCutFlow_eee_TbartChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_eee_TbartChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_eee_TbartChan->GetBinContent(i+1);
      }	
      
      if ( CutFlow_eee_TsChan!=NULL ) { 
        TabFlow1[3][4][i] = CutFlow_eee_TsChan->GetBinContent(i+1);
	TabFlow2[3][4][i] = ErrCutFlow_eee_TsChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_eee_TsChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_eee_TsChan->GetBinContent(i+1);
      }	
      if ( CutFlow_eee_TbarsChan!=NULL ) {
        TabFlow1[3][4][i] += CutFlow_eee_TbarsChan->GetBinContent(i+1);
	TabFlow2[3][4][i] += ErrCutFlow_eee_TbarsChan->GetBinContent(i+1);
        TabFlow1[4][4][i] += CutFlow_eee_TbarsChan->GetBinContent(i+1);
	TabFlow2[4][4][i] += ErrCutFlow_eee_TbarsChan->GetBinContent(i+1);
      }	
      
      
      
      if ( CutFlow_eee_WZ!=NULL ) {
        TabFlow1[3][5][i] = CutFlow_eee_WZ->GetBinContent(i+1) ;
	TabFlow2[3][5][i] = ErrCutFlow_eee_WZ->GetBinContent(i+1) ;
        TabFlow1[4][5][i] += CutFlow_eee_WZ->GetBinContent(i+1) ;
	TabFlow2[4][5][i] += ErrCutFlow_eee_WZ->GetBinContent(i+1) ;
      }	
      if ( CutFlow_eee_ZZ!=NULL ) { 
        TabFlow1[3][5][i] += CutFlow_eee_ZZ->GetBinContent(i+1) ;
	TabFlow2[3][5][i] += ErrCutFlow_eee_ZZ->GetBinContent(i+1) ;
        TabFlow1[4][5][i] += CutFlow_eee_ZZ->GetBinContent(i+1) ;
	TabFlow2[4][5][i] += ErrCutFlow_eee_ZZ->GetBinContent(i+1) ;
      }	
      if ( CutFlow_eee_WW!=NULL ) { 
        TabFlow1[3][5][i] += CutFlow_eee_WW->GetBinContent(i+1) ;
	TabFlow2[3][5][i] += ErrCutFlow_eee_WW->GetBinContent(i+1) ;
        TabFlow1[4][5][i] += CutFlow_eee_WW->GetBinContent(i+1) ;
	TabFlow2[4][5][i] += ErrCutFlow_eee_WW->GetBinContent(i+1) ;
      }	
      if ( CutFlow_eee_TZq!=NULL ) { 
        TabFlow1[3][6][i] = CutFlow_eee_TZq->GetBinContent(i+1);
	TabFlow2[3][6][i] = ErrCutFlow_eee_TZq->GetBinContent(i+1);
        TabFlow1[4][6][i] += CutFlow_eee_TZq->GetBinContent(i+1);
	TabFlow2[4][6][i] += ErrCutFlow_eee_TZq->GetBinContent(i+1);
      }	
      
      if ( CutFlow_eee_DataEG!=NULL ) {
        TabFlow1[3][100][i] = CutFlow_eee_DataEG->GetBinContent(i+1);
	TabFlow2[3][100][i] = ErrCutFlow_eee_DataEG->GetBinContent(i+1);
        TabFlow1[4][100][i] += CutFlow_eee_DataEG->GetBinContent(i+1);
	TabFlow2[4][100][i] += ErrCutFlow_eee_DataEG->GetBinContent(i+1);
      }	
      
      
      
    }
    
    
    cout<<"#########################"<<endl;
    cout<<" Fill the latex tables   "<<endl;
    cout<<"#########################"<<endl;
    
    //  string ofilenametex = string("CrossSection")+string("_")+ChannelName+string(".tex");
    //  string ofilenametex = string("CrossSection")+string(".tex");
    
    string ofilenametex = "CutFlow.tex";
    ofstream ofile(ofilenametex.c_str());
    
    ofile<<"\\documentclass[amsmath,amssymb]{revtex4}"<<endl;
    ofile<<"\\begin{document}"<<endl;
    
    ofile.setf(ios::fixed);
    ofile.precision(1);
    
    vector<string> CutName;
    CutName.push_back("Trigger");
    CutName.push_back("lepton selection");
    CutName.push_back("Z mass cut");
    CutName.push_back("jet multi.");
    CutName.push_back("b-jet multi. ");
    CutName.push_back("Mt cut ");
    
    //compute error
    for(int k0=0; k0<5; ++k0) {
      for(int k1=0; k1<101; ++k1) {
	for(int k2=0; k2<101; ++k2) {
	  //if data driven, don't usr square root
	  /*if( (k1 == 2 && k2 >= 4) || k1 == 1 ) TabFlow2[k0][k1][k2]= TabFlow2[k0][k1][k2];
	  else{
	    if (k1==0 || k1==2 || k1==4 || k1==5   ){       
	      TabFlow2[k0][k1][k2]= TabFlow2[k0][k1][k2] 
		+ pow(0.045*TabFlow1[k0][k1][k2], 2) 
		+ pow(TriggError[k0]*TabFlow1[k0][k1][k2], 2) 
		+ pow(SF_Lepton_error[k0]  *TabFlow1[k0][k1][k2], 2) ;
	      if(k2 >= 5)  TabFlow2[k0][k1][k2] +=  pow(SF_MET_error[k0], 2) ;
	      TabFlow2[k0][k1][k2]  = sqrt(TabFlow2[k0][k1][k2]);
	    }
	    else{ TabFlow2[k0][k1][k2]= sqrt(TabFlow2[k0][k1][k2] );}
	    }*/
	  TabFlow2[k0][k1][k2]= sqrt(TabFlow2[k0][k1][k2]);
	}
      }
    }  
    
    
    //**************************
    //Compute total bckgd
    for(int k0=0; k0<6; ++k0) {
      for(int k1=1; k1<10; ++k1) {
	for(int k2=0; k2<50; ++k2) {
	  TabFlow1[k0][50][k2] += TabFlow1[k0][k1][k2];
	}
      }
    } 
    
    //*********************************
    //Compute error for total bckgd
    for(int k0=0; k0<6; ++k0) {
      for(int k2=1; k2<51; ++k2) {
	for(int k1=1; k1<10; ++k1) {
	  TabFlow2[k0][50][k2] += TabFlow2[k0][k1][k2]*TabFlow2[k0][k1][k2];
	}
	TabFlow2[k0][50][k2] = sqrt(TabFlow2[k0][50][k2]);
      }
    } 
    
     
   
    
    for (int IChannel=0; IChannel<5; IChannel++) {
      
      // Summary tables
      ofile << "\\clearpage" << endl;
      ofile << "\\begin{landscape}" << endl;
      ofile << "\\begin{table}[p]" << endl;
      
      ofile << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;
      ofile << "\\hline" << endl;
      ofile << "\\hline" << endl;
      ofile << "Cut & DATA & Sum MC & signal  & Total Background & S/B \\\\" << endl;
      ofile << "\\hline" << endl;
      
      for(int ic=0; ic<CutName.size(); ++ic) {
	
	  ofile <<CutName[ic]<<" & "<<  TabFlow1[IChannel][100][ic] << " $\\pm$ "<<  TabFlow2[IChannel][100][ic] << " & " <<
	    TabFlow1[IChannel][0][ic]+TabFlow1[IChannel][50][ic] << " $\\pm$ "<< 
	    sqrt(TabFlow2[IChannel][0][ic]*TabFlow2[IChannel][0][ic]+TabFlow2[IChannel][50][ic]*TabFlow2[IChannel][50][ic]) << " & " <<
	    TabFlow1[IChannel][0][ic] << " $\\pm$ "<< TabFlow2[IChannel][0][ic] << " & " <<
	    TabFlow1[IChannel][50][ic] << " $\\pm$ "<< TabFlow2[IChannel][50][ic] << " & " <<
	    TabFlow1[IChannel][0][ic]/TabFlow1[IChannel][50][ic] <<  " \\\\" << endl;
      }
      
      // summary table
      ofile.precision(1);
      ofile << "\\hline" << endl;
      // ajout Caro
      // ofile << "$SF_{trig}$  & " << TabFlow1[IChannel][100][10] <<  " & " <<
      //     sf_trig_sum[IChannel] <<  " & " <<
      //     sf_trig_sum[IChannel] << " & " <<
      //     sf_trig_sum[IChannel] <<  " &  \\\\" << endl;
      // ofile << "$SF_{ID}$  & " << TabFlow1[IChannel][100][11] << " & " <<
      //     sf_id_sum[IChannel] << " & " <<
      //     TabFlow1[IChannel][0][11]   << " $\\pm $ " << TabFlow2[IChannel][0][11]  << " & " <<
      //     TabFlow1[IChannel][50][11]  << " $\\pm $ " << TabFlow2[IChannel][50][11] << " &  \\\\" << endl;
      // ofile << "\\hline" << endl;
      ofile << "\\hline" << endl;
      ofile << "\\end{tabular}" << endl;
      ofile << " " << endl;
      // end ajout Caro
      ofile << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;
      //ofile << "\\begin{tabular}{|l|c|c|c|c|c|c|c|}" << endl;
      ofile << "\\hline" << endl;
      ofile << "\\hline" << endl;
      ofile << "Cut & ttbar & SingleTop & DY & Diboson & tZq \\\\" << endl;
      //ofile << "Cut & TopBackg & SingleTop & DY & Wjets  & Diboson & QCD  \\\\" << endl;
      ofile << "\\hline" << endl;
      
      for(int ic=0; ic<CutName.size(); ++ic) {
	double TopBkgd    = TabFlow1[IChannel][1][ic];
	double ErrTopBkgd = TabFlow2[IChannel][1][ic];
	// to be changed
	//    double DY = DY + mll<50;
	double DY         = TabFlow1[IChannel][2][ic];
	double errDY      = TabFlow2[IChannel][2][ic];
	//double errDY      = (TabFlow2[IChannel][2][ic]*TabFlow2[IChannel][2][ic]);
	//errDY             = sqrt(errDY);
	double Wj         = TabFlow1[IChannel][3][ic];
	double errWjets   = TabFlow2[IChannel][3][ic];
	double Qcd        = TabFlow1[IChannel][8][ic];
	double errQcd     = TabFlow2[IChannel][8][ic];
	double SinglTop    = TabFlow1[IChannel][4][ic];
	double errSinglTop = sqrt(TabFlow2[IChannel][4][ic]*TabFlow2[IChannel][4][ic]);
	double Dibos      =  TabFlow1[IChannel][5][ic];
	double errDibos   = TabFlow2[IChannel][5][ic];
	double tzq      =  TabFlow1[IChannel][6][ic];
	double errtzq   = TabFlow2[IChannel][6][ic];
	
	ofile.precision(1);
	ofile <<CutName[ic]<<" & "<<
	  TopBkgd         << " $\\pm$ "<< ErrTopBkgd      << " & " <<
	  SinglTop        << " $\\pm$ "<< errSinglTop     << " & " <<
	  DY              << " $\\pm$ "<< errDY           << " & " <<
	  //Wj              << " $\\pm$ "<< errWjets        << " & " <<
	  //Dibos           << " $\\pm$ "<< errDibos        << " & " <<
	  Dibos           << " $\\pm$ "<< errDibos        << " &" << 
	  tzq           << " $\\pm$ "<< errtzq       << " \\\\" << endl;
	//Qcd             << " $\\pm$ "<< errQcd <<   " \\\\" << endl;
      }
      
      ofile.precision(1);
      // mommentane!
      ofile << "\\hline " << endl;
      ofile << "\\hline" << endl;
      ofile << "\\end{tabular}" << endl;
      
      string ChannelName;
      if (IChannel==0) ChannelName= "mumumu"; 
      else if (IChannel==1) ChannelName= "mumue"; 
      else if (IChannel==2) ChannelName= "eemu"; 
      else if (IChannel==3) ChannelName= "eee"; 
      
      if ( ChannelName == "mumumu" )  ofile << "\\caption{tri-muon cut flow, including SF for trigger efficiency, lepton and \\met selection. For the \\ttbar, diboson and single top samples, the uncertainties account for the uncertainty on the luminosity and the various selection efficiencies. The DY backgrounds is estimated from data.}" << endl;
      if ( ChannelName == "mumue" )  ofile << "\\caption{di-muon+electron cut flow, including SF for trigger efficiency, lepton and \\met selection. For the \\ttbar, diboson and single top samples, the uncertainties account for the uncertainty on the luminosity and the various selection efficiencies. The DY backgrounds is estimated from data.}" << endl;
      if ( ChannelName == "eemu" )  ofile << "\\caption{ di-electron+muon cut flow, including SF for trigger efficiency, lepton and \\met selection. For the \\ttbar, diboson and single top samples, the uncertainties account for the uncertainty on the luminosity and the various selection efficiencies. The DY backgrounds is estimated from data.}" << endl;
      if ( ChannelName == "eee" )  ofile << "\\caption{ tri electron cut flow, including SF for trigger efficiency, lepton and \\met selection. For the \\ttbar, diboson and single top samples, the uncertainties account for the uncertainty on the luminosity and the various selection efficiencies. The DY backgrounds is estimated from data.The DY backgrounds is estimated from data.  }" << endl;
      if ( ChannelName == "mumumu" ) ofile << "\\label{Table:CutFlow_mumumu}" << endl;
      if ( ChannelName == "mumue" ) ofile << "\\label{Table:CutFlow_mumue}" << endl;
      if ( ChannelName == "eemu"   ) ofile << "\\label{Table:CutFlow_eemu}" << endl;
      if ( ChannelName == "eee"  ) ofile << "\\label{Table:CutFlow_eee}" << endl;
      
      
      
      ofile << "\\end{table}" << endl;
      ofile << "\\end{landscape}" << endl;
    } // end loop IChannel
    
    
    
    ofile<<"\\end{document}"<<endl;
    //  string prodpdf = string("pdflatex CrossSection")+string("_")+ChannelName+string(".tex");
    //  string prodpdf = string("pdflatex CrossSection")+string(".tex");
    string prodpdf = string("pdflatex ")+ofilenametex;
    system(prodpdf.c_str());
    
  }
}
  
