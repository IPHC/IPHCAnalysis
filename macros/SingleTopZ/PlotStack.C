#include "TString.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "THStack.h"


#include <iostream>








void PlotStack(TString plotname, TString namechan, TString selection, bool setlogy, bool ratio){
  
  float lumi_error=0.022;
  
  
  float SF_trigger_error = 0.015;
  float SF_Lepton_error  = 0.010;
  float SF_MET_error     = 0.016;
 
  
  double frac_ee   = 0;
  double frac_mumu = 0;
  double frac_emu  = 0;
  
  
  

  TString filename;
  //filename="backup_outputProof20-12-12_10-03-23_noSFapplied_Zenriched/proof_merged_noSF_Zenriched.root"; 
  filename="backup_outputProof04-06-13_20-20/proof_merged.root";
  
  
  TString channel;
  if (namechan!="all") channel=namechan;
  else channel="mumumu";
  
  Int_t stati=0;
  Bool_t  fit=1;
  Bool_t logy=0;
  bool readStop = false;
  /*
    gStyle->SetPadRightMargin(0.13);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadBottomMargin(0.13);
  */
  /*gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    
    gStyle->SetOptDate(0);
    gStyle->SetStatColor(0);
    gStyle->SetTitleColor(1);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    
    // For the axis titles:
    
    gStyle->SetTitleColor(1, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTitleSize(0.06, "XYZ");
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYOffset(1.25);
    
    // For the axis labels:
    
    gStyle->SetLabelColor(1, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.007, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");
    
    // For the axis:
    
    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetNdivisions(510, "XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);
  */
  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0); // must be kWhite but I dunno how to do that in PyROOT
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(600); //Width of canvas
  gStyle->SetCanvasDefX(0);   //POsition on screen
  gStyle->SetCanvasDefY(0);
  
  
  // For the Pad:
  gStyle->SetPadBorderMode(0);
  // ROOT . gStyle . SetPadBorderSize(Width_t size = 1);
  gStyle->SetPadColor(0); // kWhite
  gStyle->SetPadGridX(0); //false
  gStyle->SetPadGridY(0); //false
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  
  // For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);
  
  // For the histo:
  // ROOT . gStyle . SetHistFillColor(1);
  // ROOT . gStyle . SetHistFillStyle(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);
  // ROOT . gStyle . SetLegoInnerR(Float_t rad = 0.5);
  // ROOT . gStyle . SetNumberContours(Int_t number = 20);
  
  gStyle->SetEndErrorSize(2);
  //ROOT . gStyle . SetErrorMarker(20);   /// I COMMENTED THIS OUT
  //ROOT . gStyle . SetErrorX(0.);
  
  //ROOT . gStyle . SetMarkerStyle(20);
  
  
  //For the fit/function:
  gStyle->SetOptFit(1011);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);
  
  //For the date:
  gStyle->SetOptDate(0);
  // ROOT . gStyle . SetDateX(Float_t x = 0.01);
  // ROOT . gStyle . SetDateY(Float_t y = 0.01);
  
  // For the statistics box:
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  gStyle->SetStatColor(0); // kWhite
  gStyle->SetStatFont(42);
  //ROOT . gStyle . SetStatFontSize(0.025);
  gStyle->SetStatFontSize(0.04);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.1);
  gStyle->SetStatW(0.15);
  // ROOT . gStyle . SetStatStyle(Style_t style = 1001);
  // ROOT . gStyle . SetStatX(Float_t x = 0);
  // ROOT . gStyle . SetStatY(Float_t y = 0);
  
  // Margins:
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  //ROOT . gStyle . SetPadRightMargin(0.12);
  gStyle->SetPadRightMargin(0.03);
  
  // For the Global title:
  
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
  // ROOT . gStyle . SetTitleH(0); // Set the height of the title box
  // ROOT . gStyle . SetTitleW(0); // Set the width of the title box
  // ROOT . gStyle . SetTitleX(0); // Set the position of the title box
  // ROOT . gStyle . SetTitleY(0.985); // Set the position of the title box
  // ROOT . gStyle . SetTitleStyle(Style_t style = 1001);
  // ROOT . gStyle . SetTitleBorderSize(2);
  
  // For the axis titles:
  
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  // ROOT . gStyle . SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // ROOT . gStyle . SetTitleYSize(Float_t size = 0.02);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.25);
  // ROOT . gStyle . SetTitleOffset(1.1, "Y"); // Another way to set the Offset
  
  // For the axis labels:
  
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");
  
  // For the axis:
  
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(1); // kTRUE
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
  
  // Change for log plots:
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);
  
  // Postscript options:
  gStyle->SetPaperSize(20.,20.);
  // ROOT . gStyle . SetLineScalePS(Float_t scale = 3);
  // ROOT . gStyle . SetLineStyleString(Int_t i, const char* text);
  // ROOT . gStyle . SetHeaderPS(const char* header);
  // ROOT . gStyle . SetTitlePS(const char* pstitle);
  
  // ROOT . gStyle . SetBarOffset(Float_t baroff = 0.5);
  // ROOT . gStyle . SetBarWidth(Float_t barwidth = 0.5);
  // ROOT . gStyle . SetPaintTextFormat(const char* format = "g");
  // ROOT . gStyle . SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // ROOT . gStyle . SetTimeOffset(Double_t toffset);
  // ROOT . gStyle . SetHistMinimumZero(kTRUE);
  
  
  //TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
  TCanvas *c1 = new TCanvas("c1","c1", 1000, 800);
  c1->SetBottomMargin(0.3);
  c1->SetLogy(setlogy);
  c1->cd();
  
  /*TPad *canvas_1 = new TPad("canvas_1", "canvas_1",0,0.25,1.0,0.98);
    if (ratio) {
    canvas_1->Draw();
    canvas_1->cd(); 
    canvas_1->SetLogy(setlogy);
    }
    else {
    c1->SetLogy(setlogy);
    }*/
 
 
 
  TFile * filechannel = new TFile(filename);

  TH1F * histo_Data;
  TH1F * histo_DataEG;
  TH1F * histo_DataMu;
  TH1F * histo_FCNCzut;
  TH1F * histo_TTbarBkg;
  TH1F * histo_TTbarSig;
  TH1F * histo_Zjets;
  TH1F * histo_DYToLL_M10_50;
  TH1F * histo_Wjets;
  TH1F * histo_VV ;
  TH1F * histo_WW ;
  TH1F * histo_WZ ;
  TH1F * histo_ZZ ;
  TH1F * histo_TZq;
  TH1F * histo_SingleToptW;
  TH1F * histo_TtW;
  TH1F * histo_TbartW;
  TH1F * histo_SingleToptChan;
  TH1F * histo_TtChan;
  TH1F * histo_TbartChan;
  TH1F * histo_SingleTopsChan;
  TH1F * histo_TsChan;
  TH1F * histo_TbarsChan;
  TH1F * histo_ratio;
  
  
  //cout << "plotting "  << plotname+channel+selection << endl;
  TString histo_Data_name_DataEG   = plotname+channel+selection+"_DataEG";
  TString histo_Data_name_DataMu   = plotname+channel+selection+"_DataMu";
  TString histo_Data_name_DataMuEG = plotname+channel+selection+"_DataMuEG";
  
  
  TH1F * histo_Data_DataEG   = (TH1F*)filechannel->Get(histo_Data_name_DataEG);
  TH1F * histo_Data_DataMu   = (TH1F*)filechannel->Get(histo_Data_name_DataMu);
  TH1F * histo_Data_DataMuEG = (TH1F*)filechannel->Get(histo_Data_name_DataMuEG);
  
 
  if(namechan=="mumumu"){
    histo_Data  = (TH1F*) histo_Data_DataMu->Clone()  ;
  }
  
  if(namechan=="eee"){
    histo_Data  = (TH1F*) histo_Data_DataEG->Clone()  ;
  }
  
  if(namechan=="mumue"){
    histo_Data  = (TH1F*) histo_Data_DataMu->Clone()  ;
    histo_Data->Add(histo_Data, histo_Data_DataMuEG, 1, 1);
  }
  
  if(namechan=="eemu"){
    histo_Data  = (TH1F*) histo_Data_DataEG->Clone() ;
    histo_Data->Add(histo_Data, histo_Data_DataMuEG, 1, 1);
  }
  
  
  
  
  if (namechan=="all") {
  
     TString histo_Data_name_mumumu_DataEG   = plotname+"mumumu"+selection+"_DataEG";
     TString histo_Data_name_mumumu_DataMu   = plotname+"mumumu" +selection+"_DataMu";
     TString histo_Data_name_mumumu_DataMuEG = plotname+"mumumu" +selection+"_DataMuEG";
  
     TString histo_Data_name_mumue_DataEG   = plotname+"mumue"+selection+"_DataEG";
     TString histo_Data_name_mumue_DataMu   = plotname+"mumue" +selection+"_DataMu";
     TString histo_Data_name_mumue_DataMuEG = plotname+"mumue" +selection+"_DataMuEG";
     
     TString histo_Data_name_eemu_DataEG   = plotname+"eemu"+selection+"_DataEG";
     TString histo_Data_name_eemu_DataMu   = plotname+"eemu" +selection+"_DataMu";
     TString histo_Data_name_eemu_DataMuEG = plotname+"eemu" +selection+"_DataMuEG";
     
     TString histo_Data_name_eee_DataEG   = plotname+"eee"+selection+"_DataEG";
     TString histo_Data_name_eee_DataMu   = plotname+"eee" +selection+"_DataMu";
     TString histo_Data_name_eee_DataMuEG = plotname+"eee" +selection+"_DataMuEG";
     
     
     TH1F * histo_tmp;
     
     histo_Data  = (TH1F*)filechannel->Get(histo_Data_name_mumumu_DataEG);
     histo_tmp = (TH1F*)filechannel->Get(histo_Data_name_mumumu_DataMu);
     histo_Data->Add(histo_Data, histo_tmp);
     histo_tmp = (TH1F*)filechannel->Get(histo_Data_name_mumumu_DataMuEG);
     histo_Data->Add(histo_Data, histo_tmp);
     
     
     histo_tmp  = (TH1F*)filechannel->Get(histo_Data_name_mumue_DataEG);
     histo_Data->Add(histo_Data, histo_tmp);
     histo_tmp = (TH1F*)filechannel->Get(histo_Data_name_mumue_DataMu);
     histo_Data->Add(histo_Data, histo_tmp);
     histo_tmp = (TH1F*)filechannel->Get(histo_Data_name_mumue_DataMuEG);
     histo_Data->Add(histo_Data, histo_tmp);
     
     histo_tmp  = (TH1F*)filechannel->Get(histo_Data_name_eemu_DataEG);
     histo_Data->Add(histo_Data, histo_tmp);
     histo_tmp = (TH1F*)filechannel->Get(histo_Data_name_eemu_DataMu);
     histo_Data->Add(histo_Data, histo_tmp);
     histo_tmp = (TH1F*)filechannel->Get(histo_Data_name_eemu_DataMuEG);
     histo_Data->Add(histo_Data, histo_tmp);
     
     histo_tmp  = (TH1F*)filechannel->Get(histo_Data_name_eee_DataEG);
     histo_Data->Add(histo_Data, histo_tmp);
     histo_tmp = (TH1F*)filechannel->Get(histo_Data_name_eee_DataMu);
     histo_Data->Add(histo_Data, histo_tmp);
     histo_tmp = (TH1F*)filechannel->Get(histo_Data_name_eee_DataMuEG);
     histo_Data->Add(histo_Data, histo_tmp);
  
  }
  
  
 
 
  TString histo_FCNCzut_name= plotname+channel+selection+"_FCNCzut";
  histo_FCNCzut            = (TH1F*)filechannel->Get(histo_FCNCzut_name);
  histo_FCNCzut->Scale(0.1);
  if (namechan=="all") {
    TString h_name2_FCNCzut=plotname+"mumue"+selection+"_FCNCzut";
    TH1F * h_2_FCNCzut = (TH1F*)filechannel->Get(h_name2_FCNCzut);
    h_2_FCNCzut->Scale(0.1);
    TString h_name3_FCNCzut=plotname+"eemu"+selection+"_FCNCzut";
    TH1F * h_3_FCNCzut = (TH1F*)filechannel->Get(h_name3_FCNCzut);
    h_3_FCNCzut->Scale(0.1);
    TString h_name4_FCNCzut=plotname+"eee"+selection+"_FCNCzut";
    TH1F * h_4_FCNCzut = (TH1F*)filechannel->Get(h_name4_FCNCzut);
    h_4_FCNCzut->Scale(0.1);
    histo_FCNCzut->Add(histo_FCNCzut, h_2_FCNCzut     , 1, 1);
    histo_FCNCzut->Add(histo_FCNCzut, h_3_FCNCzut     , 1, 1);
    histo_FCNCzut->Add(histo_FCNCzut, h_4_FCNCzut     , 1, 1);
  }


  TString histo_TTbarBkg_name= plotname+channel+selection+"_TTbarBkg";
  histo_TTbarBkg            = (TH1F*)filechannel->Get(histo_TTbarBkg_name);
  if (namechan=="all") {
    TString h_name2_TTbarBkg=plotname+"mumue"+selection+"_TTbarBkg";
    TH1F * h_2_TTbarBkg = (TH1F*)filechannel->Get(h_name2_TTbarBkg);
    TString h_name3_TTbarBkg=plotname+"eemu"+selection+"_TTbarBkg";
    TH1F * h_3_TTbarBkg = (TH1F*)filechannel->Get(h_name3_TTbarBkg);
    TString h_name4_TTbarBkg=plotname+"eee"+selection+"_TTbarBkg";
    TH1F * h_4_TTbarBkg = (TH1F*)filechannel->Get(h_name4_TTbarBkg);
    histo_TTbarBkg->Add(histo_TTbarBkg, h_2_TTbarBkg     , 1, 1);
    histo_TTbarBkg->Add(histo_TTbarBkg, h_3_TTbarBkg     , 1, 1);
    histo_TTbarBkg->Add(histo_TTbarBkg, h_4_TTbarBkg     , 1, 1);
  }


  TString histo_TTbarSig_name= plotname+channel+selection+"_TTbarSig";
  histo_TTbarSig            = (TH1F*)filechannel->Get(histo_TTbarSig_name);
  if (namechan=="all") {
    TString h_name2_TTbarSig=plotname+"mumue"+selection+"_TTbarSig";
    TH1F * h_2_TTbarSig = (TH1F*)filechannel->Get(h_name2_TTbarSig);
    TString h_name3_TTbarSig=plotname+"eemu"+selection+"_TTbarSig";
    TH1F * h_3_TTbarSig = (TH1F*)filechannel->Get(h_name3_TTbarSig);
    TString h_name4_TTbarSig=plotname+"eee"+selection+"_TTbarSig";
    TH1F * h_4_TTbarSig = (TH1F*)filechannel->Get(h_name4_TTbarSig);
    histo_TTbarSig->Add(histo_TTbarSig, h_2_TTbarSig     , 1, 1);
    histo_TTbarSig->Add(histo_TTbarSig, h_3_TTbarSig     , 1, 1);
    histo_TTbarSig->Add(histo_TTbarSig, h_4_TTbarSig     , 1, 1);
  }
  
  
  TString histo_Zjets_name= plotname+channel+selection+"_Zjets";
  histo_Zjets               = (TH1F*)filechannel->Get(histo_Zjets_name);
  //cout << "histo_Zjets_name " << histo_Zjets_name << endl;
  if (namechan=="all") {
    TString h_name2_Zjets=plotname+"mumue"+selection+"_Zjets";
    TH1F * h_2_Zjets = (TH1F*)filechannel->Get(h_name2_Zjets);
    TString h_name3_Zjets=plotname+"eemu"+selection+"_Zjets";
    TH1F * h_3_Zjets = (TH1F*)filechannel->Get(h_name3_Zjets);
    TString h_name4_Zjets=plotname+"eee"+selection+"_Zjets";
    TH1F * h_4_Zjets = (TH1F*)filechannel->Get(h_name4_Zjets);
    frac_ee   = histo_Zjets->Integral();
    frac_mumu = h_2_Zjets->Integral();
    frac_emu  = h_3_Zjets->Integral();
    histo_Zjets->Add(histo_Zjets, h_2_Zjets     , 1, 1);
    histo_Zjets->Add(histo_Zjets, h_3_Zjets     , 1, 1);
    histo_Zjets->Add(histo_Zjets, h_4_Zjets     , 1, 1);
  }
  

  TString histo_DYToLL_M10_50_name = plotname+channel+selection+"_DYToLL_M10-50";
  histo_DYToLL_M10_50       = (TH1F*)filechannel->Get(histo_DYToLL_M10_50_name);
  if (namechan=="all") {
    TString h_name2_DYToLL_M10_50 =plotname+"mumue"+selection+"_DYToLL_M10-50";
    TH1F * h_2_DYToLL_M10_50      = (TH1F*)filechannel->Get(h_name2_DYToLL_M10_50);
    TString h_name3_DYToLL_M10_50 =plotname+"eemu"+selection+"_DYToLL_M10-50";
    TH1F * h_3_DYToLL_M10_50      = (TH1F*)filechannel->Get(h_name3_DYToLL_M10_50);
    frac_ee   += histo_DYToLL_M10_50->Integral();
    frac_mumu += h_2_DYToLL_M10_50->Integral();
    frac_emu  += h_3_DYToLL_M10_50->Integral();
    histo_DYToLL_M10_50->Add(histo_DYToLL_M10_50, h_2_DYToLL_M10_50     , 1, 1);
    histo_DYToLL_M10_50->Add(histo_DYToLL_M10_50, h_3_DYToLL_M10_50     , 1, 1);
  }
  

  TString histo_Wjets_name = plotname+channel+selection+"_Wjets";
  histo_Wjets               = (TH1F*)filechannel->Get(histo_Wjets_name);
  if (namechan=="all") {
    TString h_name2_Wjets=plotname+"mumue"+selection+"_Wjets";
    TH1F * h_2_Wjets = (TH1F*)filechannel->Get(h_name2_Wjets);
    TString h_name3_Wjets=plotname+"eemu"+selection+"_Wjets";
    TH1F * h_3_Wjets = (TH1F*)filechannel->Get(h_name3_Wjets);
    TString h_name4_Wjets=plotname+"eee"+selection+"_Wjets";
    TH1F * h_4_Wjets = (TH1F*)filechannel->Get(h_name4_Wjets);
    histo_Wjets->Add(histo_Wjets, h_2_Wjets     , 1, 1);
    histo_Wjets->Add(histo_Wjets, h_3_Wjets     , 1, 1);
    histo_Wjets->Add(histo_Wjets, h_4_Wjets     , 1, 1);
  }
  

  TString histo_TZq_name = plotname+channel+selection+"_TZq";
  histo_TZq         = (TH1F*)filechannel->Get(histo_TZq_name);
  if (namechan=="all") {
    TString h_name2_TZq=plotname+"mumue"+selection+"_TZq";
    TH1F * h_2_TZq = (TH1F*)filechannel->Get(h_name2_TZq);
    TString h_name3_TZq=plotname+"eemu"+selection+"_TZq";
    TH1F * h_3_TZq = (TH1F*)filechannel->Get(h_name3_TZq);
    TString h_name4_TZq=plotname+"eee"+selection+"_TZq";
    TH1F * h_4_TZq = (TH1F*)filechannel->Get(h_name4_TZq);
    histo_TZq->Add(histo_TZq, h_2_TZq     , 1, 1);
    histo_TZq->Add(histo_TZq, h_3_TZq     , 1, 1);
    histo_TZq->Add(histo_TZq, h_4_TZq     , 1, 1);
  }

  TString histo_SingleToptW_name = plotname+channel+selection+"_SingleToptW";
  histo_SingleToptW         = (TH1F*)filechannel->Get(histo_SingleToptW_name);
  //     if (namechan=="all") {
  //         TString h_name2_SingleToptW=plotname+"mumu"+selection+"_SingleToptW";
  //         TH1F * h_2_SingleToptW = (TH1F*)filechannel->Get(h_name2_SingleToptW);
  //         TString h_name3_SingleToptW=plotname+"emu"+selection+"_SingleToptW";
  //         TH1F * h_3_SingleToptW = (TH1F*)filechannel->Get(h_name3_SingleToptW);
  //         histo_SingleToptW->Add(histo_SingleToptW, h_2_SingleToptW     , 1, 1);
  //         histo_SingleToptW->Add(histo_SingleToptW, h_3_SingleToptW     , 1, 1);
  //     }
  

  TString histo_TtW_name = plotname+channel+selection+"_TtW";
  histo_TtW         = (TH1F*)filechannel->Get(histo_TtW_name);
  if (namechan=="all") {
    TString h_name2_TtW=plotname+"mumue"+selection+"_TtW";
    TH1F * h_2_TtW = (TH1F*)filechannel->Get(h_name2_TtW);
    TString h_name3_TtW=plotname+"eemu"+selection+"_TtW";
    TH1F * h_3_TtW = (TH1F*)filechannel->Get(h_name3_TtW);
    TString h_name4_TtW=plotname+"eee"+selection+"_TtW";
    TH1F * h_4_TtW = (TH1F*)filechannel->Get(h_name4_TtW);
    histo_TtW->Add(histo_TtW, h_2_TtW     , 1, 1);
    histo_TtW->Add(histo_TtW, h_3_TtW     , 1, 1);
    histo_TtW->Add(histo_TtW, h_4_TtW     , 1, 1);
  }
  

  
  TString histo_TbartW_name = plotname+channel+selection+"_TbartW";
  histo_TbartW         = (TH1F*)filechannel->Get(histo_TbartW_name);
  if (namechan=="all") {
    TString h_name2_TbartW=plotname+"mumue"+selection+"_TbartW";
    TH1F * h_2_TbartW = (TH1F*)filechannel->Get(h_name2_TbartW);
    TString h_name3_TbartW=plotname+"eemu"+selection+"_TbartW";
    TH1F * h_3_TbartW = (TH1F*)filechannel->Get(h_name3_TbartW);
    TString h_name4_TbartW=plotname+"eee"+selection+"_TbartW";
    TH1F * h_4_TbartW = (TH1F*)filechannel->Get(h_name4_TbartW);
    histo_TbartW->Add(histo_TbartW, h_2_TbartW     , 1, 1);
    histo_TbartW->Add(histo_TbartW, h_3_TbartW     , 1, 1);
    histo_TbartW->Add(histo_TbartW, h_4_TbartW     , 1, 1);
  }
  
  

  TString histo_SingleToptChan_name = plotname+channel+selection+"_SingleToptChan";
  histo_SingleToptChan         = (TH1F*)filechannel->Get(histo_SingleToptChan_name);
  //     if (namechan=="all") {
  //         TString h_name2_SingleToptChan=plotname+"mumu"+selection+"_SingleToptChan";
  //         TH1F * h_2_SingleToptChan = (TH1F*)filechannel->Get(h_name2_SingleToptChan);
  //         TString h_name3_SingleToptChan=plotname+"emu"+selection+"_SingleToptChan";
  //         TH1F * h_3_SingleToptChan = (TH1F*)filechannel->Get(h_name3_SingleToptChan);
  //         histo_SingleToptChan->Add(histo_SingleToptChan, h_2_SingleToptChan     , 1, 1);
  //         histo_SingleToptChan->Add(histo_SingleToptChan, h_3_SingleToptChan     , 1, 1);
  //     }
  

  TString histo_TtChan_name = plotname+channel+selection+"_TtChan";
  histo_TtChan         = (TH1F*)filechannel->Get(histo_TtChan_name);
  if (namechan=="all") {
    TString h_name2_TtChan=plotname+"mumue"+selection+"_TtChan";
    TH1F * h_2_TtChan = (TH1F*)filechannel->Get(h_name2_TtChan);
    TString h_name3_TtChan=plotname+"eemu"+selection+"_TtChan";
    TH1F * h_3_TtChan = (TH1F*)filechannel->Get(h_name3_TtChan);
    TString h_name4_TtChan=plotname+"eee"+selection+"_TtChan";
    TH1F * h_4_TtChan = (TH1F*)filechannel->Get(h_name4_TtChan);
    histo_TtChan->Add(histo_TtChan, h_2_TtChan     , 1, 1);
    histo_TtChan->Add(histo_TtChan, h_3_TtChan     , 1, 1);
    histo_TtChan->Add(histo_TtChan, h_4_TtChan     , 1, 1);
  }
  

  
  TString histo_TbartChan_name = plotname+channel+selection+"_TbartChan";
  histo_TbartChan         = (TH1F*)filechannel->Get(histo_TbartChan_name);
  if (namechan=="all") {
    TString h_name2_TbartChan=plotname+"mumue"+selection+"_TbartChan";
    TH1F * h_2_TbartChan = (TH1F*)filechannel->Get(h_name2_TbartChan);
    TString h_name3_TbartChan=plotname+"eemu"+selection+"_TbartChan";
    TH1F * h_3_TbartChan = (TH1F*)filechannel->Get(h_name3_TbartChan);
    TString h_name4_TbartChan=plotname+"eee"+selection+"_TbartChan";
    TH1F * h_4_TbartChan = (TH1F*)filechannel->Get(h_name4_TbartChan);
    histo_TbartChan->Add(histo_TbartChan, h_2_TbartChan     , 1, 1);
    histo_TbartChan->Add(histo_TbartChan, h_3_TbartChan     , 1, 1);
    histo_TbartChan->Add(histo_TbartChan, h_4_TbartChan     , 1, 1);
  }
  
  
  TString histo_SingleTopsChan_name = plotname+channel+selection+"_SingleTopsChan";
  histo_SingleTopsChan         = (TH1F*)filechannel->Get(histo_SingleTopsChan_name);
  //     if (namechan=="all") {
  //         TString h_name2_SingleTopsChan=plotname+"mumu"+selection+"_SingleTopsChan";
  //         TH1F * h_2_SingleTopsChan = (TH1F*)filechannel->Get(h_name2_SingleTopsChan);
  //         TString h_name3_SingleTopsChan=plotname+"emu"+selection+"_SingleTopsChan";
  //         TH1F * h_3_SingleTopsChan = (TH1F*)filechannel->Get(h_name3_SingleTopsChan);
  //         histo_SingleTopsChan->Add(histo_SingleTopsChan, h_2_SingleTopsChan     , 1, 1);
  //         histo_SingleTopsChan->Add(histo_SingleTopsChan, h_3_SingleTopsChan     , 1, 1);
  //     }
  

  TString histo_TsChan_name = plotname+channel+selection+"_TsChan";
  histo_TsChan         = (TH1F*)filechannel->Get(histo_TsChan_name);
  if (namechan=="all") {
    TString h_name2_TsChan=plotname+"mumue"+selection+"_TsChan";
    TH1F * h_2_TsChan = (TH1F*)filechannel->Get(h_name2_TsChan);
    TString h_name3_TsChan=plotname+"eemu"+selection+"_TsChan";
    TH1F * h_3_TsChan = (TH1F*)filechannel->Get(h_name3_TsChan);
    TString h_name4_TsChan=plotname+"eee"+selection+"_TsChan";
    TH1F * h_4_TsChan = (TH1F*)filechannel->Get(h_name4_TsChan);
    histo_TsChan->Add(histo_TsChan, h_2_TsChan     , 1, 1);
    histo_TsChan->Add(histo_TsChan, h_3_TsChan     , 1, 1);
    histo_TsChan->Add(histo_TsChan, h_4_TsChan     , 1, 1);
  }
  

  
  TString histo_TbarsChan_name = plotname+channel+selection+"_TbarsChan";
  histo_TbarsChan         = (TH1F*)filechannel->Get(histo_TbarsChan_name);
  if (namechan=="all") {
    TString h_name2_TbarsChan=plotname+"mumue"+selection+"_TbarsChan";
    TH1F * h_2_TbarsChan = (TH1F*)filechannel->Get(h_name2_TbarsChan);
    TString h_name3_TbarsChan=plotname+"eemu"+selection+"_TbarsChan";
    TH1F * h_3_TbarsChan = (TH1F*)filechannel->Get(h_name3_TbarsChan);
    TString h_name4_TbarsChan=plotname+"eee"+selection+"_TbarsChan";
    TH1F * h_4_TbarsChan = (TH1F*)filechannel->Get(h_name4_TbarsChan);
    histo_TbarsChan->Add(histo_TbarsChan, h_2_TbarsChan     , 1, 1);
    histo_TbarsChan->Add(histo_TbarsChan, h_3_TbarsChan     , 1, 1);
    histo_TbarsChan->Add(histo_TbarsChan, h_4_TbarsChan     , 1, 1);
  }
  
  
  
  TString histo_VV_name = plotname+channel+selection+"_VV";
  histo_VV                  = (TH1F*)filechannel->Get(histo_VV_name);
  //     if (namechan=="all") {
  //         TString h_name2_VV=plotname+"mumu"+selection+"_VV";
  //         TH1F * h_2_VV = (TH1F*)filechannel->Get(h_name2_VV);
  //         TString h_name3_VV=plotname+"emu"+selection+"_VV";
  //         TH1F * h_3_VV = (TH1F*)filechannel->Get(h_name3_VV);
  //         histo_VV->Add(histo_VV, h_2_VV     , 1, 1);
  //         histo_VV->Add(histo_VV, h_3_VV     , 1, 1);
  //     }
  
  TString histo_WW_name = plotname+channel+selection+"_WW";
  histo_WW                  = (TH1F*)filechannel->Get(histo_WW_name);
  if (namechan=="all") {
    TString h_name2_WW=plotname+"mumue"+selection+"_WW";
    TH1F * h_2_WW = (TH1F*)filechannel->Get(h_name2_WW);
    TString h_name3_WW=plotname+"eemu"+selection+"_WW";
    TH1F * h_3_WW = (TH1F*)filechannel->Get(h_name3_WW);
    TString h_name4_WW=plotname+"eee"+selection+"_WW";
    TH1F * h_4_WW = (TH1F*)filechannel->Get(h_name4_WW);
    histo_WW->Add(histo_WW, h_2_WW     , 1, 1);
    histo_WW->Add(histo_WW, h_3_WW     , 1, 1);
    histo_WW->Add(histo_WW, h_4_WW     , 1, 1);
  }
  

  TString histo_WZ_name = plotname+channel+selection+"_WZ";
  histo_WZ                  = (TH1F*)filechannel->Get(histo_WZ_name);
  if (namechan=="all") {
    TString h_name2_WZ=plotname+"mumue"+selection+"_WZ";
    TH1F * h_2_WZ = (TH1F*)filechannel->Get(h_name2_WZ);
    TString h_name3_WZ=plotname+"eemu"+selection+"_WZ";
    TH1F * h_3_WZ = (TH1F*)filechannel->Get(h_name3_WZ);
    TString h_name4_WZ=plotname+"eee"+selection+"_WZ";
    TH1F * h_4_WZ = (TH1F*)filechannel->Get(h_name4_WZ);
    histo_WZ->Add(histo_WZ, h_2_WZ     , 1, 1);
    histo_WZ->Add(histo_WZ, h_3_WZ     , 1, 1);
    histo_WZ->Add(histo_WZ, h_4_WZ     , 1, 1);
  }
  

  TString histo_ZZ_name = plotname+channel+selection+"_ZZ";
  histo_ZZ                  = (TH1F*)filechannel->Get(histo_ZZ_name);
  if (namechan=="all") {
    TString h_name2_ZZ=plotname+"mumue"+selection+"_ZZ";
    TH1F * h_2_ZZ = (TH1F*)filechannel->Get(h_name2_ZZ);
    TString h_name3_ZZ=plotname+"eemu"+selection+"_ZZ";
    TH1F * h_3_ZZ = (TH1F*)filechannel->Get(h_name3_ZZ);
    TString h_name4_ZZ=plotname+"eee"+selection+"_ZZ";
    TH1F * h_4_ZZ = (TH1F*)filechannel->Get(h_name4_ZZ);
    histo_ZZ->Add(histo_ZZ, h_2_ZZ     , 1, 1);
    histo_ZZ->Add(histo_ZZ, h_3_ZZ     , 1, 1);
    histo_ZZ->Add(histo_ZZ, h_4_ZZ     , 1, 1);
  }

  

  
  if(readStop){
    
  }
  
  histo_Zjets->Add(histo_Zjets, histo_DYToLL_M10_50       , 1, 1);
  histo_TTbarBkg->Add(histo_TTbarBkg, histo_TTbarSig, 1, 1);

  
  if ( histo_SingleToptW==NULL ) { 
    histo_SingleToptW=(TH1F*)histo_TtW->Clone("SingleToptW");
    histo_SingleToptW->Add(histo_SingleToptW, histo_TbartW     , 1, 1);
    histo_SingleToptW->Add(histo_SingleToptW, histo_TtChan     , 1, 1);
    histo_SingleToptW->Add(histo_SingleToptW, histo_TbartChan  , 1, 1);
    histo_SingleToptW->Add(histo_SingleToptW, histo_TsChan     , 1, 1);
    histo_SingleToptW->Add(histo_SingleToptW, histo_TbarsChan  , 1, 1);
  }
  
  if ( histo_VV==NULL ) {
    histo_VV=(TH1F*)histo_WW->Clone("VV");
    histo_VV->Add(histo_VV, histo_WZ     , 1, 1);
    histo_VV->Add(histo_VV, histo_ZZ     , 1, 1);
  }   
  
  histo_TTbarBkg->SetFillStyle(1001);
  histo_Zjets->SetFillStyle(1001);
  histo_Wjets->SetFillStyle(1001);
  histo_VV->SetFillStyle(1001);
  histo_SingleToptW->SetFillStyle(1001);
  histo_TZq->SetFillStyle(1001);
  
  histo_TTbarBkg->SetFillColor(kRed-7);
  histo_Zjets->SetFillColor(kAzure-2);
  histo_Wjets->SetFillColor(kGreen-3);
  //    histo_VV->SetFillColor(kWhite);
  histo_VV->SetFillColor(13);
  histo_SingleToptW->SetFillColor(kMagenta);
  histo_TZq->SetFillColor(kAzure+8);
  // color code for Z/DY -> tau tau : kAzure+8
  
  TH1F * hmc= (TH1F*) histo_TTbarBkg->Clone();
  hmc->Add(hmc,histo_VV, 1., 1.);
  //hmc->Add(hmc,histo_Wjets, 1., 1.);
  hmc->Add(hmc,histo_SingleToptW, 1., 1.);
  //hmc->Add(hmc,histo_TTbarBkg, 1., 1.);
  hmc->Add(hmc,histo_Zjets, 1., 1.);
  hmc->SetName("hmc");
  
  
  
  hmc->SetMinimum(0.001);
  
  if(plotname =="InvM_ll_" ){
    int bins = 5;
    hmc->Rebin(bins);
    histo_VV->Rebin(bins);
    histo_Wjets->Rebin(bins);
    histo_SingleToptW->Rebin(bins);
    histo_TZq->Rebin(bins);
    histo_TTbarBkg->Rebin(bins);
    histo_Zjets->Rebin(bins);
    histo_Data->Rebin(bins);
    histo_FCNCzut->Rebin(bins);
  }
  
  if(plotname =="HT_" ){
    int bins = 10;
    hmc->Rebin(bins);
    histo_VV->Rebin(bins);
    histo_Wjets->Rebin(bins);
    histo_SingleToptW->Rebin(bins);
    histo_TZq->Rebin(bins);
    histo_TTbarBkg->Rebin(bins);
    histo_Zjets->Rebin(bins);
    histo_Data->Rebin(bins);
    histo_FCNCzut->Rebin(bins);
  }
  
  
  if(plotname =="MET_" ){
    int bins = 1;
    hmc->Rebin(bins);
    histo_VV->Rebin(bins);
    histo_Wjets->Rebin(bins);
    histo_SingleToptW->Rebin(bins);
    histo_TZq->Rebin(bins);
    histo_TTbarBkg->Rebin(bins);
    histo_Zjets->Rebin(bins);
    histo_Data->Rebin(bins);
    histo_FCNCzut->Rebin(bins);
  }
  
  if(plotname =="Mt_" || plotname == "mWT_" || plotname == "DijetInvM_"){
    int bins = 5;
    hmc->Rebin(bins);
    histo_VV->Rebin(bins);
    histo_Wjets->Rebin(bins);
    histo_SingleToptW->Rebin(bins);
    histo_TZq->Rebin(bins);
    histo_TTbarBkg->Rebin(bins);
    histo_Zjets->Rebin(bins);
    histo_Data->Rebin(bins);
    histo_FCNCzut->Rebin(bins);
  }
  
  if(plotname =="LeptWPt_" ){
    int bins = 5;
    hmc->Rebin(bins);
    histo_VV->Rebin(bins);
    histo_Wjets->Rebin(bins);
    histo_SingleToptW->Rebin(bins);
    histo_TZq->Rebin(bins);
    histo_TTbarBkg->Rebin(bins);
    histo_Zjets->Rebin(bins);
    histo_Data->Rebin(bins);
    histo_FCNCzut->Rebin(bins);
  }
  
  if(plotname =="LeptZPt_" ){
    int bins = 5;
    hmc->Rebin(bins);
    histo_VV->Rebin(bins);
    histo_Wjets->Rebin(bins);
    histo_SingleToptW->Rebin(bins);
    histo_TZq->Rebin(bins);
    histo_TTbarBkg->Rebin(bins);
    histo_Zjets->Rebin(bins);
    histo_Data->Rebin(bins);
    histo_FCNCzut->Rebin(bins);
  }
  
  
  if(plotname =="JetPt_" ){
    int bins = 5;
    hmc->Rebin(bins);
    histo_VV->Rebin(bins);
    histo_Wjets->Rebin(bins);
    histo_SingleToptW->Rebin(bins);
    histo_TZq->Rebin(bins);
    histo_TTbarBkg->Rebin(bins);
    histo_Zjets->Rebin(bins);
    histo_Data->Rebin(bins);
    histo_FCNCzut->Rebin(bins);
  }
  
  
  if(plotname =="Asym_" ){
    int bins = 2;
    hmc->Rebin(bins);
    histo_VV->Rebin(bins);
    histo_Wjets->Rebin(bins);
    histo_SingleToptW->Rebin(bins);
    histo_TZq->Rebin(bins);
    histo_TTbarBkg->Rebin(bins);
    histo_Zjets->Rebin(bins);
    histo_Data->Rebin(bins);
    histo_FCNCzut->Rebin(bins);
  }
  
  if(plotname =="RecoTopMass_" ){
    int bins = 10;
    hmc->Rebin(bins);
    histo_VV->Rebin(bins);
    histo_Wjets->Rebin(bins);
    histo_SingleToptW->Rebin(bins);
    histo_TZq->Rebin(bins);
    histo_TTbarBkg->Rebin(bins);
    histo_Zjets->Rebin(bins);
    histo_Data->Rebin(bins);
    histo_FCNCzut->Rebin(bins);
  }
  if(plotname =="RecoPtZ_" ){
    int bins = 10;
    hmc->Rebin(bins);
    histo_VV->Rebin(bins);
    histo_Wjets->Rebin(bins);
    histo_SingleToptW->Rebin(bins);
    histo_TZq->Rebin(bins);
    histo_TTbarBkg->Rebin(bins);
    histo_Zjets->Rebin(bins);
    histo_Data->Rebin(bins);
    histo_FCNCzut->Rebin(bins);
  }
  if(plotname =="deltaPhilb_" ){
    int bins = 20;
    hmc->Rebin(bins);
    histo_VV->Rebin(bins);
    histo_Wjets->Rebin(bins);
    histo_SingleToptW->Rebin(bins);
    histo_TZq->Rebin(bins);
    histo_TTbarBkg->Rebin(bins);
    histo_Zjets->Rebin(bins);
    histo_Data->Rebin(bins);
    histo_FCNCzut->Rebin(bins);
  }
  
  histo_ratio = (TH1F*) histo_Data->Clone();
  
  
  
  TH1F * lumiband = (TH1F*) hmc->Clone();
  frac_ee   = frac_ee/histo_Zjets->Integral();
  frac_mumu = frac_mumu/histo_Zjets->Integral() ;
  frac_emu  = frac_emu/histo_Zjets->Integral();
  
  for (int ilum=0; ilum<lumiband->GetNbinsX(); ilum++) {
    
    double error_all = 
       pow(
         (histo_VV->GetBinContent(ilum+1)+histo_SingleToptW->GetBinContent(ilum+1)+ histo_TTbarBkg->GetBinContent(ilum+1))
         *lumi_error, 
         2)+
       //*************************
       //uncertinty on trigger eff
       pow(
       	 (histo_VV->GetBinContent(ilum+1)+histo_SingleToptW->GetBinContent(ilum+1)+ histo_TTbarBkg->GetBinContent(ilum+1))
	 *SF_trigger_error
	 , 2)+
       //*************************
       //uncertinty on lepton sel
       pow(
       	 (histo_VV->GetBinContent(ilum+1)+histo_SingleToptW->GetBinContent(ilum+1)+ histo_TTbarBkg->GetBinContent(ilum+1))
	 *SF_Lepton_error
	 , 2);
     
       error_all += pow(histo_Zjets->GetBinContent(ilum+1)*0.30, 2) ;
       error_all += pow(histo_VV->GetBinContent(ilum+1)*0.20, 2) ;

     
     
     lumiband->SetBinError(ilum+1,sqrt(error_all));
     
     //modifications
     hmc->SetBinError(ilum+1,sqrt(error_all));
     
     /*cout << "------------------------------" << endl;
     cout << "bin content is Data "  << histo_Data->GetBinContent(ilum+1) << endl;
     cout << "bin content is VV   "  << histo_VV->GetBinContent(ilum+1) << endl;
     cout << "bin content is WZ   "  << histo_WZ->GetBinContent(ilum+1) << endl;
     cout << "bin content is DY   "  << histo_Zjets->GetBinContent(ilum+1) << endl;
     cout << "bin content is WZ*0.2   "  << histo_WZ->GetBinContent(ilum+1)*0.20 << endl;
     cout << "bin content is DY*0.3   "  << histo_Zjets->GetBinContent(ilum+1)*0.30<< endl;
     cout << "uncertainty MC is   "  << lumiband->GetBinError(ilum+1) << endl;
     cout << "uncertainty MC is   "  << sqrt(error_all) << endl;*/
    
  }
  TGraphErrors *thegraph = new TGraphErrors(lumiband);
  thegraph->SetFillStyle(3005);
  thegraph->SetFillColor(1);
  
  
  THStack* hs= new THStack();
  hs->Add(histo_SingleToptW);
  hs->Add(histo_TTbarBkg);
  hs->Add(histo_TZq);
  hs->Add(histo_VV);
  hs->Add(histo_Zjets);
  
  histo_Data->GetXaxis()->SetTitle("");
  histo_Data->GetYaxis()->SetTitle("");
  histo_Data->GetYaxis()->CenterTitle();
  histo_Data->GetYaxis()->SetTitle("");
  histo_Data->GetXaxis()->SetLabelSize(0);
  histo_Data->GetXaxis()->SetTitleSize(0);
  
  hmc->GetYaxis()->CenterTitle();
  hmc->GetYaxis()->SetTitle("");
  hmc->GetXaxis()->SetLabelSize(0);
  hmc->GetXaxis()->SetTitleSize(0);
  
  
  if (plotname=="NVxP_" && namechan=="all" && selection=="_aftertrigger") {
    histo_Data->Scale(1./histo_Data->Integral());
    float scaaa=1./hmc->Integral();
    hmc->Scale(scaaa);
    hmc->SetMaximum(histo_Data->GetMaximum());
  }
  
  if (plotname=="NVxP_" && namechan=="all" && selection=="_aftertrigger") {
    hmc->Draw();
  }
  else {
    hs->Draw();
    if( plotname=="Njets_" || plotname == "NBjets_") {
      hs->GetXaxis()->SetNdivisions( 505 );
      int nbb=hmc->GetNbinsX();
      for (int inbb=0; inbb<nbb; inbb++) {
	char namelast[10];
	if (inbb<nbb-1) sprintf(namelast,"%d",hs->GetXaxis()->GetBinCenter(inbb+1));
	else sprintf(namelast,"#geq %d",hs->GetXaxis()->GetBinCenter(inbb+1));
	//hs->GetXaxis()->SetBinLabel(inbb+1, namelast);
	hs->GetXaxis()->SetBinLabel(inbb+1, "");
      }
      hs->GetXaxis()->SetLabelSize(0.);
      hs->SetMaximum(histo_Data->GetMaximum()*(1.3));
      hs->SetMinimum(0.1);
    }
    thegraph->Draw("e2same");
    if (!setlogy) {
      hs->SetMinimum(0.);
      hs->SetMaximum(hs->GetMaximum()*1.4);
      if (hs->GetMaximum()<histo_Data->GetMaximum()) {
	hs->SetMaximum(histo_Data->GetMaximum()*(1.7));
      }
      if (plotname == "Met_") {
	hs->SetMaximum(hs->GetMaximum()*(1.7));
      } hs->GetXaxis()->SetLabelSize(0.);
    }
    else {
      if (hs->GetMaximum()<histo_Data->GetMaximum()) {
	hs->SetMaximum(histo_Data->GetMaximum()*(2));
      }
      hs->SetMinimum(0.1);
      //     hs->SetMinimum(0.1);
    }
    hs->GetXaxis()->SetLabelSize(0.);
    hs->SetMinimum(0.1);
    TString theXtitle = "";
    if(plotname == "Inv")theXtitle ="M_{ll} [GeV/c^{-1}]";
    else if(plotname == "Njets_")   theXtitle ="jet mult.";
    else if(plotname == "PtJet1_")   theXtitle ="Leading Jet p_{T} [GeV/c]";
    else if(plotname == "PtJet2_")   {theXtitle ="Second Leading Jet p_{T} [GeV/c]"; hs->SetMaximum(histo_Data->GetMaximum()*(1.3));}
    else if(plotname == "NBjets_")  theXtitle ="b-tagged jet mult.";
    else if(plotname == "Met_")    theXtitle ="missing E_{T} [GeV]";
    //cout << "test title  " << endl;
    if(plotname == "Inv")              histo_ratio->GetXaxis()->SetTitle("M_{ll} [GeV/c^{-1}]");
    else if(plotname == "Njets_")    histo_ratio->GetXaxis()->SetTitle("jet mult.");
    else if(plotname == "PtJet1_")   histo_ratio->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
    else if(plotname == "PtJet2_")   histo_ratio->GetXaxis()->SetTitle("Second Leading Jet p_{T} [GeV]");
    else if(plotname == "NBjets_")   histo_ratio->GetXaxis()->SetTitle("b-tagged jet mult.");
    else if(plotname == "Met_")      histo_ratio->GetXaxis()->SetTitle("missing E_{T} [GeV]");
    else if(plotname == "PtMu1_")   histo_ratio->GetXaxis()->SetTitle("Leading muon p_{T} [GeV]");
    else if(plotname == "PtMu2_")   histo_ratio->GetXaxis()->SetTitle("Second Leading muon p_{T} [GeV]");
    else if(plotname == "PtEl1_")   histo_ratio->GetXaxis()->SetTitle("Leading electron p_{T} [GeV]");
    else if(plotname == "PtEl2_")   histo_ratio->GetXaxis()->SetTitle("Second Leading electron p_{T} [GeV]");
    
    
  }
  histo_ratio->GetXaxis()->SetLabelSize(0.08);
  
  // cout << " 1 histo_ratio->GetTitle()  " << histo_ratio->GetXaxis()->GetTitle() << endl;;
  histo_Data->Draw("epsame");
  histo_Data->SetMarkerStyle(20);
  
  
  histo_Data->SetMarkerSize(1.2);
  histo_Data->Draw("epsame");
  
  histo_FCNCzut->SetLineWidth(2.0);
  histo_FCNCzut->SetLineColor(1.0);
  histo_FCNCzut->Draw("same");
  cout << histo_FCNCzut->Integral() << endl;
  
   
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->SetTextAlign(31); 
  latex->DrawLatex(0.45, 0.95, "CMS Preliminary");
  
  
  
  TLatex *latex2 = new TLatex();
  latex2->SetNDC();
  latex2->SetTextSize(0.04);
  latex2->SetTextAlign(31); 
  latex2->DrawLatex(0.87, 0.95, "4.9 fb^{-1} at #sqrt{s} = 7 TeV");
  
  TString info_data; 
  if (namechan=="eee")    info_data = "eee channel";
  if (namechan=="emumu")  info_data = "e#mu#mu channel";
  if (namechan=="mumue")  info_data = "#mu#mu e channel";
  if (namechan=="mumumu") info_data = "#mu#mu #mu channel";
  if (namechan=="all")    info_data = "eee, #mu#mu#mu, e#mu#mu, ee#mu channels";
  
  //text2 = new TLatex(0.15,0.8, info_data);
  TLatex * text2 = new TLatex(0.45,0.98, info_data);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.18);
  text2->SetY(0.92);
  //text2->SetLineWidth(2);
  text2->SetTextFont(42);
  text2->SetTextSize(0.0610687);
  //    text2->SetTextSizePixels(24);// dflt=28
  text2->Draw();
  
  histo_ratio->GetXaxis()->SetLabelSize(0.08);
  TLegend* qw = 0;
  /*if (setlogy && plotname=="Njets_" && selection=="_aftermetcut") {
    qw = new TLegend(0.17,0.70,0.40,0.98);
    }
    else {
    if (plotname!="NBjets_") {
    qw = new TLegend(0.75,0.70,0.98,0.98);
    }
    else {
    qw = new TLegend(0.75,0.70,0.98,0.98);
    }
    }*/
  qw = new TLegend(.80,.60,.95,.90);
  if(plotname == "NJet_") qw = new TLegend(.20,.50,.35,.80);
  
  
  qw->SetShadowColor(0);
  qw->SetFillColor(0);
  qw->SetLineColor(0);
  
  
  //  qw->SetHeader("CMS Preliminary, 881.8 pb^{-1}");
  //  qw->AddEntry(histo_Data,         info_data,                "p");
  qw->AddEntry(histo_Data,         "Data" ,                "ep");
  qw->AddEntry(histo_Zjets,        "DY "                  ,"f");
  qw->AddEntry(histo_VV,           "VV "                  ,"f");
  qw->AddEntry(histo_TZq,  "tZq "                  ,"f");
  qw->AddEntry(histo_TTbarBkg,     " t#bar{t} background "     ,"f");
  qw->AddEntry(histo_SingleToptW,  "single top "                  ,"f");
  //qw->AddEntry(histo_Wjets,        "W  "                  ,"f");
  //qw->AddEntry(histo_TTbarBkg,     "t#bar{t} other  "     ,"f");
  qw->AddEntry(histo_FCNCzut,     " tZ "     ,"l");
  qw->SetFillColor(0);
  qw->SetTextFont(42);
  qw->Draw();
  
  /*
    
  // caro print 
  cout << "data " << selection << "   " << histo_Data->GetMean() << "    " << histo_Data->GetMeanError() << endl;
  cout << "mc   " << selection << "   " << hmc->GetMean() << "    " << hmc->GetMeanError() << endl;
  */
  
  
  
  /*
  // caro : check histo integral with table!
  int nbinmax=histo_Data->GetNbinsX();
  cout << " selection " << selection << endl;
  cout << "data       "  << histo_Data->Integral(0,nbinmax+1) << endl;
  cout << "total bg   " << histo_TTbarBkg->Integral(0,nbinmax+1)+histo_Zjets->Integral(0,nbinmax+1) +histo_Wjets->Integral(0,nbinmax+1)+histo_SingleToptW->Integral(0,nbinmax+1)+histo_VV->Integral(0,nbinmax+1) << endl;
  cout << "ttbar bg   " << histo_TTbarBkg->Integral(0,nbinmax+1) << endl;
  cout << "Single top " << histo_SingleToptW->Integral(0,nbinmax+1) << endl;
  cout << "Zjets      " << histo_Zjets->Integral(0,nbinmax+1) << endl;
  cout << "Wjets      " << histo_Wjets->Integral(0,nbinmax+1) << endl;
  cout << "VV         " << histo_VV->Integral(0,nbinmax+1) << endl;
  */
  histo_ratio->GetXaxis()->SetLabelSize(0.08);
  
  c1->cd();
  
  if (ratio) {
    //TPad *canvas_2 = new TPad("canvas_2", "canvas_2",0,0.,1.0,0.34);
    TPad *canvas_2 = new TPad("canvas_2", "canvas_2", 0.0, 0.0, 1.0, 1.0);
    canvas_2->SetTopMargin(0.7);
    canvas_2->SetFillColor(0);
    canvas_2->SetFillStyle(0);
    canvas_2->SetGridy(1);
    canvas_2->Draw();
    canvas_2->cd(0);
    //gPad->SetBottomMargin(0.375);
    //gPad->SetGridy();
    //gPad->SetGridx();
    
    
    
    
    
    
    //histo_ratio = (TH1F*) histo_Data->Clone();
    //histo_ratio->SetName("histo_ratio");
    //cout << " 721 histo_ratio->GetTitle()  " << histo_ratio->GetXaxis()->GetTitle() << endl;;
    histo_ratio->SetTitle("");
     
     
     
     
    histo_ratio->Divide(hmc);
    //cout << " 724 histo_ratio->GetTitle()  " << histo_ratio->GetXaxis()->GetTitle() << endl;;
    if (plotname=="NVxP_" && namechan=="all" && selection=="_aftertrigger") {
      for (int iiii=0; iiii<histo_ratio->GetNbinsX(); iiii++) {
	cout << " SF_VxP1[" << iiii << "]= " << histo_ratio->GetBinContent(iiii+1) << "; " <<endl;
      }
    }
    if(plotname == "Inv")              histo_ratio->GetXaxis()->SetTitle("M_{ll} [GeV/c^{-1}]");
    else if(plotname == "Njets_")    histo_ratio->GetXaxis()->SetTitle("jet mult.");
    else if(plotname == "PtJet1_")   histo_ratio->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
    else if(plotname == "PtJet2_")   histo_ratio->GetXaxis()->SetTitle("Second Leading Jet p_{T} [GeV]");
    else if(plotname == "NBjets_")   histo_ratio->GetXaxis()->SetTitle("b-tagged jet mult.");
    else if(plotname == "Met_")      histo_ratio->GetXaxis()->SetTitle("missing E_{T} [GeV]");
    else if(plotname == "PtMu1_")   histo_ratio->GetXaxis()->SetTitle("Leading muon p_{T} [GeV]");
    else if(plotname == "PtMu2_")   histo_ratio->GetXaxis()->SetTitle("Second Leading muon p_{T} [GeV]");
    else if(plotname == "PtEl1_")   histo_ratio->GetXaxis()->SetTitle("Leading electron p_{T} [GeV]");
    else if(plotname == "PtEl2_")   histo_ratio->GetXaxis()->SetTitle("Second Leading electron p_{T} [GeV]");
    
    
    histo_ratio->GetXaxis()->SetLabelSize(0.08);
    //cout << "x label size " << histo_ratio->GetXaxis()->GetLabelSize() << endl;
    /*histo_ratio->SetMaximum( 1.5 );
      histo_ratio->SetMinimum(0.5);
      histo_ratio->SetMarkerStyle(20);
      histo_ratio->SetMarkerSize(0.7);
      histo_ratio->GetYaxis()->SetTitle("Data/MC");
      histo_ratio->GetYaxis()->SetTitleFont(42);
      histo_ratio->GetYaxis()->SetLabelFont(42);
      histo_ratio->GetXaxis()->SetTitleFont(42);
      histo_ratio->GetXaxis()->SetLabelFont(42);
      histo_ratio->GetXaxis()->SetLabelSize(0.1);*/
    
    
    if (plotname=="NVxP_" && namechan=="all" && selection=="_aftertrigger") {
      histo_ratio->GetXaxis()->SetTitle("NVxP_all_aftertrigger");
    }
    else {
      //histo_ratio->GetXaxis()->SetTitle(); 
      if( plotname=="Njets_" || plotname == "NBjets_") {
	histo_ratio->GetXaxis()->SetNdivisions( 505 );
	int nbb=histo_ratio->GetNbinsX();
	for (int inbb=0; inbb<nbb; inbb++) {
	  char namelast[10];
	  if (inbb<nbb-1) sprintf(namelast,"%d",histo_ratio->GetXaxis()->GetBinCenter(inbb+1));
	  else sprintf(namelast,"#geq %d",histo_ratio->GetXaxis()->GetBinCenter(inbb+1));
	  histo_ratio->GetXaxis()->SetBinLabel(inbb+1, namelast);
	}
	histo_ratio->GetXaxis()->SetLabelSize(0.08);
      }
    }//histo_ratio->GetXaxis()->SetLabelSize(0.15);
    /*histo_ratio->GetYaxis()->SetTitleOffset( 0.4 );
      histo_ratio->GetYaxis()->SetTitleSize( 0.1 );
      histo_ratio->GetYaxis()->SetLabelSize(0.1);
      histo_ratio->GetYaxis()->SetNdivisions( 505 );
      histo_ratio->GetXaxis()->SetTitleSize( 0.15 );*/
    histo_ratio->SetMarkerStyle(20);
    histo_ratio->SetMarkerSize(1.2);
    histo_ratio->SetMaximum( 1.5 );
    histo_ratio->SetMinimum(0.5);
    histo_ratio->GetYaxis()->SetTitle("");
    histo_ratio->GetXaxis()->SetLabelSize(0.04);
    histo_ratio->GetYaxis()->SetLabelSize(0.04);
    histo_ratio->GetYaxis()->SetNdivisions(6);
    
    histo_ratio->GetYaxis()->SetTitleSize(0.03);
    histo_ratio->SetMarkerSize(1.2);
    //histo_ratio->GetYaxis()->SetNdivisions(5);
    //ratio.Draw("e")
    
    histo_ratio->SetMinimum(0.0);
    histo_ratio->SetMaximum(2.0);
    histo_ratio->Draw("E1X0");
    
    c1->cd();
    
    
  }
  
  TString end_name;
  if(setlogy) end_name="_Logy.gif";
  else end_name=".gif"; 
  TString ratname;
  if (ratio) ratname="_r";
  else ratname="_r";
  TString outputname= "plots/"+plotname+namechan+selection+ratname+end_name;
  
  
  
  
  c1->SaveAs(outputname.Data());
  
  /*
    TFile* fout  = new TFile("interne1.root","RECREATE");
    hmc->Write();
    histo_Data->SetName("hdata");
    histo_Data->Write();
    fout->Close();
  */
  histo_Data->Delete();
  histo_TTbarBkg->Delete();
  histo_Zjets->Delete();
  histo_DYToLL_M10_50->Delete();
  histo_Wjets->Delete();
  histo_VV ->Delete();
  histo_SingleToptW->Delete();
  hs->Delete();
  hmc->Delete();
  if (ratio) {
    histo_ratio->Delete();   
  }
  thegraph->Delete();
  
  
  filechannel->Close(); 
}
