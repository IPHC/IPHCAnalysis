
#include "TString.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include <iostream>

static bool doSyst= false;
static bool doStat= true;

int rebin = 2;

bool savePlot = false;
bool doBtagUncerOnly   = false;
bool doPDFUncerOnly    = false;
bool doJESUncerOnly    = false;
bool doMatchUncerOnly  = false;
bool doScaleUncerOnly  = false;

bool doFitZ = false;

static bool plotFCNC= false;
static int chan= -1;

void PlotBDToutput(TString theVertex, TString theVariable, TString theFile){
  
  
  
  if(doSyst){
  
    doBtagUncerOnly   = false;
    doPDFUncerOnly    = false;
    doJESUncerOnly    = false;
    doMatchUncerOnly  = false;
    doScaleUncerOnly  = false;
  
  }
  
  
  
  
  double WZscale = 1.;
  double TZqscale = 0.27;
  double ZZscale = 0.06/0.08;
  
  
  Int_t stati=0;
  Bool_t  fit=1;
  Bool_t logy=0;
  
  bool setlogy = 0;
  
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
  
  // For the histo:rebin
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
  
  
  
   TFile * filechannel ;
   
   filechannel = new TFile(theFile.Data());

   std::cout << "file " << theFile << std::endl;
   
   
   TH1F * histBdt_Data      = (TH1F*)((TH1F*)filechannel->Get( (theVariable+"_Data"     ).Data()))->Clone();
   TH1F * histBdt_WZ        = (TH1F*)((TH1F*)filechannel->Get( (theVariable+"_WZ"       ).Data()))->Clone();
   TH1F * histBdt_ZZ        = (TH1F*)((TH1F*)filechannel->Get( (theVariable+"_ZZ"       ).Data()))->Clone();
   TH1F * histBdt_Zjets     = (TH1F*)((TH1F*)filechannel->Get( (theVariable+"_Zjets"    ).Data()))->Clone();
   TH1F * histBdt_DataZjets = (TH1F*)((TH1F*)filechannel->Get( (theVariable+"_DataZjets").Data()))->Clone();  
   //TH1F * histBdt_DataZjets = (TH1F*)filechannel->Get( (theVariable+"_Zjets").Data());  
   TH1F * histBdt_TTbar     = (TH1F*)((TH1F*)filechannel->Get( (theVariable+"_TTbarSig" ).Data()))->Clone();
   TH1F * histBdt_TTbarZjets= (TH1F*)((TH1F*)filechannel->Get( (theVariable+"_TTbarSigZjets" ).Data()))->Clone();
   TH1F * histBdt_TZq       = (TH1F*)((TH1F*)filechannel->Get( (theVariable+"_TZq"      ).Data()))->Clone();
   
   
   TH1F * histBdt_FCNC ;     
   if(theVariable == "MVA_BDT"){
     histBdt_FCNC = (TH1F*)filechannel->Get((theVariable+"_FCNC_"+theVertex).Data());	 
   }
   else{
     histBdt_FCNC = (TH1F*)filechannel->Get((theVariable+"_Signal").Data());	 
   } 
  
  histBdt_ZZ->Scale(ZZscale);
  
  //histBdt_FCNC->Scale(10);
  
  histBdt_Data->SetLineColor(1);
  if(theVariable!="NJets" && theVariable!="NBJets"){
    histBdt_Data->Rebin(rebin);  
    histBdt_WZ->Rebin(rebin);	     
    histBdt_ZZ->Rebin(rebin);	     
    histBdt_TTbar->Rebin(rebin);     
    histBdt_TTbarZjets->Rebin(rebin);	     
    histBdt_Zjets->Rebin(rebin);	     
    histBdt_FCNC->Rebin(rebin);  
    histBdt_DataZjets->Rebin(rebin);
    histBdt_TZq->Rebin(rebin);
  }
  
  //histBdt_WZ->Scale(0.905713);
  
  histBdt_WZ->Scale(WZscale);
  
  
  /*
  histBdt_WZ->Scale(1+2.75079*0.1);     
  histBdt_ZZ->Scale(1-0.228660*0.3);	
  histBdt_TZq->Scale(1-0.0967*0.5);
  histBdt_TTbar->Scale(1-3*0.17996);
  histBdt_DataZjets->Scale(1-0.3*0.36251);
  */
  
  
  
  histBdt_WZ->Add(histBdt_ZZ);
  //histBdt_DataZjets->Add(histBdt_DataZjets, histBdt_TTbar, 1, -0.05*histBdt_DataZjets->Integral()/histBdt_TTbar->Integral());
  if(theVariable != "MVA_BDT") {
    histBdt_DataZjets->Add(histBdt_DataZjets, histBdt_TTbarZjets, 1, -1);
    histBdt_DataZjets->Scale(histBdt_Zjets->Integral()/histBdt_DataZjets->Integral());
  }
  
  histBdt_TZq->Scale(TZqscale);
  
  cout << histBdt_DataZjets->GetMean() << endl;
  
  histBdt_WZ->SetFillColor(13);
  histBdt_DataZjets->SetFillColor(kAzure-2);
  histBdt_TZq->SetFillColor(kAzure+8);
  
  THStack* hs= new THStack();
  
  histBdt_WZ->GetXaxis()->SetLabelSize(0.);
  
  
  
  
  
  histBdt_DataZjets->GetXaxis()->SetLabelSize(0.);
  
  histBdt_TTbar->SetFillColor(kRed-7);
  
  
  hs->Add(histBdt_WZ);
  if(!doFitZ){
    hs->Add(histBdt_TTbar);
    hs->Add(histBdt_TZq);
    hs->Add(histBdt_DataZjets);
  }else{
    histBdt_Data->Add(histBdt_TTbar,	  -1.);
    histBdt_Data->Add(histBdt_TZq,	  -1.);
    histBdt_Data->Add(histBdt_DataZjets,  -1.);
  
  }


  
  
  hs->Draw("histo");
  hs->GetXaxis()->SetLabelSize(0.);
  hs->SetMaximum(histBdt_Data->GetMaximum()*1.7);
  hs->SetMinimum(0.1);
  hs->GetYaxis()->SetLabelSize(0.04);
  
  
  
  histBdt_Data->SetMarkerStyle(20);
  histBdt_Data->SetMarkerSize(1.2);
  histBdt_Data->GetXaxis()->SetLabelSize(0.);
  histBdt_Data->Draw("epsame");
  
  
  histBdt_FCNC->SetLineWidth(2.0);
  histBdt_FCNC->SetLineColor(1.);
  //histBdt_FCNC->Scale(0.01);
  histBdt_FCNC->GetXaxis()->SetLabelSize(0.);
  if(plotFCNC) histBdt_FCNC->Draw("histosame");
  
  double theWZscale = 1.0;
   if(theVertex == "zut") theWZscale = 0.72;
   //if(theVertex == "zct") theWZscale = 0.64; 
   if(theVertex == "zct") theWZscale = 0.64;
   if(theVertex == "kut") theWZscale = 0.82; 
   if(theVertex == "kct") theWZscale = 0.90; 
   
   double theWZscale_err = 1.;
   if(theVertex == "zut") theWZscale_err = 0.24;
   if(theVertex == "zct") theWZscale_err = 0.30; 
   if(theVertex == "kut") theWZscale_err = 0.23; 
   if(theVertex == "kct") theWZscale_err = 0.21; 
   
  
  TH1F* herrorband =  (TH1F*) histBdt_WZ->Clone();
  if(!doFitZ){
    herrorband->Add(histBdt_DataZjets);
    herrorband->Add(histBdt_TTbar);
    herrorband->Add(histBdt_TZq);
  }
  
  TH1F* herrorband_btagup ;
  TH1F* herrorband_btagdown   ;
  TH1F* herrorband_norw   ;
  
  
  TH1F* herrorband_PDFup ;
  TH1F* herrorband_PDFdown ;
  
  TH1F* herrorband_JESup ;
  TH1F* herrorband_JESdown ;
  
  TH1F* herrorband_Matchup ;
  TH1F* herrorband_Matchdown ;
  
  TH1F* herrorband_Scaleup ;
  TH1F* herrorband_Scaledown ;
  
  
  
  cout << "histBdt_WZ->Integral()        " << histBdt_WZ->Integral() << endl;
  cout << "histBdt_DataZjets->Integral() " << histBdt_DataZjets->Integral() << endl;
  cout << "histBdt_TTbar->Integral()     " << histBdt_TTbar->Integral() << endl;
  cout << "histBdt_TZq->Integral()       " << histBdt_TZq->Integral() << endl;
  
  
  TFile * filebta_btagup;   
  TFile * filebta_btagdown ;
  TFile * filebta_norw;
  
  TFile * filebta_PDFup ; 
  TFile * filebta_PDFdown;
  
  TFile * filebta_JESup  ;
  TFile * filebta_JESdown;
  
  TFile * filebta_Matchup ;  
  TFile * filebta_Matchdown; 
  
  TFile * filebta_Scaleup  ; 
  TFile * filebta_Scaledown ;
  
  TFile * filebta_DYup  ; 
  TFile * filebta_DYdown ;

  
  if(theVertex == "zut" && doSyst) {
    filebta_btagup    = new TFile("HistoBDToutput/TMVApp_zut_btagUp.root");
    filebta_btagdown  = new TFile("HistoBDToutput/TMVApp_zut_btagdown.root");
    //filebta_norw      = new TFile("HistoBDToutput/TMVApp_zut_noreweight.root");
  
    filebta_PDFup     = new TFile("HistoBDToutput/TMVApp_zut_PDFup.root");
    filebta_PDFdown   = new TFile("HistoBDToutput/TMVApp_zut_PDFdown.root");
  
    filebta_JESup     = new TFile("HistoBDToutput/TMVApp_zut_JESup.root");
    filebta_JESdown   = new TFile("HistoBDToutput/TMVApp_zut_JESdown.root");
  
    filebta_Matchup   = new TFile("HistoBDToutput/TMVApp_zut_Matchup.root");
    filebta_Matchdown = new TFile("HistoBDToutput/TMVApp_zut_Matchdown.root");
  
    filebta_Scaleup   = new TFile("HistoBDToutput/TMVApp_zut_Scaleup.root");
    filebta_Scaledown = new TFile("HistoBDToutput/TMVApp_zut_Scaledown.root");
  }
  
  if(theVertex == "zct"&& doSyst){
  
    filebta_btagup    = new TFile("HistoBDToutput/TMVApp_zct_btagUp.root");
    filebta_btagdown  = new TFile("HistoBDToutput/TMVApp_zct_btagdown.root");
    filebta_norw      = new TFile("HistoBDToutput/TMVApp_zct_noreweight.root");
  
    filebta_PDFup     = new TFile("HistoBDToutput/TMVApp_zct_PDFup.root");
    filebta_PDFdown   = new TFile("HistoBDToutput/TMVApp_zct_PDFdown.root");
  
    filebta_JESup     = new TFile("HistoBDToutput/TMVApp_zct_JESup.root");
    filebta_JESdown   = new TFile("HistoBDToutput/TMVApp_zct_JESdown.root");
  
    filebta_Matchup   = new TFile("HistoBDToutput/TMVApp_zct_Matchup.root");
    filebta_Matchdown = new TFile("HistoBDToutput/TMVApp_zct_Matchdown.root");
  
    filebta_Scaleup   = new TFile("HistoBDToutput/TMVApp_zct_Scaleup.root");
    filebta_Scaledown = new TFile("HistoBDToutput/TMVApp_zct_Scaledown.root");
    
    filebta_DYup   = new TFile("HistoBDToutput/TMVApp_zct_DYUp.root");
    filebta_DYdown = new TFile("HistoBDToutput/TMVApp_zct_DYDown.root");
  }
  
  if(theVertex == "kut"&& doSyst){
  
    filebta_btagup    = new TFile("HistoBDToutput/TMVApp_kut_btagUp.root");
    filebta_btagdown  = new TFile("HistoBDToutput/TMVApp_kut_btagdown.root");
    filebta_norw      = new TFile("HistoBDToutput/TMVApp_kut_noreweight.root");
  
    filebta_PDFup     = new TFile("HistoBDToutput/TMVApp_kut_PDFup.root");
    filebta_PDFdown   = new TFile("HistoBDToutput/TMVApp_kut_PDFdown.root");
  
    filebta_JESup     = new TFile("HistoBDToutput/TMVApp_kut_JESup.root");
    filebta_JESdown   = new TFile("HistoBDToutput/TMVApp_kut_JESdown.root");
  
    filebta_Matchup   = new TFile("HistoBDToutput/TMVApp_kut_Matchup.root");
    filebta_Matchdown = new TFile("HistoBDToutput/TMVApp_kut_Matchdown.root");
  
    filebta_Scaleup   = new TFile("HistoBDToutput/TMVApp_kut_Scaleup.root");
    filebta_Scaledown = new TFile("HistoBDToutput/TMVApp_kut_Scaledown.root");
  }


  
  if(theVertex == "kct"&& doSyst){
  
    filebta_btagup    = new TFile("HistoBDToutput/TMVApp_kct_btagUp.root");
    filebta_btagdown  = new TFile("HistoBDToutput/TMVApp_kct_btagdown.root");
    filebta_norw      = new TFile("HistoBDToutput/TMVApp_kct_noreweight.root");
  
    filebta_PDFup     = new TFile("HistoBDToutput/TMVApp_kct_PDFup.root");
    filebta_PDFdown   = new TFile("HistoBDToutput/TMVApp_kct_PDFdown.root");
  
    filebta_JESup     = new TFile("HistoBDToutput/TMVApp_kct_JESup.root");
    filebta_JESdown   = new TFile("HistoBDToutput/TMVApp_kct_JESdown.root");
  
    filebta_Matchup   = new TFile("HistoBDToutput/TMVApp_kct_Matchup.root");
    filebta_Matchdown = new TFile("HistoBDToutput/TMVApp_kct_Matchdown.root");
  
    filebta_Scaleup   = new TFile("HistoBDToutput/TMVApp_kct_Scaleup.root");
    filebta_Scaledown = new TFile("HistoBDToutput/TMVApp_kct_Scaledown.root");
  }

  cout << "4324 " << endl;
   
  if(doStat){
    for (int ierr=1; ierr<=herrorband->GetNbinsX(); ierr++) {
    double error_all = pow(histBdt_WZ->GetBinError(ierr), 2) // stat
                        + pow(histBdt_WZ->GetBinContent(ierr)*theWZscale_err/theWZscale, 2) // syst SF 
                        + pow(histBdt_DataZjets->GetBinError(ierr), 2) // stat
                        + pow(histBdt_DataZjets->GetBinContent(ierr)*0.30, 2) // syst SF
                        + pow(histBdt_ZZ->GetBinContent(ierr)*0.30, 2) // syst SF 
                        + pow(histBdt_TTbar->GetBinContent(ierr)*2.60, 2) // syst SF 
                        + pow(histBdt_TZq->GetBinContent(ierr)*0.50, 2) // syst SF 
                        + pow(
				(histBdt_ZZ->GetBinContent(ierr)+histBdt_TTbar->GetBinContent(ierr)+histBdt_TZq->GetBinContent(ierr))
				
				*0.025, 2); // syst lumi
				//cout << "error_all 1 " << error_all << endl;
				//cout << "theVariable " << theVariable << endl;
			if(theVariable == "NBJets" || theVariable == "NBJets_BDTcut" )	{
			   cout << "in " << endl;
			  if(ierr == 1) error_all=error_all+pow(histBdt_WZ->GetBinContent(ierr)*0.11, 2); 
			  if(ierr == 2) error_all=error_all+pow(histBdt_WZ->GetBinContent(ierr)*0.41, 2); 
			
			
			}
			
				//cout << "error_all 2 " << error_all << endl;	
    herrorband->SetBinError(ierr, sqrt(error_all));
    
    }
  }
  if(doSyst){
    cout << "starts do syst" << endl;
    //------------------------------------------
    //adding btag systematics
    //------------------------------------------
     
    TH1F * histBdt_WZ_btagup        = (TH1F*)((TH1F*)filebta_btagup->Get( (theVariable+"_WZ"       ).Data()))->Clone();
    TH1F * histBdt_ZZ_btagup        = (TH1F*)((TH1F*)filebta_btagup->Get( (theVariable+"_ZZ"       ).Data()))->Clone();
    TH1F * histBdt_Zjets_btagup     = (TH1F*)((TH1F*)filebta_btagup->Get( (theVariable+"_Zjets"    ).Data()))->Clone(); 
    TH1F * histBdt_TTbar_btagup     = (TH1F*)((TH1F*)filebta_btagup->Get( (theVariable+"_TTbarSig" ).Data()))->Clone();
    TH1F * histBdt_TZq_btagup       = (TH1F*)((TH1F*)filebta_btagup->Get( (theVariable+"_TZq"      ).Data()))->Clone();
    
    histBdt_WZ_btagup->Rebin(rebin);   
    histBdt_ZZ_btagup->Rebin(rebin);   
    histBdt_Zjets_btagup->Rebin(rebin);
    histBdt_TTbar_btagup->Rebin(rebin);
    histBdt_TZq_btagup->Rebin(rebin);  
    
    histBdt_WZ_btagup->   Scale((histBdt_WZ->Integral() - histBdt_ZZ->Integral())/histBdt_WZ_btagup->Integral());
    
    
    
    cout << "******************************* " << endl;
    herrorband_btagup =  (TH1F*) histBdt_WZ_btagup->Clone();
    
    if(!doFitZ){
      herrorband_btagup->Add(histBdt_DataZjets);
      herrorband_btagup->Add(histBdt_TTbar_btagup);
      herrorband_btagup->Add(histBdt_TZq_btagup);
      herrorband_btagup->Add(histBdt_ZZ_btagup);
    }
    
    TH1F * histBdt_WZ_btagdown        = (TH1F*)((TH1F*)filebta_btagdown->Get( (theVariable+"_WZ"       ).Data()))->Clone();
    TH1F * histBdt_ZZ_btagdown        = (TH1F*)((TH1F*)filebta_btagdown->Get( (theVariable+"_ZZ"       ).Data()))->Clone();
    TH1F * histBdt_Zjets_btagdown     = (TH1F*)((TH1F*)filebta_btagdown->Get( (theVariable+"_Zjets"    ).Data()))->Clone(); 
    TH1F * histBdt_TTbar_btagdown     = (TH1F*)((TH1F*)filebta_btagdown->Get( (theVariable+"_TTbarSig" ).Data()))->Clone();
    TH1F * histBdt_TZq_btagdown       = (TH1F*)((TH1F*)filebta_btagdown->Get( (theVariable+"_TZq"      ).Data()))->Clone();
    
    histBdt_WZ_btagdown->Rebin(rebin);   
    histBdt_ZZ_btagdown->Rebin(rebin);   
    histBdt_Zjets_btagdown->Rebin(rebin);
    histBdt_TTbar_btagdown->Rebin(rebin);
    histBdt_TZq_btagdown->Rebin(rebin);  
    
    histBdt_WZ_btagdown->Scale((histBdt_WZ->Integral() - histBdt_ZZ->Integral())/histBdt_WZ_btagdown->Integral());
    //histBdt_TTbar_btagdown->Scale(histBdt_TTbar->Integral()/histBdt_TTbar_btagdown->Integral());
    //histBdt_TZq_btagdown->Scale(histBdt_TZq->Integral()/histBdt_TZq_btagdown->Integral());
    
    
    
    herrorband_btagdown =  (TH1F*) histBdt_WZ_btagdown->Clone();
    if(!doFitZ){
      herrorband_btagdown->Add(histBdt_DataZjets);
      herrorband_btagdown->Add(histBdt_TTbar_btagdown);
      herrorband_btagdown->Add(histBdt_TZq_btagdown);
      herrorband_btagdown->Add(histBdt_ZZ_btagdown);
    }
    
    
    /*TH1F * histBdt_WZ_norw        = (TH1F*)filebta_norw->Get( (theVariable+"_WZ"       ).Data());
    TH1F * histBdt_ZZ_norw        = (TH1F*)filebta_norw->Get( (theVariable+"_ZZ"       ).Data());
    TH1F * histBdt_Zjets_norw     = (TH1F*)filebta_norw->Get( (theVariable+"_Zjets"    ).Data()); 
    TH1F * histBdt_TTbar_norw     = (TH1F*)filebta_norw->Get( (theVariable+"_TTbarSig" ).Data());
    TH1F * histBdt_TZq_norw       = (TH1F*)filebta_norw->Get( (theVariable+"_TZq"      ).Data());
    
    
    herrorband_norw =  (TH1F*) histBdt_WZ_norw->Clone();
    herrorband_norw->Add(histBdt_DataZjets);
    herrorband_norw->Add(histBdt_TTbar_norw);
    herrorband_norw->Add(histBdt_TZq_norw);
    */
    
    //------------------------------------------
    //adding PDF systematics
    //------------------------------------------
     
    TH1F * histBdt_WZ_PDFup        = (TH1F*)((TH1F*)filebta_PDFup->Get( (theVariable+"_WZ"       ).Data()))->Clone();
    TH1F * histBdt_ZZ_PDFup        = (TH1F*)((TH1F*)filebta_PDFup->Get( (theVariable+"_ZZ"       ).Data()))->Clone();
    TH1F * histBdt_Zjets_PDFup     = (TH1F*)((TH1F*)filebta_PDFup->Get( (theVariable+"_Zjets"    ).Data()))->Clone(); 
    TH1F * histBdt_TTbar_PDFup     = (TH1F*)((TH1F*)filebta_PDFup->Get( (theVariable+"_TTbarSig" ).Data()))->Clone();
    TH1F * histBdt_TZq_PDFup       = (TH1F*)((TH1F*)filebta_PDFup->Get( (theVariable+"_TZq"      ).Data()))->Clone();
    
    histBdt_WZ_PDFup->Rebin(rebin);   
    histBdt_ZZ_PDFup->Rebin(rebin);   
    histBdt_Zjets_PDFup->Rebin(rebin);
    histBdt_TTbar_PDFup->Rebin(rebin);
    histBdt_TZq_PDFup->Rebin(rebin); 
    
    histBdt_WZ_PDFup->Scale((histBdt_WZ->Integral() - histBdt_ZZ->Integral())/histBdt_WZ_PDFup->Integral());
    
    
    
    herrorband_PDFup =  (TH1F*) histBdt_WZ_PDFup->Clone();
    if(!doFitZ){
      herrorband_PDFup->Add(histBdt_DataZjets);
      herrorband_PDFup->Add(histBdt_TTbar_PDFup);
      herrorband_PDFup->Add(histBdt_TZq_PDFup);
      herrorband_PDFup->Add(histBdt_ZZ_PDFup);
    }
    
    
    
    TH1F * histBdt_WZ_PDFdown        = (TH1F*)((TH1F*)filebta_PDFdown->Get( (theVariable+"_WZ"       ).Data()))->Clone();
    TH1F * histBdt_ZZ_PDFdown        = (TH1F*)((TH1F*)filebta_PDFdown->Get( (theVariable+"_ZZ"       ).Data()))->Clone();
    TH1F * histBdt_Zjets_PDFdown     = (TH1F*)((TH1F*)filebta_PDFdown->Get( (theVariable+"_Zjets"    ).Data()))->Clone(); 
    TH1F * histBdt_TTbar_PDFdown     = (TH1F*)((TH1F*)filebta_PDFdown->Get( (theVariable+"_TTbarSig" ).Data()))->Clone();
    TH1F * histBdt_TZq_PDFdown       = (TH1F*)((TH1F*)filebta_PDFdown->Get( (theVariable+"_TZq"      ).Data()))->Clone();
    
    histBdt_WZ_PDFdown->Rebin(rebin);   
    histBdt_ZZ_PDFdown->Rebin(rebin);   
    histBdt_Zjets_PDFdown->Rebin(rebin);
    histBdt_TTbar_PDFdown->Rebin(rebin);
    histBdt_TZq_PDFdown->Rebin(rebin); 
    
    histBdt_WZ_PDFdown->Scale((histBdt_WZ->Integral() - histBdt_ZZ->Integral())/histBdt_WZ_PDFdown->Integral());
    
    herrorband_PDFdown =  (TH1F*) histBdt_WZ_PDFdown->Clone();
    if(!doFitZ){
      herrorband_PDFdown->Add(histBdt_DataZjets);
      herrorband_PDFdown->Add(histBdt_TTbar_PDFdown);
      herrorband_PDFdown->Add(histBdt_TZq_PDFdown);
      herrorband_PDFdown->Add(histBdt_ZZ_PDFdown);
    }
    
    //------------------------------------------
    //adding JES systematics
    //------------------------------------------
    
     
    TH1F * histBdt_WZ_JESup        = (TH1F*)((TH1F*)filebta_JESup->Get( (theVariable+"_WZ"       ).Data()))->Clone();
    TH1F * histBdt_ZZ_JESup        = (TH1F*)((TH1F*)filebta_JESup->Get( (theVariable+"_ZZ"       ).Data()))->Clone();
    TH1F * histBdt_Zjets_JESup     = (TH1F*)((TH1F*)filebta_JESup->Get( (theVariable+"_Zjets"    ).Data()))->Clone(); 
    TH1F * histBdt_TTbar_JESup     = (TH1F*)((TH1F*)filebta_JESup->Get( (theVariable+"_TTbarSig" ).Data()))->Clone();
    TH1F * histBdt_TZq_JESup       = (TH1F*)((TH1F*)filebta_JESup->Get( (theVariable+"_TZq"      ).Data()))->Clone();
       
    histBdt_WZ_JESup->Rebin(rebin);   
    histBdt_ZZ_JESup->Rebin(rebin);   
    histBdt_Zjets_JESup->Rebin(rebin);
    histBdt_TTbar_JESup->Rebin(rebin);
    histBdt_TZq_JESup->Rebin(rebin); 

    histBdt_WZ_JESup->Scale((histBdt_WZ->Integral() - histBdt_ZZ->Integral())/histBdt_WZ_JESup->Integral());
    
    herrorband_JESup =  (TH1F*) histBdt_WZ_JESup->Clone();
    if(!doFitZ){
      herrorband_JESup->Add(histBdt_DataZjets);
      herrorband_JESup->Add(histBdt_TTbar_JESup);
      herrorband_JESup->Add(histBdt_TZq_JESup);
      herrorband_JESup->Add(histBdt_ZZ_JESup);
    }
    
    TH1F * histBdt_WZ_JESdown        = (TH1F*)((TH1F*)filebta_JESdown->Get( (theVariable+"_WZ"       ).Data()))->Clone();
    TH1F * histBdt_ZZ_JESdown        = (TH1F*)((TH1F*)filebta_JESdown->Get( (theVariable+"_ZZ"       ).Data()))->Clone();
    TH1F * histBdt_Zjets_JESdown     = (TH1F*)((TH1F*)filebta_JESdown->Get( (theVariable+"_Zjets"    ).Data()))->Clone(); 
    TH1F * histBdt_TTbar_JESdown     = (TH1F*)((TH1F*)filebta_JESdown->Get( (theVariable+"_TTbarSig" ).Data()))->Clone();
    TH1F * histBdt_TZq_JESdown       = (TH1F*)((TH1F*)filebta_JESdown->Get( (theVariable+"_TZq"      ).Data()))->Clone();
    
    histBdt_WZ_JESdown->Rebin(rebin);   
    histBdt_ZZ_JESdown->Rebin(rebin);   
    histBdt_Zjets_JESdown->Rebin(rebin);
    histBdt_TTbar_JESdown->Rebin(rebin);
    histBdt_TZq_JESdown->Rebin(rebin); 

    histBdt_WZ_JESdown->Scale((histBdt_WZ->Integral() - histBdt_ZZ->Integral())/histBdt_WZ_JESdown->Integral());
    
    herrorband_JESdown =  (TH1F*) histBdt_WZ_JESdown->Clone();
    if(!doFitZ){
      herrorband_JESdown->Add(histBdt_DataZjets);
      herrorband_JESdown->Add(histBdt_TTbar_JESdown);
      herrorband_JESdown->Add(histBdt_TZq_JESdown);
      herrorband_JESdown->Add(histBdt_ZZ_JESdown);
    }
    
    //------------------------------------------
    //adding Match systematics
    //------------------------------------------
    
     
    TH1F * histBdt_WZ_Matchup        = (TH1F*)((TH1F*)filebta_Matchup->Get( (theVariable+"_WZ"       ).Data()))->Clone();
    TH1F * histBdt_ZZ_Matchup        = (TH1F*)((TH1F*)filebta_Matchup->Get( (theVariable+"_ZZ"       ).Data()))->Clone();
    TH1F * histBdt_Zjets_Matchup     = (TH1F*)((TH1F*)filebta_Matchup->Get( (theVariable+"_Zjets"    ).Data()))->Clone(); 
    TH1F * histBdt_TTbar_Matchup     = (TH1F*)((TH1F*)filebta_Matchup->Get( (theVariable+"_TTbarSig" ).Data()))->Clone();
    TH1F * histBdt_TZq_Matchup       = (TH1F*)((TH1F*)filebta_Matchup->Get( (theVariable+"_TZq"      ).Data()))->Clone();
      
    histBdt_WZ_Matchup->Rebin(rebin);   
    histBdt_ZZ_Matchup->Rebin(rebin);   
    histBdt_Zjets_Matchup->Rebin(rebin);
    histBdt_TTbar_Matchup->Rebin(rebin);
    histBdt_TZq_Matchup->Rebin(rebin); 

    
    histBdt_WZ_Matchup->Scale((histBdt_WZ->Integral() - histBdt_ZZ->Integral())/histBdt_WZ_Matchup->Integral());
    
    herrorband_Matchup =  (TH1F*) histBdt_WZ_Matchup->Clone();
    if(!doFitZ){
      herrorband_Matchup->Add(histBdt_DataZjets);
      herrorband_Matchup->Add(histBdt_TTbar);
      herrorband_Matchup->Add(histBdt_TZq);
      herrorband_Matchup->Add(histBdt_ZZ_Matchup);
    }
    
    TH1F * histBdt_WZ_Matchdown        = (TH1F*)((TH1F*)filebta_Matchdown->Get( (theVariable+"_WZ"       ).Data()))->Clone();
    TH1F * histBdt_ZZ_Matchdown        = (TH1F*)((TH1F*)filebta_Matchdown->Get( (theVariable+"_ZZ"       ).Data()))->Clone();
    TH1F * histBdt_Zjets_Matchdown     = (TH1F*)((TH1F*)filebta_Matchdown->Get( (theVariable+"_Zjets"    ).Data()))->Clone(); 
    TH1F * histBdt_TTbar_Matchdown     = (TH1F*)((TH1F*)filebta_Matchdown->Get( (theVariable+"_TTbarSig" ).Data()))->Clone();
    TH1F * histBdt_TZq_Matchdown       = (TH1F*)((TH1F*)filebta_Matchdown->Get( (theVariable+"_TZq"      ).Data()))->Clone();
    
    histBdt_WZ_Matchdown->Rebin(rebin);   
    histBdt_ZZ_Matchdown->Rebin(rebin);   
    histBdt_Zjets_Matchdown->Rebin(rebin);
    histBdt_TTbar_Matchdown->Rebin(rebin);
    histBdt_TZq_Matchdown->Rebin(rebin); 
    
    //cout << "histBdt_WZ_Matchdown->Integral() " << histBdt_WZ_Matchdown->Integral() << endl;
    histBdt_WZ_Matchdown->Scale((histBdt_WZ->Integral() - histBdt_ZZ->Integral())/histBdt_WZ_Matchdown->Integral());
    //cout << "histBdt_WZ->Integral() " << histBdt_WZ->Integral() << endl;
    //cout << "histBdt_WZ_Matchdown->Integral() " << histBdt_WZ_Matchdown->Integral() << endl;
    
    herrorband_Matchdown =  (TH1F*) histBdt_WZ_Matchdown->Clone();
    if(!doFitZ){
      herrorband_Matchdown->Add(histBdt_DataZjets);
      herrorband_Matchdown->Add(histBdt_TTbar);
      herrorband_Matchdown->Add(histBdt_TZq);
      herrorband_Matchdown->Add(histBdt_ZZ_Matchdown);
    }  
    
    
    //------------------------------------------
    //adding Scale systematics
    //------------------------------------------
    
     
    TH1F * histBdt_WZ_Scaleup        = (TH1F*)((TH1F*)filebta_Scaleup->Get( (theVariable+"_WZ"       ).Data()))->Clone();
    TH1F * histBdt_ZZ_Scaleup        = (TH1F*)((TH1F*)filebta_Scaleup->Get( (theVariable+"_ZZ"       ).Data()))->Clone();
    TH1F * histBdt_Zjets_Scaleup     = (TH1F*)((TH1F*)filebta_Scaleup->Get( (theVariable+"_Zjets"    ).Data()))->Clone(); 
    TH1F * histBdt_TTbar_Scaleup     = (TH1F*)((TH1F*)filebta_Scaleup->Get( (theVariable+"_TTbarSig" ).Data()))->Clone();
    TH1F * histBdt_TZq_Scaleup       = (TH1F*)((TH1F*)filebta_Scaleup->Get( (theVariable+"_TZq"      ).Data()))->Clone();
    
    histBdt_WZ_Scaleup->Rebin(rebin);   
    histBdt_ZZ_Scaleup->Rebin(rebin);   
    histBdt_Zjets_Scaleup->Rebin(rebin);
    histBdt_TTbar_Scaleup->Rebin(rebin);
    histBdt_TZq_Scaleup->Rebin(rebin); 
    
    
    histBdt_WZ_Scaleup->Scale((histBdt_WZ->Integral() - histBdt_ZZ->Integral())/histBdt_WZ_Scaleup->Integral());
    
    herrorband_Scaleup =  (TH1F*) histBdt_WZ_Scaleup->Clone();
    if(!doFitZ){
      herrorband_Scaleup->Add(histBdt_DataZjets);
      herrorband_Scaleup->Add(histBdt_TTbar);
      herrorband_Scaleup->Add(histBdt_TZq);
    }
    
    TH1F * histBdt_WZ_Scaledown        = (TH1F*)((TH1F*)filebta_Scaledown->Get( (theVariable+"_WZ"       ).Data()))->Clone();
    TH1F * histBdt_ZZ_Scaledown        = (TH1F*)((TH1F*)filebta_Scaledown->Get( (theVariable+"_ZZ"       ).Data()))->Clone();
    TH1F * histBdt_Zjets_Scaledown     = (TH1F*)((TH1F*)filebta_Scaledown->Get( (theVariable+"_Zjets"    ).Data()))->Clone(); 
    TH1F * histBdt_TTbar_Scaledown     = (TH1F*)((TH1F*)filebta_Scaledown->Get( (theVariable+"_TTbarSig" ).Data()))->Clone();
    TH1F * histBdt_TZq_Scaledown       = (TH1F*)((TH1F*)filebta_Scaledown->Get( (theVariable+"_TZq"      ).Data()))->Clone();


     cout << "check nbin **********************  " << endl;
     cout << "histBdt_WZ_btagup->GetNbinsX()    " << histBdt_WZ_btagup->GetNbinsX()   << endl;
     cout << "histBdt_WZ_btagdown->GetNbinsX()  " << histBdt_WZ_btagdown->GetNbinsX() << endl;
     cout << "histBdt_WZ_PDFup ->GetNbinsX()    " << histBdt_WZ_PDFup->GetNbinsX()    << endl;
     cout << "histBdt_WZ_PDFdown->GetNbinsX()   " << histBdt_WZ_PDFdown->GetNbinsX()  << endl;
     cout << "histBdt_WZ_JESup->GetNbinsX()     " << histBdt_WZ_JESup->GetNbinsX()    << endl;
     cout << "histBdt_WZ_JESdown->GetNbinsX()   " << histBdt_WZ_JESdown->GetNbinsX()  << endl;
     cout << "histBdt_WZ_Matchup->GetNbinsX()   " << histBdt_WZ_Matchup->GetNbinsX()  << endl;
     cout << "histBdt_WZ_Matchdown->GetNbinsX() " << histBdt_WZ_Matchdown->GetNbinsX() << endl;
     cout << "histBdt_WZ_Scaleup->GetNbinsX()   " << histBdt_WZ_Scaleup->GetNbinsX()   << endl;
     cout << "histBdt_WZ_Scaledown->GetNbinsX() " << histBdt_WZ_Scaledown ->GetNbinsX()<< endl;
    
    
    
    /*histBdt_WZ_Scaledown->Rebin(rebin);   
    histBdt_ZZ_Scaledown->Rebin(rebin);   
    histBdt_Zjets_Scaledown->Rebin(rebin);
    histBdt_TTbar_Scaledown->Rebin(rebin);
    histBdt_TZq_Scaledown->Rebin(rebin); 
    */
    histBdt_WZ_Scaledown->Scale((histBdt_WZ->Integral() - histBdt_ZZ->Integral())/histBdt_WZ_Scaledown->Integral());
    
    
    
    herrorband_Scaledown =  (TH1F*) histBdt_WZ_Scaledown->Clone();
    if(!doFitZ){
      herrorband_Scaledown->Add(histBdt_DataZjets);
      herrorband_Scaledown->Add(histBdt_TTbar);
      herrorband_Scaledown->Add(histBdt_TZq);
      herrorband_Scaledown->Add(histBdt_ZZ_Scaledown);
    }
      cout << "end do syst" << endl;
      
      
     cout << "check pointer ********************* " << endl;
     cout << "histBdt_WZ_btagdown  " << histBdt_WZ_btagdown << endl;
     cout << "histBdt_WZ_PDFup     " << histBdt_WZ_PDFup << endl;
     cout << "histBdt_WZ_PDFdown   " << histBdt_WZ_PDFdown << endl;
     cout << "histBdt_WZ_JESup     " << histBdt_WZ_JESup << endl;
     cout << "histBdt_WZ_JESdown   " << histBdt_WZ_JESdown << endl;
     cout << "histBdt_WZ_Matchup   " << histBdt_WZ_Matchup << endl;
     cout << "histBdt_WZ_Matchdown " << histBdt_WZ_Matchdown << endl;
     cout << "histBdt_WZ_Scaleup   " << histBdt_WZ_Scaleup << endl;
     cout << "histBdt_WZ_Scaledown " << histBdt_WZ_Scaledown << endl;
     
     
     
    TH1F * histBdt_DataZjets_DYUp        = (TH1F*)((TH1F*)filebta_DYup->Get( (theVariable+"_DataZjets"       ).Data()))->Clone();
    TH1F * histBdt_DataZjets_DYDown      = (TH1F*)((TH1F*)filebta_DYdown->Get( (theVariable+"_DataZjets"       ).Data()))->Clone();
    cout << "808" << endl;
    histBdt_DataZjets_DYUp->Scale(histBdt_DataZjets->Integral()/histBdt_DataZjets_DYUp->Integral());  
    histBdt_DataZjets_DYDown ->Scale(histBdt_DataZjets->Integral()/histBdt_DataZjets_DYDown->Integral());
     cout << "811" << endl;
    
      for(unsigned int i=1; i<herrorband->GetNbinsX(); i++){
  
     
        double newValue = histBdt_DataZjets->GetBinContent(i)*
     	pow(histBdt_DataZjets_DYDown->GetBinContent(i)/histBdt_DataZjets->GetBinContent(i), 0.413819);
     
        histBdt_DataZjets->SetBinContent(i, newValue);
  
  }

     
     
  }
  
  
  
  
  cout << "755 " << endl;
  
 
  TGraphAsymmErrors *thegraph = new TGraphAsymmErrors(herrorband);
  thegraph->SetFillStyle(3005);
  thegraph->SetFillColor(1);
  
  cout << "***************************" << endl;
  if(doSyst) cout << "herrorband_btagup->Integral() " << herrorband_btagup->Integral() << endl;
  cout << "herrorband->Integral()        " << herrorband->Integral() << endl;
  cout << "***************************" << endl;
  
  
  for(int i=1; i<=herrorband->GetNbinsX(); i++){
    
      double theYup_error  = 0;
      double theYdown_error = 0;
    
      double thevar=0;

    
    if(doSyst){
    
      theYup_error = pow(thegraph->GetErrorY(i)/2., 2);
      theYdown_error = pow(thegraph->GetErrorY(i)/2., 2);
    
      double thevar=0;
      //cout << "**********************   i"<< i << endl;
      //cout << "thevar 1 " << thevar << endl;
      
      
      if(herrorband_btagup->GetNbinsX() != herrorband->GetNbinsX()) cout << "different binning " << endl;
      thevar = (herrorband_btagup->GetBinContent(i+1)-herrorband->GetBinContent(i+1));
      if(thevar > 0) theYup_error  +=pow(thevar, 2);
      else theYdown_error+=pow( thevar, 2);
      //cout << "thevar 2 " << thevar << endl;
    
    
      if(herrorband_PDFup->GetNbinsX() != herrorband->GetNbinsX()) cout << "different binning " << endl;
      thevar = (herrorband_PDFup->GetBinContent(i+1)-herrorband->GetBinContent(i+1));
      if(thevar > 0) theYup_error  +=pow(thevar, 2);
      else theYdown_error+=pow( thevar, 2); 
      //cout << "thevar 3 " << thevar << endl;
      
      if(herrorband_JESup->GetNbinsX() != herrorband->GetNbinsX()) cout << "different binning " << endl;
      thevar = (herrorband_JESup->GetBinContent(i+1)-herrorband->GetBinContent(i+1));
      if(thevar > 0) theYup_error  +=pow(thevar, 2);
      else theYdown_error+=pow( thevar, 2); 
      //cout << "thevar 4 " << thevar << endl;
      
      
      
      if(herrorband_Matchup->GetNbinsX() != herrorband->GetNbinsX()) cout << "different binning " << endl;
      thevar = (herrorband_Matchup->GetBinContent(i+1)-herrorband->GetBinContent(i+1));
      if(thevar > 0) theYup_error  +=pow(thevar, 2);
      else theYdown_error+=pow( thevar, 2);
      //cout << "thevar 5 " << thevar << endl;
     
      if(herrorband_Scaleup->GetNbinsX() != herrorband->GetNbinsX()) cout << "different binning " << endl;
      thevar = (herrorband_Scaleup->GetBinContent(i+1)-herrorband->GetBinContent(i+1));
      if(thevar > 0) theYup_error  +=pow(thevar, 2);
      else theYdown_error+=pow( thevar, 2);
      //cout << "thevar 6 " << thevar << endl;
  
    
    
    
    
      if(herrorband_btagdown->GetNbinsX() != herrorband->GetNbinsX()) cout << "different binning " << endl;
      thevar = (herrorband_btagdown->GetBinContent(i+1)-herrorband->GetBinContent(i+1)) ;
      if(thevar > 0) theYup_error  +=pow(thevar, 2);
      else theYdown_error+=pow( thevar, 2);
      //cout << "thevar 7 " << thevar << endl;
      
      if(herrorband_PDFdown->GetNbinsX() != herrorband->GetNbinsX()) cout << "different binning " << endl;
      thevar = (herrorband_PDFdown->GetBinContent(i+1)-herrorband->GetBinContent(i+1)) ;
      if(thevar > 0) theYup_error  +=pow(thevar, 2);
      else theYdown_error+=pow( thevar, 2); 
      //cout << "thevar 8 " << thevar << endl;
      
      if(herrorband_JESdown->GetNbinsX() != herrorband->GetNbinsX()) cout << "different binning " << endl;
      thevar = (herrorband_JESdown->GetBinContent(i+1)-herrorband->GetBinContent(i+1))  ;
      if(thevar > 0) theYup_error  +=pow(thevar, 2);
      else theYdown_error+=pow( thevar, 2);
      //cout << "thevar 9 " << thevar << endl;
      
      if(herrorband_Matchdown->GetNbinsX() != herrorband->GetNbinsX()) cout << "different binning " << endl;
      thevar = (herrorband_Matchdown->GetBinContent(i+1)-herrorband->GetBinContent(i+1));
      if(thevar > 0) theYup_error  +=pow(thevar, 2);
      else theYdown_error+=pow( thevar, 2);
      //cout << "thevar 10 " << thevar << endl;
      
      if(herrorband_Scaledown->GetNbinsX() != herrorband->GetNbinsX()) cout << "different binning " << endl;
      thevar = (herrorband_Scaledown->GetBinContent(i+1)-herrorband->GetBinContent(i+1));
      if(thevar > 0) theYup_error  +=pow(thevar, 2);
      else theYdown_error+=pow( thevar, 2);
      //cout << "thevar 11 " << thevar << endl;
    
    
     }
     if(doStat){
    
        thevar = herrorband->GetBinError(i);
        theYup_error  +=pow(thevar, 2);
        theYdown_error+=pow( thevar, 2);
      
     }
    
    
    
    
    
    
    thegraph->SetPointEYhigh( i-1, pow(theYup_error, 0.5));
    thegraph->SetPointEYlow( i-1, pow(theYdown_error, 0.5));
    
    
    //cout << "tot theYup_error   " <<  pow(theYup_error, 0.5)<< endl;
    //cout << "tot theYdown_error " << pow(theYdown_error, 0.5) << endl;

	
  }
  
  cout << "880" << endl;
  
  thegraph->Draw("e2same");
  //thegraph->Draw("lsame");
  /*if(doBtagUncerOnly) {
    herrorband_norw->SetLineWidth(2);
    herrorband_norw->Draw("hsame");
  }*/
  
  //if(doSyst) herrorband_btagup->SetLineColor(2);
  //herrorband_btagup->SetLineWidth(2);
  
  //herrorband_btagup->Draw("same");
  herrorband->SetMarkerStyle(21) ;
  herrorband->SetMarkerSize(1.2) ;
  //herrorband->Draw("epsame");
  
  
  //herrorband_Scaledown->Draw();
  //herrorband->Draw("same");
  
  
  TH1F * histo_ratio = (TH1F*) histBdt_Data->Clone();
  TH1F * histo_mc    = (TH1F*) histBdt_WZ->Clone();
  if(!doFitZ){
    histo_mc->Add(histBdt_DataZjets);
    histo_mc->Add(histBdt_TTbar);
    histo_mc->Add(histBdt_TZq);
  }
  
  
  histo_ratio->Sumw2();
  histo_mc->Sumw2();
  
  histo_ratio->Divide(histo_ratio, histo_mc, 1, 1, "b");
  
  /*for(unsigned int i=1; i<=histBdt_WZ->GetNbinsX(); i++){
    
    histo_ratio->SetBinError(i, 10000*thegraph->GetErrorY(i-1));
    
  }*/  
  for (int ierr=1; ierr<=histo_ratio->GetNbinsX(); ierr++) {
    double denom         = histo_mc ->GetBinContent(ierr);
    double num     = histBdt_Data->GetBinContent(ierr);
    if(denom <= 0 || num <= 0) histo_ratio->SetBinContent(ierr, 0);
    if(denom <= 0 || num <= 0) histo_ratio->SetBinError(ierr, 1000);
  }
  
  
  TGraphAsymmErrors *theRatio = new TGraphAsymmErrors(histo_ratio);
  
  
  cout << "926" << endl;
  
  for (int ierr=1; ierr<=histo_ratio->GetNbinsX(); ierr++) {
    
    double num     = histBdt_Data->GetBinContent(ierr);
    double num_err = histBdt_Data->GetBinError(ierr);
    
    
    //if(histBdt_Data->GetBinContent(ierr) == -1) num_err = 1000000;
    
    double denom         = histo_mc ->GetBinContent(ierr);
    double denom_err     = herrorband->GetBinError(ierr);
    double denom_errUp   = thegraph->GetErrorYhigh(ierr-1);
    double denom_errDown = thegraph->GetErrorYlow(ierr-1);
    
    
    if(histBdt_Data->GetBinContent(ierr) == 1 ) num_err = 3;
    
    double errorUp = pow(
       pow(num_err/denom, 2)+
       pow(num*denom_errUp/(denom*denom), 2)
    , 0.5);
    
    
    double errorDown = pow(
       pow(num_err/denom, 2)+
       pow(num*denom_errDown/(denom*denom), 2)
    , 0.5);
    
    //cout << "errorUp   " << errorUp << endl;
    //cout << "errorDown " <<  errorDown<< endl;
    
    
    //histo_ratio->SetBinError(ierr, error );
    theRatio->SetPointEYhigh( ierr-1, errorDown);
    theRatio->SetPointEYlow(  ierr-1, errorUp);
    
    //cout << "get error check " << theRatio->GetErrorYhigh(ierr-1) << endl;
    //cout << "get value check " << histo_ratio->GetBinContent(ierr) << endl;
    
    cout << "********* i "  << (float(ierr)/histo_ratio->GetNbinsX())*2. -1 << endl;
    cout << "num " << num  << " pm " << num_err << endl;
    cout << "denom " <<  denom << " pm " << denom_errUp << endl;
    cout << "denom " <<  denom << " pm " << denom_errDown << endl;
    cout << "ratio " <<  num/denom << endl;
    
    histo_ratio->SetBinError(ierr-1, errorUp);
    
    //if(num > 0) histo_ratio->SetBinError(ierr-1, errorUp);
    //else histo_ratio->SetBinError(ierr-1, 1000);
    
  }
  
  cout << "970" << endl;
  /*
  cout << theVariable << endl;
  cout << theVariable << endl;
  cout << theVariable << endl;
  cout << theVariable << endl;
  cout << theVariable << endl;
  cout << theVariable << endl;
  cout << theVariable << endl;
  cout << theVariable << endl;
  cout << theVariable << endl;
  cout << theVariable << endl;
  cout << theVariable << endl;
  cout << theVariable << endl;
  cout << theVariable << endl;
  cout << theVariable << endl;
  */
  if(theVariable == "MVA_BDT"       ) histo_ratio->GetXaxis()->SetTitle("BDT output");
  if(theVariable == "topMass"       ) histo_ratio->GetXaxis()->SetTitle("m_{top} [GeV]");
  if(theVariable == "totMass"       ) histo_ratio->GetXaxis()->SetTitle("m_{tot} [GeV]");
  if(theVariable == "deltaPhilb"    ) histo_ratio->GetXaxis()->SetTitle("#Delta #phi(l_{W}-b)");
  if(theVariable == "deltaRlb"      ) histo_ratio->GetXaxis()->SetTitle("#Delta R(l_{W}-b)");
  if(theVariable == "deltaRTopZ"    ) histo_ratio->GetXaxis()->SetTitle("#Delta R(Z-t)");
  if(theVariable == "asym"          ) histo_ratio->GetXaxis()->SetTitle("q|#eta|");
  if(theVariable == "Zpt"           ) histo_ratio->GetXaxis()->SetTitle("p_{T}(Z) [GeV]");
  if(theVariable == "ZEta"          ) histo_ratio->GetXaxis()->SetTitle("#eta(Z)");
  if(theVariable == "topPt"         ) histo_ratio->GetXaxis()->SetTitle("p_{T}(t) [GeV]");
  if(theVariable == "topEta"        ) histo_ratio->GetXaxis()->SetTitle("#eta(t)");
  if(theVariable == "NJets"         ) histo_ratio->GetXaxis()->SetTitle("N_{jets}");
  if(theVariable == "NBJets"        ) histo_ratio->GetXaxis()->SetTitle("N_{bjets}");
  if(theVariable == "deltaRZl"      ) histo_ratio->GetXaxis()->SetTitle("#Delta R(l_{W}-Z)");
  if(theVariable == "deltaPhiZmet"  ) histo_ratio->GetXaxis()->SetTitle("#Delta #phi(Z-met)");
  if(theVariable == "btagDiscri"    ) histo_ratio->GetXaxis()->SetTitle("CSV discr.");
  if(theVariable == "leptWPt"       ) histo_ratio->GetXaxis()->SetTitle("p_{T}(l_{W}) [GeV]");
  if(theVariable == "leptWEta"      ) histo_ratio->GetXaxis()->SetTitle("#eta(l_{W})");
  if(theVariable == "leadJetPt"     ) histo_ratio->GetXaxis()->SetTitle("Leading jet p_{T} [GeV]");
  if(theVariable == "leadJetEta"    ) histo_ratio->GetXaxis()->SetTitle("Leading jet #eta");
  if(theVariable == "deltaPhiZleptW") histo_ratio->GetXaxis()->SetTitle("#Delta #phi(Z-l_{W})");
  
  
  
  
  if(theVariable == "topMass_BDTcut"       ) histo_ratio->GetXaxis()->SetTitle("m_{top} [GeV]");
  if(theVariable == "totMass_BDTcut"       ) histo_ratio->GetXaxis()->SetTitle("m_{tot} [GeV]");
  if(theVariable == "deltaPhilb_BDTcut"    ) histo_ratio->GetXaxis()->SetTitle("#Delta #phi(l_{W}-b)");
  if(theVariable == "deltaRlb_BDTcut"      ) histo_ratio->GetXaxis()->SetTitle("#Delta R(l_{W}-b)");
  if(theVariable == "deltaRTopZ_BDTcut"    ) histo_ratio->GetXaxis()->SetTitle("#Delta R(Z-t)");
  if(theVariable == "asym_BDTcut"          ) histo_ratio->GetXaxis()->SetTitle("q|#eta|");
  if(theVariable == "Zpt_BDTcut"           ) histo_ratio->GetXaxis()->SetTitle("p_{T}(Z) [GeV]");
  if(theVariable == "ZEta_BDTcut"          ) histo_ratio->GetXaxis()->SetTitle("#eta(Z)");
  if(theVariable == "topPt_BDTcut"         ) histo_ratio->GetXaxis()->SetTitle("p_{T}(t) [GeV]");
  if(theVariable == "topEta_BDTcut"        ) histo_ratio->GetXaxis()->SetTitle("#eta(t)");
  if(theVariable == "NJets_BDTcut"         ) histo_ratio->GetXaxis()->SetTitle("N_{jets}");
  if(theVariable == "NBJets_BDTcut"        ) histo_ratio->GetXaxis()->SetTitle("N_{bjets}");
  if(theVariable == "deltaRZl_BDTcut"      ) histo_ratio->GetXaxis()->SetTitle("#Delta R(l_{W}-Z)");
  if(theVariable == "deltaPhiZmet_BDTcut"  ) histo_ratio->GetXaxis()->SetTitle("#Delta #phi(Z-met)");
  if(theVariable == "btagDiscri_BDTcut"    ) histo_ratio->GetXaxis()->SetTitle("CSV discr.");
  if(theVariable == "leptWPt_BDTcut"       ) histo_ratio->GetXaxis()->SetTitle("p_{T}(l_{W}) [GeV]");
  if(theVariable == "leptWEta_BDTcut"      ) histo_ratio->GetXaxis()->SetTitle("#eta(l_{W})");
  if(theVariable == "leadJetPt_BDTcut"     ) histo_ratio->GetXaxis()->SetTitle("Leading jet p_{T} [GeV]");
  if(theVariable == "leadJetEta_BDTcut"    ) histo_ratio->GetXaxis()->SetTitle("Leading jet #eta");
  if(theVariable == "deltaPhiZleptW_BDTcut") histo_ratio->GetXaxis()->SetTitle("#Delta #phi(Z-l_{W})");
  
  
  
  
  
  
  histo_ratio->GetXaxis()->SetLabelSize(0.07);
  
  
  
  
  
  TPad *canvas_2 = new TPad("canvas_2", "canvas_2", 0.0, 0.0, 1.0, 1.0);
  canvas_2->SetTopMargin(0.7);
  canvas_2->SetFillColor(0);
  canvas_2->SetFillStyle(0);
  canvas_2->SetGridy(1);
  canvas_2->Draw();
  canvas_2->cd(0);
  histo_ratio->SetTitle("");
  
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
  
  histo_ratio->SetMinimum(-1.0);
  histo_ratio->SetMaximum(3.0);
  histo_ratio->SetLineColor(0);
  histo_ratio->Draw("ep");
  
  
  
  
  theRatio->SetMarkerSize(1.2);
  theRatio->Draw("ep");
  
  
  if(doFitZ && theVariable == "Zpt") {
    //theRatio->Fit("pol1");
    histo_ratio->Fit("pol1");
  }
  
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

  TString info_data = "eee, #mu#mu#mu, e#mu#mu, ee#mu channels";; 
  
  if(chan == 0) info_data= "#mu#mu#mu channel";
  if(chan == 1) info_data= "e#mu#mu channel";
  if(chan == 2) info_data= " ee#mu channel";
  if(chan == 3) info_data= "eee channel";
  
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
  
  
  
  
  
  
  
  
  
  
  
  
  TLegend* qw = new TLegend(.80,.70,.95,.90);
  
  
  if(theVariable == "deltaPhilb"     ||
     theVariable == "deltaPhiZleptW" || 
     theVariable == "deltaPhiZmet"   ||
     theVariable == "deltaRlb"	     ||
     theVariable == "deltaRTopZ"     ||
     theVariable == "deltaRZl"	   ) qw = new TLegend(.20,.60,.35,.80);
   
  qw->SetShadowColor(0);
  qw->SetFillColor(0);
  qw->SetLineColor(0);
  qw->AddEntry(histBdt_Data,         "Data" ,                "ep");
  qw->AddEntry(histBdt_WZ,           "VV "                  ,"f");
  qw->AddEntry(histBdt_DataZjets,    "non-prompt Lept. "                  ,"f");
  qw->AddEntry(histBdt_TTbar    ,    "t#bar{t}, single top"                  ,"f");
  qw->AddEntry(histBdt_TZq    ,    "tZq SM"                  ,"f");
  if(!doBtagUncerOnly) qw->AddEntry(histBdt_FCNC,     " tZ 0.1 pb"     ,"l");
  else qw->AddEntry(histBdt_FCNC,     " All Background, no btag reweighting "     ,"l");

  
  qw->SetFillColor(0);
  qw->SetTextFont(42);
  qw->Draw();
  
  if(savePlot){
  c1->SaveAs( ("plots/"+theVariable+"_"+theVertex+".pdf").Data());
  c1->SaveAs( ("plots/"+theVariable+"_"+theVertex+".gif").Data());
  }
  
}
  

void PlotBDToutput(){




   TString thevertex_zut = "zut";
   TString thevertex_zct = "zct";
   TString thevertex_kut = "kut";
   TString thevertex_kct = "kct";
   TString thevertex_xut = "xut";
   TString thevertex_xct = "xct";
   
   TString theFile_zut = "HistoBDToutput/TMVApp_zut_nom.root";
   TString theFile_zct = "HistoBDToutput/TMVApp_zct_nom.root";
   TString theFile_kut = "HistoBDToutput/TMVApp_kut_nom.root";
   TString theFile_kct = "HistoBDToutput/TMVApp_kct_nom.root";
   
    
/*TString theFile_zct = "HistoBDToutput/TMVApp_zct_nom.root";
   TString theFile_kut = "HistoBDToutput/TMVApp_kut_nom.root";
   TString theFile_kct = "HistoBDToutput/TMVApp_kct_nom.root";
   TString theFile_xut = "HistoBDToutput/TMVApp_xut_nom.root";
   TString theFile_xct = "HistoBDToutput/TMVApp_xct_nom.root";*/
   
   
   //PlotBDToutput(thevertex_zut, "NJets"        , theFile_zut); 
   //PlotBDToutput(thevertex_zut, "NBJets"        , theFile_zut); 
   //PlotBDToutput(thevertex_zut, "Zpt"           , theFile_zut);
   
   
   //PlotBDToutput(thevertex_zut, "MVA_BDT"       , theFile_zut); 
   //PlotBDToutput(thevertex_zct, "MVA_BDT"       , theFile_zct);
   //PlotBDToutput(thevertex_kut, "MVA_BDT"       , theFile_kut);
   //PlotBDToutput(thevertex_kct, "MVA_BDT"       , theFile_kct);
   
   //PlotBDToutput(thevertex_zut, "topMass"       , theFile_zut);
   
   if(chan == 0) PlotBDToutput(thevertex_zut, "MVA_BDT"       , "HistoBDToutput/TMVApp_zut_mumumu_nom.root"); 
   if(chan == 1) PlotBDToutput(thevertex_zut, "MVA_BDT"       , "HistoBDToutput/TMVApp_zut_mumue_nom.root"); 
   if(chan == 2) PlotBDToutput(thevertex_zut, "MVA_BDT"       , "HistoBDToutput/TMVApp_zut_eemu_nom.root"); 
   if(chan == 3) PlotBDToutput(thevertex_zut, "MVA_BDT"       , "HistoBDToutput/TMVApp_zut_eee_nom.root"); 
   
   //PlotBDToutput(thevertex_zut, "Zpt"           , "HistoBDToutput/TMVApp_zut_bdtcutnom.root");
   //PlotBDToutput(thevertex_zut, "Zpt"           , "HistoBDToutput/TMVApp_zut_WZZptdown.root");
   //PlotBDToutput(thevertex_kct, "btagDiscri"    , theFile_kct); 
   //PlotBDToutput(thevertex_kct, "NBJets"        , theFile_zut); 
   
   //PlotBDToutput(thevertex_kct, "btagDiscri"    , "HistoBDToutput/TMVApp_zut_noom.root"); 
   //PlotBDToutput(thevertex_kct, "btagDiscri"    , "HistoBDToutput/TMVApp_zut_nobtagcorr.root"); 
   //PlotBDToutput(thevertex_kct, "NBJets"        , "HistoBDToutput/TMVApp_zut_nom.root");
   
   /*
   PlotBDToutput(thevertex_zut, "topMass"       , theFile_zut);
   //PlotBDToutput(thevertex_zut, "totMass"       , theFile_zut);
   PlotBDToutput(thevertex_zut, "deltaPhilb"    , theFile_zut);
   //PlotBDToutput(thevertex_zut, "deltaRlb"      , theFile_zut);
   //PlotBDToutput(thevertex_zut, "deltaRTopZ"    , theFile_zut);
   PlotBDToutput(thevertex_zut, "asym"          , theFile_zut);
   PlotBDToutput(thevertex_zut, "Zpt"           , theFile_zut);
   PlotBDToutput(thevertex_zut, "ZEta"          , theFile_zut);
   //PlotBDToutput(thevertex_zut, "topPt"         , theFile_zut);
   //PlotBDToutput(thevertex_zut, "topEta"        , theFile_zut);
   PlotBDToutput(thevertex_zut, "NJets"         , theFile_zut);   
   PlotBDToutput(thevertex_zut, "NBJets"        , theFile_zut);  
   PlotBDToutput(thevertex_zut, "deltaRZl"      , theFile_zut); 
   PlotBDToutput(thevertex_zut, "deltaPhiZmet"  , theFile_zut);
   PlotBDToutput(thevertex_zut, "btagDiscri"    , theFile_zut);
   //PlotBDToutput(thevertex_zut, "leptWPt"       , theFile_zut);
   //PlotBDToutput(thevertex_zut, "leptWEta"      , theFile_zut);	    
   PlotBDToutput(thevertex_zut, "leadJetPt"     , theFile_zut);	      
   PlotBDToutput(thevertex_zut, "leadJetEta"    , theFile_zut);         
   PlotBDToutput(thevertex_zut, "deltaPhiZleptW", theFile_zut);
   
   */
   /*PlotBDToutput(thevertex_zct, "topMass"       , theFile_zct);
   PlotBDToutput(thevertex_zct, "totMass"       , theFile_zct);
   PlotBDToutput(thevertex_zct, "deltaPhilb"    , theFile_zct);
   PlotBDToutput(thevertex_zct, "deltaRlb"      , theFile_zct);
   PlotBDToutput(thevertex_zct, "deltaRTopZ"    , theFile_zct);
   PlotBDToutput(thevertex_zct, "asym"          , theFile_zct);
   PlotBDToutput(thevertex_zct, "Zpt"           , theFile_zct);
   PlotBDToutput(thevertex_zct, "ZEta"          , theFile_zct);
   PlotBDToutput(thevertex_zct, "topPt"         , theFile_zct);
   PlotBDToutput(thevertex_zct, "topEta"        , theFile_zct);
   PlotBDToutput(thevertex_zct, "NJets"         , theFile_zct);   
   PlotBDToutput(thevertex_zct, "NBJets"        , theFile_zct);  
   PlotBDToutput(thevertex_zct, "deltaRZl"      , theFile_zct); 
   PlotBDToutput(thevertex_zct, "deltaPhiZmet"  , theFile_zct);
   PlotBDToutput(thevertex_zct, "btagDiscri"    , theFile_zct);
   PlotBDToutput(thevertex_zct, "leptWPt"       , theFile_zct);
   PlotBDToutput(thevertex_zct, "leptWEta"      , theFile_zct);	    
   PlotBDToutput(thevertex_zct, "leadJetPt"     , theFile_zct);	      
   PlotBDToutput(thevertex_zct, "leadJetEta"    , theFile_zct);         
   PlotBDToutput(thevertex_zct, "deltaPhiZleptW", theFile_zct);
   */


}
