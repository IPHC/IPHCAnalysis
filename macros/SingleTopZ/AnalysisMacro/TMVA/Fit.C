#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooHistPdf.h"

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
#include <vector>
#include <TLegend.h>
#include <TApplication.h>
#include <TLatex.h>

#include <RooRealVar.h>
#include <RooDataHist.h>

#include <RooGaussian.h>
#include <RooConstVar.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooWorkspace.h>
#include <RooHistPdf.h>


#include <sstream>

using namespace std;
using namespace RooFit ;


bool ExtractHisto(std::string sel, std::string obs, std::string channel, TH1F*& zMass_Data, TH1F*& totMC, TH1F*& totMC_DY, TH1F*& DYtemplate)
{
  
  std::string filename_data    = "TMVApp.root";
  std::string filename_template  = "trainingBDT_FCNC_tZ.root";
  
  TFile *f_data      = new TFile(filename_data.c_str());
  TFile *f_template  = new TFile(filename_template.c_str());
  
  std::cout << "54" << endl;
  
  f_data->cd();
  zMass_Data = dynamic_cast<TH1F*>(gROOT->FindObject("MVA_BDT"));
  
  
  
  std::cout << "61" << endl;
  f_template->cd("Method_BDT/BDT/");
  totMC      = dynamic_cast<TH1F*>(gROOT->FindObject("MVA_BDT_S") );
  totMC_DY   = dynamic_cast<TH1F*>(gROOT->FindObject("MVA_BDT_B") );
  DYtemplate = dynamic_cast<TH1F*>(gROOT->FindObject("MVA_BDT_B") );
  std::cout << "66" << endl;
  
  
  std::cout << "nbin data    " << zMass_Data->GetNbinsX() << std::endl;
  std::cout << "min data     " << zMass_Data->GetXaxis()->GetXmin() << std::endl;
  std::cout << "max data     " << zMass_Data->GetXaxis()->GetXmax() << std::endl;
  std::cout << "nbin templte " << DYtemplate->GetNbinsX() << std::endl;
  std::cout << "max templte  " << DYtemplate->GetXaxis()->GetXmin() << std::endl;
  std::cout << "min templte  " << DYtemplate->GetXaxis()->GetXmax() << std::endl;
  std::cout << "nbin templte " << totMC->GetNbinsX() << std::endl;
  std::cout << "min templte  " << totMC->GetXaxis()->GetXmin() << std::endl;
  std::cout << "max templte  " << totMC->GetXaxis()->GetXmax() << std::endl;
  
  
  // Rebinning
  UInt_t rebin=2;
  zMass_Data->Rebin(rebin);
  std::cout << "76" << endl;
  totMC->Rebin(rebin);
  std::cout << "78" << endl;
  totMC_DY->Rebin(rebin);
  std::cout << "80" << endl;
  //DYtemplate->Rebin(rebin);
  
  std::cout << "74" << endl;
  return true;
}


void LikelihoodFit(string channel, TH1F* zMass_Data, TH1F* totMC, TH1F* totMC_DY, TH1F* DYtemplate, TApplication* theApp)
{
  std::cout << "Loading RooFit libraries ..." << std::endl;
  

  //*********************************************
  // create the RooFit environement
  RooWorkspace *w = new RooWorkspace("w",kTRUE) ;
  if (w==0) { std::cout << "ERROR : Impossible to configure RooFit !" << std::endl; return; }

  //*********************************************
  // create a RooFit variable :
  // the variable that will be used for the fit
  //*********************************************

  std::cout << "Creating RooFit variables ..." << std::endl;

  RooRealVar * LL_Zmass  ;
  LL_Zmass = new RooRealVar("LL_Zmass", "M_{ee#mu}",     zMass_Data->GetXaxis()->GetXmin() , zMass_Data->GetXaxis()->GetXmax());
  
  std::cout << "line 113 " << std::endl;
  
  //create RooDataHisto for the data (to be fitted)
  RooDataHist* histoLL_Zmass = new RooDataHist("histoFit_LL", "histoFit_LL", *LL_Zmass,  zMass_Data );
  std::cout << "line 117 " << std::endl;

  //create RooDataHisto for the MC DY : will be used to create template 1
  RooDataHist* histoFit_Template_DYEM = new RooDataHist("histoFit_Template_DYEM", 
                                                        "histoFit_Template_DYEM",
                                                        *LL_Zmass, DYtemplate);

  std::cout << "line 124 " << std::endl;
  //create RooDataHisto for the MC DY : will be used to create template 2
  RooDataHist* histoFit_Template_OtherLL = new RooDataHist("histoFit_Template_OtherLL", 
                                                           "histoFit_Template_OtherLL",
                                                           *LL_Zmass, totMC);

  std::cout << "line 130 " << std::endl;
  //convert  RooDataHisto into the template 1
  RooHistPdf* histoFit_Template_DYEM_pdf = new RooHistPdf("histoFit_Template_DYEM_pdf", 
                                                          "histoFit_Template_DYEM_pdf", 
                                                          *LL_Zmass, *histoFit_Template_DYEM);

  std::cout << "line 136 " << std::endl;
  //convert  RooDataHisto into the template 2
  RooHistPdf* histoFit_Template_OtherLL_pdf = new RooHistPdf("histoFit_Template_OtherLL_pdf", 
                                                             "histoFit_Template_OtherLL_pdf",
                                                             *LL_Zmass, *histoFit_Template_OtherLL);

  std::cout << "line 142 " << std::endl;
  //define a coefficient : it is the output of the fit
  //it represents the contribution (fraction) of one of the 2 templates with resepct to the data after the fit :
  // N(event ttbar) = coeff*N(event data)
  RooRealVar coeffModel("coeffModel", "coeffModel", 0.5, 0., 1.);

  std::cout << "line 148 " << std::endl;
  //associate the templates, the data histo and the coeff.
  RooAddPdf *pdf_sum = new RooAddPdf("pdf_sum"," test pdf sum",RooArgSet(*histoFit_Template_DYEM_pdf, *histoFit_Template_OtherLL_pdf), coeffModel);

  std::cout << "line 152 " << std::endl;
  //do the fit.
  RooFitResult* myFitResults_all = pdf_sum->fitTo(*histoLL_Zmass, Save()) ;


  //print the coeff after the fit
  std::cout << "-----------------------------------------------------------------------" << std::endl;
  std::cout << "                              Initial info" << std::endl;
  std::cout << "-----------------------------------------------------------------------" << std::endl;
  std::cout << "Nevents of Data        = " << zMass_Data->Integral() << std::endl;
  std::cout << "Nevents of DY in MC    = " << totMC_DY->Integral() << std::endl;
  std::cout << "Nevents of Other       = " << totMC->Integral() << std::endl;
  std::cout << std::endl;

  std::cout << "-----------------------------------------------------------------------" << std::endl;
  std::cout << "                              Fit results" << std::endl;
  std::cout << "-----------------------------------------------------------------------" << std::endl;
 
  //calculate the various contributions
  coeffModel.Print() ;
  std::cout << "Data = ( " << coeffModel.getVal()*zMass_Data->Integral() 
            << "  +/- "    << coeffModel.getError()*zMass_Data->Integral()
            << " ) events of background"<< std::endl;
  std::cout << "     + ( " << (1-coeffModel.getVal())*zMass_Data->Integral() 
            << " +/- " << coeffModel.getError()*zMass_Data->Integral()
            << " ) events of signal"  << std::endl;



  // Draw
  TCanvas *c1 = new TCanvas("c1", "plots",400,400,800,600);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);

  c1->SetFillColor(10);
  c1->SetFillStyle(4000);
  c1->SetBorderSize(2);
  c1->SetLogy(0);

  zMass_Data->SetMarkerStyle(20);
  zMass_Data->Draw("ep");

  
  zMass_Data->Draw("epsame");
  cout << "line 323 " << endl;

  // Drawing plot

  //create a "frame" : a kind of TCanvas used to display the result of the fit
  RooPlot* frame = LL_Zmass->frame() ;

  //plot the data on the frame
  histoLL_Zmass->plotOn(frame);

  //plots the templates (after the fit) on the frame
   pdf_sum->plotOn(frame, Components(*histoFit_Template_DYEM_pdf), VisualizeError(*myFitResults_all), FillColor(kGreen), LineWidth(2) );
   pdf_sum->plotOn(frame, Components(*histoFit_Template_OtherLL_pdf), VisualizeError(*myFitResults_all), FillColor(kRed),   LineWidth(2) );
   pdf_sum->plotOn(frame, LineStyle(kDashed), VisualizeError(*myFitResults_all), FillColor(kBlue), LineWidth(2) );

  histoLL_Zmass->plotOn(frame);
  frame->Draw("hsame");

  // Drawing legend
  TH1F * histoFitDY    = new TH1F("histoFitDY",    "histoFitDY",    0, 1, 100);
  TH1F * histoFitOther = new TH1F("histoFitOther", "histoFitOther", 0, 1, 100);
  TH1F * histoFitTot   = new TH1F("histoFitTot",   "histoFitTot",   0, 1, 100);

  histoFitDY->SetFillColor(3);
  histoFitOther->SetFillColor(2);
  histoFitTot->SetFillColor(4);
  TLegend* qw = new TLegend(0.75,0.70,0.98,0.98);
  qw->AddEntry(zMass_Data,         "Data" ,                "p");
  qw->AddEntry(histoFitTot,        "result of the Fit" ,   "f");
  qw->AddEntry(histoFitDY,         "background component" ,        "f");
  qw->AddEntry(histoFitOther,      "signal component" ,    "f");

  qw->Draw();

  //cout << "line 357 " << endl;

  //std::cout << zMass_Data->GetXaxis()->GetXmax() << std::endl;
  //std::cout << zMass_Data->GetXaxis()->GetXmin() << std::endl;
  //std::cout << zMass_Data->GetMinimum() << std::endl;
  //std::cout << zMass_Data->GetMaximum() << std::endl;
  
  double aX = zMass_Data->GetXaxis()->GetXmax()/100. - zMass_Data->GetXaxis()->GetXmin();
  double bX = zMass_Data->GetXaxis()->GetXmin();

  double aY = zMass_Data->GetMaximum()/100. - zMass_Data->GetMinimum();
  double bY = zMass_Data->GetMinimum();

  std::stringstream str;
  str << "DY scale factor = " << coeffModel.getVal()*zMass_Data->Integral()/(totMC_DY->Integral()) << " #pm " <<
    coeffModel.getError()*zMass_Data->Integral()/(totMC_DY->Integral()) << "";


  TLatex *tex1 = new TLatex(50*aX+bX,75*aY+bY,str.str().c_str());
  tex1->SetTextSize(0.03);
  tex1->SetTextColor(1);

  str.str("");
  str << "OT scale factor = " << (1.-coeffModel.getVal())*zMass_Data->Integral()/(totMC->Integral()) << " #pm " <<
    coeffModel.getError()*zMass_Data->Integral()/(totMC->Integral()) << "";

  TLatex *tex2 = new TLatex(50*aX+bX,67*aY+bY,str.str().c_str());
  tex2->SetTextSize(0.03);
  tex2->SetTextColor(1);

  tex1->Draw("same");
  tex2->Draw("same");
  cout << "line 389 " << endl;
  //theApp->Run();
  
  double count = 0;
  double count_data = 0;
  for(int ibin=0; ibin<totMC_DY->GetNbinsX()+2; ibin++){
    
    //cout << totMC->GetBinContent(ibin) << endl;
    count += totMC_DY->GetBinContent(ibin);
    count_data += DYtemplate->GetBinContent(ibin);
  }
  
  cout << "counting       " << count << endl;
  cout << "counting_data  " << count_data << endl;

  std::cout << "-----------------------------------------------------------------------" << std::endl;
  std::cout << "That's all folk!" << std::endl;
  //  delete frame;
  theApp->Run();
}


int main(int argn, char **argv)
{
  TApplication theApp("QMA_plot", &argn,argv);
  TH1F* zMass_Data=0;
  TH1F* totMC=0;
  TH1F* totMC_DY=0;
  TH1F* DYtemplate=0;
  if (!ExtractHisto("leptcut", // "leptcut"
                    "mWT", // "RecoPtZ", //"RecoTopMass", //"RecoPtZ", // RecoTopMass, //mWT
                    "eee",
                    zMass_Data,totMC, totMC_DY, DYtemplate)) abort();
  LikelihoodFit("eee",zMass_Data,totMC,totMC_DY, DYtemplate,&theApp);
  return 0;
} 
