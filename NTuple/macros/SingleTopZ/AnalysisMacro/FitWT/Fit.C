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
#include <utility>

using namespace std;
using namespace RooFit ;

int systematics 	  = -1;
int systematics_TTbar	  = 0;
int systematics_SingleTop = 0;
int systematics_WW	  = 0;
int systematics_ZZ	  = 0;
int systematics_TZq	  = 0;

bool ExtractHisto(std::string sel, std::string obs, std::string channel, TH1F*& zMass_Data, TH1F*& totMC, TH1F*& totMC_DY, TH1F*& DYtemplate)
{
  // Opening the input file
  std::cout << "Opening the input file ..." << std::endl;

  std::string filename  = "../../RootFiles/proof_woBkgSF_good.root"; // for nom sample
  std::string filename2 = "../../RootFiles/proof_DataZjets_woBkgSF_good.root"; // Z enriched
  
  TFile *f_data   = new TFile(filename.c_str());
  TFile *f_data2  = new TFile(filename2.c_str());
  if (!f_data->IsOpen() ) { std::cout << "ERROR : Cannot open the file '" << filename  << "' " << std::endl; return false;}
  if (!f_data2->IsOpen()) { std::cout << "ERROR : Cannot open the file '" << filename2 << "' " << std::endl; return false;}
  f_data->cd();

  //*********************************************
  // get the histograms and prepare the sum of MC
  //*********************************************
  std::vector<TString> histoName_List;
  histoName_List.push_back("TTbarSig");
  histoName_List.push_back("TTbarBkg");
  histoName_List.push_back("TtW");
  histoName_List.push_back("TbartW");
  histoName_List.push_back("WW");
  histoName_List.push_back("ZZ");
  histoName_List.push_back("WZ");
  histoName_List.push_back("TZq");

  std::vector<TString> histoName_List_DY;
  if(channel == "mumumu"                    ) histoName_List_DY.push_back("DataMu");
  if(channel == "eemu" || channel == "mumue") histoName_List_DY.push_back("DataMuEG");
  if(channel == "eee"                       ) histoName_List_DY.push_back("DataEG");
  
  std::vector<TString> histoName_List_DY_MC;
  histoName_List_DY_MC.push_back("Zjets");
  histoName_List_DY_MC.push_back("DYToLL_M10-50");

  TString selection;
  if(sel == "leptcut" ) selection = "afterleptsel";
  else if(sel == "btag"    ) selection = "afterbjetsel";
  else {std::cout << "ERROR : sel variable '" << sel << "' is not defined" << std::endl; return false;}


  std::vector<TString> Channels;
  std::vector<TString> Datasets;
  TString theChannel;
  if(channel == "eee")         {theChannel = "eee";    Datasets.resize(1,"DataEG");   Channels.resize(1,"eee");}
  else if(channel == "mumumu") {theChannel = "mumumu"; Datasets.resize(1,"DataMu");   Channels.resize(1,"mumumu");}
  else if(channel == "eemu"  ) {theChannel = "eemu";   Datasets.resize(1,"DataMuEG"); Channels.resize(1,"eemu");}
  else if(channel == "mumue" ) {theChannel = "mumue";  Datasets.resize(1,"DataMuEG"); Channels.resize(1,"mumue");}
  else if(channel == "all"   ) 
  {
    theChannel = "all"; 
    Channels.resize(4);   Datasets.resize(4);
    Channels[0]="eee";    Datasets[0]="DataEG";
    Channels[1]="mumumu"; Datasets[1]="DataMu";
    Channels[2]="eemu";   Datasets[2]="DataMuEG";
    Channels[3]="mumue";  Datasets[3]="DataMuEG";
  }
  else {std::cout << "ERROR : channel variable '" << theChannel << "' is not defined" << std::endl; return false; }


  //  RecoTopMass_mumumu_afterbjetsel_DataMu;1 
  std::vector<TH1F*> histo_List_Data;
  for (unsigned int i=0;i<Channels.size();i++)
  {
    std::string histoname = (obs+"_"+Channels[i]+"_"+selection+"_"+Datasets[i]).Data();
    std::cout << "Loading plot for data '" << histoname << "' ... " << std::endl;
    histo_List_Data.push_back(dynamic_cast<TH1F*>( gROOT->FindObject(histoname.c_str()) ));
    if (histo_List_Data[i]==0) { std::cout << "ERROR : plot not found" << std::endl; return false; }
  }

  //*********************************************
  // get and summ histograms for MC
  //*********************************************
  //*********************
  // Fill Histo Vector
  //*********************
  std::vector<TH1F*> histo_List;
  for (unsigned int j=0; j<Channels.size();j++)
  for(unsigned int i=0; i<histoName_List.size(); i++)
  {
    std::string histoname = (obs+"_"+Channels[j]+"_"+selection+"_"+histoName_List[i]).Data();
    std::cout << "Loading plot for MC '" << histoname << "' ... " << std::endl;
    histo_List.push_back(dynamic_cast<TH1F*>( gROOT->FindObject(histoname.c_str()) ) );
    if (histo_List.back()==0) {std::cout << "ERROR : plot not found" << std::endl; return false; }
    std::cout<<"===== "<<histo_List.back()->Integral()<<std::endl;
    
    double weightfactor = 1;
    if(systematics ==  1  ) weightfactor = 1.3;
    if(systematics == -1 ) weightfactor = 0.7;
     
    if(histoName_List[i] == "TTbarSig" &&  systematics_TTbar!=0     ) histo_List.back()->Scale( weightfactor);
    if(histoName_List[i] == "TTbarBkg" &&  systematics_TTbar!=0     ) histo_List.back()->Scale( weightfactor);
    if(histoName_List[i] == "TtW"      &&  systematics_SingleTop!=0 ) histo_List.back()->Scale( weightfactor);
    if(histoName_List[i] == "TbartW"   &&  systematics_SingleTop!=0 ) histo_List.back()->Scale( weightfactor);
    if(histoName_List[i] == "WW"       &&  systematics_WW!=0        ) histo_List.back()->Scale( weightfactor);
    if(histoName_List[i] == "ZZ"       &&  systematics_ZZ!=0        ) histo_List.back()->Scale( weightfactor);
    if(histoName_List[i] == "TZq"      &&  systematics_TZq!=0       ) histo_List.back()->Scale( weightfactor);
    
    
  }

  //*********************
  // Fill Histo Vector
  //*********************
  std::vector<TH1F*> histo_List_DY;
  for (unsigned int j=0; j<Channels.size();j++)
  for(unsigned int i=0; i<histoName_List_DY.size(); i++)
  { 
    f_data2->cd();
    std::string histoname = (obs+"_"+Channels[j]+"_"+selection+"_"+histoName_List_DY[i]).Data();
    std::cout << "Loading plot DY from data '" << histoname << "' ... " << std::endl;
    histo_List_DY.push_back(dynamic_cast<TH1F*>( gROOT->FindObject(histoname.c_str()) ));
    if (histo_List.back()==0) {std::cout << "ERROR : plot not found'" << std::endl; return false; }
  }

  //*********************
  // Fill Histo Vector
  //*********************
  std::vector<TH1F*> histo_List_DY_MC;
  for (unsigned int j=0; j<Channels.size();j++)
  for(unsigned int i=0; i<histoName_List_DY_MC.size(); i++)
  { 
    f_data->cd();
    std::string histoname = (obs+"_"+Channels[j]+"_"+selection+"_"+histoName_List_DY_MC[i]).Data();
    std::cout << "Loading plot for DY from MC '" << histoname << "' ... " << std::endl;
    histo_List_DY_MC.push_back(dynamic_cast<TH1F*>( gROOT->FindObject(histoname.c_str()) ));
    if (histo_List.back()==0) {std::cout << "ERROR : plot not found'" << std::endl; return false; }
  }

  f_data->cd();
  
  
  std::cout << "Summing Data samples ..." << std::endl;
  zMass_Data = new TH1F( "totData","totData",
               histo_List_Data.front()->GetNbinsX(), 
               histo_List_Data.front()->GetXaxis()->GetXmin(), 
               histo_List_Data.front()->GetXaxis()->GetXmax());
  for(unsigned int i=0; i<histo_List_Data.size(); i++)
  {
    histo_List_Data[i]->Sumw2();
    zMass_Data->Add( zMass_Data, histo_List_Data[i], 1, 1);
  }
 
  //*********************
  // summ up MC templates
  //*********************
  std::cout << "Summing MC samples ... " << std::endl;
  totMC = new TH1F("totMC", "totMC",histo_List_Data.front()->GetNbinsX(), 
                                    histo_List_Data.front()->GetXaxis()->GetXmin(), 
                                    histo_List_Data.front()->GetXaxis()->GetXmax());
  for(unsigned int i=0; i<histo_List.size(); i++)
  {
    histo_List[i]->Sumw2();
    totMC->Add( totMC, histo_List[i], 1, 1);
  }


  //***********************
  // summ up data templates
  //***********************
  DYtemplate = new TH1F("DYtemplate", "DYtemplate", 
                                    histo_List_Data.front()->GetNbinsX(), 
                                    histo_List_Data.front()->GetXaxis()->GetXmin(), 
                                    histo_List_Data.front()->GetXaxis()->GetXmax());
  for(unsigned int i=0; i<histo_List_DY.size(); i++)
  {
    histo_List_DY[i]->Sumw2();
    DYtemplate->Add( DYtemplate, histo_List_DY[i], 1, 1);
  }
  std::cout << DYtemplate->GetEntries() << std::endl;
  
  
  //*********************
  // summ up MC templates
  //*********************
  totMC_DY = new TH1F("totMC_DY", "totMC_DY", 
                                    histo_List_Data.front()->GetNbinsX(), 
                                    histo_List_Data.front()->GetXaxis()->GetXmin(), 
                                    histo_List_Data.front()->GetXaxis()->GetXmax());
  for(unsigned int i=0; i<histo_List_DY_MC.size(); i++)
  {
    //histo_List_DY_MC[i]->Sumw2();
    totMC_DY->Add( totMC_DY, histo_List_DY_MC[i], 1, 1);
  }
  std::cout << totMC_DY->GetEntries() << std::endl;


  // Rebinning
  UInt_t rebin=4;
  totMC->Rebin(rebin);
  totMC_DY->Rebin(rebin);
  DYtemplate->Rebin(rebin);
  zMass_Data->Rebin(rebin);


  return true;
}


std::pair<double, double >  LikelihoodFit(string channel, TH1F* zMass_Data, TH1F* totMC, TH1F* totMC_DY, TH1F* DYtemplate, 
TApplication* theApp)
{
  std::cout << "Loading RooFit libraries ..." << std::endl;
  
  std::pair <double, double > thereturn;
  thereturn.first = -1;
  thereturn.second = -1;
  
  //*********************************************
  // create the RooFit environement
  RooWorkspace *w = new RooWorkspace("w",kTRUE) ;
  if (w==0) { std::cout << "ERROR : Impossible to configure RooFit !" << std::endl; return thereturn; }

  //*********************************************
  // create a RooFit variable :
  // the variable that will be used for the fit
  //*********************************************

  std::cout << "Creating RooFit variables ..." << std::endl;

  RooRealVar * LL_Zmass  ;
  if      (channel == "eemu")   LL_Zmass = new RooRealVar("LL_Zmass", "M_{T}^{W} (ee#mu)",     zMass_Data->GetXaxis()->GetXmin() , zMass_Data->GetXaxis()->GetXmax());
  else if (channel == "mumue")  LL_Zmass = new RooRealVar("LL_Zmass", "M_{T}^{W} (#mu#mue)",   zMass_Data->GetXaxis()->GetXmin() , zMass_Data->GetXaxis()->GetXmax());
  else if (channel == "mumumu") LL_Zmass = new RooRealVar("LL_Zmass", "M_{T}^{W} (#mu#mu#mu)", zMass_Data->GetXaxis()->GetXmin() , zMass_Data->GetXaxis()->GetXmax());
  else if (channel == "eee")    LL_Zmass = new RooRealVar("LL_Zmass", "M_{T}^{W} (eee)",       zMass_Data->GetXaxis()->GetXmin() , zMass_Data->GetXaxis()->GetXmax());
  else if (channel == "all")    LL_Zmass = new RooRealVar("LL_Zmass", "M_{T}^{W} (lll)",       zMass_Data->GetXaxis()->GetXmin() , zMass_Data->GetXaxis()->GetXmax());

  //create RooDataHisto for the data (to be fitted)
  RooDataHist* histoLL_Zmass = new RooDataHist("histoFit_LL", "histoFit_LL", *LL_Zmass,  zMass_Data );

  //create RooDataHisto for the MC DY : will be used to create template 1
  RooDataHist* histoFit_Template_DYEM = new RooDataHist("histoFit_Template_DYEM", 
                                                        "histoFit_Template_DYEM",
                                                        *LL_Zmass, DYtemplate);

  //create RooDataHisto for the MC DY : will be used to create template 2
  RooDataHist* histoFit_Template_OtherLL = new RooDataHist("histoFit_Template_OtherLL", 
                                                           "histoFit_Template_OtherLL",
                                                           *LL_Zmass, totMC);

  //convert  RooDataHisto into the template 1
  RooHistPdf* histoFit_Template_DYEM_pdf = new RooHistPdf("histoFit_Template_DYEM_pdf", 
                                                          "histoFit_Template_DYEM_pdf", 
                                                          *LL_Zmass, *histoFit_Template_DYEM);

  //convert  RooDataHisto into the template 2
  RooHistPdf* histoFit_Template_OtherLL_pdf = new RooHistPdf("histoFit_Template_OtherLL_pdf", 
                                                             "histoFit_Template_OtherLL_pdf",
                                                             *LL_Zmass, *histoFit_Template_OtherLL);

  //define a coefficient : it is the output of the fit
  //it represents the contribution (fraction) of one of the 2 templates with resepct to the data after the fit :
  // N(event ttbar) = coeff*N(event data)
  RooRealVar coeffModel("coeffModel", "coeffModel", 0.5, 0., 1.);

  //associate the templates, the data histo and the coeff.
  RooAddPdf *pdf_sum = new RooAddPdf("pdf_sum"," test pdf sum",RooArgSet(*histoFit_Template_DYEM_pdf, *histoFit_Template_OtherLL_pdf), coeffModel);

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
            << " ) events of DY"<< std::endl;
  std::cout << "     + ( " << (1-coeffModel.getVal())*zMass_Data->Integral() 
            << " +/- " << coeffModel.getError()*zMass_Data->Integral()
            << " ) events of other MC"  << std::endl;



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
  zMass_Data->GetXaxis()->SetTitle("m_{T}(W)");
  
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
  qw->AddEntry(histoFitDY,         "Fake lepton" ,        "f");
  qw->AddEntry(histoFitOther,      "non-Fake lepton" ,    "f");

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
  str << "SF fake = " << coeffModel.getVal()*zMass_Data->Integral()/(totMC_DY->Integral()) << " #pm " <<
    coeffModel.getError()*zMass_Data->Integral()/(totMC_DY->Integral()) << "";


  TLatex *tex1 = new TLatex(50*aX+bX,75*aY+bY,str.str().c_str());
  tex1->SetTextSize(0.03);
  tex1->SetTextColor(1);

  str.str("");
  str << "SF non fake = " << (1.-coeffModel.getVal())*zMass_Data->Integral()/(totMC->Integral()) << " #pm " <<
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
  
  if(systematics==0) c1->Print("mWT.pdf");
 
  thereturn.first  = coeffModel.getVal()*zMass_Data->Integral()/(totMC_DY->Integral());
  thereturn.second = (1.-coeffModel.getVal())*zMass_Data->Integral()/(totMC->Integral());
  
  return thereturn;
  
  //  delete frame;
}


int main(int argn, char **argv)
{
  TApplication theApp("QMA_plot", &argn,argv);
  TH1F* zMass_Data=0;
  TH1F* totMC     =0;
  TH1F* totMC_DY  =0;
  TH1F* DYtemplate=0;
  
  std::string channel = "mumumu";

  if(argn>1) channel = argv[1];
  
  
  systematics	        = 0;
  systematics_TTbar     = 0;
  systematics_SingleTop = 0;
  systematics_WW        = 0;
  systematics_ZZ        = 0;
  systematics_TZq       = 0;
  
 
  
  
  //------------- For ttbar syst values -----------------
  systematics	        = -1;
  systematics_TTbar     = 1;
  systematics_SingleTop = 0;
  systematics_WW        = 0;
  systematics_ZZ        = 0;
  systematics_TZq       = 0;
  
  
  if (!ExtractHisto("leptcut", // "leptcut"
                    "mWT", // "RecoPtZ", //"RecoTopMass", //"RecoPtZ", // RecoTopMass, //mWT
                    channel,
                    zMass_Data,totMC, totMC_DY, DYtemplate)) abort();
		    
  
  std::pair <double, double > backTTbar1= LikelihoodFit(channel,zMass_Data,totMC,totMC_DY, DYtemplate,&theApp);
  std::cout << " backTTbarfake SF " << backTTbar1.first << "  non fake sf " << backTTbar1.second << std::endl;
  
  
  systematics	        = 1;
  
  if (!ExtractHisto("leptcut", // "leptcut"
                    "mWT", // "RecoPtZ", //"RecoTopMass", //"RecoPtZ", // RecoTopMass, //mWT
                    channel,
                    zMass_Data,totMC, totMC_DY, DYtemplate)) abort();
		    
  
  std::pair <double, double > backTTbar2= LikelihoodFit(channel,zMass_Data,totMC,totMC_DY, DYtemplate,&theApp);
  std::cout << " backTTbarfake SF " << backTTbar2.first << "  non fake sf " << backTTbar2.second << std::endl;
  
  //------------- For singletop syst values -----------------
  systematics	        = 1;
  systematics_TTbar     = 0;
  systematics_SingleTop = 1;
  systematics_WW        = 0;
  systematics_ZZ        = 0;
  systematics_TZq       = 0;
  if (!ExtractHisto("leptcut", // "leptcut"
                    "mWT", // "RecoPtZ", //"RecoTopMass", //"RecoPtZ", // RecoTopMass, //mWT
                    channel,
                    zMass_Data,totMC, totMC_DY, DYtemplate)) abort();
		    
  
  std::pair <double, double > backStop1 = LikelihoodFit(channel,zMass_Data,totMC,totMC_DY, DYtemplate,&theApp);
  std::cout << " backTTbarfake SF " << backStop1.first << "  non fake sf " << backStop1.second << std::endl;
  
  systematics	        = -1;
  if (!ExtractHisto("leptcut", // "leptcut"
                    "mWT", // "RecoPtZ", //"RecoTopMass", //"RecoPtZ", // RecoTopMass, //mWT
                    channel,
                    zMass_Data,totMC, totMC_DY, DYtemplate)) abort();
		    
  
  std::pair <double, double > backStop2 = LikelihoodFit(channel,zMass_Data,totMC,totMC_DY, DYtemplate,&theApp);
  std::cout << " backTTbarfake SF " << backStop2.first << "  non fake sf " << backStop2.second << std::endl;
    
  //------------- For WW syst syst values -----------------
  systematics	        = 1;
  systematics_TTbar     = 0;
  systematics_SingleTop = 0;
  systematics_WW        = 1;
  systematics_ZZ        = 0;
  systematics_TZq       = 0;
  if (!ExtractHisto("leptcut", // "leptcut"
                    "mWT", // "RecoPtZ", //"RecoTopMass", //"RecoPtZ", // RecoTopMass, //mWT
                    channel,
                    zMass_Data,totMC, totMC_DY, DYtemplate)) abort();
		    
  
  std::pair <double, double > backWW1 = LikelihoodFit(channel,zMass_Data,totMC,totMC_DY, DYtemplate,&theApp);
  std::cout << " backTTbarfake SF " << backWW1.first << "  non fake sf " << backWW1.second << std::endl;
  
    
  systematics	        = -1;
  if (!ExtractHisto("leptcut", // "leptcut"
                    "mWT", // "RecoPtZ", //"RecoTopMass", //"RecoPtZ", // RecoTopMass, //mWT
                    channel,
                    zMass_Data,totMC, totMC_DY, DYtemplate)) abort();
		    
  
  std::pair <double, double > backWW2 = LikelihoodFit(channel,zMass_Data,totMC,totMC_DY, DYtemplate,&theApp);
  std::cout << " backTTbarfake SF " << backWW2.first << "  non fake sf " << backWW2.second << std::endl;
    
  //------------- For ZZ syst syst values -----------------
  systematics	        = 1;
  systematics_TTbar     = 0;
  systematics_SingleTop = 0;
  systematics_WW        = 0;
  systematics_ZZ        = 1;
  systematics_TZq       = 0;
  if (!ExtractHisto("leptcut", // "leptcut"
                    "mWT", // "RecoPtZ", //"RecoTopMass", //"RecoPtZ", // RecoTopMass, //mWT
                    channel,
                    zMass_Data,totMC, totMC_DY, DYtemplate)) abort();
		    
  
  std::pair <double, double > backZZ1 = LikelihoodFit(channel,zMass_Data,totMC,totMC_DY, DYtemplate,&theApp);
  std::cout << " backTTbarfake SF " << backZZ1.first << "  non fake sf " << backZZ1.second << std::endl;
  
  
  systematics	        = -1;
  if (!ExtractHisto("leptcut", // "leptcut"
                    "mWT", // "RecoPtZ", //"RecoTopMass", //"RecoPtZ", // RecoTopMass, //mWT
                    channel,
                    zMass_Data,totMC, totMC_DY, DYtemplate)) abort();
		    
  
  std::pair <double, double > backZZ2 = LikelihoodFit(channel,zMass_Data,totMC,totMC_DY, DYtemplate,&theApp);
  std::cout << " backTTbarfake SF " << backZZ2.first << "  non fake sf " << backZZ2.second << std::endl;
  
  //------------- For TZq syst syst values -----------------
  systematics	        = 1;
  systematics_TTbar     = 0;
  systematics_SingleTop = 0;
  systematics_WW        = 0;
  systematics_ZZ        = 0;
  systematics_TZq       = 1;
  if (!ExtractHisto("leptcut", // "leptcut"
                    "mWT", // "RecoPtZ", //"RecoTopMass", //"RecoPtZ", // RecoTopMass, //mWT
                    channel,
                    zMass_Data,totMC, totMC_DY, DYtemplate)) abort();
		    
  
  std::pair <double, double > backTZq1 = LikelihoodFit(channel,zMass_Data,totMC,totMC_DY, DYtemplate,&theApp);
  std::cout << " backTZq SF " << backTZq1.first << "  non fake sf " << backTZq1.second << std::endl;
  
  
  systematics	        = -1;
  if (!ExtractHisto("leptcut", // "leptcut"
                    "mWT", // "RecoPtZ", //"RecoTopMass", //"RecoPtZ", // RecoTopMass, //mWT
                    channel,
                    zMass_Data,totMC, totMC_DY, DYtemplate)) abort();
		    
  
  std::pair <double, double > backTZq2 = LikelihoodFit(channel,zMass_Data,totMC,totMC_DY, DYtemplate,&theApp);
  std::cout << " backTZq SF " << backTZq2.first << "  non fake sf " << backTZq2.second << std::endl;
  
  
 
  //------------- For nominal values -----------------
  
  systematics	        = 0;
  systematics_TTbar     = 0;
  systematics_SingleTop = 0;
  systematics_WW        = 0;
  systematics_ZZ        = 0;
  systematics_TZq       = 0;
  
  
  if (!ExtractHisto("leptcut", // "leptcut"
                    "mWT", // "RecoPtZ", //"RecoTopMass", //"RecoPtZ", // RecoTopMass, //mWT
                    channel,
                    zMass_Data,totMC, totMC_DY, DYtemplate)) abort();
		    
  
  std::pair <double, double > nom = LikelihoodFit(channel,zMass_Data,totMC,totMC_DY, DYtemplate,&theApp);
  std::cout << "nom fake SF " << nom.first << "  non fake sf " << nom.second << std::endl;
  
  
  double ttbaruncert1 = (backTTbar1.first-nom.first)/nom.first;
  double ttbaruncert2 = (backTTbar2.first-nom.first)/nom.first;
  
  double Stopuncert1 = (backStop1.first-nom.first)/nom.first;
  double Stopuncert2 = (backStop2.first-nom.first)/nom.first;
  
  double WWuncert1 = (backWW1.first-nom.first)/nom.first;
  double WWuncert2 = (backWW2.first-nom.first)/nom.first;
  
  double ZZuncert1 = (backZZ1.first-nom.first)/nom.first;
  double ZZuncert2 = (backZZ2.first-nom.first)/nom.first;
  
  double TZquncert1 = (backTZq1.first-nom.first)/nom.first;
  double TZquncert2 = (backTZq2.first-nom.first)/nom.first;
  
  
  if(channel == "mumumu") std::cout << "  $\\mu\\mu\\mu$ &  " ;
  if(channel == "mumue")  std::cout << "  $\\mu\\mu e$   &  " ;   
  if(channel == "eemu")   std::cout << "  $ee\\mu$       &  " ;   
  if(channel == "eee")    std::cout << "  $eee$          &  " ;
   
   
  std::cout <<  ttbaruncert1 << "/" <<  ttbaruncert2  << " & " ;
  std::cout <<  Stopuncert1  << "/" <<  Stopuncert2   << " & " ;
  std::cout <<  WWuncert1    << "/" <<  WWuncert2     << " & " ;
  std::cout <<  ZZuncert1    << "/" <<  ZZuncert2     << " & " ;
  std::cout <<  TZquncert1   << "/" <<  TZquncert2    << " \\\\ " << endl;
  
  
  double totsys =  (backTTbar1.first-nom.first)*(backTTbar1.first-nom.first);
  if(fabs(backTTbar1.first-nom.first) < fabs(backTTbar2.first-nom.first) )  totsys =  (backTTbar2.first-nom.first)*(backTTbar2.first-nom.first);
  
  std::cout << "totsys1 " << totsys << std::endl;
  
  if(fabs(backStop1.first-nom.first)  > fabs(backStop2.first-nom.first)) totsys += (backStop1.first-nom.first)*(backStop1.first-nom.first);
  else totsys += (backStop2.first-nom.first)*(backStop2.first-nom.first);
  
  
  std::cout << "totsys2 " << totsys << std::endl;
  
  if(fabs(backWW1.first-nom.first)  > fabs(backWW2.first-nom.first)) totsys += (backWW1.first-nom.first)*(backWW1.first-nom.first);
  else totsys += (backWW2.first-nom.first)*(backWW2.first-nom.first);
  
  std::cout << "totsys3 " << totsys << std::endl;
  
  if(fabs(backZZ1.first-nom.first)  > fabs(backZZ2.first-nom.first)) totsys += (backZZ1.first-nom.first)*(backZZ1.first-nom.first);
  else totsys += (backZZ2.first-nom.first)*(backZZ2.first-nom.first);
  
  std::cout << "totsys4 " << totsys << std::endl;
  
  if(fabs(backTZq1.first-nom.first)  > fabs(backTZq2.first-nom.first)) totsys += (backTZq1.first-nom.first)*(backTZq1.first-nom.first);
  else totsys += (backTZq2.first-nom.first)*(backTZq2.first-nom.first);
  
  std::cout << "totsys5 " << totsys << std::endl;
  
  
  
  std::cout << "result " <<  nom.first << " +\\- " << pow(totsys, 0.5) << endl;
  
  theApp.Run();
  return 0;
} 




