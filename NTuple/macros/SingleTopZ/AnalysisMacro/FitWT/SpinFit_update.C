#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <Cintex/Cintex.h>
#include <utility>
#include "TH1F.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TPaveText.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include <map>
#include <ctime>
#include <TKey.h>
#include <TRegexp.h>
#include <TLatex.h>
#include "Math/GSLRndmEngines.h"
#include "Math/Random.h"
#include "TTree.h"
#include "RooRandom.h"
#include <TLegend.h>

/*
   .L SpiFit.C++
    j = new SpiFit()
    j->fitPoint() for a single mass
    j->fitAll() for all masses
    j->findMin() to find the minimum after fitAll()
 */

using namespace std; 
using namespace RooFit;
// using namespace RooStats;

//define GAUSS

//define SYST


static const bool systematics = true;
static const bool systematicsPDF = false;
static const int  maxPdfTemplates = 41;

static const bool useFunction = false;
static const bool do_bdf_systematics = false;

static const bool useRatio = false;

static const int maxType = 3;
static const char type [maxType][10] =//{"EE_1W", "EE_2W"};
 {"eee","eemu","mumue", "mumumu"};
// {"mumu"};
//   {"EE_0B", "EM_0B", "MM_0B"};
//   {"EE_1B", "EM_1B", "MM_1B"};
//   {"EE_2B", "EM_2B", "MM_2B"};
//{"0B","1B","2B"};
// {"1B","2B"};
// {"2B"};
//   {"EE_0B", "EE_1B", "EE_2B", 
//    "EM_0B", "EM_1B", "EM_2B", 
//    "MM_0B", "MM_1B", "MM_2B"};
//{"EM"};


static const bool fixedSample = false;
static int nTotSample[maxType];// = {11,28,41}; // PF b
// static const int nTotSample[3] = {10,31,40};  //JPT b
//static const int nTotSample[3] = {15, 43, 22}; // PF
// static const int nTotSample[3] = {13, 44, 20};  //JPT
static const bool debug = true;

static char absolutePath[100] = "FinalFiles/";
static char cut[20] = "afterbtag1";

class SpinFit {
public:

  SpinFit();
  SpinFit(string noCorr, string Corr);
  ~SpinFit(){};
  void setup();
  void quietFit(bool quiet) {quietFit_ = (quiet?-1:0);}
  double fit();
  pair<double,double> findMin();
  int fixBackground();
  void fixBackground(int fix);

  void calibration(int numberPE, float lowerLimit, float upperLimit,
	int numberBins);
  void pdfSyst(int numberPE, float ratioCorrGen);
  void do_toys(int numberPE, float ratioCorrGen, int pdfToUse=0);
  int generate_toy(float ratioCorrGen, int pdfToUse=0);
  void generateComponent(RooDataHist *rooDataHist, char * hname, int nbr);

  void printFit(bool withFit = true, bool toFile = false, char* name = "");
  void printLL(bool toFile = false, char* name = "");
  void printAllValues();
  void nicePlot();
  float rTOs(float r) {return ((r+2)/5);}
  float sTOr(float s) {return (5*s-2);}

  TCanvas *c_min;
  void cmsprelim();
  RooWorkspace *w;

public:

  void assembleDatasetsHisto();
  void assembleDatasets(map<string,RooDataHist*> mapToImport);

  vector<TString> mcLabels, dataLabels;
  map<string, TH1F *> mcSignalUncorrScaled, mcSignalCorrScaled, mcBackgroundScaled, dataHistos;
  map<string, TH1F *> mcSignalUncorrScaled_gen, mcSignalCorrScaled_gen;
  RooCategory *sample;
  RooDataHist *data;
  RooRealVar *deltaPhi, *ratioCorr;
  float ilumi;
  void typeTag(char *name);
  RooFitResult* myFitResults_all;
  RooAbsPdf* pdffit;
  map<string, RooRealVar *> nbrBackgroundEvents;
  RooAbsReal* nll;
  int quietFit_;
  int fixBckg;
//   ROOT::Math::GSLRngMT * random;
  TRandom3 _random;
  map<string, int> generatedSignal, generatedBkg;
  TH1F *toyDataHisto;
  TH1F *toy_mean, *toy_bias;
  TH1F *toy_error;
  TH1F *toy_pull;
  TH1F *calibrationH;
  TGraph *gr;
  TGraphErrors * grc;
  int nFitFailed;
  int nFitTried;
  map<string, RooDataSet *> datasets;
  string fileWithCorrelation, filenoCorrelation, fileWithCorrelation_gen,
  	fileNoCorrelation_gen, fileData, fileBackground;
TH1F* h2; 
TH1F* h1;

};
 

SpinFit::SpinFit()
{
  
  fileWithCorrelation = "ttbar_mcatnlo_SC_ref.root";
  filenoCorrelation = "ttbar_mcatnlo_NoSC_ref.root";
  
  //fileWithCorrelation_gen = "ttbar_mcatnlo_SC_ref.root";
  //fileNoCorrelation_gen = "ttbar_mcatnlo_NoSC_ref.root";
  
  //fileWithCorrelation = "ttbar_powheg_SC.root";
  //filenoCorrelation = "ttbar_powheg_NoSC.root";
  
  //fileWithCorrelation_gen = "ttbar_powheg_SC.root";
  //fileNoCorrelation_gen = "ttbar_powheg_NoSC.root";
  
  
  //fileWithCorrelation_gen = "ttbar_mcatnlo_mtop125_SC.root";
  //fileNoCorrelation_gen = "ttbar_mcatnlo_mtop125_NoSC.root";
  
  //fileWithCorrelation_gen = "ttbar_mcatnlo_SC_btagup10.root";
  //fileNoCorrelation_gen = "ttbar_mcatnlo_NoSC_btagup10.root";
  
  //fileWithCorrelation_gen = "ttbar_mcatnlo_SC_m175.root";
  //fileNoCorrelation_gen = "ttbar_mcatnlo_SC_m175.root";
  
  //fileWithCorrelation_gen = "ttbar_mcatnlo_SC_m170.root";
  //fileNoCorrelation_gen = "ttbar_mcatnlo_SC_m170.root";
  
  //fileWithCorrelation_gen = "ttbar_mcatnlo_0.5-ren-0.5_withSC.root";
  //fileNoCorrelation_gen = "ttbar_mcatnlo_fac-0.5-ren-0.5_withNoSC.root";
  
  //fileWithCorrelation_gen = "ttbar_mcatnlo_fac-2-ren-2_withSC.root";
  //fileNoCorrelation_gen = "ttbar_mcatnlo_fac-2-ren-2_withNoSC.root";
  
  
  fileWithCorrelation_gen = "ttbar_mcatnlo_SC_notau.root";
  fileNoCorrelation_gen = "ttbar_mcatnlo_NoSC_notau.root";
  
  
  fileData = "background_mcatnlo_NoSC.root";
  fileBackground = "background_mcatnlo_NoSC.root";




  setup();
}

SpinFit::SpinFit(string noCorr, string Corr)
{
  fileData = "Backgrounds_TTbar_MCATNLO_SC.root";
  fileBackground = "Backgrounds_TTbar_MCATNLO_SC.root";

  fileWithCorrelation = "TT_mcatnlo_SC_ref.root";
  filenoCorrelation = "TT_mcatnlo_NoSC_ref.root";
  
//TT_mcatnlo_NoSC_JesUp.root
//TT_mcatnlo_NoSC_JesDown.root
//TT_mcatnlo_SC_JesUp.root
//TT_mcatnlo_SC_JesDown.root

  fileWithCorrelation_gen = Corr;
  fileNoCorrelation_gen = noCorr;
  setup();
}

void SpinFit::setup()
{
  RooRandom::randomGenerator()->SetSeed(123456789) ;
  ilumi = 4.98;
  quietFit_ = 0; //Default print level of the fit change with quiet()
  nFitFailed = 0;
  nFitTried = 0;
  data = 0;

  mcLabels.push_back("TTbarBkg");
  mcLabels.push_back("DYToLL_M10-50");
  mcLabels.push_back("Zjets");
  mcLabels.push_back("TtW");
  mcLabels.push_back("TbartW");
  mcLabels.push_back("WZ");
  mcLabels.push_back("ZZ");
  mcLabels.push_back("WW");

  dataLabels.push_back("DataMu");
  dataLabels.push_back("DataEG");
  dataLabels.push_back("DataMuEG");

  char histoName[20] = "DeltaPhiLLept";

  c_min = new TCanvas("c_min","", 600, 600);

  ratioCorr = new RooRealVar("ratioCorr", "ratioCorr", 0.557, 0., 1.);

  sample = new RooCategory("sample","sample") ;
  for (int itype = 0;itype!=maxType;++itype) sample->defineType(type[itype]);

  TFile* theFile;
  TH1::AddDirectory(kFALSE);
  char name[300], hname[150], tag[50];

  sprintf(name,"%s%s",absolutePath,fileData.c_str());
  if (debug) cout << "Load data from file "<< name<<endl;
  theFile = new TFile (name) ;

  // Read data

  float binWidth;
  for (unsigned int i = 0; i!=dataLabels.size();++i) {
    float scale = 1.;

    for (int itype = 0;itype!=maxType;++itype) {
      sprintf(hname,"%s_%s_%s_%s", histoName, type[itype],cut, dataLabels[i].Data());
      if (debug) cout << hname<<endl;
      TH1F* histo  = (TH1F*) gDirectory->Get(hname) ;
      if (histo==0) {
        cout << hname<<" not there -- \n";continue;
      }
//       histo->Rebin(5);
      histo->SetLineColor(2);
      if (dataHistos.find(type[itype]) == dataHistos.end()) {
	dataHistos[type[itype]] = (TH1F*) histo->Clone("dataHistos");
      } else {
	dataHistos[type[itype]]->Add(histo);
      }
      binWidth = dataHistos[type[itype]]->GetBinWidth(0);
    }
  }
//   float binWidth = dataHistos["ee"]->GetBinWidth(0);
//   cout << "binWidth "<< binWidth<<endl;
  deltaPhi = new RooRealVar("deltaPhi", "DeltaPhiLLept", 0., binWidth*10);
  deltaPhi->setBins(10);
  cout << "Got "<<dataHistos.size() << " data histos\n";

  assembleDatasetsHisto();
  theFile->Close();

  // Get the MC templates
  sprintf(name,"%s%s",absolutePath,fileWithCorrelation.c_str());
  if (debug) cout << "Load MC histos with correlation from file "<< name<<endl;
  theFile = new TFile (name) ;
  cout << "Will now retrieve the MC histos with correlation\n";

  for (int itype = 0;itype!=maxType;++itype) {
    sprintf(hname,"%s_%s_%s_%s", histoName, type[itype],cut, "TTbarSig");
    if (debug) cout << hname<<endl;
    TH1F* histo  = (TH1F*) gDirectory->Get(hname) ;
    if (histo==0) {
      cout << hname<<" not there -- \n";return;
    }
//       histo->Rebin(5);
    histo->SetLineColor(2);
    if (mcSignalCorrScaled.find(type[itype]) == mcSignalCorrScaled.end()) {
      mcSignalCorrScaled[type[itype]] = (TH1F*) histo->Clone("mcSignalCorrScaled");
    } else {
      mcSignalCorrScaled[type[itype]]->Add(histo);
    }
  }
  cout << "Got "<<mcSignalCorrScaled.size() << " MC histos  with correlation\n";
//   binWidth = mcSignalCorrScaled["ee"]->GetBinWidth(0);
//   cout << "binWidth "<< binWidth<<endl;

  theFile->Close();

  // Get the Background templates
  sprintf(name,"%s%s",absolutePath,fileBackground.c_str());
  if (debug) cout << "Load Background histos from file "<< name<<endl;
  theFile = new TFile (name) ;
  for (unsigned int i = 0; i!=mcLabels.size();++i) {
    float scale = 1.;

    for (int itype = 0;itype!=maxType;++itype) {
      sprintf(hname,"%s_%s_%s_%s", histoName, type[itype],cut, mcLabels[i].Data());
      if (debug) cout << hname<<endl;
      TH1F* histo  = (TH1F*) gDirectory->Get(hname) ;
      if (histo==0) {
        cout << hname<<" not there -- \n";return;
      }
//       histo->Rebin(5);
      histo->SetLineColor(2);
      if (mcBackgroundScaled.find(type[itype]) == mcBackgroundScaled.end()) {
	mcBackgroundScaled[type[itype]] = (TH1F*) histo->Clone("mcBackgroundScaled");
      } else {
	mcBackgroundScaled[type[itype]]->Add(histo);
      }
    }
  }
  cout << "Got "<<mcSignalCorrScaled.size() << " MC background histos\n";

  theFile->Close();
  delete theFile;

  if (systematics) {
    sprintf(name,"%s%s",absolutePath,fileWithCorrelation_gen.c_str());
    cout << "Systematics MC histos with correlation from "  << name<<endl;
    theFile = new TFile (name) ;

    for (int itype = 0;itype!=maxType;++itype) {
      sprintf(hname,"%s_%s_%s_%s", histoName, type[itype],cut, "TTbarSig");
      if (debug) cout << hname<<endl;
      TH1F* histo  = (TH1F*) gDirectory->Get(hname) ;
      if (histo==0) {
	cout << hname<<" not there -- \n";return;
      }
  //       histo->Rebin(5);
      histo->SetLineColor(2);
      if (mcSignalCorrScaled_gen.find(type[itype]) == mcSignalCorrScaled_gen.end()) {
	mcSignalCorrScaled_gen[type[itype]] = (TH1F*) histo->Clone("mcSignalCorrScaled");
      } else {
	mcSignalCorrScaled_gen[type[itype]]->Add(histo);
      }
    }

    if (systematicsPDF) {
      for (int i = 0; i!= maxPdfTemplates; ++i) {

	for (int itype = 0;itype!=maxType;++itype) {
	  sprintf(hname,"%s_%s_%s_%i", histoName, type[itype],cut, i);
	  if (debug) cout << hname<<endl;
	  TH1F* histo  = (TH1F*) gDirectory->Get(hname) ;
	  if (histo==0) {
	    cout << hname<<" not there -- \n";return;
	  }
      //       histo->Rebin(5);
	  histo->SetLineColor(2);
	  sprintf(tag,"%s_pdf%i",type[itype],i);
	  if (mcSignalCorrScaled_gen.find(tag) == mcSignalCorrScaled_gen.end()) {
	    mcSignalCorrScaled_gen[tag] = (TH1F*) histo->Clone("mcSignalCorrScaled");
	  } else {
	    mcSignalCorrScaled_gen[tag]->Add(histo);
	  }
	}


      }
    }

    theFile->Close();
    delete theFile;
  }
  sprintf(name,"%s%s",absolutePath,filenoCorrelation.c_str());
  if (debug) cout << "Load MC histos NO correlation from file "<< name<<endl;
  theFile = new TFile (name) ;


// Get the MC templates
  cout << "Will now retrieve the MC histos NO correlation\n";
  for (int itype = 0;itype!=maxType;++itype) {
    sprintf(hname,"%s_%s_%s_%s", histoName, type[itype],cut, "TTbarSig");
    if (debug) cout << hname<<endl;
    TH1F* histo  = (TH1F*) gDirectory->Get(hname) ;
    if (histo==0) {
      cout << hname<<" not there -- \n";continue;
    }
//       histo->Rebin(5);
    histo->SetLineColor(2);
    if (mcSignalUncorrScaled.find(type[itype]) == mcSignalUncorrScaled.end()) {
      mcSignalUncorrScaled[type[itype]] = (TH1F*) histo->Clone("mcSignalUncorrScaled");
    } else {
      mcSignalUncorrScaled[type[itype]]->Add(histo);
    }
    cout << histo->Integral()<< " "<<mcSignalUncorrScaled[type[itype]]->Integral()<<endl;
  }
  cout << "Got "<<mcSignalUncorrScaled.size() << " MC histos  NO correlation\n";

  theFile->Close();
  delete theFile;


  if (systematics) {
    sprintf(name,"%s%s",absolutePath,fileNoCorrelation_gen.c_str());
    cout << "Systematics MC histos with NO correlation from "  << name<<endl;
    theFile = new TFile (name) ;

    for (int itype = 0;itype!=maxType;++itype) {
      sprintf(hname,"%s_%s_%s_%s", histoName, type[itype],cut, "TTbarSig");
      if (debug) cout << hname<<endl;
      TH1F* histo  = (TH1F*) gDirectory->Get(hname) ;
      if (histo==0) {
	cout << hname<<" not there -- \n";return;
      }
  //       histo->Rebin(5);
      histo->SetLineColor(2);
      if (mcSignalUncorrScaled_gen.find(type[itype]) == mcSignalUncorrScaled_gen.end()) {
	mcSignalUncorrScaled_gen[type[itype]] = (TH1F*) histo->Clone("mcSignalCorrScaled");
      } else {
	mcSignalUncorrScaled_gen[type[itype]]->Add(histo);
      }
    }

    if (systematicsPDF) {
      for (int i = 0; i!= maxPdfTemplates; ++i) {

	for (int itype = 0;itype!=maxType;++itype) {
	  sprintf(hname,"%s_%s_%s_%i", histoName, type[itype],cut, i);
	  if (debug) cout << hname<<endl;
	  TH1F* histo  = (TH1F*) gDirectory->Get(hname) ;
	  if (histo==0) {
	    cout << hname<<" not there -- \n";return;
	  }
      //       histo->Rebin(5);
	  histo->SetLineColor(2);
	  sprintf(tag,"%s_pdf%i",type[itype],i);
	  if (mcSignalUncorrScaled_gen.find(tag) == mcSignalUncorrScaled_gen.end()) {
	    mcSignalUncorrScaled_gen[tag] = (TH1F*) histo->Clone("mcSignalCorrScaled");
	  } else {
	    mcSignalUncorrScaled_gen[tag]->Add(histo);
	  }
	}

      }
    }

    theFile->Close();
    delete theFile;
  }

  w = new RooWorkspace("w","workspace") ;
  w->import(*deltaPhi);
  w->import(*ratioCorr);
  w->import(*sample);

/// define the datasets initially

  for (int itype = 0;itype!=maxType;++itype) {

    sprintf(hname,"histo_Corr%s",type[itype]);
      if (debug) cout << hname<<endl;
      h2 = (TH1F*) mcSignalUncorrScaled[type[itype]]->Clone();
      h2->Scale(3./h2->Integral());
      h2->Add(mcSignalCorrScaled[type[itype]],-2./mcSignalCorrScaled[type[itype]]->Integral());
      h1 = (TH1F*) mcSignalCorrScaled[type[itype]]->Clone();
      h1->Scale(3./h1->Integral());
      h1->Add(mcSignalUncorrScaled[type[itype]],-2./mcSignalUncorrScaled[type[itype]]->Integral());
      
    w->import(*(new RooDataHist(hname, hname, *deltaPhi, h1)));
    sprintf(hname,"HistPdf::corr%s(deltaPhi,histo_Corr%s)", type[itype], type[itype]);
      if (debug) cout << hname<<endl;
    w->factory(hname);


    sprintf(hname,"histo_Nocor%s",type[itype]);
      if (debug) cout << hname<<endl;
      
    w->import(*(new RooDataHist(hname, hname, *deltaPhi, h2)));
    sprintf(hname,"HistPdf::nocor%s(deltaPhi,histo_Nocor%s)", type[itype], type[itype]);
      if (debug) cout << hname<<endl;
    w->factory(hname);

    if (systematics) {
   
      {TH1F* h2S = (TH1F*) mcSignalUncorrScaled_gen[type[itype]]->Clone();
      h2S->Scale(3./h2S->Integral());
      h2S->Add(mcSignalCorrScaled_gen[type[itype]],-2./mcSignalCorrScaled_gen[type[itype]]->Integral());
      TH1F* h1S = (TH1F*) mcSignalCorrScaled_gen[type[itype]]->Clone();
      h1S->Scale(3./h1S->Integral());
      h1S->Add(mcSignalUncorrScaled_gen[type[itype]],-2./mcSignalUncorrScaled_gen[type[itype]]->Integral());

      sprintf(hname,"histo_Corr_gen%s",type[itype]);
      w->import(*(new RooDataHist(hname, hname, *deltaPhi, h1S)));
      sprintf(hname,"HistPdf::corr_gen%s(deltaPhi,histo_Corr_gen%s)", type[itype], type[itype]);
      w->factory(hname);
   
   
      sprintf(hname,"histo_Nocor_gen%s",type[itype]);
      w->import(*(new RooDataHist(hname, hname, *deltaPhi, h2S)));
      sprintf(hname,"HistPdf::nocor_gen%s(deltaPhi,histo_Nocor_gen%s)", type[itype], type[itype]);
      w->factory(hname);}

      if (systematicsPDF) {
	for (int i = 0; i!= maxPdfTemplates; ++i) {
	  sprintf(tag,"%s_pdf%i",type[itype],i);
	  TH1F* h2S = (TH1F*) mcSignalUncorrScaled_gen[tag]->Clone();
	  h2S->Scale(3./h2S->Integral());
	  h2S->Add(mcSignalCorrScaled_gen[tag],-2./mcSignalCorrScaled_gen[tag]->Integral());
	  TH1F* h1S = (TH1F*) mcSignalCorrScaled_gen[tag]->Clone();
	  h1S->Scale(3./h1S->Integral());
	  h1S->Add(mcSignalUncorrScaled_gen[tag],-2./mcSignalUncorrScaled_gen[tag]->Integral());

	  sprintf(hname,"histo_Corr_gen%s",tag);
	  w->import(*(new RooDataHist(hname, hname, *deltaPhi, h1S)));
	  sprintf(hname,"HistPdf::corr_gen%s(deltaPhi,histo_Corr_gen%s)", tag, tag);
	  w->factory(hname);


	  sprintf(hname,"histo_Nocor_gen%s",tag);
	  w->import(*(new RooDataHist(hname, hname, *deltaPhi, h2S)));
	  sprintf(hname,"HistPdf::nocor_gen%s(deltaPhi,histo_Nocor_gen%s)", tag, tag);
	  w->factory(hname);
	}
      }


    }

    sprintf(hname,"histo_bck%s",type[itype]);
      if (debug) cout << hname<<endl;
    w->import(*(new RooDataHist(hname, hname, *deltaPhi, mcBackgroundScaled[type[itype]])));
    sprintf(hname,"HistPdf::background%s(deltaPhi,histo_bck%s)", type[itype], type[itype]);
      if (debug) cout << hname<<endl;
    w->factory(hname);

//     if (systematics) {
//       sprintf(hname,"histo_bck_gen%s",type[itype]);
//       w->import(*(new RooDataHist(hname, hname, *deltaPhi, mcTotalBackgroundHistoScaled_gen[type[itype]])));
//       sprintf(hname,"HistPdf::background_gen%s(deltaPhi,histo_bck%s)", type[itype], type[itype]);
//       w->factory(hname);
//     }

    sprintf(hname,"SUM::signal%s( ratioCorr*corr%s, nocor%s )",
    	    type[itype], type[itype], type[itype]);
    if (debug) cout << hname<<endl;
    w->factory(hname);


    if (useRatio) {
      sprintf(hname,"SUM::model%s( ratioSignal%s[-1,2]*signal%s, background%s )",
    	      type[itype], type[itype], type[itype], type[itype]);
    } else {
      float totalBackground = (float) mcBackgroundScaled[type[itype]]->Integral();
      sprintf(hname,"nbrBackgroundEvents%s",type[itype]);
      if (debug) cout << hname<<endl;
      nbrBackgroundEvents[type[itype]] = new RooRealVar(hname,hname,totalBackground,0,10000);
      w->import(*nbrBackgroundEvents[type[itype]]);
      w->var(hname)->setConstant(1);
      fixBckg = 3;

      sprintf(hname,"SUM::model%s( Nsig%s[0,20000]*signal%s, nbrBackgroundEvents%s*background%s )",
    	  type[itype], type[itype], type[itype], type[itype], type[itype]);

    }


    if (debug) cout << hname<<endl;
    w->factory(hname);
//             sprintf(hname,"Nsig%s", type[itype]);
//             w->var(hname)->setConstant(1);
//             w->var(hname)->setVal(mcSignalUncorrScaled[type[itype]]->Integral());
// 	    cout << hname <<" "<< (mcSignalUncorrScaled[type[itype]]->Integral())<<endl;


  }
  
  char name2[300];
  if (debug) cout << "Build simultaneous pdf"<<endl;
  sprintf(name2,"SIMUL::model(sample");
    for (int j =0; j<maxType;++j) {
      sprintf(name,"%s,%s=model%s", name2, type[j], type[j]);
      sprintf(name2,"%s", name);
    }
     sprintf(name,"%s)",name2);
     if (debug) cout << name<<endl;
    
    w->factory(name);

  pdffit = w->pdf("model") ;
  pdffit->Print();

  w->Print();
  myFitResults_all = 0;

  //random = new ROOT::Math::GSLRngMT();
//   random->Initialize();
//   random->SetSeed(23);

}


void SpinFit::assembleDatasetsHisto()
{
  if (debug) cout << "Combine "<< maxType<<" categories\n";
  if (data) delete data;
  map<string,RooDataHist*> mapToImport;
  for (int i = 0; i!=maxType; ++i) {
    mapToImport[type[i]] = new RooDataHist(type[i], type[i], *deltaPhi, dataHistos[type[i]]);
    if (fixedSample) nTotSample[i] = (int) dataHistos[type[i]]->Integral();
  }
  data = new RooDataHist("data","combined data",*deltaPhi,Index(*sample),
	Import(mapToImport) ) ;
  for (int i = 0; i!=maxType; ++i) delete mapToImport[type[i]];
  data->Print();
}

void SpinFit::assembleDatasets(map<string,RooDataHist*> mapToImport)
{
  if (data) delete data;
  if (debug) cout << "Combine "<< maxType<<" categories\n";
  data = new RooDataHist("data","combined data",*deltaPhi,Index(*sample),
	Import(mapToImport) ) ;
  for (int i = 0; i!=maxType; ++i) delete mapToImport[type[i]];
  //data->Print();
}

double SpinFit::fit()
{
  if (myFitResults_all)  {
    delete myFitResults_all;
    delete nll;
  }
cout <<"Start fit\n";
  nll = pdffit->createNLL(*data,NumCPU(2));
cout <<"Start fit\n";

  myFitResults_all = pdffit->fitTo(*data, Save(), InitialHesse(1),PrintLevel(-1),PrintEvalErrors(-1));
   //myFitResults_all = pdffit->fitTo(*data, Save());//, InitialHesse(1));
cout <<"Start fit\n";

  if (debug) cout << "chi2 " << myFitResults_all->minNll() <<endl;
  return myFitResults_all->minNll();
}

void SpinFit::printAllValues()
{
  char hname[150];
  for (int itype = 0;itype!=maxType;++itype) {
    if (useRatio) {
      sprintf(hname,"ratioSignal%s", type[itype]);
      cout << hname <<" "<<w->var(hname)->getVal() <<" for "<<type[itype]<<endl;
    } else {
      sprintf(hname,"nbrBackgroundEvents%s", type[itype]);
      cout << hname <<" "<<w->var(hname)->getVal() 
      << " +- " << w->var(hname)->getError() 
      <<" for "<<type[itype]<<", ";
      sprintf(hname,"Nsig%s", type[itype]);
      cout << hname <<" "<<w->var(hname)->getVal() 
      << " +- " << w->var(hname)->getError() 
      <<" for "<<type[itype]<<endl;
    }
  }

      cout << "ratioCorr: "<<w->var("ratioCorr")->getVal() 
      << " +- " << w->var("ratioCorr")->getError() <<" ( "<< w->var("ratioCorr")->getErrorLo()
  << " , +"<<w->var("ratioCorr")->getErrorHi()<<" )\n";
  
  
      cout << "ratioCorrverted con: "<< 5*(w->var("ratioCorr")->getVal())-2 
      << " +- " << 5*w->var("ratioCorr")->getError() <<" ( "<< 5*w->var("ratioCorr")->getErrorLo()
  << " , +"<<5*w->var("ratioCorr")->getErrorHi()<<" )\n";
  
  
    cout << "the spin corr is  " << (5*(w->var("ratioCorr")->getVal())-2)*0.326 << "  pm " << 5*0.326*w->var("ratioCorr")->getError()<< endl;
  //ratioCorr->Print();

}
void SpinFit::fixBackground(int fix)
{
  char hname[50];
   switch ( fix )
   {
     case 0: 
       for (int itype = 0;itype!=maxType;++itype) {
          sprintf(hname,"nbrBackgroundEvents%s", type[itype]);
          w->var(hname)->setConstant(0);
       }
       break;
     case 1:
       for (int itype = 0;itype!=maxType;++itype) {
	 if (!useRatio) {
            sprintf(hname,"nbrBackgroundEvents%s", type[itype]);
            w->var(hname)->setConstant(1);
            w->var(hname)->setVal(0.);
	 } else {
            sprintf(hname,"ratioSignal%s", type[itype]);
            w->var(hname)->setConstant(1);
            w->var(hname)->setVal(1.);
	 }
       }
       break;
     case 2: cout <<"Use method fixBackground(ee, em, mm) - Nothing done now\n";
       return;
       break;
     case 3:

       if (useRatio) {
	 char tag[10];
	 for (int itype = 0;itype!=maxType;++itype) {
	   sprintf(tag,"%s172.5",type[itype]);
	   sprintf(hname,"ratioSignal%s", type[itype]);
     //       cout << tag <<" "<<hname<<" "<<(mcSignalTemplHistosScaled[tag]->Integral()/mcTotalBackgroundHistoScaled[type[itype]]->Integral())
     //       <<" "<<mcSignalTemplHistosScaled[tag]->Integral()<<" "<<mcTotalBackgroundHistoScaled[type[itype]]->Integral()<<endl;
	   w->var(hname)->setVal(mcSignalCorrScaled[type[itype]]->Integral()/
      	     (mcSignalCorrScaled[type[itype]]->Integral()+mcBackgroundScaled[type[itype]]->Integral()));
	   w->var(hname)->setConstant(1);
     //        cout << w->var(hname)->getVal() <<" for "<<type[itype]<<", \n";
	}
       } else {
	 for (int itype = 0;itype!=maxType;++itype) {
            sprintf(hname,"nbrBackgroundEvents%s", type[itype]);
            w->var(hname)->setVal(mcBackgroundScaled[type[itype]]->Integral());
            w->var(hname)->setConstant(1);
	 }
       }

       break;
     default:
       cout << "input "<<fix<<" unknown!\n";
   }
  fixBckg = fix;
  fixBackground();
}

int SpinFit::fixBackground()
{
  cout << "Choice can be: \n";
  cout << "   0: Not fixed\n";
  cout << "   1: Fixed to 0.\n";
  cout << "   2: Fixed to external input\n";
  cout << "   3: Fixed to MC\n";
  cout << "Present setting: " << fixBckg<<endl;
// cout << "Needs to be modified - Exiting\n";return 0;
  if (fixBckg>=2) {
    char hname[50];
    cout << "   Fixed to: ";
    for (int itype = 0;itype!=maxType;++itype) {
       if (useRatio) sprintf(hname,"ratioSignal%s", type[itype]);
         else sprintf(hname,"nbrBackgroundEvents%s", type[itype]);
       cout << w->var(hname)->getVal() <<" for "<<type[itype]<<", ";
    }
    cout <<endl;
  }
  return fixBckg;
}



pair<double,double> SpinFit::findMin()
{
  RooRealVar* a = (RooRealVar*) myFitResults_all->floatParsFinal().find("ratioCorr");
  cout << a->getVal()<<" "<<a->getError()<<endl;
  return pair<double,double>(a->getVal(),a->getError());
}

void SpinFit::printLL(bool toFile, char* name)
{
  c_min->cd();

// Plot likelihood scan in parameter frac
  RooPlot* frame = ratioCorr->frame(Bins(20)) ;
  nll->plotOn(frame,ShiftToZero()) ;
  frame->Draw();
  frame->SetMinimum(1e-5);
  frame->SetMaximum(250);

  cmsprelim();

  if (toFile) {
    char name2[200],name3[200];
    if (strlen(name)==0) {
      sprintf(name2,"%s/fit_result/ll", absolutePath);
      typeTag(name2);
    } else {
      sprintf(name2,"%s", name);
    }
    sprintf(name3,"%s.C",name2); c_min->Print(name3);
    sprintf(name3,"%s.png",name2); c_min->Print(name3);
    sprintf(name3,"%s.pdf",name2); c_min->Print(name3);
  }
}


void SpinFit::printFit(bool withFit, bool toFile, char* name)
{
  char hname[50];
  map<string, RooPlot*> frames;
  for (int itype = 0;itype!=maxType;++itype) {
    sprintf(hname,"#Delta #phi %s", type[itype]);
    RooPlot* frame = deltaPhi->frame(Title(hname)) ;
    sprintf(hname,"sample==sample::%s", type[itype]);
    data->plotOn(frame,Cut(hname), DrawOption("ZpE") ) ;
    frames[type[itype]] = frame;
  } 

  RooPlot* frame4 = deltaPhi->frame(Title("#Delta #phi distribution")) ;
  data->plotOn(frame4, DrawOption("Zp")) ;

  if (withFit) {
    for (int itype = 0;itype!=maxType;++itype) {
      sprintf(hname,"model%s", type[itype]);
      pdffit->plotOn(frames[type[itype]],Components(hname),ProjWData(*sample,*data), LineColor(kBlack)) ; // Normalization(10,RooAbsReal::Relative)
      sprintf(hname,"signal%s", type[itype]);
      pdffit->plotOn(frames[type[itype]],Components(hname),ProjWData(*sample,*data),LineStyle(kDotted), LineColor(kBlue)) ;
      sprintf(hname,"corr%s", type[itype]);
      pdffit->plotOn(frames[type[itype]],Components(hname),ProjWData(*sample,*data),LineStyle(kDotted), LineColor(kGreen)) ;
      sprintf(hname,"nocor%s", type[itype]);
      pdffit->plotOn(frames[type[itype]],Components(hname),ProjWData(*sample,*data),LineStyle(kDotted), LineColor(kViolet)) ;
      sprintf(hname,"background%s", type[itype]);
      pdffit->plotOn(frames[type[itype]],Components(hname),ProjWData(*sample,*data),LineStyle(kDashed), LineColor(kRed)) ;
    }

    sprintf(hname,"model");
//     w->pdf(hname)->plotOn(frame4,Components(hname)) ;
//     w->pdf(hname)->Print() ;
    pdffit->plotOn(frame4,Components(hname),ProjWData(*sample,*data), LineColor(kBlack)) ;
  }

  TCanvas* c;
  if (maxType>3) {
    c = new TCanvas("mass","Top ",1200,1200) ;
    c->Divide(3,3) ;
  } else {
    c = new TCanvas("mass","Top ",800,800) ;
    c->Divide(2,2) ;
  }
  int i=1;
  for (int itype = 0;itype!=maxType;++itype) {
    frames[type[itype]]->SetMinimum(1e-5);
    frames[type[itype]]->GetYaxis()->SetTitle("Events"); 
    c->cd(i); 
    gPad->SetLeftMargin(0.15);
    sprintf(hname,"#Delta #phi ll %s", type[itype]);
    frames[type[itype]]->GetXaxis()->SetTitle(hname); 
    frames[type[itype]]->GetYaxis()->SetTitleOffset(1.1) ; 
    frames[type[itype]]->Draw() ;
    ++i;
    cmsprelim();
  }
  
  if (i<10){
    frame4->SetMinimum(1e-5);
    frame4->GetYaxis()->SetTitle("Events");
    c->cd(i) ; gPad->SetLeftMargin(0.15) ; frame4->GetXaxis()->SetTitle("#Delta #phi - All");frame4->GetYaxis()->SetTitleOffset(1.1) ; frame4->Draw() ;
    cmsprelim();
  }

  if (toFile) {
    char name2[200],name3[200];
    if (strlen(name)==0) {
      if (withFit) sprintf(name2,"%s/fit_result/deltaPhi", absolutePath);
      else sprintf(name2,"%s/fit_result/deltaPhi", absolutePath);
      typeTag(name2);
    } else {
      sprintf(name2,"%s", name);
    }
    sprintf(name3,"%s.C",name2); c->Print(name3);
    sprintf(name3,"%s.png",name2); c->Print(name3);
    sprintf(name3,"%s.pdf",name2); c->Print(name3);
  }
//   TCanvas* c2 = new TCanvas("mass","Top ",300,300) ;
//   gPad->SetLeftMargin(0.15) ; frame4->GetXaxis()->SetTitle("Reconstructed mass - All");frame4->GetYaxis()->SetTitleOffset(1.1) ; frame4->Draw() ;
}

void SpinFit::typeTag(char *name)
{
  char name2[200];
  for (int j =0; j<maxType;++j) {
    sprintf(name2,"%s%s", name, type[j]);
    if (j<(maxType-1)) sprintf(name,"%s_", name2);
    else sprintf(name,"%s", name2);
  }
  cout << name<<endl;
}
void SpinFit::cmsprelim()
{
    TLatex* text1 = new TLatex(3.570061,25,"CMS Preliminary, 5.0 fb^{-1} at #sqrt{s} = 7 TeV");
    text1->SetNDC();
    text1->SetTextAlign(13);
    text1->SetX(0.167);
    text1->SetY(0.995);
    text1->SetTextFont(42);
    text1->SetTextSizePixels(24);
    text1->Draw();
}


void SpinFit::calibration(int numberPE, float lowerLimit, float upperLimit,
	int numberBins)
{
  time_t start,end;
  time (&start);

  nFitFailed = 0;
  nFitTried = 0;
  TVectorD x(numberBins), y(numberBins), ex(numberBins), ey(numberBins);

  int i;
  float ratioCorrGen;
  char hname[50];
//   TString name = TString("calibration_");
//   name += ilumi;
//   name.ReplaceAll ( " " , "" );
//   cout <<name<<endl;
  quietFit(true);

  TFile * out = new TFile("calibration.root","RECREATE");
  TTree* pe = new TTree("PE","");
  pe->Branch("i",&i,"i/I");
  pe->Branch("ratioCorrGen",&ratioCorrGen,"ratioCorrGen/F");
  
  
  cout << "lowerLimit " << lowerLimit << endl;
  cout << "upperLimit " << upperLimit << endl;
  
  
  for (i = 0;i!=numberBins;++i)
  {
    ratioCorrGen = lowerLimit + (upperLimit-lowerLimit)/(numberBins-1)*i;
    
    do_toys(numberPE, ratioCorrGen);
    toy_mean->Fit("gaus");
    
    
    cout << "finished toy MC " << endl;
    
    x[i]=ratioCorrGen; ex[i]=0.;
    y[i]=toy_mean->GetFunction("gaus")->GetParameter(1);
    ey[i]=toy_mean->GetFunction("gaus")->GetParameter(2);
    
    
    
//     y[pts]=toy_mean->GetMean();
//     ey[pts]=toy_mean->GetRMS();

    cout << "produce plots " << endl;
    sprintf(hname,"meanMass_%i", i);
    toy_mean->Clone(hname)->Write();
    sprintf(hname,"biasMass_%i", i);
    toy_bias->Clone(hname)->Write();
    sprintf(hname,"errMass_%i", i);
    toy_error->Clone(hname)->Write();
    sprintf(hname,"pullMass_%i", i);
    toy_pull->Clone(hname)->Write();
    pe->Fill();
  }
  for (i = 0;i!=numberBins;++i)
  {
    ratioCorrGen = lowerLimit + (upperLimit-lowerLimit)/(numberBins-1)*i;
    cout << "Ratio "<< ratioCorrGen << " - Fit: " << y[i] <<" / "<<ey[i]<<endl;
  }
  if (grc!=0) {
    delete grc;
  }
  grc = new TGraphErrors(x,y,ex,ey);
  grc->SetName("grc");
  TF1* f1 = new TF1("f1","pol1",lowerLimit, upperLimit);
  grc->Fit("f1");
  TF1* f2 = new TF1("f2","pol1",lowerLimit, upperLimit);
  f2->SetParameter(0,0.);
  f2->SetParameter(1,1.);
  f2->SetLineColor(4);

  grc->Draw("a*");
  f2->Draw("same");
  grc->Write();
  f1->Write();
  f2->Write();
  pe->Write();
//   c_min->Print(name+".png");
  out->Close();
  cout << "Failed fits: "<<nFitFailed<< " / "<< nFitTried<<endl;
  time (&end);
  double dif = difftime (end,start);
  printf ("It took  %.2lf seconds to do the whole thing, %.2lf per loop.\n", dif , dif/numberPE);
}


void SpinFit::pdfSyst(int numberPE, float ratioCorrGen)
{
  time_t start,end;
  time (&start);

  nFitFailed = 0;
  nFitTried = 0;

  char hname[50];
  quietFit(true);
  TFile * out = new TFile("pdfSystematics.root","RECREATE");
  TVectorD x(maxPdfTemplates), y(maxPdfTemplates), ex(maxPdfTemplates), ey(maxPdfTemplates);
  for (int i = 0;i<maxPdfTemplates;++i)
  {
  cout << "pdf "<<i<<endl;
    do_toys(numberPE, ratioCorrGen, i);
    toy_mean->Fit("gaus");
    
    cout << "finished toy MC " << endl;
    
    x[i]=ratioCorrGen; ex[i]=0.;
    y[i]=toy_mean->GetFunction("gaus")->GetParameter(1);
    ey[i]=toy_mean->GetFunction("gaus")->GetParameter(2);
//     cout << "produce plots " << endl;
    sprintf(hname,"meanMass_%i", i);
    toy_mean->Clone(hname)->Write();
    sprintf(hname,"biasMass_%i", i);
    toy_bias->Clone(hname)->Write();
    sprintf(hname,"errMass_%i", i);
    toy_error->Clone(hname)->Write();
    sprintf(hname,"pullMass_%i", i);
    toy_pull->Clone(hname)->Write();
  }
  for (int i = 0;i<maxPdfTemplates;++i) cout << y[i]<<endl;

  double plus, minus;
  double sumPlus = 0., sumMinus = 0.;
  cout << "Default: "<< y[0]<<endl;
  for (int i = 0;i<20;++i) {
    plus  = max( max(y[2*i+1]-y[0], y[2*i+2]-y[0]), 0.);
    minus = max( max(y[0]-y[2*i+1], y[0]-y[2*i+2]), 0.);
    cout << y[2*i+1]-y[0]<<" "<<y[2*i+2]-y[0]<<endl;
    cout << i << " "<< y[2*i+1] << " "<< y[2*i+2] << " "<< minus << " "<< plus<< " "<< endl;
    sumPlus+= (plus*plus);
    sumMinus+=(minus*minus);
  }
  sumPlus = sqrt(sumPlus );
  sumMinus= sqrt(sumMinus);
  cout << "Pdf Minus " << sumMinus;
  cout << " - Plus " << sumPlus<<endl;

  out->Close();
  cout << "Failed fits: "<<nFitFailed<< " / "<< nFitTried<<endl;
  time (&end);
  double dif = difftime (end,start);
  printf ("It took  %.2lf seconds to do the whole thing, %.2lf per loop.\n", dif , dif/numberPE);

}


void SpinFit::do_toys(int numberPE, float ratioCorrGen, int pdfToUse)
{
  quietFit(true);

  cout << "do_toys: Template ratio "<< ratioCorrGen<<endl;
  
  if (toy_error!=0) {
    delete toy_mean; delete toy_error; delete toy_pull; delete toy_bias;
  }
  toy_mean   = new TH1F("mean"  ,"Ratio",300, ratioCorrGen-1., ratioCorrGen+1.);
  toy_bias   = new TH1F("bias"  ,"Ratio bias",100, -0.5, 0.5);
  toy_error  = new TH1F("error" ,"Ratio uncertainty",500, 0, 1.);
  toy_pull   = new TH1F("pull"  ,"Ratio pull",100, -5, 5);
  toy_mean->GetXaxis()->SetNdivisions(50205);
  toy_bias->GetXaxis()->SetNdivisions(50205);
  toy_mean->SetFillColor(44);
  toy_bias->SetFillColor(44);
  toy_error->SetFillColor(44);
  toy_pull->SetFillColor(44);
  int j, stat;
  for (int i=0;i<numberPE;i++){
    int countFailures = 0;
    do {
      do { j = generate_toy((ratioCorrGen+2.)/5., pdfToUse);} while (j==0);
      fit();  
      stat =0;    
//       if (stat !=0) {
// 	++countFailures;
// 	cout << "RooFit failure "<<countFailures<<endl;
//       }
    } while ((stat !=0) && countFailures<10);
    if (stat !=0) {
      cout << "Too many consecutive failures\n";
      exit(1);
    }
    nFitFailed+=countFailures;

    pair<double,double> result = findMin();

    if (result.second>0.) {
//     cout << "******************************** " <<  ratioCorrGen << "  " << result.first<<" "<<result.second<<endl;
    
      double results_f = 5*result.first-2;
      double results_s = 5*result.second;
//     cout << "******************************** " <<  ratioCorrGen << "  " << results_f << " " << results_s <<endl;
      toy_mean->Fill(results_f);
      toy_bias->Fill(results_f-ratioCorrGen);
      toy_error->Fill(results_s);
      toy_pull->Fill((ratioCorrGen-results_f)/results_s);
    } else {
      ++nFitFailed;
    }
    ++nFitTried;
  }
//   cout << "Template mass "<< mcSignalTemplMass[templateToUse]<<endl;
  cout << "Failed fits: "<<nFitFailed<< " / "<< nFitTried<<endl;
}

void SpinFit::nicePlot()
{
  c_min->cd();

  char hname[150];
  float f = sTOr(w->var("ratioCorr")->getVal());
  TH1F* sumData, *backgr, *signalNoCorr, *signalFitted, *totalFit, *signalSM;
  TH1F* totalSM, *totalNoCorr;
  for (int itype = 0;itype!=maxType;++itype) {
    sprintf(hname,"Nsig%s", type[itype]);
    float nsig = w->var(hname)->getVal();

      if (itype == 0) {
	sumData = (TH1F*) dataHistos[type[itype]]->Clone("a");
	backgr = (TH1F*) mcBackgroundScaled[type[itype]]->Clone("s");
	signalNoCorr = (TH1F*) mcSignalUncorrScaled[type[itype]]->Clone("d");
	signalNoCorr->Scale(nsig/mcSignalUncorrScaled[type[itype]]->Integral());
	signalSM = (TH1F*) mcSignalCorrScaled[type[itype]]->Clone("d");
	signalSM->Scale(nsig/mcSignalCorrScaled[type[itype]]->Integral());
	signalFitted = (TH1F*) mcSignalUncorrScaled[type[itype]]->Clone("f");
	signalFitted->Scale((1-f)*nsig/mcSignalUncorrScaled[type[itype]]->Integral());
	signalFitted->Add(mcSignalCorrScaled[type[itype]],(f)*nsig/mcSignalCorrScaled[type[itype]]->Integral());
      } else {
	sumData->Add(dataHistos[type[itype]]);
	backgr->Add(mcBackgroundScaled[type[itype]]);
	signalNoCorr->Add(mcSignalUncorrScaled[type[itype]],nsig/mcSignalUncorrScaled[type[itype]]->Integral());
	signalFitted->Add(mcSignalUncorrScaled[type[itype]],(1-f)*nsig/mcSignalUncorrScaled[type[itype]]->Integral());
	signalFitted->Add(mcSignalCorrScaled[type[itype]],(f)*nsig/mcSignalCorrScaled[type[itype]]->Integral());
	signalSM->Add(mcSignalCorrScaled[type[itype]],nsig/mcSignalCorrScaled[type[itype]]->Integral());
      }
    }
  totalFit = (TH1F*) backgr->Clone("dataHistos");
  totalFit->Add(signalFitted);
  totalSM = (TH1F*) backgr->Clone("dataHistos");
  totalSM->Add(signalSM);
  totalNoCorr = (TH1F*) backgr->Clone("dataHistos");
  totalNoCorr->Add(signalNoCorr);

  totalSM->SetLineColor(kViolet);
  totalSM->SetLineStyle(2);
  totalSM->SetLineWidth(1.5);
  totalNoCorr->SetLineColor(kBlue);
  totalNoCorr->SetLineStyle(3);
  totalNoCorr->SetLineWidth(1.5);

  sumData->SetLineColor(kBlack);
  sumData->SetMarkerStyle(22);
  sumData->GetXaxis()->SetTitle("#Delta #phi_{l+l-}");
  sumData->GetYaxis()->SetTitle("Events");
  sumData->GetYaxis()->SetTitleSize(0.045);
  sumData->GetYaxis()->SetLabelSize(0.04);
  sumData->GetYaxis()->SetTitleOffset(1.6);
  sumData->GetYaxis()->SetTitleFont(42);
  sumData->GetYaxis()->SetLabelFont(42);
  sumData->GetXaxis()->SetTitleSize(0.045);
  sumData->GetXaxis()->SetLabelSize(0.04);
  sumData->GetXaxis()->SetTitleFont(42);
  sumData->GetXaxis()->SetLabelFont(42);
  totalFit->SetLineColor(kBlack);
  backgr->SetLineColor(kRed);
  backgr->SetLineStyle(1);
  signalNoCorr->SetLineColor(kViolet);
  signalNoCorr->SetLineStyle(2);
  signalFitted->SetLineColor(kBlue);
  signalFitted->SetLineStyle(3);
    signalNoCorr->SetLineWidth(1.5);
    signalFitted->SetLineWidth(1.5);
    totalFit->SetLineWidth(1.5);
    backgr->SetLineWidth(1.5);

  TLegend * leg = new TLegend(.55,.28,.90,.53);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->AddEntry(sumData,     "Data", "p");
  leg->AddEntry(totalFit,    "Fit", "l");
  leg->AddEntry(totalNoCorr, "Backg. + t#bar{t} without correlation", "l");
  leg->AddEntry(totalSM,     "Backg. + t#bar{t} with correlation", "l");
  leg->AddEntry(backgr, "Background", "l");
//   leg->AddEntry(signalNoCorr, "t#bar{t}, no correlation", "l");
//   leg->AddEntry(signalFitted, "t#bar{t}, SM", "l");
  
  totalNoCorr->SetLineWidth(2.);
  totalSM->SetLineWidth(2.);
  
  sumData->SetMinimum(0);
  sumData->SetMaximum(2000);
  sumData->Draw("e p x0");
  totalFit->Draw("same");
  backgr->Draw("same");

//   signalNoCorr->Draw("same");
//   signalSM->Draw("same");

  totalSM->Draw("same");
  totalNoCorr->Draw("same");

  leg->Draw();

  cmsprelim();
}

int SpinFit::generate_toy(float ratioCorrGen, int pdfToUse)
{
  cout << "generate_toy: Template ratio "<< ratioCorrGen<<endl;
  if (toyDataHisto!=0) {
    delete toyDataHisto;
  }
  toyDataHisto = (TH1F*)mcSignalCorrScaled["ee"]->Clone("toyDataHisto");
  toyDataHisto->Reset();

  char hnameB[50], hnameC[50], hnameN[50];
  int n_bkg, n_sig;
  map<string,RooDataHist*> mapToImport;

  for (int itype = 0;itype!=maxType;++itype) {
//   for (int itype = 1;itype!=2;++itype) {

//     delete datasets[type[itype]];

    if (!fixedSample) {

//       if (!systematics || !systematicsPDF) sprintf(tag,"%s%g",type[itype],mcSignalTemplMass[templateToUse]);
//       else sprintf(tag,"%s_pdf%i",type[itype],templateToUse);
// 
//       if (systematics) {
// 	n_sig = _random.Poisson( (float) mcSignalCorrScaled_gen[type[itype]]->Integral());
//       } else {
// 	n_sig = _random.Poisson( (float) mcSignalCorrScaled[type[itype]]->Integral());
//       }

      // hard cooded yields. To be fine tunned.
      float moyenne =0.;
      if (itype==0) moyenne = 1724.59;
      else if (itype==1) moyenne = 6382.7;
     else if (itype==2) moyenne = 2142.39;
      
      

      	n_sig = _random.Poisson(moyenne);

      if (debug) cout << "Generate "<< n_sig << " " << type[itype]<<" signal events, ";
//       if (!systematics || !systematicsPDF) cout << "mass "<<mcSignalTemplMass[templateToUse];
//       else cout << "pdf "<<templateToUse;
      if (debug) cout << ", with Poisson mean "<< (float) mcSignalCorrScaled_gen[type[itype]]->Integral() << endl;
// 
      //n_bkg = _random.Poisson( (float) mcBackgroundScaled[type[itype]]->Integral());
      
      if(!do_bdf_systematics) n_bkg = _random.Poisson( (float) mcBackgroundScaled[type[itype]]->Integral());
      if(do_bdf_systematics)  n_bkg = _random.Poisson( 2.*(float) mcBackgroundScaled[type[itype]]->Integral());
      
    
      
      
      if (debug) cout << "Generate "<<n_bkg << " " << " background events";
      if (debug) cout << ", with Poisson mean "<< mcBackgroundScaled[type[itype]]->Integral() << endl;
    } else {
//       vector<double> probabilities;
      double bkgProb;
//       double totalProb;
      if (systematics) {
// 	if (systematicsPDF) sprintf(tag,"%s_pdf%i",type[itype],templateToUse);
//	bkgProb = mcBackgroundScaled[type[itype]]->Integral() / nTotSample[itype];
	bkgProb = mcBackgroundScaled[type[itype]]->Integral() / 
		(mcBackgroundScaled[type[itype]]->Integral() + mcSignalCorrScaled_gen[type[itype]]->Integral());
      } else {
	bkgProb = mcBackgroundScaled[type[itype]]->Integral() / 
		(mcBackgroundScaled[type[itype]]->Integral() + mcSignalCorrScaled[type[itype]]->Integral());
//		nTotSample[itype];
      }
      if (debug) cout << bkgProb<<" "<<mcBackgroundScaled[type[itype]]->Integral()<<endl;
      if (debug) cout << nTotSample[itype]<<" "<<nTotSample[itype]*bkgProb<<endl;
//       vector<unsigned int> evtsPerChannel = random->Multinomial(nTotSample[itype], probabilities);
      n_bkg = _random.Binomial(nTotSample[itype], bkgProb);
      //n_sig = _random.Binomial(nTotSample[itype], nTotSample[itype] - n_bkg);
      n_sig = nTotSample[itype] - n_bkg;
            sprintf(hnameB,"nbrBackgroundEvents%s", type[itype]);
            w->var(hnameB)->setVal(n_bkg);
            w->var(hnameB)->setConstant(1);

    }

    generatedSignal[type[itype]] = n_sig;
    generatedBkg[type[itype]] = n_bkg;
    if (debug) cout << "Generate "<< n_sig << " " << type[itype]<<" signal events  ";
//       if (systematics && systematicsPDF)  else cout << "pdf "<<templateToUse;
    if (debug) cout << "and "<< generatedBkg[type[itype]] << " " << " background events\n";
//     cout << "with binomial prob "<< 1-bkgProb << endl;

    if (generatedBkg[type[itype]] < 0.){
      cout << "ERROR: Number of background events < 0.\n";
      assert(false);
    }

    if (systematics) {
      if (!systematicsPDF) {
	sprintf(hnameC,"corr_gen%s", type[itype]);
	sprintf(hnameN,"nocor_gen%s", type[itype]);
      } else {
	sprintf(hnameC,"corr_gen%s_pdf%i", type[itype], pdfToUse);
	sprintf(hnameN,"nocor_gen%s_pdf%i", type[itype], pdfToUse);
      }
    } else {
      sprintf(hnameC,"corr%s", type[itype]);
      sprintf(hnameN,"nocor%s", type[itype]);
    }
    //cout << hnameC<<" "<<hnameN<<endl;
    sprintf(hnameB,"background%s", type[itype]);
    mapToImport[type[itype]] = new RooDataHist(type[itype], type[itype], *deltaPhi);
    int nbrCor = (int) (generatedSignal[type[itype]]*ratioCorrGen);
    // cout << "nbrCor "<< nbrCor<<" "<< ratioCorrGen<<endl;
    if (nbrCor>0) {
      generateComponent(mapToImport[type[itype]], hnameC, nbrCor);
      generateComponent(mapToImport[type[itype]], hnameN, generatedSignal[type[itype]]-nbrCor);
    } else {
      generateComponent(mapToImport[type[itype]], hnameN, generatedSignal[type[itype]]-nbrCor);
      generateComponent(mapToImport[type[itype]], hnameC, nbrCor);
    }
    generateComponent(mapToImport[type[itype]], hnameB, generatedBkg[type[itype]]);
  }

  assembleDatasets(mapToImport);

  data->fillHistogram(toyDataHisto, *deltaPhi);

  return data->numEntries();
}

void SpinFit::generateComponent(RooDataHist *rooDataHist, char * hname, int nbr)
{
  if (nbr!=0) {
  //cout << hname<<" "<<nbr<<endl;
    RooDataHist* data2 =  w->pdf(hname)->generateBinned(*deltaPhi, abs(nbr));
    //rooDataHist->Print();
    if (nbr>0) rooDataHist->add(*data2);
    else  rooDataHist->add(*data2,(const RooFormulaVar*) 0,-1.);
    //rooDataHist->Print();
    delete data2;
  }
}

int main(int argc, const char* argv[])
{

  if (argc<1) {
    cout << "Not enough arguments given: lumi and background.\nExiting\n";
    return 1;
  }

  SpinFit *j;
  if (systematics) {
   if (argc<3) {
     cout << "Not enough arguments given: File for systematics.\nExiting\n";
     return 1;
   }
    j = new SpinFit(argv[1],argv[2]);
  } else {
    j = new SpinFit();
  }   
   //j->calibration(10,-1.,2.,31);
  j->pdfSyst(1000,0.74);
}
