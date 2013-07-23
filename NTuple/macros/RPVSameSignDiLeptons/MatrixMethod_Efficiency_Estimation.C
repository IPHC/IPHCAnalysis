#include <iomanip>
#include <iostream>
#include "NTFormat/interface/NTEvent.h"

//NTupleAnalysis classes
#include "Selection/interface/Selection.h"
#include "Selection/interface/SSDiLeptonSelection.h"
#include "Tools/interface/Dataset.h"
#include "Tools/interface/AnalysisEnvironmentLoader.h"
#include "BckgdEstimation/interface/MMEstimation.h"
#include "Tools/interface/PUWeighting.h"


#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>

using namespace IPHCTree;
using namespace std;

void ComputeEfficiencyForMC(TH1F* isotight, TH1F* isoloose, TH1F* unweighted_isotight, TH1F* unweighted_isoloose, int N_bins, TH1F* &signaleff){

  for(unsigned int i=0; i<N_bins; i++){
     float num_s_w = (float)isotight->GetBinContent(i+1);
     float den_s_w = (float)isoloose->GetBinContent(i+1);
     float num_s = (float)unweighted_isotight->GetBinContent(i+1);
     float den_s = (float)unweighted_isoloose->GetBinContent(i+1);
     float eff_s = 0.;
     float err_s = 0.;
     if(den_s != 0. && den_s_w !=0.) {
      eff_s = sqrt(num_s_w/den_s_w);
      err_s = (0.5)*sqrt((1-(num_s/den_s))/den_s);
     }
     signaleff->SetBinContent(i+1, eff_s);
     signaleff->SetBinError(i+1, err_s);
  }

}

void ComputeEfficiencyForMC2D(TH2F* isotight, TH2F* isoloose, TH2F* unweighted_isotight, TH2F* unweighted_isoloose, 
 int N_binsX, int N_binsY, TH2F* &signaleff){

  for(unsigned int i=0; i<N_binsX; i++){
    for(unsigned int j=0; j<N_binsY; j++){
       float num_s_w = (float)isotight->GetBinContent(i+1, j+1);
       float den_s_w = (float)isoloose->GetBinContent(i+1, j+1);
       float num_s = (float)unweighted_isotight->GetBinContent(i+1, j+1);
       float den_s = (float)unweighted_isoloose->GetBinContent(i+1, j+1);
       float eff_s = 0.;
       float err_s = 0.;
       if(den_s != 0. && den_s_w !=0.) {
	eff_s = sqrt(num_s_w/den_s_w);
	err_s = (0.5)*sqrt((1-(num_s/den_s))/den_s);
       }
       signaleff->SetBinContent(i+1, j+1, eff_s);
       signaleff->SetBinError(i+1, j+1, err_s);
    }
  }

}

void ComputeEfficiencyForData(TH1F* isotight, TH1F* isoloose, int N_bins, TH1F* &signaleff){

  for(unsigned int i=0; i<N_bins; i++){
     float num_s = (float)isotight->GetBinContent(i+1);
     float den_s = (float)isoloose->GetBinContent(i+1);
     float eff_s = 0.;
     float err_s = 0.;
       if(den_s != 0.) {
	eff_s = num_s/den_s;
	err_s = (sqrt(num_s*(1-(num_s/den_s)))/den_s);
       }     
     signaleff->SetBinContent(i+1, eff_s);
     signaleff->SetBinError(i+1, err_s);
     //float err_s_sys = (float)(eff_s-SignalEfficiencyNJets_MC->GetBinContent(i+1));
     //if(err_s_sys < 0.) err_s_sys = (-1)*err_s_sys;
     //SignalEfficiencyNJetsSyst->SetBinContent(i+1, eff_s);
     //SignalEfficiencyNJetsSyst->SetBinError(i+1, err_s+(err_s_sys/2));
  }

}



void ComputeFakeRateForMC(TH1F* isotighteff, TH1F* isolooseeff, TH1F* isotightfake, TH1F* isoloosefake, 
 TH1F* unweighted_isotighteff, TH1F* unweighted_isolooseeff, TH1F* unweighted_isotightfake, TH1F* unweighted_isoloosefake, 
 int N_bins, TH1F* &fakerate){

  for(unsigned int i=0; i<N_bins; i++){
     float num_s_w = (float)isotighteff->GetBinContent(i+1);
     float den_s_w = (float)isolooseeff->GetBinContent(i+1);
     float num_f_w = (float)isotightfake->GetBinContent(i+1);
     float den_f_w = (float)isoloosefake->GetBinContent(i+1);
     float num_s = (float)unweighted_isotighteff->GetBinContent(i+1);
     float den_s = (float)unweighted_isolooseeff->GetBinContent(i+1);
     float num_f = (float)unweighted_isotightfake->GetBinContent(i+1);
     float den_f = (float)unweighted_isoloosefake->GetBinContent(i+1);
     float eff_f = 0.;
     float err_f = 0.;
     if(num_s != 0. && den_f != 0. && den_s !=0. && den_s_w != 0. && den_f_w != 0.) {
      eff_f = (num_f_w/den_f_w)/sqrt(num_s_w/den_s_w);
      err_f = (sqrt(den_s/num_s)*sqrt(num_f*(1-(num_f/den_f)))/den_f) + (((num_f*den_s)/(den_f*num_s))*0.5*sqrt((1-(num_s/den_s))/den_s));
     }
     fakerate->SetBinContent(i+1, eff_f);
     fakerate->SetBinError(i+1, err_f);
  }

}

void ComputeFakeRateForMC2D(TH2F* isotighteff, TH2F* isolooseeff, TH2F* isotightfake, TH2F* isoloosefake, 
 TH2F* unweighted_isotighteff, TH2F* unweighted_isolooseeff, TH2F* unweighted_isotightfake, TH2F* unweighted_isoloosefake, 
 int N_binsX, int N_binsY, TH2F* &fakerate){

  for(unsigned int i=0; i<N_binsX; i++){
    for(unsigned int j=0; j<N_binsY; j++){
       float num_s_w = (float)isotighteff->GetBinContent(i+1, j+1);
       float den_s_w = (float)isolooseeff->GetBinContent(i+1, j+1);
       float num_f_w = (float)isotightfake->GetBinContent(i+1, j+1);
       float den_f_w = (float)isoloosefake->GetBinContent(i+1, j+1);
       float num_s = (float)unweighted_isotighteff->GetBinContent(i+1, j+1);
       float den_s = (float)unweighted_isolooseeff->GetBinContent(i+1, j+1);
       float num_f = (float)unweighted_isotightfake->GetBinContent(i+1, j+1);
       float den_f = (float)unweighted_isoloosefake->GetBinContent(i+1, j+1);
       float eff_f = 0.;
       float err_f = 0.;
       if(num_s != 0. && den_f != 0. && den_s !=0. && den_s_w != 0. && den_f_w != 0.) {
	eff_f = (num_f_w/den_f_w)/sqrt(num_s_w/den_s_w);
	err_f = (sqrt(den_s/num_s)*sqrt(num_f*(1-(num_f/den_f)))/den_f) + (((num_f*den_s)/(den_f*num_s))*0.5*sqrt((1-(num_s/den_s))/den_s));
       }
       fakerate->SetBinContent(i+1, j+1, eff_f);
       fakerate->SetBinError(i+1, j+1, err_f);
    }
  }

}

void ComputeFakeRateForData(TH1F* isotight, TH1F* isoloose, int N_bins, TH1F* &fakerate){

  for(unsigned int i=0; i<N_bins; i++){
     float num_f = (float)(isotight->GetBinContent(i+1));
     float den_f = (float)(isoloose->GetBinContent(i+1));
     float eff_f = 0.;
     float err_f = 0.;
     if(den_f != 0.) {
      eff_f = num_f/den_f;
      err_f = (sqrt(num_f*(1-(num_f/den_f)))/den_f);
     }
     fakerate->SetBinContent(i+1, eff_f);
     fakerate->SetBinError(i+1, err_f/*+err_f_sys*/);
     //float err_f_sys = (float)(eff_f-FakeRateNJets_MC->GetBinContent(i+1));
     //if(err_f_sys < 0.) err_f_sys = (-1)*err_f_sys;
     //FakeRateNJetsSyst[j]->SetBinContent(i+1, eff_f);
     //FakeRateNJetsSyst[j]->SetBinError(i+1, err_f+(err_f_sys/2));
  }


}


// -----------------------------------------------------------------------------------


void Efficiency_MC_Estimation ( string xmlFileName, bool mc=true, string lepton="Mu")
{
  cout<<"#########################"<<endl;
  cout<<"Beginning of the program"<<endl;
  cout<<"#########################"<<endl;
  
  //////////////////////
  // Initialisation
  //////////////////////
  float Luminosity = 0;
  float LumiError = 0;
  string PUWeightFileName;
  int DataType = 0; 

  AnalysisEnvironmentLoader anaEL (xmlFileName);
  vector < Dataset > datasets;
  anaEL.LoadSamples (datasets); // now the list of datasets written in the xml file is known
  int verbosity = -1;

  // isolations cuts
  float looseIso = 0.8 ; float tightIsoMMmu = 0.2; float tightIsoMMe = 0.17;
  float tightIso = tightIsoMMmu;
  if(lepton=="E") tightIso = tightIsoMMe;

  SSDiLeptonSelection sel; 
  anaEL.LoadSSDiLeptonSelection (sel); // now the parameters for the selection are given to the selection
  anaEL.LoadGeneralInfo(DataType, Luminosity, LumiError, PUWeightFileName, verbosity );
  //Load for PU:
  if(mc) sel.GeneratePUWeight(PUWeightFileName);

  IPHCTree::NTEvent * event = 0;

  /*PUWeighting  thePUReweighter;
  TFile* file1  = new TFile(PUWeightFileName.c_str(),"READ"); 
  TH1D *  hPUData = 0;
  hPUData         = (TH1D*)file1->Get("pileup");
  TH1F *  hPUMC   = new TH1F("pileup_MC", "pileup_MC", hPUData->GetXaxis()->GetNbins(), hPUData->GetXaxis()->GetXmin(), hPUData->GetXaxis()->GetXmax() );
  TFile* file2  = new TFile( "../data/CrossSection_pileup.root" ,"READ");
  hPUMC           = (TH1F*)file2->Get("pileup_TTbarSig");
  // histo in data, histo in Mc, use out-of-time pu in the reweighting
  cout << "get MC histo  " << endl;
  thePUReweighter.setPUHisto( hPUData, hPUMC);
  cout << "set MC histo in thePUReweighter " << endl;
  thePUReweighter.setUseOutOfTimePU(false); // set to true to use out-of-time PU*/

  //////////////////////
  //LOOP OVER THE DATASETS
  //////////////////////
  cout<<"#########################"<<endl;
  cout<<" Loop over the datasets  "<<endl;
  cout<<"#########################"<<endl;

  float NJets_Min = -0.5;
  float NJets_Max = 10.5;
  int NJets_N_bins = 11;

  float Eta_Min = -3;
  float Eta_Max = 3;
  int Eta_N_bins = 20;

  float PT_Min = 0;
  float PT_Max = 300;
  int PT_N_bins = 11;//60
  double PT_bins[] = { 10., 15., 20., 25., 30., 40., 50., 60., 80., 100., 150., 300. }; 
  

  float NVTX_Min = -0.5;
  float NVTX_Max = 30.5;
  int NVTX_N_bins = 31;


  TH1F * IsoTightForEfficiencyNJets = new TH1F ("IsoTightForEfficiencyNJets","IsoTightForEfficiencyNJets", NJets_N_bins, NJets_Min, NJets_Max);
  TH1F * IsoLooseForEfficiencyNJets = new TH1F ("IsoLooseForEfficiencyNJets","IsoLooseForEfficiencyNJets", NJets_N_bins, NJets_Min, NJets_Max);

  TH1F * IsoTightForEfficiencyEta = new TH1F ("IsoTightForEfficiencyEta","IsoTightForEfficiencyEta", Eta_N_bins, Eta_Min, Eta_Max);
  TH1F * IsoLooseForEfficiencyEta = new TH1F ("IsoLooseForEfficiencyEta","IsoLooseForEfficiencyEta", Eta_N_bins, Eta_Min, Eta_Max);

  //TH1F * IsoTightForEfficiencyPT = new TH1F ("IsoTightForEfficiencyPT","IsoTightForEfficiencyPT", PT_N_bins, PT_Min, PT_Max);
  //TH1F * IsoLooseForEfficiencyPT = new TH1F ("IsoLooseForEfficiencyPT","IsoLooseForEfficiencyPT", PT_N_bins, PT_Min, PT_Max);
  TH1F * IsoTightForEfficiencyPT1 = new TH1F ("IsoTightForEfficiencyPT1","IsoTightForEfficiencyPT1", PT_N_bins, PT_bins);
  TH1F * IsoLooseForEfficiencyPT1 = new TH1F ("IsoLooseForEfficiencyPT1","IsoLooseForEfficiencyPT1", PT_N_bins, PT_bins);
  TH1F * IsoTightForEfficiencyPT2 = new TH1F ("IsoTightForEfficiencyPT2","IsoTightForEfficiencyPT2", PT_N_bins, PT_bins);
  TH1F * IsoLooseForEfficiencyPT2 = new TH1F ("IsoLooseForEfficiencyPT2","IsoLooseForEfficiencyPT2", PT_N_bins, PT_bins); 

  TH1F * IsoTightForEfficiencyNVTX = new TH1F ("IsoTightForEfficiencyNVTX","IsoTightForEfficiencyNVTX", NVTX_N_bins, NVTX_Min, NVTX_Max);
  TH1F * IsoLooseForEfficiencyNVTX = new TH1F ("IsoLooseForEfficiencyNVTX","IsoLooseForEfficiencyNVTX", NVTX_N_bins, NVTX_Min, NVTX_Max);

  TH2F * IsoTightForEfficiencyNJetsvsPT1 = new TH2F ("IsoTightForEfficiencyNJetsvsPT1","IsoTightForEfficiencyNJetsvsPT1",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * IsoLooseForEfficiencyNJetsvsPT1 = new TH2F ("IsoLooseForEfficiencyNJetsvsPT1","IsoLooseForEfficiencyNJetsvsPT1",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * IsoTightForEfficiencyNJetsvsPT2 = new TH2F ("IsoTightForEfficiencyNJetsvsPT2","IsoTightForEfficiencyNJetsvsPT2",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * IsoLooseForEfficiencyNJetsvsPT2 = new TH2F ("IsoLooseForEfficiencyNJetsvsPT2","IsoLooseForEfficiencyNJetsvsPT2",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);

  TH2F * IsoTightForEfficiencyPT1vsPT2 = new TH2F ("IsoTightForEfficiencyPT1vsPT2","IsoTightForEfficiencyPT1vsPT2",
   PT_N_bins, PT_bins, PT_N_bins, PT_bins);
  TH2F * IsoLooseForEfficiencyPT1vsPT2 = new TH2F ("IsoLooseForEfficiencyPT1vsPT2","IsoLooseForEfficiencyPT1vsPT2",
   PT_N_bins, PT_bins, PT_N_bins, PT_bins);


  TH1F * IsoTightForFakeRateNJets = new TH1F ("IsoTightForFakeRateNJets","IsoTightForFakeRateNJets", NJets_N_bins, NJets_Min, NJets_Max);
  TH1F * IsoLooseForFakeRateNJets = new TH1F ("IsoLooseForFakeRateNJets","IsoLooseForFakeRateNJets", NJets_N_bins, NJets_Min, NJets_Max);

  TH1F * IsoTightForFakeRateEta = new TH1F ("IsoTightForFakeRateEta","IsoTightForFakeRateEta", Eta_N_bins, Eta_Min, Eta_Max);
  TH1F * IsoLooseForFakeRateEta = new TH1F ("IsoLooseForFakeRateEta","IsoLooseForFakeRateEta", Eta_N_bins, Eta_Min, Eta_Max);

  //TH1F * IsoTightForFakeRatePT = new TH1F ("IsoTightForFakeRatePT","IsoTightForFakeRatePT", PT_N_bins, PT_Min, PT_Max);
  //TH1F * IsoLooseForFakeRatePT = new TH1F ("IsoLooseForFakeRatePT","IsoLooseForFakeRatePT", PT_N_bins, PT_Min, PT_Max);
  TH1F * IsoTightForFakeRatePT1 = new TH1F ("IsoTightForFakeRatePT1","IsoTightForFakeRatePT1", PT_N_bins, PT_bins);
  TH1F * IsoLooseForFakeRatePT1 = new TH1F ("IsoLooseForFakeRatePT1","IsoLooseForFakeRatePT1", PT_N_bins, PT_bins);
  TH1F * IsoTightForFakeRatePT2 = new TH1F ("IsoTightForFakeRatePT2","IsoTightForFakeRatePT2", PT_N_bins, PT_bins);
  TH1F * IsoLooseForFakeRatePT2 = new TH1F ("IsoLooseForFakeRatePT2","IsoLooseForFakeRatePT2", PT_N_bins, PT_bins);

  TH1F * IsoTightForFakeRateNVTX = new TH1F ("IsoTightForFakeRateNVTX","IsoTightForFakeRateNVTX", NVTX_N_bins, NVTX_Min, NVTX_Max);
  TH1F * IsoLooseForFakeRateNVTX = new TH1F ("IsoLooseForFakeRateNVTX","IsoLooseForFakeRateNVTX", NVTX_N_bins, NVTX_Min,  NVTX_Max);

  TH2F * IsoTightForFakeRateNJetsvsPT1 = new TH2F ("IsoTightForFakeRateNJetsvsPT1","IsoTightForFakeRateNJetsvsPT1",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * IsoLooseForFakeRateNJetsvsPT1 = new TH2F ("IsoLooseForFakeRateNJetsvsPT1","IsoLooseForFakeRateNJetsvsPT1",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * IsoTightForFakeRateNJetsvsPT2 = new TH2F ("IsoTightForFakeRateNJetsvsPT2","IsoTightForFakeRateNJetsvsPT2",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * IsoLooseForFakeRateNJetsvsPT2 = new TH2F ("IsoLooseForFakeRateNJetsvsPT2","IsoLooseForFakeRateNJetsvsPT2",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);

  TH2F * IsoTightForFakeRatePT1vsPT2 = new TH2F ("IsoTightForFakeRatePT1vsPT2","IsoTightForFakeRatePT1vsPT2",
   PT_N_bins, PT_bins, PT_N_bins, PT_bins);
  TH2F * IsoLooseForFakeRatePT1vsPT2 = new TH2F ("IsoLooseForFakeRatePT1vsPT2","IsoLooseForFakeRatePT1vsPT2",
   PT_N_bins, PT_bins, PT_N_bins, PT_bins);


  TH1F * Unweighted_IsoTightForEfficiencyNJets = new TH1F ("Unweighted_IsoTightForEfficiencyNJets","Unweighted_IsoTightForEfficiencyNJets", NJets_N_bins, NJets_Min, NJets_Max);
  TH1F * Unweighted_IsoLooseForEfficiencyNJets = new TH1F ("Unweighted_IsoLooseForEfficiencyNJets","Unweighted_IsoLooseForEfficiencyNJets", NJets_N_bins, NJets_Min, NJets_Max);

  TH1F * Unweighted_IsoTightForEfficiencyEta = new TH1F ("Unweighted_IsoTightForEfficiencyEta","Unweighted_IsoTightForEfficiencyEta", Eta_N_bins, Eta_Min, Eta_Max);
  TH1F * Unweighted_IsoLooseForEfficiencyEta = new TH1F ("Unweighted_IsoLooseForEfficiencyEta","Unweighted_IsoLooseForEfficiencyEta", Eta_N_bins, Eta_Min, Eta_Max);

  //TH1F * Unweighted_IsoTightForEfficiencyPT = new TH1F ("Unweighted_IsoTightForEfficiencyPT","Unweighted_IsoTightForEfficiencyPT", PT_N_bins, PT_Min, PT_Max);
  //TH1F * Unweighted_IsoLooseForEfficiencyPT = new TH1F ("Unweighted_IsoLooseForEfficiencyPT","Unweighted_IsoLooseForEfficiencyPT", PT_N_bins, PT_Min, PT_Max);
  TH1F * Unweighted_IsoTightForEfficiencyPT1 = new TH1F ("Unweighted_IsoTightForEfficiencyPT1","Unweighted_IsoTightForEfficiencyPT1", PT_N_bins, PT_bins);
  TH1F * Unweighted_IsoLooseForEfficiencyPT1 = new TH1F ("Unweighted_IsoLooseForEfficiencyPT1","Unweighted_IsoLooseForEfficiencyPT1", PT_N_bins, PT_bins);
  TH1F * Unweighted_IsoTightForEfficiencyPT2 = new TH1F ("Unweighted_IsoTightForEfficiencyPT2","Unweighted_IsoTightForEfficiencyPT2", PT_N_bins, PT_bins);
  TH1F * Unweighted_IsoLooseForEfficiencyPT2 = new TH1F ("Unweighted_IsoLooseForEfficiencyPT2","Unweighted_IsoLooseForEfficiencyPT2", PT_N_bins, PT_bins);

  TH1F * Unweighted_IsoTightForEfficiencyNVTX = new TH1F ("Unweighted_IsoTightForEfficiencyNVTX","Unweighted_IsoTightForEfficiencyNVTX", NVTX_N_bins, NVTX_Min, NVTX_Max);
  TH1F * Unweighted_IsoLooseForEfficiencyNVTX = new TH1F ("Unweighted_IsoLooseForEfficiencyNVTX","Unweighted_IsoLooseForEfficiencyNVTX", NVTX_N_bins, NVTX_Min, NVTX_Max);

  TH2F * Unweighted_IsoTightForEfficiencyNJetsvsPT1 = 
   new TH2F("Unweighted_IsoTightForEfficiencyNJetsvsPT1","Unweighted_IsoTightForEfficiencyNJetsvsPT1",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * Unweighted_IsoLooseForEfficiencyNJetsvsPT1 = 
   new TH2F("Unweighted_IsoLooseForEfficiencyNJetsvsPT1","Unweighted_IsoLooseForEfficiencyNJetsvsPT1",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * Unweighted_IsoTightForEfficiencyNJetsvsPT2 = 
   new TH2F("Unweighted_IsoTightForEfficiencyNJetsvsPT2","Unweighted_IsoTightForEfficiencyNJetsvsPT2",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * Unweighted_IsoLooseForEfficiencyNJetsvsPT2 = 
   new TH2F("Unweighted_IsoLooseForEfficiencyNJetsvsPT2","Unweighted_IsoLooseForEfficiencyNJetsvsPT2",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);

  TH2F * Unweighted_IsoTightForEfficiencyPT1vsPT2 = new TH2F ("Unweighted_IsoTightForEfficiencyPT1vsPT2","Unweighted_IsoTightForEfficiencyPT1vsPT2",
   PT_N_bins, PT_bins, PT_N_bins, PT_bins);
  TH2F * Unweighted_IsoLooseForEfficiencyPT1vsPT2 = new TH2F ("Unweighted_IsoLooseForEfficiencyPT1vsPT2","Unweighted_IsoLooseForEfficiencyPT1vsPT2",
   PT_N_bins, PT_bins, PT_N_bins, PT_bins);


  TH1F * Unweighted_IsoTightForFakeRateNJets = new TH1F ("Unweighted_IsoTightForFakeRateNJets","Unweighted_IsoTightForFakeRateNJets", NJets_N_bins, NJets_Min, NJets_Max);
  TH1F * Unweighted_IsoLooseForFakeRateNJets = new TH1F ("Unweighted_IsoLooseForFakeRateNJets","Unweighted_IsoLooseForFakeRateNJets", NJets_N_bins, NJets_Min, NJets_Max);

  TH1F * Unweighted_IsoTightForFakeRateEta = new TH1F ("Unweighted_IsoTightForFakeRateEta","Unweighted_IsoTightForFakeRateEta", Eta_N_bins, Eta_Min, Eta_Max);
  TH1F * Unweighted_IsoLooseForFakeRateEta = new TH1F ("Unweighted_IsoLooseForFakeRateEta","Unweighted_IsoLooseForFakeRateEta", Eta_N_bins, Eta_Min, Eta_Max);

  //TH1F * Unweighted_IsoTightForFakeRatePT = new TH1F ("Unweighted_IsoTightForFakeRatePT","Unweighted_IsoTightForFakeRatePT", PT_N_bins, PT_Min, PT_Max);
  //TH1F * Unweighted_IsoLooseForFakeRatePT = new TH1F ("Unweighted_IsoLooseForFakeRatePT","Unweighted_IsoLooseForFakeRatePT", PT_N_bins, PT_Min, PT_Max);
  TH1F * Unweighted_IsoTightForFakeRatePT1 = new TH1F ("Unweighted_IsoTightForFakeRatePT1","Unweighted_IsoTightForFakeRatePT1", PT_N_bins, PT_bins);
  TH1F * Unweighted_IsoLooseForFakeRatePT1 = new TH1F ("Unweighted_IsoLooseForFakeRatePT1","Unweighted_IsoLooseForFakeRatePT1", PT_N_bins, PT_bins);
  TH1F * Unweighted_IsoTightForFakeRatePT2 = new TH1F ("Unweighted_IsoTightForFakeRatePT2","Unweighted_IsoTightForFakeRatePT2", PT_N_bins, PT_bins);
  TH1F * Unweighted_IsoLooseForFakeRatePT2 = new TH1F ("Unweighted_IsoLooseForFakeRatePT2","Unweighted_IsoLooseForFakeRatePT2", PT_N_bins, PT_bins);

  TH1F * Unweighted_IsoTightForFakeRateNVTX = new TH1F ("Unweighted_IsoTightForFakeRateNVTX","Unweighted_IsoTightForFakeRateNVTX", NVTX_N_bins, NVTX_Min, NVTX_Max);
  TH1F * Unweighted_IsoLooseForFakeRateNVTX = new TH1F ("Unweighted_IsoLooseForFakeRateNVTX","Unweighted_IsoLooseForFakeRateNVTX", NVTX_N_bins, NVTX_Min, NVTX_Max);

  TH2F * Unweighted_IsoTightForFakeRateNJetsvsPT1 = 
   new TH2F("Unweighted_IsoTightForFakeRateNJetsvsPT1","Unweighted_IsoTightForFakeRateNJetsvsPT1",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * Unweighted_IsoLooseForFakeRateNJetsvsPT1 = 
   new TH2F("Unweighted_IsoLooseForFakeRateNJetsvsPT1","Unweighted_IsoLooseForFakeRateNJetsvsPT1",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * Unweighted_IsoTightForFakeRateNJetsvsPT2 = 
   new TH2F("Unweighted_IsoTightForFakeRateNJetsvsPT2","Unweighted_IsoTightForFakeRateNJetsvsPT2",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * Unweighted_IsoLooseForFakeRateNJetsvsPT2 = 
   new TH2F("Unweighted_IsoLooseForFakeRateNJetsvsPT2","Unweighted_IsoLooseForFakeRateNJetsvsPT2",
   PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);

  TH2F * Unweighted_IsoTightForFakeRatePT1vsPT2 = new TH2F ("Unweighted_IsoTightForFakeRatePT1vsPT2","Unweighted_IsoTightForFakeRatePT1vsPT2",
   PT_N_bins, PT_bins, PT_N_bins, PT_bins);
  TH2F * Unweighted_IsoLooseForFakeRatePT1vsPT2 = new TH2F ("Unweighted_IsoLooseForFakeRatePT1vsPT2","Unweighted_IsoLooseForFakeRatePT1vsPT2",
   PT_N_bins, PT_bins, PT_N_bins, PT_bins);


  TH1F * SignalEfficiencyNJets = new TH1F("SignalEfficiencyNJets", "SignalEfficiencyNJets", NJets_N_bins, NJets_Min, NJets_Max);
  TH1F * SignalEfficiencyEta = new TH1F("SignalEfficiencyEta", "SignalEfficiencyEta", Eta_N_bins, Eta_Min, Eta_Max);
  TH1F * SignalEfficiencyPT1 = new TH1F("SignalEfficiencyPT1", "SignalEfficiencyPT1", PT_N_bins, PT_bins);
  TH1F * SignalEfficiencyPT2 = new TH1F("SignalEfficiencyPT2", "SignalEfficiencyPT2", PT_N_bins, PT_bins);
  TH1F * SignalEfficiencyNVTX = new TH1F("SignalEfficiencyNVTX", "SignalEfficiencyNVTX", NVTX_N_bins, NVTX_Min, NVTX_Max);
  TH1F * SignalEfficiencyPTfromPT1vsPT2 = new TH1F("SignalEfficiencyPTfromPT1vsPT2", "SignalEfficiencyPTfromPT1vsPT2", PT_N_bins, PT_bins);

  TH2F * SignalEfficiencyNJetsvsPT1 = new TH2F("SignalEfficiencyNJetsvsPT1", "SignalEfficiencyNJetsvsPT1", 
    PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * SignalEfficiencyNJetsvsPT2 = new TH2F("SignalEfficiencyNJetsvsPT2", "SignalEfficiencyNJetsvsPT2", 
    PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * SignalEfficiencyPT1vsPT2 = new TH2F("SignalEfficiencyPT1vsPT2", "SignalEfficiencyPT1vsPT2", 
    PT_N_bins, PT_bins, PT_N_bins, PT_bins);

  TH1F * FakeRateNJets = new TH1F("FakeRateNJets", "FakeRateNJets", NJets_N_bins, NJets_Min, NJets_Max);
  TH1F * FakeRateEta = new TH1F("FakeRateEta", "FakeRateEta", Eta_N_bins, Eta_Min, Eta_Max);
  TH1F * FakeRatePT1 = new TH1F("FakeRatePT1", "FakeRatePT1", PT_N_bins, PT_bins);
  TH1F * FakeRatePT2 = new TH1F("FakeRatePT2", "FakeRatePT2", PT_N_bins, PT_bins);
  TH1F * FakeRateNVTX = new TH1F("FakeRateNVTX", "FakeRateNVTX", NVTX_N_bins, NVTX_Min, NVTX_Max);
  TH1F * FakeRatePTfromPT1vsPT2 = new TH1F("FakeRatePTfromPT1vsPT2", "FakeRatePTfromPT1vsPT2", PT_N_bins, PT_bins);

  TH2F * FakeRateNJetsvsPT1 = new TH2F("FakeRateNJetsvsPT1", "FakeRateNJetsvsPT1", 
    PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * FakeRateNJetsvsPT2 = new TH2F("FakeRateNJetsvsPT2", "FakeRateNJetsvsPT2", 
    PT_N_bins, PT_bins, NJets_N_bins, NJets_Min, NJets_Max);
  TH2F * FakeRatePT1vsPT2 = new TH2F("FakeRatePT1vsPT2", "FakeRatePT1vsPT2", 
    PT_N_bins, PT_bins, PT_N_bins, PT_bins);
  

  TH1F * DiLepInvMass = new TH1F("DiLepInvMass", "DiLepInvMass", 150, 0, 150);
  TH1F * Iso_Lep = new TH1F("Iso_Lep", "Iso_Lep", 150, 0, 15);
  TH1F * PT_Leading_Lep = new TH1F("PT_Leading_Lep", "PT_Leading_Lep", 60, PT_Min, PT_Max);
  TH1F * PT_Inclusive_Lep = new TH1F("PT_Inclusive_Lep", "PT_Inclusive_Lep", 60, PT_Min, PT_Max);
  TH1F * NJets = new TH1F("NJets", "NJets", NJets_N_bins, NJets_Min, NJets_Max);
  TH1F * PT_Leading_Jet = new TH1F("PT_Leading_Jet", "PT_Leading_Jet", 60, PT_Min, PT_Max);
  TH1F * PT_Inclusive_Jet = new TH1F("PT_Inclusive_Jet", "PT_Inclusive_Jet", 60, PT_Min, PT_Max);
  TH1F * NVTX = new TH1F("NVTX", "NNTX", NVTX_N_bins, NVTX_Min, NVTX_Max);
  

  for (unsigned int d = 0; d < datasets.size (); d++) {
   datasets[d].eventTree ()->SetBranchAddress ("NTEvent",&event);
   cout << "dataset = " << datasets[d].Name() << endl;
   unsigned int nEvents = (int) (datasets[d].eventTree ()->GetEntries ());

    cout << "NEvents = " << nEvents << endl;
    float weight_init;
    weight_init = datasets[d].NormFactor();

    //LOOP OVER THE EVENTS

    for (unsigned int ievt = 0; ievt < nEvents; ievt++){
      datasets[d].eventTree ()->GetEntry (ievt);
      IPHCTree::NTTransient::InitializeAfterReading(event); // Important line to read new format files

      //Load event for the selection
      sel.LoadEvent(event);

      if(ievt%10000 == 0) cout << "number of processed events " << ievt << endl;

      float weight = 1;
      
      if(mc){
      	
	//Manage DY samples to avoid overlaps
	double dileptInvMass = 0;
	if( (event->mc.zAndDecays).size() > 0){
          TLorentzVector dilept = (event->mc.zAndDecays)[0].p4_Lep1_gen + (event->mc.zAndDecays)[0].p4_Lep2_gen;
          dileptInvMass = dilept.M();
	}
	if(datasets[d].Name()=="Zjets" && dileptInvMass < 50) continue;
	//if(datasets[d].Name()=="DYToMuMu_M-20"	   && (dileptInvMass > 50 || dileptInvMass < 20) ) continue;
	if(datasets[d].Name()=="DYToEE_M-20"        && (dileptInvMass > 50 || dileptInvMass < 20) ) continue;
	if(datasets[d].Name()=="DYToTauTau_M-20"    && (dileptInvMass > 50 || dileptInvMass < 20) ) continue;
	if(datasets[d].Name()=="DYToMuMu_M-10To20"   &&  dileptInvMass > 20) continue;
	if(datasets[d].Name()=="DYToEE_M-10To20"    &&  dileptInvMass > 20) continue;
	if(datasets[d].Name()=="DYToTauTau_M-10To20" &&  dileptInvMass > 20) continue;


        weight = weight_init;
	
        //float weight = weight_init*sel.GetPUWeight();
        /*if(thePUReweighter.getUseOutOfTimePU()){
           weight = weight_init*thePUReweighter.weight(event->pileup.intime_npu, event->general.runNb);
         }else{
	   weight = weight_init*thePUReweighter.weight(event->pileup.intime_npu);
       }*/
	
	DiLepInvMass->Fill(dileptInvMass, weight);	
      }
  
       std::vector<NTLepton> leptons;
       if(lepton=="Mu"){
         std::vector<NTMuon> muons =  sel.GetSelectedMuonsNoIso();
	 for(unsigned int i=0; i<muons.size(); i++){
           leptons.push_back( (NTLepton) muons[i]);
	 }
       }
       if(lepton=="E"){
	 std::vector<NTElectron> electrons =  sel.GetSelectedElectronsNoIso();
	 for(unsigned int i=0; i<electrons.size(); i++){
           leptons.push_back( (NTLepton) electrons[i]);
	 }
       }

       string channel = "";
       vector<NTElectron> candElec_Loose;
       vector<NTMuon> candMuon_Loose;
       vector<NTElectron> candElec_Tight;
       vector<NTMuon> candMuon_Tight;
       sel.GetLeptonPairForMM(candMuon_Loose, candElec_Loose, channel, looseIso, looseIso);
       sel.GetLeptonPairForMM(candMuon_Tight, candElec_Tight, channel, tightIso, tightIso);
       vector<NTJet> jets_Loose;
       vector<NTJet> jets_Tight;
       jets_Loose = sel.GetSelectedJets(candMuon_Loose, candElec_Loose);
       jets_Tight = sel.GetSelectedJets(candMuon_Tight, candElec_Tight);
       vector<NTVertex>   selVertices  = sel.GetSelectedVertex();
       NTMET met = sel.GetMET();
       vector<NTJet> jets = sel.GetJets(); // use to check njets in qcd data


        for (unsigned int j = 0; j < leptons.size(); j++) {
           Iso_Lep->Fill(leptons[j].RelIso03PF());
	   if(j==0) PT_Leading_Lep->Fill(leptons[j].p4.Perp());
	   PT_Inclusive_Lep->Fill(leptons[j].p4.Perp());
        }
	NVTX->Fill(selVertices.size());
	NJets->Fill(jets_Tight.size());
        for (unsigned int j = 0; j < jets_Tight.size(); j++) {
	   if(j==0) PT_Leading_Jet->Fill(jets_Tight[j].p4.Perp());
	   PT_Inclusive_Jet->Fill(jets_Tight[j].p4.Perp());
        }
	


       //------------------//
       // MC specific part //
       //------------------//
       if(mc)
       {
	  unsigned int loose_iso_counter_eff = 0;
	  unsigned int tight_iso_counter_eff = 0;

	  unsigned int loose_iso_counter_fake = 0;
	  unsigned int tight_iso_counter_fake = 0;
	  
	  string trigger = "";
	  if(lepton=="Mu") trigger = "mumu";
	  if(lepton=="E") trigger = "ee";

          // For efficiency
	  if(
             (lepton=="Mu" && (event->mc.TMEME == 20 || event->mc.TMEME == 11010 || event->mc.TMEME == 22000)) || //mumu
             (lepton=="E" && (event->mc.TMEME == 2 || event->mc.TMEME == 10101 || event->mc.TMEME == 20200)) || //ee
	      datasets[d].Name().find("DY") != string::npos || datasets[d].Name() == "Zjets"
             ){

	  if(leptons.size() == 2 && sel.passTriggerSelection ( &datasets[d], trigger)){
           for(unsigned int i=0; i<leptons.size(); i++){
	     if(leptons[i].RelIso03PF() < looseIso) loose_iso_counter_eff += 1;
	     if(leptons[i].RelIso03PF() < tightIso) tight_iso_counter_eff += 1;
	   }
	  }

	  if(loose_iso_counter_eff == 2){
	             // weight in case several datasets are used
                     IsoLooseForEfficiencyNJets->Fill(jets_Loose.size(),weight);
                     IsoLooseForEfficiencyEta->Fill(leptons[0].p4.Eta(),weight);
                     IsoLooseForEfficiencyPT1->Fill(leptons[0].p4.Perp(),weight);
                     IsoLooseForEfficiencyPT2->Fill(leptons[1].p4.Perp(),weight);
		     IsoLooseForEfficiencyNVTX->Fill(selVertices.size(),weight);
		     IsoLooseForEfficiencyNJetsvsPT1->Fill(leptons[0].p4.Perp(), jets_Loose.size(),weight);
		     IsoLooseForEfficiencyNJetsvsPT2->Fill(leptons[1].p4.Perp(), jets_Loose.size(),weight);
		     IsoLooseForEfficiencyPT1vsPT2->Fill(leptons[1].p4.Perp(), leptons[0].p4.Perp(),weight);
                     Unweighted_IsoLooseForEfficiencyNJets->Fill(jets_Loose.size());
                     Unweighted_IsoLooseForEfficiencyEta->Fill(leptons[0].p4.Eta());
                     Unweighted_IsoLooseForEfficiencyPT1->Fill(leptons[0].p4.Perp());
                     Unweighted_IsoLooseForEfficiencyPT2->Fill(leptons[1].p4.Perp());
                     Unweighted_IsoLooseForEfficiencyNVTX->Fill(selVertices.size());
		     Unweighted_IsoLooseForEfficiencyNJetsvsPT1->Fill(leptons[0].p4.Perp(), jets_Loose.size());
		     Unweighted_IsoLooseForEfficiencyNJetsvsPT2->Fill(leptons[1].p4.Perp(), jets_Loose.size());
		     Unweighted_IsoLooseForEfficiencyPT1vsPT2->Fill(leptons[1].p4.Perp(), leptons[0].p4.Perp());
	  }
	  if(tight_iso_counter_eff == 2){
                     IsoTightForEfficiencyNJets->Fill(jets_Tight.size(),weight);
                     IsoTightForEfficiencyEta->Fill(leptons[0].p4.Eta(),weight);
                     IsoTightForEfficiencyPT1->Fill(leptons[0].p4.Perp(),weight);
                     IsoTightForEfficiencyPT2->Fill(leptons[1].p4.Perp(),weight);
		     IsoTightForEfficiencyNVTX->Fill(selVertices.size(),weight);
		     IsoTightForEfficiencyNJetsvsPT1->Fill(leptons[0].p4.Perp(), jets_Tight.size(),weight);
		     IsoTightForEfficiencyNJetsvsPT2->Fill(leptons[1].p4.Perp(), jets_Tight.size(),weight);
		     IsoTightForEfficiencyPT1vsPT2->Fill(leptons[1].p4.Perp(), leptons[0].p4.Perp(),weight);
                     Unweighted_IsoTightForEfficiencyNJets->Fill(jets_Tight.size());
                     Unweighted_IsoTightForEfficiencyEta->Fill(leptons[0].p4.Eta());
                     Unweighted_IsoTightForEfficiencyPT1->Fill(leptons[0].p4.Perp());
                     Unweighted_IsoTightForEfficiencyPT2->Fill(leptons[1].p4.Perp());
                     Unweighted_IsoTightForEfficiencyNVTX->Fill(selVertices.size());
		     Unweighted_IsoTightForEfficiencyNJetsvsPT1->Fill(leptons[0].p4.Perp(), jets_Tight.size());
		     Unweighted_IsoTightForEfficiencyNJetsvsPT2->Fill(leptons[1].p4.Perp(), jets_Tight.size());
		     Unweighted_IsoTightForEfficiencyPT1vsPT2->Fill(leptons[1].p4.Perp(), leptons[0].p4.Perp());
	  }
	  }

          // For fake rate
	  if(
             (lepton=="Mu" && (event->mc.TMEME == 10 || // l(=mu)jets
              event->mc.TMEME == 11 || event->mc.TMEME == 11001 || event->mc.TMEME == 10110 || event->mc.TMEME == 21100 || // emu
              event->mc.TMEME == 10010 || // l(=mu)tau(->had)
              event->mc.TMEME == 11000)) || // tau(->mu)jets
             (lepton=="E" && (event->mc.TMEME == 1 || // l(=e)jets
              event->mc.TMEME == 11 || event->mc.TMEME == 11001 || event->mc.TMEME == 10110 || event->mc.TMEME == 21100 || // emu
              event->mc.TMEME == 10001 || // l(=e)tau(->had)
              event->mc.TMEME == 10100)) // tau(->e)jets
	     ){
	  if(leptons.size() == 2 && sel.passTriggerSelection ( &datasets[d], trigger)){
           for(unsigned int i=0; i<leptons.size(); i++){
	     if(leptons[i].RelIso03PF() < looseIso) loose_iso_counter_fake += 1;
	     if(leptons[i].RelIso03PF() < tightIso) tight_iso_counter_fake += 1;
	   }
	  }

	  if(loose_iso_counter_fake == 2){
                     IsoLooseForFakeRateNJets->Fill(jets_Loose.size(),weight);
                     IsoLooseForFakeRateEta->Fill(leptons[0].p4.Eta(),weight);
                     IsoLooseForFakeRatePT1->Fill(leptons[0].p4.Perp(),weight);
                     IsoLooseForFakeRatePT2->Fill(leptons[1].p4.Perp(),weight);
		     IsoLooseForFakeRateNVTX->Fill(selVertices.size(),weight);
		     IsoLooseForFakeRateNJetsvsPT1->Fill(leptons[0].p4.Perp(), jets_Loose.size(),weight);
		     IsoLooseForFakeRateNJetsvsPT2->Fill(leptons[1].p4.Perp(), jets_Loose.size(),weight);
		     IsoLooseForFakeRatePT1vsPT2->Fill(leptons[1].p4.Perp(), leptons[0].p4.Perp(),weight);
                     Unweighted_IsoLooseForFakeRateNJets->Fill(jets_Loose.size());
                     Unweighted_IsoLooseForFakeRateEta->Fill(leptons[0].p4.Eta());
                     Unweighted_IsoLooseForFakeRatePT1->Fill(leptons[0].p4.Perp());
                     Unweighted_IsoLooseForFakeRatePT2->Fill(leptons[1].p4.Perp());
                     Unweighted_IsoLooseForFakeRateNVTX->Fill(selVertices.size());
		     Unweighted_IsoLooseForFakeRateNJetsvsPT1->Fill(leptons[0].p4.Perp(), jets_Loose.size());
		     Unweighted_IsoLooseForFakeRateNJetsvsPT2->Fill(leptons[1].p4.Perp(), jets_Loose.size());
		     Unweighted_IsoLooseForFakeRatePT1vsPT2->Fill(leptons[1].p4.Perp(), leptons[0].p4.Perp());
	  }
	  if(tight_iso_counter_fake == 2){
                     IsoTightForFakeRateNJets->Fill(jets_Tight.size(),weight);
                     IsoTightForFakeRateEta->Fill(leptons[0].p4.Eta(),weight);
                     IsoTightForFakeRatePT1->Fill(leptons[0].p4.Perp(),weight);
                     IsoTightForFakeRatePT2->Fill(leptons[1].p4.Perp(),weight);
		     IsoTightForFakeRateNVTX->Fill(selVertices.size(),weight);
		     IsoTightForFakeRateNJetsvsPT1->Fill(leptons[0].p4.Perp(), jets_Tight.size(),weight);
		     IsoTightForFakeRateNJetsvsPT2->Fill(leptons[1].p4.Perp(), jets_Tight.size(),weight);
		     IsoTightForFakeRatePT1vsPT2->Fill(leptons[1].p4.Perp(), leptons[0].p4.Perp(),weight);
                     Unweighted_IsoTightForFakeRateNJets->Fill(jets_Tight.size());
                     Unweighted_IsoTightForFakeRateEta->Fill(leptons[0].p4.Eta());
                     Unweighted_IsoTightForFakeRatePT1->Fill(leptons[0].p4.Perp());
                     Unweighted_IsoTightForFakeRatePT2->Fill(leptons[1].p4.Perp());
                     Unweighted_IsoTightForFakeRateNVTX->Fill(selVertices.size());
		     Unweighted_IsoTightForFakeRateNJetsvsPT1->Fill(leptons[0].p4.Perp(), jets_Tight.size());
		     Unweighted_IsoTightForFakeRateNJetsvsPT2->Fill(leptons[1].p4.Perp(), jets_Tight.size());
		     Unweighted_IsoTightForFakeRatePT1vsPT2->Fill(leptons[1].p4.Perp(), leptons[0].p4.Perp());
	  }
	  }
       } // End of mc specific part





       //--------------------//
       // DATA specific part //
       //--------------------//
 
       if(!mc)
       {
	 bool ainvmasscutW = true;
	 if(leptons.size() == 1){
            float phi_l = leptons[0].p4.Phi();
            float pT_l = leptons[0].p4.Perp();
            float phi_met = met.p2.Phi();
            float deltaphi = phi_l - phi_met;
            if(deltaphi < -1*M_PI ) deltaphi += 2*M_PI;
            if(deltaphi > M_PI    ) deltaphi -= 2*M_PI;
            float pT_met = met.p2.Mod();
            float mT = sqrt(2*pT_l*pT_met*(1-cos(deltaphi)));
	    //          if(mT < 80 && mT > 65) ainvmasscutW = false;
            if(mT > 20) ainvmasscutW = false;
	 }


          // Lot enrichi en QCD , legere contamination en W -> prend comme syst la diff avec MC
          if(leptons.size() == 1 && met.p2.Mod() < 20 && jets.size() >= 1 && ainvmasscutW && datasets[d].Name() == "Jet"){   	    
               if (leptons[0].RelIso03PF() < looseIso){ 
                    IsoLooseForFakeRateNJets->Fill(jets_Loose.size());
                    IsoLooseForFakeRateEta->Fill(leptons[0].p4.Eta());
                    IsoLooseForFakeRatePT1->Fill(leptons[0].p4.Perp());
                    IsoLooseForFakeRateNVTX->Fill(selVertices.size());
               }
   	       if (leptons[0].RelIso03PF() < tightIso){ 
                    IsoTightForFakeRateNJets->Fill(jets_Tight.size());
                    IsoTightForFakeRateEta->Fill(leptons[0].p4.Eta());
                    IsoTightForFakeRatePT1->Fill(leptons[0].p4.Perp());
                     IsoTightForFakeRateNVTX->Fill(selVertices.size());
               }
          } // End of QCD events


         // Z events for efficiencies 
	 if(leptons.size() == 2 && datasets[d].Name() == "DataMu"){
          TLorentzVector Dilepton;
          Dilepton = leptons[0].p4 + leptons[1].p4;
          bool ainvmasscut = (Dilepton.M() < 106 && Dilepton.M() > 76);
           
	  int pair_charge = leptons[0].charge * leptons[1].charge;
          for(int ilep=0; ilep<2; ilep++)
          if(leptons[ilep].charge == -1 && leptons[ilep].RelIso03PF() < 0.2){ // Tag lepton is - charged lepton only
            if(pair_charge == -1 && ainvmasscut){ // os pair
              if (leptons[ilep].RelIso03PF() < looseIso){ 
                    IsoLooseForEfficiencyNJets->Fill(jets_Loose.size());
                    IsoLooseForEfficiencyEta->Fill(leptons[ilep].p4.Eta());
                    IsoLooseForEfficiencyPT1->Fill(leptons[ilep].p4.Perp());
                    IsoLooseForEfficiencyNVTX->Fill(selVertices.size());
		    IsoLooseForEfficiencyNJetsvsPT1->Fill(leptons[0].p4.Perp(), jets_Loose.size(),weight);
              }
	      if (leptons[ilep].RelIso03PF() < tightIso){ 
                    IsoTightForEfficiencyNJets->Fill(jets_Tight.size());
                    IsoTightForEfficiencyEta->Fill(leptons[ilep].p4.Eta());
                    IsoTightForEfficiencyPT1->Fill(leptons[ilep].p4.Perp());
                    IsoTightForEfficiencyNVTX->Fill(selVertices.size());
		    IsoTightForEfficiencyNJetsvsPT1->Fill(leptons[ilep].p4.Perp(), jets_Tight.size(),weight);
              }
            }
          }

	 } // End of Z selection
	 
	 
       } // End of data specific part
       
       
    } // end of loop over evts




  } // end of loop over the datasets 
  cout<<"#########################"<<endl;
  cout<<" Loop over the datasets  "<<endl;
  cout<<"#########################"<<endl;

  ////////////////////////////
  //  Computation after loops
  ////////////////////////////

  string leptonname = "";
  if(lepton=="Mu") leptonname = "muon";
  if(lepton=="E") leptonname = "electron";
  string datatypename = "";
  if(mc) datatypename = "MC";
  else datatypename = "DATA";
  
  TFile* file_MM = new TFile(("MatrixMethod_Efficiency_"+datatypename+"_"+lepton+".root").c_str(),"RECREATE");
  file_MM->cd();


  if(mc) {
    ComputeEfficiencyForMC (IsoTightForEfficiencyNJets, IsoLooseForEfficiencyNJets,
      Unweighted_IsoTightForEfficiencyNJets, Unweighted_IsoLooseForEfficiencyNJets, NJets_N_bins, 
      SignalEfficiencyNJets);

    ComputeEfficiencyForMC (IsoTightForEfficiencyEta, IsoLooseForEfficiencyEta,
      Unweighted_IsoTightForEfficiencyEta, Unweighted_IsoLooseForEfficiencyEta, Eta_N_bins, 
      SignalEfficiencyEta);

    ComputeEfficiencyForMC (IsoTightForEfficiencyPT1, IsoLooseForEfficiencyPT1,
      Unweighted_IsoTightForEfficiencyPT1, Unweighted_IsoLooseForEfficiencyPT1, PT_N_bins, 
      SignalEfficiencyPT1);

    ComputeEfficiencyForMC (IsoTightForEfficiencyPT2, IsoLooseForEfficiencyPT2,
      Unweighted_IsoTightForEfficiencyPT2, Unweighted_IsoLooseForEfficiencyPT2, PT_N_bins, 
      SignalEfficiencyPT2);

    ComputeEfficiencyForMC (IsoTightForEfficiencyNVTX, IsoLooseForEfficiencyNVTX,
      Unweighted_IsoTightForEfficiencyNVTX, Unweighted_IsoLooseForEfficiencyNVTX, NVTX_N_bins, 
      SignalEfficiencyNVTX);

    ComputeEfficiencyForMC2D (IsoTightForEfficiencyNJetsvsPT1, IsoLooseForEfficiencyNJetsvsPT1,
      Unweighted_IsoTightForEfficiencyNJetsvsPT1, Unweighted_IsoLooseForEfficiencyNJetsvsPT1, PT_N_bins, NJets_N_bins,
      SignalEfficiencyNJetsvsPT1);

    ComputeEfficiencyForMC2D (IsoTightForEfficiencyNJetsvsPT2, IsoLooseForEfficiencyNJetsvsPT2,
      Unweighted_IsoTightForEfficiencyNJetsvsPT2, Unweighted_IsoLooseForEfficiencyNJetsvsPT2, PT_N_bins, NJets_N_bins,
      SignalEfficiencyNJetsvsPT2);

    ComputeEfficiencyForMC2D (IsoTightForEfficiencyPT1vsPT2, IsoLooseForEfficiencyPT1vsPT2,
      Unweighted_IsoTightForEfficiencyPT1vsPT2, Unweighted_IsoLooseForEfficiencyPT1vsPT2, PT_N_bins, PT_N_bins,
      SignalEfficiencyPT1vsPT2);

    for(unsigned int bin_index=1; bin_index<PT_N_bins+1; bin_index++){
      SignalEfficiencyPTfromPT1vsPT2->SetBinContent(bin_index, SignalEfficiencyPT1vsPT2->GetBinContent(bin_index, bin_index));
      SignalEfficiencyPTfromPT1vsPT2->SetBinError(bin_index, SignalEfficiencyPT1vsPT2->GetBinError(bin_index, bin_index));
    }

    ComputeFakeRateForMC (IsoTightForEfficiencyNJets, IsoLooseForEfficiencyNJets, 
      IsoTightForFakeRateNJets, IsoLooseForFakeRateNJets, 
      Unweighted_IsoTightForEfficiencyNJets, Unweighted_IsoLooseForEfficiencyNJets, 
      Unweighted_IsoTightForFakeRateNJets, Unweighted_IsoLooseForFakeRateNJets, NJets_N_bins, 
      FakeRateNJets);

    ComputeFakeRateForMC (IsoTightForEfficiencyEta, IsoLooseForEfficiencyEta, 
      IsoTightForFakeRateEta, IsoLooseForFakeRateEta, 
      Unweighted_IsoTightForEfficiencyEta, Unweighted_IsoLooseForEfficiencyEta, 
      Unweighted_IsoTightForFakeRateEta, Unweighted_IsoLooseForFakeRateEta, Eta_N_bins, 
      FakeRateEta);

    ComputeFakeRateForMC (IsoTightForEfficiencyPT1, IsoLooseForEfficiencyPT1, 
      IsoTightForFakeRatePT1, IsoLooseForFakeRatePT1, 
      Unweighted_IsoTightForEfficiencyPT1, Unweighted_IsoLooseForEfficiencyPT1, 
      Unweighted_IsoTightForFakeRatePT1, Unweighted_IsoLooseForFakeRatePT1, PT_N_bins, 
      FakeRatePT1);

    ComputeFakeRateForMC (IsoTightForEfficiencyPT2, IsoLooseForEfficiencyPT2, 
      IsoTightForFakeRatePT2, IsoLooseForFakeRatePT2, 
      Unweighted_IsoTightForEfficiencyPT2, Unweighted_IsoLooseForEfficiencyPT2, 
      Unweighted_IsoTightForFakeRatePT2, Unweighted_IsoLooseForFakeRatePT2, PT_N_bins, 
      FakeRatePT2);

    ComputeFakeRateForMC (IsoTightForEfficiencyNVTX, IsoLooseForEfficiencyNVTX, 
      IsoTightForFakeRateNVTX, IsoLooseForFakeRateNVTX, 
      Unweighted_IsoTightForEfficiencyNVTX, Unweighted_IsoLooseForEfficiencyNVTX, 
      Unweighted_IsoTightForFakeRateNVTX, Unweighted_IsoLooseForFakeRateNVTX, NVTX_N_bins, 
      FakeRateNVTX);
      
    ComputeFakeRateForMC2D (IsoTightForEfficiencyNJetsvsPT1, IsoLooseForEfficiencyNJetsvsPT1, 
      IsoTightForFakeRateNJetsvsPT1, IsoLooseForFakeRateNJetsvsPT1, 
      Unweighted_IsoTightForEfficiencyNJetsvsPT1, Unweighted_IsoLooseForEfficiencyNJetsvsPT1, 
      Unweighted_IsoTightForFakeRateNJetsvsPT1, Unweighted_IsoLooseForFakeRateNJetsvsPT1, PT_N_bins, NJets_N_bins,
      FakeRateNJetsvsPT1);
      
    ComputeFakeRateForMC2D (IsoTightForEfficiencyNJetsvsPT2, IsoLooseForEfficiencyNJetsvsPT2, 
      IsoTightForFakeRateNJetsvsPT2, IsoLooseForFakeRateNJetsvsPT2, 
      Unweighted_IsoTightForEfficiencyNJetsvsPT2, Unweighted_IsoLooseForEfficiencyNJetsvsPT2, 
      Unweighted_IsoTightForFakeRateNJetsvsPT2, Unweighted_IsoLooseForFakeRateNJetsvsPT2, PT_N_bins, NJets_N_bins,
      FakeRateNJetsvsPT2);
      
    ComputeFakeRateForMC2D (IsoTightForEfficiencyPT1vsPT2, IsoLooseForEfficiencyPT1vsPT2, 
      IsoTightForFakeRatePT1vsPT2, IsoLooseForFakeRatePT1vsPT2, 
      Unweighted_IsoTightForEfficiencyPT1vsPT2, Unweighted_IsoLooseForEfficiencyPT1vsPT2, 
      Unweighted_IsoTightForFakeRatePT1vsPT2, Unweighted_IsoLooseForFakeRatePT1vsPT2, PT_N_bins, PT_N_bins,
      FakeRatePT1vsPT2);
    
    for(unsigned int bin_index=1; bin_index<PT_N_bins+1; bin_index++){
      FakeRatePTfromPT1vsPT2->SetBinContent(bin_index, FakeRatePT1vsPT2->GetBinContent(bin_index, bin_index));
      FakeRatePTfromPT1vsPT2->SetBinError(bin_index, FakeRatePT1vsPT2->GetBinError(bin_index, bin_index));
    }
    
  }
  else { // DATA
  
    ComputeEfficiencyForData (IsoTightForEfficiencyNJets, IsoLooseForEfficiencyNJets, NJets_N_bins, SignalEfficiencyNJets);
    ComputeEfficiencyForData (IsoTightForEfficiencyEta, IsoLooseForEfficiencyEta, Eta_N_bins, SignalEfficiencyEta);
    ComputeEfficiencyForData (IsoTightForEfficiencyPT1, IsoLooseForEfficiencyPT1, PT_N_bins, SignalEfficiencyPT1);
    ComputeEfficiencyForData (IsoTightForEfficiencyNVTX, IsoLooseForEfficiencyNVTX, NVTX_N_bins, SignalEfficiencyNVTX);

    ComputeFakeRateForData (IsoTightForFakeRateNJets, IsoLooseForFakeRateNJets, NJets_N_bins, FakeRateNJets);
    ComputeFakeRateForData (IsoTightForFakeRateEta, IsoLooseForFakeRateEta, Eta_N_bins, FakeRateEta);
    ComputeFakeRateForData (IsoTightForFakeRatePT1, IsoLooseForFakeRatePT1, PT_N_bins, FakeRatePT1);
    ComputeFakeRateForData (IsoTightForFakeRateNVTX, IsoLooseForFakeRateNVTX, NVTX_N_bins, FakeRateNVTX);
  }


  cout << "\\begin{table}" << std::endl;
  cout << "\\begin{center}" << std::endl;
  cout << "\\begin{tabular}{|";
  for(unsigned int bin_index=0; bin_index<NJets_N_bins; bin_index++){
   cout << " c |";
  }
  cout << "}" << std::endl;
  cout << "\\hline" << std::endl;
  for(unsigned int bin_index=0; bin_index<NJets_N_bins; bin_index++){
   cout << " & Njets = " << bin_index;
  }
  cout << "\\\\ \\hline" << std::endl;
  cout << "$\\epsilon_{s}$ "<<leptonname<<" case";
  for(unsigned int bin_index=0; bin_index<NJets_N_bins; bin_index++){
   cout << " & " << SignalEfficiencyNJets->GetBinContent(bin_index) << "$\\pm$" << SignalEfficiencyNJets->GetBinError(bin_index);
  }
  cout << "\\\\ \\hline" << std::endl;
  cout << "\\end{tabular}" << std::endl;
  cout << "\\caption{FIX ME} \\label{FIX ME}" << std::endl;
  cout << "\\end{center}" << std::endl;
  cout << "\\end{table}" << std::endl;



  cout << "\\begin{table}" << std::endl;
  cout << "\\begin{center}" << std::endl;
  cout << "\\begin{tabular}{|";
  for(unsigned int bin_index=0; bin_index<Eta_N_bins; bin_index++){
   cout << " c |";
  }
  cout << "}" << std::endl;
  cout << "\\hline" << std::endl;
  for(unsigned int bin_index=0; bin_index<Eta_N_bins; bin_index++){
   cout << " & $\\eta$ = " << ((Eta_Max-Eta_Min)/Eta_N_bins)*(bin_index+0.5-Eta_N_bins*0.5);
  }
  cout << "\\\\ \\hline" << std::endl;
  cout << "$\\epsilon_{s}$ "<<leptonname<<" case";
  for(unsigned int bin_index=0; bin_index<Eta_N_bins; bin_index++){
   cout << " & " << SignalEfficiencyEta->GetBinContent(bin_index) << "$\\pm$" << SignalEfficiencyEta->GetBinError(bin_index);
  }
  cout << "\\\\ \\hline" << std::endl;
  cout << "\\end{tabular}" << std::endl;
  cout << "\\caption{FIX ME} \\label{FIX ME}" << std::endl;
  cout << "\\end{center}" << std::endl;
  cout << "\\end{table}" << std::endl;



  cout << "\\begin{table}" << std::endl;
  cout << "\\begin{center}" << std::endl;
  cout << "\\begin{tabular}{|";
  for(unsigned int bin_index=0; bin_index<PT_N_bins; bin_index++){
   cout << " c |";
  }
  cout << "}" << std::endl;
  cout << "\\hline" << std::endl;
  for(unsigned int bin_index=0; bin_index<PT_N_bins; bin_index++){
   cout << " & $p_T$ = " <<   ((PT_Max-PT_Min)/PT_N_bins)*(bin_index+0.5);
  }
  cout << "\\\\ \\hline" << std::endl;
  cout << "$\\epsilon_{s}$ "<<leptonname<<" case";
  for(unsigned int bin_index=0; bin_index<PT_N_bins; bin_index++){
   cout << " & " << SignalEfficiencyPT1->GetBinContent(bin_index) << "$\\pm$" << SignalEfficiencyPT1->GetBinError(bin_index);
  }
  cout << "\\\\ \\hline" << std::endl;
  cout << "\\end{tabular}" << std::endl;
  cout << "\\caption{FIX ME} \\label{FIX ME}" << std::endl;
  cout << "\\end{center}" << std::endl;
  cout << "\\end{table}" << std::endl;




  cout << "\\begin{table}" << std::endl;
  cout << "\\begin{center}" << std::endl;
  cout << "\\begin{tabular}{|";
  for(unsigned int bin_index=0; bin_index<NJets_N_bins; bin_index++){
   cout << " c |";
  }
  cout << "}" << std::endl;
  cout << "\\hline" << std::endl;
  for(unsigned int bin_index=0; bin_index<NJets_N_bins; bin_index++){
   cout << " & Njets = " << bin_index;
  }
  cout << "\\\\ \\hline" << std::endl;
  cout << "$\\epsilon_{f}$ "<<leptonname<<" case";
  for(unsigned int bin_index=0; bin_index<NJets_N_bins; bin_index++){
   cout << " & " << FakeRateNJets->GetBinContent(bin_index) << "$\\pm$" << FakeRateNJets->GetBinError(bin_index);
  }
  cout << "\\\\ \\hline" << std::endl;
  cout << "\\end{tabular}" << std::endl;
  cout << "\\caption{FIX ME} \\label{FIX ME}" << std::endl;
  cout << "\\end{center}" << std::endl;
  cout << "\\end{table}" << std::endl;



  cout << "\\begin{table}" << std::endl;
  cout << "\\begin{center}" << std::endl;
  cout << "\\begin{tabular}{|";
  for(unsigned int bin_index=0; bin_index<Eta_N_bins; bin_index++){
   cout << " c |";
  }
  cout << "}" << std::endl;
  cout << "\\hline" << std::endl;
  for(unsigned int bin_index=0; bin_index<Eta_N_bins; bin_index++){
   cout << " & $\\eta$ = " << ((Eta_Max-Eta_Min)/Eta_N_bins)*(bin_index+0.5-Eta_N_bins*0.5);
  }
  cout << "\\\\ \\hline" << std::endl;
  cout << "$\\epsilon_{f}$ "<<leptonname<<" case";
  for(unsigned int bin_index=0; bin_index<Eta_N_bins; bin_index++){
   cout << " & " << FakeRateEta->GetBinContent(bin_index) << "$\\pm$" << FakeRateEta->GetBinError(bin_index);
  }
  cout << "\\\\ \\hline" << std::endl;
  cout << "\\end{tabular}" << std::endl;
  cout << "\\caption{FIX ME} \\label{FIX ME}" << std::endl;
  cout << "\\end{center}" << std::endl;
  cout << "\\end{table}" << std::endl;



  cout << "\\begin{table}" << std::endl;
  cout << "\\begin{center}" << std::endl;
  cout << "\\begin{tabular}{|";
  for(unsigned int bin_index=0; bin_index<PT_N_bins; bin_index++){
   cout << " c |";
  }
  cout << "}" << std::endl;
  cout << "\\hline" << std::endl;
  for(unsigned int bin_index=0; bin_index<PT_N_bins; bin_index++){
   cout << " & $p_T$ = " << ((PT_Max-PT_Min)/PT_N_bins)*(bin_index+0.5);
  }
  cout << "\\\\ \\hline" << std::endl;
  cout << "$\\epsilon_{f}$ "<<leptonname<<" case";
  for(unsigned int bin_index=0; bin_index<PT_N_bins; bin_index++){
   cout << " & " << FakeRatePT1->GetBinContent(bin_index) << "$\\pm$" << FakeRatePT1->GetBinError(bin_index);
  }
  cout << "\\\\ \\hline" << std::endl;
  cout << "\\end{tabular}" << std::endl;
  cout << "\\caption{FIX ME} \\label{FIX ME}" << std::endl;
  cout << "\\end{center}" << std::endl;
  cout << "\\end{table}" << std::endl;


  DiLepInvMass->Write();
  delete DiLepInvMass;
  Iso_Lep->Write();
  delete Iso_Lep;
  PT_Leading_Lep->Write();
  delete PT_Leading_Lep;
  PT_Inclusive_Lep->Write();
  delete PT_Inclusive_Lep;
  PT_Leading_Jet->Write();
  delete PT_Leading_Jet;
  PT_Inclusive_Jet->Write();
  delete PT_Inclusive_Jet;
  NJets->Write();
  delete NJets;
  NVTX->Write();
  delete NVTX;
 
  SignalEfficiencyNJets->Write();
  delete SignalEfficiencyNJets;
  FakeRateNJets->Write();
  delete FakeRateNJets;
  SignalEfficiencyEta->Write();
  delete SignalEfficiencyEta;
  FakeRateEta->Write();
  delete FakeRateEta;
  SignalEfficiencyPT1->Write();
  delete SignalEfficiencyPT1;
  FakeRatePT1->Write();
  delete FakeRatePT1;
  SignalEfficiencyPT2->Write();
  delete SignalEfficiencyPT2;
  FakeRatePT2->Write();
  delete FakeRatePT2;
  SignalEfficiencyNVTX->Write();
  delete SignalEfficiencyNVTX;
  FakeRateNVTX->Write();
  delete FakeRateNVTX;
  SignalEfficiencyNJetsvsPT1->Write();
  delete SignalEfficiencyNJetsvsPT1;
  FakeRateNJetsvsPT1->Write();
  delete FakeRateNJetsvsPT1;
  SignalEfficiencyNJetsvsPT2->Write();
  delete SignalEfficiencyNJetsvsPT2;
  FakeRateNJetsvsPT2->Write();
  delete FakeRateNJetsvsPT2;
  SignalEfficiencyPT1vsPT2->Write();
  delete SignalEfficiencyPT1vsPT2;
  FakeRatePT1vsPT2->Write();
  delete FakeRatePT1vsPT2;
  SignalEfficiencyPTfromPT1vsPT2->Write();
  delete SignalEfficiencyPTfromPT1vsPT2;
  FakeRatePTfromPT1vsPT2->Write();
  delete FakeRatePTfromPT1vsPT2;
 

  file_MM->Write();
  file_MM->Close();
  delete file_MM;


  cout<<"#########################"<<endl;
  cout<<"    End of the program   "<<endl;
  cout<<"#########################"<<endl;

  return;
}


int main ()
{
  Efficiency_MC_Estimation("../../config/MatrixMethod_MC.xml", true, "Mu");
  //Efficiency_MC_Estimation("E");  

  return (0);
}
