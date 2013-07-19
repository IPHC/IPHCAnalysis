#include <iomanip>
#include <iostream>
#include <time.h>
#include "NTFormat/interface/NTEvent.h"

//NTupleAnalysis classes
#include "Selection/interface/TTbarMetSelection.h"
#include "Selection/interface/SelectionTable.h"
#include "Tools/interface/Dataset.h"
#include "Tools/interface/AnalysisEnvironmentLoader.h"
#include "Plots/interface/TTbarMetHistoManager.h"
#include "Tools/interface/PUWeighting.h"
#include "Tools/interface/LumiReweightingStandAlone.h"
#include "Tools/interface/JetCorrector.h"



#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

using namespace IPHCTree;
using namespace std;

double pi=acos(-1.);

struct SortByPt
{
  bool operator()( TLorentzVector j1, TLorentzVector j2 ) const
  {
    return j1.Pt() > j2.Pt() ;
  }
};



int main (int argc, char *argv[])
{
 
  cout << "#########################" << endl;
  cout << "Beginning of the program" << endl;
  cout << "#########################" << endl;

  //////////////////////
  //Global variables
  //////////////////////
  vector < Dataset > datasets;
  TTbarMetSelection sel;
/*
  cout << sel.GetbtagAlgo() << " " << sel.GetbtagDiscriCut() << " " << sel.GetNofBtagJetsCut() << endl;
  vector<string> listCaro0=sel.GetCutList ();
  cout << listCaro0.size() << endl;
  for (int ii=0; ii<listCaro0.size(); ii++) {
    cout << listCaro0[ii] << endl;
  }
*/
  float Luminosity = 0;
  float LumiError = 0.;
  // 0: MC - 1: Data - 2 Data & MC
  int DataType = 0;
  int verbosity = -1;

//  JetCorrector JEC_L2L3Residuals;


  reweight::LumiReWeighting *LumiWeights;


  //////////////////////

  //////////////////////
  // Initialisation
  //////////////////////
  string xmlFileName;
  cout<<"argc "<<argc<<" "<<argv[0]<<endl;
  if (argc>1 ) xmlFileName = string(argv[1]);
  else xmlFileName = string ("../../config/TTbarMETAnalysis.xml");
  AnalysisEnvironmentLoader anaEL (xmlFileName);
  anaEL.LoadSamples (datasets);	// now the list of datasets written in the xml file is known
  anaEL.LoadSelection (sel);	// now the parameters for the selection are given to the selection // no specific TTbarMET parameters
  anaEL.LoadGeneralInfo (DataType, Luminosity, LumiError, verbosity);
/*  // pas pour le moment! car pas de selection specifique a TTbarMET
  int flagOriginal=sel.GetFlagb();
  int methodOriginal=sel.GetMethodb();
  int systOriginal= sel.GetSystb();
  std::cout << " For btag : flag " << flagOriginal << ", method " << methodOriginal << ", syst " << systOriginal << std::endl;
*/
  IPHCTree::NTEvent * event = 0;
  //Selection table

/*
  vector<string> listCaro=sel.GetCutList ();
  cout << listCaro.size() << endl;
  for (int ii=0; ii<listCaro.size(); ii++) {
    cout << listCaro[ii] << endl;
  }
*/
  SelectionTable selTable_e (sel.GetCutList (), datasets, string ("e"));
  SelectionTable selTable_mu (sel.GetCutList (), datasets, string ("mu"));
  SelectionTable selTable_l (sel.GetCutList (), datasets, string ("lepton"));
  
  //Book keeping of standard histos
  bool doHistoManager = true;
  TTbarMetHistoManager histoManager;
  if(doHistoManager){
  	histoManager.LoadDatasets (datasets);
  	histoManager.LoadSelectionSteps (sel.GetCutList ());
  	histoManager.LoadChannels (sel.GetChannelList ());
  	histoManager.CreateHistos ();
  }


    // JEC from JL's code
//    JEC_L2L3Residuals.LoadCorrections();


   //////////////////////

  
  //////////////////////


  cout << "The verbosity mode is " << verbosity << endl;
  cout << "The luminosity is equal to " << Luminosity << endl;
  cout << "The DataType is ";
  switch (DataType) {
  case 0:
    cout << "MC" << endl;
    break;
  case 1:
    cout << "Data" << endl;
    break;
  case 2:
    cout << "Data & MC" << endl;
    break;
  default:
    cout << " unknown" << endl;
    break;
  }
  //////////////////////



  //////////////////////
  //LOOP OVER THE DATASETS
  //////////////////////
  if (verbosity > 0) {
    cout << "#########################" << endl;
    cout << " Loop over the datasets  " << endl;
    cout << "#########################" << endl;
  }

  for (unsigned int d = 0; d < datasets.size (); d++) 
  {

    if(verbosity>2) cout<<"TTbarMET> Dataset: "<<datasets[d].Name()<<endl;
    datasets[d].eventTree ()->SetBranchAddress ("NTEvent", &event);

    unsigned int nEvents = static_cast<unsigned int>(datasets[d].eventTree ()->GetEntries ());
    cout << "TTbarMET> NEvents = " << nEvents << endl;
    cout << "TTbarMET> NEvents to run over = "<<datasets[d].NofEvtsToRunOver()<<endl;


   // PU from JL's code
   if (datasets[d].isData() == false) {
//   if(datasets[d].Name() == "DY1" || datasets[d].Name() == "signal" ){
//      string datafile = "/opt/sbg/data/data1/cms/ccollard/CMSSW/CMSSW_4_2_8_patch7/src/NTuple/NTupleAnalysis/macros/data/PUdata.root";
//      string mcfile   = "/opt/sbg/data/data1/cms/ccollard/CMSSW/CMSSW_4_2_8_patch7/src/NTuple/NTupleAnalysis/macros/data/PU3DMC_Fall11.root";
//
//      LumiWeights    = new reweight::LumiReWeighting(mcfile, datafile, "histoMCPU", "pileup" );
//    }
//    else{
      string datafile = "/opt/sbg/data/data1/cms/ccollard/CMSSW/CMSSW_4_2_8_patch7/src/NTuple/NTupleAnalysis/macros/data/default73.5mb.root";
      string mcfile   = "/opt/sbg/data/data1/cms/ccollard/CMSSW/CMSSW_4_2_8_patch7/src/NTuple/NTupleAnalysis/macros/data/PU3DMC.root";

      LumiWeights    = new reweight::LumiReWeighting(mcfile, datafile, "histoMCPU", "pileup" );
      LumiWeights->weight3D_init( 1. );

//    }
    }



    //LOOP OVER THE EVENTS
    for (unsigned int ievt = 0; ievt < datasets[d].NofEvtsToRunOver(); ievt++)
    {


      datasets[d].eventTree()->GetEntry(ievt);
      IPHCTree::NTTransient::InitializeAfterReading(event);



      float weight = 1.;
      if(datasets[d].isData() == false) weight = datasets[d].NormFactor()*Luminosity; //if Data , weight = 1

      float weightpu=1.;
      if(datasets[d].isData() == false) { // MC
//        if( datasets[d].Name() == "DY1" || datasets[d].Name() == "signal")  {
//          double ave_npu = (event->pileup.before_npu+event->pileup.intime_npu+event->pileup.after_npu)/3.;
//          weightpu = LumiWeights->ITweight3BX(ave_npu);
//        }
//        else {
          weightpu = LumiWeights->weight3D(event->pileup.before_npu, event->pileup.intime_npu, event->pileup.after_npu);
//        }
//        weight *= weightpu; //if Data , weight = 1
      }
//      else { // DATA
//         JEC_L2L3Residuals.ApplyCorrections(event); // n'appliquer la correction que pour les donnees
//      }

      if (ievt==0) {
	      cout << "TTbarMET> weight : " << weight << endl;
      }
      if (verbosity > 3)
      {
        std::cout << "event " << ievt 
                  <<" - event number=" << event->general.eventNb
                  <<" - run number="   << event->general.runNb
                  << std::endl;
      }
      if (ievt % 5000 == 0)
	      cout << "TTbarMET> Progress bar : " << ievt << endl;

      //Load event for the selection
      sel.LoadEvent(event);

      // Selection
      int selLastStep = sel.doFullSelection (&(datasets[d]), string("all"), false);

      // Table
      //  1. No Selection 
      selTable_e.Fill (d, 0, weight);
      selTable_mu.Fill (d, 0, weight);
      selTable_l.Fill (d, 0, weight);
      bool trigger_mu=false;
      bool trigger_e=false;
      //  2. Trigger
      if (selLastStep>0) {
       selTable_l.Fill (d, 1, weight);
       if (selLastStep>30) {
         trigger_e=true;
         trigger_mu=true;
         selLastStep-=30;
         selTable_e.Fill (d, 1, weight);
         selTable_mu.Fill (d, 1, weight);
       }
       else if (selLastStep>20) {
         trigger_mu=true;
         selLastStep-=20;
         selTable_mu.Fill (d, 1, weight);
       }
       else if (selLastStep>10) {
         trigger_e=true;
         selLastStep-=10;
         selTable_e.Fill (d, 1, weight);
       }
      }
      //  3. Suite de la Table
      int lep = sel.GetLeptonType();
      for (unsigned int i = 2; i < (sel.GetCutList()).size() ; i++) {
        if (selLastStep >= (int) i && lep==0)  selTable_e.Fill (d, i, weight);
        if (selLastStep >= (int) i && lep==1)  selTable_mu.Fill (d, i, weight);
        if (selLastStep >= (int) i) selTable_l.Fill (d, i, weight);
      }
      

      // Histo
      int selLastStep_e=0;
      int selLastStep_mu=0;
      if (trigger_e) selLastStep_e=1;
      if (trigger_mu) selLastStep_mu=1;
      if (trigger_e && lep==0) selLastStep_e=selLastStep;
      if (trigger_mu && lep==1) selLastStep_mu=selLastStep;
      histoManager.Fill(sel, event, sel.GetMuonsForAna(), sel.GetElectronsForAna(), selLastStep_e, 0, d, weight);
      histoManager.Fill(sel, event, sel.GetMuonsForAna(), sel.GetElectronsForAna(), selLastStep_mu, 1, d, weight);
      //histoManager.Fill(sel, event, sel.GetMuonsForAna(), sel.GetElectronsForAna(), selLastStep, 2, d, weight);



    }				// end of loop over evts


  }				// end of loop over the datasets 
  cout << "#########################" << endl;
  cout << " Loop over the datasets  " << endl;
  cout << "#########################" << endl;



  ////////////////////////////
  //  Computation after loops
  ////////////////////////////

  ////////////////////////////
  //  Histograms
  ////////////////////////////

  if(doHistoManager) histoManager.Compute ();

  if (verbosity > 0) {
    cout << "#########################" << endl;
    cout << " Fill the latex tables   " << endl;
    cout << "#########################" << endl;
  }                             
  
  string ofilename;
  if(argc>3) ofilename = argv[3];
  else ofilename = string("CrossSectionTable.tex");
  ofstream ofile (ofilename.c_str());
  ofile << "\\documentclass[8pt]{article}" << endl;
  ofile << "\\usepackage{lscape}" << endl;
  ofile << "\\begin{document}" << endl;
  ofile << "\\begin{landscape}" << endl;
  
//  ofile << "\\usepackage{geometry}" << endl;
//  ofile << "\\geometry{a4paper, textwidth=19cm, textheight=23cm}" << endl;

 //Merge channels consistently
  vector < string > mergenames;
  mergenames.push_back ("W1Jet");
  mergenames.push_back ("W2Jet");
  mergenames.push_back ("W3Jet");
  mergenames.push_back ("W4Jet");
  selTable_e.MergeDatasets (mergenames, string ("W+Jets"));
  selTable_mu.MergeDatasets (mergenames, string ("W+Jets"));
  selTable_l.MergeDatasets (mergenames, string ("W+Jets"));
  mergenames.clear ();
  mergenames.push_back ("DY1");
  mergenames.push_back ("DY2");
  selTable_e.MergeDatasets (mergenames, string ("Z+Jets"));
  selTable_mu.MergeDatasets (mergenames, string ("Z+Jets"));
  selTable_l.MergeDatasets (mergenames, string ("Z+Jets"));
  mergenames.clear ();
  mergenames.push_back ("WW");
  mergenames.push_back ("WZ");
  mergenames.push_back ("ZZ");
  selTable_e.MergeDatasets (mergenames, string ("Diboson"));
  selTable_mu.MergeDatasets (mergenames, string ("Diboson"));
  selTable_l.MergeDatasets (mergenames, string ("Diboson"));
  mergenames.clear ();
  mergenames.push_back ("singleTop1");
  mergenames.push_back ("singleTop2");
  mergenames.push_back ("singleTop3");
  mergenames.push_back ("singleTop4");
  mergenames.push_back ("singleTop5");
  mergenames.push_back ("singleTop6");
  selTable_e.MergeDatasets (mergenames, string ("SingleTop"));
  selTable_mu.MergeDatasets (mergenames, string ("SingleTop"));
  selTable_l.MergeDatasets (mergenames, string ("SingleTop"));


  //Calculations
  selTable_e.TableCalculator ();
  selTable_mu.TableCalculator ();
  selTable_l.TableCalculator ();
  //Write
  selTable_e.Write (ofile);
  selTable_mu.Write (ofile);
  selTable_l.Write (ofile);
  ofile << "\\end{landscape}" << endl;
  ofile << "\\end{document}" << endl;
  system (Form("pdflatex %s", ofilename.data()));



  if (verbosity > 0) {
    cout << "#########################" << endl;
    cout << " Write output root file " << endl;
    cout << "#########################" << endl;
  }
  string orootfilename;
  if(argc>2) orootfilename = argv[2];
  else orootfilename = string("TTbarMETanalysis.root");
  TFile *fout = new TFile (orootfilename.c_str(), "RECREATE");
  if(doHistoManager) histoManager.Write (fout);
  fout->Close ();
  
  //Clear histos before deleting the TFile
  if(doHistoManager) histoManager.Clear ();

  delete fout;



  if (verbosity > 0) {
    cout << "#########################" << endl;
    cout << "    End of the program   " << endl;
    cout << "#########################" << endl;
  }

  return (0);
}


