#include <iomanip>
#include <iostream>
#include <time.h>
#include "NTFormat/interface/NTEvent.h"

//NTupleAnalysis classes
#include "Selection/interface/SSDiLeptonSelection.h"
#include "Selection/interface/SelectionTable.h"
#include "Tools/interface/Dataset.h"
#include "Tools/interface/AnalysisEnvironmentLoader.h"
#include "Plots/interface/SSDiLepAnaHistoManager.h"


#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

using namespace IPHCTree;
using namespace std;


int main (int argc, char *argv[])
{
 
  cout << "#########################" << endl;
  cout << "Beginning of the program" << endl;
  cout << "#########################" << endl;

  //////////////////////
  //Global variables
  //////////////////////
  vector < Dataset > datasets;
  SSDiLeptonSelection sel;
  float Luminosity = 0;
  float LumiError = 0.;
  // 0: MC - 1: Data - 2 Data & MC
  int DataType = 0;
  int verbosity = -1;

  int AnaStep = 6;//which defines the cuts that the events should pass to be considered as selected

  float Nobs_mumu = 0.;
  //////////////////////

  //////////////////////
  // Initialisation
  //////////////////////
  string xmlFileName;
  cout<<"argc "<<argc<<" "<<argv[0]<<endl;
  if (argc>1 ) xmlFileName = string(argv[1]);
  else xmlFileName = string ("../../config/RPVAnalysis.xml");
  AnalysisEnvironmentLoader anaEL (xmlFileName);
  anaEL.LoadSamples (datasets);	// now the list of datasets written in the xml file is known
  anaEL.LoadSSDiLeptonSelection (sel);	// now the parameters for the selection are given to the selection
  anaEL.LoadGeneralInfo (DataType, Luminosity, LumiError, verbosity);
  int flagOriginal=sel.GetFlagb();
  int methodOriginal=sel.GetMethodb();
  int systOriginal= sel.GetSystb();
  std::cout << " For btag : flag " << flagOriginal << ", method " << methodOriginal << ", syst " << systOriginal << std::endl;
  IPHCTree::NTEvent * event = 0;
  //Selection table
  SelectionTable selTable_allChannels (sel.GetCutList (), datasets, string ("*"));
  SelectionTable selTable_ee (sel.GetCutList (), datasets, string ("ee"));
  SelectionTable selTable_emu (sel.GetCutList (), datasets, string ("emu"));
  SelectionTable selTable_mumu (sel.GetCutList (), datasets, string ("mumu"));
  
  //Book keeping of standard histos
  bool doHistoManager = true;
  SSDiLepAnaHistoManager histoManager;
  if(doHistoManager){
  	histoManager.LoadDatasets (datasets);
  	histoManager.LoadSelectionSteps (sel.GetCutList ());
  	histoManager.LoadChannels (sel.GetChannelList ());
  	histoManager.CreateHistos ();
  }
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

    if(verbosity>2) cout<<"RPVAnalysis> Dataset: "<<datasets[d].Name()<<endl;
    datasets[d].eventTree ()->SetBranchAddress ("NTEvent", &event);

    unsigned int nEvents = static_cast<unsigned int>(datasets[d].eventTree ()->GetEntries ());
    cout << "RPVAnalysis> NEvents = " << nEvents << endl;
    cout << "RPVAnalysis> NEvents to run over = "<<datasets[d].NofEvtsToRunOver()<<endl;


    //LOOP OVER THE EVENTS
    for (unsigned int ievt = 0; ievt < datasets[d].NofEvtsToRunOver(); ievt++)
    {
      float weight = 1.;
      if(datasets[d].isData() == false) weight = datasets[d].NormFactor()*Luminosity; //if Data , weight = 1
      //cout<<"weight "<<weight<<" "<<datasets[d].isData()<<endl;
      datasets[d].eventTree()->GetEntry(ievt);
      IPHCTree::NTTransient::InitializeAfterReading(event);

      if (verbosity > 3)
      {
        std::cout << "event " << ievt 
                  <<" - event number=" << event->general.eventNb
                  <<" - run number="   << event->general.runNb
                  << std::endl;
      }
      if (ievt % 1000 == 0)
	      cout << "RPVAnalysis> Progress bar : " << ievt << endl;

      //Load event for the selection
      sel.LoadEvent(event);

      //Collection of selected objects
      vector < NTElectron > selElectrons = sel.GetSelectedElectrons ();
      //vector < NTMuon > selMuons = sel.GetScaledMuons ();
      vector < NTMuon > selMuons = sel.GetSelectedMuons ();
      vector < NTJet > selJets = sel.GetSelectedJets ();
      NTMET met = sel.GetSelectedMET ();	// no criteria applyied

      //Candidate pair of lepton
      string CandType="false";		// ee - emu - mumum or false
      vector < NTElectron > candElec;
      vector < NTMuon > candMuon;
      sel.GetLeptonPair (candMuon, candElec, CandType);	// fill the variables

      //integer which define the last step of the selection that the event fullfill
      int selLastStep = 0;
      int step = 0;

      //////////////////////////////////   
      //   Fill the selection table
      //////////////////////////////////   
      step = sel.FillTable (selTable_ee, &(datasets[d]), d, weight);
      if (CandType=="ee") selLastStep = step;
      step = sel.FillTable (selTable_emu, &(datasets[d]), d, weight);
      if (CandType=="emu") selLastStep = step;
      step = sel.FillTable (selTable_mumu, &(datasets[d]), d, weight);
      if (CandType=="mumu") selLastStep = step;

      //CandType="mumu";
      //cout<<selMuons.size()<<endl;
      histoManager.Fill(sel, event, selMuons, selElectrons, selLastStep, sel.GetChannel(CandType), d, weight);
      if (CandType=="mumu" && selLastStep >= AnaStep)  Nobs_mumu+=weight;

      //     selLastStep = sel.FillTable(selTable_allChannels, &(datasets[d]), d, weight);


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
  ofile << "\\begin{document}" << endl;
//  ofile << "\\usepackage{geometry}" << endl;
//  ofile << "\\geometry{a4paper, textwidth=19cm, textheight=23cm}" << endl;

 //Merge channels consistently
  vector < string > mergenames;
  mergenames.push_back ("QCD1");
  mergenames.push_back ("QCD2");
  mergenames.push_back ("QCD3");
  mergenames.push_back ("QCD4");
  mergenames.push_back ("QCD5");
  mergenames.push_back ("QCD6");
  mergenames.push_back ("QCD7");
  mergenames.push_back ("QCD8");
  selTable_ee.MergeDatasets (mergenames, string ("QCD1"));
  selTable_emu.MergeDatasets (mergenames, string ("QCD1"));
  selTable_mumu.MergeDatasets (mergenames, string ("QCD1"));
  mergenames.clear ();
  mergenames.push_back ("QCD9");
  mergenames.push_back ("QCD10");
  mergenames.push_back ("QCD11");
  mergenames.push_back ("QCD12");
  mergenames.push_back ("QCD13");
  selTable_ee.MergeDatasets (mergenames, string ("QCD2"));
  selTable_emu.MergeDatasets (mergenames, string ("QCD2"));
  selTable_mumu.MergeDatasets (mergenames, string ("QCD2"));
  mergenames.clear ();
  mergenames.push_back ("DYToEE1");
  mergenames.push_back ("DYToEE2");
  mergenames.push_back ("DYToMuMu1");
  mergenames.push_back ("DYToMuMu2");
  mergenames.push_back ("DYToTauTau1");
  mergenames.push_back ("DYToTauTau2");
  selTable_ee.MergeDatasets (mergenames, string ("DYToLL"));
  selTable_emu.MergeDatasets (mergenames, string ("DYToLL"));
  selTable_mumu.MergeDatasets (mergenames, string ("DYToLL"));

  //Define signal
  selTable_ee.DefineFirstDataset (string ("LM1"));
  selTable_emu.DefineFirstDataset (string ("LM1"));
  selTable_mumu.DefineFirstDataset (string ("LM1"));
//  selTable_allChannels.DefineFirstDataset (string ("TTbarSignal"));
  //Calculations
  selTable_ee.TableCalculator ();
  selTable_emu.TableCalculator ();
  selTable_mumu.TableCalculator ();
//  selTable_allChannels.TableCalculator ();
  //Write
  selTable_ee.Write (ofile);
  selTable_emu.Write (ofile);
  selTable_mumu.Write (ofile);
//  selTable_allChannels.Write (ofile);

  ofile << "\\end{document}" << endl;
  system (Form("pdflatex %s", ofilename.data()));


  if (verbosity > 0) {
    cout << "#########################" << endl;
    cout << " Write output root file " << endl;
    cout << "#########################" << endl;
  }

  string orootfilename;
  if(argc>2) orootfilename = argv[2];
  else orootfilename = string("RPVAnalysis.root");
  TFile *fout = new TFile (orootfilename.c_str(), "RECREATE");
  if(doHistoManager) histoManager.Write (fout);
  //fout->Write();
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
