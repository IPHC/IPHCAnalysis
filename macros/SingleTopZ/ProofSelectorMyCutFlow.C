#define ProofSelectorMyCutFlow_cxx

//////////////////////////////////////////////////////////
//
// Example of TSelector implementation to do a Monte Carlo
// generation using Pythia8.
// See tutorials/proof/runProof.C, option "pythia8", for an
// example of how to run this selector.
//
//////////////////////////////////////////////////////////
 
#include <TCanvas.h>
#include <TFrame.h>
#include <TPaveText.h>
#include <TFormula.h>
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TParameter.h>
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"

#include "ProofSelectorMyCutFlow.h"

//_____________________________________________________________________________
ProofSelectorMyCutFlow::ProofSelectorMyCutFlow()
{
  // Constructor
  
  
  
  /*
  string ofilenametex = "HighHTEvents.txt";
  ofile = ofstream(ofilenametex.c_str());

  ofile.precision(3);
  ofile.setf(ios::fixed);
  
  */
  
  fChain     = 0;
  branch     = 0;
  event      = 0;
  dataset    = 0;
  anaEL      = 0;
  verbosity  = 0;
  DataType   = 0;
  Luminosity = 0;
  //histos
  fHist      = 0;
   
  
  
  
  applyJES  = false;
  scale     = 0; // +1 or -1
  applyJER  = false;
  ResFactor = 0.1;
  
  
  
  applyFakescale   = true;
  
  SF_Fake.push_back(3.71); //mumumu
  SF_Fake.push_back(1.50); //mumue
  SF_Fake.push_back(3.69); //eemu
  SF_Fake.push_back(2.12); //eee
  

  applyWZ          = true;
   
  SF_WZ.push_back(0.62); //mumumu
  SF_WZ.push_back(0.68); //mumue
  SF_WZ.push_back(0.63); //eemu
  SF_WZ.push_back(0.70); //eee
  
  
  applyWZ_finalSel = true;
  SF_WZ_finalSel = 0.92; // 4 channels combined
  
  
  applyLeptonSF  = true;
  applyTrigger   = true;
  
  applyLeptonSFUp    = false;
  applyLeptonSFDown  = false;

  applyTriggerUp   = false;
  applyTriggerDown = false;
  
  
  IReweight		= true;
  IDYestimateWithMetCut = true;
  
  
  IReweight_puUp	= false;
  IReweight_puDown	= false;
  
  
  useNonIsoWcand = false;
  looseIso = 0.4; //0.4
  themetcut = 35;
  rand.SetSeed(102994949);
  
  
  doPDF = false;
  //pdftype =0 ;
  pdftype =1  ;
 
  doBTagCVScorr = true;
  //pdf.Initialize();
   doBTagCSV_syst = 0;

}

//_____________________________________________________________________________
ProofSelectorMyCutFlow::~ProofSelectorMyCutFlow()
{
  // Destructor
  
  //SafeDelete(fHist);
}

//_____________________________________________________________________________
void ProofSelectorMyCutFlow::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.
  
  //fChain point to the loading tree 
  fChain = tree;
  cout << "start init tree " << endl;
  // Set branch addresses
  branch = (TBranch *) tree->GetBranch("NTEvent");
  event = new IPHCTree::NTEvent();
   branch->SetAddress(&event);
   //event is now retrieved and could be used in Process
   cout << "end init tree " << endl;
}

//_____________________________________________________________________________
void ProofSelectorMyCutFlow::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  cout << "start Begin " << endl;
  TString option = GetOption();
  cout << "end  Begin" << endl;
  
  
}

//_____________________________________________________________________________
void ProofSelectorMyCutFlow::SlaveBegin(TTree * tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  
  cout << "start SlaveBegin " << endl;
  TString option = GetOption();
  //--------------------------------------//
  //       Loading the xml file
  //--------------------------------------//
  TNamed *dsname = (TNamed *) fInput->FindObject("PROOF_DATASETNAME"); 
  datasetName = dsname->GetTitle();
  cout << "dataset name " << datasetName << endl;
  cout << "dataset name " << datasetName << endl;
  cout << "dataset name " << datasetName << endl;
  cout << "dataset name " << datasetName << endl;
  cout << "dataset name " << datasetName << endl;
  cout << "dataset name " << datasetName << endl;
  
  
  
  TNamed *xfname = (TNamed *) fInput->FindObject("PROOF_XMLFILENAME");
  string xmlFileName = xfname->GetTitle();
  anaEL = new AnalysisEnvironmentLoader(xmlFileName.c_str());
  
  
  
  anaEL->LoadSamples (datasets, datasetName); // now the list of datasets written in the xml file is known
  
  
  
  //retrieve the current dataset according to its name
  for(unsigned int d=0;d<datasets.size();d++){
    cout << "datasets.size() " << datasets.size()<< "  datasets[d].Name()" << datasets[d].Name()  << " datasetName "
	 <<datasetName  << endl;
    if(datasets[d].Name()==datasetName)dataset = &datasets[d];
  }
  cout << "load datasets "  << endl;
  anaEL->LoadDiLeptonSelection (sel); // now the parameters for the selection are given to the selection
  cout << "datasets loaded "  << endl;
  //anaEL->LoadGeneralInfo(DataType, Luminosity, verbosity );
  anaEL->LoadGeneralInfo(DataType, Luminosity, LumiError, PUWeightFileName, verbosity );
  
  //Load for PU:
  sel.GeneratePUWeight(PUWeightFileName);
  
  //******************************************
  //Load Scale Factors for lepton efficiencies
  //******************************************
  sel.LoadElScaleFactors();
  sel.LoadMuScaleFactors();
  sel.InitJESUnc("Total"); 
  
  
  anaEL->LoadWeight (sel); // now the parameters for SFBweight are initialized (for b-tag!)

   //--------------------------------------//
   //   Fill cuts and channels  	
   //--------------------------------------//
  CutName.push_back("Cut1");
   
  //--------------------------------------//
  //   Initializing variables 	
  //--------------------------------------//
  
  
  
  SF_trig_mumumu = 1.012; 
  SF_trig_mumue= 1.003; 
  SF_trig_eemu = 0.971; 
  SF_trig_eee = 0.959;
  
  
  SF_trig_mumumu_error  = 0.029; 
  SF_trig_mumue_error = 0.031;
  SF_trig_eemu_error = 0.038; 
  SF_trig_eee_error  = 0.038; 
  
  
 
  
  //**************************************
  //**************************************
  //******* fakes DD estimate ************
  //**************************************
  //**************************************
  
  
  
  sumSFlept_mumumu  = 0;
  sumSFlept_mumue   = 0;
  sumSFlept_eemu    = 0;
  sumSFlept_eee     = 0;
  
  nEvents_mumumu  = 0;
  nEvents_mumue   = 0;
  nEvents_eemu    = 0;
  nEvents_eee     = 0;
  
  
  scaleElec = 1.0; // 1 to switch off
  resolElec = 0.0; // 0 to switch off
  
  
  ITypeMC     = -1;
  ICut        = -1;  
  
  
  //************************************
  //For trigger systematics 
  
  if(applyTriggerUp){
    SF_trig_mumumu += SF_trig_mumumu_error;
    SF_trig_mumue+= SF_trig_mumue_error;  
    SF_trig_eemu += SF_trig_eemu_error;
    SF_trig_eee += SF_trig_eee_error; 
  } 
  if(applyTriggerDown){
    SF_trig_mumumu  -= SF_trig_mumumu_error;
    SF_trig_mumue -= SF_trig_mumue_error;  
    SF_trig_eemu  -= SF_trig_eemu_error;
    SF_trig_eee  -= SF_trig_eee_error; 
  } 
  
  
  
  for(unsigned int d=0;d<datasets.size();d++){
    cout << "datasets.size() " << datasets.size()<< "  datasets[d].Name()" << datasets[d].Name()  << " datasetName "
	 <<datasetName  << endl;
    if(datasets[d].Name()==datasetName)dataset = &datasets[d];
  }
  
  
  
  SF_BranchingRatio_ll = (0.108*9.)*(0.108*9.);
  SF_BranchingRatio_lj = (0.108*9.)*(0.676*1.5);
  SF_BranchingRatio_had = (0.676*1.5)*(0.676*1.5);
  
  
  
  //--------------------------------------//
  //   Managing histos  	
  //--------------------------------------//
  MyhistoManager.LoadDatasets(datasets);   
  MyhistoManager.LoadSelectionSteps(CutName);
  MyhistoManager.LoadChannels(TheChannelName);
  //example
    
  nbins = 200;
  minx = 0.;
  maxx = 350;
  
  
  //***********************
  // initiate lumi reweighting
  MyhistoManager.CreateHisto(SelABjet, "SelABjet", datasetName, "Nevents","Entries", 2, -0.5, 1.5);
  MyhistoManager.CreateHisto(SelABjet_afterjetsel, "SelABjet_afterjetsel", datasetName, "Nevents","Entries", 2, -0.5, 1.5);
 
  MyhistoManager.CreateHisto(Ntrilept_mumumu, "Ntrilept_mumumu", datasetName, "Nevents","Entries", 11, -0.5, 10.5); 
  MyhistoManager.CreateHisto(Ntrileptnoniso_mumumu, "Ntrileptnoniso_mumumu", datasetName, "Nevents","Entries", 11, -0.5, 10.5); 
  
  
  MyhistoManager.CreateHisto(Nvertex, "Nvertex", datasetName, "Nvertex", "Entries", 50, 0, 50); 
  
  MyhistoManager.CreateHisto(CutFlow_mumumu,  "CutFlow_mumumu" ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_mumue,   "CutFlow_mumue"  ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_eemu,    "CutFlow_eemu"   ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_eee,     "CutFlow_eee"    ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  
  
  MyhistoManager.SetCutFlowAxisTitleFCNCMonotop(CutFlow_mumumu,   "CutFlow_mumumu"  ,datasetName);
  MyhistoManager.SetCutFlowAxisTitleFCNCMonotop(CutFlow_mumue,    "CutFlow_mumue"   ,datasetName);
  MyhistoManager.SetCutFlowAxisTitleFCNCMonotop(CutFlow_eemu,     "CutFlow_eemu"    ,datasetName);
  MyhistoManager.SetCutFlowAxisTitleFCNCMonotop(CutFlow_eee,      "CutFlow_eee"     ,datasetName);
  
  
  MyhistoManager.CreateHisto(PU_before_mumumu, "PU_before_mumumu", datasetName, "Npileup", "Entries", 30, 0, 30); 
  MyhistoManager.CreateHisto(PU_before_mumue, "PU_before_mumue", datasetName, "Npileup", "Entries", 30, 0, 30); 
  MyhistoManager.CreateHisto(PU_before_eemu, "PU_before_eemu", datasetName, "Npileup", "Entries", 30, 0, 30); 
  MyhistoManager.CreateHisto(PU_before_eee, "PU_before_eee", datasetName, "Npileup", "Entries", 30, 0, 30);   
  
  MyhistoManager.CreateHisto(PU_intime_mumumu, "PU_intime_mumumu", datasetName, "Npileup", "Entries", 30, 0, 30); 
  MyhistoManager.CreateHisto(PU_intime_mumue, "PU_intime_mumue", datasetName, "Npileup", "Entries", 30, 0, 30); 
  MyhistoManager.CreateHisto(PU_intime_eemu, "PU_intime_eemu", datasetName, "Npileup", "Entries", 30, 0, 30); 
  MyhistoManager.CreateHisto(PU_intime_eee, "PU_intime_eee", datasetName, "Npileup", "Entries", 30, 0, 30);   
  
  MyhistoManager.CreateHisto(PU_after_mumumu, "PU_after_mumumu", datasetName, "Npileup", "Entries", 30, 0, 30); 
  MyhistoManager.CreateHisto(PU_after_mumue, "PU_after_mumue", datasetName, "Npileup", "Entries", 30, 0, 30); 
  MyhistoManager.CreateHisto(PU_after_eemu, "PU_after_eemu", datasetName, "Npileup", "Entries", 30, 0, 30); 
  MyhistoManager.CreateHisto(PU_after_eee, "PU_after_eee", datasetName, "Npileup", "Entries", 30, 0, 30);   
  

  MyhistoManager.CreateHisto(NVtx_mumumu, "NVtx_mumumu", datasetName, "Nvertex", "Entries", 31, 0, 30); 
  MyhistoManager.CreateHisto(NVtx_mumue,  "NVtx_mumue",  datasetName, "Nvertex", "Entries", 31, 0, 30); 
  MyhistoManager.CreateHisto(NVtx_eemu,   "NVtx_eemu",   datasetName, "Nvertex", "Entries", 31, 0, 30); 
  MyhistoManager.CreateHisto(NVtx_eee,    "NVtx_eee",    datasetName, "Nvertex", "Entries", 31, 0, 30); 
  
  MyhistoManager.CreateHisto(Nvtx_mumumu_afterleptsel, "Nvtx_mumumu_afterleptsel", datasetName , "NVtx", "Entries", 50, 0, 49) ;
  MyhistoManager.CreateHisto(Nvtx_mumue_afterleptsel,  "Nvtx_mumue_afterleptsel",  datasetName , "NVtx", "Entries", 50, 0, 49);
  MyhistoManager.CreateHisto(Nvtx_eemu_afterleptsel,   "Nvtx_eemu_afterleptsel",   datasetName , "NVtx", "Entries", 50, 0, 49);
  MyhistoManager.CreateHisto(Nvtx_eee_afterleptsel,    "Nvtx_eee_afterleptsel",    datasetName , "NVtx", "Entries", 50, 0, 49);
  
  MyhistoManager.CreateHisto(NVtx_mumumu_aftertrigsel, "NVtx_mumumu_aftertrigsel", datasetName, "Nvertex", "Entries", 31, 0, 30); 
  MyhistoManager.CreateHisto(NVtx_mumue_aftertrigsel,  "NVtx_mumue_aftertrigsel",  datasetName, "Nvertex", "Entries", 31, 0, 30); 
  MyhistoManager.CreateHisto(NVtx_eemu_aftertrigsel,   "NVtx_eemu_aftertrigsel",   datasetName, "Nvertex", "Entries", 31, 0, 30); 
  MyhistoManager.CreateHisto(NVtx_eee_aftertrigsel,    "NVtx_eee_aftertrigsel",    datasetName, "Nvertex", "Entries", 31, 0, 30); 
  
  MyhistoManager.CreateHisto(NVtx_mumumu_afterleptsel, "NVtx_mumumu_afterleptsel", datasetName, "Nvertex", "Entries", 31, 0, 30); 
  MyhistoManager.CreateHisto(NVtx_mumue_afterleptsel,  "NVtx_mumue_afterleptsel",  datasetName, "Nvertex", "Entries", 31, 0, 30); 
  MyhistoManager.CreateHisto(NVtx_eemu_afterleptsel,   "NVtx_eemu_afterleptsel",   datasetName, "Nvertex", "Entries", 31, 0, 30); 
  MyhistoManager.CreateHisto(NVtx_eee_afterleptsel,    "NVtx_eee_afterleptsel",    datasetName, "Nvertex", "Entries", 31, 0, 30); 

  
  
  
  
  MyhistoManager.CreateHisto(ErrCutFlow_mumumu,  "ErrCutFlow_mumumu"  ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_mumue,   "ErrCutFlow_mumue"   ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_eemu,    "ErrCutFlow_eemu"    ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_eee,     "ErrCutFlow_eee"     ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  
  
  
  MyhistoManager.CreateHisto(Mt_mumumu_afterbjetsel, "Mt_mumumu_afterbjetsel", datasetName,"Mt","Entries", 100, 0, 500); 
  MyhistoManager.CreateHisto(Mt_mumue_afterbjetsel , "Mt_mumue_afterbjetsel" , datasetName,"Mt","Entries", 100, 0, 500);
  MyhistoManager.CreateHisto(Mt_eemu_afterbjetsel  , "Mt_eemu_afterbjetsel"  , datasetName,"Mt","Entries", 100, 0, 500);
  MyhistoManager.CreateHisto(Mt_eee_afterbjetsel   , "Mt_eee_afterbjetsel"   , datasetName,"Mt","Entries", 100, 0, 500);
   
  MyhistoManager.CreateHisto(Mt_mumumu_afterbjetveto, "Mt_mumumu_afterbjetveto", datasetName,"Mt","Entries", 100, 0, 500); 
  MyhistoManager.CreateHisto(Mt_mumue_afterbjetveto , "Mt_mumue_afterbjetveto" , datasetName,"Mt","Entries", 100, 0, 500);
  MyhistoManager.CreateHisto(Mt_eemu_afterbjetveto  , "Mt_eemu_afterbjetveto"  , datasetName,"Mt","Entries", 100, 0, 500);
  MyhistoManager.CreateHisto(Mt_eee_afterbjetveto   , "Mt_eee_afterbjetveto"   , datasetName,"Mt","Entries", 100, 0, 500);
   
  
  MyhistoManager.CreateHisto(NJet_mumumu_afterZsel, "NJet_mumumu_afterZsel", datasetName,"Njets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NJet_mumue_afterZsel , "NJet_mumue_afterZsel" , datasetName,"Njets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NJet_eemu_afterZsel  , "NJet_eemu_afterZsel"  , datasetName,"Njets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NJet_eee_afterZsel   , "NJet_eee_afterZsel"   , datasetName,"Njets", "Entries",5,-0.5,4.5);
  
  MyhistoManager.CreateHisto(NJet_mumumu_afterbsel, "NJet_mumumu_afterbsel", datasetName,"Njets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NJet_mumue_afterbsel , "NJet_mumue_afterbsel" , datasetName,"Njets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NJet_eemu_afterbsel  , "NJet_eemu_afterbsel"  , datasetName,"Njets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NJet_eee_afterbsel   , "NJet_eee_afterbsel"   , datasetName,"Njets", "Entries",5,-0.5,4.5);
  
  MyhistoManager.CreateHisto(NJet_mumumu_afterleptsel_mWT110, "NJet_mumumu_afterleptsel_mWT110", datasetName,"NBjets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NJet_mumue_afterleptsel_mWT110,  "NJet_mumue_afterleptsel_mWT110",  datasetName,"NBjets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NJet_eemu_afterleptsel_mWT110,   "NJet_eemu_afterleptsel_mWT110",   datasetName,"NBjets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NJet_eee_afterleptsel_mWT110,    "NJet_eee_afterleptsel_mWT110",	 datasetName,"NBjets", "Entries",5,-0.5,4.5);
  
  MyhistoManager.CreateHisto(NLept_mumumu_afterbsel, "NLept_mumumu_afterbsel", datasetName,"NLepts", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NLept_mumue_afterbsel , "NLept_mumue_afterbsel" , datasetName,"NLepts", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NLept_eemu_afterbsel  , "NLept_eemu_afterbsel"  , datasetName,"NLepts", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NLept_eee_afterbsel   , "NLept_eee_afterbsel"   , datasetName,"NLepts", "Entries",5,-0.5,4.5);
  
  
  
  MyhistoManager.CreateHisto(NBJet_mumumu_afterZsel, "NBJet_mumumu_afterZsel", datasetName,"NBjets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NBJet_mumue_afterZsel , "NBJet_mumue_afterZsel" , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NBJet_eemu_afterZsel  , "NBJet_eemu_afterZsel"  , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NBJet_eee_afterZsel   , "NBJet_eee_afterZsel"   , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  
  MyhistoManager.CreateHisto(NBJet_mumumu_afterjetsel, "NBJet_mumumu_afterjetsel", datasetName,"NBjets", "Entries",2,-0.5,1.5);
  MyhistoManager.CreateHisto(NBJet_mumue_afterjetsel , "NBJet_mumue_afterjetsel" , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  MyhistoManager.CreateHisto(NBJet_eemu_afterjetsel  , "NBJet_eemu_afterjetsel"  , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  MyhistoManager.CreateHisto(NBJet_eee_afterjetsel   , "NBJet_eee_afterjetsel"   , datasetName,"NBjets", "Entries",2,-0.5,1.5);	
  
  
  
  
  
  MyhistoManager.CreateHisto(NBJet_mumumu_afterjetsel_bjets, "NBJet_mumumu_afterjetsel_bjets", datasetName,"NBjets", "Entries",2,-0.5,1.5);
  MyhistoManager.CreateHisto(NBJet_mumue_afterjetsel_bjets , "NBJet_mumue_afterjetsel_bjets" , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  MyhistoManager.CreateHisto(NBJet_eemu_afterjetsel_bjets  , "NBJet_eemu_afterjetsel_bjets"  , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  MyhistoManager.CreateHisto(NBJet_eee_afterjetsel_bjets   , "NBJet_eee_afterjetsel_bjets"   , datasetName,"NBjets", "Entries",2,-0.5,1.5);	
  
  
  MyhistoManager.CreateHisto(NBJet_mumumu_afterjetsel_cjets, "NBJet_mumumu_afterjetsel_cjets", datasetName,"NBjets", "Entries",2,-0.5,1.5);
  MyhistoManager.CreateHisto(NBJet_mumue_afterjetsel_cjets , "NBJet_mumue_afterjetsel_cjets" , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  MyhistoManager.CreateHisto(NBJet_eemu_afterjetsel_cjets  , "NBJet_eemu_afterjetsel_cjets"  , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  MyhistoManager.CreateHisto(NBJet_eee_afterjetsel_cjets   , "NBJet_eee_afterjetsel_cjets"   , datasetName,"NBjets", "Entries",2,-0.5,1.5);	
  
  
  MyhistoManager.CreateHisto(NBJet_mumumu_afterjetsel_ljets, "NBJet_mumumu_afterjetsel_ljets", datasetName,"NBjets", "Entries",2,-0.5,1.5);
  MyhistoManager.CreateHisto(NBJet_mumue_afterjetsel_ljets , "NBJet_mumue_afterjetsel_ljets" , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  MyhistoManager.CreateHisto(NBJet_eemu_afterjetsel_ljets  , "NBJet_eemu_afterjetsel_ljets"  , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  MyhistoManager.CreateHisto(NBJet_eee_afterjetsel_ljets   , "NBJet_eee_afterjetsel_ljets"   , datasetName,"NBjets", "Entries",2,-0.5,1.5);	
  
  
  
  
  MyhistoManager.CreateHisto(BJetDiscri_mumumu_afterjetsel_bjets, "BJetDiscri_mumumu_afterjetsel_bjets", datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  MyhistoManager.CreateHisto(BJetDiscri_mumue_afterjetsel_bjets , "BJetDiscri_mumue_afterjetsel_bjets" , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  MyhistoManager.CreateHisto(BJetDiscri_eemu_afterjetsel_bjets  , "BJetDiscri_eemu_afterjetsel_bjets"  , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  MyhistoManager.CreateHisto(BJetDiscri_eee_afterjetsel_bjets   , "BJetDiscri_eee_afterjetsel_bjets"   , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  
  
  
  MyhistoManager.CreateHisto(BJetDiscri_mumumu_afterjetsel_cjets, "BJetDiscri_mumumu_afterjetsel_cjets", datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  MyhistoManager.CreateHisto(BJetDiscri_mumue_afterjetsel_cjets , "BJetDiscri_mumue_afterjetsel_cjets" , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  MyhistoManager.CreateHisto(BJetDiscri_eemu_afterjetsel_cjets  , "BJetDiscri_eemu_afterjetsel_cjets"  , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  MyhistoManager.CreateHisto(BJetDiscri_eee_afterjetsel_cjets   , "BJetDiscri_eee_afterjetsel_cjets"   , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
 
  
  MyhistoManager.CreateHisto(BJetDiscri_mumumu_afterjetsel_ljets, "BJetDiscri_mumumu_afterjetsel_ljets", datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  MyhistoManager.CreateHisto(BJetDiscri_mumue_afterjetsel_ljets , "BJetDiscri_mumue_afterjetsel_ljets" , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  MyhistoManager.CreateHisto(BJetDiscri_eemu_afterjetsel_ljets  , "BJetDiscri_eemu_afterjetsel_ljets"  , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  MyhistoManager.CreateHisto(BJetDiscri_eee_afterjetsel_ljets   , "BJetDiscri_eee_afterjetsel_ljets"   , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
 
  
  
  
  MyhistoManager.CreateHisto(NBJet_mumumu_afterleptsel_mWT110, "NBJet_mumumu_afterleptsel_mWT110", datasetName,"NBjets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NBJet_mumue_afterleptsel_mWT110,  "NBJet_mumue_afterleptsel_mWT110",  datasetName,"NBjets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NBJet_eemu_afterleptsel_mWT110,   "NBJet_eemu_afterleptsel_mWT110",   datasetName,"NBjets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NBJet_eee_afterleptsel_mWT110,    "NBJet_eee_afterleptsel_mWT110",    datasetName,"NBjets", "Entries",5,-0.5,4.5);
  
  
  
  
  
  /*MyhistoManager.CreateHisto(NBJet_mumumu_afterjetsel, "NBJet_mumumu_afterjetsel", datasetName,"NBjets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NBJet_mumue_afterjetsel , "NBJet_mumue_afterjetsel" , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NBJet_eemu_afterjetsel  , "NBJet_eemu_afterjetsel"  , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  MyhistoManager.CreateHisto(NBJet_eee_afterjetsel   , "NBJet_eee_afterjetsel"   , datasetName,"NBjets", "Entries",5,-0.5,4.5); 
  */
   
  
  MyhistoManager.CreateHisto(InvM_ll_mumumu_afterleptsel_highSumPt, "InvM_ll_mumumu_afterleptsel_highSumPt" , datasetName,"Minv", "Entries",100,0.,1000);                     
  MyhistoManager.CreateHisto(InvM_ll_mumue_afterleptsel_highSumPt,  "InvM_ll_mumue_afterleptsel_highSumPt"  , datasetName,"Minv", "Entries",100,0.,1000);
  MyhistoManager.CreateHisto(InvM_ll_eemu_afterleptsel_highSumPt,   "InvM_ll_eemu_afterleptsel_highSumPt"   , datasetName,"Minv", "Entries",100,0.,1000);
  MyhistoManager.CreateHisto(InvM_ll_eee_afterleptsel_highSumPt,    "InvM_ll_eee_afterleptsel_highSumPt"    , datasetName,"Minv", "Entries",100,0.,1000);

 
  MyhistoManager.CreateHisto(InvM_ll_mumumu_afterleptsel, "InvM_ll_mumumu_afterleptsel" , datasetName,"Minv", "Entries",350,0.,1000);
  MyhistoManager.CreateHisto(InvM_ll_mumue_afterleptsel,  "InvM_ll_mumue_afterleptsel"  , datasetName,"Minv", "Entries",350,0.,1000);
  MyhistoManager.CreateHisto(InvM_ll_eemu_afterleptsel,   "InvM_ll_eemu_afterleptsel"   , datasetName,"Minv", "Entries",350,0.,1000);
  MyhistoManager.CreateHisto(InvM_ll_eee_afterleptsel,    "InvM_ll_eee_afterleptsel"    , datasetName,"Minv", "Entries",350,0.,1000);
  
  MyhistoManager.CreateHisto(InvM_ll_mumumu_afterleptsel_mWT110, "InvM_ll_mumumu_afterleptsel_mWT110" , datasetName,"Minv", "Entries",350,0.,1000);
  MyhistoManager.CreateHisto(InvM_ll_mumue_afterleptsel_mWT110,  "InvM_ll_mumue_afterleptsel_mWT110"  , datasetName,"Minv", "Entries",350,0.,1000);
  MyhistoManager.CreateHisto(InvM_ll_eemu_afterleptsel_mWT110,   "InvM_ll_eemu_afterleptsel_mWT110"   , datasetName,"Minv", "Entries",350,0.,1000);
  MyhistoManager.CreateHisto(InvM_ll_eee_afterleptsel_mWT110,    "InvM_ll_eee_afterleptsel_mWT110"    , datasetName,"Minv", "Entries",350,0.,1000);
 
  MyhistoManager.CreateHisto(InvM_ll_mumumu_afterleptsel_lowbin, "InvM_ll_mumumu_afterleptsel_lowbin" , datasetName,"Minv", "Entries",100,0.,200);
  MyhistoManager.CreateHisto(InvM_ll_mumue_afterleptsel_lowbin,  "InvM_ll_mumue_afterleptsel_lowbin"  , datasetName,"Minv", "Entries",100,0.,200);
  MyhistoManager.CreateHisto(InvM_ll_eemu_afterleptsel_lowbin,   "InvM_ll_eemu_afterleptsel_lowbin"   , datasetName,"Minv", "Entries",100,0.,200);
  MyhistoManager.CreateHisto(InvM_ll_eee_afterleptsel_lowbin,    "InvM_ll_eee_afterleptsel_lowbin"    , datasetName,"Minv", "Entries",100,0.,200);

  MyhistoManager.CreateHisto(InvM_ll_mumumu_afterjetsel, "InvM_ll_mumumu_afterjetsel" , datasetName,"Minv", "Entries",350, 0., 1000);
  MyhistoManager.CreateHisto(InvM_ll_mumue_afterjetsel,  "InvM_ll_mumue_afterjetsel"  , datasetName,"Minv", "Entries",350, 0., 1000);
  MyhistoManager.CreateHisto(InvM_ll_eemu_afterjetsel,   "InvM_ll_eemu_afterjetsel"   , datasetName,"Minv", "Entries",350, 0., 1000);
  MyhistoManager.CreateHisto(InvM_ll_eee_afterjetsel,    "InvM_ll_eee_afterjetsel"    , datasetName,"Minv", "Entries",350, 0., 1000);
  
  MyhistoManager.CreateHisto(InvM_ll_mumumu_afterbjetsel, "InvM_ll_mumumu_afterbjetsel" , datasetName,"Minv", "Entries",350,0.,1000);
  MyhistoManager.CreateHisto(InvM_ll_mumue_afterbjetsel,  "InvM_ll_mumue_afterbjetsel"  , datasetName,"Minv", "Entries",350,0.,1000);
  MyhistoManager.CreateHisto(InvM_ll_eemu_afterbjetsel,   "InvM_ll_eemu_afterbjetsel"   , datasetName,"Minv", "Entries",350,0.,1000);
  MyhistoManager.CreateHisto(InvM_ll_eee_afterbjetsel,    "InvM_ll_eee_afterbjetsel"    , datasetName,"Minv", "Entries",350,0.,1000);

    
  
  MyhistoManager.CreateHisto(LeptPt_mumumu_afterleptsel, "LeptPt_mumumu_afterleptsel", datasetName, "LeptPt", "Entries",350,0., 1000) ;
  MyhistoManager.CreateHisto(LeptPt_mumue_afterleptsel,  "LeptPt_mumue_afterleptsel",  datasetName, "LeptPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptPt_eemu_afterleptsel,   "LeptPt_eemu_afterleptsel",   datasetName, "LeptPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptPt_eee_afterleptsel,    "LeptPt_eee_afterleptsel",    datasetName, "LeptPt", "Entries",350,0., 1000);
  
  MyhistoManager.CreateHisto(LeptPt_mumumu_afterjetsel, "LeptPt_mumumu_afterjetsel", datasetName, "LeptPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptPt_mumue_afterjetsel,  "LeptPt_mumue_afterjetsel",  datasetName, "LeptPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptPt_eemu_afterjetsel,   "LeptPt_eemu_afterjetsel",   datasetName, "LeptPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptPt_eee_afterjetsel,    "LeptPt_eee_afterjetsel",    datasetName, "LeptPt", "Entries",350,0., 1000);
  
  MyhistoManager.CreateHisto(LeptPt_mumumu_afterbjetsel, "LeptPt_mumumu_afterbjetsel", datasetName, "LeptPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptPt_mumue_afterbjetsel,  "LeptPt_mumue_afterbjetsel",  datasetName, "LeptPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptPt_eemu_afterbjetsel,   "LeptPt_eemu_afterbjetsel",   datasetName, "LeptPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptPt_eee_afterbjetsel,    "LeptPt_eee_afterbjetsel",    datasetName, "LeptPt", "Entries",350,0., 1000);
    
  MyhistoManager.CreateHisto(LeptZPt_mumumu_afterleptsel, "LeptZPt_mumumu_afterleptsel", datasetName, "LeptZPt", "Entries",350,0., 1000) ;
  MyhistoManager.CreateHisto(LeptZPt_mumue_afterleptsel,  "LeptZPt_mumue_afterleptsel",  datasetName, "LeptZPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptZPt_eemu_afterleptsel,   "LeptZPt_eemu_afterleptsel",   datasetName, "LeptZPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptZPt_eee_afterleptsel,    "LeptZPt_eee_afterleptsel",    datasetName, "LeptZPt", "Entries",350,0., 1000);
  
  MyhistoManager.CreateHisto(LeptZPt_mumumu_afterjetsel, "LeptZPt_mumumu_afterjetsel", datasetName, "LeptZPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptZPt_mumue_afterjetsel,  "LeptZPt_mumue_afterjetsel",  datasetName, "LeptZPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptZPt_eemu_afterjetsel,   "LeptZPt_eemu_afterjetsel",   datasetName, "LeptZPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptZPt_eee_afterjetsel,    "LeptZPt_eee_afterjetsel",    datasetName, "LeptZPt", "Entries",350,0., 1000);
  
  MyhistoManager.CreateHisto(LeptZPt_mumumu_afterbjetsel, "LeptZPt_mumumu_afterbjetsel", datasetName, "LeptZPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptZPt_mumue_afterbjetsel,  "LeptZPt_mumue_afterbjetsel",  datasetName, "LeptZPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptZPt_eemu_afterbjetsel,   "LeptZPt_eemu_afterbjetsel",   datasetName, "LeptZPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptZPt_eee_afterbjetsel,    "LeptZPt_eee_afterbjetsel",    datasetName, "LeptZPt", "Entries",350,0., 1000);
    
  MyhistoManager.CreateHisto(LeptWPt_mumumu_afterleptsel, "LeptWPt_mumumu_afterleptsel", datasetName, "LeptWPt", "Entries",350,0., 1000) ;
  MyhistoManager.CreateHisto(LeptWPt_mumue_afterleptsel,  "LeptWPt_mumue_afterleptsel",  datasetName, "LeptWPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptWPt_eemu_afterleptsel,   "LeptWPt_eemu_afterleptsel",   datasetName, "LeptWPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptWPt_eee_afterleptsel,    "LeptWPt_eee_afterleptsel",    datasetName, "LeptWPt", "Entries",350,0., 1000);
  
  MyhistoManager.CreateHisto(LeptWPt_mumumu_afterjetsel, "LeptWPt_mumumu_afterjetsel", datasetName, "LeptWPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptWPt_mumue_afterjetsel,  "LeptWPt_mumue_afterjetsel",  datasetName, "LeptWPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptWPt_eemu_afterjetsel,   "LeptWPt_eemu_afterjetsel",   datasetName, "LeptWPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptWPt_eee_afterjetsel,    "LeptWPt_eee_afterjetsel",    datasetName, "LeptWPt", "Entries",350,0., 1000);
  
  MyhistoManager.CreateHisto(LeptWPt_mumumu_afterbjetsel, "LeptWPt_mumumu_afterbjetsel", datasetName, "LeptWPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptWPt_mumue_afterbjetsel,  "LeptWPt_mumue_afterbjetsel",  datasetName, "LeptWPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptWPt_eemu_afterbjetsel,   "LeptWPt_eemu_afterbjetsel",   datasetName, "LeptWPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptWPt_eee_afterbjetsel,    "LeptWPt_eee_afterbjetsel",    datasetName, "LeptWPt", "Entries",350,0., 1000);

  MyhistoManager.CreateHisto(LeptWPt_mumumu_afterbjetveto, "LeptWPt_mumumu_afterbjetveto", datasetName, "LeptWPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptWPt_mumue_afterbjetveto,  "LeptWPt_mumue_afterbjetveto",  datasetName, "LeptWPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptWPt_eemu_afterbjetveto,   "LeptWPt_eemu_afterbjetveto",   datasetName, "LeptWPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptWPt_eee_afterbjetveto,    "LeptWPt_eee_afterbjetveto",    datasetName, "LeptWPt", "Entries",350,0., 1000);

      
  MyhistoManager.CreateHisto(LeptWPt_mumumu_afterleptsel_mWT110, "LeptWPt_mumumu_afterleptsel_mWT110", datasetName, "LeptWPt", "Entries",350,0., 1000) ;
  MyhistoManager.CreateHisto(LeptWPt_mumue_afterleptsel_mWT110,  "LeptWPt_mumue_afterleptsel_mWT110",  datasetName, "LeptWPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptWPt_eemu_afterleptsel_mWT110,   "LeptWPt_eemu_afterleptsel_mWT110",   datasetName, "LeptWPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(LeptWPt_eee_afterleptsel_mWT110,    "LeptWPt_eee_afterleptsel_mWT110",    datasetName, "LeptWPt", "Entries",350,0., 1000);
  

  
  
  
  MyhistoManager.CreateHisto(JetPt_mumumu_afterleptsel, "JetPt_mumumu_afterleptsel", datasetName, "JetPt", "Entries",350,0., 1000) ;
  MyhistoManager.CreateHisto(JetPt_mumue_afterleptsel,  "JetPt_mumue_afterleptsel",  datasetName, "JetPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(JetPt_eemu_afterleptsel,   "JetPt_eemu_afterleptsel",   datasetName, "JetPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(JetPt_eee_afterleptsel,    "JetPt_eee_afterleptsel",    datasetName, "JetPt", "Entries",350,0., 1000);
  
  MyhistoManager.CreateHisto(JetPt_mumumu_afterjetsel, "JetPt_mumumu_afterjetsel", datasetName, "JetPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(JetPt_mumue_afterjetsel,  "JetPt_mumue_afterjetsel",  datasetName, "JetPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(JetPt_eemu_afterjetsel,   "JetPt_eemu_afterjetsel",   datasetName, "JetPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(JetPt_eee_afterjetsel,    "JetPt_eee_afterjetsel",    datasetName, "JetPt", "Entries",350,0., 1000);
  
  MyhistoManager.CreateHisto(JetPt_mumumu_afterbjetsel, "JetPt_mumumu_afterbjetsel", datasetName, "JetPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(JetPt_mumue_afterbjetsel,  "JetPt_mumue_afterbjetsel",  datasetName, "JetPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(JetPt_eemu_afterbjetsel,   "JetPt_eemu_afterbjetsel",   datasetName, "JetPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(JetPt_eee_afterbjetsel,    "JetPt_eee_afterbjetsel",    datasetName, "JetPt", "Entries",350,0., 1000);
  
  MyhistoManager.CreateHisto(JetPt_mumumu_afterbjetveto, "JetPt_mumumu_afterbjetveto", datasetName, "JetPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(JetPt_mumue_afterbjetveto,  "JetPt_mumue_afterbjetveto",  datasetName, "JetPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(JetPt_eemu_afterbjetveto,   "JetPt_eemu_afterbjetveto",   datasetName, "JetPt", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(JetPt_eee_afterbjetveto,    "JetPt_eee_afterbjetveto",    datasetName, "JetPt", "Entries",350,0., 1000);
  
    
  MyhistoManager.CreateHisto(JetEta_mumumu_afterleptsel, "JetEta_mumumu_afterleptsel", datasetName, "JetEta", "Entries",26, -2.5, 2.5) ;
  MyhistoManager.CreateHisto(JetEta_mumue_afterleptsel,  "JetEta_mumue_afterleptsel",  datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  MyhistoManager.CreateHisto(JetEta_eemu_afterleptsel,   "JetEta_eemu_afterleptsel",   datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  MyhistoManager.CreateHisto(JetEta_eee_afterleptsel,    "JetEta_eee_afterleptsel",    datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  
  MyhistoManager.CreateHisto(JetEta_mumumu_afterjetsel, "JetEta_mumumu_afterjetsel", datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  MyhistoManager.CreateHisto(JetEta_mumue_afterjetsel,  "JetEta_mumue_afterjetsel",  datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  MyhistoManager.CreateHisto(JetEta_eemu_afterjetsel,   "JetEta_eemu_afterjetsel",   datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  MyhistoManager.CreateHisto(JetEta_eee_afterjetsel,    "JetEta_eee_afterjetsel",    datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  
  MyhistoManager.CreateHisto(JetEta_mumumu_afterbjetsel, "JetEta_mumumu_afterbjetsel", datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  MyhistoManager.CreateHisto(JetEta_mumue_afterbjetsel,  "JetEta_mumue_afterbjetsel",  datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  MyhistoManager.CreateHisto(JetEta_eemu_afterbjetsel,   "JetEta_eemu_afterbjetsel",   datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  MyhistoManager.CreateHisto(JetEta_eee_afterbjetsel,    "JetEta_eee_afterbjetsel",    datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  
  MyhistoManager.CreateHisto(JetEta_mumumu_afterbjetveto, "JetEta_mumumu_afterbjetveto", datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  MyhistoManager.CreateHisto(JetEta_mumue_afterbjetveto,  "JetEta_mumue_afterbjetveto",  datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  MyhistoManager.CreateHisto(JetEta_eemu_afterbjetveto,   "JetEta_eemu_afterbjetveto",   datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  MyhistoManager.CreateHisto(JetEta_eee_afterbjetveto,    "JetEta_eee_afterbjetveto",    datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  
 
  
  MyhistoManager.CreateHisto(HT_mumumu_afterleptsel, "HT_mumumu_afterleptsel", datasetName, "HT", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(HT_mumue_afterleptsel,  "HT_mumue_afterleptsel",  datasetName, "HT", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(HT_eemu_afterleptsel,   "HT_eemu_afterleptsel",   datasetName, "HT", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(HT_eee_afterleptsel,    "HT_eee_afterleptsel",    datasetName, "HT", "Entries",350,0., 1000);
  
  
  MyhistoManager.CreateHisto(HT_mumumu_afterjetsel, "HT_mumumu_afterjetsel", datasetName, "HT", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(HT_mumue_afterjetsel,  "HT_mumue_afterjetsel",  datasetName, "HT", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(HT_eemu_afterjetsel,   "HT_eemu_afterjetsel",   datasetName, "HT", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(HT_eee_afterjetsel,    "HT_eee_afterjetsel",    datasetName, "HT", "Entries",350,0., 1000);
  
  MyhistoManager.CreateHisto(HT_mumumu_afterbjetsel, "HT_mumumu_afterbjetsel", datasetName, "HT", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(HT_mumue_afterbjetsel,  "HT_mumue_afterbjetsel",  datasetName, "HT", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(HT_eemu_afterbjetsel,   "HT_eemu_afterbjetsel",   datasetName, "HT", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(HT_eee_afterbjetsel,    "HT_eee_afterbjetsel",    datasetName, "HT", "Entries",350,0., 1000);
  
  MyhistoManager.CreateHisto(HT_mumumu_afterbjetveto, "HT_mumumu_afterbjetveto", datasetName, "HT", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(HT_mumue_afterbjetveto,  "HT_mumue_afterbjetveto",  datasetName, "HT", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(HT_eemu_afterbjetveto,   "HT_eemu_afterbjetveto",   datasetName, "HT", "Entries",350,0., 1000);
  MyhistoManager.CreateHisto(HT_eee_afterbjetveto,    "HT_eee_afterbjetveto",    datasetName, "HT", "Entries",350,0., 1000);
  
  
  
  MyhistoManager.CreateHisto(MET_mumumu_afterleptsel, "MET_mumumu_afterleptsel", datasetName, "MET", "Entries",100,0., 500);
  MyhistoManager.CreateHisto(MET_mumue_afterleptsel,  "MET_mumue_afterleptsel",  datasetName, "MET", "Entries",100,0., 500);
  MyhistoManager.CreateHisto(MET_eemu_afterleptsel,   "MET_eemu_afterleptsel",   datasetName, "MET", "Entries",100,0., 500);
  MyhistoManager.CreateHisto(MET_eee_afterleptsel,    "MET_eee_afterleptsel",    datasetName, "MET", "Entries",100,0., 500);
  
  MyhistoManager.CreateHisto(MET_mumumu_afterleptsel_mWT110, "MET_mumumu_afterleptsel_mWT110", datasetName, "MET", "Entries",100,0., 500);
  MyhistoManager.CreateHisto(MET_mumue_afterleptsel_mWT110,  "MET_mumue_afterleptsel_mWT110",  datasetName, "MET", "Entries",100,0., 500);
  MyhistoManager.CreateHisto(MET_eemu_afterleptsel_mWT110,   "MET_eemu_afterleptsel_mWT110",   datasetName, "MET", "Entries",100,0., 500);
  MyhistoManager.CreateHisto(MET_eee_afterleptsel_mWT110,    "MET_eee_afterleptsel_mWT110",    datasetName, "MET", "Entries",100,0., 500);
  
  
  
  
  MyhistoManager.CreateHisto(MET_mumumu_afterjetsel, "MET_mumumu_afterjetsel", datasetName, "MET", "Entries",100,0., 500);
  MyhistoManager.CreateHisto(MET_mumue_afterjetsel,  "MET_mumue_afterjetsel",  datasetName, "MET", "Entries",100,0., 500);
  MyhistoManager.CreateHisto(MET_eemu_afterjetsel,   "MET_eemu_afterjetsel",   datasetName, "MET", "Entries",100,0., 500);
  MyhistoManager.CreateHisto(MET_eee_afterjetsel,    "MET_eee_afterjetsel",    datasetName, "MET", "Entries",100,0., 500);
  
  MyhistoManager.CreateHisto(MET_mumumu_afterbjetsel, "MET_mumumu_afterbjetsel", datasetName, "MET", "Entries",100,0., 500);
  MyhistoManager.CreateHisto(MET_mumue_afterbjetsel,  "MET_mumue_afterbjetsel",  datasetName, "MET", "Entries",100,0., 500);
  MyhistoManager.CreateHisto(MET_eemu_afterbjetsel,   "MET_eemu_afterbjetsel",   datasetName, "MET", "Entries",100,0., 500);
  MyhistoManager.CreateHisto(MET_eee_afterbjetsel,    "MET_eee_afterbjetsel",    datasetName, "MET", "Entries",100,0., 500);
  
   
  MyhistoManager.CreateHisto(Asym_mumumu_afterbjetsel, "Asym_mumumu_afterbjetsel", datasetName, "Asym", "Entries", 20, -3.2, 3.2);
  MyhistoManager.CreateHisto(Asym_mumue_afterbjetsel , "Asym_mumue_afterbjetsel" , datasetName, "Asym", "Entries", 20, -3.2, 3.2);
  MyhistoManager.CreateHisto(Asym_eemu_afterbjetsel  , "Asym_eemu_afterbjetsel"  , datasetName, "Asym", "Entries", 20, -3.2, 3.2);
  MyhistoManager.CreateHisto(Asym_eee_afterbjetsel   , "Asym_eee_afterbjetsel"   , datasetName, "Asym", "Entries", 20, -3.2, 3.2);
  
  MyhistoManager.CreateHisto(RecoPtZ_mumumu_afterbjetsel, "RecoPtZ_mumumu_afterbjetsel", datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto(RecoPtZ_mumue_afterbjetsel , "RecoPtZ_mumue_afterbjetsel" , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto(RecoPtZ_eemu_afterbjetsel  , "RecoPtZ_eemu_afterbjetsel"  , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto(RecoPtZ_eee_afterbjetsel   , "RecoPtZ_eee_afterbjetsel"   , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  
  MyhistoManager.CreateHisto(RecoPtZ_mumumu_afterbjetveto, "RecoPtZ_mumumu_afterbjetveto", datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto(RecoPtZ_mumue_afterbjetveto , "RecoPtZ_mumue_afterbjetveto" , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto(RecoPtZ_eemu_afterbjetveto  , "RecoPtZ_eemu_afterbjetveto"  , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto(RecoPtZ_eee_afterbjetveto   , "RecoPtZ_eee_afterbjetveto"   , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
   
  MyhistoManager.CreateHisto(RecoPtZ_mumumu_afterleptsel, "RecoPtZ_mumumu_afterleptsel", datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto(RecoPtZ_mumue_afterleptsel , "RecoPtZ_mumue_afterleptsel" , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto(RecoPtZ_eemu_afterleptsel  , "RecoPtZ_eemu_afterleptsel"  , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto(RecoPtZ_eee_afterleptsel   , "RecoPtZ_eee_afterleptsel"   , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
   
  MyhistoManager.CreateHisto(RecoPtZ_mumumu_afterleptsel_nojet, "RecoPtZ_mumumu_afterleptsel_nojet", datasetName, "RecoPtZ", "Entries", 200, 0, 100);
  MyhistoManager.CreateHisto(RecoPtZ_mumue_afterleptsel_nojet , "RecoPtZ_mumue_afterleptsel_nojet" , datasetName, "RecoPtZ", "Entries", 200, 0, 100);
  MyhistoManager.CreateHisto(RecoPtZ_eemu_afterleptsel_nojet  , "RecoPtZ_eemu_afterleptsel_nojet"  , datasetName, "RecoPtZ", "Entries", 200, 0, 100);
  MyhistoManager.CreateHisto(RecoPtZ_eee_afterleptsel_nojet   , "RecoPtZ_eee_afterleptsel_nojet"   , datasetName, "RecoPtZ", "Entries", 200, 0, 100);

  MyhistoManager.CreateHisto( RecoTopMass_mumumu_afterbjetsel, "RecoTopMass_mumumu_afterbjetsel", datasetName, "TopMass", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto( RecoTopMass_mumue_afterbjetsel,  "RecoTopMass_mumue_afterbjetsel" , datasetName, "TopMass", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto( RecoTopMass_eemu_afterbjetsel,   "RecoTopMass_eemu_afterbjetsel"  , datasetName, "TopMass", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto( RecoTopMass_eee_afterbjetsel,    "RecoTopMass_eee_afterbjetsel"   , datasetName, "TopMass", "Entries", 200, 0, 500);
  
  MyhistoManager.CreateHisto( RecoTopMass_mumumu_afterbjetveto, "RecoTopMass_mumumu_afterbjetveto", datasetName, "TopMass", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto( RecoTopMass_mumue_afterbjetveto,  "RecoTopMass_mumue_afterbjetveto" , datasetName, "TopMass", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto( RecoTopMass_eemu_afterbjetveto,   "RecoTopMass_eemu_afterbjetveto"  , datasetName, "TopMass", "Entries", 200, 0, 500);
  MyhistoManager.CreateHisto( RecoTopMass_eee_afterbjetveto,    "RecoTopMass_eee_afterbjetveto"   , datasetName, "TopMass", "Entries", 200, 0, 500);
  
  
  MyhistoManager.CreateHisto(deltaPhilb_mumumu_afterbjetsel ,"deltaPhilb_mumumu_afterbjetsel", datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  MyhistoManager.CreateHisto(deltaPhilb_mumue_afterbjetsel  ,"deltaPhilb_mumue_afterbjetsel",  datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  MyhistoManager.CreateHisto(deltaPhilb_eemu_afterbjetsel   ,"deltaPhilb_eemu_afterbjetsel",   datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  MyhistoManager.CreateHisto(deltaPhilb_eee_afterbjetsel    ,"deltaPhilb_eee_afterbjetsel",    datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  
  MyhistoManager.CreateHisto(deltaPhilj_mumumu_afterbjetveto ,"deltaPhilj_mumumu_afterbjetveto", datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  MyhistoManager.CreateHisto(deltaPhilj_mumue_afterbjetveto  ,"deltaPhilj_mumue_afterbjetveto",  datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  MyhistoManager.CreateHisto(deltaPhilj_eemu_afterbjetveto   ,"deltaPhilj_eemu_afterbjetveto",   datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  MyhistoManager.CreateHisto(deltaPhilj_eee_afterbjetveto    ,"deltaPhilj_eee_afterbjetveto",    datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  
  
  
  
  MyhistoManager.CreateHisto(deltaR_mumumu_afterleptsel, "deltaR_mumumu_afterleptsel", datasetName, "deltaRLept", "Entries", 20, 0, 3.2);
  MyhistoManager.CreateHisto(deltaR_mumue_afterleptsel,  "deltaR_mumue_afterleptsel",  datasetName, "deltaRLept", "Entries", 20, 0, 3.2);
  MyhistoManager.CreateHisto(deltaR_eemu_afterleptsel,   "deltaR_eemu_afterleptsel",   datasetName, "deltaRLept", "Entries", 20, 0, 3.2);
  MyhistoManager.CreateHisto(deltaR_eee_afterleptsel,    "deltaR_eee_afterleptsel",    datasetName, "deltaRLept", "Entries", 20, 0, 3.2);
  
  
  
  MyhistoManager.CreateHisto(deltaRLeptJet_mumumu_afterleptsel_mWT110,"deltaRLeptJet_mumumu_afterleptsel_mWT110", datasetName, "deltaR", "Entries", 20, 0, 3.2);
  MyhistoManager.CreateHisto(deltaRLeptJet_mumue_afterleptsel_mWT110 ,"deltaRLeptJet_mumue_afterleptsel_mWT110",  datasetName, "deltaR", "Entries", 20, 0, 3.2);
  MyhistoManager.CreateHisto(deltaRLeptJet_eemu_afterleptsel_mWT110  ,"deltaRLeptJet_eemu_afterleptsel_mWT110",   datasetName, "deltaR", "Entries", 20, 0, 3.2);
  MyhistoManager.CreateHisto(deltaRLeptJet_eee_afterleptsel_mWT110   ,"deltaRLeptJet_eee_afterleptsel_mWT110",    datasetName, "deltaR", "Entries", 20, 0, 3.2);
  
  MyhistoManager.CreateHisto(deltaRLeptMet_mumumu_afterleptsel_mWT110,"deltaRLeptMet_mumumu_afterleptsel_mWT110", datasetName, "deltaR", "Entries", 20, 0, 3.2);
  MyhistoManager.CreateHisto(deltaRLeptMet_mumue_afterleptsel_mWT110 ,"deltaRLeptMet_mumue_afterleptsel_mWT110",  datasetName, "deltaR", "Entries", 20, 0, 3.2);
  MyhistoManager.CreateHisto(deltaRLeptMet_eemu_afterleptsel_mWT110  ,"deltaRLeptMet_eemu_afterleptsel_mWT110",   datasetName, "deltaR", "Entries", 20, 0, 3.2);
  MyhistoManager.CreateHisto(deltaRLeptMet_eee_afterleptsel_mWT110   ,"deltaRLeptMet_eee_afterleptsel_mWT110",    datasetName, "deltaR", "Entries", 20, 0, 3.2);
  
  
  
  MyhistoManager.CreateHisto(WmissAssing_mumumu_afterleptsel, "WmissAssing_mumumu_afterleptsel", datasetName, "MissAs", "Entries", 3, -0.5, 1.5);
  MyhistoManager.CreateHisto(WmissAssing_mumue_afterleptsel,  "WmissAssing_mumue_afterleptsel",  datasetName, "MissAs", "Entries", 3, -0.5, 1.5);
  MyhistoManager.CreateHisto(WmissAssing_eemu_afterleptsel,   "WmissAssing_eemu_afterleptsel",   datasetName, "MissAs", "Entries", 3, -0.5, 1.5);
  MyhistoManager.CreateHisto(WmissAssing_eee_afterleptsel,    "WmissAssing_eee_afterleptsel",    datasetName, "MissAs", "Entries", 3, -0.5, 1.5);
  
  
  
  MyhistoManager.CreateHisto(mWT_mumumu_afterbjetveto, "mWT_mumumu_afterbjetveto", datasetName, "MWT", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(mWT_mumue_afterbjetveto,  "mWT_mumue_afterbjetveto" , datasetName, "MWT", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(mWT_eemu_afterbjetveto,   "mWT_eemu_afterbjetveto"  , datasetName, "MWT", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(mWT_eee_afterbjetveto,    "mWT_eee_afterbjetveto"   , datasetName, "MWT", "Entries", 100, 0, 200);
  
  MyhistoManager.CreateHisto(mWT_mumumu_afterbjetsel, "mWT_mumumu_afterbjetsel", datasetName, "MWT", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(mWT_mumue_afterbjetsel,  "mWT_mumue_afterbjetsel" , datasetName, "MWT", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(mWT_eemu_afterbjetsel,   "mWT_eemu_afterbjetsel"  , datasetName, "MWT", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(mWT_eee_afterbjetsel,    "mWT_eee_afterbjetsel"   , datasetName, "MWT", "Entries", 100, 0, 200);
  
  MyhistoManager.CreateHisto(mWT_mumumu_afterjetsel, "mWT_mumumu_afterjetsel", datasetName, "MWT", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(mWT_mumue_afterjetsel,  "mWT_mumue_afterjetsel" , datasetName, "MWT", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(mWT_eemu_afterjetsel,   "mWT_eemu_afterjetsel"  , datasetName, "MWT", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(mWT_eee_afterjetsel,    "mWT_eee_afterjetsel"   , datasetName, "MWT", "Entries", 100, 0, 200);
  
  MyhistoManager.CreateHisto(mWT_mumumu_afterleptsel, "mWT_mumumu_afterleptsel", datasetName, "MWT", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(mWT_mumue_afterleptsel,  "mWT_mumue_afterleptsel" , datasetName, "MWT", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(mWT_eemu_afterleptsel,   "mWT_eemu_afterleptsel"  , datasetName, "MWT", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(mWT_eee_afterleptsel,    "mWT_eee_afterleptsel"   , datasetName, "MWT", "Entries", 100, 0, 200);
  
  
  
  MyhistoManager.CreateHisto(Charge_mumumu_afterleptsel, "Charge_mumumu_afterleptsel", datasetName, "Charge", "Entries", 11, -5, 5);
  MyhistoManager.CreateHisto(Charge_mumue_afterleptsel,  "Charge_mumue_afterleptsel",  datasetName, "Charge", "Entries", 11, -5, 5);
  MyhistoManager.CreateHisto(Charge_eemu_afterleptsel,   "Charge_eemu_afterleptsel",   datasetName, "Charge", "Entries", 11, -5, 5);
  MyhistoManager.CreateHisto(Charge_eee_afterleptsel,    "Charge_eee_afterleptsel",    datasetName, "Charge", "Entries", 11, -5, 5);
  
  MyhistoManager.CreateHisto(Charge_mumumu_afterleptsel_mWT110, "Charge_mumumu_afterleptsel_mWT110", datasetName, "Charge", "Entries", 11, -5, 5);
  MyhistoManager.CreateHisto(Charge_mumue_afterleptsel_mWT110,  "Charge_mumue_afterleptsel_mWT110",  datasetName, "Charge", "Entries", 11, -5, 5);
  MyhistoManager.CreateHisto(Charge_eemu_afterleptsel_mWT110,   "Charge_eemu_afterleptsel_mWT110",   datasetName, "Charge", "Entries", 11, -5, 5);
  MyhistoManager.CreateHisto(Charge_eee_afterleptsel_mWT110,    "Charge_eee_afterleptsel_mWT110",    datasetName, "Charge", "Entries", 11, -5, 5);
  
  
  
  
  MyhistoManager.CreateHisto(DijetInvM_mumumu_afterleptsel_inZpeak, "DijetInvM_mumumu_afterleptsel_inZpeak", datasetName, "DiJet", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(DijetInvM_mumue_afterleptsel_inZpeak,  "DijetInvM_mumue_afterleptsel_inZpeak" , datasetName, "DiJet", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(DijetInvM_eemu_afterleptsel_inZpeak,   "DijetInvM_eemu_afterleptsel_inZpeak"  , datasetName, "DiJet", "Entries", 100, 0, 200);
  MyhistoManager.CreateHisto(DijetInvM_eee_afterleptsel_inZpeak,    "DijetInvM_eee_afterleptsel_inZpeak"   , datasetName, "DiJet", "Entries", 100, 0, 200);
  
  
  
  
  //**********************
  //initiate 2D histograms
  //**********************
  
  MyhistoManager.CreateHisto2D(HT_vs_MET_mumumu_afterleptsel, "HT_vs_MET_mumumu_afterleptsel", datasetName, "HT",100,0., 1000., "MET",  100,0., 1000.)  ;
  MyhistoManager.CreateHisto2D(HT_vs_MET_mumue_afterleptsel , "HT_vs_MET_mumue_afterleptsel", datasetName, "HT", 100,0., 1000., "MET",  100,0., 1000.);
  MyhistoManager.CreateHisto2D(HT_vs_MET_eemu_afterleptsel  , "HT_vs_MET_eemu_afterleptsel",  datasetName, "HT", 100,0., 1000., "MET",  100,0., 1000.);
  MyhistoManager.CreateHisto2D(HT_vs_MET_eee_afterleptsel   , "HT_vs_MET_eee_afterleptsel",   datasetName, "HT", 100,0., 1000., "MET",  100,0., 1000.);
  
  MyhistoManager.CreateHisto2D(HT_vs_NJet_mumumu_afterleptsel, "HT_vs_NJet_mumumu_afterleptsel", datasetName, "HT", 100,0., 1000,"Njets", 5,-0.5,4.5);
  MyhistoManager.CreateHisto2D(HT_vs_NJet_mumue_afterleptsel , "HT_vs_NJet_mumue_afterleptsel",  datasetName, "HT", 100,0., 1000,"Njets", 5,-0.5,4.5);
  MyhistoManager.CreateHisto2D(HT_vs_NJet_eemu_afterleptsel  , "HT_vs_NJet_eemu_afterleptsel",   datasetName, "HT", 100,0., 1000,"Njets", 5,-0.5,4.5);
  MyhistoManager.CreateHisto2D(HT_vs_NJet_eee_afterleptsel   , "HT_vs_NJet_eee_afterleptsel",    datasetName, "HT", 100,0., 1000,"Njets", 5,-0.5,4.5);
  
  MyhistoManager.CreateHisto2D(HT_vs_NBJet_mumumu_afterleptsel, "HT_vs_NBJet_mumumu_afterleptsel", datasetName, "HT", 100,0., 1000,"NBjets", 5,-0.5,4.5);
  MyhistoManager.CreateHisto2D(HT_vs_NBJet_mumue_afterleptsel , "HT_vs_NBJet_mumue_afterleptsel",  datasetName, "HT", 100,0., 1000,"NBjets", 5,-0.5,4.5);
  MyhistoManager.CreateHisto2D(HT_vs_NBJet_eemu_afterleptsel  , "HT_vs_NBJet_eemu_afterleptsel",   datasetName, "HT", 100,0., 1000,"NBjets", 5,-0.5,4.5);
  MyhistoManager.CreateHisto2D(HT_vs_NBJet_eee_afterleptsel   , "HT_vs_NBJet_eee_afterleptsel",    datasetName, "HT", 100,0., 1000,"NBjets", 5,-0.5,4.5);
  
  MyhistoManager.CreateHisto2D(HT_vs_LeptPt_mumumu_afterleptsel, "HT_vs_LeptPt_mumumu_afterleptsel", datasetName, "HT", 100,0., 1000, "LeptPt",100,0., 1000.);
  MyhistoManager.CreateHisto2D(HT_vs_LeptPt_mumue_afterleptsel , "HT_vs_LeptPt_mumue_afterleptsel",  datasetName, "HT", 100,0., 1000, "LeptPt",100,0., 1000.);
  MyhistoManager.CreateHisto2D(HT_vs_LeptPt_eemu_afterleptsel  , "HT_vs_LeptPt_eemu_afterleptsel",   datasetName, "HT", 100,0., 1000, "LeptPt",100,0., 1000.);
  MyhistoManager.CreateHisto2D(HT_vs_LeptPt_eee_afterleptsel   , "HT_vs_LeptPt_eee_afterleptsel",    datasetName, "HT", 100,0., 1000, "LeptPt",100,0., 1000.);
  
  MyhistoManager.CreateHisto2D(HT_vs_JetPt_mumumu_afterleptsel, "HT_vs_JetPt_mumumu_afterleptsel", datasetName, "HT", 100,0., 1000, "JetPt",100,0., 1000.);
  MyhistoManager.CreateHisto2D(HT_vs_JetPt_mumue_afterleptsel  , "HT_vs_JetPt_mumue_afterleptsel",   datasetName, "HT", 100,0., 1000, "JetPt",100,0., 1000.);
  MyhistoManager.CreateHisto2D(HT_vs_JetPt_eemu_afterleptsel   , "HT_vs_JetPt_eemu_afterleptsel",    datasetName, "HT", 100,0., 1000, "JetPt",100,0., 1000.);
  MyhistoManager.CreateHisto2D(HT_vs_JetPt_eee_afterleptsel    , "HT_vs_JetPt_eee_afterleptsel",     datasetName, "HT", 100,0., 1000, "JetPt",100,0., 1000.);
  
  
  MyhistoManager.CreateHisto2D(HT_vs_Mll_mumumu_afterleptsel, "HT_vs_Mll_mumumu_afterleptsel",   datasetName, "HT", 100,0., 1000, "Mll",100,0., 1000.);
  MyhistoManager.CreateHisto2D(HT_vs_Mll_mumue_afterleptsel  , "HT_vs_Mll_mumue_afterleptsel",   datasetName, "HT", 100,0., 1000, "Mll",100,0., 1000.);
  MyhistoManager.CreateHisto2D(HT_vs_Mll_eemu_afterleptsel   , "HT_vs_Mll_eemu_afterleptsel",    datasetName, "HT", 100,0., 1000, "Mll",100,0., 1000.);
  MyhistoManager.CreateHisto2D(HT_vs_Mll_eee_afterleptsel    , "HT_vs_Mll_eee_afterleptsel",     datasetName, "HT", 100,0., 1000, "Mll",100,0., 1000.);
  
  
  MyhistoManager.CreateHisto2D(InvM_ll_vs_mWT_mumumu_afterleptsel, "InvM_ll_vs_mWT_mumumu_afterleptsel", datasetName, "Mll", 20,0., 200, "mWT", 20,0., 200.);
  MyhistoManager.CreateHisto2D(InvM_ll_vs_mWT_mumue_afterleptsel,  "InvM_ll_vs_mWT_mumue_afterleptsel" , datasetName, "Mll", 20,0., 200, "mWT", 20,0., 200.);
  MyhistoManager.CreateHisto2D(InvM_ll_vs_mWT_eemu_afterleptsel,   "InvM_ll_vs_mWT_eemu_afterleptsel"  , datasetName, "Mll", 20,0., 200, "mWT", 20,0., 200.);
  MyhistoManager.CreateHisto2D(InvM_ll_vs_mWT_eee_afterleptsel,    "InvM_ll_vs_mWT_eee_afterleptsel"   , datasetName, "Mll", 20,0., 200, "mWT", 20,0., 200.);
 
  
  TString treename = "Ttree_"+datasetName;
  
  TheTree = new TTree(treename.Data(),treename.Data());
  TheTree->Branch("tree_topMass",     &tree_topMass,     "tree_topMass/F"   );
  TheTree->Branch("tree_totMass",     &tree_totMass,     "tree_totMass/F"   );
  TheTree->Branch("tree_deltaPhilb",  &tree_deltaPhilb,  "tree_deltaPhilb/F");
  TheTree->Branch("tree_deltaRlb",    &tree_deltaRlb,    "tree_deltaRlb/F"  );
  TheTree->Branch("tree_deltaRTopZ",  &tree_deltaRTopZ,  "tree_deltaRTopZ/F");
  TheTree->Branch("tree_asym",        &tree_asym,        "tree_asym/F"      );
  TheTree->Branch("tree_Zpt",         &tree_Zpt,         "tree_Zpt/F"       );
  TheTree->Branch("tree_ZEta",        &tree_ZEta,        "tree_ZEta/F"      );
  TheTree->Branch("tree_topPt",       &tree_topPt,       "tree_topPt/F"     );
  TheTree->Branch("tree_topEta",      &tree_topEta,      "tree_topEta/F"    );
  TheTree->Branch("tree_NJets",       &tree_NJets,       "tree_NJets/F"     );
  TheTree->Branch("tree_NBJets",      &tree_NBJets,      "tree_NBJets/F"    );
  TheTree->Branch("tree_deltaRZl",    &tree_deltaRZl,    "tree_deltaRZl/F"     );
  TheTree->Branch("tree_deltaPhiZmet",&tree_deltaPhiZmet,"tree_deltaPhiZmet/F" );
  TheTree->Branch("tree_btagDiscri",  &tree_btagDiscri,  "tree_btagDiscri/F"   );
  
  TheTree->Branch("tree_totPt",      &tree_totPt,      "tree_totPt/F"   );
  TheTree->Branch("tree_totEta",     &tree_totEta,     "tree_totEta/F"   );
  
  
  TheTree->Branch("tree_leptWPt",        &tree_leptWPt        , "tree_leptWPt/F"        );
  TheTree->Branch("tree_leptWEta",       &tree_leptWEta       , "tree_leptWEta/F"       );
  TheTree->Branch("tree_leadJetPt",      &tree_leadJetPt      , "tree_leadJetPt/F"      ); 
  TheTree->Branch("tree_leadJetEta",     &tree_leadJetEta     , "tree_leadJetEta/F"     );
  TheTree->Branch("tree_deltaRZleptW",   &tree_deltaRZleptW   , "tree_deltaRZleptW/F"   );
  TheTree->Branch("tree_deltaPhiZleptW", &tree_deltaPhiZleptW , "tree_deltaPhiZleptW/F" );
  
  
  
  TheTree->Branch("tree_EvtWeight",   &tree_EvtWeight,   "tree_EvtWeight/F" );
  TheTree->Branch("tree_SampleType",  &tree_SampleType,  "tree_SampleType/I");
  TheTree->Branch("tree_Channel",     &tree_Channel,     "tree_Channel/I"   );
  
  
  
  
  if (IReweight ) {

   
    string mcfile;
    if( datasetName == "FCNCkut"  || datasetName == "FCNCkct"  || datasetName == "FCNCzut"  || datasetName == "FCNCzct") // FastSim, in-time PU only
      mcfile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PUMC_InTime_Fall11.root");
     else
       mcfile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PU3DMC_Fall11_JLA.root");
    fexists(mcfile, true);
    
    string datafile;
    if( datasetName == "FCNCkut"|| datasetName == "FCNCkct"  || datasetName == "FCNCzut"  || datasetName == "FCNCzct" ) {
      if( !IReweight_puDown && !IReweight_puUp ) datafile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PUData2011_observed_68mb.root");
      if( IReweight_puDown ) datafile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PUData2011_observed_64.6mb.root");
      if( IReweight_puUp ) datafile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PUData2011_observed_71.4mb.root");
    }
    else {
      if( !IReweight_puDown && !IReweight_puUp ) datafile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PUData2011_68mb.root");
      if( IReweight_puDown ) datafile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PUData2011_64.6mb.root");
      if( IReweight_puUp ) datafile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PUData2011_71.4mb.root");
    }
    fexists(datafile, true);
 
 

 
 
    
    LumiWeights = new reweight::LumiReWeighting(mcfile, datafile, "histoMCPU", "pileup" );
    LumiWeights->weight3D_init( 1. );
    
  }
  
  JEC_L2L3Residuals.LoadCorrections();

  if(doPDF) pdf.Initialize();
  //************************************
  
  //cout << "618 " <<  endl;
  //--------------------------------------//
  //   Output file 	
  //--------------------------------------//
  //retrieve info from the input:
  TNamed *out = (TNamed *) fInput->FindObject("PROOF_OUTPUTFILE");
  //this file will be THE file which will contains all the histograms
  fProofFile = new TProofOutputFile(out->GetTitle());
  // Open the file
  TDirectory *savedir = gDirectory;
  fFile = fProofFile->OpenFile("UPDATE");
  if (fFile && fFile->IsZombie()) SafeDelete(fFile);
  savedir->cd();
  
  
  //this file is very important !!!
  fFile->Write();
  //It is required to add in fOutput the histos you want to feedback
  //fOutput->Add(fHist);
  fOutput->Add(fFile);
  cout << "end SlaveBegin " << endl;
}

//_____________________________________________________________________________
Bool_t ProofSelectorMyCutFlow::Process(Long64_t entry)
{
  
  //---------------------------------------------------//
  // Main event loop: get entry
  //---------------------------------------------------//
  fChain->GetTree()->GetEntry(entry); 
  branch->GetEntry(entry);
  
  IPHCTree::NTTransient::InitializeAfterReading(event);

  bool isData = false;    
  if(datasetName=="DataEG" || datasetName=="DataMu" || 
     datasetName=="DataMuEG" || datasetName=="DataEGMu" ||
     datasetName=="MET1" || datasetName=="MET2") isData = true;
  if(isData) JEC_L2L3Residuals.ApplyCorrections(event); // n'appliquer la correction que pour les donnees
    
  // cout << "line 646 " << endl;
  
  //---------------------------------------------------//
  //         Doing the analysis event by event
  //---------------------------------------------------//
  int debugcc=1000;
  int maxdebugcc=10;
  //cout<<"Entry "<<entry<<endl;
  sel.LoadEvent(event);
  // cout << "Evt loaded " << endl;
  
  //Collection of selected objects
  vector<NTVertex>   selVertices  = sel.GetSelectedVertex();
  // cout << "Vertices loaded " << endl;
  // cout<<sel.GetPointer2Electrons()<<endl;
  vector<NTElectron> selElectrons; 
  // cout << "Electrons loaded " << endl;
  vector<NTMuon>     selMuons    ;
  // cout << "Muons loaded " << endl;
  
  MyhistoManager.FillHisto(Ntrilept_mumumu, "Ntrilept_mumumu", selMuons.size(), datasetName, 1, 1);
  
  vector<NTElectron> selElectronsNonIso = sel.GetSelectedElectronsNoIso();
  // cout << "No iso electrons loaded " << endl;
  vector<NTMuon>     selMuonsNonIso     = sel.GetSelectedMuonsNoIso();
  // cout << "No iso muons loaded " << endl;
  
 
  
  selElectrons = sel.GetSelectedElectrons();
  selMuons     = sel.GetSelectedMuons();
  
  MyhistoManager.FillHisto(Ntrileptnoniso_mumumu, "Ntrileptnoniso_mumumu", selMuonsNonIso.size(), datasetName, 1, 1);
  
  vector<NTElectron> ZeeCand; 
  vector<NTMuon>     ZmumuCand; 
  
  vector<NTElectron> WeCand; 
  vector<NTMuon>     WmuCand; 
 
  //NTMET met			   = sel.GetMET(); 
  NTMET met			   = sel.GetSelectedMET(applyJES, scale, applyJER, ResFactor);  
  // cout << "MET loaded " << endl;

  tree_topMass	  = -100000;
  tree_totMass	  = -100000;
  tree_deltaPhilb = -100000;
  tree_deltaRlb   = -100000;
  tree_deltaRTopZ = -100000;
  tree_asym	  = -100000;
  tree_Zpt	  = -100000;
  tree_ZEta	  = -100000;
  tree_topPt	  = -100000;
  tree_topEta	  = -100000;
  tree_SampleType = -100000;
  tree_Channel    = -100000;
  tree_EvtWeight  = -100000;

  
  
 
  double Dweight[101];
  for(int k1=0; k1<101; k1++) {
    Dweight[k1] = 0.;
  }   
  
   
  double weightITypeMC_save = Luminosity*dataset->Xsection()/dataset->getNSkimmedEvent();
  double weightITypeMC=0;
   
  //*****************************************************************
  // Loop over the datasets (lepton pairs which triggered the events)
  //*****************************************************************
  
  double EventYieldWeightError = 0;
  
  // cout << "line 686 " << endl;
  for (int IChannel=0; IChannel<4; IChannel++) {
    //      for (int IChannel=1; IChannel<2; IChannel++) {
    string ChannelName = "";
    if      (IChannel==0) ChannelName= "mumumu"; 
    else if (IChannel==1) ChannelName= "mumue" ; 
    else if (IChannel==2) ChannelName= "eemu"  ;  
    else if (IChannel==3) ChannelName= "eee"   ; 
    
    if (IChannel==0 && (datasetName=="DataEG" || datasetName=="DataMuEG")) continue;
    if (IChannel==1 && (datasetName=="DataEG" || datasetName=="DataMu")  ) continue;
    if (IChannel==2 && (datasetName=="DataMu" || datasetName=="DataEG")  ) continue;
    if (IChannel==3 && (datasetName=="DataMu" || datasetName=="DataMuEG")) continue;
   
    
    
    //*****************************************************************
    // calcul the MC weights
    //*****************************************************************    
    if ( datasetName!="DataEG"   && datasetName!="DataMu" && 
	 datasetName!="DataMuEG" && datasetName!="DataEGMu" 
	 && datasetName!="MET1"  && datasetName!="MET2") {
      
      
      if(IReweight ){
	
        //if( datasets[d].Name() == "FCNCkut" ) // FastSim, in-time PU only
        //  weightITypeMC = weightITypeMC_save*LumiWeights->ITweight(event->pileup.intime_npu);
	//else
	  weightITypeMC = weightITypeMC_save*LumiWeights->weight3D(event->pileup.before_npu, event->pileup.intime_npu, event->pileup.after_npu);

      }
      else weightITypeMC = weightITypeMC_save;
    }
    else weightITypeMC = 1;

    //double reweightPrivateProd[30] = {  1.08203 , 0.965983 , 0.948305 , 0.9541 , 0.979744 , 1.00502 , 1.014 , 1.01024 , 1.01328 , 1.01665 , 1.01683 , 1.01021 , 1.01111 , 1.00089 , 0.999908 , 0.995804 , 0.996283 , 0.99994 , 0.988876 , 1.00544 , 1.00754 , 1.0149 , 1.02353 , 1.06242 , 1.0624 , 1.07799 , 1.09178 , 1.10063 , 1.0843 , 1.10092 };
    //double reweightPrivateProd[30] = {  1.13238 , 0.949412 , 0.909075 , 0.914348 , 0.959149 , 0.990995 , 0.99766 , 1.00601 , 1.02268 , 1.03553 , 1.04722 , 1.05685 , 1.05455 , 1.03764 , 1.04671 , 1.04368 , 1.05231 , 1.05121 , 1.03173 , 1.07179 , 1.06419 , 1.07846 , 1.04479 , 1.14572 , 1.15027 , 1.19357 , 1.25374 , 1.31405 , 1.24932 , 1.30296 } ;

    
    
    
    
    //if(selVertices.size() == 0) cout << "nooooooo vertex " << endl;
    //cout << "weightITypeMC 1 " << weightITypeMC  << " selVertices.size " << selVertices.size()  << endl;
    
    //if(WZprivate){ 
    
      std::vector<double> reweightPrivateProd = GetNvertexWeight(datasetName);

      if(selVertices.size() <=31 ) weightITypeMC = weightITypeMC*reweightPrivateProd[ selVertices.size()];
      else weightITypeMC = weightITypeMC*reweightPrivateProd[31];
    
    
    //}
    //cout << "weightITypeMC 2 " << weightITypeMC  <<  endl;
    
    
    //*****************************************************************
    // determine top decay channel
    //*****************************************************************    
    
    bool IsTTbarDilept = false;
    bool IsSignal = false;
    double WeightForBranchingRatio = 1.;
    bool IsLJ = false;
    
    
    //cout << "line 737 " << endl;
    if ( datasetName=="TTbar" ) {
      
      if ( IChannel==0) { // "mumu" 
	if ( event->mc.TMEME==20 || event->mc.TMEME==11010 || event->mc.TMEME==22000 )    IsTTbarDilept = true;
	if ( !(event->mc.TMEME==20 || event->mc.TMEME==11010 || event->mc.TMEME==22000) ) IsTTbarDilept = false;
      }      
      else if ( IChannel==1) {  // "ee" 
	if ( event->mc.TMEME==2 || event->mc.TMEME==10101 || event->mc.TMEME==20200 )     IsTTbarDilept = true;
	if ( !(event->mc.TMEME==2 || event->mc.TMEME==10101 || event->mc.TMEME==20200) )  IsTTbarDilept = false;
      }      
      else if ( IChannel==2) { // "emu" 
	if ( event->mc.TMEME==11 || event->mc.TMEME==21100 || event->mc.TMEME==11001 || event->mc.TMEME==10110 )     IsTTbarDilept = true;
	if ( !(event->mc.TMEME==11 || event->mc.TMEME==21100 || event->mc.TMEME==11001 || event->mc.TMEME==10110) )  IsTTbarDilept = false;
      }      
      if ( !IsTTbarDilept && event->mc.TMEME!=0 )     IsLJ = true;
      if ( IsTTbarDilept ) {
        WeightForBranchingRatio = SF_BranchingRatio_ll;
      } else {
        if ( event->mc.TMEME==0 ){
	  WeightForBranchingRatio = SF_BranchingRatio_had;
	} else{
	  WeightForBranchingRatio = SF_BranchingRatio_lj;
	} 
      } 
      
    }
    
    
    //cout << "line 766 " << endl;
    //*****************************************************************
    // fill cutflow before any selection
    //*****************************************************************    
    
    
    
    if ( datasetName=="TTbar" ) { 
      //cout << " Dweight[ITypeMC] 689 " <<weightITypeMC_save << "  ITypeMC " << ITypeMC << endl;
      ITypeMC = 1; 
      Dweight[ITypeMC]= weightITypeMC * WeightForBranchingRatio;
      EventYieldWeightError = Dweight[ITypeMC]*Dweight[ITypeMC];


      
      //TabFlow1[IChannel][ITypeMC][0]+=Dweight[ITypeMC];
      //TabFlow2[IChannel][ITypeMC][0]+=Dweight[ITypeMC]*Dweight[ITypeMC];
      //cout << " weightITypeMC " << weightITypeMC << " Dweight[ITypeMC] 693 " << Dweight[ITypeMC]<< "  ITypeMC " << ITypeMC << endl;
    }
    
    else if ( datasetName=="Zjets" || datasetName=="DYToLL_M10-50"
	      ) { 
      ITypeMC = 2; IsSignal = false; Dweight[ITypeMC]= weightITypeMC;  
      
      EventYieldWeightError = Dweight[ITypeMC]*Dweight[ITypeMC];
      
      

      
      //TabFlow1[IChannel][ITypeMC][0]+=Dweight[ITypeMC];
      //TabFlow2[IChannel][ITypeMC][0]+=Dweight[ITypeMC]*Dweight[ITypeMC];
    }
    else if (  datasetName=="Wjets"  ) { 
      ITypeMC = 3; IsSignal = false;  Dweight[ITypeMC]= weightITypeMC; 
      EventYieldWeightError = Dweight[ITypeMC]*Dweight[ITypeMC];
      
      
 
      
      
      //TabFlow1[IChannel][ITypeMC][0]+=Dweight[ITypeMC];
      //TabFlow2[IChannel][ITypeMC][0]+=Dweight[ITypeMC]*Dweight[ITypeMC];
    }
    else if ( datasetName=="SingleToptW" || datasetName=="TtW" || datasetName=="TbartW"
	      || datasetName=="TtWScaleUp" || datasetName=="TtWScaleDown"
	      || datasetName=="TbartWScaleUp" || datasetName=="TbartWScaleDown"
	      || datasetName=="TZq"
	      || datasetName=="TZq_matchup"|| datasetName=="TZq_matchdown"
	      || datasetName=="TZq_scaleup"|| datasetName=="TZq_scaledown"
	      || datasetName=="TZq_topup"|| datasetName=="TZq_topdown") { 
	      
      ITypeMC = 4; IsSignal = false;  Dweight[ITypeMC]= weightITypeMC; 
      EventYieldWeightError = Dweight[ITypeMC]*Dweight[ITypeMC];
      
      
 
      
      //TabFlow1[IChannel][ITypeMC][0]+=Dweight[ITypeMC];
      //TabFlow2[IChannel][ITypeMC][0]+=Dweight[ITypeMC]*Dweight[ITypeMC];
    }
    else if ( datasetName=="WZ" || datasetName=="WW" || datasetName=="ZZ"  || datasetName=="VV"
     || datasetName=="WZprivate" || datasetName=="WZprivate_scaleup"|| datasetName=="WZprivate_scaledown"
     || datasetName=="WZprivate_matchup" || datasetName=="WZprivate_matchdown" 
    
    ) { 
      ITypeMC = 5; IsSignal = false;  Dweight[ITypeMC]= weightITypeMC; 
      EventYieldWeightError = Dweight[ITypeMC]*Dweight[ITypeMC];

      //TabFlow1[IChannel][ITypeMC][0]+=Dweight[ITypeMC];
      //TabFlow2[IChannel][ITypeMC][0]+=Dweight[ITypeMC]*Dweight[ITypeMC];
     
    }  else if ( 
         datasetName=="FCNCkut" || datasetName=="FCNCkct" 
       ||datasetName=="FCNCxut" || datasetName=="FCNCxct" 
       ||datasetName=="FCNCzut" || datasetName=="FCNCzct" 
       
       ||datasetName=="FCNCkutFullSim" || datasetName=="FCNCkctFullSim" 
       ||datasetName=="FCNCxutFullSim" || datasetName=="FCNCxctFullSim" 
       ||datasetName=="FCNCzutFullSim" || datasetName=="FCNCzctFullSim" 
       
       ||datasetName=="FCNCkut_matchup" || datasetName=="FCNCkct_matchup" 
       ||datasetName=="FCNCxut_matchup" || datasetName=="FCNCxct_matchup" 
       ||datasetName=="FCNCzut_matchup" || datasetName=="FCNCzct_matchup" 
       
       ||datasetName=="FCNCkut_matchdown" || datasetName=="FCNCkct_matchdown" 
       ||datasetName=="FCNCxut_matchdown" || datasetName=="FCNCxct_matchdown" 
       ||datasetName=="FCNCzut_matchdown" || datasetName=="FCNCzct_matchdown" 
       
       ||datasetName=="FCNCkut_scaleup" || datasetName=="FCNCkct_scaleup" 
       ||datasetName=="FCNCxut_scaleup" || datasetName=="FCNCxct_scaleup" 
       ||datasetName=="FCNCzut_scaleup" || datasetName=="FCNCzct_scaleup" 
       
       ||datasetName=="FCNCkut_scaledown" || datasetName=="FCNCkct_scaledown" 
       ||datasetName=="FCNCxut_scaledown" || datasetName=="FCNCxct_scaledown" 
       ||datasetName=="FCNCzut_scaledown" || datasetName=="FCNCzct_scaledown" 
       
       ||datasetName=="FCNCkut_topup" || datasetName=="FCNCkct_topup" 
       ||datasetName=="FCNCxut_topup" || datasetName=="FCNCxct_topup" 
       ||datasetName=="FCNCzut_topup" || datasetName=="FCNCzct_topup" 
       
       ||datasetName=="FCNCkut_topdown" || datasetName=="FCNCkct_topdown" 
       ||datasetName=="FCNCxut_topdown" || datasetName=="FCNCxct_topdown" 
       ||datasetName=="FCNCzut_topdown" || datasetName=="FCNCzct_topdown" 
       
       ) { 
      ITypeMC = 6; IsSignal = false;  Dweight[ITypeMC]= weightITypeMC; 
      EventYieldWeightError = Dweight[ITypeMC]*Dweight[ITypeMC];

      //TabFlow1[IChannel][ITypeMC][0]+=Dweight[ITypeMC];
      //TabFlow2[IChannel][ITypeMC][0]+=Dweight[ITypeMC]*Dweight[ITypeMC];
     
    }
    
    //cout << "834 " << endl;
    if ( datasetName=="DataEG" || datasetName=="DataMu" || 
	 datasetName=="DataMuEG" || datasetName=="DataEGMu" ||
	 datasetName=="MET1" || datasetName=="MET2") { 
      ITypeMC = 100;  Dweight[ITypeMC]= weightITypeMC; 
      EventYieldWeightError = Dweight[ITypeMC]*Dweight[ITypeMC];   
      
      
    }
   

    if(IChannel == 0) MyhistoManager.FillHisto(PU_before_mumumu, "PU_before_mumumu", event->pileup.before_npu, datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 1) MyhistoManager.FillHisto(PU_before_mumue, "PU_before_mumue", event->pileup.before_npu, datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 2) MyhistoManager.FillHisto(PU_before_eemu, "PU_before_eemu", event->pileup.before_npu, datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 3) MyhistoManager.FillHisto(PU_before_eee, "PU_before_eee", event->pileup.before_npu, datasetName, IsSignal, Dweight[ITypeMC]);

    if(IChannel == 0) MyhistoManager.FillHisto(PU_intime_mumumu, "PU_intime_mumumu", event->pileup.intime_npu, datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 1) MyhistoManager.FillHisto(PU_intime_mumue, "PU_intime_mumue", event->pileup.intime_npu, datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 2) MyhistoManager.FillHisto(PU_intime_eemu, "PU_intime_eemu", event->pileup.intime_npu, datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 3) MyhistoManager.FillHisto(PU_intime_eee, "PU_intime_eee", event->pileup.intime_npu, datasetName, IsSignal, Dweight[ITypeMC]);

    if(IChannel == 0) MyhistoManager.FillHisto(PU_after_mumumu, "PU_after_mumumu", event->pileup.after_npu, datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 1) MyhistoManager.FillHisto(PU_after_mumue, "PU_after_mumue", event->pileup.after_npu, datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 2) MyhistoManager.FillHisto(PU_after_eemu, "PU_after_eemu", event->pileup.after_npu, datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 3) MyhistoManager.FillHisto(PU_after_eee, "PU_after_eee", event->pileup.after_npu, datasetName, IsSignal, Dweight[ITypeMC]);

    
    
    
    //if(IChannel == 0) MyhistoManager.FillHisto(NVtx_mumumu, "NVtx_mumumu", selVertices.size(), datasetName, IsSignal, Dweight[ITypeMC]);
    //if(IChannel == 1) MyhistoManager.FillHisto(NVtx_mumue,  "NVtx_mumue",  selVertices.size(), datasetName, IsSignal, Dweight[ITypeMC]);
    //if(IChannel == 2) MyhistoManager.FillHisto(NVtx_eemu,   "NVtx_eemu",   selVertices.size(), datasetName, IsSignal, Dweight[ITypeMC]);
    //if(IChannel == 3) MyhistoManager.FillHisto(NVtx_eee,    "NVtx_eee",    selVertices.size(), datasetName, IsSignal, Dweight[ITypeMC]);

    if(IChannel == 0) MyhistoManager.FillHisto(NVtx_mumumu, "NVtx_mumumu", selVertices.size(), datasetName, IsSignal, 1);
    if(IChannel == 1) MyhistoManager.FillHisto(NVtx_mumue,  "NVtx_mumue",  selVertices.size(), datasetName, IsSignal, 1);
    if(IChannel == 2) MyhistoManager.FillHisto(NVtx_eemu,   "NVtx_eemu",   selVertices.size(), datasetName, IsSignal, 1);
    if(IChannel == 3) MyhistoManager.FillHisto(NVtx_eee,    "NVtx_eee",    selVertices.size(), datasetName, IsSignal, 1);
    
    
    
    MyhistoManager.FillHisto(Nvertex, "Nvertex", selVertices.size(), datasetName, IsSignal, Dweight[ITypeMC]);
 
 

    //*****************************************************************
    // pass trigger selection
    //*****************************************************************   
    
    bool passtrigger = false;
    
    if(ChannelName == "mumumu" ) passtrigger = sel.passTriggerSelection ( dataset, "mumu");
    if(ChannelName == "mumue"  ) passtrigger = sel.passTriggerSelection ( dataset, "emu" );
    if(ChannelName == "eemu"   ) passtrigger = sel.passTriggerSelection ( dataset, "emu" );
    if(ChannelName == "eee"    ) passtrigger = sel.passTriggerSelection ( dataset, "ee"  );

    
    if (   passtrigger   ) {
      //cout << "line 859, pass trigger selection " << endl;
      if(IChannel == 0) MyhistoManager.FillHisto(CutFlow_mumumu, "CutFlow_mumumu", 0, datasetName, IsSignal, Dweight[ITypeMC]);
      if(IChannel == 1) MyhistoManager.FillHisto(CutFlow_mumue,  "CutFlow_mumue" , 0, datasetName, IsSignal, Dweight[ITypeMC]);
      if(IChannel == 2) MyhistoManager.FillHisto(CutFlow_eemu,	 "CutFlow_eemu"  , 0, datasetName, IsSignal, Dweight[ITypeMC]);
      if(IChannel == 3) MyhistoManager.FillHisto(CutFlow_eee,    "CutFlow_eee"   , 0, datasetName, IsSignal, Dweight[ITypeMC]);
      if(IChannel == 0) MyhistoManager.FillHisto(ErrCutFlow_mumumu,   "ErrCutFlow_mumumu" , 0, datasetName, IsSignal, EventYieldWeightError);
      if(IChannel == 1) MyhistoManager.FillHisto(ErrCutFlow_mumue,    "ErrCutFlow_mumue"  , 0, datasetName, IsSignal, EventYieldWeightError);
      if(IChannel == 2) MyhistoManager.FillHisto(ErrCutFlow_eemu,     "ErrCutFlow_eemu"   , 0, datasetName, IsSignal, EventYieldWeightError);
      if(IChannel == 3) MyhistoManager.FillHisto(ErrCutFlow_eee,      "ErrCutFlow_eee"    , 0, datasetName, IsSignal, EventYieldWeightError);
      
      if(IChannel == 0) MyhistoManager.FillHisto(NVtx_mumumu_aftertrigsel, "NVtx_mumumu_aftertrigsel", selVertices.size(), datasetName, IsSignal, Dweight[ITypeMC]);
      if(IChannel == 1) MyhistoManager.FillHisto(NVtx_mumue_aftertrigsel, "NVtx_mumue_aftertrigsel", selVertices.size(), datasetName, IsSignal, Dweight[ITypeMC]);
      if(IChannel == 2) MyhistoManager.FillHisto(NVtx_eemu_aftertrigsel, "NVtx_eemu_aftertrigsel", selVertices.size(), datasetName, IsSignal, Dweight[ITypeMC]);
      if(IChannel == 3) MyhistoManager.FillHisto(NVtx_eee_aftertrigsel, "NVtx_eee_aftertrigsel", selVertices.size(), datasetName, IsSignal, Dweight[ITypeMC]);


      
      
      //cout << "lepton size " << selMuons.size() << endl;
      
      
       //*****************************************************************
      // select Z->ee candidate
      //*****************************************************************    
     
      int leptonFlavor = 0;
      int wcharge      = 0;
      ZeeCand.clear(); 
      ZmumuCand.clear(); 
  
      WeCand.clear(); 
      WmuCand.clear(); 

      if(selElectrons.size() >=2 ) {
        int theel1 = -1;
	int theel2 = -1;
	double mInv = 1000000;
	for(unsigned int iel1 = 0; iel1 < selElectrons.size(); iel1++){
	  for(unsigned int iel2 = 0; iel2 < selElectrons.size(); iel2++){
	    if(iel1 == iel2) continue;
	    if(selElectrons[iel1].charge == selElectrons[iel2].charge) continue;
	    TLorentzVector zeecand = selElectrons[iel1].p4 + selElectrons[iel2].p4;
	    if( fabs(zeecand.M() - 91) < fabs(mInv-91) ){
	    //if( fabs(zeecand.M() - 200) < fabs(mInv-200) ){
	      theel1 = iel1;
	      theel2 = iel2;
	      mInv = zeecand.M();
	    }
	  }
	}
	
	if(theel1>=0 && theel2>=0){ //JLA
	  ZeeCand.push_back(selElectrons[theel1]);
	  ZeeCand.push_back(selElectrons[theel2]);
	  double invM = (ZeeCand[0].p4+ZeeCand[1].p4).M();
	  //cout << "lepton origin of Zee cand " << selElectrons[theel1].LeptonOrigin << "  " << selElectrons[theel2].LeptonOrigin << "  invmass " <<  invM << endl;
	}
      }
         
      //*****************************************************************
      // select W->enu candidate
      //*****************************************************************    
      //cout << "get lepton cand W->enu " << endl;
      if(!useNonIsoWcand){
        for(unsigned int iel1 = 0; iel1 < selElectrons.size(); iel1++){
	  bool matchElec=false;
          for(unsigned int iel2 = 0; iel2 < ZeeCand.size(); iel2++){
             
	     if(fabs(selElectrons[iel1].p4.Pt() - ZeeCand[iel2].p4.Pt()) <  0.0001)  matchElec=true;
	   }
	  if(!matchElec){
	    WeCand.push_back(selElectrons[iel1]);
	    wcharge = selElectrons[iel1].charge;
	    if(selElectrons[iel1].LeptonOrigin == 10) leptonFlavor = 1;
	  }
        }
      }else{
        for(unsigned int iel1 = 0; iel1 < selElectronsNonIso.size(); iel1++){
	  bool matchElec=false;
          for(unsigned int iel2 = 0; iel2 < ZeeCand.size(); iel2++){
	     if(fabs(selElectronsNonIso[iel1].p4.Pt() - ZeeCand[iel2].p4.Pt()) <  0.0001)  matchElec=true;
	     else if(selElectrons[iel1].LeptonOrigin == 10) leptonFlavor = 1;
          }
	  if(!matchElec && selElectronsNonIso[iel1].RelIso03PF() > looseIso){ 
	    WeCand.push_back(selElectronsNonIso[iel1]);
	    wcharge = selElectronsNonIso[iel1].charge; 
	    if(selElectronsNonIso[iel1].LeptonOrigin == 10) leptonFlavor = 1;
	    //cout << "    lepton origin of We cand " << selElectronsNonIso[iel1].LeptonOrigin << endl;
	    //if(selElectronsNonIso[iel1].LeptonOrigin == 10) cout << " nofake wenu" << endl;
	  }
        }
      }
      
      
      //*****************************************************************
      // select Z->mumu candidate
      //*****************************************************************  
      if(selMuons.size() >=2 ) {
        int themu1 = -1;
	int themu2 = -1;
	double mInv = 1000000;
	for(unsigned int imu1 = 0; imu1 < selMuons.size(); imu1++){
	  for(unsigned int imu2 = 0; imu2 < selMuons.size(); imu2++){
	    if(imu1 == imu2) continue;
	    if(selMuons[imu1].charge == selMuons[imu2].charge) continue;
	    TLorentzVector zmumucand = selMuons[imu1].p4 + selMuons[imu2].p4;
	    if( fabs(zmumucand.M() - 91) < fabs(mInv-91) ){
	    //if( fabs(zmumucand.M() - 200) < fabs(mInv-200) ){
	      themu1 = imu1;
	      themu2 = imu2;
	      mInv = zmumucand.M();
	    }
	  }
	}
	if(themu1>=0 && themu2>=0){ //JLA  
	  ZmumuCand.push_back(selMuons[themu1]);
	  ZmumuCand.push_back(selMuons[themu2]);
	}
         
      }
      
       
       
      //*****************************************************************
      // select W->munu candidate
      //*****************************************************************    
      
      //cout << "get lepton cand W->munu " << endl;
      if(!useNonIsoWcand){
        //cout << "in sel W " << endl;
        for(unsigned int imu1 = 0; imu1 < selMuons.size(); imu1++){
	  bool matchMuon = false;
          for(unsigned int imu2 = 0; imu2 < ZmumuCand.size(); imu2++){
             
	     if(fabs(selMuons[imu1].p4.Pt() - ZmumuCand[imu2].p4.Pt()) <  0.0001) matchMuon = true;
	     
          } 
	  if(!matchMuon){
	    WmuCand.push_back(selMuons[imu1]);
	    wcharge = selMuons[imu1].charge;
	    if(selMuons[imu1].LeptonOrigin == 10) leptonFlavor = 1;
	  }
        }
      }else{
        for(unsigned int imu1 = 0; imu1 < selMuonsNonIso.size(); imu1++){
	  bool matchMuon = false;
          for(unsigned int imu2 = 0; imu2 < ZmumuCand.size(); imu2++){
             
	     if(fabs(selMuonsNonIso[imu1].p4.Pt() - ZmumuCand[imu2].p4.Pt())  < 0.0001) matchMuon = true;
	 
          }
	  if(!matchMuon && selMuonsNonIso[imu1].RelIso03PF() > looseIso){
	    WmuCand.push_back(selMuonsNonIso[imu1]);
	    wcharge = selMuonsNonIso[imu1].charge;
	    if(selMuonsNonIso[imu1].LeptonOrigin == 10) leptonFlavor = 1;
	  }
        }
      }
      
      //redefine the lepton coll for jet cleaning
      if(useNonIsoWcand){
         
        vector<NTElectron>  tmpElectrons = sel.GetSelectedElectronsNoIso();
        vector<NTMuon>      tmpMuons     = sel.GetSelectedMuonsNoIso();
        
        for(unsigned int iel=0; iel<tmpElectrons.size(); iel++){
          if(tmpElectrons[iel].RelIso03PF() > 0.4) selElectrons.push_back(tmpElectrons[iel]);
        }
      
        for(unsigned int imu=0; imu<tmpMuons.size(); imu++){
          if(tmpMuons[imu].RelIso03PF() > 0.4) selMuons.push_back(tmpMuons[imu]);
        }
      
      }else{
     
        selElectrons = sel.GetSelectedElectrons();
       selMuons     = sel.GetSelectedMuons();
      }
     
      //cout << "1192" << endl;
      //*****************************************************************
      // apply lepton selection
      //*****************************************************************    
      
      if( ZeeCand.size() >= 2 &&  WeCand.size() >= 1 ){  
        /*cout << "found3lept" << endl;
        cout << "Zee size " << ZeeCand.size() << endl;
        cout << "We size  " << WeCand.size() << endl;*/
	
	double detlaR1 = ZeeCand[0].p4.DeltaR(ZeeCand[1].p4);
        //cout << "deltaR 1 " << detlaR1 << endl;
	double detlaR2 = ZeeCand[0].p4.DeltaR(WeCand[0].p4);
        //cout << "deltaR 2 " << detlaR2 << endl;
	double detlaR3 = ZeeCand[1].p4.DeltaR(WeCand[0].p4);
        //cout << "deltaR 3 " << detlaR3 << endl;
        
	
        MyhistoManager.FillHisto(deltaR_eee_afterleptsel, "deltaR_eee_afterleptsel",detlaR1 , datasetName, IsSignal, Dweight[ITypeMC]);
        MyhistoManager.FillHisto(deltaR_eee_afterleptsel, "deltaR_eee_afterleptsel",detlaR2 , datasetName, IsSignal, Dweight[ITypeMC]);
        MyhistoManager.FillHisto(deltaR_eee_afterleptsel, "deltaR_eee_afterleptsel",detlaR3 , datasetName, IsSignal, Dweight[ITypeMC]);
	
	
      }
      
      
      int nbtag = 0;
      vector<NTJet>	 selJetstmp = sel.GetSelectedJets(selMuons, selElectrons, applyJES, scale, applyJER, ResFactor);
      int AlgoBtagtmp = sel.GetbtagAlgo();
      float btagDiscriCuttmp = sel.GetbtagDiscriCut ();
      
      
      for(unsigned int ijet = 0; ijet < selJetstmp.size(); ijet++){
	  
	  if ( AlgoBtagtmp==6 &&  selJetstmp[ijet].bTag["combinedSecondaryVertexBJetTags"]>= btagDiscriCuttmp) nbtag++;
      } 
	
      
      
      
      if(ZeeCand.size() ==2 && WmuCand.size() == 0 && WeCand.size() == 0 && met.p2.Mod() < 20 
        && fabs( (ZeeCand[0].p4+ZeeCand[1].p4).M()-15 ) < 91 ){
        
	vector<NTJet>      selJetsZ = sel.GetSelectedJets(selMuons, selElectrons, applyJES, scale, applyJER, ResFactor);
	if(selJetsZ.size() ==2){
	  TLorentzVector jet1 = selJetsZ[0].p4;
	  TLorentzVector jet2 = selJetsZ[1].p4;
	  double invMDijet = (jet1+jet2).M();
	  double deltaRDiJetsZ =  (jet1+jet2).DeltaR(ZeeCand[0].p4+ZeeCand[1].p4);
	  double detlaRDijets = jet1.DeltaR(jet2);
          if(deltaRDiJetsZ > 3.0 && detlaRDijets < 1.5) MyhistoManager.FillHisto(DijetInvM_eee_afterleptsel_inZpeak, "DijetInvM_eee_afterleptsel_inZpeak", invMDijet, datasetName, IsSignal, Dweight[ITypeMC]);
	}
      }
      
      
      if(ZmumuCand.size() ==2 && WmuCand.size() == 0 && WeCand.size() == 0 && met.p2.Mod() < 20 
        && fabs( (ZmumuCand[0].p4+ZmumuCand[1].p4).M()-15 ) < 91  ){
        
	vector<NTJet>      selJetsZ = sel.GetSelectedJets(selMuons, selElectrons, applyJES, scale, applyJER, ResFactor);
	if(selJetsZ.size() ==2){
	  TLorentzVector jet1 = selJetsZ[0].p4;
	  TLorentzVector jet2 = selJetsZ[1].p4;
	  double invMDijet = (jet1+jet2).M();
	  double deltaRDiJetsZ =  (jet1+jet2).DeltaR(ZmumuCand[0].p4+ZmumuCand[1].p4);
	  double detlaRDijets = jet1.DeltaR(jet2);
          if(deltaRDiJetsZ > 3.0 && detlaRDijets < 1.5) MyhistoManager.FillHisto(DijetInvM_mumumu_afterleptsel_inZpeak, "DijetInvM_mumumu_afterleptsel_inZpeak", invMDijet, datasetName, IsSignal, Dweight[ITypeMC]);
	}
      }
      
      //cout << "ZmumuCand.size() " << ZmumuCand.size() << "  WmuCand.size() " << WmuCand.size() << endl;
      //if( (WmuCand.size()+ZmumuCand.size()+WeCand.size()+ZeeCand.size()) == 3 && met.p2.Mod() > 25 ) {
      //if( (WmuCand.size()+ZmumuCand.size()+WeCand.size()+ZeeCand.size()) == 3 && nbtag == 0) {
      
      
      if( (WmuCand.size()+ZmumuCand.size()+WeCand.size()+ZeeCand.size()) == 3) {
       
	string cand3leptonChannel = "";
	if( ZmumuCand.size() == 2 ) {
	  if(WmuCand.size() == 1 ) cand3leptonChannel = "mumumu";
	  if(WeCand.size()  == 1 ) cand3leptonChannel = "mumue";
	}
	
	if( ZeeCand.size() == 2 ) {
	  if(WmuCand.size() == 1 ) cand3leptonChannel = "eemu";
	  if(WeCand.size()  == 1 ) cand3leptonChannel = "eee";
	}
		
	//cout << "line 921 , " << cand3leptonChannel << endl;
 	if(IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(CutFlow_mumumu, "CutFlow_mumumu", 1, datasetName, IsSignal, Dweight[ITypeMC]);
 	if(IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(CutFlow_mumue,  "CutFlow_mumue" , 1, datasetName, IsSignal, Dweight[ITypeMC]);
 	if(IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(CutFlow_eemu,   "CutFlow_eemu"  , 1, datasetName, IsSignal, Dweight[ITypeMC]);
 	if(IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(CutFlow_eee,    "CutFlow_eee"   , 1, datasetName, IsSignal, Dweight[ITypeMC]);
 	if(IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(ErrCutFlow_mumumu,   "ErrCutFlow_mumumu"  , 1, datasetName, IsSignal, EventYieldWeightError);
 	if(IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(ErrCutFlow_mumue,    "ErrCutFlow_mumue"   , 1, datasetName, IsSignal, EventYieldWeightError);
 	if(IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(ErrCutFlow_eemu,     "ErrCutFlow_eemu"    , 1, datasetName, IsSignal, EventYieldWeightError);
 	if(IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(ErrCutFlow_eee,      "ErrCutFlow_eee"     , 1, datasetName, IsSignal, EventYieldWeightError);
	
	 
	//*****************************************************************
        // select Z candidate
        //*****************************************************************    
	
	
	
	double lept3_Charge = 0;
	
	TLorentzVector dilept;
	TLorentzVector lept1, lept2, lept3;
	
	string theleptpair = "";
	
	if(cand3leptonChannel == "mumue") {
	  lept1 = ZmumuCand[0].p4;
	  lept2 = ZmumuCand[1].p4;
	  lept3 = WeCand[0].p4;
	  lept3_Charge= WeCand[0].charge;
	  dilept = lept1+lept2;
	}
	
	
	if(cand3leptonChannel == "eemu") {
	  lept1 = ZeeCand[0].p4;
	  lept2 = ZeeCand[1].p4;
	  lept3 = WmuCand[0].p4;
	  lept3_Charge= WmuCand[0].charge;
	  dilept = lept1+lept2;
	}
	
	if(cand3leptonChannel == "mumumu") {
	  
	  lept1 = ZmumuCand[0].p4;
	  lept2 = ZmumuCand[1].p4;
	  lept3 = WmuCand[0].p4;
	  lept3_Charge= WmuCand[0].charge;
	  dilept = lept1+lept2;
	  
	}
	
	if(cand3leptonChannel == "eee") {
	  
	  lept1 = ZeeCand[0].p4;
	  lept2 = ZeeCand[1].p4;
	  lept3 = WeCand[0].p4;
	  lept3_Charge= WeCand[0].charge;
	  dilept = lept1+lept2;
	  
	}
	
	
	 
	 
        //*****************************************************************
        // apply lepton scale factors
        //*****************************************************************    
        
    
        double LeptonSF      = 0.;
        double LeptonSFError = 0.;
 
        if(applyLeptonSF && !isData){

	  if(IChannel == 0 && cand3leptonChannel == "mumumu"){
            sumSFlept_mumumu += Dweight[ITypeMC];
            LeptonSF = sel.getLeptonScaleFactor( lept1.Pt(), lept1.Eta(), "mu") 
                      * sel.getLeptonScaleFactor( lept2.Pt(), lept2.Eta(), "mu") 
                      * sel.getLeptonScaleFactor( lept3.Pt(), lept3.Eta(), "mu");
	    if(applyLeptonSFUp)
            LeptonSF = (sel.getLeptonScaleFactor( lept1.Pt(), lept1.Eta(), "mu")+sel.getLeptonScaleFactorError( lept1.Pt(), lept1.Eta(), "mu"))
                      * (sel.getLeptonScaleFactor( lept2.Pt(), lept2.Eta(), "mu")+sel.getLeptonScaleFactorError( lept2.Pt(), lept2.Eta(), "mu"))
                      * (sel.getLeptonScaleFactor( lept3.Pt(), lept3.Eta(), "mu")+sel.getLeptonScaleFactorError( lept3.Pt(), lept3.Eta(), "mu"));		  
	    if(applyLeptonSFDown)
            LeptonSF = (sel.getLeptonScaleFactor( lept1.Pt(), lept1.Eta(), "mu")-sel.getLeptonScaleFactorError( lept1.Pt(), lept1.Eta(), "mu"))
                      * (sel.getLeptonScaleFactor( lept2.Pt(), lept2.Eta(), "mu")-sel.getLeptonScaleFactorError( lept2.Pt(), lept2.Eta(), "mu")) 
                      * (sel.getLeptonScaleFactor( lept3.Pt(), lept3.Eta(), "mu")-sel.getLeptonScaleFactorError( lept3.Pt(), lept3.Eta(), "mu"));		  
 	    nEvents_mumumu += Dweight[ITypeMC]*LeptonSF;
          }
	  
	  if(IChannel == 1 && cand3leptonChannel == "mumue"){
            sumSFlept_mumue += Dweight[ITypeMC];
            LeptonSF = sel.getLeptonScaleFactor( lept1.Pt(), lept1.Eta(), "mu")
                      * sel.getLeptonScaleFactor( lept2.Pt(), lept2.Eta(), "mu")
                      * sel.getLeptonScaleFactor( lept3.Pt(), lept3.Eta(), "e");
	    if(applyLeptonSFUp)
            LeptonSF = (sel.getLeptonScaleFactor( lept1.Pt(), lept1.Eta(), "mu")+sel.getLeptonScaleFactorError( lept1.Pt(), lept1.Eta(), "mu"))
                      * (sel.getLeptonScaleFactor( lept2.Pt(), lept2.Eta(), "mu")+sel.getLeptonScaleFactorError( lept2.Pt(), lept2.Eta(), "mu")) 
                      * (sel.getLeptonScaleFactor( lept3.Pt(), lept3.Eta(), "e")+sel.getLeptonScaleFactorError( lept3.Pt(), lept3.Eta(), "e"));		  
	    if(applyLeptonSFDown)
            LeptonSF = (sel.getLeptonScaleFactor( lept1.Pt(), lept1.Eta(), "mu")-sel.getLeptonScaleFactorError( lept1.Pt(), lept1.Eta(), "mu"))
                      * (sel.getLeptonScaleFactor( lept2.Pt(), lept2.Eta(), "mu")-sel.getLeptonScaleFactorError( lept2.Pt(), lept2.Eta(), "mu")) 
                      * (sel.getLeptonScaleFactor( lept3.Pt(), lept3.Eta(), "e")-sel.getLeptonScaleFactorError( lept3.Pt(), lept3.Eta(), "e"));		  
	    nEvents_mumue += Dweight[ITypeMC]*LeptonSF;
	  }
	  
          if(IChannel == 2 && cand3leptonChannel == "eemu"){
            sumSFlept_eemu += Dweight[ITypeMC];
            LeptonSF = sel.getLeptonScaleFactor( lept1.Pt(), lept1.Eta(), "e")
                      * sel.getLeptonScaleFactor( lept2.Pt(), lept2.Eta(), "e")
                      * sel.getLeptonScaleFactor( lept3.Pt(), lept3.Eta(), "mu");
 	    if(applyLeptonSFUp)
            LeptonSF = (sel.getLeptonScaleFactor( lept1.Pt(), lept1.Eta(), "e")+sel.getLeptonScaleFactorError( lept1.Pt(), lept1.Eta(), "e"))
                      * (sel.getLeptonScaleFactor( lept2.Pt(), lept2.Eta(), "e")+sel.getLeptonScaleFactorError( lept2.Pt(), lept2.Eta(), "e")) 
                      * (sel.getLeptonScaleFactor( lept3.Pt(), lept3.Eta(), "mu")+sel.getLeptonScaleFactorError( lept3.Pt(), lept3.Eta(), "mu"));		  
	    if(applyLeptonSFDown)
            LeptonSF = (sel.getLeptonScaleFactor( lept1.Pt(), lept1.Eta(), "e")-sel.getLeptonScaleFactorError( lept1.Pt(), lept1.Eta(), "e"))
                      * (sel.getLeptonScaleFactor( lept2.Pt(), lept2.Eta(), "e")-sel.getLeptonScaleFactorError( lept2.Pt(), lept2.Eta(), "e")) 
                      * (sel.getLeptonScaleFactor( lept3.Pt(), lept3.Eta(), "mu")-sel.getLeptonScaleFactorError( lept3.Pt(), lept3.Eta(), "mu"));		  
	    nEvents_eemu += Dweight[ITypeMC]*LeptonSF;
          }
	  
          if(IChannel == 3 && cand3leptonChannel == "eee"){
            sumSFlept_eee += Dweight[ITypeMC];
            LeptonSF = sel.getLeptonScaleFactor( lept1.Pt(), lept1.Eta(), "e")
                      * sel.getLeptonScaleFactor( lept2.Pt(), lept2.Eta(), "e")
                      * sel.getLeptonScaleFactor( lept3.Pt(), lept3.Eta(), "e");
	    if(applyLeptonSFUp)
            LeptonSF = (sel.getLeptonScaleFactor( lept1.Pt(), lept1.Eta(), "e")+sel.getLeptonScaleFactorError( lept1.Pt(), lept1.Eta(), "e"))
                      * (sel.getLeptonScaleFactor( lept2.Pt(), lept2.Eta(), "e")+sel.getLeptonScaleFactorError( lept2.Pt(), lept2.Eta(), "e")) 
                      * (sel.getLeptonScaleFactor( lept3.Pt(), lept3.Eta(), "e")+sel.getLeptonScaleFactorError( lept3.Pt(), lept3.Eta(), "e"));		  
	    if(applyLeptonSFDown)
            LeptonSF = (sel.getLeptonScaleFactor( lept1.Pt(), lept1.Eta(), "e")-sel.getLeptonScaleFactorError( lept1.Pt(), lept1.Eta(), "e"))
                      * (sel.getLeptonScaleFactor( lept2.Pt(), lept2.Eta(), "e")-sel.getLeptonScaleFactorError( lept2.Pt(), lept2.Eta(), "e")) 
                      * (sel.getLeptonScaleFactor( lept3.Pt(), lept3.Eta(), "e")-sel.getLeptonScaleFactorError( lept3.Pt(), lept3.Eta(), "e"));		  
	    nEvents_eee += Dweight[ITypeMC]*LeptonSF;
          }

          //if(event->general.eventNb%10==0) cout<< IChannel << "  " << cand3leptonChannel<<" SF "<<LeptonSF<<endl; 
	  Dweight[ITypeMC]*=LeptonSF;

        }
	
      
	 
	

        //*****************************************************************
        // apply trigger scale factors
        //*****************************************************************  
	
	if(applyTrigger  &&  !isData ){	
	
	  if(IChannel == 0 && cand3leptonChannel == "mumumu")Dweight[ITypeMC]*=SF_trig_mumumu;
          if(IChannel == 1 && cand3leptonChannel == "mumue" )Dweight[ITypeMC]*=SF_trig_mumue;
	  if(IChannel == 2 && cand3leptonChannel == "eemu"  )Dweight[ITypeMC]*=SF_trig_eemu;
	  if(IChannel == 3 && cand3leptonChannel == "eee"   )Dweight[ITypeMC]*=SF_trig_eee;
	
	}
	
	
	
        //*****************************************************************
        // apply DY scale factors
        //***************************************************************** 

	if(datasetName=="Zjets" && applyFakescale ){	
	
	  if(IChannel == 0 && cand3leptonChannel == "mumumu") Dweight[ITypeMC] = SF_Fake[0]*Dweight[ITypeMC];
          if(IChannel == 1 && cand3leptonChannel == "mumue" ) Dweight[ITypeMC] = SF_Fake[1]*Dweight[ITypeMC];
	  if(IChannel == 2 && cand3leptonChannel == "eemu"  ) Dweight[ITypeMC] = SF_Fake[2]*Dweight[ITypeMC];
	  if(IChannel == 3 && cand3leptonChannel == "eee"   ) Dweight[ITypeMC] = SF_Fake[3]*Dweight[ITypeMC];
	
	}
	
	
        //*****************************************************************
        // apply WZ scale factors
        //***************************************************************** 
	
	if(datasetName=="WZ"  && applyWZ ){	
	
	  if(IChannel == 0 && cand3leptonChannel == "mumumu") Dweight[ITypeMC] = SF_WZ[0]*Dweight[ITypeMC];
          if(IChannel == 1 && cand3leptonChannel == "mumue" ) Dweight[ITypeMC] = SF_WZ[1]*Dweight[ITypeMC];
	  if(IChannel == 2 && cand3leptonChannel == "eemu"  ) Dweight[ITypeMC] = SF_WZ[2]*Dweight[ITypeMC];
	  if(IChannel == 3 && cand3leptonChannel == "eee"   ) Dweight[ITypeMC] = SF_WZ[3]*Dweight[ITypeMC];
	
	}
	
	
        //cout << "line 1025" << endl;
	vector<NTJet>        selJets = sel.GetSelectedJets(selMuons, selElectrons, applyJES, scale, applyJER, ResFactor);
   	//vector<NTJet>      selJets = sel.GetSelectedJets();
	double dileptonIvM = dilept.M();
	
	double sumPtLeptJet = lept1.Pt() + lept2.Pt() + lept3.Pt();
	for(unsigned int i=0 ; i< selJets.size(); i++) sumPtLeptJet += selJets[i].p4.Pt();
	double theMET = met.p2.Mod();
   		
	
	 //cout << "-------------------------------------------------------" << endl;    
	
  
  
   	int NBtaggedJets = 0;
   	int idxBtag      = 0;
   	int AlgoBtag = sel.GetbtagAlgo();
   	float btagDiscriCut = sel.GetbtagDiscriCut ();
	
	bool foundASelBjet = 0;
   	//cout << "***********************" << endl;
   	for(unsigned int ijet = 0; ijet < selJets.size(); ijet++){
    	 //cout << "jet flavor " << int(selJets[ijet].partonFlavour) << endl;
    	 if(abs(selJets[ijet].partonFlavour)==5 ){
	  // cout << "found b jet " << int(selJets[ijet].partonFlavour) << endl;
	    foundASelBjet = true;
	 }
    	 if ( AlgoBtag==0 &&  selJets[ijet].bTag["trackCountingHighEffBJetTags"]	 >= btagDiscriCut) NBtaggedJets++;
    	 if ( AlgoBtag==1 &&  selJets[ijet].bTag["simpleSecondaryVertexHighEffBJetTags"] >= btagDiscriCut) NBtaggedJets++;
    	 if ( AlgoBtag==2 &&  selJets[ijet].bTag["trackCountingHighPurBJetTags"]	 >= btagDiscriCut) NBtaggedJets++;
    	 if ( AlgoBtag==3 &&  selJets[ijet].bTag["simpleSecondaryVertexHighPurBJetTags"] >= btagDiscriCut) NBtaggedJets++;
    	 if ( AlgoBtag==4 &&  selJets[ijet].bTag["jetProbabilityBJetTags"]  	         >= btagDiscriCut) NBtaggedJets++;
    	 if ( AlgoBtag==5 &&  selJets[ijet].bTag["jetBProbabilityBJetTags"] 	         >= btagDiscriCut) NBtaggedJets++;
    	 if ( AlgoBtag==6 &&  selJets[ijet].bTag["combinedSecondaryVertexBJetTags"]      >= btagDiscriCut){
   		NBtaggedJets++;
   		idxBtag = ijet;
    	 }
   	}  
    	//cout << "eventNumber " << event->general.eventNb << endl;
   	
	if(foundASelBjet) MyhistoManager.FillHisto(SelABjet, "SelABjet", 1.  , datasetName, IsSignal, 1.);
	else              MyhistoManager.FillHisto(SelABjet, "SelABjet", 0.  , datasetName, IsSignal, 1. );
        
   	//int selLastStep = sel.doFullSelection (dataset, string("all"), false);
   
   
   	 // initialisation of weightb
   	vector < float >weightb;
   	weightb.push_back (1.);
  	weightb.push_back (0.);
   	weightb.push_back (0.);
   	weightb.push_back (0.);
   	weightb.push_back (0.);
   	if (sel.GetFlagb() == 1) {	// check if the weightb computation is needed or not (depending on "flag" in xml)
    	 if (!isData) {	 // to be applied only for MC 
      	 vector < float >weight_temp = sel.GetSFBweight().GetWeigth4BSel (sel.GetMethodb(), sel.GetSystb(),selJets);
      	 weightb[0] = weight_temp[0];  //weight of the event (depending on "NofBtagJets" in xml)
      		 weightb[1] = weight_temp[1];  //proba 0 jet
      	 weightb[2] = weight_temp[2];  //proba 1 jet
      	 weightb[3] = weight_temp[3];  //proba 2 jets
      	 weightb[4] = weight_temp[4];  //proba at least 3 jets	
      	 //cout << " weight btag " << weightb[0] << " " << " p0jet " << weightb[1] << " p1jet " << weightb[2] << " p2jet " << weightb[3]
   	//	   << " prest " << weightb[4] << " check1=" << weightb[1]+weightb[2]+weightb[3]+weightb[4] << endl;
 	
     	}
   	}
   	/*cout << "nseljets " << selJets.size() << endl;
   	cout << "weightb[0] " << weightb[0] << endl;
   	cout << "weightb[1] " << weightb[1] << endl;
   	cout << "weightb[2] " << weightb[2] << endl;
   	cout << "weightb[3] " << weightb[3] << endl;
   	cout << "weightb[4] " << weightb[4] << endl;
*/

	
	int nvertex = selVertices.size();
	  
	if(IChannel == 0 && cand3leptonChannel == "mumumu"){
	
	MyhistoManager.FillHisto(Nvtx_mumumu_afterleptsel, "Nvtx_mumumu_afterleptsel", nvertex, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(InvM_ll_mumumu_afterleptsel, "InvM_ll_mumumu_afterleptsel", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(InvM_ll_mumumu_afterleptsel_lowbin, "InvM_ll_mumumu_afterleptsel_lowbin", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(LeptZPt_mumumu_afterleptsel, "LeptZPt_mumumu_afterleptsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(LeptZPt_mumumu_afterleptsel, "LeptZPt_mumumu_afterleptsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(LeptWPt_mumumu_afterleptsel, "LeptWPt_mumumu_afterleptsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	
	  
	  MyhistoManager.FillHisto(HT_mumumu_afterleptsel, "HT_mumumu_afterleptsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(MET_mumumu_afterleptsel, "MET_mumumu_afterleptsel", theMET, datasetName, IsSignal, Dweight[ITypeMC]);
	  for(unsigned int i=0 ; i< selJets.size(); i++){ 
	    MyhistoManager.FillHisto(JetPt_mumumu_afterleptsel, "JetPt_mumumu_afterleptsel",selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(JetEta_mumumu_afterleptsel, "JetEta_mumumu_afterleptsel",selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto2D(HT_vs_JetPt_mumumu_afterleptsel  ,"HT_vs_JetPt_mumumu_afterleptsel"  , sumPtLeptJet, selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  }
	  MyhistoManager.FillHisto2D(HT_vs_MET_mumumu_afterleptsel    ,"HT_vs_MET_mumumu_afterleptsel"    , sumPtLeptJet, theMET, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_NJet_mumumu_afterleptsel   ,"HT_vs_NJet_mumumu_afterleptsel"   , sumPtLeptJet, selJets.size(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_NBJet_mumumu_afterleptsel  ,"HT_vs_NBJet_mumumu_afterleptsel"  , sumPtLeptJet, NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_LeptPt_mumumu_afterleptsel ,"HT_vs_LeptPt_mumumu_afterleptsel" , sumPtLeptJet, lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_LeptPt_mumumu_afterleptsel ,"HT_vs_LeptPt_mumumu_afterleptsel" , sumPtLeptJet, lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_LeptPt_mumumu_afterleptsel ,"HT_vs_LeptPt_mumumu_afterleptsel" , sumPtLeptJet, lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  if(selJets.size() == 1 && dilept.M() > 20) MyhistoManager.FillHisto2D(HT_vs_Mll_mumumu_afterleptsel	,"HT_vs_Mll_mumumu_afterleptsel"   , sumPtLeptJet, dilept.M(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(Charge_mumumu_afterleptsel,    "Charge_mumumu_afterleptsel"   , wcharge, datasetName, IsSignal, Dweight[ITypeMC]);
	 
	}
        if(IChannel == 1 && cand3leptonChannel == "mumue" ){
	
	MyhistoManager.FillHisto(Nvtx_mumue_afterleptsel,  "Nvtx_mumue_afterleptsel", nvertex, datasetName, IsSignal, Dweight[ITypeMC]); 
	  MyhistoManager.FillHisto(InvM_ll_mumue_afterleptsel,  "InvM_ll_mumue_afterleptsel",  dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);	  
	  MyhistoManager.FillHisto(InvM_ll_mumue_afterleptsel_lowbin,  "InvM_ll_mumue_afterleptsel_lowbin",  dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);	  
	 
	  MyhistoManager.FillHisto(LeptZPt_mumue_afterleptsel, "LeptZPt_mumue_afterleptsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(LeptZPt_mumue_afterleptsel, "LeptZPt_mumue_afterleptsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(LeptWPt_mumue_afterleptsel, "LeptWPt_mumue_afterleptsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  
	  MyhistoManager.FillHisto(HT_mumue_afterleptsel, "HT_mumue_afterleptsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(MET_mumue_afterleptsel, "MET_mumue_afterleptsel", theMET, datasetName, IsSignal, Dweight[ITypeMC]);
	  for(unsigned int i=0 ; i< selJets.size(); i++){
	   MyhistoManager.FillHisto(JetPt_mumue_afterleptsel, "JetPt_mumue_afterleptsel",
	  selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	   MyhistoManager.FillHisto(JetEta_mumue_afterleptsel, "JetEta_mumue_afterleptsel",
	  selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto2D(HT_vs_JetPt_mumue_afterleptsel  ,"HT_vs_JetPt_mumue_afterleptsel"  , sumPtLeptJet, selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  
	  }
	  
	  MyhistoManager.FillHisto2D(HT_vs_MET_mumue_afterleptsel    ,"HT_vs_MET_mumue_afterleptsel"    , sumPtLeptJet, theMET, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_NJet_mumue_afterleptsel   ,"HT_vs_NJet_mumue_afterleptsel"   , sumPtLeptJet, selJets.size(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_NBJet_mumue_afterleptsel  ,"HT_vs_NBJet_mumue_afterleptsel"  , sumPtLeptJet, NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_LeptPt_mumue_afterleptsel ,"HT_vs_LeptPt_mumue_afterleptsel" , sumPtLeptJet, lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_LeptPt_mumue_afterleptsel ,"HT_vs_LeptPt_mumue_afterleptsel" , sumPtLeptJet, lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_LeptPt_mumue_afterleptsel ,"HT_vs_LeptPt_mumue_afterleptsel" , sumPtLeptJet, lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  
	  if(selJets.size() == 1 && dilept.M() > 20) MyhistoManager.FillHisto2D(HT_vs_Mll_mumue_afterleptsel   ,"HT_vs_Mll_mumue_afterleptsel"   , sumPtLeptJet, dilept.M(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(Charge_mumue_afterleptsel,    "Charge_mumue_afterleptsel"   , wcharge, datasetName, IsSignal, Dweight[ITypeMC]);
	
	  
	}
        if(IChannel == 2 && cand3leptonChannel == "eemu"  ){ 
	MyhistoManager.FillHisto(Nvtx_eemu_afterleptsel,   "Nvtx_eemu_afterleptsel", nvertex, datasetName, IsSignal, Dweight[ITypeMC]);
	
	  MyhistoManager.FillHisto(InvM_ll_eemu_afterleptsel,   "InvM_ll_eemu_afterleptsel",   dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);  
	  MyhistoManager.FillHisto(InvM_ll_eemu_afterleptsel_lowbin,   "InvM_ll_eemu_afterleptsel_lowbin",   dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);  
	  
	  MyhistoManager.FillHisto(LeptZPt_eemu_afterleptsel, "LeptZPt_eemu_afterleptsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(LeptZPt_eemu_afterleptsel, "LeptZPt_eemu_afterleptsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(LeptWPt_eemu_afterleptsel, "LeptWPt_eemu_afterleptsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  
	  MyhistoManager.FillHisto(HT_eemu_afterleptsel, "HT_eemu_afterleptsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(MET_eemu_afterleptsel, "MET_eemu_afterleptsel", theMET, datasetName, IsSignal, Dweight[ITypeMC]);
	  for(unsigned int i=0 ; i< selJets.size(); i++){
	   MyhistoManager.FillHisto(JetPt_eemu_afterleptsel, "JetPt_eemu_afterleptsel",
	  selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	   MyhistoManager.FillHisto(JetEta_eemu_afterleptsel, "JetEta_eemu_afterleptsel",
	  selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
           MyhistoManager.FillHisto2D(HT_vs_JetPt_eemu_afterleptsel  ,"HT_vs_JetPt_eemu_afterleptsel"  , sumPtLeptJet, selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
  
	  }
          
	  MyhistoManager.FillHisto2D(HT_vs_MET_eemu_afterleptsel    ,"HT_vs_MET_eemu_afterleptsel"    , sumPtLeptJet, theMET, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_NJet_eemu_afterleptsel   ,"HT_vs_NJet_eemu_afterleptsel"   , sumPtLeptJet, selJets.size(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_NBJet_eemu_afterleptsel  ,"HT_vs_NBJet_eemu_afterleptsel"  , sumPtLeptJet, NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_LeptPt_eemu_afterleptsel ,"HT_vs_LeptPt_eemu_afterleptsel" , sumPtLeptJet, lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_LeptPt_eemu_afterleptsel ,"HT_vs_LeptPt_eemu_afterleptsel" , sumPtLeptJet, lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_LeptPt_eemu_afterleptsel ,"HT_vs_LeptPt_eemu_afterleptsel" , sumPtLeptJet, lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  
	  if(selJets.size() == 1 && dilept.M() > 20) MyhistoManager.FillHisto2D(HT_vs_Mll_eemu_afterleptsel   ,"HT_vs_Mll_eemu_afterleptsel"   , sumPtLeptJet, dilept.M(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(Charge_eemu_afterleptsel,    "Charge_eemu_afterleptsel"   , wcharge, datasetName, IsSignal, Dweight[ITypeMC]);
	
	}
	if(IChannel == 3 && cand3leptonChannel == "eee"   ){ 
	
	MyhistoManager.FillHisto(Nvtx_eee_afterleptsel,    "Nvtx_eee_afterleptsel", nvertex, datasetName, IsSignal, Dweight[ITypeMC]);	
	  MyhistoManager.FillHisto(InvM_ll_eee_afterleptsel,           "InvM_ll_eee_afterleptsel",    dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(InvM_ll_eee_afterleptsel_lowbin,    "InvM_ll_eee_afterleptsel_lowbin",    dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	  
	  MyhistoManager.FillHisto(LeptZPt_eee_afterleptsel, "LeptZPt_eee_afterleptsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(LeptZPt_eee_afterleptsel, "LeptZPt_eee_afterleptsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(LeptWPt_eee_afterleptsel, "LeptWPt_eee_afterleptsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  
	  
	  MyhistoManager.FillHisto(HT_eee_afterleptsel, "HT_eee_afterleptsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(MET_eee_afterleptsel, "MET_eee_afterleptsel", theMET, datasetName, IsSignal, Dweight[ITypeMC]);
	  for(unsigned int i=0 ; i< selJets.size(); i++){
	   MyhistoManager.FillHisto(JetPt_eee_afterleptsel, "JetPt_eee_afterleptsel",
	  selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	   MyhistoManager.FillHisto(JetEta_eee_afterleptsel, "JetEta_eee_afterleptsel",
	  selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
	   MyhistoManager.FillHisto2D(HT_vs_JetPt_eee_afterleptsel  ,"HT_vs_JetPt_eee_afterleptsel"  , sumPtLeptJet, selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
 
	  }
	  
	  MyhistoManager.FillHisto2D(HT_vs_MET_eee_afterleptsel    ,"HT_vs_MET_eee_afterleptsel"    , sumPtLeptJet, theMET, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_NJet_eee_afterleptsel   ,"HT_vs_NJet_eee_afterleptsel"   , sumPtLeptJet, selJets.size(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_NBJet_eee_afterleptsel  ,"HT_vs_NBJet_eee_afterleptsel"  , sumPtLeptJet, NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_LeptPt_eee_afterleptsel ,"HT_vs_LeptPt_eee_afterleptsel" , sumPtLeptJet, lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_LeptPt_eee_afterleptsel ,"HT_vs_LeptPt_eee_afterleptsel" , sumPtLeptJet, lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto2D(HT_vs_LeptPt_eee_afterleptsel ,"HT_vs_LeptPt_eee_afterleptsel" , sumPtLeptJet, lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  
	  if(selJets.size() == 1 && dilept.M() > 20) MyhistoManager.FillHisto2D(HT_vs_Mll_eee_afterleptsel   ,"HT_vs_Mll_eee_afterleptsel"   , sumPtLeptJet, dilept.M(), datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(Charge_eee_afterleptsel,    "Charge_eee_afterleptsel"   , wcharge, datasetName, IsSignal, Dweight[ITypeMC]);
	
	}
	
	
	
	//*****************************************************************
        // pass Z mass
        //*****************************************************************  
	
	
	
	
	if(sumPtLeptJet  > 300 && fabs(dilept.M()-91) > 15){
	  
	  if(isData){
	    /*cout << "******************************************************************* "<< endl;
	    cout << "lepton 1 pt/eta/phi " << lept1.Pt() << " " << lept1.Eta() << " " <<  lept1.Phi() << endl;
	    cout << "lepton 2 pt/eta/phi " << lept2.Pt() << " " << lept2.Eta() << " " <<  lept2.Phi() << endl;
	    cout << "lepton 3 pt/eta/phi " << lept3.Pt() << " " << lept3.Eta() << " " <<  lept3.Phi() << endl;
	    cout << "dilepton Z cand     " <<  dileptonIvM << endl; 
	    cout << "sum pT		   " << sumPtLeptJet  << endl; 
	    cout << "MET		   " <<  theMET << endl; 
	    cout << "Njets		   " <<  selJets.size() << endl; 
	    for(unsigned int i=0 ; i< selJets.size(); i++) cout << " jet " << i << " pt/eta/phi " << selJets[i].p4.Pt() << " " << selJets[i].p4.Eta() << " " <<  selJets[i].p4.Phi()
	     << "  btag discri " << selJets[i].GetDiscri(string("combinedSecondaryVertexBJetTags"))     << endl;
	    */
          }
		  
	  if( IChannel == 0 && cand3leptonChannel == "mumumu")  MyhistoManager.FillHisto(InvM_ll_mumumu_afterleptsel_highSumPt, "InvM_ll_mumumu_afterleptsel_highSumPt", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
          if( IChannel == 1 && cand3leptonChannel == "mumue" )  MyhistoManager.FillHisto(InvM_ll_mumue_afterleptsel_highSumPt,  "InvM_ll_mumue_afterleptsel_highSumPt",  dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);	  
          if( IChannel == 2 && cand3leptonChannel == "eemu"  )  MyhistoManager.FillHisto(InvM_ll_eemu_afterleptsel_highSumPt,   "InvM_ll_eemu_afterleptsel_highSumPt",   dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);  
	  if( IChannel == 3 && cand3leptonChannel == "eee"   )  MyhistoManager.FillHisto(InvM_ll_eee_afterleptsel_highSumPt,    "InvM_ll_eee_afterleptsel_highSumPt",    dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	}  
	
	//*****************************************************************
	// pass 3 lept selection
	//*****************************************************************  
	
	//fixme
	if( fabs(dilept.M()-91) < 15){
	//if( fabs(dilept.M()-91) < 10000000000){
	   //cout << " weight btag pass Z " << weightb[0] << " " << " p0jet " << weightb[1] << " p1jet " << weightb[2] << " p2jet " << weightb[3]
   //		   << " prest " << weightb[4] << " check1=" << weightb[1]+weightb[2]+weightb[3]+weightb[4] << endl;;
 	
	  if(dilept.M() > 300){
	    /*cout << "dileptonIvM " << dileptonIvM
	    << "  dilept.Pt() "    << dilept.Pt() 
	    << "  lept1.Pt() "     << lept1.Pt() 
	    << "  lept2.Pt() "     << lept2.Pt()
	    << "  lept3.Pt() "     << lept3.Pt()
	    << "  Njet " << selJets.size();*/
	    for(unsigned int ijet=0; ijet<selJets.size(); ijet++) cout <<  "  jet " << ijet << " pt " << selJets[ijet].p4.Pt() << endl;
	    //cout << "  met       " << met.p2.Mod() << endl;
	    
	    /*if(selJets.size()>0) cout << "  deltaphi lept-jet 1 " << lept1.DeltaR(selJets[0].p4) << 
	                                 "  deltaphi lept-jet 2 " << lept2.DeltaR(selJets[0].p4) <<
	                                 "  deltaphi lept-jet 3 " << lept3.DeltaR(selJets[0].p4) << endl;
					 
	    cout << "  deltaphi lept-met 1 " << lept1.Phi() - met.p2.Phi()  << 
	            "  deltaphi lept-met 2 " << lept2.Phi() - met.p2.Phi()  <<
	            "  deltaphi lept-met 3 " << lept3.Phi() - met.p2.Phi()  << endl;  
  */
	  }
	  
	  
	  
 	  if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(CutFlow_mumumu, "CutFlow_mumumu", 2, datasetName, IsSignal, Dweight[ITypeMC]);
 	  if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(CutFlow_mumue,  "CutFlow_mumue" , 2, datasetName, IsSignal, Dweight[ITypeMC]);
 	  if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(CutFlow_eemu,   "CutFlow_eemu"  , 2, datasetName, IsSignal, Dweight[ITypeMC]);
 	  if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(CutFlow_eee,    "CutFlow_eee"   , 2, datasetName, IsSignal, Dweight[ITypeMC]);
 	  if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(ErrCutFlow_mumumu,   "ErrCutFlow_mumumu"  , 2, datasetName, IsSignal, EventYieldWeightError);
 	  if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(ErrCutFlow_mumue,    "ErrCutFlow_mumue"   , 2, datasetName, IsSignal, EventYieldWeightError);
 	  if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(ErrCutFlow_eemu,     "ErrCutFlow_eemu"    , 2, datasetName, IsSignal, EventYieldWeightError);
 	  if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(ErrCutFlow_eee,      "ErrCutFlow_eee"     , 2, datasetName, IsSignal, EventYieldWeightError);
	
	  if(IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(WmissAssing_mumumu_afterleptsel , "WmissAssing_mumumu_afterleptsel", leptonFlavor, datasetName, IsSignal, Dweight[ITypeMC]);
	  if(IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(WmissAssing_mumue_afterleptsel  , "WmissAssing_mumue_afterleptsel" , leptonFlavor, datasetName, IsSignal, Dweight[ITypeMC]);
	  if(IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(WmissAssing_eemu_afterleptsel   , "WmissAssing_eemu_afterleptsel"  , leptonFlavor, datasetName, IsSignal, Dweight[ITypeMC]);
	  if(IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(WmissAssing_eee_afterleptsel    , "WmissAssing_eee_afterleptsel"   , leptonFlavor, datasetName, IsSignal, Dweight[ITypeMC]);

	  
	  
	  
	  
	  double mTW = pow(
	  			2*lept3.Pt()*met.p2.Mod()*(1-cos(lept3.Phi() - met.p2.Phi()))
	   			,0.5);
	  
	  if(IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto( mWT_mumumu_afterleptsel, "mWT_mumumu_afterleptsel", mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	  if(IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto( mWT_mumue_afterleptsel,  "mWT_mumue_afterleptsel" , mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	  if(IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto( mWT_eemu_afterleptsel,   "mWT_eemu_afterleptsel"  , mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	  if(IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto( mWT_eee_afterleptsel,    "mWT_eee_afterleptsel"   , mTW, datasetName, IsSignal, Dweight[ITypeMC]);
  
	  if( IChannel == 0 && cand3leptonChannel == "mumumu")MyhistoManager.FillHisto2D(InvM_ll_vs_mWT_mumumu_afterleptsel, "InvM_ll_vs_mWT_mumumu_afterleptsel",dilept.M(), mTW, datasetName, IsSignal, Dweight[ITypeMC]);
          if( IChannel == 1 && cand3leptonChannel == "mumue") MyhistoManager.FillHisto2D(InvM_ll_vs_mWT_mumue_afterleptsel,  "InvM_ll_vs_mWT_mumue_afterleptsel" ,dilept.M(), mTW, datasetName, IsSignal, Dweight[ITypeMC]);
          if( IChannel == 2 && cand3leptonChannel == "eemu")  MyhistoManager.FillHisto2D(InvM_ll_vs_mWT_eemu_afterleptsel,   "InvM_ll_vs_mWT_eemu_afterleptsel"  ,dilept.M(), mTW, datasetName, IsSignal, Dweight[ITypeMC]);
          if( IChannel == 3 && cand3leptonChannel == "eee")   MyhistoManager.FillHisto2D(InvM_ll_vs_mWT_eee_afterleptsel,    "InvM_ll_vs_mWT_eee_afterleptsel"   ,dilept.M(), mTW, datasetName, IsSignal, Dweight[ITypeMC]);
 
	
	  
	  if(mTW > 110){
	    /*cout << "interesting events !!! " << std::endl;
	    cout<<"RUN "<<event->general.runNb<<" EVT "<<event->general.eventNb<<endl;
	    cout << "mTW " << mTW << endl;
	    cout << "dileptonIvM " << dileptonIvM
	    << "  dilept.Pt() "    << dilept.Pt() 
	    << "  lept1.Pt() "     << lept1.Pt() 
	    << "  lept2.Pt() "     << lept2.Pt()
	    << "  lept3.Pt() "     << lept3.Pt()
	    << "  Njet " << selJets.size();*/
	    for(unsigned int ijet=0; ijet<selJets.size(); ijet++) cout <<  "  jet " << ijet << " pt " << selJets[ijet].p4.Pt() << endl;
	    //cout << "  met       " << met.p2.Mod() << endl;
	    
	    /*if(selJets.size()>0) cout << "  deltaphi lept-jet 1 " << lept1.DeltaR(selJets[0].p4) << 
	                                 "  deltaphi lept-jet 2 " << lept2.DeltaR(selJets[0].p4) <<
	                                 "  deltaphi lept-jet 3 " << lept3.DeltaR(selJets[0].p4) << endl;
					 
	    cout << "  deltaphi lept-met 1 " << lept1.Phi() - met.p2.Phi()  << 
	            "  deltaphi lept-met 2 " << lept2.Phi() - met.p2.Phi()  <<
	            "  deltaphi lept-met 3 " << lept3.Phi() - met.p2.Phi()  << endl;
	   */
	   TLorentzVector tmpmet;
	   tmpmet.SetPxPyPzE(met.p2.Px(),met.p2.Py(), 0, 0 );
	   
	   int njetselec = selJets.size();
	   if(njetselec>=4) njetselec = 4;
	   
	   double deltaPhi1 = lept1.DeltaR(tmpmet);
	   double deltaPhi2 = lept2.DeltaR(tmpmet);
	   double deltaPhi3 = lept3.DeltaR(tmpmet);
	   
	   
	   
	   
	   if(IChannel == 0 && cand3leptonChannel == "mumumu"){
	     MyhistoManager.FillHisto(LeptWPt_mumumu_afterleptsel_mWT110, "LeptWPt_mumumu_afterleptsel_mWT110", lept3.Pt(),  datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(MET_mumumu_afterleptsel_mWT110,     "MET_mumumu_afterleptsel_mWT110",     theMET,      datasetName, IsSignal, Dweight[ITypeMC]);  
	     MyhistoManager.FillHisto(InvM_ll_mumumu_afterleptsel_mWT110, "InvM_ll_mumumu_afterleptsel_mWT110", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	     
	     for(unsigned int ijet=0; ijet<selJets.size(); ijet++){
	       MyhistoManager.FillHisto(deltaRLeptJet_mumumu_afterleptsel_mWT110, "deltaRLeptJet_mumumu_afterleptsel_mWT110", lept1.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	       MyhistoManager.FillHisto(deltaRLeptJet_mumumu_afterleptsel_mWT110, "deltaRLeptJet_mumumu_afterleptsel_mWT110", lept2.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	       MyhistoManager.FillHisto(deltaRLeptJet_mumumu_afterleptsel_mWT110, "deltaRLeptJet_mumumu_afterleptsel_mWT110", lept3.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	     }
	     
	      MyhistoManager.FillHisto(Charge_mumumu_afterleptsel_mWT110, "Charge_mumumu_afterleptsel_mWT110", wcharge, datasetName, IsSignal, Dweight[ITypeMC]);
	  
	     MyhistoManager.FillHisto(deltaRLeptMet_mumumu_afterleptsel_mWT110, "deltaRLeptMet_mumumu_afterleptsel_mWT110",deltaPhi1 , datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(deltaRLeptMet_mumumu_afterleptsel_mWT110, "deltaRLeptMet_mumumu_afterleptsel_mWT110",deltaPhi2 , datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(deltaRLeptMet_mumumu_afterleptsel_mWT110, "deltaRLeptMet_mumumu_afterleptsel_mWT110",deltaPhi3 , datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(NJet_mumumu_afterleptsel_mWT110  , "NJet_mumumu_afterleptsel_mWT110" , njetselec,    datasetName, IsSignal, Dweight[ITypeMC]);    
             MyhistoManager.FillHisto(NBJet_mumumu_afterleptsel_mWT110 , "NBJet_mumumu_afterleptsel_mWT110", NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]); 
	     
	     
	   }
	   if(IChannel == 1 && cand3leptonChannel == "mumue" ){ 
	     MyhistoManager.FillHisto(LeptWPt_mumue_afterleptsel_mWT110,"LeptWPt_mumue_afterleptsel_mWT110",  lept3.Pt(),  datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(MET_mumue_afterleptsel_mWT110,    "MET_mumue_afterleptsel_mWT110",      theMET,      datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(InvM_ll_mumue_afterleptsel_mWT110,"InvM_ll_mumue_afterleptsel_mWT110",  dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	     
	     
	     for(unsigned int ijet=0; ijet<selJets.size(); ijet++){
	       MyhistoManager.FillHisto(deltaRLeptJet_mumue_afterleptsel_mWT110, "deltaRLeptJet_mumue_afterleptsel_mWT110", lept1.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	       MyhistoManager.FillHisto(deltaRLeptJet_mumue_afterleptsel_mWT110, "deltaRLeptJet_mumue_afterleptsel_mWT110", lept2.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	       MyhistoManager.FillHisto(deltaRLeptJet_mumue_afterleptsel_mWT110, "deltaRLeptJet_mumue_afterleptsel_mWT110", lept3.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	     }
	        MyhistoManager.FillHisto(Charge_mumue_afterleptsel_mWT110,  "Charge_mumue_afterleptsel_mWT110" , wcharge, datasetName, IsSignal, Dweight[ITypeMC]);
	
	     
	     MyhistoManager.FillHisto(deltaRLeptMet_mumue_afterleptsel_mWT110, "deltaRLeptMet_mumue_afterleptsel_mWT110",deltaPhi1 , datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(deltaRLeptMet_mumue_afterleptsel_mWT110, "deltaRLeptMet_mumue_afterleptsel_mWT110",deltaPhi2 , datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(deltaRLeptMet_mumue_afterleptsel_mWT110, "deltaRLeptMet_mumue_afterleptsel_mWT110",deltaPhi3 , datasetName, IsSignal, Dweight[ITypeMC]);
	
	     MyhistoManager.FillHisto(NJet_mumue_afterleptsel_mWT110  , "NJet_mumue_afterleptsel_mWT110" , njetselec,    datasetName, IsSignal, Dweight[ITypeMC]);    
             MyhistoManager.FillHisto(NBJet_mumue_afterleptsel_mWT110 , "NBJet_mumue_afterleptsel_mWT110", NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]); 
	     
	   }
           if(IChannel == 2 && cand3leptonChannel == "eemu"   ){ 
	     MyhistoManager.FillHisto(LeptWPt_eemu_afterleptsel_mWT110, "LeptWPt_eemu_afterleptsel_mWT110",   lept3.Pt(),  datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(MET_eemu_afterleptsel_mWT110,     "MET_eemu_afterleptsel_mWT110",       theMET,      datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(InvM_ll_eemu_afterleptsel_mWT110, "InvM_ll_eemu_afterleptsel_mWT110",   dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	     
	     for(unsigned int ijet=0; ijet<selJets.size(); ijet++){
	       MyhistoManager.FillHisto(deltaRLeptJet_eemu_afterleptsel_mWT110, "deltaRLeptJet_eemu_afterleptsel_mWT110", lept1.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	       MyhistoManager.FillHisto(deltaRLeptJet_eemu_afterleptsel_mWT110, "deltaRLeptJet_eemu_afterleptsel_mWT110", lept2.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	       MyhistoManager.FillHisto(deltaRLeptJet_eemu_afterleptsel_mWT110, "deltaRLeptJet_eemu_afterleptsel_mWT110", lept3.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	     }
	     
	     
	     MyhistoManager.FillHisto(Charge_eemu_afterleptsel_mWT110,   "Charge_eemu_afterleptsel_mWT110"  , wcharge, datasetName, IsSignal, Dweight[ITypeMC]);
	 
	      
	     MyhistoManager.FillHisto(deltaRLeptMet_eemu_afterleptsel_mWT110, "deltaRLeptMet_eemu_afterleptsel_mWT110",deltaPhi1 , datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(deltaRLeptMet_eemu_afterleptsel_mWT110, "deltaRLeptMet_eemu_afterleptsel_mWT110",deltaPhi2 , datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(deltaRLeptMet_eemu_afterleptsel_mWT110, "deltaRLeptMet_eemu_afterleptsel_mWT110",deltaPhi3 , datasetName, IsSignal, Dweight[ITypeMC]);
	
	     MyhistoManager.FillHisto(NJet_eemu_afterleptsel_mWT110  , "NJet_eemu_afterleptsel_mWT110" , njetselec,    datasetName, IsSignal, Dweight[ITypeMC]);    
             MyhistoManager.FillHisto(NBJet_eemu_afterleptsel_mWT110 , "NBJet_eemu_afterleptsel_mWT110", NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]); 
	     
	   }
           if(IChannel == 3 && cand3leptonChannel == "eee"   ) {
	     MyhistoManager.FillHisto(LeptWPt_eee_afterleptsel_mWT110,  "LeptWPt_eee_afterleptsel_mWT110",    lept3.Pt(),  datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(MET_eee_afterleptsel_mWT110,      "MET_eee_afterleptsel_mWT110",        theMET,      datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(InvM_ll_eee_afterleptsel_mWT110,  "InvM_ll_eee_afterleptsel_mWT110",    dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	     
	     for(unsigned int ijet=0; ijet<selJets.size(); ijet++){
	       MyhistoManager.FillHisto(deltaRLeptJet_eee_afterleptsel_mWT110, "deltaRLeptJet_eee_afterleptsel_mWT110", lept1.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	       MyhistoManager.FillHisto(deltaRLeptJet_eee_afterleptsel_mWT110, "deltaRLeptJet_eee_afterleptsel_mWT110", lept2.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	       MyhistoManager.FillHisto(deltaRLeptJet_eee_afterleptsel_mWT110, "deltaRLeptJet_eee_afterleptsel_mWT110", lept3.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	     }
	     
	     
	     MyhistoManager.FillHisto(Charge_eee_afterleptsel_mWT110,    "Charge_eee_afterleptsel_mWT110"   , wcharge, datasetName, IsSignal, Dweight[ITypeMC]);
	 
	      
	     MyhistoManager.FillHisto(deltaRLeptMet_eee_afterleptsel_mWT110, "deltaRLeptMet_eee_afterleptsel_mWT110",deltaPhi1 , datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(deltaRLeptMet_eee_afterleptsel_mWT110, "deltaRLeptMet_eee_afterleptsel_mWT110",deltaPhi2 , datasetName, IsSignal, Dweight[ITypeMC]);
	     MyhistoManager.FillHisto(deltaRLeptMet_eee_afterleptsel_mWT110, "deltaRLeptMet_eee_afterleptsel_mWT110",deltaPhi3 , datasetName, IsSignal, Dweight[ITypeMC]);
	
	     MyhistoManager.FillHisto(NJet_eee_afterleptsel_mWT110  , "NJet_eee_afterleptsel_mWT110" , njetselec,    datasetName, IsSignal, Dweight[ITypeMC]);    
             MyhistoManager.FillHisto(NBJet_eee_afterleptsel_mWT110 , "NBJet_eee_afterleptsel_mWT110", NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]); 
	     
	   }

	    
	
	  }
	  
	  
	  
	  //*****************************************************************
	  // pass jet selection
	  //*****************************************************************  
	  int NSeljets = selJets.size() ;
	  if(NSeljets>4) NSeljets = 4;

	    
	  if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(NJet_mumumu_afterZsel , "NJet_mumumu_afterZsel",NSeljets, datasetName, IsSignal, Dweight[ITypeMC]);
          if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(NJet_mumue_afterZsel  , "NJet_mumue_afterZsel" ,NSeljets, datasetName, IsSignal, Dweight[ITypeMC]);
          if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(NJet_eemu_afterZsel   , "NJet_eemu_afterZsel"  ,NSeljets, datasetName, IsSignal, Dweight[ITypeMC]);
          if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(NJet_eee_afterZsel    , "NJet_eee_afterZsel"   ,NSeljets, datasetName, IsSignal, Dweight[ITypeMC]);
	  
	  if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(NBJet_mumumu_afterZsel , "NBJet_mumumu_afterZsel",NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
          if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(NBJet_mumue_afterZsel  , "NBJet_mumue_afterZsel" ,NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
          if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(NBJet_eemu_afterZsel   , "NBJet_eemu_afterZsel"  ,NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
          if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(NBJet_eee_afterZsel    , "NBJet_eee_afterZsel"   ,NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
	   

	  double Zpt = (lept1+lept2).Pt();
	  if( IChannel == 0 && cand3leptonChannel == "mumumu") {
	    MyhistoManager.FillHisto(RecoPtZ_mumumu_afterleptsel,     "RecoPtZ_mumumu_afterleptsel",	Zpt ,	  datasetName, IsSignal, Dweight[ITypeMC]);

	
	  }
	  if( IChannel == 1 && cand3leptonChannel == "mumue" ) {
	    MyhistoManager.FillHisto(RecoPtZ_mumue_afterleptsel,      "RecoPtZ_mumue_afterleptsel",	Zpt  ,     datasetName, IsSignal, Dweight[ITypeMC]);
 
	
	  }
	  if( IChannel == 2 && cand3leptonChannel == "eemu"  ) {
	    MyhistoManager.FillHisto(RecoPtZ_eemu_afterleptsel,       "RecoPtZ_eemu_afterleptsel",	Zpt  ,     datasetName, IsSignal, Dweight[ITypeMC]);

	  }
	  if( IChannel == 3 && cand3leptonChannel == "eee"   ) {
	    MyhistoManager.FillHisto(RecoPtZ_eee_afterleptsel,        "RecoPtZ_eee_afterleptsel",	 Zpt ,     datasetName, IsSignal, Dweight[ITypeMC]);	       
	
	  }
	  
	  
	  if( NSeljets== 0 ){
	    if( IChannel == 0 && cand3leptonChannel == "mumumu") {
	      MyhistoManager.FillHisto(RecoPtZ_mumumu_afterleptsel_nojet,	"RecoPtZ_mumumu_afterleptsel_nojet",    Zpt ,     datasetName, IsSignal, Dweight[ITypeMC]);

	
	    }
	    if( IChannel == 1 && cand3leptonChannel == "mumue" ) {
	      MyhistoManager.FillHisto(RecoPtZ_mumue_afterleptsel_nojet,	"RecoPtZ_mumue_afterleptsel_nojet",	  Zpt  ,     datasetName, IsSignal, Dweight[ITypeMC]);
 
	
	    }
	    if( IChannel == 2 && cand3leptonChannel == "eemu"  ) {
	      MyhistoManager.FillHisto(RecoPtZ_eemu_afterleptsel_nojet,	"RecoPtZ_eemu_afterleptsel_nojet",	  Zpt  ,     datasetName, IsSignal, Dweight[ITypeMC]);

	    }
	    if( IChannel == 3 && cand3leptonChannel == "eee"   ) {
	      MyhistoManager.FillHisto(RecoPtZ_eee_afterleptsel_nojet,	"RecoPtZ_eee_afterleptsel_nojet",	   Zpt ,     datasetName, IsSignal, Dweight[ITypeMC]);  	 
	
	    }
	  }
	  
	  
	  if(selJets.size() >= 1 ){
         //cout << " weight btag jet sel " << weightb[0] << " " << " p0jet " << weightb[1] << " p1jet " << weightb[2] << " p2jet " << weightb[3]
   	//	   << " prest " << weightb[4] << " check1=" << weightb[1]+weightb[2]+weightb[3]+weightb[4] << endl;;
 	   
	   
	   if(foundASelBjet) MyhistoManager.FillHisto(SelABjet_afterjetsel, "SelABjet_afterjetsel", 1.  , datasetName, IsSignal, 1.);
	   else              MyhistoManager.FillHisto(SelABjet_afterjetsel, "SelABjet_afterjetsel", 0.  , datasetName, IsSignal, 1. );
        
   
        //if( IChannel == 0 && cand3leptonChannel == "mumumu")
        // cout<<"RUN "<<event->general.runNb<<" EVT "<<event->general.eventNb<<endl;
	    	  	
 	    if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(CutFlow_mumumu, "CutFlow_mumumu", 3, datasetName, IsSignal, Dweight[ITypeMC]);
 	    if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(CutFlow_mumue,  "CutFlow_mumue" , 3, datasetName, IsSignal, Dweight[ITypeMC]);
 	    if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(CutFlow_eemu,   "CutFlow_eemu"  , 3, datasetName, IsSignal, Dweight[ITypeMC]);
 	    if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(CutFlow_eee,    "CutFlow_eee"   , 3, datasetName, IsSignal, Dweight[ITypeMC]);
 	    if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(ErrCutFlow_mumumu,   "ErrCutFlow_mumumu"  , 3, datasetName, IsSignal, EventYieldWeightError);
 	    if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(ErrCutFlow_mumue,    "ErrCutFlow_mumue"   , 3, datasetName, IsSignal, EventYieldWeightError);
 	    if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(ErrCutFlow_eemu,     "ErrCutFlow_eemu"    , 3, datasetName, IsSignal, EventYieldWeightError);
 	    if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(ErrCutFlow_eee,      "ErrCutFlow_eee"     , 3, datasetName, IsSignal, EventYieldWeightError);
	
		   
		 
	    if(IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto( mWT_mumumu_afterjetsel, "mWT_mumumu_afterjetsel", mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	    if(IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto( mWT_mumue_afterjetsel,  "mWT_mumue_afterjetsel" , mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	    if(IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto( mWT_eemu_afterjetsel,   "mWT_eemu_afterjetsel"  , mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	    if(IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto( mWT_eee_afterjetsel,    "mWT_eee_afterjetsel"   , mTW, datasetName, IsSignal, Dweight[ITypeMC]);

	   
		   
		   
	    if( IChannel == 0 && cand3leptonChannel == "mumumu"){
	    
	      MyhistoManager.FillHisto(InvM_ll_mumumu_afterjetsel, "InvM_ll_mumumu_afterjetsel", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	      //MyhistoManager.FillHisto(NBJet_mumumu_afterjetsel ,  "NBJet_mumumu_afterjetsel",NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
              
	      
	      if(isData) MyhistoManager.FillHisto(NBJet_mumumu_afterjetsel ,  "NBJet_mumumu_afterjetsel",NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
              else{
	        MyhistoManager.FillHisto(NBJet_mumumu_afterjetsel ,  "NBJet_mumumu_afterjetsel",0, datasetName, IsSignal, weightb[1]*Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(NBJet_mumumu_afterjetsel ,  "NBJet_mumumu_afterjetsel",1, datasetName, IsSignal, weightb[2]*Dweight[ITypeMC]);	
		
		       
	      }
	      
	      MyhistoManager.FillHisto(LeptZPt_mumumu_afterjetsel,  "LeptZPt_mumumu_afterjetsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(LeptZPt_mumumu_afterjetsel,  "LeptZPt_mumumu_afterjetsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(LeptWPt_mumumu_afterjetsel,  "LeptWPt_mumumu_afterjetsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      
	      MyhistoManager.FillHisto(HT_mumumu_afterjetsel,      "HT_mumumu_afterjetsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
 	      MyhistoManager.FillHisto(MET_mumumu_afterjetsel, "MET_mumumu_afterjetsel",theMET , datasetName, IsSignal, Dweight[ITypeMC]);
              for(unsigned int i=0 ; i< selJets.size(); i++){ 
	        MyhistoManager.FillHisto(JetPt_mumumu_afterjetsel, "JetPt_mumumu_afterjetsel",
	        selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(JetEta_mumumu_afterjetsel, "JetEta_mumumu_afterjetsel",
	        selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
	      }
	
	    }
	    if( IChannel == 1 && cand3leptonChannel == "mumue" ){
	      
	      MyhistoManager.FillHisto(InvM_ll_mumue_afterjetsel, "InvM_ll_mumue_afterjetsel", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	      
	      //MyhistoManager.FillHisto(NBJet_mumue_afterjetsel ,  "NBJet_mumue_afterjetsel",NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
              
	      if(isData) MyhistoManager.FillHisto(NBJet_mumue_afterjetsel ,  "NBJet_mumue_afterjetsel",NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
              else{
	        MyhistoManager.FillHisto(NBJet_mumue_afterjetsel ,  "NBJet_mumue_afterjetsel",0, datasetName, IsSignal, weightb[1]*Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(NBJet_mumue_afterjetsel ,  "NBJet_mumue_afterjetsel",1, datasetName, IsSignal, weightb[2]*Dweight[ITypeMC]);	        
	      }
	      
	      
	      MyhistoManager.FillHisto(LeptPt_mumue_afterjetsel, "LeptPt_mumue_afterjetsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(LeptPt_mumue_afterjetsel, "LeptPt_mumue_afterjetsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(LeptPt_mumue_afterjetsel, "LeptPt_mumue_afterjetsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(HT_mumue_afterjetsel, "HT_mumue_afterjetsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(MET_mumue_afterjetsel, "MET_mumue_afterjetsel", theMET, datasetName, IsSignal, Dweight[ITypeMC]);
              for(unsigned int i=0 ; i< selJets.size(); i++){
	         MyhistoManager.FillHisto(JetPt_mumue_afterjetsel, "JetPt_mumue_afterjetsel",
	      selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	         MyhistoManager.FillHisto(JetEta_mumue_afterjetsel, "JetEta_mumue_afterjetsel",
	      selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
	      }

	    }
	    if( IChannel == 2 && cand3leptonChannel == "eemu"  ){
	    
	      MyhistoManager.FillHisto(InvM_ll_eemu_afterjetsel, "InvM_ll_eemu_afterjetsel", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
              //MyhistoManager.FillHisto(NBJet_eemu_afterjetsel ,  "NBJet_eemu_afterjetsel",NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
             
	      if(isData) MyhistoManager.FillHisto(NBJet_eemu_afterjetsel ,  "NBJet_eemu_afterjetsel",NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
              else{
	        MyhistoManager.FillHisto(NBJet_eemu_afterjetsel ,  "NBJet_eemu_afterjetsel",0, datasetName, IsSignal, weightb[1]*Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(NBJet_eemu_afterjetsel ,  "NBJet_eemu_afterjetsel",1, datasetName, IsSignal, weightb[2]*Dweight[ITypeMC]);	        
	      }
	      
	      MyhistoManager.FillHisto(LeptZPt_eemu_afterjetsel, "LeptZPt_eemu_afterjetsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(LeptZPt_eemu_afterjetsel, "LeptZPt_eemu_afterjetsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(LeptWPt_eemu_afterjetsel, "LeptWPt_eemu_afterjetsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      
	      MyhistoManager.FillHisto(HT_eemu_afterjetsel, "HT_eemu_afterjetsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(MET_eemu_afterjetsel, "MET_eemu_afterjetsel", theMET, datasetName, IsSignal, Dweight[ITypeMC]);
              for(unsigned int i=0 ; i< selJets.size(); i++){ 
	        MyhistoManager.FillHisto(JetPt_eemu_afterjetsel, "JetPt_eemu_afterjetsel",
	        selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(JetEta_eemu_afterjetsel, "JetEta_eemu_afterjetsel",
	        selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
	      }

	    }
	    
	    if( IChannel == 3 && cand3leptonChannel == "eee"   ){
	    
	      MyhistoManager.FillHisto(InvM_ll_eee_afterjetsel, "InvM_ll_eee_afterjetsel", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	      
	      //MyhistoManager.FillHisto(NBJet_eee_afterjetsel ,  "NBJet_eee_afterjetsel",NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
              
	      if(isData) MyhistoManager.FillHisto(NBJet_eee_afterjetsel ,  "NBJet_eee_afterjetsel",NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
              else{
	        MyhistoManager.FillHisto(NBJet_eee_afterjetsel ,  "NBJet_eee_afterjetsel",0, datasetName, IsSignal, weightb[1]*Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(NBJet_eee_afterjetsel ,  "NBJet_eee_afterjetsel",1, datasetName, IsSignal, weightb[2]*Dweight[ITypeMC]);	        
	      }
	      
	      MyhistoManager.FillHisto(LeptZPt_eee_afterjetsel, "LeptZPt_eee_afterjetsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(LeptZPt_eee_afterjetsel, "LeptZPt_eee_afterjetsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(LeptWPt_eee_afterjetsel, "LeptWPt_eee_afterjetsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(HT_eee_afterjetsel, "HT_eee_afterjetsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(MET_eee_afterjetsel, "MET_eee_afterjetsel", theMET, datasetName, IsSignal, Dweight[ITypeMC]);
              for(unsigned int i=0 ; i< selJets.size(); i++){
	        MyhistoManager.FillHisto(JetPt_eee_afterjetsel, "JetPt_eee_afterjetsel",
	        selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(JetEta_eee_afterjetsel, "JetEta_eee_afterjetsel",
	        selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
	      }

	    }
	    
	    
	    double Zpt = (lept1+lept2).Pt();
	    double deltaPhilj = fabs(lept3.DeltaPhi(selJets[0].p4));
	    TLorentzVector metP4(met.p2.Px(), met.p2.Py(), 0, sqrt(met.p2.Px()*met.p2.Px() + met.p2.Py()*met.p2.Py()));
	    TLorentzVector transTop = lept3 + selJets[0].p4 + metP4;

            // Top mass computation
	     //lept3 + selJets[0].p4 + met.p4;

	    double term1 = lept3.Pz()*(lept3.Px()*met.p2.Px() + lept3.Py()*met.p2.Py() + (80.399)*(80.399)/2.);
	    double det = lept3.Px()*met.p2.Px() + lept3.Py()*met.p2.Py() + (80.399)*(80.399)/2.
			      - met.met()*met.met()*(lept3.E()*lept3.E() - lept3.Pz()*lept3.Pz());
	    if(det<0) det=0;
	    double term2 = lept3.E()*pow(det, 0.5);
	    double denom = lept3.E()*lept3.E() - lept3.Pz()*lept3.Pz();


	    double sol1 = (term1 - term2)/denom;
	    //double sol2 = (term1 + term2)/denom;

	    double neutrE = pow( pow(met.p2.Px(),2) + pow(met.p2.Py(),2) + sol1*sol1, 0.5);//neglecting neut mass

	    TLorentzVector neutrino;
	    neutrino.SetPxPyPzE( met.p2.Px(), met.p2.Py(), sol1, neutrE);

	    TLorentzVector topCand = neutrino + lept3 + selJets[0].p4 ;
	    

	    if(NBtaggedJets == 0 ){ // bjet veto
	      if( IChannel == 0 && cand3leptonChannel == "mumumu") {
		MyhistoManager.FillHisto(HT_mumumu_afterbjetveto, "HT_mumumu_afterbjetveto", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoPtZ_mumumu_afterbjetveto,     "RecoPtZ_mumumu_afterbjetveto",    Zpt ,     datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto( mWT_mumumu_afterbjetveto, "mWT_mumumu_afterbjetveto", mTW, datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(Mt_mumumu_afterbjetveto, "Mt_mumumu_afterbjetveto" , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(deltaPhilj_mumumu_afterbjetveto , "deltaPhilj_mumumu_afterbjetveto",  deltaPhilj ,    datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(LeptWPt_mumumu_afterbjetveto, "LeptWPt_mumumu_afterbjetveto", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoTopMass_mumumu_afterbjetveto, "RecoTopMass_mumumu_afterbjetveto", topCand.M() ,datasetName, IsSignal, Dweight[ITypeMC]);
		for(unsigned int i=0 ; i< selJets.size(); i++){
		  MyhistoManager.FillHisto(JetPt_mumumu_afterbjetveto, "JetPt_mumumu_afterbjetveto",
	          selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		  MyhistoManager.FillHisto(JetEta_mumumu_afterbjetveto, "JetEta_mumumu_afterbjetveto",
	          selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
		}
	      }
	      if( IChannel == 1 && cand3leptonChannel == "mumue" ) {
		MyhistoManager.FillHisto(HT_mumue_afterbjetveto, "HT_mumue_afterbjetveto", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoPtZ_mumue_afterbjetveto,      "RecoPtZ_mumue_afterbjetveto",     Zpt  ,     datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto( mWT_mumue_afterbjetveto, "mWT_mumue_afterbjetveto", mTW, datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(Mt_mumue_afterbjetveto, "Mt_mumue_afterbjetveto" , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(deltaPhilj_mumue_afterbjetveto , "deltaPhilj_mumue_afterbjetveto",  deltaPhilj ,    datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(LeptWPt_mumue_afterbjetveto, "LeptWPt_mumue_afterbjetveto", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoTopMass_mumue_afterbjetveto, "RecoTopMass_mumue_afterbjetveto", topCand.M() ,datasetName, IsSignal, Dweight[ITypeMC]);
                for(unsigned int i=0 ; i< selJets.size(); i++){
		  MyhistoManager.FillHisto(JetPt_mumue_afterbjetveto, "JetPt_mumue_afterbjetveto",
	          selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		  MyhistoManager.FillHisto(JetEta_mumue_afterbjetveto, "JetEta_mumue_afterbjetveto",
	          selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
		}
	      }
	      if( IChannel == 2 && cand3leptonChannel == "eemu"  ) {
		MyhistoManager.FillHisto(HT_eemu_afterbjetveto, "HT_eemu_afterbjetveto", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoPtZ_eemu_afterbjetveto,       "RecoPtZ_eemu_afterbjetveto",      Zpt  ,     datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto( mWT_eemu_afterbjetveto, "mWT_eemu_afterbjetveto", mTW, datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(Mt_eemu_afterbjetveto, "Mt_eemu_afterbjetveto" , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(deltaPhilj_eemu_afterbjetveto , "deltaPhilj_eemu_afterbjetveto",  deltaPhilj ,    datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(LeptWPt_eemu_afterbjetveto, "LeptWPt_eemu_afterbjetveto", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoTopMass_eemu_afterbjetveto, "RecoTopMass_eemu_afterbjetveto", topCand.M() ,datasetName, IsSignal, Dweight[ITypeMC]);
                for(unsigned int i=0 ; i< selJets.size(); i++){
		  MyhistoManager.FillHisto(JetPt_eemu_afterbjetveto, "JetPt_eemu_afterbjetveto",
	          selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		  MyhistoManager.FillHisto(JetEta_eemu_afterbjetveto, "JetEta_eemu_afterbjetveto",
	          selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
		} 
	      }
	      if( IChannel == 3 && cand3leptonChannel == "eee"   ) {
		MyhistoManager.FillHisto(HT_eee_afterbjetveto, "HT_eee_afterbjetveto", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoPtZ_eee_afterbjetveto,        "RecoPtZ_eee_afterbjetveto",        Zpt ,     datasetName, IsSignal, Dweight[ITypeMC]);   
		MyhistoManager.FillHisto( mWT_eee_afterbjetveto, "mWT_eee_afterbjetveto", mTW, datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(Mt_eee_afterbjetveto, "Mt_eee_afterbjetveto" , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(deltaPhilj_eee_afterbjetveto , "deltaPhilj_eee_afterbjetveto",  deltaPhilj ,    datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(LeptWPt_eee_afterbjetveto, "LeptWPt_eee_afterbjetveto", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoTopMass_eee_afterbjetveto, "RecoTopMass_eee_afterbjetveto", topCand.M() ,datasetName, IsSignal, Dweight[ITypeMC]);
                for(unsigned int i=0 ; i< selJets.size(); i++){
		  MyhistoManager.FillHisto(JetPt_eee_afterbjetveto, "JetPt_eee_afterbjetveto",
	          selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		  MyhistoManager.FillHisto(JetEta_eee_afterbjetveto, "JetEta_eee_afterbjetveto",
	          selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
		}
	      }
	    }
            
	    
	    
	    //------------------------------------------
	    //------------------------------------------
	    //--------- starts b-tag Data --------------
	    //------------------------------------------
	    //------------------------------------------
	      double btagdiscri = selJets[0].bTag["combinedSecondaryVertexBJetTags"];
	      if(btagdiscri<0) btagdiscri = 0;
	    if( abs(selJets[0].partonFlavour) == 5){
	      if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(BJetDiscri_mumumu_afterjetsel_bjets, "BJetDiscri_mumumu_afterjetsel_bjets" ,  btagdiscri,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(BJetDiscri_mumue_afterjetsel_bjets , "BJetDiscri_mumue_afterjetsel_bjets"  ,  btagdiscri,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(BJetDiscri_eemu_afterjetsel_bjets  , "BJetDiscri_eemu_afterjetsel_bjets"   ,  btagdiscri,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(BJetDiscri_eee_afterjetsel_bjets   , "BJetDiscri_eee_afterjetsel_bjets"    ,  btagdiscri,datasetName, IsSignal, Dweight[ITypeMC]);
	     }
	    if( abs(selJets[0].partonFlavour) ==4 ){
	      if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(BJetDiscri_mumumu_afterjetsel_cjets, "BJetDiscri_mumumu_afterjetsel_cjets" ,  btagdiscri,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(BJetDiscri_mumue_afterjetsel_cjets , "BJetDiscri_mumue_afterjetsel_cjets"  ,  btagdiscri,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(BJetDiscri_eemu_afterjetsel_cjets  , "BJetDiscri_eemu_afterjetsel_cjets"   ,  btagdiscri,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(BJetDiscri_eee_afterjetsel_cjets   , "BJetDiscri_eee_afterjetsel_cjets"    ,  btagdiscri,datasetName, IsSignal, Dweight[ITypeMC]);
	     }
	     
	     if( abs(selJets[0].partonFlavour) <4 ){
	      if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(BJetDiscri_mumumu_afterjetsel_ljets, "BJetDiscri_mumumu_afterjetsel_ljets" ,  btagdiscri,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(BJetDiscri_mumue_afterjetsel_ljets , "BJetDiscri_mumue_afterjetsel_ljets"  ,  btagdiscri,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(BJetDiscri_eemu_afterjetsel_ljets  , "BJetDiscri_eemu_afterjetsel_ljets"   ,  btagdiscri,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(BJetDiscri_eee_afterjetsel_ljets   , "BJetDiscri_eee_afterjetsel_ljets"    ,  btagdiscri,datasetName, IsSignal, Dweight[ITypeMC]);
	     }
	     
	     
	     
	     
	     
	    if(
	      (NBtaggedJets >=0  &&  NBtaggedJets <=1 && isData && !useNonIsoWcand ) ||
	      (NBtaggedJets >=0  &&  NBtaggedJets <=1 && isData && useNonIsoWcand && theMET < themetcut)
	      
	      
	      //&&
              //fabs(lept3.DeltaPhi(selJets[idxBtag].p4))  > 0.5 &&
              //lept3.DeltaR(  selJets[idxBtag].p4) > 0.5
	      ){
	    
	    
	    
	   // if(!isData)  Dweight[ITypeMC]*= (1-weightb[3]-weightb[4]);
	    
	    
	    
	    	 
	    if(IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto( mWT_mumumu_afterbjetsel, "mWT_mumumu_afterbjetsel", mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	    if(IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto( mWT_mumue_afterbjetsel,  "mWT_mumue_afterbjetsel" , mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	    if(IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto( mWT_eemu_afterbjetsel,   "mWT_eemu_afterbjetsel"  , mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	    if(IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto( mWT_eee_afterbjetsel,    "mWT_eee_afterbjetsel"   , mTW, datasetName, IsSignal, Dweight[ITypeMC]);

	
	    double nlept = selMuons.size()+selElectrons.size();
	    if(nlept>=4) nlept = 4;   
	      	  	   
	      if( IChannel == 0 && cand3leptonChannel == "mumumu"){
  		MyhistoManager.FillHisto(NJet_mumumu_afterbsel , "NJet_mumumu_afterbsel",NSeljets, datasetName, IsSignal, Dweight[ITypeMC]);
  		MyhistoManager.FillHisto(NLept_mumumu_afterbsel , "NLept_mumumu_afterbsel",nlept, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(InvM_ll_mumumu_afterbjetsel, "InvM_ll_mumumu_afterbjetsel", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(LeptZPt_mumumu_afterbjetsel, "LeptZPt_mumumu_afterbjetsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptZPt_mumumu_afterbjetsel, "LeptZPt_mumumu_afterbjetsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptWPt_mumumu_afterbjetsel, "LeptWPt_mumumu_afterbjetsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(HT_mumumu_afterbjetsel, "HT_mumumu_afterbjetsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(MET_mumumu_afterbjetsel, "MET_mumumu_afterbjetsel",theMET , datasetName, IsSignal, Dweight[ITypeMC]);
                for(unsigned int i=0 ; i< selJets.size(); i++){
		  MyhistoManager.FillHisto(JetPt_mumumu_afterbjetsel, "JetPt_mumumu_afterbjetsel",
	          selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		  MyhistoManager.FillHisto(JetEta_mumumu_afterbjetsel, "JetEta_mumumu_afterbjetsel",
	          selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
		}

	  
	      }
	      if( IChannel == 1 && cand3leptonChannel == "mumue" ){
	        MyhistoManager.FillHisto(NJet_mumue_afterbsel  , "NJet_mumue_afterbsel" ,NSeljets, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(NLept_mumue_afterbsel  , "NLept_mumue_afterbsel" ,nlept, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(InvM_ll_mumue_afterbjetsel, "InvM_ll_mumue_afterbjetsel", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(LeptZPt_mumue_afterbjetsel, "LeptZPt_mumue_afterbjetsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptZPt_mumue_afterbjetsel, "LeptZPt_mumue_afterbjetsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptWPt_mumue_afterbjetsel, "LeptWPt_mumue_afterbjetsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(HT_mumue_afterbjetsel, "HT_mumue_afterbjetsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(MET_mumue_afterbjetsel, "MET_mumue_afterbjetsel", theMET, datasetName, IsSignal, Dweight[ITypeMC]);
                for(unsigned int i=0 ; i< selJets.size(); i++){
		  MyhistoManager.FillHisto(JetPt_mumue_afterbjetsel, "JetPt_mumue_afterbjetsel",
	          selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		  MyhistoManager.FillHisto(JetEta_mumue_afterbjetsel, "JetEta_mumue_afterbjetsel",
	          selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
		}


	      }
	      if( IChannel == 2 && cand3leptonChannel == "eemu"  ){
	        MyhistoManager.FillHisto(NJet_eemu_afterbsel   , "NJet_eemu_afterbsel"  ,NSeljets, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(NLept_eemu_afterbsel   , "NLept_eemu_afterbsel"  ,nlept, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(InvM_ll_eemu_afterbjetsel, "InvM_ll_eemu_afterbjetsel", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(LeptZPt_eemu_afterbjetsel, "LeptZPt_eemu_afterbjetsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptZPt_eemu_afterbjetsel, "LeptZPt_eemu_afterbjetsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptWPt_eemu_afterbjetsel, "LeptWPt_eemu_afterbjetsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(HT_eemu_afterbjetsel, "HT_eemu_afterbjetsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(MET_eemu_afterbjetsel, "MET_eemu_afterbjetsel", theMET, datasetName, IsSignal, Dweight[ITypeMC]);
                for(unsigned int i=0 ; i< selJets.size(); i++){
		  MyhistoManager.FillHisto(JetPt_eemu_afterbjetsel, "JetPt_eemu_afterbjetsel",
	          selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		  MyhistoManager.FillHisto(JetEta_eemu_afterbjetsel, "JetEta_eemu_afterbjetsel",
	          selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
		}


	      }
	    
	      if( IChannel == 3 && cand3leptonChannel == "eee"   ){
		MyhistoManager.FillHisto(NJet_eee_afterbsel    , "NJet_eee_afterbsel"   ,NSeljets, datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(NLept_eee_afterbsel    , "NLept_eee_afterbsel"   ,nlept, datasetName, IsSignal, Dweight[ITypeMC]);
  	        MyhistoManager.FillHisto(InvM_ll_eee_afterbjetsel, "InvM_ll_eee_afterbjetsel", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);	    
	        MyhistoManager.FillHisto(LeptZPt_eee_afterbjetsel, "LeptZPt_eee_afterbjetsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptZPt_eee_afterbjetsel, "LeptZPt_eee_afterbjetsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptWPt_eee_afterbjetsel, "LeptWPt_eee_afterbjetsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(HT_eee_afterbjetsel, "HT_eee_afterbjetsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(MET_eee_afterbjetsel, "MET_eee_afterbjetsel", theMET, datasetName, IsSignal, Dweight[ITypeMC]);
                for(unsigned int i=0 ; i< selJets.size(); i++){
		  MyhistoManager.FillHisto(JetPt_eee_afterbjetsel, "JetPt_eee_afterbjetsel",
	          selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		  MyhistoManager.FillHisto(JetEta_eee_afterbjetsel, "JetEta_eee_afterbjetsel",
	          selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
		}


	      }  	  	

              //if( IChannel == 0 && cand3leptonChannel == "mumumu")      
                //   cout<<"RUN "<<event->general.runNb<<" EVT "<<event->general.eventNb<<endl;
	      
 	      if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(CutFlow_mumumu, "CutFlow_mumumu", 4, datasetName, IsSignal, Dweight[ITypeMC]);
 	      if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(CutFlow_mumue,  "CutFlow_mumue" , 4, datasetName, IsSignal, Dweight[ITypeMC]);
 	      if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(CutFlow_eemu,   "CutFlow_eemu"  , 4, datasetName, IsSignal, Dweight[ITypeMC]);
 	      if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(CutFlow_eee,    "CutFlow_eee"   , 4, datasetName, IsSignal, Dweight[ITypeMC]);
 	      if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(ErrCutFlow_mumumu,   "ErrCutFlow_mumumu"  , 4, datasetName, IsSignal, EventYieldWeightError);
 	      if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(ErrCutFlow_mumue,    "ErrCutFlow_mumue"   , 4, datasetName, IsSignal, EventYieldWeightError);
 	      if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(ErrCutFlow_eemu,     "ErrCutFlow_eemu"    , 4, datasetName, IsSignal, EventYieldWeightError);
 	      if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(ErrCutFlow_eee,      "ErrCutFlow_eee"     , 4, datasetName, IsSignal, EventYieldWeightError);
	

	      TLorentzVector metP4(met.p2.Px(), met.p2.Py(), 0, sqrt(met.p2.Px()*met.p2.Px() + met.p2.Py()*met.p2.Py()));
	      TLorentzVector transTop = lept3 + selJets[0].p4 + metP4;
	      
	      if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(Mt_mumumu_afterbjetsel, "Mt_mumumu_afterbjetsel" , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(Mt_mumue_afterbjetsel , "Mt_mumue_afterbjetsel"  , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(Mt_eemu_afterbjetsel  , "Mt_eemu_afterbjetsel"   , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(Mt_eee_afterbjetsel   , "Mt_eee_afterbjetsel"    , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
	       //cout << "line 1103" << endl;



	      // Use b jet for top mass computation
	      topCand = neutrino + lept3 + selJets[idxBtag].p4 ;
	      
	      if(dilept.M() > 300){
	        /*cout << " topcand mass " << topCand.M()  << endl;
	        cout << " topcand pT   " << topCand.Pt() << endl;
		cout << " total mass   " << (topCand+dilept).M() << endl;
		cout << " channel   " <<  cand3leptonChannel  << "  IChannel " << IChannel << endl;*/
	      }
	      
	      tree_EvtWeight  = Dweight[ITypeMC];
	      
	      
              tree_topMass    = topCand.M();
              tree_totMass    = (topCand + (lept1+lept2)).M();
              tree_deltaPhilb = fabs(lept3.DeltaPhi(selJets[idxBtag].p4));
              tree_deltaRlb   = lept3.DeltaR(  selJets[idxBtag].p4);
              tree_deltaRTopZ = (lept1+lept2).DeltaR(topCand);
              tree_asym       = lept3_Charge*fabs(lept3.Eta());
              tree_Zpt        = (lept1+lept2).Pt();
              tree_ZEta       = (lept1+lept2).Eta();
              tree_topPt      = topCand.Pt();
              tree_topEta     = topCand.Eta();
	      
	      tree_totPt        = (topCand + (lept1+lept2)).Pt();
	      tree_totEta       = (topCand + (lept1+lept2)).Eta();

  	      tree_deltaRZl     = (lept1+lept2).DeltaR(lept3);
  	      tree_deltaPhiZmet = (lept1+lept2).DeltaPhi(metP4);
  	      tree_btagDiscri   = btagdiscri;
  	      tree_NJets        = float(selJets.size());
	      
	      
	      /*if(!isData){
	        double thrand = rand.Uniform();
		NBtaggedJets = 0;
	        if( thrand > weightb[1] ) NBtaggedJets = 1;
		
	      }*/
	      
	      tree_NBJets       = float(NBtaggedJets);
	     
  	      tree_leptWPt        = lept3.Pt(); 
	      tree_leptWEta       = lept3.Eta();
  	      tree_leadJetPt      = selJets[0].p4.Pt(); 
	      tree_leadJetEta     = selJets[0].p4.Eta();
              tree_deltaRZleptW   = (lept1+lept2).DeltaR(lept3); 
	      tree_deltaPhiZleptW = (lept1+lept2).DeltaPhi(lept3);
	      
	      
	      
	      if(datasetName=="DataMu" || datasetName=="DataEG" || datasetName=="DataMuEG")  tree_SampleType   = 0;
	      if(datasetName=="FCNCkut" )        tree_SampleType   = 1;
	      if(datasetName=="TTbar" )          tree_SampleType   = 2;
	      if(datasetName=="DYToLL_M10-50" )  tree_SampleType   = 3;
	      if(datasetName=="Zjets" )          tree_SampleType   = 4;
	      if(datasetName=="Wjets" )          tree_SampleType   = 5;
	      if(datasetName=="TtW" )            tree_SampleType   = 6;
	      if(datasetName=="TbartW" )         tree_SampleType   = 7;
	      if(datasetName=="WZ" )             tree_SampleType   = 8;
	      if(datasetName=="ZZ" )             tree_SampleType   = 9;
	      if(datasetName=="WW" )             tree_SampleType   = 10;
	      if(datasetName=="FCNCkut" )        tree_SampleType   = 11;
	      if(datasetName=="FCNCkct" )        tree_SampleType   = 12;
	      if(datasetName=="FCNCxut" )        tree_SampleType   = 13;
	      if(datasetName=="FCNCxct" )        tree_SampleType   = 14;
	      if(datasetName=="FCNCzut" )        tree_SampleType   = 15;
	      if(datasetName=="FCNCzct" )        tree_SampleType   = 16;
	      
	      if( IChannel == 0 && cand3leptonChannel == "mumumu") {
	        tree_Channel = 0; TheTree->Fill();
		MyhistoManager.FillHisto(Asym_mumumu_afterbjetsel,       "Asym_mumumu_afterbjetsel",         tree_asym,    datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoPtZ_mumumu_afterbjetsel,     "RecoPtZ_mumumu_afterbjetsel",     tree_Zpt,     datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoTopMass_mumumu_afterbjetsel, "RecoTopMass_mumumu_afterbjetsel", tree_topMass ,datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(deltaPhilb_mumumu_afterbjetsel , "deltaPhilb_mumumu_afterbjetsel",  tree_deltaPhilb ,    datasetName, IsSignal, Dweight[ITypeMC]);

		
	      }
	      if( IChannel == 1 && cand3leptonChannel == "mumue" ) {
	        tree_Channel = 1; TheTree->Fill();
		MyhistoManager.FillHisto(Asym_mumue_afterbjetsel,         "Asym_mumue_afterbjetsel",          tree_asym,    datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoPtZ_mumue_afterbjetsel,      "RecoPtZ_mumue_afterbjetsel",       tree_Zpt,     datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoTopMass_mumue_afterbjetsel,  "RecoTopMass_mumue_afterbjetsel" ,  tree_topMass ,datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(deltaPhilb_mumue_afterbjetsel ,  "deltaPhilb_mumue_afterbjetsel",    tree_deltaPhilb ,    datasetName, IsSignal, Dweight[ITypeMC]);
 
		
	      }
	      if( IChannel == 2 && cand3leptonChannel == "eemu"  ) {
	        tree_Channel = 2; TheTree->Fill();
		MyhistoManager.FillHisto(Asym_eemu_afterbjetsel,          "Asym_eemu_afterbjetsel",           tree_asym,    datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoPtZ_eemu_afterbjetsel,       "RecoPtZ_eemu_afterbjetsel",        tree_Zpt,     datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoTopMass_eemu_afterbjetsel,   "RecoTopMass_eemu_afterbjetsel" ,   tree_topMass ,datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(deltaPhilb_eemu_afterbjetsel ,   "deltaPhilb_eemu_afterbjetsel",     tree_deltaPhilb ,    datasetName, IsSignal, Dweight[ITypeMC]);

	      }
	      if( IChannel == 3 && cand3leptonChannel == "eee"   ) {
	        tree_Channel = 3; TheTree->Fill();
		MyhistoManager.FillHisto(Asym_eee_afterbjetsel,           "Asym_eee_afterbjetsel",            tree_asym,    datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoPtZ_eee_afterbjetsel,        "RecoPtZ_eee_afterbjetsel",         tree_Zpt,     datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(RecoTopMass_eee_afterbjetsel,    "RecoTopMass_eee_afterbjetsel" ,    tree_topMass ,datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(deltaPhilb_eee_afterbjetsel ,    "deltaPhilb_eee_afterbjetsel",      tree_deltaPhilb ,    datasetName, IsSignal, Dweight[ITypeMC]);
	   
		
	      }
	      
  
	      
	      
	      
	      
	      
	      
	      if(transTop.Mt() > 150){

                if( IChannel == 0 && cand3leptonChannel == "mumumu")
                   //cout<<"RUN "<<event->general.runNb<<" EVT "<<event->general.eventNb<<endl;

		if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(CutFlow_mumumu, "CutFlow_mumumu", 5, datasetName, IsSignal, Dweight[ITypeMC]);
 	        if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(CutFlow_mumue,  "CutFlow_mumue" , 5, datasetName, IsSignal, Dweight[ITypeMC]);
 	        if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(CutFlow_eemu,   "CutFlow_eemu"  , 5, datasetName, IsSignal, Dweight[ITypeMC]);
 	        if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(CutFlow_eee,    "CutFlow_eee"   , 5, datasetName, IsSignal, Dweight[ITypeMC]);
 	        if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(ErrCutFlow_mumumu,   "ErrCutFlow_mumumu"  , 5, datasetName, IsSignal, EventYieldWeightError);
 	        if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(ErrCutFlow_mumue,    "ErrCutFlow_mumue"   , 5, datasetName, IsSignal, EventYieldWeightError);
 	        if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(ErrCutFlow_eemu,     "ErrCutFlow_eemu"    , 5, datasetName, IsSignal, EventYieldWeightError);
 	        if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(ErrCutFlow_eee,      "ErrCutFlow_eee"     , 5, datasetName, IsSignal, EventYieldWeightError);
	      }
	    } // end selection btag Data
	    
	    
	    
	     
	    
	    //----------------------------------------
	    //----------------------------------------
	    //--------- starts b-tag MC --------------
	    //----------------------------------------
	    //----------------------------------------
	    
	    
	      //cout << " weight btag btag sel " << weightb[0] << " " << " p0jet " << weightb[1] << " p1jet " << weightb[2] << " p2jet " << weightb[3]
   //		   << " prest " << weightb[4] << " check1=" << weightb[1]+weightb[2]+weightb[3]+weightb[4]  << endl;
 	
	    
	    if(!isData){
	    
	      double sumweight = weightb[1];
	      double therand = rand.Uniform();
	      
		   double btagCorrDiscri = 0;
	      //fixme
	      if(!doBTagCVScorr){
	        if( weightb[1] > 0 && therand  <  weightb[1]                                     ) {NBtaggedJets = 0;}
	        sumweight+=weightb[2];
	        if( weightb[2] > 0 && therand  >= (sumweight-weightb[2]) && therand < sumweight  ) {NBtaggedJets = 1;}
	        sumweight+=weightb[3];
	        if( weightb[3] > 0 && therand  >= (sumweight-weightb[3]) && therand < sumweight  ) {NBtaggedJets = 2;}
	        sumweight+=weightb[4];
	        if( weightb[4] > 0 && therand  >= (sumweight-weightb[4])                         )  NBtaggedJets = 3;
		
		btagCorrDiscri = functionDiscriLjet(selJets[0].bTag["combinedSecondaryVertexBJetTags"]);
		
	      }else{
	         NBtaggedJets = 0;
		 for(unsigned int ijet = 0; ijet < selJets.size(); ijet++){
    	         //cout << "jet flavor " << int(selJets[ijet].partonFlavour) << endl;
    	           if(abs(selJets[ijet].partonFlavour)==5 ){
	           // cout << "found b jet " << int(selJets[ijet].partonFlavour) << endl;
	             foundASelBjet = true;
	           }
		   
		   
    	           if ( AlgoBtag==0 &&  selJets[ijet].bTag["trackCountingHighEffBJetTags"]	  >= btagDiscriCut)  NBtaggedJets++;
    	           if ( AlgoBtag==1 &&  selJets[ijet].bTag["simpleSecondaryVertexHighEffBJetTags"] >= btagDiscriCut) NBtaggedJets++;
    	           if ( AlgoBtag==2 &&  selJets[ijet].bTag["trackCountingHighPurBJetTags"]	  >= btagDiscriCut)  NBtaggedJets++;
    	           if ( AlgoBtag==3 &&  selJets[ijet].bTag["simpleSecondaryVertexHighPurBJetTags"] >= btagDiscriCut) NBtaggedJets++;
    	           if ( AlgoBtag==4 &&  selJets[ijet].bTag["jetProbabilityBJetTags"]  	          >= btagDiscriCut)  NBtaggedJets++;
    	           if ( AlgoBtag==5 &&  selJets[ijet].bTag["jetBProbabilityBJetTags"] 	          >= btagDiscriCut)  NBtaggedJets++;
		   
		   
		   if(abs(selJets[ijet].partonFlavour) == 5){ 
		      if(doBTagCSV_syst == 0) btagCorrDiscri = functionDiscriBjet(     selJets[ijet].bTag["combinedSecondaryVertexBJetTags"]);
		      if(doBTagCSV_syst == 1) btagCorrDiscri = functionDiscriBjet_up(  selJets[ijet].bTag["combinedSecondaryVertexBJetTags"]);
		      if(doBTagCSV_syst == 2) btagCorrDiscri = functionDiscriBjet_down(selJets[ijet].bTag["combinedSecondaryVertexBJetTags"]);
				}
		   if(abs(selJets[ijet].partonFlavour) == 4) {
		      if(doBTagCSV_syst == 0) btagCorrDiscri = functionDiscriCjet(     selJets[ijet].bTag["combinedSecondaryVertexBJetTags"]);
		      if(doBTagCSV_syst == 1) btagCorrDiscri = functionDiscriCjet_up(  selJets[ijet].bTag["combinedSecondaryVertexBJetTags"]);
		      if(doBTagCSV_syst == 2) btagCorrDiscri = functionDiscriCjet_down(selJets[ijet].bTag["combinedSecondaryVertexBJetTags"]);
				}
		   if(abs(selJets[ijet].partonFlavour) < 4) {
		     if(doBTagCSV_syst == 0) btagCorrDiscri = functionDiscriLjet(     selJets[ijet].bTag["combinedSecondaryVertexBJetTags"]);
		     if(doBTagCSV_syst == 1) btagCorrDiscri = functionDiscriLjet_up(  selJets[ijet].bTag["combinedSecondaryVertexBJetTags"]);
		     if(doBTagCSV_syst == 2) btagCorrDiscri = functionDiscriLjet_down(selJets[ijet].bTag["combinedSecondaryVertexBJetTags"]);
		   }
		 		
		 		
    	           if ( AlgoBtag==6 &&  btagCorrDiscri      >= btagDiscriCut){
   		    NBtaggedJets++;
   		    idxBtag = ijet;
    	           }
   	         } 
	       }
	      
	      
	      
	      if(abs(selJets[0].partonFlavour) == 5   && selJets.size() == 1){
	          MyhistoManager.FillHisto(NBJet_mumumu_afterjetsel_bjets ,  "NBJet_mumumu_afterjetsel_bjets",0, datasetName, IsSignal, weightb[1]*Dweight[ITypeMC]);
	          MyhistoManager.FillHisto(NBJet_mumumu_afterjetsel_bjets ,  "NBJet_mumumu_afterjetsel_bjets",1, datasetName, IsSignal, weightb[2]*Dweight[ITypeMC]);	
		}
		if(abs(selJets[0].partonFlavour) == 4 && selJets.size() == 1){
	          MyhistoManager.FillHisto(NBJet_mumumu_afterjetsel_cjets ,  "NBJet_mumumu_afterjetsel_cjets",0, datasetName, IsSignal, weightb[1]*Dweight[ITypeMC]);
	          MyhistoManager.FillHisto(NBJet_mumumu_afterjetsel_cjets ,  "NBJet_mumumu_afterjetsel_cjets",1, datasetName, IsSignal, weightb[2]*Dweight[ITypeMC]);	
		}
		if(abs(selJets[0].partonFlavour) <=3  && selJets.size() == 1){
	          MyhistoManager.FillHisto(NBJet_mumumu_afterjetsel_ljets ,  "NBJet_mumumu_afterjetsel_ljets",0, datasetName, IsSignal, weightb[1]*Dweight[ITypeMC]);
	          MyhistoManager.FillHisto(NBJet_mumumu_afterjetsel_ljets ,  "NBJet_mumumu_afterjetsel_ljets",1, datasetName, IsSignal, weightb[2]*Dweight[ITypeMC]);	
		}
	   
	      //fixme
	      if(
	        (NBtaggedJets >=0  &&  NBtaggedJets <=1 && !useNonIsoWcand ) ||
		(NBtaggedJets >=0  &&  NBtaggedJets <=1 &&  useNonIsoWcand && theMET < themetcut)
		
		){
	      //if(NBtaggedJets >=0  &&  NBtaggedJets <=100000){
	    
	    
	    //if(!isData)  Dweight[ITypeMC]*= (1-weightb[3]-weightb[4]);
	    
	    
	    if(datasetName=="WZ" && applyWZ_finalSel) Dweight[ITypeMC]*=SF_WZ_finalSel;
	    	 
	    if(IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto( mWT_mumumu_afterbjetsel, "mWT_mumumu_afterbjetsel", mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	    if(IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto( mWT_mumue_afterbjetsel,  "mWT_mumue_afterbjetsel" , mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	    if(IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto( mWT_eemu_afterbjetsel,   "mWT_eemu_afterbjetsel"  , mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	    if(IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto( mWT_eee_afterbjetsel,    "mWT_eee_afterbjetsel"   , mTW, datasetName, IsSignal, Dweight[ITypeMC]);
            
	
	    double nlept = selMuons.size()+selElectrons.size();
	    if(nlept>=4) nlept = 4;   
	      	  	   
	      if( IChannel == 0 && cand3leptonChannel == "mumumu"){
  		MyhistoManager.FillHisto(NJet_mumumu_afterbsel , "NJet_mumumu_afterbsel",NSeljets, datasetName, IsSignal, Dweight[ITypeMC]);
  		MyhistoManager.FillHisto(NLept_mumumu_afterbsel , "NLept_mumumu_afterbsel",nlept, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(InvM_ll_mumumu_afterbjetsel, "InvM_ll_mumumu_afterbjetsel", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(LeptZPt_mumumu_afterbjetsel, "LeptZPt_mumumu_afterbjetsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptZPt_mumumu_afterbjetsel, "LeptZPt_mumumu_afterbjetsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptWPt_mumumu_afterbjetsel, "LeptWPt_mumumu_afterbjetsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(HT_mumumu_afterbjetsel, "HT_mumumu_afterbjetsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(MET_mumumu_afterbjetsel, "MET_mumumu_afterbjetsel",theMET , datasetName, IsSignal, Dweight[ITypeMC]);
                for(unsigned int i=0 ; i< selJets.size(); i++){
		  MyhistoManager.FillHisto(JetPt_mumumu_afterbjetsel, "JetPt_mumumu_afterbjetsel",
	          selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		  MyhistoManager.FillHisto(JetEta_mumumu_afterbjetsel, "JetEta_mumumu_afterbjetsel",
	          selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
		}

	  
	      }
	      if( IChannel == 1 && cand3leptonChannel == "mumue" ){
	        MyhistoManager.FillHisto(NJet_mumue_afterbsel  , "NJet_mumue_afterbsel" ,NSeljets, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(NLept_mumue_afterbsel  , "NLept_mumue_afterbsel" ,nlept, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(InvM_ll_mumue_afterbjetsel, "InvM_ll_mumue_afterbjetsel", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(LeptZPt_mumue_afterbjetsel, "LeptZPt_mumue_afterbjetsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptZPt_mumue_afterbjetsel, "LeptZPt_mumue_afterbjetsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptWPt_mumue_afterbjetsel, "LeptWPt_mumue_afterbjetsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(HT_mumue_afterbjetsel, "HT_mumue_afterbjetsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(MET_mumue_afterbjetsel, "MET_mumue_afterbjetsel", theMET, datasetName, IsSignal, Dweight[ITypeMC]);
                for(unsigned int i=0 ; i< selJets.size(); i++){
		  MyhistoManager.FillHisto(JetPt_mumue_afterbjetsel, "JetPt_mumue_afterbjetsel",
	          selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		  MyhistoManager.FillHisto(JetEta_mumue_afterbjetsel, "JetEta_mumue_afterbjetsel",
	          selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
		}


	      }
	      if( IChannel == 2 && cand3leptonChannel == "eemu"  ){
	        MyhistoManager.FillHisto(NJet_eemu_afterbsel   , "NJet_eemu_afterbsel"  ,NSeljets, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(NLept_eemu_afterbsel   , "NLept_eemu_afterbsel"  ,nlept, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(InvM_ll_eemu_afterbjetsel, "InvM_ll_eemu_afterbjetsel", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(LeptZPt_eemu_afterbjetsel, "LeptZPt_eemu_afterbjetsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptZPt_eemu_afterbjetsel, "LeptZPt_eemu_afterbjetsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptWPt_eemu_afterbjetsel, "LeptWPt_eemu_afterbjetsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(HT_eemu_afterbjetsel, "HT_eemu_afterbjetsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(MET_eemu_afterbjetsel, "MET_eemu_afterbjetsel", theMET, datasetName, IsSignal, Dweight[ITypeMC]);
                for(unsigned int i=0 ; i< selJets.size(); i++){
		  MyhistoManager.FillHisto(JetPt_eemu_afterbjetsel, "JetPt_eemu_afterbjetsel",
	          selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		  MyhistoManager.FillHisto(JetEta_eemu_afterbjetsel, "JetEta_eemu_afterbjetsel",
	          selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
		}


	      }
	    
	      if( IChannel == 3 && cand3leptonChannel == "eee"   ){
		MyhistoManager.FillHisto(NJet_eee_afterbsel    , "NJet_eee_afterbsel"   ,NSeljets, datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(NLept_eee_afterbsel    , "NLept_eee_afterbsel"   ,nlept, datasetName, IsSignal, Dweight[ITypeMC]);
  	        MyhistoManager.FillHisto(InvM_ll_eee_afterbjetsel, "InvM_ll_eee_afterbjetsel", dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);	    
	        MyhistoManager.FillHisto(LeptZPt_eee_afterbjetsel, "LeptZPt_eee_afterbjetsel", lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptZPt_eee_afterbjetsel, "LeptZPt_eee_afterbjetsel", lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(LeptWPt_eee_afterbjetsel, "LeptWPt_eee_afterbjetsel", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(HT_eee_afterbjetsel, "HT_eee_afterbjetsel", sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(MET_eee_afterbjetsel, "MET_eee_afterbjetsel", theMET, datasetName, IsSignal, Dweight[ITypeMC]);
                for(unsigned int i=0 ; i< selJets.size(); i++){
		  MyhistoManager.FillHisto(JetPt_eee_afterbjetsel, "JetPt_eee_afterbjetsel",
	          selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		  MyhistoManager.FillHisto(JetEta_eee_afterbjetsel, "JetEta_eee_afterbjetsel",
	          selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
		}


	      }  	  	

              //if( IChannel == 0 && cand3leptonChannel == "mumumu")      
                //   cout<<"RUN "<<event->general.runNb<<" EVT "<<event->general.eventNb<<endl;
		
	      if(doPDF && !isData){
	      
    	         UncertaintyType result = pdf.Calculate(event->mc); 
    	         std::cout << "weight : mean=" << result.Mean << " - max=" << result.Max << " - min=" << result.Min << std::endl; 

    	         /*if(pdftype == 0)  Dweight[ITypeMC]= Dweight[ITypeMC]*(result.Min/result.Mean);
    	         if(pdftype == 1)  Dweight[ITypeMC]= Dweight[ITypeMC]*(result.Max/result.Mean);*/
		 
    	         if(pdftype == 0 && result.Mean!=0)  Dweight[ITypeMC]= Dweight[ITypeMC]*(result.Min/result.Mean);
    	         if(pdftype == 1 && result.Mean!=0)  Dweight[ITypeMC]= Dweight[ITypeMC]*(result.Max/result.Mean);
		 
 	      }  

 	      if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(CutFlow_mumumu, "CutFlow_mumumu", 4, datasetName, IsSignal, Dweight[ITypeMC]);
 	      if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(CutFlow_mumue,  "CutFlow_mumue" , 4, datasetName, IsSignal, Dweight[ITypeMC]);
 	      if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(CutFlow_eemu,   "CutFlow_eemu"  , 4, datasetName, IsSignal, Dweight[ITypeMC]);
 	      if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(CutFlow_eee,    "CutFlow_eee"   , 4, datasetName, IsSignal, Dweight[ITypeMC]);
 	      if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(ErrCutFlow_mumumu,   "ErrCutFlow_mumumu"  , 4, datasetName, IsSignal, EventYieldWeightError);
 	      if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(ErrCutFlow_mumue,    "ErrCutFlow_mumue"   , 4, datasetName, IsSignal, EventYieldWeightError);
 	      if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(ErrCutFlow_eemu,     "ErrCutFlow_eemu"    , 4, datasetName, IsSignal, EventYieldWeightError);
 	      if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(ErrCutFlow_eee,      "ErrCutFlow_eee"     , 4, datasetName, IsSignal, EventYieldWeightError);
	

	      TLorentzVector metP4(met.p2.Px(), met.p2.Py(), 0, sqrt(met.p2.Px()*met.p2.Px() + met.p2.Py()*met.p2.Py()));
	      TLorentzVector transTop = lept3 + selJets[0].p4 + metP4;
	      
	      if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(Mt_mumumu_afterbjetsel, "Mt_mumumu_afterbjetsel" , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(Mt_mumue_afterbjetsel , "Mt_mumue_afterbjetsel"  , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(Mt_eemu_afterbjetsel  , "Mt_eemu_afterbjetsel"   , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
              if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(Mt_eee_afterbjetsel   , "Mt_eee_afterbjetsel"    , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
	       //cout << "line 1103" << endl;
	       
	  
	      // Use b jet for top mass computation
	      topCand = neutrino + lept3 + selJets[idxBtag].p4 ;
	      std::cout << "*******************  evt nb " << event->general.eventNb  <<  cand3leptonChannel  << "  IChannel " << IChannel << std::endl;
	      std::cout << "x1 " <<  event->mc.x.first <<   "  x2 " <<  event->mc.x.second  << " Q_scale " << event->mc.Q_scale << endl;
		
	      //cout << "Dweight[ITypeMC] 1 " << Dweight[ITypeMC] << endl;
	      tree_EvtWeight  = Dweight[ITypeMC];
	      
	      //cout << "Dweight[ITypeMC] 2 " << Dweight[ITypeMC] << endl;
	      
              tree_topMass    = topCand.M();
              tree_totMass    = (topCand + (lept1+lept2)).M();
              tree_deltaPhilb = fabs(lept3.DeltaPhi(selJets[idxBtag].p4));
              tree_deltaRlb   = lept3.DeltaR(  selJets[idxBtag].p4);
              tree_deltaRTopZ = (lept1+lept2).DeltaR(topCand);
              tree_asym       = lept3_Charge*fabs(lept3.Eta());
              tree_Zpt        = (lept1+lept2).Pt();
              tree_ZEta       = (lept1+lept2).Eta();
              tree_topPt      = topCand.Pt();
              tree_topEta     = topCand.Eta();
	      
	      tree_totPt        = (topCand + (lept1+lept2)).Pt();
	      tree_totEta       = (topCand + (lept1+lept2)).Eta();

  	      tree_deltaRZl     = (lept1+lept2).DeltaR(lept3);
  	      tree_deltaPhiZmet = (lept1+lept2).DeltaPhi(metP4);
  	      tree_btagDiscri   = btagCorrDiscri;//btagdiscri; //selJets[0].bTag["combinedSecondaryVertexBJetTags"];
  	      tree_NJets        = float(selJets.size());
	      
	      
	      /*if(!isData){
	        double thrand = rand.Uniform();
		NBtaggedJets = 0;
	        if( thrand > weightb[1] ) NBtaggedJets = 1;
		
	      }*/
	      
	      tree_NBJets       = float(NBtaggedJets);
	     
  	      tree_leptWPt        = lept3.Pt(); 
	      tree_leptWEta       = lept3.Eta();
  	      tree_leadJetPt      = selJets[0].p4.Pt(); 
	      tree_leadJetEta     = selJets[0].p4.Eta();
              tree_deltaRZleptW   = (lept1+lept2).DeltaR(lept3); 
	      tree_deltaPhiZleptW = (lept1+lept2).DeltaPhi(lept3);
	      
	      
	      
	      if(datasetName=="DataMu" || datasetName=="DataEG" || datasetName=="DataMuEG")  tree_SampleType   = 0;
	      if(datasetName=="TTbar" )          tree_SampleType   = 2;
	      if(datasetName=="DYToLL_M10-50" )  tree_SampleType   = 3;
	      if(datasetName=="Zjets" )          tree_SampleType   = 4;
	      if(datasetName=="Wjets" )          tree_SampleType   = 5;
	      if(datasetName=="TtW" )            tree_SampleType   = 6;
	      if(datasetName=="TbartW" )         tree_SampleType   = 7;
	      if(datasetName=="WZ" )             tree_SampleType   = 8;
	      if(datasetName=="ZZ" )             tree_SampleType   = 9;
	      if(datasetName=="WW" )             tree_SampleType   = 10;
	      if(datasetName=="FCNCkut" )        tree_SampleType   = 11;
	      if(datasetName=="FCNCkct" )        tree_SampleType   = 12;
	      if(datasetName=="FCNCxut" )        tree_SampleType   = 13;
	      if(datasetName=="FCNCxct" )        tree_SampleType   = 14;
	      if(datasetName=="FCNCzut" )        tree_SampleType   = 15;
	      if(datasetName=="FCNCzct" )        tree_SampleType   = 16;
	      
	      if( IChannel == 0 && cand3leptonChannel == "mumumu") {
	        tree_Channel = 0; TheTree->Fill();
		MyhistoManager.FillHisto(Asym_mumumu_afterbjetsel,       "Asym_mumumu_afterbjetsel",         tree_asym,    datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoPtZ_mumumu_afterbjetsel,     "RecoPtZ_mumumu_afterbjetsel",     tree_Zpt,     datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoTopMass_mumumu_afterbjetsel, "RecoTopMass_mumumu_afterbjetsel", tree_topMass ,datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(deltaPhilb_mumumu_afterbjetsel , "deltaPhilb_mumumu_afterbjetsel",  tree_deltaPhilb ,    datasetName, IsSignal, Dweight[ITypeMC]);

		
	      }
	      if( IChannel == 1 && cand3leptonChannel == "mumue" ) {
	        tree_Channel = 1; TheTree->Fill();
		MyhistoManager.FillHisto(Asym_mumue_afterbjetsel,         "Asym_mumue_afterbjetsel",          tree_asym,    datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoPtZ_mumue_afterbjetsel,      "RecoPtZ_mumue_afterbjetsel",       tree_Zpt,     datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoTopMass_mumue_afterbjetsel,  "RecoTopMass_mumue_afterbjetsel" ,  tree_topMass ,datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(deltaPhilb_mumue_afterbjetsel ,  "deltaPhilb_mumue_afterbjetsel",    tree_deltaPhilb ,    datasetName, IsSignal, Dweight[ITypeMC]);
 
		
	      }
	      if( IChannel == 2 && cand3leptonChannel == "eemu"  ) {
	        tree_Channel = 2; TheTree->Fill();
		MyhistoManager.FillHisto(Asym_eemu_afterbjetsel,          "Asym_eemu_afterbjetsel",           tree_asym,    datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoPtZ_eemu_afterbjetsel,       "RecoPtZ_eemu_afterbjetsel",        tree_Zpt,     datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoTopMass_eemu_afterbjetsel,   "RecoTopMass_eemu_afterbjetsel" ,   tree_topMass ,datasetName, IsSignal, Dweight[ITypeMC]);
                MyhistoManager.FillHisto(deltaPhilb_eemu_afterbjetsel ,   "deltaPhilb_eemu_afterbjetsel",     tree_deltaPhilb ,    datasetName, IsSignal, Dweight[ITypeMC]);

	      }
	      if( IChannel == 3 && cand3leptonChannel == "eee"   ) {
	        tree_Channel = 3; TheTree->Fill();
		MyhistoManager.FillHisto(Asym_eee_afterbjetsel,           "Asym_eee_afterbjetsel",            tree_asym,    datasetName, IsSignal, Dweight[ITypeMC]);
		MyhistoManager.FillHisto(RecoPtZ_eee_afterbjetsel,        "RecoPtZ_eee_afterbjetsel",         tree_Zpt,     datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(RecoTopMass_eee_afterbjetsel,    "RecoTopMass_eee_afterbjetsel" ,    tree_topMass ,datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(deltaPhilb_eee_afterbjetsel ,    "deltaPhilb_eee_afterbjetsel",      tree_deltaPhilb ,    datasetName, IsSignal, Dweight[ITypeMC]);
	   
		
	      }
	      
  
	      
	      
	      
	      
	      
	      
	      if(transTop.Mt() > 150){

                if( IChannel == 0 && cand3leptonChannel == "mumumu")
                   //cout<<"RUN "<<event->general.runNb<<" EVT "<<event->general.eventNb<<endl;

		if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(CutFlow_mumumu, "CutFlow_mumumu", 5, datasetName, IsSignal, Dweight[ITypeMC]);
 	        if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(CutFlow_mumue,  "CutFlow_mumue" , 5, datasetName, IsSignal, Dweight[ITypeMC]);
 	        if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(CutFlow_eemu,   "CutFlow_eemu"  , 5, datasetName, IsSignal, Dweight[ITypeMC]);
 	        if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(CutFlow_eee,    "CutFlow_eee"   , 5, datasetName, IsSignal, Dweight[ITypeMC]);
 	        if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(ErrCutFlow_mumumu,   "ErrCutFlow_mumumu"  , 5, datasetName, IsSignal, EventYieldWeightError);
 	        if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(ErrCutFlow_mumue,    "ErrCutFlow_mumue"   , 5, datasetName, IsSignal, EventYieldWeightError);
 	        if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(ErrCutFlow_eemu,     "ErrCutFlow_eemu"    , 5, datasetName, IsSignal, EventYieldWeightError);
 	        if( IChannel == 3 && cand3leptonChannel == "eee"   ) MyhistoManager.FillHisto(ErrCutFlow_eee,      "ErrCutFlow_eee"     , 5, datasetName, IsSignal, EventYieldWeightError);
	      }
	      }
	    } // end selection btag Data
	    
	    
	    
	    
	    
	  } // end selection njet
	} //end selection Z cand    
      } // end selection 3 leptons
    } // end selection trigger
  } //end loops on datasets
  
  
  return kTRUE;
}

//_____________________________________________________________________________
void ProofSelectorMyCutFlow::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
  if(fProofFile) fProofFile->Print();
  if (fFile) {
    Bool_t cleanup = kFALSE;
    TDirectory *savedir = gDirectory;
    fFile->cd();
    
  MyhistoManager.WriteMyHisto(SelABjet, "all");
  MyhistoManager.WriteMyHisto(SelABjet_afterjetsel, "all");
  MyhistoManager.WriteMyHisto(Ntrilept_mumumu, "all");
  MyhistoManager.WriteMyHisto(Ntrileptnoniso_mumumu, "all");
  
   
  MyhistoManager.WriteMyHisto(Nvertex, "all");  
  
  MyhistoManager.WriteMyHisto(CutFlow_mumumu, "all" );
  MyhistoManager.WriteMyHisto(CutFlow_mumue,  "all" );
  MyhistoManager.WriteMyHisto(CutFlow_eemu,   "all" );
  MyhistoManager.WriteMyHisto(CutFlow_eee,    "all" );
  
  
  
  MyhistoManager.WriteMyHisto(PU_before_mumumu, "all");  
  MyhistoManager.WriteMyHisto(PU_before_mumue, "all");  
  MyhistoManager.WriteMyHisto(PU_before_eemu, "all");  
  MyhistoManager.WriteMyHisto(PU_before_eee, "all");  

  MyhistoManager.WriteMyHisto(PU_intime_mumumu, "all");  
  MyhistoManager.WriteMyHisto(PU_intime_mumue, "all");  
  MyhistoManager.WriteMyHisto(PU_intime_eemu, "all");  
  MyhistoManager.WriteMyHisto(PU_intime_eee, "all");  

  MyhistoManager.WriteMyHisto(PU_after_mumumu, "all");  
  MyhistoManager.WriteMyHisto(PU_after_mumue, "all");  
  MyhistoManager.WriteMyHisto(PU_after_eemu, "all");  
  MyhistoManager.WriteMyHisto(PU_after_eee, "all");  
  
  
  MyhistoManager.WriteMyHisto(NVtx_mumumu, "all");  
  MyhistoManager.WriteMyHisto(NVtx_mumue, "all");  
  MyhistoManager.WriteMyHisto(NVtx_eemu, "all");  
  MyhistoManager.WriteMyHisto(NVtx_eee, "all");  
  
  MyhistoManager.WriteMyHisto(Nvtx_mumumu_afterleptsel, "all" );
  MyhistoManager.WriteMyHisto(Nvtx_mumue_afterleptsel , "all" );
  MyhistoManager.WriteMyHisto(Nvtx_eemu_afterleptsel  , "all" );
  MyhistoManager.WriteMyHisto(Nvtx_eee_afterleptsel   , "all" );
  
  
  MyhistoManager.WriteMyHisto(ErrCutFlow_mumumu,  "all");
  MyhistoManager.WriteMyHisto(ErrCutFlow_mumue,   "all");
  MyhistoManager.WriteMyHisto(ErrCutFlow_eemu,    "all");
  MyhistoManager.WriteMyHisto(ErrCutFlow_eee,     "all");
  
  
  
  MyhistoManager.WriteMyHisto(Mt_mumumu_afterbjetsel, "all"); 
  MyhistoManager.WriteMyHisto(Mt_mumue_afterbjetsel , "all");
  MyhistoManager.WriteMyHisto(Mt_eemu_afterbjetsel  , "all");
  MyhistoManager.WriteMyHisto(Mt_eee_afterbjetsel   , "all");
    
  MyhistoManager.WriteMyHisto(Mt_mumumu_afterbjetveto, "all"); 
  MyhistoManager.WriteMyHisto(Mt_mumue_afterbjetveto , "all");
  MyhistoManager.WriteMyHisto(Mt_eemu_afterbjetveto  , "all");
  MyhistoManager.WriteMyHisto(Mt_eee_afterbjetveto   , "all");
    
    
  
  MyhistoManager.WriteMyHisto(NJet_mumumu_afterZsel,"all");
  MyhistoManager.WriteMyHisto(NJet_mumue_afterZsel, "all");
  MyhistoManager.WriteMyHisto(NJet_eemu_afterZsel,  "all");
  MyhistoManager.WriteMyHisto(NJet_eee_afterZsel,   "all");
  
  
  MyhistoManager.WriteMyHisto(NJet_mumumu_afterbsel,"all");
  MyhistoManager.WriteMyHisto(NJet_mumue_afterbsel, "all");
  MyhistoManager.WriteMyHisto(NJet_eemu_afterbsel,  "all");
  MyhistoManager.WriteMyHisto(NJet_eee_afterbsel,   "all");
  
  
  MyhistoManager.WriteMyHisto(NLept_mumumu_afterbsel,"all");
  MyhistoManager.WriteMyHisto(NLept_mumue_afterbsel, "all");
  MyhistoManager.WriteMyHisto(NLept_eemu_afterbsel,  "all");
  MyhistoManager.WriteMyHisto(NLept_eee_afterbsel,   "all");
  
  MyhistoManager.WriteMyHisto(NBJet_mumumu_afterZsel, "all");
  MyhistoManager.WriteMyHisto(NBJet_mumue_afterZsel,  "all");
  MyhistoManager.WriteMyHisto(NBJet_eemu_afterZsel,   "all");
  MyhistoManager.WriteMyHisto(NBJet_eee_afterZsel,    "all");
  
  
  MyhistoManager.WriteMyHisto(NBJet_mumumu_afterjetsel, "all");
  MyhistoManager.WriteMyHisto(NBJet_mumue_afterjetsel,  "all");
  MyhistoManager.WriteMyHisto(NBJet_eemu_afterjetsel,   "all");
  MyhistoManager.WriteMyHisto(NBJet_eee_afterjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(InvM_ll_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto(InvM_ll_mumue_afterleptsel,  "all");
  MyhistoManager.WriteMyHisto(InvM_ll_eemu_afterleptsel,   "all");
  MyhistoManager.WriteMyHisto(InvM_ll_eee_afterleptsel,    "all");
  
  MyhistoManager.WriteMyHisto(InvM_ll_mumumu_afterleptsel_mWT110, "all");
  MyhistoManager.WriteMyHisto(InvM_ll_mumue_afterleptsel_mWT110,  "all");
  MyhistoManager.WriteMyHisto(InvM_ll_eemu_afterleptsel_mWT110,   "all");
  MyhistoManager.WriteMyHisto(InvM_ll_eee_afterleptsel_mWT110,    "all");
  
  MyhistoManager.WriteMyHisto(InvM_ll_mumumu_afterleptsel_lowbin, "all");
  MyhistoManager.WriteMyHisto(InvM_ll_mumue_afterleptsel_lowbin,  "all");
  MyhistoManager.WriteMyHisto(InvM_ll_eemu_afterleptsel_lowbin,   "all");
  MyhistoManager.WriteMyHisto(InvM_ll_eee_afterleptsel_lowbin,    "all");
    
  
  MyhistoManager.WriteMyHisto(InvM_ll_mumumu_afterjetsel, "all");
  MyhistoManager.WriteMyHisto(InvM_ll_mumue_afterjetsel,  "all");
  MyhistoManager.WriteMyHisto(InvM_ll_eemu_afterjetsel,   "all");
  MyhistoManager.WriteMyHisto(InvM_ll_eee_afterjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(InvM_ll_mumumu_afterbjetsel, "all");
  MyhistoManager.WriteMyHisto(InvM_ll_mumue_afterbjetsel,  "all");
  MyhistoManager.WriteMyHisto(InvM_ll_eemu_afterbjetsel,   "all");
  MyhistoManager.WriteMyHisto(InvM_ll_eee_afterbjetsel,    "all");
    
    
  MyhistoManager.WriteMyHisto(LeptPt_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto(LeptPt_mumue_afterleptsel,  "all");
  MyhistoManager.WriteMyHisto(LeptPt_eemu_afterleptsel,   "all");
  MyhistoManager.WriteMyHisto(LeptPt_eee_afterleptsel,    "all");
  
  MyhistoManager.WriteMyHisto(LeptPt_mumumu_afterjetsel, "all");
  MyhistoManager.WriteMyHisto(LeptPt_mumue_afterjetsel,  "all");
  MyhistoManager.WriteMyHisto(LeptPt_eemu_afterjetsel,   "all");
  MyhistoManager.WriteMyHisto(LeptPt_eee_afterjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(LeptPt_mumumu_afterbjetsel, "all");
  MyhistoManager.WriteMyHisto(LeptPt_mumue_afterbjetsel,  "all");
  MyhistoManager.WriteMyHisto(LeptPt_eemu_afterbjetsel,   "all");
  MyhistoManager.WriteMyHisto(LeptPt_eee_afterbjetsel,    "all");
  
    
    
  MyhistoManager.WriteMyHisto(LeptZPt_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto(LeptZPt_mumue_afterleptsel,  "all");
  MyhistoManager.WriteMyHisto(LeptZPt_eemu_afterleptsel,   "all");
  MyhistoManager.WriteMyHisto(LeptZPt_eee_afterleptsel,    "all");
  
  MyhistoManager.WriteMyHisto(LeptZPt_mumumu_afterjetsel, "all");
  MyhistoManager.WriteMyHisto(LeptZPt_mumue_afterjetsel,  "all");
  MyhistoManager.WriteMyHisto(LeptZPt_eemu_afterjetsel,   "all");
  MyhistoManager.WriteMyHisto(LeptZPt_eee_afterjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(LeptZPt_mumumu_afterbjetsel, "all");
  MyhistoManager.WriteMyHisto(LeptZPt_mumue_afterbjetsel,  "all");
  MyhistoManager.WriteMyHisto(LeptZPt_eemu_afterbjetsel,   "all");
  MyhistoManager.WriteMyHisto(LeptZPt_eee_afterbjetsel,    "all");
  
    
    
  MyhistoManager.WriteMyHisto(LeptWPt_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto(LeptWPt_mumue_afterleptsel,  "all");
  MyhistoManager.WriteMyHisto(LeptWPt_eemu_afterleptsel,   "all");
  MyhistoManager.WriteMyHisto(LeptWPt_eee_afterleptsel,    "all");
  
  MyhistoManager.WriteMyHisto(LeptWPt_mumumu_afterjetsel, "all");
  MyhistoManager.WriteMyHisto(LeptWPt_mumue_afterjetsel,  "all");
  MyhistoManager.WriteMyHisto(LeptWPt_eemu_afterjetsel,   "all");
  MyhistoManager.WriteMyHisto(LeptWPt_eee_afterjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(LeptWPt_mumumu_afterbjetsel, "all");
  MyhistoManager.WriteMyHisto(LeptWPt_mumue_afterbjetsel,  "all");
  MyhistoManager.WriteMyHisto(LeptWPt_eemu_afterbjetsel,   "all");
  MyhistoManager.WriteMyHisto(LeptWPt_eee_afterbjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(LeptWPt_mumumu_afterbjetveto, "all");
  MyhistoManager.WriteMyHisto(LeptWPt_mumue_afterbjetveto,  "all");
  MyhistoManager.WriteMyHisto(LeptWPt_eemu_afterbjetveto,   "all");
  MyhistoManager.WriteMyHisto(LeptWPt_eee_afterbjetveto,    "all");
  
    
    
  MyhistoManager.WriteMyHisto(LeptWPt_mumumu_afterleptsel_mWT110, "all");
  MyhistoManager.WriteMyHisto(LeptWPt_mumue_afterleptsel_mWT110,  "all");
  MyhistoManager.WriteMyHisto(LeptWPt_eemu_afterleptsel_mWT110,   "all");
  MyhistoManager.WriteMyHisto(LeptWPt_eee_afterleptsel_mWT110,    "all");
  
    
  MyhistoManager.WriteMyHisto(JetPt_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto(JetPt_mumue_afterleptsel,  "all");
  MyhistoManager.WriteMyHisto(JetPt_eemu_afterleptsel,   "all");
  MyhistoManager.WriteMyHisto(JetPt_eee_afterleptsel,    "all");
  
  MyhistoManager.WriteMyHisto(JetPt_mumumu_afterjetsel, "all");
  MyhistoManager.WriteMyHisto(JetPt_mumue_afterjetsel,  "all");
  MyhistoManager.WriteMyHisto(JetPt_eemu_afterjetsel,   "all");
  MyhistoManager.WriteMyHisto(JetPt_eee_afterjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(JetPt_mumumu_afterbjetsel, "all");
  MyhistoManager.WriteMyHisto(JetPt_mumue_afterbjetsel,  "all");
  MyhistoManager.WriteMyHisto(JetPt_eemu_afterbjetsel,   "all");
  MyhistoManager.WriteMyHisto(JetPt_eee_afterbjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(JetPt_mumumu_afterbjetveto, "all");
  MyhistoManager.WriteMyHisto(JetPt_mumue_afterbjetveto,  "all");
  MyhistoManager.WriteMyHisto(JetPt_eemu_afterbjetveto,   "all");
  MyhistoManager.WriteMyHisto(JetPt_eee_afterbjetveto,    "all");
  
    
    
  MyhistoManager.WriteMyHisto(JetEta_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto(JetEta_mumue_afterleptsel,  "all");
  MyhistoManager.WriteMyHisto(JetEta_eemu_afterleptsel,   "all");
  MyhistoManager.WriteMyHisto(JetEta_eee_afterleptsel,    "all");
  
  MyhistoManager.WriteMyHisto(JetEta_mumumu_afterjetsel, "all");
  MyhistoManager.WriteMyHisto(JetEta_mumue_afterjetsel,  "all");
  MyhistoManager.WriteMyHisto(JetEta_eemu_afterjetsel,   "all");
  MyhistoManager.WriteMyHisto(JetEta_eee_afterjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(JetEta_mumumu_afterbjetsel, "all");
  MyhistoManager.WriteMyHisto(JetEta_mumue_afterbjetsel,  "all");
  MyhistoManager.WriteMyHisto(JetEta_eemu_afterbjetsel,   "all");
  MyhistoManager.WriteMyHisto(JetEta_eee_afterbjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(JetEta_mumumu_afterbjetveto, "all");
  MyhistoManager.WriteMyHisto(JetEta_mumue_afterbjetveto,  "all");
  MyhistoManager.WriteMyHisto(JetEta_eemu_afterbjetveto,   "all");
  MyhistoManager.WriteMyHisto(JetEta_eee_afterbjetveto,    "all");
  
  
  MyhistoManager.WriteMyHisto(HT_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto(HT_mumue_afterleptsel,  "all");
  MyhistoManager.WriteMyHisto(HT_eemu_afterleptsel,   "all");
  MyhistoManager.WriteMyHisto(HT_eee_afterleptsel,    "all");
  
  
  MyhistoManager.WriteMyHisto(HT_mumumu_afterjetsel, "all");
  MyhistoManager.WriteMyHisto(HT_mumue_afterjetsel,  "all");
  MyhistoManager.WriteMyHisto(HT_eemu_afterjetsel,   "all");
  MyhistoManager.WriteMyHisto(HT_eee_afterjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(HT_mumumu_afterbjetsel, "all");
  MyhistoManager.WriteMyHisto(HT_mumue_afterbjetsel,  "all");
  MyhistoManager.WriteMyHisto(HT_eemu_afterbjetsel,   "all");
  MyhistoManager.WriteMyHisto(HT_eee_afterbjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(HT_mumumu_afterbjetveto, "all");
  MyhistoManager.WriteMyHisto(HT_mumue_afterbjetveto,  "all");
  MyhistoManager.WriteMyHisto(HT_eemu_afterbjetveto,   "all");
  MyhistoManager.WriteMyHisto(HT_eee_afterbjetveto,    "all");
  
  
  
  
  
  MyhistoManager.WriteMyHisto(MET_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto(MET_mumue_afterleptsel,  "all");
  MyhistoManager.WriteMyHisto(MET_eemu_afterleptsel,   "all");
  MyhistoManager.WriteMyHisto(MET_eee_afterleptsel,    "all");
  
  MyhistoManager.WriteMyHisto(MET_mumumu_afterleptsel_mWT110, "all");
  MyhistoManager.WriteMyHisto(MET_mumue_afterleptsel_mWT110,  "all");
  MyhistoManager.WriteMyHisto(MET_eemu_afterleptsel_mWT110,   "all");
  MyhistoManager.WriteMyHisto(MET_eee_afterleptsel_mWT110,    "all");
  
  
  MyhistoManager.WriteMyHisto(MET_mumumu_afterjetsel, "all");
  MyhistoManager.WriteMyHisto(MET_mumue_afterjetsel,  "all");
  MyhistoManager.WriteMyHisto(MET_eemu_afterjetsel,   "all");
  MyhistoManager.WriteMyHisto(MET_eee_afterjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(MET_mumumu_afterbjetsel, "all");
  MyhistoManager.WriteMyHisto(MET_mumue_afterbjetsel,  "all");
  MyhistoManager.WriteMyHisto(MET_eemu_afterbjetsel,   "all");
  MyhistoManager.WriteMyHisto(MET_eee_afterbjetsel,    "all");
  
  
  
  MyhistoManager.WriteMyHisto(InvM_ll_mumumu_afterleptsel_highSumPt, "all");
  MyhistoManager.WriteMyHisto(InvM_ll_mumue_afterleptsel_highSumPt, "all");
  MyhistoManager.WriteMyHisto(InvM_ll_eemu_afterleptsel_highSumPt, "all");
  MyhistoManager.WriteMyHisto(InvM_ll_eee_afterleptsel_highSumPt, "all");
  
  
  
  
  MyhistoManager.WriteMyHisto(Asym_mumumu_afterbjetsel, "all");
  MyhistoManager.WriteMyHisto(Asym_mumue_afterbjetsel,  "all");
  MyhistoManager.WriteMyHisto(Asym_eemu_afterbjetsel,   "all");
  MyhistoManager.WriteMyHisto(Asym_eee_afterbjetsel,    "all");
  
  
  
  MyhistoManager.WriteMyHisto(RecoPtZ_mumumu_afterbjetsel, "all");
  MyhistoManager.WriteMyHisto(RecoPtZ_mumue_afterbjetsel,  "all");
  MyhistoManager.WriteMyHisto(RecoPtZ_eemu_afterbjetsel,   "all");
  MyhistoManager.WriteMyHisto(RecoPtZ_eee_afterbjetsel,    "all");
  
  
  
  MyhistoManager.WriteMyHisto(RecoPtZ_mumumu_afterbjetveto, "all");
  MyhistoManager.WriteMyHisto(RecoPtZ_mumue_afterbjetveto,  "all");
  MyhistoManager.WriteMyHisto(RecoPtZ_eemu_afterbjetveto,   "all");
  MyhistoManager.WriteMyHisto(RecoPtZ_eee_afterbjetveto,    "all");
  
  MyhistoManager.WriteMyHisto(RecoPtZ_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto(RecoPtZ_mumue_afterleptsel,  "all");
  MyhistoManager.WriteMyHisto(RecoPtZ_eemu_afterleptsel,   "all");
  MyhistoManager.WriteMyHisto(RecoPtZ_eee_afterleptsel,    "all");
  
  
  MyhistoManager.WriteMyHisto(RecoPtZ_mumumu_afterleptsel_nojet, "all");
  MyhistoManager.WriteMyHisto(RecoPtZ_mumue_afterleptsel_nojet,  "all");
  MyhistoManager.WriteMyHisto(RecoPtZ_eemu_afterleptsel_nojet,   "all");
  MyhistoManager.WriteMyHisto(RecoPtZ_eee_afterleptsel_nojet,    "all");
  
  MyhistoManager.WriteMyHisto(RecoTopMass_mumumu_afterbjetsel  , "all");
  MyhistoManager.WriteMyHisto(RecoTopMass_mumue_afterbjetsel,   "all");
  MyhistoManager.WriteMyHisto(RecoTopMass_eemu_afterbjetsel,    "all");
  MyhistoManager.WriteMyHisto(RecoTopMass_eee_afterbjetsel,     "all");
  
  MyhistoManager.WriteMyHisto(RecoTopMass_mumumu_afterbjetveto  , "all");
  MyhistoManager.WriteMyHisto(RecoTopMass_mumue_afterbjetveto,   "all");
  MyhistoManager.WriteMyHisto(RecoTopMass_eemu_afterbjetveto,    "all");
  MyhistoManager.WriteMyHisto(RecoTopMass_eee_afterbjetveto,     "all");
  
  
  
  MyhistoManager.WriteMyHisto(deltaPhilb_mumumu_afterbjetsel, "all");
  MyhistoManager.WriteMyHisto(deltaPhilb_mumue_afterbjetsel,  "all");
  MyhistoManager.WriteMyHisto(deltaPhilb_eemu_afterbjetsel,   "all");
  MyhistoManager.WriteMyHisto(deltaPhilb_eee_afterbjetsel,    "all");
  
  MyhistoManager.WriteMyHisto(deltaPhilj_mumumu_afterbjetveto, "all");
  MyhistoManager.WriteMyHisto(deltaPhilj_mumue_afterbjetveto,  "all");
  MyhistoManager.WriteMyHisto(deltaPhilj_eemu_afterbjetveto,   "all");
  MyhistoManager.WriteMyHisto(deltaPhilj_eee_afterbjetveto,    "all");
  
  MyhistoManager.WriteMyHisto(deltaR_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto(deltaR_mumue_afterleptsel,  "all");
  MyhistoManager.WriteMyHisto(deltaR_eemu_afterleptsel,   "all");
  MyhistoManager.WriteMyHisto(deltaR_eee_afterleptsel,    "all");
  
    
  MyhistoManager.WriteMyHisto(WmissAssing_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto(WmissAssing_mumue_afterleptsel,  "all");
  MyhistoManager.WriteMyHisto(WmissAssing_eemu_afterleptsel,   "all");
  MyhistoManager.WriteMyHisto(WmissAssing_eee_afterleptsel,    "all");
  
  
  MyhistoManager.WriteMyHisto(DijetInvM_mumumu_afterleptsel_inZpeak,    "all");
  MyhistoManager.WriteMyHisto(DijetInvM_mumue_afterleptsel_inZpeak,     "all");
  MyhistoManager.WriteMyHisto(DijetInvM_eemu_afterleptsel_inZpeak,      "all");
  MyhistoManager.WriteMyHisto(DijetInvM_eee_afterleptsel_inZpeak,       "all");
  
  
  
  
  
  
  MyhistoManager.WriteMyHisto(mWT_mumumu_afterbjetveto, "all");
  MyhistoManager.WriteMyHisto(mWT_mumue_afterbjetveto,  "all");
  MyhistoManager.WriteMyHisto(mWT_eemu_afterbjetveto,   "all");
  MyhistoManager.WriteMyHisto(mWT_eee_afterbjetveto,    "all");

  MyhistoManager.WriteMyHisto(mWT_mumumu_afterbjetsel, "all");
  MyhistoManager.WriteMyHisto(mWT_mumue_afterbjetsel,  "all");
  MyhistoManager.WriteMyHisto(mWT_eemu_afterbjetsel,   "all");
  MyhistoManager.WriteMyHisto(mWT_eee_afterbjetsel,    "all");

  MyhistoManager.WriteMyHisto(mWT_mumumu_afterjetsel, "all");
  MyhistoManager.WriteMyHisto(mWT_mumue_afterjetsel,  "all");
  MyhistoManager.WriteMyHisto(mWT_eemu_afterjetsel,   "all");
  MyhistoManager.WriteMyHisto(mWT_eee_afterjetsel,    "all");

  MyhistoManager.WriteMyHisto(mWT_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto(mWT_mumue_afterleptsel,  "all");
  MyhistoManager.WriteMyHisto(mWT_eemu_afterleptsel,   "all");
  MyhistoManager.WriteMyHisto(mWT_eee_afterleptsel,    "all");

  MyhistoManager.WriteMyHisto(Charge_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto(Charge_mumue_afterleptsel,  "all");
  MyhistoManager.WriteMyHisto(Charge_eemu_afterleptsel,   "all");
  MyhistoManager.WriteMyHisto(Charge_eee_afterleptsel,    "all");
  
  MyhistoManager.WriteMyHisto(Charge_mumumu_afterleptsel_mWT110, "all");
  MyhistoManager.WriteMyHisto(Charge_mumue_afterleptsel_mWT110,  "all");
  MyhistoManager.WriteMyHisto(Charge_eemu_afterleptsel_mWT110,   "all");
  MyhistoManager.WriteMyHisto(Charge_eee_afterleptsel_mWT110,    "all");
  
  
  
  
  
  MyhistoManager.WriteMyHisto(deltaRLeptJet_mumumu_afterleptsel_mWT110, "all");
  MyhistoManager.WriteMyHisto(deltaRLeptJet_mumue_afterleptsel_mWT110,  "all");
  MyhistoManager.WriteMyHisto(deltaRLeptJet_eemu_afterleptsel_mWT110,   "all");
  MyhistoManager.WriteMyHisto(deltaRLeptJet_eee_afterleptsel_mWT110,    "all");
  
  MyhistoManager.WriteMyHisto(deltaRLeptMet_mumumu_afterleptsel_mWT110, "all");
  MyhistoManager.WriteMyHisto(deltaRLeptMet_mumue_afterleptsel_mWT110,  "all");
  MyhistoManager.WriteMyHisto(deltaRLeptMet_eemu_afterleptsel_mWT110,   "all");
  MyhistoManager.WriteMyHisto(deltaRLeptMet_eee_afterleptsel_mWT110,    "all");
  
  
  MyhistoManager.WriteMyHisto(NJet_mumumu_afterleptsel_mWT110, "all");
  MyhistoManager.WriteMyHisto(NJet_mumue_afterleptsel_mWT110,  "all");
  MyhistoManager.WriteMyHisto(NJet_eemu_afterleptsel_mWT110,   "all");
  MyhistoManager.WriteMyHisto(NJet_eee_afterleptsel_mWT110,    "all");
  
  
  
  
  MyhistoManager.WriteMyHisto(NBJet_mumumu_afterleptsel_mWT110, "all");
  MyhistoManager.WriteMyHisto(NBJet_mumue_afterleptsel_mWT110,  "all");
  MyhistoManager.WriteMyHisto(NBJet_eemu_afterleptsel_mWT110,   "all");
  MyhistoManager.WriteMyHisto(NBJet_eee_afterleptsel_mWT110,    "all");
  
  
  MyhistoManager.WriteMyHisto(NBJet_mumumu_afterjetsel_bjets, "all");
  MyhistoManager.WriteMyHisto(NBJet_mumue_afterjetsel_bjets,  "all"); 
  MyhistoManager.WriteMyHisto(NBJet_eemu_afterjetsel_bjets,   "all");  
  MyhistoManager.WriteMyHisto(NBJet_eee_afterjetsel_bjets,    "all");   


  MyhistoManager.WriteMyHisto(NBJet_mumumu_afterjetsel_cjets, "all");
  MyhistoManager.WriteMyHisto(NBJet_mumue_afterjetsel_cjets,  "all"); 
  MyhistoManager.WriteMyHisto(NBJet_eemu_afterjetsel_cjets,   "all");  
  MyhistoManager.WriteMyHisto(NBJet_eee_afterjetsel_cjets,    "all");   


  MyhistoManager.WriteMyHisto(NBJet_mumumu_afterjetsel_ljets, "all");
  MyhistoManager.WriteMyHisto(NBJet_mumue_afterjetsel_ljets,  "all"); 
  MyhistoManager.WriteMyHisto(NBJet_eemu_afterjetsel_ljets,   "all");  
  MyhistoManager.WriteMyHisto(NBJet_eee_afterjetsel_ljets ,   "all");  

  
  
  MyhistoManager.WriteMyHisto(BJetDiscri_mumumu_afterjetsel_bjets,    "all");
  MyhistoManager.WriteMyHisto(BJetDiscri_mumue_afterjetsel_bjets,     "all");
  MyhistoManager.WriteMyHisto(BJetDiscri_eemu_afterjetsel_bjets,      "all");
  MyhistoManager.WriteMyHisto(BJetDiscri_eee_afterjetsel_bjets,       "all");
  
  
  
  MyhistoManager.WriteMyHisto(BJetDiscri_mumumu_afterjetsel_cjets,    "all");
  MyhistoManager.WriteMyHisto(BJetDiscri_mumue_afterjetsel_cjets,     "all");
  MyhistoManager.WriteMyHisto(BJetDiscri_eemu_afterjetsel_cjets,      "all");
  MyhistoManager.WriteMyHisto(BJetDiscri_eee_afterjetsel_cjets,       "all");
  
  
  MyhistoManager.WriteMyHisto(BJetDiscri_mumumu_afterjetsel_ljets,    "all");
  MyhistoManager.WriteMyHisto(BJetDiscri_mumue_afterjetsel_ljets,     "all");
  MyhistoManager.WriteMyHisto(BJetDiscri_eemu_afterjetsel_ljets,      "all");
  MyhistoManager.WriteMyHisto(BJetDiscri_eee_afterjetsel_ljets,       "all");
  
  
  
  
  MyhistoManager.WriteMyHisto2D(HT_vs_MET_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_MET_mumue_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_MET_eemu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_MET_eee_afterleptsel, "all");
 
  MyhistoManager.WriteMyHisto2D(HT_vs_NJet_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_NJet_mumue_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_NJet_eemu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_NJet_eee_afterleptsel, "all");
 
  MyhistoManager.WriteMyHisto2D(HT_vs_NBJet_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_NBJet_mumue_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_NBJet_eemu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_NBJet_eee_afterleptsel, "all");
 
  MyhistoManager.WriteMyHisto2D(HT_vs_LeptPt_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_LeptPt_mumue_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_LeptPt_eemu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_LeptPt_eee_afterleptsel, "all");
 
  MyhistoManager.WriteMyHisto2D(HT_vs_JetPt_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_JetPt_mumue_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_JetPt_eemu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_JetPt_eee_afterleptsel, "all");
  
  
  
  MyhistoManager.WriteMyHisto2D(HT_vs_Mll_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_Mll_mumue_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_Mll_eemu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(HT_vs_Mll_eee_afterleptsel, "all");
  
  MyhistoManager.WriteMyHisto2D(InvM_ll_vs_mWT_mumumu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(InvM_ll_vs_mWT_mumue_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(InvM_ll_vs_mWT_eemu_afterleptsel, "all");
  MyhistoManager.WriteMyHisto2D(InvM_ll_vs_mWT_eee_afterleptsel, "all");
  
  
  
  TheTree->Write();
  
  
  
   //The following line is mandatory to copy everythin in a common RootFile
    fOutput->Add(fProofFile);
    
  CutFlow_mumumu.clear();
  CutFlow_mumue.clear();
  CutFlow_eemu.clear();
  CutFlow_eee.clear();
  ErrCutFlow_mumumu.clear();
  ErrCutFlow_mumue.clear();
  ErrCutFlow_eemu.clear();
  ErrCutFlow_eee.clear();
  
  
  
  Ntrilept_mumumu.clear();
  Ntrileptnoniso_mumumu.clear();
 
  
  
  DijetInvM_mumumu_afterleptsel_inZpeak.clear();
  DijetInvM_mumue_afterleptsel_inZpeak.clear();
  DijetInvM_eemu_afterleptsel_inZpeak.clear();
  DijetInvM_eee_afterleptsel_inZpeak.clear();
  
  
  
  Mt_mumumu_afterbjetsel.clear();
  Mt_mumue_afterbjetsel.clear();
  Mt_eemu_afterbjetsel.clear();
  Mt_eee_afterbjetsel.clear();
   
  Mt_mumumu_afterbjetveto.clear();
  Mt_mumue_afterbjetveto.clear();
  Mt_eemu_afterbjetveto.clear();
  Mt_eee_afterbjetveto.clear();
   
  
  NJet_mumumu_afterZsel.clear();
  NJet_mumue_afterZsel.clear();
  NJet_eemu_afterZsel.clear();
  NJet_eee_afterZsel.clear();
  
  
  NJet_mumumu_afterbsel.clear();
  NJet_mumue_afterbsel.clear();
  NJet_eemu_afterbsel.clear();
  NJet_eee_afterbsel.clear();
  
  NJet_mumumu_afterleptsel_mWT110.clear();
  NJet_mumue_afterleptsel_mWT110.clear();
  NJet_eemu_afterleptsel_mWT110.clear();
  NJet_eee_afterleptsel_mWT110.clear();
  
  
  NLept_mumumu_afterbsel.clear();
  NLept_mumue_afterbsel.clear();
  NLept_eemu_afterbsel.clear();
  NLept_eee_afterbsel.clear();
  
  
  NBJet_mumumu_afterZsel.clear();
  NBJet_mumue_afterZsel.clear();
  NBJet_eemu_afterZsel.clear();
  NBJet_eee_afterZsel.clear();
  
  
  NBJet_mumumu_afterjetsel.clear();
  NBJet_mumue_afterjetsel.clear();
  NBJet_eemu_afterjetsel.clear();
  NBJet_eee_afterjetsel.clear();
  
  
  
  NBJet_mumumu_afterleptsel_mWT110.clear();
  NBJet_mumue_afterleptsel_mWT110.clear();
  NBJet_eemu_afterleptsel_mWT110.clear();
  NBJet_eee_afterleptsel_mWT110.clear();
  
  //to be filled
  
  Nvtx_mumumu_afterleptsel.clear();
  Nvtx_mumue_afterleptsel.clear();
  Nvtx_eemu_afterleptsel.clear();
  Nvtx_eee_afterleptsel.clear();
  
  InvM_ll_mumumu_afterleptsel.clear();
  InvM_ll_mumue_afterleptsel.clear();
  InvM_ll_eemu_afterleptsel.clear();
  InvM_ll_eee_afterleptsel.clear();
  
  InvM_ll_mumumu_afterleptsel_mWT110.clear();
  InvM_ll_mumue_afterleptsel_mWT110.clear();
  InvM_ll_eemu_afterleptsel_mWT110.clear();
  InvM_ll_eee_afterleptsel_mWT110.clear();
  
  InvM_ll_mumumu_afterleptsel_lowbin.clear();
  InvM_ll_mumue_afterleptsel_lowbin.clear();
  InvM_ll_eemu_afterleptsel_lowbin.clear();
  InvM_ll_eee_afterleptsel_lowbin.clear();
  
  InvM_ll_mumumu_afterleptsel_highSumPt.clear();
  InvM_ll_mumue_afterleptsel_highSumPt.clear();
  InvM_ll_eemu_afterleptsel_highSumPt.clear();
  InvM_ll_eee_afterleptsel_highSumPt.clear();
  
  InvM_ll_mumumu_afterjetsel.clear();
  InvM_ll_mumue_afterjetsel.clear();
  InvM_ll_eemu_afterjetsel.clear();
  InvM_ll_eee_afterjetsel.clear();
  
  InvM_ll_mumumu_afterbjetsel.clear();
  InvM_ll_mumue_afterbjetsel.clear();
  InvM_ll_eemu_afterbjetsel.clear();
  InvM_ll_eee_afterbjetsel.clear();
  
  LeptPt_mumumu_afterleptsel.clear();
  LeptPt_mumue_afterleptsel.clear();
  LeptPt_eemu_afterleptsel.clear();
  LeptPt_eee_afterleptsel.clear();
  
  LeptPt_mumumu_afterjetsel.clear();
  LeptPt_mumue_afterjetsel.clear();
  LeptPt_eemu_afterjetsel.clear();
  LeptPt_eee_afterjetsel.clear();
  
  LeptPt_mumumu_afterbjetsel.clear();
  LeptPt_mumue_afterbjetsel.clear();
  LeptPt_eemu_afterbjetsel.clear();
  LeptPt_eee_afterbjetsel.clear();
  
  
  
  LeptZPt_mumumu_afterleptsel.clear();
  LeptZPt_mumue_afterleptsel.clear();
  LeptZPt_eemu_afterleptsel.clear();
  LeptZPt_eee_afterleptsel.clear();
  
  LeptZPt_mumumu_afterjetsel.clear();
  LeptZPt_mumue_afterjetsel.clear();
  LeptZPt_eemu_afterjetsel.clear();
  LeptZPt_eee_afterjetsel.clear();
  
  LeptZPt_mumumu_afterbjetsel.clear();
  LeptZPt_mumue_afterbjetsel.clear();
  LeptZPt_eemu_afterbjetsel.clear();
  LeptZPt_eee_afterbjetsel.clear();
  
  
  
  LeptWPt_mumumu_afterleptsel.clear();
  LeptWPt_mumue_afterleptsel.clear();
  LeptWPt_eemu_afterleptsel.clear();
  LeptWPt_eee_afterleptsel.clear();
  
  LeptWPt_mumumu_afterjetsel.clear();
  LeptWPt_mumue_afterjetsel.clear();
  LeptWPt_eemu_afterjetsel.clear();
  LeptWPt_eee_afterjetsel.clear();
  
  LeptWPt_mumumu_afterbjetsel.clear();
  LeptWPt_mumue_afterbjetsel.clear();
  LeptWPt_eemu_afterbjetsel.clear();
  LeptWPt_eee_afterbjetsel.clear();
  
  LeptWPt_mumumu_afterbjetveto.clear();
  LeptWPt_mumue_afterbjetveto.clear();
  LeptWPt_eemu_afterbjetveto.clear();
  LeptWPt_eee_afterbjetveto.clear();
  
  
  LeptWPt_mumumu_afterleptsel_mWT110.clear();
  LeptWPt_mumue_afterleptsel_mWT110.clear();
  LeptWPt_eemu_afterleptsel_mWT110.clear();
  LeptWPt_eee_afterleptsel_mWT110.clear();
  
  
  JetPt_mumumu_afterleptsel.clear();
  JetPt_mumue_afterleptsel.clear();
  JetPt_eemu_afterleptsel.clear();
  JetPt_eee_afterleptsel.clear();
  
  JetPt_mumumu_afterjetsel.clear();
  JetPt_mumue_afterjetsel.clear();
  JetPt_eemu_afterjetsel.clear();
  JetPt_eee_afterjetsel.clear();
  
  JetPt_mumumu_afterbjetsel.clear();
  JetPt_mumue_afterbjetsel.clear();
  JetPt_eemu_afterbjetsel.clear();
  JetPt_eee_afterbjetsel.clear();
  
  JetPt_mumumu_afterbjetveto.clear();
  JetPt_mumue_afterbjetveto.clear();
  JetPt_eemu_afterbjetveto.clear();
  JetPt_eee_afterbjetveto.clear();
  
  JetEta_mumumu_afterleptsel.clear();
  JetEta_mumue_afterleptsel.clear();
  JetEta_eemu_afterleptsel.clear();
  JetEta_eee_afterleptsel.clear();
  
  JetEta_mumumu_afterjetsel.clear();
  JetEta_mumue_afterjetsel.clear();
  JetEta_eemu_afterjetsel.clear();
  JetEta_eee_afterjetsel.clear();
  
  JetEta_mumumu_afterbjetsel.clear();
  JetEta_mumue_afterbjetsel.clear();
  JetEta_eemu_afterbjetsel.clear();
  JetEta_eee_afterbjetsel.clear();
  
  JetEta_mumumu_afterbjetveto.clear();
  JetEta_mumue_afterbjetveto.clear();
  JetEta_eemu_afterbjetveto.clear();
  JetEta_eee_afterbjetveto.clear();
  
  HT_mumumu_afterleptsel.clear();
  HT_mumue_afterleptsel.clear();
  HT_eemu_afterleptsel.clear();
  HT_eee_afterleptsel.clear();
  
  
  HT_mumumu_afterjetsel.clear();
  HT_mumue_afterjetsel.clear();
  HT_eemu_afterjetsel.clear();
  HT_eee_afterjetsel.clear();
  
  HT_mumumu_afterbjetsel.clear();
  HT_mumue_afterbjetsel.clear();
  HT_eemu_afterbjetsel.clear();
  HT_eee_afterbjetsel.clear();
  
  HT_mumumu_afterbjetveto.clear();
  HT_mumue_afterbjetveto.clear();
  HT_eemu_afterbjetveto.clear();
  HT_eee_afterbjetveto.clear();
  
  
  
  
  
  MET_mumumu_afterleptsel.clear();
  MET_mumue_afterleptsel.clear();
  MET_eemu_afterleptsel.clear();
  MET_eee_afterleptsel.clear();
  
  
  
  MET_mumumu_afterleptsel_mWT110.clear();
  MET_mumue_afterleptsel_mWT110.clear();
  MET_eemu_afterleptsel_mWT110.clear();
  MET_eee_afterleptsel_mWT110.clear();
  
  
  MET_mumumu_afterjetsel.clear();
  MET_mumue_afterjetsel.clear();
  MET_eemu_afterjetsel.clear();
  MET_eee_afterjetsel.clear();
  
  MET_mumumu_afterbjetsel.clear();
  MET_mumue_afterbjetsel.clear();
  MET_eemu_afterbjetsel.clear();
  MET_eee_afterbjetsel.clear();
  
  Asym_mumumu_afterbjetsel.clear();
  Asym_mumue_afterbjetsel.clear();
  Asym_eemu_afterbjetsel.clear();
  Asym_eee_afterbjetsel.clear();
  
  
  
  mWT_mumumu_afterjetsel.clear();
  mWT_mumue_afterjetsel.clear();
  mWT_eemu_afterjetsel.clear();
  mWT_eee_afterjetsel.clear();
  
  
  
  RecoPtZ_mumumu_afterbjetsel.clear();
  RecoPtZ_mumue_afterbjetsel.clear();
  RecoPtZ_eemu_afterbjetsel.clear();
  RecoPtZ_eee_afterbjetsel.clear();
  
  RecoPtZ_mumumu_afterbjetveto.clear();
  RecoPtZ_mumue_afterbjetveto.clear();
  RecoPtZ_eemu_afterbjetveto.clear();
  RecoPtZ_eee_afterbjetveto.clear();
  
  
  RecoPtZ_mumumu_afterleptsel.clear();
  RecoPtZ_mumue_afterleptsel.clear();
  RecoPtZ_eemu_afterleptsel.clear();
  RecoPtZ_eee_afterleptsel.clear();
  
  
  RecoPtZ_mumumu_afterleptsel_nojet.clear();
  RecoPtZ_mumue_afterleptsel_nojet.clear();
  RecoPtZ_eemu_afterleptsel_nojet.clear();
  RecoPtZ_eee_afterleptsel_nojet.clear();
  
  
  RecoTopMass_mumumu_afterbjetsel.clear();
  RecoTopMass_mumue_afterbjetsel.clear();
  RecoTopMass_eemu_afterbjetsel.clear();
  RecoTopMass_eee_afterbjetsel.clear();
  
  RecoTopMass_mumumu_afterbjetveto.clear();
  RecoTopMass_mumue_afterbjetveto.clear();
  RecoTopMass_eemu_afterbjetveto.clear();
  RecoTopMass_eee_afterbjetveto.clear();
  
  
  deltaPhilb_mumumu_afterbjetsel.clear();
  deltaPhilb_mumue_afterbjetsel.clear();
  deltaPhilb_eemu_afterbjetsel.clear();
  deltaPhilb_eee_afterbjetsel.clear();
  
  deltaPhilj_mumumu_afterbjetveto.clear();
  deltaPhilj_mumue_afterbjetveto.clear();
  deltaPhilj_eemu_afterbjetveto.clear();
  deltaPhilj_eee_afterbjetveto.clear();
  
  
  
  
  deltaR_mumumu_afterleptsel.clear();
  deltaR_mumue_afterleptsel.clear();
  deltaR_eemu_afterleptsel.clear();
  deltaR_eee_afterleptsel.clear();
  
  
  
  deltaRLeptJet_mumumu_afterleptsel_mWT110.clear();
  deltaRLeptJet_mumue_afterleptsel_mWT110.clear();
  deltaRLeptJet_eemu_afterleptsel_mWT110.clear();
  deltaRLeptJet_eee_afterleptsel_mWT110.clear();
  
  deltaRLeptMet_mumumu_afterleptsel_mWT110.clear();
  deltaRLeptMet_mumue_afterleptsel_mWT110.clear();
  deltaRLeptMet_eemu_afterleptsel_mWT110.clear();
  deltaRLeptMet_eee_afterleptsel_mWT110.clear();
  
  
  
  WmissAssing_mumumu_afterleptsel.clear();
  WmissAssing_mumue_afterleptsel.clear();
  WmissAssing_eemu_afterleptsel.clear();
  WmissAssing_eee_afterleptsel.clear();
  
  
  mWT_mumumu_afterleptsel.clear();
  mWT_mumue_afterleptsel.clear();
  mWT_eemu_afterleptsel.clear();
  mWT_eee_afterleptsel.clear();
  
  
  mWT_mumumu_afterbjetsel.clear();
  mWT_mumue_afterbjetsel.clear();
  mWT_eemu_afterbjetsel.clear();
  mWT_eee_afterbjetsel.clear();
  
  mWT_mumumu_afterbjetveto.clear();
  mWT_mumue_afterbjetveto.clear();
  mWT_eemu_afterbjetveto.clear();
  mWT_eee_afterbjetveto.clear();
  
  
  Charge_mumumu_afterleptsel.clear();
  Charge_mumue_afterleptsel.clear();
  Charge_eemu_afterleptsel.clear();
  Charge_eee_afterleptsel.clear();
  
  Charge_mumumu_afterleptsel_mWT110.clear();
  Charge_mumue_afterleptsel_mWT110.clear();
  Charge_eemu_afterleptsel_mWT110.clear();
  Charge_eee_afterleptsel_mWT110.clear();
  
 
  Nvertex.clear();
 
  InvM_ll_vs_mWT_mumumu_afterleptsel.clear();
  InvM_ll_vs_mWT_mumue_afterleptsel.clear();
  InvM_ll_vs_mWT_eemu_afterleptsel.clear();
  InvM_ll_vs_mWT_eee_afterleptsel.clear();
  
  HT_vs_MET_mumumu_afterleptsel.clear();
  HT_vs_MET_mumue_afterleptsel.clear();
  HT_vs_MET_eemu_afterleptsel.clear();
  HT_vs_MET_eee_afterleptsel.clear();
  
  HT_vs_NJet_mumumu_afterleptsel.clear();
  HT_vs_NJet_mumue_afterleptsel.clear();
  HT_vs_NJet_eemu_afterleptsel.clear();
  HT_vs_NJet_eee_afterleptsel.clear();
  
  HT_vs_NBJet_mumumu_afterleptsel.clear();
  HT_vs_NBJet_mumue_afterleptsel.clear();
  HT_vs_NBJet_eemu_afterleptsel.clear();
  HT_vs_NBJet_eee_afterleptsel.clear();
  
  HT_vs_LeptPt_mumumu_afterleptsel.clear();
  HT_vs_LeptPt_mumue_afterleptsel.clear();
  HT_vs_LeptPt_eemu_afterleptsel.clear();
  HT_vs_LeptPt_eee_afterleptsel.clear();
  
  HT_vs_JetPt_mumumu_afterleptsel.clear();
  HT_vs_JetPt_mumue_afterleptsel.clear();
  HT_vs_JetPt_eemu_afterleptsel.clear();
  HT_vs_JetPt_eee_afterleptsel.clear();
  
  
  
  HT_vs_Mll_mumumu_afterleptsel.clear();
  HT_vs_Mll_mumue_afterleptsel.clear();
  HT_vs_Mll_eemu_afterleptsel.clear();
  HT_vs_Mll_eee_afterleptsel.clear();
  
  
  
  
  
  
    cout << "global lepton SF mumumu : " << nEvents_mumumu/sumSFlept_mumumu << endl;
    cout << "global lepton SF mumue  : " << nEvents_mumue/sumSFlept_mumue << endl;
    cout << "global lepton SF eemu   : " << nEvents_eemu/sumSFlept_eemu << endl;
    cout << "global lepton SF eee    : " << nEvents_eee/sumSFlept_eee << endl;
  
  
  
  
  
  
    //delete file1  ;
    //delete hPUMC ;  
    //delete file2 ; 
    delete TheTree;
    delete anaEL;
    delete LumiWeights;
  }
}

//_____________________________________________________________________________
void ProofSelectorMyCutFlow::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  //Possibility to retrieve information from the merged file and perform some calculation or plotting tasks
  delete event ;
  //
  // Create canvas
  //
  //TList* list = fOutput->GetOutputList() ;
  /*
    TIter next_object((TList*) fOutput);
    TObject* obj ;
    cout << "-- Retrieved objects:" << endl ;
    while ((obj = next_object())) { TString objname = obj->GetName() ; cout << " " << objname << endl ; }
    
    if ((fi = dynamic_cast<TFile *>(fOutput->FindObject("blabla.root")))) {
    cout<<"Warning"<<endl;
    fi->Write("toto.root");
    cout<<"Warning"<<endl;
    }
    // Final update
    c1->cd();
    c1->Update();
  */
}


std::vector<double> ProofSelectorMyCutFlow::GetNvertexWeight(TString datasetName){

    
    std::vector<double> reweightPrivateProd;
    
    // for kut
    if(datasetName == "FCNCkut"){
     double reweightPrivateProd_tmp[31] = {  0 , 0.907961 , 0.874912 , 0.852147 , 0.935751 , 0.99628 , 1.03062 , 1.07897 , 1.07771 , 1.04004 , 1.04899 , 1.00515 , 0.998906 , 1.07023 , 1.04179 , 1.0139 , 1.03238 , 0.970094 , 1.01184 , 0.982466 , 0.961056 , 0.929841 , 1.04474 , 1.07386 , 1.36141 , 1.66165 , 1.89831 , 1.2499 , 1.40118 , 1.40688 , 0  };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }
    
    // for kct
    else if(datasetName == "FCNCkct"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.712725 , 0.656884 , 0.753967 , 0.960977 , 1.08034 , 1.13014 , 1.13592 , 1.02528 , 0.950872 , 0.965542 , 0.982825 , 1.03358 , 1.13427 , 1.18196 , 1.19021 , 1.19707 , 1.15671 , 1.20606 , 1.1747 , 1.19182 , 1.17983 , 1.35037 , 1.48941 , 1.68852 , 1.79627 , 1.63421 , 4.39171 , 3.72483 , 3.31721 , 0  };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

 
    // for zut
    else if(datasetName == "FCNCzut"){
     double reweightPrivateProd_tmp[31]= {  14.2631 , 0.697375 , 0.809106 , 0.869676 , 0.899343 , 0.956851 , 1.00431 , 1.05732 , 1.06992 , 1.04496 , 1.04483 , 1.03554 , 1.11835 , 1.09721 , 1.12549 , 1.18966 , 1.13856 , 1.11547 , 1.10133 , 1.09732 , 1.01627 , 0.990315 , 0.956436 , 1.01717 , 1.04933 , 1.49887 , 1.9157 , 1.67661 , 1.39616 , 1.51968 , 0  };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

 
    // for zct
    else if(datasetName == "FCNCzct"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.692144 , 0.817213 , 0.921795 , 0.909188 , 0.920288 , 0.945239 , 1.04917 , 1.10059 , 1.08028 , 1.06175 , 0.964076 , 0.960537 , 0.996581 , 1.00754 , 1.07115 , 1.18039 , 1.32376 , 1.31678 , 1.43264 , 1.29996 , 1.44647 , 1.68392 , 2.07507 , 2.74374 , 4.04146 , 4.44285 , 4.11708 , 12.5709 , 18.6586 , 0  };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for kut
    else if(datasetName == "FCNCkut_scaleup"){
     double reweightPrivateProd_tmp[31]= {  15.1379 , 0.551451 , 0.627648 , 0.837559 , 1.04268 , 1.11849 , 1.06767 , 1.03912 , 0.997338 , 0.930993 , 0.980354 , 0.991143 , 1.08498 , 1.16188 , 1.18731 , 1.19847 , 1.23108 , 1.3609 , 1.19203 , 1.29347 , 1.29243 , 1.20005 , 1.43753 , 1.28594 , 1.39211 , 1.59079 , 1.92023 , 1.60149 , 1.48179 , 1.05188 , 0 };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

 
 
 
 
 
    // for zut
    else if(datasetName == "FCNCzut_scaleup"){
     double reweightPrivateProd_tmp[31]= {  40.0456 , 0.795185 , 0.952177 , 0.97578 , 0.965051 , 0.945924 , 0.947971 , 0.957591 , 0.906844 , 0.890566 , 0.947077 , 1.00407 , 1.08736 , 1.12073 , 1.19486 , 1.24646 , 1.29959 , 1.34313 , 1.33116 , 1.26193 , 1.23704 , 1.26307 , 1.27553 , 1.44578 , 1.5832 , 1.63655 , 1.75838 , 2.54194 , 1.65842 , 2.13335 , 0  };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

 
    // for zct
    else if(datasetName == "FCNCzct_scaleup"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.938284 , 0.923777 , 0.991521 , 1.0126 , 0.916387 , 0.833552 , 0.854684 , 0.846158 , 0.826356 , 0.870978 , 0.990614 , 1.07197 , 1.15189 , 1.26719 , 1.569 , 1.68764 , 1.8092 , 2.11854 , 2.23978 , 2.3613 , 3.14149 , 2.98911 , 4.45639 , 4.47219 , 6.46073 , 13.035 , 9.0594 , 7.68374 , 22.8096 , 0  };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for kut
    else if(datasetName == "FCNCkut_scaledown"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.745082 , 0.809829 , 0.898396 , 0.969705 , 0.984589 , 1.00252 , 1.01246 , 1.03578 , 1.0233 , 1.02326 , 1.02118 , 0.998246 , 1.02082 , 1.05408 , 1.06153 , 1.11223 , 1.14877 , 1.09703 , 1.12478 , 1.2645 , 1.19921 , 1.38483 , 1.55163 , 2.08558 , 2.16232 , 2.47029 , 2.37721 , 2.09688 , 3.89044 , 0  };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }
    // for kct
    else if(datasetName == "FCNCkct_scaledown"){
     double reweightPrivateProd_tmp[31]= { 0 , 0.703597 , 0.727591 , 1.28895 , 1.61283 , 1.42395 , 1.18867 , 0.922161 , 0.834933 , 0.76834 , 0.754498 , 0.779636 , 0.804281 , 0.871038 , 1.03117 , 1.07804 , 1.14148 , 1.20386 , 1.18028 , 1.35359 , 1.46127 , 1.69478 , 1.86092 , 2.46661 , 2.07899 , 3.19344 , 4.24813 , 5.16684 , 9.34882 , 6.93812 , 0  } ;
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

 
 
 
 
 
 
 
    // for zut
    else if(datasetName == "FCNCzut_scaledown"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.910875 , 1.00395 , 0.979439 , 0.965395 , 0.981437 , 1.02583 , 1.06932 , 1.11198 , 1.01895 , 0.958027 , 0.91515 , 0.8691 , 0.925065 , 0.942151 , 0.91838 , 1.03337 , 0.996567 , 1.09232 , 1.1754 , 1.18656 , 1.36198 , 1.39598 , 1.50597 , 1.86258 , 1.68986 , 2.82751 , 3.69909 , 3.87895 , 3.95824 , 0  };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

 
    // for zct
    else if(datasetName == "FCNCzct_scaledown"){
     double reweightPrivateProd_tmp[31]= {  36.327 , 0.910983 , 0.695497 , 0.748991 , 0.964121 , 1.11941 , 1.12051 , 1.02468 , 1.0327 , 1.04627 , 1.08918 , 1.118 , 1.10832 , 1.07402 , 1.05069 , 1.05511 , 1.03379 , 1.00201 , 0.948706 , 0.982789 , 1.01898 , 0.946641 , 0.935019 , 1.02865 , 1.29136 , 0.998972 , 1.20211 , 1.22654 , 1.11757 , 2.23299 , 0  };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for kut
    else if(datasetName == "FCNCkut_matchup"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.762897 , 0.764029 , 0.884739 , 0.844435 , 0.76688 , 0.794838 , 0.817928 , 0.880097 , 0.977448 , 1.1276 , 1.23209 , 1.38688 , 1.4933 , 1.49885 , 1.74964 , 1.88273 , 2.009 , 2.23029 , 2.43118 , 2.99231 , 2.75765 , 3.4293 , 2.87149 , 6.48376 , 9.01983 , 5.99938 , 14.5937 , 6.6014 , 0 , 0  };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }




    // for kct
    else if(datasetName == "FCNCkct_matchup"){
     double reweightPrivateProd_tmp[31]= { 0 , 0.698017 , 0.719593 , 0.743384 , 0.747518 , 0.872695 , 1.15409 , 1.25711 , 1.21695 , 1.08732 , 1.03162 , 1.05062 , 1.09648 , 1.12552 , 1.21891 , 1.23022 , 1.24022 , 1.11835 , 1.17996 , 1.16561 , 1.31478 , 1.29517 , 1.3183 , 1.74909 , 1.5515 , 1.59628 , 1.85804 , 2.5827 , 2.54897 , 2.97265 , 0  } ;
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }







 
    // for zut
    else if(datasetName == "FCNCzut_matchup"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.841243 , 1.44438 , 1.45269 , 1.21106 , 1.00853 , 0.944302 , 0.886085 , 0.894557 , 0.850107 , 0.839717 , 0.849464 , 0.91427 , 0.941008 , 1.00291 , 1.05408 , 1.12743 , 1.12569 , 1.09497 , 1.07089 , 1.16715 , 1.08432 , 1.1993 , 1.22656 , 2.17072 , 1.57554 , 2.24958 , 2.0572 , 5.30425 , 2.4603 , 0 };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

 
    // for zct
    else if(datasetName == "FCNCzct_matchup"){
     double reweightPrivateProd_tmp[31]= {  17.5605 , 0.596561 , 0.688948 , 0.832575 , 1.11862 , 1.20398 , 1.13964 , 1.0901 , 0.986408 , 0.966402 , 0.967177 , 1.00395 , 0.965944 , 1.0281 , 1.07703 , 1.07344 , 1.05414 , 1.1098 , 1.05357 , 1.06254 , 1.03471 , 0.948667 , 1.14736 , 1.0849 , 1.23809 , 1.1112 , 1.27288 , 1.79786 , 1.71893 , 3.11834 , 0 };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for kut
    else if(datasetName == "FCNCkut_matchdown"){
     double reweightPrivateProd_tmp[31]= {  50.4888 , 0.707003 , 0.747977 , 0.818406 , 0.853953 , 0.949151 , 1.03146 , 1.02546 , 1.0191 , 0.981608 , 1.04264 , 1.05989 , 1.07686 , 1.16623 , 1.2507 , 1.26121 , 1.23155 , 1.25147 , 1.31989 , 1.35846 , 1.33442 , 1.33587 , 1.46197 , 2.1134 , 1.64292 , 1.74777 , 2.45278 , 2.28917 , 2.36364 , 2.52159 , 0  };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for kct
    else if(datasetName == "FCNCkct_matchdown"){
     double reweightPrivateProd_tmp[31]= { 0 , 0.663241 , 0.633681 , 0.901797 , 1.28888 , 1.18535 , 1.05819 , 0.918223 , 0.866225 , 0.86927 , 0.930297 , 0.967001 , 0.985765 , 1.06207 , 1.09894 , 1.13441 , 1.34766 , 1.41431 , 1.27042 , 1.46696 , 1.70405 , 1.54996 , 1.58161 , 1.93062 , 2.62567 , 1.91766 , 2.12583 , 5.17113 , 3.50872 , 0 , 0  } ; 
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }









 
    // for zut
    else if(datasetName == "FCNCzut_matchdown"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.722843 , 1.0616 , 1.05239 , 0.966535 , 0.911082 , 0.898049 , 0.885217 , 0.851363 , 0.887488 , 0.928911 , 1.0243 , 1.00966 , 1.17086 , 1.23454 , 1.32743 , 1.35913 , 1.35151 , 1.45697 , 1.6665 , 1.83476 , 1.97466 , 2.24081 , 2.87341 , 2.89334 , 3.01879 , 6.24678 , 6.51233 , 3.535 , 0 , 0 };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

 
    // for zct
    else if(datasetName == "FCNCzct_matchdown"){
     double reweightPrivateProd_tmp[31]= {  20.1409 , 0.846448 , 0.805092 , 0.784806 , 0.849833 , 0.985978 , 1.03221 , 1.05127 , 1.00932 , 0.991493 , 1.01662 , 1.04498 , 1.07889 , 1.13574 , 1.17261 , 1.1563 , 1.27704 , 1.23133 , 1.28447 , 1.13989 , 1.24961 , 1.26282 , 1.20286 , 1.23771 , 1.7388 , 1.85198 , 2.55487 , 6.39236 , 2.40964 , 2.01182 , 0  };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for kut
    else if(datasetName == "FCNCkut_topup"){
     double reweightPrivateProd_tmp[31]= {  0 , 1.03004 , 0.960673 , 1.08199 , 1.01634 , 0.966311 , 0.986895 , 0.979055 , 0.990711 , 0.975136 , 0.969891 , 0.960413 , 0.982091 , 1.02613 , 1.01741 , 1.00422 , 1.03575 , 0.986282 , 0.98827 , 1.00417 , 1.03863 , 1.01223 , 1.09841 , 1.31879 , 1.48158 , 1.664 , 2.55901 , 2.12651 , 3.16059 , 2.46288 , 0  };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }


    // for kct
    else if(datasetName == "FCNCkct_topup"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.578898 , 0.612173 , 0.783089 , 0.944304 , 1.08885 , 1.23533 , 1.3013 , 1.18122 , 1.0511 , 0.993603 , 0.956608 , 0.965738 , 0.980513 , 1.04356 , 1.02931 , 1.05583 , 1.07781 , 1.10004 , 1.08476 , 1.10286 , 1.17233 , 1.24038 , 1.48819 , 1.73382 , 1.97577 , 2.3626 , 3.03526 , 2.71852 , 3.15237 , 0 } ; 
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }






 
    // for zut
    else if(datasetName == "FCNCzut_topup"){
     double reweightPrivateProd_tmp[31]= {  0 , 1.35821 , 1.05842 , 0.962271 , 1.05702 , 1.07218 , 0.951359 , 0.925066 , 0.951599 , 0.933318 , 0.938659 , 0.906345 , 0.933653 , 0.948235 , 0.99131 , 1.0065 , 1.09457 , 1.12502 , 1.23387 , 1.17505 , 1.23879 , 1.24366 , 1.2617 , 1.58252 , 1.68672 , 1.85151 , 2.48519 , 3.55041 , 3.33557 , 2.47546 , 0  };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

     
    // for zct
    else if(datasetName == "FCNCzct_topup"){
     double reweightPrivateProd_tmp[31]= {  41.5534 , 0.764022 , 0.718226 , 0.860688 , 0.876271 , 0.893366 , 0.992805 , 1.09479 , 1.12378 , 1.07315 , 1.08742 , 1.10158 , 1.06466 , 1.09712 , 1.06265 , 1.14782 , 1.17845 , 1.13282 , 1.14177 , 1.13039 , 1.14367 , 1.18426 , 1.31807 , 1.48169 , 1.47715 , 1.94077 , 2.15634 , 2.27384 , 1.44331 , 2.55425 , 0  };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for kut
    else if(datasetName == "FCNCkut_topdown"){
     double reweightPrivateProd_tmp[31]= {  9.125 , 0.980759 , 1.14075 , 0.959074 , 0.875701 , 0.884457 , 0.897158 , 0.940597 , 0.952344 , 0.979812 , 1.03483 , 1.11915 , 1.12613 , 1.17311 , 1.12708 , 1.10673 , 1.16769 , 1.15769 , 1.16511 , 1.18459 , 1.10422 , 1.12766 , 1.23306 , 1.36658 , 1.40733 , 1.4026 , 1.71583 , 1.71803 , 1.6775 , 1.96317 , 0  };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for kct
    else if(datasetName == "FCNCkct_topdown"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.60354 , 0.6397 , 0.986734 , 1.26481 , 1.18199 , 1.0277 , 0.902555 , 0.858423 , 0.838187 , 0.954761 , 0.993142 , 1.06956 , 1.17508 , 1.13329 , 1.14773 , 1.24981 , 1.17382 , 1.39394 , 1.23581 , 1.29021 , 1.41998 , 1.54785 , 1.81228 , 2.04698 , 2.54562 , 5.0144 , 5.03396 , 3.07408 , 3.25913 , 0 } ; 
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }





 
    // for zut
    else if(datasetName == "FCNCzut_topdown"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.608823 , 0.733166 , 0.784112 , 0.911624 , 0.957957 , 1.00966 , 1.07282 , 1.20318 , 1.24168 , 1.15363 , 1.14219 , 1.03602 , 1.04011 , 1.08254 , 1.08701 , 1.05019 , 1.06686 , 1.05638 , 1.03762 , 0.938022 , 1.0018 , 1.00685 , 1.30297 , 1.47545 , 1.20477 , 1.72019 , 1.66048 , 1.22909 , 1.50506 , 0  };
      reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

 
    // for zct
    else if(datasetName == "FCNCzct_topdown"){
     double reweightPrivateProd_tmp[31]= {  28.6416 , 0.932728 , 0.772171 , 0.824685 , 1.00605 , 0.995949 , 0.940026 , 0.908997 , 0.980909 , 1.00116 , 1.02356 , 1.09983 , 1.08679 , 1.14295 , 1.22058 , 1.14861 , 1.22009 , 1.24146 , 1.27108 , 1.15123 , 1.31055 , 1.16199 , 1.07328 , 1.19028 , 1.40885 , 1.53229 , 1.21106 , 2.27258 , 1.92749 , 2.86093 , 0  };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for WZ private
    else if(datasetName == "WZprivate"){
     double reweightPrivateProd_tmp[31]= {  5.60045 , 1.05978 , 0.79694 , 0.781762 , 0.883931 , 0.968911 , 0.981463 , 1.00412 , 1.07303 , 1.1176 , 1.14079 , 1.16191 , 1.13931 , 1.06695 , 1.10416 , 1.0904 , 1.11693 , 1.07321 , 1.01076 , 0.984842 , 0.968539 , 0.921287 , 0.883743 , 1.04358 , 1.16515 , 1.24631 , 1.27875 , 1.48123 , 1.50757 , 1.49177 , 0 };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for WZ private matchup
    else if(datasetName == "WZprivate_matchup"){
     double reweightPrivateProd_tmp[31]= {  3.97022 , 0.641981 , 0.735831 , 0.851583 , 0.850335 , 0.942161 , 1.07982 , 1.07316 , 0.989848 , 0.940942 , 0.949548 , 0.990726 , 1.04251 , 1.14313 , 1.21043 , 1.32572 , 1.37966 , 1.38649 , 1.50986 , 1.57969 , 1.60366 , 1.60547 , 1.61868 , 2.59631 , 2.72351 , 3.26013 , 4.72966 , 4.72528 , 4.27493 , 9.51777 , 0 };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for WZ private matchdown
    else if(datasetName == "WZprivate_matchdown"){
     double reweightPrivateProd_tmp[31]= {  4.97165 , 0.943345 , 0.995265 , 1.11425 , 1.1736 , 0.966863 , 0.86165 , 0.806637 , 0.863769 , 0.927898 , 1.02027 , 1.06075 , 1.05789 , 1.14124 , 1.133 , 1.1397 , 1.13974 , 1.16056 , 1.16869 , 1.11429 , 1.14865 , 1.18922 , 1.21618 , 1.40092 , 1.31446 , 1.69609 , 1.56576 , 2.25416 , 1.86199 , 3.97283 , 0  };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for WZ private scaleup
    else if(datasetName == "WZprivate_scaleup"){
     double reweightPrivateProd_tmp[31]= {  3.8893 , 0.758531 , 0.851 , 0.914643 , 1.01246 , 0.991989 , 0.903277 , 0.87064 , 0.890912 , 0.93049 , 0.979772 , 1.03688 , 1.14387 , 1.2145 , 1.2883 , 1.31274 , 1.31185 , 1.35892 , 1.45586 , 1.33054 , 1.35821 , 1.30395 , 1.35144 , 1.4297 , 1.41312 , 1.72276 , 2.26047 , 2.16019 , 2.66497 , 2.29005 , 0 };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for WZ private scaledown
    else if(datasetName == "WZprivate_scaledown"){
     double reweightPrivateProd_tmp[31]= {  3.8805 , 0.673692 , 0.823117 , 0.879461 , 0.944958 , 0.992174 , 0.943534 , 0.906341 , 0.929292 , 0.934526 , 0.989993 , 1.07523 , 1.10036 , 1.20668 , 1.25037 , 1.31577 , 1.36481 , 1.35275 , 1.39318 , 1.42934 , 1.43312 , 1.41591 , 1.40472 , 1.69364 , 1.93794 , 2.33865 , 2.35352 , 2.22736 , 1.91977 , 2.92862 , 0  };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for  TZq
    else if(datasetName == "TZq"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.973806 , 0.887764 , 0.822885 , 0.872453 , 0.939798 , 0.940942 , 0.99365 , 1.03431 , 1.05627 , 1.05617 , 1.14854 , 1.05693 , 1.12571 , 1.08706 , 1.05765 , 1.12804 , 1.15374 , 1.24244 , 1.21423 , 1.19548 , 1.11615 , 1.25506 , 1.3633 , 1.38904 , 1.48147 , 3.832 , 2.99618 , 2.03297 , 1.50874 , 0 };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for  TZq matchup
    else if(datasetName == "TZq_matchup"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.826837 , 0.812389 , 0.84089 , 0.865891 , 0.893724 , 0.955301 , 1.00185 , 1.03626 , 1.03305 , 1.01178 , 1.00955 , 1.05831 , 1.11531 , 1.15173 , 1.19961 , 1.31234 , 1.32658 , 1.39132 , 1.44248 , 1.62267 , 1.58019 , 1.82971 , 2.03203 , 2.33843 , 2.62749 , 7.57307 , 5.26334 , 8.33299 , 18.5527 , 0  };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

    // for  TZq matchdown
    else if(datasetName == "TZq_matchdown"){
     double reweightPrivateProd_tmp[31]= {  0 , 0.88192 , 0.741805 , 0.887731 , 0.999574 , 1.02477 , 1.1351 , 0.971376 , 1.00494 , 0.946958 , 0.988118 , 1.0501 , 1.02462 , 1.02386 , 1.02603 , 0.975517 , 1.0667 , 1.18406 , 1.15494 , 1.11834 , 1.00729 , 1.19381 , 2.10356 , 1.58955 , 2.16182 , 1.5037 , 1.63359 , 5.6768 , 1.28394 , 1.42929 , 0 };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }
    
    else
    {
     double reweightPrivateProd_tmp[31]= {  1 ,1 , 1 , 1 , 1 , 1 ,1 , 1 , 1.  , 1 ,1 , 1 ,1 , 1, 1 , 1 ,1 ,1 , 1  , 1 , 1.  , 1.  ,1, 1 , 1 , 1 , 1 , 1 , 1  , 1  , 1 };
     reweightPrivateProd = TableToVector(reweightPrivateProd_tmp, 31);
    }

 

    return reweightPrivateProd;

}


std::vector<double> ProofSelectorMyCutFlow::TableToVector(double *theTable, int size){



   std::vector<double> thevector;
   
   
   for(int i=0; i<size; i++){
   
     thevector.push_back(theTable[i]);
   } 



  return thevector;



}













