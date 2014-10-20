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
#include "ProofSelectorMyCutFlowTools.h"

//_____________________________________________________________________________
ProofSelectorMyCutFlow::ProofSelectorMyCutFlow()
{
  // Constructor
  
  

  cout << "start proof constructor " << endl;
  fChain     = 0;
  branch     = 0;
  event      = 0;
  dataset    = 0;
  anaEL      = 0;
  verbosity  = 0;
  DataType   = 0;
  Luminosity = 0;
  //histos
  //fHist      = 0;
  
  nselevents_mumumu = 0;
  isTTbarMCatNLO = false;
  
  //------------------------// 
  //initialize the variables
  //------------------------// 
  //------------------------// 
  //for PU
  IReweight             = true;
  IReweight_puUp	= false;
  IReweight_puDown	= false;
  IReweight_pu	        = true;
  //------------------------// 
  //for JES uncertainties
  doJESuncert = false;
  upOrDown    = false;
  applyJES    = false; 
  applyJER    = false;
  ResFactor   = 0;
  //------------------------// 
  //for PDF unceratinties
  doPDF = false;
  //pdftype =0 ;
  pdftype =1  ;
  //------------------------// 
  //loose iso on W
  //for backgr studies
  useNonIsoWcand = false;
  looseIso = 0.35; //0.4
  themetcut = 1000000;
  //------------------------// 
  //tight iso for W cand.
  tightIso_e  = 0.15;//0.15;
  tightIso_mu = 0.20;//0.20;
  elecIso_default = 0.15;
  //------------------------// 
  applyTrigger     = true;
  applyTriggerUp   = false;
  applyTriggerDown = false;
  //------------------------// 
  applyLeptonSF      = true;
  applyLeptonSFUp    = false;
  applyLeptonSFDown  = false;
  //------------------------// 
  


  //------------------------// 
  //initialize scale factors
  //------------------------// 
  //trigger SF
  /*SF_trig_mumumu = 0.967; 
  SF_trig_mumue  = 0.953; 
  SF_trig_eemu   = 0.953; 
  SF_trig_eee    = 0.974;
  SF_trig_mumumu_error  = 0.012; 
  SF_trig_mumue_error   = 0.010;
  SF_trig_eemu_error    = 0.010; 
  SF_trig_eee_error     = 0.010; */
  SF_trig_mumumu = 1; 
  SF_trig_mumue  = 1 ; 
  SF_trig_eemu   = 1 ; 
  SF_trig_eee    = 1 ;
  SF_trig_mumumu_error  = 0.012; 
  SF_trig_mumue_error   = 0.010;
  SF_trig_eemu_error    = 0.010; 
  SF_trig_eee_error     = 0.010; 
  //------------------------// 
  //for WZ+jets
  applyWZ          = true;
  SF_WZ.push_back(1.); //mumumu
  SF_WZ.push_back(1.); //mumue
  SF_WZ.push_back(1.); //eemu
  SF_WZ.push_back(1.); //eee
  //------------------------// 
  //for fake leptons
  applyFakescale   = true;
  SF_Fake.push_back(1.); //mumumu
  SF_Fake.push_back(1.); //mumue
  SF_Fake.push_back(1.); //eemu
  SF_Fake.push_back(1.); //eee
  
  
  cout << "end proof constructor " << endl;

}

//_____________________________________________________________________________
ProofSelectorMyCutFlow::~ProofSelectorMyCutFlow()
{
  // Destructor
  cout << "in destructor " << endl;
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
  //tree->cd("MyModule);
  
  if(!tree) cout << "tree with null pointer " << endl;
  cout << "	GetNbranches " << tree->GetNbranches() << endl;
  branch = (TBranch *) tree->GetBranch("NTEvent");
  //branch = (TBranch *) tree->GetBranch("NTSampleInfo");
  
  event = new IPHCTree::NTEvent();
  
  if(!branch) cout << "get branch failed " << endl;
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
  anaEL->LoadGeneralInfo(DataType, Luminosity, LumiError, PUWeightFileName, verbosity );
  
   //******************************************
  //Load Scale Factors for lepton efficiencies
  //******************************************


   cout << "load muon SF   " << endl,
   sel.LoadMuIDScaleFactors();
   sel.LoadMuIsoScaleFactors12();
   sel.LoadMuIsoScaleFactors20();
   cout << "muon SF loaded " << endl,


  
  
  anaEL->LoadWeight (sel); // now the parameters for SFBweight are initialized (for b-tag!)
  
  
  //--------------------------------------//
  //   retrive datasets	
  //--------------------------------------//
  
  
  for(unsigned int d=0;d<datasets.size();d++){
    cout << "datasets.size() " << datasets.size()<< "  datasets[d].Name()" << datasets[d].Name()  << " datasetName "
	 <<datasetName  << endl;
    if(datasets[d].Name()==datasetName)dataset = &datasets[d];
  }
  
  
  //--------------------------------------//
  //   Fill cuts and channels  	
  //--------------------------------------//
  CutName.push_back("Cut1");
  
  
  //--------------------------------------//
  //   Determine is dataset or MC samples 	
  //--------------------------------------//
  bool isData_sample = false;    
  if(datasetName=="DataEG" || datasetName=="DataMu" || 
     datasetName=="DataMuEG" || datasetName=="DataEGMu" ||
     datasetName=="DataMuZenriched" || datasetName=="DataMuEGZenriched" || 
     datasetName=="DataEGZenriched" ||
     datasetName=="MET1" || datasetName=="MET2") isData_sample = true;
  
    
  
  //--------------------------------------//
  //   Managing histos  	
  //--------------------------------------//
  cout << "load  MyhistoManager " << endl;
  cout << "dataset " <<  datasetName<< endl;
  cout << "the channel " << TheChannelName.size() << endl;
  MyhistoManager.LoadDatasets(datasets);   
  MyhistoManager.LoadSelectionSteps(CutName);
  MyhistoManager.LoadChannels(TheChannelName);
  cout << " MyhistoManager loaded" << endl;
  
  ITypeMC = -1 ;
  
  
  
  createTheHisto(&MyhistoManager);
  
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
  //fFile = fProofFile->OpenFile("NEW");
  if (fFile && fFile->IsZombie()) SafeDelete(fFile);
  savedir->cd();
  fFile->cd();

  
  //--------------------------------------//
  //   Output TTree 	
  //--------------------------------------//
  TString treename = "Ttree_"+datasetName;
  TheTree = new TTree(treename.Data(),treename.Data());
  TheTree->Branch("tree_cosThetaStar",     &tree_cosThetaStar,     "tree_cosThetaStar/F"   );
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
  
  
  TheTree->Branch("tree_met", &tree_met , "tree_met/F" );
  TheTree->Branch("tree_mTW", &tree_mTW , "tree_mTW/F" );
  
  
  TheTree->Branch("tree_EvtWeight",   &tree_EvtWeight,   "tree_EvtWeight/F" );
  TheTree->Branch("tree_SampleType",  &tree_SampleType,  "tree_SampleType/I");
  TheTree->Branch("tree_Channel",     &tree_Channel,     "tree_Channel/I"   );

   
  
   //--------------------------------------//
  //   Output TTree 	
  //--------------------------------------//
  TString smalltreename = "SmallTree_"+datasetName;
  SmallTree = new TTree(smalltreename.Data(),smalltreename.Data());
   
  SmallTree->Branch("smalltree_nlepton",    &smalltree_nlepton,     "smalltree_nlepton/I");
  SmallTree->Branch("smalltree_lept_pt",     smalltree_lept_pt,     "smalltree_lept_pt[smalltree_nlepton]/F");
  SmallTree->Branch("smalltree_lept_eta",    smalltree_lept_eta,    "smalltree_lept_eta[smalltree_nlepton]/F");
  SmallTree->Branch("smalltree_lept_phi",    smalltree_lept_phi,    "smalltree_lept_phi[smalltree_nlepton]/F");
  SmallTree->Branch("smalltree_lept_iso",    smalltree_lept_iso,    "smalltree_lept_iso[smalltree_nlepton]/F");
  SmallTree->Branch("smalltree_lept_flav",   smalltree_lept_flav,   "smalltree_lept_flav[smalltree_nlepton]/I");
   
  SmallTree->Branch("smalltree_njets",          &smalltree_njets,         "smalltree_njets/I");
  SmallTree->Branch("smalltree_jet_pt",         smalltree_jet_pt,         "smalltree_jet_pt[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_eta",        smalltree_jet_eta,        "smalltree_jet_eta[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_phi",        smalltree_jet_phi,        "smalltree_jet_phi[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_btagdiscri", smalltree_jet_btagdiscri, "smalltree_jet_btagdiscri[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_btagdiscri_up",   smalltree_jet_btagdiscri_up,   "smalltree_jet_btagdiscri_up[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_btagdiscri_down", smalltree_jet_btagdiscri_down, "smalltree_jet_btagdiscri_down[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_flav",       smalltree_jet_flav,       "smalltree_jet_flav[smalltree_njets]/I");
  
  
  SmallTree->Branch("smalltree_jesup_njets",          &smalltree_jesup_njets,        "smalltree_jesup_njets/I");
  SmallTree->Branch("smalltree_jet_jesup_pt",         smalltree_jet_jesup_pt,        "smalltree_jet_jesup_pt[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jesup_eta",        smalltree_jet_jesup_eta,       "smalltree_jet_jesup_eta[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jesup_phi",        smalltree_jet_jesup_phi,       "smalltree_jet_jesup_phi[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jesup_btagdiscri", smalltree_jet_jesup_btagdiscri,"smalltree_jet_jesup_btagdiscri[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jesup_flav",       smalltree_jet_jesup_flav,      "smalltree_jet_jesup_flav[smalltree_njets]/I");
   
  SmallTree->Branch("smalltree_jesdown_njets",          &smalltree_jesdown_njets,        "smalltree_jesdown_njets/I");
  SmallTree->Branch("smalltree_jet_jesdown_pt",         smalltree_jet_jesdown_pt,        "smalltree_jet_jesdown_pt[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jesdown_eta",        smalltree_jet_jesdown_eta,       "smalltree_jet_jesdown_eta[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jesdown_phi",        smalltree_jet_jesdown_phi,       "smalltree_jet_jesdown_phi[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jesdown_btagdiscri", smalltree_jet_jesdown_btagdiscri,"smalltree_jet_jesdown_btagdiscri[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jesdown_flav",       smalltree_jet_jesdown_flav,      "smalltree_jet_jesdown_flav[smalltree_njets]/I");
 
    
  
  SmallTree->Branch("smalltree_jerup_njets",          &smalltree_jerup_njets,           "smalltree_jerup_njets/I");
  SmallTree->Branch("smalltree_jet_jerup_pt",         smalltree_jet_jerup_pt,           "smalltree_jet_jerup_pt[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jerup_eta",        smalltree_jet_jerup_eta,          "smalltree_jet_jerup_eta[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jerup_phi",        smalltree_jet_jerup_phi,          "smalltree_jet_jerup_phi[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jerup_btagdiscri", smalltree_jet_jerup_btagdiscri,   "smalltree_jet_jerup_btagdiscri[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jerup_flav",       smalltree_jet_jerup_flav,         "smalltree_jet_jerup_flav[smalltree_njets]/I");
   
  SmallTree->Branch("smalltree_jerdown_njets",          &smalltree_jerdown_njets,         "smalltree_jerdown_njets/I");
  SmallTree->Branch("smalltree_jet_jerdown_pt",         smalltree_jet_jerdown_pt,         "smalltree_jet_jerdown_pt[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jerdown_eta",        smalltree_jet_jerdown_eta,        "smalltree_jet_jerdown_eta[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jerdown_phi",        smalltree_jet_jerdown_phi,        "smalltree_jet_jerdown_phi[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jerdown_btagdiscri", smalltree_jet_jerdown_btagdiscri, "smalltree_jet_jerdown_btagdiscri[smalltree_njets]/F");
  SmallTree->Branch("smalltree_jet_jerdown_flav",       smalltree_jet_jerdown_flav,       "smalltree_jet_jerdown_flav[smalltree_njets]/I");


  SmallTree->Branch("smalltree_met_jesup_pt",   &smalltree_met_jesup_pt,    "smalltree_met_jesup_pt/F");
  SmallTree->Branch("smalltree_met_jesup_phi",  &smalltree_met_jesup_phi,   "smalltree_met_jesup_phi/F");
  
  SmallTree->Branch("smalltree_met_jesdown_pt",   &smalltree_met_jesdown_pt,    "smalltree_met_jesdown_pt/F");
  SmallTree->Branch("smalltree_met_jesdown_phi",  &smalltree_met_jesdown_phi,   "smalltree_met_jesdown_phi/F");
  
  SmallTree->Branch("smalltree_met_jerup_pt",   &smalltree_met_jerup_pt,    "smalltree_met_jerup_pt/F");
  SmallTree->Branch("smalltree_met_jerup_phi",  &smalltree_met_jerup_phi,   "smalltree_met_jerup_phi/F");
  
  SmallTree->Branch("smalltree_met_jerdown_pt",   &smalltree_met_jerdown_pt,    "smalltree_met_jerdown_pt/F");
  SmallTree->Branch("smalltree_met_jerdown_phi",  &smalltree_met_jerdown_phi,   "smalltree_met_jerdown_phi/F");
  
  
  
  SmallTree->Branch("smalltree_met_unclsup_pt",   &smalltree_met_unclsup_pt,    "smalltree_met_unclsup_pt/F");
  SmallTree->Branch("smalltree_met_unclsup_phi",  &smalltree_met_unclsup_phi,   "smalltree_met_unclsup_phi/F");
  
  SmallTree->Branch("smalltree_met_unclsdown_pt",   &smalltree_met_unclsdown_pt,    "smalltree_met_unclsdown_pt/F");
  SmallTree->Branch("smalltree_met_unclsdown_phi",  &smalltree_met_unclsdown_phi,   "smalltree_met_unclsdown_phi/F");
  
  
  SmallTree->Branch("smalltree_met_pt",   &smalltree_met_pt,    "smalltree_met_pt/F");
  SmallTree->Branch("smalltree_met_phi",  &smalltree_met_phi,   "smalltree_met_phi/F");
  
  
  SmallTree->Branch("smalltree_weight_trigup",   &smalltree_weight_trigup,   "smalltree_weight_trigup/F");
  SmallTree->Branch("smalltree_weight_trigdown", &smalltree_weight_trigdown, "smalltree_weight_trigdown/F");
  
  
  SmallTree->Branch("smalltree_weight_leptup",     &smalltree_weight_leptup,     "smalltree_weight_leptup/F");
  SmallTree->Branch("smalltree_weight_leptdown",   &smalltree_weight_leptdown,   "smalltree_weight_leptdown/F");
  
  SmallTree->Branch("smalltree_weight_PDFup",     &smalltree_weight_PDFup,     "smalltree_weight_PDFup/F");
  SmallTree->Branch("smalltree_weight_PDFdown",   &smalltree_weight_PDFdown,   "smalltree_weight_PDFdown/F");
  
  SmallTree->Branch("smalltree_evtweight",   &smalltree_evtweight,   "smalltree_evtweight/F");
 
   
  SmallTree->Branch("smalltree_IChannel",    &smalltree_IChannel,     "smalltree_IChannel/I");
   
   //--------------------------------------//
   //   Initialize PU reweighting for 8 TeV 	
   //--------------------------------------//
  
  string PUdatafilename;
  TH1D * dataPUhisto;
  if( !isData_sample && IReweight_pu) {
    if( !IReweight_puDown && !IReweight_puUp ) PUdatafilename = getenv( "CMSSW_BASE" )+string("/src/IPHCAnalysis/NTuple/macros/data/PileupHistogram2012_MuEG.root");
    if( IReweight_puDown )                     PUdatafilename = getenv( "CMSSW_BASE" )+string("/src/IPHCAnalysis/NTuple/macros/data/");
    if( IReweight_puUp )                       PUdatafilename = getenv( "CMSSW_BASE" )+string("/src/IPHCAnalysis/NTuple/macros/data/");
  }
  else {
    if( !IReweight_puDown && !IReweight_puUp ) PUdatafilename = getenv( "CMSSW_BASE" )+string("/src/IPHCAnalysis/NTuple/macros/data/PileupHistogram2012_MuEG.root");
    if( IReweight_puDown )                     PUdatafilename = getenv( "CMSSW_BASE" )+string("/src/IPHCAnalysis/NTuple/macros/data/");
    if( IReweight_puUp )                       PUdatafilename = getenv( "CMSSW_BASE" )+string("/src/IPHCAnalysis/NTuple/macros/data/");
  }
  fexists(PUdatafilename, true);
  
  TFile* pufile = new TFile(PUdatafilename.c_str() );
  dataPUhisto = (TH1D*) pufile->Get( "pileup" );
  LumiWeights = new  PUWeighting( );
  LumiWeights->initPUSummer12_S10(&*dataPUhisto);
  
  
  
  //--------------------------------------//
  //   For JET uncertainties 	
  //--------------------------------------//
  
  
  //cout << " isData_sample  " << isData_sample << endl;
  string jesuncertcalibpath = "";
  if(isData_sample){
    jesuncertcalibpath = getenv( "CMSSW_BASE" )+string("/src/IPHCAnalysis/NTuple/macros/data/Summer13_V5_DATA_UncertaintySources_AK5PFchs.txt");
    JetCorrectorParameters *pjetparam = new JetCorrectorParameters( jesuncertcalibpath.c_str(), "Total");
    theJESuncertainty = new JetCorrectionUncertainty( *pjetparam );
  }
  else{
    jesuncertcalibpath = getenv( "CMSSW_BASE" )+string("/src/IPHCAnalysis/NTuple/macros/data/Fall12_V7_MC_Uncertainty_AK5PFchs.txt");
    theJESuncertainty = new JetCorrectionUncertainty( jesuncertcalibpath.c_str() );    
  }
  
  
  
  //--------------------------------------//
  //   for PDF uncertainties	
  //--------------------------------------//
  
  
  if(doPDF) pdf.Initialize();
  
 
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
     datasetName=="DataEGZenriched" || datasetName=="DataMuZenriched" || 
     datasetName=="DataMuEGZenriched" ||
     datasetName=="MET1" || datasetName=="MET2") isData = true;
  
  
 
  //---------------------------------
  //load events 
  //---------------------------------
  sel.LoadEvent(event);
  
  
  
  //---------------------------------
  //initialize output tree
  //---------------------------------
   
  smalltree_nlepton         = -1000;  
  smalltree_njets           = -1000;  
  smalltree_jesup_njets     = -1000;    
  smalltree_jesdown_njets   = -1000;   
  smalltree_jerup_njets     = -1000;  
  smalltree_jerdown_njets   = -1000;    
  smalltree_met_pt          = -1000;
  smalltree_met_phi         = -1000;
  smalltree_met_jesup_pt    = -1000;
  smalltree_met_jesup_phi   = -1000;
  smalltree_met_jesdown_pt  = -1000;
  smalltree_met_jesdown_phi = -1000;
  smalltree_met_jerup_pt    = -1000;
  smalltree_met_jerup_phi   = -1000;
  smalltree_met_jerdown_pt  = -1000;
  smalltree_met_jerdown_phi = -1000;
  smalltree_met_unclsup_pt   = -1000;
  smalltree_met_unclsup_phi  = -1000;
  smalltree_met_unclsdown_pt = -1000;
  smalltree_met_unclsdown_phi= -1000;
  smalltree_evtweight       = -1000;
  smalltree_weight_trigup   = -1000;
  smalltree_weight_trigdown = -1000;
  smalltree_weight_leptup   = -1000;
  smalltree_weight_leptdown = -1000;
  smalltree_weight_PDFup    = -1000;
  smalltree_weight_PDFdown  = -1000;
  smalltree_IChannel        = -1000;
  
  //---------------------------------
  //get collections from events 
  //---------------------------------
  double rho = event->pileup.rho_PUUE_dens;
  
  //Collection of selected objects
  vector<NTVertex>   selVertices        = sel.GetSelectedVertex();
  vector<NTElectron> selElectrons       = sel.GetSelectedElectronsRhoIso(20, 2.5, elecIso_default, 0, 0, 0, 0, rho); 
  vector<NTMuon>     selMuons           = sel.GetSelectedMuonsDeltaBetaIso();
  vector<NTElectron> selElectronsNonIso = sel.GetSelectedElectronsNoIsoDileptonTTbar(20, 2.5, 0,  0, 0, 0);
  vector<NTMuon>     selMuonsNonIso     = sel.GetSelectedMuonsNoIsoDileptonTTbar();
  
  
  
  //FIXME add JER
  
  
  //selection of veto leptons
  vector<NTElectron> looseElectrons ;
  vector<NTMuon>     looseMuons	    ;
  
  vector<NTElectron> theelectrons   = *(sel.GetPointer2Electrons());
  vector<NTMuon>     themuons	    = *(sel.GetPointer2Muons());
  
  
  
  for(unsigned int i=0; i<themuons.size(); i++){ 
    
    if(!themuons[i].isGlobalMuon  && !themuons[i].isTrackerMuon )  continue;
    if(themuons[i].p4.Pt() < 10 )                                  continue;
    if( fabs(themuons[i].p4.Eta()) > 2.4 )                         continue;
    if(sel.RelIso04PFDeltaBeta( themuons[i]  ) > 0.20)             continue;
    looseMuons.push_back(themuons[i]);
      
  }
  for(unsigned int i=0; i<theelectrons.size(); i++){ 
    if(theelectrons[i].p4Gsf.Pt()< 10 || fabs(theelectrons[i].p4Gsf.Eta()) > 2.5 )  continue;
    if(theelectrons[i].ID["mvaTrigV0"] < 0  ) continue;
    if(sel.EffArea03PF(theelectrons[i], rho) > 0.15)                                continue;
    looseElectrons.push_back(theelectrons[i]);  

  } 
  
 

  
   //---------------------------------
  //initiate lepton candidate collections 
  //---------------------------------
  vector<NTElectron> ZeeCand; 
  vector<NTMuon>     ZmumuCand; 
  
  vector<NTElectron> WeCand; 
  vector<NTMuon>     WmuCand; 






      

  tree_cosThetaStar= -10000;
  tree_EvtWeight  = -10000;    
  tree_topMass    = -10000;
  tree_totMass    = -10000; 
  tree_deltaPhilb = -10000; 
  tree_deltaRlb   = -10000; 
  tree_deltaRTopZ = -10000; 
  tree_asym	  = -10000; 
  tree_Zpt	  = -10000; 
  tree_ZEta	  = -10000; 
  tree_topPt	  = -10000; 
  tree_topEta	  = -10000; 
  tree_totPt	    =  -10000;
  tree_totEta	    = -10000; 
  tree_deltaRZl     = -10000; 
  tree_deltaPhiZmet = -10000; 
  tree_btagDiscri   = -10000; 
  tree_NJets	    = -10000; 
  tree_NBJets	      =  -10000;
  tree_leptWPt        = -10000; 
  tree_leptWEta       =  -10000;
  tree_leadJetPt      = -10000; 
  tree_leadJetEta     = -10000; 
  tree_deltaRZleptW   = -10000; 
  tree_deltaPhiZleptW = -10000; 
  tree_met =  -10000;
  tree_mTW =  -10000;
  tree_Channel = -10000; 
  tree_SampleType= -10000; 
  
  setNullHistoPointer();
  

  //---------------------------------
  //initialize MC weights 
  //---------------------------------
  double Dweight[101];
  for(int k1=0; k1<101; k1++) {
    Dweight[k1] = 0.;
  }
  
  double weightITypeMC_save = Luminosity*dataset->Xsection()/dataset->getNSkimmedEvent();
  double weightITypeMC=0;
  
 
  
      
  
  
  //*****************************************************************
  // calcul the MC weights
  //*****************************************************************
  //to do : use isData parameter
  if ( datasetName!="DataMu" && datasetName!="DataEG" &&
       datasetName!="DataMuEG" && datasetName!="DataMuZenriched" && datasetName!="DataEGZenriched" &&
       datasetName!="DataMuEGZenriched" ) { 
    if(IReweight ){
      weightITypeMC = weightITypeMC_save*LumiWeights->weight_Summer12_S10(event->pileup.Tnpv);
    }
    else weightITypeMC = weightITypeMC_save;
  }
  else weightITypeMC = 1;
  
  
 
  
   if(datasetName=="DataMuZenriched" || datasetName=="DataEGZenriched" ||
       datasetName=="DataMuEGZenriched")   useNonIsoWcand = true;
       else  useNonIsoWcand = false;
  
  
  //*****************************************************************
  // determine top decay channel
  //*****************************************************************    
  bool IsTTbarDilept = false;
  bool IsSignal      = false;
  double WeightForBranchingRatio = 1.;
  
  
  
  //*****************************************************************
  // determine MC evetn weight
  //*****************************************************************
  
  std::vector< double > thereweight = ProofSelectorMyCutFlowTools::determineWeights(datasetName, weightITypeMC, WeightForBranchingRatio);
  ITypeMC = thereweight[0];
  Dweight[ITypeMC] = thereweight[1];
  double EventYieldWeightError = thereweight[2];
  if(thereweight[3] > 0) IsSignal = true; else IsSignal = false;
  
  
  
  //*****************************************************************
  // pass trigger selection
  //*****************************************************************
  
  bool passtrigger_mumu = false;
  bool passtrigger_emu  = false;
  bool passtrigger_ee   = false;
  
     

 
  
  //if( Dweight[ITypeMC] != 1) cout << "norm issue" << endl;
  
  //to do : update trigger selection
  passtrigger_mumu   = sel.passTriggerSelection8TeV ( dataset, "mumu");
  passtrigger_emu    = sel.passTriggerSelection8TeV ( dataset, "emu" );
  passtrigger_ee     = sel.passTriggerSelection8TeV ( dataset, "ee" );
  
  
  int  IChannel= -1;
  
  if(isData){
    if(passtrigger_mumu  && !passtrigger_emu && !passtrigger_ee &&   (datasetName=="DataMu"   || datasetName== "DataMuZenriched"  ) )    IChannel = 0;
    if(passtrigger_emu   &&                                          (datasetName=="DataMuEG" || datasetName== "DataMuEGZenriched") )  IChannel = 12;
    if(passtrigger_ee    && !passtrigger_emu && !passtrigger_mumu && (datasetName=="DataEG"   || datasetName== "DataEGZenriched"  ) )    IChannel = 3;
  }else{
    if( passtrigger_mumu && !passtrigger_emu && !passtrigger_ee ) IChannel = 0;
    if( passtrigger_emu                                         ) IChannel = 12;
    if(!passtrigger_mumu && !passtrigger_emu &&  passtrigger_ee ) IChannel = 3;
  }
  
 
  

      
  
  if( IChannel != -1 ){
  

    
    
    
    if(IChannel == 0)  MyhistoManager.FillHisto(CutFlow_mumumu, "CutFlow_mumumu", 0, datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 12) MyhistoManager.FillHisto(CutFlow_mumue,  "CutFlow_mumue" , 0, datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 12) MyhistoManager.FillHisto(CutFlow_eemu,   "CutFlow_eemu"  , 0, datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 3)  MyhistoManager.FillHisto(CutFlow_eee,    "CutFlow_eee"   , 0, datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 0)  MyhistoManager.FillHisto(ErrCutFlow_mumumu,   "ErrCutFlow_mumumu" , 0, datasetName, IsSignal, EventYieldWeightError);
    if(IChannel == 12) MyhistoManager.FillHisto(ErrCutFlow_mumue,    "ErrCutFlow_mumue"  , 0, datasetName, IsSignal, EventYieldWeightError);
    if(IChannel == 12) MyhistoManager.FillHisto(ErrCutFlow_eemu,     "ErrCutFlow_eemu"	, 0, datasetName, IsSignal, EventYieldWeightError);
    if(IChannel == 3)  MyhistoManager.FillHisto(ErrCutFlow_eee,      "ErrCutFlow_eee"	, 0, datasetName, IsSignal, EventYieldWeightError);
    
    if(IChannel == 0)  MyhistoManager.FillHisto(NVtx_mumumu_aftertrigsel, "NVtx_mumumu_aftertrigsel", selVertices.size(), datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 12) MyhistoManager.FillHisto(NVtx_mumue_aftertrigsel,  "NVtx_mumue_aftertrigsel",  selVertices.size(), datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 12) MyhistoManager.FillHisto(NVtx_eemu_aftertrigsel,   "NVtx_eemu_aftertrigsel",   selVertices.size(), datasetName, IsSignal, Dweight[ITypeMC]);
    if(IChannel == 3)  MyhistoManager.FillHisto(NVtx_eee_aftertrigsel,    "NVtx_eee_aftertrigsel",    selVertices.size(), datasetName, IsSignal, Dweight[ITypeMC]);


      
    
    ProofSelectorMyCutFlowTools::determineLeptonCandidates(useNonIsoWcand, looseIso, rho ,
			      &selElectrons,       &selMuons, 
			      &selElectronsNonIso, &selMuonsNonIso, 
			      &ZeeCand, &ZmumuCand, 
			      &WeCand,  &WmuCand,
			      tightIso_e, tightIso_mu);
    
    

      
     
    if(ZeeCand.size()>0 ){
      if( IChannel == 12 ) MyhistoManager.FillHisto(CutFlow_eemu,   "CutFlow_eemu"  , 1, datasetName, IsSignal, Dweight[ITypeMC]);
      if( IChannel == 3 ) MyhistoManager.FillHisto(CutFlow_eee,    "CutFlow_eee"   , 1, datasetName, IsSignal, Dweight[ITypeMC]);
    }
    
    if(ZmumuCand.size()>0 ){
       if( IChannel == 0 ) MyhistoManager.FillHisto(CutFlow_mumumu, "CutFlow_mumumu", 1, datasetName, IsSignal, Dweight[ITypeMC]);
       if( IChannel == 12 ) MyhistoManager.FillHisto(CutFlow_mumue,  "CutFlow_mumue" , 1, datasetName, IsSignal, Dweight[ITypeMC]);
    }
    
      
    
    
    
    //*************************************************
    //              select 3 lepton events
    //*************************************************
    
    TString leptChannel = "";
    if( ZmumuCand.size() == 2 ) {
      if(WmuCand.size() == 1 && IChannel== 0  )  {	       leptChannel = "mumumu";}// for mumumu events
      if(WeCand.size()  == 1 && IChannel== 12 )  {IChannel= 1; leptChannel = "mumue"; }// for mumue events
    }
    
    if( ZeeCand.size() == 2 ) {
      if(WmuCand.size() == 1  && IChannel== 12  )  {IChannel= 2; leptChannel = "eemu";  }// for eemu events
      if(WeCand.size()  == 1  && IChannel== 3	)  {		leptChannel = "eee";   }// for eee events
    }
    
    if(ZmumuCand.size()>1){
      double invMassLL = (ZmumuCand[0].p4+ZmumuCand[1].p4).M();
      MyhistoManager.FillHisto(InvM_ll_mumumu_afterdileptsel, "InvM_ll_mumumu_afterdileptsel",  invMassLL, datasetName, IsSignal, Dweight[ITypeMC]);
    }
    
    
    //fixme
    //if( (WmuCand.size()+ZmumuCand.size()+WeCand.size()+ZeeCand.size()) == 3) {
   if( 
        (    (!useNonIsoWcand &&
	     ((selElectrons.size()  + selMuons.size() ) == 3 ) &&
	     ((looseElectrons.size()+looseMuons.size()) == 3 ) ) ||
	     
	     (useNonIsoWcand  &&
	     ( (WmuCand.size()+ZmumuCand.size()+WeCand.size()+ZeeCand.size()) == 3 ))
	      
	)&& 
	   ((IChannel== 0 && leptChannel == "mumumu") ||
	    (IChannel== 1 && leptChannel == "mumue") ||
	    (IChannel== 2 && leptChannel == "eemu") ||
	    (IChannel== 3 && leptChannel == "eee") )
    ) {
	
	
      
      //create pointers to histograms
      
      
      
     
      //*************************************************
      //              determine met 
      //*************************************************
    

       
      NTMET met;
      if(isData){
         met   = sel.GetSelectedMET(false, &*theJESuncertainty, upOrDown , false, 0);  
         smalltree_met_pt = met.p2.Mod();
         smalltree_met_phi = met.p2.Phi();

      }
      else{    
        NTMET met_unclsup, met_unclsdown;
        NTMET met_jesup, met_jesdown;
        NTMET met_jerup, met_jerdown;    
        NTMET metnosmear = sel.GetSelectedMET(false,    &*theJESuncertainty, upOrDown , false,    0) ;
 
    
      
        vector<NTJet>  looseJetsNoSmear = sel.GetSelectedJetsLoose(selMuons, selElectrons, false, theJESuncertainty, upOrDown, false,    0);     
        vector<NTJet>  looseJetsSmear   = sel.GetSelectedJetsLoose(selMuons, selElectrons, false, theJESuncertainty, upOrDown, true,     0);
      
        met=  SmearedMET(looseJetsSmear, looseJetsNoSmear, metnosmear);  
        
        smalltree_met_pt = met.p2.Mod();
        smalltree_met_phi = met.p2.Phi();
         
      
    
        vector<NTJet>  looseJetsNoSmear_JESup  = sel.GetSelectedJetsLoose(selMuons, selElectrons, false, theJESuncertainty, true, false,    -100);
        vector<NTJet>  looseJetsSmear_JESup    = sel.GetSelectedJetsLoose(selMuons, selElectrons, true, theJESuncertainty, true, true,     0);
        met_jesup=  SmearedMET(looseJetsSmear_JESup, looseJetsNoSmear_JESup, metnosmear);  
        smalltree_met_jesup_pt = met_jesup.p2.Mod();
        smalltree_met_jesup_phi = met_jesup.p2.Phi();
    
    
      
      
        vector<NTJet>  looseJetsNoSmear_JESdown  = sel.GetSelectedJetsLoose(selMuons, selElectrons, false, theJESuncertainty, true, false,    -100);
        vector<NTJet>  looseJetsSmear_JESdown    = sel.GetSelectedJetsLoose(selMuons, selElectrons, true, theJESuncertainty, false, true,     0);
        met_jesdown=  SmearedMET(looseJetsSmear_JESdown, looseJetsNoSmear_JESdown, metnosmear); 
        smalltree_met_jesdown_pt = met_jesdown.p2.Mod();
        smalltree_met_jesdown_phi = met_jesdown.p2.Phi();
    
    
      
      
        vector<NTJet>  looseJetsNoSmear_JERup  = sel.GetSelectedJetsLoose(selMuons, selElectrons, false, theJESuncertainty, true, false,    1);
        vector<NTJet>  looseJetsSmear_JERup    = sel.GetSelectedJetsLoose(selMuons, selElectrons, false, theJESuncertainty, true, true,     1);
        met_jerup=  SmearedMET(looseJetsSmear_JERup, looseJetsNoSmear_JERup, metnosmear); 
        smalltree_met_jerup_pt = met_jerup.p2.Mod();
        smalltree_met_jerup_phi = met_jerup.p2.Phi();
    
    
      
        vector<NTJet>  looseJetsNoSmear_JERdown  = sel.GetSelectedJetsLoose(selMuons, selElectrons, false, theJESuncertainty, true, false,    -1);
        vector<NTJet>  looseJetsSmear_JERdown    = sel.GetSelectedJetsLoose(selMuons, selElectrons, false, theJESuncertainty, true, true,     -1);
        met_jerdown=  SmearedMET(looseJetsSmear_JERdown, looseJetsNoSmear_JERdown, metnosmear); 
        smalltree_met_jerdown_pt = met_jerdown.p2.Mod();
        smalltree_met_jerdown_phi = met_jerdown.p2.Phi();
    
      
        
	NTMET met_uncl =  (sel.GetUnclusScaledMET(false, 1.1));
        met_unclsup = SmearedMET(looseJetsSmear, looseJetsNoSmear,  met_uncl);
        smalltree_met_unclsup_pt  = met_unclsup.p2.Mod();
        smalltree_met_unclsup_phi = met_unclsup.p2.Phi();
    
    
        met_uncl  = (sel.GetUnclusScaledMET(false, 0.9));
        met_unclsdown = SmearedMET(looseJetsSmear, looseJetsNoSmear, met_uncl);
        smalltree_met_unclsdown_pt  = met_unclsdown.p2.Mod();
        smalltree_met_unclsdown_phi = met_unclsdown.p2.Phi();
    
      } 
    
      
      
      /*if(leptChannel == "mumumu"){
        double pt  = selMuons[0].p4.Pt();
	double eta = selMuons[0].p4.Eta() ;
        cout <<  pt<< "  " << eta << endl;
        cout << "ID SF    " << sel.getScaleFactorMuID(pt, eta)     << endl;
        cout << "Iso12 SF " << sel.getScaleFactorMuIso12(pt, eta)  << endl;
        cout << "Iso20 SF " << sel.getScaleFactorMuIso20(pt, eta)  << endl;
      }
      cout << "start study global SF electons " << endl;
      if(leptChannel == "eee"){
        double pt  = selElectrons[0].p4.Pt();
	double eta = selElectrons[0].p4.Eta() ;
        cout <<  pt<< "  " << eta << endl;
	cout << "global SF " << sel.getSscaleFactorElectronAllID05(pt, eta)  << endl;
      
      }
      cout << "done " << endl;*/
      
      //fill the SmallTree
      
      smalltree_nlepton = 3;
      
      
      //*****************************************************************
      // define  Z/W candidate TLorentzVector
      //*****************************************************************    
      
      
    
      double lept3_Charge = 0;
      TLorentzVector dilept;
      TLorentzVector lept1, lept2, lept3;
      if(leptChannel == "mumue") {
	lept1 = ZmumuCand[0].p4;
	lept2 = ZmumuCand[1].p4;
	lept3 = WeCand[0].p4Gsf;
	lept3_Charge= WeCand[0].charge;
	dilept = lept1+lept2;
	
	smalltree_lept_pt[0]  = ZmumuCand[0].p4.Pt();
        smalltree_lept_eta[0] = ZmumuCand[0].p4.Eta();
        smalltree_lept_phi[0] = ZmumuCand[0].p4.Phi();
        smalltree_lept_iso[0] = sel.RelIso04PFDeltaBeta( ZmumuCand[0]  )  ;
        smalltree_lept_flav[0] =  (ZmumuCand[0].charge)*13 ;
	
	smalltree_lept_pt[1]  = ZmumuCand[1].p4.Pt();
        smalltree_lept_eta[1] = ZmumuCand[1].p4.Eta();
        smalltree_lept_phi[1] = ZmumuCand[1].p4.Phi();
        smalltree_lept_iso[1] = sel.RelIso04PFDeltaBeta( ZmumuCand[1] )  ;
        smalltree_lept_flav[1] =  (ZmumuCand[1].charge)*13 ;
	
	smalltree_lept_pt[2]  = WeCand[0].p4Gsf.Pt();
        smalltree_lept_eta[2] = WeCand[0].p4Gsf.Eta();
        smalltree_lept_phi[2] = WeCand[0].p4Gsf.Phi();
        smalltree_lept_iso[2] = sel.EffArea03PF(WeCand[0], rho)  ;
        smalltree_lept_flav[2] =  (WeCand[0].charge)*11 ;
	
	
	if(lept3.Pt() < 20) cout << "check pt " << lept3.Pt() << "  " <<  WeCand[0].p4.Pt() << endl;
	
      }else  if(leptChannel == "eemu") {
	lept1 = ZeeCand[0].p4Gsf;
	lept2 = ZeeCand[1].p4Gsf;
	lept3 = WmuCand[0].p4;
	lept3_Charge= WmuCand[0].charge;
	dilept = lept1+lept2;
	
	
	smalltree_lept_pt[0]  = ZeeCand[0].p4Gsf.Pt();
        smalltree_lept_eta[0] = ZeeCand[0].p4Gsf.Eta();
        smalltree_lept_phi[0] = ZeeCand[0].p4Gsf.Phi();
        smalltree_lept_iso[0] = sel.EffArea03PF(ZeeCand[0], rho)  ;
        smalltree_lept_flav[0] =  (ZeeCand[0].charge)*11 ;
	
	smalltree_lept_pt[1]  = ZeeCand[1].p4Gsf.Pt();
        smalltree_lept_eta[1] = ZeeCand[1].p4Gsf.Eta();
        smalltree_lept_phi[1] = ZeeCand[1].p4Gsf.Phi();
        smalltree_lept_iso[1] = sel.EffArea03PF(ZeeCand[1], rho)  ;
        smalltree_lept_flav[1] =  (ZeeCand[1].charge)*11 ;
	
	smalltree_lept_pt[2]  = WmuCand[0].p4.Pt();
        smalltree_lept_eta[2] = WmuCand[0].p4.Eta();
        smalltree_lept_phi[2] = WmuCand[0].p4.Phi();
        smalltree_lept_iso[2] = sel.RelIso04PFDeltaBeta( WmuCand[0] )  ;
        smalltree_lept_flav[2] =  (WmuCand[0].charge)*13 ;
	
	
      }else  if(leptChannel == "mumumu") {
	lept1 = ZmumuCand[0].p4;
	lept2 = ZmumuCand[1].p4;
	lept3 = WmuCand[0].p4;
	lept3_Charge= WmuCand[0].charge;
	dilept = lept1+lept2;
	
	smalltree_lept_pt[0]  = ZmumuCand[0].p4.Pt();
        smalltree_lept_eta[0] = ZmumuCand[0].p4.Eta();
        smalltree_lept_phi[0] = ZmumuCand[0].p4.Phi();
        smalltree_lept_iso[0] = sel.RelIso04PFDeltaBeta( ZmumuCand[0] )  ;
        smalltree_lept_flav[0] =  (ZmumuCand[0].charge)*13 ;
	
	smalltree_lept_pt[1]  = ZmumuCand[1].p4.Pt();
        smalltree_lept_eta[1] = ZmumuCand[1].p4.Eta();
        smalltree_lept_phi[1] = ZmumuCand[1].p4.Phi();
        smalltree_lept_iso[1] = sel.RelIso04PFDeltaBeta(ZmumuCand[1]  )  ;
        smalltree_lept_flav[1] =  (ZmumuCand[1].charge)*13 ;
	
	smalltree_lept_pt[2]  = WmuCand[0].p4.Pt();
        smalltree_lept_eta[2] = WmuCand[0].p4.Eta();
        smalltree_lept_phi[2] = WmuCand[0].p4.Phi();
        smalltree_lept_iso[2] = sel.RelIso04PFDeltaBeta( WmuCand[0] ) ;
        smalltree_lept_flav[2] =  (WmuCand[0].charge)*13 ;
	
	
	
      }else  if(leptChannel == "eee") {
	lept1 = ZeeCand[0].p4Gsf;
	lept2 = ZeeCand[1].p4Gsf;
	lept3 = WeCand[0].p4Gsf;
	lept3_Charge= WeCand[0].charge;
	dilept = lept1+lept2;
	
	
	
	smalltree_lept_pt[0]  = ZeeCand[0].p4Gsf.Pt();
        smalltree_lept_eta[0] = ZeeCand[0].p4Gsf.Eta();
        smalltree_lept_phi[0] = ZeeCand[0].p4Gsf.Phi();
        smalltree_lept_iso[0] = sel.EffArea03PF(ZeeCand[0], rho)  ;
        smalltree_lept_flav[0] =  (ZeeCand[0].charge)*11 ;
	
	smalltree_lept_pt[1]  = ZeeCand[1].p4Gsf.Pt();
        smalltree_lept_eta[1] = ZeeCand[1].p4Gsf.Eta();
        smalltree_lept_phi[1] = ZeeCand[1].p4Gsf.Phi();
        smalltree_lept_iso[1] = sel.EffArea03PF(ZeeCand[1], rho)  ;
        smalltree_lept_flav[1] =  (ZeeCand[1].charge)*11 ;
	
	smalltree_lept_pt[2]  = WeCand[0].p4Gsf.Pt();
        smalltree_lept_eta[2] = WeCand[0].p4Gsf.Eta();
        smalltree_lept_phi[2] = WeCand[0].p4Gsf.Phi();
        smalltree_lept_iso[2] = sel.EffArea03PF(WeCand[0], rho) ;
        smalltree_lept_flav[2] =  (WeCand[0].charge)*11 ;
	
      }
      
      
    
 
      smalltree_evtweight = Dweight[ITypeMC];
      //*****************************************************************
      // apply lepton scale factors
      //***************************************************************** 
      double leptonSF           = 1.;
      double leptonSF_scaleup   = 1.;
      double leptonSF_scaledown = 1.;
      
      if(applyLeptonSF && !isData){
        leptonSF	   = getLeptonSF(lept1, lept2, lept3,  leptChannel, 0)   ;
        leptonSF_scaleup   = getLeptonSF(lept1, lept2, lept3,  leptChannel, 1)   ;
        leptonSF_scaledown = getLeptonSF(lept1, lept2, lept3,  leptChannel, -1)   ;
        Dweight[ITypeMC]*=leptonSF;
      }
      
    
      
      
      
      
      
      int wcharge = lept3_Charge;
      
      //*****************************************************************
      // apply trigger scale factors
      //*****************************************************************  
            

      
      if(applyTrigger  &&  !isData ){ 
	if(IChannel == 0 ) Dweight[ITypeMC]*=SF_trig_mumumu;
	if(IChannel == 1 ) Dweight[ITypeMC]*=SF_trig_mumue;
	if(IChannel == 2 ) Dweight[ITypeMC]*=SF_trig_eemu;
	if(IChannel == 3 ) Dweight[ITypeMC]*=SF_trig_eee;
      }
      
            

      
      //*****************************************************************
      // apply DY scale factors
      //***************************************************************** 
      
      if(datasetName=="Zjets" && applyFakescale ){    
	if(IChannel == 0 ) Dweight[ITypeMC] = SF_Fake[0]*Dweight[ITypeMC];
	if(IChannel == 1 ) Dweight[ITypeMC] = SF_Fake[1]*Dweight[ITypeMC];
	if(IChannel == 2 ) Dweight[ITypeMC] = SF_Fake[2]*Dweight[ITypeMC];
	if(IChannel == 3 ) Dweight[ITypeMC] = SF_Fake[3]*Dweight[ITypeMC];
      }
      
      
            
 
      

      //*****************************************************************
      // apply WZ scale factors
      //***************************************************************** 
      
      
      if(datasetName=="WZ"  && applyWZ ){     
	if(IChannel == 0 ) Dweight[ITypeMC] = SF_WZ[0]*Dweight[ITypeMC];
	if(IChannel == 1 ) Dweight[ITypeMC] = SF_WZ[1]*Dweight[ITypeMC];
	if(IChannel == 2 ) Dweight[ITypeMC] = SF_WZ[2]*Dweight[ITypeMC];
	if(IChannel == 3 ) Dweight[ITypeMC] = SF_WZ[3]*Dweight[ITypeMC];
      }
      
      defineHistoPointer(IChannel);
      
     
      smalltree_evtweight        = Dweight[ITypeMC];
      smalltree_weight_leptup    = leptonSF_scaleup*Dweight[ITypeMC]/leptonSF;
      smalltree_weight_leptdown  = leptonSF_scaledown*Dweight[ITypeMC]/leptonSF;
	
	
      smalltree_evtweight	= Dweight[ITypeMC];
      smalltree_weight_trigup	= Dweight[ITypeMC];
      smalltree_weight_trigdown = Dweight[ITypeMC];
      
      
      if(pCutFlow == 0)     cout << "null pointer pCutFlow   " << endl;
      if(pErrCutFlow == 0)  cout << "null pointer pErrCutFlow" << endl;
      
      
      
      MyhistoManager.FillHisto(*pCutFlow,      ("CutFlow_"+leptChannel).Data(),    2, datasetName, IsSignal, Dweight[ITypeMC]);          
      MyhistoManager.FillHisto(*pErrCutFlow,   ("ErrCutFlow_"+leptChannel).Data(), 2, datasetName, IsSignal, EventYieldWeightError);
      
      
      
     

      double dileptonIvM = dilept.M();
      
      //***************************************************************** 
      //get selected jets 
      //*****************************************************************
      
      vector<NTElectron> emptyele; 
      vector<NTMuon>     emptymu; 
      
      
      
      
      vector<NTJet>  selJets = sel.GetSelectedJets(selMuons, selElectrons, false, theJESuncertainty, upOrDown, true, 0.);
      
      smalltree_njets = 0;
      for(unsigned int i=0 ; i< selJets.size(); i++){
      
        smalltree_jet_pt[i]   = selJets[i].p4.Pt();
        smalltree_jet_eta[i]  = selJets[i].p4.Eta();
        smalltree_jet_phi[i]  = selJets[i].p4.Phi();
        smalltree_jet_flav[i] = selJets[i].partonFlavour;
        smalltree_jet_btagdiscri[i]      = selJets[i].bTag["combinedSecondaryVertexBJetTags"];
	
        smalltree_njets++;
      } 
       

      double sumPtLeptJet = lept1.Pt() + lept2.Pt() + lept3.Pt();
      for(unsigned int i=0 ; i< selJets.size(); i++) sumPtLeptJet += selJets[i].p4.Pt();
      double theMET = met.p2.Mod();


      //***************************************************************** 
      //get b-tag info 
      //*****************************************************************
      
      /*if(leptChannel == "eemu"){
        cout << event->general.runNb << " " << event->general.eventNb << " " 
        << lept1.Pt()  << " " << lept1.Eta() << " "
        << lept2.Pt()  << " " << lept2.Eta() << " "
        << lept3.Pt()  << " " << lept3.Eta() << " ";
        for(unsigned int j=0; j< selJets.size(); j++){
          cout <<  selJets[j].p4.Pt() << " " <<  selJets[j].p4.Eta() << " ";
        }
        cout << met.p2.Mod() << endl;
      }*/

      
      int NBtaggedJets = 0;
      int idxBtag_tmp      = -10000;
      int idxBtag          = -10000;
      int AlgoBtag = sel.GetbtagAlgo();
      float btagDiscriCut = sel.GetbtagDiscriCut ();
      //cout << "-----------btagDiscriCut " << btagDiscriCut << endl;
      
      bool foundASelBjet = 0;
      for(unsigned int ijet = 0; ijet < selJets.size(); ijet++){
	if(abs(selJets[ijet].partonFlavour)==5 ){
	  foundASelBjet = true;
	}
	//cout << "  btagDiscriCut  " << btagDiscriCut << endl;
	//cout << "  bjet discri " << selJets[ijet].bTag["combinedSecondaryVertexBJetTags"] << endl;
	if ( selJets[ijet].bTag["combinedSecondaryVertexBJetTags"]      >= btagDiscriCut){
	  NBtaggedJets++;
	  //cout << "NBtaggedJets " << NBtaggedJets << endl;
	  if(idxBtag_tmp < 0) idxBtag_tmp = ijet;
	  else if(selJets[ijet].bTag["combinedSecondaryVertexBJetTags"] > selJets[idxBtag_tmp].bTag["combinedSecondaryVertexBJetTags"]) idxBtag_tmp = ijet; 
	}
      }  
      
      if(idxBtag_tmp > 0 ){ 
        
        for(unsigned int ijet = 0; ijet < selJets.size(); ijet++){
  	  if( fabs(selJets[idxBtag_tmp].p4.Pt() - selJets[ijet].p4.Pt() ) < 0.00001) idxBtag = idxBtag_tmp;
        }
      }
      
      //cout << "NBtaggedJets " << NBtaggedJets << endl;
      
      
      double btagdiscri = -10000;
      if(idxBtag>=0) btagdiscri = selJets[idxBtag].bTag["combinedSecondaryVertexBJetTags"];
      else if( selJets.size()>0)  btagdiscri = selJets[0].bTag["combinedSecondaryVertexBJetTags"];
      
            
      if(foundASelBjet) MyhistoManager.FillHisto(SelABjet, "SelABjet", 1.  , datasetName, IsSignal, 1.);
      else	        MyhistoManager.FillHisto(SelABjet, "SelABjet", 0.  , datasetName, IsSignal, 1. );
      
      
      if(idxBtag<0 && selJets.size() > 0) idxBtag = 0;
     
       
      
      int nvertex = selVertices.size();
            
         
     
      //smalltree_met_pt  = theMET;
      //smalltree_met_phi = met.p2.Phi(); 
      
      //-----------------------------------------------
      //for systematic jes up and down on jet selection
      //-----------------------------------------------
      
      vector<NTJet>  selJets_jesup   = sel.GetSelectedJets(selMuons, selElectrons, true, theJESuncertainty, true,  true, 0);  
      smalltree_jesup_njets = 0;
      for(unsigned int i=0 ; i< selJets_jesup.size(); i++){
        
        smalltree_jet_jesup_pt[i]   = selJets_jesup[i].p4.Pt();
        smalltree_jet_jesup_eta[i]  = selJets_jesup[i].p4.Eta();
        smalltree_jet_jesup_phi[i]  = selJets_jesup[i].p4.Phi();
        smalltree_jet_jesup_flav[i] = selJets_jesup[i].partonFlavour;
        smalltree_jet_jesup_btagdiscri[i]      = selJets_jesup[i].bTag["combinedSecondaryVertexBJetTags"];
	
        smalltree_jesup_njets++;
      } 
      
      
      
      vector<NTJet>  selJets_jesdown = sel.GetSelectedJets(selMuons, selElectrons, true, theJESuncertainty, false, true, 0);
      smalltree_jesdown_njets = 0;
      for(unsigned int i=0 ; i< selJets_jesdown.size(); i++){
      
        smalltree_jet_jesdown_pt[i]   = selJets_jesdown[i].p4.Pt();
        smalltree_jet_jesdown_eta[i]  = selJets_jesdown[i].p4.Eta();
        smalltree_jet_jesdown_phi[i]  = selJets_jesdown[i].p4.Phi();
        smalltree_jet_jesdown_flav[i] = selJets_jesdown[i].partonFlavour;
        smalltree_jet_jesdown_btagdiscri[i]      = selJets_jesdown[i].bTag["combinedSecondaryVertexBJetTags"];
	
        smalltree_jesdown_njets++;
      } 
      
      
      

      

      //-----------------------------------------------
      //for systematic jer up and down on jet selection
      //-----------------------------------------------
      vector<NTJet>  selJets_jerup   = sel.GetSelectedJets(selMuons, selElectrons, false, theJESuncertainty, true,  true, 1.);  
      smalltree_jerup_njets = 0;
      for(unsigned int i=0 ; i< selJets_jerup.size(); i++){
      
        smalltree_jet_jerup_pt[i]   = selJets_jerup[i].p4.Pt();
        smalltree_jet_jerup_eta[i]  = selJets_jerup[i].p4.Eta();
        smalltree_jet_jerup_phi[i]  = selJets_jerup[i].p4.Phi();
        smalltree_jet_jerup_flav[i] = selJets_jerup[i].partonFlavour;
        smalltree_jet_jerup_btagdiscri[i]      = selJets_jerup[i].bTag["combinedSecondaryVertexBJetTags"];
	
        smalltree_jerup_njets++;
      } 
      
      
      
      vector<NTJet>  selJets_jerdown = sel.GetSelectedJets(selMuons, selElectrons, false, theJESuncertainty, false, true, -1.);
      smalltree_jerdown_njets = 0;
      for(unsigned int i=0 ; i< selJets_jerdown.size(); i++){
      
        smalltree_jet_jerdown_pt[i]   = selJets_jerdown[i].p4.Pt();
        smalltree_jet_jerdown_eta[i]  = selJets_jerdown[i].p4.Eta();
        smalltree_jet_jerdown_phi[i]  = selJets_jerdown[i].p4.Phi();
        smalltree_jet_jerdown_flav[i] = selJets_jerdown[i].partonFlavour;
        smalltree_jet_jerdown_btagdiscri[i]      = selJets_jerdown[i].bTag["combinedSecondaryVertexBJetTags"];
	
        smalltree_jerdown_njets++;
      } 
      
      
      
      
      
      
      SmallTree->Fill();
      
      
      MyhistoManager.FillHisto(*pNvtx_afterleptsel,           ("Nvtx_"+leptChannel+"_afterleptsel").Data(),                nvertex, datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pInvM_ll_afterleptsel,        ("InvM_ll_"+leptChannel+"_afterleptsel").Data(),		dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pInvM_ll_afterleptsel_lowbin, ("InvM_ll_"+leptChannel+"_afterleptsel_lowbin").Data(),	dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pLeptZPt_afterleptsel,        ("LeptZPt_"+leptChannel+"_afterleptsel").Data(),		lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pLeptZPt_afterleptsel,        ("LeptZPt_"+leptChannel+"_afterleptsel").Data(),		lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pLeptWPt_afterleptsel,        ("LeptWPt_"+leptChannel+"_afterleptsel").Data(),		lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
      
      MyhistoManager.FillHisto(*pLeptPt_afterleptsel,        ("LeptPt_"+leptChannel+"_afterleptsel").Data(),		lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pLeptPt_afterleptsel,        ("LeptPt_"+leptChannel+"_afterleptsel").Data(),		lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pLeptPt_afterleptsel,        ("LeptPt_"+leptChannel+"_afterleptsel").Data(),		lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
   
      

      MyhistoManager.FillHisto(*pHT_afterleptsel,  ("HT_"+leptChannel+"_afterleptsel").Data(), sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pMET_afterleptsel, ("MET_"+leptChannel+"_afterleptsel").Data(), theMET, datasetName, IsSignal, Dweight[ITypeMC]);
      for(unsigned int i=0 ; i< selJets.size(); i++){ 
	MyhistoManager.FillHisto(*pJetPt_afterleptsel,          ("JetPt_"+leptChannel+"_afterleptsel").Data(),         selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	MyhistoManager.FillHisto(*pJetEta_afterleptsel,         ("JetEta_"+leptChannel+"_afterleptsel").Data(),        selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
	MyhistoManager.FillHisto2D(*pHT_vs_JetPt_afterleptsel  ,("HT_vs_JetPt_"+leptChannel+"_afterleptsel").Data()  , sumPtLeptJet, selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
      }
      
      
      MyhistoManager.FillHisto2D(*pHT_vs_MET_afterleptsel	  ,("HT_vs_MET_"+leptChannel+"_afterleptsel").Data()    , sumPtLeptJet, theMET, datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto2D(*pHT_vs_NJet_afterleptsel   ,("HT_vs_NJet_"+leptChannel+"_afterleptsel").Data()   , sumPtLeptJet, selJets.size(), datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto2D(*pHT_vs_NBJet_afterleptsel  ,("HT_vs_NBJet_"+leptChannel+"_afterleptsel").Data()  , sumPtLeptJet, NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto2D(*pHT_vs_LeptPt_afterleptsel ,("HT_vs_LeptPt_"+leptChannel+"_afterleptsel").Data() , sumPtLeptJet, lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto2D(*pHT_vs_LeptPt_afterleptsel ,("HT_vs_LeptPt_"+leptChannel+"_afterleptsel").Data() , sumPtLeptJet, lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto2D(*pHT_vs_LeptPt_afterleptsel ,("HT_vs_LeptPt_"+leptChannel+"_afterleptsel").Data() , sumPtLeptJet, lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
      if(selJets.size() == 1 && dilept.M() > 20) MyhistoManager.FillHisto2D(*pHT_vs_Mll_afterleptsel   ,("HT_vs_Mll_"+leptChannel+"_afterleptsel").Data()   , sumPtLeptJet, dilept.M(), datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pCharge_afterleptsel,    "Charge_"+leptChannel+"_afterleptsel"   , wcharge, datasetName, IsSignal, Dweight[ITypeMC]);
      
      


 
      
      //*****************************************************************
      // pass Z mass
      //***************************************************************** 
      
            
  
      
      //signal region
      if( fabs(dilept.M()-91) < 15){
       	
	
	MyhistoManager.FillHisto(*pCutFlow,     ("CutFlow_"+leptChannel).Data(),      3, datasetName, IsSignal, Dweight[ITypeMC]);
	MyhistoManager.FillHisto(*pErrCutFlow,	("ErrCutFlow_"+leptChannel).Data()  , 3, datasetName, IsSignal, EventYieldWeightError);
	//MyhistoManager.FillHisto(*pWmissAssing_afterZsel , ("WmissAssing_"+leptChannel+"afterZsel").Data(), leptonFlavor, datasetName, IsSignal, Dweight[ITypeMC]);
	
	
	
	double mTW = pow(
			 2*lept3.Pt()*met.p2.Mod()*(1-cos(lept3.Phi() - met.p2.Phi()))
			 ,0.5);
	
	
	
	MyhistoManager.FillHisto( *pmWT_afterZsel, ("mWT_"+leptChannel+"_afterZsel").Data(), mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	//MyhistoManager.FillHisto2D(*pInvM_ll_vs_mWT_afterZsel, ("InvM_ll_vs_mWT_"+leptChannel+"_afterZsel").Data(),dilept.M(), mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	
	
        //*****************************************************************
        // pass jet selection
        //***************************************************************** 
      	int NSeljets = selJets.size() ;
	if(NSeljets>4) NSeljets = 4;
	
	MyhistoManager.FillHisto(*pNBJet_afterZsel ,        ("NBJet_"+leptChannel+"_afterZsel").Data(), NBtaggedJets,      datasetName, IsSignal, Dweight[ITypeMC]); 
        
	MyhistoManager.FillHisto(*pNJet_afterZsel ,        ("NJet_"+leptChannel+"_afterZsel").Data(),NSeljets ,      datasetName, IsSignal, Dweight[ITypeMC]); 
        

        //cout << "before jet sel" <<  selJets.size() << endl;
	
  	if(selJets.size() >= 1  ){
	  
	  MyhistoManager.FillHisto(*pCutFlow,      ("CutFlow_"+leptChannel).Data(),    4, datasetName, IsSignal, Dweight[ITypeMC]);    
          MyhistoManager.FillHisto(*pErrCutFlow,   ("ErrCutFlow_"+leptChannel).Data(), 4, datasetName, IsSignal, EventYieldWeightError);
	  
	  MyhistoManager.FillHisto( *pmWT_afterjetsel, ("mWT_"+leptChannel+"_afterjetsel").Data(), mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	
	  MyhistoManager.FillHisto(*pInvM_ll_afterjetsel, ("InvM_ll_"+leptChannel+"_afterjetsel").Data(), dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
          
	 
 

	  
	  MyhistoManager.FillHisto(*pNBJet_afterjetsel ,  ("NBJet_"+leptChannel+"_afterjetsel").Data(),NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
          
	  
	  MyhistoManager.FillHisto(*pLeptZPt_afterjetsel,  ("LeptZPt_"+leptChannel+"_afterjetsel").Data(), lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(*pLeptZPt_afterjetsel,  ("LeptZPt_"+leptChannel+"_afterjetsel").Data(), lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(*pLeptWPt_afterjetsel,  ("LeptWPt_"+leptChannel+"_afterjetsel").Data(), lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  
	  MyhistoManager.FillHisto(*pHT_afterjetsel,      ("HT_"+leptChannel+"_afterjetsel").Data(), sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
 	  MyhistoManager.FillHisto(*pMET_afterjetsel,     ("MET_"+leptChannel+"_afterjetsel").Data(),theMET , datasetName, IsSignal, Dweight[ITypeMC]);
          for(unsigned int i=0 ; i< selJets.size(); i++){ 
	    MyhistoManager.FillHisto(*pJetPt_afterjetsel, ("JetPt_"+leptChannel+"_afterjetsel").Data(),
	    selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pJetEta_afterjetsel, ("JetEta_"+leptChannel+"_afterjetsel").Data(),
	    selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
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
	  
	  TLorentzVector theWcand = neutrino + lept3;

	  TLorentzVector topCand = neutrino + lept3 + selJets[idxBtag].p4 ;
	  
	

 
	  
	  if(NBtaggedJets == 0 ){ // bjet veto
	    MyhistoManager.FillHisto(*pHT_afterbjetveto,          ("HT_"+leptChannel+"_afterbjetveto").Data()        , sumPtLeptJet,  datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pRecoPtZ_afterbjetveto,     ("RecoPtZ_"+leptChannel+"_afterbjetveto").Data()   , Zpt ,          datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto( *pmWT_afterbjetveto,        ("mWT_"+leptChannel+"_afterbjetveto").Data()       , mTW,           datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pMt_afterbjetveto,          ("Mt_"+leptChannel+"_afterbjetveto").Data()        , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pdeltaPhilj_afterbjetveto , ("deltaPhilj_"+leptChannel+"_afterbjetveto").Data(),  deltaPhilj ,  datasetName, IsSignal, Dweight[ITypeMC]);
            MyhistoManager.FillHisto(*pLeptWPt_afterbjetveto,     ("LeptWPt_"+leptChannel+"_afterbjetveto").Data()   , lept3.Pt(),    datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pRecoTopMass_afterbjetveto, ("RecoTopMass_"+leptChannel+"_afterbjetveto").Data(), topCand.M() , datasetName, IsSignal, Dweight[ITypeMC]);
	    for(unsigned int i=0 ; i< selJets.size(); i++){
	      MyhistoManager.FillHisto(*pJetPt_afterbjetveto, ("JetPt_"+leptChannel+"_afterbjetveto").Data(),
	      selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(*pJetEta_afterbjetveto, ("JetEta_"+leptChannel+"_afterbjetveto").Data(),
	      selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
	    }
	  }
	  
	  //cout << "before b-tagging " << NBtaggedJets << endl;
	
	  if(
	    //(NBtaggedJets >=0  &&  NBtaggedJets <=1  && !useNonIsoWcand && theMET > themetcut) //||
	    //(NBtaggedJets >=0  &&  NBtaggedJets <=1  && !useNonIsoWcand ) ||
	    (NBtaggedJets <=1  && !useNonIsoWcand ) ||
	    //(NBtaggedJets ==1  && !useNonIsoWcand ) ||
	    //(NBtaggedJets >=0  &&  NBtaggedJets <=1  && useNonIsoWcand && theMET < themetcut)
	    ( NBtaggedJets <=1  &&  useNonIsoWcand )
	    //(NBtaggedJets ==1  && useNonIsoWcand && theMET < themetcut)
	    
	    
	    ){
	    
	   
 
	    
	    //cout << "after b-tagging " << endl;
            MyhistoManager.FillHisto(*pCutFlow,      ("CutFlow_"+leptChannel).Data(),	 5, datasetName, IsSignal, Dweight[ITypeMC]);	 
            MyhistoManager.FillHisto(*pErrCutFlow,   ("ErrCutFlow_"+leptChannel).Data(), 5, datasetName, IsSignal, EventYieldWeightError);
	    
	    MyhistoManager.FillHisto( *pmWT_afterbjetsel, ("mWT_"+leptChannel+"_afterbjetsel").Data(), mTW, datasetName, IsSignal, Dweight[ITypeMC]);

	    double nlept = selMuons.size()+selElectrons.size();
	    if(nlept>=4) nlept = 4;   
	      
		  	   
            MyhistoManager.FillHisto(*pNJet_afterbsel ,      ("NJet_"+leptChannel+"_afterbsel").Data(),        NSeljets, datasetName, IsSignal, Dweight[ITypeMC]);
  	    MyhistoManager.FillHisto(*pNLept_afterbsel ,     ("NLept_"+leptChannel+"_afterbsel").Data(),       nlept, datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pInvM_ll_afterbjetsel, ("InvM_ll_"+leptChannel+"_afterbjetsel").Data(),  dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pLeptZPt_afterbjetsel, ("LeptZPt_"+leptChannel+"_afterbjetsel").Data(),  lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pLeptWPt_afterbjetsel, ("LeptWPt_"+leptChannel+"_afterbjetsel").Data(),  lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pHT_afterbjetsel,      ("HT_"+leptChannel+"_afterbjetsel").Data(),       sumPtLeptJet, datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pMET_afterbjetsel,     ("MET_"+leptChannel+"_afterbjetsel").Data(),      theMET , datasetName, IsSignal, Dweight[ITypeMC]);
            for(unsigned int i=0 ; i< selJets.size(); i++){
	      MyhistoManager.FillHisto(*pJetPt_afterbjetsel, ("JetPt_"+leptChannel+"_afterbjetsel").Data(),
	      selJets[i].p4.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(*pJetEta_afterbjetsel, ("JetEta_"+leptChannel+"_afterbjetsel").Data(),
	      selJets[i].p4.Eta(), datasetName, IsSignal, Dweight[ITypeMC]);
	    }

	    
	    TLorentzVector metP4(met.p2.Px(), met.p2.Py(), 0, sqrt(met.p2.Px()*met.p2.Px() + met.p2.Py()*met.p2.Py()));
	    TLorentzVector transTop = lept3 + selJets[0].p4 + metP4;
	    
	    MyhistoManager.FillHisto(*pMt_afterbjetsel, ("Mt_"+leptChannel+"_afterbjetsel").Data() , transTop.Mt() ,datasetName, IsSignal, Dweight[ITypeMC]);
        
	    
 
	

	    // Use b jet for top mass computation
	    //topCand = neutrino + lept3 + selJets[idxBtag].p4 ;
	    
            
	    TLorentzVector theWcand_topRF = theWcand;
	    TLorentzVector letp_WRF = lept3;
	    
	    theWcand_topRF.Boost( TVector3(-1.0*topCand.BoostVector().Px(), -1.0*topCand.BoostVector().Py() ,-1.0*topCand.BoostVector().Pz()));
	    letp_WRF.Boost(       TVector3(-1.0*theWcand.BoostVector().Px(),-1.0*theWcand.BoostVector().Py(),-1.0*theWcand.BoostVector().Pz()));
	    
	    
	    double cosThetaStar = cos(letp_WRF.Vect().Angle(theWcand_topRF.Vect()));
            
	    tree_cosThetaStar = cosThetaStar;
	    tree_EvtWeight  = Dweight[ITypeMC];
	  


	  
            tree_topMass    = topCand.M();
            tree_totMass    = (topCand + (lept1+lept2)).M();
            tree_deltaPhilb = fabs(lept3.DeltaPhi(selJets[idxBtag].p4));
            tree_deltaRlb   = lept3.DeltaR(  selJets[idxBtag].p4);
            tree_deltaRTopZ = (lept1+lept2).DeltaR(topCand);
            tree_asym	  = lept3_Charge*fabs(lept3.Eta());
            tree_Zpt	  = (lept1+lept2).Pt();
            tree_ZEta	  = (lept1+lept2).Eta();
            tree_topPt	  = topCand.Pt();
            tree_topEta	  = topCand.Eta();
	  
	    //if(topCand.M() < 0) cout << "topCand.M() " << tree_topMass << endl;
	    //if(topCand.M() > 100000) cout << "topCand.M() " <<  tree_topMass<< endl;
	    
	    tree_totPt	    = (topCand + (lept1+lept2)).Pt();
	    tree_totEta	    = (topCand + (lept1+lept2)).Eta();

  	    tree_deltaRZl     = (lept1+lept2).DeltaR(lept3);
  	    tree_deltaPhiZmet = (lept1+lept2).DeltaPhi(metP4);
  	    tree_btagDiscri   = btagdiscri;
	    
 
  	    tree_NJets	    = float(selJets.size());
	  
	
 
	    tree_NBJets	    = float(NBtaggedJets);
	  
  	    tree_leptWPt        = lept3.Pt(); 
	    tree_leptWEta       = lept3.Eta();
  	    tree_leadJetPt      = selJets[idxBtag].p4.Pt(); 
	    tree_leadJetEta     = selJets[idxBtag].p4.Eta();
            tree_deltaRZleptW   = (lept1+lept2).DeltaR(lept3); 
	    tree_deltaPhiZleptW = (lept1+lept2).DeltaPhi(lept3);
	  
	    tree_met = theMET;
	    tree_mTW = mTW;
	  
	  
	  
 
	    MyhistoManager.FillHisto(*pAsym_afterbjetsel,	       ("Asym_"+leptChannel+"_afterbjetsel").Data(),         tree_asym,    datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pRecoPtZ_afterbjetsel,     ("RecoPtZ_"+leptChannel+"_afterbjetsel").Data(),     tree_Zpt,     datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pRecoTopMass_afterbjetsel, ("RecoTopMass_"+leptChannel+"_afterbjetsel").Data(), tree_topMass ,datasetName, IsSignal, Dweight[ITypeMC]);
            MyhistoManager.FillHisto(*pdeltaPhilb_afterbjetsel , ("deltaPhilb_"+leptChannel+"_afterbjetsel").Data(),  tree_deltaPhilb ,    datasetName, IsSignal, Dweight[ITypeMC]);

	    MyhistoManager.FillHisto(*pcosThetaStar_afterbjetsel, ("cosThetaStar_"+leptChannel+"_afterbjetsel").Data(), tree_cosThetaStar  ,    datasetName, IsSignal, Dweight[ITypeMC]);

	  
	  
	    if(datasetName=="DataMu" || datasetName=="DataEG" || datasetName=="DataMuEG")  tree_SampleType   = 0;
	    if(datasetName=="FCNCkut" )	     tree_SampleType   = 1;
	    if(datasetName=="TT" )	     tree_SampleType   = 2;
	    if(datasetName=="DYToLL_M10-50" )  tree_SampleType   = 3;
	    if(datasetName=="Zjets" )	     tree_SampleType   = 4;
	    if(datasetName=="Wjets" )	     tree_SampleType   = 5;
	    if(datasetName=="TtW" )	     tree_SampleType   = 6;
	    if(datasetName=="TbartW" )	     tree_SampleType   = 7;
	    if(datasetName=="WZ" )	     tree_SampleType   = 8;
	    if(datasetName=="ZZ" )	     tree_SampleType   = 9;
	    if(datasetName=="WW" )	     tree_SampleType   = 10;
	    if(datasetName=="FCNCkut" )	     tree_SampleType   = 11;
	    if(datasetName=="FCNCkct" )	     tree_SampleType   = 12;
	    if(datasetName=="FCNCxut" )	     tree_SampleType   = 13;
	    if(datasetName=="FCNCxct" )	     tree_SampleType   = 14;
	    if(datasetName=="FCNCzut" )	     tree_SampleType   = 15;
	    if(datasetName=="FCNCzct" )	     tree_SampleType   = 16;
	    if(datasetName=="tZq" )	     tree_SampleType   = 17;
	    if(leptChannel == "mumumu") tree_Channel = 0; 
	    if(leptChannel == "mumue" ) tree_Channel = 1; 
	    if(leptChannel == "eemu"  ) tree_Channel = 2; 
	    if(leptChannel == "eee"   ) tree_Channel = 3; 
	    
	    if(theMET > 20){
	      MyhistoManager.FillHisto(*pCutFlow,      ("CutFlow_"+leptChannel).Data(),	 6, datasetName, IsSignal, Dweight[ITypeMC]);	 
              MyhistoManager.FillHisto(*pErrCutFlow,   ("ErrCutFlow_"+leptChannel).Data(), 6, datasetName, IsSignal, EventYieldWeightError);
	     
	      /*if(leptChannel == "eee"){
	        cout << event->general.runNb << " " << event->general.eventNb << " " 
                << lept1.Pt()  << " " << lept1.Eta() << " "
                << lept2.Pt()  << " " << lept2.Eta() << " "
                << lept3.Pt()  << " " << lept3.Eta() << " ";
	        for(unsigned int j=0; j< selJets.size(); j++){
                  cout <<  selJets[j].p4.Pt() << " " <<  selJets[j].p4.Eta() << " ";
	        }
                cout << met.p2.Mod() << endl;
	      }*/
 
	    }
	  
	    //TheTree->Fill();
	  
	  } //bjet selection
	  

	} // at least on jet
	
      } // selection Z inv M
    
    }//3 lepton selection
    
  }//pass trigger selection
  
  
  
  
  
  
  
  return kTRUE;
}

//_____________________________________________________________________________
void ProofSelectorMyCutFlow::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
  if(fProofFile) fProofFile->Print();
  cout << "SlaveTerminate " << fFile << endl;
  if (fFile) {
    //Bool_t cleanup = kFALSE;
    //TDirectory *savedir = gDirectory;
    fFile->cd();
    
    
    
    
  /*MyhistoManager.WriteMyHisto(CutFlow_mumumu, "all" );
  MyhistoManager.WriteMyHisto(CutFlow_mumue,  "all" );
  MyhistoManager.WriteMyHisto(CutFlow_eemu,   "all" );
  MyhistoManager.WriteMyHisto(CutFlow_eee,    "all" );*/
  
    
    
    WriteTheHisto(&*fFile, &MyhistoManager);
    fFile->cd();
  
  //cout << "selected mumumu events " << nselevents_mumumu << endl;
    TheTree->Write();
    SmallTree->Write();
  
   //The following line is mandatory to copy everything in a common RootFile
    fOutput->Add(fProofFile);

    cleanHistoVector();
    delete TheTree, LumiWeights, anaEL;
    //delete fProofFile;
    fFile->Close("R");
    
  }
}





//_____________________________________________________________________________
void ProofSelectorMyCutFlow::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  cout << "stat terminate " << endl;
  //Possibility to retrieve information from the merged file and perform some calculation or plotting tasks
  delete event ;
  
  cout << "event deleted " << endl;
}



void ProofSelectorMyCutFlow::createTheHisto(HistoManager *thehistomanag){
  
  
  //create all histograms
  thehistomanag->CreateHisto(SelABjet, "SelABjet", datasetName, "Nevents","Entries", 2, -0.5, 1.5);
  thehistomanag->CreateHisto(SelABjet_afterjetsel, "SelABjet_afterjetsel", datasetName, "Nevents","Entries", 2, -0.5, 1.5);
 
  thehistomanag->CreateHisto(Ntrilept_mumumu, "Ntrilept_mumumu", datasetName, "Nevents","Entries", 11, -0.5, 10.5); 
  thehistomanag->CreateHisto(Ntrileptnoniso_mumumu, "Ntrileptnoniso_mumumu", datasetName, "Nevents","Entries", 11, -0.5, 10.5); 
  
  thehistomanag->CreateHisto(InvM_ll_mumumu_afterdileptsel, "InvM_ll_mumumu_afterdileptsel", datasetName, "DileptMll","Entries", 150, 60, 130);
  
  
  thehistomanag->CreateHisto(Nvertex, "Nvertex", datasetName, "Nvertex", "Entries", 50, 0, 50); 
  
  thehistomanag->CreateHisto(CutFlow_mumumu,  "CutFlow_mumumu" ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  thehistomanag->CreateHisto(CutFlow_mumue,   "CutFlow_mumue"  ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  thehistomanag->CreateHisto(CutFlow_eemu,    "CutFlow_eemu"   ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  thehistomanag->CreateHisto(CutFlow_eee,     "CutFlow_eee"    ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  
  
  thehistomanag->SetCutFlowAxisTitleFCNCMonotop(CutFlow_mumumu,   "CutFlow_mumumu"  ,datasetName);
  thehistomanag->SetCutFlowAxisTitleFCNCMonotop(CutFlow_mumue,    "CutFlow_mumue"   ,datasetName);
  thehistomanag->SetCutFlowAxisTitleFCNCMonotop(CutFlow_eemu,     "CutFlow_eemu"    ,datasetName);
  thehistomanag->SetCutFlowAxisTitleFCNCMonotop(CutFlow_eee,      "CutFlow_eee"     ,datasetName);
  
  
  thehistomanag->CreateHisto(NJetLight_mumu, "NJetLight_mumu", datasetName, "NJet", "Entries", 10, -0.5, 9.5);
  thehistomanag->CreateHisto(NJetLight_ee,   "NJetLight_ee",   datasetName, "NJet", "Entries", 10, -0.5, 9.5);
  
  thehistomanag->CreateHisto(NJetHeavy_mumu, "NJetHeavy_mumu", datasetName, "NJet", "Entries", 10, -0.5, 9.5);
  thehistomanag->CreateHisto(NJetHeavy_ee,   "NJetHeavy_ee"  , datasetName, "NJet", "Entries", 10, -0.5, 9.5);
  
  thehistomanag->CreateHisto(FlavComp_mumu, "FlavComp_mumu",   datasetName, "NJet", "Entries", 5,  -0.5, 4.5);
  thehistomanag->CreateHisto(FlavComp_ee,   "FlavComp_ee",     datasetName, "NJet", "Entries", 5,  -0.5, 4.5);
  
  
  thehistomanag->CreateHisto(PU_before_mumumu, "PU_before_mumumu", datasetName, "Npileup", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(PU_before_mumue,  "PU_before_mumue",  datasetName, "Npileup", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(PU_before_eemu,   "PU_before_eemu",   datasetName, "Npileup", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(PU_before_eee,    "PU_before_eee",    datasetName, "Npileup", "Entries", 60, 0, 60);   
  
  thehistomanag->CreateHisto(PU_intime_mumumu, "PU_intime_mumumu", datasetName, "Npileup", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(PU_intime_mumue,  "PU_intime_mumue",  datasetName, "Npileup", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(PU_intime_eemu,   "PU_intime_eemu",   datasetName, "Npileup", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(PU_intime_eee,    "PU_intime_eee",    datasetName, "Npileup", "Entries", 60, 0, 60);   
  
  thehistomanag->CreateHisto(PU_after_mumumu, "PU_after_mumumu", datasetName, "Npileup", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(PU_after_mumue,  "PU_after_mumue",  datasetName, "Npileup", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(PU_after_eemu,   "PU_after_eemu",   datasetName, "Npileup", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(PU_after_eee,    "PU_after_eee",    datasetName, "Npileup", "Entries", 60, 0, 60);   
  

  thehistomanag->CreateHisto(NVtx_mumumu, "NVtx_mumumu", datasetName, "Nvertex", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(NVtx_mumue,  "NVtx_mumue",  datasetName, "Nvertex", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(NVtx_eemu,   "NVtx_eemu",   datasetName, "Nvertex", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(NVtx_eee,    "NVtx_eee",    datasetName, "Nvertex", "Entries", 60, 0, 60); 
  
  thehistomanag->CreateHisto(Nvtx_mumumu_afterleptsel, "Nvtx_mumumu_afterleptsel", datasetName , "NVtx", "Entries", 60, 0, 60) ;
  thehistomanag->CreateHisto(Nvtx_mumue_afterleptsel,  "Nvtx_mumue_afterleptsel",  datasetName , "NVtx", "Entries", 60, 0, 60);
  thehistomanag->CreateHisto(Nvtx_eemu_afterleptsel,   "Nvtx_eemu_afterleptsel",   datasetName , "NVtx", "Entries", 60, 0, 60);
  thehistomanag->CreateHisto(Nvtx_eee_afterleptsel,    "Nvtx_eee_afterleptsel",    datasetName , "NVtx", "Entries", 60, 0, 60);
  
  thehistomanag->CreateHisto(NVtx_mumumu_aftertrigsel, "NVtx_mumumu_aftertrigsel", datasetName, "Nvertex", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(NVtx_mumue_aftertrigsel,  "NVtx_mumue_aftertrigsel",  datasetName, "Nvertex", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(NVtx_eemu_aftertrigsel,   "NVtx_eemu_aftertrigsel",   datasetName, "Nvertex", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(NVtx_eee_aftertrigsel,    "NVtx_eee_aftertrigsel",    datasetName, "Nvertex", "Entries", 60, 0, 60); 
  
  thehistomanag->CreateHisto(NVtx_mumumu_afterleptsel, "NVtx_mumumu_afterleptsel", datasetName, "Nvertex", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(NVtx_mumue_afterleptsel,  "NVtx_mumue_afterleptsel",  datasetName, "Nvertex", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(NVtx_eemu_afterleptsel,   "NVtx_eemu_afterleptsel",   datasetName, "Nvertex", "Entries", 60, 0, 60); 
  thehistomanag->CreateHisto(NVtx_eee_afterleptsel,    "NVtx_eee_afterleptsel",    datasetName, "Nvertex", "Entries", 60, 0, 60); 

  
  
  
  
  thehistomanag->CreateHisto(ErrCutFlow_mumumu,  "ErrCutFlow_mumumu"  ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  thehistomanag->CreateHisto(ErrCutFlow_mumue,   "ErrCutFlow_mumue"   ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  thehistomanag->CreateHisto(ErrCutFlow_eemu,    "ErrCutFlow_eemu"    ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  thehistomanag->CreateHisto(ErrCutFlow_eee,     "ErrCutFlow_eee"     ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  
  
  
  thehistomanag->CreateHisto(Mt_mumumu_afterbjetsel, "Mt_mumumu_afterbjetsel", datasetName,"Mt","Entries", 100, 0, 500); 
  thehistomanag->CreateHisto(Mt_mumue_afterbjetsel , "Mt_mumue_afterbjetsel" , datasetName,"Mt","Entries", 100, 0, 500);
  thehistomanag->CreateHisto(Mt_eemu_afterbjetsel  , "Mt_eemu_afterbjetsel"  , datasetName,"Mt","Entries", 100, 0, 500);
  thehistomanag->CreateHisto(Mt_eee_afterbjetsel   , "Mt_eee_afterbjetsel"   , datasetName,"Mt","Entries", 100, 0, 500);
   
  thehistomanag->CreateHisto(Mt_mumumu_afterbjetveto, "Mt_mumumu_afterbjetveto", datasetName,"Mt","Entries", 100, 0, 500); 
  thehistomanag->CreateHisto(Mt_mumue_afterbjetveto , "Mt_mumue_afterbjetveto" , datasetName,"Mt","Entries", 100, 0, 500);
  thehistomanag->CreateHisto(Mt_eemu_afterbjetveto  , "Mt_eemu_afterbjetveto"  , datasetName,"Mt","Entries", 100, 0, 500);
  thehistomanag->CreateHisto(Mt_eee_afterbjetveto   , "Mt_eee_afterbjetveto"   , datasetName,"Mt","Entries", 100, 0, 500);
   
  
  thehistomanag->CreateHisto(NJet_mumumu_afterZsel, "NJet_mumumu_afterZsel", datasetName,"Njets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NJet_mumue_afterZsel , "NJet_mumue_afterZsel" , datasetName,"Njets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NJet_eemu_afterZsel  , "NJet_eemu_afterZsel"  , datasetName,"Njets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NJet_eee_afterZsel   , "NJet_eee_afterZsel"   , datasetName,"Njets", "Entries",5,-0.5,4.5);
  
  thehistomanag->CreateHisto(NJet_mumumu_afterbsel, "NJet_mumumu_afterbsel", datasetName,"Njets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NJet_mumue_afterbsel , "NJet_mumue_afterbsel" , datasetName,"Njets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NJet_eemu_afterbsel  , "NJet_eemu_afterbsel"  , datasetName,"Njets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NJet_eee_afterbsel   , "NJet_eee_afterbsel"   , datasetName,"Njets", "Entries",5,-0.5,4.5);
  
  thehistomanag->CreateHisto(NJet_mumumu_afterleptsel_mWT110, "NJet_mumumu_afterleptsel_mWT110", datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NJet_mumue_afterleptsel_mWT110,  "NJet_mumue_afterleptsel_mWT110",  datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NJet_eemu_afterleptsel_mWT110,   "NJet_eemu_afterleptsel_mWT110",   datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NJet_eee_afterleptsel_mWT110,    "NJet_eee_afterleptsel_mWT110",	 datasetName,"NBjets", "Entries",5,-0.5,4.5);
  
  
  
  
  
  
  thehistomanag->CreateHisto(NLept_mumumu_afterbsel, "NLept_mumumu_afterbsel", datasetName,"NLepts", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NLept_mumue_afterbsel , "NLept_mumue_afterbsel" , datasetName,"NLepts", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NLept_eemu_afterbsel  , "NLept_eemu_afterbsel"  , datasetName,"NLepts", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NLept_eee_afterbsel   , "NLept_eee_afterbsel"   , datasetName,"NLepts", "Entries",5,-0.5,4.5);
  
  
  
  thehistomanag->CreateHisto(NBJet_mumumu_afterZsel, "NBJet_mumumu_afterZsel", datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_mumue_afterZsel , "NBJet_mumue_afterZsel" , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_eemu_afterZsel  , "NBJet_eemu_afterZsel"  , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_eee_afterZsel   , "NBJet_eee_afterZsel"   , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  
   
  thehistomanag->CreateHisto(mWT_mumumu_afterZsel, "mWT_mumumu_afterZsel", datasetName,"mWT", "Entries",100,0,200);
  thehistomanag->CreateHisto(mWT_mumue_afterZsel , "mWT_mumue_afterZsel" , datasetName,"mWT", "Entries",100,0,200);
  thehistomanag->CreateHisto(mWT_eemu_afterZsel  , "mWT_eemu_afterZsel"  , datasetName,"mWT", "Entries",100,0,200);
  thehistomanag->CreateHisto(mWT_eee_afterZsel   , "mWT_eee_afterZsel"   , datasetName,"mWT", "Entries",100,0,200);
 
  
  
  
  thehistomanag->CreateHisto(NBJet_mumumu_afterjetsel, "NBJet_mumumu_afterjetsel", datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_mumue_afterjetsel , "NBJet_mumue_afterjetsel" , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_eemu_afterjetsel  , "NBJet_eemu_afterjetsel"  , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_eee_afterjetsel   , "NBJet_eee_afterjetsel"   , datasetName,"NBjets", "Entries",5,-0.5,4.5);	
  
  
  
  
  
  thehistomanag->CreateHisto(NBJet_mumumu_afterjetsel_bjets, "NBJet_mumumu_afterjetsel_bjets", datasetName,"NBjets", "Entries",2,-0.5,1.5);
  thehistomanag->CreateHisto(NBJet_mumue_afterjetsel_bjets , "NBJet_mumue_afterjetsel_bjets" , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  thehistomanag->CreateHisto(NBJet_eemu_afterjetsel_bjets  , "NBJet_eemu_afterjetsel_bjets"  , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  thehistomanag->CreateHisto(NBJet_eee_afterjetsel_bjets   , "NBJet_eee_afterjetsel_bjets"   , datasetName,"NBjets", "Entries",2,-0.5,1.5);	
  
  
  thehistomanag->CreateHisto(NBJet_mumumu_afterjetsel_cjets, "NBJet_mumumu_afterjetsel_cjets", datasetName,"NBjets", "Entries",2,-0.5,1.5);
  thehistomanag->CreateHisto(NBJet_mumue_afterjetsel_cjets , "NBJet_mumue_afterjetsel_cjets" , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  thehistomanag->CreateHisto(NBJet_eemu_afterjetsel_cjets  , "NBJet_eemu_afterjetsel_cjets"  , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  thehistomanag->CreateHisto(NBJet_eee_afterjetsel_cjets   , "NBJet_eee_afterjetsel_cjets"   , datasetName,"NBjets", "Entries",2,-0.5,1.5);	
  
  
  thehistomanag->CreateHisto(NBJet_mumumu_afterjetsel_ljets, "NBJet_mumumu_afterjetsel_ljets", datasetName,"NBjets", "Entries",2,-0.5,1.5);
  thehistomanag->CreateHisto(NBJet_mumue_afterjetsel_ljets , "NBJet_mumue_afterjetsel_ljets" , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  thehistomanag->CreateHisto(NBJet_eemu_afterjetsel_ljets  , "NBJet_eemu_afterjetsel_ljets"  , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  thehistomanag->CreateHisto(NBJet_eee_afterjetsel_ljets   , "NBJet_eee_afterjetsel_ljets"   , datasetName,"NBjets", "Entries",2,-0.5,1.5);	
  
  
  
  
  thehistomanag->CreateHisto(BJetDiscri_mumumu_afterjetsel_bjets, "BJetDiscri_mumumu_afterjetsel_bjets", datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  thehistomanag->CreateHisto(BJetDiscri_mumue_afterjetsel_bjets , "BJetDiscri_mumue_afterjetsel_bjets" , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  thehistomanag->CreateHisto(BJetDiscri_eemu_afterjetsel_bjets  , "BJetDiscri_eemu_afterjetsel_bjets"  , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  thehistomanag->CreateHisto(BJetDiscri_eee_afterjetsel_bjets   , "BJetDiscri_eee_afterjetsel_bjets"   , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  
  
  
  thehistomanag->CreateHisto(BJetDiscri_mumumu_afterjetsel_cjets, "BJetDiscri_mumumu_afterjetsel_cjets", datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  thehistomanag->CreateHisto(BJetDiscri_mumue_afterjetsel_cjets , "BJetDiscri_mumue_afterjetsel_cjets" , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  thehistomanag->CreateHisto(BJetDiscri_eemu_afterjetsel_cjets  , "BJetDiscri_eemu_afterjetsel_cjets"  , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  thehistomanag->CreateHisto(BJetDiscri_eee_afterjetsel_cjets   , "BJetDiscri_eee_afterjetsel_cjets"   , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
 
  
  thehistomanag->CreateHisto(BJetDiscri_mumumu_afterjetsel_ljets, "BJetDiscri_mumumu_afterjetsel_ljets", datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  thehistomanag->CreateHisto(BJetDiscri_mumue_afterjetsel_ljets , "BJetDiscri_mumue_afterjetsel_ljets" , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  thehistomanag->CreateHisto(BJetDiscri_eemu_afterjetsel_ljets  , "BJetDiscri_eemu_afterjetsel_ljets"  , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
  thehistomanag->CreateHisto(BJetDiscri_eee_afterjetsel_ljets   , "BJetDiscri_eee_afterjetsel_ljets"   , datasetName,"BjetsDiscri", "Entries",100, 0, 1.);
 
  
  
  
  thehistomanag->CreateHisto(NBJet_mumumu_afterleptsel_mWT110, "NBJet_mumumu_afterleptsel_mWT110", datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_mumue_afterleptsel_mWT110,  "NBJet_mumue_afterleptsel_mWT110",  datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_eemu_afterleptsel_mWT110,   "NBJet_eemu_afterleptsel_mWT110",   datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_eee_afterleptsel_mWT110,    "NBJet_eee_afterleptsel_mWT110",    datasetName,"NBjets", "Entries",5,-0.5,4.5);
  
  
  thehistomanag->CreateHisto(NJet_mumumu_afterleptsel, "NJet_mumumu_afterleptsel", datasetName,"Njets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NJet_mumue_afterleptsel,  "NJet_mumue_afterleptsel",  datasetName,"Njets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NJet_eemu_afterleptsel,   "NJet_eemu_afterleptsel",   datasetName,"Njets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NJet_eee_afterleptsel,    "NJet_eee_afterleptsel",    datasetName,"Njets", "Entries",5,-0.5,4.5);
  
  
  
  
  
  /*thehistomanag->CreateHisto(NBJet_mumumu_afterjetsel, "NBJet_mumumu_afterjetsel", datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_mumue_afterjetsel , "NBJet_mumue_afterjetsel" , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_eemu_afterjetsel  , "NBJet_eemu_afterjetsel"  , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_eee_afterjetsel   , "NBJet_eee_afterjetsel"   , datasetName,"NBjets", "Entries",5,-0.5,4.5); 
  */
   
  
  thehistomanag->CreateHisto(InvM_ll_mumumu_afterleptsel_highSumPt, "InvM_ll_mumumu_afterleptsel_highSumPt" , datasetName,"Minv", "Entries",100,0.,1000);                     
  thehistomanag->CreateHisto(InvM_ll_mumue_afterleptsel_highSumPt,  "InvM_ll_mumue_afterleptsel_highSumPt"  , datasetName,"Minv", "Entries",100,0.,1000);
  thehistomanag->CreateHisto(InvM_ll_eemu_afterleptsel_highSumPt,   "InvM_ll_eemu_afterleptsel_highSumPt"   , datasetName,"Minv", "Entries",100,0.,1000);
  thehistomanag->CreateHisto(InvM_ll_eee_afterleptsel_highSumPt,    "InvM_ll_eee_afterleptsel_highSumPt"    , datasetName,"Minv", "Entries",100,0.,1000);

 
  thehistomanag->CreateHisto(InvM_ll_mumumu_afterleptsel, "InvM_ll_mumumu_afterleptsel" , datasetName,"Minv", "Entries",100,0.,250);
  thehistomanag->CreateHisto(InvM_ll_mumue_afterleptsel,  "InvM_ll_mumue_afterleptsel"  , datasetName,"Minv", "Entries",100,0.,250);
  thehistomanag->CreateHisto(InvM_ll_eemu_afterleptsel,   "InvM_ll_eemu_afterleptsel"   , datasetName,"Minv", "Entries",100,0.,250);
  thehistomanag->CreateHisto(InvM_ll_eee_afterleptsel,    "InvM_ll_eee_afterleptsel"    , datasetName,"Minv", "Entries",100,0.,250);
  
  thehistomanag->CreateHisto(InvM_ll_mumumu_afterleptsel_mWT110, "InvM_ll_mumumu_afterleptsel_mWT110" , datasetName,"Minv", "Entries",100,0.,250);
  thehistomanag->CreateHisto(InvM_ll_mumue_afterleptsel_mWT110,  "InvM_ll_mumue_afterleptsel_mWT110"  , datasetName,"Minv", "Entries",100,0.,250);
  thehistomanag->CreateHisto(InvM_ll_eemu_afterleptsel_mWT110,   "InvM_ll_eemu_afterleptsel_mWT110"   , datasetName,"Minv", "Entries",100,0.,250);
  thehistomanag->CreateHisto(InvM_ll_eee_afterleptsel_mWT110,    "InvM_ll_eee_afterleptsel_mWT110"    , datasetName,"Minv", "Entries",100,0.,250);
 
  thehistomanag->CreateHisto(InvM_ll_mumumu_afterleptsel_lowbin, "InvM_ll_mumumu_afterleptsel_lowbin" , datasetName,"Minv", "Entries",100,0.,200);
  thehistomanag->CreateHisto(InvM_ll_mumue_afterleptsel_lowbin,  "InvM_ll_mumue_afterleptsel_lowbin"  , datasetName,"Minv", "Entries",100,0.,200);
  thehistomanag->CreateHisto(InvM_ll_eemu_afterleptsel_lowbin,   "InvM_ll_eemu_afterleptsel_lowbin"   , datasetName,"Minv", "Entries",100,0.,200);
  thehistomanag->CreateHisto(InvM_ll_eee_afterleptsel_lowbin,    "InvM_ll_eee_afterleptsel_lowbin"    , datasetName,"Minv", "Entries",100,0.,200);

  thehistomanag->CreateHisto(InvM_ll_mumumu_afterjetsel, "InvM_ll_mumumu_afterjetsel" , datasetName,"Minv", "Entries",100,0.,200);
  thehistomanag->CreateHisto(InvM_ll_mumue_afterjetsel,  "InvM_ll_mumue_afterjetsel"  , datasetName,"Minv", "Entries",100,0.,200);
  thehistomanag->CreateHisto(InvM_ll_eemu_afterjetsel,   "InvM_ll_eemu_afterjetsel"   , datasetName,"Minv", "Entries",100,0.,200);
  thehistomanag->CreateHisto(InvM_ll_eee_afterjetsel,    "InvM_ll_eee_afterjetsel"    , datasetName,"Minv", "Entries",100,0.,200);
  
  thehistomanag->CreateHisto(InvM_ll_mumumu_afterbjetsel, "InvM_ll_mumumu_afterbjetsel" , datasetName,"Minv", "Entries",100,0.,250);
  thehistomanag->CreateHisto(InvM_ll_mumue_afterbjetsel,  "InvM_ll_mumue_afterbjetsel"  , datasetName,"Minv", "Entries",100,0.,250);
  thehistomanag->CreateHisto(InvM_ll_eemu_afterbjetsel,   "InvM_ll_eemu_afterbjetsel"   , datasetName,"Minv", "Entries",100,0.,250);
  thehistomanag->CreateHisto(InvM_ll_eee_afterbjetsel,    "InvM_ll_eee_afterbjetsel"    , datasetName,"Minv", "Entries",100,0.,250);

    
  
  thehistomanag->CreateHisto(LeptPt_mumumu_afterleptsel, "LeptPt_mumumu_afterleptsel", datasetName, "LeptPt", "Entries",100,0.,200) ;
  thehistomanag->CreateHisto(LeptPt_mumue_afterleptsel,  "LeptPt_mumue_afterleptsel",  datasetName, "LeptPt", "Entries",100,0.,200);
  thehistomanag->CreateHisto(LeptPt_eemu_afterleptsel,   "LeptPt_eemu_afterleptsel",   datasetName, "LeptPt", "Entries",100,0.,200);
  thehistomanag->CreateHisto(LeptPt_eee_afterleptsel,    "LeptPt_eee_afterleptsel",    datasetName, "LeptPt", "Entries",100,0.,200);
  
  thehistomanag->CreateHisto(LeptPt_mumumu_afterjetsel, "LeptPt_mumumu_afterjetsel", datasetName, "LeptPt", "Entries",100,0.,200);
  thehistomanag->CreateHisto(LeptPt_mumue_afterjetsel,  "LeptPt_mumue_afterjetsel",  datasetName, "LeptPt", "Entries",100,0.,200);
  thehistomanag->CreateHisto(LeptPt_eemu_afterjetsel,   "LeptPt_eemu_afterjetsel",   datasetName, "LeptPt", "Entries",100,0.,200);
  thehistomanag->CreateHisto(LeptPt_eee_afterjetsel,    "LeptPt_eee_afterjetsel",    datasetName, "LeptPt", "Entries",100,0.,200);
  
  thehistomanag->CreateHisto(LeptPt_mumumu_afterbjetsel, "LeptPt_mumumu_afterbjetsel", datasetName, "LeptPt", "Entries",100,0.,200);
  thehistomanag->CreateHisto(LeptPt_mumue_afterbjetsel,  "LeptPt_mumue_afterbjetsel",  datasetName, "LeptPt", "Entries",100,0.,200);
  thehistomanag->CreateHisto(LeptPt_eemu_afterbjetsel,   "LeptPt_eemu_afterbjetsel",   datasetName, "LeptPt", "Entries",100,0.,200);
  thehistomanag->CreateHisto(LeptPt_eee_afterbjetsel,    "LeptPt_eee_afterbjetsel",    datasetName, "LeptPt", "Entries",100,0.,200);
    
  thehistomanag->CreateHisto(LeptZPt_mumumu_afterleptsel, "LeptZPt_mumumu_afterleptsel", datasetName, "LeptZPt", "Entries",350,0., 1000) ;
  thehistomanag->CreateHisto(LeptZPt_mumue_afterleptsel,  "LeptZPt_mumue_afterleptsel",  datasetName, "LeptZPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptZPt_eemu_afterleptsel,   "LeptZPt_eemu_afterleptsel",   datasetName, "LeptZPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptZPt_eee_afterleptsel,    "LeptZPt_eee_afterleptsel",    datasetName, "LeptZPt", "Entries",350,0., 1000);
  
  thehistomanag->CreateHisto(LeptZPt_mumumu_afterjetsel, "LeptZPt_mumumu_afterjetsel", datasetName, "LeptZPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptZPt_mumue_afterjetsel,  "LeptZPt_mumue_afterjetsel",  datasetName, "LeptZPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptZPt_eemu_afterjetsel,   "LeptZPt_eemu_afterjetsel",   datasetName, "LeptZPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptZPt_eee_afterjetsel,    "LeptZPt_eee_afterjetsel",    datasetName, "LeptZPt", "Entries",350,0., 1000);
  
  thehistomanag->CreateHisto(LeptZPt_mumumu_afterbjetsel, "LeptZPt_mumumu_afterbjetsel", datasetName, "LeptZPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptZPt_mumue_afterbjetsel,  "LeptZPt_mumue_afterbjetsel",  datasetName, "LeptZPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptZPt_eemu_afterbjetsel,   "LeptZPt_eemu_afterbjetsel",   datasetName, "LeptZPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptZPt_eee_afterbjetsel,    "LeptZPt_eee_afterbjetsel",    datasetName, "LeptZPt", "Entries",350,0., 1000);
    
  thehistomanag->CreateHisto(LeptWPt_mumumu_afterleptsel, "LeptWPt_mumumu_afterleptsel", datasetName, "LeptWPt", "Entries",350,0., 1000) ;
  thehistomanag->CreateHisto(LeptWPt_mumue_afterleptsel,  "LeptWPt_mumue_afterleptsel",  datasetName, "LeptWPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptWPt_eemu_afterleptsel,   "LeptWPt_eemu_afterleptsel",   datasetName, "LeptWPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptWPt_eee_afterleptsel,    "LeptWPt_eee_afterleptsel",    datasetName, "LeptWPt", "Entries",350,0., 1000);
  
  thehistomanag->CreateHisto(LeptWPt_mumumu_afterjetsel, "LeptWPt_mumumu_afterjetsel", datasetName, "LeptWPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptWPt_mumue_afterjetsel,  "LeptWPt_mumue_afterjetsel",  datasetName, "LeptWPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptWPt_eemu_afterjetsel,   "LeptWPt_eemu_afterjetsel",   datasetName, "LeptWPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptWPt_eee_afterjetsel,    "LeptWPt_eee_afterjetsel",    datasetName, "LeptWPt", "Entries",350,0., 1000);
  
  thehistomanag->CreateHisto(LeptWPt_mumumu_afterbjetsel, "LeptWPt_mumumu_afterbjetsel", datasetName, "LeptWPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptWPt_mumue_afterbjetsel,  "LeptWPt_mumue_afterbjetsel",  datasetName, "LeptWPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptWPt_eemu_afterbjetsel,   "LeptWPt_eemu_afterbjetsel",   datasetName, "LeptWPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptWPt_eee_afterbjetsel,    "LeptWPt_eee_afterbjetsel",    datasetName, "LeptWPt", "Entries",350,0., 1000);

  thehistomanag->CreateHisto(LeptWPt_mumumu_afterbjetveto, "LeptWPt_mumumu_afterbjetveto", datasetName, "LeptWPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptWPt_mumue_afterbjetveto,  "LeptWPt_mumue_afterbjetveto",  datasetName, "LeptWPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptWPt_eemu_afterbjetveto,   "LeptWPt_eemu_afterbjetveto",   datasetName, "LeptWPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptWPt_eee_afterbjetveto,    "LeptWPt_eee_afterbjetveto",    datasetName, "LeptWPt", "Entries",350,0., 1000);

      
  thehistomanag->CreateHisto(LeptWPt_mumumu_afterleptsel_mWT110, "LeptWPt_mumumu_afterleptsel_mWT110", datasetName, "LeptWPt", "Entries",350,0., 1000) ;
  thehistomanag->CreateHisto(LeptWPt_mumue_afterleptsel_mWT110,  "LeptWPt_mumue_afterleptsel_mWT110",  datasetName, "LeptWPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptWPt_eemu_afterleptsel_mWT110,   "LeptWPt_eemu_afterleptsel_mWT110",   datasetName, "LeptWPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptWPt_eee_afterleptsel_mWT110,    "LeptWPt_eee_afterleptsel_mWT110",    datasetName, "LeptWPt", "Entries",350,0., 1000);
  

  
  
  thehistomanag->CreateHisto(JetPt_mumumu_afterleptsel, "JetPt_mumumu_afterleptsel", datasetName, "JetPt", "Entries",100,0., 300) ;
  thehistomanag->CreateHisto(JetPt_mumue_afterleptsel,  "JetPt_mumue_afterleptsel",  datasetName, "JetPt", "Entries",100,0., 300);
  thehistomanag->CreateHisto(JetPt_eemu_afterleptsel,   "JetPt_eemu_afterleptsel",   datasetName, "JetPt", "Entries",100,0., 300);
  thehistomanag->CreateHisto(JetPt_eee_afterleptsel,    "JetPt_eee_afterleptsel",    datasetName, "JetPt", "Entries",100,0., 300);
  
  thehistomanag->CreateHisto(JetPt_mumumu_afterjetsel, "JetPt_mumumu_afterjetsel", datasetName, "JetPt", "Entries",100,0., 300);
  thehistomanag->CreateHisto(JetPt_mumue_afterjetsel,  "JetPt_mumue_afterjetsel",  datasetName, "JetPt", "Entries",100,0., 300);
  thehistomanag->CreateHisto(JetPt_eemu_afterjetsel,   "JetPt_eemu_afterjetsel",   datasetName, "JetPt", "Entries",100,0., 300);
  thehistomanag->CreateHisto(JetPt_eee_afterjetsel,    "JetPt_eee_afterjetsel",    datasetName, "JetPt", "Entries",100,0., 300);
  
  thehistomanag->CreateHisto(JetPt_mumumu_afterbjetsel, "JetPt_mumumu_afterbjetsel", datasetName, "JetPt", "Entries",100,0., 300);
  thehistomanag->CreateHisto(JetPt_mumue_afterbjetsel,  "JetPt_mumue_afterbjetsel",  datasetName, "JetPt", "Entries",100,0., 300);
  thehistomanag->CreateHisto(JetPt_eemu_afterbjetsel,   "JetPt_eemu_afterbjetsel",   datasetName, "JetPt", "Entries",100,0., 300);
  thehistomanag->CreateHisto(JetPt_eee_afterbjetsel,    "JetPt_eee_afterbjetsel",    datasetName, "JetPt", "Entries",100,0., 300);
  
  thehistomanag->CreateHisto(JetPt_mumumu_afterbjetveto, "JetPt_mumumu_afterbjetveto", datasetName, "JetPt", "Entries",100,0., 300);
  thehistomanag->CreateHisto(JetPt_mumue_afterbjetveto,  "JetPt_mumue_afterbjetveto",  datasetName, "JetPt", "Entries",100,0., 300);
  thehistomanag->CreateHisto(JetPt_eemu_afterbjetveto,   "JetPt_eemu_afterbjetveto",   datasetName, "JetPt", "Entries",100,0., 300);
  thehistomanag->CreateHisto(JetPt_eee_afterbjetveto,    "JetPt_eee_afterbjetveto",    datasetName, "JetPt", "Entries",100,0., 300);
  
    
  thehistomanag->CreateHisto(JetEta_mumumu_afterleptsel, "JetEta_mumumu_afterleptsel", datasetName, "JetEta", "Entries",26, -2.5, 2.5) ;
  thehistomanag->CreateHisto(JetEta_mumue_afterleptsel,  "JetEta_mumue_afterleptsel",  datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  thehistomanag->CreateHisto(JetEta_eemu_afterleptsel,   "JetEta_eemu_afterleptsel",   datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  thehistomanag->CreateHisto(JetEta_eee_afterleptsel,    "JetEta_eee_afterleptsel",    datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  
  thehistomanag->CreateHisto(JetEta_mumumu_afterjetsel, "JetEta_mumumu_afterjetsel", datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  thehistomanag->CreateHisto(JetEta_mumue_afterjetsel,  "JetEta_mumue_afterjetsel",  datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  thehistomanag->CreateHisto(JetEta_eemu_afterjetsel,   "JetEta_eemu_afterjetsel",   datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  thehistomanag->CreateHisto(JetEta_eee_afterjetsel,    "JetEta_eee_afterjetsel",    datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  
  thehistomanag->CreateHisto(JetEta_mumumu_afterbjetsel, "JetEta_mumumu_afterbjetsel", datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  thehistomanag->CreateHisto(JetEta_mumue_afterbjetsel,  "JetEta_mumue_afterbjetsel",  datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  thehistomanag->CreateHisto(JetEta_eemu_afterbjetsel,   "JetEta_eemu_afterbjetsel",   datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  thehistomanag->CreateHisto(JetEta_eee_afterbjetsel,    "JetEta_eee_afterbjetsel",    datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  
  thehistomanag->CreateHisto(JetEta_mumumu_afterbjetveto, "JetEta_mumumu_afterbjetveto", datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  thehistomanag->CreateHisto(JetEta_mumue_afterbjetveto,  "JetEta_mumue_afterbjetveto",  datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  thehistomanag->CreateHisto(JetEta_eemu_afterbjetveto,   "JetEta_eemu_afterbjetveto",   datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  thehistomanag->CreateHisto(JetEta_eee_afterbjetveto,    "JetEta_eee_afterbjetveto",    datasetName, "JetEta", "Entries",26, -2.5, 2.5 );
  
 
  
  thehistomanag->CreateHisto(HT_mumumu_afterleptsel, "HT_mumumu_afterleptsel", datasetName, "HT", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(HT_mumue_afterleptsel,  "HT_mumue_afterleptsel",  datasetName, "HT", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(HT_eemu_afterleptsel,   "HT_eemu_afterleptsel",   datasetName, "HT", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(HT_eee_afterleptsel,    "HT_eee_afterleptsel",    datasetName, "HT", "Entries",350,0., 1000);
  
  
  thehistomanag->CreateHisto(HT_mumumu_afterjetsel, "HT_mumumu_afterjetsel", datasetName, "HT", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(HT_mumue_afterjetsel,  "HT_mumue_afterjetsel",  datasetName, "HT", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(HT_eemu_afterjetsel,   "HT_eemu_afterjetsel",   datasetName, "HT", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(HT_eee_afterjetsel,    "HT_eee_afterjetsel",    datasetName, "HT", "Entries",350,0., 1000);
  
  thehistomanag->CreateHisto(HT_mumumu_afterbjetsel, "HT_mumumu_afterbjetsel", datasetName, "HT", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(HT_mumue_afterbjetsel,  "HT_mumue_afterbjetsel",  datasetName, "HT", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(HT_eemu_afterbjetsel,   "HT_eemu_afterbjetsel",   datasetName, "HT", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(HT_eee_afterbjetsel,    "HT_eee_afterbjetsel",    datasetName, "HT", "Entries",350,0., 1000);
  
  thehistomanag->CreateHisto(HT_mumumu_afterbjetveto, "HT_mumumu_afterbjetveto", datasetName, "HT", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(HT_mumue_afterbjetveto,  "HT_mumue_afterbjetveto",  datasetName, "HT", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(HT_eemu_afterbjetveto,   "HT_eemu_afterbjetveto",   datasetName, "HT", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(HT_eee_afterbjetveto,    "HT_eee_afterbjetveto",    datasetName, "HT", "Entries",350,0., 1000);
  
  
  
  thehistomanag->CreateHisto(MET_mumumu_afterleptsel, "MET_mumumu_afterleptsel", datasetName, "MET", "Entries",50,0., 250);
  thehistomanag->CreateHisto(MET_mumue_afterleptsel,  "MET_mumue_afterleptsel",  datasetName, "MET", "Entries",50,0., 250);
  thehistomanag->CreateHisto(MET_eemu_afterleptsel,   "MET_eemu_afterleptsel",   datasetName, "MET", "Entries",50,0., 250);
  thehistomanag->CreateHisto(MET_eee_afterleptsel,    "MET_eee_afterleptsel",    datasetName, "MET", "Entries",50,0., 250);
  
  thehistomanag->CreateHisto(MET_mumumu_afterleptsel_mWT110, "MET_mumumu_afterleptsel_mWT110", datasetName, "MET", "Entries",50,0., 250);
  thehistomanag->CreateHisto(MET_mumue_afterleptsel_mWT110,  "MET_mumue_afterleptsel_mWT110",  datasetName, "MET", "Entries",50,0., 250);
  thehistomanag->CreateHisto(MET_eemu_afterleptsel_mWT110,   "MET_eemu_afterleptsel_mWT110",   datasetName, "MET", "Entries",50,0., 250);
  thehistomanag->CreateHisto(MET_eee_afterleptsel_mWT110,    "MET_eee_afterleptsel_mWT110",    datasetName, "MET", "Entries",50,0., 250);
  
  
  
  
  thehistomanag->CreateHisto(MET_mumumu_afterjetsel, "MET_mumumu_afterjetsel", datasetName, "MET", "Entries",50,0., 250);
  thehistomanag->CreateHisto(MET_mumue_afterjetsel,  "MET_mumue_afterjetsel",  datasetName, "MET", "Entries",50,0., 250);
  thehistomanag->CreateHisto(MET_eemu_afterjetsel,   "MET_eemu_afterjetsel",   datasetName, "MET", "Entries",50,0., 250);
  thehistomanag->CreateHisto(MET_eee_afterjetsel,    "MET_eee_afterjetsel",    datasetName, "MET", "Entries",50,0., 250);
  
  thehistomanag->CreateHisto(MET_mumumu_afterbjetsel, "MET_mumumu_afterbjetsel", datasetName, "MET", "Entries",50,0., 250);
  thehistomanag->CreateHisto(MET_mumue_afterbjetsel,  "MET_mumue_afterbjetsel",  datasetName, "MET", "Entries",50,0., 250);
  thehistomanag->CreateHisto(MET_eemu_afterbjetsel,   "MET_eemu_afterbjetsel",   datasetName, "MET", "Entries",50,0., 250);
  thehistomanag->CreateHisto(MET_eee_afterbjetsel,    "MET_eee_afterbjetsel",    datasetName, "MET", "Entries",50,0., 250);
    
   
  thehistomanag->CreateHisto(Asym_mumumu_afterbjetsel, "Asym_mumumu_afterbjetsel", datasetName, "Asym", "Entries", 20, -3.2, 3.2);
  thehistomanag->CreateHisto(Asym_mumue_afterbjetsel , "Asym_mumue_afterbjetsel" , datasetName, "Asym", "Entries", 20, -3.2, 3.2);
  thehistomanag->CreateHisto(Asym_eemu_afterbjetsel  , "Asym_eemu_afterbjetsel"  , datasetName, "Asym", "Entries", 20, -3.2, 3.2);
  thehistomanag->CreateHisto(Asym_eee_afterbjetsel   , "Asym_eee_afterbjetsel"   , datasetName, "Asym", "Entries", 20, -3.2, 3.2);
  
  thehistomanag->CreateHisto(RecoPtZ_mumumu_afterbjetsel, "RecoPtZ_mumumu_afterbjetsel", datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto(RecoPtZ_mumue_afterbjetsel , "RecoPtZ_mumue_afterbjetsel" , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto(RecoPtZ_eemu_afterbjetsel  , "RecoPtZ_eemu_afterbjetsel"  , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto(RecoPtZ_eee_afterbjetsel   , "RecoPtZ_eee_afterbjetsel"   , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  
  thehistomanag->CreateHisto(RecoPtZ_mumumu_afterbjetveto, "RecoPtZ_mumumu_afterbjetveto", datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto(RecoPtZ_mumue_afterbjetveto , "RecoPtZ_mumue_afterbjetveto" , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto(RecoPtZ_eemu_afterbjetveto  , "RecoPtZ_eemu_afterbjetveto"  , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto(RecoPtZ_eee_afterbjetveto   , "RecoPtZ_eee_afterbjetveto"   , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
   
  thehistomanag->CreateHisto(RecoPtZ_mumumu_afterleptsel, "RecoPtZ_mumumu_afterleptsel", datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto(RecoPtZ_mumue_afterleptsel , "RecoPtZ_mumue_afterleptsel" , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto(RecoPtZ_eemu_afterleptsel  , "RecoPtZ_eemu_afterleptsel"  , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto(RecoPtZ_eee_afterleptsel   , "RecoPtZ_eee_afterleptsel"   , datasetName, "RecoPtZ", "Entries", 200, 0, 500);
   
  thehistomanag->CreateHisto(RecoPtZ_mumumu_afterleptsel_nojet, "RecoPtZ_mumumu_afterleptsel_nojet", datasetName, "RecoPtZ", "Entries", 200, 0, 100);
  thehistomanag->CreateHisto(RecoPtZ_mumue_afterleptsel_nojet , "RecoPtZ_mumue_afterleptsel_nojet" , datasetName, "RecoPtZ", "Entries", 200, 0, 100);
  thehistomanag->CreateHisto(RecoPtZ_eemu_afterleptsel_nojet  , "RecoPtZ_eemu_afterleptsel_nojet"  , datasetName, "RecoPtZ", "Entries", 200, 0, 100);
  thehistomanag->CreateHisto(RecoPtZ_eee_afterleptsel_nojet   , "RecoPtZ_eee_afterleptsel_nojet"   , datasetName, "RecoPtZ", "Entries", 200, 0, 100);

  thehistomanag->CreateHisto( RecoTopMass_mumumu_afterbjetsel, "RecoTopMass_mumumu_afterbjetsel", datasetName, "TopMass", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto( RecoTopMass_mumue_afterbjetsel,  "RecoTopMass_mumue_afterbjetsel" , datasetName, "TopMass", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto( RecoTopMass_eemu_afterbjetsel,   "RecoTopMass_eemu_afterbjetsel"  , datasetName, "TopMass", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto( RecoTopMass_eee_afterbjetsel,    "RecoTopMass_eee_afterbjetsel"   , datasetName, "TopMass", "Entries", 200, 0, 500);
  
  thehistomanag->CreateHisto( RecoTopMass_mumumu_afterbjetveto, "RecoTopMass_mumumu_afterbjetveto", datasetName, "TopMass", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto( RecoTopMass_mumue_afterbjetveto,  "RecoTopMass_mumue_afterbjetveto" , datasetName, "TopMass", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto( RecoTopMass_eemu_afterbjetveto,   "RecoTopMass_eemu_afterbjetveto"  , datasetName, "TopMass", "Entries", 200, 0, 500);
  thehistomanag->CreateHisto( RecoTopMass_eee_afterbjetveto,    "RecoTopMass_eee_afterbjetveto"   , datasetName, "TopMass", "Entries", 200, 0, 500);
  
  
  thehistomanag->CreateHisto(deltaPhilb_mumumu_afterbjetsel ,"deltaPhilb_mumumu_afterbjetsel", datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  thehistomanag->CreateHisto(deltaPhilb_mumue_afterbjetsel  ,"deltaPhilb_mumue_afterbjetsel",  datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  thehistomanag->CreateHisto(deltaPhilb_eemu_afterbjetsel   ,"deltaPhilb_eemu_afterbjetsel",   datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  thehistomanag->CreateHisto(deltaPhilb_eee_afterbjetsel    ,"deltaPhilb_eee_afterbjetsel",    datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  
  thehistomanag->CreateHisto(deltaPhilj_mumumu_afterbjetveto ,"deltaPhilj_mumumu_afterbjetveto", datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  thehistomanag->CreateHisto(deltaPhilj_mumue_afterbjetveto  ,"deltaPhilj_mumue_afterbjetveto",  datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  thehistomanag->CreateHisto(deltaPhilj_eemu_afterbjetveto   ,"deltaPhilj_eemu_afterbjetveto",   datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  thehistomanag->CreateHisto(deltaPhilj_eee_afterbjetveto    ,"deltaPhilj_eee_afterbjetveto",    datasetName, "deltaPhilb", "Entries", 200, 0, 3.2);
  
  
  
  
  thehistomanag->CreateHisto(deltaR_mumumu_afterleptsel, "deltaR_mumumu_afterleptsel", datasetName, "deltaRLept", "Entries", 20, 0, 3.2);
  thehistomanag->CreateHisto(deltaR_mumue_afterleptsel,  "deltaR_mumue_afterleptsel",  datasetName, "deltaRLept", "Entries", 20, 0, 3.2);
  thehistomanag->CreateHisto(deltaR_eemu_afterleptsel,   "deltaR_eemu_afterleptsel",   datasetName, "deltaRLept", "Entries", 20, 0, 3.2);
  thehistomanag->CreateHisto(deltaR_eee_afterleptsel,    "deltaR_eee_afterleptsel",    datasetName, "deltaRLept", "Entries", 20, 0, 3.2);
  
  
  
  thehistomanag->CreateHisto(deltaRLeptJet_mumumu_afterleptsel_mWT110,"deltaRLeptJet_mumumu_afterleptsel_mWT110", datasetName, "deltaR", "Entries", 20, 0, 3.2);
  thehistomanag->CreateHisto(deltaRLeptJet_mumue_afterleptsel_mWT110 ,"deltaRLeptJet_mumue_afterleptsel_mWT110",  datasetName, "deltaR", "Entries", 20, 0, 3.2);
  thehistomanag->CreateHisto(deltaRLeptJet_eemu_afterleptsel_mWT110  ,"deltaRLeptJet_eemu_afterleptsel_mWT110",   datasetName, "deltaR", "Entries", 20, 0, 3.2);
  thehistomanag->CreateHisto(deltaRLeptJet_eee_afterleptsel_mWT110   ,"deltaRLeptJet_eee_afterleptsel_mWT110",    datasetName, "deltaR", "Entries", 20, 0, 3.2);
  
  thehistomanag->CreateHisto(deltaRLeptMet_mumumu_afterleptsel_mWT110,"deltaRLeptMet_mumumu_afterleptsel_mWT110", datasetName, "deltaR", "Entries", 20, 0, 3.2);
  thehistomanag->CreateHisto(deltaRLeptMet_mumue_afterleptsel_mWT110 ,"deltaRLeptMet_mumue_afterleptsel_mWT110",  datasetName, "deltaR", "Entries", 20, 0, 3.2);
  thehistomanag->CreateHisto(deltaRLeptMet_eemu_afterleptsel_mWT110  ,"deltaRLeptMet_eemu_afterleptsel_mWT110",   datasetName, "deltaR", "Entries", 20, 0, 3.2);
  thehistomanag->CreateHisto(deltaRLeptMet_eee_afterleptsel_mWT110   ,"deltaRLeptMet_eee_afterleptsel_mWT110",    datasetName, "deltaR", "Entries", 20, 0, 3.2);
  
  
  
  thehistomanag->CreateHisto(WmissAssing_mumumu_afterleptsel, "WmissAssing_mumumu_afterleptsel", datasetName, "MissAs", "Entries", 3, -0.5, 1.5);
  thehistomanag->CreateHisto(WmissAssing_mumue_afterleptsel,  "WmissAssing_mumue_afterleptsel",  datasetName, "MissAs", "Entries", 3, -0.5, 1.5);
  thehistomanag->CreateHisto(WmissAssing_eemu_afterleptsel,   "WmissAssing_eemu_afterleptsel",   datasetName, "MissAs", "Entries", 3, -0.5, 1.5);
  thehistomanag->CreateHisto(WmissAssing_eee_afterleptsel,    "WmissAssing_eee_afterleptsel",    datasetName, "MissAs", "Entries", 3, -0.5, 1.5);
  
  
  
  thehistomanag->CreateHisto(mWT_mumumu_afterbjetveto, "mWT_mumumu_afterbjetveto", datasetName, "MWT", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(mWT_mumue_afterbjetveto,  "mWT_mumue_afterbjetveto" , datasetName, "MWT", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(mWT_eemu_afterbjetveto,   "mWT_eemu_afterbjetveto"  , datasetName, "MWT", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(mWT_eee_afterbjetveto,    "mWT_eee_afterbjetveto"   , datasetName, "MWT", "Entries", 100, 0, 200);
  
  thehistomanag->CreateHisto(mWT_mumumu_afterbjetsel, "mWT_mumumu_afterbjetsel", datasetName, "MWT", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(mWT_mumue_afterbjetsel,  "mWT_mumue_afterbjetsel" , datasetName, "MWT", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(mWT_eemu_afterbjetsel,   "mWT_eemu_afterbjetsel"  , datasetName, "MWT", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(mWT_eee_afterbjetsel,    "mWT_eee_afterbjetsel"   , datasetName, "MWT", "Entries", 100, 0, 200);
  
  thehistomanag->CreateHisto(mWT_mumumu_afterjetsel, "mWT_mumumu_afterjetsel", datasetName, "MWT", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(mWT_mumue_afterjetsel,  "mWT_mumue_afterjetsel" , datasetName, "MWT", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(mWT_eemu_afterjetsel,   "mWT_eemu_afterjetsel"  , datasetName, "MWT", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(mWT_eee_afterjetsel,    "mWT_eee_afterjetsel"   , datasetName, "MWT", "Entries", 100, 0, 200);
  
  thehistomanag->CreateHisto(mWT_mumumu_afterleptsel, "mWT_mumumu_afterleptsel", datasetName, "MWT", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(mWT_mumue_afterleptsel,  "mWT_mumue_afterleptsel" , datasetName, "MWT", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(mWT_eemu_afterleptsel,   "mWT_eemu_afterleptsel"  , datasetName, "MWT", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(mWT_eee_afterleptsel,    "mWT_eee_afterleptsel"   , datasetName, "MWT", "Entries", 100, 0, 200);
  
  thehistomanag->CreateHisto(cosThetaStar_mumumu_afterbjetsel, "cosThetaStar_mumumu_afterbjetsel", datasetName, "cosThetaStar", "Entries", 20, -1, 1);
  thehistomanag->CreateHisto(cosThetaStar_mumue_afterbjetsel,  "cosThetaStar_mumue_afterbjetsel",  datasetName, "cosThetaStar", "Entries", 20, -1, 1);
  thehistomanag->CreateHisto(cosThetaStar_eemu_afterbjetsel,   "cosThetaStar_eemu_afterbjetsel",   datasetName, "cosThetaStar", "Entries", 20, -1, 1);
  thehistomanag->CreateHisto(cosThetaStar_eee_afterbjetsel,    "cosThetaStar_eee_afterbjetsel",    datasetName, "cosThetaStar", "Entries", 20, -1, 1);

  
  thehistomanag->CreateHisto(Charge_mumumu_afterleptsel, "Charge_mumumu_afterleptsel", datasetName, "Charge", "Entries", 11, -5, 5);
  thehistomanag->CreateHisto(Charge_mumue_afterleptsel,  "Charge_mumue_afterleptsel",  datasetName, "Charge", "Entries", 11, -5, 5);
  thehistomanag->CreateHisto(Charge_eemu_afterleptsel,   "Charge_eemu_afterleptsel",   datasetName, "Charge", "Entries", 11, -5, 5);
  thehistomanag->CreateHisto(Charge_eee_afterleptsel,    "Charge_eee_afterleptsel",    datasetName, "Charge", "Entries", 11, -5, 5);
  
  thehistomanag->CreateHisto(Charge_mumumu_afterleptsel_mWT110, "Charge_mumumu_afterleptsel_mWT110", datasetName, "Charge", "Entries", 11, -5, 5);
  thehistomanag->CreateHisto(Charge_mumue_afterleptsel_mWT110,  "Charge_mumue_afterleptsel_mWT110",  datasetName, "Charge", "Entries", 11, -5, 5);
  thehistomanag->CreateHisto(Charge_eemu_afterleptsel_mWT110,   "Charge_eemu_afterleptsel_mWT110",   datasetName, "Charge", "Entries", 11, -5, 5);
  thehistomanag->CreateHisto(Charge_eee_afterleptsel_mWT110,    "Charge_eee_afterleptsel_mWT110",    datasetName, "Charge", "Entries", 11, -5, 5);
  
  
  
  
  thehistomanag->CreateHisto(DijetInvM_mumumu_afterleptsel_inZpeak, "DijetInvM_mumumu_afterleptsel_inZpeak", datasetName, "DiJet", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(DijetInvM_mumue_afterleptsel_inZpeak,  "DijetInvM_mumue_afterleptsel_inZpeak" , datasetName, "DiJet", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(DijetInvM_eemu_afterleptsel_inZpeak,   "DijetInvM_eemu_afterleptsel_inZpeak"  , datasetName, "DiJet", "Entries", 100, 0, 200);
  thehistomanag->CreateHisto(DijetInvM_eee_afterleptsel_inZpeak,    "DijetInvM_eee_afterleptsel_inZpeak"   , datasetName, "DiJet", "Entries", 100, 0, 200);
  
  
  
  
  //**********************
  //initiate 2D histograms
  //**********************
  
  thehistomanag->CreateHisto2D(HT_vs_MET_mumumu_afterleptsel, "HT_vs_MET_mumumu_afterleptsel", datasetName, "HT",100,0., 1000., "MET",  100,0., 1000.)  ;
  thehistomanag->CreateHisto2D(HT_vs_MET_mumue_afterleptsel , "HT_vs_MET_mumue_afterleptsel", datasetName, "HT", 100,0., 1000., "MET",  100,0., 1000.);
  thehistomanag->CreateHisto2D(HT_vs_MET_eemu_afterleptsel  , "HT_vs_MET_eemu_afterleptsel",  datasetName, "HT", 100,0., 1000., "MET",  100,0., 1000.);
  thehistomanag->CreateHisto2D(HT_vs_MET_eee_afterleptsel   , "HT_vs_MET_eee_afterleptsel",   datasetName, "HT", 100,0., 1000., "MET",  100,0., 1000.);
  
  thehistomanag->CreateHisto2D(HT_vs_NJet_mumumu_afterleptsel, "HT_vs_NJet_mumumu_afterleptsel", datasetName, "HT", 100,0., 1000,"Njets", 5,-0.5,4.5);
  thehistomanag->CreateHisto2D(HT_vs_NJet_mumue_afterleptsel , "HT_vs_NJet_mumue_afterleptsel",  datasetName, "HT", 100,0., 1000,"Njets", 5,-0.5,4.5);
  thehistomanag->CreateHisto2D(HT_vs_NJet_eemu_afterleptsel  , "HT_vs_NJet_eemu_afterleptsel",   datasetName, "HT", 100,0., 1000,"Njets", 5,-0.5,4.5);
  thehistomanag->CreateHisto2D(HT_vs_NJet_eee_afterleptsel   , "HT_vs_NJet_eee_afterleptsel",    datasetName, "HT", 100,0., 1000,"Njets", 5,-0.5,4.5);
  
  thehistomanag->CreateHisto2D(HT_vs_NBJet_mumumu_afterleptsel, "HT_vs_NBJet_mumumu_afterleptsel", datasetName, "HT", 100,0., 1000,"NBjets", 5,-0.5,4.5);
  thehistomanag->CreateHisto2D(HT_vs_NBJet_mumue_afterleptsel , "HT_vs_NBJet_mumue_afterleptsel",  datasetName, "HT", 100,0., 1000,"NBjets", 5,-0.5,4.5);
  thehistomanag->CreateHisto2D(HT_vs_NBJet_eemu_afterleptsel  , "HT_vs_NBJet_eemu_afterleptsel",   datasetName, "HT", 100,0., 1000,"NBjets", 5,-0.5,4.5);
  thehistomanag->CreateHisto2D(HT_vs_NBJet_eee_afterleptsel   , "HT_vs_NBJet_eee_afterleptsel",    datasetName, "HT", 100,0., 1000,"NBjets", 5,-0.5,4.5);
  
  thehistomanag->CreateHisto2D(HT_vs_LeptPt_mumumu_afterleptsel, "HT_vs_LeptPt_mumumu_afterleptsel", datasetName, "HT", 100,0., 1000, "LeptPt",100,0., 1000.);
  thehistomanag->CreateHisto2D(HT_vs_LeptPt_mumue_afterleptsel , "HT_vs_LeptPt_mumue_afterleptsel",  datasetName, "HT", 100,0., 1000, "LeptPt",100,0., 1000.);
  thehistomanag->CreateHisto2D(HT_vs_LeptPt_eemu_afterleptsel  , "HT_vs_LeptPt_eemu_afterleptsel",   datasetName, "HT", 100,0., 1000, "LeptPt",100,0., 1000.);
  thehistomanag->CreateHisto2D(HT_vs_LeptPt_eee_afterleptsel   , "HT_vs_LeptPt_eee_afterleptsel",    datasetName, "HT", 100,0., 1000, "LeptPt",100,0., 1000.);
  
  thehistomanag->CreateHisto2D(HT_vs_JetPt_mumumu_afterleptsel, "HT_vs_JetPt_mumumu_afterleptsel", datasetName, "HT", 100,0., 1000, "JetPt",100,0., 1000.);
  thehistomanag->CreateHisto2D(HT_vs_JetPt_mumue_afterleptsel  , "HT_vs_JetPt_mumue_afterleptsel",   datasetName, "HT", 100,0., 1000, "JetPt",100,0., 1000.);
  thehistomanag->CreateHisto2D(HT_vs_JetPt_eemu_afterleptsel   , "HT_vs_JetPt_eemu_afterleptsel",    datasetName, "HT", 100,0., 1000, "JetPt",100,0., 1000.);
  thehistomanag->CreateHisto2D(HT_vs_JetPt_eee_afterleptsel    , "HT_vs_JetPt_eee_afterleptsel",     datasetName, "HT", 100,0., 1000, "JetPt",100,0., 1000.);
  
  
  thehistomanag->CreateHisto2D(HT_vs_Mll_mumumu_afterleptsel, "HT_vs_Mll_mumumu_afterleptsel",   datasetName, "HT", 100,0., 1000, "Mll",100,0., 1000.);
  thehistomanag->CreateHisto2D(HT_vs_Mll_mumue_afterleptsel  , "HT_vs_Mll_mumue_afterleptsel",   datasetName, "HT", 100,0., 1000, "Mll",100,0., 1000.);
  thehistomanag->CreateHisto2D(HT_vs_Mll_eemu_afterleptsel   , "HT_vs_Mll_eemu_afterleptsel",    datasetName, "HT", 100,0., 1000, "Mll",100,0., 1000.);
  thehistomanag->CreateHisto2D(HT_vs_Mll_eee_afterleptsel    , "HT_vs_Mll_eee_afterleptsel",     datasetName, "HT", 100,0., 1000, "Mll",100,0., 1000.);
  
  
  thehistomanag->CreateHisto2D(InvM_ll_vs_mWT_mumumu_afterleptsel, "InvM_ll_vs_mWT_mumumu_afterleptsel", datasetName, "Mll", 20,0., 200, "mWT", 20,0., 200.);
  thehistomanag->CreateHisto2D(InvM_ll_vs_mWT_mumue_afterleptsel,  "InvM_ll_vs_mWT_mumue_afterleptsel" , datasetName, "Mll", 20,0., 200, "mWT", 20,0., 200.);
  thehistomanag->CreateHisto2D(InvM_ll_vs_mWT_eemu_afterleptsel,   "InvM_ll_vs_mWT_eemu_afterleptsel"  , datasetName, "Mll", 20,0., 200, "mWT", 20,0., 200.);
  thehistomanag->CreateHisto2D(InvM_ll_vs_mWT_eee_afterleptsel,    "InvM_ll_vs_mWT_eee_afterleptsel"   , datasetName, "Mll", 20,0., 200, "mWT", 20,0., 200.);
  //TheTree->Write();
}
 
void ProofSelectorMyCutFlow::defineHistoPointer(int thechannel){

  //cout << "thechannel " << thechannel << endl;
  if( thechannel == 0){
    pCutFlow                        = &CutFlow_mumumu ;
    pErrCutFlow                     = &ErrCutFlow_mumumu ;
    pDijetInvM_afterleptsel_inZpeak = &DijetInvM_mumumu_afterleptsel_inZpeak;
    pMt_afterbjetsel                = &Mt_mumumu_afterbjetsel;
    pMt_afterbjetveto               = &Mt_mumumu_afterbjetveto;
    pNJet_afterZsel                 = &NJet_mumumu_afterZsel;  
    pmWT_afterZsel                 = &mWT_mumumu_afterZsel;  
    pNJet_afterbsel                 = &NJet_mumumu_afterbsel;
    pNJet_afterleptsel_mWT110       = &NJet_mumumu_afterleptsel_mWT110;
    pNLept_afterbsel                = &NLept_mumumu_afterbsel;
    pNBJet_afterZsel                = &NBJet_mumumu_afterZsel ;
    pNBJet_afterjetsel              = &NBJet_mumumu_afterjetsel;
    pNBJet_afterjetsel_bjets        = &NBJet_mumumu_afterjetsel_bjets;
    pNBJet_afterjetsel_cjets        = &NBJet_mumumu_afterjetsel_cjets;
    pNBJet_afterjetsel_ljets        = &NBJet_mumumu_afterjetsel_ljets;
    pBJetDiscri_afterjetsel_bjets   = &BJetDiscri_mumumu_afterjetsel_bjets;
    pBJetDiscri_afterjetsel_cjets   = &BJetDiscri_mumumu_afterjetsel_cjets;
    pBJetDiscri_afterjetsel_ljets   = &BJetDiscri_mumumu_afterjetsel_ljets;
    pNBJet_afterleptsel_mWT110      = &NBJet_mumumu_afterleptsel_mWT110;
    pNBJet_afterleptsel             = &NBJet_mumumu_afterleptsel;
    pNJet_afterleptsel              = &NJet_mumumu_afterleptsel;
    pNvtx_afterleptsel              = &Nvtx_mumumu_afterleptsel;
    pInvM_ll_afterleptsel           = &InvM_ll_mumumu_afterleptsel;  
    pInvM_ll_afterleptsel_mWT110    = &InvM_ll_mumumu_afterleptsel_mWT110;  
    pInvM_ll_afterleptsel_lowbin    = &InvM_ll_mumumu_afterleptsel_lowbin;  
    pInvM_ll_afterleptsel_highSumPt = &InvM_ll_mumumu_afterleptsel_highSumPt;  
    pInvM_ll_afterjetsel            = &InvM_ll_mumumu_afterjetsel;
    pInvM_ll_afterbjetsel           = &InvM_ll_mumumu_afterbjetsel;  
    pLeptPt_afterleptsel            = &LeptPt_mumumu_afterleptsel;  
    pLeptPt_afterjetsel             = &LeptPt_mumumu_afterjetsel;  
    pLeptPt_afterbjetsel            = &LeptPt_mumumu_afterbjetsel;  
    pLeptZPt_afterleptsel           = &LeptZPt_mumumu_afterleptsel;  
    pLeptZPt_afterjetsel            = &LeptZPt_mumumu_afterjetsel;  
    pLeptZPt_afterbjetsel           = &LeptZPt_mumumu_afterbjetsel;
    pLeptWPt_afterleptsel           = &LeptWPt_mumumu_afterleptsel;  
    pLeptWPt_afterjetsel            = &LeptWPt_mumumu_afterjetsel;
    pLeptWPt_afterbjetsel           = &LeptWPt_mumumu_afterbjetsel;
    pLeptWPt_afterbjetveto          = &LeptWPt_mumumu_afterbjetveto;
    pLeptWPt_afterleptsel_mWT110    = &LeptWPt_mumumu_afterleptsel_mWT110;
    pJetPt_afterleptsel             = &JetPt_mumumu_afterleptsel;
    pJetPt_afterjetsel              = &JetPt_mumumu_afterjetsel;
    pJetPt_afterbjetsel             = &JetPt_mumumu_afterbjetsel;
    pJetPt_afterbjetveto            = &JetPt_mumumu_afterbjetveto;  
    pJetEta_afterleptsel            = &JetEta_mumumu_afterleptsel;
    pJetEta_afterjetsel             = &JetEta_mumumu_afterjetsel;  
    pJetEta_afterbjetsel            = &JetEta_mumumu_afterbjetsel;  
    pJetEta_afterbjetveto           = &JetEta_mumumu_afterbjetveto;
    pHT_afterleptsel                = &HT_mumumu_afterleptsel;  
    pHT_afterjetsel                 = &HT_mumumu_afterjetsel ;
    pHT_afterbjetsel                = &HT_mumumu_afterbjetsel;  
    pHT_afterbjetveto               = &HT_mumumu_afterbjetveto;  
    pMET_afterleptsel               = &MET_mumumu_afterleptsel;
    pMET_afterleptsel_mWT110        = &MET_mumumu_afterleptsel_mWT110;
    pMET_afterjetsel                = &MET_mumumu_afterjetsel;
    pMET_afterbjetsel               = &MET_mumumu_afterbjetsel;
    pAsym_afterbjetsel              = &Asym_mumumu_afterbjetsel;  
    pmWT_afterjetsel                = &mWT_mumumu_afterjetsel;
    pRecoPtZ_afterbjetsel           = &RecoPtZ_mumumu_afterbjetsel;
    pRecoPtZ_afterbjetveto          = &RecoPtZ_mumumu_afterbjetveto;
    pRecoPtZ_afterleptsel           = &RecoPtZ_mumumu_afterleptsel;  
    pRecoPtZ_afterleptsel_nojet     = &RecoPtZ_mumumu_afterleptsel_nojet;
    pRecoTopMass_afterbjetsel       = &RecoTopMass_mumumu_afterbjetsel;
    pRecoTopMass_afterbjetveto      = &RecoTopMass_mumumu_afterbjetveto;  
    pdeltaPhilb_afterbjetsel        = &deltaPhilb_mumumu_afterbjetsel;
    pdeltaPhilj_afterbjetveto       = &deltaPhilj_mumumu_afterbjetveto;  
    pdeltaR_afterleptsel            = &deltaR_mumumu_afterleptsel;  
    pdeltaRLeptJet_afterleptsel_mWT110 = &deltaRLeptJet_mumumu_afterleptsel_mWT110;  
    pdeltaRLeptMet_afterleptsel_mWT110 = &deltaRLeptMet_mumumu_afterleptsel_mWT110;  
    pWmissAssing_afterleptsel       = &WmissAssing_mumumu_afterleptsel;  
    pmWT_afterleptsel               = &mWT_mumumu_afterleptsel; 
    pmWT_afterbjetsel               = &mWT_mumumu_afterbjetsel;  
    pmWT_afterbjetveto              = &mWT_mumumu_afterbjetveto;  
    pcosThetaStar_afterbjetsel      = &cosThetaStar_mumumu_afterbjetsel; 
    pCharge_afterleptsel            = &Charge_mumumu_afterleptsel;  
    pCharge_afterleptsel_mWT110     = &Charge_mumumu_afterleptsel_mWT110;     
    pNvertex                        = &Nvertex ; 
    pInvM_ll_vs_mWT_afterleptsel    = &InvM_ll_vs_mWT_mumumu_afterleptsel;  
    pHT_vs_MET_afterleptsel         = &HT_vs_MET_mumumu_afterleptsel;  
    pHT_vs_NJet_afterleptsel        = &HT_vs_NJet_mumumu_afterleptsel;  
    pHT_vs_NBJet_afterleptsel       = &HT_vs_NBJet_mumumu_afterleptsel;  
    pHT_vs_LeptPt_afterleptsel      = &HT_vs_LeptPt_mumumu_afterleptsel;
    pHT_vs_JetPt_afterleptsel       = &HT_vs_JetPt_mumumu_afterleptsel;
    pHT_vs_Mll_afterleptsel         = &HT_vs_Mll_mumumu_afterleptsel;
  }else if(thechannel == 1){ 
    pCutFlow                        = &CutFlow_mumue ;
    pErrCutFlow                     = &ErrCutFlow_mumue ;   
    pDijetInvM_afterleptsel_inZpeak = &DijetInvM_mumue_afterleptsel_inZpeak;
    pMt_afterbjetsel                = &Mt_mumue_afterbjetsel;
    pMt_afterbjetveto               = &Mt_mumue_afterbjetveto;
    pNJet_afterZsel                 = &NJet_mumue_afterZsel;  
    pmWT_afterZsel                  = &mWT_mumue_afterZsel;  
    pNJet_afterbsel                 = &NJet_mumue_afterbsel;
    pNJet_afterleptsel_mWT110       = &NJet_mumue_afterleptsel_mWT110;
    pNLept_afterbsel                = &NLept_mumue_afterbsel;
    pNBJet_afterZsel                = &NBJet_mumue_afterZsel ;
    pNBJet_afterjetsel              = &NBJet_mumue_afterjetsel;
    pNBJet_afterjetsel_bjets        = &NBJet_mumue_afterjetsel_bjets;
    pNBJet_afterjetsel_cjets        = &NBJet_mumue_afterjetsel_cjets;
    pNBJet_afterjetsel_ljets        = &NBJet_mumue_afterjetsel_ljets;
    pBJetDiscri_afterjetsel_bjets   = &BJetDiscri_mumue_afterjetsel_bjets;
    pBJetDiscri_afterjetsel_cjets   = &BJetDiscri_mumue_afterjetsel_cjets;
    pBJetDiscri_afterjetsel_ljets   = &BJetDiscri_mumue_afterjetsel_ljets;
    pNBJet_afterleptsel             = &NBJet_mumue_afterleptsel ;
    pNJet_afterleptsel              = &NJet_mumue_afterleptsel ;
    pNvtx_afterleptsel              = &Nvtx_mumue_afterleptsel;
    pInvM_ll_afterleptsel           = &InvM_ll_mumue_afterleptsel;  
    pInvM_ll_afterleptsel_mWT110    = &InvM_ll_mumue_afterleptsel_mWT110;  
    pInvM_ll_afterleptsel_lowbin    = &InvM_ll_mumue_afterleptsel_lowbin;  
    pInvM_ll_afterleptsel_highSumPt = &InvM_ll_mumue_afterleptsel_highSumPt;  
    pInvM_ll_afterjetsel            = &InvM_ll_mumue_afterjetsel;
    pInvM_ll_afterbjetsel           = &InvM_ll_mumue_afterbjetsel;  
    pLeptPt_afterleptsel            = &LeptPt_mumue_afterleptsel;  
    pLeptPt_afterjetsel             = &LeptPt_mumue_afterjetsel;  
    pLeptPt_afterbjetsel            = &LeptPt_mumue_afterbjetsel;  
    pLeptZPt_afterleptsel           = &LeptZPt_mumue_afterleptsel;  
    pLeptZPt_afterjetsel            = &LeptZPt_mumue_afterjetsel;  
    pLeptZPt_afterbjetsel           = &LeptZPt_mumue_afterbjetsel;
    pLeptWPt_afterleptsel           = &LeptWPt_mumue_afterleptsel;  
    pLeptWPt_afterjetsel            = &LeptWPt_mumue_afterjetsel;
    pLeptWPt_afterbjetsel           = &LeptWPt_mumue_afterbjetsel;
    pLeptWPt_afterbjetveto          = &LeptWPt_mumue_afterbjetveto;
    pLeptWPt_afterleptsel_mWT110    = &LeptWPt_mumue_afterleptsel_mWT110;
    pJetPt_afterleptsel             = &JetPt_mumue_afterleptsel;
    pJetPt_afterjetsel              = &JetPt_mumue_afterjetsel;
    pJetPt_afterbjetsel             = &JetPt_mumue_afterbjetsel;
    pJetPt_afterbjetveto            = &JetPt_mumue_afterbjetveto;  
    pJetEta_afterleptsel            = &JetEta_mumue_afterleptsel;
    pJetEta_afterjetsel             = &JetEta_mumue_afterjetsel;  
    pJetEta_afterbjetsel            = &JetEta_mumue_afterbjetsel;  
    pJetEta_afterbjetveto           = &JetEta_mumue_afterbjetveto;
    pHT_afterleptsel                = &HT_mumue_afterleptsel;  
    pHT_afterjetsel                 = &HT_mumue_afterjetsel ;
    pHT_afterbjetsel                = &HT_mumue_afterbjetsel;  
    pHT_afterbjetveto               = &HT_mumue_afterbjetveto;  
    pMET_afterleptsel               = &MET_mumue_afterleptsel;
    pMET_afterleptsel_mWT110        = &MET_mumue_afterleptsel_mWT110;
    pMET_afterjetsel                = &MET_mumue_afterjetsel;
    pMET_afterbjetsel                = &MET_mumue_afterbjetsel;
    pAsym_afterbjetsel              = &Asym_mumue_afterbjetsel;  
    pmWT_afterjetsel                = &mWT_mumue_afterjetsel;
    pRecoPtZ_afterbjetsel           = &RecoPtZ_mumue_afterbjetsel;
    pRecoPtZ_afterbjetveto          = &RecoPtZ_mumue_afterbjetveto;
    pRecoPtZ_afterleptsel           = &RecoPtZ_mumue_afterleptsel;  
    pRecoPtZ_afterleptsel_nojet     = &RecoPtZ_mumue_afterleptsel_nojet;
    pRecoTopMass_afterbjetsel       = &RecoTopMass_mumue_afterbjetsel;
    pRecoTopMass_afterbjetveto      = &RecoTopMass_mumue_afterbjetveto;  
    pdeltaPhilb_afterbjetsel        = &deltaPhilb_mumue_afterbjetsel;
    pdeltaPhilj_afterbjetveto       = &deltaPhilj_mumue_afterbjetveto;  
    pdeltaR_afterleptsel            = &deltaR_mumue_afterleptsel;  
    pdeltaRLeptJet_afterleptsel_mWT110 = &deltaRLeptJet_mumue_afterleptsel_mWT110;  
    pdeltaRLeptMet_afterleptsel_mWT110 = &deltaRLeptMet_mumue_afterleptsel_mWT110;  
    pWmissAssing_afterleptsel       = &WmissAssing_mumue_afterleptsel;  
    pmWT_afterleptsel               = &mWT_mumue_afterleptsel; 
    pmWT_afterbjetsel               = &mWT_mumue_afterbjetsel;  
    pmWT_afterbjetveto              = &mWT_mumue_afterbjetveto; 
    pcosThetaStar_afterbjetsel      = &cosThetaStar_mumue_afterbjetsel;  
    pCharge_afterleptsel            = &Charge_mumue_afterleptsel;  
    pCharge_afterleptsel_mWT110     = &Charge_mumue_afterleptsel_mWT110;     
    pNvertex                        = &Nvertex ; 
    pInvM_ll_vs_mWT_afterleptsel    = &InvM_ll_vs_mWT_mumue_afterleptsel;  
    pHT_vs_MET_afterleptsel         = &HT_vs_MET_mumue_afterleptsel;  
    pHT_vs_NJet_afterleptsel        = &HT_vs_NJet_mumue_afterleptsel;  
    pHT_vs_NBJet_afterleptsel       = &HT_vs_NBJet_mumue_afterleptsel;  
    pHT_vs_LeptPt_afterleptsel      = &HT_vs_LeptPt_mumue_afterleptsel;
    pHT_vs_JetPt_afterleptsel       = &HT_vs_JetPt_mumue_afterleptsel;
    pHT_vs_Mll_afterleptsel         = &HT_vs_Mll_mumue_afterleptsel;
  }else if(thechannel == 2){   
    pCutFlow                        = &CutFlow_eemu ;
    pErrCutFlow                     = &ErrCutFlow_eemu ;   
    pDijetInvM_afterleptsel_inZpeak = &DijetInvM_eemu_afterleptsel_inZpeak;
    pMt_afterbjetsel                = &Mt_eemu_afterbjetsel;
    pMt_afterbjetveto               = &Mt_eemu_afterbjetveto;
    pNJet_afterZsel                 = &NJet_eemu_afterZsel;  
    pmWT_afterZsel                  = &mWT_eemu_afterZsel;  
    pNJet_afterbsel                 = &NJet_eemu_afterbsel;
    pNJet_afterleptsel_mWT110       = &NJet_eemu_afterleptsel_mWT110;
    pNLept_afterbsel                = &NLept_eemu_afterbsel;
    pNBJet_afterZsel                = &NBJet_eemu_afterZsel ;
    pNBJet_afterjetsel              = &NBJet_eemu_afterjetsel;
    pNBJet_afterjetsel_bjets        = &NBJet_eemu_afterjetsel_bjets;
    pNBJet_afterjetsel_cjets        = &NBJet_eemu_afterjetsel_cjets;
    pNBJet_afterjetsel_ljets        = &NBJet_eemu_afterjetsel_ljets;
    pBJetDiscri_afterjetsel_bjets   = &BJetDiscri_eemu_afterjetsel_bjets;
    pBJetDiscri_afterjetsel_cjets   = &BJetDiscri_eemu_afterjetsel_cjets;
    pBJetDiscri_afterjetsel_ljets   = &BJetDiscri_eemu_afterjetsel_ljets;
    pNBJet_afterleptsel             = &NBJet_eemu_afterleptsel ;
    pNJet_afterleptsel              = &NJet_eemu_afterleptsel ;
    pNvtx_afterleptsel              = &Nvtx_eemu_afterleptsel;
    pInvM_ll_afterleptsel           = &InvM_ll_eemu_afterleptsel;  
    pInvM_ll_afterleptsel_mWT110    = &InvM_ll_eemu_afterleptsel_mWT110;  
    pInvM_ll_afterleptsel_lowbin    = &InvM_ll_eemu_afterleptsel_lowbin;  
    pInvM_ll_afterleptsel_highSumPt = &InvM_ll_eemu_afterleptsel_highSumPt;  
    pInvM_ll_afterjetsel            = &InvM_ll_eemu_afterjetsel;
    pInvM_ll_afterbjetsel           = &InvM_ll_eemu_afterbjetsel;  
    pLeptPt_afterleptsel            = &LeptPt_eemu_afterleptsel;  
    pLeptPt_afterjetsel             = &LeptPt_eemu_afterjetsel;  
    pLeptPt_afterbjetsel            = &LeptPt_eemu_afterbjetsel;  
    pLeptZPt_afterleptsel           = &LeptZPt_eemu_afterleptsel;  
    pLeptZPt_afterjetsel            = &LeptZPt_eemu_afterjetsel;  
    pLeptZPt_afterbjetsel           = &LeptZPt_eemu_afterbjetsel;
    pLeptWPt_afterleptsel           = &LeptWPt_eemu_afterleptsel;  
    pLeptWPt_afterjetsel            = &LeptWPt_eemu_afterjetsel;
    pLeptWPt_afterbjetsel           = &LeptWPt_eemu_afterbjetsel;
    pLeptWPt_afterbjetveto          = &LeptWPt_eemu_afterbjetveto;
    pLeptWPt_afterleptsel_mWT110    = &LeptWPt_eemu_afterleptsel_mWT110;
    pJetPt_afterleptsel             = &JetPt_eemu_afterleptsel;
    pJetPt_afterjetsel              = &JetPt_eemu_afterjetsel;
    pJetPt_afterbjetsel             = &JetPt_eemu_afterbjetsel;
    pJetPt_afterbjetveto            = &JetPt_eemu_afterbjetveto;  
    pJetEta_afterleptsel            = &JetEta_eemu_afterleptsel;
    pJetEta_afterjetsel             = &JetEta_eemu_afterjetsel;  
    pJetEta_afterbjetsel            = &JetEta_eemu_afterbjetsel;  
    pJetEta_afterbjetveto           = &JetEta_eemu_afterbjetveto;
    pHT_afterleptsel                = &HT_eemu_afterleptsel;  
    pHT_afterjetsel                 = &HT_eemu_afterjetsel ;
    pHT_afterbjetsel                = &HT_eemu_afterbjetsel;  
    pHT_afterbjetveto               = &HT_eemu_afterbjetveto;  
    pMET_afterleptsel               = &MET_eemu_afterleptsel;
    pMET_afterleptsel_mWT110        = &MET_eemu_afterleptsel_mWT110;
    pMET_afterjetsel                = &MET_eemu_afterjetsel;
    pMET_afterbjetsel                = &MET_eemu_afterbjetsel;
    pAsym_afterbjetsel              = &Asym_eemu_afterbjetsel;  
    pmWT_afterjetsel                = &mWT_eemu_afterjetsel;
    pRecoPtZ_afterbjetsel           = &RecoPtZ_eemu_afterbjetsel;
    pRecoPtZ_afterbjetveto          = &RecoPtZ_eemu_afterbjetveto;
    pRecoPtZ_afterleptsel           = &RecoPtZ_eemu_afterleptsel;  
    pRecoPtZ_afterleptsel_nojet     = &RecoPtZ_eemu_afterleptsel_nojet;
    pRecoTopMass_afterbjetsel       = &RecoTopMass_eemu_afterbjetsel;
    pRecoTopMass_afterbjetveto      = &RecoTopMass_eemu_afterbjetveto;  
    pdeltaPhilb_afterbjetsel        = &deltaPhilb_eemu_afterbjetsel;
    pdeltaPhilj_afterbjetveto       = &deltaPhilj_eemu_afterbjetveto;  
    pdeltaR_afterleptsel            = &deltaR_eemu_afterleptsel;  
    pdeltaRLeptJet_afterleptsel_mWT110 = &deltaRLeptJet_eemu_afterleptsel_mWT110;  
    pdeltaRLeptMet_afterleptsel_mWT110 = &deltaRLeptMet_eemu_afterleptsel_mWT110;  
    pWmissAssing_afterleptsel       = &WmissAssing_eemu_afterleptsel;  
    pmWT_afterleptsel               = &mWT_eemu_afterleptsel; 
    pmWT_afterbjetsel               = &mWT_eemu_afterbjetsel;  
    pmWT_afterbjetveto              = &mWT_eemu_afterbjetveto;  
    pcosThetaStar_afterbjetsel      = &cosThetaStar_eemu_afterbjetsel;  
    pCharge_afterleptsel            = &Charge_eemu_afterleptsel;  
    pCharge_afterleptsel_mWT110     = &Charge_eemu_afterleptsel_mWT110;     
    pNvertex                        = &Nvertex ; 
    pInvM_ll_vs_mWT_afterleptsel    = &InvM_ll_vs_mWT_eemu_afterleptsel;  
    pHT_vs_MET_afterleptsel         = &HT_vs_MET_eemu_afterleptsel;  
    pHT_vs_NJet_afterleptsel        = &HT_vs_NJet_eemu_afterleptsel;  
    pHT_vs_NBJet_afterleptsel       = &HT_vs_NBJet_eemu_afterleptsel;  
    pHT_vs_LeptPt_afterleptsel      = &HT_vs_LeptPt_eemu_afterleptsel;
    pHT_vs_JetPt_afterleptsel       = &HT_vs_JetPt_eemu_afterleptsel;
    pHT_vs_Mll_afterleptsel         = &HT_vs_Mll_eemu_afterleptsel;
  }else if(thechannel == 3){    
    pCutFlow                        = &CutFlow_eee ;
    pErrCutFlow                     = &ErrCutFlow_eee ;   
    pDijetInvM_afterleptsel_inZpeak = &DijetInvM_eee_afterleptsel_inZpeak;
    pMt_afterbjetsel                = &Mt_eee_afterbjetsel;
    pMt_afterbjetveto               = &Mt_eee_afterbjetveto;
    pNJet_afterZsel                 = &NJet_eee_afterZsel;  
    pmWT_afterZsel                  = &mWT_eee_afterZsel;  
    pNJet_afterbsel                 = &NJet_eee_afterbsel;
    pNJet_afterleptsel_mWT110       = &NJet_eee_afterleptsel_mWT110;
    pNLept_afterbsel                = &NLept_eee_afterbsel;
    pNBJet_afterZsel                = &NBJet_eee_afterZsel ;
    pNBJet_afterjetsel              = &NBJet_eee_afterjetsel;
    pNBJet_afterjetsel_bjets        = &NBJet_eee_afterjetsel_bjets;
    pNBJet_afterjetsel_cjets        = &NBJet_eee_afterjetsel_cjets;
    pNBJet_afterjetsel_ljets        = &NBJet_eee_afterjetsel_ljets;
    pBJetDiscri_afterjetsel_bjets   = &BJetDiscri_eee_afterjetsel_bjets;
    pBJetDiscri_afterjetsel_cjets   = &BJetDiscri_eee_afterjetsel_cjets;
    pBJetDiscri_afterjetsel_ljets   = &BJetDiscri_eee_afterjetsel_ljets;
    pNBJet_afterleptsel             = &NBJet_eee_afterleptsel ;
    pNJet_afterleptsel              = &NJet_eee_afterleptsel ;
    pNvtx_afterleptsel              = &Nvtx_eee_afterleptsel;
    pInvM_ll_afterleptsel           = &InvM_ll_eee_afterleptsel;  
    pInvM_ll_afterleptsel_mWT110    = &InvM_ll_eee_afterleptsel_mWT110;  
    pInvM_ll_afterleptsel_lowbin    = &InvM_ll_eee_afterleptsel_lowbin;  
    pInvM_ll_afterleptsel_highSumPt = &InvM_ll_eee_afterleptsel_highSumPt;  
    pInvM_ll_afterjetsel            = &InvM_ll_eee_afterjetsel;
    pInvM_ll_afterbjetsel           = &InvM_ll_eee_afterbjetsel;  
    pLeptPt_afterleptsel            = &LeptPt_eee_afterleptsel;  
    pLeptPt_afterjetsel             = &LeptPt_eee_afterjetsel;  
    pLeptPt_afterbjetsel            = &LeptPt_eee_afterbjetsel;  
    pLeptZPt_afterleptsel           = &LeptZPt_eee_afterleptsel;  
    pLeptZPt_afterjetsel            = &LeptZPt_eee_afterjetsel;  
    pLeptZPt_afterbjetsel           = &LeptZPt_eee_afterbjetsel;
    pLeptWPt_afterleptsel           = &LeptWPt_eee_afterleptsel;  
    pLeptWPt_afterjetsel            = &LeptWPt_eee_afterjetsel;
    pLeptWPt_afterbjetsel           = &LeptWPt_eee_afterbjetsel;
    pLeptWPt_afterbjetveto          = &LeptWPt_eee_afterbjetveto;
    pLeptWPt_afterleptsel_mWT110    = &LeptWPt_eee_afterleptsel_mWT110;
    pJetPt_afterleptsel             = &JetPt_eee_afterleptsel;
    pJetPt_afterjetsel              = &JetPt_eee_afterjetsel;
    pJetPt_afterbjetsel             = &JetPt_eee_afterbjetsel;
    pJetPt_afterbjetveto            = &JetPt_eee_afterbjetveto;  
    pJetEta_afterleptsel            = &JetEta_eee_afterleptsel;
    pJetEta_afterjetsel             = &JetEta_eee_afterjetsel;  
    pJetEta_afterbjetsel            = &JetEta_eee_afterbjetsel;  
    pJetEta_afterbjetveto           = &JetEta_eee_afterbjetveto;
    pHT_afterleptsel                = &HT_eee_afterleptsel;  
    pHT_afterjetsel                 = &HT_eee_afterjetsel ;
    pHT_afterbjetsel                = &HT_eee_afterbjetsel;  
    pHT_afterbjetveto               = &HT_eee_afterbjetveto;  
    pMET_afterleptsel               = &MET_eee_afterleptsel;
    pMET_afterleptsel_mWT110        = &MET_eee_afterleptsel_mWT110;
    pMET_afterjetsel                = &MET_eee_afterjetsel;
    pMET_afterbjetsel                = &MET_eee_afterbjetsel;
    pAsym_afterbjetsel              = &Asym_eee_afterbjetsel;  
    pmWT_afterjetsel                = &mWT_eee_afterjetsel;
    pRecoPtZ_afterbjetsel           = &RecoPtZ_eee_afterbjetsel;
    pRecoPtZ_afterbjetveto          = &RecoPtZ_eee_afterbjetveto;
    pRecoPtZ_afterleptsel           = &RecoPtZ_eee_afterleptsel;  
    pRecoPtZ_afterleptsel_nojet     = &RecoPtZ_eee_afterleptsel_nojet;
    pRecoTopMass_afterbjetsel       = &RecoTopMass_eee_afterbjetsel;
    pRecoTopMass_afterbjetveto      = &RecoTopMass_eee_afterbjetveto;  
    pdeltaPhilb_afterbjetsel        = &deltaPhilb_eee_afterbjetsel;
    pdeltaPhilj_afterbjetveto       = &deltaPhilj_eee_afterbjetveto;  
    pdeltaR_afterleptsel            = &deltaR_eee_afterleptsel;  
    pdeltaRLeptJet_afterleptsel_mWT110 = &deltaRLeptJet_eee_afterleptsel_mWT110;  
    pdeltaRLeptMet_afterleptsel_mWT110 = &deltaRLeptMet_eee_afterleptsel_mWT110;  
    pWmissAssing_afterleptsel       = &WmissAssing_eee_afterleptsel;  
    pmWT_afterleptsel               = &mWT_eee_afterleptsel; 
    pmWT_afterbjetsel               = &mWT_eee_afterbjetsel;  
    pmWT_afterbjetveto              = &mWT_eee_afterbjetveto;
    pcosThetaStar_afterbjetsel      = &cosThetaStar_eee_afterbjetsel;   
    pCharge_afterleptsel            = &Charge_eee_afterleptsel;  
    pCharge_afterleptsel_mWT110     = &Charge_eee_afterleptsel_mWT110;     
    pNvertex                        = &Nvertex ; 
    pInvM_ll_vs_mWT_afterleptsel    = &InvM_ll_vs_mWT_eee_afterleptsel;  
    pHT_vs_MET_afterleptsel         = &HT_vs_MET_eee_afterleptsel;  
    pHT_vs_NJet_afterleptsel        = &HT_vs_NJet_eee_afterleptsel;  
    pHT_vs_NBJet_afterleptsel       = &HT_vs_NBJet_eee_afterleptsel;  
    pHT_vs_LeptPt_afterleptsel      = &HT_vs_LeptPt_eee_afterleptsel;
    pHT_vs_JetPt_afterleptsel       = &HT_vs_JetPt_eee_afterleptsel;
    pHT_vs_Mll_afterleptsel         = &HT_vs_Mll_eee_afterleptsel;
  }


}



void ProofSelectorMyCutFlow::setNullHistoPointer(){

    pCutFlow                        = 0;
    pErrCutFlow                     = 0;
    pDijetInvM_afterleptsel_inZpeak = 0;
    pMt_afterbjetsel                = 0;
    pMt_afterbjetveto               = 0;
    pNJet_afterZsel                 = 0;
    pmWT_afterZsel                  = 0;
    pNJet_afterbsel                 = 0;
    pNJet_afterleptsel_mWT110       = 0;
    pNLept_afterbsel                = 0;
    pNBJet_afterZsel                = 0;
    pNBJet_afterjetsel              = 0;
    pNBJet_afterjetsel_bjets        = 0;
    pNBJet_afterjetsel_cjets        = 0;
    pNBJet_afterjetsel_ljets        = 0;
    pBJetDiscri_afterjetsel_bjets   = 0;
    pBJetDiscri_afterjetsel_cjets   = 0;
    pBJetDiscri_afterjetsel_ljets   = 0;
    pNBJet_afterleptsel             = 0;
    pNJet_afterleptsel              = 0;
    pNvtx_afterleptsel              = 0;
    pInvM_ll_afterleptsel           = 0;
    pInvM_ll_afterleptsel_mWT110    = 0;
    pInvM_ll_afterleptsel_lowbin    = 0;
    pInvM_ll_afterleptsel_highSumPt = 0;  
    pInvM_ll_afterjetsel            = 0;
    pInvM_ll_afterbjetsel           = 0;
    pLeptPt_afterleptsel            = 0;
    pLeptPt_afterjetsel             = 0;
    pLeptPt_afterbjetsel            = 0;
    pLeptZPt_afterleptsel           = 0;
    pLeptZPt_afterjetsel            = 0;
    pLeptZPt_afterbjetsel           = 0;
    pLeptWPt_afterleptsel           = 0;
    pLeptWPt_afterjetsel            = 0;
    pLeptWPt_afterbjetsel           = 0;
    pLeptWPt_afterbjetveto          = 0;
    pLeptWPt_afterleptsel_mWT110    = 0;
    pJetPt_afterleptsel             = 0;
    pJetPt_afterjetsel              = 0;
    pJetPt_afterbjetsel             = 0;
    pJetPt_afterbjetveto            = 0;
    pJetEta_afterleptsel            = 0;
    pJetEta_afterjetsel             = 0;
    pJetEta_afterbjetsel            = 0;
    pJetEta_afterbjetveto           = 0;
    pHT_afterleptsel                = 0;
    pHT_afterjetsel                 = 0;
    pHT_afterbjetsel                = 0;
    pHT_afterbjetveto               = 0;
    pMET_afterleptsel               = 0;
    pMET_afterleptsel_mWT110        = 0;
    pMET_afterjetsel                = 0;
    pMET_afterbjetsel                =
    pAsym_afterbjetsel              = 0;
    pmWT_afterjetsel                = 0;
    pRecoPtZ_afterbjetsel           = 0;
    pRecoPtZ_afterbjetveto          = 0;
    pRecoPtZ_afterleptsel           = 0;
    pRecoPtZ_afterleptsel_nojet     = 0;
    pRecoTopMass_afterbjetsel       = 0;
    pRecoTopMass_afterbjetveto      = 0;
    pdeltaPhilb_afterbjetsel        = 0;
    pdeltaPhilj_afterbjetveto       = 0;
    pdeltaR_afterleptsel            = 0;
    pdeltaRLeptJet_afterleptsel_mWT110 = 0;
    pdeltaRLeptMet_afterleptsel_mWT110 = 0;
    pWmissAssing_afterleptsel       = 0;
    pmWT_afterleptsel               = 0;
    pmWT_afterbjetsel               = 0;
    pmWT_afterbjetveto              = 0;
    pcosThetaStar_afterbjetsel      = 0;
    pCharge_afterleptsel            = 0;
    pCharge_afterleptsel_mWT110     = 0;
    pNvertex                        = 0;
    pInvM_ll_vs_mWT_afterleptsel    = 0;
    pHT_vs_MET_afterleptsel         = 0;
    pHT_vs_NJet_afterleptsel        = 0;
    pHT_vs_NBJet_afterleptsel       = 0;
    pHT_vs_LeptPt_afterleptsel      = 0;
    pHT_vs_JetPt_afterleptsel       = 0;
    pHT_vs_Mll_afterleptsel         = 0;
 


}




void  ProofSelectorMyCutFlow::WriteTheHisto(TFile* theoutputfile, HistoManager *thehistomanag){
  
  
  theoutputfile->cd();
  //TheTree->Write();
  thehistomanag->WriteMyHisto(SelABjet, "all");
  thehistomanag->WriteMyHisto(SelABjet_afterjetsel, "all");
  thehistomanag->WriteMyHisto(Ntrilept_mumumu, "all");
  thehistomanag->WriteMyHisto(Ntrileptnoniso_mumumu, "all");
  
   
  thehistomanag->WriteMyHisto(Nvertex, "all");  
  
  
  thehistomanag->WriteMyHisto(InvM_ll_mumumu_afterdileptsel, "all");  
  
  
  thehistomanag->WriteMyHisto(CutFlow_mumumu, "all" );
  thehistomanag->WriteMyHisto(CutFlow_mumue,  "all" );
  thehistomanag->WriteMyHisto(CutFlow_eemu,   "all" );
  thehistomanag->WriteMyHisto(CutFlow_eee,    "all" );
  
  thehistomanag->WriteMyHisto(NJetLight_mumu, "all" );
  thehistomanag->WriteMyHisto(NJetLight_ee, "all" ); 
  thehistomanag->WriteMyHisto(NJetHeavy_mumu, "all" );
  thehistomanag->WriteMyHisto(NJetHeavy_ee, "all" );  
  thehistomanag->WriteMyHisto(FlavComp_mumu, "all" );
  thehistomanag->WriteMyHisto(FlavComp_ee, "all" );
  
  
  thehistomanag->WriteMyHisto(PU_before_mumumu, "all");  
  thehistomanag->WriteMyHisto(PU_before_mumue, "all");  
  thehistomanag->WriteMyHisto(PU_before_eemu, "all");  
  thehistomanag->WriteMyHisto(PU_before_eee, "all");  

  thehistomanag->WriteMyHisto(PU_intime_mumumu, "all");  
  thehistomanag->WriteMyHisto(PU_intime_mumue, "all");  
  thehistomanag->WriteMyHisto(PU_intime_eemu, "all");  
  thehistomanag->WriteMyHisto(PU_intime_eee, "all");  

  thehistomanag->WriteMyHisto(PU_after_mumumu, "all");  
  thehistomanag->WriteMyHisto(PU_after_mumue, "all");  
  thehistomanag->WriteMyHisto(PU_after_eemu, "all");  
  thehistomanag->WriteMyHisto(PU_after_eee, "all");  
  
  
  thehistomanag->WriteMyHisto(NVtx_mumumu, "all");  
  thehistomanag->WriteMyHisto(NVtx_mumue, "all");  
  thehistomanag->WriteMyHisto(NVtx_eemu, "all");  
  thehistomanag->WriteMyHisto(NVtx_eee, "all");  
  
  thehistomanag->WriteMyHisto(Nvtx_mumumu_afterleptsel, "all" );
  thehistomanag->WriteMyHisto(Nvtx_mumue_afterleptsel , "all" );
  thehistomanag->WriteMyHisto(Nvtx_eemu_afterleptsel  , "all" );
  thehistomanag->WriteMyHisto(Nvtx_eee_afterleptsel   , "all" );
  
  
  thehistomanag->WriteMyHisto(ErrCutFlow_mumumu,  "all");
  thehistomanag->WriteMyHisto(ErrCutFlow_mumue,   "all");
  thehistomanag->WriteMyHisto(ErrCutFlow_eemu,    "all");
  thehistomanag->WriteMyHisto(ErrCutFlow_eee,     "all");
  
  
  
  thehistomanag->WriteMyHisto(Mt_mumumu_afterbjetsel, "all"); 
  thehistomanag->WriteMyHisto(Mt_mumue_afterbjetsel , "all");
  thehistomanag->WriteMyHisto(Mt_eemu_afterbjetsel  , "all");
  thehistomanag->WriteMyHisto(Mt_eee_afterbjetsel   , "all");
    
  thehistomanag->WriteMyHisto(Mt_mumumu_afterbjetveto, "all"); 
  thehistomanag->WriteMyHisto(Mt_mumue_afterbjetveto , "all");
  thehistomanag->WriteMyHisto(Mt_eemu_afterbjetveto  , "all");
  thehistomanag->WriteMyHisto(Mt_eee_afterbjetveto   , "all");
    
    
  
  thehistomanag->WriteMyHisto(NJet_mumumu_afterZsel,"all");
  thehistomanag->WriteMyHisto(NJet_mumue_afterZsel, "all");
  thehistomanag->WriteMyHisto(NJet_eemu_afterZsel,  "all");
  thehistomanag->WriteMyHisto(NJet_eee_afterZsel,   "all");
  
  
  thehistomanag->WriteMyHisto(mWT_mumumu_afterZsel,"all");
  thehistomanag->WriteMyHisto(mWT_mumue_afterZsel, "all");
  thehistomanag->WriteMyHisto(mWT_eemu_afterZsel,  "all");
  thehistomanag->WriteMyHisto(mWT_eee_afterZsel,   "all");
  
  thehistomanag->WriteMyHisto(NJet_mumumu_afterbsel,"all");
  thehistomanag->WriteMyHisto(NJet_mumue_afterbsel, "all");
  thehistomanag->WriteMyHisto(NJet_eemu_afterbsel,  "all");
  thehistomanag->WriteMyHisto(NJet_eee_afterbsel,   "all");
  
  
  thehistomanag->WriteMyHisto(NLept_mumumu_afterbsel,"all");
  thehistomanag->WriteMyHisto(NLept_mumue_afterbsel, "all");
  thehistomanag->WriteMyHisto(NLept_eemu_afterbsel,  "all");
  thehistomanag->WriteMyHisto(NLept_eee_afterbsel,   "all");
  
  thehistomanag->WriteMyHisto(NBJet_mumumu_afterZsel, "all");
  thehistomanag->WriteMyHisto(NBJet_mumue_afterZsel,  "all");
  thehistomanag->WriteMyHisto(NBJet_eemu_afterZsel,   "all");
  thehistomanag->WriteMyHisto(NBJet_eee_afterZsel,    "all");
  
  
  thehistomanag->WriteMyHisto(NBJet_mumumu_afterjetsel, "all");
  thehistomanag->WriteMyHisto(NBJet_mumue_afterjetsel,  "all");
  thehistomanag->WriteMyHisto(NBJet_eemu_afterjetsel,   "all");
  thehistomanag->WriteMyHisto(NBJet_eee_afterjetsel,    "all");
  
  thehistomanag->WriteMyHisto(InvM_ll_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(InvM_ll_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(InvM_ll_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(InvM_ll_eee_afterleptsel,    "all");
  
  thehistomanag->WriteMyHisto(InvM_ll_mumumu_afterleptsel_mWT110, "all");
  thehistomanag->WriteMyHisto(InvM_ll_mumue_afterleptsel_mWT110,  "all");
  thehistomanag->WriteMyHisto(InvM_ll_eemu_afterleptsel_mWT110,   "all");
  thehistomanag->WriteMyHisto(InvM_ll_eee_afterleptsel_mWT110,    "all");
  
  thehistomanag->WriteMyHisto(InvM_ll_mumumu_afterleptsel_lowbin, "all");
  thehistomanag->WriteMyHisto(InvM_ll_mumue_afterleptsel_lowbin,  "all");
  thehistomanag->WriteMyHisto(InvM_ll_eemu_afterleptsel_lowbin,   "all");
  thehistomanag->WriteMyHisto(InvM_ll_eee_afterleptsel_lowbin,    "all");
    
  
  thehistomanag->WriteMyHisto(InvM_ll_mumumu_afterjetsel, "all");
  thehistomanag->WriteMyHisto(InvM_ll_mumue_afterjetsel,  "all");
  thehistomanag->WriteMyHisto(InvM_ll_eemu_afterjetsel,   "all");
  thehistomanag->WriteMyHisto(InvM_ll_eee_afterjetsel,    "all");
  
  thehistomanag->WriteMyHisto(InvM_ll_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(InvM_ll_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(InvM_ll_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(InvM_ll_eee_afterbjetsel,    "all");
    
    
  thehistomanag->WriteMyHisto(LeptPt_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(LeptPt_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(LeptPt_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(LeptPt_eee_afterleptsel,    "all");
  
  thehistomanag->WriteMyHisto(LeptPt_mumumu_afterjetsel, "all");
  thehistomanag->WriteMyHisto(LeptPt_mumue_afterjetsel,  "all");
  thehistomanag->WriteMyHisto(LeptPt_eemu_afterjetsel,   "all");
  thehistomanag->WriteMyHisto(LeptPt_eee_afterjetsel,    "all");
  
  thehistomanag->WriteMyHisto(LeptPt_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(LeptPt_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(LeptPt_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(LeptPt_eee_afterbjetsel,    "all");
  
    
    
  thehistomanag->WriteMyHisto(LeptZPt_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(LeptZPt_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(LeptZPt_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(LeptZPt_eee_afterleptsel,    "all");
  
  thehistomanag->WriteMyHisto(LeptZPt_mumumu_afterjetsel, "all");
  thehistomanag->WriteMyHisto(LeptZPt_mumue_afterjetsel,  "all");
  thehistomanag->WriteMyHisto(LeptZPt_eemu_afterjetsel,   "all");
  thehistomanag->WriteMyHisto(LeptZPt_eee_afterjetsel,    "all");
  
  thehistomanag->WriteMyHisto(LeptZPt_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(LeptZPt_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(LeptZPt_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(LeptZPt_eee_afterbjetsel,    "all");
  
    
    
  thehistomanag->WriteMyHisto(LeptWPt_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(LeptWPt_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(LeptWPt_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(LeptWPt_eee_afterleptsel,    "all");
  
  thehistomanag->WriteMyHisto(LeptWPt_mumumu_afterjetsel, "all");
  thehistomanag->WriteMyHisto(LeptWPt_mumue_afterjetsel,  "all");
  thehistomanag->WriteMyHisto(LeptWPt_eemu_afterjetsel,   "all");
  thehistomanag->WriteMyHisto(LeptWPt_eee_afterjetsel,    "all");
  
  thehistomanag->WriteMyHisto(LeptWPt_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(LeptWPt_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(LeptWPt_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(LeptWPt_eee_afterbjetsel,    "all");
  
  thehistomanag->WriteMyHisto(LeptWPt_mumumu_afterbjetveto, "all");
  thehistomanag->WriteMyHisto(LeptWPt_mumue_afterbjetveto,  "all");
  thehistomanag->WriteMyHisto(LeptWPt_eemu_afterbjetveto,   "all");
  thehistomanag->WriteMyHisto(LeptWPt_eee_afterbjetveto,    "all");
  
    
    
  thehistomanag->WriteMyHisto(LeptWPt_mumumu_afterleptsel_mWT110, "all");
  thehistomanag->WriteMyHisto(LeptWPt_mumue_afterleptsel_mWT110,  "all");
  thehistomanag->WriteMyHisto(LeptWPt_eemu_afterleptsel_mWT110,   "all");
  thehistomanag->WriteMyHisto(LeptWPt_eee_afterleptsel_mWT110,    "all");
  
    
  thehistomanag->WriteMyHisto(JetPt_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(JetPt_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(JetPt_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(JetPt_eee_afterleptsel,    "all");
  
  thehistomanag->WriteMyHisto(JetPt_mumumu_afterjetsel, "all");
  thehistomanag->WriteMyHisto(JetPt_mumue_afterjetsel,  "all");
  thehistomanag->WriteMyHisto(JetPt_eemu_afterjetsel,   "all");
  thehistomanag->WriteMyHisto(JetPt_eee_afterjetsel,    "all");
  
  thehistomanag->WriteMyHisto(JetPt_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(JetPt_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(JetPt_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(JetPt_eee_afterbjetsel,    "all");
  
  thehistomanag->WriteMyHisto(JetPt_mumumu_afterbjetveto, "all");
  thehistomanag->WriteMyHisto(JetPt_mumue_afterbjetveto,  "all");
  thehistomanag->WriteMyHisto(JetPt_eemu_afterbjetveto,   "all");
  thehistomanag->WriteMyHisto(JetPt_eee_afterbjetveto,    "all");
  
    
    
  thehistomanag->WriteMyHisto(JetEta_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(JetEta_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(JetEta_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(JetEta_eee_afterleptsel,    "all");
  
  thehistomanag->WriteMyHisto(JetEta_mumumu_afterjetsel, "all");
  thehistomanag->WriteMyHisto(JetEta_mumue_afterjetsel,  "all");
  thehistomanag->WriteMyHisto(JetEta_eemu_afterjetsel,   "all");
  thehistomanag->WriteMyHisto(JetEta_eee_afterjetsel,    "all");
  
  thehistomanag->WriteMyHisto(JetEta_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(JetEta_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(JetEta_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(JetEta_eee_afterbjetsel,    "all");
  
  thehistomanag->WriteMyHisto(JetEta_mumumu_afterbjetveto, "all");
  thehistomanag->WriteMyHisto(JetEta_mumue_afterbjetveto,  "all");
  thehistomanag->WriteMyHisto(JetEta_eemu_afterbjetveto,   "all");
  thehistomanag->WriteMyHisto(JetEta_eee_afterbjetveto,    "all");
  
  
  thehistomanag->WriteMyHisto(HT_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(HT_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(HT_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(HT_eee_afterleptsel,    "all");
  
  
  thehistomanag->WriteMyHisto(HT_mumumu_afterjetsel, "all");
  thehistomanag->WriteMyHisto(HT_mumue_afterjetsel,  "all");
  thehistomanag->WriteMyHisto(HT_eemu_afterjetsel,   "all");
  thehistomanag->WriteMyHisto(HT_eee_afterjetsel,    "all");
  
  thehistomanag->WriteMyHisto(HT_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(HT_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(HT_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(HT_eee_afterbjetsel,    "all");
  
  thehistomanag->WriteMyHisto(HT_mumumu_afterbjetveto, "all");
  thehistomanag->WriteMyHisto(HT_mumue_afterbjetveto,  "all");
  thehistomanag->WriteMyHisto(HT_eemu_afterbjetveto,   "all");
  thehistomanag->WriteMyHisto(HT_eee_afterbjetveto,    "all");
  
  
  
  
  
  thehistomanag->WriteMyHisto(MET_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(MET_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(MET_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(MET_eee_afterleptsel,    "all");
  
  thehistomanag->WriteMyHisto(MET_mumumu_afterleptsel_mWT110, "all");
  thehistomanag->WriteMyHisto(MET_mumue_afterleptsel_mWT110,  "all");
  thehistomanag->WriteMyHisto(MET_eemu_afterleptsel_mWT110,   "all");
  thehistomanag->WriteMyHisto(MET_eee_afterleptsel_mWT110,    "all");
  
  
  thehistomanag->WriteMyHisto(MET_mumumu_afterjetsel, "all");
  thehistomanag->WriteMyHisto(MET_mumue_afterjetsel,  "all");
  thehistomanag->WriteMyHisto(MET_eemu_afterjetsel,   "all");
  thehistomanag->WriteMyHisto(MET_eee_afterjetsel,    "all");
  
  thehistomanag->WriteMyHisto(MET_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(MET_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(MET_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(MET_eee_afterbjetsel,    "all");
  
  thehistomanag->WriteMyHisto(MET_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(MET_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(MET_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(MET_eee_afterbjetsel,    "all");
  
  
  
  thehistomanag->WriteMyHisto(InvM_ll_mumumu_afterleptsel_highSumPt, "all");
  thehistomanag->WriteMyHisto(InvM_ll_mumue_afterleptsel_highSumPt, "all");
  thehistomanag->WriteMyHisto(InvM_ll_eemu_afterleptsel_highSumPt, "all");
  thehistomanag->WriteMyHisto(InvM_ll_eee_afterleptsel_highSumPt, "all");
  
  
  
  
  thehistomanag->WriteMyHisto(Asym_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(Asym_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(Asym_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(Asym_eee_afterbjetsel,    "all");
  
  
  
  thehistomanag->WriteMyHisto(RecoPtZ_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(RecoPtZ_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(RecoPtZ_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(RecoPtZ_eee_afterbjetsel,    "all");
  
  
  
  thehistomanag->WriteMyHisto(RecoPtZ_mumumu_afterbjetveto, "all");
  thehistomanag->WriteMyHisto(RecoPtZ_mumue_afterbjetveto,  "all");
  thehistomanag->WriteMyHisto(RecoPtZ_eemu_afterbjetveto,   "all");
  thehistomanag->WriteMyHisto(RecoPtZ_eee_afterbjetveto,    "all");
  
  thehistomanag->WriteMyHisto(RecoPtZ_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(RecoPtZ_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(RecoPtZ_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(RecoPtZ_eee_afterleptsel,    "all");
  
  
  thehistomanag->WriteMyHisto(RecoPtZ_mumumu_afterleptsel_nojet, "all");
  thehistomanag->WriteMyHisto(RecoPtZ_mumue_afterleptsel_nojet,  "all");
  thehistomanag->WriteMyHisto(RecoPtZ_eemu_afterleptsel_nojet,   "all");
  thehistomanag->WriteMyHisto(RecoPtZ_eee_afterleptsel_nojet,    "all");
  
  thehistomanag->WriteMyHisto(RecoTopMass_mumumu_afterbjetsel  , "all");
  thehistomanag->WriteMyHisto(RecoTopMass_mumue_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(RecoTopMass_eemu_afterbjetsel,    "all");
  thehistomanag->WriteMyHisto(RecoTopMass_eee_afterbjetsel,     "all");
  
  thehistomanag->WriteMyHisto(RecoTopMass_mumumu_afterbjetveto  , "all");
  thehistomanag->WriteMyHisto(RecoTopMass_mumue_afterbjetveto,   "all");
  thehistomanag->WriteMyHisto(RecoTopMass_eemu_afterbjetveto,    "all");
  thehistomanag->WriteMyHisto(RecoTopMass_eee_afterbjetveto,     "all");
  
  
  
  thehistomanag->WriteMyHisto(deltaPhilb_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(deltaPhilb_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(deltaPhilb_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(deltaPhilb_eee_afterbjetsel,    "all");
  
  thehistomanag->WriteMyHisto(deltaPhilj_mumumu_afterbjetveto, "all");
  thehistomanag->WriteMyHisto(deltaPhilj_mumue_afterbjetveto,  "all");
  thehistomanag->WriteMyHisto(deltaPhilj_eemu_afterbjetveto,   "all");
  thehistomanag->WriteMyHisto(deltaPhilj_eee_afterbjetveto,    "all");
  
  thehistomanag->WriteMyHisto(deltaR_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(deltaR_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(deltaR_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(deltaR_eee_afterleptsel,    "all");
  
    
  thehistomanag->WriteMyHisto(WmissAssing_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(WmissAssing_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(WmissAssing_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(WmissAssing_eee_afterleptsel,    "all");
  
  
  thehistomanag->WriteMyHisto(DijetInvM_mumumu_afterleptsel_inZpeak,    "all");
  thehistomanag->WriteMyHisto(DijetInvM_mumue_afterleptsel_inZpeak,     "all");
  thehistomanag->WriteMyHisto(DijetInvM_eemu_afterleptsel_inZpeak,      "all");
  thehistomanag->WriteMyHisto(DijetInvM_eee_afterleptsel_inZpeak,       "all");
  
  
  
  
  
  
  thehistomanag->WriteMyHisto(mWT_mumumu_afterbjetveto, "all");
  thehistomanag->WriteMyHisto(mWT_mumue_afterbjetveto,  "all");
  thehistomanag->WriteMyHisto(mWT_eemu_afterbjetveto,   "all");
  thehistomanag->WriteMyHisto(mWT_eee_afterbjetveto,    "all");

  thehistomanag->WriteMyHisto(mWT_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(mWT_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(mWT_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(mWT_eee_afterbjetsel,    "all");

  thehistomanag->WriteMyHisto(mWT_mumumu_afterjetsel, "all");
  thehistomanag->WriteMyHisto(mWT_mumue_afterjetsel,  "all");
  thehistomanag->WriteMyHisto(mWT_eemu_afterjetsel,   "all");
  thehistomanag->WriteMyHisto(mWT_eee_afterjetsel,    "all");

  thehistomanag->WriteMyHisto(mWT_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(mWT_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(mWT_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(mWT_eee_afterleptsel,    "all");

  thehistomanag->WriteMyHisto(cosThetaStar_mumumu_afterbjetsel, "all");
  thehistomanag->WriteMyHisto(cosThetaStar_mumue_afterbjetsel,  "all");
  thehistomanag->WriteMyHisto(cosThetaStar_eemu_afterbjetsel,   "all");
  thehistomanag->WriteMyHisto(cosThetaStar_eee_afterbjetsel,    "all");
  

  thehistomanag->WriteMyHisto(Charge_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(Charge_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(Charge_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(Charge_eee_afterleptsel,    "all");
  
  thehistomanag->WriteMyHisto(Charge_mumumu_afterleptsel_mWT110, "all");
  thehistomanag->WriteMyHisto(Charge_mumue_afterleptsel_mWT110,  "all");
  thehistomanag->WriteMyHisto(Charge_eemu_afterleptsel_mWT110,   "all");
  thehistomanag->WriteMyHisto(Charge_eee_afterleptsel_mWT110,    "all");
  
  
  
  
  
  thehistomanag->WriteMyHisto(deltaRLeptJet_mumumu_afterleptsel_mWT110, "all");
  thehistomanag->WriteMyHisto(deltaRLeptJet_mumue_afterleptsel_mWT110,  "all");
  thehistomanag->WriteMyHisto(deltaRLeptJet_eemu_afterleptsel_mWT110,   "all");
  thehistomanag->WriteMyHisto(deltaRLeptJet_eee_afterleptsel_mWT110,    "all");
  
  thehistomanag->WriteMyHisto(deltaRLeptMet_mumumu_afterleptsel_mWT110, "all");
  thehistomanag->WriteMyHisto(deltaRLeptMet_mumue_afterleptsel_mWT110,  "all");
  thehistomanag->WriteMyHisto(deltaRLeptMet_eemu_afterleptsel_mWT110,   "all");
  thehistomanag->WriteMyHisto(deltaRLeptMet_eee_afterleptsel_mWT110,    "all");
  
  
  thehistomanag->WriteMyHisto(NJet_mumumu_afterleptsel_mWT110, "all");
  thehistomanag->WriteMyHisto(NJet_mumue_afterleptsel_mWT110,  "all");
  thehistomanag->WriteMyHisto(NJet_eemu_afterleptsel_mWT110,   "all");
  thehistomanag->WriteMyHisto(NJet_eee_afterleptsel_mWT110,    "all");
  
  thehistomanag->WriteMyHisto(NJet_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(NJet_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(NJet_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(NJet_eee_afterleptsel,    "all");
  
  
  
  
  thehistomanag->WriteMyHisto(NBJet_mumumu_afterleptsel_mWT110, "all");
  thehistomanag->WriteMyHisto(NBJet_mumue_afterleptsel_mWT110,  "all");
  thehistomanag->WriteMyHisto(NBJet_eemu_afterleptsel_mWT110,   "all");
  thehistomanag->WriteMyHisto(NBJet_eee_afterleptsel_mWT110,    "all");
  
  thehistomanag->WriteMyHisto(NBJet_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto(NBJet_mumue_afterleptsel,  "all");
  thehistomanag->WriteMyHisto(NBJet_eemu_afterleptsel,   "all");
  thehistomanag->WriteMyHisto(NBJet_eee_afterleptsel,    "all");
  
  thehistomanag->WriteMyHisto(NBJet_mumumu_afterjetsel_bjets, "all");
  thehistomanag->WriteMyHisto(NBJet_mumue_afterjetsel_bjets,  "all"); 
  thehistomanag->WriteMyHisto(NBJet_eemu_afterjetsel_bjets,   "all");  
  thehistomanag->WriteMyHisto(NBJet_eee_afterjetsel_bjets,    "all");   


  thehistomanag->WriteMyHisto(NBJet_mumumu_afterjetsel_cjets, "all");
  thehistomanag->WriteMyHisto(NBJet_mumue_afterjetsel_cjets,  "all"); 
  thehistomanag->WriteMyHisto(NBJet_eemu_afterjetsel_cjets,   "all");  
  thehistomanag->WriteMyHisto(NBJet_eee_afterjetsel_cjets,    "all");   


  thehistomanag->WriteMyHisto(NBJet_mumumu_afterjetsel_ljets, "all");
  thehistomanag->WriteMyHisto(NBJet_mumue_afterjetsel_ljets,  "all"); 
  thehistomanag->WriteMyHisto(NBJet_eemu_afterjetsel_ljets,   "all");  
  thehistomanag->WriteMyHisto(NBJet_eee_afterjetsel_ljets ,   "all");  

  
  
  thehistomanag->WriteMyHisto(BJetDiscri_mumumu_afterjetsel_bjets,    "all");
  thehistomanag->WriteMyHisto(BJetDiscri_mumue_afterjetsel_bjets,     "all");
  thehistomanag->WriteMyHisto(BJetDiscri_eemu_afterjetsel_bjets,      "all");
  thehistomanag->WriteMyHisto(BJetDiscri_eee_afterjetsel_bjets,       "all");
  
  
  
  thehistomanag->WriteMyHisto(BJetDiscri_mumumu_afterjetsel_cjets,    "all");
  thehistomanag->WriteMyHisto(BJetDiscri_mumue_afterjetsel_cjets,     "all");
  thehistomanag->WriteMyHisto(BJetDiscri_eemu_afterjetsel_cjets,      "all");
  thehistomanag->WriteMyHisto(BJetDiscri_eee_afterjetsel_cjets,       "all");
  
  
  thehistomanag->WriteMyHisto(BJetDiscri_mumumu_afterjetsel_ljets,    "all");
  thehistomanag->WriteMyHisto(BJetDiscri_mumue_afterjetsel_ljets,     "all");
  thehistomanag->WriteMyHisto(BJetDiscri_eemu_afterjetsel_ljets,      "all");
  thehistomanag->WriteMyHisto(BJetDiscri_eee_afterjetsel_ljets,       "all");
  
  
  
  
  thehistomanag->WriteMyHisto2D(HT_vs_MET_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_MET_mumue_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_MET_eemu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_MET_eee_afterleptsel, "all");
 
  thehistomanag->WriteMyHisto2D(HT_vs_NJet_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_NJet_mumue_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_NJet_eemu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_NJet_eee_afterleptsel, "all");
 
  thehistomanag->WriteMyHisto2D(HT_vs_NBJet_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_NBJet_mumue_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_NBJet_eemu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_NBJet_eee_afterleptsel, "all");
 
  thehistomanag->WriteMyHisto2D(HT_vs_LeptPt_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_LeptPt_mumue_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_LeptPt_eemu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_LeptPt_eee_afterleptsel, "all");
 
  thehistomanag->WriteMyHisto2D(HT_vs_JetPt_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_JetPt_mumue_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_JetPt_eemu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_JetPt_eee_afterleptsel, "all");
  
  
  
  thehistomanag->WriteMyHisto2D(HT_vs_Mll_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_Mll_mumue_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_Mll_eemu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(HT_vs_Mll_eee_afterleptsel, "all");
  
  thehistomanag->WriteMyHisto2D(InvM_ll_vs_mWT_mumumu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(InvM_ll_vs_mWT_mumue_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(InvM_ll_vs_mWT_eemu_afterleptsel, "all");
  thehistomanag->WriteMyHisto2D(InvM_ll_vs_mWT_eee_afterleptsel, "all");


}



void ProofSelectorMyCutFlow::cleanHistoVector(){
  
  
  
  
  SelABjet.clear();
  SelABjet_afterjetsel.clear();
  
  Ntrilept_mumumu.clear();
  Ntrileptnoniso_mumumu.clear();
  
  
  CutFlow_mumumu.clear();
  CutFlow_mumue.clear();
  CutFlow_eemu.clear();
  CutFlow_eee.clear();
  
  
  ErrCutFlow_mumumu.clear();
  ErrCutFlow_mumue.clear();
  ErrCutFlow_eemu.clear();
  ErrCutFlow_eee.clear();
  
  PU_before_mumumu.clear();
  PU_before_mumue.clear();
  PU_before_eemu.clear();
  PU_before_eee.clear();

  PU_intime_mumumu.clear();
  PU_intime_mumue.clear();
  PU_intime_eemu.clear();
  PU_intime_eee.clear();
  
  PU_after_mumumu.clear();
  PU_after_mumue.clear();
  PU_after_eemu.clear();
  PU_after_eee.clear(); 
  
  
  NVtx_mumumu.clear();
  NVtx_mumue.clear();
  NVtx_eemu.clear();
  NVtx_eee.clear(); 
  
  NVtx_mumumu_aftertrigsel.clear();
  NVtx_mumue_aftertrigsel.clear();
  NVtx_eemu_aftertrigsel.clear();
  NVtx_eee_aftertrigsel.clear(); 
  
  NVtx_mumumu_afterleptsel.clear();
  NVtx_mumue_afterleptsel.clear();
  NVtx_eemu_afterleptsel.clear();
  NVtx_eee_afterleptsel.clear(); 
  
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
  
  
  
  
  NBJet_mumumu_afterjetsel_bjets.clear();
  NBJet_mumue_afterjetsel_bjets .clear();
  NBJet_eemu_afterjetsel_bjets  .clear();
  NBJet_eee_afterjetsel_bjets  .clear();
  
  
  NBJet_mumumu_afterjetsel_cjets.clear();
  NBJet_mumue_afterjetsel_cjets .clear();
  NBJet_eemu_afterjetsel_cjets  .clear();
  NBJet_eee_afterjetsel_cjets  .clear();
  
  
  NBJet_mumumu_afterjetsel_ljets.clear();
  NBJet_mumue_afterjetsel_ljets .clear();
  NBJet_eemu_afterjetsel_ljets  .clear();
  NBJet_eee_afterjetsel_ljets  .clear();
  
  
  
  BJetDiscri_mumumu_afterjetsel_bjets.clear();
  BJetDiscri_mumue_afterjetsel_bjets.clear();
  BJetDiscri_eemu_afterjetsel_bjets.clear();
  BJetDiscri_eee_afterjetsel_bjets.clear();
  
  BJetDiscri_mumumu_afterjetsel_cjets.clear();
  BJetDiscri_mumue_afterjetsel_cjets.clear();
  BJetDiscri_eemu_afterjetsel_cjets.clear();
  BJetDiscri_eee_afterjetsel_cjets.clear();
  
  BJetDiscri_mumumu_afterjetsel_ljets.clear();
  BJetDiscri_mumue_afterjetsel_ljets.clear();
  BJetDiscri_eemu_afterjetsel_ljets.clear();
  BJetDiscri_eee_afterjetsel_ljets.clear();
  
  NBJet_mumumu_afterleptsel_mWT110.clear();
  NBJet_mumue_afterleptsel_mWT110.clear();
  NBJet_eemu_afterleptsel_mWT110.clear();
  NBJet_eee_afterleptsel_mWT110.clear();
  
  
  NBJet_mumumu_afterleptsel.clear();
  NBJet_mumue_afterleptsel.clear();
  NBJet_eemu_afterleptsel.clear();
  NBJet_eee_afterleptsel.clear();
  
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
  
  
  cosThetaStar_mumumu_afterbjetsel.clear();
  cosThetaStar_mumue_afterbjetsel.clear();
  cosThetaStar_eemu_afterbjetsel.clear();
  cosThetaStar_eee_afterbjetsel.clear();
  
  
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
  
  setNullHistoPointer();
}






bool ProofSelectorMyCutFlow::selectedEventForSynch(int event_nbr){


  if( event_nbr == 49782|| event_nbr ==49787|| event_nbr ==49787|| event_nbr ==49793|| event_nbr ==49796|| event_nbr ==49801|| 
  event_nbr ==49805|| event_nbr ==49811|| event_nbr ==49819|| event_nbr ==49819|| event_nbr ==49824|| event_nbr ==49826|| 
  event_nbr ==49826|| event_nbr ==49828|| event_nbr ==49828|| event_nbr ==49832|| event_nbr ==49836|| event_nbr ==49836|| 
  event_nbr ==49846|| event_nbr ==49849|| event_nbr ==49861|| event_nbr ==49867|| event_nbr ==49867|| event_nbr ==49872|| 
  event_nbr ==49872|| event_nbr ==49876|| event_nbr ==49882|| event_nbr ==49883|| event_nbr ==49885|| event_nbr ==49886|| 
  event_nbr ==49889|| event_nbr ==49893|| event_nbr ==49894|| event_nbr ==49904|| event_nbr ==49910|| event_nbr ==49929|| 
  event_nbr ==49929|| event_nbr ==49933|| event_nbr ==49933|| event_nbr ==49941|| event_nbr ==49944|| event_nbr ==49947|| 
  event_nbr ==49949|| event_nbr ==49951|| event_nbr ==49951|| event_nbr ==49953|| event_nbr ==49964|| event_nbr ==49964|| 
  event_nbr ==49967|| event_nbr ==49967|| event_nbr ==49968|| event_nbr ==49971|| event_nbr ==49973|| event_nbr ==49979|| 
  event_nbr ==49989|| event_nbr ==49999|| event_nbr ==50009|| event_nbr ==50012|| event_nbr ==50018|| event_nbr ==50023|| 
  event_nbr ==50023|| event_nbr ==50024|| event_nbr ==50024|| event_nbr ==50027|| event_nbr ==50028|| event_nbr ==50030|| 
  event_nbr ==50031|| event_nbr ==50032|| event_nbr ==50032|| event_nbr ==50042|| event_nbr ==50052|| event_nbr ==50056|| 
  event_nbr ==50056|| event_nbr ==50066|| event_nbr ==50067|| event_nbr ==50067|| event_nbr ==50070|| event_nbr ==50075|| 
  event_nbr ==50077|| event_nbr ==50082|| event_nbr ==50082|| event_nbr ==50084|| event_nbr ==50084|| event_nbr ==50086|| 
  event_nbr ==50089|| event_nbr ==50092|| event_nbr ==50098|| event_nbr ==50106|| event_nbr ==50115|| event_nbr ==50117|| 
  event_nbr ==50117|| event_nbr ==50122|| event_nbr ==50131|| event_nbr ==50131|| event_nbr ==50138|| event_nbr ==50142|| 
  event_nbr ==50142|| event_nbr ==50149|| event_nbr ==50162|| event_nbr ==50168|| event_nbr ==50178) return true;
  else return false;


}


bool ProofSelectorMyCutFlow::selectedEventForSynch_step0(int event_nbr){ 
  
  if( 
    event_nbr == 153483||
    event_nbr == 153496||
    event_nbr == 153529||
    event_nbr == 153544||
    event_nbr == 153553||
    event_nbr == 153579||
    event_nbr == 153583||
    event_nbr == 153513||
    event_nbr == 153594||
    event_nbr == 153606||
    event_nbr == 153513||
    event_nbr == 153594||
    event_nbr == 153606||
    event_nbr == 153550 ) return true;
    else return false;
  
  
}



double ProofSelectorMyCutFlow::getLeptonSF(TLorentzVector lept1, TLorentzVector lept2, TLorentzVector lept3, TString channel, int syst){
  
  double theSF = 0;
  
    
  if(channel == "mumumu"){
    
    std::vector<double> SF_muID1  = sel.getScaleFactorMuID(   lept1.Pt(), lept1.Eta());
    std::vector<double> SF_muIso1 = sel.getScaleFactorMuIso20(lept1.Pt(), lept1.Eta());
    std::vector<double> SF_muID2  = sel.getScaleFactorMuID(   lept2.Pt(), lept2.Eta());
    std::vector<double> SF_muIso2 = sel.getScaleFactorMuIso20(lept2.Pt(), lept2.Eta());
    std::vector<double> SF_muID3  = sel.getScaleFactorMuID(   lept3.Pt(), lept3.Eta());
    std::vector<double> SF_muIso3 = sel.getScaleFactorMuIso20(lept3.Pt(), lept3.Eta());
    
    if(syst == 0) theSF  =  SF_muID1[0]*SF_muIso1[0]*SF_muID2[0]*SF_muIso2[0]*SF_muID3[0]*SF_muIso3[0];    
    if(syst == 1) theSF  = (SF_muID1[0]+SF_muID1[1])*(SF_muIso1[0]+SF_muIso1[1])*
    			   (SF_muID2[0]+SF_muID2[1])*(SF_muIso2[0]+SF_muIso2[1])*
			   (SF_muID3[0]+SF_muID3[1])*(SF_muIso3[0]+SF_muIso3[1]);			   
    if(syst == -1)theSF  = (SF_muID1[0]-SF_muID1[2])*(SF_muIso1[0]-SF_muIso1[2])*
    			   (SF_muID2[0]-SF_muID2[2])*(SF_muIso2[0]-SF_muIso2[2])*
		           (SF_muID3[0]-SF_muID3[2])*(SF_muIso3[0]-SF_muIso3[2]);	   
    
  }
  else if(channel == "mumue"){
    
    std::vector<double> SF_muID1  = sel.getScaleFactorMuID(   lept1.Pt(), lept1.Eta());
    std::vector<double> SF_muIso1 = sel.getScaleFactorMuIso20(lept1.Pt(), lept1.Eta());
    std::vector<double> SF_muID2  = sel.getScaleFactorMuID(   lept2.Pt(), lept2.Eta());
    std::vector<double> SF_muIso2 = sel.getScaleFactorMuIso20(lept2.Pt(), lept2.Eta());
    std::vector<double> SF_elAll3 = sel.getSscaleFactorElectronAllID05(lept3.Pt(), lept3.Eta());
    
    if(syst == 0) theSF  =  SF_muID1[0]*SF_muIso1[0]*SF_muID2[0]*SF_muIso2[0]*SF_elAll3[0];    
    if(syst == 1) theSF  = (SF_muID1[0]+SF_muID1[1])*(SF_muIso1[0]+SF_muIso1[1])*
    			   (SF_muID2[0]+SF_muID2[1])*(SF_muIso2[0]+SF_muIso2[1])*
			   (SF_elAll3[0]+SF_elAll3[1]);			   
    if(syst == -1)theSF  = (SF_muID1[0]-SF_muID1[2])*(SF_muIso1[0]-SF_muIso1[2])*
    			   (SF_muID2[0]-SF_muID2[2])*(SF_muIso2[0]-SF_muIso2[2])*
			   (SF_elAll3[0]-SF_elAll3[2]);
    
  }
  else if(channel == "eemu"){
    
    std::vector<double> SF_elAll1 = sel.getSscaleFactorElectronAllID05(lept1.Pt(), lept1.Eta());
    std::vector<double> SF_elAll2 = sel.getSscaleFactorElectronAllID05(lept2.Pt(), lept2.Eta());
    std::vector<double> SF_muID3  = sel.getScaleFactorMuID(lept3.Pt(), lept3.Eta());
    std::vector<double> SF_muIso3 = sel.getScaleFactorMuIso20(lept3.Pt(), lept3.Eta());
    
    
    if(syst == 0) theSF  =  SF_elAll1[0]*SF_elAll2[0]*SF_muID3[0]*SF_muIso3[0];    
    if(syst == 1) theSF  = (SF_elAll1[0]+SF_elAll1[1])*
			   (SF_elAll2[0]+SF_elAll2[1])*
			   (SF_muID3[0]+SF_muID3[1])*(SF_muIso3[0]+SF_muIso3[1]);			   
    if(syst == -1)theSF  = (SF_elAll1[0]-SF_elAll1[2])*
			   (SF_elAll2[0]-SF_elAll2[2])*
			   (SF_muID3[0]-SF_muID3[2])*(SF_muIso3[0]+SF_muIso3[2]) ;   
    
  }
  else if(channel == "eee"){
    std::vector<double> SF_elAll1 = sel.getSscaleFactorElectronAllID05(lept1.Pt(), lept1.Eta());    
    std::vector<double> SF_elAll2 = sel.getSscaleFactorElectronAllID05(lept2.Pt(), lept2.Eta()); 
    std::vector<double> SF_elAll3 = sel.getSscaleFactorElectronAllID05(lept3.Pt(), lept3.Eta());
    
    
    if(syst == 0) theSF  =  SF_elAll1[0]*SF_elAll2[0]*SF_elAll3[0];    
    if(syst == 1) theSF  = (SF_elAll1[0]+SF_elAll1[1])*
			   (SF_elAll2[0]+SF_elAll2[1])*
			   (SF_elAll3[0]+SF_elAll3[1]);			   
    if(syst == -1)theSF  = (SF_elAll1[0]-SF_elAll1[2])*
			   (SF_elAll2[0]-SF_elAll2[2])*
			   (SF_elAll3[0]-SF_elAll3[2]);
			       
    
    cout << "nom value electron " << SF_elAll1[0]  << " + " << SF_elAll1[1]  << " " << SF_elAll1[2]  << endl;
    
  }
    
  
  return theSF;


}


double ProofSelectorMyCutFlow::getTriggerSF(int i){

  return 0;
}




NTMET ProofSelectorMyCutFlow::SmearedMET(vector<NTJet> &selJets, vector<NTJet> &selJetsNoSmear, NTMET &met){

  IPHCTree::NTMET newMET;

    
  TVector2  themet(met.p2.Px(),met.p2.Py() );
  
    
  for (unsigned int i=0; i<selJetsNoSmear.size(); i++)
  {
    TVector2 jetDivector( (selJetsNoSmear[i]).p4.Px(), (selJetsNoSmear[i]).p4.Py());
    themet = themet + jetDivector;
  } 
  
  //cout << "---------------" << endl;
  
  for (unsigned int i=0; i<selJets.size(); i++)
  {
    TVector2 jetDivector((selJets[i]).p4.Px(), (selJets[i]).p4.Py());
    themet = themet - jetDivector;
  } 
  
    
  
  newMET = met;
  newMET.p2.Set(themet.Px(),
                themet.Py());
  
  
    

  return newMET;

}










