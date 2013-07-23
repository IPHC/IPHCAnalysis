#define ProofSelectorMatrixMethod_cxx

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

#include "ProofSelectorMatrixMethod.h"

//_____________________________________________________________________________
ProofSelectorMatrixMethod::ProofSelectorMatrixMethod()
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
  
  
  applyDYscale = true ;
  applyFakescale = true  ;
  
  IReweight		= true;
  IDYestimateWithMetCut = true;
  IReweight_puUp	= false;
  IReweight_puDown	= false;
  
  
  useNonIsoWcand = false;
  
}

//_____________________________________________________________________________
ProofSelectorMatrixMethod::~ProofSelectorMatrixMethod()
{
  // Destructor
  
  //SafeDelete(fHist);
}

//_____________________________________________________________________________
void ProofSelectorMatrixMethod::Init(TTree *tree)
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
void ProofSelectorMatrixMethod::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  cout << "start Begin " << endl;
  TString option = GetOption();
  cout << "end  Begin" << endl;
  
  
}

//_____________________________________________________________________________
void ProofSelectorMatrixMethod::SlaveBegin(TTree * tree)
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
  anaEL->LoadWeight (sel); // now the parameters for SFBweight are initialized (for b-tag!)
  
  //Load for PU:
  sel.GeneratePUWeight(PUWeightFileName);
  
  //******************************************
  //Load Scale Factors for lepton efficiencies
  //********d**********************************
  sel.LoadElScaleFactors();
  sel.LoadMuScaleFactors();
  
   //--------------------------------------//
   //   Fill cuts and channels  	
   //--------------------------------------//
  CutName.push_back("Cut1");
   
  //--------------------------------------//
  //   Initializing variables 	
  //--------------------------------------//
  
    
  
  
  SF_trig_mu = 0.977; 
  SF_trig_emu= 1.008; 
  SF_trig_ee = 0.962; 
  
  SF_trig_mu_error  = 0.014; 
  SF_trig_emu_error = 0.007; 
  SF_trig_ee_error  = 0.016; 
  
  SF_e  =1.0; 
  SF_mu =1.0; 
  
  /*SF_met_mumu = 1.00237;
    SF_met_emu  = 1.0;
    SF_met_ee   = 1.00237;
    
    
    SF_met_mumu_error = 0.00629333;
    SF_met_emu_error  = 0;
    SF_met_ee_error   = 0.00629333;*/
  
  
  //**************************************
  //**************************************
  //******* DY DD estimate ***************
  //**************************************
  //**************************************
  
  
  // for PLR plot
  SF_DY_ee = 1.;
  SF_DY_mm = 1.;
  SF_DY_em = 1.;
  
  double lumisf = 1./1.037;
  
  if(applyDYscale){ 
    
    
    //after dilept invMass
    vSF_DY_ee.push_back(1) ;
    vSF_DY_ee_error.push_back(0) ;
    //after jet cut
    vSF_DY_ee.push_back(1.169*lumisf) ;
    vSF_DY_ee_error.push_back(0.077*lumisf) ;
    //after met cut
    vSF_DY_ee.push_back( 1.9281*lumisf/vSF_DY_ee[1]) ;
    vSF_DY_ee_error.push_back(0.404145*lumisf) ;
    //for 1 btag 
    vSF_DY_ee.push_back(2.11816*lumisf/(1.9281*lumisf)) ;
    vSF_DY_ee_error.push_back(0.53244*lumisf) ;
    //for E btag 
    vSF_DY_ee.push_back(2.1314*lumisf/(1.9281*lumisf)) ;
    vSF_DY_ee_error.push_back(1.06501*lumisf) ;
    
    //after dilept invMass
    vSF_DY_mm.push_back(1) ;
    vSF_DY_mm_error.push_back(0) ;
    //after jet cut
    vSF_DY_mm.push_back(1.257*lumisf) ;
    vSF_DY_mm_error.push_back(0.078*lumisf) ;
    //after met cut
    vSF_DY_mm.push_back(1.82219*lumisf/vSF_DY_mm[1]) ;
    vSF_DY_mm_error.push_back(0.365828*lumisf) ;
    //for 1 btag 
    vSF_DY_mm.push_back(1.73089*lumisf/(1.82219*lumisf)) ;
    vSF_DY_mm_error.push_back(0.428215*lumisf) ;
    //for E btag 
    vSF_DY_mm.push_back(0.819327*lumisf/(1.82219*lumisf)) ;
    vSF_DY_mm_error.push_back(0.667288*lumisf) ;
    
    
    
    
    //after dilept invMass
    vSF_DY_em.push_back(1) ;
    vSF_DY_em_error.push_back(0) ;
    //after jet cut
    vSF_DY_em.push_back(1.38872*lumisf) ;
    vSF_DY_em_error.push_back(0.296232*lumisf) ;
    //after met cut
    vSF_DY_em.push_back(1) ;
    vSF_DY_em_error.push_back(0.296232*lumisf) ;
    //for 1 btag 
    vSF_DY_em.push_back(1) ;
    vSF_DY_em_error.push_back(0.296232*lumisf) ;
    //for E btag 
    vSF_DY_em.push_back(1) ;
    vSF_DY_em_error.push_back(0.296232*lumisf) ;
    
    
    // for PLR plot
    SF_DY_ee = 1.9281*lumisf;
    SF_DY_mm = 1.82219*lumisf;
    SF_DY_em = 1.38872*lumisf;
    
    
  }else{
    
    
    for(unsigned int i=0; i< 5; i++){
      vSF_DY_em.push_back(1.) ;
      vSF_DY_em_error.push_back(0.) ;
      vSF_DY_ee.push_back(1.) ;
      vSF_DY_ee_error.push_back(0.) ;
      vSF_DY_mm.push_back(1.) ;
      vSF_DY_mm_error.push_back(0.) ;
    }
    
    
    
  }
  SF_BranchingRatio_ll = (0.108*9.)*(0.108*9.);
  SF_BranchingRatio_lj = (0.108*9.)*(0.676*1.5);
  SF_BranchingRatio_had = (0.676*1.5)*(0.676*1.5);
  
  
  
  
  
  //**************************************
  //**************************************
  //******* fakes DD estimate ************
  //**************************************
  //**************************************
  
  
  if(applyFakescale){  
    //after MET selection
    vSF_FakeBck_ee.push_back(0.966*lumisf); 
    vSF_FakeBck_ee_error.push_back(0.674*lumisf); 
    //after 1 b jet
    vSF_FakeBck_ee.push_back(0.35*lumisf/vSF_FakeBck_ee[0]); 
    vSF_FakeBck_ee_error.push_back(0.615*lumisf); 
    //after 2 b jet
    vSF_FakeBck_ee.push_back(0.35*lumisf/vSF_FakeBck_ee[0]); 
    vSF_FakeBck_ee_error.push_back(0.615*lumisf); 
    
    
    
    //after MET selection
    vSF_FakeBck_mm.push_back(4.9*lumisf); 
    vSF_FakeBck_mm_error.push_back(2.3*lumisf); 
   //after 1 b jet
    vSF_FakeBck_mm.push_back(4.0*lumisf/vSF_FakeBck_mm[0]); 
    vSF_FakeBck_mm_error.push_back(2.2*lumisf); 
    //after 2 b jet
    vSF_FakeBck_mm.push_back(4.0*lumisf/vSF_FakeBck_mm[0]); 
    vSF_FakeBck_mm_error.push_back(2.2*lumisf);
    
   
    
    //after MET selection
    vSF_FakeBck_em.push_back(3.8*lumisf); 
    vSF_FakeBck_em_error.push_back(0.9*lumisf); 
    //after 1 b jet
    vSF_FakeBck_em.push_back(2.5*lumisf/vSF_FakeBck_em[0]); 
    vSF_FakeBck_em_error.push_back(0.25*lumisf); 
    //after 2 b jet
    vSF_FakeBck_em.push_back(2.5*lumisf/vSF_FakeBck_em[0]); 
    vSF_FakeBck_em_error.push_back(0.25*lumisf); 
    
    
    
    
    
    SF_Wjets_ee = 0.966*lumisf ;
    SF_Wjets_mm = 4.9*lumisf ;
    SF_Wjets_em = 3.8*lumisf;
    
    SF_QCD_ee = 1.;
    SF_QCD_mm = 1.;
    SF_QCD_em = 1.;
    
  }else{
    
    
    for(unsigned int i=0; i< 5; i++){
      vSF_FakeBck_em.push_back(1.) ;
      vSF_FakeBck_em_error.push_back(0.) ;
      vSF_FakeBck_ee.push_back(1.) ;
      vSF_FakeBck_ee_error.push_back(0.) ;
      vSF_FakeBck_mm.push_back(1.) ;
      vSF_FakeBck_mm_error.push_back(0.) ;
    }
  }
  
  
  
  //**************************************
  //**************************************
  //******* MET DD estimate **************
  //**************************************
  //**************************************
  
  //SF_met_mumu = 1.00237;
  //SF_met_emu  = 1.0;
  //SF_met_ee   = 1.00237;
  //SF_met_mumu_error = 0.00629333;
  //SF_met_emu_error  = 0;
  //SF_met_ee_error   = 0.00629333;
   
  
  vSF_met_ee.push_back(1.0029);
  vSF_met_ee_error.push_back(0.0167596);
  vSF_met_ee.push_back(1.00839/vSF_met_ee[0]);
  vSF_met_ee_error.push_back(0.0159795);
  vSF_met_ee.push_back(1.01456/vSF_met_ee[0]);
  vSF_met_ee_error.push_back(0.017556);
  
  vSF_met_mumu.push_back(1.0029);
  vSF_met_mumu_error.push_back(0.0135555);
  vSF_met_mumu.push_back(1.00839/vSF_met_mumu[0]);
  vSF_met_mumu_error.push_back(0.012474);
  vSF_met_mumu.push_back(1.01456/vSF_met_mumu[0]);
  vSF_met_mumu_error.push_back(0.0179107);
  
  vSF_met_emu.push_back(1.);
  vSF_met_emu_error.push_back(0.);
  vSF_met_emu.push_back(1.);
  vSF_met_emu_error.push_back(0.);
  vSF_met_emu.push_back(1.);
  vSF_met_emu_error.push_back(0.);
  
  sumSFlept_ee    = 0;
  sumSFlept_mumu  = 0;
  sumSFlept_emu   = 0;
  
  nEvents_ee   = 0;
  nEvents_mumu   = 0;
  nEvents_emu   = 0;
  
  
  scaleElec = 1.0; // 1 to switch off
  resolElec = 0.0; // 0 to switch off
  
  ApplyLeptonSF = true;
  
  ITypeMC     = -1;
  ICut        = -1;  
  
  
  //************************************
  //For trigger systematics 
  
  if(datasetName=="TTbarTriggerUp"){
    SF_trig_mu += SF_trig_mu_error;
    SF_trig_emu+= SF_trig_emu_error;  
    SF_trig_ee += SF_trig_ee_error; 
  } 
  if(datasetName=="TTbarTriggerDown"){
    SF_trig_mu  -= SF_trig_mu_error;
    SF_trig_emu -= SF_trig_emu_error;  
    SF_trig_ee  -= SF_trig_ee_error; 
  } 
  
  //************************************
  
  //************************************
  //For MET systematics 
  
  
  /*if(datasetName=="TTbarMETUp"){
    SF_met_mumu  += SF_met_mumu_error;
    SF_met_emu += SF_met_emu_error;  
    SF_met_ee  += SF_met_ee_error; 
    } 
    if(datasetName=="TTbarMETDown"){
    SF_met_mumu  -= SF_met_mumu_error;
     SF_met_emu -= SF_met_emu_error;  
     SF_met_ee  -= SF_met_ee_error; 
     } */
  
  
  
  for(unsigned int d=0;d<datasets.size();d++){
    cout << "datasets.size() " << datasets.size()<< "  datasets[d].Name()" << datasets[d].Name()  << " datasetName "
	 <<datasetName  << endl;
    if(datasets[d].Name()==datasetName)dataset = &datasets[d];
  }
  
  
  
  
  
  
  //--------------------------------------//
  //   Managing histos  	
  //--------------------------------------//
  MyhistoManager.LoadDatasets(datasets);   
  MyhistoManager.LoadSelectionSteps(CutName);
  MyhistoManager.LoadChannels(ChannelName);
  //example
  
  nbins = 200;
  minx = 0.;
  maxx = 350;
  
  
  //***********************
  // initiate lumi reweighting
  
  
  
  MyhistoManager.CreateHisto(CutFlow_mumumu_tight,  "CutFlow_mumumu_tight" ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_mumue_tight,   "CutFlow_mumue_tight"  ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_eemu_tight,    "CutFlow_eemu_tight"   ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_eee_tight,     "CutFlow_eee_tight"    ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  
  //MyhistoManager.SetCutFlowAxisTitleFCNCMonotop(CutFlow_mumumu_tight,   "CutFlow_mumumu_tight"  ,datasetName);
  //MyhistoManager.SetCutFlowAxisTitleFCNCMonotop(CutFlow_mumue_tight,    "CutFlow_mumue_tight"   ,datasetName);
  //MyhistoManager.SetCutFlowAxisTitleFCNCMonotop(CutFlow_eemu_tight,     "CutFlow_eemu_tight"    ,datasetName);
  //MyhistoManager.SetCutFlowAxisTitleFCNCMonotop(CutFlow_eee_tight,      "CutFlow_eee_tight"     ,datasetName);
  
  MyhistoManager.CreateHisto(ErrCutFlow_mumumu_tight,  "ErrCutFlow_mumumu_tight"  ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_mumue_tight,   "ErrCutFlow_mumue_tight"   ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_eemu_tight,    "ErrCutFlow_eemu_tight"    ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_eee_tight,     "ErrCutFlow_eee_tight"     ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);

  MyhistoManager.CreateHisto(CutFlow_mumumu_loose,  "CutFlow_mumumu_loose" ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_mumue_loose,   "CutFlow_mumue_loose"  ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_eemu_loose,    "CutFlow_eemu_loose"   ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_eee_loose,     "CutFlow_eee_loose"    ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  
  MyhistoManager.CreateHisto(ErrCutFlow_mumumu_loose,  "ErrCutFlow_mumumu_loose"  ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_mumue_loose,   "ErrCutFlow_mumue_loose"   ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_eemu_loose,    "ErrCutFlow_eemu_loose"    ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_eee_loose,     "ErrCutFlow_eee_loose"     ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);

  
  MyhistoManager.CreateHisto(ThirdLepPt_mumumu_tight_cut1,  "ThirdLepPt_mumumu_tight_cut1" ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_mumue_tight_cut1,   "ThirdLepPt_mumue_tight_cut1"  ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_eemu_tight_cut1,    "ThirdLepPt_eemu_tight_cut1"   ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_eee_tight_cut1,     "ThirdLepPt_eee_tight_cut1"    ,datasetName,"ThirdLepPt","Entries",100,0,200);
  
  MyhistoManager.CreateHisto(ThirdLepPt_mumumu_tight_cut2,  "ThirdLepPt_mumumu_tight_cut2" ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_mumue_tight_cut2,   "ThirdLepPt_mumue_tight_cut2"  ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_eemu_tight_cut2,    "ThirdLepPt_eemu_tight_cut2"   ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_eee_tight_cut2,     "ThirdLepPt_eee_tight_cut2"    ,datasetName,"ThirdLepPt","Entries",100,0,200);
  
  MyhistoManager.CreateHisto(ThirdLepPt_mumumu_tight_cut3,  "ThirdLepPt_mumumu_tight_cut3" ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_mumue_tight_cut3,   "ThirdLepPt_mumue_tight_cut3"  ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_eemu_tight_cut3,    "ThirdLepPt_eemu_tight_cut3"   ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_eee_tight_cut3,     "ThirdLepPt_eee_tight_cut3"    ,datasetName,"ThirdLepPt","Entries",100,0,200);
  
  MyhistoManager.CreateHisto(ThirdLepPt_mumumu_tight_cut4,  "ThirdLepPt_mumumu_tight_cut4" ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_mumue_tight_cut4,   "ThirdLepPt_mumue_tight_cut4"  ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_eemu_tight_cut4,    "ThirdLepPt_eemu_tight_cut4"   ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_eee_tight_cut4,     "ThirdLepPt_eee_tight_cut4"    ,datasetName,"ThirdLepPt","Entries",100,0,200);
  
  MyhistoManager.CreateHisto(ThirdLepPt_mumumu_tight_cut5,  "ThirdLepPt_mumumu_tight_cut5" ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_mumue_tight_cut5,   "ThirdLepPt_mumue_tight_cut5"  ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_eemu_tight_cut5,    "ThirdLepPt_eemu_tight_cut5"   ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_eee_tight_cut5,     "ThirdLepPt_eee_tight_cut5"    ,datasetName,"ThirdLepPt","Entries",100,0,200);
  
  MyhistoManager.CreateHisto(ThirdLepPt_mumumu_tight_cut6,  "ThirdLepPt_mumumu_tight_cut6" ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_mumue_tight_cut6,   "ThirdLepPt_mumue_tight_cut6"  ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_eemu_tight_cut6,    "ThirdLepPt_eemu_tight_cut6"   ,datasetName,"ThirdLepPt","Entries",100,0,200);
  MyhistoManager.CreateHisto(ThirdLepPt_eee_tight_cut6,     "ThirdLepPt_eee_tight_cut6"    ,datasetName,"ThirdLepPt","Entries",100,0,200);
  
  
  if (IReweight ) {
    
    if(datasetName == "DYToLL_M10-50" || datasetName == "FCNCkut" ){
      string datafile = "/opt/sbg/data/data1/cms/jandrea/TopIPHC_2012_01_25/CMSSW_4_2_8_patch7/src/MiniTreeAnalysis/NTupleAnalysis/macros/data/PUdata.root";
      string mcfile   = "/opt/sbg/data/data1/cms/jandrea/TopIPHC_2012_01_25/CMSSW_4_2_8_patch7/src/MiniTreeAnalysis/NTupleAnalysis/macros/data/PU3DMC_Fall11.root";
      
      LumiWeights    = new reweight::LumiReWeighting(mcfile, datafile, "histoMCPU", "pileup" );
    }
    else{
      string datafile = "";
      if(!IReweight_puUp && !IReweight_puDown){
	datafile = "/opt/sbg/data/data1/cms/jandrea/TopIPHC_2012_01_25/CMSSW_4_2_8_patch7/src/MiniTreeAnalysis/NTupleAnalysis/macros/data/default73.5mb.root";
      }
      
      if( IReweight_puUp)   datafile = "/opt/sbg/data/data1/cms/jandrea/TopIPHC_2012_01_25/CMSSW_4_2_8_patch7/src/MiniTreeAnalysis/NTupleAnalysis/macros/data/default73.5mbUp.root";
      if( IReweight_puDown) datafile = "/opt/sbg/data/data1/cms/jandrea/TopIPHC_2012_01_25/CMSSW_4_2_8_patch7/src/MiniTreeAnalysis/NTupleAnalysis/macros/data/default73.5mbDown.root";
      
      string mcfile   = "/opt/sbg/data/data1/cms/jandrea/TopIPHC_2012_01_25/CMSSW_4_2_8_patch7/src/MiniTreeAnalysis/NTupleAnalysis/macros/data/PU3DMC.root";
      
      LumiWeights    = new reweight::LumiReWeighting(mcfile, datafile, "histoMCPU", "pileup" );
      if(!IReweight_puUp && !IReweight_puDown)  LumiWeights->weight3D_init( 1. );
      if( IReweight_puDown                   )  LumiWeights->weight3D_init( 1. );
      if( IReweight_puUp                     )  LumiWeights->weight3D_init( 1. );
      
    }
  }
  
  JEC_L2L3Residuals.LoadCorrections();

  //************************************
  
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
  fOutput->Add(fFile);
  cout << "end SlaveBegin " << endl;
}




//_____________________________________________________________________________
void ProofSelectorMatrixMethod::FillCutFlowHistos(int IChannel, string cand3leptonChannel, bool IsSignal,
 double weight, double EventYieldWeightError, int iCut, bool useLooseIso)
{

  if(!useLooseIso) {
    if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(CutFlow_mumumu_tight, "CutFlow_mumumu_tight", iCut, datasetName, IsSignal, weight);
    if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(CutFlow_mumue_tight,  "CutFlow_mumue_tight" , iCut, datasetName, IsSignal, weight);
    if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(CutFlow_eemu_tight,   "CutFlow_eemu_tight"  , iCut, datasetName, IsSignal, weight);
    if( IChannel == 3 && cand3leptonChannel == "eee"	) MyhistoManager.FillHisto(CutFlow_eee_tight,   "CutFlow_eee_tight"   , iCut, datasetName, IsSignal, weight);
    if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(ErrCutFlow_mumumu_tight,	"ErrCutFlow_mumumu_tight"  , iCut, datasetName, IsSignal, EventYieldWeightError);
    if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(ErrCutFlow_mumue_tight,	"ErrCutFlow_mumue_tight"   , iCut, datasetName, IsSignal, EventYieldWeightError);
    if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(ErrCutFlow_eemu_tight,	"ErrCutFlow_eemu_tight"    , iCut, datasetName, IsSignal, EventYieldWeightError);
    if( IChannel == 3 && cand3leptonChannel == "eee"	) MyhistoManager.FillHisto(ErrCutFlow_eee_tight,	"ErrCutFlow_eee_tight"     , iCut, datasetName, IsSignal, EventYieldWeightError);
  }
  else {
    if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(CutFlow_mumumu_loose, "CutFlow_mumumu_loose", iCut, datasetName, IsSignal, weight);
    if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(CutFlow_mumue_loose,  "CutFlow_mumue_loose" , iCut, datasetName, IsSignal, weight);
    if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(CutFlow_eemu_loose,   "CutFlow_eemu_loose"  , iCut, datasetName, IsSignal, weight);
    if( IChannel == 3 && cand3leptonChannel == "eee"	) MyhistoManager.FillHisto(CutFlow_eee_loose,    "CutFlow_eee_loose"  , iCut, datasetName, IsSignal, weight);
    if( IChannel == 0 && cand3leptonChannel == "mumumu") MyhistoManager.FillHisto(ErrCutFlow_mumumu_loose,	"ErrCutFlow_mumumu_loose"  , iCut, datasetName, IsSignal, EventYieldWeightError);
    if( IChannel == 1 && cand3leptonChannel == "mumue" ) MyhistoManager.FillHisto(ErrCutFlow_mumue_loose,	"ErrCutFlow_mumue_loose"   , iCut, datasetName, IsSignal, EventYieldWeightError);
    if( IChannel == 2 && cand3leptonChannel == "eemu"  ) MyhistoManager.FillHisto(ErrCutFlow_eemu_loose,	"ErrCutFlow_eemu_loose"    , iCut, datasetName, IsSignal, EventYieldWeightError);
    if( IChannel == 3 && cand3leptonChannel == "eee"	) MyhistoManager.FillHisto(ErrCutFlow_eee_loose,	"ErrCutFlow_eee_loose"     , iCut, datasetName, IsSignal, EventYieldWeightError);
  }
  
}

//_____________________________________________________________________________
void ProofSelectorMatrixMethod::ApplySelectionAndFillHistos(IPHCTree::NTEvent *event, bool useLooseIso)
{
  if(!useLooseIso) useNonIsoWcand = false;
  else useNonIsoWcand = true;
  
  //---------------------------------------------------//
  //         Doing the analysis event by event
  //---------------------------------------------------//
  int debugcc=1000;
  int maxdebugcc=10;
  //cout<<"Entry "<<entry<<endl;
  sel.LoadEvent(event);
  
  //Collection of selected objects
  vector<NTVertex>   selVertices  = sel.GetSelectedVertex();
  vector<IPHCTree::NTElectron> selElectrons = sel.GetSelectedElectrons();
  vector<NTMuon>     selMuons     = sel.GetSelectedMuons();
  
  vector<NTElectron> selElectronsNonIso = sel.GetSelectedElectronsNoIso();
  vector<NTMuon>     selMuonsNonIso     = sel.GetSelectedMuonsNoIso();
  
  vector<NTElectron> selElectronsLoose;
  vector<NTMuon>     selMuonsLoose;
  
  double ElectronLooseIso = 0.8; //0.8 ; 999
  double MuonLooseIso = 0.8; //0.8 ; 999
  
  vector<NTElectron> ZeeCand; 
  vector<NTMuon>     ZmumuCand; 
  
  vector<NTElectron> WeCand; 
  vector<NTMuon>     WmuCand; 
 
  NTMET met			   = sel.GetMET(); 

  
 
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
   
    bool isData = false;
    
    if(datasetName=="DataEG" || datasetName=="DataMuEG" || datasetName=="DataMu") isData = true;
    
    
    //*****************************************************************
    // calcul the MC weights
    //*****************************************************************    
    if ( datasetName!="DataEG"   && datasetName!="DataMu" && 
	 datasetName!="DataMuEG" && datasetName!="DataEGMu" 
	 && datasetName!="MET1"  && datasetName!="MET2") {
      
      
      if(IReweight ){
	
	
	if(datasetName != "DYToLL_M10-50" && datasetName != "FCNCkut" ){
	  weightITypeMC = weightITypeMC_save*LumiWeights->weight3D(event->pileup.before_npu, event->pileup.intime_npu, event->pileup.after_npu);
	}               
	if( datasetName == "DYToLL_M10-50" || datasetName == "FCNCkut")  {
	  double ave_npu = (event->pileup.before_npu+event->pileup.intime_npu+event->pileup.after_npu)/3.;
	  weightITypeMC = weightITypeMC_save*LumiWeights->ITweight3BX(ave_npu);
	}
	
	
      }
      else weightITypeMC = weightITypeMC_save;
    }
    else weightITypeMC = 1;
    
    
    //*****************************************************************
    // determine top decay channel
    //*****************************************************************    
    
    bool IsTTbarDilept = false;
    bool IsSignal = false;
    bool IsData = false;
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
	      ) { 
      ITypeMC = 4; IsSignal = false;  Dweight[ITypeMC]= weightITypeMC; 
      EventYieldWeightError = Dweight[ITypeMC]*Dweight[ITypeMC];
      
      
 
      
      //TabFlow1[IChannel][ITypeMC][0]+=Dweight[ITypeMC];
      //TabFlow2[IChannel][ITypeMC][0]+=Dweight[ITypeMC]*Dweight[ITypeMC];
    }
    else if ( datasetName=="WZ" || datasetName=="WW" || datasetName=="ZZ"  || datasetName=="VV") { 
      ITypeMC = 5; IsSignal = false;  Dweight[ITypeMC]= weightITypeMC; 
      EventYieldWeightError = Dweight[ITypeMC]*Dweight[ITypeMC];

      //TabFlow1[IChannel][ITypeMC][0]+=Dweight[ITypeMC];
      //TabFlow2[IChannel][ITypeMC][0]+=Dweight[ITypeMC]*Dweight[ITypeMC];
     
    }  else if ( datasetName=="FCNCkut") { 
      ITypeMC = 6; IsSignal = false;  Dweight[ITypeMC]= weightITypeMC; 
      EventYieldWeightError = Dweight[ITypeMC]*Dweight[ITypeMC];

      //TabFlow1[IChannel][ITypeMC][0]+=Dweight[ITypeMC];
      //TabFlow2[IChannel][ITypeMC][0]+=Dweight[ITypeMC]*Dweight[ITypeMC];
     
    }
    
    //cout << "834 " << endl;
    if ( datasetName=="DataEG" || datasetName=="DataMu" || 
	 datasetName=="DataMuEG" || datasetName=="DataEGMu" ||
	 datasetName=="MET1" || datasetName=="MET2") { 
      ITypeMC = 100; IsData = true;  Dweight[ITypeMC]= weightITypeMC; 
      EventYieldWeightError = Dweight[ITypeMC]*Dweight[ITypeMC];   
      
      

      
      
    }
    
    //cout << "line 648 " << endl;
    //*****************************************************************
    // determine lepton scale factors
    //*****************************************************************    
    
    
    if(ApplyLeptonSF){
      //to be updated for 3 leptons selection
    }
   
	string cand3leptonChannel = "";
    //cout << "line 854  " << endl;
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
      if(IChannel==0) FillCutFlowHistos( IChannel, "mumumu", IsSignal,  Dweight[ITypeMC], EventYieldWeightError, 0, useNonIsoWcand);
      if(IChannel==1) FillCutFlowHistos( IChannel, "mumue", IsSignal,  Dweight[ITypeMC], EventYieldWeightError, 0, useNonIsoWcand);
      if(IChannel==2) FillCutFlowHistos( IChannel, "eemu", IsSignal,  Dweight[ITypeMC], EventYieldWeightError, 0, useNonIsoWcand);
      if(IChannel==3) FillCutFlowHistos( IChannel, "eee", IsSignal,  Dweight[ITypeMC], EventYieldWeightError, 0, useNonIsoWcand);


      //*****************************************************************
      // pass 3 lepton requirements
      //*****************************************************************    
      
      
      //*****************************************************************
      // select Z->ee candidate
      //*****************************************************************    
      
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
	      theel1 = iel1;
	      theel2 = iel2;
	      mInv = zeecand.M();
	    }
	  }
	}
	
	if(theel1>=0 && theel2>=0){ //JLA
	  ZeeCand.push_back(selElectrons[theel1]);
	  ZeeCand.push_back(selElectrons[theel2]);
	}
      }
         
      //*****************************************************************
      // select W->enu candidate
      //*****************************************************************    
      
      if(!useNonIsoWcand){
        for(unsigned int iel1 = 0; iel1 < selElectrons.size(); iel1++){
	  bool matchElec=false;
          for(unsigned int iel2 = 0; iel2 < ZeeCand.size(); iel2++){
             
	     if(fabs(selElectrons[iel1].p4.Pt() - ZeeCand[iel2].p4.Pt()) <  0.0001)  matchElec=true;
	   
          }
	  if(!matchElec && selElectrons[iel1].p4.Pt()>30) WeCand.push_back(selElectrons[iel1]);
        }
      }else{
      
        selElectronsLoose.clear();
	for(unsigned int iel1 = 0; iel1 < selElectronsNonIso.size(); iel1++)
          if(selElectronsNonIso[iel1].RelIso03PF()<ElectronLooseIso) selElectronsLoose.push_back(selElectronsNonIso[iel1]);

        for(unsigned int iel1 = 0; iel1 < selElectronsLoose.size(); iel1++){
	  bool matchElec=false;
          for(unsigned int iel2 = 0; iel2 < ZeeCand.size(); iel2++){
             
	     if(fabs(selElectronsLoose[iel1].p4.Pt() - ZeeCand[iel2].p4.Pt()) <  0.0001)  matchElec=true;
	   
          }
	  if(!matchElec && selElectronsLoose[iel1].p4.Pt()>30) WeCand.push_back(selElectronsLoose[iel1]);
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
      
      if(!useNonIsoWcand){
        //cout << "in sel W " << endl;
        for(unsigned int imu1 = 0; imu1 < selMuons.size(); imu1++){
	  bool matchMuon = false;
          for(unsigned int imu2 = 0; imu2 < ZmumuCand.size(); imu2++){
             
	     if(fabs(selMuons[imu1].p4.Pt() - ZmumuCand[imu2].p4.Pt()) <  0.0001) matchMuon = true;
	     
          } 
	  if(!matchMuon && selMuons[imu1].p4.Pt()>30) WmuCand.push_back(selMuons[imu1]);
        }
      }else{

        selMuonsLoose.clear();
        for(unsigned int imu1 = 0; imu1 < selMuonsNonIso.size(); imu1++)
          if(selMuonsNonIso[imu1].RelIso03PF()<MuonLooseIso) selMuonsLoose.push_back(selMuonsNonIso[imu1]);

        for(unsigned int imu1 = 0; imu1 < selMuonsLoose.size(); imu1++){
	  bool matchMuon = false;
          for(unsigned int imu2 = 0; imu2 < ZmumuCand.size(); imu2++){
             
	     if(fabs(selMuonsLoose[imu1].p4.Pt() - ZmumuCand[imu2].p4.Pt())  < 0.0001) matchMuon = true;
	 
          }
	  if(!matchMuon && selMuonsLoose[imu1].p4.Pt()>30) WmuCand.push_back(selMuonsLoose[imu1]);
        }
      }
      
      
      //*****************************************************************
      // apply lepton selection
      //*****************************************************************    
      
      if( (WmuCand.size()+ZmumuCand.size()+WeCand.size()+ZeeCand.size()) == 3 && met.p2.Pt() > 25) {
      
        //cout<<"Wmu "<<WmuCand.size()<<" Zmu "<<ZmumuCand.size()<<"  We "<<WeCand.size()<<" Ze "<<ZeeCand.size()<<endl;
       
	string cand3leptonChannel = "";
	if( ZmumuCand.size() == 2 ) {
	  if(WmuCand.size() == 1 ) cand3leptonChannel = "mumumu";
	  if(WeCand.size()  == 1 ) cand3leptonChannel = "mumue";
	}
	
	if( ZeeCand.size() == 2 ) {
	  if(WmuCand.size() == 1 ) cand3leptonChannel = "eemu";
	  if(WeCand.size() == 1 ) cand3leptonChannel = "eee";
	}
		
        FillCutFlowHistos( IChannel, cand3leptonChannel, IsSignal,  Dweight[ITypeMC], EventYieldWeightError, 1, useNonIsoWcand);
 
 			      
	//*****************************************************************
        // select Z candidate
        //*****************************************************************    
	
	 
	
	TLorentzVector dilept;
	TLorentzVector lept1, lept2, lept3;
	
	string theleptpair = "";
	
	if(cand3leptonChannel == "mumue") {
	  lept1 = ZmumuCand[0].p4;
	  lept2 = ZmumuCand[1].p4;
	  lept3 = WeCand[0].p4;
	  dilept = lept1+lept2;
	}
	
	
	if(cand3leptonChannel == "eemu") {
	  lept1 = ZeeCand[0].p4;
	  lept2 = ZeeCand[1].p4;
	  lept3 = WmuCand[0].p4;
	  dilept = lept1+lept2;
	}
	
	if(cand3leptonChannel == "mumumu") {
	  
	  lept1 = ZmumuCand[0].p4;
	  lept2 = ZmumuCand[1].p4;
	  lept3 = WmuCand[0].p4;
	  dilept = lept1+lept2;
	  
	}
	
	if(cand3leptonChannel == "eee") {
	  
	  lept1 = ZeeCand[0].p4;
	  lept2 = ZeeCand[1].p4;
	  lept3 = WeCand[0].p4;
	  dilept = lept1+lept2;
	  
	}
	
	
       if( IChannel == 0 && cand3leptonChannel == "mumumu") 
        MyhistoManager.FillHisto(ThirdLepPt_mumumu_tight_cut1, "ThirdLepPt_mumumu_tight_cut1", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
       if( IChannel == 1 && cand3leptonChannel == "mumue" ) 
        MyhistoManager.FillHisto(ThirdLepPt_mumue_tight_cut1, "ThirdLepPt_mumue_tight_cut1", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
       if( IChannel == 2 && cand3leptonChannel == "eemu"  ) 
        MyhistoManager.FillHisto(ThirdLepPt_eemu_tight_cut1, "ThirdLepPt_eemu_tight_cut1", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
       if( IChannel == 3 && cand3leptonChannel == "eee"   ) 
        MyhistoManager.FillHisto(ThirdLepPt_eee_tight_cut1, "ThirdLepPt_eee_tight_cut1", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	
	
	
	
        //cout << "line 1025" << endl;
	vector<NTJet>      selJets = sel.GetSelectedJets(selMuons, selElectrons);
	double dileptonIvM = dilept.M();
	
	double sumPtLeptJet = lept1.Pt() + lept2.Pt() + lept3.Pt();
	for(unsigned int i=0 ; i< selJets.size(); i++) sumPtLeptJet += selJets[i].p4.Pt();
	double theMET = met.p2.Mod();
		
	int NBtaggedJets = 0;
	int AlgoBtag = sel.GetbtagAlgo();
	float btagDiscriCut = sel.GetbtagDiscriCut ();
	//cout << "***********************" << endl;
	for(unsigned int ijet = 0; ijet < selJets.size(); ijet++){
	  
	  if ( AlgoBtag==0 &&  selJets[ijet].bTag["trackCountingHighEffBJetTags"]	  >= btagDiscriCut) NBtaggedJets++;
	  if ( AlgoBtag==1 &&  selJets[ijet].bTag["simpleSecondaryVertexHighEffBJetTags"]    >= btagDiscriCut) NBtaggedJets++;
	  if ( AlgoBtag==2 &&  selJets[ijet].bTag["trackCountingHighPurBJetTags"]	  >= btagDiscriCut) NBtaggedJets++;
	  if ( AlgoBtag==3 &&  selJets[ijet].bTag["simpleSecondaryVertexHighPurBJetTags"]    >= btagDiscriCut) NBtaggedJets++;
	  if ( AlgoBtag==4 &&  selJets[ijet].bTag["jetProbabilityBJetTags"]		  >= btagDiscriCut) NBtaggedJets++;
	  if ( AlgoBtag==5 &&  selJets[ijet].bTag["jetBProbabilityBJetTags"]		  >= btagDiscriCut) NBtaggedJets++;
	  if ( AlgoBtag==6 &&  selJets[ijet].bTag["combinedSecondaryVertexBJetTags"]	  >= btagDiscriCut) NBtaggedJets++;
	
	  //cout << "jet " << ijet << " has a flavor " << selJets[ijet].partonFlavour << endl;
	  //cout << "NBtaggedJets " << NBtaggedJets << endl;
	}  
	
		
	
	//*****************************************************************
        // pass Z mass
        //*****************************************************************  
			
	if( fabs(dilept.M()-91) < 15){
	
          FillCutFlowHistos( IChannel, cand3leptonChannel, IsSignal,  Dweight[ITypeMC], EventYieldWeightError, 2, useNonIsoWcand);
	  if( IChannel == 0 && cand3leptonChannel == "mumumu") 
           MyhistoManager.FillHisto(ThirdLepPt_mumumu_tight_cut2, "ThirdLepPt_mumumu_tight_cut2", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  if( IChannel == 1 && cand3leptonChannel == "mumue" ) 
           MyhistoManager.FillHisto(ThirdLepPt_mumue_tight_cut2, "ThirdLepPt_mumue_tight_cut2", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  if( IChannel == 2 && cand3leptonChannel == "eemu"  ) 
           MyhistoManager.FillHisto(ThirdLepPt_eemu_tight_cut2, "ThirdLepPt_eemu_tight_cut2", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	  if( IChannel == 3 && cand3leptonChannel == "eee"   ) 
           MyhistoManager.FillHisto(ThirdLepPt_eee_tight_cut2, "ThirdLepPt_eee_tight_cut2", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);

	  
	  //*****************************************************************
	  // pass jet selection
	  //*****************************************************************  
	  int NSeljets = selJets.size() ;
	  if(NSeljets>4) NSeljets = 4;

	  if(selJets.size() >= 1 ){
	  	
            FillCutFlowHistos( IChannel, cand3leptonChannel, IsSignal,  Dweight[ITypeMC], EventYieldWeightError, 3, useNonIsoWcand);
	    if( IChannel == 0 && cand3leptonChannel == "mumumu") 
             MyhistoManager.FillHisto(ThirdLepPt_mumumu_tight_cut3, "ThirdLepPt_mumumu_tight_cut3", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	    if( IChannel == 1 && cand3leptonChannel == "mumue" ) 
             MyhistoManager.FillHisto(ThirdLepPt_mumue_tight_cut3, "ThirdLepPt_mumue_tight_cut3", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	    if( IChannel == 2 && cand3leptonChannel == "eemu"  ) 
             MyhistoManager.FillHisto(ThirdLepPt_eemu_tight_cut3, "ThirdLepPt_eemu_tight_cut3", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	    if( IChannel == 3 && cand3leptonChannel == "eee"   ) 
             MyhistoManager.FillHisto(ThirdLepPt_eee_tight_cut3, "ThirdLepPt_eee_tight_cut3", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);

	    //*****************************************************************
	    // pass btag selection
	    //*****************************************************************  

	    if(NBtaggedJets == 1 ){
	      	  	   	      
              FillCutFlowHistos( IChannel, cand3leptonChannel, IsSignal,  Dweight[ITypeMC], EventYieldWeightError, 4, useNonIsoWcand);
	      if( IChannel == 0 && cand3leptonChannel == "mumumu") 
               MyhistoManager.FillHisto(ThirdLepPt_mumumu_tight_cut4, "ThirdLepPt_mumumu_tight_cut4", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      if( IChannel == 1 && cand3leptonChannel == "mumue" ) 
               MyhistoManager.FillHisto(ThirdLepPt_mumue_tight_cut4, "ThirdLepPt_mumue_tight_cut4", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      if( IChannel == 2 && cand3leptonChannel == "eemu"  ) 
               MyhistoManager.FillHisto(ThirdLepPt_eemu_tight_cut4, "ThirdLepPt_eemu_tight_cut4", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
	      if( IChannel == 3 && cand3leptonChannel == "eee"   ) 
               MyhistoManager.FillHisto(ThirdLepPt_eee_tight_cut4, "ThirdLepPt_eee_tight_cut4", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);

    
   	      //*****************************************************************
	      // pass Mtop selection
	      //*****************************************************************  

	      TLorentzVector metP4(met.p2.Px(), met.p2.Py(), 0, sqrt(met.p2.Px()*met.p2.Px() + met.p2.Py()*met.p2.Py()));
	      TLorentzVector transTop = lept3 + selJets[0].p4 + metP4;

	      if(transTop.Mt() > 150){

                FillCutFlowHistos( IChannel, cand3leptonChannel, IsSignal,  Dweight[ITypeMC], EventYieldWeightError, 5, useNonIsoWcand);
		if( IChannel == 0 && cand3leptonChannel == "mumumu") 
        	 MyhistoManager.FillHisto(ThirdLepPt_mumumu_tight_cut5, "ThirdLepPt_mumumu_tight_cut5", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		if( IChannel == 1 && cand3leptonChannel == "mumue" ) 
        	 MyhistoManager.FillHisto(ThirdLepPt_mumue_tight_cut5, "ThirdLepPt_mumue_tight_cut5", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		if( IChannel == 2 && cand3leptonChannel == "eemu"  ) 
        	 MyhistoManager.FillHisto(ThirdLepPt_eemu_tight_cut5, "ThirdLepPt_eemu_tight_cut5", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
		if( IChannel == 3 && cand3leptonChannel == "eee"   ) 
        	 MyhistoManager.FillHisto(ThirdLepPt_eee_tight_cut5, "ThirdLepPt_eee_tight_cut5", lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);

	      }
	    } // end selection btag 
	  } // end selection njet
	} //end selection Z cand    
      } // end selection 3 leptons
    } // end selection trigger
  } //end loops on datasets
  
  
}

//_____________________________________________________________________________
Bool_t ProofSelectorMatrixMethod::Process(Long64_t entry)
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

  ApplySelectionAndFillHistos(event, false); //tight selection
  ApplySelectionAndFillHistos(event, true); //loose selection
  
}

//_____________________________________________________________________________
void ProofSelectorMatrixMethod::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
  if(fProofFile) fProofFile->Print();
  if (fFile) {
    Bool_t cleanup = kFALSE;
    TDirectory *savedir = gDirectory;
    fFile->cd();
    
   
  MyhistoManager.WriteMyHisto(CutFlow_mumumu_tight, "all" );
  MyhistoManager.WriteMyHisto(CutFlow_mumue_tight,  "all" );
  MyhistoManager.WriteMyHisto(CutFlow_eemu_tight,   "all" );
  MyhistoManager.WriteMyHisto(CutFlow_eee_tight,    "all" );
  
  MyhistoManager.WriteMyHisto(ErrCutFlow_mumumu_tight,  "all");
  MyhistoManager.WriteMyHisto(ErrCutFlow_mumue_tight,   "all");
  MyhistoManager.WriteMyHisto(ErrCutFlow_eemu_tight,    "all");
  MyhistoManager.WriteMyHisto(ErrCutFlow_eee_tight,     "all");
  
  MyhistoManager.WriteMyHisto(CutFlow_mumumu_loose, "all" );
  MyhistoManager.WriteMyHisto(CutFlow_mumue_loose,  "all" );
  MyhistoManager.WriteMyHisto(CutFlow_eemu_loose,   "all" );
  MyhistoManager.WriteMyHisto(CutFlow_eee_loose,    "all" );
  
  MyhistoManager.WriteMyHisto(ErrCutFlow_mumumu_loose,  "all");
  MyhistoManager.WriteMyHisto(ErrCutFlow_mumue_loose,   "all");
  MyhistoManager.WriteMyHisto(ErrCutFlow_eemu_loose,    "all");
  MyhistoManager.WriteMyHisto(ErrCutFlow_eee_loose,     "all");
  
  MyhistoManager.WriteMyHisto(ThirdLepPt_mumumu_tight_cut1, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_mumue_tight_cut1, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_eemu_tight_cut1, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_eee_tight_cut1, "all" );
   
  MyhistoManager.WriteMyHisto(ThirdLepPt_mumumu_tight_cut2, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_mumue_tight_cut2, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_eemu_tight_cut2, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_eee_tight_cut2, "all" );

  MyhistoManager.WriteMyHisto(ThirdLepPt_mumumu_tight_cut3, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_mumue_tight_cut3, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_eemu_tight_cut3, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_eee_tight_cut3, "all" );

  MyhistoManager.WriteMyHisto(ThirdLepPt_mumumu_tight_cut4, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_mumue_tight_cut4, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_eemu_tight_cut4, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_eee_tight_cut4, "all" );

  MyhistoManager.WriteMyHisto(ThirdLepPt_mumumu_tight_cut5, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_mumue_tight_cut5, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_eemu_tight_cut5, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_eee_tight_cut5, "all" );

  MyhistoManager.WriteMyHisto(ThirdLepPt_mumumu_tight_cut6, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_mumue_tight_cut6, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_eemu_tight_cut6, "all" );
  MyhistoManager.WriteMyHisto(ThirdLepPt_eee_tight_cut6, "all" );

   //The following line is mandatory to copy everythin in a common RootFile
    fOutput->Add(fProofFile);
    
  CutFlow_mumumu_tight.clear();
  CutFlow_mumue_tight.clear();
  CutFlow_eemu_tight.clear();
  CutFlow_eee_tight.clear();
  
  ErrCutFlow_mumumu_tight.clear();
  ErrCutFlow_mumue_tight.clear();
  ErrCutFlow_eemu_tight.clear();
  ErrCutFlow_eee_tight.clear();
  
  CutFlow_mumumu_loose.clear();
  CutFlow_mumue_loose.clear();
  CutFlow_eemu_loose.clear();
  CutFlow_eee_loose.clear();
  
  ErrCutFlow_mumumu_loose.clear();
  ErrCutFlow_mumue_loose.clear();
  ErrCutFlow_eemu_loose.clear();
  ErrCutFlow_eee_loose.clear();
  
  ThirdLepPt_mumumu_tight_cut1.clear();
  ThirdLepPt_mumue_tight_cut1.clear();
  ThirdLepPt_eemu_tight_cut1.clear();
  ThirdLepPt_eee_tight_cut1.clear();
  
  ThirdLepPt_mumumu_tight_cut2.clear();
  ThirdLepPt_mumue_tight_cut2.clear();
  ThirdLepPt_eemu_tight_cut2.clear();
  ThirdLepPt_eee_tight_cut2.clear();
  
  ThirdLepPt_mumumu_tight_cut3.clear();
  ThirdLepPt_mumue_tight_cut3.clear();
  ThirdLepPt_eemu_tight_cut3.clear();
  ThirdLepPt_eee_tight_cut3.clear();
  
  ThirdLepPt_mumumu_tight_cut4.clear();
  ThirdLepPt_mumue_tight_cut4.clear();
  ThirdLepPt_eemu_tight_cut4.clear();
  ThirdLepPt_eee_tight_cut4.clear();
  
  ThirdLepPt_mumumu_tight_cut5.clear();
  ThirdLepPt_mumue_tight_cut5.clear();
  ThirdLepPt_eemu_tight_cut5.clear();
  ThirdLepPt_eee_tight_cut5.clear();
  
  ThirdLepPt_mumumu_tight_cut6.clear();
  ThirdLepPt_mumue_tight_cut6.clear();
  ThirdLepPt_eemu_tight_cut6.clear();
  ThirdLepPt_eee_tight_cut6.clear();
  
    delete anaEL;
    delete LumiWeights;
  }
}

//_____________________________________________________________________________
void ProofSelectorMatrixMethod::Terminate()
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
