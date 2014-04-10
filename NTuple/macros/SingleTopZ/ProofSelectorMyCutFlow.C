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
  
  
  //------------------------// 
  //initialize the variables
  //------------------------// 
  //------------------------// 
  //for PU
  IReweight             = true;
  IReweight_puUp	= false;
  IReweight_puDown	= false;
  IReweight_pu	        = false;
  //------------------------// 
  //for JES uncertainties
  doJESuncert = false;
  upOrDown    = false;
  applyJES    = false; 
  applyJER    = false;
  ResFactor   = 0.1;
  //------------------------// 
  //for PDF unceratinties
  doPDF = false;
  //pdftype =0 ;
  pdftype =1  ;
  //------------------------// 
  //loose iso on W
  //for backgr studies
  useNonIsoWcand = false;
  looseIso = 0.4; //0.4
  themetcut = 35;
  //------------------------// 
  applyTrigger     = false;
  applyTriggerUp   = false;
  applyTriggerDown = false;
  //------------------------// 
  applyLeptonSF      = false;
  applyLeptonSFUp    = false;
  applyLeptonSFDown  = false;
  //------------------------// 
  
  
  
  //------------------------// 
  //initialize scale factors
  //------------------------// 
  //trigger SF
  SF_trig_mumumu = 1.0; 
  SF_trig_mumue  = 1.0; 
  SF_trig_eemu   = 1.0; 
  SF_trig_eee    = 1.0;
  SF_trig_mumumu_error  = 0.0; 
  SF_trig_mumue_error   = 0.0;
  SF_trig_eemu_error    = 0.0; 
  SF_trig_eee_error     = 0.0; 
  
  //------------------------// 
  //for WZ+jets
  applyWZ          = false;
  SF_WZ.push_back(1.0); //mumumu
  SF_WZ.push_back(1.0); //mumue
  SF_WZ.push_back(1.0); //eemu
  SF_WZ.push_back(1.0); //eee
  //------------------------// 
  //for fake leptons
  applyFakescale   = false;
  SF_Fake.push_back(1.0); //mumumu
  SF_Fake.push_back(1.0); //mumue
  SF_Fake.push_back(1.0); //eemu
  SF_Fake.push_back(1.0); //eee
  
  
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
  sel.LoadElScaleFactors();
  sel.LoadMuScaleFactors();
  
  
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
     datasetName=="MET1" || datasetName=="MET2") isData_sample = true;
  
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
  
  
  
  //--------------------------------------//
  //   Managing histos  	
  //--------------------------------------//
  MyhistoManager.LoadDatasets(datasets);   
  MyhistoManager.LoadSelectionSteps(CutName);
  MyhistoManager.LoadChannels(TheChannelName);
  
  ITypeMC = -1 ;
  
  
  
  createTheHisto(&MyhistoManager);
  
  
  //--------------------------------------//
  //   Output TTree 	
  //--------------------------------------//
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
  
  
  TheTree->Branch("tree_met", &tree_met , "tree_met/F" );
  TheTree->Branch("tree_mTW", &tree_mTW , "tree_mTW/F" );
  
  
  TheTree->Branch("tree_EvtWeight",   &tree_EvtWeight,   "tree_EvtWeight/F" );
  TheTree->Branch("tree_SampleType",  &tree_SampleType,  "tree_SampleType/I");
  TheTree->Branch("tree_Channel",     &tree_Channel,     "tree_Channel/I"   );


  
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
  
  
  
  //---------------------------------
  //load events 
  //---------------------------------
  sel.LoadEvent(event);
  
  
  
  
  //cout<<"RUN "<<event->general.runNb<<" EVT "<<event->general.eventNb<<endl;
  
  //std::pair<int, int> thepair(event->general.runNb, event->general.eventNb) ;
  
  //if(the_run_evt_map;
  
  //---------------------------------
  //get collections from events 
  //---------------------------------
  double rho = event->pileup.rho_PUUE_dens;
  
  //Collection of selected objects
  vector<NTVertex>   selVertices        = sel.GetSelectedVertex();
  vector<NTElectron> selElectrons       = sel.GetSelectedElectronsRhoIso(20, 2.5, 0.15, 0, 0, 0, 0, rho); 
  vector<NTMuon>     selMuons           = sel.GetSelectedMuonsDeltaBetaIso();
  vector<NTElectron> selElectronsNonIso = sel.GetSelectedElectronsNoIso();
  vector<NTMuon>     selMuonsNonIso     = sel.GetSelectedMuonsNoIso();
  NTMET met   = sel.GetSelectedMET(applyJES, &*theJESuncertainty, upOrDown , applyJER, ResFactor);  
  
  
  //---------------------------------
  //initiate lepton candidate collections 
  //---------------------------------
  vector<NTElectron> ZeeCand; 
  vector<NTMuon>     ZmumuCand; 
  
  vector<NTElectron> WeCand; 
  vector<NTMuon>     WmuCand; 
 

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
       datasetName!="DataMuEG" ) { 
    if(IReweight ){
      weightITypeMC = weightITypeMC_save*LumiWeights->weight_Summer12_S10(event->pileup.Tnpv);
    }
    else weightITypeMC = weightITypeMC_save;
  }
  else weightITypeMC = 1;
  
  
  
  
  //*****************************************************************
  // determine top decay channel
  //*****************************************************************    
  bool IsTTbarDilept = false;
  bool IsSignal      = false;
  double WeightForBranchingRatio = 1.;
  
  
  
  //*****************************************************************
  // determine MC evetn weight
  //*****************************************************************
  
  std::vector< double > thereweight = determineWeights(datasetName, weightITypeMC, WeightForBranchingRatio);
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
    if(passtrigger_mumu  && !passtrigger_emu && !passtrigger_ee && datasetName=="DataMu"  )   IChannel = 0;
    if(passtrigger_emu   &&                                        datasetName=="DataMuEG")   IChannel = 12;
    if(passtrigger_ee    && !passtrigger_emu &&  passtrigger_ee && datasetName=="DataEG"  )   IChannel = 3;
  }else{
    if( passtrigger_mumu && !passtrigger_emu && !passtrigger_ee ) IChannel = 0;
    if( passtrigger_emu                                         ) IChannel = 12;
    if(!passtrigger_mumu && !passtrigger_emu &&  passtrigger_ee ) IChannel = 3;
  }
  /*cout << "------------------------------------------" << endl;  
    cout << "passtrigger_mumu " << passtrigger_mumu << endl;
  cout << "passtrigger_emu  " << passtrigger_emu  << endl;
  cout << "passtrigger_ee   " << passtrigger_ee   << endl;
  cout << "selElectrons.size() " << selElectrons.size() << endl;
  cout << "selMuons.size()     " << selMuons.size() << endl;
  cout << "IChannel " << IChannel << endl;*/
  
  
  
  if( IChannel != -1 ){
    
    //if(IChannel != 0 && passtrigger_mumu == 1){
    //  cout << "problemmmmmmm "  <<  IChannel << endl;
    //  cout << " datasetName " <<  datasetName << endl;
    //}
    
    for(unsigned int i=0; i<selMuons.size(); i++){
      
      if(selMuons[i].p4.Pt() < 20 || fabs(selMuons[i].p4.Eta()) > 2.5 ) cout << "unselected leptons with pt =" <<  selMuons[i].p4.Pt()  << " and eta = " << selMuons[i].p4.Eta() << endl;
      
    }
    
    
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


    
    
    
    
    determineLeptonCandidates(useNonIsoWcand, looseIso, rho ,
			      &selElectrons,       &selMuons, 
			      &selElectronsNonIso, &selMuonsNonIso, 
			      &ZeeCand, &ZmumuCand, 
			      &WeCand,  &WmuCand);
    
    
    if(ZeeCand.size()>0 ){
      if( IChannel == 12 ) MyhistoManager.FillHisto(CutFlow_eemu,   "CutFlow_eemu"  , 1, datasetName, IsSignal, Dweight[ITypeMC]);
      if( IChannel == 3 ) MyhistoManager.FillHisto(CutFlow_eee,    "CutFlow_eee"   , 1, datasetName, IsSignal, Dweight[ITypeMC]);
    }
    
    if(ZmumuCand.size()>0 ){
       if( IChannel == 0 ) MyhistoManager.FillHisto(CutFlow_mumumu, "CutFlow_mumumu", 1, datasetName, IsSignal, Dweight[ITypeMC]);
       if( IChannel == 12 ) MyhistoManager.FillHisto(CutFlow_mumue,  "CutFlow_mumue" , 1, datasetName, IsSignal, Dweight[ITypeMC]);
    }
    
    
    
    /*
      cout << "ZeeCand   " << ZeeCand.size()    << endl;
      cout << "ZmumuCand " << ZmumuCand.size()  << endl;
      cout << "WeCand    " << WeCand.size()     << endl;
      cout << "WmuCand   " << WmuCand.size()    << endl;
    */
    //*************************************************
    //              select 3 lepton events
    //*************************************************
    
    TString leptChannel = "";
    //if(isData){
    if( ZmumuCand.size() == 2 ) {
      if(WmuCand.size() == 1 && IChannel== 0  )  {	       leptChannel = "mumumu";}// for mumumu events
      if(WeCand.size()  == 1 && IChannel== 12 )  {IChannel= 1; leptChannel = "mumue"; }// for mumue events
    }
    
    if( ZeeCand.size() == 2 ) {
      if(WmuCand.size() == 1  && IChannel== 12  )  {IChannel= 2; leptChannel = "eemu";  }// for eemu events
      if(WeCand.size()  == 1  && IChannel== 3	)  {		leptChannel = "eee";   }// for eee events
    }
    
    
    //if( (WmuCand.size()+ZmumuCand.size()+WeCand.size()+ZeeCand.size()) == 3) {
    if( ( (selElectrons.size() + selMuons.size()) == 3 ) 
       && (
	   (IChannel== 0 && leptChannel == "mumumu") ||
	   (IChannel== 1 && leptChannel == "mumue") ||
	   (IChannel== 2 && leptChannel == "eemu") ||
	   (IChannel== 3 && leptChannel == "eee") )
	
	) {
      
      
      //cout << "leptChannel " <<leptChannel  << endl;
      
      
      if(IChannel == 12){
	
	cout << "ZmumuCand.size()  " << ZmumuCand.size() << endl;
	cout << "ZeeCand.size()    " << ZeeCand.size()   << endl;
	cout << "WmuCand.size()    " << WmuCand.size()   << endl;
	cout << "WeCand.size()     " << WeCand.size()    << endl;
	
	if(WeCand.size()> 1) {
	  for(unsigned int i=0;i< WeCand.size(); i++){cout << "WeCand["<<i<<"].Charge " << WeCand[i].charge << endl;}
	}
	if(WmuCand.size()> 1) {
	  for(unsigned int i=0;i< WmuCand.size(); i++){cout << "WmuCand["<<i<<"].Charge " << WmuCand[i].charge << endl;}
	}
      } 
      
      //create pointers to histograms
      
      
      defineHistoPointer(IChannel);
      cout << "IChannel " << IChannel << "  leptChannel " <<  leptChannel << endl;
      //if(leptChannel != "mumumu") continue;
      nselevents_mumumu++;
      //if( leptChannel == "mumumu"){
      /* if( fabs(ZmumuCand[0].p4.Pt() - ZmumuCand[1].p4.Pt() ) < 0.0001) cout << "same lept in Z cand " << endl;
	 if( fabs(ZmumuCand[0].p4.Pt() - WmuCand[0].p4.Pt()   ) < 0.0001) cout << "Z[0] cand match W cand " << endl;
	 if( fabs(ZmumuCand[1].p4.Pt() - WmuCand[0].p4.Pt()   ) < 0.0001) cout << "Z[1] cand match W cand" << endl;
	 
	 if( sel.RelIso03PFDeltaBeta(ZmumuCand[0]) > 0.15 ) cout <<        sel.RelIso03PFDeltaBeta(ZmumuCand[0] ) << endl;
	 if( sel.RelIso03PFDeltaBeta(ZmumuCand[1]) > 0.15 ) cout <<        sel.RelIso03PFDeltaBeta(ZmumuCand[1] ) << endl;
	 if( sel.RelIso03PFDeltaBeta(WmuCand[0]  ) > 0.15 ) cout <<        sel.RelIso03PFDeltaBeta(WmuCand[0]   ) << endl;
	 
	 
	 if(ZmumuCand[0].p4.Pt() < 20)  cout << "found low pt lepton from Z1 " << endl;
	 if(ZmumuCand[1].p4.Pt() < 20)  cout << "found low pt lepton fromZ2  " << endl;
	 if(WmuCand[0].p4.Pt()   < 20)  cout << "found low pt lepton fromW   " << endl;
	 
	 if(fabs(ZmumuCand[0].p4.Eta()) > 2.5)  cout << "found low pt lepton from Z1 " << endl;
	 if(fabs(ZmumuCand[1].p4.Eta()) > 2.5)  cout << "found low pt lepton fromZ2  " << endl;
	 if(fabs(WmuCand[0].p4.Eta()  ) > 2.5)  cout << "found low pt lepton fromW   " << endl;*/
        
      //} 
      
      
      
      if(leptChannel == "mumumu"){
      
      cout << "---------------------------------------------------------" << endl;
      cout << "Z[0] pt " << ZmumuCand[0].p4.Pt() << ", eta " <<  ZmumuCand[0].p4.Eta()  <<  
        ", iso   "  << sel.RelIso03PFDeltaBeta(ZmumuCand[0])  << 
        ", GM    "  << ZmumuCand[0].isGlobalMuon << 
        ", TM    "  << ZmumuCand[0].isTrackerMuon<< 
        ", C2    "  << ZmumuCand[0].Chi2	   << 
        ", NHit  "  << ZmumuCand[0].NTrValidHits  << 
        ", NVhit "  << ZmumuCand[0].NValidHits    << 
        ", DO    "  <<  fabs( ZmumuCand[0].D0Inner) << endl; 
      cout << "Z[1] pt " << ZmumuCand[1].p4.Pt() << ", eta " <<  ZmumuCand[1].p4.Eta()  <<  
        ", iso   "  << sel.RelIso03PFDeltaBeta(ZmumuCand[1])  << 
        ", GM    "  << ZmumuCand[1].isGlobalMuon << 
        ", TM    "  << ZmumuCand[1].isTrackerMuon<< 
        ", C2    "  << ZmumuCand[1].Chi2	   << 
        ", NHit  "  << ZmumuCand[1].NTrValidHits  << 
        ", NVhit "  << ZmumuCand[1].NValidHits    << 
        ", DO    "  << fabs(ZmumuCand[1].D0Inner) << endl; 
      cout << "W    pt " << WmuCand[0].p4.Pt() << ", eta " <<  WmuCand[0].p4.Eta()  <<  
        ", iso   "  << sel.RelIso03PFDeltaBeta(WmuCand[0])  << 
        ", GM    "  << WmuCand[0].isGlobalMuon << 
        ", TM    "  << WmuCand[0].isTrackerMuon<< 
        ", C2    "  << WmuCand[0].Chi2	   << 
        ", NHit  "  << WmuCand[0].NTrValidHits  << 
        ", NVhit "  << WmuCand[0].NValidHits    << 
        ", DO    "  << fabs(WmuCand[0].D0Inner) << endl; 
	cout << "Dweight[ITypeMC] " << Dweight[ITypeMC]  << endl;
      
      }
      
      
      
      if(pCutFlow == 0)     cout << "null pointer pCutFlow   " << endl;
      if(pErrCutFlow == 0)  cout << "null pointer pErrCutFlow" << endl;
      
      
      
      MyhistoManager.FillHisto(*pCutFlow,      ("CutFlow_"+leptChannel).Data(),    2, datasetName, IsSignal, Dweight[ITypeMC]);   
       
      MyhistoManager.FillHisto(*pErrCutFlow,   ("ErrCutFlow_"+leptChannel).Data(), 2, datasetName, IsSignal, EventYieldWeightError);
      
      
      
      
      
      
      //*****************************************************************
      // define  Z/W candidate TLorentzVector
      //*****************************************************************    
      
      double lept3_Charge = 0;
      TLorentzVector dilept;
      TLorentzVector lept1, lept2, lept3;
      if(leptChannel == "mumue") {
	lept1 = ZmumuCand[0].p4;
	lept2 = ZmumuCand[1].p4;
	lept3 = WeCand[0].p4;
	lept3_Charge= WeCand[0].charge;
	dilept = lept1+lept2;
      }else  if(leptChannel == "eemu") {
	lept1 = ZeeCand[0].p4;
	lept2 = ZeeCand[1].p4;
	lept3 = WmuCand[0].p4;
	lept3_Charge= WmuCand[0].charge;
	dilept = lept1+lept2;
      }else  if(leptChannel == "mumumu") {
	lept1 = ZmumuCand[0].p4;
	lept2 = ZmumuCand[1].p4;
	lept3 = WmuCand[0].p4;
	lept3_Charge= WmuCand[0].charge;
	dilept = lept1+lept2;
      }else  if(leptChannel == "eee") {
	lept1 = ZeeCand[0].p4;
	lept2 = ZeeCand[1].p4;
	lept3 = WeCand[0].p4;
	lept3_Charge= WeCand[0].charge;
	dilept = lept1+lept2;
      }
      
      
            

      int wcharge = lept3_Charge;
      
      //*****************************************************************
      // apply trigger scale factors
      //*****************************************************************  
            

      
      if(applyTrigger  &&  !isData ){ 
	if(IChannel == 0 )Dweight[ITypeMC]*=SF_trig_mumumu;
	if(IChannel == 1 )Dweight[ITypeMC]*=SF_trig_mumue;
	if(IChannel == 2 )Dweight[ITypeMC]*=SF_trig_eemu;
	if(IChannel == 3 )Dweight[ITypeMC]*=SF_trig_eee;
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
            

      double dileptonIvM = dilept.M();
      
      //***************************************************************** 
      //get selected jets 
      //*****************************************************************
      vector<NTJet>	  selJets = sel.GetSelectedJets(selMuons, selElectrons, applyJES, theJESuncertainty, upOrDown, applyJER, ResFactor);
      
            

      double sumPtLeptJet = lept1.Pt() + lept2.Pt() + lept3.Pt();
      for(unsigned int i=0 ; i< selJets.size(); i++) sumPtLeptJet += selJets[i].p4.Pt();
      double theMET = met.p2.Mod();
      
            

      //***************************************************************** 
      //get b-tag info 
      //*****************************************************************
      
            

      int NBtaggedJets = 0;
      int idxBtag      = -10000;
      int AlgoBtag = sel.GetbtagAlgo();
      float btagDiscriCut = sel.GetbtagDiscriCut ();
      
      
      bool foundASelBjet = 0;
      for(unsigned int ijet = 0; ijet < selJets.size(); ijet++){
	if(abs(selJets[ijet].partonFlavour)==5 ){
	  foundASelBjet = true;
	}
	if ( AlgoBtag==0 &&  selJets[ijet].bTag["trackCountingHighEffBJetTags"]	      >= btagDiscriCut) NBtaggedJets++;
	if ( AlgoBtag==1 &&  selJets[ijet].bTag["simpleSecondaryVertexHighEffBJetTags"] >= btagDiscriCut) NBtaggedJets++;
	if ( AlgoBtag==2 &&  selJets[ijet].bTag["trackCountingHighPurBJetTags"]	      >= btagDiscriCut) NBtaggedJets++;
	if ( AlgoBtag==3 &&  selJets[ijet].bTag["simpleSecondaryVertexHighPurBJetTags"] >= btagDiscriCut) NBtaggedJets++;
	if ( AlgoBtag==4 &&  selJets[ijet].bTag["jetProbabilityBJetTags"] 	      >= btagDiscriCut) NBtaggedJets++;
	if ( AlgoBtag==5 &&  selJets[ijet].bTag["jetBProbabilityBJetTags"]	      >= btagDiscriCut) NBtaggedJets++;
	if ( AlgoBtag==6 &&  selJets[ijet].bTag["combinedSecondaryVertexBJetTags"]      >= btagDiscriCut){
	  NBtaggedJets++;
	  if(idxBtag < 0) idxBtag = ijet;
	}
      }  
            cout << __LINE__ << idxBtag  << endl;

      
      double btagdiscri = -10000;
      if(idxBtag>=0) btagdiscri = selJets[idxBtag].bTag["combinedSecondaryVertexBJetTags"];
      
            
      if(foundASelBjet) MyhistoManager.FillHisto(SelABjet, "SelABjet", 1.  , datasetName, IsSignal, 1.);
      else	       MyhistoManager.FillHisto(SelABjet, "SelABjet", 0.  , datasetName, IsSignal, 1. );
      
      
      if(idxBtag<0 && selJets.size() > 0) idxBtag = 0;
     
            
      
      int nvertex = selVertices.size();
            

      
      MyhistoManager.FillHisto(*pNvtx_afterleptsel,           ("Nvtx_"+leptChannel+"_afterleptsel").Data(),                nvertex, datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pInvM_ll_afterleptsel,        ("InvM_ll_"+leptChannel+"_afterleptsel").Data(),		dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pInvM_ll_afterleptsel_lowbin, ("InvM_ll_"+leptChannel+"_afterleptsel_lowbin").Data(),	dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pLeptZPt_afterleptsel,        ("LeptZPt_"+leptChannel+"_afterleptsel").Data(),		lept1.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pLeptZPt_afterleptsel,        ("LeptZPt_"+leptChannel+"_afterleptsel").Data(),		lept2.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(*pLeptWPt_afterleptsel,        ("LeptWPt_"+leptChannel+"_afterleptsel").Data(),		lept3.Pt(), datasetName, IsSignal, Dweight[ITypeMC]);
      
      
            

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
      
            

      //test for events with high ST
      if(sumPtLeptJet  > 300 && fabs(dilept.M()-91) > 15){
	
	if(isData){
	  /*cout << "******************************************************************* "<< endl;
	    cout << "lepton 1 pt/eta/phi " << lept1.Pt() << " " << lept1.Eta() << " " <<  lept1.Phi() << endl;
	    cout << "lepton 2 pt/eta/phi " << lept2.Pt() << " " << lept2.Eta() << " " <<  lept2.Phi() << endl;
	    cout << "lepton 3 pt/eta/phi " << lept3.Pt() << " " << lept3.Eta() << " " <<  lept3.Phi() << endl;
	    cout << "dilepton Z cand     " <<  dileptonIvM << endl; 
	    cout << "sum pT 	       " << sumPtLeptJet  << endl; 
	    cout << "MET		       " <<  theMET << endl; 
	    cout << "Njets  	       " <<  selJets.size() << endl; 
	    for(unsigned int i=0 ; i< selJets.size(); i++) cout << " jet " << i << " pt/eta/phi " << selJets[i].p4.Pt() << " " << selJets[i].p4.Eta() << " " <<  selJets[i].p4.Phi()
	    << "  btag discri " << selJets[i].GetDiscri(string("combinedSecondaryVertexBJetTags"))     << endl;
	  */
	}
	
	
        
	MyhistoManager.FillHisto(*pInvM_ll_afterleptsel_highSumPt, ("InvM_ll_"+leptChannel+"_afterleptsel_highSumPt").Data(), dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
      }  
      
      //signal region
      if( fabs(dilept.M()-91) < 15){
	
        //test for high Invarian masses
    	if(dilept.M() > 300){
	  /*cout << "dileptonIvM " << dileptonIvM
	    << "  dilept.Pt() "	 << dilept.Pt() 
	    << "  lept1.Pt() "	 << lept1.Pt() 
	    << "  lept2.Pt() "	 << lept2.Pt()
	    << "  lept3.Pt() "	 << lept3.Pt()
	    << "  Njet " << selJets.size();*/
	  for(unsigned int ijet=0; ijet<selJets.size(); ijet++) cout <<  "  jet " << ijet << " pt " << selJets[ijet].p4.Pt() << endl;
	  //cout << "  met	 " << met.p2.Mod() << endl;
	  
	  /*if(selJets.size()>0) cout << "  deltaphi lept-jet 1 " << lept1.DeltaR(selJets[0].p4) << 
	    "  deltaphi lept-jet 2 " << lept2.DeltaR(selJets[0].p4) <<
	    "  deltaphi lept-jet 3 " << lept3.DeltaR(selJets[0].p4) << endl;
	    
	    cout << "  deltaphi lept-met 1 " << lept1.Phi() - met.p2.Phi()  << 
	    "  deltaphi lept-met 2 " << lept2.Phi() - met.p2.Phi()  <<
	    "  deltaphi lept-met 3 " << lept3.Phi() - met.p2.Phi()  << endl;  
	  */
	}
	
	
	
      
      /*if(leptChannel == "mumumu"){
	
      cout << "------------------afterZmass cut--------------------------" << endl;
      cout << "Z[0] pt " << ZmumuCand[0].p4.Pt() << ", eta " <<  ZmumuCand[0].p4.Eta()  <<  
        ", iso   "  << sel.RelIso03PFDeltaBeta(ZmumuCand[0])  << 
        ", GM    "  << ZmumuCand[0].isGlobalMuon << 
        ", TM    "  << ZmumuCand[0].isTrackerMuon<< 
        ", C2    "  << ZmumuCand[0].Chi2	   << 
        ", NHit  "  << ZmumuCand[0].NTrValidHits  << 
        ", NVhit "  << ZmumuCand[0].NValidHits    << 
        ", DO    "  <<  fabs( ZmumuCand[0].D0Inner) << endl; 
      cout << "Z[1] pt " << ZmumuCand[1].p4.Pt() << ", eta " <<  ZmumuCand[1].p4.Eta()  <<  
        ", iso   "  << sel.RelIso03PFDeltaBeta(ZmumuCand[1])  << 
        ", GM    "  << ZmumuCand[1].isGlobalMuon << 
        ", TM    "  << ZmumuCand[1].isTrackerMuon<< 
        ", C2    "  << ZmumuCand[1].Chi2	   << 
        ", NHit  "  << ZmumuCand[1].NTrValidHits  << 
        ", NVhit "  << ZmumuCand[1].NValidHits    << 
        ", DO    "  << fabs(ZmumuCand[1].D0Inner) << endl; 
      cout << "W    pt " << WmuCand[0].p4.Pt() << ", eta " <<  WmuCand[0].p4.Eta()  <<  
        ", iso   "  << sel.RelIso03PFDeltaBeta(WmuCand[0])  << 
        ", GM    "  << WmuCand[0].isGlobalMuon << 
        ", TM    "  << WmuCand[0].isTrackerMuon<< 
        ", C2    "  << WmuCand[0].Chi2	   << 
        ", NHit  "  << WmuCand[0].NTrValidHits  << 
        ", NVhit "  << WmuCand[0].NValidHits    << 
        ", DO    "  << fabs(WmuCand[0].D0Inner) << endl; 
	cout << "Dweight[ITypeMC] " << Dweight[ITypeMC]  << endl;
	
      }*/
	
      
	
	
	MyhistoManager.FillHisto(*pCutFlow,     ("CutFlow_"+leptChannel).Data(),      3, datasetName, IsSignal, Dweight[ITypeMC]);
	MyhistoManager.FillHisto(*pErrCutFlow,	("ErrCutFlow_"+leptChannel).Data()  , 3, datasetName, IsSignal, EventYieldWeightError);
	//MyhistoManager.FillHisto(*pWmissAssing_afterleptsel , ("WmissAssing_"+leptChannel+"afterleptsel").Data(), leptonFlavor, datasetName, IsSignal, Dweight[ITypeMC]);
	
	
	
	double mTW = pow(
			 2*lept3.Pt()*met.p2.Mod()*(1-cos(lept3.Phi() - met.p2.Phi()))
			 ,0.5);
	
	
	
	MyhistoManager.FillHisto( *pmWT_afterleptsel, ("mWT_"+leptChannel+"_afterleptsel").Data(), mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	MyhistoManager.FillHisto2D(*pInvM_ll_vs_mWT_afterleptsel, ("InvM_ll_vs_mWT_"+leptChannel+"_afterleptsel").Data(),dilept.M(), mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	
	
	
	
	
	//test for events in the MWT tail
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
	  
	  MyhistoManager.FillHisto(*pLeptWPt_afterleptsel_mWT110, ("LeptWPt_"+leptChannel+"_afterleptsel_mWT110").Data(), lept3.Pt(),  datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(*pMET_afterleptsel_mWT110,	   ("MET_"+leptChannel+"_afterleptsel_mWT110").Data(),     theMET,	   datasetName, IsSignal, Dweight[ITypeMC]);  
	  MyhistoManager.FillHisto(*pInvM_ll_afterleptsel_mWT110, ("InvM_ll_"+leptChannel+"_afterleptsel_mWT110").Data(), dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
	  
	  for(unsigned int ijet=0; ijet<selJets.size(); ijet++){
	    MyhistoManager.FillHisto(*pdeltaRLeptJet_afterleptsel_mWT110, ("deltaRLeptJet_"+leptChannel+"_afterleptsel_mWT110").Data(), lept1.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pdeltaRLeptJet_afterleptsel_mWT110, ("deltaRLeptJet_"+leptChannel+"_afterleptsel_mWT110").Data(), lept2.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pdeltaRLeptJet_afterleptsel_mWT110, ("deltaRLeptJet_"+leptChannel+"_afterleptsel_mWT110").Data(), lept3.DeltaR(selJets[ijet].p4), datasetName, IsSignal, Dweight[ITypeMC]);
	  }
	  
	
	
	  MyhistoManager.FillHisto(*pCharge_afterleptsel_mWT110, ("Charge_"+leptChannel+"_afterleptsel_mWT110").Data(), wcharge, datasetName, IsSignal, Dweight[ITypeMC]);
	  
	  MyhistoManager.FillHisto(*pdeltaRLeptMet_afterleptsel_mWT110, ("deltaRLeptMet_"+leptChannel+"_afterleptsel_mWT110").Data(),deltaPhi1 , datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(*pdeltaRLeptMet_afterleptsel_mWT110, ("deltaRLeptMet_"+leptChannel+"_afterleptsel_mWT110").Data(),deltaPhi2 , datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(*pdeltaRLeptMet_afterleptsel_mWT110, ("deltaRLeptMet_"+leptChannel+"_afterleptsel_mWT110").Data(),deltaPhi3 , datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(*pNJet_afterleptsel_mWT110  ,        ("NJet_"+leptChannel+"_afterleptsel_mWT110").Data() , njetselec,	 datasetName, IsSignal, Dweight[ITypeMC]);    
	  MyhistoManager.FillHisto(*pNBJet_afterleptsel_mWT110 ,        ("NBJet_"+leptChannel+"_afterleptsel_mWT110").Data(), NBtaggedJets,      datasetName, IsSignal, Dweight[ITypeMC]); 
	  
	} // outmTW > 110
	
	
	
	
	  
        //*****************************************************************
        // pass jet selection
        //***************************************************************** 
      	int NSeljets = selJets.size() ;
	if(NSeljets>4) NSeljets = 4;

	
  	if(selJets.size() >= 1 ){
	  
	  MyhistoManager.FillHisto(*pCutFlow,      ("CutFlow_"+leptChannel).Data(),    4, datasetName, IsSignal, Dweight[ITypeMC]);    
          MyhistoManager.FillHisto(*pErrCutFlow,   ("ErrCutFlow_"+leptChannel).Data(), 4, datasetName, IsSignal, EventYieldWeightError);
	  
	  MyhistoManager.FillHisto( *pmWT_afterjetsel, ("mWT_"+leptChannel+"_afterjetsel").Data(), mTW, datasetName, IsSignal, Dweight[ITypeMC]);
	
	  MyhistoManager.FillHisto(*pInvM_ll_afterjetsel, ("InvM_ll_"+leptChannel+"_afterjetsel").Data(), dileptonIvM, datasetName, IsSignal, Dweight[ITypeMC]);
          
	  
	  if(isData) MyhistoManager.FillHisto(*pNBJet_afterjetsel ,  ("NBJet_"+leptChannel+"_afterjetsel").Data(),NBtaggedJets, datasetName, IsSignal, Dweight[ITypeMC]);
          else{
	    MyhistoManager.FillHisto(*pNBJet_afterjetsel ,  ("NBJet_"+leptChannel+"_afterjetsel").Data(),0, datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(*pNBJet_afterjetsel ,  ("NBJet_"+leptChannel+"_afterjetsel").Data(),1, datasetName, IsSignal, Dweight[ITypeMC]); 
	
	    	   
	  }
	  
	
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

	  TLorentzVector topCand = neutrino + lept3 + selJets[0].p4 ;
	  
	
	  
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
	  
	    
	
	  if(
	    (NBtaggedJets >=0  &&  NBtaggedJets <=1 && isData && !useNonIsoWcand ) ||
	    (NBtaggedJets >=0  &&  NBtaggedJets <=1 && isData && useNonIsoWcand && theMET < themetcut)
	    
	    
	    ){

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
	    topCand = neutrino + lept3 + selJets[idxBtag].p4 ;
	    

	    
	  } //bjet selection
	  
	  
	  

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
	  
	    
	  tree_totPt	    = (topCand + (lept1+lept2)).Pt();
	  tree_totEta	    = (topCand + (lept1+lept2)).Eta();

  	  tree_deltaRZl     = (lept1+lept2).DeltaR(lept3);
  	  tree_deltaPhiZmet = (lept1+lept2).DeltaPhi(metP4);
  	  tree_btagDiscri   = btagdiscri;
	    
  	  tree_NJets	    = float(selJets.size());
	  
	
	  
	  /*if(!isData){
	    double thrand = rand.Uniform();
	    NBtaggedJets = 0;
	    if( thrand > weightb[1] ) NBtaggedJets = 1;
	
	  }*/
	  
	  tree_NBJets	    = float(NBtaggedJets);
	  
  	  tree_leptWPt        = lept3.Pt(); 
	  tree_leptWEta       = lept3.Eta();
  	  tree_leadJetPt      = selJets[idxBtag].p4.Pt(); 
	  tree_leadJetEta     = selJets[idxBtag].p4.Eta();
          tree_deltaRZleptW   = (lept1+lept2).DeltaR(lept3); 
	  tree_deltaPhiZleptW = (lept1+lept2).DeltaPhi(lept3);
	  
	  tree_met = theMET;
	  tree_mTW = mTW;
	  
	  
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
	  tree_Channel = 0; TheTree->Fill();
	  
	
	  
	  MyhistoManager.FillHisto(*pAsym_afterbjetsel,	       ("Asym_"+leptChannel+"_afterbjetsel").Data(),         tree_asym,    datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(*pRecoPtZ_afterbjetsel,     ("RecoPtZ_"+leptChannel+"_afterbjetsel").Data(),     tree_Zpt,     datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(*pRecoTopMass_afterbjetsel, ("RecoTopMass_"+leptChannel+"_afterbjetsel").Data(), tree_topMass ,datasetName, IsSignal, Dweight[ITypeMC]);
          MyhistoManager.FillHisto(*pdeltaPhilb_afterbjetsel , ("deltaPhilb_"+leptChannel+"_afterbjetsel").Data(),  tree_deltaPhilb ,    datasetName, IsSignal, Dweight[ITypeMC]);

	  
	  
	  
	  
	
	  
	  
	  

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
    
  cout << "selected mumumu events " << nselevents_mumumu << endl;

  
   //The following line is mandatory to copy everything in a common RootFile
    fOutput->Add(fProofFile);

    cleanHistoVector();
    
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





//_____________________________________________________________________________
std::vector<double> ProofSelectorMyCutFlow::determineWeights(TString thedatasetName, double weightITypeMC, double WeightForBranchingRatio){
  
  
  
  //------------------------------------
  //to calculate the event weight
  //returns a vector of double,
  //containing the infor for reweighting
  //------------------------------------

      
    double ITypeMC = 0;
    double Dweight = 0 ;
    double EventYieldWeightError = 0 ;
    double IsSignal = -1 ;
    
    Dweight= weightITypeMC * WeightForBranchingRatio;
    EventYieldWeightError = Dweight*Dweight;
    
      
    if ( thedatasetName=="TTbar" ) {
      ITypeMC = -1;
      Dweight= weightITypeMC * WeightForBranchingRatio;
      EventYieldWeightError = Dweight*Dweight;
    }
    else {
      IsSignal = -1;
      Dweight= weightITypeMC;
      EventYieldWeightError = Dweight*Dweight;
      if ( thedatasetName=="Zjets" || thedatasetName=="DYToLL_M10-50") ITypeMC = 2;
      else if ( thedatasetName=="Wjets" ) ITypeMC = 3;
      else if ( thedatasetName=="SingleToptW" || thedatasetName=="TtW" || thedatasetName=="TbartW"
		|| thedatasetName=="TtWScaleUp" || thedatasetName=="TtWScaleDown"
		|| thedatasetName=="TbartWScaleUp" || thedatasetName=="TbartWScaleDown") ITypeMC = 4;
      else if ( thedatasetName=="WZ" || thedatasetName=="WW" || thedatasetName=="ZZ" 
             || thedatasetName=="WZ_scaleup"|| thedatasetName=="WZ_scaledown"  
             || thedatasetName=="WZ_matchup"|| thedatasetName=="WZ_matchdown" ) ITypeMC = 5;
      else if ( 
       		  thedatasetName=="FCNCkut" || thedatasetName=="FCNCkct" 
     		||thedatasetName=="FCNCxut" || thedatasetName=="FCNCxct" 
     		||thedatasetName=="FCNCzut" || thedatasetName=="FCNCzct" 
     
     		||thedatasetName=="FCNCkutFullSim" || thedatasetName=="FCNCkctFullSim" 
     		||thedatasetName=="FCNCxutFullSim" || thedatasetName=="FCNCxctFullSim" 
    		||thedatasetName=="FCNCzutFullSim" || thedatasetName=="FCNCzctFullSim" 
     
     		||thedatasetName=="FCNCkut_matchup" || thedatasetName=="FCNCkct_matchup" 
     		||thedatasetName=="FCNCxut_matchup" || thedatasetName=="FCNCxct_matchup" 
     		||thedatasetName=="FCNCzut_matchup" || thedatasetName=="FCNCzct_matchup" 
     
     		||thedatasetName=="FCNCkut_matchdown" || thedatasetName=="FCNCkct_matchdown" 
     		||thedatasetName=="FCNCxut_matchdown" || thedatasetName=="FCNCxct_matchdown" 
     		||thedatasetName=="FCNCzut_matchdown" || thedatasetName=="FCNCzct_matchdown" 
     
     		||thedatasetName=="FCNCkut_scaleup" || thedatasetName=="FCNCkct_scaleup" 
    		||thedatasetName=="FCNCxut_scaleup" || thedatasetName=="FCNCxct_scaleup" 
     		||thedatasetName=="FCNCzut_scaleup" || thedatasetName=="FCNCzct_scaleup" 
     
    		||thedatasetName=="FCNCkut_scaledown" || thedatasetName=="FCNCkct_scaledown" 
     		||thedatasetName=="FCNCxut_scaledown" || thedatasetName=="FCNCxct_scaledown" 
     		||thedatasetName=="FCNCzut_scaledown" || thedatasetName=="FCNCzct_scaledown" 
     
     		||thedatasetName=="FCNCkut_topup" || thedatasetName=="FCNCkct_topup" 
     		||thedatasetName=="FCNCxut_topup" || thedatasetName=="FCNCxct_topup" 
     		||thedatasetName=="FCNCzut_topup" || thedatasetName=="FCNCzct_topup" 
     
     		||thedatasetName=="FCNCkut_topdown" || thedatasetName=="FCNCkct_topdown" 
     		||thedatasetName=="FCNCxut_topdown" || thedatasetName=="FCNCxct_topdown" 
   
      )  ITypeMC = 6; 
      else if(thedatasetName=="tZq") ITypeMC = 7;
      else if(thedatasetName=="TT" )ITypeMC = 8;
   }
  
  if ( thedatasetName=="DataDiEG" || thedatasetName=="DataDiMu" ||
       thedatasetName=="DataDiMuEG" || thedatasetName=="DataDiEGMu") {
     ITypeMC = 100;
     Dweight= weightITypeMC;
     EventYieldWeightError = Dweight*Dweight;
  }




  
  std::vector< double > thereturn;
  thereturn.push_back(ITypeMC);
  thereturn.push_back(Dweight);
  thereturn.push_back(EventYieldWeightError);
  thereturn.push_back(IsSignal);
  
  
  return thereturn;
}






 
void  ProofSelectorMyCutFlow::determineLeptonCandidates(
  		bool UseLooseWcand, float looseIsoCut, double rhocorr,
  		std::vector<NTElectron> *selE,        std::vector<NTMuon> *selM, 
  		std::vector<NTElectron> *selENonIso,  std::vector<NTMuon> *selMNonIso, 
		std::vector<NTElectron> *theZeeCand,  std::vector<NTMuon> *theZmumuCand, 
		std::vector<NTElectron> *theWeCand,   std::vector<NTMuon> *theWmuCand
		){
 
      //*****************************************************************
      // select Z->ee candidate
      //*****************************************************************    
     
  int leptonFlavor = 0;
  int wcharge	   = 0;
  
  
  theZeeCand->clear(); 
  theZmumuCand->clear(); 
  
  theWeCand->clear(); 
  theWmuCand->clear(); 

  if(selE->size() >=2 ) {
    int theel1 = -1;
    int theel2 = -1;
    double mInv = 1000000;
    for(unsigned int iel1 = 0; iel1 < selE->size(); iel1++){
      for(unsigned int iel2 = 0; iel2 < selE->size(); iel2++){
    	if(iel1 == iel2) continue;
    	if((*selE)[iel1].charge == (*selE)[iel2].charge) continue;
    	TLorentzVector theZee = (*selE)[iel1].p4 + (*selE)[iel2].p4;
    	if( fabs(theZee.M() - 91) < fabs(mInv-91) ){
    	//if( fabs(theZeeCand.M() - 200) < fabs(mInv-200) ){
    	  theel1 = iel1;
    	  theel2 = iel2;
    	  mInv = theZee.M();
    	}
      }
    }

    if(theel1>=0 && theel2>=0){ //JLA
      theZeeCand->push_back((*selE)[theel1]);
      theZeeCand->push_back((*selE)[theel2]);
      //double invM = (theZeeCand[0].p4+theZeeCand[1].p4).M();
      //cout << "lepton origin of Zee cand " << selE[theel1].LeptonOrigin << "  " << selE[theel2].LeptonOrigin << "  invmass " <<  invM << endl;
    }
  }
     
  //*****************************************************************
  // select W->enu candidate
  //*****************************************************************	 
  //cout << "get lepton cand W->enu " << endl;
  if(!UseLooseWcand){
    for(unsigned int iel1 = 0; iel1 < selE->size(); iel1++){
      bool matchElec=false;
      for(unsigned int iel2 = 0; iel2 < theZeeCand->size(); iel2++){
  	 
    	 if(fabs((*selE)[iel1].p4.Pt() - (*theZeeCand)[iel2].p4.Pt()) <  0.0001)  matchElec=true;
       }
      if(!matchElec){
    	theWeCand->push_back((*selE)[iel1]);
    	wcharge = (*selE)[iel1].charge;
    	if( (*selE)[iel1].LeptonOrigin == 10) leptonFlavor = 1;
      }
    }
  }else{
    for(unsigned int iel1 = 0; iel1 < selENonIso->size(); iel1++){
      bool matchElec=false;
      for(unsigned int iel2 = 0; iel2 < theZeeCand->size(); iel2++){
    	 if(fabs( (*selENonIso)[iel1].p4.Pt() - (*theZeeCand)[iel2].p4.Pt()) <  0.0001)  matchElec=true;
    	 else if( (*selE)[iel1].LeptonOrigin == 10) leptonFlavor = 1;
      }
      //if(!matchElec && selENonIso[iel1].RelIso03PF() > looseIsoCut){ 
      if(!matchElec && sel.EffArea03PF( (*selENonIso)[iel1], rhocorr) > looseIsoCut){ 
    	theWeCand->push_back((*selENonIso)[iel1]);
    	wcharge = (*selENonIso)[iel1].charge; 
    	if( (*selENonIso)[iel1].LeptonOrigin == 10) leptonFlavor = 1;
    	//cout << "    lepton origin of We cand " << selENonIso[iel1].LeptonOrigin << endl;
    	//if(selENonIso[iel1].LeptonOrigin == 10) cout << " nofake wenu" << endl;
      }
    }
  }
  
  
  //*****************************************************************
  // select Z->mumu candidate
  //*****************************************************************  
  if(selM->size() >=2 ) {
    int themu1 = -1;
    int themu2 = -1;
    double mInv = 1000000;
    for(unsigned int imu1 = 0; imu1 < selM->size(); imu1++){
      for(unsigned int imu2 = 0; imu2 < selM->size(); imu2++){
    	if(imu1 == imu2) continue;
    	if((*selM)[imu1].charge == (*selM)[imu2].charge) continue;
    	TLorentzVector theZmumu = (*selM)[imu1].p4 + (*selM)[imu2].p4;
    	if( fabs(theZmumu.M() - 91) < fabs(mInv-91) ){
    	//if( fabs(theZmumuCand.M() - 200) < fabs(mInv-200) ){
    	  themu1 = imu1;
    	  themu2 = imu2;
    	  mInv = theZmumu.M();
    	}
      }
    }
    if(themu1>=0 && themu2>=0){ //JLA  
      theZmumuCand->push_back((*selM)[themu1]);
      theZmumuCand->push_back((*selM)[themu2]);
    }
     
  }
  
   
   
  //*****************************************************************
  // select W->munu candidate
  //*****************************************************************	 
  
  //cout << "get lepton cand W->munu " << endl;
  if(!UseLooseWcand){
    //cout << "in sel W " << endl;
    for(unsigned int imu1 = 0; imu1 < selM->size(); imu1++){
      bool matchMuon = false;
      for(unsigned int imu2 = 0; imu2 < theZmumuCand->size(); imu2++){
  	 
    	 if(fabs((*selM)[imu1].p4.Pt() -(* theZmumuCand)[imu2].p4.Pt()) <  0.0001) matchMuon = true;
    	 
      } 
      if(!matchMuon){
    	theWmuCand->push_back( (*selM)[imu1]);
    	wcharge = (*selM)[imu1].charge;
    	if( (*selM)[imu1].LeptonOrigin == 10) leptonFlavor = 1;
      }
    }
  }else{
    for(unsigned int imu1 = 0; imu1 < selMNonIso->size(); imu1++){
      bool matchMuon = false;
      for(unsigned int imu2 = 0; imu2 < theZmumuCand->size(); imu2++){
  	 
    	 if(fabs( (*selMNonIso)[imu1].p4.Pt() - (*theZmumuCand)[imu2].p4.Pt())  < 0.0001) matchMuon = true;
     
      }
      //if(!matchMuon && selMNonIso[imu1].RelIso03PF() > looseIsoCut){
      if(!matchMuon && sel.RelIso03PFDeltaBeta( (*selMNonIso)[imu1]) > looseIsoCut){
    	theWmuCand->push_back( (*selMNonIso)[imu1]);
    	wcharge = (*selMNonIso)[imu1].charge;
    	if( (*selMNonIso)[imu1].LeptonOrigin == 10) leptonFlavor = 1;
      }
    }
  }
  
  
  /*
  //redefine the lepton coll for jet cleaning
  if(UseLooseWcand){
     
    vector<NTElectron>  tmpElectrons = sel.GetSelectedElectronsNoIso();
    vector<NTMuon>	tmpMuons     = sel.GetSelectedMuonsNoIso();
    
    for(unsigned int iel=0; iel<tmpElectrons.size(); iel++){
      if(sel.EffArea03PF(tmpElectrons[iel], rho ) > 0.4) selE->push_back(tmpElectrons[iel]);
    }
  
    for(unsigned int imu=0; imu<tmpMuons.size(); imu++){
      if(sel.RelIso03PFDeltaBeta(tmpMuons[imu]) > 0.4) selM->push_back(tmpMuons[imu]);
    }
  
  }else{
  
    selE = sel.GetSelectedElectrons();
    selM	= sel.GetSelectedMuons();
  }*/
  
}
















void ProofSelectorMyCutFlow::createTheHisto(HistoManager *thehistomanag){
  
  
  //create all histograms
  thehistomanag->CreateHisto(SelABjet, "SelABjet", datasetName, "Nevents","Entries", 2, -0.5, 1.5);
  thehistomanag->CreateHisto(SelABjet_afterjetsel, "SelABjet_afterjetsel", datasetName, "Nevents","Entries", 2, -0.5, 1.5);
 
  thehistomanag->CreateHisto(Ntrilept_mumumu, "Ntrilept_mumumu", datasetName, "Nevents","Entries", 11, -0.5, 10.5); 
  thehistomanag->CreateHisto(Ntrileptnoniso_mumumu, "Ntrileptnoniso_mumumu", datasetName, "Nevents","Entries", 11, -0.5, 10.5); 
  
  
  thehistomanag->CreateHisto(Nvertex, "Nvertex", datasetName, "Nvertex", "Entries", 50, 0, 50); 
  
  thehistomanag->CreateHisto(CutFlow_mumumu,  "CutFlow_mumumu" ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  thehistomanag->CreateHisto(CutFlow_mumue,   "CutFlow_mumue"  ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  thehistomanag->CreateHisto(CutFlow_eemu,    "CutFlow_eemu"   ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  thehistomanag->CreateHisto(CutFlow_eee,     "CutFlow_eee"    ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  
  
  thehistomanag->SetCutFlowAxisTitleFCNCMonotop(CutFlow_mumumu,   "CutFlow_mumumu"  ,datasetName);
  thehistomanag->SetCutFlowAxisTitleFCNCMonotop(CutFlow_mumue,    "CutFlow_mumue"   ,datasetName);
  thehistomanag->SetCutFlowAxisTitleFCNCMonotop(CutFlow_eemu,     "CutFlow_eemu"    ,datasetName);
  thehistomanag->SetCutFlowAxisTitleFCNCMonotop(CutFlow_eee,      "CutFlow_eee"     ,datasetName);
  
  
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
  
  thehistomanag->CreateHisto(NBJet_mumumu_afterjetsel, "NBJet_mumumu_afterjetsel", datasetName,"NBjets", "Entries",2,-0.5,1.5);
  thehistomanag->CreateHisto(NBJet_mumue_afterjetsel , "NBJet_mumue_afterjetsel" , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  thehistomanag->CreateHisto(NBJet_eemu_afterjetsel  , "NBJet_eemu_afterjetsel"  , datasetName,"NBjets", "Entries",2,-0.5,1.5);
  thehistomanag->CreateHisto(NBJet_eee_afterjetsel   , "NBJet_eee_afterjetsel"   , datasetName,"NBjets", "Entries",2,-0.5,1.5);	
  
  
  
  
  
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
  
  
  
  
  
  /*thehistomanag->CreateHisto(NBJet_mumumu_afterjetsel, "NBJet_mumumu_afterjetsel", datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_mumue_afterjetsel , "NBJet_mumue_afterjetsel" , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_eemu_afterjetsel  , "NBJet_eemu_afterjetsel"  , datasetName,"NBjets", "Entries",5,-0.5,4.5);
  thehistomanag->CreateHisto(NBJet_eee_afterjetsel   , "NBJet_eee_afterjetsel"   , datasetName,"NBjets", "Entries",5,-0.5,4.5); 
  */
   
  
  thehistomanag->CreateHisto(InvM_ll_mumumu_afterleptsel_highSumPt, "InvM_ll_mumumu_afterleptsel_highSumPt" , datasetName,"Minv", "Entries",100,0.,1000);                     
  thehistomanag->CreateHisto(InvM_ll_mumue_afterleptsel_highSumPt,  "InvM_ll_mumue_afterleptsel_highSumPt"  , datasetName,"Minv", "Entries",100,0.,1000);
  thehistomanag->CreateHisto(InvM_ll_eemu_afterleptsel_highSumPt,   "InvM_ll_eemu_afterleptsel_highSumPt"   , datasetName,"Minv", "Entries",100,0.,1000);
  thehistomanag->CreateHisto(InvM_ll_eee_afterleptsel_highSumPt,    "InvM_ll_eee_afterleptsel_highSumPt"    , datasetName,"Minv", "Entries",100,0.,1000);

 
  thehistomanag->CreateHisto(InvM_ll_mumumu_afterleptsel, "InvM_ll_mumumu_afterleptsel" , datasetName,"Minv", "Entries",350,0.,1000);
  thehistomanag->CreateHisto(InvM_ll_mumue_afterleptsel,  "InvM_ll_mumue_afterleptsel"  , datasetName,"Minv", "Entries",350,0.,1000);
  thehistomanag->CreateHisto(InvM_ll_eemu_afterleptsel,   "InvM_ll_eemu_afterleptsel"   , datasetName,"Minv", "Entries",350,0.,1000);
  thehistomanag->CreateHisto(InvM_ll_eee_afterleptsel,    "InvM_ll_eee_afterleptsel"    , datasetName,"Minv", "Entries",350,0.,1000);
  
  thehistomanag->CreateHisto(InvM_ll_mumumu_afterleptsel_mWT110, "InvM_ll_mumumu_afterleptsel_mWT110" , datasetName,"Minv", "Entries",350,0.,1000);
  thehistomanag->CreateHisto(InvM_ll_mumue_afterleptsel_mWT110,  "InvM_ll_mumue_afterleptsel_mWT110"  , datasetName,"Minv", "Entries",350,0.,1000);
  thehistomanag->CreateHisto(InvM_ll_eemu_afterleptsel_mWT110,   "InvM_ll_eemu_afterleptsel_mWT110"   , datasetName,"Minv", "Entries",350,0.,1000);
  thehistomanag->CreateHisto(InvM_ll_eee_afterleptsel_mWT110,    "InvM_ll_eee_afterleptsel_mWT110"    , datasetName,"Minv", "Entries",350,0.,1000);
 
  thehistomanag->CreateHisto(InvM_ll_mumumu_afterleptsel_lowbin, "InvM_ll_mumumu_afterleptsel_lowbin" , datasetName,"Minv", "Entries",100,0.,200);
  thehistomanag->CreateHisto(InvM_ll_mumue_afterleptsel_lowbin,  "InvM_ll_mumue_afterleptsel_lowbin"  , datasetName,"Minv", "Entries",100,0.,200);
  thehistomanag->CreateHisto(InvM_ll_eemu_afterleptsel_lowbin,   "InvM_ll_eemu_afterleptsel_lowbin"   , datasetName,"Minv", "Entries",100,0.,200);
  thehistomanag->CreateHisto(InvM_ll_eee_afterleptsel_lowbin,    "InvM_ll_eee_afterleptsel_lowbin"    , datasetName,"Minv", "Entries",100,0.,200);

  thehistomanag->CreateHisto(InvM_ll_mumumu_afterjetsel, "InvM_ll_mumumu_afterjetsel" , datasetName,"Minv", "Entries",350, 0., 1000);
  thehistomanag->CreateHisto(InvM_ll_mumue_afterjetsel,  "InvM_ll_mumue_afterjetsel"  , datasetName,"Minv", "Entries",350, 0., 1000);
  thehistomanag->CreateHisto(InvM_ll_eemu_afterjetsel,   "InvM_ll_eemu_afterjetsel"   , datasetName,"Minv", "Entries",350, 0., 1000);
  thehistomanag->CreateHisto(InvM_ll_eee_afterjetsel,    "InvM_ll_eee_afterjetsel"    , datasetName,"Minv", "Entries",350, 0., 1000);
  
  thehistomanag->CreateHisto(InvM_ll_mumumu_afterbjetsel, "InvM_ll_mumumu_afterbjetsel" , datasetName,"Minv", "Entries",350,0.,1000);
  thehistomanag->CreateHisto(InvM_ll_mumue_afterbjetsel,  "InvM_ll_mumue_afterbjetsel"  , datasetName,"Minv", "Entries",350,0.,1000);
  thehistomanag->CreateHisto(InvM_ll_eemu_afterbjetsel,   "InvM_ll_eemu_afterbjetsel"   , datasetName,"Minv", "Entries",350,0.,1000);
  thehistomanag->CreateHisto(InvM_ll_eee_afterbjetsel,    "InvM_ll_eee_afterbjetsel"    , datasetName,"Minv", "Entries",350,0.,1000);

    
  
  thehistomanag->CreateHisto(LeptPt_mumumu_afterleptsel, "LeptPt_mumumu_afterleptsel", datasetName, "LeptPt", "Entries",350,0., 1000) ;
  thehistomanag->CreateHisto(LeptPt_mumue_afterleptsel,  "LeptPt_mumue_afterleptsel",  datasetName, "LeptPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptPt_eemu_afterleptsel,   "LeptPt_eemu_afterleptsel",   datasetName, "LeptPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptPt_eee_afterleptsel,    "LeptPt_eee_afterleptsel",    datasetName, "LeptPt", "Entries",350,0., 1000);
  
  thehistomanag->CreateHisto(LeptPt_mumumu_afterjetsel, "LeptPt_mumumu_afterjetsel", datasetName, "LeptPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptPt_mumue_afterjetsel,  "LeptPt_mumue_afterjetsel",  datasetName, "LeptPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptPt_eemu_afterjetsel,   "LeptPt_eemu_afterjetsel",   datasetName, "LeptPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptPt_eee_afterjetsel,    "LeptPt_eee_afterjetsel",    datasetName, "LeptPt", "Entries",350,0., 1000);
  
  thehistomanag->CreateHisto(LeptPt_mumumu_afterbjetsel, "LeptPt_mumumu_afterbjetsel", datasetName, "LeptPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptPt_mumue_afterbjetsel,  "LeptPt_mumue_afterbjetsel",  datasetName, "LeptPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptPt_eemu_afterbjetsel,   "LeptPt_eemu_afterbjetsel",   datasetName, "LeptPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(LeptPt_eee_afterbjetsel,    "LeptPt_eee_afterbjetsel",    datasetName, "LeptPt", "Entries",350,0., 1000);
    
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
  

  
  
  
  thehistomanag->CreateHisto(JetPt_mumumu_afterleptsel, "JetPt_mumumu_afterleptsel", datasetName, "JetPt", "Entries",350,0., 1000) ;
  thehistomanag->CreateHisto(JetPt_mumue_afterleptsel,  "JetPt_mumue_afterleptsel",  datasetName, "JetPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(JetPt_eemu_afterleptsel,   "JetPt_eemu_afterleptsel",   datasetName, "JetPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(JetPt_eee_afterleptsel,    "JetPt_eee_afterleptsel",    datasetName, "JetPt", "Entries",350,0., 1000);
  
  thehistomanag->CreateHisto(JetPt_mumumu_afterjetsel, "JetPt_mumumu_afterjetsel", datasetName, "JetPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(JetPt_mumue_afterjetsel,  "JetPt_mumue_afterjetsel",  datasetName, "JetPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(JetPt_eemu_afterjetsel,   "JetPt_eemu_afterjetsel",   datasetName, "JetPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(JetPt_eee_afterjetsel,    "JetPt_eee_afterjetsel",    datasetName, "JetPt", "Entries",350,0., 1000);
  
  thehistomanag->CreateHisto(JetPt_mumumu_afterbjetsel, "JetPt_mumumu_afterbjetsel", datasetName, "JetPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(JetPt_mumue_afterbjetsel,  "JetPt_mumue_afterbjetsel",  datasetName, "JetPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(JetPt_eemu_afterbjetsel,   "JetPt_eemu_afterbjetsel",   datasetName, "JetPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(JetPt_eee_afterbjetsel,    "JetPt_eee_afterbjetsel",    datasetName, "JetPt", "Entries",350,0., 1000);
  
  thehistomanag->CreateHisto(JetPt_mumumu_afterbjetveto, "JetPt_mumumu_afterbjetveto", datasetName, "JetPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(JetPt_mumue_afterbjetveto,  "JetPt_mumue_afterbjetveto",  datasetName, "JetPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(JetPt_eemu_afterbjetveto,   "JetPt_eemu_afterbjetveto",   datasetName, "JetPt", "Entries",350,0., 1000);
  thehistomanag->CreateHisto(JetPt_eee_afterbjetveto,    "JetPt_eee_afterbjetveto",    datasetName, "JetPt", "Entries",350,0., 1000);
  
    
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
  
  
  
  thehistomanag->CreateHisto(MET_mumumu_afterleptsel, "MET_mumumu_afterleptsel", datasetName, "MET", "Entries",100,0., 500);
  thehistomanag->CreateHisto(MET_mumue_afterleptsel,  "MET_mumue_afterleptsel",  datasetName, "MET", "Entries",100,0., 500);
  thehistomanag->CreateHisto(MET_eemu_afterleptsel,   "MET_eemu_afterleptsel",   datasetName, "MET", "Entries",100,0., 500);
  thehistomanag->CreateHisto(MET_eee_afterleptsel,    "MET_eee_afterleptsel",    datasetName, "MET", "Entries",100,0., 500);
  
  thehistomanag->CreateHisto(MET_mumumu_afterleptsel_mWT110, "MET_mumumu_afterleptsel_mWT110", datasetName, "MET", "Entries",100,0., 500);
  thehistomanag->CreateHisto(MET_mumue_afterleptsel_mWT110,  "MET_mumue_afterleptsel_mWT110",  datasetName, "MET", "Entries",100,0., 500);
  thehistomanag->CreateHisto(MET_eemu_afterleptsel_mWT110,   "MET_eemu_afterleptsel_mWT110",   datasetName, "MET", "Entries",100,0., 500);
  thehistomanag->CreateHisto(MET_eee_afterleptsel_mWT110,    "MET_eee_afterleptsel_mWT110",    datasetName, "MET", "Entries",100,0., 500);
  
  
  
  
  thehistomanag->CreateHisto(MET_mumumu_afterjetsel, "MET_mumumu_afterjetsel", datasetName, "MET", "Entries",100,0., 500);
  thehistomanag->CreateHisto(MET_mumue_afterjetsel,  "MET_mumue_afterjetsel",  datasetName, "MET", "Entries",100,0., 500);
  thehistomanag->CreateHisto(MET_eemu_afterjetsel,   "MET_eemu_afterjetsel",   datasetName, "MET", "Entries",100,0., 500);
  thehistomanag->CreateHisto(MET_eee_afterjetsel,    "MET_eee_afterjetsel",    datasetName, "MET", "Entries",100,0., 500);
  
  thehistomanag->CreateHisto(MET_mumumu_afterbjetsel, "MET_mumumu_afterbjetsel", datasetName, "MET", "Entries",100,0., 500);
  thehistomanag->CreateHisto(MET_mumue_afterbjetsel,  "MET_mumue_afterbjetsel",  datasetName, "MET", "Entries",100,0., 500);
  thehistomanag->CreateHisto(MET_eemu_afterbjetsel,   "MET_eemu_afterbjetsel",   datasetName, "MET", "Entries",100,0., 500);
  thehistomanag->CreateHisto(MET_eee_afterbjetsel,    "MET_eee_afterbjetsel",    datasetName, "MET", "Entries",100,0., 500);
    
   
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
    pNBJet_afterleptsel_mWT110      = &NBJet_mumue_afterleptsel_mWT110;
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
    pNBJet_afterleptsel_mWT110      = &NBJet_eemu_afterleptsel_mWT110;
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
    pNBJet_afterleptsel_mWT110      = &NBJet_eee_afterleptsel_mWT110;
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





void  ProofSelectorMyCutFlow::WriteTheHisto(TFile* theoutputfile, HistoManager *thehistomanag){
  
  
  theoutputfile->cd();
  
  thehistomanag->WriteMyHisto(SelABjet, "all");
  thehistomanag->WriteMyHisto(SelABjet_afterjetsel, "all");
  thehistomanag->WriteMyHisto(Ntrilept_mumumu, "all");
  thehistomanag->WriteMyHisto(Ntrileptnoniso_mumumu, "all");
  
   
  thehistomanag->WriteMyHisto(Nvertex, "all");  
  
  thehistomanag->WriteMyHisto(CutFlow_mumumu, "all" );
  thehistomanag->WriteMyHisto(CutFlow_mumue,  "all" );
  thehistomanag->WriteMyHisto(CutFlow_eemu,   "all" );
  thehistomanag->WriteMyHisto(CutFlow_eee,    "all" );
  
  
  
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
  
  
  
  
  thehistomanag->WriteMyHisto(NBJet_mumumu_afterleptsel_mWT110, "all");
  thehistomanag->WriteMyHisto(NBJet_mumue_afterleptsel_mWT110,  "all");
  thehistomanag->WriteMyHisto(NBJet_eemu_afterleptsel_mWT110,   "all");
  thehistomanag->WriteMyHisto(NBJet_eee_afterleptsel_mWT110,    "all");
  
  
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
  
  
}








