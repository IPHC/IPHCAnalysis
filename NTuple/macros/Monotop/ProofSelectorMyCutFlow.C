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
   
  
  
  applyFakescale = true;
  applyLeptonSF  = true;
  applyWZ        = true;
  applyTrigger   = true;
  
  
  IReweight		= true;
  IDYestimateWithMetCut = true;
  IReweight_puUp	= false;
  IReweight_puDown	= false;
  
  
  useNonIsoWcand = false;
  
  rand.SetSeed(102994949221);

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
  if(tree == 0) cout << "TTree null pointer " << endl; 
  // Set branch addresses
  cout << "GetNbranches() " << tree->GetNbranches() << endl;
  //branch = (TBranch *) tree->GetBranch("NTEvent");
  branch = (TBranch *) tree->FindBranch("NTEvent");
  
  cout << "test number of branches " << tree->GetListOfBranches()->GetEntries()  << endl;
  cout << "test branch list        " << tree->GetListOfBranches()->First()->GetName() << endl;
  cout << "file directory is " << tree->GetDirectory()->GetPath() << endl;
  
  if(branch == 0) cout << "TBranch null pointer " << endl; 
  cout << "construct NTevent pointer " << endl;
  event = new IPHCTree::NTEvent();
  cout << "set  NTevent pointer " << endl;
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
  
  
  
  
  
 
  
  //**************************************
  //**************************************
  //******* fakes DD estimate ************
  //**************************************
  //**************************************
  
  
  
  
  scaleElec = 1.0; // 1 to switch off
  resolElec = 0.0; // 0 to switch off
  
  
  ITypeMC     = -1;
  ICut        = -1;  
  
  
  //************************************
  //For trigger systematics 
  
  /*if(datasetName=="TTbarTriggerUp"){
    SF_trig_mumu += SF_trig_mumu_error;
    SF_trig_emu+= SF_trig_emu_error;  
    SF_trig_ee += SF_trig_ee_error; 
  } 
  if(datasetName=="TTbarTriggerDown"){
    SF_trig_mumu  -= SF_trig_mumu_error;
    SF_trig_emu -= SF_trig_emu_error;  
    SF_trig_ee  -= SF_trig_ee_error; 
  }*/ 
  
  
  
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
  
   
  MyhistoManager.CreateHisto(CutFlow_mumumu,  "CutFlow_mumumu" ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_mumue,   "CutFlow_mumue"  ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_eemu,    "CutFlow_eemu"   ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_eee,     "CutFlow_eee"    ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  
  
  MyhistoManager.SetCutFlowAxisTitleFCNCMonotop(CutFlow_mumumu,   "CutFlow_mumumu"  ,datasetName);
  MyhistoManager.SetCutFlowAxisTitleFCNCMonotop(CutFlow_mumue,    "CutFlow_mumue"   ,datasetName);
  MyhistoManager.SetCutFlowAxisTitleFCNCMonotop(CutFlow_eemu,     "CutFlow_eemu"    ,datasetName);
  MyhistoManager.SetCutFlowAxisTitleFCNCMonotop(CutFlow_eee,      "CutFlow_eee"     ,datasetName);    
  
  MyhistoManager.CreateHisto(ErrCutFlow_mumumu,  "ErrCutFlow_mumumu"  ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_mumue,   "ErrCutFlow_mumue"   ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_eemu,    "ErrCutFlow_eemu"    ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_eee,     "ErrCutFlow_eee"     ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  
 
  /*MyhistoManager.CreateHisto2D(InvM_ll_vs_mWT_mumumu_afterleptsel, "InvM_ll_vs_mWT_mumumu_afterleptsel", datasetName, "Mll", 20,0., 200, "mWT", 20,0., 200.);
  MyhistoManager.CreateHisto2D(InvM_ll_vs_mWT_mumue_afterleptsel,  "InvM_ll_vs_mWT_mumue_afterleptsel" , datasetName, "Mll", 20,0., 200, "mWT", 20,0., 200.);
  MyhistoManager.CreateHisto2D(InvM_ll_vs_mWT_eemu_afterleptsel,   "InvM_ll_vs_mWT_eemu_afterleptsel"  , datasetName, "Mll", 20,0., 200, "mWT", 20,0., 200.);
  MyhistoManager.CreateHisto2D(InvM_ll_vs_mWT_eee_afterleptsel,    "InvM_ll_vs_mWT_eee_afterleptsel"   , datasetName, "Mll", 20,0., 200, "mWT", 20,0., 200.);
  */
  
  TString treename = "Ttree_"+datasetName;
  
  TheTree = new TTree(treename.Data(),treename.Data());
  TheTree->Branch("tree_topMass",     &tree_topMass,     "tree_topMass/F"   );
  TheTree->Branch("tree_deltaPhilb",  &tree_deltaPhilb,  "tree_deltaPhilb/F");
  TheTree->Branch("tree_asym",        &tree_asym,        "tree_asym/F"      );
  TheTree->Branch("tree_Zpt",         &tree_Zpt,         "tree_Zpt/F"       );
  TheTree->Branch("tree_EvtWeight",   &tree_EvtWeight,   "tree_EvtWeight/F" );
  TheTree->Branch("tree_SampleType",  &tree_SampleType,  "tree_SampleType/I");
  TheTree->Branch("tree_Channel",     &tree_Channel,     "tree_Channel/I"   );
  
  
  
  
  if (IReweight ) {
    
    string mcfile;
    if( datasetName == "FCNCkut" ) // FastSim, in-time PU only
      mcfile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PUMC_InTime_Fall11.root");
    else
      mcfile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PU3DMC_Fall11_JLA.root");
    fexists(mcfile, true);

    string datafile;
    if( datasetName == "FCNCkut" ) {
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
  vector<NTElectron> selElectrons = sel.GetSelectedElectrons();
  // cout << "Electrons loaded " << endl;
  vector<NTMuon>     selMuons     = sel.GetSelectedMuons();
  // cout << "Muons loaded " << endl;
  
  vector<NTElectron> selElectronsNonIso = sel.GetSelectedElectronsNoIso();
  // cout << "No iso electrons loaded " << endl;
  vector<NTMuon>     selMuonsNonIso     = sel.GetSelectedMuonsNoIso();
  // cout << "No iso muons loaded " << endl;
  
  vector<NTElectron> ZeeCand; 
  vector<NTMuon>     ZmumuCand; 
  
  vector<NTElectron> WeCand; 
  vector<NTMuon>     WmuCand; 
  
  NTMET met			   = sel.GetMET(); 
  // cout << "MET loaded " << endl;
  
  tree_topMass	  = -100000;
  tree_deltaPhilb = -100000;
  tree_asym	  = -100000;
  tree_Zpt	  = -100000;
  tree_SampleType = -100000;
  tree_Channel    = -100000;
  tree_EvtWeight  = -100000;
  
  
  //cout << "***************************" << endl;
  //cout << "***************************" << endl;
  //cout << "***************************" << endl;
  
  
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
  
  
  //*****************************************************************
  // calcul the MC weights
  //*****************************************************************    
  if ( datasetName!="DataEG"   && datasetName!="DataMu" && 
       datasetName!="DataMuEG" && datasetName!="DataEGMu" 
       && datasetName!="MET1"  && datasetName!="MET2") {
    
    
    if(IReweight ){
      
      if( datasetName == "FCNCkut" ) // FastSim, in-time PU only
	weightITypeMC = weightITypeMC_save*LumiWeights->ITweight(event->pileup.intime_npu);
      else
	weightITypeMC = weightITypeMC_save*LumiWeights->weight3D(event->pileup.before_npu, event->pileup.intime_npu, event->pileup.after_npu);
      
    }
    else weightITypeMC = weightITypeMC_save;
  }
  else weightITypeMC = 1;
  
  
  //*****************************************************************
  // determine top decay channel
  //*****************************************************************    
  
  bool IsTTbarDilept = false;
  bool IsSignal = false;
  double WeightForBranchingRatio = 1.;
  bool IsLJ = false;
  
  
  
  
  
  
  
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
    ITypeMC = 100;  Dweight[ITypeMC]= weightITypeMC; 
    EventYieldWeightError = Dweight[ITypeMC]*Dweight[ITypeMC];   
    
    
  }
  
  
  
  //*****************************************************************
  // pass trigger selection
  //*****************************************************************   
  
  bool passtrigger = false;
  
  passtrigger = sel.passTriggerSelection ( dataset, "mumu");
  
  
  if (   passtrigger   ) {
    
    
  } // end selection trigger

  return kTRUE;
  
} //end loops on datasets


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
    
    
      
  MyhistoManager.WriteMyHisto(CutFlow_mumumu, "all" );
  MyhistoManager.WriteMyHisto(CutFlow_mumue,  "all" );
  MyhistoManager.WriteMyHisto(CutFlow_eemu,   "all" );
  MyhistoManager.WriteMyHisto(CutFlow_eee,    "all" );
  
  MyhistoManager.WriteMyHisto(ErrCutFlow_mumumu,  "all");
  MyhistoManager.WriteMyHisto(ErrCutFlow_mumue,   "all");
  MyhistoManager.WriteMyHisto(ErrCutFlow_eemu,    "all");
  MyhistoManager.WriteMyHisto(ErrCutFlow_eee,     "all");
  
  
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
  
  
    //delete file1  ;
    //delete hPUMC ;  
    //delete file2 ; 
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
