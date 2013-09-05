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
  
  
  
  applyLeptonSF  = true;
  applyTrigger   = true;
  
  
  
  IReweight		= true;
  
  doBTagCVScorr = true;
  //pdf.Initialize();
  doBTagCSV_syst = 0;
  
  
  
  
  applyJES  = false;
  scale     = -1; // +1 or -1
  applyJER  = false;
  ResFactor = 0.1;
  
  applyLeptonSFUp    = false;
  applyLeptonSFDown  = false;
  
  applyTriggerUp   = false;
  applyTriggerDown = false;
  
  //to do rename
  IReweight_puUp	= false;
  IReweight_puDown	= false;
  
  
  
  
  doPDF = false;
  //pdftype =0 ;
  pdftype =1  ;
  
  PDFmode=-1; // PDF4LHC recipe
  //PDFmode=0 ;// enveloppe avec les 4 NNPDF2.1
  //PDFmode=1 ;// NNPDF21
  //PDFmode=2 ;// NNPDF21_mc15
  //PDFmode=3 ;// NNPDF21_mc16
  //PDFmode=4 ;// NNPDF21_mc17
  
  
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
  tree->ls();
  cout << "start init tree " << endl;
  // Set branch addresses
  branch = (TBranch *) tree->GetBranch("NTEvent");
  cout << "branch adress retrieved " << endl;
  event = new IPHCTree::NTEvent();
  cout << "create event pointer " << endl;
  if(tree==0)   cout << "link to the tree is null" << endl;
  if(branch==0) cout << "branch pointer is a null pointer 1" << endl;
  if((TBranch *) tree->GetBranch("NTEvent") == 0) cout << "branch pointer is a null pointer 2" << endl;
  if((TBranch *) tree->GetBranch("NTEvent") == 0) cout << "branch pointer is a null pointer 2" << endl;
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
  
  
 
  
  //**************************************
  //**************************************
  //******* fakes DD estimate ************
  //**************************************
  //**************************************
  
  
  
  ITypeMC     = -1;
  ICut        = -1;  
  
  
  //************************************
  //For trigger systematics 
  
  if(applyTriggerUp){
    
    
    
  } 
  if(applyTriggerDown){
    
    
    
  } 
  
  
  
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
  MyhistoManager.LoadChannels(TheChannelName);
  //example
  
  nbins = 200;
  minx = 0.;
  maxx = 350;
  
  
  //***********************
  // initiate histograms
  
  MyhistoManager.CreateHisto(CutFlow_mumu,  "CutFlow_mumu" ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_emu,   "CutFlow_emu"  ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(CutFlow_ee,    "CutFlow_ee"   ,datasetName,"CutFlow","Entries",15,-0.5,14.5);
  
  MyhistoManager.CreateHisto(ErrCutFlow_mumu,  "ErrCutFlow_mumu" ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_emu,   "ErrCutFlow_emu"  ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  MyhistoManager.CreateHisto(ErrCutFlow_ee,    "ErrCutFlow_ee"   ,datasetName,"ErrCutFlow","Entries",15,-0.5,14.5);
  
  
  
  
  
  
  if (IReweight ) {
    
    //to do : update PU reweighting input file
    string mcfile = getenv( "CMSSW_BASE" )+string("/src/IPHCAnalysis/NTuple/macros/data/PU3DMC_Fall11_JLA.root");
    fexists(mcfile, true);
    
    string datafile;
    if( !IReweight_puDown && !IReweight_puUp ) datafile = getenv( "CMSSW_BASE" )+string("/src/IPHCAnalysis/NTuple/macros/data/PUData2011_68mb.root");
    if( IReweight_puDown ) datafile = getenv( "CMSSW_BASE" )+string("/src/IPHCAnalysis/NTuple/macros/data/PUData2011_64.6mb.root");
    if( IReweight_puUp ) datafile = getenv( "CMSSW_BASE" )+string("/src/IPHCAnalysis/NTuple/macros/data/PUData2011_71.4mb.root");
    
    fexists(datafile, true);
    
    
    //to do : update PU for the last recipe
    LumiWeights = new reweight::LumiReWeighting(mcfile, datafile, "histoMCPU", "pileup" );
    LumiWeights->weight3D_init( 1. );
    
  }
  
  
  //if(doPDF) pdf.Initialize();
  
  if (doPDF)
    { 
      if (PDFmode==-1) pdf.Initialize();
      else pdf2.Initialize();
    }
  
  
  
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
  
  
  cout << "line 347 " << endl;
  
  
  //---------------------------------------------------//
  // Main event loop: get entry
  //---------------------------------------------------//
  fChain->GetTree()->GetEntry(entry); 
  branch->GetEntry(entry);
  
  IPHCTree::NTTransient::InitializeAfterReading(event);
  
  
  //to do : add "isData" from the xml config file
  bool isData = false;    
  if(datasetName=="DataDiEG"   || datasetName=="DataDiMu" || 
     datasetName=="DataDiMuEG" || datasetName=="DataDiEGMu" ) isData = true;
  
  
  //---------------------------------------------------//
  //         Doing the analysis event by event
  //---------------------------------------------------//
  int debugcc=1000;
  int maxdebugcc=10;
  sel.LoadEvent(event);
  
  //Collection of selected objects
  vector<NTVertex>   selVertices  = sel.GetSelectedVertex();
  vector<NTElectron> selElectrons = sel.GetSelectedElectrons(); 
  vector<NTMuon>     selMuons     = sel.GetSelectedMuons();
  //NTMET met			   = sel.GetMET(); 
  
  NTMET met			   = sel.GetSelectedMET(applyJES, scale, applyJER, ResFactor);  
  // cout << "MET loaded " << endl;
  
  
  
  
  double Dweight[101];
  for(int k1=0; k1<101; k1++) {
    Dweight[k1] = 0.;
  }   
  
  
  double weightITypeMC_save = Luminosity*dataset->Xsection()/dataset->getNSkimmedEvent();
  double weightITypeMC=0;
  
  //*****************************************************************
  // Loop over the channels (lepton pairs which triggered the events)
  //*****************************************************************
  
  double EventYieldWeightError = 0;
  
  
  //to do : add DY overlap
  
  //for (int IChannel=0; IChannel<3; IChannel++) {
    
    int IChannel = -1;
    
    if( datasetName=="DataDiMu"   ) IChannel = 0;
    if( datasetName=="DataDiMuEG" ) IChannel = 1;
    if( datasetName=="DataDiEG"   ) IChannel = 2;
    
    
    TString ChannelName = "";
    if      (IChannel==0) ChannelName= "mumu"; 
    else if (IChannel==1) ChannelName= "emu" ; 
    else if (IChannel==2) ChannelName= "ee"  ;  
    
    //if (IChannel==0 && datasetName!="DataDiMu"   ) continue;
    //if (IChannel==1 && datasetName!="DataDiMuEG" ) continue;
    //if (IChannel==2 && datasetName!="DataDiEG"   ) continue;
    
    
    
    //*****************************************************************
    // calcul the MC weights
    //*****************************************************************   
    //to do : use isData parameter 
    if ( datasetName!="DataDiEG"   && datasetName!="DataDiMu" && 
	 datasetName!="DataDiMuEG" && datasetName!="DataDiEGMu" ) {
      
      
      if(IReweight ){
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
    //bool IsLJ = false;
    
    
    
    
    //to do : fix "is signal" in filling histograms
    //*****************************************************************
    // determine top decay channel
    //*****************************************************************    
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
      
    }
    
    
    
    
    //to do : fix IsSignal
    
    //*****************************************************************
    // determine MC evetn weight
    //*****************************************************************    
    
    std::vector< double > thereweight = determineWeights(datasetName, weightITypeMC, WeightForBranchingRatio);
    ITypeMC = thereweight[0];
    Dweight[ITypeMC] =   thereweight[1];
    EventYieldWeightError = thereweight[2];
    if(thereweight[3] > 0) IsSignal = true; else IsSignal = false;
    
    
    
    //*****************************************************************
    // fill cutflow before any selection
    //*****************************************************************   
    
    
     
    if(IChannel == 0){
      MyhistoManager.FillHisto(CutFlow_mumu,      "CutFlow_mumu",    0, datasetName, IsSignal, Dweight[ITypeMC]);    
      MyhistoManager.FillHisto(ErrCutFlow_mumu,   "ErrCutFlow_mumu", 0, datasetName, IsSignal, EventYieldWeightError);
    }
    else if(IChannel == 1){
      MyhistoManager.FillHisto(CutFlow_emu,       "CutFlow_emu",    0, datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(ErrCutFlow_emu,    "ErrCutFlow_emu", 0, datasetName, IsSignal, EventYieldWeightError);
    }
    else if(IChannel == 2) {
      MyhistoManager.FillHisto(CutFlow_ee,        "CutFlow_ee",    0, datasetName, IsSignal, Dweight[ITypeMC]);
      MyhistoManager.FillHisto(ErrCutFlow_ee,     "ErrCutFlow_ee", 0, datasetName, IsSignal, EventYieldWeightError);
    }
    
  cout << "line 497 " << endl;
    
    //*****************************************************************
    // pass trigger selection
    //*****************************************************************   
    
    bool passtrigger = false;
    
    
    //to do : update trigger selection
    if(ChannelName == "mumu" ) passtrigger = sel.passTriggerSelection ( dataset, "mumu");
    if(ChannelName == "emu"  ) passtrigger = sel.passTriggerSelection ( dataset, "emu" );
    if(ChannelName == "ee"   ) passtrigger = sel.passTriggerSelection ( dataset, "ee" );
    
    
    if (   passtrigger   ) {
      
      //*****************************************************************
      // fill cutflow after trigger selection
      //*****************************************************************   
      
      if(     IChannel == 0){
        MyhistoManager.FillHisto(CutFlow_mumu,      "CutFlow_mumu",    1, datasetName, IsSignal, Dweight[ITypeMC]);    
        MyhistoManager.FillHisto(ErrCutFlow_mumu,   "ErrCutFlow_mumu", 1, datasetName, IsSignal, EventYieldWeightError);
      }
      else if(IChannel == 1){
        MyhistoManager.FillHisto(CutFlow_emu,       "CutFlow_emu",    1, datasetName, IsSignal, Dweight[ITypeMC]);
        MyhistoManager.FillHisto(ErrCutFlow_emu,    "ErrCutFlow_emu", 1, datasetName, IsSignal, EventYieldWeightError);
      }
      else if(IChannel == 2 ) {
        MyhistoManager.FillHisto(CutFlow_ee,        "CutFlow_ee",    1, datasetName, IsSignal, Dweight[ITypeMC]);
        MyhistoManager.FillHisto(ErrCutFlow_ee,     "ErrCutFlow_ee", 1, datasetName, IsSignal, EventYieldWeightError);
      }
      
      
      
      //*****************************************************************
      // determine decay channel
      //***************************************************************** 
      string decayChannel_tmp = "";
      sel.GetLeptonPair(selMuons, selElectrons, decayChannel_tmp ); 
      TString decayChannel = decayChannel_tmp;
      //to do : makes pat iso configurable
      
      
      
      std::vector<TH1F>* pCutFlow;
      std::vector<TH1F>* pErrCutFlow;
      
      if(IChannel == 0 && decayChannel == "mumu") pCutFlow = &CutFlow_mumu;
      if(IChannel == 1 && decayChannel == "emu" ) pCutFlow = &CutFlow_emu;
      if(IChannel == 2 && decayChannel == "ee"  ) pCutFlow = &CutFlow_ee;
      
      
      if(IChannel == 0 && decayChannel == "mumu") pErrCutFlow = &ErrCutFlow_mumu;
      if(IChannel == 1 && decayChannel == "emu" ) pErrCutFlow = &ErrCutFlow_emu;
      if(IChannel == 2 && decayChannel == "ee"  ) pErrCutFlow = &ErrCutFlow_ee;
      
       
      
      //*****************************************************************
      // apply dilepton selection
      //*****************************************************************   
   
      cout << "line 545 " << endl;
      if(decayChannel != ""){
      
	//*****************************************************************
	// fill cutflow after lepton selection
	//*****************************************************************   
	
	 MyhistoManager.FillHisto(*pCutFlow,      ("CutFlow_"+decayChannel).Data(),    2, datasetName, IsSignal, Dweight[ITypeMC]);    
	 MyhistoManager.FillHisto(*pErrCutFlow,   ("ErrCutFlow_"+decayChannel).Data(), 2, datasetName, IsSignal, EventYieldWeightError);
	
	
	
	/*
	if(     IChannel == 0 && decayChannel == "mumu"){
	  MyhistoManager.FillHisto(CutFlow_mumu,      "CutFlow_mumu",    2, datasetName, IsSignal, Dweight[ITypeMC]);    
	  MyhistoManager.FillHisto(ErrCutFlow_mumu,   "ErrCutFlow_mumu", 2, datasetName, IsSignal, EventYieldWeightError);
	}
	else if(IChannel == 1 && decayChannel == "emu" ){
	  MyhistoManager.FillHisto(CutFlow_emu,       "CutFlow_emu",    2, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(ErrCutFlow_emu,    "ErrCutFlow_emu", 2, datasetName, IsSignal, EventYieldWeightError);
	}
	else if(IChannel == 2 && decayChannel == "ee"  ) {
	  MyhistoManager.FillHisto(CutFlow_ee,        "CutFlow_ee",    2, datasetName, IsSignal, Dweight[ITypeMC]);
	  MyhistoManager.FillHisto(ErrCutFlow_ee,     "ErrCutFlow_ee", 2, datasetName, IsSignal, EventYieldWeightError);
	}
	*/
	//*****************************************************************
	// apply trigger scale factors
	//*****************************************************************  
	//to do update to pT eta dependent SF ?
	if(applyTrigger  &&  !isData ){	
	  
	  //if(IChannel == 0 &&  decayChannel == "mumu") Dweight[ITypeMC]*=SF_trig_mumu;
	  //if(IChannel == 1 &&  decayChannel == "emu" ) Dweight[ITypeMC]*=SF_trig_emu;
	  //if(IChannel == 2 &&  decayChannel == "ee"  ) Dweight[ITypeMC]*=SF_trig_ee;
	  
	}
	
	//*****************************************************************
	// apply lepton scale factors
	//*****************************************************************    
	
	
	double LeptonSF      = 0.;
	
	if(applyLeptonSF && !isData){
	  //to do : update to last SF values 
	  if(IChannel == 0 &&  decayChannel== "mumu"){
	    LeptonSF = sel.getLeptonScaleFactor( selMuons[0].p4.Pt(), selMuons[0].p4.Eta(), "mu") 
	      * sel.getLeptonScaleFactor( selMuons[1].p4.Pt(), selMuons[1].p4.Eta(), "mu");
	    if(applyLeptonSFUp)
	      LeptonSF = (sel.getLeptonScaleFactor( selMuons[0].p4.Pt(), selMuons[0].p4.Eta(), "mu")+sel.getLeptonScaleFactorError( selMuons[0].p4.Pt(), selMuons[0].p4.Eta(), "mu"))
		* (sel.getLeptonScaleFactor( selMuons[1].p4.Pt(), selMuons[1].p4.Eta(), "mu")+sel.getLeptonScaleFactorError( selMuons[1].p4.Pt(), selMuons[1].p4.Eta(), "mu"));		  
	    if(applyLeptonSFDown)
	      LeptonSF = (sel.getLeptonScaleFactor( selMuons[0].p4.Pt(), selMuons[0].p4.Eta(), "mu")-sel.getLeptonScaleFactorError( selMuons[0].p4.Pt(), selMuons[0].p4.Eta(), "mu"))
		* (sel.getLeptonScaleFactor( selMuons[1].p4.Pt(), selMuons[1].p4.Eta(), "mu")-sel.getLeptonScaleFactorError( selMuons[1].p4.Pt(), selMuons[1].p4.Eta(), "mu"));		  
	  }
	  if(IChannel == 1 &&  decayChannel== "emu"){
	    LeptonSF = sel.getLeptonScaleFactor( selMuons[0].p4.Pt()    , selMuons[0].p4.Eta()    , "mu") 
	      * sel.getLeptonScaleFactor( selElectrons[0].p4.Pt(), selElectrons[0].p4.Eta(), "e");
	    if(applyLeptonSFUp)
	      LeptonSF = (sel.getLeptonScaleFactor( selMuons[0].p4.Pt()    , selMuons[0].p4.Eta(),     "mu")+sel.getLeptonScaleFactorError( selMuons[0].p4.Pt(), selMuons[0].p4.Eta(), "mu"))
		* (sel.getLeptonScaleFactor( selElectrons[0].p4.Pt(), selElectrons[0].p4.Eta(), "mu")+sel.getLeptonScaleFactorError( selElectrons[0].p4.Pt(), selElectrons[0].p4.Eta(), "e"));		  
	    if(applyLeptonSFDown)
	      LeptonSF = (sel.getLeptonScaleFactor( selMuons[0].p4.Pt()    , selMuons[0].p4.Eta(),     "mu")-sel.getLeptonScaleFactorError( selMuons[0].p4.Pt(), selMuons[0].p4.Eta(), "mu"))
		* (sel.getLeptonScaleFactor( selElectrons[0].p4.Pt(), selElectrons[0].p4.Eta(), "mu")-sel.getLeptonScaleFactorError( selElectrons[0].p4.Pt(), selElectrons[0].p4.Eta(), "e"));		  
	  }
	  if(IChannel == 2 &&  decayChannel== "ee"){
	    LeptonSF = sel.getLeptonScaleFactor( selElectrons[0].p4.Pt(), selElectrons[0].p4.Eta(), "e") 
	      * sel.getLeptonScaleFactor( selElectrons[1].p4.Pt(), selElectrons[1].p4.Eta(), "e");
	    if(applyLeptonSFUp)
	      LeptonSF = (sel.getLeptonScaleFactor( selElectrons[0].p4.Pt(), selElectrons[0].p4.Eta(), "e")+sel.getLeptonScaleFactorError( selElectrons[0].p4.Pt(), selElectrons[0].p4.Eta(), "e"))
		* (sel.getLeptonScaleFactor( selElectrons[1].p4.Pt(), selElectrons[1].p4.Eta(), "e")+sel.getLeptonScaleFactorError( selElectrons[1].p4.Pt(), selElectrons[1].p4.Eta(), "e"));		  
	    if(applyLeptonSFDown)
	      LeptonSF = (sel.getLeptonScaleFactor( selElectrons[0].p4.Pt(), selElectrons[0].p4.Eta(), "e")-sel.getLeptonScaleFactorError( selElectrons[0].p4.Pt(), selElectrons[0].p4.Eta(), "e"))
		* (sel.getLeptonScaleFactor( selElectrons[1].p4.Pt(), selElectrons[1].p4.Eta(), "e")-sel.getLeptonScaleFactorError( selElectrons[1].p4.Pt(), selElectrons[1].p4.Eta(), "e"));		  
	  }
	  //to do for synch, leptonSF = 1
	  LeptonSF=1;
	  Dweight[ITypeMC]*=LeptonSF;
	}
	
	
	
	
  cout << "line 624 " << endl;
	
	
	  
        float InvDilMass = 0;
        if (ChannelName=="mumu")  InvDilMass = (selMuons[0].p4     +selMuons[1].p4    ).M();
        if (ChannelName=="emu")   InvDilMass = (selMuons[0].p4     +selElectrons[0].p4).M();
        if (ChannelName=="ee")    InvDilMass = (selElectrons[0].p4 +selElectrons[1].p4).M();
	
	
	
	//*****************************************************************
	// apply lepton invM selection
	//*****************************************************************   
	
	if( fabs(InvDilMass-91) < 15) {
	  
	  
	  //*****************************************************************
	  // fill cutflow lepton invM selection
	  //*****************************************************************   
	
	  if(     IChannel == 0 && decayChannel == "mumu"){
	    MyhistoManager.FillHisto(CutFlow_mumu,      "CutFlow_mumu",    3, datasetName, IsSignal, Dweight[ITypeMC]);    
	    MyhistoManager.FillHisto(ErrCutFlow_mumu,   "ErrCutFlow_mumu", 3, datasetName, IsSignal, EventYieldWeightError);
	  }
	  else if(IChannel == 1 && decayChannel == "emu" ){
	    MyhistoManager.FillHisto(CutFlow_emu,       "CutFlow_emu",    3, datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(ErrCutFlow_emu,    "ErrCutFlow_emu", 3, datasetName, IsSignal, EventYieldWeightError);
	  }
	  else if(IChannel == 2 && decayChannel == "ee"  ) {
	    MyhistoManager.FillHisto(CutFlow_ee,        "CutFlow_ee",    3, datasetName, IsSignal, Dweight[ITypeMC]);
	    MyhistoManager.FillHisto(ErrCutFlow_ee,     "ErrCutFlow_ee", 3, datasetName, IsSignal, EventYieldWeightError);
	  }
	  
	  
  cout << "line 660 " << endl;
	  
	  //*****************************************************************
	  // apply jets mutli. selection
	  //*****************************************************************   
		
	  vector<NTJet>  selJets = sel.GetSelectedJets(selMuons, selElectrons, applyJES, scale, applyJER, ResFactor);
	  if(selJets.size() >= 2){
	    
	    
	    //*****************************************************************
	    // fill cutflow jets mutli. selection
	    //*****************************************************************   
	
	    if(     IChannel == 0 && decayChannel == "mumu"){
	      MyhistoManager.FillHisto(CutFlow_mumu,      "CutFlow_mumu",    3, datasetName, IsSignal, Dweight[ITypeMC]);    
	      MyhistoManager.FillHisto(ErrCutFlow_mumu,   "ErrCutFlow_mumu", 3, datasetName, IsSignal, EventYieldWeightError);
	    }
	    else if(IChannel == 1 && decayChannel == "emu" ){
	      MyhistoManager.FillHisto(CutFlow_emu,       "CutFlow_emu",    3, datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(ErrCutFlow_emu,    "ErrCutFlow_emu", 3, datasetName, IsSignal, EventYieldWeightError);
	    }
	    else if(IChannel == 2 && decayChannel == "ee"  ) {
	      MyhistoManager.FillHisto(CutFlow_ee,        "CutFlow_ee",    3, datasetName, IsSignal, Dweight[ITypeMC]);
	      MyhistoManager.FillHisto(ErrCutFlow_ee,     "ErrCutFlow_ee", 3, datasetName, IsSignal, EventYieldWeightError);
	    }
	  
	    
	    //*****************************************************************
	    // apply met selection
	    //***************************************************************** 
	    double theMET = met.p2.Mod();
	    
  cout << "line 693 " << endl;
	    if(decayChannel == "emu" || theMET > 40){
	      
	       
	      //*****************************************************************
	      // fill cutflow jets mutli. selection
	      //*****************************************************************   
	
	      if(     IChannel == 0 && decayChannel == "mumu"){
	        MyhistoManager.FillHisto(CutFlow_mumu,      "CutFlow_mumu",    4, datasetName, IsSignal, Dweight[ITypeMC]);    
	        MyhistoManager.FillHisto(ErrCutFlow_mumu,   "ErrCutFlow_mumu", 4, datasetName, IsSignal, EventYieldWeightError);
	      }
	      else if(IChannel == 1 && decayChannel == "emu" ){
	        MyhistoManager.FillHisto(CutFlow_emu,       "CutFlow_emu",    4, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(ErrCutFlow_emu,    "ErrCutFlow_emu", 4, datasetName, IsSignal, EventYieldWeightError);
	      }
	      else if(IChannel == 2 && decayChannel == "ee"  ) {
	        MyhistoManager.FillHisto(CutFlow_ee,        "CutFlow_ee",    4, datasetName, IsSignal, Dweight[ITypeMC]);
	        MyhistoManager.FillHisto(ErrCutFlow_ee,     "ErrCutFlow_ee", 4, datasetName, IsSignal, EventYieldWeightError);
	      }
	     
	     
	     
	     
  cout << "line 717 " << endl;
	      //*****************************************************************
	      // calculate btagging SF
	      //***************************************************************** 
	      
	      int NBtaggedJets    = 0;
	      int idxBtag         = 0;
	      int AlgoBtag        = sel.GetbtagAlgo();
	      float btagDiscriCut = sel.GetbtagDiscriCut ();
	
	      bool foundASelBjet = 0;
	      for(unsigned int ijet = 0; ijet < selJets.size(); ijet++){
	        if(abs(selJets[ijet].partonFlavour)==5 ){
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
	      
	      
	      // initialisation of weightb
	      vector < float >weightb;
	      weightb.push_back (1.);
	      weightb.push_back (1.);
	      weightb.push_back (1.);
	      weightb.push_back (1.);
	      weightb.push_back (1.);
	
	      if (sel.GetFlagb() == 1) {	// check if the weightb computation is needed or not (depending on "flag" in xml)
	        if (!isData) {	 // to be applied only for MC 
	          vector < float >weight_temp = sel.GetSFBweight().GetWeigth4BSel (sel.GetMethodb(), sel.GetSystb(),selJets);
	          weightb[0] = weight_temp[0];  //weight of the event (depending on "NofBtagJets" in xml)
	          weightb[1] = weight_temp[1];  //proba 0 jet
	          weightb[2] = weight_temp[2];  //proba 1 jet
	          weightb[3] = weight_temp[3];  //proba 2 jets
	          weightb[4] = weight_temp[4];  //proba at least 3 jets		
	        }
	      }
	      
	      
	      double btagweight_1selBtag = 1;
	      double btagweight_2selBtag = 1;
	      if(!isData) {
	        btagweight_1selBtag = 1-weightb[1];
	        btagweight_2selBtag = weightb[0];
	      }
	      
	      
	      
	      
	      //*****************************************************************
	      // apply btag selection
	      //***************************************************************** 
	      
  cout << "line 779 " << endl;

	      if( (isData && NBtaggedJets >=1) || !isData){
	        
	        //*****************************************************************
	        // fill cutflow jets mutli. selection
	        //*****************************************************************   
	        
  cout << "line 785 " << endl;
		
		Dweight[ITypeMC] *= btagweight_1selBtag;
		
	        if(     IChannel == 0 && decayChannel == "mumu"){
	          MyhistoManager.FillHisto(CutFlow_mumu,      "CutFlow_mumu",    5, datasetName, IsSignal, Dweight[ITypeMC]);    
	          MyhistoManager.FillHisto(ErrCutFlow_mumu,   "ErrCutFlow_mumu", 5, datasetName, IsSignal, EventYieldWeightError);
	        }
	        else if(IChannel == 1 && decayChannel == "emu" ){
	          MyhistoManager.FillHisto(CutFlow_emu,       "CutFlow_emu",    5, datasetName, IsSignal, Dweight[ITypeMC]);
	          MyhistoManager.FillHisto(ErrCutFlow_emu,    "ErrCutFlow_emu", 5, datasetName, IsSignal, EventYieldWeightError);
	        }
	        else if(IChannel == 2 && decayChannel == "ee"  ) {
	          MyhistoManager.FillHisto(CutFlow_ee,        "CutFlow_ee",    5, datasetName, IsSignal, Dweight[ITypeMC]);
	          MyhistoManager.FillHisto(ErrCutFlow_ee,     "ErrCutFlow_ee", 5, datasetName, IsSignal, EventYieldWeightError);
	        }
		
	        if( (isData && NBtaggedJets >=2) || !isData){
	        
	          //*****************************************************************
	          // fill cutflow jets mutli. selection
	          //*****************************************************************   
	  
		  Dweight[ITypeMC] = btagweight_2selBtag*Dweight[ITypeMC]/btagweight_1selBtag;
		  
  cout << "line 812 " << endl;
		  
	          if(     IChannel == 0 && decayChannel == "mumu"){
	            MyhistoManager.FillHisto(CutFlow_mumu,      "CutFlow_mumu",    6, datasetName, IsSignal, Dweight[ITypeMC]);    
	            MyhistoManager.FillHisto(ErrCutFlow_mumu,   "ErrCutFlow_mumu", 6, datasetName, IsSignal, EventYieldWeightError);
	          }
	          else if(IChannel == 1 && decayChannel == "emu" ){
	            MyhistoManager.FillHisto(CutFlow_emu,       "CutFlow_emu",    6, datasetName, IsSignal, Dweight[ITypeMC]);
	            MyhistoManager.FillHisto(ErrCutFlow_emu,    "ErrCutFlow_emu", 6, datasetName, IsSignal, EventYieldWeightError);
	          }
	          else if(IChannel == 2 && decayChannel == "ee"  ) {
	            MyhistoManager.FillHisto(CutFlow_ee,        "CutFlow_ee",    6, datasetName, IsSignal, Dweight[ITypeMC]);
	            MyhistoManager.FillHisto(ErrCutFlow_ee,     "ErrCutFlow_ee", 6, datasetName, IsSignal, EventYieldWeightError);
	          }
		  
		  
	        } // end selection 2bjet sel
	      } // end selection 1bjet sel
	    } // end selection met sel
	  } // end selection NJet sel
	} // end selection Minv
      } // end selection lepton pair
    } // end selection trigger
  //} //end loops on datasets
  
  
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
    
  MyhistoManager.WriteMyHisto(CutFlow_mumu, "all" );
  MyhistoManager.WriteMyHisto(CutFlow_emu,  "all" );
  MyhistoManager.WriteMyHisto(CutFlow_ee,   "all" );
  
    
  MyhistoManager.WriteMyHisto(ErrCutFlow_mumu, "all" );
  MyhistoManager.WriteMyHisto(ErrCutFlow_emu,  "all" );
  MyhistoManager.WriteMyHisto(ErrCutFlow_ee,   "all" );
  
  
   //The following line is mandatory to copy everythin in a common RootFile
    fOutput->Add(fProofFile);
    
  
    
    
   // delete TheTree;
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






std::vector<double> ProofSelectorMyCutFlow::determineWeights(TString datasetName, double weightITypeMC, double WeightForBranchingRatio){

    
    double ITypeMC = 0;
    double Dweight = 0 ;
    double EventYieldWeightError = 0 ;
    double IsSignal = 1 ;
    
    Dweight= weightITypeMC * WeightForBranchingRatio;
    EventYieldWeightError = Dweight*Dweight;
    
      
    if ( datasetName=="TTbar" ) { 
      ITypeMC = 1; 
      Dweight= weightITypeMC * WeightForBranchingRatio;
      EventYieldWeightError = Dweight*Dweight;
    }
    else {   
      IsSignal = -1; 
      Dweight= weightITypeMC; 
      EventYieldWeightError = Dweight*Dweight;
      if ( datasetName=="Zjets" || datasetName=="DYToLL_M10-50"
	      ) ITypeMC = 2; 
      else if (  datasetName=="Wjets"  ) ITypeMC = 3;
      else if ( datasetName=="SingleToptW" || datasetName=="TtW" || datasetName=="TbartW"
	      || datasetName=="TtWScaleUp" || datasetName=="TtWScaleDown"
	      || datasetName=="TbartWScaleUp" || datasetName=="TbartWScaleDown") ITypeMC = 4;  
      else if ( datasetName=="WZ" || datasetName=="WW" || datasetName=="ZZ"  )   ITypeMC = 5; 
    } 
    
    
    
    if ( datasetName=="DataDiEG" || datasetName=="DataDiMu" || 
	 datasetName=="DataDiMuEG" || datasetName=="DataDiEGMu") { 
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





