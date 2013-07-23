
//////////////////////////////////////////////////////////
//
// Example of TSelector implementation to do a Monte Carlo
// generation using Pythia8.
// See tutorials/proof/runProof.C, option "pythia8", for an
// example of how to run this selector.
//
//////////////////////////////////////////////////////////

#ifndef ProofSelectorMyCutFlow_h
#define ProofSelectorMyCutFlow_h

#include <TSelector.h>
#include <TTree.h>
#include <TFile.h>
#include <TProofOutputFile.h>
#include <TRandom.h>

#include "NTFormat/interface/NTEvent.h"
#include "Plots/interface/HistoManager.h"



#include "Tools/interface/Dataset.h"
#include "Tools/interface/AnalysisEnvironmentLoader.h"
#include "Selection/interface/DiLeptonSelection.h"
#include "Plots/interface/DiLepAnaHistoManager.h"
#include "Tools/interface/PUWeighting.h"
#include "Tools/interface/LumiReweightingStandAlone.h"
#include "Tools/interface/JetCorrector.h"


#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h> 
#include <TH3.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TLorentzVector.h>
//#include <iostream>

class TH1F;
class TBranch;
class TTree;

class AnalysisEnvironmentLoader;
class DiLeptonSelection;

class ProofSelectorMyCutFlow : public TSelector {
 public :
  
  // Specific members
  //Access to the tree and outputs
  TTree* fChain;
  TBranch* branch;
  IPHCTree::NTEvent* event;
  TFile            *fFile;
  TProofOutputFile *fProofFile; // For optimized merging of the ntuple
  //Pointer on results from xml file  reading
  AnalysisEnvironmentLoader* anaEL; 
  //Minimimal info
  vector<Dataset> datasets;
  Dataset* dataset;
  vector<string> CutName;
  vector<string> ChannelName;
  vector<string> VecChannelName;
  DiLeptonSelection sel; 
  float Luminosity;
  int verbosity;
  // 0: MC - 1: Data - 2 Data & MC
  int DataType;
  //Info analysis macro specific 
  

  JetCorrector JEC_L2L3Residuals;
  
  
  reweight::LumiReWeighting *LumiWeights;
  float LumiError ;
  string PUWeightFileName;
  reweight::PoissonMeanShifter PShiftUp_;
  reweight::PoissonMeanShifter PShiftDown_;

  
  // Here define the studied channel (ee/mumu/emu)
  //  string ChannelName  = "ee";  // "mumu", "ee", "emu"
  // on va tourner sur les 3 canaux en meme temps!!!
  
  bool IReweight             ;
  bool IReweight_puUp        ;
  bool IReweight_puDown      ;
  bool IDYestimateWithMetCut ;
  
  bool useNonIsoWcand;
  bool applyTrigger ;
  
  bool applyWZ ;
  bool applyFakescale ;
  bool applyLeptonSF;
  
  //Here define Scale Factors
  //SF_trigger applied for mumu
  double SF_trig_mumu ;  
  double SF_trig_emu;  
  double SF_trig_ee ;  
  
  double SF_trig_mumu_error ;  
  double SF_trig_emu_error;  
  double SF_trig_ee_error ;  
  
  double SF_BranchingRatio_ll ; 
  double SF_BranchingRatio_lj ; 
  double SF_BranchingRatio_had ; 
    
   
   
  
  bool PUup;
  bool PUdown;
  
  TRandom rand;
  
  
  
 
  
  
  double scaleElec; // 1 to switch off
  double resolElec; // 0 to switch off
  
  
  
  
  int ITypeMC ;
  int ICut    ;  
  
    
  
  TH1F* fHist;
  
  std::vector<TH1F> CutFlow_mumumu;
  std::vector<TH1F> CutFlow_mumue;
  std::vector<TH1F> CutFlow_eemu;
  std::vector<TH1F> CutFlow_eee;
  std::vector<TH1F> ErrCutFlow_mumumu;
  std::vector<TH1F> ErrCutFlow_mumue;
  std::vector<TH1F> ErrCutFlow_eemu;
  std::vector<TH1F> ErrCutFlow_eee;
  


  
  
  
  //TTree and banches used for BDT
  TTree *TheTree;
  
  
  float tree_topMass;
  float tree_deltaPhilb;
  float tree_asym;
  float tree_Zpt;
  int   tree_SampleType;
  int   tree_Channel;
  float tree_EvtWeight;
  
  ofstream ofile;
  
  //------------------------------------
  //Definition of the various histograms
  //------------------------------------
  int nbins ;
  float minx;
  float maxx;
  HistoManager MyhistoManager;
  
  
  std::map<TString, int> initMCevents;
  std::map<TString, int> skimMCevents;
 
  //------------------------------------
  // for PileUP reweighting
  //------------------------------------
  PUWeighting  thePUReweighter;
  
  
  string datasetName ;
  
  
  //------------------------------------
  //definition of member functions
  //------------------------------------
  ProofSelectorMyCutFlow();
  virtual ~ProofSelectorMyCutFlow();
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual Bool_t  Process(Long64_t entry);
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();
  
  ClassDef(ProofSelectorMyCutFlow,0);
};

#endif
