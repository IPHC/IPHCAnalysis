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


#include <TSelector.h>
#include <TProofOutputFile.h>
#include <TRandom.h>

class TH1F;
class TBranch;
class TTree;

class AnalysisEnvironmentLoader;
class TTbarMetSelection;

double pi=acos(-1.);


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
  int idataset;
  TTbarMetSelection sel; 
  float Luminosity_e;
  float Luminosity_mu;
  int verbosity;
  // 0: MC - 1: Data - 2 Data & MC
  int DataType;
  //Info analysis macro specific 
  
   JetCorrector JEC_L2L3Residuals;


  reweight::LumiReWeighting *LumiWeights;

  
  //------------------------------------
  //Definition of the various histograms
  //------------------------------------
  bool doHistoManager;
  TTbarMetHistoManager* histoManager;
  
  
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
