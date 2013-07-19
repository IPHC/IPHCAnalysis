
//////////////////////////////////////////////////////////
//
// Example of TSelector implementation to do a Monte Carlo
// generation using Pythia8.
// See tutorials/proof/runProof.C, option "pythia8", for an
// example of how to run this selector.
//
//////////////////////////////////////////////////////////

#ifndef ProofSelectorRPVAnalysis_h
#define ProofSelectorRPVAnalysis_h

#include <TSelector.h>
#include <TTree.h>
#include <TFile.h>
#include <TProofOutputFile.h>

#include "NTFormat/interface/NTEvent.h"
#include "Plots/interface/SSDiLepAnaHistoManager.h"
#include "Tools/interface/Dataset.h"
#include "Tools/interface/AnalysisEnvironmentLoader.h"
#include "Selection/interface/SSDiLeptonSelection.h"

class TH1F;
class TBranch;
class TTree;

class AnalysisEnvironmentLoader;
class SSDiLeptonSelection;

class ProofSelectorRPVAnalysis : public TSelector {
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
   vector<string> CutName;
   vector<string> ChannelName;
   SSDiLeptonSelection sel; 
   float Luminosity;
   int verbosity;
   // 0: MC - 1: Data - 2 Data & MC
   int DataType;
   //Info analysis macro specific 
   TH1F* fHist;
   
   SSDiLepAnaHistoManager* histoManager;

  //------------------------------------
  // for PileUP reweighting
  //------------------------------------
  TH1D *  hPUData ;
  TH1F *  hPUMC   ; 
  
  float LumiError ;
  string PUWeightFileName;
  //PUWeighting  thePUReweighter;
  
  
  string datasetName ;

   ProofSelectorRPVAnalysis();
   virtual ~ProofSelectorRPVAnalysis();
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

   ClassDef(ProofSelectorRPVAnalysis,0);
};

#endif
