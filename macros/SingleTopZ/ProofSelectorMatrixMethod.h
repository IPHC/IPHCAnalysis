
//////////////////////////////////////////////////////////
//
// Example of TSelector implementation to do a Monte Carlo
// generation using Pythia8.
// See tutorials/proof/runProof.C, option "pythia8", for an
// example of how to run this selector.
//
//////////////////////////////////////////////////////////

#ifndef ProofSelectorMatrixMethod_h
#define ProofSelectorMatrixMethod_h

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
#include <TCanvas.h>
#include <TRandom.h>
#include <TLorentzVector.h>
//#include <iostream>

class TH1F;
class TBranch;
class TTree;

class AnalysisEnvironmentLoader;
class DiLeptonSelection;

class ProofSelectorMatrixMethod : public TSelector {
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
  
  bool applyDYscale ;
  bool applyFakescale ;
  
  //Here define Scale Factors
  //SF_trigger applied for mumu
  double SF_trig_mu ;  
  double SF_trig_emu;  
  double SF_trig_ee ;  
  
  double SF_trig_mu_error ;  
  double SF_trig_emu_error;  
  double SF_trig_ee_error ;  
  
  double SF_e  ;  
  double SF_mu ; 
    
  double SF_BranchingRatio_ll ; 
  double SF_BranchingRatio_lj ; 
  double SF_BranchingRatio_had ; 
    
  double SF_DY_ee ;
  double SF_DY_mm ;
  double SF_DY_em ;
  
  double  SF_DY_Njets_ee ;
  double  SF_DY_Njets_mm ;
  double  SF_DY_Njets_em ;
  
  /*double SF_met_mumu;
  double SF_met_emu;
  double SF_met_ee;

  double SF_met_mumu_error;
  double SF_met_emu_error;
  double SF_met_ee_error;*/
   
  
  std::vector<double> vSF_met_ee ;
  std::vector<double> vSF_met_mumu ;
  std::vector<double> vSF_met_emu ;
  
  std::vector<double> vSF_met_ee_error ;
  std::vector<double> vSF_met_mumu_error ;
  std::vector<double> vSF_met_emu_error ;
  
  
  bool PUup;
  bool PUdown;
  
  TRandom rand;
  
   
  std::vector<double> vSF_DY_ee ;
  std::vector<double> vSF_DY_mm ;
  std::vector<double> vSF_DY_em ;
   
  std::vector<double> vSF_DY_ee_error ;
  std::vector<double> vSF_DY_mm_error ;
  std::vector<double> vSF_DY_em_error ;
  
  std::vector<double> vSF_FakeBck_ee; 
  std::vector<double> vSF_FakeBck_mm; 
  std::vector<double> vSF_FakeBck_em; 
  
  
  std::vector<double> vSF_FakeBck_ee_error; 
  std::vector<double> vSF_FakeBck_mm_error; 
  std::vector<double> vSF_FakeBck_em_error; 
  
  double sumSFlept_ee;
  double sumSFlept_mumu;
  double sumSFlept_emu;
  
  double nEvents_ee;
  double nEvents_mumu;
  double nEvents_emu;
  
   
   
  double SF_Wjets_ee;
  double SF_Wjets_mm;
  double SF_Wjets_em;
  
  double SF_QCD_ee;
  double SF_QCD_mm;
  double SF_QCD_em;
  
  double scaleElec; // 1 to switch off
  double resolElec; // 0 to switch off
  
  
  bool ApplyLeptonSF;
  
  
  int ITypeMC ;
  int ICut    ;  
  
    
  std::vector<TH1F> CutFlow_mumumu_tight;
  std::vector<TH1F> CutFlow_mumue_tight;
  std::vector<TH1F> CutFlow_eemu_tight;
  std::vector<TH1F> CutFlow_eee_tight;
  std::vector<TH1F> ErrCutFlow_mumumu_tight;
  std::vector<TH1F> ErrCutFlow_mumue_tight;
  std::vector<TH1F> ErrCutFlow_eemu_tight;
  std::vector<TH1F> ErrCutFlow_eee_tight;
  
  std::vector<TH1F> CutFlow_mumumu_loose;
  std::vector<TH1F> CutFlow_mumue_loose;
  std::vector<TH1F> CutFlow_eemu_loose;
  std::vector<TH1F> CutFlow_eee_loose;
  std::vector<TH1F> ErrCutFlow_mumumu_loose;
  std::vector<TH1F> ErrCutFlow_mumue_loose;
  std::vector<TH1F> ErrCutFlow_eemu_loose;
  std::vector<TH1F> ErrCutFlow_eee_loose;
  
  std::vector<TH1F> ThirdLepPt_mumumu_tight_cut1;
  std::vector<TH1F> ThirdLepPt_mumue_tight_cut1;
  std::vector<TH1F> ThirdLepPt_eemu_tight_cut1;
  std::vector<TH1F> ThirdLepPt_eee_tight_cut1;

  std::vector<TH1F> ThirdLepPt_mumumu_tight_cut2;
  std::vector<TH1F> ThirdLepPt_mumue_tight_cut2;
  std::vector<TH1F> ThirdLepPt_eemu_tight_cut2;
  std::vector<TH1F> ThirdLepPt_eee_tight_cut2;

  std::vector<TH1F> ThirdLepPt_mumumu_tight_cut3;
  std::vector<TH1F> ThirdLepPt_mumue_tight_cut3;
  std::vector<TH1F> ThirdLepPt_eemu_tight_cut3;
  std::vector<TH1F> ThirdLepPt_eee_tight_cut3;
  
  std::vector<TH1F> ThirdLepPt_mumumu_tight_cut4;
  std::vector<TH1F> ThirdLepPt_mumue_tight_cut4;
  std::vector<TH1F> ThirdLepPt_eemu_tight_cut4;
  std::vector<TH1F> ThirdLepPt_eee_tight_cut4;

  std::vector<TH1F> ThirdLepPt_mumumu_tight_cut5;
  std::vector<TH1F> ThirdLepPt_mumue_tight_cut5;
  std::vector<TH1F> ThirdLepPt_eemu_tight_cut5;
  std::vector<TH1F> ThirdLepPt_eee_tight_cut5;

  std::vector<TH1F> ThirdLepPt_mumumu_tight_cut6;
  std::vector<TH1F> ThirdLepPt_mumue_tight_cut6;
  std::vector<TH1F> ThirdLepPt_eemu_tight_cut6;
  std::vector<TH1F> ThirdLepPt_eee_tight_cut6;

  ofstream ofile;
  
  //------------------------------------
  //Definition of the various histograms
  //------------------------------------
  int nbins ;
  float minx;
  float maxx;
  HistoManager MyhistoManager;
  
  
  //------------------------------------
  // for PileUP reweighting
  //------------------------------------
  PUWeighting  thePUReweighter;
  
  
  string datasetName ;
  
  
  //------------------------------------
  //definition of member functions
  //------------------------------------
  ProofSelectorMatrixMethod();
  virtual ~ProofSelectorMatrixMethod();
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    FillCutFlowHistos(int IChannel, string cand3leptonChannel, bool IsSignal,
                                     double weight, double EventYieldWeightError, int iCut, bool useLooseIso);
  virtual void    ApplySelectionAndFillHistos(IPHCTree::NTEvent *event, bool useLooseIso);
  virtual Bool_t  Process(Long64_t entry);
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();
  
  ClassDef(ProofSelectorMatrixMethod,0);
};

#endif
