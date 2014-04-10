
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
#include "Tools/interface/BtagHardcodedConditions.h"


#include "Tools/interface/PDFReweighter.h"
#include "Tools/interface/PUWeighting.h"


#include "JetCorrections/interface/JetCorrectionUncertainty.h"


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
  vector<string> TheChannelName;
  vector<string> VecChannelName;
  DiLeptonSelection sel; 
  float Luminosity;
  int verbosity;
  // 0: MC - 1: Data - 2 Data & MC
  int DataType;
  //Info analysis macro specific 
  
  
  int nselevents_mumumu;
  
  string datasetName ;
  //------------------------------------
  //For JES uncerainties
  //------------------------------------
  
  JetCorrectionUncertainty *theJESuncertainty;
  bool  doJESuncert;
  bool  upOrDown;
  bool  applyJES;
  bool  applyJER;
  float ResFactor;
  
  //------------------------------------
  //For PU reweighting
  //------------------------------------
  PUWeighting *LumiWeights;
  float LumiError ;
  string PUWeightFileName;
  bool IReweight;
  bool IReweight_puUp  ;
  bool IReweight_puDown;
  bool IReweight_pu;
  
   
   
  //------------------------------------
  //For PDF reweighting
  //------------------------------------
  bool doPDF ;
  int pdftype ;
  PDFReweighter pdf; 
  
 
  //------------------------------------
  //Definition of the various histograms
  //------------------------------------
  HistoManager MyhistoManager;
  int ITypeMC ;
  
    
  //------------------------------------
  //loose isolation on W cand
  //for background estimations
  //------------------------------------
  bool useNonIsoWcand;
  float looseIso;
  float themetcut;
   
 
  //------------------------------------
  //Trigger scale factors
  //------------------------------------
  
  bool applyTrigger ;
  bool applyTriggerUp ;
  bool applyTriggerDown ;
  
  double SF_trig_mumumu ;  
  double SF_trig_mumue;  
  double SF_trig_eemu ;
  double SF_trig_eee ;  
  
  double SF_trig_mumumu_error ;  
  double SF_trig_mumue_error;  
  double SF_trig_eemu_error ;  
  double SF_trig_eee_error ;
  
  //------------------------------------
  //Lepton scale factors
  //------------------------------------
  
  bool applyLeptonSF;
  bool applyLeptonSFUp;
  bool applyLeptonSFDown;
  
  
  //------------------------------------
  //background scale factors
  //------------------------------------
  bool applyWZ ;
  bool applyWZ_finalSel ;
  bool applyFakescale ;
  std::vector <double > SF_WZ;
  std::vector <double > SF_Fake;
  
  
  std::vector<TH1F> SelABjet;
  std::vector<TH1F> SelABjet_afterjetsel;
  
  std::vector<TH1F> Ntrilept_mumumu;
  std::vector<TH1F> Ntrileptnoniso_mumumu;
  
  
  std::vector<TH1F> CutFlow_mumumu;
  std::vector<TH1F> CutFlow_mumue;
  std::vector<TH1F> CutFlow_eemu;
  std::vector<TH1F> CutFlow_eee;
  std::vector<TH1F> ErrCutFlow_mumumu;
  std::vector<TH1F> ErrCutFlow_mumue;
  std::vector<TH1F> ErrCutFlow_eemu;
  std::vector<TH1F> ErrCutFlow_eee;
  
  std::vector<TH1F> PU_before_mumumu;
  std::vector<TH1F> PU_before_mumue;
  std::vector<TH1F> PU_before_eemu;
  std::vector<TH1F> PU_before_eee;

  std::vector<TH1F> PU_intime_mumumu;
  std::vector<TH1F> PU_intime_mumue;
  std::vector<TH1F> PU_intime_eemu;
  std::vector<TH1F> PU_intime_eee;
  
  std::vector<TH1F> PU_after_mumumu;
  std::vector<TH1F> PU_after_mumue;
  std::vector<TH1F> PU_after_eemu;
  std::vector<TH1F> PU_after_eee; 
  
  
  std::vector<TH1F> NVtx_mumumu;
  std::vector<TH1F> NVtx_mumue;
  std::vector<TH1F> NVtx_eemu;
  std::vector<TH1F> NVtx_eee; 
  
  std::vector<TH1F> NVtx_mumumu_aftertrigsel;
  std::vector<TH1F> NVtx_mumue_aftertrigsel;
  std::vector<TH1F> NVtx_eemu_aftertrigsel;
  std::vector<TH1F> NVtx_eee_aftertrigsel; 
  
  std::vector<TH1F> NVtx_mumumu_afterleptsel;
  std::vector<TH1F> NVtx_mumue_afterleptsel;
  std::vector<TH1F> NVtx_eemu_afterleptsel;
  std::vector<TH1F> NVtx_eee_afterleptsel; 
  
  std::vector<TH1F> DijetInvM_mumumu_afterleptsel_inZpeak;
  std::vector<TH1F> DijetInvM_mumue_afterleptsel_inZpeak;
  std::vector<TH1F> DijetInvM_eemu_afterleptsel_inZpeak;
  std::vector<TH1F> DijetInvM_eee_afterleptsel_inZpeak;
  
  
  
  std::vector<TH1F> Mt_mumumu_afterbjetsel;
  std::vector<TH1F> Mt_mumue_afterbjetsel;
  std::vector<TH1F> Mt_eemu_afterbjetsel;
  std::vector<TH1F> Mt_eee_afterbjetsel;
   
  std::vector<TH1F> Mt_mumumu_afterbjetveto;
  std::vector<TH1F> Mt_mumue_afterbjetveto;
  std::vector<TH1F> Mt_eemu_afterbjetveto;
  std::vector<TH1F> Mt_eee_afterbjetveto;
   
  
  std::vector<TH1F> NJet_mumumu_afterZsel;
  std::vector<TH1F> NJet_mumue_afterZsel;
  std::vector<TH1F> NJet_eemu_afterZsel;
  std::vector<TH1F> NJet_eee_afterZsel;
  
  
  std::vector<TH1F> NJet_mumumu_afterbsel;
  std::vector<TH1F> NJet_mumue_afterbsel;
  std::vector<TH1F> NJet_eemu_afterbsel;
  std::vector<TH1F> NJet_eee_afterbsel;
  
  std::vector<TH1F> NJet_mumumu_afterleptsel_mWT110;
  std::vector<TH1F> NJet_mumue_afterleptsel_mWT110;
  std::vector<TH1F> NJet_eemu_afterleptsel_mWT110;
  std::vector<TH1F> NJet_eee_afterleptsel_mWT110;
  
  
  std::vector<TH1F> NLept_mumumu_afterbsel;
  std::vector<TH1F> NLept_mumue_afterbsel;
  std::vector<TH1F> NLept_eemu_afterbsel;
  std::vector<TH1F> NLept_eee_afterbsel;
  
  
  std::vector<TH1F> NBJet_mumumu_afterZsel;
  std::vector<TH1F> NBJet_mumue_afterZsel;
  std::vector<TH1F> NBJet_eemu_afterZsel;
  std::vector<TH1F> NBJet_eee_afterZsel;
  
  
  std::vector<TH1F> NBJet_mumumu_afterjetsel;
  std::vector<TH1F> NBJet_mumue_afterjetsel;
  std::vector<TH1F> NBJet_eemu_afterjetsel;
  std::vector<TH1F> NBJet_eee_afterjetsel;
  
  
  
  
  std::vector<TH1F> NBJet_mumumu_afterjetsel_bjets;
  std::vector<TH1F> NBJet_mumue_afterjetsel_bjets ;
  std::vector<TH1F> NBJet_eemu_afterjetsel_bjets  ;
  std::vector<TH1F> NBJet_eee_afterjetsel_bjets	 ;
  
  
  std::vector<TH1F> NBJet_mumumu_afterjetsel_cjets;
  std::vector<TH1F> NBJet_mumue_afterjetsel_cjets ;
  std::vector<TH1F> NBJet_eemu_afterjetsel_cjets  ;
  std::vector<TH1F> NBJet_eee_afterjetsel_cjets	 ;
  
  
  std::vector<TH1F> NBJet_mumumu_afterjetsel_ljets;
  std::vector<TH1F> NBJet_mumue_afterjetsel_ljets ;
  std::vector<TH1F> NBJet_eemu_afterjetsel_ljets  ;
  std::vector<TH1F> NBJet_eee_afterjetsel_ljets	 ;
  
  
  
 
  
  
  
  
  
  
  
  
  std::vector<TH1F> BJetDiscri_mumumu_afterjetsel_bjets;
  std::vector<TH1F> BJetDiscri_mumue_afterjetsel_bjets;
  std::vector<TH1F> BJetDiscri_eemu_afterjetsel_bjets;
  std::vector<TH1F> BJetDiscri_eee_afterjetsel_bjets;
  
  std::vector<TH1F> BJetDiscri_mumumu_afterjetsel_cjets;
  std::vector<TH1F> BJetDiscri_mumue_afterjetsel_cjets;
  std::vector<TH1F> BJetDiscri_eemu_afterjetsel_cjets;
  std::vector<TH1F> BJetDiscri_eee_afterjetsel_cjets;
  
  std::vector<TH1F> BJetDiscri_mumumu_afterjetsel_ljets;
  std::vector<TH1F> BJetDiscri_mumue_afterjetsel_ljets;
  std::vector<TH1F> BJetDiscri_eemu_afterjetsel_ljets;
  std::vector<TH1F> BJetDiscri_eee_afterjetsel_ljets;
  
  std::vector<TH1F> NBJet_mumumu_afterleptsel_mWT110;
  std::vector<TH1F> NBJet_mumue_afterleptsel_mWT110;
  std::vector<TH1F> NBJet_eemu_afterleptsel_mWT110;
  std::vector<TH1F> NBJet_eee_afterleptsel_mWT110;
  
  //to be filled
  
  std::vector<TH1F> Nvtx_mumumu_afterleptsel;
  std::vector<TH1F> Nvtx_mumue_afterleptsel;
  std::vector<TH1F> Nvtx_eemu_afterleptsel;
  std::vector<TH1F> Nvtx_eee_afterleptsel;
  
  std::vector<TH1F> InvM_ll_mumumu_afterleptsel;
  std::vector<TH1F> InvM_ll_mumue_afterleptsel;
  std::vector<TH1F> InvM_ll_eemu_afterleptsel;
  std::vector<TH1F> InvM_ll_eee_afterleptsel;
  
  std::vector<TH1F> InvM_ll_mumumu_afterleptsel_mWT110;
  std::vector<TH1F> InvM_ll_mumue_afterleptsel_mWT110;
  std::vector<TH1F> InvM_ll_eemu_afterleptsel_mWT110;
  std::vector<TH1F> InvM_ll_eee_afterleptsel_mWT110;
  
  std::vector<TH1F> InvM_ll_mumumu_afterleptsel_lowbin;
  std::vector<TH1F> InvM_ll_mumue_afterleptsel_lowbin;
  std::vector<TH1F> InvM_ll_eemu_afterleptsel_lowbin;
  std::vector<TH1F> InvM_ll_eee_afterleptsel_lowbin;
  
  std::vector<TH1F> InvM_ll_mumumu_afterleptsel_highSumPt;
  std::vector<TH1F> InvM_ll_mumue_afterleptsel_highSumPt;
  std::vector<TH1F> InvM_ll_eemu_afterleptsel_highSumPt;
  std::vector<TH1F> InvM_ll_eee_afterleptsel_highSumPt;
  
  std::vector<TH1F> InvM_ll_mumumu_afterjetsel;
  std::vector<TH1F> InvM_ll_mumue_afterjetsel;
  std::vector<TH1F> InvM_ll_eemu_afterjetsel;
  std::vector<TH1F> InvM_ll_eee_afterjetsel;
  
  std::vector<TH1F> InvM_ll_mumumu_afterbjetsel;
  std::vector<TH1F> InvM_ll_mumue_afterbjetsel;
  std::vector<TH1F> InvM_ll_eemu_afterbjetsel;
  std::vector<TH1F> InvM_ll_eee_afterbjetsel;
  
  std::vector<TH1F> LeptPt_mumumu_afterleptsel;
  std::vector<TH1F> LeptPt_mumue_afterleptsel;
  std::vector<TH1F> LeptPt_eemu_afterleptsel;
  std::vector<TH1F> LeptPt_eee_afterleptsel;
  
  std::vector<TH1F> LeptPt_mumumu_afterjetsel;
  std::vector<TH1F> LeptPt_mumue_afterjetsel;
  std::vector<TH1F> LeptPt_eemu_afterjetsel;
  std::vector<TH1F> LeptPt_eee_afterjetsel;
  
  std::vector<TH1F> LeptPt_mumumu_afterbjetsel;
  std::vector<TH1F> LeptPt_mumue_afterbjetsel;
  std::vector<TH1F> LeptPt_eemu_afterbjetsel;
  std::vector<TH1F> LeptPt_eee_afterbjetsel;
  
  
  
  std::vector<TH1F> LeptZPt_mumumu_afterleptsel;
  std::vector<TH1F> LeptZPt_mumue_afterleptsel;
  std::vector<TH1F> LeptZPt_eemu_afterleptsel;
  std::vector<TH1F> LeptZPt_eee_afterleptsel;
  
  std::vector<TH1F> LeptZPt_mumumu_afterjetsel;
  std::vector<TH1F> LeptZPt_mumue_afterjetsel;
  std::vector<TH1F> LeptZPt_eemu_afterjetsel;
  std::vector<TH1F> LeptZPt_eee_afterjetsel;
  
  std::vector<TH1F> LeptZPt_mumumu_afterbjetsel;
  std::vector<TH1F> LeptZPt_mumue_afterbjetsel;
  std::vector<TH1F> LeptZPt_eemu_afterbjetsel;
  std::vector<TH1F> LeptZPt_eee_afterbjetsel;
  
  
  
  std::vector<TH1F> LeptWPt_mumumu_afterleptsel;
  std::vector<TH1F> LeptWPt_mumue_afterleptsel;
  std::vector<TH1F> LeptWPt_eemu_afterleptsel;
  std::vector<TH1F> LeptWPt_eee_afterleptsel;
  
  std::vector<TH1F> LeptWPt_mumumu_afterjetsel;
  std::vector<TH1F> LeptWPt_mumue_afterjetsel;
  std::vector<TH1F> LeptWPt_eemu_afterjetsel;
  std::vector<TH1F> LeptWPt_eee_afterjetsel;
  
  std::vector<TH1F> LeptWPt_mumumu_afterbjetsel;
  std::vector<TH1F> LeptWPt_mumue_afterbjetsel;
  std::vector<TH1F> LeptWPt_eemu_afterbjetsel;
  std::vector<TH1F> LeptWPt_eee_afterbjetsel;
  
  std::vector<TH1F> LeptWPt_mumumu_afterbjetveto;
  std::vector<TH1F> LeptWPt_mumue_afterbjetveto;
  std::vector<TH1F> LeptWPt_eemu_afterbjetveto;
  std::vector<TH1F> LeptWPt_eee_afterbjetveto;
  
  
  std::vector<TH1F> LeptWPt_mumumu_afterleptsel_mWT110;
  std::vector<TH1F> LeptWPt_mumue_afterleptsel_mWT110;
  std::vector<TH1F> LeptWPt_eemu_afterleptsel_mWT110;
  std::vector<TH1F> LeptWPt_eee_afterleptsel_mWT110;
  
  
  std::vector<TH1F> JetPt_mumumu_afterleptsel;
  std::vector<TH1F> JetPt_mumue_afterleptsel;
  std::vector<TH1F> JetPt_eemu_afterleptsel;
  std::vector<TH1F> JetPt_eee_afterleptsel;
  
  std::vector<TH1F> JetPt_mumumu_afterjetsel;
  std::vector<TH1F> JetPt_mumue_afterjetsel;
  std::vector<TH1F> JetPt_eemu_afterjetsel;
  std::vector<TH1F> JetPt_eee_afterjetsel;
  
  std::vector<TH1F> JetPt_mumumu_afterbjetsel;
  std::vector<TH1F> JetPt_mumue_afterbjetsel;
  std::vector<TH1F> JetPt_eemu_afterbjetsel;
  std::vector<TH1F> JetPt_eee_afterbjetsel;
  
  std::vector<TH1F> JetPt_mumumu_afterbjetveto;
  std::vector<TH1F> JetPt_mumue_afterbjetveto;
  std::vector<TH1F> JetPt_eemu_afterbjetveto;
  std::vector<TH1F> JetPt_eee_afterbjetveto;
  
  std::vector<TH1F> JetEta_mumumu_afterleptsel;
  std::vector<TH1F> JetEta_mumue_afterleptsel;
  std::vector<TH1F> JetEta_eemu_afterleptsel;
  std::vector<TH1F> JetEta_eee_afterleptsel;
  
  std::vector<TH1F> JetEta_mumumu_afterjetsel;
  std::vector<TH1F> JetEta_mumue_afterjetsel;
  std::vector<TH1F> JetEta_eemu_afterjetsel;
  std::vector<TH1F> JetEta_eee_afterjetsel;
  
  std::vector<TH1F> JetEta_mumumu_afterbjetsel;
  std::vector<TH1F> JetEta_mumue_afterbjetsel;
  std::vector<TH1F> JetEta_eemu_afterbjetsel;
  std::vector<TH1F> JetEta_eee_afterbjetsel;
  
  std::vector<TH1F> JetEta_mumumu_afterbjetveto;
  std::vector<TH1F> JetEta_mumue_afterbjetveto;
  std::vector<TH1F> JetEta_eemu_afterbjetveto;
  std::vector<TH1F> JetEta_eee_afterbjetveto;
  
  std::vector<TH1F> HT_mumumu_afterleptsel;
  std::vector<TH1F> HT_mumue_afterleptsel;
  std::vector<TH1F> HT_eemu_afterleptsel;
  std::vector<TH1F> HT_eee_afterleptsel;
  
  
  std::vector<TH1F> HT_mumumu_afterjetsel;
  std::vector<TH1F> HT_mumue_afterjetsel;
  std::vector<TH1F> HT_eemu_afterjetsel;
  std::vector<TH1F> HT_eee_afterjetsel;
  
  std::vector<TH1F> HT_mumumu_afterbjetsel;
  std::vector<TH1F> HT_mumue_afterbjetsel;
  std::vector<TH1F> HT_eemu_afterbjetsel;
  std::vector<TH1F> HT_eee_afterbjetsel;
  
  std::vector<TH1F> HT_mumumu_afterbjetveto;
  std::vector<TH1F> HT_mumue_afterbjetveto;
  std::vector<TH1F> HT_eemu_afterbjetveto;
  std::vector<TH1F> HT_eee_afterbjetveto;
  
  
  
  
  
  std::vector<TH1F> MET_mumumu_afterleptsel;
  std::vector<TH1F> MET_mumue_afterleptsel;
  std::vector<TH1F> MET_eemu_afterleptsel;
  std::vector<TH1F> MET_eee_afterleptsel;
  
  
  
  std::vector<TH1F> MET_mumumu_afterleptsel_mWT110;
  std::vector<TH1F> MET_mumue_afterleptsel_mWT110;
  std::vector<TH1F> MET_eemu_afterleptsel_mWT110;
  std::vector<TH1F> MET_eee_afterleptsel_mWT110;
  
  
  std::vector<TH1F> MET_mumumu_afterjetsel;
  std::vector<TH1F> MET_mumue_afterjetsel;
  std::vector<TH1F> MET_eemu_afterjetsel;
  std::vector<TH1F> MET_eee_afterjetsel;
  
  std::vector<TH1F> MET_mumumu_afterbjetsel;
  std::vector<TH1F> MET_mumue_afterbjetsel;
  std::vector<TH1F> MET_eemu_afterbjetsel;
  std::vector<TH1F> MET_eee_afterbjetsel;
  
  std::vector<TH1F> Asym_mumumu_afterbjetsel;
  std::vector<TH1F> Asym_mumue_afterbjetsel;
  std::vector<TH1F> Asym_eemu_afterbjetsel;
  std::vector<TH1F> Asym_eee_afterbjetsel;
  
  
  
  std::vector<TH1F> mWT_mumumu_afterjetsel;
  std::vector<TH1F> mWT_mumue_afterjetsel;
  std::vector<TH1F> mWT_eemu_afterjetsel;
  std::vector<TH1F> mWT_eee_afterjetsel;
  
  
  
  std::vector<TH1F> RecoPtZ_mumumu_afterbjetsel;
  std::vector<TH1F> RecoPtZ_mumue_afterbjetsel;
  std::vector<TH1F> RecoPtZ_eemu_afterbjetsel;
  std::vector<TH1F> RecoPtZ_eee_afterbjetsel;
  
  std::vector<TH1F> RecoPtZ_mumumu_afterbjetveto;
  std::vector<TH1F> RecoPtZ_mumue_afterbjetveto;
  std::vector<TH1F> RecoPtZ_eemu_afterbjetveto;
  std::vector<TH1F> RecoPtZ_eee_afterbjetveto;
  
  
  std::vector<TH1F> RecoPtZ_mumumu_afterleptsel;
  std::vector<TH1F> RecoPtZ_mumue_afterleptsel;
  std::vector<TH1F> RecoPtZ_eemu_afterleptsel;
  std::vector<TH1F> RecoPtZ_eee_afterleptsel;
  
  
  std::vector<TH1F> RecoPtZ_mumumu_afterleptsel_nojet;
  std::vector<TH1F> RecoPtZ_mumue_afterleptsel_nojet;
  std::vector<TH1F> RecoPtZ_eemu_afterleptsel_nojet;
  std::vector<TH1F> RecoPtZ_eee_afterleptsel_nojet;
  
  
  std::vector<TH1F> RecoTopMass_mumumu_afterbjetsel;
  std::vector<TH1F> RecoTopMass_mumue_afterbjetsel;
  std::vector<TH1F> RecoTopMass_eemu_afterbjetsel;
  std::vector<TH1F> RecoTopMass_eee_afterbjetsel;
  
  std::vector<TH1F> RecoTopMass_mumumu_afterbjetveto;
  std::vector<TH1F> RecoTopMass_mumue_afterbjetveto;
  std::vector<TH1F> RecoTopMass_eemu_afterbjetveto;
  std::vector<TH1F> RecoTopMass_eee_afterbjetveto;
  
  
  std::vector<TH1F> deltaPhilb_mumumu_afterbjetsel;
  std::vector<TH1F> deltaPhilb_mumue_afterbjetsel;
  std::vector<TH1F> deltaPhilb_eemu_afterbjetsel;
  std::vector<TH1F> deltaPhilb_eee_afterbjetsel;
  
  std::vector<TH1F> deltaPhilj_mumumu_afterbjetveto;
  std::vector<TH1F> deltaPhilj_mumue_afterbjetveto;
  std::vector<TH1F> deltaPhilj_eemu_afterbjetveto;
  std::vector<TH1F> deltaPhilj_eee_afterbjetveto;
  
  
  
  
  std::vector<TH1F> deltaR_mumumu_afterleptsel;
  std::vector<TH1F> deltaR_mumue_afterleptsel;
  std::vector<TH1F> deltaR_eemu_afterleptsel;
  std::vector<TH1F> deltaR_eee_afterleptsel;
  
  
  
  std::vector<TH1F> deltaRLeptJet_mumumu_afterleptsel_mWT110;
  std::vector<TH1F> deltaRLeptJet_mumue_afterleptsel_mWT110;
  std::vector<TH1F> deltaRLeptJet_eemu_afterleptsel_mWT110;
  std::vector<TH1F> deltaRLeptJet_eee_afterleptsel_mWT110;
  
  std::vector<TH1F> deltaRLeptMet_mumumu_afterleptsel_mWT110;
  std::vector<TH1F> deltaRLeptMet_mumue_afterleptsel_mWT110;
  std::vector<TH1F> deltaRLeptMet_eemu_afterleptsel_mWT110;
  std::vector<TH1F> deltaRLeptMet_eee_afterleptsel_mWT110;
  
  
  
  std::vector<TH1F> WmissAssing_mumumu_afterleptsel;
  std::vector<TH1F> WmissAssing_mumue_afterleptsel;
  std::vector<TH1F> WmissAssing_eemu_afterleptsel;
  std::vector<TH1F> WmissAssing_eee_afterleptsel;
  
  
  std::vector<TH1F> mWT_mumumu_afterleptsel;
  std::vector<TH1F> mWT_mumue_afterleptsel;
  std::vector<TH1F> mWT_eemu_afterleptsel;
  std::vector<TH1F> mWT_eee_afterleptsel;
  
  
  std::vector<TH1F> mWT_mumumu_afterbjetsel;
  std::vector<TH1F> mWT_mumue_afterbjetsel;
  std::vector<TH1F> mWT_eemu_afterbjetsel;
  std::vector<TH1F> mWT_eee_afterbjetsel;
  
  std::vector<TH1F> mWT_mumumu_afterbjetveto;
  std::vector<TH1F> mWT_mumue_afterbjetveto;
  std::vector<TH1F> mWT_eemu_afterbjetveto;
  std::vector<TH1F> mWT_eee_afterbjetveto;
  
  
  std::vector<TH1F> Charge_mumumu_afterleptsel;
  std::vector<TH1F> Charge_mumue_afterleptsel;
  std::vector<TH1F> Charge_eemu_afterleptsel;
  std::vector<TH1F> Charge_eee_afterleptsel;
  
  std::vector<TH1F> Charge_mumumu_afterleptsel_mWT110;
  std::vector<TH1F> Charge_mumue_afterleptsel_mWT110;
  std::vector<TH1F> Charge_eemu_afterleptsel_mWT110;
  std::vector<TH1F> Charge_eee_afterleptsel_mWT110;
  
   
  std::vector<TH1F> Nvertex;
 
  std::vector<TH2D> InvM_ll_vs_mWT_mumumu_afterleptsel;
  std::vector<TH2D> InvM_ll_vs_mWT_mumue_afterleptsel;
  std::vector<TH2D> InvM_ll_vs_mWT_eemu_afterleptsel;
  std::vector<TH2D> InvM_ll_vs_mWT_eee_afterleptsel;
  
  std::vector<TH2D> HT_vs_MET_mumumu_afterleptsel;
  std::vector<TH2D> HT_vs_MET_mumue_afterleptsel;
  std::vector<TH2D> HT_vs_MET_eemu_afterleptsel;
  std::vector<TH2D> HT_vs_MET_eee_afterleptsel;
  
  std::vector<TH2D> HT_vs_NJet_mumumu_afterleptsel;
  std::vector<TH2D> HT_vs_NJet_mumue_afterleptsel;
  std::vector<TH2D> HT_vs_NJet_eemu_afterleptsel;
  std::vector<TH2D> HT_vs_NJet_eee_afterleptsel;
  
  std::vector<TH2D> HT_vs_NBJet_mumumu_afterleptsel;
  std::vector<TH2D> HT_vs_NBJet_mumue_afterleptsel;
  std::vector<TH2D> HT_vs_NBJet_eemu_afterleptsel;
  std::vector<TH2D> HT_vs_NBJet_eee_afterleptsel;
  
  std::vector<TH2D> HT_vs_LeptPt_mumumu_afterleptsel;
  std::vector<TH2D> HT_vs_LeptPt_mumue_afterleptsel;
  std::vector<TH2D> HT_vs_LeptPt_eemu_afterleptsel;
  std::vector<TH2D> HT_vs_LeptPt_eee_afterleptsel;
  
  std::vector<TH2D> HT_vs_JetPt_mumumu_afterleptsel;
  std::vector<TH2D> HT_vs_JetPt_mumue_afterleptsel;
  std::vector<TH2D> HT_vs_JetPt_eemu_afterleptsel;
  std::vector<TH2D> HT_vs_JetPt_eee_afterleptsel;
  
  
  
  std::vector<TH2D> HT_vs_Mll_mumumu_afterleptsel;
  std::vector<TH2D> HT_vs_Mll_mumue_afterleptsel;
  std::vector<TH2D> HT_vs_Mll_eemu_afterleptsel;
  std::vector<TH2D> HT_vs_Mll_eee_afterleptsel;
  
   
  
  std::vector<TH1F> *pCutFlow;
  std::vector<TH1F> *pErrCutFlow;
       
  std::vector<TH1F> *pDijetInvM_afterleptsel_inZpeak;
  std::vector<TH1F> *pMt_afterbjetsel;
  std::vector<TH1F> *pMt_afterbjetveto;
  std::vector<TH1F> *pNJet_afterZsel;  
  std::vector<TH1F> *pNJet_afterbsel;
  std::vector<TH1F> *pNJet_afterleptsel_mWT110;
  std::vector<TH1F> *pNLept_afterbsel;
  std::vector<TH1F> *pNBJet_afterZsel;
  std::vector<TH1F> *pNBJet_afterjetsel;
  std::vector<TH1F> *pNBJet_afterjetsel_bjets;
  std::vector<TH1F> *pNBJet_afterjetsel_cjets;
  std::vector<TH1F> *pNBJet_afterjetsel_ljets;
  std::vector<TH1F> *pBJetDiscri_afterjetsel_bjets;
  std::vector<TH1F> *pBJetDiscri_afterjetsel_cjets;
  std::vector<TH1F> *pBJetDiscri_afterjetsel_ljets;
  std::vector<TH1F> *pNBJet_afterleptsel_mWT110;
  std::vector<TH1F> *pNvtx_afterleptsel;
  std::vector<TH1F> *pInvM_ll_afterleptsel;  
  std::vector<TH1F> *pInvM_ll_afterleptsel_mWT110;  
  std::vector<TH1F> *pInvM_ll_afterleptsel_lowbin;  
  std::vector<TH1F> *pInvM_ll_afterleptsel_highSumPt;  
  std::vector<TH1F> *pInvM_ll_afterjetsel;
  std::vector<TH1F> *pInvM_ll_afterbjetsel;  
  std::vector<TH1F> *pLeptPt_afterleptsel;  
  std::vector<TH1F> *pLeptPt_afterjetsel;  
  std::vector<TH1F> *pLeptPt_afterbjetsel;  
  std::vector<TH1F> *pLeptZPt_afterleptsel;  
  std::vector<TH1F> *pLeptZPt_afterjetsel;  
  std::vector<TH1F> *pLeptZPt_afterbjetsel;
  std::vector<TH1F> *pLeptWPt_afterleptsel;  
  std::vector<TH1F> *pLeptWPt_afterjetsel;
  std::vector<TH1F> *pLeptWPt_afterbjetsel;
  std::vector<TH1F> *pLeptWPt_afterbjetveto;
  std::vector<TH1F> *pLeptWPt_afterleptsel_mWT110;
  std::vector<TH1F> *pJetPt_afterleptsel;
  std::vector<TH1F> *pJetPt_afterjetsel;
  std::vector<TH1F> *pJetPt_afterbjetsel;
  std::vector<TH1F> *pJetPt_afterbjetveto;  
  std::vector<TH1F> *pJetEta_afterleptsel;
  std::vector<TH1F> *pJetEta_afterjetsel;  
  std::vector<TH1F> *pJetEta_afterbjetsel;  
  std::vector<TH1F> *pJetEta_afterbjetveto;
  std::vector<TH1F> *pHT_afterleptsel;  
  std::vector<TH1F> *pHT_afterjetsel;
  std::vector<TH1F> *pHT_afterbjetsel;  
  std::vector<TH1F> *pHT_afterbjetveto;  
  std::vector<TH1F> *pMET_afterleptsel;
  std::vector<TH1F> *pMET_afterleptsel_mWT110;
  std::vector<TH1F> *pMET_afterjetsel;
  std::vector<TH1F> *pMET_afterbjetsel;
  std::vector<TH1F> *pAsym_afterbjetsel;  
  std::vector<TH1F> *pmWT_afterjetsel;
  std::vector<TH1F> *pRecoPtZ_afterbjetsel;
  std::vector<TH1F> *pRecoPtZ_afterbjetveto;
  std::vector<TH1F> *pRecoPtZ_afterleptsel;  
  std::vector<TH1F> *pRecoPtZ_afterleptsel_nojet;
  std::vector<TH1F> *pRecoTopMass_afterbjetsel;
  std::vector<TH1F> *pRecoTopMass_afterbjetveto;  
  std::vector<TH1F> *pdeltaPhilb_afterbjetsel;
  std::vector<TH1F> *pdeltaPhilj_afterbjetveto;  
  std::vector<TH1F> *pdeltaR_afterleptsel;  
  std::vector<TH1F> *pdeltaRLeptJet_afterleptsel_mWT110;  
  std::vector<TH1F> *pdeltaRLeptMet_afterleptsel_mWT110;  
  std::vector<TH1F> *pWmissAssing_afterleptsel;  
  std::vector<TH1F> *pmWT_afterleptsel; 
  std::vector<TH1F> *pmWT_afterbjetsel;  
  std::vector<TH1F> *pmWT_afterbjetveto;  
  std::vector<TH1F> *pCharge_afterleptsel;  
  std::vector<TH1F> *pCharge_afterleptsel_mWT110;   
  std::vector<TH1F> *pNvertex; 
  std::vector<TH2D> *pInvM_ll_vs_mWT_afterleptsel;  
  std::vector<TH2D> *pHT_vs_MET_afterleptsel;  
  std::vector<TH2D> *pHT_vs_NJet_afterleptsel;  
  std::vector<TH2D> *pHT_vs_NBJet_afterleptsel;  
  std::vector<TH2D> *pHT_vs_LeptPt_afterleptsel;
  std::vector<TH2D> *pHT_vs_JetPt_afterleptsel;
  std::vector<TH2D> *pHT_vs_Mll_afterleptsel;
  
  
  
  
  
  
  
  //------------------------------------
  //BTag scale factors
  //------------------------------------
  
  BtagHardcodedConditions BTagSFManager;
  
   //------------------------------------
  //TTree and banches used for BDT
  //------------------------------------
  TTree *TheTree;
  
  
  float tree_topMass;
  float tree_totMass;
  float tree_deltaPhilb;
  float tree_deltaRlb;
  float tree_deltaRTopZ;
  float tree_asym;
  float tree_Zpt;
  float tree_ZEta;
  float tree_topPt;
  float tree_topEta;
  
  float tree_deltaRZl;
  float tree_deltaPhiZmet;
  float tree_btagDiscri;
  float tree_NJets;
  float tree_NBJets;
  float tree_totPt;
  float tree_totEta;
  float tree_met;
  float tree_mTW;
  
  float tree_leptWPt, tree_leptWEta;
  float tree_leadJetPt, tree_leadJetEta;
  float tree_deltaRZleptW, tree_deltaPhiZleptW;
   
  int   tree_SampleType;
  int   tree_Channel;
  
  float tree_EvtWeight;
  
  
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
  
  //to determine the MC event weight
  std::vector< double > determineWeights(TString, double , double);
  
  void createTheHisto(HistoManager *thehistomanag);
  void WriteTheHisto(TFile* theoutputfile, HistoManager *thehistomanag);
  void cleanHistoVector();
  
  
  void determineLeptonCandidates(
  		bool UseLooseWcand, float looseIsoCut, double rhocorr,
  		std::vector<NTElectron> *selE,        std::vector<NTMuon> *selM, 
  		std::vector<NTElectron> *selENonIso,  std::vector<NTMuon> *selMNonIso, 
		std::vector<NTElectron> *theZeeCand,  std::vector<NTMuon> *theZmumuCand, 
		std::vector<NTElectron> *theWeCand,   std::vector<NTMuon> *theWmuCand
		); 
  
  void defineHistoPointer(int channel);
  
};

#endif
