
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
  bool applyTriggerUp ;
  bool applyTriggerDown ;
  
  bool applyWZ ;
  bool applyWZ_finalSel ;
  bool applyFakescale ;

  bool applyLeptonSF;
  bool applyLeptonSFUp;
  bool applyLeptonSFDown;
  
  
  double SF_WZ_finalSel;
  double looseIso;
  
  std::vector <double > SF_WZ;
  std::vector <double > SF_Fake;
  
  
  bool doBTagCVScorr ;
  int doBTagCSV_syst;
  
  //Here define Scale Factors
  //SF_trigger applied for mumu
  double SF_trig_mumumu ;  
  double SF_trig_mumue;  
  double SF_trig_eemu ;
  double SF_trig_eee ;  
  
  double SF_trig_mumumu_error ;  
  double SF_trig_mumue_error;  
  double SF_trig_eemu_error ;  
  double SF_trig_eee_error ;
  
  double SF_BranchingRatio_ll ; 
  double SF_BranchingRatio_lj ; 
  double SF_BranchingRatio_had ; 
    
   
  bool  applyJES;
  float scale;
  bool  applyJER;
  float ResFactor;
   
  
  bool PUup;
  bool PUdown;
  
  TRandom rand;
  
  double themetcut;
  
  bool doPDF ;
  int pdftype ;
  PDFReweighter pdf; 
  
  
  
  double sumSFlept_mumumu;
  double sumSFlept_mumue;
  double sumSFlept_eemu;
  double sumSFlept_eee;
  
  double nEvents_mumumu;
  double nEvents_mumue;
  double nEvents_eemu;
  double nEvents_eee;
  
   
   
  double SF_Wjets_ee;
  double SF_Wjets_mm;
  double SF_Wjets_em;
  
  double SF_QCD_ee;
  double SF_QCD_mm;
  double SF_QCD_em;
  
  double scaleElec; // 1 to switch off
  double resolElec; // 0 to switch off
  
  
  
  
  int ITypeMC ;
  int ICut    ;  
  
    
  
  TH1F* fHist;
  
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
  
  
  //TTree and banches used for BDT
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
  
  float tree_leptWPt, tree_leptWEta;
  float tree_leadJetPt, tree_leadJetEta;
  float tree_deltaRZleptW, tree_deltaPhiZleptW;
   
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
  //BTag scale factors
  //------------------------------------
  
  BtagHardcodedConditions BTagSFManager;
  
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
  
  
  std::vector<double>  GetNvertexWeight(TString datasetName);
  
  std::vector<double>  TableToVector(double * theTable, int size);
  
  
  
double functionDiscriLjet(double CVSvalue){
  
  double a = 0.958359;
  double b = 1.06131;
  double c = -1.9045;
  double d = 0.884829;
  
  double corrected_cvs = a*CVSvalue +b*pow(CVSvalue,2)  +c*pow(CVSvalue,3) +d*pow(CVSvalue,4) ;
      return corrected_cvs;
}

double functionDiscriCjet(double CVSvalue){
  double a = 1.11389;
  double b = -0.538485;
  double c = 0.428048;
  double d = 0;
  double corrected_cvs = a*CVSvalue +b*pow(CVSvalue,2)  +c*pow(CVSvalue,3) +d*pow(CVSvalue,4) ;
      return corrected_cvs;
}

double functionDiscriBjet(double CVSvalue){
  double a = 0.942535;
  double b = -0.264954;
  double c = 0.328049;
  double d = 0;
  double corrected_cvs = a*CVSvalue +b*pow(CVSvalue,2)  +c*pow(CVSvalue,3) +d*pow(CVSvalue,4) ;
      return corrected_cvs;
}

  
    
double functionDiscriLjet_up(double CVSvalue){
  
  double a = 1.01955;
double b = 1.06079;
double c = -2.02085;
double d = 0.940512;

  double corrected_cvs = a*CVSvalue +b*pow(CVSvalue,2)  +c*pow(CVSvalue,3) +d*pow(CVSvalue,4) ;
      return corrected_cvs;
}

double functionDiscriCjet_up(double CVSvalue){
double a = 1.12734;
double b = -0.374471;
double c = 0.24919;
double d = 0;

  double corrected_cvs = a*CVSvalue +b*pow(CVSvalue,2)  +c*pow(CVSvalue,3) +d*pow(CVSvalue,4) ;
      return corrected_cvs;
}

double functionDiscriBjet_up(double CVSvalue){
double a = 1.12734;
double b = -0.374471;
double c = 0.24919;
double d = 0;

  double corrected_cvs = a*CVSvalue +b*pow(CVSvalue,2)  +c*pow(CVSvalue,3) +d*pow(CVSvalue,4) ;
      return corrected_cvs;
}

  
  
  
    
double functionDiscriLjet_down(double CVSvalue){
  double a = 0.979725;
double b = 0.528246;
double c = -0.838917;
double d = 0.330947;

  
  double corrected_cvs = a*CVSvalue +b*pow(CVSvalue,2)  +c*pow(CVSvalue,3) +d*pow(CVSvalue,4) ;
      return corrected_cvs;
}

double functionDiscriCjet_down(double CVSvalue){
double a = 1.01008;
double b = -0.422312;
double c = 0.418221;
double d = 0;


  double corrected_cvs = a*CVSvalue +b*pow(CVSvalue,2)  +c*pow(CVSvalue,3) +d*pow(CVSvalue,4) ;
      return corrected_cvs;
}

double functionDiscriBjet_down(double CVSvalue){

double a = 0.777746;
double b = -0.209792;
double c = 0.440256;
double d = 0;

  double corrected_cvs = a*CVSvalue +b*pow(CVSvalue,2)  +c*pow(CVSvalue,3) +d*pow(CVSvalue,4) ;
      return corrected_cvs;
}

  
  
  
  
};

#endif
