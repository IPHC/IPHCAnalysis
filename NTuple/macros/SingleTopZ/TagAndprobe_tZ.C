#include <iomanip>
#include <iostream>
#include "NTFormat/interface/NTEvent.h"

//NTupleAnalysis classes
#include "Selection/interface/DiLeptonSelection.h"
#include "Tools/interface/Dataset.h"
#include "Tools/interface/AnalysisEnvironmentLoader.h"
#include "Plots/interface/DiLepAnaHistoManager.h"
#include "BckgdEstimation/interface/DYEstimation.h"
#include "BckgdEstimation/interface/MMEstimation.h"
#include "BckgdEstimation/interface/FakeRate_LeptEff.h"
#include "Tools/interface/TagAndProbe.h"
#include "Tools/interface/PUWeighting.h"
#include "Tools/interface/LumiReweightingStandAlone.h"
#include "Tools/interface/JetCorrector.h"


#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

using namespace IPHCTree;
using namespace std;

int main ()
{
  cout<<"#########################"<<endl;
  cout<<"Beginning of the program"<<endl;
  cout<<"#########################"<<endl;
  
  //////////////////////
  //Global variables
  //////////////////////
  vector < Dataset > datasets;
  DiLeptonSelection sel; 
  float Luminosity = 0;
  float LumiError = 0;
  // 0: MC - 1: Data - 2 Data & MC
  int DataType = 0; 
  int verbosity = -1;

  //////////////////////

  string PUWeightFileName;
  //////////////////////
  // Initialisation
  //////////////////////
  string xmlFileName = string ("../../config/TagAndProb_2011AB_FCNCkut.xml");
  AnalysisEnvironmentLoader anaEL (xmlFileName);
  anaEL.LoadSamples (datasets); // now the list of datasets written in the xml file is known
  anaEL.LoadDiLeptonSelection (sel); // now the parameters for the selection are given to the selection
  anaEL.LoadGeneralInfo(DataType, Luminosity, LumiError, PUWeightFileName, verbosity );
  
    
  IPHCTree::NTEvent * event = 0;
  //Selection table
  // 4 tables: ee - emu - mumu - allChannels
  
   
  //CutFlow tables
  double TabFlow1[5][101][101]; //central values [channel][typemc][cuttype]
  double TabFlow2[5][101][101]; // errors
  //dataset weights
  double Dweight[101]; 
  for(int k0=0; k0<5; ++k0) {
    for(int k1=0; k1<101; ++k1) {
      Dweight[k1] = 0.;
      for(int k2=0; k2<101; ++k2) {
        TabFlow1[k0][k1][k2]=0.;
        TabFlow2[k0][k1][k2]=0.;
      }
    }
  } 
  vector<string> CutName;
  vector<string> VecChannelName;
  
  // Here define the studied channel (ee/mumu/emu)
  string ChannelName  = "";
  //
  
  bool IReweightWithVxP = true;
  bool IReweight        = true;

  
  
  //////////////////////
 
  //////////////////////
  
  
  
    
  cout << "loading datasets " << endl;
  //Book keeping of my histos
  HistoManager MyhistoManager; 
  MyhistoManager.LoadDatasets(datasets);    
  cout << "datasets loaded " << endl;
  
  
  //modifdg
  std::vector<TH1F> MyHistos;
  std::vector<TH1F> MyHistos1;
  std::vector<TH1F> HInvM_ll_pair;
  std::vector<TH1F> HInvM_emu_pair;
  std::vector<TH1F> HInvM_emu_pair_aftermetcut;
  std::vector<TH1F> HNjets_aftermetcut;
  std::vector<TH1F> HMet;
  std::vector<TH1F> HMet_aftermetcut;
  std::vector<TH1F> HTriggerEff_afterDileptSel_pT;
  std::vector<TH1F> HTriggerEff_afterDileptSel_eta;
  std::vector<TH1F> HTriggerEff_afterDileptSel_phi;
  std::vector<TH1F> HTriggerEff_afterDileptSel_selTrig_pT;
  std::vector<TH1F> HTriggerEff_afterDileptSel_selTrig_eta;
  std::vector<TH1F> HTriggerEff_afterDileptSel_selTrig_phi;
  
  std::vector<TH1F> HNjets;
  std::vector<TH1F> HNjets_Loose;
  std::vector<TH1F> HNjets_Medium;
  
  std::vector<TH1F> HNjets_afterLeptSel;
  std::vector<TH1F> HNjets_Loose_afterLeptSel;
  std::vector<TH1F> HNjets_Medium_afterLeptSel;
  
  std::vector<TH1F> deltaR_Jet_Lepton_FullSel;
  std::vector<TH1F> deltaR_Jet_MET_FullSel;
  std::vector<TH1F> deltaR_Muon_MET_FullSel;
  std::vector<TH1F> deltaR_Elec_MET_FullSel;
  std::vector<TH1F> deltaR_Muon_Electron_FullSel;
  
  std::vector<TH1F> deltaPhi_Jet_Lepton_FullSel;
  std::vector<TH1F> deltaPhi_Jet_MET_FullSel;
  std::vector<TH1F> deltaPhi_Muon_MET_FullSel;
  std::vector<TH1F> deltaPhi_Elec_MET_FullSel;
  std::vector<TH1F> deltaPhi_DiLept_MET_FullSel;
  


  std::vector<TH1F> HNmuons_ElecSel;
  std::vector<TH1F> HNelec_ElecSel;
  std::vector<TH1F> HNmuons_MuonSel;
  std::vector<TH1F> HNelec_MuonSel;
  std::vector<TH1F> HNmuons_DiMuonElecSel;
  std::vector<TH1F> HNelec_DiMuonElecSel;
  std::vector<TH1F> HNmuons_DiElecMuonSel;
  std::vector<TH1F> HNelec_DiElecMuonSel;

  
  std::vector<TH1F> HTriggerEff_ElecSel_pT;
  std::vector<TH1F> HTriggerEff_ElecSel_eta;
  std::vector<TH1F> HTriggerEff_ElecSel_phi;
  std::vector<TH1F> HTriggerEff_ElecSel_selTrig_pT;
  std::vector<TH1F> HTriggerEff_ElecSel_selTrig_eta;
  std::vector<TH1F> HTriggerEff_ElecSel_selTrig_phi; 
  
  std::vector<TH1F> HTriggerEff_MuonSel_pT;
  std::vector<TH1F> HTriggerEff_MuonSel_eta;
  std::vector<TH1F> HTriggerEff_MuonSel_phi;
  std::vector<TH1F> HTriggerEff_MuonSel_selTrig_pT;
  std::vector<TH1F> HTriggerEff_MuonSel_selTrig_eta;
  std::vector<TH1F> HTriggerEff_MuonSel_selTrig_phi;
  
  std::vector<TH1F> HTriggerEff_DiMuonElecSel_pT;
  std::vector<TH1F> HTriggerEff_DiMuonElecSel_eta;
  std::vector<TH1F> HTriggerEff_DiMuonElecSel_phi;
  std::vector<TH1F> HTriggerEff_DiMuonElecSel_selTrig_pT;
  std::vector<TH1F> HTriggerEff_DiMuonElecSel_selTrig_eta;
  std::vector<TH1F> HTriggerEff_DiMuonElecSel_selTrig_phi;
  
  std::vector<TH1F> HTriggerEff_DiElecMuonSel_pT;
  std::vector<TH1F> HTriggerEff_DiElecMuonSel_eta;
  std::vector<TH1F> HTriggerEff_DiElecMuonSel_phi;
  std::vector<TH1F> HTriggerEff_DiElecMuonSel_selTrig_pT;
  std::vector<TH1F> HTriggerEff_DiElecMuonSel_selTrig_eta;
  std::vector<TH1F> HTriggerEff_DiElecMuonSel_selTrig_phi;
  
  
  std::vector<TH2D> correlationBTAGMET_mumumu;
  std::vector<TH2D> correlationMET_mumumu;
  
  std::vector<TH2D> correlationBTAGMET_eee;
  std::vector<TH2D> correlationMET_eee;
  
  std::vector<TH2D> correlationBTAGMET_mumue;
  std::vector<TH2D> correlationMET_mumue;
  
  std::vector<TH2D> correlationBTAGMET_eemu;
  std::vector<TH2D> correlationMET_eemu;
  
  //*************************
  //*************************
  //*************************
  // for MC closure test
  
  
  std::vector<TH1F> genEl_Tight_pt;
  std::vector<TH1F> genEl_Tight_eta;
  std::vector<TH1F> genEl_Tight_phi;
  std::vector<TH1F> genEl_Tight_njet;
  
  std::vector<TH1F> genEl_LooseID_pt;
  std::vector<TH1F> genEl_LooseID_eta;
  std::vector<TH1F> genEl_LooseID_phi;
  std::vector<TH1F> genEl_LooseID_njet;
  
  std::vector<TH1F> genEl_LooseIso_pt;
  std::vector<TH1F> genEl_LooseIso_eta;
  std::vector<TH1F> genEl_LooseIso_phi;
  std::vector<TH1F> genEl_LooseIso_njet;
  
  
  std::vector<TH1F> genMu_Tight_pt;
  std::vector<TH1F> genMu_Tight_eta;
  std::vector<TH1F> genMu_Tight_phi;
  std::vector<TH1F> genMu_Tight_njet;
  
  std::vector<TH1F> genMu_LooseID_pt;
  std::vector<TH1F> genMu_LooseID_eta;
  std::vector<TH1F> genMu_LooseID_phi;
  std::vector<TH1F> genMu_LooseID_njet;
  
  std::vector<TH1F> genMu_LooseIso_pt;
  std::vector<TH1F> genMu_LooseIso_eta;
  std::vector<TH1F> genMu_LooseIso_phi;
  std::vector<TH1F> genMu_LooseIso_njet;
  
  //endmodifdg
  
  
  int ITypeMC	  = -1;
  int ICut	  = -1;  
  int IChannel    = -1;  
 
    
  cout<<"The verbosity mode is "<<verbosity <<endl;
  cout<<"The luminosity is equal to "<< Luminosity<<endl;

  
  JetCorrector JEC_L2L3Residuals;
  JEC_L2L3Residuals.LoadCorrections();

  
  //////////////////////
  //LOOP OVER THE DATASETS
  //////////////////////
  if(verbosity>0) {
  	cout<<"#########################"<<endl;
  	cout<<" Loop over the datasets  "<<endl;
	cout<<"#########################"<<endl;
  }

  for (unsigned int d = 0; d < datasets.size (); d++) {
    

    datasets[d].eventTree ()->SetBranchAddress ("NTEvent",&event);

    unsigned int nEvents = (int) (datasets[d].eventTree ()->GetEntries ());
    cout << "NEvents = " << nEvents << endl;
    
   
    
    TagAndProbe theEffPlots;
    theEffPlots.CreateHistos("muons",     datasets[d].Name());
    theEffPlots.CreateHistos("electrons", datasets[d].Name());
    
    MyhistoManager.CreateHisto(HTriggerEff_afterDileptSel_pT, "HTriggerEff_afterDileptSel_pT" ,datasets[d].Name(),"trig eff Pt"  ,"Entries", 60, 20, 300);
    MyhistoManager.CreateHisto(HTriggerEff_afterDileptSel_eta,"HTriggerEff_afterDileptSel_eta" ,datasets[d].Name(),"trig eff eta","Entries", 21,-2.5,2.5);   
    MyhistoManager.CreateHisto(HTriggerEff_afterDileptSel_phi,"HTriggerEff_afterDileptSel_phi" ,datasets[d].Name(),"trig eff phi","Entries", 20, -3.2, 3.2);
    
    
    MyhistoManager.CreateHisto(HTriggerEff_afterDileptSel_selTrig_pT, "HTriggerEff_afterDileptSel_selTrig_pT" ,datasets[d].Name(),"trig eff Pt"  ,"Entries", 60, 20, 300);
    MyhistoManager.CreateHisto(HTriggerEff_afterDileptSel_selTrig_eta,"HTriggerEff_afterDileptSel_selTrig_eta" ,datasets[d].Name(),"trig eff eta","Entries", 21,-2.5,2.5);
    MyhistoManager.CreateHisto(HTriggerEff_afterDileptSel_selTrig_phi,"HTriggerEff_afterDileptSel_selTrig_phi" ,datasets[d].Name(),"trig eff phi","Entries", 20, -3.2, 3.2);
    
    
    MyhistoManager.CreateHisto(HNjets,        "Njets"        ,datasets[d].Name(),"Njets","Entries",4,-0.5,3.5);
    MyhistoManager.CreateHisto(HNjets_Loose,  "Njets_Loose"  ,datasets[d].Name(),"Njets_Loose","Entries",4,-0.5,3.5);
    MyhistoManager.CreateHisto(HNjets_Medium, "Njets_Medium" ,datasets[d].Name(),"Njets_Medium","Entries",4,-0.5,3.5);
    
    MyhistoManager.CreateHisto(HNjets_afterLeptSel,        "Njets_afterLeptSel"        ,datasets[d].Name(),"Njets_afterLeptSel"       ,"Entries",4,-0.5,3.5);
    MyhistoManager.CreateHisto(HNjets_Loose_afterLeptSel,  "Njets_Loose_afterLeptSel"  ,datasets[d].Name(),"Njets_Loose_afterLeptSel" ,"Entries",4,-0.5,3.5);
    MyhistoManager.CreateHisto(HNjets_Medium_afterLeptSel, "Njets_Medium_afterLeptSel" ,datasets[d].Name(),"Njets_Medium_afterLeptSel","Entries",4,-0.5,3.5);
    
    MyhistoManager.CreateHisto(deltaR_Jet_Lepton_FullSel   , "deltaR_Jet_Lepton_FullSel"   , datasets[d].Name(), "delta R", "Entries", 20, 0, 4 );         
    MyhistoManager.CreateHisto(deltaR_Jet_MET_FullSel      , "deltaR_Jet_MET_FullSel"      , datasets[d].Name(), "delta R", "Entries", 20, 0, 4 );         
    MyhistoManager.CreateHisto(deltaR_Muon_MET_FullSel     , "deltaR_Muon_MET_FullSel"     , datasets[d].Name(), "delta R", "Entries", 20, 0, 4 );           
    MyhistoManager.CreateHisto(deltaR_Elec_MET_FullSel     , "deltaR_Elec_MET_FullSel"     , datasets[d].Name(), "delta R", "Entries", 20, 0, 4 );       
    //MyhistoManager.CreateHisto(deltaR_Muon_Electron_FullSel, "deltaR_Muon_Electron_FullSel", datasets[d].Name(), "delta R", "Entries", 20, 0, 4 );         
    MyhistoManager.CreateHisto(deltaPhi_Jet_Lepton_FullSel   , "deltaPhi_Jet_Lepton_FullSel"   , datasets[d].Name(), "delta R", "Entries", 20, 0, 4 );         
    MyhistoManager.CreateHisto(deltaPhi_Jet_MET_FullSel      , "deltaPhi_Jet_MET_FullSel"      , datasets[d].Name(), "delta R", "Entries", 20, 0, 4 );         
    MyhistoManager.CreateHisto(deltaPhi_Muon_MET_FullSel     , "deltaPhi_Muon_MET_FullSel"     , datasets[d].Name(), "delta R", "Entries", 20, 0, 4 );           
    MyhistoManager.CreateHisto(deltaPhi_Elec_MET_FullSel     , "deltaPhi_Elec_MET_FullSel"     , datasets[d].Name(), "delta R", "Entries", 20, 0, 4 );       
              
    MyhistoManager.CreateHisto(deltaPhi_DiLept_MET_FullSel     , "deltaPhi_DiLept_MET_FullSel"     , datasets[d].Name(), "delta R", "Entries", 20, 0, 4 ); 
    
    
    
    MyhistoManager.CreateHisto(HNmuons_ElecSel,  "Nmuons_ElecSel"  ,datasets[d].Name(),"Nmuons","Entries",5,-0.5,4.5);
    MyhistoManager.CreateHisto(HNelec_ElecSel,   "Nelec_ElecSel"   ,datasets[d].Name(),"Nelec" ,"Entries",5,-0.5,4.5);
    MyhistoManager.CreateHisto(HNmuons_MuonSel,  "Nmuons_MuonSel"  ,datasets[d].Name(),"Nmuons","Entries",5,-0.5,4.5);
    MyhistoManager.CreateHisto(HNelec_MuonSel,   "Nelec_MuonSel"   ,datasets[d].Name(),"Nelec" ,"Entries",5,-0.5,4.5);
    MyhistoManager.CreateHisto(HNmuons_DiMuonElecSel,  "Nmuons_DiMuonElecSel"  ,datasets[d].Name(),"Nmuons","Entries",5,-0.5,4.5);
    MyhistoManager.CreateHisto(HNelec_DiMuonElecSel,   "Nelec_DiMuonElecSel"   ,datasets[d].Name(),"Nelec" ,"Entries",5,-0.5,4.5);
    MyhistoManager.CreateHisto(HNmuons_DiElecMuonSel,  "Nmuons_DiElecMuonSel"  ,datasets[d].Name(),"Nmuons","Entries",5,-0.5,4.5);
    MyhistoManager.CreateHisto(HNelec_DiElecMuonSel,   "Nelec_DiElecMuonSel"   ,datasets[d].Name(),"Nelec" ,"Entries",5,-0.5,4.5);
    
    
    MyhistoManager.CreateHisto(HTriggerEff_ElecSel_pT,  "HTriggerEff_ElecSel_pT",  datasets[d].Name(),"trig eff el Pt"  ,"Entries", 60, 20, 300);
    MyhistoManager.CreateHisto(HTriggerEff_ElecSel_eta, "HTriggerEff_ElecSel_eta", datasets[d].Name(),"trig eff el eta" ,"Entries", 21,-2.5,2.5);
    MyhistoManager.CreateHisto(HTriggerEff_ElecSel_phi, "HTriggerEff_ElecSel_phi", datasets[d].Name(),"trig eff el phi" ,"Entries", 20, -3.2, 3.2);
    MyhistoManager.CreateHisto(HTriggerEff_ElecSel_selTrig_pT,  "HTriggerEff_ElecSel_selTrig_pT",  datasets[d].Name(),"trig eff sel el Pt"  ,"Entries", 60, 20, 300);
    MyhistoManager.CreateHisto(HTriggerEff_ElecSel_selTrig_eta, "HTriggerEff_ElecSel_selTrig_eta", datasets[d].Name(),"trig eff sel el eta" ,"Entries", 21,-2.5,2.5);
    MyhistoManager.CreateHisto(HTriggerEff_ElecSel_selTrig_phi, "HTriggerEff_ElecSel_selTrig_phi", datasets[d].Name(),"trig eff sel el phi" ,"Entries", 20, -3.2, 3.2);
  
    MyhistoManager.CreateHisto(HTriggerEff_MuonSel_pT,  "HTriggerEff_MuonSel_pT",  datasets[d].Name(),"trig eff el Pt"  ,"Entries", 60, 20, 300);
    MyhistoManager.CreateHisto(HTriggerEff_MuonSel_eta, "HTriggerEff_MuonSel_eta", datasets[d].Name(),"trig eff el eta" ,"Entries", 21,-2.5,2.5);
    MyhistoManager.CreateHisto(HTriggerEff_MuonSel_phi, "HTriggerEff_MuonSel_phi", datasets[d].Name(),"trig eff el phi" ,"Entries", 20, -3.2, 3.2);
    MyhistoManager.CreateHisto(HTriggerEff_MuonSel_selTrig_pT,  "HTriggerEff_MuonSel_selTrig_pT",  datasets[d].Name(),"trig eff sel el Pt"  ,"Entries", 60, 20, 300);
    MyhistoManager.CreateHisto(HTriggerEff_MuonSel_selTrig_eta, "HTriggerEff_MuonSel_selTrig_eta", datasets[d].Name(),"trig eff sel el eta" ,"Entries", 21,-2.5,2.5);
    MyhistoManager.CreateHisto(HTriggerEff_MuonSel_selTrig_phi, "HTriggerEff_MuonSel_selTrig_phi", datasets[d].Name(),"trig eff sel el phi" ,"Entries", 20, -3.2, 3.2);
  
    MyhistoManager.CreateHisto(HTriggerEff_DiMuonElecSel_pT,  "HTriggerEff_DiMuonElecSel_pT",  datasets[d].Name(),"trig eff el Pt"  ,"Entries", 60, 20, 300);
    MyhistoManager.CreateHisto(HTriggerEff_DiMuonElecSel_eta, "HTriggerEff_DiMuonElecSel_eta", datasets[d].Name(),"trig eff el eta" ,"Entries", 21,-2.5,2.5);
    MyhistoManager.CreateHisto(HTriggerEff_DiMuonElecSel_phi, "HTriggerEff_DiMuonElecSel_phi", datasets[d].Name(),"trig eff el phi" ,"Entries", 20, -3.2, 3.2);
    MyhistoManager.CreateHisto(HTriggerEff_DiMuonElecSel_selTrig_pT,  "HTriggerEff_DiMuonElecSel_selTrig_pT",  datasets[d].Name(),"trig eff sel el Pt"  ,"Entries", 60, 20, 300);
    MyhistoManager.CreateHisto(HTriggerEff_DiMuonElecSel_selTrig_eta, "HTriggerEff_DiMuonElecSel_selTrig_eta", datasets[d].Name(),"trig eff sel el eta" ,"Entries", 21,-2.5,2.5);
    MyhistoManager.CreateHisto(HTriggerEff_DiMuonElecSel_selTrig_phi, "HTriggerEff_DiMuonElecSel_selTrig_phi", datasets[d].Name(),"trig eff sel el phi" ,"Entries", 20, -3.2, 3.2);
  
    MyhistoManager.CreateHisto(HTriggerEff_DiElecMuonSel_pT,  "HTriggerEff_DiElecMuonSel_pT",  datasets[d].Name(),"trig eff el Pt"  ,"Entries", 60, 20, 300);
    MyhistoManager.CreateHisto(HTriggerEff_DiElecMuonSel_eta, "HTriggerEff_DiElecMuonSel_eta", datasets[d].Name(),"trig eff el eta" ,"Entries", 21,-2.5,2.5);
    MyhistoManager.CreateHisto(HTriggerEff_DiElecMuonSel_phi, "HTriggerEff_DiElecMuonSel_phi", datasets[d].Name(),"trig eff el phi" ,"Entries", 20, -3.2, 3.2);
    MyhistoManager.CreateHisto(HTriggerEff_DiElecMuonSel_selTrig_pT,  "HTriggerEff_DiElecMuonSel_selTrig_pT",  datasets[d].Name(),"trig eff sel el Pt"  ,"Entries", 60, 20, 300);
    MyhistoManager.CreateHisto(HTriggerEff_DiElecMuonSel_selTrig_eta, "HTriggerEff_DiElecMuonSel_selTrig_eta", datasets[d].Name(),"trig eff sel el eta" ,"Entries", 21,-2.5,2.5);
    MyhistoManager.CreateHisto(HTriggerEff_DiElecMuonSel_selTrig_phi, "HTriggerEff_DiElecMuonSel_selTrig_phi", datasets[d].Name(),"trig eff sel el phi" ,"Entries", 20, -3.2, 3.2);

    //float ptRange[] = {20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 65, 80, 300};
    float ptRange[] = {20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 100};
   
    MyhistoManager.CreateHisto(genEl_Tight_pt,   "genEl_Tight_pt",   datasets[d].Name(),"eff gen el "  ,"Entries", 10, ptRange);
    MyhistoManager.CreateHisto(genEl_Tight_eta,  "genEl_Tight_eta",  datasets[d].Name(),"eff gen el "  ,"Entries", 40, -2.5 , 2.5);
    MyhistoManager.CreateHisto(genEl_Tight_phi,  "genEl_Tight_phi",  datasets[d].Name(),"eff gen el el "  ,"Entries", 40, -3.14, 3.14);
    MyhistoManager.CreateHisto(genEl_Tight_njet, "genEl_Tight_njet", datasets[d].Name(),"eff gen "  ,"Entries", 7,  -0.5 , 6.5);
  
    MyhistoManager.CreateHisto(genEl_LooseID_pt,   "genEl_LooseID_pt",   datasets[d].Name(),"eff gen el "  ,"Entries", 10, ptRange);
    MyhistoManager.CreateHisto(genEl_LooseID_eta,  "genEl_LooseID_eta",  datasets[d].Name(),"eff gen el "  ,"Entries", 40, -2.5 , 2.5);
    MyhistoManager.CreateHisto(genEl_LooseID_phi,  "genEl_LooseID_phi",  datasets[d].Name(),"eff gen el "  ,"Entries", 40, -3.14, 3.14);
    MyhistoManager.CreateHisto(genEl_LooseID_njet, "genEl_LooseID_njet", datasets[d].Name(),"eff gen el "  ,"Entries", 7,  -0.5 , 6.5);
  
    MyhistoManager.CreateHisto(genEl_LooseIso_pt,   "genEl_LooseIso_pt",   datasets[d].Name(),"eff gen el "  ,"Entries", 10, ptRange);
    MyhistoManager.CreateHisto(genEl_LooseIso_eta,  "genEl_LooseIso_eta",  datasets[d].Name(),"eff gen el "  ,"Entries", 40, -2.5 , 2.5);
    MyhistoManager.CreateHisto(genEl_LooseIso_phi,  "genEl_LooseIso_phi",  datasets[d].Name(),"eff gen el "  ,"Entries", 40, -3.14, 3.14);
    MyhistoManager.CreateHisto(genEl_LooseIso_njet, "genEl_LooseIso_njet", datasets[d].Name(),"eff gen el "  ,"Entries", 7,  -0.5 , 6.5);
  
  
  
    MyhistoManager.CreateHisto(genMu_Tight_pt,   "genMu_Tight_pt",   datasets[d].Name(),"eff gen mu "  ,"Entries", 10, ptRange);
    MyhistoManager.CreateHisto(genMu_Tight_eta,  "genMu_Tight_eta",  datasets[d].Name(),"eff gen mu "  ,"Entries", 40, -2.5 , 2.5);
    MyhistoManager.CreateHisto(genMu_Tight_phi,  "genMu_Tight_phi",  datasets[d].Name(),"eff gen mu "  ,"Entries", 40, -3.14, 3.14);
    MyhistoManager.CreateHisto(genMu_Tight_njet, "genMu_Tight_njet", datasets[d].Name(),"eff gen mu "  ,"Entries", 7,  -0.5 , 6.5);
  
    MyhistoManager.CreateHisto(genMu_LooseID_pt,   "genMu_LooseID_pt",   datasets[d].Name(),"eff gen mu "  ,"Entries", 10, ptRange);
    MyhistoManager.CreateHisto(genMu_LooseID_eta,  "genMu_LooseID_eta",  datasets[d].Name(),"eff gen mu "  ,"Entries", 40, -2.5 , 2.5);
    MyhistoManager.CreateHisto(genMu_LooseID_phi,  "genMu_LooseID_phi",  datasets[d].Name(),"eff gen mu "  ,"Entries", 40, -3.14, 3.14);
    MyhistoManager.CreateHisto(genMu_LooseID_njet, "genMu_LooseID_njet", datasets[d].Name(),"eff gen mu "  ,"Entries", 7,  -0.5 , 6.5);
  
    MyhistoManager.CreateHisto(genMu_LooseIso_pt,   "genMu_LooseIso_pt",   datasets[d].Name(),"eff gen mu "  ,"Entries", 10, ptRange);
    MyhistoManager.CreateHisto(genMu_LooseIso_eta,  "genMu_LooseIso_eta",  datasets[d].Name(),"eff gen mu "  ,"Entries", 40, -2.5 , 2.5);
    MyhistoManager.CreateHisto(genMu_LooseIso_phi,  "genMu_LooseIso_phi",  datasets[d].Name(),"eff gen mu "  ,"Entries", 40, -3.14, 3.14);
    MyhistoManager.CreateHisto(genMu_LooseIso_njet, "genMu_LooseIso_njet", datasets[d].Name(),"eff gen mu "  ,"Entries", 7,  -0.5 , 6.5);
  
  
 
    MyhistoManager.CreateHisto2D(correlationBTAGMET_mumumu, "correlationBTAGMET_mumumu", datasets[d].Name(),"NTrigBTAGMET "  , 2, -0.5 , 1.5, "NTrigLept", 2, -0.5, 1.5);
    MyhistoManager.CreateHisto2D(correlationMET_mumumu, "correlationMET_mumumu", datasets[d].Name(),"NTrigMET "   , 2, -0.5 , 1.5, "NTrigLept", 2, -0.5, 1.5);
  
    MyhistoManager.CreateHisto2D(correlationBTAGMET_eee, "correlationBTAGMET_eee", datasets[d].Name(),"NTrigBTAGMET "  , 2, -0.5 , 1.5, "NTrigLept", 2, -0.5, 1.5);
    MyhistoManager.CreateHisto2D(correlationMET_eee, "correlationMET_eee", datasets[d].Name(),"NTrigMET "   , 2, -0.5 , 1.5, "NTrigLept", 2, -0.5, 1.5);
  
    MyhistoManager.CreateHisto2D(correlationBTAGMET_mumue, "correlationBTAGMET_mumue", datasets[d].Name(),"NTrigBTAGMET "  , 2, -0.5 , 1.5, "NTrigLept", 2, -0.5, 1.5);
    MyhistoManager.CreateHisto2D(correlationMET_mumue, "correlationMET_mumue", datasets[d].Name(),"NTrigMET "   , 2, -0.5 , 1.5, "NTrigLept", 2, -0.5, 1.5);
  
    MyhistoManager.CreateHisto2D(correlationBTAGMET_eemu, "correlationBTAGMET_eemu", datasets[d].Name(),"NTrigBTAGMET "  , 2, -0.5 , 1.5, "NTrigLept", 2, -0.5, 1.5);
    MyhistoManager.CreateHisto2D(correlationMET_eemu, "correlationMET_eemu", datasets[d].Name(),"NTrigMET "   , 2, -0.5 , 1.5, "NTrigLept", 2, -0.5, 1.5);
  
 
    //////////////////////
    // PU Reweighting
    //////////////////////

    reweight::LumiReWeighting *LumiWeights;
    
    if (IReweight ) {

      string mcfile;
      if( datasets[d].Name() == "FCNCkut" ) // FastSim, in-time PU only
	mcfile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PUMC_InTime_Fall11.root");
      else
	mcfile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PU3DMC_Fall11_JLA.root");
      fexists(mcfile, true);

      string datafile;
      if( datasets[d].Name() == "FCNCkut" )
	datafile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PUData2011_observed_68mb.root");
      else
	datafile = getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/PUData2011_68mb.root");
      fexists(datafile, true);

      LumiWeights = new reweight::LumiReWeighting(mcfile, datafile, "histoMCPU", "pileup" );
      LumiWeights->weight3D_init( 1. );

    }


    //LOOP OVER THE EVENTS
    for (unsigned int ievt = 0; ievt < nEvents; ievt++) { //nEvents
      float weight = 1.;
      datasets[d].eventTree ()->GetEntry (ievt);
      IPHCTree::NTTransient::InitializeAfterReading(event);
      if (verbosity > 3)
	//cout << "run number=" << event->runNb << endl;    
     cout << "event number=" << event->general.eventNb << endl;
     if(ievt%10000 == 0) cout << "number of processed events " << ievt << endl;
     //Load event for the selection
     //if(ievt == 3000000) break;
     
     
     bool isData = false;    
     if(datasets[d].Name()=="DataEG" || datasets[d].Name()=="DataMu" || 
	datasets[d].Name()=="DataMuEG" || datasets[d].Name()=="DataEGMu" ||
	datasets[d].Name()=="MET1" || datasets[d].Name()=="MET2") isData = true;
     if(isData) JEC_L2L3Residuals.ApplyCorrections(event); // n'appliquer la correction que pour les donnees


     //********************************
     //**** check trigger paths
     //********************************
     
      //Dweight[ITypeMC] = Luminosity*datasets[d].Xsection()/datasets[d].NofEvtsToRunOver();
      //Dweight[ITypeMC] = Luminosity*datasets[d].Xsection()/datasets[d].nEvents;
      //if ( datasets[d].Name()=="Data" || datasets[d].Name()=="DATA") Dweight[ITypeMC] = 1;
     
     
     
    
     //cout << " pass trigger TriggerPassed_mumu_DATA " << TriggerPassed_mumu_DATA << endl;
     //cout << " pass trigger TriggerPassed_ee_DATA " << TriggerPassed_ee_DATA << endl;
     //std::ostream & ostrigger:
     //if(ievt == 1) event->PrintTriggerList (std::cout);
     //if(ievt == 1) event->  PrintTriggerPassed(std::cout);
      
      bool IsSignal = false;
      int  LastStep = 0;
     
     
      //if(!IsData) cout  << TriggerPassed_ee_MC << "  sel.eleHLTMatch " << sel.eleHLTMatch << endl;
    
     sel.LoadEvent(event);

     //Collection of selected objects
     vector<NTElectron> selElectrons = sel.GetSelectedElectrons();
     vector<NTMuon>     selMuons = sel.GetSelectedMuons();     
     vector<NTVertex>   selVertices  = sel.GetSelectedVertex();
     vector<NTJet>      selJets = sel.GetSelectedJets();
     NTMET met = sel.GetMET(); // no criteria applyied
     
     // Set collections for TagAndProbe
     vector<NTMuon>      selMuonsNonIso = sel.GetSelectedMuonsNoIso();
     vector<NTElectron>  selElecsNonIso = sel.GetSelectedElectronsNoIso();
     theEffPlots.SetLooseIsoMuonCollection(     selMuonsNonIso );
     theEffPlots.SetLooseIsoElectronCollection( selElecsNonIso );
     
     vector<NTMuon>      selMuonsNonID  = sel.GetSelectedMuonIsoNonID();
     vector<NTElectron>  selElecsNonID  = sel.GetSelectedElectronsIsoNonID();
     theEffPlots.SetLooseIDMuonCollection(      selMuonsNonID  );
     theEffPlots.SetLooseIDElectronCollection(  selElecsNonID  );

     theEffPlots.SetTightMuonCollection( selMuons );
     theEffPlots.SetTightElectronCollection( selElectrons );
       
     
     
     
     
     //Manage DY samples to avoid overlaps
      double dileptInvMass = 0;
      if( (event->mc.zAndDecays).size() > 0){
        TLorentzVector dilept = (event->mc.zAndDecays)[0].p4_Lep1_gen + (event->mc.zAndDecays)[0].p4_Lep2_gen;
	dileptInvMass = dilept.M();
      }
      
      /*if(datasets[d].Name()=="Zjets" && dileptInvMass < 50) continue;      
      if(datasets[d].Name()=="DYToMuMu_M-20"	   && (dileptInvMass > 50 || dileptInvMass < 20) ) continue;
      if(datasets[d].Name()=="DYToEE_M-20"	   && (dileptInvMass > 50 || dileptInvMass < 20) ) continue;
      if(datasets[d].Name()=="DYToTauTau_M-20"	   && (dileptInvMass > 50 || dileptInvMass < 20) ) continue;
      if(datasets[d].Name()=="DYToMuMu_M-10To20"   &&  dileptInvMass > 20) continue;
      if(datasets[d].Name()=="DYToEE_M-10To20"	   &&  dileptInvMass > 20) continue;
      if(datasets[d].Name()=="DYToTauTau_M-10To20" &&  dileptInvMass > 20) continue;*/
      
      
      
     
      double weightITypeMC_save = Luminosity*datasets[d].Xsection()/datasets[d].getNSkimmedEvent();
      double weightITypeMC=0;
 
     
     
     
     for (int IChannel=0; IChannel<3; IChannel++) {
       //      for (int IChannel=1; IChannel<2; IChannel++) {
       string ChannelName;
       if (IChannel==0) ChannelName= "mumu"; 
       else if (IChannel==1) ChannelName= "ee"; 
       else if (IChannel==2) ChannelName= "emu"; 
       
       
       if (IChannel==0 && (datasets[d].Name()=="DataEG" || datasets[d].Name()=="DataMuEG")) continue;
          if (IChannel==1 && (datasets[d].Name()=="DataMu" || datasets[d].Name()=="DataMuEG")) continue;
          if (IChannel==2 && (datasets[d].Name()=="DataMu" || datasets[d].Name()=="DataEG")) continue;


       //*****************************************************************
       // calcul MC weights
       //*****************************************************************    
       if (!isData) {

	 if(IReweight ){

           if( datasets[d].Name() == "FCNCkut" ) // FastSim, in-time PU only
             weightITypeMC = weightITypeMC_save*LumiWeights->ITweight(event->pileup.intime_npu);
	   else
	     weightITypeMC = weightITypeMC_save*LumiWeights->weight3D(event->pileup.before_npu, event->pileup.intime_npu, event->pileup.after_npu);

	 }
	 else weightITypeMC = weightITypeMC_save;
       }
       else weightITypeMC = 1;


      
       if ( datasets[d].Name()=="TTbar" ) {
            if ( IChannel==0) { // "mumu" 
	      if ( event->mc.TMEME==20 || event->mc.TMEME==11010 || event->mc.TMEME==22000 )    IsSignal = true;
	      if ( !(event->mc.TMEME==20 || event->mc.TMEME==11010 || event->mc.TMEME==22000) ) IsSignal = false;
	    }      
            else if ( IChannel==1) {  // "ee" 
	      if ( event->mc.TMEME==2 || event->mc.TMEME==10101 || event->mc.TMEME==20200 )     IsSignal = true;
	      if ( !(event->mc.TMEME==2 || event->mc.TMEME==10101 || event->mc.TMEME==20200) )  IsSignal = false;
	    }      
            else if ( IChannel==2) { // "emu" 
	      if ( event->mc.TMEME==11 || event->mc.TMEME==21100 || event->mc.TMEME==11001 || event->mc.TMEME==10110 )     IsSignal = true;
	      if ( !(event->mc.TMEME==11 || event->mc.TMEME==21100 || event->mc.TMEME==11001 || event->mc.TMEME==10110) )  IsSignal = false;
	    }      
          }
	  
	  
	  
       bool passTrigger = sel.passTriggerSelection ( &datasets[d], ChannelName);
       
       if(passTrigger) theEffPlots.FillHistos(selJets, "muons",     datasets[d].Name(), weightITypeMC);
       if(passTrigger) theEffPlots.FillHistos(selJets, "electrons", datasets[d].Name(), weightITypeMC);
              
    
      
       //*************************************
       //get sel MC truth lept sel efficiency
       //for muons
       
       if( passTrigger) {
	 for(unsigned int imu = 0; imu<selMuons.size(); imu++){
	   if(selMuons[imu].LeptonOrigin != 10 && selMuons[imu].LeptonOrigin != 1 ) continue;
	   MyhistoManager.FillHisto(genMu_Tight_pt,   "genMu_Tight_pt",   selMuons[imu].p4.Pt() , datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(genMu_Tight_eta,  "genMu_Tight_eta",  selMuons[imu].p4.Eta(), datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(genMu_Tight_phi,  "genMu_Tight_phi",  selMuons[imu].p4.Phi(), datasets[d].Name(), IsSignal, weightITypeMC);
	   if(imu == 0) MyhistoManager.FillHisto(genMu_Tight_njet, "genMu_Tight_njet", selJets.size(), datasets[d].Name(), IsSignal, weightITypeMC);
	 }
	 for(unsigned int imu = 0; imu<selMuonsNonID.size(); imu++){
	   if(selMuonsNonID[imu].LeptonOrigin != 10 && selMuonsNonID[imu].LeptonOrigin != 1 ) continue;
	   MyhistoManager.FillHisto(genMu_LooseID_pt,   "genMu_LooseID_pt",  selMuonsNonID[imu].p4.Pt() , datasets[d].Name(), IsSignal, weightITypeMC);  
	   MyhistoManager.FillHisto(genMu_LooseID_eta,  "genMu_LooseID_eta", selMuonsNonID[imu].p4.Eta(), datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(genMu_LooseID_phi,  "genMu_LooseID_phi", selMuonsNonID[imu].p4.Phi(), datasets[d].Name(), IsSignal, weightITypeMC);
	   if(imu == 0) MyhistoManager.FillHisto(genMu_LooseID_njet, "genMu_LooseID_njet",selJets.size(), datasets[d].Name(), IsSignal, weightITypeMC);
	 }
	 
	 for(unsigned int imu = 0; imu<selMuonsNonIso.size(); imu++){
	   if(selMuonsNonIso[imu].LeptonOrigin != 10 && selMuonsNonIso[imu].LeptonOrigin != 1 ) continue;
	   MyhistoManager.FillHisto(genMu_LooseIso_pt,   "genMu_LooseIso_pt",  selMuonsNonIso[imu].p4.Pt() , datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(genMu_LooseIso_eta,  "genMu_LooseIso_eta", selMuonsNonIso[imu].p4.Eta(), datasets[d].Name(), IsSignal, weightITypeMC); 
	   MyhistoManager.FillHisto(genMu_LooseIso_phi,  "genMu_LooseIso_phi", selMuonsNonIso[imu].p4.Phi(), datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(genMu_LooseIso_njet, "genMu_LooseIso_njet",selJets.size(), datasets[d].Name(), IsSignal, 1);
	 }
       }
       
       
       
       //*************************************
       //get sel MC truth lept sel efficiency
       //for electrons
       if( passTrigger ){
	 for(unsigned int iel = 0; iel<selElectrons.size(); iel++){
	   if(selElectrons[iel].LeptonOrigin != 10 && selElectrons[iel].LeptonOrigin != 1 ) continue;
	   MyhistoManager.FillHisto(genEl_Tight_pt,   "genEl_Tight_pt",   selElectrons[iel].p4.Pt() , datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(genEl_Tight_eta,  "genEl_Tight_eta",  selElectrons[iel].p4.Eta(), datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(genEl_Tight_phi,  "genEl_Tight_phi",  selElectrons[iel].p4.Phi(), datasets[d].Name(), IsSignal, weightITypeMC);
	   if(iel == 0) MyhistoManager.FillHisto(genEl_Tight_njet, "genEl_Tight_njet", selJets.size(), datasets[d].Name(), IsSignal, weightITypeMC);
	 }
	 for(unsigned int iel = 0; iel<selElecsNonID.size(); iel++){
	   if(selElecsNonID[iel].LeptonOrigin != 10 && selElecsNonID[iel].LeptonOrigin != 1 ) continue;
	   MyhistoManager.FillHisto(genEl_LooseID_pt,   "genEl_LooseID_pt",  selElecsNonID[iel].p4.Pt() , datasets[d].Name(), IsSignal, weightITypeMC);  
	   MyhistoManager.FillHisto(genEl_LooseID_eta,  "genEl_LooseID_eta", selElecsNonID[iel].p4.Eta(), datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(genEl_LooseID_phi,  "genEl_LooseID_phi", selElecsNonID[iel].p4.Phi(), datasets[d].Name(), IsSignal, weightITypeMC);
	   if(iel == 0) MyhistoManager.FillHisto(genEl_LooseID_njet, "genEl_LooseID_njet",selJets.size(), datasets[d].Name(), IsSignal, weightITypeMC);
	 }
	 
	 for(unsigned int iel = 0; iel<selElecsNonIso.size(); iel++){
	   if(selElecsNonIso[iel].LeptonOrigin != 10 && selElecsNonIso[iel].LeptonOrigin != 1 ) continue;
	   MyhistoManager.FillHisto(genEl_LooseIso_pt,   "genEl_LooseIso_pt",  selElecsNonIso[iel].p4.Pt() , datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(genEl_LooseIso_eta,  "genEl_LooseIso_eta", selElecsNonIso[iel].p4.Eta(), datasets[d].Name(), IsSignal, weightITypeMC); 
	   MyhistoManager.FillHisto(genEl_LooseIso_phi,  "genEl_LooseIso_phi", selElecsNonIso[iel].p4.Phi(), datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(genEl_LooseIso_njet, "genEl_LooseIso_njet",selJets.size(), datasets[d].Name(), IsSignal, 1);
	 }
       }

       

       //Candidate pair of lepton
       string CandType; // ee - emu - mumu or false
       vector<NTElectron> candElec;
       vector<NTMuon> candMuon;
       //sel.GetLeptonPair(candMuon, candElec, CandType ); // fill the variables       
       // JLA : REPLACE BY 3 LEPTONS SELECTION AND CHANNELS
       
             
      //------------------------ 
      // Do 3 leptons selection 
      //------------------------
   
      vector<NTElectron> ZeeCand; 
      vector<NTMuon>     ZmumuCand; 

      vector<NTElectron> WeCand; 
      vector<NTMuon>     WmuCand; 
      

       //*****************************************************************
      // select Z->ee candidate
      //*****************************************************************    
     
      int leptonFlavor = 0;
      
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
	    //if( fabs(zeecand.M() - 200) < fabs(mInv-200) ){
	      theel1 = iel1;
	      theel2 = iel2;
	      mInv = zeecand.M();
	    }
	  }
	}
	
	if(theel1>=0 && theel2>=0){ //JLA
	  ZeeCand.push_back(selElectrons[theel1]);
	  ZeeCand.push_back(selElectrons[theel2]);
	  double invM = (ZeeCand[0].p4+ZeeCand[1].p4).M();
	  //cout << "lepton origin of Zee cand " << selElectrons[theel1].LeptonOrigin << "  " << selElectrons[theel2].LeptonOrigin << "  invmass " <<  invM << endl;
	}
      }
         
      //*****************************************************************
      // select W->enu candidate
      //*****************************************************************    
      
        for(unsigned int iel1 = 0; iel1 < selElectrons.size(); iel1++){
	  bool matchElec=false;
          for(unsigned int iel2 = 0; iel2 < ZeeCand.size(); iel2++){
             
	     if(fabs(selElectrons[iel1].p4.Pt() - ZeeCand[iel2].p4.Pt()) <  0.0001)  matchElec=true;
	   }
	  if(!matchElec){
	    WeCand.push_back(selElectrons[iel1]);
	    if(selElectrons[iel1].LeptonOrigin == 10) leptonFlavor = 1;
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
	    //if( fabs(zmumucand.M() - 200) < fabs(mInv-200) ){
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
      
        for(unsigned int imu1 = 0; imu1 < selMuons.size(); imu1++){
	  bool matchMuon = false;
          for(unsigned int imu2 = 0; imu2 < ZmumuCand.size(); imu2++){
             
	     if(fabs(selMuons[imu1].p4.Pt() - ZmumuCand[imu2].p4.Pt()) <  0.0001) matchMuon = true;
	     
          } 
	  if(!matchMuon){
	    WmuCand.push_back(selMuons[imu1]);
	    if(selMuons[imu1].LeptonOrigin == 10) leptonFlavor = 1;
	  }
        }
           

      // CandType setting
      if( (WmuCand.size()+ZmumuCand.size()+WeCand.size()+ZeeCand.size()) != 3 ) continue;

      if(ZmumuCand.size()==2 && WmuCand.size()==1) CandType="mumumu";
      if(ZeeCand.size()==2 && WeCand.size()==1) CandType="eee";
      if(ZmumuCand.size()==2 && WeCand.size()==1) CandType="mumue";
      if(ZeeCand.size()==2 && WmuCand.size()==1) CandType="eemu";

      // try with >= 2 leptons for Zjets
      //if( selMuons.size()+selElectrons.size() <2 ) continue; // 2lept incl
      //if( selMuons.size()+selElectrons.size() <3 ) continue; // 3lept incl
/*      if( selMuons.size()+selElectrons.size() !=3 ) continue; // 3lept excl 

      if(ChannelName=="mumu" && ZmumuCand.size()==2) CandType="mumumu";
      if(ChannelName=="ee" && ZeeCand.size()==2) CandType="eee";
      if(ChannelName=="emu" && selMuons.size()>=1 && selElectrons.size()>=1) CandType="mumue";
      if(ChannelName=="emu" && selMuons.size()>=1 && selElectrons.size()>=1) CandType="eemu";
*/      //std::cout<<ChannelName<<"  mu "<<selMuons.size()<<" mumu "<<ZmumuCand.size()<<" e "<<selElectrons.size()<<" ee "<<ZeeCand.size()<<"  "<<CandType<<endl; 

      // try >= 3 leptons for Data eee
/*      if( selElectrons.size()+selMuons.size() < 3 ) continue;

      if(ZmumuCand.size()==2 && WmuCand.size()>=1 && ZeeCand.size()==0 && WeCand.size()==0) CandType="mumumu";
      if(ZeeCand.size()==2 && WeCand.size()>=1) CandType="eee";
      if(ZmumuCand.size()==2 && WeCand.size()>=1 && ZeeCand.size()==0 && WmuCand.size()==0) CandType="mumue";
      if(ZeeCand.size()==2 && WmuCand.size()>=1 && ZmumuCand.size()==0 && WeCand.size()==0) CandType="eemu";
      if(CandType=="eee") std::cout<<" mu "<<selMuons.size()<<" e "<<selElectrons.size()<<std::endl;
*/

      //if(datasets[d].Name()=="Zjets") std::cout<< "Zjets : cand "<<CandType<<" trig "<<ChannelName<<" "<<passTrigger<<"  mumu "<<ZmumuCand.size()<<" mu "<<WmuCand.size()<<" ee "<<ZeeCand.size()<<" e "<<WeCand.size()<<std::endl;

      
      if(WeCand.size()) candElec = WeCand;
      else candElec = selElectrons;
      if(WmuCand.size()) candMuon = WmuCand;
      else candMuon = selMuons;


      
       
       //***********************************************
       //check correlations of leptons and Btag triggers
       //***********************************************   
     
      double passTriggerBTAGMET = 0;
      double passTriggerMET     = 0;
      const NTTrigger* triggers = sel.GetPointer2Trigger();
      //triggers->Dump();
      

      std::pair<bool, bool> trigBTagMET_MET = sel.passMETTriggerSelection (&datasets[d]);
      if(trigBTagMET_MET.first  == true) passTriggerBTAGMET = 1;
      if(trigBTagMET_MET.second == true) passTriggerMET     = 1;
      
      if (ChannelName=="mumu" && CandType=="mumumu"){ //mumunu
        MyhistoManager.FillHisto2D (correlationBTAGMET_mumumu, "correlationBTAGMET_mumumu", passTriggerBTAGMET,passTrigger , datasets[d].Name (), IsSignal, weightITypeMC );
        MyhistoManager.FillHisto2D (correlationMET_mumumu,     "correlationMET_mumumu",     passTriggerMET,    passTrigger , datasets[d].Name (), IsSignal, weightITypeMC );
      }
      if (ChannelName=="ee" && CandType=="eee"){ //eee
        MyhistoManager.FillHisto2D (correlationBTAGMET_eee, "correlationBTAGMET_eee", passTriggerBTAGMET,passTrigger , datasets[d].Name (), IsSignal, weightITypeMC );
        MyhistoManager.FillHisto2D (correlationMET_eee,     "correlationMET_eee",     passTriggerMET,    passTrigger , datasets[d].Name (), IsSignal, weightITypeMC );
      }
      if (ChannelName=="emu" && CandType=="mumue"){ //mumue
        MyhistoManager.FillHisto2D (correlationBTAGMET_mumue, "correlationBTAGMET_mumue", passTriggerBTAGMET,passTrigger , datasets[d].Name (), IsSignal, weightITypeMC );
        MyhistoManager.FillHisto2D (correlationMET_mumue,     "correlationMET_mumue",     passTriggerMET,    passTrigger , datasets[d].Name (), IsSignal, weightITypeMC );
      }
      if (ChannelName=="emu" && CandType=="eemu"){ //eemu
        MyhistoManager.FillHisto2D (correlationBTAGMET_eemu, "correlationBTAGMET_eemu", passTriggerBTAGMET,passTrigger , datasets[d].Name (), IsSignal, weightITypeMC );
        MyhistoManager.FillHisto2D (correlationMET_eemu,     "correlationMET_eemu",     passTriggerMET,    passTrigger , datasets[d].Name (), IsSignal, weightITypeMC );
      }
       
       
       
       


       if( (ChannelName=="ee" && CandType=="eee" && isData) ||
           (ChannelName=="ee" && CandType=="eee" && !isData && passTriggerMET==1)){
         //cout << "    in ee  " << CandType << " candidate passing trigger " <<  passTrigger << endl;
	 MyhistoManager.FillHisto(HTriggerEff_ElecSel_pT,  "HTriggerEff_ElecSel_pT",  candElec[0].p4.Pt(), datasets[d].Name(), IsSignal, weightITypeMC);
	 MyhistoManager.FillHisto(HTriggerEff_ElecSel_eta, "HTriggerEff_ElecSel_eta", candElec[0].p4.Eta(),datasets[d].Name(), IsSignal, weightITypeMC);
	 MyhistoManager.FillHisto(HTriggerEff_ElecSel_phi, "HTriggerEff_ElecSel_phi", candElec[0].p4.Phi(),datasets[d].Name(), IsSignal, weightITypeMC);
	 
         MyhistoManager.FillHisto(HNmuons_ElecSel,  "HNmuons_ElecSel",  selMuons.size(), datasets[d].Name(), IsSignal, weightITypeMC); 
         MyhistoManager.FillHisto(HNelec_ElecSel,  "HNelec_ElecSel",  selElectrons.size(), datasets[d].Name(), IsSignal, weightITypeMC);

	 if( passTrigger ){
	   MyhistoManager.FillHisto(HTriggerEff_ElecSel_selTrig_pT,  "HTriggerEff_ElecSel_selTrig_pT",  candElec[0].p4.Pt(), datasets[d].Name(),
				    IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(HTriggerEff_ElecSel_selTrig_eta, "HTriggerEff_ElecSel_selTrig_eta", candElec[0].p4.Eta(),datasets[d].Name(),
				    IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(HTriggerEff_ElecSel_selTrig_phi, "HTriggerEff_ElecSel_selTrig_phi", candElec[0].p4.Phi(),datasets[d].Name(),
				    IsSignal, weightITypeMC);
	 }
       }
       
       
       
       if( (ChannelName=="mumu" && CandType=="mumumu" && isData) ||
           (ChannelName=="mumu" && CandType=="mumumu" && !isData && passTriggerMET==1)){
        // cout << "    in mumu  " << CandType << " candidate passing trigger " <<  passTrigger << endl;
	 MyhistoManager.FillHisto(HTriggerEff_MuonSel_pT,  "HTriggerEff_MuonSel_pT",  candMuon[0].p4.Pt(), datasets[d].Name(), IsSignal, weightITypeMC);
	 MyhistoManager.FillHisto(HTriggerEff_MuonSel_eta, "HTriggerEff_MuonSel_eta", candMuon[0].p4.Eta(),datasets[d].Name(), IsSignal, weightITypeMC);
	 MyhistoManager.FillHisto(HTriggerEff_MuonSel_phi, "HTriggerEff_MuonSel_phi", candMuon[0].p4.Phi(),datasets[d].Name(), IsSignal, weightITypeMC);
	 
         MyhistoManager.FillHisto(HNmuons_MuonSel,  "HNmuons_MuonSel",  selMuons.size(), datasets[d].Name(), IsSignal, weightITypeMC);
         MyhistoManager.FillHisto(HNelec_MuonSel,  "HNelec_MuonSel",  selElectrons.size(), datasets[d].Name(), IsSignal, weightITypeMC);

	 if(   passTrigger   ){
	   MyhistoManager.FillHisto(HTriggerEff_MuonSel_selTrig_pT,  "HTriggerEff_MuonSel_selTrig_pT",  candMuon[0].p4.Pt(), datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(HTriggerEff_MuonSel_selTrig_eta, "HTriggerEff_MuonSel_selTrig_eta", candMuon[0].p4.Eta(),datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(HTriggerEff_MuonSel_selTrig_phi, "HTriggerEff_MuonSel_selTrig_phi", candMuon[0].p4.Phi(),datasets[d].Name(), IsSignal, weightITypeMC);
	 }
       }
       
       
       if( (ChannelName=="emu" && CandType=="emu"  && isData) ||
           (ChannelName=="emu" && CandType=="emu"  && !isData && passTriggerMET==1)){
	 
          // cout << "    in emu  " << CandType << " candidate passing trigger " <<  passTrigger << endl;
	   //MyhistoManager.FillHisto(HInvM_ll_pair, "InvDilMassPair" ,InvDilMass,datasets[d].Name(), IsSignal, 1); 
	   MyhistoManager.FillHisto(HTriggerEff_afterDileptSel_pT, "HTriggerEff_afterDileptSel_pT" , candMuon[0].p4.Pt() , datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(HTriggerEff_afterDileptSel_eta,"HTriggerEff_afterDileptSel_eta" , candMuon[0].p4.Eta(), datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(HTriggerEff_afterDileptSel_phi,"HTriggerEff_afterDileptSel_phi" , candMuon[0].p4.Phi(), datasets[d].Name(), IsSignal, weightITypeMC);
	   
	   if ( passTrigger && CandType=="emu"){
	     MyhistoManager.FillHisto(HTriggerEff_afterDileptSel_selTrig_pT,  "HTriggerEff_afterDileptSel_selTrig_pT" ,  candMuon[0].p4.Pt(), datasets[d].Name(), IsSignal, weightITypeMC);
	     MyhistoManager.FillHisto(HTriggerEff_afterDileptSel_selTrig_eta, "HTriggerEff_afterDileptSel_selTrig_eta" , candMuon[0].p4.Eta(), datasets[d].Name(), IsSignal, weightITypeMC);
	     MyhistoManager.FillHisto(HTriggerEff_afterDileptSel_selTrig_phi, "HTriggerEff_afterDileptSel_selTrig_phi" , candMuon[0].p4.Phi(), datasets[d].Name(), IsSignal, weightITypeMC);
	   }
       }
     
       if( (ChannelName=="emu" && CandType=="mumue"  && isData) ||
           (ChannelName=="emu" && CandType=="mumue"  && !isData && passTriggerMET==1)){
	 
	   MyhistoManager.FillHisto(HTriggerEff_DiMuonElecSel_pT, "HTriggerEff_DiMuonElecSel_pT" , candElec[0].p4.Pt() , datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(HTriggerEff_DiMuonElecSel_eta,"HTriggerEff_DiMuonElecSel_eta" , candElec[0].p4.Eta(), datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(HTriggerEff_DiMuonElecSel_phi,"HTriggerEff_DiMuonElecSel_phi" , candElec[0].p4.Phi(), datasets[d].Name(), IsSignal, weightITypeMC);
	   
           MyhistoManager.FillHisto(HNmuons_DiMuonElecSel,  "HNmuons_DiMuonElecSel",  selMuons.size(), datasets[d].Name(), IsSignal, weightITypeMC); 
           MyhistoManager.FillHisto(HNelec_DiMuonElecSel,  "HNelec_DiMuonElecSel",  selElectrons.size(), datasets[d].Name(), IsSignal, weightITypeMC);

	   if ( passTrigger && CandType=="mumue"){
	     MyhistoManager.FillHisto(HTriggerEff_DiMuonElecSel_selTrig_pT,  "HTriggerEff_DiMuonElecSel_selTrig_pT" , candElec[0].p4.Pt(), datasets[d].Name(), IsSignal, weightITypeMC);
	     MyhistoManager.FillHisto(HTriggerEff_DiMuonElecSel_selTrig_eta, "HTriggerEff_DiMuonElecSel_selTrig_eta" , candElec[0].p4.Eta(), datasets[d].Name(), IsSignal, weightITypeMC);
	     MyhistoManager.FillHisto(HTriggerEff_DiMuonElecSel_selTrig_phi, "HTriggerEff_DiMuonElecSel_selTrig_phi" , candElec[0].p4.Phi(), datasets[d].Name(), IsSignal, weightITypeMC);
	   }
       }
     
       if( (ChannelName=="emu" && CandType=="eemu"  && isData) ||
           (ChannelName=="emu" && CandType=="eemu"  && !isData && passTriggerMET==1)){
	 
	   MyhistoManager.FillHisto(HTriggerEff_DiElecMuonSel_pT, "HTriggerEff_DiElecMuonSel_pT" , candMuon[0].p4.Pt() , datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(HTriggerEff_DiElecMuonSel_eta,"HTriggerEff_DiElecMuonSel_eta" , candMuon[0].p4.Eta(), datasets[d].Name(), IsSignal, weightITypeMC);
	   MyhistoManager.FillHisto(HTriggerEff_DiElecMuonSel_phi,"HTriggerEff_DiElecMuonSel_phi" , candMuon[0].p4.Phi(), datasets[d].Name(), IsSignal, weightITypeMC);
	   
           MyhistoManager.FillHisto(HNmuons_DiElecMuonSel,  "HNmuons_DiElecMuonSel",  selMuons.size(), datasets[d].Name(), IsSignal, weightITypeMC);
           MyhistoManager.FillHisto(HNelec_DiElecMuonSel,  "HNelec_DiElecMuonSel",  selElectrons.size(), datasets[d].Name(), IsSignal, weightITypeMC);

	   if ( passTrigger && CandType=="eemu"){
	     MyhistoManager.FillHisto(HTriggerEff_DiElecMuonSel_selTrig_pT,  "HTriggerEff_DiElecMuonSel_selTrig_pT" ,  candMuon[0].p4.Pt(), datasets[d].Name(), IsSignal, weightITypeMC);
	     MyhistoManager.FillHisto(HTriggerEff_DiElecMuonSel_selTrig_eta, "HTriggerEff_DiElecMuonSel_selTrig_eta" , candMuon[0].p4.Eta(), datasets[d].Name(), IsSignal, weightITypeMC);
	     MyhistoManager.FillHisto(HTriggerEff_DiElecMuonSel_selTrig_phi, "HTriggerEff_DiElecMuonSel_selTrig_phi" , candMuon[0].p4.Phi(), datasets[d].Name(), IsSignal, weightITypeMC);
	   }
       }
     }
     
     // cout some values based on objects
     //if(ievt == 10000) break ;
    }  // end of loop over evts
    TString outputFileName = "tagAndProbeHisto_"+datasets[d].Name()+".root";
    theEffPlots.Savehistos(outputFileName);
    
    
    
    }		
		// end of loop over the datasets 
  cout<<"#########################"<<endl;
  cout<<" Loop over the datasets  "<<endl;
  cout<<"#########################"<<endl;
 
  TFile* fout  = new TFile("TriggerPlots.root","RECREATE");


  MyhistoManager.WriteMyHisto(HTriggerEff_afterDileptSel_pT ,"all");
  MyhistoManager.WriteMyHisto(HTriggerEff_afterDileptSel_eta,"all");
  MyhistoManager.WriteMyHisto(HTriggerEff_afterDileptSel_phi,"all");
  MyhistoManager.WriteMyHisto(HTriggerEff_afterDileptSel_selTrig_pT ,"all");
  MyhistoManager.WriteMyHisto(HTriggerEff_afterDileptSel_selTrig_eta,"all");
  MyhistoManager.WriteMyHisto(HTriggerEff_afterDileptSel_selTrig_phi,"all");
   
  MyhistoManager.WriteMyHisto(deltaR_Jet_Lepton_FullSel,"all"  );         
  MyhistoManager.WriteMyHisto(deltaR_Jet_MET_FullSel   ,"all"  );	    
  MyhistoManager.WriteMyHisto(deltaR_Muon_MET_FullSel  ,"all"  );	      
  MyhistoManager.WriteMyHisto(deltaR_Elec_MET_FullSel  ,"all"  );	  
  
  MyhistoManager.WriteMyHisto(deltaPhi_Jet_Lepton_FullSel,"all"  );         
  MyhistoManager.WriteMyHisto(deltaPhi_Jet_MET_FullSel   ,"all"  );	    
  MyhistoManager.WriteMyHisto(deltaPhi_Muon_MET_FullSel  ,"all"  );	      
  MyhistoManager.WriteMyHisto(deltaPhi_Elec_MET_FullSel  ,"all"  );	      
  MyhistoManager.WriteMyHisto(deltaPhi_DiLept_MET_FullSel  ,"all"  );	
  
  MyhistoManager.WriteMyHisto(HNmuons_ElecSel,  "all");
  MyhistoManager.WriteMyHisto(HNelec_ElecSel,  "all");
  MyhistoManager.WriteMyHisto(HNmuons_MuonSel,  "all");
  MyhistoManager.WriteMyHisto(HNelec_MuonSel,  "all");
  MyhistoManager.WriteMyHisto(HNmuons_DiMuonElecSel,  "all");
  MyhistoManager.WriteMyHisto(HNelec_DiMuonElecSel,  "all");
  MyhistoManager.WriteMyHisto(HNmuons_DiElecMuonSel,  "all");
  MyhistoManager.WriteMyHisto(HNelec_DiElecMuonSel,  "all");

  MyhistoManager.WriteMyHisto(HTriggerEff_ElecSel_pT,  "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_ElecSel_eta, "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_ElecSel_phi, "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_ElecSel_selTrig_pT,  "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_ElecSel_selTrig_eta, "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_ElecSel_selTrig_phi, "all");
  
  MyhistoManager.WriteMyHisto(HTriggerEff_MuonSel_pT,  "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_MuonSel_eta, "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_MuonSel_phi, "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_MuonSel_selTrig_pT,  "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_MuonSel_selTrig_eta, "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_MuonSel_selTrig_phi, "all");
   
  MyhistoManager.WriteMyHisto(HTriggerEff_DiMuonElecSel_pT,  "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_DiMuonElecSel_eta, "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_DiMuonElecSel_phi, "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_DiMuonElecSel_selTrig_pT,  "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_DiMuonElecSel_selTrig_eta, "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_DiMuonElecSel_selTrig_phi, "all");
   
  MyhistoManager.WriteMyHisto(HTriggerEff_DiElecMuonSel_pT,  "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_DiElecMuonSel_eta, "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_DiElecMuonSel_phi, "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_DiElecMuonSel_selTrig_pT,  "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_DiElecMuonSel_selTrig_eta, "all");
  MyhistoManager.WriteMyHisto(HTriggerEff_DiElecMuonSel_selTrig_phi, "all");
   
  MyhistoManager.WriteMyHisto(genEl_Tight_pt,  "all");
  MyhistoManager.WriteMyHisto(genEl_Tight_eta, "all");
  MyhistoManager.WriteMyHisto(genEl_Tight_phi, "all");
  MyhistoManager.WriteMyHisto(genEl_Tight_njet,"all");
  
  MyhistoManager.WriteMyHisto(genEl_LooseID_pt,  "all");
  MyhistoManager.WriteMyHisto(genEl_LooseID_eta, "all");
  MyhistoManager.WriteMyHisto(genEl_LooseID_phi, "all");
  MyhistoManager.WriteMyHisto(genEl_LooseID_njet,"all");
  
  MyhistoManager.WriteMyHisto(genEl_LooseIso_pt,  "all");
  MyhistoManager.WriteMyHisto(genEl_LooseIso_eta, "all");
  MyhistoManager.WriteMyHisto(genEl_LooseIso_phi, "all");
  MyhistoManager.WriteMyHisto(genEl_LooseIso_njet,"all");
  
  
  MyhistoManager.WriteMyHisto(genMu_Tight_pt,  "all");
  MyhistoManager.WriteMyHisto(genMu_Tight_eta, "all");
  MyhistoManager.WriteMyHisto(genMu_Tight_phi, "all");
  MyhistoManager.WriteMyHisto(genMu_Tight_njet,"all");
  
  MyhistoManager.WriteMyHisto(genMu_LooseID_pt,  "all");
  MyhistoManager.WriteMyHisto(genMu_LooseID_eta, "all");
  MyhistoManager.WriteMyHisto(genMu_LooseID_phi, "all");
  MyhistoManager.WriteMyHisto(genMu_LooseID_njet,"all");
  
  MyhistoManager.WriteMyHisto(genMu_LooseIso_pt,  "all");
  MyhistoManager.WriteMyHisto(genMu_LooseIso_eta, "all");
  MyhistoManager.WriteMyHisto(genMu_LooseIso_phi, "all");
  MyhistoManager.WriteMyHisto(genMu_LooseIso_njet,"all");
  MyhistoManager.WriteMyHisto2D(correlationBTAGMET_mumumu,"all");
  MyhistoManager.WriteMyHisto2D(correlationMET_mumumu,"all");
  MyhistoManager.WriteMyHisto2D(correlationBTAGMET_eee,"all");
  MyhistoManager.WriteMyHisto2D(correlationMET_eee,"all");
  MyhistoManager.WriteMyHisto2D(correlationBTAGMET_mumue,"all");
  MyhistoManager.WriteMyHisto2D(correlationMET_mumue,"all");
  MyhistoManager.WriteMyHisto2D(correlationBTAGMET_eemu,"all");
  MyhistoManager.WriteMyHisto2D(correlationMET_eemu,"all");
 
  ////////////////////////////
  //  Computation after loops
  ////////////////////////////





  ///////////////


  return (0);
}
