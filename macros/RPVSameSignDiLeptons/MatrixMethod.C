#include <iomanip>
#include <iostream>
#include "NTFormat/interface/NTEvent.h"

//NTupleAnalysis classes
#include "Selection/interface/SSDiLeptonSelection.h"
#include "Tools/interface/Dataset.h"
#include "Tools/interface/AnalysisEnvironmentLoader.h"
#include "BckgdEstimation/interface/MMEstimation.h"
#include "Tools/interface/PUWeighting.h"


#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

using namespace IPHCTree;
using namespace std;

void MatrixMethod (string xmlFileName, bool mc=true, string channel="") // MuMu or EE for data, do both for MC
{
  cout<<"#########################"<<endl;
  cout<<"Beginning of the program"<<endl;
  cout<<"#########################"<<endl;
  
  //////////////////////
  // Initialisation
  //////////////////////
  float Luminosity = 0;
  float LumiError = 0;
  string PUWeightFileName;
  int DataType = 0; 

  AnalysisEnvironmentLoader anaEL (xmlFileName);
  vector < Dataset > datasets;
  anaEL.LoadSamples (datasets); // now the list of datasets written in the xml file is known
  int verbosity = -1;

  // Matrix Method
  int MMselStepCut = 3;
  float looseIsoMM = 0.8; float tightIsoMMmu = 0.2; float tightIsoMMe = 0.17;
  unsigned int nbinsMM = 11; float lowEdgeMM = -0.5; float highEdgeMM = 10.5; // Histo nb events en fonction #jets
  MMEstimation MMestEE(datasets, looseIsoMM, tightIsoMMe, nbinsMM, lowEdgeMM, highEdgeMM, "EE"); 
  MMEstimation MMestMuMu(datasets,looseIsoMM, tightIsoMMmu, nbinsMM, lowEdgeMM, highEdgeMM, "MuMu"); 

  SSDiLeptonSelection sel;
  anaEL.LoadSSDiLeptonSelection (sel); // now the parameters for the selection are given to the selection
  anaEL.LoadGeneralInfo(DataType, Luminosity, LumiError, PUWeightFileName, verbosity );
  //Load for PU:
  if(mc) sel.GeneratePUWeight(PUWeightFileName);

  IPHCTree::NTEvent * event = 0;

  /*PUWeighting  thePUReweighter;
  TFile* file1  = new TFile(PUWeightFileName.c_str(),"READ"); 
  TH1D *  hPUData = 0;
  hPUData         = (TH1D*)file1->Get("pileup");
  TH1F *  hPUMC   = new TH1F("pileup_MC", "pileup_MC", hPUData->GetXaxis()->GetNbins(), hPUData->GetXaxis()->GetXmin(), hPUData->GetXaxis()->GetXmax() );
  TFile* file2  = new TFile( "../data/CrossSection_pileup.root" ,"READ");
  hPUMC           = (TH1F*)file2->Get("pileup_TTbarSig");
  // histo in data, histo in Mc, use out-of-time pu in the reweighting
  cout << "get MC histo  " << endl;
  thePUReweighter.setPUHisto( hPUData, hPUMC);
  cout << "set MC histo in thePUReweighter " << endl;
  thePUReweighter.setUseOutOfTimePU(false); // set to true to use out-of-time PU
  */
  
  //////////////////////
  //LOOP OVER THE DATASETS
  //////////////////////
  cout<<"#########################"<<endl;
  cout<<" Loop over the datasets  "<<endl;
  cout<<"#########################"<<endl;

  for (unsigned int d = 0; d < datasets.size (); d++) {
   datasets[d].eventTree ()->SetBranchAddress ("NTEvent",&event);
   cout << "dataset = " << datasets[d].Name() << endl;
   unsigned int nEvents = (int) (datasets[d].eventTree ()->GetEntries ());
    cout << "NEvents = " << nEvents << endl;
    float weight_init;
    weight_init = datasets[d].NormFactor()*Luminosity;
    cout << "weight_init = " << weight_init << endl;
    //LOOP OVER THE EVENTS
    for (unsigned int ievt = 0; ievt < nEvents; ievt++) {
      datasets[d].eventTree ()->GetEntry (ievt);
      IPHCTree::NTTransient::InitializeAfterReading(event); // Important line to read new format files
      //Load event for the selection
      sel.LoadEvent(event);
      if(ievt%10000 == 0) cout << "number of processed events " << ievt << endl;
      
      float weight = 1;
      
      if(mc){
      	
	//Manage DY samples to avoid overlaps
	double dileptInvMass = 0;
	if( (event->mc.zAndDecays).size() > 0){
          TLorentzVector dilept = (event->mc.zAndDecays)[0].p4_Lep1_gen + (event->mc.zAndDecays)[0].p4_Lep2_gen;
          dileptInvMass = dilept.M();
	}
	if(datasets[d].Name()=="Zjets" && dileptInvMass < 50) continue;
	if(datasets[d].Name()=="DYToMuMu_M-20"	   && (dileptInvMass > 50 || dileptInvMass < 20) ) continue;
	if(datasets[d].Name()=="DYToEE_M-20"        && (dileptInvMass > 50 || dileptInvMass < 20) ) continue;
	if(datasets[d].Name()=="DYToTauTau_M-20"    && (dileptInvMass > 50 || dileptInvMass < 20) ) continue;
	if(datasets[d].Name()=="DYToMuMu_M-10To20"   &&  dileptInvMass > 20) continue;
	if(datasets[d].Name()=="DYToEE_M-10To20"    &&  dileptInvMass > 20) continue;
	if(datasets[d].Name()=="DYToTauTau_M-10To20" &&  dileptInvMass > 20) continue;


        weight = weight_init;
        //float weight = weight_init*sel.GetPUWeight();
        /*if(thePUReweighter.getUseOutOfTimePU()){
           weight = weight_init*thePUReweighter.weight(event->pileup.intime_npu, event->general.runNb);
         }else{
	   weight = weight_init*thePUReweighter.weight(event->pileup.intime_npu);
       }*/
      }
      
      if(mc || (!mc && channel=="EE")) MMestEE.CountNSel(sel, datasets[d], "ee_ss", weight, MMselStepCut, &(event->mc));
      if(mc || (!mc && channel=="MuMu")) MMestMuMu.CountNSel(sel, datasets[d], "mumu_ss", weight, MMselStepCut, &(event->mc));
    } // end of loop over evts
  } // end of loop over the datasets 
  cout<<"#########################"<<endl;
  cout<<" Loop over the datasets  "<<endl;
  cout<<"#########################"<<endl;
  
  ////////////////////////////
  //  Computation after loops
  ////////////////////////////
  // Matrix Method estimation for ee and mumu cases 
  vector<struct MMEpsilons> valMMEpsilons; 
  for(unsigned int bin_index = 0; bin_index < nbinsMM; bin_index++){
    struct MMEpsilons valMMEpsilonsTmp;
    valMMEpsilonsTmp.EpsilonSignal = 0.99;
    valMMEpsilonsTmp.EpsilonSignalErr = 0.05;
    valMMEpsilonsTmp.EpsilonFake = 0.20;
    valMMEpsilonsTmp.EpsilonFakeErr = 0.10;
    valMMEpsilons.push_back(valMMEpsilonsTmp);
  }
  unsigned int NbIterations = 10000;
  bool doStatistical = true; bool doSystematic = true; bool doCorrections = false;
  string filename = "MatrixMethod_OutPut_";
  string mcname = "MC";
  if(!mc) mcname = "DATA";
  
  if(mc || (!mc && channel=="EE")) {
    MMestEE.RunTheMatrixMethod(valMMEpsilons, NbIterations, doStatistical, doSystematic, doCorrections);
    MMestEE.PrintMMEstimated();
    MMestEE.WriteMMFile((filename+"EE_"+mcname+".root").c_str());
  }
  if(mc || (!mc && channel=="MuMu")) {
    MMestMuMu.RunTheMatrixMethod(valMMEpsilons,  NbIterations, doStatistical, doSystematic, doCorrections);
    MMestMuMu.PrintMMEstimated();
    MMestMuMu.WriteMMFile((filename+"MuMu_"+mcname+".root").c_str());
  }

  cout<<"#########################"<<endl;
  cout<<"    End of the program   "<<endl;
  cout<<"#########################"<<endl;

  return;
}

int main ()
{
  MatrixMethod("../../config/MatrixMethod_MC.xml", true);
  // Ancienne macro MC_NoWeight <-> Tourner sur les MC en mode data
  
  return (0);
}
