#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "test/TMVAGui.C"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"


void BDT_FCNC_tZ(){
  
  
  //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();
  
  
  
  //---------------------------------------------------------------
  //  Output file
  TString outfileName( "bdtTMVA_FCNC_tZ.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  
  //---------------------------------------------------------------
  // Initialize the TMVA factory
  TMVA::Factory *factory = new TMVA::Factory( "doBDT_FCNC_tZ", outputFile,
					      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
  
  
  
  //---------------------------------------------------------------
  //Events weight, to filled with the cutflow (nbr of signal and bckg events)
  Double_t signalWeight     = 1;
  Double_t backgroundWeight = 1;
  
  //---------------------------------------------------------------
  //  Input files which contain the TTrees for signal and bckg
  TFile *input_sig = TFile::Open( "proof.root" );
  TFile *input_wz = TFile::Open( "proof.root" );
  
  //---------------------------------------------------------------
  //  Get the TTrees for signal and bckg
  TTree *signal     = (TTree*)input_sig->Get("Ttree_FCNCkut");
  TTree *background_WZ = (TTree*)input_wz->Get("Ttree_WZ");
  
  //---------------------------------------------------------------
  //Add the TTree to the TMVA::factory
  factory->AddSignalTree    ( signal,            signalWeight     );
  factory->AddBackgroundTree( background_WZ,     backgroundWeight );
  
  
  //---------------------------------------------------------------
  // Define the variables to take into account for the trainning
  // the variables should be branches of the trees (same name for signal and background)
  factory->AddVariable("tree_topMass",    'F');
  factory->AddVariable("tree_deltaPhilb", 'F');
  factory->AddVariable("tree_asym",       'F');
  factory->AddVariable("tree_Zpt",        'F');
  
  
  
  
  
  
  
  
  //--------------------------------------------------------------- 
  // Apply additional cuts on the signal and background samples (can be different)
  TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";
  
  factory->PrepareTrainingAndTestTree( mycuts, mycutb,
				       "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
  
  
  
  factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );
  
  
  
  
  // Train MVAs using the set of training events
  factory->TrainAllMethods();
  
  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();
  
  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();
  
  // --------------------------------------------------------------
  
  // Save the output
  outputFile->Close();
  
  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;
  
  delete factory;
  
  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVAGui( outfileName );
  
  
}
