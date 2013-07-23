



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


void testBDT(){


   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();
/*
  TString outfileName( "TMVA.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "testBDT", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );


 
   // global event weights per tree (see below for setting event-wise weights)
   //Double_t signalWeight     = 0.003582;
   //Double_t backgroundWeight = 0.0269;
   
   Double_t signalWeight     = 1;
   Double_t backgroundWeight = 1;
   
   TFile *input_sig = TFile::Open( "signal_exclusif.root" );
   TFile *input_wz = TFile::Open( "bruit_w_z.root" );
   
   TTree *signal     = (TTree*)input_sig->Get("tree");
   TTree *background = (TTree*)input_wz->Get("tree");

   // You can add an arbitrary number of signal or background trees
   factory->AddSignalTree    ( signal,     signalWeight     );
   factory->AddBackgroundTree( background, backgroundWeight );
   
   
   factory->AddVariable("PT_z"   , 'F');
   factory->AddVariable("ASYM"    , 'F');
   factory->AddVariable("PHI_lw_b", 'F');
   factory->AddVariable("M_top", 'F');
   */
   
   
   
   TString outfileName( "bdtTMVA_FCNC_tZ.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "doBDT_FCNC_tZ", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );


 
   // global event weights per tree (see below for setting event-wise weights)
   //Double_t signalWeight     = 0.003582;
   //Double_t backgroundWeight = 0.0269;
   
   Double_t signalWeight     = 1;
   Double_t backgroundWeight = 1;
   
   TFile *input_sig = TFile::Open( "proof.root" );
   TFile *input_wz = TFile::Open( "proof.root" );
   
   TTree *signal     = (TTree*)input_sig->Get("Ttree_FCNCkut");
   
   
   TTree *background_WZ = (TTree*)input_wz->Get("Ttree_WZ");
   /*TTree *background_ZZ = (TTree*)input_wz->Get("Ttree_ZZ");
   TTree *background_WW = (TTree*)input_wz->Get("Ttree_WW");
   
   TTree *background_TTbar  = (TTree*)input_wz->Get("Ttree_TTbar");
   TTree *background_Zjets  = (TTree*)input_wz->Get("Ttree_Zjets");
   TTree *background_Wjets  = (TTree*)input_wz->Get("Ttree_Wjets");
   TTree *background_TtW    = (TTree*)input_wz->Get("Ttree_TtW");
   TTree *background_TbartW = (TTree*)input_wz->Get("Ttree_TbartW");*/

   // You can add an arbitrary number of signal or background trees
   factory->AddSignalTree    ( signal,            signalWeight     );
   factory->AddBackgroundTree( background_WZ,     backgroundWeight );
   /*factory->AddBackgroundTree( background_ZZ,     backgroundWeight );
   factory->AddBackgroundTree( background_WW,     backgroundWeight );
   factory->AddBackgroundTree( background_TTbar,  backgroundWeight );
   factory->AddBackgroundTree( background_Zjets,  backgroundWeight );
   factory->AddBackgroundTree( background_Wjets,  backgroundWeight );
   factory->AddBackgroundTree( background_TtW,    backgroundWeight );
   factory->AddBackgroundTree( background_TbartW, backgroundWeight );*/
   
   
   factory->AddVariable("tree_topMass",    'F');
   factory->AddVariable("tree_deltaPhilb", 'F');
   factory->AddVariable("tree_asym",       'F');
   factory->AddVariable("tree_Zpt",        'F');
   
   
   
   
   
   
   
   
   
   // to set weights. The variable must exist in the tree
   //    for signal    : factory->SetSignalWeightExpression    ("weight1*weight2");
   //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
   
   
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
