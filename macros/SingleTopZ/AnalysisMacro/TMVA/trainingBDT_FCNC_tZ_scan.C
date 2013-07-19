#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <sstream>

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


void doBDT_FCNC_tZ_scan(TString systtype, double Zut, double Zct){


   std::ostringstream ossZut;
   ossZut << Zut ;
   TString sZut = ossZut.str();
  
   std::ostringstream ossZct;
   ossZct << Zct ;
   TString sZct = ossZct.str();
  
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   TString outfileName( ("ScanTraining"+systtype+"/trainingBDT_FCNC_tZ_"+sZut+"_"+sZct+"_"+systtype+".root").Data() );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::Factory *factory = new TMVA::Factory( "BDT_trainning_zut", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );


 
   // global event weights per tree (see below for setting event-wise weights)
   //Double_t signalWeight     = 0.003582;
   //Double_t backgroundWeight = 0.0269;
   
   Double_t signalWeight     = 1;
   Double_t backgroundWeight = 1;
   
   TFile *input_sig   = TFile::Open( ("FilesBenchmark/TreeMVA_benchmark_"+sZut+"_"+sZct+"_"+systtype+".root").Data());
   TFile *input_wz    = TFile::Open( "../RootFiles/backup_outputProof02-01-13_18-48-30_AllSFinclWZ/proof_merged.root" );
   
   input_sig->cd();
   TTree *signal= (TTree*)input_sig->Get(("TreeMVA_benchmark_"+sZut+"_"+sZct).Data() );
   
              
   input_wz->cd();
   TTree *background_WZ     = (TTree*)input_wz->Get("Ttree_WZ");
   
   
    // You can add an arbitrary number of signal or background trees
   factory->AddSignalTree    ( signal,            signalWeight     );
   factory->AddBackgroundTree( background_WZ,     backgroundWeight );
     
   
   factory->AddVariable("tree_topMass",    'F'); 
   factory->AddVariable("tree_totMass",    'F');  
   factory->AddVariable("tree_deltaPhilb", 'F');
   factory->AddVariable("tree_deltaRlb",   'F');
   factory->AddVariable("tree_deltaRTopZ", 'F');
   factory->AddVariable("tree_asym",       'F');
   factory->AddVariable("tree_Zpt",        'F');
   factory->AddVariable("tree_ZEta",       'F');
   factory->AddVariable("tree_topPt",      'F');
   factory->AddVariable("tree_topEta",     'F');
   factory->AddVariable("tree_NJets",        'F');
   factory->AddVariable("tree_NBJets",         'F');
   factory->AddVariable("tree_deltaRZl ",     'F');   
   factory->AddVariable("tree_deltaPhiZmet",  'F');
   factory->AddVariable("tree_btagDiscri",    'F');  	    
   factory->AddVariable("tree_leptWPt",       'F');  		   
   factory->AddVariable("tree_leptWEta",      'F');  		   
   factory->AddVariable("tree_leadJetPt",     'F');  	     
   factory->AddVariable("tree_leadJetEta",    'F');  	    
   //factory->AddVariable("tree_deltaRZleptW",  'F');  	  
   factory->AddVariable("tree_deltaPhiZleptW",'F');  	
   
   
   
         
   
   // to set weights. The variable must exist in the tree
   //    for signal    : factory->SetSignalWeightExpression    ("weight1*weight2");
   //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
   
   
   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   
   
   
   //factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );
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
   //if (!gROOT->IsBatch()) TMVAGui( outfileName );


}



void trainingBDT_FCNC_tZ_scan(){

  std::vector<double> bench_Zut, bench_Zct;
 
  for(int bench=0; bench<6;   bench++){
    bench_Zut.push_back((bench)*0.00004);
  }
  for(int bench=0; bench<6; bench++){
    bench_Zct.push_back((bench)*0.0004);
  }


 for(unsigned int ibench = 0; ibench<bench_Zut.size(); ibench++){
   for(unsigned int jbench = 0; jbench<bench_Zct.size(); jbench++){
     
     if(bench_Zut[ibench] == 0 & bench_Zct[jbench]==0 ) continue;
     doBDT_FCNC_tZ_scan ( "nom", bench_Zut[ibench], bench_Zct[jbench]);
   
   } 
 }
 
   
   
}



