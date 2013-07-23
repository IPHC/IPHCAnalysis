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


void doBDT_FCNC_tZ(TString thevertex, bool doScan,  double Vut=0, double Vct=0, TString syst = ""){

   
   
   
   
   
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

  TString outfileName( "trainingBDT_FCNC_" );
  outfileName+=thevertex+".root";
  TString sVut;
  TString sVct;
  if(doScan){
  
   std::ostringstream ossVut;
   ossVut << Vut ;
   sVut = ossVut.str();
  
   std::ostringstream ossVct;
   ossVct << Vct ;
   sVct = ossVct.str();
  
   outfileName =  "ScanRootFiles/ScanTrainingBDT_FCNC_"+thevertex+"_"+sVut+"_"+sVct+".root";
  }
  
  
  
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  
  
  
  
  
  TMVA::Factory *factory ;
  
  
  if(!doScan){ 
  
    cout << "no scan "<< endl;
    cout << "no scan "<< endl;
    cout << "no scan "<< endl;
    cout << "no scan "<< endl;
    if(thevertex == "zut")      factory = new TMVA::Factory( "BDT_trainning_zut", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    else if(thevertex == "zct") factory = new TMVA::Factory( "BDT_trainning_zct", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    else if(thevertex == "kut") factory = new TMVA::Factory( "BDT_trainning_kut", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    else if(thevertex == "kct") factory = new TMVA::Factory( "BDT_trainning_kct", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    else if(thevertex == "xut") factory = new TMVA::Factory( "BDT_trainning_xut", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    else if(thevertex == "xct") factory = new TMVA::Factory( "BDT_trainning_xct", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    //else if(thevertex == "tZq") factory = new TMVA::Factory( "BDT_trainning_tZq", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
  }else{
    if(thevertex == "zut")      factory = new TMVA::Factory( TString("BDT_trainning_zut")+"_"+sVut+"_"+sVct+"_"+syst, outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    else if(thevertex == "zct") factory = new TMVA::Factory( TString("BDT_trainning_zct")+"_"+sVut+"_"+sVct+"_"+syst, outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    else if(thevertex == "kut") factory = new TMVA::Factory( TString("BDT_trainning_kut")+"_"+sVut+"_"+sVct+"_"+syst, outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    else if(thevertex == "kct") factory = new TMVA::Factory( TString("BDT_trainning_kct")+"_"+sVut+"_"+sVct+"_"+syst, outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    else if(thevertex == "xut") factory = new TMVA::Factory( TString("BDT_trainning_xut")+"_"+sVut+"_"+sVct+"_"+syst, outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    else if(thevertex == "xct") factory = new TMVA::Factory( TString("BDT_trainning_xct")+"_"+sVut+"_"+sVct+"_"+syst, outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

  }

   
   /*TFile *input_sig1   = TFile::Open( "../../SystProofFiles_btadDiscri_update/proof_nom_btagcorr2.root" );
   TFile *input_wz   = TFile::Open( "../../SystProofFiles_btadDiscri_update/proof_nom_btagcorr2.root" );
   TFile *input_DY = TFile::Open( "../../SystProofFiles_btadDiscri_update/proof_nom_btagcorr2.root" );
 */
 
   
   TFile *input_sig1   = TFile::Open( "../../SystRootFiles_2013_07_01/proof_nom.root" );
   TFile *input_wz   = TFile::Open(   "../../SystRootFiles_2013_07_01/proof_nom.root" );
   TFile *input_DY = TFile::Open(     "../../SystRootFiles_2013_07_01/proof_DYenriched.root" );
 
 
 
   cout << "line 74 " <<  "FilesBenchmark/TreeMVA_benchmark_"+sVut+"_"+sVct+"_"+syst+"_"+"zqt.root"   << endl;
   
   TTree *signal1;
   if(thevertex == "zut") signal1 = (TTree*)input_sig1->Get("Ttree_FCNCzut");
   else if(thevertex == "zct") signal1 = (TTree*)input_sig1->Get("Ttree_FCNCzct");
   else if(thevertex == "kut") signal1 = (TTree*)input_sig1->Get("Ttree_FCNCkut");
   else if(thevertex == "kct") signal1 = (TTree*)input_sig1->Get("Ttree_FCNCkct");
   else if(thevertex == "xut") signal1 = (TTree*)input_sig1->Get("Ttree_FCNCxut");
   else if(thevertex == "xct") signal1 = (TTree*)input_sig1->Get("Ttree_FCNCxct");
   //else if(thevertex == "tZq") signal1 = (TTree*)input_sig1->Get("Ttree_TZq");
   if(doScan){
     if(thevertex == "zut") signal1 = (TTree*)input_sig1->Get( ("TreeMVA_benchmark_"+sVut+"_"+sVct+"_zqt").Data());
     if(thevertex == "kut") signal1 = (TTree*)input_sig1->Get( ("TreeMVA_benchmark_"+sVut+"_"+sVct+"_kqt").Data());
   }
   
  if ((TTree*)input_sig1->Get("Ttree_FCNCzut") == 0 )  cout << "no signal " << endl;;
  if(signal1 == 0) cout << "no signal " << endl;
  
  
   
  
  
  
  cout << " thevertex " << thevertex << endl;
  cout << " thevertex " << thevertex << endl;
  cout << " thevertex " << thevertex << endl;
  cout << " thevertex " << thevertex << endl;
  cout << " thevertex " << thevertex << endl;
  cout << " thevertex " << thevertex << endl;
  cout << " thevertex " << thevertex << endl;
   
   TTree *background_WZ     = (TTree*)input_wz->Get("Ttree_WZ");
   //TTree *background_DY     = (TTree*)input_wz->Get("Ttree_DataZJets");
   
   
   Double_t signalWeight1     = 1;
   Double_t signalWeight2     = 1;
   
   Double_t backgroundWeight1  = 1;
   Double_t backgroundWeight2  = 1;
   
   double sumWeight1 = 0.;
   double sumWeight2 = 0.;
  
    // You can add an arbitrary number of signal or background trees
   factory->AddSignalTree               ( signal1,            signalWeight1     );
   //if(doScan) factory->AddSignalTree    ( signal2,            signalWeight2     );
   factory->AddBackgroundTree( background_WZ,     backgroundWeight1 );
     
    
   cout << "line 169 " << endl;
   
   
   factory->AddVariable("tree_topMass",    'F'); 
   factory->AddVariable("tree_deltaPhilb", 'F');
   factory->AddVariable("tree_asym",       'F');
   factory->AddVariable("tree_Zpt",        'F');
   factory->AddVariable("tree_ZEta",       'F');
   ////factory->AddVariable("tree_topEta",	  'F');
   factory->AddVariable("tree_NJets",	    'F');
   factory->AddVariable("tree_NBJets",         'F');
   ////factory->AddVariable("tree_deltaRZl",     'F');   
   factory->AddVariable("tree_deltaPhiZmet",  'F');
   //if(thevertex != "zut" && thevertex != "kut" ) factory->AddVariable("tree_btagDiscri",    'F');  
   factory->AddVariable("tree_btagDiscri",    'F');  
   
	     	     
   factory->AddVariable("tree_leadJetEta",    'F');  	   
   factory->AddVariable("tree_deltaPhiZleptW",'F');  	
   
   
   
   //factory->AddVariable("tree_topPt",      'F');
   //factory->AddVariable("tree_leptWPt",       'F');  		   
   //factory->AddVariable("tree_leptWEta",      'F');  		   
   ////factory->AddVariable("tree_leadJetPt",     'F');    
   //factory->AddVariable("tree_deltaRZleptW",  'F');  
   //factory->AddVariable("tree_deltaRlb",   'F');
   //factory->AddVariable("tree_deltaRTopZ", 'F');
   //factory->AddVariable("tree_totMass",    'F');  	
   
        
   cout << "line 202 " << endl; 
   
   // to set weights. The variable must exist in the tree
   factory->SetSignalWeightExpression	 ("tree_EvtWeight");
   factory->SetBackgroundWeightExpression("tree_EvtWeight");
   
   
   // Apply additional cuts on the signal and background samples (can be different)
   //TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   //TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";
   TCut mycuts = "tree_mTW>20"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = "tree_mTW>20"; // for example: TCut mycutb = "abs(var1)<0.5";

   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   
   
   
   cout << "line 218" << endl; 
   //if(thevertex == "zut")      factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=60:nEventsMin=1200:MaxDepth=10:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );
   /*if(thevertex == "zut")      factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=60:nEventsMin=1000:MaxDepth=10:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );
   else if(thevertex == "zct") factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=500:MaxDepth=6:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );
   else if(thevertex == "kut") factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=60:nEventsMin=1200:MaxDepth=10:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );
   else if(thevertex == "kct") factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=500:MaxDepth=4:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );
*/
   if(thevertex == "zut")      factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );
   else if(thevertex == "zct") factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );
   else if(thevertex == "kut") factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );
   else if(thevertex == "kct") factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );



   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   cout << "line 238" << endl; 
   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch() && !doScan) TMVAGui( outfileName );


   cout << "signalWeight1 " << signalWeight1 << " sumWeight1 " << sumWeight1 << endl;
   cout << "signalWeight2 " << signalWeight2 << " sumWeight2 " << sumWeight2 << endl;
}



void trainingBDT_FCNC_tZ(){



   TString thevertex_zut = "zut";
   TString thevertex_zct = "zct";
   TString thevertex_kut = "kut";
   TString thevertex_kct = "kct";
  // TString thevertex_xut = "xut";
   //TString thevertex_xct = "xct";

   doBDT_FCNC_tZ (thevertex_zut, false);
   doBDT_FCNC_tZ (thevertex_zct, false);
   doBDT_FCNC_tZ (thevertex_kut, false);
   doBDT_FCNC_tZ (thevertex_kct, false);
   
   //doBDT_FCNC_tZ ("tZq", false);
   
   //doBDT_FCNC_tZ ( thevertex_zct,  false);
   //doBDT_FCNC_tZ ( thevertex_zut,  true, 0., 1., "nom");
   
   
   
   /*cout << "WARNING btag discri removed !!!!" << endl;
   cout << "WARNING btag discri removed !!!!" << endl;
   cout << "WARNING btag discri removed !!!!" << endl;
   cout << "WARNING btag discri removed !!!!" << endl;
   cout << "WARNING btag discri removed !!!!" << endl;
   cout << "WARNING btag discri removed !!!!" << endl;*/
   

}



