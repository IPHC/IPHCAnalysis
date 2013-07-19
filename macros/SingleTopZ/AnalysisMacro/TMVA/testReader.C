


#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "test/TMVAGui.C"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"



void testReader(){
   // This loads the library
   TMVA::Tools::Instance();

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    
   
   
   
   Float_t topMass,deltaPhilb,asym,Zpt ;
   reader->AddVariable("tree_topMass",    &topMass    );
   reader->AddVariable("tree_deltaPhilb", &deltaPhilb );
   reader->AddVariable("tree_asym",	  &asym       );
   reader->AddVariable("tree_Zpt",	  &Zpt        );
   
   
   
   UInt_t nbin = 40;
   TH1F* histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.18, 0.25 );
 
   reader->BookMVA( "BDT", "weights/doBDT_FCNC_tZ_BDT.weights.xml" ); 
   
   TFile *input = new TFile("proof.root", "read");
   
   TTree* theTree_DataMu = (TTree*)input->Get("Ttree_DataMu");
   theTree_DataMu->SetBranchAddress( "tree_topMass",    &topMass );
   theTree_DataMu->SetBranchAddress( "tree_deltaPhilb", &deltaPhilb );
   theTree_DataMu->SetBranchAddress( "tree_asym",       &asym );
   theTree_DataMu->SetBranchAddress( "tree_Zpt",        &Zpt );
   
   
   for (Long64_t ievt=0; ievt<theTree_DataMu->GetEntries();ievt++) {    
     if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DataMu->GetEntry(ievt);
     histBdt->Fill( reader->EvaluateMVA( "BDT"           ) );
   }
   
   
   
   TTree* theTree_DataEG = (TTree*)input->Get("Ttree_DataEG");
   theTree_DataEG->SetBranchAddress( "tree_topMass",    &topMass );
   theTree_DataEG->SetBranchAddress( "tree_deltaPhilb", &deltaPhilb );
   theTree_DataEG->SetBranchAddress( "tree_asym",       &asym );
   theTree_DataEG->SetBranchAddress( "tree_Zpt",        &Zpt );
   
   
   for (Long64_t ievt=0; ievt<theTree_DataEG->GetEntries();ievt++) {
     if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DataEG->GetEntry(ievt);
     histBdt    ->Fill( reader->EvaluateMVA( "BDT"           ) );
   }
   
   
   
   /*
   TTree* theTree_DataMuEG = (TTree*)input->Get("Ttree_DataMuEG");
   theTree_DataMuEG->SetBranchAddress( "tree_topMass",    &topMass );
   theTree_DataMuEG->SetBranchAddress( "tree_deltaPhilb", &deltaPhilb );
   theTree_DataMuEG->SetBranchAddress( "tree_asym",       &asym );
   theTree_DataMuEG->SetBranchAddress( "tree_Zpt",        &Zpt );
   
   
   for (Long64_t ievt=0; ievt<theTree_DataMuEG->GetEntries();ievt++) {  
     if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DataMuEG->GetEntry(ievt);
     histBdt    ->Fill( reader->EvaluateMVA( "BDT"           ) );
   }*/
   
   
   
   
   TFile *target  = new TFile( "TMVApp.root","RECREATE" );
   target->cd();
   histBdt->Write();
   target->Close();
   
}
