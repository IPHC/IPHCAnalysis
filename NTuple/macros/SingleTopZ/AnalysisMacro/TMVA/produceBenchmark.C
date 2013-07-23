


#include <cstdlib>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>


#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TRandom.h"



/*std::pair<double, double>  calculCrossSection(double Zut, double vct){
  
  std::pair<double, double> xspair;
  xspair.first = 2682000*Zut*Zut;
  xspair.second = 195000*vct*vct;
  
  return xspair;
}*/



void produceBenchmark(TTree* thevutTree, TTree* thevctTree, double vut, double vct, TString syst, TString sign){
  std::cout << "-----------------------------------" << vct << std::endl;
  
  std::cout << "producing becnhmark for vut = " << vut << "  , vct = " << vct << std::endl;
  
  TRandom r;
  r.SetSeed(int( (vut+vct)*10000)*123 );
  
  //std::pair<double, double>  theCrossSections = calculCrossSection(vut, vct);
  
  std::ostringstream ossvut;
  ossvut << vut ;
  TString svut = ossvut.str();
  
  std::ostringstream ossvct;
  ossvct << vct ;
  TString svct = ossvct.str();
  
  
  
   double sumWeight1 = 0.;
   double sumWeight2 = 0.;
   
   
   Float_t weight1, weight2;
   
   
   thevutTree->SetBranchAddress( "tree_EvtWeight",  &weight1); 
   thevctTree->SetBranchAddress( "tree_EvtWeight",  &weight2); 
   
   for(int ievent = 0; ievent < thevutTree->GetEntries(); ievent++){
      thevutTree->GetEntry(ievent);
      sumWeight1+=weight1;
   }
     
   for(int ievent = 0; ievent < thevctTree->GetEntries(); ievent++){
      thevctTree->GetEntry(ievent);
      sumWeight2+=weight2;
   }
     
  
   double integral = vut*sumWeight1+vct*sumWeight2;
   
   
   double signalWeight1  = vut*sumWeight1/integral;
   double signalWeight2  = vct*sumWeight2/integral;
  
  cout << "signalWeight1 " << signalWeight1 << endl;
  cout << "signalWeight2 " << signalWeight2 << endl;
  
  
  
  //-----------------------------------------------------------------------
  //defined the variables used by the BDT reader
  // tree to fill 
  //-----------------------------------------------------------------------
  
  TFile * output = new TFile( ("FilesBenchmark/TreeMVA_benchmark_"+svut+"_"+svct+"_"+syst+"_"+sign+".root").Data(), "recreate");
  
  output->cd();
  TTree *outputTree = new TTree(("TreeMVA_benchmark_"+svut+"_"+svct+"_"+sign).Data(), ("TreeMVA_benchmark_"+svut+"_"+svct+"_"+sign).Data());
  
  Float_t topMass,deltaPhilb,asym,Zpt, deltaRlb, deltaRTopZ, ZEta, topPt, topEta, totPt, totEta;
  Float_t totMass, deltaRZl, deltaPhiZmet, btagDiscri, NJets, NBJets, EvtWeight;
  Float_t leptWPt,leptWEta,leadJetPt,leadJetEta,deltaRZleptW,deltaPhiZleptW;
  
  outputTree->Branch("tree_topMass",     &topMass,     "tree_topMass/F"   );
  outputTree->Branch("tree_totMass",     &totMass,     "tree_totMass/F"   );
  outputTree->Branch("tree_deltaPhilb",  &deltaPhilb,  "tree_deltaPhilb/F");
  outputTree->Branch("tree_deltaRlb",    &deltaRlb,    "tree_deltaRlb/F"  );
  outputTree->Branch("tree_deltaRTopZ",  &deltaRTopZ,  "tree_deltaRTopZ/F");
  outputTree->Branch("tree_asym",        &asym,        "tree_asym/F"	  );
  outputTree->Branch("tree_Zpt",         &Zpt,         "tree_Zpt/F"	  );
  outputTree->Branch("tree_ZEta",        &ZEta,        "tree_ZEta/F"	  );
  outputTree->Branch("tree_topPt",       &topPt,       "tree_topPt/F"	  );
  outputTree->Branch("tree_topEta",      &topEta,      "tree_topEta/F"    );
  outputTree->Branch("tree_NJets",       &NJets,       "tree_NJets/F"	  );
  outputTree->Branch("tree_NBJets",      &NBJets,      "tree_NBJets/F"    );
  outputTree->Branch("tree_deltaRZl",    &deltaRZl,    "tree_deltaRZl/F"     );
  outputTree->Branch("tree_deltaPhiZmet",&deltaPhiZmet,"tree_deltaPhiZmet/F" );
  outputTree->Branch("tree_btagDiscri",  &btagDiscri,  "tree_btagDiscri/F"   );
  
  outputTree->Branch("tree_totPt",      &totPt,      "tree_totPt/F"   );
  outputTree->Branch("tree_totEta",     &totEta,     "tree_totEta/F"   );
  
  
  outputTree->Branch("tree_leptWPt",        &leptWPt	    , "tree_leptWPt/F"        );
  outputTree->Branch("tree_leptWEta",       &leptWEta	    , "tree_leptWEta/F"       );
  outputTree->Branch("tree_leadJetPt",      &leadJetPt      , "tree_leadJetPt/F"      ); 
  outputTree->Branch("tree_leadJetEta",     &leadJetEta     , "tree_leadJetEta/F"     );
  outputTree->Branch("tree_deltaRZleptW",   &deltaRZleptW   , "tree_deltaRZleptW/F"   );
  outputTree->Branch("tree_deltaPhiZleptW", &deltaPhiZleptW , "tree_deltaPhiZleptW/F" );
  
  
  
  outputTree->Branch("tree_EvtWeight",   &EvtWeight,   "tree_EvtWeight/F" );
  
  
  //-----------------------------------------------------------------------
  //defined the variables used by the BDT reader
  // tree to read 
  //-----------------------------------------------------------------------
  Float_t topMass_tmp,deltaPhilb_tmp,asym_tmp,Zpt_tmp, deltaRlb_tmp, deltaRTopZ_tmp, ZEta_tmp, topPt_tmp, topEta_tmp;
  Float_t totMass_tmp, deltaRZl_tmp, deltaPhiZmet_tmp, btagDiscri_tmp, NJets_tmp, NBJets_tmp, EvtWeight_tmp;
  Float_t leptWPt_tmp,leptWEta_tmp,leadJetPt_tmp,leadJetEta_tmp,deltaRZleptW_tmp,deltaPhiZleptW_tmp, totPt_tmp, totEta_tmp;

  
  
  thevutTree->SetBranchAddress( "tree_topMass",	       &topMass_tmp);
  thevutTree->SetBranchAddress( "tree_totMass",	       &totMass_tmp );
  thevutTree->SetBranchAddress( "tree_deltaPhilb",     &deltaPhilb_tmp );
  thevutTree->SetBranchAddress( "tree_deltaRlb",       &deltaRlb_tmp);
  thevutTree->SetBranchAddress( "tree_deltaRTopZ",     &deltaRTopZ_tmp);
  thevutTree->SetBranchAddress( "tree_asym",	       &asym_tmp );
  thevutTree->SetBranchAddress( "tree_Zpt", 	       &Zpt_tmp );
  thevutTree->SetBranchAddress( "tree_ZEta",	       &ZEta_tmp);
  thevutTree->SetBranchAddress( "tree_topPt",	       &topPt_tmp);
  thevutTree->SetBranchAddress( "tree_topEta",	       &topEta_tmp);
  thevutTree->SetBranchAddress( "tree_NJets",	       &NJets_tmp);
  thevutTree->SetBranchAddress( "tree_NBJets",	       &NBJets_tmp);	      
  thevutTree->SetBranchAddress( "tree_deltaRZl",       &deltaRZl_tmp);   
  thevutTree->SetBranchAddress( "tree_deltaPhiZmet",   &deltaPhiZmet_tmp);
  thevutTree->SetBranchAddress( "tree_btagDiscri",     &btagDiscri_tmp);
  
  thevutTree->SetBranchAddress( "tree_leptWPt",	       &leptWPt_tmp);	       
  thevutTree->SetBranchAddress( "tree_leptWEta",       &leptWEta_tmp);	       
  thevutTree->SetBranchAddress( "tree_leadJetPt",      &leadJetPt_tmp);		 
  thevutTree->SetBranchAddress( "tree_leadJetEta",     &leadJetEta_tmp); 	
  thevutTree->SetBranchAddress( "tree_deltaPhiZleptW", &deltaPhiZleptW_tmp);
  thevutTree->SetBranchAddress( "tree_EvtWeight",  &EvtWeight_tmp);
  
  
  
   
  for (Long64_t ievt=0; ievt<thevutTree->GetEntries();ievt++) { 
    
    if( r.Uniform() >  signalWeight1 ) continue;
    //cout << "in tree fill vut" << endl;
    thevutTree->GetEntry(ievt); 
    
    topMass        = topMass_tmp;
    deltaPhilb     = deltaPhilb_tmp;
    asym           = asym_tmp;
    Zpt            = Zpt_tmp; 
    deltaRlb       = deltaRlb_tmp;
    deltaRTopZ     = deltaRTopZ_tmp;
    ZEta           = ZEta_tmp; 
    topPt          = topPt_tmp;
    topEta         = topEta_tmp;
    totMass        = totMass_tmp;
    deltaRZl       = deltaRZl_tmp; 
    deltaPhiZmet   = deltaPhiZmet_tmp; 
    btagDiscri     = btagDiscri_tmp; 
    NJets          = NJets_tmp; 
    NBJets         = NBJets_tmp; 
    EvtWeight      = EvtWeight_tmp;
    leptWPt        = leptWPt_tmp;
    leptWEta       = leptWEta_tmp;
    leadJetPt      = leadJetPt_tmp;
    leadJetEta     = leadJetEta_tmp;
    deltaRZleptW   = deltaRZleptW_tmp;
    deltaPhiZleptW = deltaPhiZleptW_tmp;
    EvtWeight      = EvtWeight_tmp ;
    
    outputTree->Fill();
  }
  
  
  
  
  //-----------------------------------------------------------------------
  //defined the variables used by the BDT reader
  // tree to read 
  //-----------------------------------------------------------------------
   
  thevctTree->SetBranchAddress( "tree_topMass",	      &topMass_tmp);
  thevctTree->SetBranchAddress( "tree_totMass",	      &totMass_tmp );
  thevctTree->SetBranchAddress( "tree_deltaPhilb",    &deltaPhilb_tmp );
  thevctTree->SetBranchAddress( "tree_deltaRlb",      &deltaRlb_tmp);
  thevctTree->SetBranchAddress( "tree_deltaRTopZ",    &deltaRTopZ_tmp);
  thevctTree->SetBranchAddress( "tree_asym",	      &asym_tmp );
  thevctTree->SetBranchAddress( "tree_Zpt", 	      &Zpt_tmp );
  thevctTree->SetBranchAddress( "tree_ZEta",	      &ZEta_tmp);
  thevctTree->SetBranchAddress( "tree_topPt",	      &topPt_tmp);
  thevctTree->SetBranchAddress( "tree_topEta",	      &topEta_tmp);
  thevctTree->SetBranchAddress( "tree_NJets",	      &NJets_tmp);
  thevctTree->SetBranchAddress( "tree_NBJets",	      &NBJets_tmp);		
  thevctTree->SetBranchAddress( "tree_deltaRZl",      &deltaRZl_tmp);   
  thevctTree->SetBranchAddress( "tree_deltaPhiZmet",  &deltaPhiZmet_tmp);
  thevctTree->SetBranchAddress( "tree_btagDiscri",    &btagDiscri_tmp);
  
  thevctTree->SetBranchAddress("tree_leptWPt",	      &leptWPt_tmp);	       
  thevctTree->SetBranchAddress("tree_leptWEta",	      &leptWEta_tmp);	       
  thevctTree->SetBranchAddress("tree_leadJetPt",      &leadJetPt_tmp);		 
  thevctTree->SetBranchAddress("tree_leadJetEta",     &leadJetEta_tmp); 	
  thevctTree->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW_tmp);
  thevctTree->SetBranchAddress("tree_EvtWeight",      &EvtWeight_tmp);
  
  
  
  
  for (Long64_t ievt=0; ievt<thevctTree->GetEntries();ievt++) { 
    
    
    if( r.Uniform() <  signalWeight1 ) continue;
    
    //cout << "in tree fill vct" << endl;
    thevctTree->GetEntry(ievt);
    
    topMass        = topMass_tmp;
    deltaPhilb     = deltaPhilb_tmp;
    asym           = asym_tmp;
    Zpt            = Zpt_tmp; 
    deltaRlb       = deltaRlb_tmp;
    deltaRTopZ     = deltaRTopZ_tmp;
    ZEta           = ZEta_tmp; 
    topPt          = topPt_tmp;
    topEta         = topEta_tmp;
    totMass        = totMass_tmp;
    deltaRZl       = deltaRZl_tmp; 
    deltaPhiZmet   = deltaPhiZmet_tmp; 
    btagDiscri     = btagDiscri_tmp; 
    NJets          = NJets_tmp; 
    NBJets         = NBJets_tmp; 
    EvtWeight      = EvtWeight_tmp;
    leptWPt        = leptWPt_tmp;
    leptWEta       = leptWEta_tmp;
    leadJetPt      = leadJetPt_tmp;
    leadJetEta     = leadJetEta_tmp;
    deltaRZleptW   = deltaRZleptW_tmp;
    deltaPhiZleptW = deltaPhiZleptW_tmp;
    EvtWeight      = EvtWeight_tmp ;
    outputTree->Fill();
  }
  
  
  
  
  outputTree->Write();
  
  std::cout << "production done "<< std::endl;
  
  
  std::cout << "created file " << ("FilesBenchmark/TreeMVA_benchmark_"+svut+"_"+svct+"_"+syst+".root").Data() << std::endl;
  std::cout << "containing the tree " << ("TreeMVA_benchmark_"+svut+"_"+svct).Data() << std::endl;
  delete outputTree;
  delete  output;
}





void produceBenchmark(TFile *input_FCNC, TString syst, TString sign){

  
 TTree* theTree_vut;
 TTree* theTree_vct;
 
 if(sign == "zqt"){
   theTree_vut = (TTree*)input_FCNC->Get("Ttree_FCNCzut");
   theTree_vct = (TTree*)input_FCNC->Get("Ttree_FCNCzct");
 }else{
   if(sign == "kqt"){
   theTree_vut = (TTree*)input_FCNC->Get("Ttree_FCNCkut");
   theTree_vct = (TTree*)input_FCNC->Get("Ttree_FCNCkct");
   }
 }
 
 
 std::vector<double> bench_vut, bench_vct;
 
 for(int bench=0; bench<6;   bench++){
   bench_vut.push_back((bench)*0.00004);
 }
 for(int bench=0; bench<6; bench++){
   bench_vct.push_back((bench)*0.0004);
 }
 
 for(unsigned int ibench = 0; ibench<bench_vut.size(); ibench++){
   for(unsigned int jbench = 0; jbench<bench_vct.size(); jbench++){
     
     if(bench_vut[ibench] == 0 & bench_vct[jbench]==0 ) continue;
     
     produceBenchmark(theTree_vut, theTree_vct, bench_vut[ibench], bench_vct[jbench], syst, sign); 
   } 
 }
 
 //produceBenchmark(theTree_vut, theTree_vct, 0.0, 1.0, syst, sign); 
  
 
 
 
}



void produceBenchmark(){

 TFile *input_FCNC    = new TFile("../../RootFiles/proof_woWZSF.root", "read");
 produceBenchmark(input_FCNC, "nom", "zqt");
 produceBenchmark(input_FCNC, "nom", "kqt");
 
 
}



