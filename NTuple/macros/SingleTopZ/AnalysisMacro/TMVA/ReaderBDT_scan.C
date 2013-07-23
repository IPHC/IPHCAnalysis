


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





std::map<TString,std::vector<TH1F*> > theHistoMap;

std::pair<double, double>  calculCrossSection(double Zut, double Zct){
  
  std::pair<double, double> xspair;
  xspair.first = 2682000*Zut*Zut;
  xspair.second = 195000*Zct*Zct;
  
  return xspair;
}

void fillHisto(TString sample, std::vector<double> theVar, int ievt, double weight=1){
   
   
   
   std::vector<TH1F*> histovect = theHistoMap[sample];
   
   histovect[0]->Fill(theVar[0],         weight); //topMass
   histovect[1]->Fill(theVar[1],         weight); //totMass
   histovect[2]->Fill(fabs(theVar[2]),   weight); //deltaPhilb
   histovect[3]->Fill(theVar[3],         weight); //deltaRlb
   histovect[4]->Fill(theVar[4],         weight); //deltaRTopZ
   histovect[5]->Fill(theVar[5],         weight); //asym
   histovect[6]->Fill(theVar[6],         weight); //Zpt
   histovect[7]->Fill(theVar[7],         weight); //ZEta
   histovect[8]->Fill(theVar[8],         weight); //topPt
   histovect[9]->Fill(theVar[9],         weight); //topEta
   histovect[10]->Fill(theVar[10],       weight); //NJets
   histovect[11]->Fill(theVar[11],       weight); //NBJets
   histovect[12]->Fill(theVar[12],       weight); //deltaRZl
   histovect[13]->Fill(fabs(theVar[13]), weight); //deltaPhiZmet
   histovect[14]->Fill(theVar[14],       weight); //btagDiscri
   histovect[15]->Fill(theVar[15],       weight); //leptWPt
   histovect[16]->Fill(theVar[16],       weight); //leptWEta
   histovect[17]->Fill(theVar[17],       weight); //leadJetPt
   histovect[18]->Fill(theVar[18],       weight); //leadJetEta
   histovect[19]->Fill(fabs(theVar[19]), weight); //deltaPhiZleptW
   

 
}


void writeHisto(TString sample){

   std::vector<TH1F*> histovect = theHistoMap[sample];
   
   histovect[0]->Write(); //topMass
   histovect[1]->Write(); //totMass
   histovect[2]->Write(); //deltaPhilb
   histovect[3]->Write(); //deltaRlb
   histovect[4]->Write(); //deltaRTopZ
   histovect[5]->Write(); //asym
   histovect[6]->Write(); //Zpt
   histovect[7]->Write(); //ZEta
   histovect[8]->Write(); //topPt
   histovect[9]->Write(); //topEta
   histovect[10]->Write(); //NJets
   histovect[11]->Write(); //NBJets
   histovect[12]->Write(); //deltaRZl
   histovect[13]->Write(); //deltaPhiZmet
   histovect[14]->Write(); //btagDiscri
   histovect[15]->Write(); //leptWPt
   histovect[16]->Write(); //leptWEta
   histovect[17]->Write(); //leadJetPt
   histovect[18]->Write(); //leadJetEta
   histovect[19]->Write(); //deltaPhiZleptW

}


void scaleHisto(TString sample, double thescale){

   std::vector<TH1F*> histovect = theHistoMap[sample];
   
   histovect[0]->Scale(thescale); //topMass
   histovect[1]->Scale(thescale); //totMass
   histovect[2]->Scale(thescale); //deltaPhilb
   histovect[3]->Scale(thescale); //deltaRlb
   histovect[4]->Scale(thescale); //deltaRTopZ
   histovect[5]->Scale(thescale); //asym
   histovect[6]->Scale(thescale); //Zpt
   histovect[7]->Scale(thescale); //ZEta
   histovect[8]->Scale(thescale); //topPt
   histovect[9]->Scale(thescale); //topEta
   histovect[10]->Scale(thescale); //NJets
   histovect[11]->Scale(thescale); //NBJets
   histovect[12]->Scale(thescale); //deltaRZl
   histovect[13]->Scale(thescale); //deltaPhiZmet
   histovect[14]->Scale(thescale); //btagDiscri
   histovect[15]->Scale(thescale); //leptWPt
   histovect[16]->Scale(thescale); //leptWEta
   histovect[17]->Scale(thescale); //leadJetPt
   histovect[18]->Scale(thescale); //leadJetEta
   histovect[19]->Scale(thescale); //deltaPhiZleptW

}


void ReaderBDT_scan(TString stringinput, TString syst, double Zut, double Zct){


   
  //double WZscale = 0.91;
  double WZscale = 1.;
  float bdt_BCK_cut = -0.20;
  r.SetSeed(int( (Zut+Zct)*10000)*123 );
  
  
  std::ostringstream ossZut;
  ossZut << Zut ;
  TString sZut = ossZut.str();
  
  std::ostringstream ossZct;
  ossZct << Zct ;
  TString sZct = ossZct.str();

   
   
   // This loads the library
   TMVA::Tools::Instance();
   
   //create the BDT reader
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );      
   
   TFile *input         = new TFile(stringinput, "read");
   TFile *input_FCNC    = new TFile( ("FilesBenchmark/TreeMVA_benchmark_"+sZut+"_"+sZct+"_"+syst+".root").Data(), "read");
   
   
   TFile *input_zjets   = new TFile("../RootFiles/backup_outputProof20-12-12_10-03-23_noSFapplied_Zenriched/proof_merged_noSF_Zenriched.root", "read");
   
   
   //*************************************************
   //*************************************************
   //declaration of histograms
   //*************************************************
   //*************************************************
   
   
   TFile *target  ;
   
   if(thevertex == "zut") target= new TFile( ("HistoBDToutput/TMVApp_zut_"+syst+".root").Data(),"RECREATE" );
   if(thevertex == "zct") target= new TFile( ("HistoBDToutput/TMVApp_zct_"+syst+".root").Data(),"RECREATE" );
   if(thevertex == "kut") target= new TFile( ("HistoBDToutput/TMVApp_kut_"+syst+".root").Data(),"RECREATE" );
   if(thevertex == "kct") target= new TFile( ("HistoBDToutput/TMVApp_kct_"+syst+".root").Data(),"RECREATE" );
   if(thevertex == "xut") target= new TFile( ("HistoBDToutput/TMVApp_xut_"+syst+".root").Data(),"RECREATE" );
   if(thevertex == "xct") target= new TFile( ("HistoBDToutput/TMVApp_xct_"+syst+".root").Data(),"RECREATE" );
   
   
   std::vector<TString> samples;
   
   samples.push_back("Data");	 
   samples.push_back("Zjets");	 
   samples.push_back("WZ");	 
   samples.push_back("DataZjets");
   samples.push_back("TTbarSig"); 
   samples.push_back("Wjets");	 
   samples.push_back("TtW");	 
   samples.push_back("TbartW");   
   samples.push_back("TtChan");   
   samples.push_back("TbartChan");
   samples.push_back("TsChan");   
   samples.push_back("TbarsChan");
   samples.push_back("ZZ");	 
   samples.push_back("WW");	 
   samples.push_back("Signal");	 
   
   
   
   std::vector<TString> variables;
   

   variables.push_back("topMass");
   variables.push_back("totMass");
   variables.push_back("deltaPhilb");
   variables.push_back("deltaRlb");
   variables.push_back("deltaRTopZ");
   variables.push_back("asym");
   variables.push_back("Zpt");
   variables.push_back("ZEta");
   variables.push_back("topPt");
   variables.push_back("topEta");
   variables.push_back("NJets");   
   variables.push_back("NBJets");  
   variables.push_back("deltaRZl"); 
   variables.push_back("deltaPhiZmet");
   variables.push_back("btagDiscri");
   variables.push_back("leptWPt");
   variables.push_back("leptWEta");	   
   variables.push_back("leadJetPt");	     
   variables.push_back("leadJetEta");	      
   variables.push_back("deltaPhiZleptW");
   
   
   for(unsigned int i=0; i<samples.size(); i++){
     int nbins = 0; double xmin = 0; double xmax=0;
     
     std::vector<TH1F*> histovect;
     
     for(unsigned int j=0; j<variables.size(); j++){
     
       if(variables[j]=="topMass")       {nbins = 20; xmin = 0;    xmax = 500;};
       if(variables[j]=="totMass")       {nbins = 20; xmin = 100;  xmax = 700;};
       if(variables[j]=="deltaPhilb")    {nbins = 10; xmin = 0;    xmax = 3.2;};
       if(variables[j]=="deltaRlb")      {nbins = 10; xmin = 0;    xmax = 3.2;};
       if(variables[j]=="deltaRTopZ")    {nbins = 10; xmin = 0;    xmax = 3.2;};
       if(variables[j]=="asym")          {nbins = 10; xmin = -2.5; xmax = 2.5;};
       if(variables[j]=="Zpt")           {nbins = 20; xmin = 0;    xmax = 500;};
       if(variables[j]=="ZEta")          {nbins = 10; xmin = -2.5; xmax = 2.5;};
       if(variables[j]=="topPt")         {nbins = 20; xmin = 0;    xmax = 500;};
       if(variables[j]=="topEta")        {nbins = 10; xmin = -2.5; xmax = 2.5;};
       if(variables[j]=="NJets")         {nbins = 4;  xmin = 0.5;  xmax = 4.5;};   
       if(variables[j]=="NBJets")        {nbins = 2;  xmin = -0.5; xmax = 1.5;};  
       if(variables[j]=="deltaRZl")      {nbins = 10; xmin = 0;    xmax = 3.2;}; 
       if(variables[j]=="deltaPhiZmet")  {nbins = 10; xmin = 0;    xmax = 3.2;};
       if(variables[j]=="btagDiscri")    {nbins = 20; xmin = -1;   xmax = 1;};
       if(variables[j]=="leptWPt")       {nbins = 20; xmin = 0;    xmax = 350;};
       if(variables[j]=="leptWEta")      {nbins = 10; xmin = -2.5; xmax = 2.5;};        
       if(variables[j]=="leadJetPt")     {nbins = 20; xmin = 0;    xmax = 500;};	 
       if(variables[j]=="leadJetEta")    {nbins = 20; xmin = -2.5; xmax = 2.5;};	  
       if(variables[j]=="deltaPhiZleptW"){nbins = 10; xmin = 0;    xmax = 3.2;};
       
       //cout << "(variables[j]) "  << variables[j] << " nbins " << nbins << " xmin " << xmin << " xmax " << xmax << endl;
       
       TH1F * histo = new TH1F((variables[j]+"_"+samples[i]).Data(), (variables[j]+"_"+samples[i]).Data(), nbins, xmin, xmax);
       histovect.push_back(histo);
     }
     
     theHistoMap[samples[i]] = histovect;
   }
   
   
   
   
   
   
   //defined the variables used by the BDT reader
   Float_t topMass,deltaPhilb,asym,Zpt, deltaRlb, deltaRTopZ, ZEta, topPt, topEta;
   Float_t totMass, deltaRZl, deltaPhiZmet, btagDiscri, NJets, NBJets, EvtWeight;
   Float_t leptWPt,leptWEta,leadJetPt,leadJetEta,deltaRZleptW,deltaPhiZleptW;
   
   
   Int_t Channel;

   reader->AddVariable("tree_topMass",    &topMass    );
   reader->AddVariable("tree_totMass",    &totMass    );
   reader->AddVariable("tree_deltaPhilb", &deltaPhilb );
   reader->AddVariable("tree_deltaRlb",     &deltaRlb);
   reader->AddVariable("tree_deltaRTopZ",   &deltaRTopZ);
   reader->AddVariable("tree_asym",	  &asym       );
   reader->AddVariable("tree_Zpt",	  &Zpt        );
   
   reader->AddVariable("tree_ZEta",	     &ZEta);
   reader->AddVariable("tree_topPt",         &topPt);
   reader->AddVariable("tree_topEta",        &topEta); 
   reader->AddVariable("tree_NJets",         &NJets);	  
   reader->AddVariable("tree_NBJets",        &NBJets);	    
   reader->AddVariable("tree_deltaRZl",      &deltaRZl);   
   reader->AddVariable("tree_deltaPhiZmet",  &deltaPhiZmet);
   reader->AddVariable("tree_btagDiscri",    &btagDiscri);     
    		    
   reader->AddVariable("tree_leptWPt",        &leptWPt);	    
   reader->AddVariable("tree_leptWEta",       &leptWEta);	    
   reader->AddVariable("tree_leadJetPt",      &leadJetPt);	      
   reader->AddVariable("tree_leadJetEta",     &leadJetEta);	     
   //reader->AddVariable("tree_deltaRZleptW",   &deltaRZleptW);	   
   reader->AddVariable("tree_deltaPhiZleptW", &deltaPhiZleptW);  
   
   
   UInt_t nbin = 40;
   TH1F* histBdt_Data           = new TH1F( "MVA_BDT_Data",      "MVA_BDT_Data",       nbin, -0.5, 0.8 );
   TH1F* histBdt_Zjets          = new TH1F( "MVA_BDT_Zjets",     "MVA_BDT_Zjets",      nbin, -0.5, 0.8 );
   TH1F* histBdt_WZ             = new TH1F( "MVA_BDT_WZ",        "MVA_BDT_WZ",	       nbin, -0.5, 0.8 );
   TH1F* histBdt_DataZjets      = new TH1F( "MVA_BDT_DataZjets", "MVA_BDT_DataZjets",  nbin, -0.5, 0.8 );
   
   TH1F* histBdt_TTbarSig       = new TH1F( "MVA_BDT_TTbarSig",	 "MVA_BDT_TTbarSig",   nbin, -0.5, 0.8 );
   TH1F* histBdt_Wjets          = new TH1F( "MVA_BDT_Wjets",	 "MVA_BDT_Wjets",      nbin, -0.5, 0.8 );
   TH1F* histBdt_TtW            = new TH1F( "MVA_BDT_TtW",	 "MVA_BDT_TtW",        nbin, -0.5, 0.8 );
   TH1F* histBdt_TbartW         = new TH1F( "MVA_BDT_TbartW",	 "MVA_BDT_TbartW",     nbin, -0.5, 0.8 );
   TH1F* histBdt_TtChan         = new TH1F( "MVA_BDT_TtChan",	 "MVA_BDT_TtChan",     nbin, -0.5, 0.8 );
   TH1F* histBdt_TbartChan      = new TH1F( "MVA_BDT_TbartChan", "MVA_BDT_TbartChan",  nbin, -0.5, 0.8 );
   TH1F* histBdt_TsChan         = new TH1F( "MVA_BDT_TsChan",	 "MVA_BDT_TsChan",     nbin, -0.5, 0.8 );
   TH1F* histBdt_TbarsChan      = new TH1F( "MVA_BDT_TbarsChan", "MVA_BDT_TbarsChan",  nbin, -0.5, 0.8 );
   TH1F* histBdt_ZZ             = new TH1F( "MVA_BDT_ZZ",	 "MVA_BDT_ZZ",         nbin, -0.5, 0.8 );
   TH1F* histBdt_WW             = new TH1F( "MVA_BDT_WW",	 "MVA_BDT_WW",         nbin, -0.5, 0.8 );
   
   TH1F* histBdt_FCNC;
   
   if(thevertex == "zut") histBdt_FCNC = new TH1F( "MVA_BDT_FCNC_zut",   "MVA_BDT_FCNC_zut", nbin, -0.5, 0.8 );
   if(thevertex == "zct") histBdt_FCNC = new TH1F( "MVA_BDT_FCNC_zct",   "MVA_BDT_FCNC_zct", nbin, -0.5, 0.8 );
   if(thevertex == "kut") histBdt_FCNC = new TH1F( "MVA_BDT_FCNC_kut",   "MVA_BDT_FCNC_kut", nbin, -0.5, 0.8 );
   if(thevertex == "kct") histBdt_FCNC = new TH1F( "MVA_BDT_FCNC_kct",   "MVA_BDT_FCNC_kct", nbin, -0.5, 0.8 );
   if(thevertex == "xut") histBdt_FCNC = new TH1F( "MVA_BDT_FCNC_xut",   "MVA_BDT_FCNC_xut", nbin, -0.5, 0.8 );
   if(thevertex == "xct") histBdt_FCNC = new TH1F( "MVA_BDT_FCNC_xct",   "MVA_BDT_FCNC_xct", nbin, -0.5, 0.8 );
   
   //Get the BDT trainning through an xml file, created during the trainning phase
   if(thevertex == "zut") reader->BookMVA( "BDT", "weights/BDT_trainning_zut_BDT.weights.xml" ); 
   if(thevertex == "zct") reader->BookMVA( "BDT", "weights/BDT_trainning_zct_BDT.weights.xml" ); 
   if(thevertex == "kut") reader->BookMVA( "BDT", "weights/BDT_trainning_kut_BDT.weights.xml" ); 
   if(thevertex == "kct") reader->BookMVA( "BDT", "weights/BDT_trainning_kct_BDT.weights.xml" ); 
   if(thevertex == "xut") reader->BookMVA( "BDT", "weights/BDT_trainning_xut_BDT.weights.xml" ); 
   if(thevertex == "xct") reader->BookMVA( "BDT", "weights/BDT_trainning_xct_BDT.weights.xml" ); 
   
   float ndataEvents_BCKenriched = 0;
   
   //-----------------------------------------------------
   //for Data
   //-----------------------------------------------------
   input->cd();
   //define the tree to read
   TTree* theTree_DataMu = (TTree*)input->Get("Ttree_DataMu");
   theTree_DataMu->SetBranchAddress( "tree_topMass",	  &topMass);
   theTree_DataMu->SetBranchAddress( "tree_totMass",      &totMass );
   theTree_DataMu->SetBranchAddress( "tree_deltaPhilb",   &deltaPhilb );
   theTree_DataMu->SetBranchAddress( "tree_deltaRlb",	  &deltaRlb);
   theTree_DataMu->SetBranchAddress( "tree_deltaRTopZ",   &deltaRTopZ);
   theTree_DataMu->SetBranchAddress( "tree_asym",         &asym );
   theTree_DataMu->SetBranchAddress( "tree_Zpt",          &Zpt );
   theTree_DataMu->SetBranchAddress( "tree_ZEta",	  &ZEta);
   theTree_DataMu->SetBranchAddress( "tree_topPt",	  &topPt);
   theTree_DataMu->SetBranchAddress( "tree_topEta",	  &topEta);
   theTree_DataMu->SetBranchAddress( "tree_NJets",	  &NJets);
   theTree_DataMu->SetBranchAddress( "tree_NBJets",	  &NBJets);		 
   theTree_DataMu->SetBranchAddress( "tree_deltaRZl",	  &deltaRZl);   
   theTree_DataMu->SetBranchAddress( "tree_deltaPhiZmet", &deltaPhiZmet);
   theTree_DataMu->SetBranchAddress( "tree_btagDiscri",   &btagDiscri);
   
   theTree_DataMu->SetBranchAddress("tree_leptWPt",	  &leptWPt);		
   theTree_DataMu->SetBranchAddress("tree_leptWEta",	  &leptWEta);		
   theTree_DataMu->SetBranchAddress("tree_leadJetPt",	  &leadJetPt);  	  
   theTree_DataMu->SetBranchAddress("tree_leadJetEta",     &leadJetEta);	 
   //theTree_DataMu->SetBranchAddress("tree_deltaRZleptW",   &deltaRZleptW);     
   theTree_DataMu->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);	
    	
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DataMu->GetEntries();ievt++) {    
     if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     
     theTree_DataMu->GetEntry(ievt);
     histBdt_Data->Fill( reader->EvaluateMVA( "BDT"           ) );
     
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     fillHisto("Data", theVar, ievt);
     
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut) ndataEvents_BCKenriched++;
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   //define the tree to read
   TTree* theTree_DataEG = (TTree*)input->Get("Ttree_DataEG");
   theTree_DataEG->SetBranchAddress( "tree_topMass",	  &topMass);
   theTree_DataEG->SetBranchAddress( "tree_totMass",      &totMass );
   theTree_DataEG->SetBranchAddress( "tree_deltaPhilb",   &deltaPhilb );
   theTree_DataEG->SetBranchAddress( "tree_deltaRlb",	  &deltaRlb);
   theTree_DataEG->SetBranchAddress( "tree_deltaRTopZ",   &deltaRTopZ);
   theTree_DataEG->SetBranchAddress( "tree_asym",         &asym );
   theTree_DataEG->SetBranchAddress( "tree_Zpt",          &Zpt );
   theTree_DataEG->SetBranchAddress( "tree_ZEta",	  &ZEta);
   theTree_DataEG->SetBranchAddress( "tree_topPt",	  &topPt);
   theTree_DataEG->SetBranchAddress( "tree_topEta",	  &topEta);
   theTree_DataEG->SetBranchAddress( "tree_NJets",	  &NJets);
   theTree_DataEG->SetBranchAddress( "tree_NBJets",	  &NBJets);	 
   theTree_DataEG->SetBranchAddress( "tree_deltaRZl",	  &deltaRZl);   
   theTree_DataEG->SetBranchAddress( "tree_deltaPhiZmet", &deltaPhiZmet);
   theTree_DataEG->SetBranchAddress( "tree_btagDiscri",   &btagDiscri);  
   
   theTree_DataEG->SetBranchAddress("tree_leptWPt",	  &leptWPt);		
   theTree_DataEG->SetBranchAddress("tree_leptWEta",	  &leptWEta);		
   theTree_DataEG->SetBranchAddress("tree_leadJetPt",	  &leadJetPt);  	  
   theTree_DataEG->SetBranchAddress("tree_leadJetEta",     &leadJetEta);	 
   //theTree_DataEG->SetBranchAddress("tree_deltaRZleptW",   &deltaRZleptW);     
   theTree_DataEG->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DataEG->GetEntries();ievt++) {    
     if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DataEG->GetEntry(ievt);
     histBdt_Data->Fill( reader->EvaluateMVA( "BDT"           ) );
     
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     fillHisto("Data", theVar, ievt);
     
     
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut) ndataEvents_BCKenriched++;
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   //define the tree to read
   TTree* theTree_DataMuEG = (TTree*)input->Get("Ttree_DataMuEG");
   theTree_DataMuEG->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_DataMuEG->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_DataMuEG->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_DataMuEG->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_DataMuEG->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_DataMuEG->SetBranchAddress( "tree_asym",        &asym );
   theTree_DataMuEG->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_DataMuEG->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_DataMuEG->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_DataMuEG->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_DataMuEG->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_DataMuEG->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_DataMuEG->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_DataMuEG->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_DataMuEG->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   theTree_DataMuEG->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_DataMuEG->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_DataMuEG->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_DataMuEG->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_DataMuEG->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_DataMuEG->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);

   
   
   
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DataMuEG->GetEntries();ievt++) {    
     if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DataMuEG->GetEntry(ievt);
     histBdt_Data->Fill( reader->EvaluateMVA( "BDT"           ) );
     
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     fillHisto("Data", theVar, ievt);
     
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut) ndataEvents_BCKenriched++;
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
  
   
   cout << "total number of events " << theTree_DataMu->GetEntries()+theTree_DataEG->GetEntries()+theTree_DataMuEG->GetEntries() << endl;
   
   
   //-----------------------------------------------------
   //for WZ jets from MC 
   //-----------------------------------------------------
   input->cd();
   
   //define the tree to read
   TTree* theTree_WZ = (TTree*)input->Get("Ttree_WZ");
   theTree_WZ->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_WZ->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_WZ->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_WZ->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_WZ->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_WZ->SetBranchAddress( "tree_asym",        &asym );
   theTree_WZ->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_WZ->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_WZ->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_WZ->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_WZ->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_WZ->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_WZ->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_WZ->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_WZ->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_WZ->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_WZ->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_WZ->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_WZ->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_WZ->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_WZ->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   
   theTree_WZ->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   float nwzEvents_BCKenriched = 0;
   float nwzEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_WZ->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_WZ->GetEntry(ievt);
     
     EvtWeight*=WZscale;
     
     
     histBdt_WZ->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
     
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     fillHisto("WZ", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       nwzEvents_BCKenriched = nwzEvents_BCKenriched+EvtWeight;
       nwzEvents_BCKenriched_unw++;
     }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   
   //-----------------------------------------------------
   //for TTbarSig jets from MC 
   //-----------------------------------------------------
   //define the tree to read
   TTree* theTree_TTbarSig = (TTree*)input->Get("Ttree_TTbar");
   theTree_TTbarSig->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_TTbarSig->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_TTbarSig->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_TTbarSig->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_TTbarSig->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_TTbarSig->SetBranchAddress( "tree_asym",        &asym );
   theTree_TTbarSig->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_TTbarSig->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_TTbarSig->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_TTbarSig->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_TTbarSig->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_TTbarSig->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_TTbarSig->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_TTbarSig->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_TTbarSig->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_TTbarSig->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_TTbarSig->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_TTbarSig->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_TTbarSig->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_TTbarSig->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_TTbarSig->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   
   theTree_TTbarSig->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   float nTTbarSigEvents_BCKenriched = 0;
   float nTTbarSigEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TTbarSig->GetEntries();ievt++) {    
     if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TTbarSig->GetEntry(ievt);
     histBdt_TTbarSig->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
     
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     fillHisto("TTbarSig", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       nTTbarSigEvents_BCKenriched = nTTbarSigEvents_BCKenriched+EvtWeight;
       nTTbarSigEvents_BCKenriched_unw++;
     }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   //-----------------------------------------------------
   //for Wjets jets from MC 
   //-----------------------------------------------------
   
   //define the tree to read
   TTree* theTree_Wjets = (TTree*)input->Get("Ttree_Wjets");
   theTree_Wjets->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_Wjets->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_Wjets->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_Wjets->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_Wjets->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_Wjets->SetBranchAddress( "tree_asym",        &asym );
   theTree_Wjets->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_Wjets->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_Wjets->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_Wjets->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_Wjets->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_Wjets->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_Wjets->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_Wjets->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_Wjets->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_Wjets->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_Wjets->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_Wjets->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_Wjets->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_Wjets->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_Wjets->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   
   theTree_Wjets->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   float nWjetsEvents_BCKenriched = 0;
   float nWjetsEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_Wjets->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_Wjets->GetEntry(ievt);
     histBdt_Wjets->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
     
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     
     fillHisto("Wjets", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       nWjetsEvents_BCKenriched = nWjetsEvents_BCKenriched+EvtWeight;
       nWjetsEvents_BCKenriched_unw++;
     }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   
   
   //-----------------------------------------------------
   //for TtW jets from MC 
   //-----------------------------------------------------
   //define the tree to read
   TTree* theTree_TtW = (TTree*)input->Get("Ttree_TtW");
   theTree_TtW->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_TtW->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_TtW->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_TtW->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_TtW->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_TtW->SetBranchAddress( "tree_asym",        &asym );
   theTree_TtW->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_TtW->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_TtW->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_TtW->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_TtW->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_TtW->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_TtW->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_TtW->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_TtW->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_TtW->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_TtW->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_TtW->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_TtW->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_TtW->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_TtW->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   
   theTree_TtW->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   float nTtWEvents_BCKenriched = 0;
   float nTtWEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TtW->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TtW->GetEntry(ievt);
     histBdt_TtW->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight); 
       
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     
 
     fillHisto("TtW", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       nTtWEvents_BCKenriched = nTtWEvents_BCKenriched+EvtWeight;
       nTtWEvents_BCKenriched_unw++;
      }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   
   
   
   
   //-----------------------------------------------------
   //for TbartW jets from MC 
   //-----------------------------------------------------
   //define the tree to read
   TTree* theTree_TbartW = (TTree*)input->Get("Ttree_TbartW");
   theTree_TbartW->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_TbartW->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_TbartW->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_TbartW->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_TbartW->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_TbartW->SetBranchAddress( "tree_asym",        &asym );
   theTree_TbartW->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_TbartW->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_TbartW->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_TbartW->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_TbartW->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_TbartW->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_TbartW->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_TbartW->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_TbartW->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_TbartW->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_TbartW->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_TbartW->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_TbartW->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_TbartW->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_TbartW->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   
   theTree_TbartW->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   float nTbartWEvents_BCKenriched = 0;
   float nTbartWEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TbartW->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TbartW->GetEntry(ievt);
     histBdt_TbartW->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
       
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     fillHisto("TbartW", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
      nTbartWEvents_BCKenriched = nTbartWEvents_BCKenriched+EvtWeight;
      nTbartWEvents_BCKenriched_unw++;
     }
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   
   
   
   
   
   
   
   
   
   
   //-----------------------------------------------------
   //for TtChan jets from MC 
   //-----------------------------------------------------
   //define the tree to read
   TTree* theTree_TtChan = (TTree*)input->Get("Ttree_TtChan");
   theTree_TtChan->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_TtChan->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_TtChan->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_TtChan->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_TtChan->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_TtChan->SetBranchAddress( "tree_asym",        &asym );
   theTree_TtChan->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_TtChan->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_TtChan->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_TtChan->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_TtChan->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_TtChan->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_TtChan->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_TtChan->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_TtChan->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_TtChan->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_TtChan->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_TtChan->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_TtChan->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_TtChan->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_TtChan->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   
   theTree_TtChan->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   float nTtChanEvents_BCKenriched = 0;
   float nTtChanEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TtChan->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TtChan->GetEntry(ievt);
     histBdt_TtChan->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight); 
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     

     fillHisto("TtChan", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
      nTtChanEvents_BCKenriched = nTtChanEvents_BCKenriched+EvtWeight;
      nTtChanEvents_BCKenriched_unw++;
     }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   
   //-----------------------------------------------------
   //for TbartChan jets from MC 
   //-----------------------------------------------------
   //define the tree to read
   TTree* theTree_TbartChan = (TTree*)input->Get("Ttree_TbartChan");
   theTree_TbartChan->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_TbartChan->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_TbartChan->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_TbartChan->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_TbartChan->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_TbartChan->SetBranchAddress( "tree_asym",        &asym );
   theTree_TbartChan->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_TbartChan->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_TbartChan->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_TbartChan->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_TbartChan->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_TbartChan->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_TbartChan->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_TbartChan->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_TbartChan->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_TbartChan->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_TbartChan->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_TbartChan->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_TbartChan->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_TbartChan->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_TbartChan->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   
   theTree_TbartChan->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   float nTbartChanEvents_BCKenriched = 0;
   float nTbartChanEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TbartChan->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TbartChan->GetEntry(ievt);
     histBdt_TbartChan->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     

     fillHisto("TbartChan", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       nTbartChanEvents_BCKenriched = nTbartChanEvents_BCKenriched+EvtWeight;
       nTbartChanEvents_BCKenriched_unw++;
      }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   
   
   
   
   
   
   //-----------------------------------------------------
   //for TsChan jets from MC 
   //-----------------------------------------------------
   //define the tree to read
   TTree* theTree_TsChan = (TTree*)input->Get("Ttree_TsChan");
   theTree_TsChan->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_TsChan->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_TsChan->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_TsChan->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_TsChan->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_TsChan->SetBranchAddress( "tree_asym",        &asym );
   theTree_TsChan->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_TsChan->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_TsChan->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_TsChan->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_TsChan->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_TsChan->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_TsChan->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_TsChan->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_TsChan->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_TsChan->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_TsChan->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_TsChan->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_TsChan->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_TsChan->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_TsChan->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   
   theTree_TsChan->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   float nTsChanEvents_BCKenriched = 0;
   float nTsChanEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TsChan->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TsChan->GetEntry(ievt);
     histBdt_TsChan->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     

     fillHisto("TsChan", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       nTsChanEvents_BCKenriched = nTsChanEvents_BCKenriched+EvtWeight;
       nTsChanEvents_BCKenriched_unw++;
     }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   
   
   
   //-----------------------------------------------------
   //for TbarsChan jets from MC 
   //-----------------------------------------------------
   //define the tree to read
   TTree* theTree_TbarsChan = (TTree*)input->Get("Ttree_TbarsChan");
   theTree_TbarsChan->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_TbarsChan->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_TbarsChan->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_TbarsChan->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_TbarsChan->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_TbarsChan->SetBranchAddress( "tree_asym",        &asym );
   theTree_TbarsChan->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_TbarsChan->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_TbarsChan->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_TbarsChan->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_TbarsChan->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_TbarsChan->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_TbarsChan->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_TbarsChan->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_TbarsChan->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_TbarsChan->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_TbarsChan->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_TbarsChan->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_TbarsChan->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_TbarsChan->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_TbarsChan->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   
   theTree_TbarsChan->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   float nTbarsChanEvents_BCKenriched = 0;
   float nTbarsChanEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TbarsChan->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TbarsChan->GetEntry(ievt);
     histBdt_TbarsChan->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     fillHisto("TbarsChan", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       nTbarsChanEvents_BCKenriched = nTbarsChanEvents_BCKenriched+EvtWeight;
       nTbarsChanEvents_BCKenriched_unw++;
      }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   
    
   
   
   
   //-----------------------------------------------------
   //for ZZ jets from MC 
   //-----------------------------------------------------
   //define the tree to read
   TTree* theTree_ZZ = (TTree*)input->Get("Ttree_ZZ");
   theTree_ZZ->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_ZZ->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_ZZ->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_ZZ->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_ZZ->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_ZZ->SetBranchAddress( "tree_asym",        &asym );
   theTree_ZZ->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_ZZ->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_ZZ->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_ZZ->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_ZZ->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_ZZ->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_ZZ->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_ZZ->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_ZZ->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_ZZ->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_ZZ->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_ZZ->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_ZZ->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_ZZ->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_ZZ->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   
   theTree_ZZ->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   float nZZEvents_BCKenriched = 0;
   float nZZEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_ZZ->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_ZZ->GetEntry(ievt);
     histBdt_ZZ->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     fillHisto("ZZ", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       nZZEvents_BCKenriched = nZZEvents_BCKenriched+EvtWeight;
       nZZEvents_BCKenriched_unw++;
     }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   
   
   
    
   
   
   
   //-----------------------------------------------------
   //for WW jets from MC 
   //-----------------------------------------------------
   //define the tree to read
   TTree* theTree_WW = (TTree*)input->Get("Ttree_WW");
   theTree_WW->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_WW->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_WW->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_WW->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_WW->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_WW->SetBranchAddress( "tree_asym",        &asym );
   theTree_WW->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_WW->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_WW->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_WW->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_WW->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_WW->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_WW->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_WW->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_WW->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_WW->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_WW->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_WW->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_WW->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_WW->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_WW->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   
   theTree_WW->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   float nWWEvents_BCKenriched = 0;
   float nWWEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_WW->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_WW->GetEntry(ievt);
     histBdt_WW->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     fillHisto("WW", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       nWWEvents_BCKenriched = nWWEvents_BCKenriched+EvtWeight;
       nWWEvents_BCKenriched_unw++;
     }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   
   
   
   
   
   //-----------------------------------------------------
   //for Z jets from MC = > used for normalization
   //-----------------------------------------------------
     //define the tree to read
   TTree* theTree_Zjets = (TTree*)input->Get("Ttree_Zjets");
   theTree_Zjets->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_Zjets->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_Zjets->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_Zjets->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_Zjets->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_Zjets->SetBranchAddress( "tree_asym",        &asym );
   theTree_Zjets->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_Zjets->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_Zjets->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_Zjets->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_Zjets->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_Zjets->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_Zjets->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_Zjets->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_Zjets->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   theTree_Zjets->SetBranchAddress("tree_leptWPt",       &leptWPt);	 
   theTree_Zjets->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_Zjets->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_Zjets->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_Zjets->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_Zjets->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);

   
   theTree_Zjets->SetBranchAddress( "tree_EvtWeight",  &EvtWeight); 
   theTree_Zjets->SetBranchAddress( "tree_Channel",  &Channel);    
   
   
   float mcexpectedZ = 0;
   float mcexpectedZ_unw = 0;
   float mcexpectedZ_mumumu = 0;
   float mcexpectedZ_mumue  = 0;
   float mcexpectedZ_eemu   = 0;
   float mcexpectedZ_eee    = 0;
   
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_Zjets->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_Zjets->GetEntry(ievt);
     histBdt_Zjets ->Fill( reader->EvaluateMVA( "BDT") ,  EvtWeight);
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     fillHisto("Zjets", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       mcexpectedZ = mcexpectedZ+EvtWeight;
       mcexpectedZ_unw++;
       
       if(Channel == 0) mcexpectedZ_mumumu += EvtWeight;
       if(Channel == 1) mcexpectedZ_mumue  += EvtWeight;
       if(Channel == 2) mcexpectedZ_eemu   += EvtWeight;
       if(Channel == 3) mcexpectedZ_eee    += EvtWeight;
       
       
      }
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   
   //-----------------------------------------------------
   //for Z jets from MC = > used for normalization
   //-----------------------------------------------------
     //define the tree to read
   TTree* theTree_DYToLL_M10_50 = (TTree*)input->Get("Ttree_DYToLL_M10-50");
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_asym",        &asym );
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   theTree_DYToLL_M10_50->SetBranchAddress("tree_leptWPt",       &leptWPt);	 
   theTree_DYToLL_M10_50->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_DYToLL_M10_50->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_DYToLL_M10_50->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_DYToLL_M10_50->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_DYToLL_M10_50->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);

   
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_Channel",  &Channel);   
   
   
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DYToLL_M10_50->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DYToLL_M10_50->GetEntry(ievt);
     histBdt_Zjets ->Fill( reader->EvaluateMVA( "BDT") ,  EvtWeight);
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     fillHisto("Zjets", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       mcexpectedZ = mcexpectedZ+EvtWeight;
       mcexpectedZ++;
       if(Channel == 0) mcexpectedZ_mumumu += EvtWeight;
       if(Channel == 1) mcexpectedZ_mumue  += EvtWeight;
       if(Channel == 2) mcexpectedZ_eemu   += EvtWeight;
       if(Channel == 3) mcexpectedZ_eee    += EvtWeight;
     }
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
  
  
  
   
   
   
   
   input_FCNC->cd();
   //define the tree to read
   TTree* theTree_FCNC ;
   if(thevertex == "zut") theTree_FCNC = (TTree*)input->Get("Ttree_FCNCzut");
   if(thevertex == "zct") theTree_FCNC = (TTree*)input->Get("Ttree_FCNCzct");
   if(thevertex == "kut") theTree_FCNC = (TTree*)input->Get("Ttree_FCNCkut");
   if(thevertex == "kct") theTree_FCNC = (TTree*)input->Get("Ttree_FCNCkct");
   if(thevertex == "xut") theTree_FCNC = (TTree*)input->Get("Ttree_FCNCxut");
   if(thevertex == "xct") theTree_FCNC = (TTree*)input->Get("Ttree_FCNCxct");
   
   
   theTree_FCNC->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_FCNC->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_FCNC->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_FCNC->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_FCNC->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_FCNC->SetBranchAddress( "tree_asym",        &asym );
   theTree_FCNC->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_FCNC->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_FCNC->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_FCNC->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_FCNC->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_FCNC->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_FCNC->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_FCNC->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_FCNC->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   theTree_FCNC->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_FCNC->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_FCNC->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_FCNC->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_FCNC->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_FCNC->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);

   
   theTree_FCNC->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   float nsignalEvents_BCKenriched=0;
   float nsignalEvents_BCKenriched_unw=0;
   
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_FCNC->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_FCNC->GetEntry(ievt);
     histBdt_FCNC ->Fill( reader->EvaluateMVA( "BDT") , EvtWeight);
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     fillHisto("Signal", theVar, ievt, EvtWeight*0.1);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       nsignalEvents_BCKenriched = nsignalEvents_BCKenriched + EvtWeight;
       nsignalEvents_BCKenriched_unw++;
     }
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
  
  
  
  
  
   //-----------------------------------------------------
   //for Z jets from data
   //-----------------------------------------------------
   input_zjets->cd();

  
     //define the tree to read
   TTree* theTree_DataMuZjets = (TTree*)input_zjets->Get("Ttree_DataMu");
   theTree_DataMuZjets->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_DataMuZjets->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_DataMuZjets->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_DataMuZjets->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_DataMuZjets->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_DataMuZjets->SetBranchAddress( "tree_asym",        &asym );
   theTree_DataMuZjets->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_DataMuZjets->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_DataMuZjets->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_DataMuZjets->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_DataMuZjets->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_DataMuZjets->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_DataMuZjets->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_DataMuZjets->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_DataMuZjets->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   theTree_DataMuZjets->SetBranchAddress("tree_leptWPt",       &leptWPt);	 
   theTree_DataMuZjets->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_DataMuZjets->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_DataMuZjets->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_DataMuZjets->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_DataMuZjets->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);

   
   theTree_DataMuZjets->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   theTree_DataMuZjets->SetBranchAddress( "tree_Channel",  &Channel);  
   
   float ndataEvents_Ztemplate = 0;
   float ndataEvents_Ztemplate_mumumu = 0;
   
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DataMuZjets->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DataMuZjets->GetEntry(ievt);
     histBdt_DataZjets ->Fill( reader->EvaluateMVA( "BDT") );
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     fillHisto("DataZjets", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
      ndataEvents_Ztemplate++;
      if(Channel == 0) ndataEvents_Ztemplate_mumumu++;
      
      }
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
  
  
  
  
  
  
     //define the tree to read
   TTree* theTree_DataEGZjets = (TTree*)input_zjets->Get("Ttree_DataEG");
   theTree_DataEGZjets->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_DataEGZjets->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_DataEGZjets->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_DataEGZjets->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_DataEGZjets->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_DataEGZjets->SetBranchAddress( "tree_asym",        &asym );
   theTree_DataEGZjets->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_DataEGZjets->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_DataEGZjets->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_DataEGZjets->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_DataEGZjets->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_DataEGZjets->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_DataEGZjets->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_DataEGZjets->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_DataEGZjets->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   theTree_DataEGZjets->SetBranchAddress("tree_leptWPt",       &leptWPt);	 
   theTree_DataEGZjets->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_DataEGZjets->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_DataEGZjets->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_DataEGZjets->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_DataEGZjets->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);

   
   theTree_DataEGZjets->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);
   theTree_DataEGZjets->SetBranchAddress( "tree_Channel",  &Channel);    
   
   
   float ndataEvents_Ztemplate_eee = 0;
   
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DataEGZjets->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DataEGZjets->GetEntry(ievt);
     histBdt_DataZjets ->Fill( reader->EvaluateMVA( "BDT") );
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     fillHisto("DataZjets",theVar , ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
        ndataEvents_Ztemplate++;
     
        if(Channel == 3) ndataEvents_Ztemplate_eee++;
      
      }
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
  
  
  
  
  
  
     //define the tree to read
   TTree* theTree_DataMuEGZjets = (TTree*)input_zjets->Get("Ttree_DataMuEG");
   theTree_DataMuEGZjets->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_DataMuEGZjets->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_DataMuEGZjets->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_DataMuEGZjets->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_DataMuEGZjets->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_DataMuEGZjets->SetBranchAddress( "tree_asym",        &asym );
   theTree_DataMuEGZjets->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_DataMuEGZjets->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_DataMuEGZjets->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_DataMuEGZjets->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_DataMuEGZjets->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_DataMuEGZjets->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_DataMuEGZjets->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_DataMuEGZjets->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_DataMuEGZjets->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   theTree_DataMuEGZjets->SetBranchAddress("tree_leptWPt",       &leptWPt);	 
   theTree_DataMuEGZjets->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_DataMuEGZjets->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_DataMuEGZjets->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_DataMuEGZjets->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_DataMuEGZjets->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);

   
   theTree_DataMuEGZjets->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   theTree_DataMuEGZjets->SetBranchAddress( "tree_Channel",  &Channel);    
   
   
   float ndataEvents_Ztemplate_eemu = 0;
   float ndataEvents_Ztemplate_mumue = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DataMuEGZjets->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DataMuEGZjets->GetEntry(ievt);
     histBdt_DataZjets ->Fill( reader->EvaluateMVA( "BDT") );
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     theVar.push_back(deltaRlb);
     theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     theVar.push_back(leptWPt);  	  
     theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     fillHisto("DataZjets", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
     
       ndataEvents_Ztemplate++;
       
        if(Channel == 1) ndataEvents_Ztemplate_mumue++;
        if(Channel == 2) ndataEvents_Ztemplate_eemu++;
       
     }
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
  
  
  
  
  
  
  
   cout << "************************** " << endl;
   cout << "number of data events with BDT <  " << bdt_BCK_cut << " is " <<  ndataEvents_BCKenriched << endl;
   cout << "number of WZ events with BDT   <  " << bdt_BCK_cut << " is " <<  nwzEvents_BCKenriched << endl;
   cout << "************************** " << endl;
   
   //calculate Z contributionn
   
   double fractionZpassCut = ndataEvents_Ztemplate/
   		(theTree_DataMuZjets->GetEntries()+theTree_DataEGZjets->GetEntries()+theTree_DataMuEGZjets->GetEntries());
   
   double nexpectedZ = fractionZpassCut*histBdt_Zjets->Integral();
   
   cout << "fractionZpassCut " << fractionZpassCut << endl;
   cout << "histBdt_Zjets->Integral() " << histBdt_Zjets->Integral() << endl;
   
   
   
   double nOtherExpected = 
   	  nexpectedZ + 
          nTTbarSigEvents_BCKenriched + 
	  nWjetsEvents_BCKenriched + 
	  nTtWEvents_BCKenriched + 
	  nTbartWEvents_BCKenriched + 
	  nTtChanEvents_BCKenriched + 
	  nTbartChanEvents_BCKenriched + 
	  nTsChanEvents_BCKenriched + 
	  nTbarsChanEvents_BCKenriched + 
	  nZZEvents_BCKenriched + 
	  nWWEvents_BCKenriched;
   
   cout << "number of Z+jets events with BDT   <  " << bdt_BCK_cut << " is " << nexpectedZ << endl;
   cout << "nTTbarSigEvents_BCKenriched  " << nTTbarSigEvents_BCKenriched  << endl;
   cout << "nWjetsEvents_BCKenriched     " << nWjetsEvents_BCKenriched     << endl;
   cout << "nTtWEvents_BCKenriched       " << nTtWEvents_BCKenriched       << endl;
   cout << "nTbartWEvents_BCKenriched    " << nTbartWEvents_BCKenriched    << endl;
   cout << "nTtChanEvents_BCKenriched    " << nTtChanEvents_BCKenriched    << endl;
   cout << "nTbartChanEvents_BCKenriched " << nTbartChanEvents_BCKenriched << endl;
   cout << "nTsChanEvents_BCKenriched    " << nTsChanEvents_BCKenriched    << endl;
   cout << "nTbarsChanEvents_BCKenriched " << nTbarsChanEvents_BCKenriched << endl;
   cout << "nZZEvents_BCKenriched        " << nZZEvents_BCKenriched        << endl;
   cout << "nWWEvents_BCKenriched        " << nWWEvents_BCKenriched        << endl;
   
   
   ndataEvents_BCKenriched -= nOtherExpected;
   
   double WZscaleFactor       = ndataEvents_BCKenriched/nwzEvents_BCKenriched;
   /*double WZscaleFactor_error = pow(  
   				   pow(  pow(ndataEvents_BCKenriched,0.5)/nwzEvents_BCKenriched  , 2)+
				   pow( ndataEvents_BCKenriched*pow(nwzEvents_BCKenriched,0.5)/(nwzEvents_BCKenriched*nwzEvents_BCKenriched)  , 2)
				   ,0.5);
   */
   
   double error_data = pow(  pow(ndataEvents_BCKenriched,0.5)/nwzEvents_BCKenriched, 2);
   
   double error_MC_stat = 0;
   
   if(nTTbarSigEvents_BCKenriched_unw  > 0 ) error_MC_stat +=  pow( (nTTbarSigEvents_BCKenriched/nTTbarSigEvents_BCKenriched_unw)  *pow(nTTbarSigEvents_BCKenriched_unw ,0.5)/nwzEvents_BCKenriched, 2);
   if(nWjetsEvents_BCKenriched_unw     > 0 ) error_MC_stat +=  pow( (nWjetsEvents_BCKenriched/nWjetsEvents_BCKenriched_unw)        *pow(nWjetsEvents_BCKenriched_unw    ,0.5)/nwzEvents_BCKenriched, 2);
   if(nTtWEvents_BCKenriched_unw       > 0 ) error_MC_stat +=  pow( (nTtWEvents_BCKenriched/nTtWEvents_BCKenriched_unw)            *pow(nTtWEvents_BCKenriched_unw      ,0.5)/nwzEvents_BCKenriched, 2);
   if(nTbartWEvents_BCKenriched_unw    > 0 ) error_MC_stat +=  pow( (nTbartWEvents_BCKenriched/nTbartWEvents_BCKenriched_unw)      *pow(nTbartWEvents_BCKenriched_unw   ,0.5)/nwzEvents_BCKenriched, 2);
   if(nTtChanEvents_BCKenriched_unw    > 0 ) error_MC_stat +=  pow( (nTtChanEvents_BCKenriched/nTtChanEvents_BCKenriched_unw)      *pow(nTtChanEvents_BCKenriched_unw   ,0.5)/nwzEvents_BCKenriched, 2);
   if(nTbartChanEvents_BCKenriched_unw > 0 ) error_MC_stat +=  pow( (nTbartChanEvents_BCKenriched/nTbartChanEvents_BCKenriched_unw)*pow(nTbartChanEvents_BCKenriched_unw,0.5)/nwzEvents_BCKenriched, 2);
   if(nTsChanEvents_BCKenriched_unw    > 0 ) error_MC_stat +=  pow( (nTsChanEvents_BCKenriched/nTsChanEvents_BCKenriched_unw)      *pow(nTsChanEvents_BCKenriched_unw   ,0.5)/nwzEvents_BCKenriched, 2);
   if(nTbarsChanEvents_BCKenriched_unw > 0 ) error_MC_stat +=  pow( (nTbarsChanEvents_BCKenriched/nTbarsChanEvents_BCKenriched_unw)*pow(nTbarsChanEvents_BCKenriched_unw,0.5)/nwzEvents_BCKenriched, 2);
   if(nZZEvents_BCKenriched_unw        > 0 ) error_MC_stat +=  pow( (nZZEvents_BCKenriched/nZZEvents_BCKenriched_unw)              *pow(nZZEvents_BCKenriched_unw       ,0.5)/nwzEvents_BCKenriched, 2);
   if(nWWEvents_BCKenriched_unw        > 0 ) error_MC_stat +=  pow( (nWWEvents_BCKenriched/nWWEvents_BCKenriched_unw)              *pow(nWWEvents_BCKenriched_unw       ,0.5)/nwzEvents_BCKenriched, 2);
  
   double normuncert = 0.30;
   double error_MC_weight =  
    pow( nTTbarSigEvents_BCKenriched_unw*normuncert/nwzEvents_BCKenriched, 2);
    
    error_MC_weight += pow( nWjetsEvents_BCKenriched    *normuncert/nwzEvents_BCKenriched, 2);
    error_MC_weight += pow( nTtWEvents_BCKenriched      *normuncert/nwzEvents_BCKenriched, 2);
    error_MC_weight += pow( nTbartWEvents_BCKenriched   *normuncert/nwzEvents_BCKenriched, 2);
    error_MC_weight += pow( nTtChanEvents_BCKenriched   *normuncert/nwzEvents_BCKenriched, 2);
    error_MC_weight += pow( nTsChanEvents_BCKenriched   *normuncert/nwzEvents_BCKenriched, 2);
    error_MC_weight += pow( nTbarsChanEvents_BCKenriched*normuncert/nwzEvents_BCKenriched, 2);
    error_MC_weight += pow( nZZEvents_BCKenriched       *normuncert/nwzEvents_BCKenriched, 2);
    error_MC_weight += pow( nWWEvents_BCKenriched       *normuncert/nwzEvents_BCKenriched, 2);
    
   
   cout << "---------------------------------------" << endl;
   double meanZjetsSF = 
	(3.7687*mcexpectedZ_mumumu +
	1.51623*mcexpectedZ_mumue + 
	3.64732*mcexpectedZ_eemu  +
	2.23334*mcexpectedZ_eee )/
	(mcexpectedZ_mumumu+mcexpectedZ_mumue+mcexpectedZ_eemu+mcexpectedZ_eee)  ;

   double meanZjetsSFError = 
	(1.94*mcexpectedZ_mumumu +
	0.25*mcexpectedZ_mumue + 
	2.63*mcexpectedZ_eemu  +
	0.31*mcexpectedZ_eee)/
	(mcexpectedZ_mumumu+mcexpectedZ_mumue+mcexpectedZ_eemu+mcexpectedZ_eee)
	   ;
   
    
  
   
   
   cout << "ndataEvents_Ztemplate_mumumu " << ndataEvents_Ztemplate_mumumu << endl;
   cout << "ndataEvents_Ztemplate_mumue  " << ndataEvents_Ztemplate_mumue  << endl;
   cout << "ndataEvents_Ztemplate_eemu   " << ndataEvents_Ztemplate_eemu   << endl;
   cout << "ndataEvents_Ztemplate_eee    " << ndataEvents_Ztemplate_eee    << endl;
   cout << "ndataEvents_Ztemplate        " << ndataEvents_Ztemplate        << endl;
   
   cout << "meanZjetsSF                  " << meanZjetsSF  << endl;
   cout << "meanZjetsSFError             " << meanZjetsSFError << endl;
   
   cout << "  " << endl;
   
   
   //uncertainty from DY SF
   double error_Zjets  = pow( meanZjetsSFError*mcexpectedZ/nwzEvents_BCKenriched , 2);
   
   //uncertainty from MC stat ==> not needed
   //cout << "error_Zjets 1 " << error_Zjets << endl;
         // error_Zjets += pow( (mcexpectedZ/mcexpectedZ_unw)*pow(mcexpectedZ_unw, 0.5)/nwzEvents_BCKenriched  , 2);
   
   //cout << "error_Zjets 2 " << error_Zjets << endl;
   
   //uncertainty from fraction of DY events passing the BDT selection < -0.2
   double ntotdata = theTree_DataMuZjets->GetEntries()+theTree_DataEGZjets->GetEntries()+theTree_DataMuEGZjets->GetEntries(); 
   
   
   double fracError = (ntotdata-ndataEvents_Ztemplate)*pow(ndataEvents_Ztemplate, 0.5)/(ntotdata*ntotdata);
   		
 
   
   cout << "   fracError " << fracError << endl;
          error_Zjets += pow( mcexpectedZ*fracError/nwzEvents_BCKenriched, 2);
   cout << "error_Zjets 3 " << error_Zjets << endl;
   
   cout << "---------------------------------------" << endl;
  
   double WZscaleFactor_error = pow(error_data + error_MC_stat + error_MC_weight + error_Zjets, 0.5);
   cout << "uncertainty from Data          " << pow(error_data,      0.5) << endl;
   cout << "uncertainty from Zjets         " << pow(error_Zjets,     0.5) << endl;
   cout << "uncertainty from MC stat       " << pow(error_MC_stat,   0.5) << endl;
   cout << "uncertainty from normalisation " << pow(error_MC_weight, 0.5) << endl;
   
   
   cout << "---------------------------------------" << endl;
   
   cout << "the overall WZ scale factor is " << WZscaleFactor << " pm " << WZscaleFactor_error << endl;
  
   cout << "signal contamination for       <  " << bdt_BCK_cut << " is " <<  nsignalEvents_BCKenriched 
   
     << "  (" <<  nsignalEvents_BCKenriched/ndataEvents_BCKenriched<< "%) " << endl;
   
   
   
   std::cout << "WZ integral 1 " << histBdt_WZ->Integral() << endl;
   histBdt_WZ->Scale(WZscaleFactor);
   std::cout << "WZ integral 2 " << histBdt_WZ->Integral() << endl;
   
   
   
   
   target->cd();
   
   for(unsigned int i=0; i<samples.size(); i++){
   
     scaleHisto(samples[i],  WZscaleFactor);
     //writeHisto(samples[i] );
     
   }
   
  
   target->Write();
   
   
   histBdt_Data->Write();
   histBdt_WZ->Write();
   
   
   histBdt_Zjets->Write();
   histBdt_FCNC->Write();      
   histBdt_DataZjets->Write();
   
   histBdt_TTbarSig->Write();	 
   histBdt_Wjets->Write();	 
   histBdt_TtW->Write();	 
   histBdt_TbartW->Write();   
   histBdt_TtChan ->Write();  
   histBdt_TbartChan->Write();
   histBdt_TsChan->Write();   
   histBdt_TbarsChan->Write();
   histBdt_ZZ->Write();	 
   histBdt_WW->Write();
   
   
   
   
   
   target->Close();
   
}


void ReaderBDT(){

   TString thevertex_zut = "zut";
   TString thevertex_zct = "zct";
   TString thevertex_kut = "kut";
   TString thevertex_kct = "kct";
   TString thevertex_xut = "xut";
   TString thevertex_xct = "xct";
   
   
   
   
   //******************************************
   // for Zut ac
   //******************************************
   
   ReaderBDT(thevertex_zut, 
   	"../RootFiles/backup_outputProof02-01-13_18-48-30_AllSFinclWZ/proof_merged.root",
   	"../RootFiles/backup_outputProof02-01-13_18-48-30_AllSFinclWZ/proof_merged.root",
	 "nom");
   /*
   //for Jes Up
   ReaderBDT(thevertex_zut, 
   	"../RootFiles/backup_outputProof03-01-13_11-04-06_JESup/proof_merged_JESup.root",
   	"../RootFiles/backup_outputProof03-01-13_11-04-06_JESup/proof_merged_JESup.root",
	 "JESup");
	
   //for Jes down
   ReaderBDT(thevertex_zut, 
   	"../RootFiles/backup_outputProof03-01-13_12-02-35_JESdown/proof_merged_JESdown.root",
   	"../RootFiles/backup_outputProof03-01-13_12-02-35_JESdown/proof_merged_JESdown.root",
	 "JESdown");
	 
   //for JER 
   ReaderBDT(thevertex_zut, 
   	"../RootFiles/backup_outputProof03-01-13_12-52-49_JER/proof_merged_JER.root",
   	"../RootFiles/backup_outputProof03-01-13_12-52-49_JER/proof_merged_JER.root",
	 "JER");
	 
   //for LeptSF 
   ReaderBDT(thevertex_zut, 
   	"../RootFiles/backup_outputProof03-01-13_13-56-13_leptonSFoff/proof_merged_leptSFOff.root",
   	"../RootFiles/backup_outputProof03-01-13_13-56-13_leptonSFoff/proof_merged_leptSFOff.root",
	 "LeptSF");
	 
   //for TrigSF 
   ReaderBDT(thevertex_zut, 
   	"../RootFiles/backup_outputProof03-01-13_14-41-10_triggerOff/proof_merged_trigSFOff.root",
   	"../RootFiles/backup_outputProof03-01-13_14-41-10_triggerOff/proof_merged_trigSFOff.root",
	 "TrigSF");
   
   
   //for PU up
   ReaderBDT(thevertex_zut, 
   	"../RootFiles/backup_outputProof03-01-13_15-42-38_PUup/proof_merged_PUup.root",
   	"../RootFiles/backup_outputProof03-01-13_15-42-38_PUup/proof_merged_PUup.root",
	 "PUup");
   
	 
   //for PU down
   ReaderBDT(thevertex_zut, 
   	"../RootFiles/backup_outputProof03-01-13_16-40-33_PUdown/proof_merged_PUdown.root",
   	"../RootFiles/backup_outputProof03-01-13_16-40-33_PUdown/proof_merged_PUdown.root",
	 "PUdown");
   
   //for btag up
   ReaderBDT(thevertex_zut, 
   	"../RootFiles/backup_outputProof03-01-13_18-03-26_btagup/proof_merged_btagup.root",
   	"../RootFiles/backup_outputProof03-01-13_18-03-26_btagup/proof_merged_btagup.root",
	 "btagup");
   
   
   
   //for btag down
   ReaderBDT(thevertex_zut, 
   	"../RootFiles/backup_outputProof03-01-13_19-01-07_btagdown/proof_merged_btagdown.root",
   	"../RootFiles/backup_outputProof03-01-13_19-01-07_btagdown/proof_merged_btagdown.root",
	 "btagdown");
   
      */
   //******************************************
   // for Zct ac
   //******************************************
   
   ReaderBDT(thevertex_zct, 
   	"../RootFiles/backup_outputProof02-01-13_18-48-30_AllSFinclWZ/proof_merged.root",
   	"../RootFiles/backup_outputProof02-01-13_18-48-30_AllSFinclWZ/proof_merged.root",
	 "nom");
   /*
   //for Jes Up
   ReaderBDT(thevertex_zct, 
   	"../RootFiles/backup_outputProof03-01-13_11-04-06_JESup/proof_merged_JESup.root",
   	"../RootFiles/backup_outputProof03-01-13_11-04-06_JESup/proof_merged_JESup.root",
	 "JESup");
	
   //for Jes down
   ReaderBDT(thevertex_zct, 
   	"../RootFiles/backup_outputProof03-01-13_12-02-35_JESdown/proof_merged_JESdown.root",
   	"../RootFiles/backup_outputProof03-01-13_12-02-35_JESdown/proof_merged_JESdown.root",
	 "JESdown");
	 
   //for JER 
   ReaderBDT(thevertex_zct, 
   	"../RootFiles/backup_outputProof03-01-13_12-52-49_JER/proof_merged_JER.root",
   	"../RootFiles/backup_outputProof03-01-13_12-52-49_JER/proof_merged_JER.root",
	 "JER");
	 
   //for LeptSF 
   ReaderBDT(thevertex_zct, 
   	"../RootFiles/backup_outputProof03-01-13_13-56-13_leptonSFoff/proof_merged_leptSFOff.root",
   	"../RootFiles/backup_outputProof03-01-13_13-56-13_leptonSFoff/proof_merged_leptSFOff.root",
	 "LeptSF");
	 
   //for TrigSF 
   ReaderBDT(thevertex_zct, 
   	"../RootFiles/backup_outputProof03-01-13_14-41-10_triggerOff/proof_merged_trigSFOff.root",
   	"../RootFiles/backup_outputProof03-01-13_14-41-10_triggerOff/proof_merged_trigSFOff.root",
	 "TrigSF");
   
   
   //for PU up
   ReaderBDT(thevertex_zct, 
   	"../RootFiles/backup_outputProof03-01-13_15-42-38_PUup/proof_merged_PUup.root",
   	"../RootFiles/backup_outputProof03-01-13_15-42-38_PUup/proof_merged_PUup.root",
	 "PUup");
   
	 
   //for PU down
   ReaderBDT(thevertex_zct, 
   	"../RootFiles/backup_outputProof03-01-13_16-40-33_PUdown/proof_merged_PUdown.root",
   	"../RootFiles/backup_outputProof03-01-13_16-40-33_PUdown/proof_merged_PUdown.root",
	 "PUdown");
   
   //for btag up
   ReaderBDT(thevertex_zct, 
   	"../RootFiles/backup_outputProof03-01-13_18-03-26_btagup/proof_merged_btagup.root",
   	"../RootFiles/backup_outputProof03-01-13_18-03-26_btagup/proof_merged_btagup.root",
	 "btagup");
   
   
   
   //for btag down
   ReaderBDT(thevertex_zct, 
   	"../RootFiles/backup_outputProof03-01-13_19-01-07_btagdown/proof_merged_btagdown.root",
   	"../RootFiles/backup_outputProof03-01-13_19-01-07_btagdown/proof_merged_btagdown.root",
	 "btagdown");
   */
   
}
