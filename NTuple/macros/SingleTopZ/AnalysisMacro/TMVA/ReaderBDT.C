


#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


bool Private = false;

double WZscale = 1.0;
//static   double WZscale_zut = 1.;
/*
static   double WZscale_zut = 0.72;
static   double WZscale_zct = 0.64;
static   double WZscale_kut = 0.82;
static   double WZscale_kct = 0.90;
*/

static   double WZscale_zut = 0.927/0.92;
static   double WZscale_zct = 0.851/0.92;
static   double WZscale_kut = 1.059/0.92 ;
static   double WZscale_kct = 1.051/0.92 ;
 
 
 
static   float bdt_BCK_cut    = 0.1;
static   float bdt_BCK_cut_up = 0.1;
 
static bool   apply_BCK_cut    = true;
static bool   apply_BCK_cut_up = false;

static bool   applyMWTcut = true;
static double mWTcut = 20;
std::map<TString,std::vector<TH1F*> > theHistoMap;

static int chan = -1;

void fillHisto(TString sample, std::vector<double> theVar, int ievt, double weight=1){
   
   
   
   std::vector<TH1F*> histovect = theHistoMap[sample];
   
   histovect[0]->Fill(theVar[0],         weight); //topMass
   //histovect[1]->Fill(theVar[1],         weight); //totMass
   histovect[1]->Fill(fabs(theVar[1]),   weight); //deltaPhilb
   //histovect[2]->Fill(theVar[2],         weight); //deltaRlb
   //histovect[4]->Fill(theVar[4],         weight); //deltaRTopZ
   histovect[2]->Fill(theVar[2],         weight); //asym
   histovect[3]->Fill(theVar[3],         weight); //Zpt
   histovect[4]->Fill(theVar[4],         weight); //ZEta
   //histovect[8]->Fill(theVar[8],         weight); //topPt
   histovect[5]->Fill(theVar[5],         weight); //topEta
   histovect[6]->Fill(theVar[6],       weight); //NJets
   histovect[7]->Fill(theVar[7],       weight); //NBJets
   histovect[8]->Fill(theVar[8],       weight); //deltaRZl
   histovect[9]->Fill(fabs(theVar[9]), weight); //deltaPhiZmet
   histovect[10]->Fill(theVar[10],       weight); //btagDiscri
   //histovect[15]->Fill(theVar[15],       weight); //leptWPt
   //histovect[16]->Fill(theVar[16],       weight); //leptWEta
   histovect[11]->Fill(theVar[11],       weight); //leadJetPt
   histovect[12]->Fill(theVar[12],       weight); //leadJetEta
   histovect[13]->Fill(fabs(theVar[13]), weight); //deltaPhiZleptW
   

 
}


void writeHisto(TString sample){

   std::vector<TH1F*> histovect = theHistoMap[sample];
   
   histovect[0]->Write(); //topMass
   //histovect[1]->Write(); //totMass
   histovect[1]->Write(); //deltaPhilb
   //histovect[2]->Write(); //deltaRlb
   //histovect[4]->Write(); //deltaRTopZ
   histovect[2]->Write(); //asym
   histovect[3]->Write(); //Zpt
   histovect[4]->Write(); //ZEta
   //histovect[8]->Write(); //topPt
   histovect[5]->Write(); //topEta
   histovect[6]->Write(); //NJets
   histovect[7]->Write(); //NBJets
   histovect[8]->Write(); //deltaRZl
   histovect[9]->Write(); //deltaPhiZmet
   histovect[10]->Write(); //btagDiscri
   //histovect[15]->Write(); //leptWPt
   //histovect[16]->Write(); //leptWEta
   histovect[11]->Write(); //leadJetPt
   histovect[12]->Write(); //leadJetEta
   histovect[13]->Write(); //deltaPhiZleptW

}


void scaleHisto(TString sample, double thescale){

   std::vector<TH1F*> histovect = theHistoMap[sample];
   
   histovect[0]->Scale(thescale); //topMass
   //histovect[1]->Scale(thescale); //totMass
   histovect[1]->Scale(thescale); //deltaPhilb
   //histovect[2]->Scale(thescale); //deltaRlb
   //histovect[4]->Scale(thescale); //deltaRTopZ
   histovect[2]->Scale(thescale); //asym
   histovect[3]->Scale(thescale); //Zpt
   histovect[4]->Scale(thescale); //ZEta
   //histovect[8]->Scale(thescale); //topPt
   histovect[5]->Scale(thescale); //topEta
   histovect[6]->Scale(thescale); //NJets
   histovect[7]->Scale(thescale); //NBJets
   histovect[8]->Scale(thescale); //deltaRZl
   histovect[9]->Scale(thescale); //deltaPhiZmet
   histovect[10]->Scale(thescale); //btagDiscri
   //histovect[15]->Scale(thescale); //leptWPt
   //histovect[16]->Scale(thescale); //leptWEta
   histovect[11]->Scale(thescale); //leadJetPt
   histovect[12]->Scale(thescale); //leadJetEta
   histovect[13]->Scale(thescale); //deltaPhiZleptW

}


void ReaderBDT(TString thevertex, TString stringinput, TString stringinput_FCNC , TString syst){

   
   
   
   
   TH2F * pTZ_vs_pTJets_MC   = new TH2F("pTZ_vs_pTJets_MC",   "pTZ_vs_pTJets_MC",   10, 0, 200, 10, 0, 200);
   TH2F * pTZ_vs_pTJets_DY   = new TH2F("pTZ_vs_pTJets_DY",   "pTZ_vs_pTJets_DY",   10, 0, 200, 10, 0, 200);
   TH2F * pTZ_vs_pTJets_Data = new TH2F("pTZ_vs_pTJets_Data", "pTZ_vs_pTJets_Data", 10, 0, 200, 10, 0, 200);
   
   
  //double WZscale = 0.91;
   
   
   // This loads the library
   TMVA::Tools::Instance();
   
   //create the BDT reader
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    
   
    
   TFile *input         = new TFile(stringinput, "read");
   TFile *input_FCNC    = new TFile(stringinput_FCNC, "read");
   
   
   
   TFile *input_data    = new TFile("../../SystRootFiles_2013_07_01/proof_nom.root", "read");
   TFile *input_zjets   = new TFile("../../SystRootFiles_2013_07_01/proof_DYenriched.root", "read");
   
   if(syst=="DYenriched03") input_zjets   = new TFile("../../SystRootFiles_2013_07_01/proof_DYenriched03.root", "read");
   if(syst=="DYenriched05") input_zjets   = new TFile("../../SystRootFiles_2013_07_01/proof_DYenriched05.root", "read");
   
   //*************************************************
   //*************************************************
   //declaration of histograms
   //*************************************************
   //*************************************************
   
   
   TFile *target  ;
   
   if(thevertex == "zut") WZscale = WZscale_zut;
   if(thevertex == "zct") WZscale = WZscale_zct; 
   if(thevertex == "kut") WZscale = WZscale_kut; 
   if(thevertex == "kct") WZscale = WZscale_kct; 
   
   cout << "WZscale " << WZscale << endl;
   
   if(!apply_BCK_cut && !apply_BCK_cut_up && chan == -1){
   
    if(thevertex == "zut") target= new TFile( ("HistoBDToutput/TMVApp_zut_"+syst+".root").Data(),"RECREATE" );
    if(thevertex == "zct") target= new TFile( ("HistoBDToutput/TMVApp_zct_"+syst+".root").Data(),"RECREATE" );
    if(thevertex == "kut") target= new TFile( ("HistoBDToutput/TMVApp_kut_"+syst+".root").Data(),"RECREATE" );
    if(thevertex == "kct") target= new TFile( ("HistoBDToutput/TMVApp_kct_"+syst+".root").Data(),"RECREATE" );
    if(thevertex == "xut") target= new TFile( ("HistoBDToutput/TMVApp_xut_"+syst+".root").Data(),"RECREATE" );
    if(thevertex == "xct") target= new TFile( ("HistoBDToutput/TMVApp_xct_"+syst+".root").Data(),"RECREATE" );
    if(thevertex == "tzq") target= new TFile( ("HistoBDToutput/TMVApp_tzq_"+syst+".root").Data(),"RECREATE" );
   }
   else if(apply_BCK_cut) {
   
     if(thevertex == "zut") target= new TFile( ("HistoBDToutput/TMVApp_zut_bdtcut"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "zct") target= new TFile( ("HistoBDToutput/TMVApp_zct_bdtcut"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "kut") target= new TFile( ("HistoBDToutput/TMVApp_kut_bdtcut"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "kct") target= new TFile( ("HistoBDToutput/TMVApp_kct_bdtcut"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "xut") target= new TFile( ("HistoBDToutput/TMVApp_xut_bdtcut"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "xct") target= new TFile( ("HistoBDToutput/TMVApp_xct_bdtcut"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "tzq") target= new TFile( ("HistoBDToutput/TMVApp_tzq_bdtcut"+syst+".root").Data(),"RECREATE" );
   
   }
   else if(apply_BCK_cut_up) {
   
     if(thevertex == "zut") target= new TFile( ("HistoBDToutput/TMVApp_zut_bdtcutup_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "zct") target= new TFile( ("HistoBDToutput/TMVApp_zct_bdtcutup_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "kut") target= new TFile( ("HistoBDToutput/TMVApp_kut_bdtcutup_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "kct") target= new TFile( ("HistoBDToutput/TMVApp_kct_bdtcutup_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "xut") target= new TFile( ("HistoBDToutput/TMVApp_xut_bdtcutup_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "xct") target= new TFile( ("HistoBDToutput/TMVApp_xct_bdtcutup_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "tzq") target= new TFile( ("HistoBDToutput/TMVApp_tzq_bdtcutup_"+syst+".root").Data(),"RECREATE" );
   
   }
   else if(chan == 0){
     if(thevertex == "zut") target= new TFile( ("HistoBDToutput/TMVApp_zut_mumumu_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "zct") target= new TFile( ("HistoBDToutput/TMVApp_zct_mumumu_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "kut") target= new TFile( ("HistoBDToutput/TMVApp_kut_mumumu_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "kct") target= new TFile( ("HistoBDToutput/TMVApp_kct_mumumu_"+syst+".root").Data(),"RECREATE" );   
   }
   else if(chan == 1){
     if(thevertex == "zut") target= new TFile( ("HistoBDToutput/TMVApp_zut_mumue_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "zct") target= new TFile( ("HistoBDToutput/TMVApp_zct_mumue_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "kut") target= new TFile( ("HistoBDToutput/TMVApp_kut_mumue_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "kct") target= new TFile( ("HistoBDToutput/TMVApp_kct_mumue_"+syst+".root").Data(),"RECREATE" );   
   }
   else if(chan == 2){
     if(thevertex == "zut") target= new TFile( ("HistoBDToutput/TMVApp_zut_eemu_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "zct") target= new TFile( ("HistoBDToutput/TMVApp_zct_eemu_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "kut") target= new TFile( ("HistoBDToutput/TMVApp_kut_eemu_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "kct") target= new TFile( ("HistoBDToutput/TMVApp_kct_eemu_"+syst+".root").Data(),"RECREATE" );   
   }
   else if(chan == 3){
     if(thevertex == "zut") target= new TFile( ("HistoBDToutput/TMVApp_zut_eee_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "zct") target= new TFile( ("HistoBDToutput/TMVApp_zct_eee_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "kut") target= new TFile( ("HistoBDToutput/TMVApp_kut_eee_"+syst+".root").Data(),"RECREATE" );
     if(thevertex == "kct") target= new TFile( ("HistoBDToutput/TMVApp_kct_eee_"+syst+".root").Data(),"RECREATE" );   
   }
   
   std::vector<TString> samples;
   
   samples.push_back("Data");	 
   samples.push_back("Zjets");	 
   samples.push_back("WZ");	 
   samples.push_back("DataZjets");
   samples.push_back("TTbarSig"); 
   samples.push_back("TTbarSigZjets"); 
   samples.push_back("Wjets");	 
   samples.push_back("TtW");	 
   samples.push_back("TbartW");   
   samples.push_back("TtChan");   
   samples.push_back("TbartChan");
   samples.push_back("TsChan");   
   samples.push_back("TbarsChan");
   samples.push_back("TZq");
   samples.push_back("ZZ");	 
   samples.push_back("WW");	 
   samples.push_back("Signal");	 
   
   // Histos for control plots : training without the var. checked
   samples.push_back("BDTcut_Data");	 
   samples.push_back("BDTcut_Zjets");
   samples.push_back("BDTcut_WZ");
   samples.push_back("BDTcut_DataZjets");
   samples.push_back("BDTcut_TTbarSig");
   samples.push_back("BDTcut_TTbarSigZjets");  
   samples.push_back("BDTcut_Wjets");	 
   samples.push_back("BDTcut_TtW");	 
   samples.push_back("BDTcut_TbartW");   
   samples.push_back("BDTcut_TtChan");   
   samples.push_back("BDTcut_TbartChan");
   samples.push_back("BDTcut_TsChan");   
   samples.push_back("BDTcut_TbarsChan");
   samples.push_back("BDTcut_TZq");
   samples.push_back("BDTcut_ZZ");	 
   samples.push_back("BDTcut_WW");	 
   samples.push_back("BDTcut_Signal");	 
   
   
   std::vector<TString> variables;
   

   variables.push_back("topMass");
   //variables.push_back("totMass");
   variables.push_back("deltaPhilb");
   //variables.push_back("deltaRlb");
   //variables.push_back("deltaRTopZ");
   variables.push_back("asym");
   variables.push_back("Zpt");
   variables.push_back("ZEta");
   //variables.push_back("topPt");
   variables.push_back("topEta");
   variables.push_back("NJets");   
   variables.push_back("NBJets");  
   variables.push_back("deltaRZl"); 
   variables.push_back("deltaPhiZmet");
   variables.push_back("btagDiscri");
   //variables.push_back("leptWPt");
   //variables.push_back("leptWEta");	   
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
   Float_t met, mWT;
   
   Int_t Channel;

   reader->AddVariable("tree_topMass",    &topMass    ); 
   reader->AddVariable("tree_deltaPhilb", &deltaPhilb );
   reader->AddVariable("tree_asym",	  &asym       );
   reader->AddVariable("tree_Zpt",	  &Zpt        );
   reader->AddVariable("tree_ZEta",	     &ZEta);
   ////reader->AddVariable("tree_topEta",        &topEta); 
   reader->AddVariable("tree_NJets",         &NJets);	  
   reader->AddVariable("tree_NBJets",        &NBJets);	    
   ////reader->AddVariable("tree_deltaRZl",      &deltaRZl);   
   reader->AddVariable("tree_deltaPhiZmet",  &deltaPhiZmet); 
   //if(thevertex != "zut" && thevertex != "kut" ) reader->AddVariable("tree_btagDiscri",    &btagDiscri);   
   reader->AddVariable("tree_btagDiscri",    &btagDiscri);    
   reader->AddVariable("tree_leadJetEta",     &leadJetEta);
   reader->AddVariable("tree_deltaPhiZleptW", &deltaPhiZleptW); 	 
 
      
    
   //reader->AddVariable("tree_totMass",    &totMass    );
   //reader->AddVariable("tree_topPt",         &topPt);		    
   //reader->AddVariable("tree_leptWPt",        &leptWPt);	    
   //reader->AddVariable("tree_leptWEta",       &leptWEta);	    
   ////reader->AddVariable("tree_leadJetPt",      &leadJetPt);	     
   //reader->AddVariable("tree_deltaRZleptW",   &deltaRZleptW);	
   //reader->AddVariable("tree_deltaRlb",     &deltaRlb);
   //reader->AddVariable("tree_deltaRTopZ",   &deltaRTopZ);
     
   
   UInt_t nbin = 40;
   if(apply_BCK_cut_up) nbin = 2;
   TH1F* histBdt_Data           = new TH1F( "MVA_BDT_Data",      "MVA_BDT_Data",       nbin, -1, 1 );
   TH1F* histBdt_Zjets          = new TH1F( "MVA_BDT_Zjets",     "MVA_BDT_Zjets",      nbin, -1, 1 );
   TH1F* histBdt_WZ             = new TH1F( "MVA_BDT_WZ",        "MVA_BDT_WZ",	       nbin, -1, 1 );
   TH1F* histBdt_DataZjets      = new TH1F( "MVA_BDT_DataZjets", "MVA_BDT_DataZjets",  nbin, -1, 1 );
   
   TH1F* histBdt_TTbarSig       = new TH1F( "MVA_BDT_TTbarSig",	 "MVA_BDT_TTbarSig",   nbin, -1, 1 );
   TH1F* histBdt_TTbarSigZjets  = new TH1F( "MVA_BDT_TTbarSigZjets",	 "MVA_BDT_TTbarSigZjets",   nbin, -1, 1 );
   TH1F* histBdt_Wjets          = new TH1F( "MVA_BDT_Wjets",	 "MVA_BDT_Wjets",      nbin, -1, 1 );
   TH1F* histBdt_TtW            = new TH1F( "MVA_BDT_TtW",	 "MVA_BDT_TtW",        nbin, -1, 1 );
   TH1F* histBdt_TbartW         = new TH1F( "MVA_BDT_TbartW",	 "MVA_BDT_TbartW",     nbin, -1, 1 );
   TH1F* histBdt_TtChan         = new TH1F( "MVA_BDT_TtChan",	 "MVA_BDT_TtChan",     nbin, -1, 1 );
   TH1F* histBdt_TbartChan      = new TH1F( "MVA_BDT_TbartChan", "MVA_BDT_TbartChan",  nbin, -1, 1 );
   TH1F* histBdt_TsChan         = new TH1F( "MVA_BDT_TsChan",	 "MVA_BDT_TsChan",     nbin, -1, 1 );
   TH1F* histBdt_TbarsChan      = new TH1F( "MVA_BDT_TbarsChan", "MVA_BDT_TbarsChan",  nbin, -1, 1 );
   TH1F* histBdt_TZq            = new TH1F( "MVA_BDT_TZq",       "MVA_BDT_TZq",        nbin, -1, 1 );
   TH1F* histBdt_ZZ             = new TH1F( "MVA_BDT_ZZ",	 "MVA_BDT_ZZ",         nbin, -1, 1 );
   TH1F* histBdt_WW             = new TH1F( "MVA_BDT_WW",	 "MVA_BDT_WW",         nbin, -1, 1 );
   
   TH1F* histBdt_FCNC;
   
   if(thevertex == "zut") histBdt_FCNC = new TH1F( "MVA_BDT_FCNC_zut",   "MVA_BDT_FCNC_zut", nbin, -1, 1 );
   if(thevertex == "zct") histBdt_FCNC = new TH1F( "MVA_BDT_FCNC_zct",   "MVA_BDT_FCNC_zct", nbin, -1, 1 );
   if(thevertex == "kut") histBdt_FCNC = new TH1F( "MVA_BDT_FCNC_kut",   "MVA_BDT_FCNC_kut", nbin, -1, 1 );
   if(thevertex == "kct") histBdt_FCNC = new TH1F( "MVA_BDT_FCNC_kct",   "MVA_BDT_FCNC_kct", nbin, -1, 1 );
   if(thevertex == "xut") histBdt_FCNC = new TH1F( "MVA_BDT_FCNC_xut",   "MVA_BDT_FCNC_xut", nbin, -1, 1 );
   if(thevertex == "xct") histBdt_FCNC = new TH1F( "MVA_BDT_FCNC_xct",   "MVA_BDT_FCNC_xct", nbin, -1, 1 );
   //if(thevertex == "xct") histBdt_FCNC = new TH1F( "MVA_BDT_FCNC_xct",   "MVA_BDT_FCNC_xct", nbin, -1, 1 );
   
   //Get the BDT trainning through an xml file, created during the trainning phase
   if(thevertex == "zut") reader->BookMVA( "BDT", "weights/BDT_trainning_zut_BDT.weights.xml" ); 
   if(thevertex == "zct") reader->BookMVA( "BDT", "weights/BDT_trainning_zct_BDT.weights.xml" ); 
   if(thevertex == "kut") reader->BookMVA( "BDT", "weights/BDT_trainning_kut_BDT.weights.xml" ); 
   if(thevertex == "kct") reader->BookMVA( "BDT", "weights/BDT_trainning_kct_BDT.weights.xml" ); 
   if(thevertex == "xut") reader->BookMVA( "BDT", "weights/BDT_trainning_xut_BDT.weights.xml" ); 
   if(thevertex == "xct") reader->BookMVA( "BDT", "weights/BDT_trainning_xct_BDT.weights.xml" ); 
   //if(thevertex == "tzq") reader->BookMVA( "BDT", "weights/BDT_trainning_tZq_BDT.weights.xml" ); 
   
   
   histBdt_Data->Sumw2()        ;
   histBdt_Zjets->Sumw2()       ;
   histBdt_WZ->Sumw2()	       ;
   histBdt_DataZjets->Sumw2()   ;
   
   histBdt_TTbarSig->Sumw2()    ;
   histBdt_TTbarSigZjets->Sumw2()    ;
   histBdt_Wjets->Sumw2()       ;
   histBdt_TtW->Sumw2()         ;
   histBdt_TbartW->Sumw2()      ;
   histBdt_TtChan->Sumw2()      ;
   histBdt_TbartChan->Sumw2()   ;
   histBdt_TsChan->Sumw2()      ;
   histBdt_TbarsChan->Sumw2()   ;
   histBdt_TZq->Sumw2()        ;
   histBdt_ZZ->Sumw2()	       ;
   histBdt_WW->Sumw2()	       ;
 
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   float ndataEvents_BCKenriched = 0;
   
   //-----------------------------------------------------
   //for Data
   //-----------------------------------------------------
   cout << "start to read data " << endl;
   input_data->cd();
   //define the tree to read
   TTree* theTree_DataMu = (TTree*)input_data->Get("Ttree_DataMu");
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
   theTree_DataMu->SetBranchAddress("tree_met", &met);	
   if(applyMWTcut)theTree_DataMu->SetBranchAddress("tree_mTW", &mWT);	
   theTree_DataMu->SetBranchAddress("tree_Channel", &Channel);
    
    
    cout << "line 392 " << endl;	
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DataMu->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     
     
     theTree_DataMu->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut && applyMWTcut) continue;
     
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     histBdt_Data->Fill( reader->EvaluateMVA( "BDT"           ) );
     
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     ////theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     pTZ_vs_pTJets_Data->Fill(Zpt, leadJetPt);
     fillHisto("Data", theVar, ievt);
     
     
     
     if(reader->EvaluateMVA( "BDT") > 0.4){
       cout << "*********** " << reader->EvaluateMVA( "BDT") << endl;
       cout << "found large event, pt Z  " << Zpt << endl;
       cout << "found large event, topMass " << topMass << endl;
       cout << "found large event, deltaPhilb " << deltaPhilb << endl;
       cout << "found large event, asym " <<  asym<< endl;
       cout << "found large event, ZEta " <<  ZEta<< endl;
       cout << "found large event, topEta " << topEta << endl;
       cout << "found large event, NJets " << NJets << endl;
       cout << "found large event, NBJets " << NBJets << endl;
       cout << "found large event, deltaRZl " << deltaRZl << endl;
       cout << "found large event, deltaPhiZmet " << deltaPhiZmet << endl;
       cout << "found large event, btagDiscri " <<  btagDiscri<< endl;
        
       cout << "found large event, leptWPt  " <<  leptWPt<< endl;		
       cout << "found large event, leptWEta  " <<  leptWEta<< endl;		
       cout << "found large event, leadJetPt  " << leadJetPt << endl;  	  
       cout << "found large event, leadJetEta  " << leadJetEta << endl;
       cout << "found large event, met  " << met << endl;
       cout << "found large event, mWT  " << mWT << endl;
        
     }
     
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut)
     {
       fillHisto("BDTcut_Data", theVar, ievt);
       ndataEvents_BCKenriched++;
     }
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   //define the tree to read
   TTree* theTree_DataEG = (TTree*)input_data->Get("Ttree_DataEG");
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
   theTree_DataEG->SetBranchAddress("tree_met", &met);	
   if(applyMWTcut)theTree_DataEG->SetBranchAddress("tree_mTW", &mWT);
   theTree_DataEG->SetBranchAddress("tree_Channel", &Channel);		
   
   cout << "line 645 " <<endl;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DataEG->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DataEG->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     histBdt_Data->Fill( reader->EvaluateMVA( "BDT"           ) );
     
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     ////theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     pTZ_vs_pTJets_Data->Fill(Zpt, leadJetPt);
     fillHisto("Data", theVar, ievt);
     
     
     if(reader->EvaluateMVA( "BDT") > 0.4){
       cout << "*********** " << reader->EvaluateMVA( "BDT") << endl;
       cout << "found large event, pt Z  " << Zpt << endl;
       cout << "found large event, topMass " << topMass << endl;
       cout << "found large event, deltaPhilb " << deltaPhilb << endl;
       cout << "found large event, asym " <<  asym<< endl;
       cout << "found large event, ZEta " <<  ZEta<< endl;
       cout << "found large event, topEta " << topEta << endl;
       cout << "found large event, NJets " << NJets << endl;
       cout << "found large event, NBJets " << NBJets << endl;
       cout << "found large event, deltaRZl " << deltaRZl << endl;
       cout << "found large event, deltaPhiZmet " << deltaPhiZmet << endl;
       cout << "found large event, btagDiscri " <<  btagDiscri<< endl;
       cout << "found large event, leptWPt  " <<  leptWPt<< endl;		
       cout << "found large event, leptWEta  " <<  leptWEta<< endl;		
       cout << "found large event, leadJetPt  " << leadJetPt << endl;  	  
       cout << "found large event, leadJetEta  " << leadJetEta << endl;
       cout << "found large event, met  " << met << endl;
       cout << "found large event, mWT  " << mWT << endl;
        
        
     }
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut)
     {
       fillHisto("BDTcut_Data", theVar, ievt);
       ndataEvents_BCKenriched++;
     }
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   //define the tree to read
   TTree* theTree_DataMuEG = (TTree*)input_data->Get("Ttree_DataMuEG");
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
   theTree_DataMuEG->SetBranchAddress("tree_met", &met);	
   if(applyMWTcut)theTree_DataMuEG->SetBranchAddress("tree_mTW", &mWT);	
   theTree_DataMuEG->SetBranchAddress("tree_Channel", &Channel);			

   
   
   
   cout << "line 539 " <<endl;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DataMuEG->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DataMuEG->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     histBdt_Data->Fill( reader->EvaluateMVA( "BDT"           ) );
     
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     ////theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     
     if(reader->EvaluateMVA( "BDT") > 0.4){
       cout << "*********** " << reader->EvaluateMVA( "BDT") << endl;
       cout << "found large event, pt Z  " << Zpt << endl;
       cout << "found large event, topMass " << topMass << endl;
       cout << "found large event, deltaPhilb " << deltaPhilb << endl;
       cout << "found large event, asym " <<  asym<< endl;
       cout << "found large event, ZEta " <<  ZEta<< endl;
       cout << "found large event, topEta " << topEta << endl;
       cout << "found large event, NJets " << NJets << endl;
       cout << "found large event, NBJets " << NBJets << endl;
       cout << "found large event, deltaRZl " << deltaRZl << endl;
       cout << "found large event, deltaPhiZmet " << deltaPhiZmet << endl;
       cout << "found large event, btagDiscri " <<  btagDiscri<< endl;
       cout << "found large event, leptWPt  " <<  leptWPt<< endl;		
       cout << "found large event, leptWEta  " <<  leptWEta<< endl;		
       cout << "found large event, leadJetPt  " << leadJetPt << endl;  	  
       cout << "found large event, leadJetEta  " << leadJetEta << endl;
       cout << "found large event, met  " << met << endl;
       cout << "found large event, mWT  " << mWT << endl;
        
        
        
     }
     pTZ_vs_pTJets_Data->Fill(Zpt, leadJetPt);
     fillHisto("Data", theVar, ievt);
     
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut)
     {
       fillHisto("BDTcut_Data", theVar, ievt);
       ndataEvents_BCKenriched++;
     }
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
  
   
   cout << "total number of events " << theTree_DataMu->GetEntries()+theTree_DataEG->GetEntries()+theTree_DataMuEG->GetEntries() << endl;
   
   
   //-----------------------------------------------------
   //for WZ jets from MC 
   //-----------------------------------------------------
   cout << "start to read WZ " << endl;
   input->cd();
   
   //define the tree to read
   TTree* theTree_WZ ;
   
   if(syst == "Scaleup")        { theTree_WZ = (TTree*)input->Get("Ttree_WZprivate_scaleup");  }
   else if(syst == "Scaledown") { theTree_WZ = (TTree*)input->Get("Ttree_WZprivate_scaledown");}
   else if(syst == "Matchup")   { theTree_WZ = (TTree*)input->Get("Ttree_WZprivate_matchup");  }
   else if(syst == "Matchdown") { theTree_WZ = (TTree*)input->Get("Ttree_WZprivate_matchdown"); cout << "in matchdown " << endl;}
   else if(Private)             { theTree_WZ = (TTree*)input->Get("Ttree_WZprivate");}
   else                         { theTree_WZ = (TTree*)input->Get("Ttree_WZ");}
   
  
   if(theTree_WZ == 0) cout << "null pointer for theTree_WZ  " << endl;
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
   if(applyMWTcut)theTree_WZ->SetBranchAddress("tree_mTW", &mWT);
   theTree_WZ->SetBranchAddress("tree_Channel", &Channel);		
   
   
    theTree_WZ->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   double sumWZ_nocut = 0;
   
   
   float nwzEvents_BCKenriched = 0;
   float nwzEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_WZ->GetEntries();ievt++) {    
     if (ievt%1000 == 0) std::cout << "--- ... WZ Processing event: " << ievt << std::endl;
     theTree_WZ->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     sumWZ_nocut +=  EvtWeight*WZscale;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     
     //cout << reader->EvaluateMVA( "BDT")  <<   "  " << EvtWeight << endl;
     histBdt_WZ->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight*WZscale);
     
     //if(EvtWeight == nan) cout << "******** found diverging value " << endl;
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     ////theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     pTZ_vs_pTJets_MC->Fill(Zpt, leadJetPt); 
     
     
     fillHisto("WZ", theVar, ievt, EvtWeight*WZscale);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_WZ", theVar, ievt, EvtWeight);
       nwzEvents_BCKenriched = nwzEvents_BCKenriched+EvtWeight;
       nwzEvents_BCKenriched_unw++;
     }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   cout << "712 " << endl;
   //to get the proper normalization
   TTree* theTree_WZ_nom = (TTree*)input_data->Get("Ttree_WZ");
   
   theTree_WZ_nom->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_WZ_nom->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_WZ_nom->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_WZ_nom->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_WZ_nom->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_WZ_nom->SetBranchAddress( "tree_asym",        &asym );
   theTree_WZ_nom->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_WZ_nom->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_WZ_nom->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_WZ_nom->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_WZ_nom->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_WZ_nom->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_WZ_nom->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_WZ_nom->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_WZ_nom->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_WZ_nom->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_WZ_nom->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_WZ_nom->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_WZ_nom->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_WZ_nom->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_WZ_nom->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   	
   if(applyMWTcut)theTree_WZ_nom->SetBranchAddress("tree_mTW", &mWT);	
   
   theTree_WZ_nom->SetBranchAddress( "tree_EvtWeight",  &EvtWeight); 
   theTree_WZ_nom->SetBranchAddress("tree_Channel", &Channel);		 
     
   double sumWeight_WZnom=0;
   
   for (Long64_t ievt=0; ievt<theTree_WZ_nom->GetEntries();ievt++) { 
     theTree_WZ_nom->GetEntry(ievt); 
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     //if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut   ) continue;
     //if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     sumWeight_WZnom +=EvtWeight*WZscale;
   }  
 
   //cout << "sumWeight_WZnom " << sumWeight_WZnom << endl;
   //cout << "histBdt_WZ->Integral() " << histBdt_WZ->Integral() << endl;
   //histBdt_WZ->Scale(sumWeight_WZnom/histBdt_WZ->Integral());
   histBdt_WZ->Scale(sumWeight_WZnom/sumWZ_nocut);
   
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
   
   if(applyMWTcut)theTree_TTbarSig->SetBranchAddress("tree_mTW", &mWT);	
   
   theTree_TTbarSig->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   theTree_TTbarSig->SetBranchAddress("tree_Channel", &Channel);		 
   
   float nTTbarSigEvents_BCKenriched = 0;
   float nTTbarSigEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TTbarSig->GetEntries();ievt++) {    
     if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TTbarSig->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_TTbarSig->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
     
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     pTZ_vs_pTJets_MC->Fill(Zpt, leadJetPt); 
     fillHisto("TTbarSig", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_TTbarSig", theVar, ievt, EvtWeight);
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
   
   
   if(applyMWTcut)theTree_Wjets->SetBranchAddress("tree_mTW", &mWT);
   theTree_Wjets->SetBranchAddress( "tree_EvtWeight",  &EvtWeight); 
   theTree_Wjets->SetBranchAddress("tree_Channel", &Channel);	 
   
   float nWjetsEvents_BCKenriched = 0;
   float nWjetsEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_Wjets->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_Wjets->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_Wjets->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
     
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     
     pTZ_vs_pTJets_MC->Fill(Zpt, leadJetPt); 
     fillHisto("Wjets", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_Wjets", theVar, ievt, EvtWeight);
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
   
   
   if(applyMWTcut)theTree_TtW->SetBranchAddress("tree_mTW", &mWT);
   theTree_TtW->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   theTree_TtW->SetBranchAddress("tree_Channel", &Channel);	 
   
   float nTtWEvents_BCKenriched = 0;
   float nTtWEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TtW->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TtW->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_TtW->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight); 
       
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     
     pTZ_vs_pTJets_MC->Fill(Zpt, leadJetPt); 
 
     fillHisto("TtW", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_TtW", theVar, ievt, EvtWeight);
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
   
   if(applyMWTcut)theTree_TbartW->SetBranchAddress("tree_mTW", &mWT);
   
   theTree_TbartW->SetBranchAddress( "tree_EvtWeight",  &EvtWeight); 
   theTree_TbartW->SetBranchAddress("tree_Channel", &Channel);	  
   
   float nTbartWEvents_BCKenriched = 0;
   float nTbartWEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TbartW->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TbartW->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut)    continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_TbartW->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
       
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     pTZ_vs_pTJets_MC->Fill(Zpt, leadJetPt); 
     fillHisto("TbartW", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
      fillHisto("BDTcut_TbartW", theVar, ievt, EvtWeight);
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
   
   if(applyMWTcut)theTree_TtChan->SetBranchAddress("tree_mTW", &mWT);
   
   theTree_TtChan->SetBranchAddress( "tree_EvtWeight",  &EvtWeight); 
   theTree_TtChan->SetBranchAddress("tree_Channel", &Channel);	   
   
   float nTtChanEvents_BCKenriched = 0;
   float nTtChanEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TtChan->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TtChan->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_TtChan->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight); 
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     pTZ_vs_pTJets_MC->Fill(Zpt, leadJetPt); 

     fillHisto("TtChan", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
      fillHisto("BDTcut_TtChan", theVar, ievt, EvtWeight);
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
   
   if(applyMWTcut)theTree_TbartChan->SetBranchAddress("tree_mTW", &mWT);
   
   theTree_TbartChan->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);
   theTree_TbartChan->SetBranchAddress("tree_Channel", &Channel);	     
   
   float nTbartChanEvents_BCKenriched = 0;
   float nTbartChanEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TbartChan->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TbartChan->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_TbartChan->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     

     pTZ_vs_pTJets_MC->Fill(Zpt, leadJetPt); 
     fillHisto("TbartChan", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_TbartChan", theVar, ievt, EvtWeight);
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
   
   if(applyMWTcut)theTree_TsChan->SetBranchAddress("tree_mTW", &mWT);
   
   theTree_TsChan->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);
   theTree_TsChan->SetBranchAddress("tree_Channel", &Channel);	     
   
   float nTsChanEvents_BCKenriched = 0;
   float nTsChanEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TsChan->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TsChan->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_TsChan->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     pTZ_vs_pTJets_MC->Fill(Zpt, leadJetPt); 

     fillHisto("TsChan", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_TsChan", theVar, ievt, EvtWeight);
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
   
   
   if(applyMWTcut)theTree_TbarsChan->SetBranchAddress("tree_mTW", &mWT);
   theTree_TbarsChan->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);
   theTree_TbarsChan->SetBranchAddress("tree_Channel", &Channel);	       
   
   float nTbarsChanEvents_BCKenriched = 0;
   float nTbarsChanEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TbarsChan->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TbarsChan->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_TbarsChan->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     ////theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     pTZ_vs_pTJets_MC->Fill(Zpt, leadJetPt); 
     fillHisto("TbarsChan", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_TbarsChan", theVar, ievt, EvtWeight);
       nTbarsChanEvents_BCKenriched = nTbarsChanEvents_BCKenriched+EvtWeight;
       nTbarsChanEvents_BCKenriched_unw++;
      }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   
    
   
   cout << "1323  " << endl;
  //-----------------------------------------------------
   //for TZq jets from MC 
   //-----------------------------------------------------
   //define the tree to read
   TTree* theTree_TZq ;
   if     (syst == "Matchup" )   theTree_TZq= (TTree*)input->Get("Ttree_TZq_matchup");
   else if(syst == "Matchdown" ) theTree_TZq= (TTree*)input->Get("Ttree_TZq_matchdown");
   else                        theTree_TZq= (TTree*)input->Get("Ttree_TZq");
   theTree_TZq->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_TZq->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_TZq->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_TZq->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_TZq->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_TZq->SetBranchAddress( "tree_asym",        &asym );
   theTree_TZq->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_TZq->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_TZq->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_TZq->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_TZq->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_TZq->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_TZq->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_TZq->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_TZq->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_TZq->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_TZq->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_TZq->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_TZq->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_TZq->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_TZq->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   
   if(applyMWTcut)theTree_TZq->SetBranchAddress("tree_mTW", &mWT);
   theTree_TZq->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);
   theTree_TZq->SetBranchAddress("tree_Channel", &Channel);	    
   
   float nTZqEvents_BCKenriched = 0;
   float nTZqEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TZq->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TZq->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_TZq->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
    // //theVar.push_back(deltaRlb);
    // //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     pTZ_vs_pTJets_MC->Fill(Zpt, leadJetPt); 
     fillHisto("TZq", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_TZq", theVar, ievt, EvtWeight);
       nTZqEvents_BCKenriched = nTZqEvents_BCKenriched+EvtWeight;
       nTZqEvents_BCKenriched_unw++;
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
   
   
   if(applyMWTcut)theTree_ZZ->SetBranchAddress("tree_mTW", &mWT);
   theTree_ZZ->SetBranchAddress( "tree_EvtWeight",  &EvtWeight); 
   theTree_ZZ->SetBranchAddress("tree_Channel", &Channel);	     
   
   float nZZEvents_BCKenriched = 0;
   float nZZEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_ZZ->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_ZZ->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_ZZ->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ///theVar.push_back(deltaRlb);
     ////theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     pTZ_vs_pTJets_MC->Fill(Zpt, leadJetPt); 
     fillHisto("ZZ", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_ZZ", theVar, ievt, EvtWeight);
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
   
   
   if(applyMWTcut) theTree_WW->SetBranchAddress("tree_mTW", &mWT);
   
   theTree_WW->SetBranchAddress( "tree_EvtWeight",  &EvtWeight); 
   theTree_WW->SetBranchAddress("tree_Channel", &Channel);	     
   
   float nWWEvents_BCKenriched = 0;
   float nWWEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_WW->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_WW->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_WW->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     //theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     //theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     //theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     //theVar.push_back(leptWPt);  	  
     //theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     pTZ_vs_pTJets_MC->Fill(Zpt, leadJetPt); 
     fillHisto("WW", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_WW", theVar, ievt, EvtWeight);
       nWWEvents_BCKenriched = nWWEvents_BCKenriched+EvtWeight;
       nWWEvents_BCKenriched_unw++;
     }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
   
   
   
   
   
   
   cout << "1567  " << endl;
   
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

   
   if(applyMWTcut) theTree_Zjets->SetBranchAddress("tree_mTW", &mWT);
   theTree_Zjets->SetBranchAddress( "tree_EvtWeight",  &EvtWeight); 
   theTree_Zjets->SetBranchAddress( "tree_Channel",  &Channel);  
   
   
   float mcexpectedZ = 0;
   float mcexpectedZ_unw = 0;
   float mcexpectedZ_mumumu = 0;
   float mcexpectedZ_mumue  = 0;
   float mcexpectedZ_eemu   = 0;
   float mcexpectedZ_eee    = 0;
   
   double totexpectedZjets = 0;
   
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_Zjets->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_Zjets->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     totexpectedZjets += EvtWeight;
     if(apply_BCK_cut && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_Zjets ->Fill( reader->EvaluateMVA( "BDT") ,  EvtWeight);
         
     std::vector<double > theVar;
     theVar.push_back(topMass);
     //theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     //theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     //theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     //theVar.push_back(leptWPt);  	  
     //theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     fillHisto("Zjets", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_Zjets", theVar, ievt, EvtWeight);
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
   
   if(applyMWTcut) theTree_DYToLL_M10_50->SetBranchAddress("tree_leptWPt",       &leptWPt);	 
   theTree_DYToLL_M10_50->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_DYToLL_M10_50->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_DYToLL_M10_50->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_DYToLL_M10_50->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_DYToLL_M10_50->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);

   
   theTree_DYToLL_M10_50->SetBranchAddress("tree_mTW", &mWT);
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   theTree_DYToLL_M10_50->SetBranchAddress( "tree_Channel",  &Channel);   
   
   
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DYToLL_M10_50->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DYToLL_M10_50->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     totexpectedZjets += EvtWeight;
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_Zjets ->Fill( reader->EvaluateMVA( "BDT") ,  EvtWeight);
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     //theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     //theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     //theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     //theVar.push_back(leptWPt);  	  
     //theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     fillHisto("Zjets", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_Zjets", theVar, ievt, EvtWeight);
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
  
  
  
   
   cout << "1737  " << endl;
   
   
   input_FCNC->cd();
   //define the tree to read
   TTree* theTree_FCNC ;
   
   if(syst == "Matchup"){
     if(thevertex == "zut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzut_matchup");
     if(thevertex == "zct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzct_matchup");
     if(thevertex == "kut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkut_matchup");
     if(thevertex == "kct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkct_matchup");
     if(thevertex == "xut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxut_matchup");
     if(thevertex == "xct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxct_matchup");
   }else if(syst == "Matchdown"){
     if(thevertex == "zut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzut_matchdown");
     if(thevertex == "zct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzct_matchdown");
     if(thevertex == "kut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkut_matchdown");
     if(thevertex == "kct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkct_matchdown");
     if(thevertex == "xut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxut_matchdown");
     if(thevertex == "xct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxct_matchdown");
   }else if(syst == "Scaleup"){
     if(thevertex == "zut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzut_scaleup");
     if(thevertex == "zct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzct_scaleup");
     if(thevertex == "kut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkut_scaleup");
     if(thevertex == "kct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkct_scaleup");
     if(thevertex == "xut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxut_scaleup");
     if(thevertex == "xct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxct_scaleup");
   }else if(syst == "Scaledown"){
     if(thevertex == "zut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzut_scaledown");
     if(thevertex == "zct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzct_scaledown");
     if(thevertex == "kut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkut_scaledown");
     if(thevertex == "kct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkct_scaledown");
     if(thevertex == "xut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxut_scaledown");
     if(thevertex == "xct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxct_scaledown");
   }else if(syst == "Mtopup"){
     if(thevertex == "zut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzut_topup");
     if(thevertex == "zct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzct_topup");
     if(thevertex == "kut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkut_topup");
     if(thevertex == "kct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkct_topup");
     if(thevertex == "xut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxut_topup");
     if(thevertex == "xct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxct_topup");
   }else if(syst == "Mtopdown"){
     if(thevertex == "zut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzut_topdown");
     if(thevertex == "zct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzct_topdown");
     if(thevertex == "kut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkut_topdown");
     if(thevertex == "kct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkct_topdown");
     if(thevertex == "xut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxut_topdown");
     if(thevertex == "xct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxct_topdown");
   }else {
     if(thevertex == "zut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzutFullSim");
     if(thevertex == "zct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzctFullSim");
     if(thevertex == "kut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkutFullSim");
     if(thevertex == "kct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkctFullSim");
     if(thevertex == "xut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxutFullSim");
     if(thevertex == "xct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxctFullSim");
     
     /*if(thevertex == "zut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzut");
     if(thevertex == "zct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCzct");
     if(thevertex == "kut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkut");
     if(thevertex == "kct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCkct");
     if(thevertex == "xut") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxut");
     if(thevertex == "xct") theTree_FCNC = (TTree*)input_FCNC->Get("Ttree_FCNCxct");*/
   }
   
   float nsignalEvents_BCKenriched=0;
   float nsignalEvents_BCKenriched_unw=0;
   if(thevertex != "tzq"){
   cout << "start reading FCNC " << endl;
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

   
   if(applyMWTcut) theTree_FCNC->SetBranchAddress("tree_mTW", &mWT);
   theTree_FCNC->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   theTree_FCNC->SetBranchAddress( "tree_Channel",  &Channel);  
   
   
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_FCNC->GetEntries();ievt++) {    
     if (ievt%100 == 0) std::cout << "--- ... FCNC Processing event: " << ievt << std::endl;
     theTree_FCNC->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     if(apply_BCK_cut && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     
     
     histBdt_FCNC ->Fill( reader->EvaluateMVA( "BDT") , EvtWeight);
     
     //cout << reader->EvaluateMVA( "BDT") << "  " << EvtWeight << endl;
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     //theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     //theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     //theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     //theVar.push_back(leptWPt);  	  
     //theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     fillHisto("Signal", theVar, ievt, EvtWeight*0.1);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_Signal", theVar, ievt, EvtWeight*0.1);
       nsignalEvents_BCKenriched = nsignalEvents_BCKenriched + EvtWeight;
       nsignalEvents_BCKenriched_unw++;
     }
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
  
  
  
  
  }
   //-----------------------------------------------------
   //for Z jets from data
   //-----------------------------------------------------
   input_zjets->cd();
   
   
   double totexpectedZjets_data = 0;
  
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

   
   if(applyMWTcut) theTree_DataMuZjets->SetBranchAddress("tree_mTW", &mWT);
   theTree_DataMuZjets->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   theTree_DataMuZjets->SetBranchAddress( "tree_Channel",  &Channel);  
   
   float ndataEvents_Ztemplate = 0;
   float ndataEvents_Ztemplate_mumumu = 0;
   
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DataMuZjets->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DataMuZjets->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     totexpectedZjets_data++;
     if(apply_BCK_cut && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     histBdt_DataZjets ->Fill( reader->EvaluateMVA( "BDT") );
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     //theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     //theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     //theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     //theVar.push_back(leptWPt);  	  
     //theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     pTZ_vs_pTJets_DY->Fill(Zpt, leadJetPt);
     
      
     fillHisto("DataZjets", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
      fillHisto("BDTcut_DataZjets", theVar, ievt, EvtWeight);
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

   
   if(applyMWTcut) theTree_DataEGZjets->SetBranchAddress("tree_mTW", &mWT);
   theTree_DataEGZjets->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);
   theTree_DataEGZjets->SetBranchAddress( "tree_Channel",  &Channel);    
   
   
   float ndataEvents_Ztemplate_eee = 0;
   
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DataEGZjets->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DataEGZjets->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     totexpectedZjets_data++;
     if(apply_BCK_cut && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     histBdt_DataZjets ->Fill( reader->EvaluateMVA( "BDT") );
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     //theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     //theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     //theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     //theVar.push_back(leptWPt);  	  
     //theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     pTZ_vs_pTJets_DY->Fill(Zpt, leadJetPt);
     fillHisto("DataZjets",theVar , ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
        fillHisto("BDTcut_DataZjets",theVar , ievt, EvtWeight);
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

   
   if(applyMWTcut) theTree_DataMuEGZjets->SetBranchAddress("tree_mTW", &mWT);
   theTree_DataMuEGZjets->SetBranchAddress( "tree_EvtWeight",  &EvtWeight);  
   
   theTree_DataMuEGZjets->SetBranchAddress( "tree_Channel",  &Channel);    
   
   
   float ndataEvents_Ztemplate_eemu = 0;
   float ndataEvents_Ztemplate_mumue = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_DataMuEGZjets->GetEntries();ievt++) {    
     //if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_DataMuEGZjets->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     totexpectedZjets_data++;
     if(apply_BCK_cut && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     histBdt_DataZjets ->Fill( reader->EvaluateMVA( "BDT") );
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     //theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     //theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     //theVar.push_back(leptWPt);  	  
     //theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     pTZ_vs_pTJets_DY->Fill(Zpt, leadJetPt);
     fillHisto("DataZjets", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
     
       fillHisto("BDTcut_DataZjets", theVar, ievt, EvtWeight);
       ndataEvents_Ztemplate++;
       
        if(Channel == 1) ndataEvents_Ztemplate_mumue++;
        if(Channel == 2) ndataEvents_Ztemplate_eemu++;
       
     }
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
  
  
   //if(syst == "DYDown"){
     
     //histBdt_DataZjets->Add(histBdt_DataZjets, histBdt_TTbarSig, 1, -0.05*histBdt_DataZjets->Integral()/histBdt_TTbarSig->Integral());
     
     
   //}
  
  
   double totexpectedTTbarZjets = 0;
  //-----------------------------------------------------
   //for TTbarSig jets from MC 
   //-----------------------------------------------------
   //define the tree to read
   TTree* theTree_TTbarSigZjets = (TTree*)input_zjets->Get("Ttree_TTbar");
   theTree_TTbarSigZjets->SetBranchAddress( "tree_totMass",	   &totMass);	
   theTree_TTbarSigZjets->SetBranchAddress( "tree_topMass",     &topMass );
   theTree_TTbarSigZjets->SetBranchAddress( "tree_deltaRlb",	   &deltaRlb);
   theTree_TTbarSigZjets->SetBranchAddress( "tree_deltaRTopZ",  &deltaRTopZ);
   theTree_TTbarSigZjets->SetBranchAddress( "tree_deltaPhilb",  &deltaPhilb );
   theTree_TTbarSigZjets->SetBranchAddress( "tree_asym",        &asym );
   theTree_TTbarSigZjets->SetBranchAddress( "tree_Zpt",         &Zpt );
   theTree_TTbarSigZjets->SetBranchAddress( "tree_ZEta",	   &ZEta);
   theTree_TTbarSigZjets->SetBranchAddress( "tree_topPt",	   &topPt);
   theTree_TTbarSigZjets->SetBranchAddress( "tree_topEta",	   &topEta); 
   theTree_TTbarSigZjets->SetBranchAddress( "tree_NJets",	   &NJets);
   theTree_TTbarSigZjets->SetBranchAddress( "tree_NBJets",	   &NBJets);
   theTree_TTbarSigZjets->SetBranchAddress( "tree_deltaRZl",   &deltaRZl);   
   theTree_TTbarSigZjets->SetBranchAddress( "tree_deltaPhiZmet",&deltaPhiZmet);
   theTree_TTbarSigZjets->SetBranchAddress( "tree_btagDiscri",  &btagDiscri);  
   
   
   theTree_TTbarSigZjets->SetBranchAddress("tree_leptWPt",       &leptWPt);	     
   theTree_TTbarSigZjets->SetBranchAddress("tree_leptWEta",      &leptWEta);	     
   theTree_TTbarSigZjets->SetBranchAddress("tree_leadJetPt",     &leadJetPt);	      
   theTree_TTbarSigZjets->SetBranchAddress("tree_leadJetEta",	&leadJetEta);	      
   //theTree_TTbarSigZjets->SetBranchAddress("tree_deltaRZleptW",	&deltaRZleptW);     
   theTree_TTbarSigZjets->SetBranchAddress("tree_deltaPhiZleptW", &deltaPhiZleptW);
   
   if(applyMWTcut)theTree_TTbarSigZjets->SetBranchAddress("tree_mTW", &mWT);	
   
   theTree_TTbarSigZjets->SetBranchAddress( "tree_EvtWeight",  &EvtWeight); 
   theTree_TTbarSigZjets->SetBranchAddress( "tree_Channel",  &Channel);    
   
   float nTTbarSigZjetsEvents_BCKenriched = 0;
   float nTTbarSigZjetsEvents_BCKenriched_unw = 0;
   //loop on the events and calculate the BDT output
   for (Long64_t ievt=0; ievt<theTree_TTbarSigZjets->GetEntries();ievt++) {    
     if (ievt%10 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree_TTbarSigZjets->GetEntry(ievt);
     
     if(chan != Channel && chan != -1) continue;
     if(mWT < mWTcut&& applyMWTcut) continue;
     totexpectedTTbarZjets+=EvtWeight;
     
     if(apply_BCK_cut    && reader->EvaluateMVA( "BDT") > bdt_BCK_cut) continue;
     if(apply_BCK_cut_up && reader->EvaluateMVA( "BDT") < bdt_BCK_cut_up) continue;
     if(EvtWeight !=EvtWeight)  continue; // check for nan
     histBdt_TTbarSigZjets->Fill( reader->EvaluateMVA( "BDT"),  EvtWeight);
     
     
     std::vector<double > theVar;
     theVar.push_back(topMass);
     ////theVar.push_back(totMass );
     theVar.push_back(deltaPhilb );
     ////theVar.push_back(deltaRlb);
     //theVar.push_back(deltaRTopZ);
     theVar.push_back(asym );
     theVar.push_back(Zpt );
     theVar.push_back(ZEta);
     ////theVar.push_back(topPt);
     theVar.push_back(topEta);
     theVar.push_back(NJets);
     theVar.push_back(NBJets);		   
     theVar.push_back(deltaRZl);   
     theVar.push_back(deltaPhiZmet);
     theVar.push_back(btagDiscri);
   
     ////theVar.push_back(leptWPt);  	  
     ////theVar.push_back(leptWEta); 	  
     theVar.push_back(leadJetPt);	    
     theVar.push_back(leadJetEta);	     
     theVar.push_back(deltaPhiZleptW);
     
     
     fillHisto("TTbarSigZjets", theVar, ievt, EvtWeight);
     if(reader->EvaluateMVA( "BDT") < bdt_BCK_cut){
       fillHisto("BDTcut_TTbarSigZjets", theVar, ievt, EvtWeight);
     }
     
     //to calculate the BDT output
     //reader->EvaluateMVA( "BDT"           ) 
     
   }
 
  
    
   
   target->cd();
   target->Write();
   
   
   histBdt_Data->Write();
   
   
  
   
   
   
   histBdt_WZ->Write();
   
   cout << "totexpectedZjets_data " << totexpectedZjets_data << endl;
   cout << "totexpectedZjets      " << totexpectedZjets      << endl;
   
   
   histBdt_Zjets->Write();
  
  
   histBdt_FCNC->Scale(0.1); 
   
   histBdt_FCNC->Write();
   cout << "number of signal events " << histBdt_FCNC->Integral() << endl;
   
   
   if(syst == "DYUp"){
      totexpectedZjets_data -= 2*totexpectedTTbarZjets;
      histBdt_DataZjets->Add(histBdt_DataZjets, histBdt_TTbarSigZjets, 1, -2);
      histBdt_DataZjets->Scale( totexpectedZjets/totexpectedZjets_data);
   }
   if(syst != "DYUp" && syst != "DYDown") {
      totexpectedZjets_data -= totexpectedTTbarZjets;
      histBdt_DataZjets->Add(histBdt_DataZjets, histBdt_TTbarSigZjets, 1, -1);
      histBdt_DataZjets->Scale( totexpectedZjets/totexpectedZjets_data);

   
   }
   
   histBdt_DataZjets->Write();
   
   
   
   
   pTZ_vs_pTJets_DY->Scale( totexpectedZjets/totexpectedZjets_data);
   pTZ_vs_pTJets_MC->Add(pTZ_vs_pTJets_MC, pTZ_vs_pTJets_DY, 1, 1);
     
   histBdt_TTbarSig->Write();	 
   histBdt_Wjets->Write();	 
   histBdt_TtW->Write();	 
   histBdt_TbartW->Write();   
   histBdt_TtChan ->Write();  
   histBdt_TbartChan->Write();
   histBdt_TsChan->Write();   
   histBdt_TbarsChan->Write();
   histBdt_TZq->Write();	 
   histBdt_ZZ->Write();	 
   histBdt_WW->Write();
   
   
   pTZ_vs_pTJets_MC->Write();  
   pTZ_vs_pTJets_Data->Write();
   
   
   target->Close();
   
}


void ReaderBDT(){

   TString thevertex_zut = "zut";
   TString thevertex_zct = "zct";
   TString thevertex_kut = "kut";
   TString thevertex_kct = "kct";
   TString thevertex_xut = "xut";
   TString thevertex_xct = "xct";
   TString thevertex_tzq = "tzq";
   
   
   /*
  ReaderBDT(thevertex_zut, 
    "../../SystRootFiles_2013_07_01/proof_nobtagcorr.root",
    "../../SystRootFiles_2013_07_01/proof_nobtagcorr.root",
     "nobtagcorr");
   */
   
   
	

  ReaderBDT(thevertex_zut, 
    "../../SystRootFiles_2013_07_01/proof_nom.root",
    "../../SystRootFiles_2013_07_01/proof_nom.root",
     "nom");

  ReaderBDT(thevertex_zct, 
    "../../SystRootFiles_2013_07_01/proof_nom.root",
    "../../SystRootFiles_2013_07_01/proof_nom.root",
     "nom");
    
  ReaderBDT(thevertex_kut, 
    "../../SystRootFiles_2013_07_01/proof_nom.root",
    "../../SystRootFiles_2013_07_01/proof_nom.root",
     "nom");
      
  ReaderBDT(thevertex_kct, 
    "../../SystRootFiles_2013_07_01/proof_nom.root",
    "../../SystRootFiles_2013_07_01/proof_nom.root",
     "nom");

	
   
   
   //******************************************
   // for Zut ac
   //******************************************
 
   
   //for Jes Up
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_JESup.root",
   	"../../SystRootFiles_2013_07_01/proof_JESup.root",
   	 "JESup");
	
   //for Jes down
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_JESdown.root",
   	"../../SystRootFiles_2013_07_01/proof_JESdown.root",
	 "JESdown");
	 
   //for JER 
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_JER.root",
   	"../../SystRootFiles_2013_07_01/proof_JER.root",
	 "JER");
	 
   //for LeptSF 
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_LeptSFup.root",
   	"../../SystRootFiles_2013_07_01/proof_LeptSFup.root",
	 "LeptSFup");
	 
   //for LeptSF 
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_LeptSFdown.root",
   	"../../SystRootFiles_2013_07_01/proof_LeptSFdown.root",
	 "LeptSFdown");
	 

   
   
   //for PU up
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_PUup.root",
   	"../../SystRootFiles_2013_07_01/proof_PUup.root",
	 "PUup");
   
	 
   //for PU down
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_PUdown.root",
   	"../../SystRootFiles_2013_07_01/proof_PUdown.root",
	 "PUdown");
   
   
   
   //for btag up
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_btagup.root",
   	"../../SystRootFiles_2013_07_01/proof_btagup.root",
	 "btagUp");
   
   
   
   //for btag down
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_btagdown.root",
   	"../../SystRootFiles_2013_07_01/proof_btagdown.root",
	 "btagdown");
  
      
  
   //for match up
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Matchup");
   
   
   //for match down
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Matchdown");
   
   
   //for scale up
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Scaleup");
   
   
   //for scale down
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Scaledown");
   
   //for top mass up
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Mtopup");
   
   
   //for top mass down
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Mtopdown");

   
   
	 
    ReaderBDT(thevertex_zut, 
   	//"../../SystRootFiles_2013_07_01/proof_PDFup.root",
   	//"../../SystRootFiles_2013_07_01/proof_PDFup.root",
   	"../../SystProofFiles_btadDiscri_update/proof_PDFup.root",
   	"../../SystProofFiles_btadDiscri_update/proof_PDFup.root",
	 "PDFup"); 
    
    ReaderBDT(thevertex_zut, 
   	//"../../SystRootFiles_2013_07_01/proof_PDFdown.root",
   	//"../../SystRootFiles_2013_07_01/proof_PDFdown.root",
   	"../../SystProofFiles_btadDiscri_update/proof_PDFdown.root",
   	"../../SystProofFiles_btadDiscri_update/proof_PDFdown.root",
	 "PDFdown");
	
	 
	 
   //for DY up
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYUp");
	 
   //for DY down
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYDown");
	
      //for DY up
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYenriched03");
	 
   //for DY down
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYenriched05");  
	 
   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "WZZptup"); 

   ReaderBDT(thevertex_zut, 
   	"../../SystRootFiles_2013_07_01/proof_ZptReweight.root",
   	"../../SystRootFiles_2013_07_01/proof_ZptReweight.root",
	 "WZZptdown"); 
	 

   //******************************************
   // for Zct ac
   //******************************************
 
   
   //for Jes Up
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_JESup.root",
   	"../../SystRootFiles_2013_07_01/proof_JESup.root",
   	 "JESup");
	
   //for Jes down
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_JESdown.root",
   	"../../SystRootFiles_2013_07_01/proof_JESdown.root",
	 "JESdown");
	 
   //for JER 
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_JER.root",
   	"../../SystRootFiles_2013_07_01/proof_JER.root",
	 "JER");
	 
   //for LeptSF 
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_LeptSFup.root",
   	"../../SystRootFiles_2013_07_01/proof_LeptSFup.root",
	 "LeptSFup");
	 
   //for LeptSF 
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_LeptSFdown.root",
   	"../../SystRootFiles_2013_07_01/proof_LeptSFdown.root",
	 "LeptSFdown");
	 

   
   
   //for PU up
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_PUup.root",
   	"../../SystRootFiles_2013_07_01/proof_PUup.root",
	 "PUup");
   
	 
   //for PU down
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_PUdown.root",
   	"../../SystRootFiles_2013_07_01/proof_PUdown.root",
	 "PUdown");
   
   
   
   //for btag up
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_btagup.root",
   	"../../SystRootFiles_2013_07_01/proof_btagup.root",
	 "btagUp");
   
   
   
   //for btag down
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_btagdown.root",
   	"../../SystRootFiles_2013_07_01/proof_btagdown.root",
	 "btagdown");
  
      
   
   //for match up
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Matchup");
   
   
   //for match down
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Matchdown");
   
   
   //for scale up
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Scaleup");
   
   
   //for scale down
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Scaledown");
   
   //for top mass up
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Mtopup");
   
   
   //for top mass down
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Mtopdown");

    
   
   //for DY up
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYUp");
	 
   //for DY down
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYDown");
	
      //for DY up
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYenriched03");
	 
   //for DY down
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYenriched05"); 
 
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_PDFup.root",
   	"../../SystRootFiles_2013_07_01/proof_PDFup.root",
	 "PDFup");
	 
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_PDFdown.root",
   	"../../SystRootFiles_2013_07_01/proof_PDFdown.root",
	 "PDFdown"); 
		 
   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "WZZptup"); 

   ReaderBDT(thevertex_zct, 
   	"../../SystRootFiles_2013_07_01/proof_ZptReweight.root",
   	"../../SystRootFiles_2013_07_01/proof_ZptReweight.root",
	 "WZZptdown"); 
	 

   //******************************************
   // for kut ac
   //******************************************
 
   
   //for Jes Up
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_JESup.root",
   	"../../SystRootFiles_2013_07_01/proof_JESup.root",
   	 "JESup");
	
   //for Jes down
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_JESdown.root",
   	"../../SystRootFiles_2013_07_01/proof_JESdown.root",
	 "JESdown");
	 
   //for JER 
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_JER.root",
   	"../../SystRootFiles_2013_07_01/proof_JER.root",
	 "JER");
	 
   //for LeptSF 
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_LeptSFup.root",
   	"../../SystRootFiles_2013_07_01/proof_LeptSFup.root",
	 "LeptSFup");
	 
   //for LeptSF 
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_LeptSFdown.root",
   	"../../SystRootFiles_2013_07_01/proof_LeptSFdown.root",
	 "LeptSFdown");
	 

   
   
   //for PU up
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_PUup.root",
   	"../../SystRootFiles_2013_07_01/proof_PUup.root",
	 "PUup");
   
	 
   //for PU down
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_PUdown.root",
   	"../../SystRootFiles_2013_07_01/proof_PUdown.root",
	 "PUdown");
   
   
   
   //for btag up
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_btagup.root",
   	"../../SystRootFiles_2013_07_01/proof_btagup.root",
	 "btagUp");
   
   
   
   //for btag down
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_btagdown.root",
   	"../../SystRootFiles_2013_07_01/proof_btagdown.root",
	 "btagdown");
  
      
   
   //for match up
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Matchup");
   
   
   //for match down
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Matchdown");
   
   
   //for scale up
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Scaleup");
   
   
   //for scale down
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Scaledown");
   
   //for top mass up
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Mtopup");
   
   
   //for top mass down
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Mtopdown");

   
   //for DY up
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYUp");
	 
   //for DY down
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYDown");
   	
   //for DY up
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYenriched03");
	 
   //for DY down
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYenriched05"); 
 	
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_PDFup.root",
   	"../../SystRootFiles_2013_07_01/proof_PDFup.root",
	 "PDFup");
	 
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_PDFdown.root",
   	"../../SystRootFiles_2013_07_01/proof_PDFdown.root",
	 "PDFdown");
	 
	
   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "WZZptup"); 

   ReaderBDT(thevertex_kut, 
   	"../../SystRootFiles_2013_07_01/proof_ZptReweight.root",
   	"../../SystRootFiles_2013_07_01/proof_ZptReweight.root",
	 "WZZptdown"); 
  
   //******************************************
   // for kct ac
   //******************************************
 
   
   //for Jes Up
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_JESup.root",
   	"../../SystRootFiles_2013_07_01/proof_JESup.root",
   	 "JESup");
	
   //for Jes down
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_JESdown.root",
   	"../../SystRootFiles_2013_07_01/proof_JESdown.root",
	 "JESdown");
	 
   //for JER 
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_JER.root",
   	"../../SystRootFiles_2013_07_01/proof_JER.root",
	 "JER");
	 
   //for LeptSF 
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_LeptSFup.root",
   	"../../SystRootFiles_2013_07_01/proof_LeptSFup.root",
	 "LeptSFup");
	 
   //for LeptSF 
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_LeptSFdown.root",
   	"../../SystRootFiles_2013_07_01/proof_LeptSFdown.root",
	 "LeptSFdown");
	 

   
   
   //for PU up
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_PUup.root",
   	"../../SystRootFiles_2013_07_01/proof_PUup.root",
	 "PUup");
   
	 
   //for PU down
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_PUdown.root",
   	"../../SystRootFiles_2013_07_01/proof_PUdown.root",
	 "PUdown");
   
   
   
   //for btag up
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_btagup.root",
   	"../../SystRootFiles_2013_07_01/proof_btagup.root",
	 "btagUp");
   
   
   
   //for btag down
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_btagdown.root",
   	"../../SystRootFiles_2013_07_01/proof_btagdown.root",
	 "btagdown");
  
      
   
   //for match up
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Matchup");
   
   
   //for match down
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Matchdown");
   
   
   //for scale up
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Scaleup");
   
   
   //for scale down
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Scaledown");
   
   //for top mass up
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Mtopup");
   
   
   //for top mass down
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "Mtopdown");

   
   
   //for DY up
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYUp");
	 
   //for DY down
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYDown");
	 
	
   
   //for DY up
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYenriched03");
	 
   //for DY down
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "DYenriched05"); 
	 
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_PDFup.root",
   	"../../SystRootFiles_2013_07_01///proof_PDFup.root",
	 "PDFup");
	 
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_PDFdown.root",
   	"../../SystRootFiles_2013_07_01/proof_PDFdown.root",
	 "PDFdown"); 
	
   
   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
   	"../../SystRootFiles_2013_07_01/proof_nom.root",
	 "WZZptup"); 

   ReaderBDT(thevertex_kct, 
   	"../../SystRootFiles_2013_07_01/proof_ZptReweight.root",
   	"../../SystRootFiles_2013_07_01/proof_ZptReweight.root",
	 "WZZptdown"); 
	 
	 

   
}
