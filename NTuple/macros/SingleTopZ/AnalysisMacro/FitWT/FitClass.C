#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <Cintex/Cintex.h>
#include <utility>
#include "TH1F.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TPaveText.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include <map>
#include <ctime>
#include <TKey.h>
#include <TRegexp.h>
#include <TLatex.h>
#include "Math/GSLRndmEngines.h"
#include "Math/Random.h"
#include "TTree.h"
#include "RooRandom.h"
#include <TLegend.h>


static const bool debug = true;


class TemplateFit {
public:

  TemplateFit();
  ~TemplateFit(){};
  


public:


  std::vector<std::vector<TString> > histoName_List_backgrounds;
  std::vector<TString >              histoName_List_signal;
  std::vector<TString >              histoName_List_data;
  
  std::vector< TString >             channel_List;
  std::vector< TString >             rootFileName_List_backgrounds;
  std::vector< TString >             rootFileName_List_signal;
  std::vector< TString >             rootFileName_List_data;
  
  std::vector< std::vector< TH1F *> > templates_backgrounds;
  std::vector< TH1F *>                templates_signal;
  std::vector< TH1F *>                histo_data;
  
  std::vector< bool >                isFloatingBckg;
  
  TString observable;
  TString selection;
  
  
  
  RooWorkspace *w;
  
  RooRealVar *variable, *fraction;
  RooCategory *sample;
  
  void configure();
  bool getInputs();
  void doFit();
  void displayResults();

};






TemplateFit::TemplateFit(){
  
  
  observable = "0";
  selection  = "0";
  
}





 

void TemplateFit::configure()
{
  
  
  //*********************************************************************
  // define 1st template for background 1
  std::vector<TString> histoName_List_backgrounds1;
  histoName_List_backgrounds1.push_back("DataMu");
  
  
  //*********************************************************************
  // define 2nd template for background 2
  std::vector<TString> histoName_List_backgrounds2;
  histoName_List_backgrounds2.push_back("TTbarSig");
  histoName_List_backgrounds2.push_back("TTbarBkg");
  histoName_List_backgrounds2.push_back("TtW");
  histoName_List_backgrounds2.push_back("TbartW");
  histoName_List_backgrounds2.push_back("WW");
  histoName_List_backgrounds2.push_back("ZZ");
  
  
  //*********************************************************************
  //push back of background template and define status (floating-fixed)
  histoName_List_backgrounds.push_back(histoName_List_backgrounds1);
  isFloatingBckg.push_back(true);
  histoName_List_backgrounds.push_back(histoName_List_backgrounds2);
  isFloatingBckg.push_back(false);
  
  
  //*********************************************************************
  // define signal template
  histoName_List_signal.push_back("WZ");
  //histoName_List_signal.push_back("WZ");
  
  
  //*********************************************************************
  // define data distribution
  histoName_List_data.push_back("DataMu");
  //histoName_List_data.push_back("DataMuEG");
  //histoName_List_data.push_back("DataMuEG");
  //histoName_List_data.push_back("DataEG");
  
  
  channel_List.push_back("mumumu");
  //channel_List.push_back("mumue");
  //channel_List.push_back("eemu");
  //channel_List.push_back("eee");
  
  
  
  
  rootFileName_List_backgrounds.push_back("proof_nom.root");
  rootFileName_List_backgrounds.push_back("proof_ZjetsTemplate.root");
  rootFileName_List_signal.push_back("proof_nom.root");
  rootFileName_List_data.push_back("proof_nom.root");
  
  
  
  observable = "mWT";
  selection  = "afterleptsel";
  
  
}


//######################################################################################################
//######################################################################################################
//######################################################################################################
//######################################################################################################
//######################################################################################################
//######################################################################################################

bool TemplateFit::getInputs(){
  
  if(histoName_List_data.size() != channel_List.size()){ 
    std::cout << "mismatch in the number of channels and the number of datasets. Please correct." << std::endl;
    return false;
  }
  
  
  if(rootFileName_List_backgrounds.size() != histoName_List_backgrounds.size()){
    std::cout << "mistmatch in the number of backgrounds template and the number of input files. There should be 1 file for 1 template " << std::endl;
    return false;
  }
  
  
  /*if(rootFileName_List_signal.size() != histoName_List_signal.size()){
    std::cout << "mistmatch in the number of signal template and the number of input files. There should be 1 file for 1 template " << std::endl;
    return false;
  }*/
  
  if(histoName_List_signal.size()  != channel_List.size()){    
    std::cout << "mistmatch in the number of channels and the number of signal samples. Please correct. " << std::endl;
    return false;
  }
  
  if(observable == "0"){
    std::cout << "observable not define " << endl; return false;
  }
  if(selection == "0"){
    std::cout << "selection not define " << endl; return false;  
  }
  
  for(unsigned int ichan = 0; ichan<channel_List.size(); ichan++){
    
    std::cout << "get histograms for channel " << channel_List[ichan] << " with the corresponding dataset being " << histoName_List_data[ichan] << endl;
    
    //*********************************************************************
    // get data histograms
    //*********************************************************************
    if(debug) cout << "roof file for data " << rootFileName_List_data[0] << endl;
    TFile *f_data   = new TFile(rootFileName_List_data[0].Data());
    f_data->cd();
    std::string histoname = (observable+"_"+channel_List[ichan]+"_"+selection+"_"+histoName_List_data[ichan]).Data();
    std::cout << "Loading histo for Data '" << histoname << "' ... " << std::endl;
    histo_data.push_back(dynamic_cast<TH1F*>( gROOT->FindObject(histoname.c_str()) ) );
    if(histo_data.back()==0) {std::cout << "ERROR : plot not found'" << std::endl; return false; }
    
    histo_data[ichan]->Draw();
    
    
    
      
    //*********************************************************************
    // get signal histograms
    //*********************************************************************
    if(debug) cout << "root file for signal " << rootFileName_List_signal[0] << endl;
    TFile *f_signal   = new TFile(rootFileName_List_signal[0].Data());
    cout << "line 231 " << endl;
    f_signal->cd();
    cout << "line 233 " << endl;
    histoname = (observable+"_"+channel_List[ichan]+"_"+selection+"_"+histoName_List_signal[ichan]).Data();
    cout << "line 235 " << endl;
    std::cout << "Loading template for signal '" << histoname << "' ... " << std::endl;
    templates_signal.push_back(dynamic_cast<TH1F*>( gROOT->FindObject(histoname.c_str()) ) );
    if(templates_signal.back()==0) {std::cout << "ERROR : plot not found'" << std::endl; return false; }
    
    templates_signal[ichan]->Draw("same");
    
    
  
    
    
    
    
  
    
    
    //*********************************************************************
    // get background templates
    //*********************************************************************
    
    std::vector< TH1F* > tmptemplates_backgrounds;
    
    for(unsigned itemp=0; itemp<histoName_List_backgrounds.size(); itemp++){
       
       TFile *f_temp   = new TFile(rootFileName_List_backgrounds[0].Data());
       f_temp->cd();
       
       
       std::vector<TH1F*> histo_List;

       for(unsigned int iback=0; iback<histoName_List_backgrounds[itemp].size(); iback++){
         
	   histoname = (observable+"_"+channel_List[ichan]+"_"+selection+"_"+(histoName_List_backgrounds[itemp])[iback]).Data();
           std::cout << "Loading template for background nb " << itemp << " with name " << histoname << std::endl;
           histo_List.push_back(dynamic_cast<TH1F*>( gROOT->FindObject(histoname.c_str()) ));
           if (histo_List.back()==0) {std::cout << "ERROR : plot not found'" << std::endl; return false; }
	   
       }
       
       TH1F *sumHisto = new TH1F( ("sumHisto_"+channel_List[ichan]+"_"+histoName_List_backgrounds[itemp].front()).Data() ,
                                  ("sumHisto_"+channel_List[ichan]+"_"+histoName_List_backgrounds[itemp].front()).Data() ,
                                    histo_List.front()->GetNbinsX(), 
                                    histo_List.front()->GetXaxis()->GetXmin(), 
                                    histo_List.front()->GetXaxis()->GetXmax() );
       
      for(unsigned int i=0; i<histo_List.size(); i++)
       {
         histo_List[i]->Sumw2();
         sumHisto->Add( sumHisto, histo_List[i], 1, 1);
       }
       tmptemplates_backgrounds.push_back(sumHisto);
       sumHisto->SetLineColor(itemp);
       //sumHisto->Draw("lsame");
    }
    templates_backgrounds.push_back(tmptemplates_backgrounds);
  }//end of loop on channel
  
  return true;
}






//######################################################################################################
//######################################################################################################
//######################################################################################################
//######################################################################################################
//######################################################################################################
//######################################################################################################


void TemplateFit::doFit(){
  
  //---------------------------
  //define basics for workspace
  //---------------------------
  w = new RooWorkspace("w","workspace") ;
  variable = new RooRealVar("variable", "variable", histo_data.front()->GetXaxis()->GetXmin(), histo_data.front()->GetXaxis()->GetXmax());
  sample = new RooCategory("sample","sample") ;
  for (unsigned int itype = 0;itype<channel_List.size();itype++) sample->defineType(channel_List[itype]);
  
  w->import(*variable); 
  
  
  w->import(*sample);
  
  //--------------------------------
  //get templates and construct pdfs
  //--------------------------------
  for (unsigned int itype = 0;itype<channel_List.size();itype++){
      if(debug) cout << "************for the channel " << channel_List[itype] << endl;
      
      
      //***************
      //for backgrounds
      for(unsigned int ibackg=0; ibackg<templates_backgrounds.size(); ibackg++){
        
	
        if(debug) cout << "   ********* for the bckg template " << (histoName_List_backgrounds[ibackg])[0] << endl;
	
	//get 
	TString hname = "histo_bkg_"+(histoName_List_backgrounds[ibackg])[0]+"_"+channel_List[itype];
	if (debug) cout << "hname  " << hname  << endl;
        w->import(*(new RooDataHist(hname.Data(), hname.Data(), *variable, (templates_backgrounds[itype])[ibackg] )));
	TString hname2 = "HistPdf::bkg_"+channel_List[itype]+"_"+(histoName_List_backgrounds[ibackg])[0]+"(variable,"+ hname+")";
        w->factory(hname2.Data());
	if (debug) cout << "hname2 " << hname2 << endl;
	
	
      }
      
      
      
      //***************
      //for signal 
      if(debug) cout << "   ********* for the signal template " << histoName_List_signal[itype] << endl;
	
      TString hname = "histo_sig_"+histoName_List_signal[itype]+"_"+channel_List[itype];
      if (debug) cout << "hname  " << hname  << endl;
      w->import(*(new RooDataHist(hname.Data(), hname.Data(), *variable,  histoName_List_signal[itype] )));
      TString hname2 = "HistPdf::sig_"+channel_List[itype]+"_"+histoName_List_signal[itype]+"(variable,"+ hname+")";
      w->factory(hname2.Data());
      if (debug) cout << "hname2 " << hname2 << endl;
      
      
      
      //***********************
      //define model for signal
      /*for(unsigned int ibackg=1; ibackg<templates_backgrounds.size(); ibackg++){
        
	TString hname = "";
	
	
      }*/
    
      
      
      //*****************************
      //define sum pdf for bakgrounds
      for(unsigned int ibackg=0; ibackg<templates_backgrounds[itype].size()-1; ibackg++){
        
	
	
	if(ibackg == 0){
          //***************************************
	  //sum pdf for signal and first background
	  TString factionname = "fraction_"+channel_List[itype]+"_"+(histoName_List_backgrounds[ibackg])[0];
	  if(debug) cout << "the fraction name for iteration " << ibackg << " is " << factionname << endl;
          fraction = new RooRealVar(factionname.Data(), factionname.Data(), 0.5, 0., 1.);
          w->import(*fraction);
	  TString hnamebck = "";
	  TString back1  = "bkg_"+channel_List[itype]+"_"+(histoName_List_backgrounds[ibackg])[0];
	  TString signal = "sig_"+channel_List[itype]+"_"+histoName_List_signal[itype];
          hnamebck = "SUM::thesum0_"+channel_List[itype]+"("+factionname.Data()+"*"+signal+", "+back1+" )";
	  if(debug) cout << "the SUM pdf name for iteration " << ibackg << " is " << hnamebck << endl;
	  w->factory(hnamebck.Data());
	  
	}else{
          //******************************************
	  //sum pdf for SUM::pdf-1 and next background
	  
	  char frac[50] = "";
	  sprintf(frac, "fraction%i_", ibackg);
	  
	  TString factionname = frac+channel_List[itype]+"_"+(histoName_List_backgrounds[ibackg])[0];
	  if(debug) cout << "the fraction name for iteration " << ibackg << " is " << factionname << endl;
          fraction = new RooRealVar(factionname.Data(), factionname.Data(), 0.5, 0., 1.);
          w->import(*fraction);
	  
	   char back[50]        = "";
	   char sumpdf[50]      = "";
	   char sumpdf_prev[50] = "";
	  
	  sprintf(back,        "bkg%i_",    ibackg  );
	  sprintf(sumpdf,      "thesum%i_", ibackg  );
	  sprintf(sumpdf_prev, "thesum%i_", ibackg-1);
	  
	  
	  TString hnamebck     = "";
	  TString hnamebck_tmp = "";
	  TString back1  = back+channel_List[itype]+"_"+(histoName_List_backgrounds[ibackg])[0];
	  hnamebck_tmp = sumpdf+channel_List[itype]+"("+factionname+"*"+sumpdf_prev+", "+back1+" )";
          hnamebck = "SUM::"+hnamebck_tmp;
	  
	  if(debug) cout << "the SUM pdf name for iteration " << ibackg << " is " << hnamebck << endl;
	  w->factory(hnamebck.Data());

	  
	  
	}
      }
      
      //**************************
      //define model 
      if( (isFloatingBckg.back()) ){
        //if last background is floating
	// add another ratio
        int ibackg = (templates_backgrounds[itype].size()-1);
	char frac[50] = "";
	sprintf(frac, "fraction%i_", ibackg);
	
	TString factionname = frac+channel_List[itype]+"_"+(histoName_List_backgrounds[ibackg])[0];
	if(debug) cout << "the fraction name for iteration " << ibackg << " is " << factionname << endl;
        fraction = new RooRealVar(factionname.Data(), factionname.Data(), 0.5, 0., 1.);
        w->import(*fraction);
	
	 char back[50]        = "";
	 char sumpdf[50]      = "";
	 char sumpdf_prev[50] = "";
	
	sprintf(back,	     "bkg%i_",    ibackg  );
	sprintf(sumpdf,      "thesum%i_", ibackg  );
	sprintf(sumpdf_prev, "thesum%i_", ibackg-1);
	
	
	TString hnamebck     = "";
	TString hnamebck_tmp = "";
	TString back1  = back+channel_List[itype]+"_"+(histoName_List_backgrounds[ibackg])[0];
	hnamebck_tmp = sumpdf+channel_List[itype]+"("+factionname+"*"+sumpdf_prev+", "+back1+" )";
        hnamebck = "SUM::"+hnamebck_tmp;
	
	if(debug) cout << "the SUM pdf name for iteration " << ibackg << " is " << hnamebck << endl;
	w->factory(hnamebck.Data());
	
      }else{
        //if last background is fixed
	//fit the number of events
	
        //if last background is floating
	// add another ratio
        int ibackg = (templates_backgrounds[itype].size()-1);
	char frac[50] = "";
	sprintf(frac, "fraction%i_", ibackg);
	
	TString factionname = frac+channel_List[itype]+"_"+(histoName_List_backgrounds[ibackg])[0];
	if(debug) cout << "the fraction name for iteration " << ibackg << " is " << factionname << endl;
        fraction = new RooRealVar(factionname.Data(), factionname.Data(), 0.5, 0., 1.);
        w->import(*fraction);
	
	 char back[50]        = "";
	 char sumpdf[50]      = "";
	 char sumpdf_prev[50] = "";
	
	sprintf(back,	     "bkg%i_",    ibackg  );
	sprintf(sumpdf,      "thesum%i_", ibackg  );
	sprintf(sumpdf_prev, "thesum%i_", ibackg-1);
	
	
	TString hnamebck     = "";
	TString hnamebck_tmp = "";
	TString back1  = back+channel_List[itype]+"_"+(histoName_List_backgrounds[ibackg])[0];
	hnamebck_tmp = sumpdf+channel_List[itype]+"("+factionname+"*"+sumpdf_prev+", "+back1+" )";
        hnamebck = "SUM::model_"channel_List[itype]+"("+;
	
	if(debug) cout << "the SUM pdf name for iteration " << ibackg << " is " << hnamebck << endl;
	w->factory(hnamebck.Data());
      }
      
      
  }
  

  
  
  //--------------------------------
  //define models
  //--------------------------------
  
   
  
  
  
  
}

void TemplateFit::displayResults(){

}

