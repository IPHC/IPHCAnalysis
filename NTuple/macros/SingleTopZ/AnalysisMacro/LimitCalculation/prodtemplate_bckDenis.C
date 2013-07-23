#include "TH1F.h"
#include "TFile.h"

#include "TString.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include <iostream>

void prod(){

// using namespace RooFit;
// using namespace RooStats;


//////////////////////////////////////////////////////////////////////////////////////////////////

  
   
  TString filename;
  filename="TMVApp.root"; //with add plots
  

 
  TFile * filechannel = new TFile(filename,"update");

  TH1F * histo_Data;
  TH1F * histo_Zjets;
  TH1F * histo_ZjetsMC;
  TH1F * histo_WZ ;
  TH1F * histo_Sig;


  std::cout<<"passage Data "<<std::endl;
  histo_Data = (TH1F*)filechannel->Get("MVA_BDT_Data");

  std::cout<<"passage Zjets "<<std::endl;
  histo_Zjets               = (TH1F*)filechannel->Get("MVA_BDT_DataZjets");
  histo_ZjetsMC             = (TH1F*)filechannel->Get("MVA_BDT_Zjets");
  
  
  histo_Zjets->Scale(histo_ZjetsMC->Integral()/histo_Zjets->Integral());
  
  
  std::cout<<"passage WZ "<<std::endl;
  histo_WZ                  = (TH1F*)filechannel->Get("MVA_BDT_WZ");
  histo_WZ->Scale(0.905713);  
  
  
  
  std::cout<<"passage Sig "<<std::endl;
  histo_Sig                  = (TH1F*)filechannel->Get("MVA_BDT_FCNCkut");


// Needed by Theta: correct empty bins of the model corresponding to non empty bins for the data (if not, Theta crashs)
      for (unsigned int k = 1; k < histo_Sig->GetNbinsX()+1; k++) {
        double b_s = histo_Sig->GetBinContent(k);
        double b_d = histo_Data->GetBinContent(k);
        double b_b1 = histo_Zjets->GetBinContent(k);
        double b_b2 = histo_WZ->GetBinContent(k);
	if ( b_s<=0 ) {
	    b_s = 0.000001;
            histo_Sig->SetBinContent(k,b_s);
	}
	if ( b_d<=0 ) {
	    b_d = 0.000001;
            histo_Data->SetBinContent(k,b_d);
	}
	if ( b_b1<=0 ) {
	    b_b1 = 0.000001;
            histo_Zjets->SetBinContent(k,b_b1);
	}
	if ( b_b2<=0 ) {
	    b_b2 = 0.000001;
            histo_WZ->SetBinContent(k,b_b2);
	}
      }





  TH1F *  histo_Sig_MinusOneSigma   = new TH1F ("histo_Sig_MinusOneSigma", "histo_Sig_MinusOneSigma", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );
  TH1F *  histo_Sig_PlusOneSigma   = new TH1F ("histo_Sig_PlusOneSigma", "histo_Sig_PlusOneSigma", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );

  TH1F *  histo_Zjets_MinusOneSigma   = new TH1F ("histo_Zjets_MinusOneSigma", "histo_Zjets_MinusOneSigma", histo_Zjets->GetXaxis()->GetNbins(),histo_Zjets->GetXaxis()->GetXmin(),       histo_Zjets->GetXaxis()->GetXmax() );
  TH1F *  histo_Zjets_PlusOneSigma   = new TH1F ("histo_Zjets_PlusOneSigma", "histo_Zjets_PlusOneSigma", histo_Zjets->GetXaxis()->GetNbins(),histo_Zjets->GetXaxis()->GetXmin(),       histo_Zjets->GetXaxis()->GetXmax() );

  TH1F *  histo_WZ_MinusOneSigma   = new TH1F ("histo_WZ_MinusOneSigma", "histo_WZ_MinusOneSigma", histo_WZ->GetXaxis()->GetNbins(),histo_WZ->GetXaxis()->GetXmin(),       histo_WZ->GetXaxis()->GetXmax() );
  TH1F *  histo_WZ_PlusOneSigma   = new TH1F ("histo_WZ_PlusOneSigma", "histo_WZ_PlusOneSigma", histo_WZ->GetXaxis()->GetNbins(),histo_WZ->GetXaxis()->GetXmin(),       histo_WZ->GetXaxis()->GetXmax() );


  
   

      
// Here I build an arbitrary distorded histo corresponding to a +/- one sigma systematic effect affecting the shape & normaliastion

      for (unsigned int k = 1; k < histo_Sig->GetNbinsX()+1; k++) {
        double shiftedbin_minus = 0;
	shiftedbin_minus = 0.999999*histo_Sig->GetBinContent(k);
        double shiftedbin_plus = 0;
	shiftedbin_plus =  1.000001*histo_Sig->GetBinContent(k);
        histo_Sig_MinusOneSigma->SetBinContent(k,shiftedbin_minus);
        histo_Sig_PlusOneSigma->SetBinContent(k,shiftedbin_plus);
      }
      for (unsigned int k = 1; k < histo_Zjets->GetNbinsX()+1; k++) {
        double shiftedbin_minus = 0;
	shiftedbin_minus = 0.999999*histo_Zjets->GetBinContent(k);
        double shiftedbin_plus = 0;
	shiftedbin_plus =  1.000001*histo_Zjets->GetBinContent(k);
        histo_Zjets_MinusOneSigma->SetBinContent(k,shiftedbin_minus);
        histo_Zjets_PlusOneSigma->SetBinContent(k,shiftedbin_plus);
      }
      for (unsigned int k = 1; k < histo_WZ->GetNbinsX()+1; k++) {
        double shiftedbin_minus = 0;
	shiftedbin_minus = 0.999999*histo_WZ->GetBinContent(k);
        double shiftedbin_plus = 0;
	shiftedbin_plus =  1.000001*histo_WZ->GetBinContent(k);
        histo_WZ_MinusOneSigma->SetBinContent(k,shiftedbin_minus);
        histo_WZ_PlusOneSigma->SetBinContent(k,shiftedbin_plus);
      }

       std::cout<<"norm Sig: "<<histo_Sig->Integral()<<" "<<histo_Sig_MinusOneSigma->Integral()<<" "<<histo_Sig_PlusOneSigma->Integral()<<std::endl;
       std::cout<<"norm Zjets: "<<histo_Zjets->Integral()<<" "<<histo_Zjets_MinusOneSigma->Integral()<<" "<<histo_Zjets_PlusOneSigma->Integral()<<std::endl;
       std::cout<<"norm WZ: "<<histo_WZ->Integral()<<" "<<histo_WZ_MinusOneSigma->Integral()<<" "<<histo_WZ_PlusOneSigma->Integral()<<std::endl;
       
       histo_Sig->SetName("MVABDT__FCNCkut62");
       histo_Data->SetName("MVABDT__DATA");
       histo_Zjets->SetName("MVABDT__DataZjets");
       histo_WZ->SetName("MVABDT__WZ");
       histo_Sig_MinusOneSigma->SetName("MVABDT__FCNCkut62__DumSignal__minus");
       histo_Sig_PlusOneSigma->SetName("MVABDT__FCNCkut62__DumSignal__plus");
       histo_Zjets_MinusOneSigma->SetName("MVABDT__DataZjets__DumZjets__minus");
       histo_Zjets_PlusOneSigma->SetName("MVABDT__DataZjets__DumZjets__plus");
       histo_WZ_MinusOneSigma->SetName("MVABDT__WZ__DumWZ__minus");
       histo_WZ_PlusOneSigma->SetName("MVABDT__WZ__DumWZ__plus");

// Add 2 additionnal points  
  
  
  TH1F *  histo_Sig_05   = new TH1F ("histo_Sig_05", "histo_Sig_05", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );
  TH1F *  histo_Sig_05_MinusOneSigma   = new TH1F ("histo_Sig_05_MinusOneSigma", "histo_Sig_05_MinusOneSigma", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );
  TH1F *  histo_Sig_05_PlusOneSigma   = new TH1F ("histo_Sig_05_PlusOneSigma", "histo_Sig_05_PlusOneSigma", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );
  
  
  TH1F *  histo_Sig_10   = new TH1F ("histo_Sig_10", "histo_Sig_10", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );
  TH1F *  histo_Sig_10_MinusOneSigma   = new TH1F ("histo_Sig_10_MinusOneSigma", "histo_Sig_10_MinusOneSigma", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );
  TH1F *  histo_Sig_10_PlusOneSigma   = new TH1F ("histo_Sig_10_PlusOneSigma", "histo_Sig_10_PlusOneSigma", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );

  
  TH1F *  histo_Sig_20   = new TH1F ("histo_Sig_20", "histo_Sig_20", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );
  TH1F *  histo_Sig_20_MinusOneSigma   = new TH1F ("histo_Sig_20_MinusOneSigma", "histo_Sig_20_MinusOneSigma", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );
  TH1F *  histo_Sig_20_PlusOneSigma   = new TH1F ("histo_Sig_20_PlusOneSigma", "histo_Sig_20_PlusOneSigma", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );


  TH1F *  histo_Sig_30   = new TH1F ("histo_Sig_30", "histo_Sig_30", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );
  TH1F *  histo_Sig_30_MinusOneSigma   = new TH1F ("histo_Sig_30_MinusOneSigma", "histo_Sig_30_MinusOneSigma", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );
  TH1F *  histo_Sig_30_PlusOneSigma   = new TH1F ("histo_Sig_30_PlusOneSigma", "histo_Sig_30_PlusOneSigma", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );

  TH1F *  histo_Sig_50   = new TH1F ("histo_Sig_50", "histo_Sig_50", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );
  TH1F *  histo_Sig_50_MinusOneSigma   = new TH1F ("histo_Sig_50_MinusOneSigma", "histo_Sig_50_MinusOneSigma", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );
  TH1F *  histo_Sig_50_PlusOneSigma   = new TH1F ("histo_Sig_50_PlusOneSigma", "histo_Sig_50_PlusOneSigma", histo_Sig->GetXaxis()->GetNbins(),histo_Sig->GetXaxis()->GetXmin(),       histo_Sig->GetXaxis()->GetXmax() );
  
      for (unsigned int k = 1; k < histo_Sig->GetNbinsX()+1; k++) {
        histo_Sig_05->SetBinContent(k,05.*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
        histo_Sig_05_MinusOneSigma->SetBinContent(k,05*0.99999*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
        histo_Sig_05_PlusOneSigma->SetBinContent(k,05*1.00001*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
	
        histo_Sig_10->SetBinContent(k,10.*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
        histo_Sig_10_MinusOneSigma->SetBinContent(k,10*0.99999*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
        histo_Sig_10_PlusOneSigma->SetBinContent(k,10*1.00001*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
	
        histo_Sig_20->SetBinContent(k,20.*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
        histo_Sig_20_MinusOneSigma->SetBinContent(k,20*0.99999*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
        histo_Sig_20_PlusOneSigma->SetBinContent(k,20*1.00001*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
	
        histo_Sig_30->SetBinContent(k,30.*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
        histo_Sig_30_MinusOneSigma->SetBinContent(k,30*0.99999*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
        histo_Sig_30_PlusOneSigma->SetBinContent(k,30*1.00001*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
	
        histo_Sig_50->SetBinContent(k,50.*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
        histo_Sig_50_MinusOneSigma->SetBinContent(k,50*0.99999*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
        histo_Sig_50_PlusOneSigma->SetBinContent(k,50*1.00001*histo_Sig->GetBinContent(k)/histo_Sig->Integral());
      }
  
       histo_Sig_05->SetName("MVABDT__FCNCkut05");
       histo_Sig_05_MinusOneSigma->SetName("MVABDT__FCNCkut05__DumSignal__minus");
       histo_Sig_05_PlusOneSigma->SetName("MVABDT__FCNCkut05__DumSignal__plus");
       
       histo_Sig_10->SetName("MVABDT__FCNCkut10");
       histo_Sig_10_MinusOneSigma->SetName("MVABDT__FCNCkut10__DumSignal__minus");
       histo_Sig_10_PlusOneSigma->SetName("MVABDT__FCNCkut10__DumSignal__plus");
       
       histo_Sig_20->SetName("MVABDT__FCNCkut20");
       histo_Sig_20_MinusOneSigma->SetName("MVABDT__FCNCkut20__DumSignal__minus");
       histo_Sig_20_PlusOneSigma->SetName("MVABDT__FCNCkut20__DumSignal__plus");
       
       histo_Sig_30->SetName("MVABDT__FCNCkut30");
       histo_Sig_30_MinusOneSigma->SetName("MVABDT__FCNCkut30__DumSignal__minus");
       histo_Sig_30_PlusOneSigma->SetName("MVABDT__FCNCkut30__DumSignal__plus");
       
       
       histo_Sig_50->SetName("MVABDT__FCNCkut50");
       histo_Sig_50_MinusOneSigma->SetName("MVABDT__FCNCkut50__DumSignal__minus");
       histo_Sig_50_PlusOneSigma->SetName("MVABDT__FCNCkut50__DumSignal__plus");
  

  TFile * filechannel1 = new TFile("NewFileToBeUsedForThetaWithAutoNamingConvention.root","new");
       histo_Sig->Write();
       histo_Data->Write();
       histo_Zjets->Write();
       histo_WZ->Write();
       histo_Sig_MinusOneSigma->Write();
       histo_Sig_PlusOneSigma->Write();
       histo_Zjets_MinusOneSigma->Write();
       histo_Zjets_PlusOneSigma->Write();
       histo_WZ_MinusOneSigma->Write();
       histo_WZ_PlusOneSigma->Write();

// add histo for making bands

       histo_Sig_05->Write();
       histo_Sig_05_MinusOneSigma->Write();
       histo_Sig_05_PlusOneSigma->Write();
       
       histo_Sig_10->Write();
       histo_Sig_10_MinusOneSigma->Write();
       histo_Sig_10_PlusOneSigma->Write();
       
       
       
       histo_Sig_20->Write();
       histo_Sig_20_MinusOneSigma->Write();
       histo_Sig_20_PlusOneSigma->Write();
       
       
       histo_Sig_30->Write();
       histo_Sig_30_MinusOneSigma->Write();
       histo_Sig_30_PlusOneSigma->Write();
       
       
       histo_Sig_50->Write();
       histo_Sig_50_MinusOneSigma->Write();
       histo_Sig_50_PlusOneSigma->Write();
       
       filechannel1->Close();  
       
       filechannel->Close();  

} 
