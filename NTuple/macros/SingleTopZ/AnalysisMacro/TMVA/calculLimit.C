#include "TFile.h"
#include "TH1F.h"
#include "TLimitDataSource.h"
#include "TLimit.h"
#include "TConfidenceLevel.h"
#include <iostream>





void calculLimit(){

  TFile* infile=new TFile("TMVA.root","READ");
  //infile->cd("Method_BDT/BDT");
  
  
  TH1F* sh=(TH1F*)infile->Get("Method_BDT/BDT/MVA_BDT_S_high");
  TH1F* bh=(TH1F*)infile->Get("Method_BDT/BDT/MVA_BDT_B_high");
 // TH1F* dh=(TH1F*)infile->Get("Method_BDT/BDT/MVA_BDT_S_high");
  
  
  
  
  
  //sh->Rebin(100);
  //bh->Rebin(100);
  //dh->Rebin(1);
  
  TH1F * dh = new TH1F("dh", "dh", sh->GetXaxis()->GetNbins(),sh->GetXaxis()->GetXmin() ,  sh->GetXaxis()->GetXmax());
  dh->Add(sh, bh, 1, 1);
  
  sh->SetLineColor(1);
  bh->SetLineColor(2);
  dh->SetLineColor(4);
  
  
  sh->Rebin(100);
  bh->Rebin(100);
  dh->Rebin(100);
  
  
  dh->Draw("h");
  sh->Draw("hsame");
  bh->Draw("hsame");
  
  
  TLimitDataSource* mydatasource = new TLimitDataSource(sh,bh,dh);
  TConfidenceLevel *myconfidence = TLimit::ComputeLimit(mydatasource,50000);
  cout << "  CLs    : " << myconfidence->CLs()  << endl;
  cout << "  CLsb   : " << myconfidence->CLsb() << endl;
  cout << "  CLb    : " << myconfidence->CLb()  << endl;
  cout << "< CLs >  : " << myconfidence->GetExpectedCLs_b()  << endl;
  cout << "< CLsb > : " << myconfidence->GetExpectedCLsb_b() << endl;
  cout << "< CLb >  : " << myconfidence->GetExpectedCLb_b()  << endl;
  delete myconfidence;
  delete mydatasource;



}
