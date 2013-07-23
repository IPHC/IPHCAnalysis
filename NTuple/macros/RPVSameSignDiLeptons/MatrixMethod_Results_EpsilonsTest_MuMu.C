MatrixMethod_Results_EpsilonsTest_MuMu(string selection, string eventType){


  float NjetsMin = -0.5;
  float NjetsMax = 5.5;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  TFile * file_Expected   = new TFile("MatrixMethod_OutPut_MuMu_MC_Fast_PullTest.root");
  TFile * file_Estimated  = new TFile("MatrixMethod_OutPut_MuMu_MC_Fast_EpsilonsTest.root");

  file_Expected->cd();
  TH1F* tight_TTDilept  = (TH1F*)gROOT->FindObject("MMExpected_TightMuMu_TTbarSignal");
  TH1F* medium_TTDilept = (TH1F*)gROOT->FindObject("MMExpected_MediumMuMu_TTbarSignal");
  TH1F* loose_TTDilept  = (TH1F*)gROOT->FindObject("MMExpected_LooseMuMu_TTbarSignal");

  tight_TTDilept->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  medium_TTDilept->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  loose_TTDilept->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);

  TH1F* tight_TTSemi  = (TH1F*)gROOT->FindObject("MMExpected_TightMuMu_TTbarSemileptonic");
  TH1F* medium_TTSemi = (TH1F*)gROOT->FindObject("MMExpected_MediumMuMu_TTbarSemileptonic");
  TH1F* loose_TTSemi  = (TH1F*)gROOT->FindObject("MMExpected_LooseMuMu_TTbarSemileptonic");

  tight_TTSemi->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  medium_TTSemi->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  loose_TTSemi->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);

  TH1F* tight_TTBckg  = (TH1F*)gROOT->FindObject("MMExpected_TightMuMu_TTbar");
  TH1F* medium_TTBckg = (TH1F*)gROOT->FindObject("MMExpected_MediumMuMu_TTbar");
  TH1F* loose_TTBckg  = (TH1F*)gROOT->FindObject("MMExpected_LooseMuMu_TTbar");

  tight_TTBckg->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  medium_TTBckg->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  loose_TTBckg->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);

  tight_TTDilept->SetFillColor(kRed+1);
  medium_TTDilept->SetFillColor(kRed+1);
  loose_TTDilept->SetFillColor(kRed+1);

  tight_TTSemi->SetFillColor(kRed-7);
  medium_TTSemi->SetFillColor(kRed-7);
  loose_TTSemi->SetFillColor(kRed-7);

  tight_TTBckg->SetFillColor(kRed-7);
  medium_TTBckg->SetFillColor(kRed-7);
  loose_TTBckg->SetFillColor(kRed-7);
  
  THStack* tightSignalLike = new THStack("tightSignalLike","tight signal like");
  tightSignalLike->Add(tight_TTDilept);
  THStack* mediumSignalLike = new THStack("mediumSignalLike","medium signal like");
  mediumSignalLike->Add(medium_TTDilept);
  THStack* looseSignalLike = new THStack( "looseSignalLike", "loose signal like");
  looseSignalLike->Add(loose_TTDilept);

  THStack* tightWLike = new THStack( "tightWLike"," tight W-like" );
  tightWLike->Add(tight_TTSemi);
  THStack* mediumWLike = new THStack("mediumWLike","medium W-like");
  mediumWLike->Add(medium_TTSemi);
  THStack* looseWLike = new THStack( "looseWLike", "loose W-like" );
  looseWLike->Add(loose_TTSemi);
  
  THStack* tightQCDLike = new THStack("tightQCDLike",  "tight QCD-like" );
  tightQCDLike->Add(tight_TTBckg);
  THStack* mediumQCDLike = new THStack("mediumQCDLike","medium QCD-like");
  mediumQCDLike->Add(medium_TTBckg);
  THStack* looseQCDLike = new THStack("looseQCDLike",  "loose QCD-like" );
  looseQCDLike->Add(loose_TTBckg);

  
  file_Estimated->cd();
  
  TH1F* estimate_Loose_QCD  = (TH1F*)gROOT->FindObject("MMEstimated_LooseMuMu_QCD" );
  TH1F* estimate_Tight_QCD  = (TH1F*)gROOT->FindObject("MMEstimated_TightMuMu_QCD" );
  TH1F* estimate_Medium_QCD = (TH1F*)gROOT->FindObject("MMEstimated_MediumMuMu_QCD");

  estimate_Loose_QCD->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  estimate_Tight_QCD->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  estimate_Medium_QCD->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);

  TH1F* estimate_Loose_WJets  = (TH1F*)gROOT->FindObject("MMEstimated_LooseMuMu_WJets" );
  TH1F* estimate_Tight_WJets  = (TH1F*)gROOT->FindObject("MMEstimated_TightMuMu_WJets" );
  TH1F* estimate_Medium_WJets = (TH1F*)gROOT->FindObject("MMEstimated_MediumMuMu_WJets");

  estimate_Loose_WJets->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  estimate_Tight_WJets->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  estimate_Medium_WJets->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);

  TH1F* estimate_Loose_Signal  = (TH1F*)gROOT->FindObject("MMEstimated_LooseMuMu_Signal" );
  TH1F* estimate_Tight_Signal  = (TH1F*)gROOT->FindObject("MMEstimated_TightMuMu_Signal" );
  TH1F* estimate_Medium_Signal = (TH1F*)gROOT->FindObject("MMEstimated_MediumMuMu_Signal");

  estimate_Loose_Signal->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  estimate_Tight_Signal->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
  estimate_Medium_Signal->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);

  title = new TLatex(10.,10.,"CMS preliminary");
  title->SetNDC();
  title->SetTextAlign(12);
  title->SetX(0.60);
  title->SetY(0.70);
  title->SetTextFont(42);
  title->SetTextSize(0.04);
  title->SetTextSizePixels(24);

    
  TLegend* qw = 0;
  qw = new TLegend(0.70,0.75,0.95,0.95);
  qw->SetFillColor(kWhite);

  if( selection == "Tight"){
    
    
    if(eventType == "Signal"){
      TCanvas * canv = new TCanvas(("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), ("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), 300, 400);
      estimate_Tight_Signal->SetMarkerStyle(20);
      estimate_Tight_Signal->GetXaxis()->SetTitle("Jet multiplicity");
      estimate_Tight_Signal->GetYaxis()->SetTitle("Number of events");
      estimate_Tight_Signal->GetYaxis()->SetLabelSize(0.025);
      estimate_Tight_Signal->SetTitle("Signal-like #mu#mu (tight)");
      estimate_Tight_Signal->Draw();
      tightSignalLike->Draw("hesame");
      estimate_Tight_Signal->Draw("same");
      TH1F * tmp  = new TH1F("tmp",  "tmp",  11, -0.5,  10.5);
      tmp->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
      tmp->Add(tmp, tight_TTDilept);
      TGraphErrors *thegraph = new TGraphErrors(tmp);
      thegraph->SetFillStyle(3005);
      thegraph->SetFillColor(1);
      thegraph->Draw("e2same");
      qw->AddEntry(tight_TTDilept, " t#bar{t} signal ","f");
      qw->AddEntry(estimate_Tight_Signal,    " Matrix Method  ","p");
      qw->Draw("same");
      title->Draw("same");
      canv->SaveAs("TightSignalMuMu_EpsilonsTest.pdf");
    }
    
    
    if(eventType == "Wlike"){
      TCanvas * canv = new TCanvas(("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), ("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), 300, 400);
      estimate_Tight_WJets->SetMarkerStyle(20);
      estimate_Tight_WJets->GetXaxis()->SetTitle("Jet multiplicity");
      estimate_Tight_WJets->GetYaxis()->SetTitle("Number of events");
      estimate_Tight_WJets->GetYaxis()->SetLabelSize(0.025);
      estimate_Tight_WJets->SetTitle("W+jets-like #mu#mu (tight)");
      estimate_Tight_WJets->Draw();
      tightWLike->Draw("hesame");
      estimate_Tight_WJets->Draw("same");
      TH1F * tmp  = new TH1F("tmp",  "tmp",  11, -0.5,  10.5);
      tmp->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
      tmp->Add(tmp, tight_TTSemi);
      TGraphErrors *thegraph = new TGraphErrors(tmp);
      thegraph->SetFillStyle(3005);
      thegraph->SetFillColor(1);
      thegraph->Draw("e2same");
      qw->AddEntry(tight_TTSemi, " t#bar{t} semileptonic ","f");
      qw->AddEntry(estimate_Tight_WJets,    " Matrix Method  ","p");
      qw->Draw("same");
      title->Draw("same");
      canv->SaveAs("TightWlikeMuMu_EpsilonsTest.pdf");
    }
    
    
    if(eventType == "QCDlike"){
      TCanvas * canv = new TCanvas(("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), ("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), 300, 400);
      estimate_Tight_QCD->SetMarkerStyle(20);
      estimate_Tight_QCD->GetXaxis()->SetTitle("Jet multiplicity");
      estimate_Tight_QCD->GetYaxis()->SetTitle("Number of events");
      estimate_Tight_QCD->GetYaxis()->SetLabelSize(0.025);
      estimate_Tight_QCD->SetTitle("QCD-like #mu#mu (tight)");
      estimate_Tight_QCD->Draw();
      tightQCDLike->Draw("hesame");
      estimate_Tight_QCD->Draw("same");
      TH1F * tmp  = new TH1F("tmp",  "tmp",  11, -0.5,  10.5);
      tmp->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
      tmp->Add(tmp, tight_TTBckg);
      TGraphErrors *thegraph = new TGraphErrors(tmp);
      thegraph->SetFillStyle(3005);
      thegraph->SetFillColor(1);
      estimate_Tight_QCD->Draw("same");
      thegraph->Draw("e2same");
      qw->AddEntry(tight_TTBckg, " t#bar{t} hadronic ","f");
      qw->AddEntry(estimate_Tight_QCD, " Matrix Method  ","p");
      qw->Draw("same");
      title->Draw("same");
      canv->SaveAs("TightQCDMuMu_EpsilonsTest.pdf");      
    }
    
    
  }
  
  if( selection == "Medium"){
    
    
    if(eventType == "Signal"){
      TCanvas * canv = new TCanvas(("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), ("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), 300, 400);
      estimate_Medium_Signal->SetMarkerStyle(20);
      estimate_Medium_Signal->GetXaxis()->SetTitle("Jet multiplicity");
      estimate_Medium_Signal->GetYaxis()->SetTitle("Number of events");
      estimate_Medium_Signal->GetYaxis()->SetLabelSize(0.025);
      estimate_Medium_Signal->SetTitle("Signal-like #mu#mu (medium)");
      estimate_Medium_Signal->Draw();
      mediumSignalLike->Draw("hesame");
      estimate_Medium_Signal->Draw("same");
      TH1F * tmp  = new TH1F("tmp",  "tmp",  11, -0.5,  10.5);
      tmp->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
      tmp->Add(tmp, medium_TTDilept);
      TGraphErrors *thegraph = new TGraphErrors(tmp);
      thegraph->SetFillStyle(3005);
      thegraph->SetFillColor(1);
      thegraph->Draw("e2same");
      qw->AddEntry(medium_TTDilept, " t#bar{t} signal ","f");
      qw->AddEntry(estimate_Medium_Signal,    " Matrix Method  ","p");
      qw->Draw("same");
      title->Draw("same");
      canv->SaveAs("MediumSignalMuMu_EpsilonsTest.pdf");
    }
    
    
    if(eventType == "Wlike"){
      TCanvas * canv = new TCanvas(("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), ("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), 300, 400);
      estimate_Medium_WJets->SetMarkerStyle(20);
      estimate_Medium_WJets->GetXaxis()->SetTitle("Jet multiplicity");
      estimate_Medium_WJets->GetYaxis()->SetTitle("Number of events");
      estimate_Medium_WJets->GetYaxis()->SetLabelSize(0.025);
      estimate_Medium_WJets->SetTitle("W+jets-like #mu#mu (medium)");
      estimate_Medium_WJets->Draw();
      mediumWLike->Draw("hesame");
      estimate_Medium_WJets->Draw("same");
      TH1F * tmp  = new TH1F("tmp",  "tmp",  11, -0.5,  10.5);
      tmp->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
      tmp->Add(tmp, medium_TTSemi);
      TGraphErrors *thegraph = new TGraphErrors(tmp);
      thegraph->SetFillStyle(3005);
      thegraph->SetFillColor(1);
      thegraph->Draw("e2same");
      qw->AddEntry(medium_TTSemi, " t#bar{t} semileptonic ","f");
      qw->AddEntry(estimate_Medium_WJets,    " Matrix Method  ","p");
      qw->Draw("same");
      title->Draw("same");
      canv->SaveAs("MediumWlikeMuMu_EpsilonsTest.pdf");      
    }
    
    
    if(eventType == "QCDlike"){
      TCanvas * canv = new TCanvas(("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), ("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), 300, 400);
      estimate_Medium_QCD->SetMarkerStyle(20);
      estimate_Medium_QCD->GetXaxis()->SetTitle("Jet multiplicity");
      estimate_Medium_QCD->GetYaxis()->SetTitle("Number of events");
      estimate_Medium_QCD->GetYaxis()->SetLabelSize(0.025);
      estimate_Medium_QCD->SetTitle("QCD-like #mu#mu (medium)");
      estimate_Medium_QCD->Draw();
      mediumQCDLike->Draw("hesame");
      estimate_Medium_QCD->Draw("same"); 
      TH1F * tmp  = new TH1F("tmp",  "tmp",  11, -0.5,  10.5);
      tmp->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
      tmp->Add(tmp, medium_TTBckg);
      TGraphErrors *thegraph = new TGraphErrors(tmp);
      thegraph->SetFillStyle(3005);
      thegraph->SetFillColor(1);
      thegraph->Draw("e2same");
      qw->AddEntry(medium_TTBckg, " t#bar{t} hadronic ","f");
      qw->AddEntry(estimate_Medium_QCD, " Matrix Method  ","p");
      qw->Draw("same");
      title->Draw("same");
      canv->SaveAs("MediumQCDMuMu_EpsilonsTest.pdf");      
    }
    
    
  }
  
  
  
  
  
  if( selection == "Loose"){
    
    
    if(eventType == "Signal"){
      TCanvas * canv = new TCanvas(("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), ("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), 300, 400);
      estimate_Loose_Signal->SetMarkerStyle(20);
      estimate_Loose_Signal->GetXaxis()->SetTitle("Jet multiplicity");
      estimate_Loose_Signal->GetYaxis()->SetTitle("Number of events");
      estimate_Loose_Signal->GetYaxis()->SetLabelSize(0.025);
      estimate_Loose_Signal->SetTitle("Signal-like #mu#mu (loose)");
      estimate_Loose_Signal->Draw();
      looseSignalLike->Draw("hesame");
      estimate_Loose_Signal->Draw("same");
      TH1F * tmp  = new TH1F("tmp",  "tmp",  11, -0.5,  10.5);
      tmp->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
      tmp->Add(tmp, loose_TTDilept);
      TGraphErrors *thegraph = new TGraphErrors(tmp);
      thegraph->SetFillStyle(3005);
      thegraph->SetFillColor(1);
      thegraph->Draw("e2same");
      qw->AddEntry(loose_TTDilept, " t#bar{t} signal ","f");
      qw->AddEntry(estimate_Loose_Signal,    " Matrix Method  ","p");
      qw->Draw("same");
      title->Draw("same");
      canv->SaveAs("LooseSignalMuMu_EpsilonsTest.pdf");      
    }
    
    
    if(eventType == "Wlike"){
      TCanvas * canv = new TCanvas(("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), ("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), 300, 400);
      estimate_Loose_WJets->SetMarkerStyle(20);
      estimate_Loose_WJets->GetXaxis()->SetTitle("Jet multiplicity");
      estimate_Loose_WJets->GetYaxis()->SetTitle("Number of events");
      estimate_Loose_WJets->GetYaxis()->SetLabelSize(0.025);
      estimate_Loose_WJets->SetTitle("W+jets-like #mu#mu (loose)");
      estimate_Loose_WJets->Draw();
      looseWLike->Draw("hesame");
      estimate_Loose_WJets->Draw("same");
      TH1F * tmp  = new TH1F("tmp",  "tmp",  11, -0.5,  10.5);
      tmp->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
      tmp->Add(tmp, loose_TTSemi);
      TGraphErrors *thegraph = new TGraphErrors(tmp);
      thegraph->SetFillStyle(3005);
      thegraph->SetFillColor(1);
      thegraph->Draw("e2same");
      qw->AddEntry(loose_TTSemi, " t#bar{t} semileptonic ","f");
      qw->AddEntry(estimate_Loose_WJets,    " Matrix Method  ","p");
      qw->Draw("same");
      title->Draw("same");
      canv->SaveAs("LooseWlikeMuMu_EpsilonsTest.pdf");      
    }
    
    
    if(eventType == "QCDlike"){
      TCanvas * canv = new TCanvas(("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), ("canvEpsilonsTest_MuMu_"+selection+eventType).c_str(), 300,400);
      estimate_Loose_QCD->SetMarkerStyle(20);
      estimate_Loose_QCD->GetXaxis()->SetTitle("Jet multiplicity");
      estimate_Loose_QCD->GetYaxis()->SetTitle("Number of events");
      estimate_Loose_QCD->GetYaxis()->SetLabelSize(0.025);
      estimate_Loose_QCD->SetTitle("QCD-like #mu#mu (loose)");
      estimate_Loose_QCD->Draw(); 
      looseQCDLike->Draw("hesame");
      estimate_Loose_QCD->Draw("same"); 
      TH1F * tmp  = new TH1F("tmp",  "tmp",  11, -0.5,  10.5);
      tmp->GetXaxis()->SetRangeUser(NjetsMin,NjetsMax);
      tmp->Add(tmp, loose_TTBckg);
      TGraphErrors *thegraph = new TGraphErrors(tmp);
      thegraph->SetFillStyle(3005);
      thegraph->SetFillColor(1);
      thegraph->Draw("e2same");
      qw->AddEntry(loose_TTBckg, " t#bar{t} hadroni ","f");
      qw->AddEntry(estimate_Loose_QCD, " Matrix Method  ","p");
      qw->Draw("same");
      title->Draw("same");
      canv->SaveAs("LooseQCDMuMu_EpsilonsTest.pdf");      
    }
    
    
  }
   
 
 
}
