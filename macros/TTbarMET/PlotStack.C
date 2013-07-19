PlotStack(TString categorie, TString plotname, TString namechan, TString selection, bool setlogy, bool norma, bool ratio){
  
  TString channel;
  if (namechan!="all") channel=namechan;
  else channel="merged";

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  gStyle->SetOptDate(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleColor(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

// For the axis titles:

  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.25);

// For the axis labels:

  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
  
  
  TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
  c1->cd();
   
  TPad *canvas_1 = new TPad("canvas_1", "canvas_1",0,0.25,1.0,0.98);
  if (ratio) {
  canvas_1->Draw();
  canvas_1->cd(); 
  canvas_1->SetLogy(setlogy);
  }
  else {
  c1->SetLogy(setlogy);
  }
  
 
  TString file_tmp;
  file_tmp="../backup_outputProof16-07-12_12_data/proof.root";
  TFile * filedata = new TFile(file_tmp);

  TString filename;
  filename="../backup_outputProof13-07-12_fullmc/proof.root";
  TFile * filechannel = new TFile(filename);

  TH1F * histo_Data;
  TH1F * histo_DataEG;
  TH1F * histo_DataMu;

  TH1F * histo_Signal;
  TH1F * histo_TTbar;
  TH1F * histo_W1J;
  TH1F * histo_W2J;
  TH1F * histo_W3J;
  TH1F * histo_W4J;
  TH1F * histo_WJ;
  TH1F * histo_DY1;
  TH1F * histo_DY2;
  TH1F * histo_ZJ;
  TH1F * histo_WW;
  TH1F * histo_WZ;
  TH1F * histo_ZZ;
  TH1F * histo_VV ;
  TH1F * histo_SingleTop1;
  TH1F * histo_SingleTop2;
  TH1F * histo_SingleTop3;
  TH1F * histo_SingleTop4;
  TH1F * histo_SingleTop5;
  TH1F * histo_SingleTop6;
  TH1F * histo_STop;
  TH1F * histo_ratio;
  
  cout << "plotting "  << categorie << " " << channel << " " << selection << " " << plotname << endl;

  if (channel=="e"){
   TString histo_Data_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/DataEG/"+plotname;
   histo_Data              = (TH1F*)filedata->Get(histo_Data_name);
   cout << " histo_Data_name " << histo_Data_name << endl;
  }
  else if (channel=="mu") {
   TString histo_Data_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/DataMu/"+plotname;
   histo_Data              = (TH1F*)filedata->Get(histo_Data_name);
   cout << " histo_Data_name " << histo_Data_name << endl;
  }
  else if (channel=="merged") {
   TString histo_Data_name1 = "TTbarMetAnaPlots/"+categorie+"/e/"+selection+"/DataEG/"+plotname;
   histo_DataEG              = (TH1F*)filedata->Get(histo_Data_name1);
   cout << " histo_Data_name " << histo_Data_name1 << endl;
   TString histo_Data_name2 = "TTbarMetAnaPlots/"+categorie+"/mu/"+selection+"/DataMu/"+plotname;
   histo_DataMu              = (TH1F*)filedata->Get(histo_Data_name2);
   cout << " histo_Data_name " << histo_Data_name2 << endl;
   histo_Data               = (TH1F*)filedata->Get(histo_Data_name1);
   histo_Data->Add(histo_Data, histo_DataMu     , 1, 1);
  }

   TString histo_Signal_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/signal/"+plotname;
   histo_Signal              = (TH1F*)filechannel->Get(histo_Signal_name);
   cout << " histo_Signal_name " << histo_Signal_name << endl;

   TString histo_TTbar_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/ttbar/"+plotname;
   histo_TTbar              = (TH1F*)filechannel->Get(histo_TTbar_name);

   TString histo_W1J_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/W1Jet/"+plotname;
   histo_W1J              = (TH1F*)filechannel->Get(histo_W1J_name);
   TString histo_W2J_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/W2Jet/"+plotname;
   histo_W2J              = (TH1F*)filechannel->Get(histo_W2J_name);
   TString histo_W3J_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/W3Jet/"+plotname;
   histo_W3J              = (TH1F*)filechannel->Get(histo_W3J_name);
   TString histo_W4J_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/W4Jet/"+plotname;
   histo_W4J              = (TH1F*)filechannel->Get(histo_W4J_name);

   histo_WJ               = (TH1F*)filechannel->Get(histo_W1J_name);
   histo_WJ->Add(histo_WJ, histo_W2J     , 1, 1);
   histo_WJ->Add(histo_WJ, histo_W3J     , 1, 1);
   histo_WJ->Add(histo_WJ, histo_W4J     , 1, 1);


   TString histo_DY1_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/DY1/"+plotname;
   histo_DY1              = (TH1F*)filechannel->Get(histo_DY1_name);
   TString histo_DY2_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/DY2/"+plotname;
   histo_DY2              = (TH1F*)filechannel->Get(histo_DY2_name);

   histo_ZJ              = (TH1F*)filechannel->Get(histo_DY1_name);
   histo_ZJ->Add(histo_ZJ, histo_DY2     , 1, 1);

   TString histo_WW_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/WW/"+plotname;
   histo_WW              = (TH1F*)filechannel->Get(histo_WW_name);
   TString histo_WZ_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/WZ/"+plotname;
   histo_WZ              = (TH1F*)filechannel->Get(histo_WZ_name);
   TString histo_ZZ_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/ZZ/"+plotname;
   histo_ZZ              = (TH1F*)filechannel->Get(histo_ZZ_name);

   histo_VV              = (TH1F*)filechannel->Get(histo_WW_name);
   histo_VV->Add(histo_VV, histo_WZ     , 1, 1);
   histo_VV->Add(histo_VV, histo_ZZ     , 1, 1);

   TString histo_ST1_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/singleTop1/"+plotname;
   histo_SingleTop1              = (TH1F*)filechannel->Get(histo_ST1_name);
   TString histo_ST2_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/singleTop2/"+plotname;
   histo_SingleTop2              = (TH1F*)filechannel->Get(histo_ST2_name);
   TString histo_ST3_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/singleTop3/"+plotname;
   histo_SingleTop3              = (TH1F*)filechannel->Get(histo_ST3_name);
   TString histo_ST4_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/singleTop4/"+plotname;
   histo_SingleTop4              = (TH1F*)filechannel->Get(histo_ST4_name);
   TString histo_ST5_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/singleTop5/"+plotname;
   histo_SingleTop5              = (TH1F*)filechannel->Get(histo_ST5_name);
   TString histo_ST6_name = "TTbarMetAnaPlots/"+categorie+"/"+channel+"/"+selection+"/singleTop6/"+plotname;
   histo_SingleTop6              = (TH1F*)filechannel->Get(histo_ST6_name);

   histo_STop              = (TH1F*)filechannel->Get(histo_ST1_name);
   histo_STop->Add(histo_STop, histo_SingleTop2     , 1, 1);
   histo_STop->Add(histo_STop, histo_SingleTop3     , 1, 1);
   histo_STop->Add(histo_STop, histo_SingleTop4     , 1, 1);
   histo_STop->Add(histo_STop, histo_SingleTop5     , 1, 1);
   histo_STop->Add(histo_STop, histo_SingleTop6     , 1, 1);

    // hmc sans le signal mc
    TH1F * hmc= (TH1F*) histo_TTbar->Clone();
    hmc->Add(hmc,histo_WJ, 1., 1.);
    hmc->Add(hmc,histo_ZJ, 1., 1.);
    hmc->Add(hmc,histo_VV, 1., 1.);
    hmc->Add(hmc,histo_STop, 1., 1.);
    hmc->SetName("hmc");

    TH1F * hlinesig = (TH1F*) histo_Signal->Clone();

    histo_Signal->SetFillStyle(1001);
    histo_TTbar->SetFillStyle(1001);
    histo_WJ->SetFillStyle(1001);
    histo_ZJ->SetFillStyle(1001);
    histo_VV->SetFillStyle(1001);
    histo_STop->SetFillStyle(1001);
    
    histo_Signal->SetFillColor(kRed+1);
    histo_TTbar->SetFillColor(kRed-7);
    histo_WJ->SetFillColor(kGreen-3);
    histo_ZJ->SetFillColor(kAzure-2);
    histo_VV->SetFillColor(13);
    histo_STop->SetFillColor(kMagenta);
  


/*
    TH1F * lumiband = (TH1F*) hmc->Clone();
    for (int ilum=0; ilum<lumiband->GetNbinsX(); ilum++) {
     lumiband->SetBinError(ilum+1,lumiband->GetBinContent(ilum+1)*lumi_error);
    }
    TGraphErrors *thegraph = new TGraphErrors(lumiband);
    thegraph->SetFillStyle(3005);
    thegraph->SetFillColor(1);
*/

  
    THStack* hs= new THStack();
    hs->Add(histo_WJ);
    hs->Add(histo_ZJ);
    hs->Add(histo_VV);
    hs->Add(histo_STop);
    hs->Add(histo_TTbar);
    hs->Add(histo_Signal);
  
  
   if (!norma) {
    if (histo_Data->GetMaximum() > hs->GetMaximum()) hs->SetMaximum(histo_Data->GetMaximum()*1.2);
    if (plotname=="e_Eta") hs->SetMaximum(hs->GetMaximum()*1.2);
    hs->Draw("hist");
    hs->GetXaxis()->SetTitle(plotname);
    histo_Data->GetXaxis()->SetTitle("");
    histo_Data->GetYaxis()->SetTitle("");
    histo_Data->SetMarkerStyle(20);
    histo_Data->SetMarkerSize(0.7);
    histo_Data->Draw("epsame");
   }
   else {
    float integralmc=hmc->Integral();
    float integralsig=hlinesig->Integral();
    hmc->Scale(1./integralmc);
    hlinesig->Scale(1./integralsig);
    float maxi = hmc->GetMaximum();
    float max2 = hlinesig->GetMaximum();
    if (max2>maxi) hmc->SetMaximum(max2*1.2);
    hmc->SetLineWidth(2);
    hlinesig->SetLineColor(kRed+1);
    hlinesig->SetLineWidth(2);
    hmc->Draw("hist");
    hlinesig->Draw("same,hist");
   }

  
    text1 = new TLatex(0.35,0.95,"CMS Simu, 2.3 fb^{-1}");
    text1->SetNDC();
    text1->SetTextAlign(13);
    text1->SetX(0.35);
    text1->SetY(0.92);
    text1->SetTextFont(42);
    text1->SetTextSize(0.0610687);
    text1->Draw();

  TString info_data; 
  if (namechan=="e")   info_data = "e+jets channel";
  if (namechan=="mu") info_data = "#mu+jets channel";
  if (namechan=="all")  info_data = "l(e, #mu)+jets channel";

    text2 = new TLatex(0.15,0.8, info_data);
    text2->SetNDC();
    text2->SetTextAlign(13);
    text2->SetX(0.35);
    text2->SetY(0.85);
    text2->SetTextFont(42);
    text2->SetTextSize(0.0610687);
    text2->Draw();


  TLegend* qw = 0;
  if (!norma) {
  qw = new TLegend(0.75,0.70,0.98,0.98);
  qw->AddEntry(histo_Data,      "Data" ,                "p");
  qw->AddEntry(histo_WJ,        "W+jets  "                  ,"f");
  qw->AddEntry(histo_ZJ,        "Z+jets "                  ,"f");
  qw->AddEntry(histo_VV,        "diboson "                  ,"f");
  qw->AddEntry(histo_STop,       "single t"                  ,"f");
  qw->AddEntry(histo_TTbar,     "t#bar{t}  ","f");
  qw->AddEntry(histo_Signal,     "t#bar{t}+MET signal "     ,"f");
  }
  else {
  qw = new TLegend(0.75,0.75,0.98,0.98);
  qw->AddEntry(hlinesig,     "t#bar{t}+MET signal "     ,"l");
  qw->AddEntry(hmc,     "#sum backgrounds  ","l");
  }
  qw->SetFillColor(0);
  qw->SetTextFont(42);
  qw->Draw();

  c1->cd();

  if (ratio) {
  TPad *canvas_2 = new TPad("canvas_2", "canvas_2",0,0.,1.0,0.34);
  canvas_2->Draw();
  canvas_2->cd();
  gPad->SetBottomMargin(0.375);
  gPad->SetGridy();
  gPad->SetGridx();

  histo_ratio = (TH1F*) histo_Data->Clone();
  histo_ratio->SetName("histo_ratio");
  histo_ratio->SetTitle("");
  histo_ratio->Divide(hmc);

  histo_ratio->SetMarkerStyle(20);
  histo_ratio->SetMarkerSize(0.7);
  histo_ratio->GetXaxis()->SetTitle(plotname);
  histo_ratio->GetYaxis()->SetTitle("Data/SM MC");
  histo_ratio->GetYaxis()->SetTitleFont(42);
  histo_ratio->GetYaxis()->SetLabelFont(42);
  histo_ratio->GetXaxis()->SetTitleFont(42);
  histo_ratio->GetXaxis()->SetLabelFont(42);
  histo_ratio->GetXaxis()->SetLabelSize(0.1);

  histo_ratio->GetYaxis()->SetTitleOffset( 0.4 );
  histo_ratio->GetYaxis()->SetTitleSize( 0.1 );
  histo_ratio->GetYaxis()->SetLabelSize(0.1);
  histo_ratio->GetYaxis()->SetNdivisions( 505 );
  histo_ratio->GetXaxis()->SetTitleSize( 0.15 );

  histo_ratio->SetMinimum(0.);
  histo_ratio->SetMaximum(2.);
  histo_ratio->Draw("E1X0");

   c1->cd();


  }
 
  TString end_name;
  if(setlogy) end_name="_Logy.gif";
  else end_name=".gif"; 
  TString ratname;
  if (ratio) ratname="_r";
  else ratname="_r";
  TString aa;
  if (!norma)  aa="0";
  else aa="1";
  TString outputname= "plots/"+plotname+"_"+namechan+"_"+selection+aa+ratname+end_name;
  
  c1->SaveAs(outputname.Data());

/*
  histo_Signal->Delete();
  histo_TTbar->Delete();
  histo_ZJ->Delete();
  histo_DY1->Delete();
  histo_DY2->Delete();
  histo_WJ->Delete();
  histo_W1J->Delete();
  histo_W2J->Delete();
  histo_W3J->Delete();
  histo_W4J->Delete();
  histo_VV ->Delete();
  histo_WW ->Delete();
  histo_WZ ->Delete();
  histo_ZZ ->Delete();
  histo_SingleTop1->Delete();
  histo_SingleTop2->Delete();
  histo_SingleTop3->Delete();
  histo_SingleTop4->Delete();
  histo_SingleTop5->Delete();
  histo_SingleTop6->Delete();
  histo_STop->Delete();
  hs->Delete();
  hmc->Delete();
  if (ratio) {
   histo_ratio->Delete();   
  }
*/
/*
  thegraph->Delete();
*/


  filechannel->Close(); 
  filedata->Close(); 
}
