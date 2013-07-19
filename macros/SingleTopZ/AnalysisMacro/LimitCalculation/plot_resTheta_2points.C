#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TText.h>
#include <TString.h>

#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "Riostream.h"

////
////
//// Usage .L plot_resTheta.C+
//// plot()
//// voir les inputs a changer
//// pour les files (observed.txt et expected.txt) issus de theta, OTER LA PREMIERE LIGNE 
////


// square
inline double sqr (double arg) {return arg*arg;}



void plot(string Type, bool isZut)
{

// Type = "events" ("") or "param"

// Inputs: Number of points && range of plots in x,y 
 enum { Npoints = 3 };
 float Xmin = 0.;
 float Xmax = 60.;
 float Ymin = 0.;
 float Ymax = 6.0;
 if(isZut ) Ymax = 4.0;
 float sigma_init = 0.10; // section efficace initiale du MC
 float A = 0. ; 
 float Ninit = 0; 
 
 if(isZut){
 //for zut
  A = 2678000./(1000*1000);///(pow(2, 0.5)*1000.) ; // parametre du lagrangien du modele
  //A =  1.34; ///(pow(2, 0.5)*1000.) ; // parametre du lagrangien du modele
  Ninit = 57.0*4; // Nombre d'events signal MC initial =  MVABDT__FCNCkut57.Integral()
 }
 else{
   //for zct
    //A = 194400./(pow(2, 0.5)*1000.); // parametre du lagrangien du modele
    A =  194400; // parametre du lagrangien du modele
    Ninit = 56.6; // Nombre d'events signal MC initial =  MVABDT__FCNCkut57.Integral()
 }
 
 
 int stati=0;
 bool fit= 0;
 bool logx=0;
 
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//TCanvas *c1 = new TCanvas("c1", "plots",200,0,700,700);
TCanvas *c1 = new TCanvas("c1", "plots",200,0,800,700);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);
c1->SetGrid(1, 1);



TPad* pad1 = new TPad("pad1","This is pad1",0.02,0.02,0.98,0.98,21);

pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogx(logx);
// pad1->SetGrid();
   pad1->SetTopMargin(0.05);
   pad1->SetBottomMargin(0.15);
   pad1->SetRightMargin(0.04);
   pad1->SetLeftMargin(0.16);
pad1->SetGrid(1, 1);

gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleFont(42);
gStyle->SetTitleColor(1);
gStyle->SetTitleTextColor(1);
gStyle->SetTitleFillColor(10);
gStyle->SetTitleFontSize(0.05);
gStyle->SetTitleW(0.4);
gStyle->SetTitleH(0.09);
// gStyle->SetTitleX(0); // Set the position of the title box
// gStyle->SetTitleY(0.985); // Set the position of the title box
// gStyle->SetTitleStyle(Style_t style = 1001);
// gStyle->SetTitleBorderSize(2);
gStyle->SetOptStat(stati);
gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
// gStyle->SetPadGridX(true); gStyle->SetPadGridY(true);
gStyle->SetPadGridX(false); gStyle->SetPadGridY(false);

if (fit) {
gStyle->SetStatW(0.2);
gStyle->SetStatH(0.1);
gStyle->SetOptFit(111);
} else {
gStyle->SetStatW(0.3);
gStyle->SetStatH(0.2);
gStyle->SetOptFit(0);
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


 float x[Npoints], ex[Npoints], ex_low[Npoints], ex_high[Npoints];
 float resobs[Npoints], resobs_err[Npoints];
 float resexp[Npoints], resexp_err[Npoints], resexp_2sigma[Npoints], resexp_1sigma[Npoints], resexp_2sigma_low_err[Npoints], resexp_2sigma_high_err[Npoints],resexp_1sigma_low_err[Npoints], resexp_1sigma_high_err[Npoints];
 
 ifstream in_obs;
 if(isZut) in_obs.open("cls_limits_observed.txt");
 else      in_obs.open("cls_limits_observed.txt");
 
 
 int inc = 0;
 
  while (!in_obs.eof()) {
    in_obs >> x[inc] >> resobs[inc];
    ex[inc]  =  0., ex_low[inc]  =  0., ex_high[inc]  =  0., resobs_err[inc] = 0.;
    cout<<"Lecture : x,resobs "<<x[inc]<<" "<<resobs[inc]<<" "<<endl;
    inc++;
  }
  in_obs.close();
  
 ifstream in_exp;
 if(isZut) in_exp.open("cls_limits_expected.txt");
 else      in_exp.open("cls_limits_expected.txt");

 inc = 0;
 float dum;
 
  while (!in_exp.eof()) {
    in_exp >> dum >> resexp[inc] >> resexp_2sigma_low_err[inc] >> resexp_2sigma_high_err[inc] >> resexp_1sigma_low_err[inc] >> resexp_1sigma_high_err[inc];
    resexp_err[inc] = 0., resexp_2sigma[inc] = resexp[inc], resexp_1sigma[inc] = resexp[inc];
//    cout<<"Lecture : x,resexp "<<x[inc]<<" "<<resexp[inc]<<" "<<endl;
    inc++;
  }
  in_obs.close();
  
 
 for (int k=0; k<Npoints; k++) {
    resexp_2sigma_low_err[k] = resexp[k] - resexp_2sigma_low_err[k];
    resexp_2sigma_high_err[k] = resexp_2sigma_high_err[k] - resexp[k];
    resexp_1sigma_low_err[k] = resexp[k] - resexp_1sigma_low_err[k];
    resexp_1sigma_high_err[k] = resexp_1sigma_high_err[k] - resexp[k];
    cout<<" resexp_2sigma "<<resexp_2sigma_low_err[k]<<" "<<resexp_2sigma_high_err[k]<<endl;
 }
 
// Transformation nevents-> lambda (eventually) 
 if (Type=="param") {
    for (int k=0; k<Npoints; k++) {
       if ( k==0 ) Xmin = x[k]-0.1*x[k];
       x[k] = sqrt(  x[k]*sigma_init/(A*Ninit)); 
       if ( k==0 ) {
         
         Xmin = sqrt(  Xmin*sigma_init/(A*Ninit)); 
         Xmax = sqrt(  Xmax*sigma_init/(A*Ninit)); 
       }	 
       cout<<"k,x,xmin,xmax "<<k<<" "<<x[k]<<" "<<Xmin<<" "<<Xmax<<endl;
    }
 }
 
 
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  pad1->cd();

  TH2F* histo = new TH2F("histo","",100,Xmin,Xmax,100,Ymin,Ymax);

  histo->Draw(""); 
 if (Type=="param") {
  histo->SetLabelSize(0.022, "XYZ");
  } else {
  histo->SetLabelSize(0.04, "XYZ");
  }
 if (Type!="param") {
  histo->SetTitleSize(0.05, "XYZ"); 
  histo->GetXaxis()->SetTitle("Nevents");
  } else {
    histo->SetTitleSize(0.04, "XYZ"); 
    if(isZut) histo->GetXaxis()->SetTitle("#kappa_{gut}/#Lambda");
    else histo->GetXaxis()->SetTitle("#kappa_{gct}/#Lambda");
    //histo->GetXaxis()->SetTitle("Coupling");
  }
  histo->SetLabelFont(42, "XYZ"); 
  histo->SetTitleFont(42, "XYZ");
  histo->GetYaxis()->SetTitle("95% CL Limit on #sigma/#sigma_{pred.}");
  histo->SetTitleOffset(1.0,"X");
  histo->SetTitleOffset(0.6,"Y");
//   histo->SetNdivisions(509, "XYZ");
  histo->GetXaxis()->SetMoreLogLabels();
 if (Type!="param") {
   histo->GetXaxis()->SetNoExponent();
 }  

  TLine* Line1 = new TLine(Xmin,1,Xmax,1);

  //include the official CMS label
//   TPaveText* pt = new TPaveText(0.55,0.21,0.90,0.31,"brNDC");   
  TPaveText* pt = new TPaveText(0.55,0.85,0.90,0.95,"brNDC");   
  pt->SetBorderSize(0);   
  pt->SetFillStyle(0);  
  pt->SetTextAlign(13);   
  pt->SetTextFont(22);   
  pt->SetTextSize(0.03);   
  pt->AddText("CMS prelim. at 7 TeV, 4.9 fb^{-1}");   
  pt->Draw(""); 

 TGraphErrors* UL_obs = new TGraphErrors(Npoints,x,resobs,ex,resobs_err);
//   sf1stat->Draw("P"); 
//   sf1stat->SetFillStyle(3005);
//   sf1stat->SetFillColor(kBlue);
//   sf1stat->SetMarkerColor(kBlue);
//   sf1stat->SetLineColor(kBlue);
//   sf1stat->SetLineStyle(1);
   UL_obs->SetMarkerStyle(20);
//   sf1stat->SetMarkerSize(1.3);
   UL_obs->SetLineWidth(3);

 TGraphErrors* UL_exp = new TGraphErrors(Npoints,x,resexp,ex,resexp_err);
   UL_exp->SetMarkerStyle(21);
   UL_exp->SetLineWidth(3);
   //UL_exp->SetLineSize(0);
   UL_exp->SetMarkerSize(0);
   UL_exp->SetLineStyle(2);
  
  TGraphAsymmErrors* UL_exp_2sigma = new TGraphAsymmErrors(Npoints, x, resexp_2sigma, ex_low, ex_high, resexp_2sigma_low_err, resexp_2sigma_high_err);
  UL_exp_2sigma->SetMarkerStyle(21);

  TGraphAsymmErrors* UL_exp_1sigma = new TGraphAsymmErrors(Npoints, x, resexp_1sigma, ex_low, ex_high, resexp_1sigma_low_err, resexp_1sigma_high_err);
  UL_exp_1sigma->SetMarkerStyle(21);

  UL_exp_2sigma->Draw("3"); // for drawing bands
  UL_exp_2sigma->SetFillColor(kYellow);
  UL_exp_2sigma->SetLineStyle(2);
  UL_exp_1sigma->Draw("3same"); // for drawing bands
  UL_exp_1sigma->SetFillColor(kGreen);
  UL_exp_1sigma->SetLineStyle(2);
  UL_obs->Draw("PLsame"); 
  UL_exp->Draw("PLsame"); 


  TLegend* leg = new TLegend(0.60,0.55,0.87,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);   
//  leg->SetHeader("");
  leg->AddEntry(UL_obs,"Observed","PL");
  //leg->AddEntry(UL_exp,"Expeted","PL");
  leg->AddEntry(UL_exp_1sigma,"Expected #pm 1#sigma","FL");
  leg->AddEntry(UL_exp_2sigma,"Expected #pm 2#sigma","FL");
  leg->Draw();
  
   Line1->SetLineWidth(2);
   //Line1->SetLineColor(2);
  Line1->Draw();
  
  float UppLimitobs = 0.;
  float UppLimitexp = 0.;

 // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 //// look for the UL(observed/expected) in the values list
 if ( resobs[0]>1 && resobs[Npoints-1]<1){
    for (int k=0; k<Npoints-1; k++) {
      if ( resobs[k]>1 && resobs[k+1]<1){
         float a = (resobs[k+1]-resobs[k])/(x[k+1]-x[k]);
	 float b = resobs[k]- a*x[k];
	 UppLimitobs = (1.-b)/a;
	 cout<<"OBSERVED UPPER LIMIT VALUE = "<<UppLimitobs<<endl;
      }
    }
    
 } else {std::cout<<"Pas de limite observee!, changer le range de recherche "<<std::endl;}
  
 if ( resexp[0]>1 && resexp[Npoints-1]<1){
    for (int k=0; k<Npoints-1; k++) {
      if ( resexp[k]>1 && resexp[k+1]<1){
         float a = (resexp[k+1]-resexp[k])/(x[k+1]-x[k]);
	 float b = resexp[k]- a*x[k];
	 UppLimitexp = (1.-b)/a;
	 cout<<"EXPECTED UPPER LIMIT VALUE = "<<UppLimitexp<<endl;
      }
    }
    
 } else {std::cout<<"Pas de limite attendue!, changer le range de recherche "<<std::endl;}
  
  TLine* Line2 = new TLine(UppLimitobs,Ymin,UppLimitobs,Ymax*0.7);
  //Line2->Draw();
  Line2->SetLineWidth(3);
  Line2->SetLineColor(2);
  float gr1x[2],gr1y[2],egr1x[2],egr1y[2] ;
  gr1x[0] = UppLimitobs;
  gr1y[0] = Ymin;
  gr1x[1] = UppLimitobs;
  gr1y[1] = Ymax;
  egr1x[0] = 0.;
  egr1y[0] = 0.;
  egr1x[1] = 0.;
  egr1y[1] = 0.;
  TGraphErrors* dummy = new TGraphErrors(2,gr1x,gr1y,egr1x,egr1y);
  dummy->SetLineWidth(-2002);
  dummy->SetLineColor(2);
  dummy->SetFillStyle(3004);
  //dummy->Draw();

// to add code Higgs results (to be commented if not needed)
 /*dum = 0;
 float resobs_codehiggs[Npoints], resobs_err_codehiggs[Npoints];
 ifstream in_obs_codehiggs;
 in_obs_codehiggs.open("cls_limits_observed_codehiggs.txt");

 inc = 0;
  while (!in_obs_codehiggs.eof()) {
    in_obs_codehiggs >> dum >> resobs_codehiggs[inc];
     resobs_err_codehiggs[inc] = 0.;
    inc++;
  }
  in_obs_codehiggs.close();
 TGraphErrors* UL_obs_codehiggs = new TGraphErrors(Npoints,x,resobs_codehiggs,ex,resobs_err_codehiggs);
   UL_obs_codehiggs->SetMarkerStyle(32);
   UL_obs_codehiggs->SetMarkerColor(4);
   UL_obs_codehiggs->SetLineWidth(3);
   UL_obs_codehiggs->SetLineColor(4);
   UL_obs_codehiggs->SetLineStyle(2);
   UL_obs_codehiggs->Draw("PLsame"); 
   leg->AddEntry(UL_obs_codehiggs,"observed limit (code higgs)","PL");

// end to add code Higgs results 
*/

  c1->Update();

//$$  TString oname(title);
//$$  oname.Append(".pdf");
//$$  c1->SaveAs(oname);
}

