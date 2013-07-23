
#include "TLimitDataSource.h"
#include "TConfidenceLevel.h"
#include "TLimit.h"
#include "TObjString.h"
#include "TVectorD.h"
#include "TH1D.h"
#include <iostream>

using namespace std;

void TLimitsExp()
{

  // Data no b-tag :

  //double fac = 1000./3.1;

  //  double _nmc   = 1.2*fac;
  //  double _ndata = 1*fac;  
  //  double _nsig  = 150;

  double _nmc   = 210;//.17.7;
  double _nsig  = 170;
  
  double _ndata = 390;  

  TH1D sig1("sig1","expected signal" ,1,0,1);
  TH1D bck1("bck1","expected bkg"    ,1,0,1);
  TH1D dat1("dat1","data"            ,1,0,1);

  TVectorD errorb(1);
  TVectorD errors(1);

  TObjArray names(1);
  TObjString name1("Total error");

  names.AddLast(&name1);

  errorb[0] = 0.1;
  errors[0] = 0.1;

  sig1.SetBinContent(1,_nsig);
  bck1.SetBinContent(1,_nmc);
  dat1.SetBinContent(1,_ndata);

  TLimitDataSource* mynewdatasource  = new TLimitDataSource();
  mynewdatasource->AddChannel(&sig1,&bck1,&dat1,&errors,&errorb,&names);
  //mynewdatasource->AddChannel(&sig1,&bck1,&dat1);
  TConfidenceLevel *mynewconfidence = TLimit::ComputeLimit(mynewdatasource,200000,false);

  double _clb      = (Float_t)mynewconfidence->CLb();
  double _cls_exp  = (Float_t)mynewconfidence->GetExpectedCLs_b();
  double _cls_obt  = (Float_t)mynewconfidence->CLs();
  double _cls_obt2 = (Float_t)mynewconfidence->CLs(kTRUE);

  cout << "CLb = " << _clb << endl;
  cout << "Expect    1-CLs = " << _cls_exp << endl;
  cout << "Observed  1-CLs = " << _cls_obt << endl;
  cout << "Observed2 1-CLs = " << _cls_obt2 << endl;


}
