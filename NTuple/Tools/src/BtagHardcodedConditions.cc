// -*- C++ -*-
//
// Helper class, provides jet tagging eff, scale factors, etc.
//
// Normally, you would want to get this via official conditions
// mechanisms (event setup). At the moment, there is a memory
// leak somewhere in core software when one gets some conditions
// from fwlite event setup from a root file.
//
// This class is hopefully a temporary solution for the memory leak.
//
// Gena Kukartsev, November 2012
//


#include <cmath>
#include "../interface/BtagHardcodedConditions.h"

using namespace std;
BtagHardcodedConditions::BtagHardcodedConditions() {
  float SFb_TCHPT_temp[14] = { 0.0543376, 0.0534339, 0.0266156, 0.0271337, 0.0276364, 0.0308838, 0.0381656, 0.0336979, 0.0336773, 0.0347688, 0.0376865, 0.0556052, 0.0598105, 0.0861122 };
  float SFb_CSVL_temp[14] = { 0.0188743, 0.0161816, 0.0139824, 0.0152644, 0.0161226, 0.0157396, 0.0161619, 0.0168747, 0.0257175, 0.026424, 0.0264928, 0.0315127, 0.030734, 0.0438259 };
  float SFb_CSVM_temp[14] = { 0.0295675, 0.0295095, 0.0210867, 0.0219349, 0.0227033, 0.0204062, 0.0185857, 0.0256242, 0.0383341, 0.0409675, 0.0420284, 0.0541299, 0.0578761, 0.0655432 };
  float SFb_CSVT_temp[14] = { 0.0364717, 0.0362281, 0.0232876, 0.0249618, 0.0261482, 0.0290466, 0.0300033, 0.0453252, 0.0685143, 0.0653621, 0.0712586, 0.094589, 0.0777011, 0.0866563 };
  float SFb_JPL_temp[14] = { 0.0250319, 0.0250197, 0.0212994, 0.0225867, 0.0239025, 0.026476, 0.0264219, 0.0156582, 0.0222798, 0.0223169, 0.0225454, 0.0405975, 0.0405668, 0.0415829 };
  float SFb_JPM_temp[14] = { 0.0352594, 0.0353008, 0.0299008, 0.0276606, 0.0292312, 0.0336607, 0.0284701, 0.029544, 0.0358872, 0.0367869, 0.0375048, 0.0597367, 0.0653152, 0.074242 };
  float SFb_JPT_temp[14] = { 0.0475813, 0.0472359, 0.0378328, 0.0334787, 0.034681, 0.0398312, 0.0481646, 0.0392262, 0.0463086, 0.0534565, 0.0545823, 0.102505, 0.113198, 0.138116 };
  fillArray(SFb_TCHPT_error, SFb_TCHPT_temp,14);
  fillArray(SFb_CSVL_error, SFb_CSVL_temp,14);
  fillArray(SFb_CSVM_error, SFb_CSVM_temp,14);
  fillArray(SFb_CSVT_error, SFb_CSVT_temp,14);
  fillArray(SFb_JPL_error, SFb_JPL_temp,14);
  fillArray(SFb_JPM_error, SFb_JPM_temp,14);
  fillArray(SFb_JPT_error, SFb_JPT_temp,14);
  float ptminT[14] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
  for (int i=0;i<14;++i) ptRange.push_back(ptminT[i]);
}


BtagHardcodedConditions::~BtagHardcodedConditions() {
}

std::string BtagHardcodedConditions::getAlgoName(const std::string & op){
  if( op == "JPL")   	  return "jetProbabilityBJetTags";
  else if( op == "JPM")   return "jetProbabilityBJetTags";
  else if( op == "JPT")   return "jetProbabilityBJetTags";
  else if( op == "JP")    return "jetProbabilityBJetTags";
  else if( op == "TCHPT") return "trackCountingHighPurBJetTags";
  else if( op == "TCHP")  return "trackCountingHighPurBJetTags";
  else if( op == "CSVL")  return "combinedSecondaryVertexBJetTags";
  else if( op == "CSVM")  return "combinedSecondaryVertexBJetTags";
  else if( op == "CSVT")  return "combinedSecondaryVertexBJetTags";
  else if( op == "CSV")   return "combinedSecondaryVertexBJetTags";
  else{
     std::cout  << "Unknown tagger/operating point: "<< op << std::endl;
     return "not defined";  
  }
}

float BtagHardcodedConditions::getDiscriminant(const std::string & op){
  if( op == "JPL")   	  return 0.275;
  else if( op == "JPM")   return 0.545;
  else if( op == "JPT")   return 0.79;
  else if( op == "TCHPT") return 3.41;
  else if( op == "CSVL")  return 0.244;
  else if( op == "CSVM")  return 0.679;
  else if( op == "CSVT")  return 0.898;
  else {
    std::cout << "Unknown operating point: "<< op << std::endl;
    return 10000;  
  }
   
}


double BtagHardcodedConditions::GetBtagEfficiency(double pt, double eta,
					      std::string tagger)
{

  // tag eff, x - discriminant
  // from https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/eff_b_c-ttbar_payload.txt
  float d = getDiscriminant(tagger);
  if ((tagger=="CSV") || (getAlgoTag(tagger)=="CSV")){
    return -4.46932526854*d*d*d*d+7.32781975653*d*d*d-3.78459588569*d*d+0.221027515486*d+0.970299300468;
  } else if ((tagger=="JP") || (getAlgoTag(tagger)=="JP")){
    return -1.3411375438*d*d*d*d+1.86566177901*d*d*d-0.59348240368*d*d-0.893938089125*d+1.22839928411;
  } else if ((tagger=="TCHP") || (getAlgoTag(tagger)=="TCHP")){
    return 9.83842428415e-06*d*d*d*d +  -0.000556835427293*d*d*d +  0.0123891144567*d*d +  -0.141658673059*d +  0.804455651041;
  }

  // unknown tagger, return default
  return -100.0;
}



double BtagHardcodedConditions::GetBtagScaleFactor(double pt, double eta,
					      std::string tagger){
// // THis is 2011 ttbar SF
//   if ((tagger=="CSV") || (getAlgoTag(tagger)=="CSV")){
//     return -0.113472343605*discriminant +  1.04926963159;
//   } else if ((tagger=="JP") || (getAlgoTag(tagger)=="JP")){
//     return -0.193012703295*discriminant+1.05615935852;
//   } else if ((tagger=="TCHP") || (getAlgoTag(tagger)=="TCHP")){
//     return 0.00066342906116*discriminant +  0.93334455507;
//   }

// This is 2011 muon-in-jet
  if (pt>670) pt=670;
  else if (pt<30) pt=30;

  double SFb=0;
  if( tagger=="JPL")   SFb = 0.969851*((1.+(-6.06362e-05*pt))/(1.+(-0.000156638*pt)));
  else if( tagger=="JPM")   SFb = 0.90806*((1.+(0.000236997*pt))/(1.+(5.49455e-05*pt)));
  else if( tagger=="JPT")    SFb = 0.835882*((1.+(0.00167826*pt))/(1.+(0.00120221*pt)));
  else if( tagger=="TCHPT")  SFb = 0.895596*((1.+(9.43219e-05*pt))/(1.+(-4.63927e-05*pt)));
  else if( tagger=="CSVL")  SFb = 1.02658*((1.+(0.0195388*pt))/(1.+(0.0209145*pt)));
  else if( tagger=="CSVM")  SFb = 0.6981*((1.+(0.414063*pt))/(1.+(0.300155*pt)));
  else if( tagger=="CSVT")  SFb = 0.901615*((1.+(0.552628*pt))/(1.+(0.547195*pt)));
  return SFb;
}


double BtagHardcodedConditions::GetBtagSFUncertainty2011(double pt, double eta,
					      std::string tagger)
{
  if (pt<30) return 0.12;
  int bin = findBin(pt);
  float err = -1;

  if( tagger=="JPL")   	   err =  SFb_JPL_error[bin];
  else if( tagger=="JPM")  err =  SFb_JPM_error[bin];
  else if( tagger=="JPT")  err =  SFb_JPT_error[bin];
  else if( tagger=="TCHPT")err =  SFb_TCHPT_error[bin];
  else if( tagger=="CSVL") err =  SFb_CSVL_error[bin];
  else if( tagger=="CSVM") err =  SFb_CSVM_error[bin];
  else if( tagger=="CSVT") err =  SFb_CSVT_error[bin];

  if (pt>670) err*=2.0;
  return err;
}

double BtagHardcodedConditions::GetBtagSFUncertUp(double pt, double eta,
					      std::string tagger)
{
  return GetBtagSFUncertainty2011(pt, eta, tagger)*1.5;
}

double BtagHardcodedConditions::GetBtagSFUncertDown(double pt, double eta,
					      std::string tagger)
{
  return GetBtagSFUncertainty2011(pt, eta, tagger)*1.5;
}


double BtagHardcodedConditions::GetMistagRate(double pt, double eta,
					      std::string tagger){

  // mistag, x-pT
  // from https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/MistagFuncs.C

  if (pt>670) pt=670;
  else if (pt<20) pt=20;
  double _absEta = abs(eta);

  if( tagger == "CSVL" && _absEta>=0.0 && _absEta<0.5)
    return 242534*(((1+(0.0182863*pt))+(4.50105e-05*(pt*pt)))/(1+(108569*pt)));
  if( tagger == "CSVL" && _absEta>=0.5 && _absEta<1.0)
    return 129.938*(((1+(0.0197657*pt))+(4.73472e-05*(pt*pt)))/(1+(55.2415*pt)));
  if( tagger == "CSVL" && _absEta>=1.0 && _absEta<1.5)
    return 592.214*(((1+(0.00671207*pt))+(6.46109e-05*(pt*pt)))/(1+(134.318*pt)));
  if( tagger == "CSVL" && _absEta>=1.5 && _absEta<2.4)
    return 93329*(((1+(0.0219705*pt))+(3.76566e-05*(pt*pt)))/(1+(18245.1*pt)));
  if( tagger == "CSVM" && _absEta>=0.0 && _absEta<0.8)
    return (0.00967751+(2.54564e-05*pt))+(-6.92256e-10*(pt*pt));
  if( tagger == "CSVM" && _absEta>=0.8 && _absEta<1.6)
    return (0.00974141+(5.09503e-05*pt))+(2.0641e-08*(pt*pt));
  if( tagger == "CSVM" && _absEta>=1.6 && _absEta<2.4)
    return (0.013595+(0.000104538*pt))+(-1.36087e-08*(pt*pt));
  if( tagger == "CSVT" && _absEta>=0.0 && _absEta<2.4)
    return 0.00315116*(((1+(-0.00769281*pt))+(2.58066e-05*(pt*pt)))+(-2.02149e-08*(pt*(pt*pt))));
  if( tagger == "JBPL" && _absEta>=0.0 && _absEta<0.5)
    return (0.0277261+(0.000808207*pt))+(-6.44146e-07*(pt*pt));
  if( tagger == "JBPL" && _absEta>=0.5 && _absEta<1.0)
    return (0.0278926+(0.000827697*pt))+(-7.01497e-07*(pt*pt));
  if( tagger == "JBPL" && _absEta>=1.0 && _absEta<1.5)
    return (0.0221411+(0.000900444*pt))+(-6.52873e-07*(pt*pt));
  if( tagger == "JBPL" && _absEta>=1.5 && _absEta<2.4)
    return (0.0227045+(0.000808122*pt))+(-5.67134e-07*(pt*pt));
  if( tagger == "JBPM" && _absEta>=0.0 && _absEta<0.8)
    return (((0.00206106+(0.000105851*pt))+(2.691e-08*(pt*pt)))+(-4.34651e-11*(pt*(pt*pt))))+(-6.73107e-14*(pt*(pt*(pt*pt))));
  if( tagger == "JBPM" && _absEta>=0.8 && _absEta<1.6)
    return (((0.00318438+(4.40327e-05*pt))+(3.46922e-07*(pt*pt)))+(-3.93396e-10*(pt*(pt*pt))))+(3.94283e-14*(pt*(pt*(pt*pt))));
  if( tagger == "JBPM" && _absEta>=1.6 && _absEta<2.4)
    return (((0.00209833+(4.27753e-05*pt))+(1.96076e-07*(pt*pt)))+(6.19275e-11*(pt*(pt*pt))))+(-2.63318e-13*(pt*(pt*(pt*pt))));
  if( tagger == "JBPT" && _absEta>=0.0 && _absEta<2.4)
    return (-3.36681e-05+(1.37292e-05*pt))+(1.78479e-08*(pt*pt));
  if( tagger == "JPL" && _absEta>=0.0 && _absEta<0.5)
    return (0.060001+(0.000332202*pt))+(-2.36709e-07*(pt*pt));
  if( tagger == "JPL" && _absEta>=0.5 && _absEta<1.0)
    return (0.0597675+(0.000370979*pt))+(-2.94673e-07*(pt*pt));
  if( tagger == "JPL" && _absEta>=1.0 && _absEta<1.5)
    return (0.0483728+(0.000528418*pt))+(-3.17825e-07*(pt*pt));
  if( tagger == "JPL" && _absEta>=1.5 && _absEta<2.4)
    return (0.0463159+(0.000546644*pt))+(-3.40486e-07*(pt*pt));
  if( tagger == "JPM" && _absEta>=0.0 && _absEta<0.8)
    return (0.00727084+(4.48901e-05*pt))+(-4.42894e-09*(pt*pt));
  if( tagger == "JPM" && _absEta>=0.8 && _absEta<1.6)
    return (0.00389156+(6.35508e-05*pt))+(1.54183e-08*(pt*pt));
  if( tagger == "JPM" && _absEta>=1.6 && _absEta<2.4)
    return (0.0032816+(4.18867e-05*pt))+(7.44912e-08*(pt*pt));
  if( tagger == "JPT" && _absEta>=0.0 && _absEta<2.4)
    return (0.000379966+(8.30969e-06*pt))+(1.10364e-08*(pt*pt));
  if( tagger == "SSVHEM" && _absEta>=0.0 && _absEta<0.8)
    return (((0.000547883+(0.00023023*pt))+(-7.31792e-07*(pt*pt)))+(1.15659e-09*(pt*(pt*pt))))+(-7.00641e-13*(pt*(pt*(pt*pt))));
  if( tagger == "SSVHEM" && _absEta>=0.8 && _absEta<1.6)
    return (((0.000615562+(0.000240254*pt))+(-7.00237e-07*(pt*pt)))+(1.2566e-09*(pt*(pt*pt))))+(-8.59011e-13*(pt*(pt*(pt*pt))));
  if( tagger == "SSVHEM" && _absEta>=1.6 && _absEta<2.4)
    return (((0.000372388+(0.000309735*pt))+(-4.35952e-07*(pt*pt)))+(3.63763e-10*(pt*(pt*pt))))+(-2.11993e-13*(pt*(pt*(pt*pt))));
  if( tagger == "SSVHPT" && _absEta>=0.0 && _absEta<2.4)
    return (-2.9605e-05+(2.35624e-05*pt))+(-1.77552e-08*(pt*pt));
  if( tagger == "TCHEL" && _absEta>=0.0 && _absEta<0.5)
    return (((-0.0235318+(0.00268868*pt))+(-6.47688e-06*(pt*pt)))+(7.92087e-09*(pt*(pt*pt))))+(-4.06519e-12*(pt*(pt*(pt*pt))));
  if( tagger == "TCHEL" && _absEta>=0.5 && _absEta<1.0)
    return (((-0.0257274+(0.00289337*pt))+(-7.48879e-06*(pt*pt)))+(9.84928e-09*(pt*(pt*pt))))+(-5.40844e-12*(pt*(pt*(pt*pt))));
  if( tagger == "TCHEL" && _absEta>=1.0 && _absEta<1.5)
    return (((-0.0310046+(0.00307803*pt))+(-7.94145e-06*(pt*pt)))+(1.06889e-08*(pt*(pt*pt))))+(-6.08971e-12*(pt*(pt*(pt*pt))));
  if( tagger == "TCHEL" && _absEta>=1.5 && _absEta<2.4)
    return (((-0.0274561+(0.00301096*pt))+(-8.89588e-06*(pt*pt)))+(1.40142e-08*(pt*(pt*pt))))+(-8.95723e-12*(pt*(pt*(pt*pt))));
  if( tagger == "TCHEM" && _absEta>=0.0 && _absEta<0.8)
    return (0.000919586+(0.00026266*pt))+(-1.75723e-07*(pt*pt));
  if( tagger == "TCHEM" && _absEta>=0.8 && _absEta<1.6)
    return (-0.00364137+(0.000350371*pt))+(-1.89967e-07*(pt*pt));
  if( tagger == "TCHEM" && _absEta>=1.6 && _absEta<2.4)
    return (-0.00483904+(0.000367751*pt))+(-1.36152e-07*(pt*pt));
  if( tagger == "TCHPM" && _absEta>=0.0 && _absEta<0.8)
    return (((-0.00464673+(0.000247485*pt))+(9.13236e-07*(pt*pt)))+(-2.49994e-09*(pt*(pt*pt))))+(1.65678e-12*(pt*(pt*(pt*pt))));
  if( tagger == "TCHPM" && _absEta>=0.8 && _absEta<1.6)
    return (((-0.0060878+(0.000297422*pt))+(1.13369e-06*(pt*pt)))+(-2.84945e-09*(pt*(pt*pt))))+(1.64721e-12*(pt*(pt*(pt*pt))));
  if( tagger == "TCHPM" && _absEta>=1.6 && _absEta<2.4)
    return (((-0.00836219+(0.000391889*pt))+(2.78156e-07*(pt*pt)))+(-6.14017e-10*(pt*(pt*pt))))+(-1.30592e-13*(pt*(pt*(pt*pt))));
  if( tagger == "TCHPT" && _absEta>=0.0 && _absEta<2.4)
    return (-0.00101+(4.70405e-05*pt))+(8.3338e-09*(pt*pt));

  // unknown tagger, return default
  return -100.0;
}


double BtagHardcodedConditions::GetMistagScaleFactor(double pt, double eta,
						     std::string tagger){
  return GetMistagSF2011(pt, eta, tagger, "mean") * GetMistagCorr2012(pt, eta, tagger);
}

double BtagHardcodedConditions::GetMistagSFUncertDown(double pt, double eta,
						     std::string tagger){
  return GetMistagCorr2012(pt, eta, tagger) * (pt>670?2.0:1.0) *
	(GetMistagSF2011(pt, eta, tagger, "mean")- GetMistagSF2011(pt, eta, tagger, "min"));
}

double BtagHardcodedConditions::GetMistagSFUncertUp(double pt, double eta,
						     std::string tagger){
  return GetMistagCorr2012(pt, eta, tagger) * (pt>670?2.0:1.0) *
	(GetMistagSF2011(pt, eta, tagger, "max")- GetMistagSF2011(pt, eta, tagger, "mean"));
}

double BtagHardcodedConditions::GetMistagCorr2012(double pt, double eta,
						     std::string tagger){

  // 2012 correction factors from https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlighCorrFunc_2012_06_06.txt

  if (pt>670) pt=670;
  else if (pt<20) pt=20;
  double corrFact=0;

  if( tagger=="JPL")   corrFact= 1.01744 + 0.000813491*pt - 6.01592e-07*pt*pt;
  else if( tagger=="JPM")   corrFact= 0.964487 + 0.00134038*pt - 1.43995e-06*pt*pt;
  else if( tagger=="JPT")   corrFact= 1.06452 + 0.00080354*pt - 1.34109e-06*pt*pt;
  else if( tagger=="TCHPT") corrFact= 1.08376 - 0.000666189*pt + 1.01272e-06*pt*pt;
  else if( tagger=="CSVL")  corrFact= 0.979396 + 0.000205898*pt + 2.49868e-07*pt*pt;
  else if( tagger=="CSVM")  corrFact= 1.10422 - 0.000523856*pt + 1.14251e-06*pt*pt;
  else if( tagger=="CSVT")  corrFact= 1.19275 - 0.00191042*pt + 2.92205e-06*pt*pt;
  return corrFact;
}

double BtagHardcodedConditions::GetMistagSF2011(double pt, double eta,
					     std::string tagger, std::string meanminmax)
{

  double _absEta = abs(eta);
  double sf = -1;
  if ((pt<670)||(tagger[tagger.length()-1]=='T')) {
    if (pt<20) pt=20;
    else if (pt>670) pt=670;
    if( tagger=="CSVL"&& _absEta>=0.0 && _absEta<0.5) {
      if( meanminmax == "mean" ) sf = ((1.07536+(0.000175506*pt))+(-8.63317e-07*(pt*pt)))+(3.27516e-10*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.994425+(-8.66392e-05*pt))+(-3.03813e-08*(pt*pt)))+(-3.52151e-10*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.15628+(0.000437668*pt))+(-1.69625e-06*(pt*pt)))+(1.00718e-09*(pt*(pt*pt)));
    }
    else if( tagger=="CSVL" && _absEta>=0.5 && _absEta<1.0) {
      if( meanminmax == "mean" ) sf = ((1.07846+(0.00032458*pt))+(-1.30258e-06*(pt*pt)))+(8.50608e-10*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.998088+(6.94916e-05*pt))+(-4.82731e-07*(pt*pt)))+(1.63506e-10*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.15882+(0.000579711*pt))+(-2.12243e-06*(pt*pt)))+(1.53771e-09*(pt*(pt*pt)));
    }
    else if( tagger=="CSVL" && _absEta>=1.0 && _absEta<1.5) {
      if( meanminmax == "mean" ) sf = ((1.08294+(0.000474818*pt))+(-1.43857e-06*(pt*pt)))+(1.13308e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((1.00294+(0.000289844*pt))+(-7.9845e-07*(pt*pt)))+(5.38525e-10*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.16292+(0.000659848*pt))+(-2.07868e-06*(pt*pt)))+(1.72763e-09*(pt*(pt*pt)));
    }
    else if( tagger=="CSVL" && _absEta>=1.5 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = ((1.0617+(0.000173654*pt))+(-5.29009e-07*(pt*pt)))+(5.55931e-10*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.979816+(0.000138797*pt))+(-3.14503e-07*(pt*pt)))+(2.38124e-10*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.14357+(0.00020854*pt))+(-7.43519e-07*(pt*pt)))+(8.73742e-10*(pt*(pt*pt)));
    }
    else if( tagger=="CSVM" && _absEta>=0.0 && _absEta<0.8) {
      if( meanminmax == "mean" ) sf = ((1.06182+(0.000617034*pt))+(-1.5732e-06*(pt*pt)))+(3.02909e-10*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.972455+(7.51396e-06*pt))+(4.91857e-07*(pt*pt)))+(-1.47661e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.15116+(0.00122657*pt))+(-3.63826e-06*(pt*pt)))+(2.08242e-09*(pt*(pt*pt)));
    }
    else if( tagger=="CSVM" && _absEta>=0.8 && _absEta<1.6) {
      if( meanminmax == "mean" ) sf = ((1.111+(-9.64191e-06*pt))+(1.80811e-07*(pt*pt)))+(-5.44868e-10*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((1.02055+(-0.000378856*pt))+(1.49029e-06*(pt*pt)))+(-1.74966e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.20146+(0.000359543*pt))+(-1.12866e-06*(pt*pt)))+(6.59918e-10*(pt*(pt*pt)));
    }
    else if( tagger=="CSVM" && _absEta>=1.6 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = ((1.08498+(-0.000701422*pt))+(3.43612e-06*(pt*pt)))+(-4.11794e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.983476+(-0.000607242*pt))+(3.17997e-06*(pt*pt)))+(-4.01242e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.18654+(-0.000795808*pt))+(3.69226e-06*(pt*pt)))+(-4.22347e-09*(pt*(pt*pt)));
    }
    else if( tagger=="CSVT" && _absEta>=0.0 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = ((0.948463+(0.00288102*pt))+(-7.98091e-06*(pt*pt)))+(5.50157e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.899715+(0.00102278*pt))+(-2.46335e-06*(pt*pt)))+(9.71143e-10*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((0.997077+(0.00473953*pt))+(-1.34985e-05*(pt*pt)))+(1.0032e-08*(pt*(pt*pt)));
    }
    else if( tagger=="JBPL" && _absEta>=0.0 && _absEta<0.5) {
      if( meanminmax == "mean" ) sf = ((0.996303+(-0.00049586*pt))+(1.48662e-06*(pt*pt)))+(-1.60955e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.909313+(-0.000483037*pt))+(1.48507e-06*(pt*pt)))+(-1.60327e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.08332+(-0.000508763*pt))+(1.48816e-06*(pt*pt)))+(-1.61583e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JBPL" && _absEta>=0.5 && _absEta<1.0) {
      if( meanminmax == "mean" ) sf = ((1.01607+(-0.000958122*pt))+(3.12318e-06*(pt*pt)))+(-3.13777e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.925793+(-0.000877501*pt))+(2.88538e-06*(pt*pt)))+(-2.9089e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.10639+(-0.0010389*pt))+(3.36098e-06*(pt*pt)))+(-3.36665e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JBPL" && _absEta>=1.0 && _absEta<1.5) {
      if( meanminmax == "mean" ) sf = ((1.04234+(-0.00109152*pt))+(3.71686e-06*(pt*pt)))+(-3.57219e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.947786+(-0.000985917*pt))+(3.39659e-06*(pt*pt)))+(-3.28635e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.13696+(-0.00119731*pt))+(4.03713e-06*(pt*pt)))+(-3.85803e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JBPL" && _absEta>=1.5 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = ((0.960685+(-0.000514241*pt))+(2.69297e-06*(pt*pt)))+(-3.12123e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.875356+(-0.000455763*pt))+(2.42337e-06*(pt*pt)))+(-2.83637e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.04606+(-0.000572874*pt))+(2.96257e-06*(pt*pt)))+(-3.40609e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JBPM" && _absEta>=0.0 && _absEta<0.8) {
      if( meanminmax == "mean" ) sf = ((0.932447+(0.000285676*pt))+(-1.03771e-06*(pt*pt)))+(4.52275e-10*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.822274+(0.000138316*pt))+(-4.14616e-07*(pt*pt)))+(-9.7638e-11*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.0426+(0.000433059*pt))+(-1.6608e-06*(pt*pt)))+(1.00219e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JBPM" && _absEta>=0.8 && _absEta<1.6) {
      if( meanminmax == "mean" ) sf = ((0.924959+(0.000170347*pt))+(-1.56056e-07*(pt*pt)))+(-2.06751e-10*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.822394+(-2.61379e-05*pt))+(6.08356e-07*(pt*pt)))+(-9.28476e-10*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.02752+(0.000366822*pt))+(-9.20467e-07*(pt*pt)))+(5.14974e-10*(pt*(pt*pt)));
    }
    else if( tagger=="JBPM" && _absEta>=1.6 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = ((0.846053+(0.000224848*pt))+(2.87503e-07*(pt*pt)))+(-5.93182e-10*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.714511+(0.000568422*pt))+(-7.56289e-07*(pt*pt)))+(2.61634e-10*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((0.977599+(-0.000118755*pt))+(1.3313e-06*(pt*pt)))+(-1.448e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JBPT" && _absEta>=0.0 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = ((0.771257+(0.00238891*pt))+(-6.2112e-06*(pt*pt)))+(4.33595e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.666101+(0.00163462*pt))+(-3.92728e-06*(pt*pt)))+(2.48081e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((0.87631+(0.00314343*pt))+(-8.49513e-06*(pt*pt)))+(6.19109e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JPL" && _absEta>=0.0 && _absEta<0.5) {
      if( meanminmax == "mean" ) sf = ((1.02571+(-0.000391686*pt))+(1.01948e-06*(pt*pt)))+(-1.16475e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.931859+(-0.00045457*pt))+(1.25431e-06*(pt*pt)))+(-1.36433e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.11958+(-0.00032886*pt))+(7.84649e-07*(pt*pt)))+(-9.65161e-10*(pt*(pt*pt)));
    }
    else if( tagger=="JPL" && _absEta>=0.5 && _absEta<1.0) {
      if( meanminmax == "mean" ) sf = ((1.03375+(-0.00068776*pt))+(2.13443e-06*(pt*pt)))+(-2.24163e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.936905+(-0.000681017*pt))+(2.13885e-06*(pt*pt)))+(-2.22607e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.13063+(-0.000694616*pt))+(2.13001e-06*(pt*pt)))+(-2.25719e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JPL" && _absEta>=1.0 && _absEta<1.5) {
      if( meanminmax == "mean" ) sf = ((1.03597+(-0.000778058*pt))+(3.02129e-06*(pt*pt)))+(-3.0478e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.938438+(-0.00074623*pt))+(2.89732e-06*(pt*pt)))+(-2.92483e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.13355+(-0.000810039*pt))+(3.14525e-06*(pt*pt)))+(-3.17077e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JPL" && _absEta>=1.5 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = ((0.95897+(-0.000111286*pt))+(1.6091e-06*(pt*pt)))+(-2.18387e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.867768+(-9.92078e-05*pt))+(1.46903e-06*(pt*pt)))+(-2.02118e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.0502+(-0.000123474*pt))+(1.74917e-06*(pt*pt)))+(-2.34655e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JPM" && _absEta>=0.0 && _absEta<0.8) {
      if( meanminmax == "mean" ) sf = ((0.970028+(0.00118179*pt))+(-4.23119e-06*(pt*pt)))+(3.61065e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.840326+(0.000626372*pt))+(-2.08293e-06*(pt*pt)))+(1.57604e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.09966+(0.00173739*pt))+(-6.37946e-06*(pt*pt)))+(5.64527e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JPM" && _absEta>=0.8 && _absEta<1.6) {
      if( meanminmax == "mean" ) sf = ((0.918387+(0.000898595*pt))+(-2.00643e-06*(pt*pt)))+(1.26486e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.790843+(0.000548016*pt))+(-6.70941e-07*(pt*pt)))+(1.90355e-11*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.0459+(0.00124924*pt))+(-3.34192e-06*(pt*pt)))+(2.51068e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JPM" && _absEta>=1.6 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = ((0.790103+(0.00117865*pt))+(-2.07334e-06*(pt*pt)))+(1.42608e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.667144+(0.00105593*pt))+(-1.43608e-06*(pt*pt)))+(5.24039e-10*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((0.913027+(0.00130143*pt))+(-2.71061e-06*(pt*pt)))+(2.32812e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JPT" && _absEta>=0.0 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = ((0.831392+(0.00269525*pt))+(-7.33391e-06*(pt*pt)))+(5.73942e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.671888+(0.0020106*pt))+(-5.03177e-06*(pt*pt)))+(3.74225e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((0.990774+(0.00338018*pt))+(-9.63606e-06*(pt*pt)))+(7.73659e-09*(pt*(pt*pt)));
    }
    else if( tagger=="SSVHEM" && _absEta>=0.0 && _absEta<0.8) {
      if( meanminmax == "mean" ) sf = ((0.86318+(0.000801639*pt))+(-1.64119e-06*(pt*pt)))+(2.59121e-10*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.790364+(0.000463086*pt))+(-4.35934e-07*(pt*pt)))+(-9.08296e-10*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((0.935969+(0.0011402*pt))+(-2.84645e-06*(pt*pt)))+(1.42654e-09*(pt*(pt*pt)));
    }
    else if( tagger=="SSVHEM" && _absEta>=0.8 && _absEta<1.6) {
      if( meanminmax == "mean" ) sf = ((0.958973+(-0.000269555*pt))+(1.381e-06*(pt*pt)))+(-1.87744e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.865771+(-0.000279908*pt))+(1.34144e-06*(pt*pt)))+(-1.75588e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.0522+(-0.000259296*pt))+(1.42056e-06*(pt*pt)))+(-1.999e-09*(pt*(pt*pt)));
    }
    else if( tagger=="SSVHEM" && _absEta>=1.6 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = ((0.923033+(-0.000898227*pt))+(4.74565e-06*(pt*pt)))+(-6.11053e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.828021+(-0.000731926*pt))+(4.19613e-06*(pt*pt)))+(-5.81379e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.01812+(-0.00106483*pt))+(5.29518e-06*(pt*pt)))+(-6.40728e-09*(pt*(pt*pt)));
    }
    else if( tagger=="SSVHPT" && _absEta>=0.0 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = ((0.97409+(0.000646241*pt))+(-2.86294e-06*(pt*pt)))+(2.79484e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.807222+(0.00103676*pt))+(-3.6243e-06*(pt*pt)))+(3.17368e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.14091+(0.00025586*pt))+(-2.10157e-06*(pt*pt)))+(2.41599e-09*(pt*(pt*pt)));
    }
    else if( tagger=="TCHEL" && _absEta>=0.0 && _absEta<0.5) {
      if( meanminmax == "mean" ) sf = (1.13615*((1+(-0.00119852*pt))+(1.17888e-05*(pt*pt))))+(-9.8581e-08*(pt*(pt*(pt/(1+(0.00689317*pt))))));
      else if( meanminmax == "min" ) sf = (1.0369*((1+(-0.000945578*pt))+(7.73273e-06*(pt*pt))))+(-4.47791e-08*(pt*(pt*(pt/(1+(0.00499343*pt))))));
      else if( meanminmax == "max" ) sf = (1.22179*((1+(-0.000946228*pt))+(7.37821e-06*(pt*pt))))+(-4.8451e-08*(pt*(pt*(pt/(1+(0.0047976*pt))))));
    }
    else if( tagger=="TCHEL" && _absEta>=0.5 && _absEta<1.0) {
      if( meanminmax == "mean" ) sf = (1.13277*((1+(-0.00084146*pt))+(3.80313e-06*(pt*pt))))+(-8.75061e-09*(pt*(pt*(pt/(1+(0.00118695*pt))))));
      else if( meanminmax == "min" ) sf = (0.983748*((1+(7.13613e-05*pt))+(-1.08648e-05*(pt*pt))))+(2.96162e-06*(pt*(pt*(pt/(1+(0.282104*pt))))));
      else if( meanminmax == "max" ) sf = (1.22714*((1+(-0.00085562*pt))+(3.74425e-06*(pt*pt))))+(-8.91028e-09*(pt*(pt*(pt/(1+(0.00109346*pt))))));
    }
    else if( tagger=="TCHEL" && _absEta>=1.0 && _absEta<1.5) {
      if( meanminmax == "mean" ) sf = (1.17163*((1+(-0.000828475*pt))+(3.0769e-06*(pt*pt))))+(-4.692e-09*(pt*(pt*(pt/(1+(0.000337759*pt))))));
      else if( meanminmax == "min" ) sf = (1.0698*((1+(-0.000731877*pt))+(2.56922e-06*(pt*pt))))+(-3.0318e-09*(pt*(pt*(pt/(1+(5.04118e-05*pt))))));
      else if( meanminmax == "max" ) sf = (1.27351*((1+(-0.000911891*pt))+(3.5465e-06*(pt*pt))))+(-6.69625e-09*(pt*(pt*(pt/(1+(0.000590847*pt))))));
    }
    else if( tagger=="TCHEL" && _absEta>=1.5 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = (1.14554*((1+(-0.000128043*pt))+(4.10899e-07*(pt*pt))))+(-2.07565e-10*(pt*(pt*(pt/(1+(-0.00118618*pt))))));
      else if( meanminmax == "min" ) sf = (1.04766*((1+(-6.87499e-05*pt))+(2.2454e-07*(pt*pt))))+(-1.18395e-10*(pt*(pt*(pt/(1+(-0.00128734*pt))))));
      else if( meanminmax == "max" ) sf = (1.24367*((1+(-0.000182494*pt))+(5.92637e-07*(pt*pt))))+(-3.3745e-10*(pt*(pt*(pt/(1+(-0.00107694*pt))))));
    }
    else if( tagger=="TCHEM" && _absEta>=0.0 && _absEta<0.8) {
      if( meanminmax == "mean" ) sf = (1.2875*((1+(-0.000356371*pt))+(1.08081e-07*(pt*pt))))+(-6.89998e-11*(pt*(pt*(pt/(1+(-0.0012139*pt))))));
      else if( meanminmax == "min" ) sf = (1.11418*((1+(-0.000442274*pt))+(1.53463e-06*(pt*pt))))+(-4.93683e-09*(pt*(pt*(pt/(1+(0.00152436*pt))))));
      else if( meanminmax == "max" ) sf = (1.47515*((1+(-0.000484868*pt))+(2.36817e-07*(pt*pt))))+(-2.05073e-11*(pt*(pt*(pt/(1+(-0.00142819*pt))))));
    }
    else if( tagger=="TCHEM" && _absEta>=0.8 && _absEta<1.6) {
      if( meanminmax == "mean" ) sf = (1.24986*((1+(-0.00039734*pt))+(5.37486e-07*(pt*pt))))+(-1.74023e-10*(pt*(pt*(pt/(1+(-0.00112954*pt))))));
      else if( meanminmax == "min" ) sf = (1.08828*((1+(-0.000208737*pt))+(1.50487e-07*(pt*pt))))+(-2.54249e-11*(pt*(pt*(pt/(1+(-0.00141477*pt))))));
      else if( meanminmax == "max" ) sf = (1.41211*((1+(-0.000559603*pt))+(9.50754e-07*(pt*pt))))+(-5.81148e-10*(pt*(pt*(pt/(1+(-0.000787359*pt))))));
    }
    else if( tagger=="TCHEM" && _absEta>=1.6 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = (1.10763*((1+(-0.000105805*pt))+(7.11718e-07*(pt*pt))))+(-5.3001e-10*(pt*(pt*(pt/(1+(-0.000821215*pt))))));
      else if( meanminmax == "min" ) sf = (0.958079*((1+(0.000327804*pt))+(-4.09511e-07*(pt*pt))))+(-1.95933e-11*(pt*(pt*(pt/(1+(-0.00143323*pt))))));
      else if( meanminmax == "max" ) sf = (1.26236*((1+(-0.000524055*pt))+(2.08863e-06*(pt*pt))))+(-2.29473e-09*(pt*(pt*(pt/(1+(-0.000276268*pt))))));
    }
    else if( tagger=="TCHPM" && _absEta>=0.0 && _absEta<0.8) {
      if( meanminmax == "mean" ) sf = ((1.27011+(-0.000869141*pt))+(2.49796e-06*(pt*pt)))+(-2.62962e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((1.12949+(-0.000678492*pt))+(2.02219e-06*(pt*pt)))+(-2.21675e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.41077+(-0.00105992*pt))+(2.97373e-06*(pt*pt)))+(-3.0425e-09*(pt*(pt*pt)));
    }
    else if( tagger=="TCHPM" && _absEta>=0.8 && _absEta<1.6) {
      if( meanminmax == "mean" ) sf = ((1.36167+(-0.00153237*pt))+(4.54567e-06*(pt*pt)))+(-4.38874e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((1.21289+(-0.00126411*pt))+(3.81676e-06*(pt*pt)))+(-3.75847e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.51053+(-0.00180085*pt))+(5.27457e-06*(pt*pt)))+(-5.01901e-09*(pt*(pt*pt)));
    }
    else if( tagger=="TCHPM" && _absEta>=1.6 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = ((1.22696+(0.000249231*pt))+(9.55279e-08*(pt*pt)))+(-1.04034e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((1.07572+(0.00055366*pt))+(-9.55796e-07*(pt*pt)))+(-3.73943e-11*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.3782+(-5.52498e-05*pt))+(1.14685e-06*(pt*pt)))+(-2.04329e-09*(pt*(pt*pt)));
    }
    else if( tagger=="TCHPT" && _absEta>=0.0 && _absEta<2.4) {
      if( meanminmax == "mean" ) sf = ((1.20711+(0.000681067*pt))+(-1.57062e-06*(pt*pt)))+(2.83138e-10*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((1.03418+(0.000428273*pt))+(-5.43024e-07*(pt*pt)))+(-6.18061e-10*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.38002+(0.000933875*pt))+(-2.59821e-06*(pt*pt)))+(1.18434e-09*(pt*(pt*pt)));
    }
  } else {
    pt=670;
    if( tagger=="CSVM" && _absEta>=0.0 && _absEta<2.4)
    {
      if( meanminmax == "mean" ) sf = ((1.04318+(0.000848162*pt))+(-2.5795e-06*(pt*pt)))+(1.64156e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.962627+(0.000448344*pt))+(-1.25579e-06*(pt*pt)))+(4.82283e-10*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.12368+(0.00124806*pt))+(-3.9032e-06*(pt*pt)))+(2.80083e-09*(pt*(pt*pt)));
    }
    else if( tagger=="CSVL" && _absEta>=0.0 && _absEta<2.4)
    {
      if( meanminmax == "mean" ) sf = ((1.0344+(0.000962994*pt))+(-3.65392e-06*(pt*pt)))+(3.23525e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.956023+(0.000825106*pt))+(-3.18828e-06*(pt*pt)))+(2.81787e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.11272+(0.00110104*pt))+(-4.11956e-06*(pt*pt)))+(3.65263e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JBPL" && _absEta>=0.0 && _absEta<2.4)
    {
      if( meanminmax == "mean" ) sf = ((0.974968+(-0.000192541*pt))+(7.08162e-07*(pt*pt)))+(-9.7623e-10*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.891217+(-0.000201377*pt))+(7.41513e-07*(pt*pt)))+(-9.86349e-10*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.05873+(-0.000183755*pt))+(6.74811e-07*(pt*pt)))+(-9.6611e-10*(pt*(pt*pt)));
    }
    else if( tagger=="JBPM" && _absEta>=0.0 && _absEta<2.4)
    {
      if( meanminmax == "mean" ) sf = ((0.859232+(0.00117533*pt))+(-3.51857e-06*(pt*pt)))+(2.63162e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.769143+(0.000865281*pt))+(-2.47018e-06*(pt*pt)))+(1.72476e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((0.949262+(0.00148551*pt))+(-4.56695e-06*(pt*pt)))+(3.53847e-09*(pt*(pt*pt)));
    }
    else if( tagger=="JPL" && _absEta>=0.0 && _absEta<2.4)
    {
      if( meanminmax == "mean" ) sf = ((0.978061+(0.000211142*pt))+(-6.67003e-07*(pt*pt)))+(2.94232e-10*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.889182+(0.000124259*pt))+(-3.83838e-07*(pt*pt)))+(5.99164e-11*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.06693+(0.00029804*pt))+(-9.50169e-07*(pt*pt)))+(5.28547e-10*(pt*(pt*pt)));
    }
    else if( tagger=="JPM" && _absEta>=0.0 && _absEta<2.4)
    {
      if( meanminmax == "mean" ) sf = ((0.871294+(0.00215201*pt))+(-6.77675e-06*(pt*pt)))+(5.79197e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.7654+(0.00149792*pt))+(-4.47192e-06*(pt*pt)))+(3.67664e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((0.977076+(0.00280638*pt))+(-9.08158e-06*(pt*pt)))+(7.9073e-09*(pt*(pt*pt)));
    }
    else if( tagger=="SSVHEM" && _absEta>=0.0 && _absEta<2.4)
    {
      if( meanminmax == "mean" ) sf = ((0.890254+(0.000553319*pt))+(-1.29993e-06*(pt*pt)))+(4.19294e-10*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((0.817099+(0.000421567*pt))+(-9.46432e-07*(pt*pt)))+(1.62339e-10*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((0.963387+(0.000685092*pt))+(-1.65343e-06*(pt*pt)))+(6.76249e-10*(pt*(pt*pt)));
    }
    else if( tagger=="TCHEL" && _absEta>=0.0 && _absEta<2.4)
    {
      if( meanminmax == "mean" ) sf = (1.10649*((1+(-9.00297e-05*pt))+(2.32185e-07*(pt*pt))))+(-4.04925e-10*(pt*(pt*(pt/(1+(-0.00051036*pt))))));
      else if( meanminmax == "min" ) sf = (1.01541*((1+(-6.04627e-05*pt))+(1.38195e-07*(pt*pt))))+(-2.83043e-10*(pt*(pt*(pt/(1+(-0.000633609*pt))))));
      else if( meanminmax == "max" ) sf = (1.19751*((1+(-0.000114197*pt))+(3.08558e-07*(pt*pt))))+(-5.27598e-10*(pt*(pt*(pt/(1+(-0.000422372*pt))))));
    }
    else if( tagger=="TCHEM" && _absEta>=0.0 && _absEta<2.4)
    {
      if( meanminmax == "mean" ) sf = (1.06268*((1+(0.00390509*pt))+(-5.85405e-05*(pt*pt))))+(7.87135e-07*(pt*(pt*(pt/(1+(0.01259*pt))))));
      else if( meanminmax == "min" ) sf = (0.967092*((1+(0.00201431*pt))+(-1.49359e-05*(pt*pt))))+(6.94324e-08*(pt*(pt*(pt/(1+(0.00459787*pt))))));
      else if( meanminmax == "max" ) sf = (1.22691*((1+(0.00211682*pt))+(-2.07959e-05*(pt*pt))))+(1.72938e-07*(pt*(pt*(pt/(1+(0.00658853*pt))))));
    }
    else if( tagger=="TCHPM" && _absEta>=0.0 && _absEta<2.4)
    {
      if( meanminmax == "mean" ) sf = ((1.27417+(-0.000449095*pt))+(1.0719e-06*(pt*pt)))+(-1.35208e-09*(pt*(pt*pt)));
      else if( meanminmax == "min" ) sf = ((1.14205+(-0.000350151*pt))+(8.43333e-07*(pt*pt)))+(-1.14104e-09*(pt*(pt*pt)));
      else if( meanminmax == "max" ) sf = ((1.4063+(-0.000548107*pt))+(1.30047e-06*(pt*pt)))+(-1.56311e-09*(pt*(pt*pt)));
    }
  }
  
  return sf;
}


	 
