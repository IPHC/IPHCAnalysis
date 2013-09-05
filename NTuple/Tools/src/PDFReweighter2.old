#include "Tools/interface/PDFReweighter2.h"

//#ifdef PDFREWEIGHT_ENABLE

using namespace IPHCTree;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// -----------------------------------------------------------------------------
bool PDFReweighter2::Initialize()
{
  std::cout << "===============================================================" << std::endl;
  std::cout << "===============================================================" << std::endl;
  std::cout << "                 INITIALIZATION PDF-REWEIGHTER " << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::string tab="   ";
  std::cout << tab << "       _==/          i     i          \\==_         " << std::endl; 
  std::cout << tab << "     /XX/            |\\___/|            \\XX\\       " << std::endl; 
  std::cout << tab << "   /XXXX\\            |XXXXX|            /XXXX\\     " << std::endl; 
  std::cout << tab << "  |XXXXXX\\_         _XXXXXXX_         _/XXXXXX|    " << std::endl; 
  std::cout << tab << " XXXXXXXXXXXxxxxxxxXXXXXXXXXXXxxxxxxxXXXXXXXXXXX   " << std::endl; 
  std::cout << tab << "|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX|  " << std::endl; 
  std::cout << tab << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  " << std::endl; 
  std::cout << tab << "|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX|  " << std::endl; 
  std::cout << tab << " XXXXXX/^^^^\"\\XXXXXXXXXXXXXXXXXXXXX/^^^^^\\XXXXXX   " << std::endl; 
  std::cout << tab << "  |XXX|       \\XXX/^^\\XXXXX/^^\\XXX/       |XXX|    " << std::endl; 
  std::cout << tab << "    \\XX\\       \\X/    \\XXX/    \\X/       /XX/      " << std::endl; 
  std::cout << tab << "       \"\\       \"      \\X/      \"       /\"        " << std::endl; 
  std::cout << std::endl;
  std::cout << std::endl;

  // CTEQ uncertainties PDF
  NNPDF_pdf_  = std::make_pair("NNPDF21_100.LHgrid", 0);
  LHAPDF::initPDFSet(1,NNPDF_pdf_.first);
  NNPDF_pdf_.second = LHAPDF::numberPDF(1);
  std::cout << " => PDF '" << NNPDF_pdf_.first << "' is loaded with " 
            << NNPDF_pdf_.second << " subsets" << std::endl << std::endl;

  // CTEQ uncertainties PDF
  NNPDF15_pdf_  = std::make_pair("NNPDF21_mc_15_100.LHgrid", 0);
  LHAPDF::initPDFSet(2,NNPDF15_pdf_.first);
  NNPDF15_pdf_.second = LHAPDF::numberPDF(2);
  std::cout << " => PDF '" << NNPDF15_pdf_.first << "' is loaded with " 
            << NNPDF15_pdf_.second << " subsets" << std::endl << std::endl;

  // CTEQ uncertainties PDF
  NNPDF16_pdf_  = std::make_pair("NNPDF21_mc_16_100.LHgrid", 0);
  LHAPDF::initPDFSet(3,NNPDF16_pdf_.first);
  NNPDF16_pdf_.second = LHAPDF::numberPDF(3);
  std::cout << " => PDF '" << NNPDF16_pdf_.first << "' is loaded with " 
            << NNPDF16_pdf_.second << " subsets" << std::endl << std::endl;

  // CTEQ uncertainties PDF
  NNPDF17_pdf_  = std::make_pair("NNPDF21_mc_17_100.LHgrid", 0);
  LHAPDF::initPDFSet(4,NNPDF17_pdf_.first);
  NNPDF17_pdf_.second = LHAPDF::numberPDF(4);
  std::cout << " => PDF '" << NNPDF17_pdf_.first << "' is loaded with " 
            << NNPDF17_pdf_.second << " subsets" << std::endl << std::endl;

  // Initial PDF
  REF_pdf_ = std::make_pair("cteq6l.LHpdf", 0);
  //REF_pdf_ = std::make_pair("cteq6ll.LHpdf", 0);
  LHAPDF::initPDFSet(5,REF_pdf_.first);
  REF_pdf_.second = LHAPDF::numberPDF(5);
  if (REF_pdf_.second == 1) REF_pdf_.second=0;
  std::cout << " => PDF '" << REF_pdf_.first << "' is loaded with Best Fit + " 
            << REF_pdf_.second << " subsets" << std::endl << std::endl;
  if (REF_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  std::cout << "===============================================================" << std::endl;
  std::cout << "===============================================================" << std::endl;
  std::cout << std::flush;

  return true;
}


// -----------------------------------------------------------------------------
// Finalize
// -----------------------------------------------------------------------------
bool PDFReweighter2::Finalize()
{
  std::cout << "ERIC FINALIZE" << std::endl;

  return true;
}


// -----------------------------------------------------------------------------
// Calculate
// -----------------------------------------------------------------------------
UncertaintyType PDFReweighter2::Calculate(const NTMonteCarlo& mc,
                                          unsigned int mode) const
{
  Double_t pdf1=0.;
  Double_t pdf2=0.;
  UncertaintyType result;

  // Calculate initial pdf values
  CalculateRefWeight(mc,pdf1,pdf2);
  if (pdf1<1e-20 || pdf2<1e-20)
  {
    std::cout << "WARNING: initial pdf are equal to zero";
    std::cout << "         return empty weights" << std::endl;
    return UncertaintyType();
  }
  double ipdf=pdf1*pdf2;

  // Calculate NNPDF weights
  UncertaintyType NNPDFvalue   = CalculateNNPDF(1,mc,pdf1,pdf2);
  UncertaintyType NNPDF15value = CalculateNNPDF(2,mc,pdf1,pdf2);
  UncertaintyType NNPDF16value = CalculateNNPDF(3,mc,pdf1,pdf2);
  UncertaintyType NNPDF17value = CalculateNNPDF(4,mc,pdf1,pdf2);

  // Calculating envelope upper limit
  Double_t NNPDFup   = NNPDFvalue.Mean   + NNPDFvalue.Max;
  Double_t NNPDF15up = NNPDF15value.Mean + NNPDFvalue.Max;
  Double_t NNPDF16up = NNPDF16value.Mean + NNPDFvalue.Max;
  Double_t NNPDF17up = NNPDF17value.Mean + NNPDFvalue.Max;

  // Calculating envelope lower limit
  Double_t NNPDFdown   = NNPDFvalue.Mean   - NNPDFvalue.Min;
  Double_t NNPDF15down = NNPDF15value.Mean - NNPDFvalue.Min;
  Double_t NNPDF16down = NNPDF16value.Mean - NNPDFvalue.Min;
  Double_t NNPDF17down = NNPDF17value.Mean - NNPDFvalue.Min;

  // Calculating the result
  if (mode==0)
  {
    result.Max = std::max(NNPDFup,std::max(NNPDF15up,
                                  std::max(NNPDF16up,NNPDF17up)));
    result.Min = std::min(NNPDFdown,std::min(NNPDF15down,
                                    std::min(NNPDF16down,NNPDF17down)));
  }
  else if (mode==1)
  {
    result.Max = NNPDFup;
    result.Min = NNPDFdown;
  }
  else if (mode==2)
  {
    result.Max = NNPDF15up;
    result.Min = NNPDF15down;
  }
  else if (mode==3)
  {
    result.Max = NNPDF16up;
    result.Min = NNPDF16down;
  }
  else if (mode==4)
  {
    result.Max = NNPDF17up;
    result.Min = NNPDF17down;
  }

  /*  // Good boundaries ?
  if (result.Min<0.01 || result.Max>10)
  {
    std::cout << "WARNING: weight boundaries are excessive ! "
              << "Mean:" << result.Mean 
              << " Max:" << result.Max 
              << " Min:" << result.Min << std::endl;
    result.Min=0; result.Max=0;
    }*/

  // Calculating envelope mid-point
  result.Mean = (result.Max+result.Min)/2.;

  return result;
}


// -----------------------------------------------------------------------------
// CalculateNNPDF
// -----------------------------------------------------------------------------
UncertaintyType PDFReweighter2::CalculateNNPDF(unsigned int number, const NTMonteCarlo& mc, 
                               Double_t& pdf1, Double_t& pdf2) const
{
  // Computing best fit
  Double_t bestfit = CalculateWeight(number,0,mc,pdf1,pdf2);

  // Computing global weight related to uncertainties
  std::vector<Double_t> weights;

  if      (number==1) weights.resize(NNPDF_pdf_.second,0.);
  else if (number==2) weights.resize(NNPDF15_pdf_.second,0.);
  else if (number==3) weights.resize(NNPDF16_pdf_.second,0.);
  else if (number==4) weights.resize(NNPDF17_pdf_.second,0.);

  for (unsigned int subset=0;subset<weights.size();subset++)
  {
    weights[subset] = CalculateWeight(number,subset+1,mc,pdf1,pdf2);
  }
  
  // Master equations
  Double_t Sigma = 0;
  for (unsigned int i=0;i<weights.size();i+=2)
  {
    Sigma += pow( weights[i] - bestfit, 2);
  }
  Sigma /= (weights.size()-1);
  Sigma = sqrt(Sigma);
  UncertaintyType uncert_value;
  uncert_value.Max=Sigma;
  uncert_value.Min=Sigma;
  uncert_value.Mean=bestfit;
  
  return uncert_value;
}


// -----------------------------------------------------------------------------
// CalculateNNPDF
// -----------------------------------------------------------------------------
UncertaintyType PDFReweighter2::CalculateNNPDFalone(unsigned int number,
                               Int_t parton, Double_t x, Double_t Q2, 
                               Double_t& pdf) const
{
  // Computing best fit
  Double_t bestfit = CalculateWeightalone(number,0,parton,x,Q2,pdf);

  // Computing global weight related to uncertainties
  std::vector<Double_t> weights;

  if      (number==1) weights.resize(NNPDF_pdf_.second,0.);
  else if (number==2) weights.resize(NNPDF15_pdf_.second,0.);
  else if (number==3) weights.resize(NNPDF16_pdf_.second,0.);
  else if (number==4) weights.resize(NNPDF17_pdf_.second,0.);

  for (unsigned int subset=0;subset<weights.size();subset++)
  {
    weights[subset] = CalculateWeightalone(number,subset+1,parton,x,Q2,pdf);
  }
  
  // Master equations
  Double_t Sigma = 0;
  for (unsigned int i=0;i<weights.size();i+=2)
  {
    Sigma += pow( weights[i] - bestfit, 2);
  }
  Sigma /= (weights.size()-1);
  Sigma = sqrt(Sigma);
  UncertaintyType uncert_value;
  uncert_value.Max=Sigma;
  uncert_value.Min=Sigma;
  uncert_value.Mean=bestfit;
  
  return uncert_value;
}


// -----------------------------------------------------------------------------
// CalculateWeight
// -----------------------------------------------------------------------------
Double_t PDFReweighter2::CalculateWeight(unsigned int pdf, unsigned int subset,
                         const NTMonteCarlo& mc, 
                         const Double_t& pdf1, 
                         const Double_t& pdf2) const
{
  Int_t parton1 = static_cast<Int_t>(mc.partonFlavor.first);
  if (parton1==21) parton1=0;
  Int_t parton2 = static_cast<Int_t>(mc.partonFlavor.second);
  if (parton2==21) parton2=0;

  LHAPDF::usePDFMember(pdf,subset);
  Double_t newpdf1 = LHAPDF::xfx(pdf,mc.x.first, 
                                 mc.Q_scale, 
                                 parton1) / mc.x.first;
  Double_t newpdf2 = LHAPDF::xfx(pdf,mc.x.second, 
                                 mc.Q_scale, 
                                 parton2) / mc.x.second;
  return newpdf1/pdf1*newpdf2/pdf2;
}


// -----------------------------------------------------------------------------
// CalculateRefWeight
// -----------------------------------------------------------------------------
void PDFReweighter2::CalculateRefWeight(const NTMonteCarlo& mc, Double_t& pdf1, Double_t& pdf2) const
{
  Int_t parton1 = static_cast<Int_t>(mc.partonFlavor.first);
  if (parton1==21) parton1=0;
  Int_t parton2 = static_cast<Int_t>(mc.partonFlavor.second);
  if (parton2==21) parton2=0;

  LHAPDF::usePDFMember(5,0);
  pdf1 = LHAPDF::xfx(5,mc.x.first, 
                     mc.Q_scale, 
                     parton1)/mc.x.first;
  pdf2 = LHAPDF::xfx(5,mc.x.second, 
                     mc.Q_scale, 
                     parton2)/mc.x.second;
}

// -----------------------------------------------------------------------------
// CalculateWeight
// -----------------------------------------------------------------------------
Double_t PDFReweighter2::CalculateWeightalone(unsigned int pdf, 
                         unsigned int subset,
                         Int_t parton, Double_t x, Double_t Q2, 
                         const Double_t& oldpdf) const
{
  if (parton==21) parton=0;

  LHAPDF::usePDFMember(pdf,subset);
  Double_t newpdf = LHAPDF::xfx(pdf,x, 
                         Q2, 
                         parton) / x;
  return newpdf/oldpdf;
}


// -----------------------------------------------------------------------------
// CalculateRefWeight
// -----------------------------------------------------------------------------
Double_t PDFReweighter2::CalculateRefWeightalone(
                          Int_t parton, Double_t x, Double_t Q2, 
                                            Double_t& pdf) const
{
  if (parton==21) parton=0;

  LHAPDF::usePDFMember(5,0);
  Double_t newpdf = LHAPDF::xfx(5,x,
                         Q2, 
                         parton)/x;
  pdf = newpdf;
  return newpdf;
}


//#endif



// -----------------------------------------------------------------------------
// Calculate
// -----------------------------------------------------------------------------
void PDFReweighter2::CalculateAlone(Int_t parton, double x, double Q2, 
                                   double& ipdf, double& mean, double& up, double& down)const
{
  // Calculate initial pdf values
  CalculateRefWeightalone(parton,x,Q2,ipdf);
  std::cout << "ipdf=" << ipdf << std::endl;
  UncertaintyType result;

  // Calculate NNPDF weights
  std::cout << "NNPDF" << std::endl;
  UncertaintyType NNPDFvalue = CalculateNNPDFalone(1,parton,x,Q2,ipdf);
  UncertaintyType NNPDF15value = CalculateNNPDFalone(2,parton,x,Q2,ipdf);
  UncertaintyType NNPDF16value = CalculateNNPDFalone(3,parton,x,Q2,ipdf);
  UncertaintyType NNPDF17value = CalculateNNPDFalone(4,parton,x,Q2,ipdf);

  // Calculating envelope upper limit
  Double_t NNPDFup   = NNPDFvalue.Mean   + NNPDFvalue.Max;
  Double_t NNPDF15up = NNPDF15value.Mean + NNPDFvalue.Max;
  Double_t NNPDF16up = NNPDF16value.Mean + NNPDFvalue.Max;
  Double_t NNPDF17up = NNPDF17value.Mean + NNPDFvalue.Max;
  result.Max = std::max(NNPDFup,std::max(NNPDF15up,std::max(NNPDF16up,NNPDF17up)));

  // Calculating envelope lower limit
  Double_t NNPDFdown   = NNPDFvalue.Mean   - NNPDFvalue.Min;
  Double_t NNPDF15down = NNPDF15value.Mean - NNPDFvalue.Min;
  Double_t NNPDF16down = NNPDF16value.Mean - NNPDFvalue.Min;
  Double_t NNPDF17down = NNPDF17value.Mean - NNPDFvalue.Min;
  result.Min = std::min(NNPDFdown,std::min(NNPDF15down,std::min(NNPDF16down,NNPDF17down)));

  // Calculating envelope mid-point
  result.Mean = (result.Max+result.Min)/2.;

  mean = result.Mean;
  down = result.Min;
  up   = result.Max;
}
