#include "Tools/interface/PDFReweighter.h"

#ifdef PDFREWEIGHT_ENABLE

using namespace IPHCTree;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// -----------------------------------------------------------------------------
bool PDFReweighter::Initialize()
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

  // Initial PDF
  REF_pdf_ = std::make_pair("cteq6l.LHpdf", 0);
  LHAPDF::initPDFSet(1,REF_pdf_.first);
  REF_pdf_.second = LHAPDF::numberPDF(1);
  if (REF_pdf_.second == 1) REF_pdf_.second=0;
  std::cout << " => PDF '" << REF_pdf_.first << "' is loaded with Best Fit + " 
            << REF_pdf_.second << " subsets" << std::endl << std::endl;
  if (REF_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // CTEQ uncertainties PDF
  CTEQ_uncertainties_pdf_  = std::make_pair("cteq66.LHgrid", 0);
  LHAPDF::initPDFSet(2,CTEQ_uncertainties_pdf_.first);
  CTEQ_uncertainties_pdf_.second = LHAPDF::numberPDF(2);
  std::cout << " => PDF '" << CTEQ_uncertainties_pdf_.first << "' is loaded with Best Fit + " 
            << CTEQ_uncertainties_pdf_.second << " subsets" << std::endl << std::endl;
  if (CTEQ_uncertainties_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // CTEQ alphas PDF
  CTEQ_alphas_pdf_ = std::make_pair("cteq66alphas.LHgrid", 0);
  LHAPDF::initPDFSet(3,CTEQ_alphas_pdf_.first);
  CTEQ_alphas_pdf_.second = LHAPDF::numberPDF(3);
  std::cout << " => PDF '" << CTEQ_alphas_pdf_.first << "' is loaded with "
            << CTEQ_alphas_pdf_.second+1 << " members" << std::endl << std::endl;
  if (CTEQ_alphas_pdf_.second!=4)
  {
    std::cout << "ERROR: wrong number of subsets for CTEQ uncertainties" << std::endl;
    return false;
  }

  // MSTW uncertainties PDF
  MSTW_uncertainties_pdf_ = std::make_pair("MSTW2008nlo68cl.LHgrid", 0);
  LHAPDF::initPDFSet(4,MSTW_uncertainties_pdf_.first);
  MSTW_uncertainties_pdf_.second = LHAPDF::numberPDF(4);
  std::cout << " => PDF '" << MSTW_uncertainties_pdf_.first << "' is loaded with Best Fit + " 
            << MSTW_uncertainties_pdf_.second << " subsets" << std::endl << std::endl;
  if (MSTW_uncertainties_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // MSTW alphas up PDF
  MSTW_alphasP_pdf_ = std::make_pair("MSTW2008nlo68cl_asmz+68cl.LHgrid", 0);
  LHAPDF::initPDFSet(5,MSTW_alphasP_pdf_.first);
  MSTW_alphasP_pdf_.second = LHAPDF::numberPDF(5);
  std::cout << " => PDF '" << MSTW_alphasP_pdf_.first << "' is loaded with Best Fit + " 
            << MSTW_alphasP_pdf_.second << " subsets" << std::endl << std::endl;
  if (MSTW_alphasP_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // MSTW alphas down PDF
  MSTW_alphasM_pdf_ = std::make_pair("MSTW2008nlo68cl_asmz-68cl.LHgrid", 0);
  LHAPDF::initPDFSet(6,MSTW_alphasM_pdf_.first);
  MSTW_alphasM_pdf_.second = LHAPDF::numberPDF(6);
  std::cout << " => PDF '" << MSTW_alphasM_pdf_.first << "' is loaded with Best Fit + " 
            << MSTW_alphasM_pdf_.second << " subsets" << std::endl << std::endl;
  if (MSTW_alphasM_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // NNPDF uncertainties PDF
  NNPDF_uncertainties_pdf_ = std::make_pair("NNPDF20_100.LHgrid", 0);
  LHAPDF::initPDFSet(7,NNPDF_uncertainties_pdf_.first);
  NNPDF_uncertainties_pdf_.second = LHAPDF::numberPDF(7);
  std::cout << " => PDF '" << NNPDF_uncertainties_pdf_.first << "' is loaded with Best Fit + " 
            << NNPDF_uncertainties_pdf_.second << " subsets" << std::endl << std::endl;
  if (NNPDF_uncertainties_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // NNPDF alphas up PDF
  NNPDF_alphasP_pdf_ = std::make_pair("NNPDF20_as_0120_100.LHgrid", 0);
  LHAPDF::initPDFSet(8,NNPDF_alphasP_pdf_.first);
  NNPDF_alphasP_pdf_.second = LHAPDF::numberPDF(8);
  std::cout << " => PDF '" << NNPDF_alphasP_pdf_.first << "' is loaded with Best Fit + " 
            << NNPDF_alphasP_pdf_.second << " subsets" << std::endl << std::endl;
  if (NNPDF_alphasP_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // NNPDF alphas down PDF
  NNPDF_alphasM_pdf_ = std::make_pair("NNPDF20_as_0118_100.LHgrid", 0);
  LHAPDF::initPDFSet(9,NNPDF_alphasM_pdf_.first);
  NNPDF_alphasM_pdf_.second = LHAPDF::numberPDF(9);
  std::cout << " => PDF '" << NNPDF_alphasM_pdf_.first << "' is loaded with Best Fit + " 
            << NNPDF_alphasM_pdf_.second << " subsets" << std::endl << std::endl;
  if (NNPDF_alphasM_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  std::cout << "===============================================================" << std::endl;
  std::cout << "===============================================================" << std::endl;

  return true;
}


// -----------------------------------------------------------------------------
// Finalize
// -----------------------------------------------------------------------------
bool PDFReweighter::Finalize()
{
  return true;
}


// -----------------------------------------------------------------------------
// Calculate
// -----------------------------------------------------------------------------
UncertaintyType PDFReweighter::Calculate(const NTMonteCarlo& mc) const
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

  // Calculate CTEQ weights
  UncertaintyType CTEQvalue  = CalculateCTEQ(mc,pdf1,pdf2);
  UncertaintyType MSTWvalue  = CalculateMSTW(mc,pdf1,pdf2);
  UncertaintyType NNPDFvalue = CalculateNNPDF(mc,pdf1,pdf2);

  // Calculating envelope upper limit
  Double_t CTEQup  = CTEQvalue.Mean  + CTEQvalue.Max;
  Double_t MSTWup  = MSTWvalue.Mean  + MSTWvalue.Max;
  Double_t NNPDFup = NNPDFvalue.Mean + NNPDFvalue.Max;
  result.Max = std::max(CTEQup,std::max(MSTWup,NNPDFup));

  // Calculating envelope upper limit
  Double_t CTEQdown  = CTEQvalue.Mean  - CTEQvalue.Min;
  Double_t MSTWdown  = MSTWvalue.Mean  - MSTWvalue.Min;
  Double_t NNPDFdown = NNPDFvalue.Mean - NNPDFvalue.Min;
  result.Min = std::min(CTEQdown,std::min(MSTWdown,NNPDFdown));

  // Good boundaries ?
  if (result.Min<0.01 || result.Max>10)
  {
    std::cout << "WARNING: weight boundaries are excessive ! "
              << "Mean:" << result.Mean 
              << " Max:" << result.Max 
              << " Min:" << result.Min << std::endl;
    result.Min=0; result.Max=0;
  }

  // Calculating envelope mid-point
  result.Mean = (result.Max+result.Min)/2.;
  return result;
}


// -----------------------------------------------------------------------------
// CalculateCTEQ
// -----------------------------------------------------------------------------
UncertaintyType PDFReweighter::CalculateCTEQ(const NTMonteCarlo& mc, 
                                             Double_t& pdf1, 
                                             Double_t& pdf2) const
{
  // Computing best fit
  Double_t bestfit = CalculateWeight(2,0,mc,pdf1,pdf2);

  // Computing weights related to uncertainties
  std::vector<Double_t> weights(CTEQ_uncertainties_pdf_.second,0.);
  for (unsigned int subset=0;subset<weights.size();subset++)
  {
    weights[subset] = CalculateWeight(2,subset+1,mc,pdf1,pdf2);
  } 

  // Computing global weight related to uncertainties
  // Master equations
  Double_t DeltaP = 0;
  Double_t DeltaM = 0;
  for (unsigned int i=0;i<weights.size();i+=2)
  {
    DeltaP += pow(std::max(0., std::max(weights[i]   - bestfit,
                                        weights[i+1] - bestfit)),2);
    DeltaM += pow(std::max(0., std::max(bestfit - weights[i],
                                        bestfit - weights[i+1])),2);
  }
  DeltaP = 1./1.64485 * sqrt(DeltaP);
  DeltaM = 1./1.64485 * sqrt(DeltaM);
  UncertaintyType uncert_value;
  uncert_value.Max=DeltaP;
  uncert_value.Min=DeltaM;

  // Computing bound related to alpha_s
  UncertaintyType alphas_value;
  double myweight = CalculateWeight(3,3,mc,pdf1,pdf2);
  alphas_value.Max = (myweight - bestfit)/(5./6.);
  myweight = CalculateWeight(3,1,mc,pdf1,pdf2);
  alphas_value.Min = (myweight - bestfit)/(5./6.);

  // Combining alphas and uncertainties
  UncertaintyType combination;
  combination.Mean = bestfit;
  combination.Max = sqrt( pow(uncert_value.Max, 2) + 
                          pow(alphas_value.Max, 2) );
  combination.Min = sqrt( pow(uncert_value.Min,2) + 
                          pow(alphas_value.Min,2) );
  return combination;
}


// -----------------------------------------------------------------------------
// CalculateMSTW
// -----------------------------------------------------------------------------
UncertaintyType PDFReweighter::CalculateMSTW(const NTMonteCarlo& mc, 
                                             Double_t& pdf1, 
                                             Double_t& pdf2) const
{
  // Computing best fit
  Double_t bestfit = CalculateWeight(4,0,mc,pdf1,pdf2);

  // Computing weights related to uncertainties
  std::vector<Double_t> weights(MSTW_uncertainties_pdf_.second,0.);
  for (unsigned int subset=0;subset<weights.size();subset++)
  {
    weights[subset] = CalculateWeight(4,subset+1,mc,pdf1,pdf2);
  } 

  // Computing global weight related to uncertainties
  // Master equations
  Double_t DeltaP = 0;
  Double_t DeltaM = 0;
  for (unsigned int i=0;i<weights.size();i+=2)
  {
    DeltaP += pow(std::max(0., std::max(weights[i]   - bestfit,
                                        weights[i+1] - bestfit)),2);
    DeltaM += pow(std::max(0., std::max(bestfit - weights[i],
                                        bestfit - weights[i+1])),2);
  }
  DeltaP = sqrt(DeltaP);
  DeltaM = sqrt(DeltaM);
  UncertaintyType uncert_value;
  uncert_value.Max=DeltaP;
  uncert_value.Min=DeltaM;

  // Computing bound related to alpha_s
  UncertaintyType alphas_value;
  double myweight = CalculateWeight(5,0,mc,pdf1,pdf2);
  alphas_value.Max = (myweight - bestfit)/(5./4.);
  myweight = CalculateWeight(6,0,mc,pdf1,pdf2);
  alphas_value.Min = (myweight - bestfit)/(5./4.);

  // Combining alphas and uncertainties
  UncertaintyType combination;
  combination.Mean = bestfit;
  combination.Max  = sqrt( pow(uncert_value.Max, 2) + 
                           pow(alphas_value.Max, 2) );
  combination.Min  = sqrt( pow(uncert_value.Min,2) + 
                           pow(alphas_value.Min,2) );
  return combination;
}


// -----------------------------------------------------------------------------
// CalculateNNPDF
// -----------------------------------------------------------------------------
UncertaintyType PDFReweighter::CalculateNNPDF(const NTMonteCarlo& mc, 
                               Double_t& pdf1, Double_t& pdf2) const
{
  // Computing best fit
  Double_t bestfit = CalculateWeight(7,0,mc,pdf1,pdf2);

  // Computing weights related to uncertainties
  std::vector<Double_t> weights(NNPDF_uncertainties_pdf_.second,0.);
  for (unsigned int subset=0;subset<weights.size();subset++)
  {
    weights[subset] = CalculateWeight(7,subset+1,mc,pdf1,pdf2);
  } 

  // Computing global weight related to uncertainties
  // Master equations
  Double_t DeltaP = 0;
  Double_t DeltaM = 0;
  for (unsigned int i=0;i<weights.size();i+=2)
  {
    DeltaP += pow(std::max(0., std::max(weights[i]   - bestfit,
                                        weights[i+1] - bestfit)),2);
    DeltaM += pow(std::max(0., std::max(bestfit - weights[i],
                                        bestfit - weights[i+1])),2);
  }
  DeltaP = sqrt(DeltaP);
  DeltaM = sqrt(DeltaM);
  UncertaintyType uncert_value;
  uncert_value.Max=DeltaP;
  uncert_value.Min=DeltaM;

  // Computing bound related to alpha_s
  UncertaintyType alphas_value;
  double myweight = CalculateWeight(8,0,mc,pdf1,pdf2);
  alphas_value.Max = (myweight - bestfit)/(5./4.);
  myweight = CalculateWeight(9,0,mc,pdf1,pdf2);
  alphas_value.Min = (myweight - bestfit)/(5./4.);

  // Combining alphas and uncertainties
  UncertaintyType combination;
  combination.Mean = bestfit;
  combination.Max  = sqrt( pow(uncert_value.Max, 2) + 
                           pow(alphas_value.Max, 2) );
  combination.Min  = sqrt( pow(uncert_value.Min,2) + 
                           pow(alphas_value.Min,2) );
  return combination;
}


// -----------------------------------------------------------------------------
// CalculateWeight
// -----------------------------------------------------------------------------
Double_t PDFReweighter::CalculateWeight(unsigned int pdf, unsigned int subset,
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
void PDFReweighter::CalculateRefWeight(const NTMonteCarlo& mc, Double_t& pdf1, Double_t& pdf2) const
{
  Int_t parton1 = static_cast<Int_t>(mc.partonFlavor.first);
  if (parton1==21) parton1=0;
  Int_t parton2 = static_cast<Int_t>(mc.partonFlavor.second);
  if (parton2==21) parton2=0;

  LHAPDF::usePDFMember(1,0);
  pdf1 = LHAPDF::xfx(1,mc.x.first, 
                     mc.Q_scale, 
                     parton1)/mc.x.first;
  pdf2 = LHAPDF::xfx(1,mc.x.second, 
                     mc.Q_scale, 
                     parton2)/mc.x.second;
}

#endif
