#include "Tools/interface/PDFReweighter.h"

//#ifdef PDFREWEIGHT_ENABLE

using namespace IPHCTree;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// -----------------------------------------------------------------------------
bool PDFReweighter::Initialize()
{
  cteq  = new TH1F("cteq","cteq",100,0,5);
  mstw  = new TH1F("mstw","mstw",100,0,5);
  nnpdf = new TH1F("nnpdf","nnpdf",100,0,5);
  cteq2  = new TH1F("cteq2","cteq2",100,0,5);
  mstw2  = new TH1F("mstw2","mstw2",100,0,5);
  nnpdf2 = new TH1F("nnpdf2","nnpdf2",100,0,5);
  cteq0  = new TH1F("cteq0","cteq0",100,0,5);
  mstw0  = new TH1F("mstw0","mstw0",100,0,5);
  nnpdf0 = new TH1F("nnpdf0","nnpdf0",100,0,5);
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
  CTEQ_uncertainties_pdf_  = std::make_pair("cteq66.LHgrid", 0);
  LHAPDF::initPDFSet(1,CTEQ_uncertainties_pdf_.first);
  CTEQ_uncertainties_pdf_.second = LHAPDF::numberPDF(1);
  std::cout << " => PDF '" << CTEQ_uncertainties_pdf_.first << "' is loaded with Best Fit + " 
            << CTEQ_uncertainties_pdf_.second << " subsets" << std::endl << std::endl;
  if (CTEQ_uncertainties_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // CTEQ alphas PDF
  CTEQ_alphas_pdf_ = std::make_pair("cteq66alphas.LHgrid", 0);
  LHAPDF::initPDFSet(2,CTEQ_alphas_pdf_.first);
  CTEQ_alphas_pdf_.second = LHAPDF::numberPDF(2);
  std::cout << " => PDF '" << CTEQ_alphas_pdf_.first << "' is loaded with "
            << CTEQ_alphas_pdf_.second+1 << " members" << std::endl << std::endl;
  if (CTEQ_alphas_pdf_.second!=4)
  {
    std::cout << "ERROR: wrong number of subsets for CTEQ uncertainties" << std::endl;
     return false;
  }

  // MSTW uncertainties PDF
  MSTW_uncertainties_pdf_ = std::make_pair("MSTW2008nlo68cl.LHgrid", 0);
  LHAPDF::initPDFSet(3,MSTW_uncertainties_pdf_.first);
  MSTW_uncertainties_pdf_.second = LHAPDF::numberPDF(3);
  std::cout << " => PDF '" << MSTW_uncertainties_pdf_.first << "' is loaded with Best Fit + " 
            << MSTW_uncertainties_pdf_.second << " subsets" << std::endl << std::endl;
  if (MSTW_uncertainties_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // MSTW alphas up PDF
  MSTW_alphasP_pdf_ = std::make_pair("MSTW2008nlo68cl_asmz+68cl.LHgrid", 0);
  LHAPDF::initPDFSet(4,MSTW_alphasP_pdf_.first);
  MSTW_alphasP_pdf_.second = LHAPDF::numberPDF(4);
  std::cout << " => PDF '" << MSTW_alphasP_pdf_.first << "' is loaded with Best Fit + " 
            << MSTW_alphasP_pdf_.second << " subsets" << std::endl << std::endl;
  if (MSTW_alphasP_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // MSTW alphas down PDF
  MSTW_alphasM_pdf_ = std::make_pair("MSTW2008nlo68cl_asmz-68cl.LHgrid", 0);
  LHAPDF::initPDFSet(5,MSTW_alphasM_pdf_.first);
  MSTW_alphasM_pdf_.second = LHAPDF::numberPDF(5);
  std::cout << " => PDF '" << MSTW_alphasM_pdf_.first << "' is loaded with Best Fit + " 
            << MSTW_alphasM_pdf_.second << " subsets" << std::endl << std::endl;
  if (MSTW_alphasM_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // NNPDF uncertainties PDF
  NNPDF_uncertainties_pdf_ = std::make_pair("NNPDF20_100.LHgrid", 0);
  LHAPDF::initPDFSet(6,NNPDF_uncertainties_pdf_.first);
  NNPDF_uncertainties_pdf_.second = LHAPDF::numberPDF(6);
  std::cout << " => PDF '" << NNPDF_uncertainties_pdf_.first << "' is loaded with Best Fit + " 
            << NNPDF_uncertainties_pdf_.second << " subsets" << std::endl << std::endl;
  if (NNPDF_uncertainties_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // NNPDF alphas up PDF
  NNPDF_alphasP_pdf_ = std::make_pair("NNPDF20_as_0120_100.LHgrid", 0);
  LHAPDF::initPDFSet(7,NNPDF_alphasP_pdf_.first);
  NNPDF_alphasP_pdf_.second = LHAPDF::numberPDF(7);
  std::cout << " => PDF '" << NNPDF_alphasP_pdf_.first << "' is loaded with Best Fit + " 
            << NNPDF_alphasP_pdf_.second << " subsets" << std::endl << std::endl;
  if (NNPDF_alphasP_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // NNPDF alphas down PDF
  NNPDF_alphasM_pdf_ = std::make_pair("NNPDF20_as_0118_100.LHgrid", 0);
  LHAPDF::initPDFSet(8,NNPDF_alphasM_pdf_.first);
  NNPDF_alphasM_pdf_.second = LHAPDF::numberPDF(8);
  std::cout << " => PDF '" << NNPDF_alphasM_pdf_.first << "' is loaded with Best Fit + " 
            << NNPDF_alphasM_pdf_.second << " subsets" << std::endl << std::endl;
  if (NNPDF_alphasM_pdf_.second%2!=0)
  {
    std::cout << "ERROR: odd-number of subsets" << std::endl;
    return false;
  }

  // Initial PDF
  REF_pdf_ = std::make_pair("cteq6l.LHpdf", 0);
  //REF_pdf_ = std::make_pair("cteq6ll.LHpdf", 0);
  LHAPDF::initPDFSet(9,REF_pdf_.first);
  REF_pdf_.second = LHAPDF::numberPDF(9);
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
bool PDFReweighter::Finalize()
{
  std::cout << "ERIC FINALIZE" << std::endl;

  cteq0->SetLineColor(kRed);
  mstw0->SetLineColor(kBlue);
  nnpdf0->SetLineColor(kGreen);
  cteq->SetLineColor(kRed);
  mstw->SetLineColor(kBlue);
  nnpdf->SetLineColor(kGreen);
  cteq2->SetLineColor(kRed);
  mstw2->SetLineColor(kBlue);
  nnpdf2->SetLineColor(kGreen);

  TCanvas * eric1 = new TCanvas ("eric_max","eric_max");
  THStack * stack1 = new THStack();
  stack1->Add(cteq);
  stack1->Add(mstw);
  stack1->Add(nnpdf);
  stack1->Draw("nostack");
  eric1->Write();

  TCanvas * eric2 = new TCanvas ("eric_min","eric_min");
  THStack * stack2 = new THStack();
  stack2->Add(cteq2);
  stack2->Add(mstw2);
  stack2->Add(nnpdf2);
  stack2->Draw("nostack");
  eric2->Write();

  TCanvas * eric3 = new TCanvas ("eric_mean","eric_mean");
  THStack * stack3 = new THStack();
  stack3->Add(cteq0);
  stack3->Add(mstw0);
  stack3->Add(nnpdf0);
  stack3->Draw("nostack");
  eric3->Write();

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
  std::cout << "CTEQ" << std::endl;
  UncertaintyType CTEQvalue = CalculateCTEQ(mc,pdf1,pdf2);
  cteq ->Fill(CTEQvalue.Max);
  cteq2->Fill(CTEQvalue.Min);
  cteq0->Fill(CTEQvalue.Mean);
  
  // Calculate MSTW weights
  std::cout << "LSTW" << std::endl;
  UncertaintyType MSTWvalue = CalculateMSTW(mc,pdf1,pdf2);
  mstw ->Fill(MSTWvalue.Max);
  mstw2->Fill(MSTWvalue.Min);
  mstw0->Fill(MSTWvalue.Mean);
  
  // Calculate NNPDF weights
  std::cout << "NNPDF" << std::endl;
  UncertaintyType NNPDFvalue = CalculateNNPDF(mc,pdf1,pdf2);
  nnpdf ->Fill(NNPDFvalue.Max);
  nnpdf2->Fill(NNPDFvalue.Min);
  nnpdf0->Fill(NNPDFvalue.Mean);

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
  Double_t bestfit = CalculateWeight(1,0,mc,pdf1,pdf2);

  // Computing weights related to uncertainties
  std::vector<Double_t> weights(CTEQ_uncertainties_pdf_.second,0.);
  for (unsigned int subset=0;subset<weights.size();subset++)
  {
    weights[subset] = CalculateWeight(1,subset+1,mc,pdf1,pdf2);
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
  uncert_value.Mean=bestfit;

  // Computing bound related to alpha_s
  UncertaintyType alphas_value;
  double myweight = CalculateWeight(2,3,mc,pdf1,pdf2);
  alphas_value.Max = (myweight - bestfit)/(5./6.);
  myweight = CalculateWeight(2,1,mc,pdf1,pdf2);
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
  Double_t bestfit = CalculateWeight(3,0,mc,pdf1,pdf2);

  // Computing weights related to uncertainties
  std::cout << "size=" << MSTW_uncertainties_pdf_.second << std::endl;
  std::vector<Double_t> weights(MSTW_uncertainties_pdf_.second,0.);
  for (unsigned int subset=0;subset<weights.size();subset++)
  {
    weights[subset] = CalculateWeight(3,subset+1,mc,pdf1,pdf2);
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
  std::cout << "4" << std::endl;
  
  DeltaP = sqrt(DeltaP);
  DeltaM = sqrt(DeltaM);
  UncertaintyType uncert_value;
  uncert_value.Max=DeltaP;
  uncert_value.Min=DeltaM;
  uncert_value.Mean=bestfit;

  // Computing bound related to alpha_s
  UncertaintyType alphas_value;
  double myweight = CalculateWeight(4,0,mc,pdf1,pdf2);
  alphas_value.Max = (myweight - bestfit)/(5./4.);
  myweight = CalculateWeight(5,0,mc,pdf1,pdf2);
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
  Double_t bestfit = CalculateWeight(6,0,mc,pdf1,pdf2);

  // Computing global weight related to uncertainties
  std::vector<Double_t> weights(NNPDF_uncertainties_pdf_.second,0.);
  for (unsigned int subset=0;subset<weights.size();subset++)
  {
    weights[subset] = CalculateWeight(6,subset+1,mc,pdf1,pdf2);
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
  
  // Computing bound related to alpha_s
  UncertaintyType alphas_value;
  double myweight = CalculateWeight(7,0,mc,pdf1,pdf2);
  alphas_value.Max = (myweight - bestfit);
  myweight = CalculateWeight(8,0,mc,pdf1,pdf2);
  alphas_value.Min = (myweight - bestfit);

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

  LHAPDF::usePDFMember(9,0);
  pdf1 = LHAPDF::xfx(9,mc.x.first, 
                     mc.Q_scale, 
                     parton1)/mc.x.first;
  pdf2 = LHAPDF::xfx(9,mc.x.second, 
                     mc.Q_scale, 
                     parton2)/mc.x.second;
}

//#endif
