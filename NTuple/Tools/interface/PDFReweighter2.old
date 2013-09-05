#ifndef PDF_REWEIGHTER2_h
#define PDF_REWEIGHTER2_h

//#ifdef PDFREWEIGHT_ENABLE

// STL headers
#include <iomanip>
#include <iostream>
#include <time.h>

// NTupleAnalysis headers
#include "NTFormat/interface/NTEvent.h"
#include "Selection/interface/SSDiLeptonSelection.h"
#include "Selection/interface/SelectionTable.h"
#include "Plots/interface/SSDiLepAnaHistoManager.h"
#include "Tools/interface/JetCorrector.h"
#include "Tools/interface/Dataset.h"
#include "Tools/interface/AnalysisEnvironmentLoader.h"
#include "Tools/interface/PDFReweighter.h"

// ROOT headers
#include <Rtypes.h>

// LHAPDF headers
#include "LHAPDF/include/LHAPDF/LHAPDF.h"


using namespace IPHCTree;
using namespace std;



class PDFReweighter2
{
  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------

private:

  std::pair<std::string,unsigned int> REF_pdf_;
  std::pair<std::string,unsigned int> NNPDF_pdf_;
  std::pair<std::string,unsigned int> NNPDF15_pdf_;
  std::pair<std::string,unsigned int> NNPDF16_pdf_;
  std::pair<std::string,unsigned int> NNPDF17_pdf_;

  // -------------------------------------------------------------
  //                    public method members
  // -------------------------------------------------------------
public:

  /// Constructor without argument
  PDFReweighter2()
  {}

  /// Destructor
  ~PDFReweighter2()
  {}

  /// Initialize function
  bool Initialize();

  /// Finalize function
  bool Finalize();

  /// Calculate envelop related to weights
  UncertaintyType Calculate(const IPHCTree::NTMonteCarlo& mc, unsigned int mode) const;


void CalculateAlone(Int_t parton, double x, double Q2, 
                    double& ipdf, double& mean, double& up, double& down)const;


  // -------------------------------------------------------------
  //                    private method members
  // -------------------------------------------------------------
private:

  /// Calculate NNPDF weights
  UncertaintyType CalculateNNPDF(unsigned int number, const IPHCTree::NTMonteCarlo& mc, 
                                 Double_t& pdf1, Double_t& pdf2) const;

  /// Calculate NNPDF weights
  UncertaintyType CalculateNNPDFalone(unsigned int number, Int_t parton, Double_t x, Double_t Q2, 
                                Double_t& pdf) const;

  /// Generic function to calculate PDF weights
  Double_t CalculateWeight(unsigned int pdf, unsigned int subset,
                           const IPHCTree::NTMonteCarlo& mc, 
                           const Double_t& pdf1, 
                           const Double_t& pdf2) const;

  /// Calculate initial weights
  void CalculateRefWeight(const IPHCTree::NTMonteCarlo& mc, 
                          Double_t& pdf1, Double_t& pdf2) const;

  /// Generic function to calculate PDF weights
  Double_t CalculateWeightalone(unsigned int oldpdf, 
                                unsigned int subset,
                           Int_t parton, Double_t x, Double_t Q2, 
                           const Double_t& pdf) const;

  /// Calculate initial weights
  Double_t CalculateRefWeightalone(
                          Int_t parton, Double_t x, Double_t Q2, 
                          Double_t& pdf) const;


};

#endif

//#endif
