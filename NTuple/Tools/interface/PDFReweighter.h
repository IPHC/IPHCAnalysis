#ifndef PDF_REWEIGHTER_h
#define PDF_REWEIGHTER_h

// STL headers
#include <iomanip>
#include <iostream>
#include <time.h>

// NTupleAnalysis headers
#include "NTFormat/interface/NTEvent.h"
#include "Selection/interface/DiLeptonSelection.h"
#include "Selection/interface/SelectionTable.h"
#include "Plots/interface/DiLepAnaHistoManager.h"
#include "Tools/interface/JetCorrector.h"
#include "Tools/interface/Dataset.h"
#include "Tools/interface/AnalysisEnvironmentLoader.h"

// ROOT headers
#include <Rtypes.h>

// LHAPDF headers
//#include "LHAPDF/LHAPDF.h"

namespace LHAPDF
{
  class PDF;
}

struct UncertaintyType
{
  Double_t Mean;
  Double_t Min;
  Double_t Max;

  UncertaintyType()
  { Mean=0.; Min=0.; Max=0.; }
};


class PDFReweighter
{
  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------

private:

  const LHAPDF::PDF* REF_pdf_;
  const LHAPDF::PDF* CTEQ_uncertainties_pdf_;
  const LHAPDF::PDF* CTEQ_alphas_pdf_;
  const LHAPDF::PDF* MSTW_uncertainties_pdf_;
  const LHAPDF::PDF* MSTW_alphasP_pdf_;
  const LHAPDF::PDF* MSTW_alphasM_pdf_;
  const LHAPDF::PDF* NNPDF_uncertainties_pdf_;
  const LHAPDF::PDF* NNPDF_alphasP_pdf_;
  const LHAPDF::PDF* NNPDF_alphasM_pdf_;

  // -------------------------------------------------------------
  //                    public method members
  // -------------------------------------------------------------
public:

  /// Constructor without argument
  PDFReweighter()
  {}

  /// Destructor
  ~PDFReweighter()
  {}

  /// Initialize function
  bool Initialize();

  /// Finalize function
  bool Finalize();

  /// Calculate envelop related to weights
  UncertaintyType Calculate(const IPHCTree::NTMonteCarlo& mc) const;


  // -------------------------------------------------------------
  //                    private method members
  // -------------------------------------------------------------
private:

  /// Calculate CTEQ weights
  UncertaintyType CalculateCTEQ(const IPHCTree::NTMonteCarlo& mc, 
                                Double_t& pdf1, Double_t& pdf2) const;


  /// Calculate MSTW weights
  UncertaintyType CalculateMSTW(const IPHCTree::NTMonteCarlo& mc, 
                                Double_t& pdf1, Double_t& pdf2) const;

  /// Calculate NNPDF weights
  UncertaintyType CalculateNNPDF(const IPHCTree::NTMonteCarlo& mc, 
                                 Double_t& pdf1, Double_t& pdf2) const;

  /// Generic function to calculate PDF weights
  Double_t CalculateWeight(unsigned int pdf, unsigned int subset,
                           const IPHCTree::NTMonteCarlo& mc, 
                           const Double_t& pdf1, 
                           const Double_t& pdf2) const;

  /// Calculate initial weights
  void CalculateRefWeight(const IPHCTree::NTMonteCarlo& mc, 
                          Double_t& pdf1, Double_t& pdf2) const;

};

#endif

//#endif
