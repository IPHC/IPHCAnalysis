#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include <iostream>

using namespace std;


static bool nomatching_scale_singal = false;

void prodtemplate(std::string signalname, 
                  std::string outputname,
                  bool WZ_PRIVATE, 
                  bool WZ_SYS,
                  bool TZq_SYS,
                  bool doJES, bool doJER, bool doBTag, bool doPU, bool doLept, 
                  bool doTopMass, bool doPDF, bool doScale, bool doMatch, bool doDY);


int main(int argc,char *argv[])
{

  std::cout << "Welcome to ProdTemplate macro" << std::endl << std::endl;
  std::cout << "By default, all the boolean parameters are set to true." << std::endl;

  std::vector<std::string> arguments;
  if (argc>1)
  {
    for (int i=1;i<argc;i++) arguments.push_back(argv[i]);
  }

  bool WZ_PRIVATE  = false;
  bool WZ_SYS      = true;
  bool TZq_SYS     = true;

  bool doJES=true;
  bool doJER=true;
  bool doBTag=true;
  bool doPU=true;
  bool doLept=true;
  bool doTopMass=true;
  bool doPDF=true; //dfd
  bool doScale=true;
  bool doMatch=true;
  bool doDY=true;
  std::string signalname="kct";
  std::string outputname;
  

  
  if(signalname == "zut" ) outputname="ThetaFile_zut.root";
  if(signalname == "zct" ) outputname="ThetaFile_zct.root";
  if(signalname == "kut" ) outputname="ThetaFile_kut.root";
  if(signalname == "kct" ) outputname="ThetaFile_kct.root";
  if(signalname == "tzq" ) outputname="ThetaFile_tzq.root";

  for (unsigned int i=0;i<arguments.size();i++)
  {
    if (arguments[i]=="WZ_PRIVATE=true")       WZ_PRIVATE=true;
    else if (arguments[i]=="WZ_PRIVATE=false") WZ_PRIVATE=false;
    else if (arguments[i]=="WZ_SYS=true")      WZ_SYS=true;
    else if (arguments[i]=="WZ_SYS=false")     WZ_SYS=false;
    else if (arguments[i]=="TZq_SYS=true")     TZq_SYS=true;
    else if (arguments[i]=="TZq_SYS=false")    TZq_SYS=false;
    else if (arguments[i]=="doJES=false")      doJES=false;
    else if (arguments[i]=="doJER=false")      doJER=false;
    else if (arguments[i]=="doBTag=false")     doBTag=false;
    else if (arguments[i]=="doPU=false")       doPU=false;
    else if (arguments[i]=="doLept=false")     doLept=false;
    else if (arguments[i]=="doTopMass=false")  doTopMass=false;
    else if (arguments[i]=="doPDF=false")      doPDF=false;
    else if (arguments[i]=="doScale=false")    doScale=false;
    else if (arguments[i]=="doMatch=false")    doMatch=false;
    else if (arguments[i]=="doDY=false")       doDY=false;
    else if (arguments[i]=="signal=zut")       signalname="zut";
    else if (arguments[i]=="signal=zct")       signalname="zct";
    else if (arguments[i]=="signal=kut")       signalname="kut";
    else if (arguments[i]=="signal=kct")       signalname="kct";
    else if (arguments[i].find("output=")==0)
    {
      outputname=arguments[i].substr(7,std::string::npos);
    }
    else 
    {
      std::cout << "ERROR: the parameter '" << arguments[i] 
                << "' is unknown. Skip it!" << std::endl;
    }
  }

  // Display configuration
  std::cout << "--------------------------------------" << std::endl;
  std::cout << "Configuration: " << std::endl;
  std::cout << std::endl;
  std::cout << " Signal name = " << signalname << std::endl;
  std::cout << " Output name = " << outputname << std::endl;
  std::cout << std::endl;
  std::cout << " WZ_PRIVATE=" << WZ_PRIVATE << "   WZ_SYS ="  << WZ_SYS  << "   TZq_SYS  =" << TZq_SYS   << std::endl;
  std::cout << " doJES     =" << doJES      << "   doJER  ="  << doJER   << "   doBTag   =" << doBTag    << std::endl;
  std::cout << " doPU      =" << doPU       << "   doLept ="  << doLept  << "   doTopMass=" << doTopMass << std::endl;
  std::cout << " doPDF     =" << doPDF      << "   doScale=" << doScale << "   doMatch  =" << doMatch   << std::endl; 
  std::cout << "--------------------------------------" << std::endl;

  prodtemplate(signalname,outputname,WZ_PRIVATE, WZ_SYS, TZq_SYS, doJES, doJER, doBTag, doPU, doLept, doTopMass, doPDF, doScale, doMatch, doDY);


  return 0;

}


void prodtemplate(std::string signalname, std::string outputname,
                  bool WZ_PRIVATE, 
                  bool WZ_SYS,
                  bool TZq_SYS,
                  bool doJES, bool doJER, bool doBTag, bool doPU, bool doLept, 
                  bool doTopMass, bool doPDF, bool doScale, bool doMatch, bool doDY)
{
  // Initializing scale factors
  double WZscaleFactor = 0.;
  if (WZ_PRIVATE) WZscaleFactor = 1.0;
  else WZscaleFactor = 1.0;
  double TZqscaleFactor  = 0.27;
  //  double datasf          = 4.0;
  double ZjetscaleFactor = 1.0;
  double ZZscaleFactor   = 1.0;
  double TTbarscaleFactor= 1.0;

  // Initializing WZ name
  std::string WZname;
  if (WZ_PRIVATE) WZname = "WZprivate";
  else WZname = "WZ";
  std::string MVA_WZname="MVA_BDT_"+WZname;

  //----------------------------------------------------------------------------------------
  //   STEP 1 : Opening ROOT files
  //----------------------------------------------------------------------------------------
  std::cout << "Opening ROOT files ..." << std::endl;

  // Opening file  
  TFile * nomFile       = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_nom.root").c_str());
  if (!nomFile->IsOpen()) return;

  TFile * nomJES_up     = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_JESup.root").c_str());
  if (!nomJES_up->IsOpen()) return;

  TFile * nomJES_down   = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_JESdown.root").c_str()  );
  if (!nomJES_down->IsOpen()) return;

  TFile * nomJER        = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_JER.root").c_str()      );
  if (!nomJER->IsOpen()) return;

  TFile * nomLeptUp     = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_LeptSFup.root").c_str()   );
  if (!nomLeptUp->IsOpen()) return;

  TFile * nomLeptDown   = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_LeptSFdown.root").c_str()   );
  if (!nomLeptDown->IsOpen()) return;

  TFile * nomPU_up      = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_PUup.root").c_str()     );
  if (!nomPU_up->IsOpen()) return;

  TFile * nomPU_down    = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_PUdown.root").c_str()   );
  if (!nomPU_down->IsOpen()) return;

  TFile * nomBTag_up    = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_btagUp.root").c_str()   );
  if (!nomBTag_up->IsOpen()) return;

  TFile * nomBTag_down  = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_btagdown.root").c_str() );
  if (!nomBTag_down->IsOpen()) return;

  TFile * nomScale_up   = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_Scaleup.root").c_str()    );
  if (!nomScale_up->IsOpen()) return;

  TFile * nomScale_down = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_Scaledown.root").c_str()  );
  if (!nomScale_down->IsOpen()) return;

  TFile * nomMatch_up   = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_Matchup.root").c_str()     );
  if (!nomMatch_up->IsOpen()) return;

  TFile * nomMatch_down = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_Matchdown.root").c_str()   );
  if (!nomMatch_down->IsOpen()) return;

  TFile * nomMtop_up    = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_Mtopup.root").c_str()   );
  if (!nomMtop_up->IsOpen()) return;

  TFile * nomMtop_down  = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_Mtopdown.root").c_str() );
  if (!nomMtop_down->IsOpen()) return;

  TFile * nomPDF_up     = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_PDFup.root").c_str()   );
  //TFile * nomPDF_up     = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_nom.root").c_str()   );
  if (!nomPDF_up->IsOpen()) return;

  TFile * nomPDF_down   = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_PDFdown.root").c_str() );
  //TFile * nomPDF_down   = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_nom.root").c_str() );
  if (!nomPDF_down->IsOpen()) return;
  
  
  
  TFile * nomDY_up    = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_DYUp.root").c_str()   );
  if (!nomDY_up->IsOpen()) return;

  TFile * nomDY_down  = new TFile(("../TMVA/HistoBDToutput/TMVApp_"+signalname+"_DYDown.root").c_str() );
  if (!nomDY_down->IsOpen()) return;


  
  
  
  
  
  
  //nominal templates 
  TH1F histo_Data;
  TH1F histo_Zjets;
  TH1F histo_ZjetsMC;
  TH1F histo_WZ ;
  TH1F histo_TZq;
  TH1F histo_Sig;
  TH1F histo_TTbar;
  TH1F histo_ZZ;       
  
  //templates for JES
  TH1F  histo_TTbar_JES_up;
  TH1F  histo_TTbar_JES_down;
  TH1F  histo_ZZ_JES_up;
  TH1F  histo_ZZ_JES_down;
  TH1F  histo_WZ_JES_up;
  TH1F  histo_WZ_JES_down;
  TH1F  histo_TZq_JES_up;
  TH1F  histo_TZq_JES_down;
  TH1F  histo_Sig_JES_up;
  TH1F  histo_Sig_JES_down;
  
  //templates for JER
  TH1F  histo_TTbar_JER_up;
  TH1F  histo_ZZ_JER_up;
  TH1F  histo_WZ_JER_up;
  TH1F  histo_TZq_JER_up;
  TH1F  histo_Sig_JER_up;
  TH1F  histo_TTbar_JER_down;
  TH1F  histo_ZZ_JER_down;
  TH1F  histo_WZ_JER_down;
  TH1F  histo_TZq_JER_down;
  TH1F  histo_Sig_JER_down;
  
  //templates for BTag
  TH1F  histo_TTbar_BTag_up;
  TH1F  histo_TTbar_BTag_down;
  TH1F  histo_ZZ_BTag_up;
  TH1F  histo_ZZ_BTag_down;
  TH1F  histo_WZ_BTag_up;
  TH1F  histo_WZ_BTag_down;
  TH1F  histo_TZq_BTag_up;
  TH1F  histo_TZq_BTag_down;
  TH1F  histo_Sig_BTag_up;
  TH1F  histo_Sig_BTag_down;
  
  //templates for PU
  TH1F  histo_TTbar_PU_up;
  TH1F  histo_TTbar_PU_down;
  TH1F  histo_ZZ_PU_up;
  TH1F  histo_ZZ_PU_down;
  TH1F  histo_WZ_PU_up;
  TH1F  histo_WZ_PU_down;
  TH1F  histo_TZq_PU_up;
  TH1F  histo_TZq_PU_down;
  TH1F  histo_Sig_PU_up;
  TH1F  histo_Sig_PU_down;
  
  //templates for LeptSF
  TH1F  histo_TTbar_Lept_up;
  TH1F  histo_TTbar_Lept_down;
  TH1F  histo_ZZ_Lept_up;
  TH1F  histo_ZZ_Lept_down;
  TH1F  histo_WZ_Lept_up;
  TH1F  histo_WZ_Lept_down;
  TH1F  histo_TZq_Lept_up;
  TH1F  histo_TZq_Lept_down;
  TH1F  histo_Sig_Lept_up;
  TH1F  histo_Sig_Lept_down;

  //templates for mtop
  TH1F histo_Sig_Mtop_up;
  TH1F histo_Sig_Mtop_down;

  //templates for pdf
  TH1F histo_Sig_PDF_up;
  TH1F histo_Sig_PDF_down;
  TH1F histo_WZ_PDF_up;
  TH1F histo_WZ_PDF_down;
  TH1F histo_TZq_PDF_up;
  TH1F histo_TZq_PDF_down;

  //templates for scale
  TH1F histo_Sig_Scale_up;
  TH1F histo_Sig_Scale_down;
  TH1F histo_WZ_Scale_up;
  TH1F histo_WZ_Scale_down;

  //templates for match
  TH1F histo_Sig_Match_up;
  TH1F histo_Sig_Match_down;
  TH1F histo_WZ_Match_up;
  TH1F histo_WZ_Match_down;
  TH1F histo_TZq_Match_up;
  TH1F histo_TZq_Match_down;
  
  //templates for DY
  TH1F histo_Zjets_DY_up;
  TH1F histo_Zjets_DY_down;
   
  //nominal templates 
   histo_Data.Sumw2();
   histo_Zjets.Sumw2();
   histo_ZjetsMC.Sumw2();
   histo_WZ .Sumw2();
   histo_TZq.Sumw2();
   histo_Sig.Sumw2();
   histo_TTbar.Sumw2();
   histo_ZZ.Sumw2();       
  
  //templates for JES
    histo_TTbar_JES_up.Sumw2();
    histo_TTbar_JES_down.Sumw2();
    histo_ZZ_JES_up.Sumw2();
    histo_ZZ_JES_down.Sumw2();
    histo_WZ_JES_up.Sumw2();
    histo_WZ_JES_down.Sumw2();
    histo_TZq_JES_up.Sumw2();
    histo_TZq_JES_down.Sumw2();
    histo_Sig_JES_up.Sumw2();
    histo_Sig_JES_down.Sumw2();
  
  //templates for JER
    histo_TTbar_JER_up.Sumw2();
    histo_ZZ_JER_up.Sumw2();
    histo_WZ_JER_up.Sumw2();
    histo_TZq_JER_up.Sumw2();
    histo_Sig_JER_up.Sumw2();
    histo_TTbar_JER_down.Sumw2();
    histo_ZZ_JER_down.Sumw2();
    histo_WZ_JER_down.Sumw2();
    histo_TZq_JER_down.Sumw2();
    histo_Sig_JER_down.Sumw2();
  
  //templates for BTag
    histo_TTbar_BTag_up.Sumw2();
    histo_TTbar_BTag_down.Sumw2();
    histo_ZZ_BTag_up.Sumw2();
    histo_ZZ_BTag_down.Sumw2();
    histo_WZ_BTag_up.Sumw2();
    histo_WZ_BTag_down.Sumw2();
    histo_TZq_BTag_up.Sumw2();
    histo_TZq_BTag_down.Sumw2();
    histo_Sig_BTag_up.Sumw2();
    histo_Sig_BTag_down.Sumw2();
  
  //templates for PU
    histo_TTbar_PU_up.Sumw2();
    histo_TTbar_PU_down.Sumw2();
    histo_ZZ_PU_up.Sumw2();
    histo_ZZ_PU_down.Sumw2();
    histo_WZ_PU_up.Sumw2();
    histo_WZ_PU_down.Sumw2();
    histo_TZq_PU_up.Sumw2();
    histo_TZq_PU_down.Sumw2();
    histo_Sig_PU_up.Sumw2();
    histo_Sig_PU_down.Sumw2();
  
  //templates for LeptSF
    histo_TTbar_Lept_up.Sumw2();
    histo_TTbar_Lept_down.Sumw2();
    histo_ZZ_Lept_up.Sumw2();
    histo_ZZ_Lept_down.Sumw2();
    histo_WZ_Lept_up.Sumw2();
    histo_WZ_Lept_down.Sumw2();
    histo_TZq_Lept_up.Sumw2();
    histo_TZq_Lept_down.Sumw2();
    histo_Sig_Lept_up.Sumw2();
    histo_Sig_Lept_down.Sumw2();

  //templates for mtop
   histo_Sig_Mtop_up.Sumw2();
   histo_Sig_Mtop_down.Sumw2();

  //templates for pdf
   histo_Sig_PDF_up.Sumw2();
   histo_Sig_PDF_down.Sumw2();
   histo_WZ_PDF_up.Sumw2();
   histo_WZ_PDF_down.Sumw2();
   histo_TZq_PDF_up.Sumw2();
   histo_TZq_PDF_down.Sumw2();

  //templates for scale
   histo_Sig_Scale_up.Sumw2();
   histo_Sig_Scale_down.Sumw2();
   histo_WZ_Scale_up.Sumw2();
   histo_WZ_Scale_down.Sumw2();

  //templates for match
   histo_Sig_Match_up.Sumw2();
   histo_Sig_Match_down.Sumw2();
   histo_WZ_Match_up.Sumw2();
   histo_WZ_Match_down.Sumw2();
   histo_TZq_Match_up.Sumw2();
   histo_TZq_Match_down.Sumw2();
   
   
  //templates for DY
   histo_Zjets_DY_up.Sumw2();
   histo_Zjets_DY_down.Sumw2();

  //----------------------------------------------------------------------------------------
  //   STEP 2 : Getting histograms from ROOT files
  //----------------------------------------------------------------------------------------
  std::cout << "Getting histograms from ROOT files ..." << std::endl;
  std::cout << " - nominal" << std::endl;
  nomFile->cd();
  
  histo_Data = *(TH1F*)nomFile->Get("MVA_BDT_Data");
  //histo_Data.Scale(datasf);
  histo_Zjets               = *(TH1F*)nomFile->Get("MVA_BDT_DataZjets");
  histo_ZjetsMC             = *(TH1F*)nomFile->Get("MVA_BDT_Zjets");
  histo_Zjets.Scale(histo_ZjetsMC.Integral()/histo_Zjets.Integral());
  histo_Zjets.Scale(ZjetscaleFactor);
  
  histo_WZ                  = *(TH1F*)nomFile->Get(MVA_WZname.c_str());
  histo_WZ.Scale(WZscaleFactor);  
  histo_Sig                 = *(TH1F*)nomFile->Get(("MVA_BDT_FCNC_"+signalname).c_str());
  histo_TTbar               = *(TH1F*)nomFile->Get("MVA_BDT_TTbarSig");
  histo_TTbar.Scale(TTbarscaleFactor);
  histo_ZZ                  = *(TH1F*)nomFile->Get("MVA_BDT_ZZ");
  histo_ZZ.Scale(ZZscaleFactor);
  histo_TZq                 = *(TH1F*)nomFile->Get("MVA_BDT_TZq");  
  histo_TZq.Scale(TZqscaleFactor);  
  
  /*
  histo_Data = histo_WZ;
  histo_Data.Add(&histo_TTbar);
  histo_Data.Add(&histo_ZZ);
  histo_Data.Add(&histo_Zjets);
  histo_Data.Add(&histo_TZq);
  */
  
  //template for JES up
  std::cout << " - JES up" << std::endl;
  nomJES_up->cd();
  
  histo_TTbar_JES_up = *(TH1F*)nomJES_up->Get("MVA_BDT_TTbarSig");
  histo_ZZ_JES_up    = *(TH1F*)nomJES_up->Get("MVA_BDT_ZZ"   );
  histo_WZ_JES_up    = *(TH1F*)nomJES_up->Get(MVA_WZname.c_str()   );
  histo_WZ_JES_up.Scale(WZscaleFactor);  
  histo_TZq_JES_up   = *(TH1F*)nomJES_up->Get("MVA_BDT_TZq");
  histo_TZq_JES_up.Scale(TZqscaleFactor);  
  histo_Sig_JES_up   = *(TH1F*)nomJES_up->Get(("MVA_BDT_FCNC_"+signalname).c_str() );

  
  //template for JES down
  std::cout << " - JES down" << std::endl;
  nomJES_down->cd();
  
  histo_TTbar_JES_down = *(TH1F*)nomJES_down->Get("MVA_BDT_TTbarSig");
  histo_ZZ_JES_down    = *(TH1F*)nomJES_down->Get("MVA_BDT_ZZ"   );
  histo_WZ_JES_down    = *(TH1F*)nomJES_down->Get(MVA_WZname.c_str()   );
  histo_WZ_JES_down.Scale(WZscaleFactor);  
  histo_TZq_JES_down   = *(TH1F*)nomJES_down->Get("MVA_BDT_TZq");
  histo_TZq_JES_down.Scale(TZqscaleFactor);  
  histo_Sig_JES_down   = *(TH1F*)nomJES_down->Get(("MVA_BDT_FCNC_"+signalname).c_str()  );
  
  
  
  //template for JER up
  std::cout << " - JER up" << std::endl;
  nomJER->cd();
  
  histo_TTbar_JER_up = *(TH1F*)nomJER->Get("MVA_BDT_TTbarSig");
  histo_ZZ_JER_up    = *(TH1F*)nomJER->Get("MVA_BDT_ZZ"   );
  histo_WZ_JER_up    = *(TH1F*)nomJER->Get(MVA_WZname.c_str()   );
  histo_WZ_JER_up.Scale(WZscaleFactor);  
  histo_TZq_JER_up   = *(TH1F*)nomJER->Get("MVA_BDT_TZq");
  histo_TZq_JER_up.Scale(TZqscaleFactor);  
  histo_Sig_JER_up   = *(TH1F*)nomJER->Get(("MVA_BDT_FCNC_"+signalname).c_str()  );
  
  
  //template for JER down
  std::cout << " - JER down" << std::endl;
  nomFile->cd();

  histo_TTbar_JER_down = *(TH1F*)nomFile->Get("MVA_BDT_TTbarSig");
  histo_ZZ_JER_down    = *(TH1F*)nomFile->Get("MVA_BDT_ZZ"   );
  histo_WZ_JER_down    = *(TH1F*)nomFile->Get(MVA_WZname.c_str()   );
  histo_WZ_JER_down.Scale(WZscaleFactor);  
  histo_TZq_JER_down   = *(TH1F*)nomFile->Get("MVA_BDT_TZq");
  histo_TZq_JER_down.Scale(TZqscaleFactor);  
  histo_Sig_JER_down   = *(TH1F*)nomFile->Get(("MVA_BDT_FCNC_"+signalname).c_str()  );
  
  
  //template for BTag up
  std::cout << " - BTag up" << std::endl;
  nomBTag_up->cd();
  
  histo_TTbar_BTag_up = *(TH1F*)nomBTag_up->Get("MVA_BDT_TTbarSig");
  histo_ZZ_BTag_up    = *(TH1F*)nomBTag_up->Get("MVA_BDT_ZZ"   );
  histo_WZ_BTag_up    = *(TH1F*)nomBTag_up->Get(MVA_WZname.c_str()   );
  histo_WZ_BTag_up.Scale(WZscaleFactor);  
  histo_TZq_BTag_up   = *(TH1F*)nomBTag_up->Get("MVA_BDT_TZq");
  histo_TZq_BTag_up.Scale(TZqscaleFactor);  
  histo_Sig_BTag_up   = *(TH1F*)nomBTag_up->Get(("MVA_BDT_FCNC_"+signalname).c_str()  );
  
  
  //template for BTag down
  std::cout << " - BTag down" << std::endl;
  nomBTag_down->cd();
  
  histo_TTbar_BTag_down = *(TH1F*)nomBTag_down->Get("MVA_BDT_TTbarSig");
  histo_ZZ_BTag_down    = *(TH1F*)nomBTag_down->Get("MVA_BDT_ZZ"   );
  histo_WZ_BTag_down    = *(TH1F*)nomBTag_down->Get(MVA_WZname.c_str()   );
  histo_WZ_BTag_down.Scale(WZscaleFactor);  
  histo_TZq_BTag_down   = *(TH1F*)nomBTag_down->Get("MVA_BDT_TZq");
  histo_TZq_BTag_down.Scale(TZqscaleFactor);  
  histo_Sig_BTag_down   = *(TH1F*)nomBTag_down->Get(("MVA_BDT_FCNC_"+signalname).c_str()  );
  
  
  //template for PU up
  std::cout << " - PU up" << std::endl;
  nomPU_up->cd();
  
  histo_TTbar_PU_up = *(TH1F*)nomPU_up->Get("MVA_BDT_TTbarSig");
  histo_ZZ_PU_up    = *(TH1F*)nomPU_up->Get("MVA_BDT_ZZ"   );
  histo_WZ_PU_up    = *(TH1F*)nomPU_up->Get(MVA_WZname.c_str()   );
  histo_WZ_PU_up.Scale(WZscaleFactor);  
  histo_TZq_PU_up   = *(TH1F*)nomPU_up->Get("MVA_BDT_TZq");
  histo_TZq_PU_up.Scale(TZqscaleFactor);  
  histo_Sig_PU_up   = *(TH1F*)nomPU_up->Get(("MVA_BDT_FCNC_"+signalname).c_str()  );
  
  
  //template for PU down
  std::cout << " - PU down" << std::endl;
  nomPU_down->cd();
  
  histo_TTbar_PU_down = *(TH1F*)nomPU_down->Get("MVA_BDT_TTbarSig");
  histo_ZZ_PU_down    = *(TH1F*)nomPU_down->Get("MVA_BDT_ZZ"   );
  histo_WZ_PU_down    = *(TH1F*)nomPU_down->Get(MVA_WZname.c_str()   );
  histo_WZ_PU_down.Scale(WZscaleFactor);  
  histo_TZq_PU_down   = *(TH1F*)nomPU_down->Get("MVA_BDT_TZq");
  histo_TZq_PU_down.Scale(TZqscaleFactor);  
  histo_Sig_PU_down   = *(TH1F*)nomPU_down->Get(("MVA_BDT_FCNC_"+signalname).c_str()  );

  
  //template for Lept up
  std::cout << " - Lept up" << std::endl;
  nomLeptUp->cd();
  
  histo_TTbar_Lept_up = *(TH1F*)nomLeptUp->Get("MVA_BDT_TTbarSig");
  histo_ZZ_Lept_up    = *(TH1F*)nomLeptUp->Get("MVA_BDT_ZZ"   );
  histo_WZ_Lept_up    = *(TH1F*)nomLeptUp->Get(MVA_WZname.c_str()   );
  histo_WZ_Lept_up.Scale(WZscaleFactor);  
  histo_TZq_Lept_up   = *(TH1F*)nomLeptUp->Get("MVA_BDT_TZq");
  histo_TZq_Lept_up.Scale(TZqscaleFactor);  
  histo_Sig_Lept_up   = *(TH1F*)nomLeptUp->Get(("MVA_BDT_FCNC_"+signalname).c_str()  );
  
  
  //template for Lept down
  std::cout << " - Lept down" << std::endl;
  nomLeptDown->cd();
  
  histo_TTbar_Lept_down = *(TH1F*)nomLeptDown->Get("MVA_BDT_TTbarSig");
  histo_ZZ_Lept_down    = *(TH1F*)nomLeptDown->Get("MVA_BDT_ZZ"   );
  histo_WZ_Lept_down    = *(TH1F*)nomLeptDown->Get(MVA_WZname.c_str()   );
  histo_WZ_Lept_down.Scale(WZscaleFactor);  
  histo_TZq_Lept_down    = *(TH1F*)nomLeptDown->Get("MVA_BDT_TZq"   );
  histo_TZq_Lept_down.Scale(TZqscaleFactor);  
  histo_Sig_Lept_down   = *(TH1F*)nomLeptDown->Get(("MVA_BDT_FCNC_"+signalname).c_str()  );

  // template for mtop up
  std::cout << " - TopMass up" << std::endl;
  nomMtop_up->cd();
  histo_Sig_Mtop_up = *(TH1F*)nomMtop_up->Get(("MVA_BDT_FCNC_"+signalname).c_str() );
  
  // template for mtop down
  std::cout << " - TopMass down" << std::endl;
  nomMtop_down->cd();
  histo_Sig_Mtop_down = *(TH1F*)nomMtop_down->Get(("MVA_BDT_FCNC_"+signalname).c_str() );

  // template for scale up
  std::cout << " - Scale up" << std::endl;
  nomScale_up->cd();
  histo_Sig_Scale_up = *(TH1F*)nomScale_up->Get(("MVA_BDT_FCNC_"+signalname).c_str() );
  if (WZ_SYS)
  {
    histo_WZ_Scale_up  = *(TH1F*)nomScale_up->Get(MVA_WZname.c_str());
    histo_WZ_Scale_up.Scale(WZscaleFactor);  
  }

  // template for scale down
  std::cout << " - Scale down" << std::endl;
  nomScale_down->cd();
  histo_Sig_Scale_down = *(TH1F*)nomScale_down->Get(("MVA_BDT_FCNC_"+signalname).c_str() );
  if (WZ_SYS)
  {
    histo_WZ_Scale_down  = *(TH1F*)nomScale_down->Get(MVA_WZname.c_str());
    histo_WZ_Scale_down.Scale(WZscaleFactor);  
  }

  // template for match up
  std::cout << " - Match up" << std::endl;
  nomMatch_up->cd();
  histo_Sig_Match_up = *(TH1F*)nomMatch_up->Get(("MVA_BDT_FCNC_"+signalname).c_str() );
  if (WZ_SYS)
  {
    histo_WZ_Match_up  = *(TH1F*)nomMatch_up->Get(MVA_WZname.c_str());
    histo_WZ_Match_up.Scale(WZscaleFactor);  
  }
  if (TZq_SYS)
  {
    histo_TZq_Match_up  = *(TH1F*)nomMatch_up->Get("MVA_BDT_TZq");
    histo_TZq_Match_up.Scale(TZqscaleFactor);  
  }

  // template for match down
  std::cout << " - Match down" << std::endl;
  nomMatch_down->cd();
  histo_Sig_Match_down = *(TH1F*)nomMatch_down->Get(("MVA_BDT_FCNC_"+signalname).c_str() );
  if (WZ_SYS)
  {
    histo_WZ_Match_down  = *(TH1F*)nomMatch_down->Get(MVA_WZname.c_str());
    histo_WZ_Match_down.Scale(WZscaleFactor);  
  }
  if (TZq_SYS)
  {
    histo_TZq_Match_down  = *(TH1F*)nomMatch_down->Get("MVA_BDT_TZq");
    histo_TZq_Match_down.Scale(TZqscaleFactor);  
  }





  // template for DY up
  std::cout << " - DY up" << std::endl;
  nomDY_up->cd();
  histo_Zjets_DY_up = *(TH1F*)nomDY_up->Get("MVA_BDT_DataZjets" );
  
  // template for DY down
  std::cout << " - DY down" << std::endl;
  nomDY_down->cd();
  histo_Zjets_DY_down = *(TH1F*)nomDY_down->Get("MVA_BDT_DataZjets");


  histo_Zjets_DY_up.Scale(histo_Zjets.Integral()/histo_Zjets_DY_up.Integral());
  histo_Zjets_DY_down.Scale(histo_Zjets.Integral()/histo_Zjets_DY_down.Integral());





  // template for PDF up
  std::cout << " - PDF up" << std::endl;
  nomPDF_up->cd();
  histo_Sig_PDF_up = *(TH1F*)nomPDF_up->Get(("MVA_BDT_FCNC_"+signalname).c_str() );
  histo_WZ_PDF_up  = *(TH1F*)nomPDF_up->Get(MVA_WZname.c_str());
  histo_WZ_PDF_up.Scale(WZscaleFactor);  
  histo_TZq_PDF_up  = *(TH1F*)nomPDF_up->Get("MVA_BDT_TZq");
  histo_TZq_PDF_up.Scale(TZqscaleFactor);  

  // template for PDF down
  std::cout << " - PDF down" << std::endl;
  nomPDF_down->cd();
  histo_Sig_PDF_down = *(TH1F*)nomPDF_down->Get(("MVA_BDT_FCNC_"+signalname).c_str() );
  histo_WZ_PDF_down  = *(TH1F*)nomPDF_down->Get(MVA_WZname.c_str());
  histo_WZ_PDF_down.Scale(WZscaleFactor);  
  histo_TZq_PDF_down  = *(TH1F*)nomPDF_down->Get("MVA_BDT_TZq");
  histo_TZq_PDF_down.Scale(TZqscaleFactor);  
  
  
  //----------------------------------------------------------------------------------------
  //   STEP 3 : Needed by Theta: correct empty bins of the model 
  //            corresponding to non empty bins for the data (if not, Theta crashs)
  //----------------------------------------------------------------------------------------
  std::cout << "Correcting empty bins (required by THETA) ... " << std::endl;
  for (int k = 1; k < histo_Sig.GetNbinsX()+1; k++) 
  {
    //nominal templates 
    if ( histo_Sig.GetBinContent(k)  <=0 ) histo_Sig.SetBinContent(  k,1e-6);
    if ( histo_Data.GetBinContent(k) <=0 ) histo_Data.SetBinContent( k,1e-6);
    if ( histo_Zjets.GetBinContent(k)<=0 ) histo_Zjets.SetBinContent(k,1e-6);
    if ( histo_ZjetsMC.GetBinContent(k)<=0 ) histo_ZjetsMC.SetBinContent(k,1e-6);
    if ( histo_WZ.GetBinContent(k)   <=0 ) histo_WZ.SetBinContent(   k,1e-6);
    if ( histo_TZq.GetBinContent(k)  <=0 ) histo_TZq.SetBinContent(  k,1e-6);
    if ( histo_ZZ.GetBinContent(k)   <=0 ) histo_ZZ.SetBinContent(   k,1e-6);    
    if ( histo_TTbar.GetBinContent(k)<=0 ) histo_TTbar.SetBinContent(k,1e-6);
    
    //templates for JES up
    if ( histo_TTbar_JES_up.GetBinContent(k)  <=0 ) histo_TTbar_JES_up.SetBinContent(k,1e-6);
    if ( histo_ZZ_JES_up.GetBinContent(k)     <=0 ) histo_ZZ_JES_up.SetBinContent(   k,1e-6); 
    if ( histo_WZ_JES_up.GetBinContent(k)     <=0 ) histo_WZ_JES_up.SetBinContent(   k,1e-6); 
    if ( histo_TZq_JES_up.GetBinContent(k)    <=0 ) histo_TZq_JES_up.SetBinContent(  k,1e-6); 
    if ( histo_Sig_JES_up.GetBinContent(k)    <=0 ) histo_Sig_JES_up.SetBinContent(  k,1e-6);
   
    //templates for JES down
    if ( histo_TTbar_JES_down.GetBinContent(k)  <=0 ) histo_TTbar_JES_down.SetBinContent(k,1e-6);
    if ( histo_ZZ_JES_down.GetBinContent(k)     <=0 ) histo_ZZ_JES_down.SetBinContent(   k,1e-6);
    if ( histo_WZ_JES_down.GetBinContent(k)     <=0 ) histo_WZ_JES_down.SetBinContent(   k,1e-6);  
    if ( histo_TZq_JES_down.GetBinContent(k)    <=0 ) histo_TZq_JES_down.SetBinContent(  k,1e-6);  
    if ( histo_Sig_JES_down.GetBinContent(k)    <=0 ) histo_Sig_JES_down.SetBinContent(  k,1e-6);
   
    //templates for JER up
    if ( histo_TTbar_JER_up.GetBinContent(k)  <=0 ) histo_TTbar_JER_up.SetBinContent(k,1e-6);
    if ( histo_ZZ_JER_up.GetBinContent(k)     <=0 ) histo_ZZ_JER_up.SetBinContent(   k,1e-6);
    if ( histo_WZ_JER_up.GetBinContent(k)     <=0 ) histo_WZ_JER_up.SetBinContent(   k,1e-6);
    if ( histo_TZq_JER_up.GetBinContent(k)    <=0 ) histo_TZq_JER_up.SetBinContent(  k,1e-6);
    if ( histo_Sig_JER_up.GetBinContent(k)    <=0 ) histo_Sig_JER_up.SetBinContent(  k,1e-6);
    
    //templates for JER down
    if ( histo_TTbar_JER_down.GetBinContent(k)  <=0 ) histo_TTbar_JER_down.SetBinContent(k,1e-6);
    if ( histo_ZZ_JER_down.GetBinContent(k)     <=0 ) histo_ZZ_JER_down.SetBinContent(   k,1e-6);	
    if ( histo_WZ_JER_down.GetBinContent(k)     <=0 ) histo_WZ_JER_down.SetBinContent(   k,1e-6);	
    if ( histo_TZq_JER_down.GetBinContent(k)    <=0 ) histo_TZq_JER_down.SetBinContent(  k,1e-6);	
    if ( histo_Sig_JER_down.GetBinContent(k)    <=0 ) histo_Sig_JER_down.SetBinContent(  k,1e-6);	
    
    // templates for BTag up
    if ( histo_TTbar_BTag_up.GetBinContent(k)  <=0 ) histo_TTbar_BTag_up.SetBinContent(k,1e-6);
    if ( histo_ZZ_BTag_up.GetBinContent(k)     <=0 ) histo_ZZ_BTag_up.SetBinContent(   k,1e-6); 
    if ( histo_WZ_BTag_up.GetBinContent(k)     <=0 ) histo_WZ_BTag_up.SetBinContent(   k,1e-6); 
    if ( histo_TZq_BTag_up.GetBinContent(k)    <=0 ) histo_TZq_BTag_up.SetBinContent(  k,1e-6); 
    if ( histo_Sig_BTag_up.GetBinContent(k)    <=0 ) histo_Sig_BTag_up.SetBinContent(  k,1e-6);
    
    // templates for BTag down
    if ( histo_TTbar_BTag_down.GetBinContent(k)  <=0 ) histo_TTbar_BTag_down.SetBinContent(k,1e-6);
    if ( histo_ZZ_BTag_down.GetBinContent(k)     <=0 ) histo_ZZ_BTag_down.SetBinContent(   k,1e-6);
    if ( histo_WZ_BTag_down.GetBinContent(k)     <=0 ) histo_WZ_BTag_down.SetBinContent(   k,1e-6);  
    if ( histo_TZq_BTag_down.GetBinContent(k)    <=0 ) histo_TZq_BTag_down.SetBinContent(  k,1e-6);  
    if ( histo_Sig_BTag_down.GetBinContent(k)    <=0 ) histo_Sig_BTag_down.SetBinContent(  k,1e-6);
     
    if ( histo_TTbar_PU_up.GetBinContent(k)  <=0 ) histo_TTbar_PU_up.SetBinContent(k,1e-6);
    if ( histo_ZZ_PU_up.GetBinContent(k)     <=0 ) histo_ZZ_PU_up.SetBinContent(   k,1e-6); 
    if ( histo_WZ_PU_up.GetBinContent(k)     <=0 ) histo_WZ_PU_up.SetBinContent(   k,1e-6); 
    if ( histo_TZq_PU_up.GetBinContent(k)    <=0 ) histo_TZq_PU_up.SetBinContent(  k,1e-6); 
    if ( histo_Sig_PU_up.GetBinContent(k)    <=0 ) histo_Sig_PU_up.SetBinContent(  k,1e-6);
    
    if ( histo_TTbar_PU_down.GetBinContent(k)  <=0 ) histo_TTbar_PU_down.SetBinContent(k,1e-6);
    if ( histo_ZZ_PU_down.GetBinContent(k)     <=0 ) histo_ZZ_PU_down.SetBinContent(   k,1e-6);
    if ( histo_WZ_PU_down.GetBinContent(k)     <=0 ) histo_WZ_PU_down.SetBinContent(   k,1e-6);  
    if ( histo_TZq_PU_down.GetBinContent(k)    <=0 ) histo_TZq_PU_down.SetBinContent(  k,1e-6);  
    if ( histo_Sig_PU_down.GetBinContent(k)    <=0 ) histo_Sig_PU_down.SetBinContent(  k,1e-6);
    
    if ( histo_TTbar_Lept_up.GetBinContent(k)  <=0 ) histo_TTbar_Lept_up.SetBinContent(k,1e-6);
    if ( histo_ZZ_Lept_up.GetBinContent(k)     <=0 ) histo_ZZ_Lept_up.SetBinContent(   k,1e-6); 
    if ( histo_WZ_Lept_up.GetBinContent(k)     <=0 ) histo_WZ_Lept_up.SetBinContent(   k,1e-6); 
    if ( histo_TZq_Lept_up.GetBinContent(k)    <=0 ) histo_TZq_Lept_up.SetBinContent(  k,1e-6); 
    if ( histo_Sig_Lept_up.GetBinContent(k)    <=0 ) histo_Sig_Lept_up.SetBinContent(  k,1e-6);
    
    if ( histo_TTbar_Lept_down.GetBinContent(k)  <=0 ) histo_TTbar_Lept_down.SetBinContent(k,1e-6);
    if ( histo_ZZ_Lept_down.GetBinContent(k)     <=0 ) histo_ZZ_Lept_down.SetBinContent(   k,1e-6);
    if ( histo_WZ_Lept_down.GetBinContent(k)     <=0 ) histo_WZ_Lept_down.SetBinContent(   k,1e-6);  
    if ( histo_TZq_Lept_down.GetBinContent(k)    <=0 ) histo_TZq_Lept_down.SetBinContent(  k,1e-6);  
    if ( histo_Sig_Lept_down.GetBinContent(k)    <=0 ) histo_Sig_Lept_down.SetBinContent(  k,1e-6);

    if (histo_Sig_Mtop_up.GetBinContent(k)   <=0 ) histo_Sig_Mtop_up.SetBinContent(k,1e-6);
    if (histo_Sig_Mtop_down.GetBinContent(k) <=0 ) histo_Sig_Mtop_down.SetBinContent(k,1e-6);

    if (histo_Sig_PDF_up.GetBinContent(k)   <=0 ) histo_Sig_PDF_up.SetBinContent(k,1e-6);
    if (histo_Sig_PDF_down.GetBinContent(k) <=0 ) histo_Sig_PDF_down.SetBinContent(k,1e-6);
    if (histo_WZ_PDF_up.GetBinContent(k)    <=0 ) histo_WZ_PDF_up.SetBinContent(k,1e-6);
    if (histo_WZ_PDF_down.GetBinContent(k)  <=0 ) histo_WZ_PDF_down.SetBinContent(k,1e-6);
    if (histo_TZq_PDF_up.GetBinContent(k)   <=0 ) histo_TZq_PDF_up.SetBinContent(k,1e-6);
    if (histo_TZq_PDF_down.GetBinContent(k) <=0 ) histo_TZq_PDF_down.SetBinContent(k,1e-6);

    if (histo_Sig_Scale_up.GetBinContent(k)   <=0 ) histo_Sig_Scale_up.SetBinContent(k,1e-6);
    if (histo_Sig_Scale_down.GetBinContent(k) <=0 ) histo_Sig_Scale_down.SetBinContent(k,1e-6);
    if (histo_WZ_Scale_up.GetBinContent(k)    <=0 ) histo_WZ_Scale_up.SetBinContent(k,1e-6);
    if (histo_WZ_Scale_down.GetBinContent(k)  <=0 ) histo_WZ_Scale_down.SetBinContent(k,1e-6);

    if (histo_Sig_Match_up.GetBinContent(k)   <=0 ) histo_Sig_Match_up.SetBinContent(k,1e-6);
    if (histo_Sig_Match_down.GetBinContent(k) <=0 ) histo_Sig_Match_down.SetBinContent(k,1e-6);
    if (histo_WZ_Match_up.GetBinContent(k)    <=0 ) histo_WZ_Match_up.SetBinContent(k,1e-6);
    if (histo_WZ_Match_down.GetBinContent(k)  <=0 ) histo_WZ_Match_down.SetBinContent(k,1e-6);
    if (histo_TZq_Match_up.GetBinContent(k)   <=0 ) histo_TZq_Match_up.SetBinContent(k,1e-6);
    if (histo_TZq_Match_down.GetBinContent(k) <=0 ) histo_TZq_Match_down.SetBinContent(k,1e-6);

    if (histo_Zjets_DY_up.GetBinContent(k)   <=0 ) histo_Zjets_DY_up.SetBinContent(k,1e-6);
    if (histo_Zjets_DY_down.GetBinContent(k) <=0 ) histo_Zjets_DY_down.SetBinContent(k,1e-6);
  }

  //----------------------------------------------------------------------------------------
  //   STEP 4 : Setting histogram names
  //----------------------------------------------------------------------------------------
  std::cout << "Setting histogram names ... " << std::endl;
  if(signalname != "tzq") histo_Sig.SetName  (("MVABDT__FCNC_"+signalname+"74").c_str() );
  
  histo_Data.SetName ("MVABDT__DATA"      );
  histo_Zjets.SetName("MVABDT__DataZjets" );
  histo_WZ.SetName   ("MVABDT__WZ"        );
  histo_TZq.SetName  ("MVABDT__TZq"        );
  histo_ZZ.SetName   ("MVABDT__ZZ"        );
  histo_TTbar.SetName("MVABDT__TTbar"     );
  
  if(signalname != "tzq") histo_Sig_JES_up.SetName  (("MVABDT__FCNC_"+signalname+"74__JES__plus").c_str() );
  
  histo_WZ_JES_up.SetName   ("MVABDT__WZ__JES__plus"        );
  histo_TZq_JES_up.SetName   ("MVABDT__TZq__JES__plus"        );
  histo_ZZ_JES_up.SetName   ("MVABDT__ZZ__JES__plus"        );
  histo_TTbar_JES_up.SetName("MVABDT__TTbar__JES__plus"     );
  
  if(signalname != "tzq") histo_Sig_JES_down.SetName  (("MVABDT__FCNC_"+signalname+"74__JES__minus").c_str() );
  histo_WZ_JES_down.SetName   ("MVABDT__WZ__JES__minus"        );
  histo_TZq_JES_down.SetName   ("MVABDT__TZq__JES__minus"        );
  histo_ZZ_JES_down.SetName   ("MVABDT__ZZ__JES__minus"        );
  histo_TTbar_JES_down.SetName("MVABDT__TTbar__JES__minus"     );

  if(signalname != "tzq") histo_Sig_JER_up.SetName  (("MVABDT__FCNC_"+signalname+"74__JER__plus").c_str() );
  histo_WZ_JER_up.SetName   ("MVABDT__WZ__JER__plus"        );
  histo_TZq_JER_up.SetName  ("MVABDT__TZq__JER__plus"        );
  histo_ZZ_JER_up.SetName   ("MVABDT__ZZ__JER__plus"        );
  histo_TTbar_JER_up.SetName("MVABDT__TTbar__JER__plus"     );
  
  if(signalname != "tzq") histo_Sig_JER_down.SetName  (("MVABDT__FCNC_"+signalname+"74__JER__minus").c_str() );
  histo_WZ_JER_down.SetName   ("MVABDT__WZ__JER__minus"        );
  histo_TZq_JER_down.SetName  ("MVABDT__TZq__JER__minus"        );
  histo_ZZ_JER_down.SetName   ("MVABDT__ZZ__JER__minus"        );
  histo_TTbar_JER_down.SetName("MVABDT__TTbar__JER__minus"     );
 
  if(signalname != "tzq") histo_Sig_BTag_up.SetName  (("MVABDT__FCNC_"+signalname+"74__BTag__plus").c_str() );
  histo_WZ_BTag_up.SetName   ("MVABDT__WZ__BTag__plus"        );
  histo_TZq_BTag_up.SetName  ("MVABDT__TZq__BTag__plus"        );
  histo_ZZ_BTag_up.SetName   ("MVABDT__ZZ__BTag__plus"        );
  histo_TTbar_BTag_up.SetName("MVABDT__TTbar__BTag__plus"     );
  
  if(signalname != "tzq") histo_Sig_BTag_down.SetName  (("MVABDT__FCNC_"+signalname+"74__BTag__minus").c_str() );
  histo_WZ_BTag_down.SetName   ("MVABDT__WZ__BTag__minus"        );
  histo_TZq_BTag_down.SetName  ("MVABDT__TZq__BTag__minus"        );
  histo_ZZ_BTag_down.SetName   ("MVABDT__ZZ__BTag__minus"        );
  histo_TTbar_BTag_down.SetName("MVABDT__TTbar__BTag__minus"     );
  
  if(signalname != "tzq") histo_Sig_PU_up.SetName  (("MVABDT__FCNC_"+signalname+"74__PU__plus").c_str() );
  histo_WZ_PU_up.SetName   ("MVABDT__WZ__PU__plus"        );
  histo_TZq_PU_up.SetName  ("MVABDT__TZq__PU__plus"        );
  histo_ZZ_PU_up.SetName   ("MVABDT__ZZ__PU__plus"        );
  histo_TTbar_PU_up.SetName("MVABDT__TTbar__PU__plus"     );
  
  if(signalname != "tzq") histo_Sig_PU_down.SetName  (("MVABDT__FCNC_"+signalname+"74__PU__minus").c_str() );
  histo_WZ_PU_down.SetName   ("MVABDT__WZ__PU__minus"        );
  histo_TZq_PU_down.SetName  ("MVABDT__TZq__PU__minus"       );
  histo_ZZ_PU_down.SetName   ("MVABDT__ZZ__PU__minus"        );
  histo_TTbar_PU_down.SetName("MVABDT__TTbar__PU__minus"     );
  
  if(signalname != "tzq") histo_Sig_Lept_up.SetName  (("MVABDT__FCNC_"+signalname+"74__Lept__plus").c_str() );
  histo_WZ_Lept_up.SetName   ("MVABDT__WZ__Lept__plus"        );
  histo_TZq_Lept_up.SetName   ("MVABDT__TZq__Lept__plus"        );
  histo_ZZ_Lept_up.SetName   ("MVABDT__ZZ__Lept__plus"        );
  histo_TTbar_Lept_up.SetName("MVABDT__TTbar__Lept__plus"     );
  
  if(signalname != "tzq") histo_Sig_Lept_down.SetName  (("MVABDT__FCNC_"+signalname+"74__Lept__minus").c_str() );
  histo_WZ_Lept_down.SetName   ("MVABDT__WZ__Lept__minus"        );
  histo_TZq_Lept_down.SetName  ("MVABDT__TZq__Lept__minus"        );
  histo_ZZ_Lept_down.SetName   ("MVABDT__ZZ__Lept__minus"        );
  histo_TTbar_Lept_down.SetName("MVABDT__TTbar__Lept__minus"     );
  
  if(signalname != "tzq") histo_Sig_Mtop_up.SetName(("MVABDT__FCNC_"+signalname+"74__Mtop__plus").c_str());
  if(signalname != "tzq") histo_Sig_Mtop_down.SetName(("MVABDT__FCNC_"+signalname+"74__Mtop__minus").c_str());

  if(signalname != "tzq") histo_Sig_PDF_up.SetName(("MVABDT__FCNC_"+signalname+"74__PDF__plus").c_str());
  if(signalname != "tzq") histo_Sig_PDF_down.SetName(("MVABDT__FCNC_"+signalname+"74__PDF__minus").c_str());
  histo_WZ_PDF_up.SetName("MVABDT__WZ__PDF__plus");
  histo_WZ_PDF_down.SetName("MVABDT__WZ__PDF__minus");
  histo_TZq_PDF_up.SetName("MVABDT__TZq__PDF__plus");
  histo_TZq_PDF_down.SetName("MVABDT__TZq__PDF__minus");

  if(signalname != "tzq") histo_Sig_Scale_up.SetName(("MVABDT__FCNC_"+signalname+"74__Scale__plus").c_str());
  if(signalname != "tzq") histo_Sig_Scale_down.SetName(("MVABDT__FCNC_"+signalname+"74__Scale__minus").c_str());
  histo_WZ_Scale_up.SetName("MVABDT__WZ__Scale__plus");
  histo_WZ_Scale_down.SetName("MVABDT__WZ__Scale__minus");

  if(signalname != "tzq") histo_Sig_Match_up.SetName(("MVABDT__FCNC_"+signalname+"74__Match__plus").c_str());
  if(signalname != "tzq") histo_Sig_Match_down.SetName(("MVABDT__FCNC_"+signalname+"74__Match__minus").c_str());
  histo_WZ_Match_up.SetName("MVABDT__WZ__Match__plus");
  histo_WZ_Match_down.SetName("MVABDT__WZ__Match__minus");
  histo_TZq_Match_up.SetName("MVABDT__TZq__Match__plus");
  histo_TZq_Match_down.SetName("MVABDT__TZq__Match__minus");
  
  
  
  histo_Zjets_DY_up.SetName("MVABDT__DataZjets__DY__plus");
  histo_Zjets_DY_down.SetName("MVABDT__DataZjets__DY__minus");

  //----------------------------------------------------------------------------------------
  //   STEP 5 : Creating additional histograms related to signal with other xsection values
  //----------------------------------------------------------------------------------------
  std::cout << "Creating additional signal histograms with different xsection values ... " << std::endl;
  double xmin = histo_Sig.GetXaxis()->GetXmin();
  double xmax = histo_Sig.GetXaxis()->GetXmax();
  unsigned int nbins = histo_Sig.GetXaxis()->GetNbins();
  std::vector<double> xsections;
  std::vector<std::string> xnames;

  // ------------------- TO CHANGE BY THE USER
  if(signalname != "kct" ) {xsections.push_back(5.); xnames.push_back("05"); }
  if(signalname != "kct" ) {xsections.push_back(10.); xnames.push_back("10"); }
  xsections.push_back(20.); xnames.push_back("20"); 
  xsections.push_back(30.); xnames.push_back("30"); 
  xsections.push_back(50.); xnames.push_back("50"); 
  if(signalname == "kct" ) {xsections.push_back(70.); xnames.push_back("70"); }
  if(signalname == "kct" ) {xsections.push_back(80.); xnames.push_back("80"); }
  // -------------- TO CHANGE BY THE USER
  
  
  std::vector<TH1F*> histo_newSig;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i];
    histo_newSig.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }

  std::vector<TH1F*> histo_newSig_JES_up;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_JES_up";
    histo_newSig_JES_up.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }
  
  std::vector<TH1F*> histo_newSig_JES_down;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_JES_down";
    histo_newSig_JES_down.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }

  std::vector<TH1F*> histo_newSig_JER_up;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_JER_up";
    histo_newSig_JER_up.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }
  
  std::vector<TH1F*> histo_newSig_JER_down;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_JER_down";
    histo_newSig_JER_down.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }

  std::vector<TH1F*> histo_newSig_BTag_up;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_BTag_up";
    histo_newSig_BTag_up.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }
  
  std::vector<TH1F*> histo_newSig_BTag_down;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_BTag_down";
    histo_newSig_BTag_down.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }

  std::vector<TH1F*> histo_newSig_PU_up;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_PU_up";
    histo_newSig_PU_up.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }
  
  std::vector<TH1F*> histo_newSig_PU_down;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_PU_down";
    histo_newSig_PU_down.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }

  std::vector<TH1F*> histo_newSig_Lept_up;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_Lept_up";
    histo_newSig_Lept_up.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }
  
  std::vector<TH1F*> histo_newSig_Lept_down;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_Lept_down";
    histo_newSig_Lept_down.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }
  
  std::vector<TH1F*> histo_newSig_Mtop_up;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_Mtop_up";
    histo_newSig_Mtop_up.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }
  
  std::vector<TH1F*> histo_newSig_Mtop_down;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_Mtop_down";
    histo_newSig_Mtop_down.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }

  std::vector<TH1F*> histo_newSig_PDF_up;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_PDF_up";
    histo_newSig_PDF_up.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }
  
  std::vector<TH1F*> histo_newSig_PDF_down;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_PDF_down";
    histo_newSig_PDF_down.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }

  std::vector<TH1F*> histo_newSig_Scale_up;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_Scale_up";
    histo_newSig_Scale_up.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }
  
  std::vector<TH1F*> histo_newSig_Scale_down;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_Scale_down";
    histo_newSig_Scale_down.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }

  std::vector<TH1F*> histo_newSig_Match_up;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_Match_up";
    histo_newSig_Match_up.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }
  
  std::vector<TH1F*> histo_newSig_Match_down;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    std::string histoname = "histo_Sig_"+xnames[i]+"_Match_down";
    histo_newSig_Match_down.push_back(new TH1F(histoname.c_str(), histoname.c_str(), nbins, xmin, xmax) );
  }


  //----------------------------------------------------------------------------------------
  //   STEP 6 : Scaling additional signal histograms
  //----------------------------------------------------------------------------------------
  std::cout << "Scaling additional signal histograms ... " << std::endl;
  for (int k = 1; k < histo_Sig.GetNbinsX()+1; k++) 
    for (unsigned int i=0;i<xsections.size();i++)
    {
      /*histo_newSig[i]           -> SetBinContent(k,xsections[i]*histo_Sig.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_JES_up[i]    -> SetBinContent(k,xsections[i]*histo_Sig_JES_up.GetBinContent(k)/histo_Sig_JES_up.Integral());
      histo_newSig_JES_down[i]  -> SetBinContent(k,xsections[i]*histo_Sig_JES_down.GetBinContent(k)/histo_Sig_JES_down.Integral());
      histo_newSig_JER_up[i]    -> SetBinContent(k,xsections[i]*histo_Sig_JER_up.GetBinContent(k)/histo_Sig_JER_up.Integral());
      histo_newSig_JER_down[i]  -> SetBinContent(k,xsections[i]*histo_Sig_JER_down.GetBinContent(k)/histo_Sig_JER_down.Integral());
      histo_newSig_BTag_up[i]   -> SetBinContent(k,xsections[i]*histo_Sig_BTag_up.GetBinContent(k)/histo_Sig_BTag_up.Integral());
      histo_newSig_BTag_down[i] -> SetBinContent(k,xsections[i]*histo_Sig_BTag_down.GetBinContent(k)/histo_Sig_BTag_down.Integral());
      histo_newSig_PU_up[i]     -> SetBinContent(k,xsections[i]*histo_Sig_PU_up.GetBinContent(k)/histo_Sig_PU_up.Integral());
      histo_newSig_PU_down[i]   -> SetBinContent(k,xsections[i]*histo_Sig_PU_down.GetBinContent(k)/histo_Sig_PU_down.Integral());
      histo_newSig_Lept_up[i]   -> SetBinContent(k,xsections[i]*histo_Sig_Lept_up.GetBinContent(k)/histo_Sig_Lept_up.Integral());
      histo_newSig_Lept_down[i] -> SetBinContent(k,xsections[i]*histo_Sig_Lept_down.GetBinContent(k)/histo_Sig_Lept_down.Integral());
      histo_newSig_Mtop_up[i]   -> SetBinContent(k,xsections[i]*histo_Sig_Mtop_up.GetBinContent(k)/histo_Sig_Mtop_up.Integral());
      histo_newSig_Mtop_down[i] -> SetBinContent(k,xsections[i]*histo_Sig_Mtop_down.GetBinContent(k)/histo_Sig_Mtop_down.Integral());
      histo_newSig_PDF_up[i]    -> SetBinContent(k,xsections[i]*histo_Sig_PDF_up.GetBinContent(k)/histo_Sig_PDF_up.Integral());
      histo_newSig_PDF_down[i]  -> SetBinContent(k,xsections[i]*histo_Sig_PDF_down.GetBinContent(k)/histo_Sig_PDF_down.Integral());
      histo_newSig_Scale_up[i]  -> SetBinContent(k,xsections[i]*histo_Sig_Scale_up.GetBinContent(k)/histo_Sig_Scale_up.Integral());
      histo_newSig_Scale_down[i]-> SetBinContent(k,xsections[i]*histo_Sig_Scale_down.GetBinContent(k)/histo_Sig_Scale_down.Integral());
      histo_newSig_Match_up[i]  -> SetBinContent(k,xsections[i]*histo_Sig_Match_up.GetBinContent(k)/histo_Sig_Match_up.Integral());
      histo_newSig_Match_down[i]-> SetBinContent(k,xsections[i]*histo_Sig_Match_down.GetBinContent(k)/histo_Sig_Match_down.Integral());
      */
      
      
      histo_newSig[i]           -> SetBinContent(k,xsections[i]*histo_Sig.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_JES_up[i]    -> SetBinContent(k,xsections[i]*histo_Sig_JES_up.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_JES_down[i]  -> SetBinContent(k,xsections[i]*histo_Sig_JES_down.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_JER_up[i]    -> SetBinContent(k,xsections[i]*histo_Sig_JER_up.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_JER_down[i]  -> SetBinContent(k,xsections[i]*histo_Sig_JER_down.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_BTag_up[i]   -> SetBinContent(k,xsections[i]*histo_Sig_BTag_up.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_BTag_down[i] -> SetBinContent(k,xsections[i]*histo_Sig_BTag_down.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_PU_up[i]     -> SetBinContent(k,xsections[i]*histo_Sig_PU_up.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_PU_down[i]   -> SetBinContent(k,xsections[i]*histo_Sig_PU_down.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_Lept_up[i]   -> SetBinContent(k,xsections[i]*histo_Sig_Lept_up.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_Lept_down[i] -> SetBinContent(k,xsections[i]*histo_Sig_Lept_down.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_Mtop_up[i]   -> SetBinContent(k,xsections[i]*histo_Sig_Mtop_up.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_Mtop_down[i] -> SetBinContent(k,xsections[i]*histo_Sig_Mtop_down.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_PDF_up[i]    -> SetBinContent(k,xsections[i]*histo_Sig_PDF_up.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_PDF_down[i]  -> SetBinContent(k,xsections[i]*histo_Sig_PDF_down.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_Scale_up[i]  -> SetBinContent(k,xsections[i]*histo_Sig_Scale_up.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_Scale_down[i]-> SetBinContent(k,xsections[i]*histo_Sig_Scale_down.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_Match_up[i]  -> SetBinContent(k,xsections[i]*histo_Sig_Match_up.GetBinContent(k)/histo_Sig.Integral());
      histo_newSig_Match_down[i]-> SetBinContent(k,xsections[i]*histo_Sig_Match_down.GetBinContent(k)/histo_Sig.Integral());
      
      
      
      
      
  }

  //----------------------------------------------------------------------------------------
  //   STEP 7 : Setting name to additional histograms 
  //----------------------------------------------------------------------------------------
  std::cout << "Setting name to additional signal histograms ... " << std::endl;
  for (unsigned int i=0;i<xsections.size();i++)
  {
    if(signalname != "tzq")  {
      histo_newSig[i]           -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]).c_str());
      histo_newSig_JES_up[i]    -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__JES__plus"   ).c_str());
      histo_newSig_JES_down[i]  -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__JES__minus"  ).c_str());
      histo_newSig_JER_up[i]    -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__JER__plus"   ).c_str());
      histo_newSig_JER_down[i]  -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__JER__minus"  ).c_str());
      histo_newSig_BTag_up[i]   -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__BTag__plus"  ).c_str());
      histo_newSig_BTag_down[i] -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__BTag__minus" ).c_str());
      histo_newSig_PU_up[i]     -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__PU__plus"    ).c_str());
      histo_newSig_PU_down[i]   -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__PU__minus"   ).c_str());
      histo_newSig_Lept_up[i]   -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__Lept__plus"  ).c_str());
      histo_newSig_Lept_down[i] -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__Lept__minus" ).c_str());
      histo_newSig_Mtop_up[i]   -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__Mtop__plus"  ).c_str());
      histo_newSig_Mtop_down[i] -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__Mtop__minus" ).c_str());
      histo_newSig_PDF_up[i]    -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__PDF__plus"   ).c_str());
      histo_newSig_PDF_down[i]  -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__PDF__minus"  ).c_str());
      histo_newSig_Scale_up[i]  -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__Scale__plus" ).c_str());
      histo_newSig_Scale_down[i]-> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__Scale__minus").c_str());
      histo_newSig_Match_up[i]  -> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__Match__plus" ).c_str());
      histo_newSig_Match_down[i]-> SetName(("MVABDT__FCNC_"+signalname+xnames[i]+"__Match__minus").c_str());
    }
    
    
    
    
  }
  



  //----------------------------------------------------------------------------------------
  //   STEP 7.5 : Rebinning histogram names
  //----------------------------------------------------------------------------------------
  std::cout << "Rebinning histograms ... " << std::endl;
  unsigned int rebin = 1;
  histo_Sig.Rebin(rebin);
  histo_Data.Rebin(rebin);
  histo_Zjets.Rebin(rebin);
  histo_WZ.Rebin(rebin);
  histo_TZq.Rebin(rebin);
  histo_ZZ.Rebin(rebin);
  histo_TTbar.Rebin(rebin);
  
  histo_Sig_JES_up.Rebin(rebin);
  histo_WZ_JES_up.Rebin(rebin);
  histo_TZq_JES_up.Rebin(rebin);
  histo_ZZ_JES_up.Rebin(rebin);
  histo_TTbar_JES_up.Rebin(rebin);
  
  histo_Sig_JES_down.Rebin(rebin);
  histo_WZ_JES_down.Rebin(rebin);
  histo_TZq_JES_down.Rebin(rebin);
  histo_ZZ_JES_down.Rebin(rebin);
  histo_TTbar_JES_down.Rebin(rebin);

  histo_Sig_JER_up.Rebin(rebin);
  histo_WZ_JER_up.Rebin(rebin);
  histo_TZq_JER_up.Rebin(rebin);
  histo_ZZ_JER_up.Rebin(rebin);
  histo_TTbar_JER_up.Rebin(rebin);
  
  histo_Sig_JER_down.Rebin(rebin);
  histo_WZ_JER_down.Rebin(rebin);
  histo_TZq_JER_down.Rebin(rebin);
  histo_ZZ_JER_down.Rebin(rebin);
  histo_TTbar_JER_down.Rebin(rebin);
 
  histo_Sig_BTag_up.Rebin(rebin);
  histo_WZ_BTag_up.Rebin(rebin);
  histo_TZq_BTag_up.Rebin(rebin);
  histo_ZZ_BTag_up.Rebin(rebin);
  histo_TTbar_BTag_up.Rebin(rebin);
  
  histo_Sig_BTag_down.Rebin(rebin);
  histo_WZ_BTag_down.Rebin(rebin);
  histo_TZq_BTag_down.Rebin(rebin);
  histo_ZZ_BTag_down.Rebin(rebin);
  histo_TTbar_BTag_down.Rebin(rebin);
  
  histo_Sig_PU_up.Rebin(rebin);
  histo_WZ_PU_up.Rebin(rebin);
  histo_TZq_PU_up.Rebin(rebin);
  histo_ZZ_PU_up.Rebin(rebin);
  histo_TTbar_PU_up.Rebin(rebin);
  
  histo_Sig_PU_down.Rebin(rebin);
  histo_WZ_PU_down.Rebin(rebin);
  histo_TZq_PU_down.Rebin(rebin);
  histo_ZZ_PU_down.Rebin(rebin);
  histo_TTbar_PU_down.Rebin(rebin);
  
  histo_Sig_Lept_up.Rebin(rebin);
  histo_WZ_Lept_up.Rebin(rebin);
  histo_TZq_Lept_up.Rebin(rebin);
  histo_ZZ_Lept_up.Rebin(rebin);
  histo_TTbar_Lept_up.Rebin(rebin);
  
  histo_Sig_Lept_down.Rebin(rebin);
  histo_WZ_Lept_down.Rebin(rebin);
  histo_TZq_Lept_down.Rebin(rebin);
  histo_ZZ_Lept_down.Rebin(rebin);
  histo_TTbar_Lept_down.Rebin(rebin);
  
  histo_Sig_Mtop_up.Rebin(rebin);
  histo_Sig_Mtop_down.Rebin(rebin);

  histo_Sig_PDF_up.Rebin(rebin);
  histo_Sig_PDF_down.Rebin(rebin);
  histo_WZ_PDF_up.Rebin(rebin);
  histo_WZ_PDF_down.Rebin(rebin);

  histo_Sig_Scale_up.Rebin(rebin);
  histo_Sig_Scale_down.Rebin(rebin);
  histo_WZ_Scale_up.Rebin(rebin);
  histo_WZ_Scale_down.Rebin(rebin);

  histo_Sig_Match_up.Rebin(rebin);
  histo_Sig_Match_down.Rebin(rebin);
  histo_WZ_Match_up.Rebin(rebin);
  histo_WZ_Match_down.Rebin(rebin);
  histo_TZq_Match_up.Rebin(rebin);
  histo_TZq_Match_down.Rebin(rebin);
  
  histo_Zjets_DY_up.Rebin(rebin);
  histo_Zjets_DY_down.Rebin(rebin);

  
  

  for (unsigned int i=0;i<xsections.size();i++)
  {
    histo_newSig[i]           -> Rebin(rebin); 
    histo_newSig_JES_up[i]    -> Rebin(rebin); 
    histo_newSig_JES_down[i]  -> Rebin(rebin); 
    histo_newSig_JER_up[i]    -> Rebin(rebin); 
    histo_newSig_JER_down[i]  -> Rebin(rebin); 
    histo_newSig_BTag_up[i]   -> Rebin(rebin); 
    histo_newSig_BTag_down[i] -> Rebin(rebin); 
    histo_newSig_PU_up[i]     -> Rebin(rebin); 
    histo_newSig_PU_down[i]   -> Rebin(rebin); 
    histo_newSig_Lept_up[i]   -> Rebin(rebin); 
    histo_newSig_Lept_down[i] -> Rebin(rebin); 
    histo_newSig_Mtop_up[i]   -> Rebin(rebin); 
    histo_newSig_Mtop_down[i] -> Rebin(rebin); 
    histo_newSig_PDF_up[i]    -> Rebin(rebin); 
    histo_newSig_PDF_down[i]  -> Rebin(rebin); 
    histo_newSig_Scale_up[i]  -> Rebin(rebin); 
    histo_newSig_Scale_down[i]-> Rebin(rebin); 
    histo_newSig_Match_up[i]  -> Rebin(rebin); 
    histo_newSig_Match_down[i]-> Rebin(rebin); 
  }


  //----------------------------------------------------------------------------------------
  //   STEP 8 : Saving histograms
  //----------------------------------------------------------------------------------------
  std::cout << "Saving histograms ... " << std::endl;
  //TFile * outputfile = new TFile("NewFileToBeUsedForThetaWithAutoNamingConvention.root","new");
  TFile * outputfile = new TFile(outputname.c_str(),"RECREATE");
  outputfile->cd();
  
  /*
    // add histo for making bands
    histo_Sig_05->Rebin(4);
    histo_Sig_10->Rebin(4);
    histo_Sig_20->Rebin(4);
    histo_Sig_30->Rebin(4);
    histo_Sig_50->Rebin(4);
    histo_Sig.Rebin(4);
  
    histo_Data.Rebin(4);
    histo_Zjets.Rebin(4);
    histo_WZ.Rebin(4);
    histo_ZZ.Rebin(4);
    histo_TTbar.Rebin(4);
  */
  
  // add histo for making bands
  std::cout << " - nominal" << std::endl;
  for (unsigned int i=0;i<xsections.size();i++) histo_newSig[i]->Write();
  histo_Sig.Write();
  histo_Data.Write();
  histo_Zjets.Write();
  histo_WZ.Write();
  histo_TZq.Write();
  histo_ZZ.Write();
  histo_TTbar.Write();

  if (doJES)
  {
    std::cout << " - JES up" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_JES_up[i]->Write();
    histo_Sig_JES_up.Write();
    histo_WZ_JES_up.Write();
    histo_TZq_JES_up.Write();
    histo_ZZ_JES_up.Write();
    histo_TTbar_JES_up.Write();
  }

  if (doJES)
  {
    std::cout << " - JES down" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_JES_down[i]->Write();
    histo_Sig_JES_down.Write();
    histo_WZ_JES_down.Write();
    histo_TZq_JES_down.Write();
    histo_ZZ_JES_down.Write();
    histo_TTbar_JES_down.Write();
  }

  if (doJER)
  {
    std::cout << " - JER up" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_JER_up[i]->Write();
    histo_Sig_JER_up.Write();
    histo_WZ_JER_up.Write();
    histo_TZq_JER_up.Write();
    histo_ZZ_JER_up.Write();
    histo_TTbar_JER_up.Write();
  }

  if (doJER)
  {
    std::cout << " - JER down" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_JER_down[i]->Write();
    histo_Sig_JER_down.Write();
    histo_WZ_JER_down.Write();
    histo_TZq_JER_down.Write();
    histo_ZZ_JER_down.Write();
    histo_TTbar_JER_down.Write();
  }

  if (doBTag)
  {
    std::cout << " - BTag up" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_BTag_up[i]->Write();
    histo_Sig_BTag_up.Write();
    histo_WZ_BTag_up.Write();
    histo_TZq_BTag_up.Write();
    histo_ZZ_BTag_up.Write();
    histo_TTbar_BTag_up.Write();
  }

  if (doBTag)
  {
    std::cout << " - BTag down" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_BTag_down[i]->Write();
    histo_Sig_BTag_down.Write();
    histo_WZ_BTag_down.Write();
    histo_TZq_BTag_down.Write();
    histo_ZZ_BTag_down.Write();
    histo_TTbar_BTag_down.Write();
  }

  if (doPU)
  {
    std::cout << " - PU up" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_PU_up[i]->Write();
    histo_Sig_PU_up.Write();
    histo_WZ_PU_up.Write();
    histo_TZq_PU_up.Write();
    histo_ZZ_PU_up.Write();
    histo_TTbar_PU_up.Write();
  }

  if (doPU)
  {
    std::cout << " - PU down" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_PU_down[i]->Write();
    histo_Sig_PU_down.Write();
    histo_WZ_PU_down.Write();
    histo_TZq_PU_down.Write();
    histo_ZZ_PU_down.Write();
    histo_TTbar_PU_down.Write();
  }

  if (doLept)
  {
    std::cout << " - Lept up" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_Lept_up[i]->Write();
    histo_Sig_Lept_up.Write();
    histo_WZ_Lept_up.Write();
    histo_TZq_Lept_up.Write();
    histo_ZZ_Lept_up.Write();
    histo_TTbar_Lept_up.Write();
  }

  if (doLept)
  {
    std::cout << " - Lept down" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_Lept_down[i]->Write();
    histo_Sig_Lept_down.Write();
    histo_WZ_Lept_down.Write();
    histo_TZq_Lept_down.Write();
    histo_ZZ_Lept_down.Write();
    histo_TTbar_Lept_down.Write();
  }

  if (doTopMass)
  {
    std::cout << " - TopMass up" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_Mtop_up[i]->Write();
    histo_Sig_Mtop_up.Write();
  }

  if (doTopMass)
  {
    std::cout << " - TopMass down" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_Mtop_down[i]->Write();
    histo_Sig_Mtop_down.Write();
  }

  if (doPDF)
  {
    std::cout << " - PDF up" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_PDF_up[i]->Write();
    histo_Sig_PDF_up.Write();
    histo_WZ_PDF_up.Write();
    histo_TZq_PDF_up.Write();
  }

  if (doPDF)
  {
    std::cout << " - PDF down" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_PDF_down[i]->Write();
    histo_Sig_PDF_down.Write();
    histo_WZ_PDF_down.Write();
    histo_TZq_PDF_down.Write();
  }

  if (doScale)
  {
    std::cout << " - Scale up" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_Scale_up[i]->Write();
    if(!nomatching_scale_singal) histo_Sig_Scale_up.Write();
    if (WZ_SYS) histo_WZ_Scale_up.Write();
  }

  if (doScale)
  {
    std::cout << " - Scale down" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_Scale_down[i]->Write();
    if(!nomatching_scale_singal)histo_Sig_Scale_down.Write();
    if (WZ_SYS) histo_WZ_Scale_down.Write();
  }

  if (doMatch)
  {
    std::cout << " - Match up" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_Match_up[i]->Write();
    if(!nomatching_scale_singal)histo_Sig_Match_up.Write();
    if (WZ_SYS)  histo_WZ_Match_up.Write();
    if (TZq_SYS) histo_TZq_Match_up.Write();
  }

  if (doMatch)
  {
    std::cout << " - Match down" << std::endl;
    for (unsigned int i=0;i<xsections.size();i++) histo_newSig_Match_down[i]->Write();
    if(!nomatching_scale_singal)histo_Sig_Match_down.Write();
    if (WZ_SYS)  histo_WZ_Match_down.Write();
    if (TZq_SYS) histo_TZq_Match_down.Write();
  }

  if (doDY)
  {
    std::cout << " - DY up" << std::endl;
    histo_Zjets_DY_up.Write();
    std::cout << " - DY down" << std::endl;
    histo_Zjets_DY_down.Write();
  }
  
  
  std::cout << "Closing output file ... " << std::endl;
  outputfile->Close();  

  std::cout << "End!" << std::endl;  
} 
