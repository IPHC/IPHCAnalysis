#include "Tools/interface/Dataset.h"
#include "Tools/interface/AnalysisEnvironmentLoader.h"
#include "BckgdEstimation/interface/MMEstimation.h"

void FastAnalysis(string datatype="MC", string channel="MuMu") // MC or DATA, MuMu or EE
{
  cout<<"#########################"<<endl;
  cout<<"Beginning of the program"<<endl;
  cout<<"#########################"<<endl;
  
  //////////////////////
  // Initialisation
  //////////////////////
  vector < Dataset > datasets;
  int verbosity = -1;

  ////////////////////////////////
  // Matrix Method declaration
  ///////////////////////////////
  float looseIso = 0.8 ; float tightIsoMMmu = 0.2; float tightIsoMMe = 0.17;
  float tightIso = tightIsoMMmu;
  if(channel=="EE") tightIso = tightIsoMMe;
  unsigned int nbinsMM = 11; float lowEdgeMM = -0.5; float highEdgeMM = 10.5;
  MMEstimation MMest(datasets,looseIso, tightIso, nbinsMM, lowEdgeMM, highEdgeMM, channel); 
  //////////////////////

  //////////////////////////////////////////
  // Running the Matrix Method in fast way
  //////////////////////////////////////////

  MMest.ReadMMFile(("MatrixMethod_OutPut_"+channel+"_"+datatype+".root").c_str());

  TH1F * hSignalEfficiencyNJets;
  TH1F * hFakeRateNJets;
  
  string onelept = channel;
  onelept.replace(0, channel.length()/2, "");
  TFile* file = new TFile(("MatrixMethod_Efficiency_"+datatype+"_"+onelept+".root").c_str());
  file->cd();
  hSignalEfficiencyNJets = (TH1F*)gDirectory->Get("SignalEfficiencyNJets");
  hFakeRateNJets = (TH1F*)gDirectory->Get("FakeRateNJets");
  
  vector<struct MMEpsilons> valMMEpsilons;
  for(unsigned int bin_index = 0; bin_index < nbinsMM; bin_index++){
    struct MMEpsilons valMMEpsilonsTmp;
    valMMEpsilonsTmp.EpsilonSignal = hSignalEfficiencyNJets->GetBinContent(bin_index+1);
    valMMEpsilonsTmp.EpsilonSignalErr = hSignalEfficiencyNJets->GetBinError(bin_index+1);
    valMMEpsilonsTmp.EpsilonFake = hFakeRateNJets->GetBinContent(bin_index+1);
    valMMEpsilonsTmp.EpsilonFakeErr = hFakeRateNJets->GetBinError(bin_index+1);
    valMMEpsilons.push_back(valMMEpsilonsTmp);
  }

  file->Close();
  delete file;

  unsigned int NbIterations = 100000;
  bool doStatistical = true; bool doSystematic = true; bool doCorrections = false;
  MMest.RunTheMatrixMethod(valMMEpsilons,  NbIterations, doStatistical, doSystematic, doCorrections);
  MMest.PrintMMEstimated();
  MMest.WriteMMFileFast(("MatrixMethod_OutPut_"+channel+"_"+datatype+"_Fast.root").c_str());

  cout<<"#########################"<<endl;
  cout<<"    End of the program   "<<endl;
  cout<<"#########################"<<endl;

  return;

}

int main ()
{

  //FastAnalysis("MC", "EE");
  FastAnalysis("MC", "MuMu");

  return (0);
}
