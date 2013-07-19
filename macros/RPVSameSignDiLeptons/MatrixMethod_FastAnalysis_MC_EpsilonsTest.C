#include "Tools/interface/Dataset.h"
#include "Tools/interface/AnalysisEnvironmentLoader.h"
#include "BckgdEstimation/interface/MMEstimation.h"

int main ()
{
  cout<<"#########################"<<endl;
  cout<<"Beginning of the program"<<endl;
  cout<<"#########################"<<endl;
  
  //////////////////////
  // Initialisation
  //////////////////////
  float Luminosity = 0;
  float LumiError = 0;
  string PUWeightFileName;
  int DataType =0;

  string xmlFileName = string ("../../config/MatrixMethod_MC.xml");
  AnalysisEnvironmentLoader anaEl (xmlFileName);
  vector < Dataset > datasets;
  anaEl.LoadSamples(datasets);
  int verbosity = -1;
  anaEl.LoadGeneralInfo(DataType, Luminosity, LumiError, PUWeightFileName, verbosity);

  ////////////////////////////////
  // Matrix Method declaration
  ///////////////////////////////
  float looseIsoMM = 0.8; float tightIsoMMe = 0.17; float tightIsoMMmu = 0.2;
  unsigned int nbinsMM = 11; float lowEdgeMM = -0.5; float highEdgeMM = 10.5;
  MMEstimation MMestEE(datasets, looseIsoMM, tightIsoMMe, nbinsMM, lowEdgeMM, highEdgeMM, "EE");
  MMEstimation MMestMuMu(datasets, looseIsoMM, tightIsoMMe, nbinsMM, lowEdgeMM, highEdgeMM, "MuMu");

  //////////////////////

  //////////////////////////////////////////
  // Running the Matrix Method in fast way
  //////////////////////////////////////////
  //MMestEE.ReadMMFileForPullTest("MatrixMethod_OutPut_EE_MC_forNoWeight.root", 0.95, 0.2);
  MMestMuMu.ReadMMFileForPullTest("MatrixMethod_OutPut_MuMu_MC_NoWeight.root", 0.95, 0.2);

  /*vector<struct MMEpsilons> valMMEpsilons_EE;
  for(unsigned int bin_index = 0; bin_index < nbinsMM; bin_index++){
    struct MMEpsilons valMMEpsilonsTmp_EE;
    valMMEpsilonsTmp_EE.EpsilonSignal = 0.95;
    valMMEpsilonsTmp_EE.EpsilonSignalErr = 0.1;
    valMMEpsilonsTmp_EE.EpsilonFake = 0.2;
    valMMEpsilonsTmp_EE.EpsilonFakeErr = 0.1;
    valMMEpsilons_EE.push_back(valMMEpsilonsTmp_EE);
  }*/

  vector<struct MMEpsilons> valMMEpsilons_MuMu;
  for(unsigned int bin_index = 0; bin_index < nbinsMM; bin_index++){
    struct MMEpsilons valMMEpsilonsTmp_MuMu;
    valMMEpsilonsTmp_MuMu.EpsilonSignal = 0.95;
    valMMEpsilonsTmp_MuMu.EpsilonSignalErr = 0.1;
    valMMEpsilonsTmp_MuMu.EpsilonFake = 0.2;
    valMMEpsilonsTmp_MuMu.EpsilonFakeErr = 0.1;
    valMMEpsilons_MuMu.push_back(valMMEpsilonsTmp_MuMu);
  }


  float weight = datasets[0].NormFactor()*Luminosity;

  unsigned int NbIterations = 100000;
  bool doStatistical = true; bool doSystematic = false;
  
  /*MMestEE.RunTheMatrixMethodForEpsilonsTest(valMMEpsilons_EE, NbIterations, doStatistical, doSystematic, weight);
  MMestEE.PrintMMEstimated();
  MMestEE.WriteMMFileFastForPullTest("MatrixMethod_OutPut_EE_MC_forNoWeight.root","MatrixMethod_OutPut_EE_MC_Fast_EpsilonsTest.root", 0.95, 0.2);*/


  MMestMuMu.RunTheMatrixMethodForEpsilonsTest(valMMEpsilons_MuMu, NbIterations, doStatistical, doSystematic, weight);
  MMestMuMu.PrintMMEstimated();
  MMestMuMu.WriteMMFileFastForPullTest("MatrixMethod_OutPut_MuMu_MC_NoWeight.root","MatrixMethod_OutPut_MuMu_MC_Fast_EpsilonsTest.root", 0.95, 0.2);


  cout<<"#########################"<<endl;
  cout<<"    End of the program   "<<endl;
  cout<<"#########################"<<endl;

  return (0);
}
