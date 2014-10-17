void SETUP()
{
// Load the libraries
gSystem->Load(TString(gSystem->ExpandPathName("$NTUPLEDATAFORMAT_PATH")+TString("/src/libNTuple.so")));
gSystem->Load(TString(gSystem->ExpandPathName("$NTUPLEANA_PATH")+TString("/../.lib/libNTupleAna.so")));
// Set the include paths
/*gROOT->ProcessLine(TString(".include ")+TString(gSystem->ExpandPathName("$NTUPLEANA_PATH")+TString("/../BckgdEstimation/interface")));
gROOT->ProcessLine(TString(".include ")+TString(gSystem->ExpandPathName("$NTUPLEANA_PATH")+TString("/../EffEstimation/interface")));
gROOT->ProcessLine(TString(".include ")+TString(gSystem->ExpandPathName("$NTUPLEANA_PATH")+TString("/../Measurements/interface")));
gROOT->ProcessLine(TString(".include ")+TString(gSystem->ExpandPathName("$NTUPLEANA_PATH")+TString("/../Plots/interface")));
gROOT->ProcessLine(TString(".include ")+TString(gSystem->ExpandPathName("$NTUPLEANA_PATH")+TString("/../Selection/interface")));
gROOT->ProcessLine(TString(".include ")+TString(gSystem->ExpandPathName("$NTUPLEANA_PATH")+TString("/../Tools/interface")));
gROOT->ProcessLine(TString(".include ")+TString(gSystem->ExpandPathName("$NTUPLEANA_PATH")+TString("/../tinyxml")));
gROOT->ProcessLine(TString(".include ")+TString(gSystem->ExpandPathName("$NTUPLEDATAFORMAT_PATH")+TString("/../interface")));*/
gROOT->ProcessLine(TString(".include ")+TString(gSystem->ExpandPathName("$NTUPLEDATAFORMAT_PATH")+TString("/../interface")));
gROOT->ProcessLine(TString(".include ")+TString(gSystem->ExpandPathName("$NTUPLEDATAFORMAT_PATH")+TString("/../")));
gROOT->ProcessLine(TString(".include ")+TString(gSystem->ExpandPathName("$NTUPLEANA_PATH")+TString("/../")));
}

