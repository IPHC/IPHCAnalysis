#ifndef JetHistoManager_h
#define JetHistoManager_h

#include "NTFormat/interface/NTJet.h"

#include "Plots/interface/HistoManager.h"


using namespace IPHCTree;

class JetHistoManager: public HistoManager{

  public:
	JetHistoManager();
	~JetHistoManager();

	//Initialisation methods

	void CreateHistos(); /** Create a bunch of standard histos */

	//Fill methods
	void Fill(const vector<NTJet>& jets, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight);
	void FillSelStep(const vector<NTJet>& jets, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight);


};

#endif
