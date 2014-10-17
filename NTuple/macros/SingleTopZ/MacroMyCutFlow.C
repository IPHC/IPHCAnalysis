#include <TSelector.h>
#include <TChain.h>
#include <vector>

#include "../../Tools/interface/Dataset.h"
#include "../../Tools/interface/AnalysisEnvironmentLoader.h"
#include "TSystem.h"
#include "SETUP.C"

int main(int argc, char* argv[]){

  SETUP(); // load libraries and files paths

  string macroName = "ProofSelectorMyCutFlow.C+"; //"+" should be put at the end to use ACLIC complication - This macro should inherit from TSelector 
  //In order to allow the node to access the xml, the name should be given with the full path
  string xmlFileName = getenv( "CMSSW_BASE" )+string("/src/IPHCAnalysis/NTuple/config/MyCutFlow_2012_FCNCkut.xml");
  string outputFileName = "proof.root";
    
  //---------------------------------------//
  // Xml Loading  & Dataset registration
  //---------------------------------------//
  
  vector < Dataset > datasets;
  AnalysisEnvironmentLoader anaEL (xmlFileName);
  anaEL.LoadSamples (datasets); // now the list of datasets written in the xml file is known
    
  //---------------------------------------//
  // 	Processing of the datasets
  //---------------------------------------//
  
  string outputFileNameModif = outputFileName.substr(0,outputFileName.size()-5);
  
  for(unsigned int i=0;i<datasets.size();i++){
    
    cout<<"################################################################"<<endl;
    cout<<"########### Processing the dataset "<<datasets[i].Name()<<endl;
    cout<<"################################################################"<<endl;

    TChain *tree = new TChain("tree");
    for(unsigned int j=0;j<datasets[i].Filenames().size();j++){
      string treePath = datasets[i].Filenames()[j]+"/MyModule/Event";
      tree->Add(treePath.c_str());
    }
    
    TSelector *selector = TSelector::GetSelector(macroName.c_str());
    TList *inputs = new TList();
    inputs->Add(new TNamed("PROOF_DATASETNAME", datasets[i].Name()));
    inputs->Add(new TNamed("PROOF_XMLFILENAME", xmlFileName));
    string outputnameSample = "proof_"+datasets[i].Name();
    inputs->Add(new TNamed("PROOF_OUTPUTFILE", outputnameSample));
    selector->SetInputList( inputs );
    tree->Process(selector);
    
  }
  
  cout << "start backuping proof root files " << endl;
  system("mkdir backup_outputProof`date +\"%d-%m-%y_%H-%M\"`;mv proof*.root  backup_outputProof`date +\"%d-%m-%y_%H-%M\"`/.");
  
  cout<<"###############################################################"<<endl;
  cout<<"################ 	   May your job 	##############"<<endl;
  cout<<"################      Live long and prosper	##############"<<endl;
  cout<<"###############################################################"<<endl;
  cout<< "  							      "<< endl;
  cout<< "  			     _  			      "<< endl;
  cout<< "  			  .-T |   _			      "<< endl;
  cout<< "  			  | | |  / |			      "<< endl;
  cout<< "  			  | | | / /`|			      "<< endl;
  cout<< "  		       _  | | |/ / /			      "<< endl;
  cout<< "  		       \\`\\| \'.\' / / 		      "<< endl;
  cout<< "  			\\ \\`-. \'--|  		      "<< endl;
  cout<< "  			 \\    \'   |			      "<< endl;
  cout<< "  			  \\ \\  .` /			      "<< endl;
  cout<< "  			    |	 |			      "<< endl;
  cout<< "  							      "<< endl;
  cout<< "  							      "<< endl;
  cout<<"###############################################################"<<endl;
  cout<<"###############################################################"<<endl;
    
  
  return (0);

}
