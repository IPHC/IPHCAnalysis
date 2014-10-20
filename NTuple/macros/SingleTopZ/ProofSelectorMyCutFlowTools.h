#ifndef ProofSelectorMyCutFlowTools_h
#define ProofSelectorMyCutFlowTools_h

#include "NTFormat/interface/NTEvent.h"

class ProofSelectorMyCutFlowTools{
 public :
  static std::vector< double > determineWeights(TString, double , double);

  static void determineLeptonCandidates(
                bool UseLooseWcand, float looseIsoCut, double rhocorr,
                std::vector<NTElectron> *selE,        std::vector<NTMuon> *selM,
                std::vector<NTElectron> *selENonIso,  std::vector<NTMuon> *selMNonIso,
                std::vector<NTElectron> *theZeeCand,  std::vector<NTMuon> *theZmumuCand,
                std::vector<NTElectron> *theWeCand,   std::vector<NTMuon> *theWmuCand,
		float tightIso_e, float tightIso_mu
                );

};




//_____________________________________________________________________________
std::vector<double> ProofSelectorMyCutFlowTools::determineWeights(TString thedatasetName, double weightITypeMC, double WeightForBranchingRatio){
  
  
  
  //------------------------------------
  //to calculate the event weight
  //returns a vector of double,
  //containing the infor for reweighting
  //------------------------------------

      
    double ITypeMC = 0;
    double Dweight = 0 ;
    double EventYieldWeightError = 0 ;
    double IsSignal = -1 ;
    
    Dweight= weightITypeMC * WeightForBranchingRatio;
    EventYieldWeightError = Dweight*Dweight;
    
      
    if ( thedatasetName=="TTbar" ) {
      ITypeMC = -1;
      Dweight= weightITypeMC * WeightForBranchingRatio;
      EventYieldWeightError = Dweight*Dweight;
    }
    else {
      IsSignal = -1;
      Dweight= weightITypeMC;
      EventYieldWeightError = Dweight*Dweight;
      if ( thedatasetName=="Zjets" || thedatasetName=="DYToLL_M10-50") ITypeMC = 2;
      else if ( thedatasetName=="Wjets" ) ITypeMC = 3;
      else if ( thedatasetName=="SingleToptW" || thedatasetName=="TtW" || thedatasetName=="TbartW"
		|| thedatasetName=="TtWScaleUp" || thedatasetName=="TtWScaleDown"
		|| thedatasetName=="TbartWScaleUp" || thedatasetName=="TbartWScaleDown") ITypeMC = 4;
      else if ( thedatasetName=="WZ" || thedatasetName=="WW" || thedatasetName=="ZZ" 
             || thedatasetName=="WZ_scaleup"|| thedatasetName=="WZ_scaledown"  
             || thedatasetName=="WZ_matchup"|| thedatasetName=="WZ_matchdown" ) ITypeMC = 5;
      else if ( 
       		  thedatasetName=="FCNCkut" || thedatasetName=="FCNCkct" 
     		||thedatasetName=="FCNCxut" || thedatasetName=="FCNCxct" 
     		||thedatasetName=="FCNCzut" || thedatasetName=="FCNCzct" 
     
     		||thedatasetName=="FCNCkutFullSim" || thedatasetName=="FCNCkctFullSim" 
     		||thedatasetName=="FCNCxutFullSim" || thedatasetName=="FCNCxctFullSim" 
    		||thedatasetName=="FCNCzutFullSim" || thedatasetName=="FCNCzctFullSim" 
     
     		||thedatasetName=="FCNCkut_matchup" || thedatasetName=="FCNCkct_matchup" 
     		||thedatasetName=="FCNCxut_matchup" || thedatasetName=="FCNCxct_matchup" 
     		||thedatasetName=="FCNCzut_matchup" || thedatasetName=="FCNCzct_matchup" 
     
     		||thedatasetName=="FCNCkut_matchdown" || thedatasetName=="FCNCkct_matchdown" 
     		||thedatasetName=="FCNCxut_matchdown" || thedatasetName=="FCNCxct_matchdown" 
     		||thedatasetName=="FCNCzut_matchdown" || thedatasetName=="FCNCzct_matchdown" 
     
     		||thedatasetName=="FCNCkut_scaleup" || thedatasetName=="FCNCkct_scaleup" 
    		||thedatasetName=="FCNCxut_scaleup" || thedatasetName=="FCNCxct_scaleup" 
     		||thedatasetName=="FCNCzut_scaleup" || thedatasetName=="FCNCzct_scaleup" 
     
    		||thedatasetName=="FCNCkut_scaledown" || thedatasetName=="FCNCkct_scaledown" 
     		||thedatasetName=="FCNCxut_scaledown" || thedatasetName=="FCNCxct_scaledown" 
     		||thedatasetName=="FCNCzut_scaledown" || thedatasetName=="FCNCzct_scaledown" 
     
     		||thedatasetName=="FCNCkut_topup" || thedatasetName=="FCNCkct_topup" 
     		||thedatasetName=="FCNCxut_topup" || thedatasetName=="FCNCxct_topup" 
     		||thedatasetName=="FCNCzut_topup" || thedatasetName=="FCNCzct_topup" 
     
     		||thedatasetName=="FCNCkut_topdown" || thedatasetName=="FCNCkct_topdown" 
     		||thedatasetName=="FCNCxut_topdown" || thedatasetName=="FCNCxct_topdown" 
   
      )  ITypeMC = 6; 
      else if(thedatasetName=="tZq") ITypeMC = 7;
      else if(thedatasetName=="TT" )ITypeMC = 8;
   }
  
  if ( thedatasetName=="DataDiEG" || thedatasetName=="DataDiMu" ||
       thedatasetName=="DataDiMuEG" || thedatasetName=="DataDiEGMu") {
     ITypeMC = 100;
     Dweight= weightITypeMC;
     EventYieldWeightError = Dweight*Dweight;
  }




  
  std::vector< double > thereturn;
  thereturn.push_back(ITypeMC);
  thereturn.push_back(Dweight);
  thereturn.push_back(EventYieldWeightError);
  thereturn.push_back(IsSignal);
  
  
  return thereturn;
}






 
void  ProofSelectorMyCutFlowTools::determineLeptonCandidates(
  		bool UseLooseWcand, float looseIsoCut, double rhocorr,
  		std::vector<NTElectron> *selE,        std::vector<NTMuon> *selM, 
  		std::vector<NTElectron> *selENonIso,  std::vector<NTMuon> *selMNonIso, 
		std::vector<NTElectron> *theZeeCand,  std::vector<NTMuon> *theZmumuCand, 
		std::vector<NTElectron> *theWeCand,   std::vector<NTMuon> *theWmuCand,
		float tightIso_e, float tightIso_mu
		){
 
      //*****************************************************************
      // select Z->ee candidate
      //*****************************************************************    
     
  int leptonFlavor = 0;
  int wcharge	   = 0;
  
  
  theZeeCand->clear(); 
  theZmumuCand->clear(); 
  
  theWeCand->clear(); 
  theWmuCand->clear(); 

  if(selE->size() >=2 ) {
    int theel1 = -1;
    int theel2 = -1;
    double mInv = 1000000;
    for(unsigned int iel1 = 0; iel1 < selE->size(); iel1++){
      for(unsigned int iel2 = 0; iel2 < selE->size(); iel2++){
    	if(iel1 == iel2) continue;
    	if((*selE)[iel1].charge == (*selE)[iel2].charge) continue;
    	TLorentzVector theZee = (*selE)[iel1].p4Gsf + (*selE)[iel2].p4Gsf;
    	if( fabs(theZee.M() - 91) < fabs(mInv-91) ){
    	//if( fabs(theZeeCand.M() - 200) < fabs(mInv-200) ){
    	  theel1 = iel1;
    	  theel2 = iel2;
    	  mInv = theZee.M();
    	}
      }
    }

    if(theel1>=0 && theel2>=0){ //JLA
      theZeeCand->push_back((*selE)[theel1]);
      theZeeCand->push_back((*selE)[theel2]);
      //double invM = (theZeeCand[0].p4Gsf+theZeeCand[1].p4Gsf).M();
      //cout << "lepton origin of Zee cand " << selE[theel1].LeptonOrigin << "  " << selE[theel2].LeptonOrigin << "  invmass " <<  invM << endl;
    }
  }
     
  //*****************************************************************
  // select W->enu candidate
  //*****************************************************************	 
  //cout << "get lepton cand W->enu " << endl;
  if(!UseLooseWcand){
    for(unsigned int iel1 = 0; iel1 < selE->size(); iel1++){
      bool matchElec=false;
      for(unsigned int iel2 = 0; iel2 < theZeeCand->size(); iel2++){
  	 
    	 if(fabs((*selE)[iel1].p4Gsf.Pt() - (*theZeeCand)[iel2].p4Gsf.Pt()) <  0.0001)  matchElec=true;
       }
      if(!matchElec && Selection::EffArea03PF( (*selE)[iel1], rhocorr) < tightIso_e){
    	theWeCand->push_back((*selE)[iel1]);
    	wcharge = (*selE)[iel1].charge;
    	if( (*selE)[iel1].LeptonOrigin == 10) leptonFlavor = 1;
      }
    }
  }else{
    for(unsigned int iel1 = 0; iel1 < selENonIso->size(); iel1++){
      bool matchElec=false;
      for(unsigned int iel2 = 0; iel2 < theZeeCand->size(); iel2++){
    	 if(fabs( (*selENonIso)[iel1].p4Gsf.Pt() - (*theZeeCand)[iel2].p4Gsf.Pt()) <  0.0001)  matchElec=true;
    	 else if( (*selE)[iel1].LeptonOrigin == 10) leptonFlavor = 1;
      }
      //if(!matchElec && selENonIso[iel1].RelIso03PF() > looseIsoCut){ 
      //if(!matchElec) cout << " eleciso " << Selection::EffArea03PF( (*selENonIso)[iel1], rhocorr) << " looseIsoCut " <<looseIsoCut << endl;
      if(!matchElec && Selection::EffArea03PF( (*selENonIso)[iel1], rhocorr) > looseIsoCut){ 
    	theWeCand->push_back((*selENonIso)[iel1]);
    	wcharge = (*selENonIso)[iel1].charge; 
    	if( (*selENonIso)[iel1].LeptonOrigin == 10) leptonFlavor = 1;
    	//cout << "    lepton origin of We cand " << selENonIso[iel1].LeptonOrigin << endl;
    	//if(selENonIso[iel1].LeptonOrigin == 10) cout << " nofake wenu" << endl;
      }
    }
  }
  
  
  //*****************************************************************
  // select Z->mumu candidate
  //*****************************************************************  
  if(selM->size() >=2 ) {
    int themu1 = -1;
    int themu2 = -1;
    double mInv = 1000000;
    for(unsigned int imu1 = 0; imu1 < selM->size(); imu1++){
      for(unsigned int imu2 = 0; imu2 < selM->size(); imu2++){
    	if(imu1 == imu2) continue;
    	if((*selM)[imu1].charge == (*selM)[imu2].charge) continue;
    	TLorentzVector theZmumu = (*selM)[imu1].p4 + (*selM)[imu2].p4;
    	if( fabs(theZmumu.M() - 91) < fabs(mInv-91) ){
    	//if( fabs(theZmumuCand.M() - 200) < fabs(mInv-200) ){
    	  themu1 = imu1;
    	  themu2 = imu2;
    	  mInv = theZmumu.M();
    	}
      }
    }
    if(themu1>=0 && themu2>=0){ //JLA  
      theZmumuCand->push_back((*selM)[themu1]);
      theZmumuCand->push_back((*selM)[themu2]);
    }
     
  }
  
   
   
  //*****************************************************************
  // select W->munu candidate
  //*****************************************************************	 
  
  //cout << "get lepton cand W->munu " << endl;
  if(!UseLooseWcand){
    //cout << "in sel W " << endl;
    for(unsigned int imu1 = 0; imu1 < selM->size(); imu1++){
      bool matchMuon = false;
      for(unsigned int imu2 = 0; imu2 < theZmumuCand->size(); imu2++){
  	 
    	 if(fabs((*selM)[imu1].p4.Pt() -(* theZmumuCand)[imu2].p4.Pt()) <  0.0001) matchMuon = true;
    	 
      } 
      if(!matchMuon &&  Selection::RelIso04PFDeltaBeta( (*selM)[imu1]) < tightIso_mu) {
    	theWmuCand->push_back( (*selM)[imu1]);
    	wcharge = (*selM)[imu1].charge;
    	if( (*selM)[imu1].LeptonOrigin == 10) leptonFlavor = 1;
      }
    }
  }else{
    for(unsigned int imu1 = 0; imu1 < selMNonIso->size(); imu1++){
      bool matchMuon = false;
      for(unsigned int imu2 = 0; imu2 < theZmumuCand->size(); imu2++){
  	 
    	 if(fabs( (*selMNonIso)[imu1].p4.Pt() - (*theZmumuCand)[imu2].p4.Pt())  < 0.0001) matchMuon = true;
     
      }
      //if(!matchMuon && selMNonIso[imu1].RelIso03PF() > looseIsoCut){
      if(!matchMuon && Selection::RelIso04PFDeltaBeta( (*selMNonIso)[imu1]) > looseIsoCut){
    	theWmuCand->push_back( (*selMNonIso)[imu1]);
    	wcharge = (*selMNonIso)[imu1].charge;
    	if( (*selMNonIso)[imu1].LeptonOrigin == 10) leptonFlavor = 1;
      }
    }
  }
  
  
  /*
  //redefine the lepton coll for jet cleaning
  if(UseLooseWcand){
     
    vector<NTElectron>  tmpElectrons = sel.GetSelectedElectronsNoIso();
    vector<NTMuon>	tmpMuons     = sel.GetSelectedMuonsNoIso();
    
    for(unsigned int iel=0; iel<tmpElectrons.size(); iel++){
      if(sel.EffArea03PF(tmpElectrons[iel], rho ) > 0.4) selE->push_back(tmpElectrons[iel]);
    }
  
    for(unsigned int imu=0; imu<tmpMuons.size(); imu++){
      if(sel.RelIso03PFDeltaBeta(tmpMuons[imu]) > 0.4) selM->push_back(tmpMuons[imu]);
    }
  
  }else{
  
    selE = sel.GetSelectedElectrons();
    selM	= sel.GetSelectedMuons();
  }*/
  
}



#endif
