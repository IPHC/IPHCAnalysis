// Revision du 25 Janv 2012

#include "../interface/DiLeptonSelection.h"


DiLeptonSelection::DiLeptonSelection ()
{
  selCode_ = -1;
  MinMassCut_ = 0.;
  ZMassWindow_ = pair < float, float >(9999., 0.);
  METCuts_ = pair < float, float >(9999., 0.);
  btagAlgo_ = -1;
  btagDiscriCut_ = -999.;
  NofBtagJets_ = 0;
  //Fill the table of cuts
  cuts_.push_back ("All");
  cuts_.push_back ("All dilept");
  cuts_.push_back ("Trigger");
  cuts_.push_back ("DiLepton pair");
  cuts_.push_back ("Z mass veto");
  cuts_.push_back ("NJets cut");
  cuts_.push_back ("MET cuts");
  cuts_.push_back ("NbtagJets cut");
  //Fill Channels
  channels_.push_back (string ("ee"));
  channels_.push_back (string ("emu"));
  channels_.push_back (string ("mumu"));

}


int DiLeptonSelection::GetChannel (string & CandPairType)
{
  for (unsigned int i = 0; i < channels_.size (); i++) {
    if (channels_[i] == CandPairType)
      return (int) i;
  }
  return -999;
}

DiLeptonSelection::DiLeptonSelection (const DiLeptonSelection & s):Selection (s)
{
  selCode_ = s.selCode_;
  channel_ = s.channel_;	//
  MinMassCut_ = s.MinMassCut_;
  METCuts_ = s.METCuts_;
  ZMassWindow_ = s.ZMassWindow_;
  cuts_ = s.cuts_;
  channels_ = s.channels_;	//
  btagAlgo_ = s.btagAlgo_;
  btagDiscriCut_ = s.btagDiscriCut_;
  NofBtagJets_ = s.NofBtagJets_;

}

DiLeptonSelection::~DiLeptonSelection ()
{
}


void DiLeptonSelection::SetParameters (float MinValue, pair < float, float >METCuts, pair < float, float >ZMassWindow, int btagAlgo, float btagDiscriCut, int NofBtagJets)
{
  MinMassCut_ = MinValue;
  METCuts_ = METCuts;
  ZMassWindow_ = ZMassWindow;
  btagAlgo_ = btagAlgo;
  btagDiscriCut_ = btagDiscriCut;
  NofBtagJets_ = NofBtagJets;
}



void DiLeptonSelection::LoadEvent (const NTEvent * event)
{
  Event::LoadEvent (event);
  channel_.LoadEvent (event);
}

bool DiLeptonSelection::GetLeptonPair (std::vector < NTMuon > muon_in, std::vector < IPHCTree::NTElectron > elec_in, std::vector < NTMuon > &muon_out, std::vector < NTElectron > &elec_out,
				       string & CandPairType, bool isForMM, float iso1_in, float iso2_in)
{
  
  
  
  
  //important: reset the out collections
  muon_out.clear ();
  elec_out.clear ();
  
  float sum_pT_ee = 0.;
  bool pass_elec = false;
  int ie1 = -1;
  int ie2 = -1;
  if (elec_in.size () >= 2) {
    for (unsigned int i = 0; i < elec_in.size (); i++) {
      for (unsigned int j = i + 1; j < elec_in.size (); j++) {
	if (pass_elec)
	  continue;
	if ( elec_in[i].charge != elec_in[j].charge) {
	  pass_elec = true;
	  sum_pT_ee = elec_in[i].p4.Pt () + elec_in[j].p4.Pt ();
	  ie1 = i;
	  ie2 = j;
	}
      }
    }
  }

  float sum_pT_mumu = 0.;
  bool pass_muon = false;
  int imu1 = -1;
  int imu2 = -1;
  if (muon_in.size () >= 2) {
    for (unsigned int i = 0; i < muon_in.size (); i++) {
      for (unsigned int j = i + 1; j < muon_in.size (); j++) {
	if (pass_muon)
	  continue;
	if (muon_in[i].charge != muon_in[j].charge) {	  pass_muon = true;
	  sum_pT_mumu = muon_in[i].p4.Pt () + muon_in[j].p4.Pt ();
	  imu1 = i;
	  imu2 = j;
	}
      }
    }
  }


  float sum_pT_emu_start = 0.;
  float sum_pT_emu = 0.;
  int je1 = -1;
  int jmu2 = -1;
  if (muon_in.size () >= 1 && elec_in.size () >= 1) {
    for (unsigned int i = 0; i < muon_in.size (); i++) {
      for (unsigned int j = 0; j < elec_in.size (); j++) {
	if (muon_in[i].charge != elec_in[j].charge){
	  sum_pT_emu = muon_in[i].p4.Pt () + elec_in[j].p4.Pt ();
	  if (sum_pT_emu > sum_pT_emu_start) {
	    sum_pT_emu_start = sum_pT_emu;
	    je1 = j;
	    jmu2 = i;
	  }
	}
      }
    }
  }


  float sum[3] = { sum_pT_ee, sum_pT_mumu, sum_pT_emu };
  int sortedIndices[3];
  TMath::Sort (3, sum, sortedIndices);
  if (sortedIndices[0] == 0 && sum_pT_ee != 0.) {
    elec_out.push_back (elec_in[ie1]);
    elec_out.push_back (elec_in[ie2]);
  }
  else if (sortedIndices[0] == 1 && sum_pT_mumu != 0.) {
    muon_out.push_back (muon_in[imu1]);
    muon_out.push_back (muon_in[imu2]);
  }
  else if (sortedIndices[0] == 2 && sum_pT_emu != 0.) {
    elec_out.push_back (elec_in[je1]);
    muon_out.push_back (muon_in[jmu2]);
  }





  if (elec_out.size () + muon_out.size () == 2) {
    if (muon_out.size () == 2) CandPairType = "mumu";
    if (elec_out.size () == 2) CandPairType = "ee";
    if (muon_out.size () == 1 && elec_out.size () == 1) CandPairType = "emu";
  }
  if (CandPairType == "ee" || CandPairType == "emu" || CandPairType == "mumu")
    return true;
  else
    return false;

}






bool DiLeptonSelection::GetLeptonPairForMM(float iso1_in , float iso2_in, std::vector<NTMuon> muon_in,std::vector<IPHCTree::NTElectron> elec_in, std::vector<NTMuon>& muon_out,std::vector<NTElectron>& elec_out, string& CandPairType){
  
  muon_out.clear();
  elec_out.clear();
  
  float sum_pT_ee = 0.;
  bool pass_elec = false;
  int ie1 = -1;
  int ie2 = -1;
  if ( elec_in.size()>=2 ) {
    for(unsigned int i=0; i<elec_in.size(); i++){
      for(unsigned int j=i+1; j<elec_in.size(); j++){
        if ( pass_elec ) continue;
        if ( (elec_in[i].charge != elec_in[j].charge) 
        && ( (RelIso03PF(elec_in[i])<iso1_in && RelIso03PF(elec_in[j])<iso2_in)  || (RelIso03PF(elec_in[i])<iso2_in && RelIso03PF(elec_in[j])<iso1_in) )  
        ) {
            pass_elec = true;
            sum_pT_ee = elec_in[i].p4.Pt()+elec_in[j].p4.Pt();
            ie1 = i; ie2 = j;
          }
        }
      }
   }
  
  
  
  
  float sum_pT_mumu = 0.;
  bool pass_muon = false;
  int imu1 = -1;
  int imu2 = -1;
  if ( muon_in.size()>=2 ) {
    for(unsigned int i=0; i<muon_in.size(); i++){
      for(unsigned int j=i+1; j<muon_in.size(); j++){
        if ( pass_muon ) continue;
        if ( (muon_in[i].charge != muon_in[j].charge) 
        && ( (RelIso03PFDeltaBeta(muon_in[i])<iso1_in && RelIso03PFDeltaBeta(muon_in[j])<iso2_in)  || (RelIso03PFDeltaBeta(muon_in[i])<iso2_in && RelIso03PFDeltaBeta(muon_in[j])<iso1_in) )                
           ) {
            pass_muon = true;
            sum_pT_mumu = muon_in[i].p4.Pt()+muon_in[j].p4.Pt();
            imu1 = i; imu2 = j;
          }
       }
     }
   }

  
  float sum_pT_emu_start = 0.;
  float sum_pT_emu       = 0.;
  int je1  = -1;
  int jmu2 = -1;
  if ( muon_in.size()>=1 && elec_in.size()>=1 ) {
    for(unsigned int i=0; i<muon_in.size(); i++){
      for(unsigned int j=0; j<elec_in.size(); j++){
	if ( (muon_in[i].charge != elec_in[j].charge) 
	   && ( (RelIso03PFDeltaBeta(muon_in[i])<iso1_in && RelIso03PF(elec_in[j])<iso2_in) )                      ){
         sum_pT_emu = muon_in[i].p4.Pt()+elec_in[j].p4.Pt();
         if (sum_pT_emu>sum_pT_emu_start)
         {
           sum_pT_emu_start = sum_pT_emu;
           je1 = j; jmu2 = i;
         }
      }
    }
   }
 }
 
 
 float sum[3] = {sum_pT_ee, sum_pT_mumu, sum_pT_emu};
 int sortedIndices[3];
 TMath::Sort(3,sum,sortedIndices);
 if      (sortedIndices[0]==0 && sum_pT_ee   !=0.){ elec_out.push_back(elec_in[ie1]);  elec_out.push_back(elec_in[ie2]); }
 else if (sortedIndices[0]==1 && sum_pT_mumu !=0.){ muon_out.push_back(muon_in[imu1]); muon_out.push_back(muon_in[imu2]);}
 else if (sortedIndices[0]==2 && sum_pT_emu  !=0.){ elec_out.push_back(elec_in[je1]);  muon_out.push_back(muon_in[jmu2]);}
 
 if(elec_out.size()+muon_out.size()==2){
 	if(muon_out.size()==2){
		if(fabs(muon_out[0].p4.Eta())<2.1 || fabs(muon_out[1].p4.Eta())<2.1){ CandPairType = "mumu"; }
		else{
			CandPairType = "false";
		}	
	}
 	if(elec_out.size()==2){ CandPairType = "ee";	}
 	if(muon_out.size()==1 && elec_out.size()==1){
// necessary ??????? to be verified !!    
//		if(elec_out[0].ET_SC<17 && fabs(muon_out[0].p4.Eta())>2.1){
		if( fabs(muon_out[0].p4.Eta())>2.1){
			CandPairType = "false";
		}
		else {CandPairType = "emu";}	
	}
 }
 if(CandPairType == "ee" || CandPairType == "emu" || CandPairType == "mumu" ) return true; 
 else return false; 

}


bool DiLeptonSelection::GetLeptonPairForMM(float iso1_in , float iso2_in, std::vector<NTMuon>& muon_out,std::vector<NTElectron>& elec_out, string& CandPairType){
  return GetLeptonPairForMM(iso1_in , iso2_in, GetSelectedMuonsNoIso(),  GetSelectedElectronsNoIso(),  muon_out, elec_out, CandPairType);
}





bool DiLeptonSelection::GetLeptonPair (std::vector < NTMuon > &muon_out, std::vector < NTElectron > &elec_out, string & CandPairType, bool isForMM, float iso1_in, float iso2_in, float rho)
{ 
  
  return GetLeptonPair (GetSelectedMuonsDileptonTTbar (), GetSelectedElectronsDileptonTTbar (), muon_out, elec_out, CandPairType, isForMM, iso1_in, iso2_in);
}


bool DiLeptonSelection::GetLeptonPairElectronScaled (std::vector < NTMuon > &muon_out, std::vector < NTElectron > &elec_out, string & CandPairType, float scale, bool isForMM, float iso1_in, float iso2_in, float rho)
{
  
  return GetLeptonPair (GetSelectedMuonsDileptonTTbar (), GetSelectedElectronsDileptonTTbar (true, scale, false, 1, rho), muon_out, elec_out, CandPairType, isForMM, iso1_in, iso2_in);
}

bool DiLeptonSelection::GetLeptonPairElectronSmeared (std::vector < NTMuon > &muon_out, std::vector < NTElectron > &elec_out, string & CandPairType,  float resol, bool isForMM, float iso1_in, float iso2_in, float rho)
{
  return GetLeptonPair (GetSelectedMuonsDileptonTTbar (), GetSelectedElectronsDileptonTTbar (false, 1, true, resol, rho), muon_out, elec_out, CandPairType, isForMM, iso1_in, iso2_in);
}



bool DiLeptonSelection::TestIsolationOfPair (float iso1_in, float iso2_in, std::vector < NTMuon > muon_in, std::vector < IPHCTree::NTElectron > elec_in)
{

  bool pass_elec = false;
  if (elec_in.size () == 2) {
    for (unsigned int i = 0; i < elec_in.size (); i++) {
      for (unsigned int j = i + 1; j < elec_in.size (); j++) {
	if (pass_elec)
	  continue;
	if ((RelIso03PF(elec_in[i]) < iso1_in && RelIso03PF(elec_in[j]) < iso2_in) || (RelIso03PF(elec_in[i]) < iso2_in && RelIso03PF(elec_in[j]) < iso1_in)) {
	  pass_elec = true;
	}
      }
    }
  }



  bool pass_muon = false;
  if (muon_in.size () == 2) {
    for (unsigned int i = 0; i < muon_in.size (); i++) {
      for (unsigned int j = i + 1; j < muon_in.size (); j++) {
	if (pass_muon)
	  continue;
	if ((RelIso03PFDeltaBeta(muon_in[i]) < iso1_in && RelIso03PFDeltaBeta(muon_in[j]) < iso2_in) || (RelIso03PFDeltaBeta(muon_in[i]) < iso2_in && RelIso03PFDeltaBeta(muon_in[j]) < iso1_in)) {
	  pass_muon = true;
	}
      }
    }
  }


  bool pass_emu = false;
  if ((muon_in.size () + elec_in.size ()) == 2) {
    for (unsigned int i = 0; i < muon_in.size (); i++) {
      for (unsigned int j = 0; j < elec_in.size (); j++) {
	if (pass_emu)
	  continue;
	if (RelIso03PFDeltaBeta(muon_in[i]) < iso1_in && RelIso03PF(elec_in[j]) < iso2_in) {
	  pass_emu = true;
	}
      }
    }
  }

  if (pass_elec || pass_muon || pass_emu)
    return true;
  else
    return false;

}


TLorentzVector DiLeptonSelection::DiLeptonCand (const std::vector < NTMuon > &muons_cand, const std::vector < NTElectron > &electrons_cand)
{
  TLorentzVector DiLepton;
  if (muons_cand.size () + electrons_cand.size () != 2)
    return DiLepton;
  for (unsigned int i = 0; i < muons_cand.size (); i++)
    DiLepton += muons_cand[i].p4;
  for (unsigned int i = 0; i < electrons_cand.size (); i++)
    DiLepton += electrons_cand[i].p4;
  return DiLepton;
}

float DiLeptonSelection::DiLeptonMass (const std::vector < NTMuon > &muons_cand, const std::vector < NTElectron > &electrons_cand)
{
  return DiLeptonCand (muons_cand, electrons_cand).M ();
}

float DiLeptonSelection::DiLeptonMT (const std::vector < NTMuon > &muons_cand, const std::vector < NTElectron > &electrons_cand)
{
  return DiLeptonCand (muons_cand, electrons_cand).Mt ();
}

bool DiLeptonSelection::DiLeptonMassCut (float MinValue, pair < float, float >ZMassWindow, const std::vector < NTMuon > &muons_cand, const std::vector < NTElectron > &electrons_cand,
					 string channelName)
{


  string leptonpairname;
  if (channelName != string ("") && channelName != string ("*") && channelName != string ("allChannels")) {
    leptonpairname = channelName;
  }
  else {
    if (electrons_cand.size () == 2)
      leptonpairname = "ee";
    if (muons_cand.size () == 2)
      leptonpairname = "mumu";
    if (muons_cand.size () == 1 && electrons_cand.size () == 1)
      leptonpairname = "emu";
  }

  bool iresult = true;
  float mass = DiLeptonMass (muons_cand, electrons_cand);

  // reject low DY mass;
  if (mass < MinValue)
    iresult = false;
  // Z veto;
  if (leptonpairname != "emu" && ZMassWindow.first < mass && mass < ZMassWindow.second)
    iresult = false;

  return iresult;

}

bool DiLeptonSelection::DiLeptonMassCut (const std::vector < NTMuon > &muons_cand, const std::vector < NTElectron > &electrons_cand, string channelName)
{
  return DiLeptonMassCut (MinMassCut_, ZMassWindow_, muons_cand, electrons_cand, channelName);
}

int DiLeptonSelection::doFullSelection (Dataset * dataset, string channelName,
	bool applyTTbarMCCut, bool print, bool isForMM, bool MMInverted,
	float iso1_in, float iso2_in, bool applyJES, float JESParam,
	bool applyEEScale, float EEScaleParam, bool applyEEResol, float EEResolParam,
	bool applyMEScale, float MEScaleParam,
	bool applyJER, float JERFactor,	bool applyMETS, float METScale)
{
  vector < float >weightb;
  weightb.push_back (1.);
  weightb.push_back (0.);
  weightb.push_back (0.);
  weightb.push_back (0.);
  weightb.push_back (0.);
  int FinalStep = doFullSelection (dataset, weightb, channelName, applyTTbarMCCut, print,
	isForMM, MMInverted, iso1_in, iso2_in, applyJES, JESParam, 
	applyEEScale, EEScaleParam, applyEEResol, EEResolParam,
	applyMEScale, MEScaleParam,
	applyJER, JERFactor, applyMETS, METScale);
  return FinalStep;
}

int DiLeptonSelection::doFullSelection (Dataset * dataset, vector < float >&weightb,
	string channelName, bool applyTTbarMCCut, bool print, bool isForMM,
	bool MMInverted, float iso1_in, float iso2_in, bool applyJES, float JESParam,
	bool applyEEScale, float EEScaleParam, bool applyEEResol, float EEResolParam,
	bool applyMEScale, float MEScaleParam,
	bool applyJER, float JERFactor, bool applyMETS, float METScale)
{
 
  //clear object collections
  jetsAna.clear();
  bjetsAna.clear();
  electronsAna.clear();
  muonsAna.clear();

  bool applyChannel = false;
  if (channelName != string ("") && channelName != string ("*") && channelName != string ("allChannels"))
    applyChannel = true;
  bool filterTTbarMC = true;
  if (applyTTbarMCCut && !channel_.isSignal (channelName))
    filterTTbarMC = false;
  //cout << "applyTTbarMCCut " << applyTTbarMCCut <<  "  " << channel_.isSignal (channelName)  << endl;
  //boolean for the selection step: true = pass the cut
  bool step_trigger = false;
  bool step_pair = false;
  bool step_Zveto = false;
  bool step_jets = false;
  bool step_met = false;
  bool step_bjets = false;

  TString dump;
  ostringstream runNumber_oss;
  ostringstream eventNumber_oss;
  runNumber_oss << getRunNumber();
  eventNumber_oss << getEventNumber();
  dump = channel_.ChannelName () + string (" | ") + runNumber_oss.str () + string (" | ") + eventNumber_oss.str () + string (" | ");
  //dump+=string("")+runNumberS;  

  double METEMu = METCuts_.first;
  double METLL = METCuts_.second;
  std::vector < NTMuon > muon_cand;
  std::vector < NTElectron > elec_cand;
  pairType_ = "";
  dimass_=0.;

  int FinalStep = 0;
  //Step 1        MC-match
  //cout << "filterTTbarMC " << filterTTbarMC << endl;
  if (filterTTbarMC) {
    FinalStep++;
    //Step 2        Trigger
    //cout << "test trigger " << FinalStep << endl;
    if (passTriggerSelection (dataset, channelName)) {
      step_trigger = true;
      //FinalStep++;
    }
    //Step 3        Dilepton pair choice
    if (((!isForMM && 
		GetLeptonPair (GetSelectedMuons(applyMEScale, MEScaleParam),
			GetSelectedElectrons(applyEEScale, EEScaleParam, applyEEResol, EEResolParam),
			muon_cand, elec_cand, pairType_) == true) || 
	( isForMM && !MMInverted && GetLeptonPair (GetSelectedMuonsNoIso (), GetSelectedElectronsNoIso (), muon_cand, elec_cand, pairType_, true, iso1_in, iso2_in) == true ) || 
	( isForMM && MMInverted && GetLeptonPair (GetSelectedMuonsNoIso (), GetSelectedElectronsNoIso (), muon_cand, elec_cand, pairType_, true, 100000, 100000) == true && TestIsolationOfPair (iso1_in, iso2_in, muon_cand, elec_cand)) ) 
	&& (!applyChannel || (applyChannel && pairType_ == channelName))) {
      if (!applyChannel && step_trigger) step_trigger = passTriggerSelection(dataset,pairType_);
      step_pair = true;
      //FinalStep++;
      muonsAna = muon_cand;
      electronsAna = elec_cand;

      float lep1PtxCharge = 0;
      float lep2PtxCharge = 0;
      float lep1RelIso = 0;
      float lep2RelIso = 0;
      ostringstream lep1PtxCharge_oss;
      ostringstream lep2PtxCharge_oss;
      ostringstream lep1RelIso_oss;
      ostringstream lep2RelIso_oss;
      if (pairType_ == "mumu") {
	lep1PtxCharge = muon_cand[0].p4.Pt () * muon_cand[0].charge;
	lep2PtxCharge = muon_cand[1].p4.Pt () * muon_cand[1].charge;
	lep1RelIso = RelIso03PFDeltaBeta(muon_cand[0]);
	lep2RelIso = RelIso03PFDeltaBeta(muon_cand[1]);
      }
      if (pairType_ == "ee") {
	lep1PtxCharge = elec_cand[0].p4.Pt () * elec_cand[0].charge;
	lep2PtxCharge = elec_cand[1].p4.Pt () * elec_cand[1].charge;
	lep1RelIso = RelIso03PF(elec_cand[0]);
	lep2RelIso = RelIso03PF(elec_cand[1]);
      }
      if (pairType_ == "emu") {
	if (elec_cand[0].p4.Pt () > muon_cand[0].p4.Pt ()) {

	  lep1PtxCharge = elec_cand[0].p4.Pt () * elec_cand[0].charge;
	  lep2PtxCharge = muon_cand[0].p4.Pt () * muon_cand[0].charge;
	  lep1RelIso = RelIso03PF(elec_cand[0]);
	  lep2RelIso = RelIso03PFDeltaBeta(muon_cand[0]);
	}
	else {
	  lep2PtxCharge = elec_cand[0].p4.Pt () * elec_cand[0].charge;
	  lep1PtxCharge = muon_cand[0].p4.Pt () * muon_cand[0].charge;
	  lep2RelIso = RelIso03PF(elec_cand[0]);
	  lep1RelIso = RelIso03PFDeltaBeta(muon_cand[0]);
	}
      }
      dimass_ = DiLeptonMass (muon_cand, elec_cand);
      ostringstream dimass_oss;
      dimass_oss << dimass_;
      lep1PtxCharge_oss << lep1PtxCharge;
      lep2PtxCharge_oss << lep2PtxCharge;
      lep1RelIso_oss << lep1RelIso;
      lep2RelIso_oss << lep2RelIso;
      dump += lep1PtxCharge_oss.str () + "," + lep2PtxCharge_oss.str () + " | " + lep1RelIso_oss.str () + "," + lep2RelIso_oss.str () + " | " + dimass_oss.str () + " | ";
      //Step 4     Z mass veto 
      if (DiLeptonMassCut (muon_cand, elec_cand, channelName)) {
	//FinalStep++;
	step_Zveto = true;
      }
      //Step 5    Minimal jet multiplicity 
      //vector < NTJet > SelectedJets0 = GetSelectedJetsResCorr (GetSelectedMuons (), GetSelectedElectrons ()); //-> ResCorr: has to be reintegrated correctly
      vector < NTJet > SelectedJets = GetSelectedJets (muon_cand, elec_cand, applyJES, JESParam, applyJER, JERFactor);
      // compute the weightb for b-tag
      if (GetFlagb() == 1) {
	if (!dataset->isData ()) {	//MC
	  vector < float >weight_temp = GetSFBweight().GetWeigth4BSel (GetMethodb(), GetSystb(), SelectedJets);
	  weightb[0] *= weight_temp[0];	//weight of the event
	  weightb[1] = weight_temp[1];	//proba 0 jet
	  weightb[2] = weight_temp[2];	//proba 1 jet
	  weightb[3] = weight_temp[3];	//proba 2 jets
	  weightb[4] = weight_temp[4];	//proba at least 3 jets
	}
      }
      jetsAna = SelectedJets;

      if (SelectedJets.size () >= 2) {
	//FinalStep++;
	step_jets = true;
      }
      //Step 6  MET cuts
      if (pairType_ == string ("emu")) {
	if (GetSelectedMET (applyJES, JESParam, applyMETS, METScale).p2.Mod () > METEMu) {
	  //FinalStep++;
	  step_met = true;
	}
      }
      else {
	if (GetSelectedMET (applyJES, JESParam, applyMETS, METScale).p2.Mod () > METLL) {
	  //FinalStep++;
	  step_met = true;
	}
      }

      //Step 7 btagging
      vector < NTJet > btagjets;
      vector < float >btagDiscri;
      ostringstream jet1pt_oss;
      ostringstream jet2pt_oss;
      ostringstream jet1bdisc_oss;
      ostringstream jet2bdisc_oss;
      ostringstream met_oss;
      if(step_jets){
      	jet1pt_oss << SelectedJets[0].p4.Pt ();
      	jet2pt_oss << SelectedJets[1].p4.Pt ();
      	jet1bdisc_oss << SelectedJets[0].bTag["trackCountingHighEffBJetTags"];
      	jet2bdisc_oss << SelectedJets[1].bTag["trackCountingHighEffBJetTags"];
     	 met_oss << GetMET ().p2.Mod ();
      	dump += jet1pt_oss.str () + " , " + jet2pt_oss.str () + " | " + met_oss.str () + " | " + jet1bdisc_oss.str () + " , " + jet2bdisc_oss.str () + " | " + pairType_;
      }

      for (unsigned int j = 0; j < SelectedJets.size (); j++) {
	if (passBtagSelection(SelectedJets[j])) {
	  btagjets.push_back (SelectedJets[j]);
	  btagDiscri.push_back (sfb_dilep_.getBtagDiscr(SelectedJets[j], btagAlgo_));
	}
      }
      bjetsAna = btagjets;
      if ((int) btagjets.size () >= NofBtagJets_) {
	//FinalStep++;
	step_bjets = true;
      }
      if (GetFlagb() == 1) {
	if (dataset->isData ()) {	//Data
	  if (step_bjets) {
              weightb[0] *= 1;
          }
          else  {
             weightb[0]=0;
          }
          for (int ibj=0; ibj<4; ibj++) {
	    if ((int) btagjets.size ()== ibj) weightb[ibj+1] = 1;	
            else weightb[ibj] = 0;
          }
          if ((int) btagjets.size ()>=3) weightb[4] = 1;
	}
      }
    }
  }

  //Compute a collection of jets and b-jets and a weight even if the selection of pair fails
  if(!step_pair){
      vector < NTJet > SelectedJets = GetSelectedJets (GetSelectedMuons(), GetSelectedElectrons(), applyJES, JESParam);
      // compute the weightb for b-tag
      if (GetFlagb() == 1) {
	if (!dataset->isData ()) {	//MC
	  vector < float >weight_temp = GetSFBweight().GetWeigth4BSel (GetMethodb(), GetSystb(), SelectedJets);
	  weightb[0] *= weight_temp[0];	//weight of the event
	  weightb[1] = weight_temp[1];	//proba 0 jet
	  weightb[2] = weight_temp[2];	//proba 1 jet
	  weightb[3] = weight_temp[3];	//proba 2 jets
	  weightb[4] = weight_temp[4];	//proba at least 3 jets
	}
      }
      jetsAna = SelectedJets;
      vector < NTJet > btagjets;
      vector < float >btagDiscri;
      for (unsigned int j = 0; j < SelectedJets.size (); j++) {
	if (passBtagSelection(SelectedJets[j])) {
	  btagjets.push_back (SelectedJets[j]);
	  btagDiscri.push_back (sfb_dilep_.getBtagDiscr(SelectedJets[j], btagAlgo_));
	}
      }
      bjetsAna = btagjets;
      if (GetFlagb() == 1) {
	if (dataset->isData ()) {	//Data
          if ((int) btagjets.size () >= NofBtagJets_) {
              weightb[0] *= 1;
          }
          else  {
             weightb[0]=0;
          }
          for (int ibj=0; ibj<4; ibj++) {
	    if ((int) btagjets.size ()== ibj) weightb[ibj+1] = 1;	
            else weightb[ibj] = 0;
          }
          if ((int) btagjets.size ()>=3) weightb[4] = 1;
	}
      }
  }


  //compute FinalStep
  if (step_trigger) {
    FinalStep++;
    if (step_pair) {
      FinalStep++;
      if (step_Zveto) {
	FinalStep++;
	if (step_jets) {
	  FinalStep++;
	  if (step_met) {
	    FinalStep++;
	    if (step_bjets) {
	      FinalStep++;
	      //the event is selected
	      if (print)
		cout << dump << endl;
	    }
	  }
	}
      }
    }
  }
  selCode_ = step_trigger*maskTriggerCut + step_pair*maskPairCut + 
	step_Zveto * maskZvetoCut + step_jets * maskJetCut + 
	step_met * maskMETCut + step_bjets * maskBjetsCut;

  return FinalStep;
}













int DiLeptonSelection::doFullSelectionForMM (float iso1_in, float iso2_in, Dataset * dataset, string channelName, bool applyTTbarMCCut, bool print, bool isForMM, bool MMInverted, bool applyJES, float JESParam, bool applyEES, float EESParam, bool applyMES, float MESParam, bool applyJER, float JERFactor, bool applyMETS, float METScale)
{
  vector < float >weightb;
  weightb.push_back (1.);
  weightb.push_back (0.);
  weightb.push_back (0.);
  weightb.push_back (0.);
  weightb.push_back (0.);
  int FinalStep = doFullSelectionForMM (iso1_in, iso2_in, dataset, weightb, channelName, applyTTbarMCCut, print, isForMM, MMInverted, applyJES, JESParam, applyEES, EESParam, applyMES, MESParam, applyJER, JERFactor, applyMETS, METScale);
  return FinalStep;
}

int DiLeptonSelection::doFullSelectionForMM (float iso1_in, float iso2_in, Dataset * dataset, vector < float >&weightb, string channelName, bool applyTTbarMCCut, bool print, bool isForMM, bool MMInverted,  bool applyJES, float JESParam, bool applyEES, float EESParam, bool applyMES, float MESParam, bool applyJER, float JERFactor, bool applyMETS, float METScale)
{
 
  //clear object collections
  jetsAna.clear();
  bjetsAna.clear();
  electronsAna.clear();
  muonsAna.clear();

  bool applyChannel = false;
  if (channelName != string ("") && channelName != string ("*") && channelName != string ("allChannels"))
    applyChannel = true;
  bool filterTTbarMC = true;
  if (applyTTbarMCCut && !channel_.isSignal (channelName))
    filterTTbarMC = false;
  //cout << "applyTTbarMCCut " << applyTTbarMCCut <<  "  " << channel_.isSignal (channelName)  << endl;
  //boolean for the selection step: true = pass the cut
  bool step_trigger = false;
  bool step_pair = false;
  bool step_Zveto = false;
  bool step_jets = false;
  bool step_met = false;
  bool step_bjets = false;

  TString dump;
  ostringstream runNumber_oss;
  ostringstream eventNumber_oss;
  runNumber_oss << getRunNumber();
  eventNumber_oss << getEventNumber();
  dump = channel_.ChannelName () + string (" | ") + runNumber_oss.str () + string (" | ") + eventNumber_oss.str () + string (" | ");
  //dump+=string("")+runNumberS;  

  double METEMu = METCuts_.first;
  double METLL = METCuts_.second;
  std::vector < NTMuon > muon_cand;
  std::vector < NTElectron > elec_cand;
  pairType_ = "";
  dimass_=0.;

  int FinalStep = 0;
  //Step 1        MC-match
  //cout << "filterTTbarMC " << filterTTbarMC << endl;
  if (filterTTbarMC) {
    FinalStep++;
    //Step 2        Trigger
    //cout << "test trigger " << FinalStep << endl;
    if (passTriggerSelection (dataset, channelName)) {
      step_trigger = true;
      //FinalStep++;
    }
    //Step 3        Dilepton pair choice
    if (GetLeptonPairForMM(0.8, 0.8, GetSelectedMuonsNoIso(), GetSelectedElectronsNoIso(), muon_cand, elec_cand, pairType_) == true && (!applyChannel || (applyChannel && pairType_==channelName))){
      if (TestIsolationOfPair(iso1_in, iso2_in, muon_cand, elec_cand)){

      step_pair = true;
      //FinalStep++;
      muonsAna = muon_cand;
      electronsAna = elec_cand;

      float lep1PtxCharge = 0;
      float lep2PtxCharge = 0;
      float lep1RelIso = 0;
      float lep2RelIso = 0;
      ostringstream lep1PtxCharge_oss;
      ostringstream lep2PtxCharge_oss;
      ostringstream lep1RelIso_oss;
      ostringstream lep2RelIso_oss;
      if (pairType_ == "mumu") {
	lep1PtxCharge = muon_cand[0].p4.Pt () * muon_cand[0].charge;
	lep2PtxCharge = muon_cand[1].p4.Pt () * muon_cand[1].charge;
	lep1RelIso = RelIso03PFDeltaBeta(muon_cand[0]);
	lep2RelIso = RelIso03PFDeltaBeta(muon_cand[1]);
      }
      if (pairType_ == "ee") {
	lep1PtxCharge = elec_cand[0].p4.Pt () * elec_cand[0].charge;
	lep2PtxCharge = elec_cand[1].p4.Pt () * elec_cand[1].charge;
	lep1RelIso = RelIso03PF(elec_cand[0]);
	lep2RelIso = RelIso03PF(elec_cand[1]);
      }
      if (pairType_ == "emu") {
	if (elec_cand[0].p4.Pt () > muon_cand[0].p4.Pt ()) {

	  lep1PtxCharge = elec_cand[0].p4.Pt () * elec_cand[0].charge;
	  lep2PtxCharge = muon_cand[0].p4.Pt () * muon_cand[0].charge;
	  lep1RelIso = RelIso03PF(elec_cand[0]);
	  lep2RelIso = RelIso03PFDeltaBeta(muon_cand[0]);
	}
	else {
	  lep2PtxCharge = elec_cand[0].p4.Pt () * elec_cand[0].charge;
	  lep1PtxCharge = muon_cand[0].p4.Pt () * muon_cand[0].charge;
	  lep2RelIso = RelIso03PF(elec_cand[0]);
	  lep1RelIso = RelIso03PFDeltaBeta(muon_cand[0]);
	}
      }
      dimass_ = DiLeptonMass (muon_cand, elec_cand);
      ostringstream dimass_oss;
      dimass_oss << dimass_;
      lep1PtxCharge_oss << lep1PtxCharge;
      lep2PtxCharge_oss << lep2PtxCharge;
      lep1RelIso_oss << lep1RelIso;
      lep2RelIso_oss << lep2RelIso;
      dump += lep1PtxCharge_oss.str () + "," + lep2PtxCharge_oss.str () + " | " + lep1RelIso_oss.str () + "," + lep2RelIso_oss.str () + " | " + dimass_oss.str () + " | ";
      //Step 4     Z mass veto 
      if (DiLeptonMassCut (muon_cand, elec_cand, channelName)) {
	//FinalStep++;
	step_Zveto = true;
      }
      //Step 5    Minimal jet multiplicity 
      //vector < NTJet > SelectedJets0 = GetSelectedJetsResCorr (GetSelectedMuons (), GetSelectedElectrons ()); //-> ResCorr: has to be reintegrated correctly
      vector < NTJet > SelectedJets = GetSelectedJets (muon_cand, elec_cand, applyJES, JESParam, applyJER, JERFactor);
      // compute the weightb for b-tag
      if (GetFlagb() == 1) {
	if (!dataset->isData ()) {	//MC
	  vector < float >weight_temp = GetSFBweight().GetWeigth4BSel (GetMethodb(), GetSystb(), SelectedJets);
	  weightb[0] *= weight_temp[0];	//weight of the event
	  weightb[1] = weight_temp[1];	//proba 0 jet
	  weightb[2] = weight_temp[2];	//proba 1 jet
	  weightb[3] = weight_temp[3];	//proba 2 jets
	  weightb[4] = weight_temp[4];	//proba at least 3 jets
	}
      }
      jetsAna = SelectedJets;

      if (SelectedJets.size () >= 2) {
	//FinalStep++;
	step_jets = true;
      }
      //Step 6  MET cuts
      if (pairType_ == string ("emu")) {
	if (GetSelectedMET (applyJES, JESParam, applyMETS, METScale).p2.Mod () > METEMu) {
	  //FinalStep++;
	  step_met = true;
	}
      }
      else {
	if (GetSelectedMET (applyJES, JESParam, applyMETS, METScale).p2.Mod () > METLL) {
	  //FinalStep++;
	  step_met = true;
	}
      }

      //Step 7 btagging
      vector < NTJet > btagjets;
      vector < float >btagDiscri;
      ostringstream jet1pt_oss;
      ostringstream jet2pt_oss;
      ostringstream jet1bdisc_oss;
      ostringstream jet2bdisc_oss;
      ostringstream met_oss;
      if(step_jets){
      	jet1pt_oss << SelectedJets[0].p4.Pt ();
      	jet2pt_oss << SelectedJets[1].p4.Pt ();
      	jet1bdisc_oss << SelectedJets[0].bTag["trackCountingHighEffBJetTags"];
      	jet2bdisc_oss << SelectedJets[1].bTag["trackCountingHighEffBJetTags"];
     	 met_oss << GetMET ().p2.Mod ();
      	dump += jet1pt_oss.str () + " , " + jet2pt_oss.str () + " | " + met_oss.str () + " | " + jet1bdisc_oss.str () + " , " + jet2bdisc_oss.str () + " | " + pairType_;
      }

      for (unsigned int j = 0; j < SelectedJets.size (); j++) {
	if (passBtagSelection(SelectedJets[j])) {
	  btagjets.push_back (SelectedJets[j]);
	  btagDiscri.push_back (sfb_dilep_.getBtagDiscr(SelectedJets[j], btagAlgo_));
	}
      }
      bjetsAna = btagjets;
      if ((int) btagjets.size () >= NofBtagJets_) {
	//FinalStep++;
	step_bjets = true;
      }
      if (GetFlagb() == 1) {
	if (dataset->isData ()) {	//Data
	  if (step_bjets) {
              weightb[0] *= 1;
          }
          else  {
             weightb[0]=0;
          }
          for (int ibj=0; ibj<4; ibj++) {
	    if ((int) btagjets.size ()== ibj) weightb[ibj+1] = 1;	
            else weightb[ibj] = 0;
          }
          if ((int) btagjets.size ()>=3) weightb[4] = 1;
	}
      }
     }
    }
  }

  //Compute a collection of jets and b-jets and a weight even if the selection of pair fails
  if(!step_pair){
      vector < NTJet > SelectedJets = GetSelectedJets (GetSelectedMuons(), GetSelectedElectrons(), applyJES, JESParam);
      // compute the weightb for b-tag
      if (GetFlagb() == 1) {
	if (!dataset->isData ()) {	//MC
	  vector < float >weight_temp = GetSFBweight().GetWeigth4BSel (GetMethodb(), GetSystb(), SelectedJets);
	  weightb[0] *= weight_temp[0];	//weight of the event
	  weightb[1] = weight_temp[1];	//proba 0 jet
	  weightb[2] = weight_temp[2];	//proba 1 jet
	  weightb[3] = weight_temp[3];	//proba 2 jets
	  weightb[4] = weight_temp[4];	//proba at least 3 jets
	}
      }
      jetsAna = SelectedJets;
      vector < NTJet > btagjets;
      vector < float >btagDiscri;
      for (unsigned int j = 0; j < SelectedJets.size (); j++) {
	if (passBtagSelection(SelectedJets[j])) {
	  btagjets.push_back (SelectedJets[j]);
	  btagDiscri.push_back (sfb_dilep_.getBtagDiscr(SelectedJets[j], btagAlgo_));
	}
      }
      bjetsAna = btagjets;
      if (GetFlagb() == 1) {
	if (dataset->isData ()) {	//Data
          if ((int) btagjets.size () >= NofBtagJets_) {
              weightb[0] *= 1;
          }
          else  {
             weightb[0]=0;
          }
          for (int ibj=0; ibj<4; ibj++) {
	    if ((int) btagjets.size ()== ibj) weightb[ibj+1] = 1;	
            else weightb[ibj] = 0;
          }
          if ((int) btagjets.size ()>=3) weightb[4] = 1;
	}
      }
  }


  //compute FinalStep
  if (step_trigger) {
    FinalStep++;
    if (step_pair) {
      FinalStep++;
      if (step_Zveto) {
	FinalStep++;
	if (step_jets) {
	  FinalStep++;
	  if (step_met) {
	    FinalStep++;
	    if (step_bjets) {
	      FinalStep++;
	      //the event is selected
	      if (print)
		cout << dump << endl;
	    }
	  }
	}
      }
    }
  }
  selCode_ = step_trigger*maskTriggerCut + step_pair*maskPairCut + 
	step_Zveto * maskZvetoCut + step_jets * maskJetCut + 
	step_met * maskMETCut + step_bjets * maskBjetsCut;

  return FinalStep;
}


int DiLeptonSelection::FillTable (SelectionTable & selTable, Dataset * dataset, int idataset, float weight)
{

  int sel = doFullSelection (dataset, selTable.Channel (), false);	// true-> has to be modified !!
  for (unsigned int i = 0; i < cuts_.size () + 1; i++)
    if (sel >= (int) i)
      selTable.Fill (idataset, i, weight);
  return sel;
}

     
int DiLeptonSelection::FillTableForMM(float iso1_in, float iso2_in, SelectionTable& selTable, Dataset* dataset, int idataset, float weight){
      
  int sel = doFullSelectionForMM(iso1_in, iso2_in, dataset, selTable.Channel(), false); // true-> has to be modified !!
  /*
        for(unsigned int i=0;i<cuts_.size()+1;i++)
                if(sel >=(int)i ) selTable.Fill(idataset,i, weight);
      
  */
    	return sel;
      
}
  



//******************************************************
//******** Dilepton trigger selection ******************
//******************************************************
//bool DiLeptonSelection::passTriggerSelection(string datasetName, string channelName){
bool DiLeptonSelection::passTriggerSelection (Dataset * dataset, string channelName)
{

  //bool match_HLT_Ele10_LW_L1R_recoEl = eleHLTMatch;
//cout << "start TS\n";
  ////cout << eleHLTMatch << endl;
  bool passEl = false;
  bool passMu = false;
  bool passElMu = false;
  int skim = -1;

  string datasetName = dataset->Name ();
//   if (datasetName == "MuData")
//     skim = 0;
//   if (datasetName == "EGData")
//     skim = 1;
// to be compatible with MyCutFlow and others
  if (datasetName.compare(0,6,"DataMu")==0)
    skim = 0;
  if (datasetName.compare(0,6,"DataEG")==0)
    skim = 1;
//  if (datasetName == "DataMuEG" || datasetName == "DataEGMu")
  if ( datasetName.compare(0,8,"DataMuEG")==0 || datasetName.compare(0,8,"DataEGMu")==0 )
    skim = 2;
  
  
  //cout << " datasetName " << datasetName << endl;
  
  const NTTrigger* triggers = GetPointer2Trigger();

  if (!dataset->isData ()) {	//MC
    //cout << " test trigger list --------------------------------" << endl;
    
    
passMu = (triggers->IsFired("HLT_DoubleMu6_v1") // Summer11
               || triggers->IsFired("HLT_DoubleMu6_v8")
              || triggers->IsFired("HLT_DoubleMu7_v12")); // Fall11

passEl = (triggers->IsFired("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2") // Summer11
               || triggers->IsFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8")
               || triggers->IsFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10"));

passElMu = (triggers->IsFired("HLT_Mu8_Ele17_CaloIdL_v2") || triggers->IsFired("HLT_Mu10_Ele10_CaloIdL_v3") // Summer11
                || triggers->IsFired("HLT_Mu8_Ele17_CaloIdL_v9") ||
               triggers->IsFired("HLT_Mu8_Ele17_CaloIdL_v13")); // Fall11 
    
    

  }
  else {
    // DATA --> Taken from TopDileptonRefAnalysis2010Pass6
    // WARNING : COULD SOMEONE CHECK THE < AND <= FOR THE RUN NUMBER????

    int runNumber = getRunNumber();

// Caro on basis of Ana's mail 10 october 2011A
      // a priori : no change since previous version for mumu, but change for ee and emu!


///////////////////////double muon channel  ///////////////////////

    if ( runNumber <= 165208 )
      passMu = (triggers->IsFired("HLT_DoubleMu7_v1") || triggers->IsFired("HLT_DoubleMu7_v2"));          

    if (  runNumber  >= 165209 && runNumber  <=178419)
      passMu = (    triggers->IsFired("HLT_Mu13_Mu8_v2") 
	         || triggers->IsFired("HLT_Mu13_Mu8_v3")
	         || triggers->IsFired("HLT_Mu13_Mu8_v4")
	         || triggers->IsFired("HLT_Mu13_Mu8_v5")
	         || triggers->IsFired("HLT_Mu13_Mu8_v6")
	         || triggers->IsFired("HLT_Mu13_Mu8_v7")   );

     if (  runNumber  >= 178420 && runNumber  <=999999 )
       passMu = (triggers->IsFired("HLT_Mu17_Mu8_v10") || triggers->IsFired("HLT_Mu17_Mu8_v11"));          

     if (  runNumber  >= 178420 && runNumber  <=999999)
       passMu = (triggers->IsFired("HLT_Mu17_TkMu8_v3") || triggers->IsFired("HLT_Mu17_TkMu8_v4"));          


///////////////////////double electron channel  ///////////////////////

    if ( runNumber <= 170901 )
      passEl = (    triggers->IsFired("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1") 
	         || triggers->IsFired("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2")
	         || triggers->IsFired("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3")
	         || triggers->IsFired("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4")
	         || triggers->IsFired("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5")
	         || triggers->IsFired("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6")
	         || triggers->IsFired("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v7")
	         || triggers->IsFired("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v8")   );

    else if ( runNumber>=  170902 )   
      passEl = (    triggers->IsFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6") 
	         || triggers->IsFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7")
	         || triggers->IsFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8")
	         || triggers->IsFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9")
	         || triggers->IsFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10")   );


///////////////////////e mu channel  ///////////////////////

    if ( runNumber  <=175972 ) // boundary of 2e33/3e33 menu
      passElMu = (   triggers->IsFired("HLT_Mu17_Ele8_CaloIdL_v1") 
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdL_v2")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdL_v3")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdL_v4")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdL_v5")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdL_v6")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdL_v7")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdL_v8")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdL_v9")   );


    if ( runNumber  >= 175973 )

      passElMu = (   triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4") 
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8")   );

    if (  runNumber  <= 167913 )
      passElMu = (   triggers->IsFired("HLT_Mu8_Ele17_CaloIdL_v1") 
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdL_v2")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdL_v3")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdL_v4")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdL_v5")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdL_v6")   );

    if ( runNumber  >=  167914 )
      passElMu = (   triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3") 
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8")   );
		  
    if ( runNumber  >=  190000 )
      passElMu = (   triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1") 
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9")
                  || triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10")
		   );	  
    if ( runNumber  >=  190000 )
      passElMu = (   triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1") 
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9")
                  || triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10")
		   );
       
  } // end DATA
  

  if (channelName == string ("ee") && ( (skim==-1) || (skim==1) ) ) {
    if (passEl)
      return true;
    else
      return false;
  }
  if (channelName == string ("mumu") && ( (skim==-1) || (skim==0) ) ) {
    if (passMu)
      return true;
    else
      return false;
  }
  if (channelName == string ("emu")) {
     bool thereturn = false;
//     if (  (passEl && skim==1 ) || ( passMu && skim==0 && !passEl ) ) thereturn = true;
//     if ( skim == -1 && (passEl  || passMu ) ) thereturn = true;
    if ( skim == -1 &&  passElMu) return true;
    if ( skim == 2 && (passElMu) )
      thereturn = true;
    return thereturn;
  }
  if (channelName == string ("") || channelName == string ("*") || channelName == string ("allChannels")) {
    if (passEl || passMu || passElMu)
      return true;
    return false;
  }

  return false;

}



//******************************************************
//******** Dilepton trigger selection ******************
//******************************************************
//bool DiLeptonSelection::passTriggerSelection(string datasetName, string channelName){
bool DiLeptonSelection::passTriggerSelection8TeV (Dataset * dataset, string channelName)
{


  bool passEl = false;
  bool passMu = false;
  bool passElMu = false;
 
  string datasetName = dataset->Name ();
  
  
  //cout << " datasetName " << datasetName << endl;
  
  const NTTrigger* triggers = GetPointer2Trigger();

  if (!dataset->isData ()) {	//MC
    //cout << " test trigger list --------------------------------" << endl;
    
    
    passMu = (triggers->IsFired("HLT_Mu17_Mu8_v17") || triggers->IsFired("HLT_Mu17_TkMu8_v10") ); // Summer12 53X

    passEl = (triggers->IsFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17")  );

    passElMu = (triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7") || triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7") );  
    
    

  }
  else {
    // DATA --> Taken from TopDileptonRefAnalysis2010Pass6
    // WARNING : COULD SOMEONE CHECK THE < AND <= FOR THE RUN NUMBER????

    //    int runNumber = getRunNumber();

//  
   /* 
   std::vector<IPHCTree::NTTriggerPathType> subtable_mu;
    
   triggers->GetSubTable("HLT_Mu17_Mu8_v*", subtable_mu);
    
   std::vector<IPHCTree::NTTriggerPathType> subtable_el;
    
   triggers->GetSubTable("HLT_Ele*", subtable_el);
   cout << "******************************" << endl;
   for(unsigned int i=0; i<subtable_mu.size(); i++){
    
      cout << subtable_mu[i].name << endl;
    
   }
   for(unsigned int i=0; i<subtable_el.size(); i++){
    
      cout << subtable_el[i].name << endl;
    
   }
   */
    
    int runNumber = getRunNumber();
    //if(runNumber >= 190456 && runNumber <= 193805){
      
      
      std::vector<IPHCTree::NTTriggerPathType> subtable_dimu1, subtable_dimu2 ;
      
      
      triggers->GetSubTable("HLT_Mu17_Mu8_v*", subtable_dimu1);
      triggers->GetSubTable("HLT_Mu17_TkMu8_v*", subtable_dimu2);
      for(unsigned int i=0; i<subtable_dimu1.size(); i++){
        if(subtable_dimu1[i].fired == 1) passMu = true;
      }
      for(unsigned int i=0; i<subtable_dimu2.size(); i++){
        if(subtable_dimu2[i].fired == 1) passMu = true;
      }
      
      
      
      std::vector<IPHCTree::NTTriggerPathType> subtable_diel ;
      triggers->GetSubTable("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", subtable_diel);
      for(unsigned int i=0; i< subtable_diel.size(); i++){
        if(subtable_diel[i].fired == 1)  passEl= true;
      }
    
      
      
      std::vector<IPHCTree::NTTriggerPathType> subtable_elmu1, subtable_elmu2 ;
      
      triggers->GetSubTable("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", subtable_elmu1);
      triggers->GetSubTable("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", subtable_elmu2);
      for(unsigned int i=0; i<subtable_elmu1.size(); i++){
        if(subtable_elmu1[i].fired == 1)  passElMu = true;
      }
		
		
    /*
    }else if(runNumber >= 193806 && runNumber<=198021) {
    
      
      passMu = (triggers->IsFired("HLT_Mu17_Mu8_v*") || triggers->IsFired("HLT_Mu17_TkMu8_v*"));
    
      
      passEl = (triggers->IsFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"));
    
      
      passElMu = (triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*") || 
    	        triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"));
    }else if(runNumber >= 198022 ) {
    
      
      passMu = (triggers->IsFired("HLT_Mu17_Mu8_v21") || triggers->IsFired("HLT_Mu17_TkMu8_v13"));
    
      
      passEl = (triggers->IsFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"));
    
      
      passElMu = (triggers->IsFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*") || 
    	        triggers->IsFired("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"));
      
    }*/
       
  } // end DATA
  

  if (channelName == string ("ee") ){
   if( passEl ) return true;
    else return false;
  }
  if (channelName == string ("mumu") ){
    if( passMu ) return true;
    else return false;
  }
  
  if (channelName == string ("emu") ){
    if( passElMu ) return true;
    else return false;
  }

}

//******************************************************
//************* MET trigger selection ******************
//******************************************************
//bool DiLeptonSelection::passTriggerSelection(string datasetName, string channelName){
std::pair<bool,bool> DiLeptonSelection::passMETTriggerSelection (Dataset * dataset)
{
  
 
  bool passTriggerBTAGMET = 0;
  bool passTriggerMET     = 0;
 
  const NTTrigger* triggers = GetPointer2Trigger();
 
  if (!dataset->isData ()) {	//MC
    //cout << " test trigger list --------------------------------" << endl;
    
    passTriggerMET = (    triggers->IsFired("HLT_CentralJet80_MET65_v1") // Summer 11
                      || triggers->IsFired("HLT_CentralJet80_MET80_v1")
                      || triggers->IsFired("HLT_CentralJet80_MET100_v1")
                      || triggers->IsFired("HLT_CentralJet80_MET160_v1")
                      || triggers->IsFired("HLT_DiJet60_MET45_v1")
                      || triggers->IsFired("HLT_PFMHT150_v2")
                      || triggers->IsFired("HLT_MET100_v1")
                      || triggers->IsFired("HLT_MET120_v1")
                      || triggers->IsFired("HLT_MET200_v1")

                      || triggers->IsFired("HLT_CentralJet80_MET65_v7") // Fall 11
                      || triggers->IsFired("HLT_CentralJet80_MET80_v6")
                      || triggers->IsFired("HLT_CentralJet80_MET100_v7")
                      || triggers->IsFired("HLT_CentralJet80_MET160_v7")
                      || triggers->IsFired("HLT_DiJet60_MET45_v7")
                      || triggers->IsFired("HLT_PFMHT150_v12")
                      || triggers->IsFired("HLT_MET65_v7")
                      || triggers->IsFired("HLT_MET100_v7")
                      || triggers->IsFired("HLT_MET120_v7")
                      || triggers->IsFired("HLT_MET200_v7")
                      || triggers->IsFired("HLT_MET400_v7")   );
   
    if( passTriggerMET ) 
      passTriggerBTAGMET = (   triggers->IsFired("HLT_BTagMu_DiJet20_Mu5_v2") // Summer 11
                	    || triggers->IsFired("HLT_BTagMu_DiJet60_Mu7_v2")
                	    || triggers->IsFired("HLT_BTagMu_DiJet80_Mu9_v2")
                	    || triggers->IsFired("HLT_BTagMu_DiJet100_Mu9_v2")

                            || triggers->IsFired("HLT_BTagMu_DiJet20_Mu5_v10") // Fall 11
                	    || triggers->IsFired("HLT_BTagMu_DiJet60_Mu5_v10")
                	    || triggers->IsFired("HLT_BTagMu_DiJet80_Mu5_v10")
                	    || triggers->IsFired("HLT_BTagMu_DiJet100_Mu5_v10")   );
     
   }
  else {
    //always return true is data. To be used on MET/BTagMET SD only
     passTriggerBTAGMET = 1; passTriggerMET = 1;
 }
 
 
 std::pair<bool, bool> theReturn;
 
 theReturn.first  = passTriggerBTAGMET;
 theReturn.second = passTriggerMET;
 return theReturn;
}


float DiLeptonSelection::GetMinValueMassCut ()
{
  return MinMassCut_;
}

pair < float, float >DiLeptonSelection::GetMETCut ()
{
  return METCuts_;
}

pair < float, float >DiLeptonSelection::GetZmassWindowCut ()
{
  return ZMassWindow_;
}

     int DiLeptonSelection::GetbtagAlgo () const
     {
       return btagAlgo_;
     }

     float DiLeptonSelection::GetbtagDiscriCut () const
     {
       return btagDiscriCut_;
     }

     int DiLeptonSelection::GetNofBtagJetsCut () const
     {
       return NofBtagJets_;
     }

     int DiLeptonSelection::FillTablewBweight (SelectionTable & selTable, Dataset * dataset, int idataset, float weight, int method_b, vector < float >&weightb)
{

  ///// USEFUL IF weightb COMPUTED OUTSIDE doFullSelection
  ///// OTHERWISE USE THE OTHER FillTablewBweight, SEE BELOW!

  //int sel = doFullSelection (dataset, selTable.Channel (), true);
  int sel = doFullSelection (dataset, selTable.Channel (), false);
  // For cuts up to b-tag:
  for (unsigned int i = 0; i < cuts_.size () - 1; i++) {
    if (sel >= (int) i)
      selTable.Fill (idataset, i, weight);
    //        std::cout << "i " << i << " sel " << sel << " for cut " << cuts_[i] << std::endl;
  }
  // For btag cut : cuts_.size()==8
  int cut_num_btag = (int) cuts_.size () - 1;
  if (method_b == 0) {
    // Discri cut + SF
    if (sel >= cut_num_btag)
      selTable.Fill (idataset, cut_num_btag, weightb[0]);
    //      std::cout << "cut_num_btag0 " << cut_num_btag << " sel " << sel << " for cut " << cuts_[cut_num_btag] << std::endl;
  }
  else if (method_b >= 1) {
    // No discri cut + Eff
    if (sel >= cut_num_btag - 1)
      selTable.Fill (idataset, cut_num_btag, weightb[0]);
    //      std::cout << "cut_num_btag1 " << cut_num_btag << " sel " << sel << " for cut " << cuts_[cut_num_btag] << std::endl;
  }

  return sel;
}


int DiLeptonSelection::FillTablewBweight (SelectionTable & selTable, Dataset * dataset, int idataset, float weight, vector < float >&weightb)
{
  int sel = 0;
  // need to initialize weightb each time we inter in FillTablewBweight
  weightb[0] = weight;		// weight of the event
  weightb[1] = 0.;		// Proba of 0 jet
  weightb[2] = 0.;		// Proba of 1 jet;
  weightb[3] = 0.;		// Proba of 2 jets;
  weightb[4] = 0.;		// Proba of at least 3 jets;
  //cout << " go through FillTablewBweight " << endl;
  if (GetFlagb() == 0) {
    sel = FillTable (selTable, dataset, idataset, weight);
  }
  else if (GetFlagb() == 1) {
    sel = doFullSelection(dataset, weightb, selTable.Channel(), false);  // not all the possible input variables defined here!
    // For cuts up to b-tag:
    for (unsigned int i = 0; i < cuts_.size () - 1; i++) {
      if (sel >= (int) i)
	selTable.Fill (idataset, i, weight);
      //        std::cout << "i " << i << " sel " << sel << " for cut " << cuts_[i] << std::endl;
    }
    // For btag cut : cuts_.size()==8
    int cut_num_btag = (int) cuts_.size () - 1;
    if (GetMethodb() == 0) {
      // Discri cut + SF
      if (sel >= cut_num_btag)
	selTable.Fill (idataset, cut_num_btag, weightb[0]);
      //      std::cout << "cut_num_btag0 " << cut_num_btag << " sel " << sel << " for cut " << cuts_[cut_num_btag] << std::endl;
    }
    else if (GetMethodb() >= 1) {
      // No discri cut + Eff
      if (sel >= cut_num_btag - 1)
	selTable.Fill (idataset, cut_num_btag, weightb[0]);
      //      std::cout << "cut_num_btag1 " << cut_num_btag << " sel " << sel << " for cut " << cuts_[cut_num_btag] << std::endl;
    }
  }

  return sel;
}

bool DiLeptonSelection::passBtagSelection(const NTJet & jet) const
{
  return Selection::passBtagSelection(jet, btagAlgo_, btagDiscriCut_);
}


double DiLeptonSelection::getLeptonScaleFactor(double pt, double eta, string lepton){

  double the_getScaleFactor = 0;

  if(pt > 100) pt = 99;
  if(fabs(eta) > 2.5) eta = 2.4;

  if(lepton == "e"){
     int binx = getScaleFactEl()->GetXaxis()->FindBin( pt );
     int biny = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta) );
     the_getScaleFactor = getScaleFactEl()->GetBinContent( binx, biny );
  }
  if(lepton == "mu"){
     int binx = getScaleFactMu()->GetXaxis()->FindBin( pt );
     int biny = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta) );
     the_getScaleFactor = getScaleFactMu()->GetBinContent( binx, biny );
  }

  return the_getScaleFactor;
}

double DiLeptonSelection::getLeptonScaleFactorError(double pt, double eta, string lepton){

  double the_getScaleFactor = 0;

  if(pt > 100) pt = 99;
  if(fabs(eta) > 2.5) eta = 2.4;

  if(lepton == "e"){
     int binx = getScaleFactEl()->GetXaxis()->FindBin( pt );
     int biny = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta) );

     the_getScaleFactor = getScaleFactEl()->GetBinError( binx, biny );
  }
  if(lepton == "mu"){
     int binx = getScaleFactMu()->GetXaxis()->FindBin( pt );
     int biny = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta) );

     the_getScaleFactor = getScaleFactMu()->GetBinError( binx, biny );
  }

  return the_getScaleFactor;
}

double DiLeptonSelection::getLeptonScaleFactor(double pt1, double eta1, double pt2, double eta2, string channel){
  double the_getScaleFactor = 0;
  
  if(pt1 > 100) pt1 = 99;
  if(fabs(eta1) > 2.5) eta1 = 2.4;
  
  if(pt2 > 100) pt2 = 99;
  if(fabs(eta2) > 2.5) eta2 = 2.4;
//  cout << " debug " << pt1 << " " << eta1 << " " << pt2 << " " << eta2 << endl;
  if(channel == "ee"){
     int binx1 = getScaleFactEl()->GetXaxis()->FindBin( pt1 );
     int biny1 = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta1) );
     
     int binx2 = getScaleFactEl()->GetXaxis()->FindBin( pt2 );
     int biny2 = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta2) );
     
     the_getScaleFactor = getScaleFactEl()->GetBinContent( binx1, biny1 )*getScaleFactEl()->GetBinContent( binx2, biny2 );

  }
  
  if(channel == "mumu"){
     int binx1 = getScaleFactMu()->GetXaxis()->FindBin( pt1 );
     int biny1 = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta1) );
     int binx2 = getScaleFactMu()->GetXaxis()->FindBin( pt2 );
     int biny2 = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta2) );
     
     the_getScaleFactor = getScaleFactMu()->GetBinContent( binx1, biny1 )*getScaleFactMu()->GetBinContent( binx2, biny2 );
  }

  if(channel == "emu"){
  
     int binx1 = getScaleFactEl()->GetXaxis()->FindBin( pt1 );
     int biny1 = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta1) );
     
     int binx2 = getScaleFactMu()->GetXaxis()->FindBin( pt2 );
     int biny2 = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta2) );
     
     the_getScaleFactor = getScaleFactEl()->GetBinContent( binx1, biny1 )*getScaleFactMu()->GetBinContent( binx2, biny2 );
  }
  
//  cout << " debug the_getScaleFactor " << the_getScaleFactor << endl;
  
  return the_getScaleFactor;
}





double DiLeptonSelection::getLeptonScaleFactorError(double pt1, double eta1, double pt2, double eta2, string channel){
  double the_getScaleFactor = 0;
  
  if(pt1 > 100) pt1 = 99;
  if(fabs(eta1) > 2.5) eta1 = 2.4;
  
  if(pt2 > 100) pt2 = 99;
  if(fabs(eta2) > 2.5) eta2 = 2.4;
  
  if(channel == "ee"){
     int binx1 = getScaleFactEl()->GetXaxis()->FindBin( pt1 );
     int biny1 = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta1) );
     
     int binx2 = getScaleFactEl()->GetXaxis()->FindBin( pt2 );
     int biny2 = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta2) );
     
     the_getScaleFactor = getScaleFactEl()->GetBinError( binx1, biny1 )*getScaleFactEl()->GetBinContent( binx2, biny2 );

  }
  
  if(channel == "mumu"){
  
     int binx1 = getScaleFactMu()->GetXaxis()->FindBin( pt1 );
     int biny1 = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta1) );
     
     int binx2 = getScaleFactMu()->GetXaxis()->FindBin( pt2 );
     int biny2 = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta2) );
     
     the_getScaleFactor = getScaleFactMu()->GetBinError( binx1, biny1 )*getScaleFactMu()->GetBinContent( binx2, biny2 );
  }

  if(channel == "emu"){
  
     int binx1 = getScaleFactEl()->GetXaxis()->FindBin( pt1 );
     int biny1 = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta1) );
     
     int binx2 = getScaleFactMu()->GetXaxis()->FindBin( pt2 );
     int biny2 = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta2) );
     
     the_getScaleFactor = getScaleFactEl()->GetBinError( binx1, biny1 )*getScaleFactMu()->GetBinContent( binx2, biny2 );
  }
  
  
  return the_getScaleFactor;
}












