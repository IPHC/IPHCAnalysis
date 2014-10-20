#ifndef Selection_h
#define Selection_h

// IPHC headers
#include "NTFormat/interface/NTEvent.h"
#include "Selection/interface/Event.h"
#include "Selection/interface/Requirement.h"
#include "EffEstimation/interface/SFBweight.h"

//#include "../../../../JR_Standalone/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "../../../../JR_Standalone/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetCorrections/interface/JetCorrectionUncertainty.h"
#include "JetCorrections/interface/JetCorrectorParameters.h"
#include "JetCorrections/interface/JetCorrectionUncertainty.h"



// ROOT headers
#include "TH2F.h"
#include "TGraphAsymmErrors.h"

// STL headers
#include <memory>
#include <vector>


struct HighestPt
{
  bool operator()( NTJet j1, NTJet j2 ) const
  {
    return j1.p4.Pt() > j2.p4.Pt() ;
  }
  bool operator()( NTLepton j1, NTLepton j2 ) const
  {
    return j1.p4.Pt() > j2.p4.Pt() ;
  }
};

void InitVectorOfWeight(vector<float>& weightb);

void LoadBWeight(Selection& sel, vector<float>& weightb, const vector<NTJet>& theselJets);	


//! \class Selection
//! History for generated event
class Selection : public Event
{

  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 public:

  //! Constructor without argument
  Selection();

  //! Copy constructor
  Selection(const Selection &);

  //! Destructor
  ~Selection()
  {  
    delete histo_jesunc_;
    delete scaleFactEl;
    delete scaleFactMu;
  }


  // -------------- accessor to vertex collections -------------

  //! Get selected vertices
  std::vector<IPHCTree::NTVertex> GetSelectedVertex() const; 

  // -------------- accessor to electron collections -------------

  //! Get scaled electrons
  std::vector<IPHCTree::NTElectron> GetScaledElectrons(float scale = 1.) const;

  //! Get smeared electrons
  std::vector<IPHCTree::NTElectron> GetSmearedElectrons(float resol = 1.)const;

  //! Get selected electrons
  std::vector<IPHCTree::NTElectron> GetSelectedElectrons(
                          bool applyLES = false, float scale = 1.,
                          bool applyLER = false, float resol = 1.) const;

  //! Get selected electrons
  std::vector<IPHCTree::NTElectron> GetSelectedElectrons(
                          float PtThr, float EtaThr,
                          float ElectronRelIso, bool applyLES = false,
                          float scale = 1., bool applyLER = false,
                          float resol = 1.) const;   

  //! Get selected electrons no iso
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsNoIso(
                          bool applyLES = false, float scale = 1.,
                          bool applyLER = false, float resol = 1.) const;

  //! Get selected electrons no iso
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsNoIso(
                          float PtThr, float EtaThr, 
                          bool applyLES = false, float scale = 1., 
                          bool applyLER = false, float resol = 1.) const;

  //! Get selected electrons iso non id
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsIsoNonID(
                          bool applyLES = false, float scale = 1.,
                          bool applyLER = false, float resol = 1.) const;

  //! Get selected electrons non iso non id
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsNoIsoNonID(
                          bool applyLES = false, float scale = 1.,
                          bool applyLER = false, float resol = 1.) const;
			  
			  
  // -------------- accessor to electron collections for ttdilep -------------
  
    //! Get scaled electrons
  std::vector<IPHCTree::NTElectron> GetScaledElectronsDileptonTTbar(float scale = 1.) const;

  //! Get smeared electrons
  std::vector<IPHCTree::NTElectron> GetSmearedElectronsDileptonTTbar(float resol = 1.)const;

  //! Get selected electrons
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsDileptonTTbar(
                          bool applyLES = false, float scale = 1.,
                          bool applyLER = false, float resol = 1., float rho=0.) const;
			  
  //! Get selected electrons
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsRhoIso(
                          bool applyLES = false, float scale = 1.,
                          bool applyLER = false, float resol = 1., float rho=0.) const;

  //! Get selected electrons
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsDileptonTTbar(
                          float PtThr, float EtaThr,
                          float ElectronRelIso, bool applyLES = false,
                          float scale = 1., bool applyLER = false,
                          float resol = 1., float rho=0) const;   

  //! Get selected electrons
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsRhoIso(
                          float PtThr, float EtaThr,
                          float ElectronRelIso, bool applyLES = false,
                          float scale = 1., bool applyLER = false,
                          float resol = 1., float rho=0) const;   

  //! Get selected electrons no iso
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsNoIsoDileptonTTbar(
                          bool applyLES = false, float scale = 1.,
                          bool applyLER = false, float resol = 1.) const;

  //! Get selected electrons no iso
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsNoIsoDileptonTTbar(
                          float PtThr, float EtaThr, 
                          bool applyLES = false, float scale = 1., 
                          bool applyLER = false, float resol = 1.) const;

  //! Get selected electrons iso non id
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsIsoNonIDDileptonTTbar(
                          bool applyLES = false, float scale = 1.,
                          bool applyLER = false, float resol = 1.) const;

  //! Get selected electrons non iso non id
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsNoIsoNonIDDileptonTTbar(
                          bool applyLES = false, float scale = 1.,
                          bool applyLER = false, float resol = 1.) const;
			  
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsLooseDileptonTTbar(float rho)  const;
  
  // -------------- accessor to muon collections -------------

  //! Get scaled muons
  std::vector<IPHCTree::NTMuon> GetScaledMuons(float scale = 1.) const;

  //! Get selected muons no iso
  std::vector<IPHCTree::NTMuon> GetSelectedMuonsNoIso(
                                           bool applyLES = false,
                                           float scale = 1.) const;

  //! Get selected muons
  std::vector<IPHCTree::NTMuon> GetSelectedMuons(
                                           bool applyLES = false,
                                           float scale = 1.) const;

  //! Get selected muons no Iso
  std::vector<IPHCTree::NTMuon> GetSelectedMuonsNoIso(
                                           float PtThr, float EtaThr,
                                           bool applyLES = false, 
                                           float scale = 1.) const;

  //! Get selected muons
  std::vector<IPHCTree::NTMuon> GetSelectedMuons(
                                           float PtThr, float EtaThr,
                                           float MuonRelIso, 
                                           bool applyLES = false,
                                           float scale = 1.) const;

  //! Get selected muons
  std::vector<IPHCTree::NTMuon> GetSelectedMuonIsoNonID(
                                           bool applyLES = false,
                                           float scale = 1.) const;




  // -------------- accessor to muon collections for dilepton -------------

  //! Get scaled muons
  std::vector<IPHCTree::NTMuon> GetScaledMuonsDileptonTTbar(float scale = 1.) const;

  //! Get selected muons no iso
  std::vector<IPHCTree::NTMuon> GetSelectedMuonsNoIsoDileptonTTbar(
                                           bool applyLES = false,
                                           float scale = 1.) const;

  //! Get selected muons
  std::vector<IPHCTree::NTMuon> GetSelectedMuonsDileptonTTbar(
                                           bool applyLES = false,
                                           float scale = 1.) const;
					   
  //! Get selected muons
  std::vector<IPHCTree::NTMuon> GetSelectedMuonsDeltaBetaIso(
                                           bool applyLES = false,
                                           float scale = 1.) const;

  //! Get selected muons no Iso
  std::vector<IPHCTree::NTMuon> GetSelectedMuonsNoIsoDileptonTTbar(
                                           float PtThr, float EtaThr,
                                           bool applyLES = false, 
                                           float scale = 1.) const;

  //! Get selected muons
  std::vector<IPHCTree::NTMuon> GetSelectedMuonsDileptonTTbar(
                                           float PtThr, float EtaThr,
                                           float MuonRelIso, 
                                           bool applyLES = false,
                                           float scale = 1.) const;

  //! Get selected muons
  std::vector<IPHCTree::NTMuon> GetSelectedMuonsDeltaBetaIso(
                                           float PtThr, float EtaThr,
                                           float MuonRelIso, 
                                           bool applyLES = false,
                                           float scale = 1.) const;

  //! Get selected muons
  std::vector<IPHCTree::NTMuon> GetSelectedMuonIsoNonIDDileptonTTbar(
                                           bool applyLES = false,
                                           float scale = 1.) const;









  // -------- accessor to lepton + jet collections ------------
      
  //! GetSelected muons for l+jets
  std::vector<IPHCTree::NTMuon> GetSelectedMuonsForLJets(
                                    const std::vector<IPHCTree::NTJet>& SelectedJets,
                                    bool applyLES = false, 
                                    float scale = 1.) const;

  std::vector<IPHCTree::NTMuon> GetSelectedMuonsNoIsoForLJets(
                                    const std::vector<IPHCTree::NTJet>& SelectedJets,
                               float PtThr, float EtaThr,
                               bool applyLES = false, float scale = 1.) const;

  std::vector<IPHCTree::NTMuon> GetSelectedMuonsForLJets(
                                    const std::vector<IPHCTree::NTJet>& SelectedJets,
                               float PtThr, float EtaThr, 
                               float MuonRelIso, bool applyLES = false,
                               float scale = 1.) const;

  std::vector<IPHCTree::NTMuon> GetVetoMuonsForLJets(
                              bool applyLES = false, float scale = 1.) const;

  //! GetSelected electrons for l+jets
  std::vector<IPHCTree::NTElectron> GetSelectedElectronsForLJets(
                                    const std::vector<IPHCTree::NTJet>& SelectedJets,
                                bool applyLES = false, float scale = 1.,
                                bool applyLER = false, float resol = 1.) const;

  std::vector<IPHCTree::NTElectron> GetSelectedElectronsNoIsoForLJets(
                                    const std::vector<IPHCTree::NTJet>& SelectedJets,
                                           float PtThr, float EtaThr,
                                           bool applyLES = false, float scale = 1.,
                                           bool applyLER = false, float resol = 1.) const;

  std::vector<IPHCTree::NTElectron> GetSelectedElectronsForLJets(
                                    const std::vector<IPHCTree::NTJet>& SelectedJets,
                               float PtThr, float EtaThr,
                               float ElectronRelIso,
                               bool applyLES = false, float scale = 1.,
                               bool applyLER = false, float resol = 1.) const;

 std::vector<IPHCTree::NTElectron> GetVetoElectronsForLJets(
                               bool applyLES = false, float scale = 1.,
                               bool applyLER = false, float resol = 1.) const;




  //! GetSelected loose muons for jets
  std::vector<IPHCTree::NTMuon> GetSelectedLooseMuonsForMuJets(
                               bool applyLES = false, float scale = 1.) const;

  //! GetSelected loose electrons for jets
  std::vector<IPHCTree::NTElectron> GetSelectedLooseElectronsForMuJets(
                               bool applyLES = false,
                               float scale = 1.) const;


  // -------------- accessor to jet collections -------------

  //! Get clean Jets
  std::vector<IPHCTree::NTJet> GetCleanJets(
                         std::vector<IPHCTree::NTJet> jet_cand, 
			 std::vector<IPHCTree::NTMuon> muon_cand, 
			 std::vector<IPHCTree::NTElectron> elec_cand) const;

  //! Get scaled jets
  std::vector<IPHCTree::NTJet> GetScaledJets(float factor = 1.) const; 

  //! Get scaled jets
  std::vector<IPHCTree::NTJet> GetScaledJets( JetCorrectionUncertainty* theJESuncertainty, bool upOrDown , int jetResFactor = 0 ) const; 
  
  
  //! Get smeared jets
  std::vector<IPHCTree::NTJet> GetSmearedJets(
                          const vector<IPHCTree::NTJet>& injets,
                          int jetResFactor = 0) const; 

  //! Get selected jets
  std::vector<IPHCTree::NTJet> GetSelectedJets(bool applyJES = false, 
                                               float scale = 1., 
                                               bool applyJER = false,
                                               int ResFactor = 0) const;

  //! Get selected jets
  std::vector<IPHCTree::NTJet> GetSelectedJets(
                          const std::vector<IPHCTree::NTMuon>& muon_cand,
                          const std::vector<IPHCTree::NTElectron>& elec_cand,
                          bool applyJES = false, float scale = 1.,
                          bool applyJER = false, int ResFactor = 0) const;



  //! Get selected jets for 8TeV uncertainties
  std::vector<IPHCTree::NTJet> GetSelectedJets(
                          const std::vector<IPHCTree::NTMuon>& muon_cand,
                          const std::vector<IPHCTree::NTElectron>& elec_cand,
                          bool applyJES = false, JetCorrectionUncertainty* theJESuncertainty = 0, bool upOrDown = true, 
                          bool applyJER = false, int ResFactor = 0) const;
			  
			  

  //! Get selected jets for 8TeV uncertainties
  std::vector<IPHCTree::NTJet> GetSelectedJetsLoose(
                          const std::vector<IPHCTree::NTMuon>& muon_cand,
                          const std::vector<IPHCTree::NTElectron>& elec_cand,
                          bool applyJES = false, JetCorrectionUncertainty* theJESuncertainty = 0, bool upOrDown = true, 
                          bool applyJER = false, int ResFactor = 0) const;
			  
			  

  //! Get selected jets
  std::vector<IPHCTree::NTJet> GetSelectedJets(
                          const std::vector<IPHCTree::NTMuon>& muon_cand,
                          const std::vector<IPHCTree::NTElectron>& elec_cand,
                          const std::vector<IPHCTree::NTTau>& tau_cand,
                          bool applyJES = false, float scale = 1.,
                          bool applyJER = false, int ResFactor = 0) const;
 
  //! Get selected jets 
  std::vector<IPHCTree::NTJet> GetSelectedJets(
                          const std::vector<NTTau>& tau_cand, 
			  bool applyJES = false, float scale = 1., 
			  bool applyJER = false, int ResFactor = 0) const;
 
  //! Get selected jets
  std::vector<IPHCTree::NTJet> GetSelectedJets(
                          float PtThr, float EtaThr,
                          bool applyJES = false, float scale = 1.,
                          bool applyJER = false, int ResFactor = 0) const;

  //! Get selected jets For L+Jets
  std::vector<IPHCTree::NTJet> GetSelectedJetsForLJets(
                          bool applyJES = false, float scale = 1.,
                          bool applyJER = false, int ResFactor = 0) const;
			  
  std::vector<IPHCTree::NTJet> GetSelectedJetsForLJets(
                          const std::vector<IPHCTree::NTMuon>& muon_cand,
                          const std::vector<IPHCTree::NTElectron>& elec_cand,
                          bool applyJES = false, float scale = 1.,
                          bool applyJER = false, int ResFactor = 0) const;


  //! Get selected B jets
  std::vector<IPHCTree::NTJet> GetSelectedBJets(
                          const std::vector<IPHCTree::NTJet>& SelectedJets,
                          const int& algo, 
                          const float & discricut) const;


  // -------------- accessor to tau collections -------------

  //! Get scaled taus
  std::vector<IPHCTree::NTTau> GetScaledTaus(float scale = 1.) const;

  //! Get selected taus
  std::vector<IPHCTree::NTTau> GetSelectedTaus(
                          float PtThr, float EtaThr,
                          bool applyLES = false, float scale = 1.,
                          int isoLevel = 0., bool antiLep = true) const;

  //! Get selected taus
  std::vector<IPHCTree::NTTau> GetSelectedTaus(
                          bool applyLES = false, float scale = 1.,
                          int isoLevel = 0., bool antiLep = true) const;

  //! Get selected taus 
  std::vector<IPHCTree::NTTau> GetSelectedTaus(
                          const std::vector<IPHCTree::NTMuon>& muon_cand,
                          const std::vector<IPHCTree::NTElectron>& elec_cand,
                          float PtThr, float EtaThr,
                          bool applyLES = false , float scale = 1.,
                          int isoLevel = 0., bool antiLep = true) const;

  //! Get selected taus
  std::vector<IPHCTree::NTTau> GetSelectedTaus(
                          const std::vector<IPHCTree::NTMuon>& muon_cand,
                          const std::vector<IPHCTree::NTElectron>& elec_cand,
                          bool applyLES = false , float scale = 1.,
                          int isoLevel = 0., bool antiLep = true) const;   

  //! Get taus
  std::vector<IPHCTree::NTTau> GetTaus(
                     const std::vector<IPHCTree::NTMuon>& muon_cand,
                     const std::vector<IPHCTree::NTElectron>& elec_cand) const;

  // -------------- accessor to MET collections -------------

  //! Get scaled MET
  IPHCTree::NTMET GetScaledMET(float scale = 1.) const; 
  
  //! Get scaled MET
  IPHCTree::NTMET GetScaledMET(JetCorrectionUncertainty* theJESuncertainty, bool upOrDown) const; 

  //! Get smeared MET
  IPHCTree::NTMET GetSmearedMET(const std::vector<IPHCTree::NTJet>& injets,
                                int jetResFactor = 1) const;

  //! Get unclusmeared MET
  IPHCTree::NTMET GetUnclusScaledMET(
                              bool applyUnclusScale,
                              float scale) const;
  //! Get smeared MET
  IPHCTree::NTMET GetSelectedMET(bool applyJES = false, 
                         float scale = 1., 
                         bool applyJER = false,
                         int ResFactor = 0) const;
      
  //! Get smeared MET for 8TeV uncertainties
  IPHCTree::NTMET GetSelectedMET(bool applyJES = false, 
  		         JetCorrectionUncertainty* theJESuncertainty = 0, 
			 bool upOrDown = true,
                         bool applyJER = false,
                         int ResFactor = 0) const;
      
  //! tools for Type1 MET calculated ourselves, on the flight.  
  IPHCTree::NTMET GetScaledType1METWithJER(vector<IPHCTree::NTJet> injets, 
                                           bool applyJES = false, float scale = 1., 
					   bool applyJER = false, int ResFactor = 0) const;
					   
  IPHCTree::NTMET GetType1MET(vector<IPHCTree::NTJet> injets, 
                              bool applyJES = false, float scale = 1., 
			      bool applyJER = false, int ResFactor = 0) const;  
			      
  IPHCTree::NTMET GetScaledType1MET(IPHCTree::NTMET &themet, float scale = 1.) const; 
  
  IPHCTree::NTMET GetSmearedType1MET(IPHCTree::NTMET &themet, vector<NTJet> injets, int jetResFactor = 0) const;

  IPHCTree::NTMET GetGenMET(std::vector<IPHCTree::WDecaysMC> &wAndDecays);

  bool isAnEventSelected(unsigned int nElectrons, unsigned int nMuons,
                         unsigned int nTaus, unsigned int nJets);

  // -------------- BTag methods -------------

  //! btag weight
  void InitSFBWeight(int flagb, int methodb,
                     int systb, int btagAlgo,
                     float btagDiscriCut, int btagNjets);

  void ReInitSFBWeight (int flagb, int methodb,
                        int systb, int btagAlgo, 
			float btagDiscriCut, int btagNjets);

  void ResetParameters4Bweight(int flagb, int methodb, int systb);

  //! Tells whether a particular jet passes the b-tagging requirement
  bool passBtagSelection(const NTJet & jet, const int& algo, const float & discricut) const;

  //! Returns the b-tagging discriminant of a particular jet according to the algorithm specified in the config file.
  double getBtagDiscr(const NTJet & jet, const int& algo) const;

  int GetFlagb () const
  { return flag_btagweight_; }

  int GetMethodb () const
  { return methodb_; }

  int GetSystb () const
  { return systb_; }

  const SFBweight& GetSFBweight() const
  { return sfb_; }
  

  // -------------- PileUp methods -------------

  void GeneratePUWeight(string PUWeightFileName);
  float GetPUWeight();
  vector<double> generate_flat10_weights(TH1D* data_npu_estimated);

  // -------------- Scale factors -------------
      
  void LoadElScaleFactors();
  void LoadMuScaleFactors();
  
  void LoadMuIDScaleFactors();
  std::vector <double > getScaleFactorMuID(double pT, double eta);
  
  void LoadMuIsoScaleFactors12();
  std::vector <double > getScaleFactorMuIso12(double pT, double eta);
  
  void LoadMuIsoScaleFactors20();
  std::vector <double > getScaleFactorMuIso20(double pT, double eta);
  
  
  std::vector <double > getSscaleFactorElectronAllID05(double pT, double eta);
  
  
  
  TH2F * getScaleFactEl() {return scaleFactEl;};
  TH2F * getScaleFactMu() {return scaleFactMu;};

  // -------------- JES -------------

  //Treatment of the JES Uncertainty
  void InitJESUnc ();
  
  //Treatment of the JES Uncertainty
  //void InitJESUnc (char* jecContrib = "Total");
  void InitJESUnc (char* jecContrib );
  void changeJESUncSource (char* jecContrib) {jecSource = jecContrib;}
  
  
  TVector2 GetUnclusMET() const;
  bool looseJetId(const NTJet & theJet) const;
  bool cleanJet(const NTJet & theJet,const std::vector<NTMuon> & muon_cand,
	       const std::vector<NTElectron> & elec_cand) const;
  
  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------

   
   
 public:

  //! requirements
  Requirement cfg;
   
static   double RelIso03PFDeltaBeta(IPHCTree::NTMuon &themuon) {
     
     return (
      ( themuon.isolation["PF03Char"] + max(0., themuon.isolation["PF03Neut"] + themuon.isolation["PF03Phot"]-0.5*themuon.isolation["PF03PU"] ) 
      )/ themuon.p4.Pt() );
      
         
     //return (
     // ( themuon.isolation["PF03Char"] + themuon.isolation["PF03Neut"] + themuon.isolation["PF03Phot"]) 
     // / themuon.p4.Pt() );
      
      
   }
   
  
static   double RelIso04PFDeltaBeta(IPHCTree::NTMuon &themuon) {
     
     return (
      ( themuon.isolation["PF04Char"] + max(0., themuon.isolation["PF04Neut"] + themuon.isolation["PF04Phot"]-0.5*themuon.isolation["PF04PU"] ) 
      )/ themuon.p4.Pt() );
   } 
      
static   double RelIso03PF(IPHCTree::NTElectron &theelec) {
     
     return (
      ( theelec.isolation["PATCharH"] +theelec.isolation["PATNeutH"] + theelec.isolation["PATPhoto"] ) 
      / theelec.p4.Pt() );
   }
   
   
   
      
static   float EffArea03PF(IPHCTree::NTElectron &theelec, double rho) {
     
      
      /*return (
      ( theelec.isolation["PATCharH"] + 
      max(theelec.isolation["PATNeutH"] + theelec.isolation["PATPhoto"]- rho*theelec.isolation["Aeff"], 0. )
      
      )/ theelec.p4.Pt() );*/
     
     	return (
      ( theelec.isolation["PATCharH"] + 
      max(theelec.isolation["PATNeutH"] + theelec.isolation["PATPhoto"]- rho*AeffDR03_2012(theelec.etaSuperCluster ), 0. )
      
      )/ theelec.p4Gsf.Pt() );
		
		
      //return isolation;
      
      
   }
   


static float AeffDR03_2012(double eta){


float theAeff = 0.;                                                                                                       
                                                                                                                            
  double abseta = fabs(eta);                                                                                                
                                                                                                                            
  if(abseta<1.0)                        theAeff = 0.13;                                                                     
  else if(abseta>1.0   && abseta<1.479) theAeff = 0.14;                                                                     
  else if(abseta>1.479 && abseta<2.0)   theAeff = 0.07;                                                                     
  else if(abseta>2.0   && abseta<2.2)   theAeff = 0.09;                                                                     
  else if(abseta>2.2   && abseta<2.3)   theAeff = 0.11;                                                                     
  else if(abseta>2.3   && abseta<2.4)   theAeff = 0.11;                                                                     
  else if(abseta>2.4)                   theAeff = 0.14 ;                                                                                                                                                                                           
                                                                                                                            
   return theAeff;                                                  


     
   
}
  
  
       
static   float EffArea04PF(IPHCTree::NTElectron &theelec, double rho) {
     
      
      return (
      ( theelec.isolation["PATCharH"] + 
      max(theelec.isolation["PATNeutH"] + theelec.isolation["PATPhoto"]- rho*AeffDR04_2012(theelec.etaSuperCluster), 0. )
      
      )/ theelec.p4.Pt() );
   }
   


static float AeffDR04_2012(double eta){

  float theAeff = 0.;
  
  double abseta = fabs(eta);
    
  if(abseta<1.0)			theAeff = 0.208;
  else if(abseta>1.0   && abseta<1.479) theAeff = 0.209;
  else if(abseta>1.479 && abseta<2.0)	theAeff = 0.115;
  else if(abseta>2.0   && abseta<2.2)	theAeff = 0.143;
  else if(abseta>2.2   && abseta<2.3)	theAeff = 0.183;
  else if(abseta>2.3   && abseta<2.4)	theAeff = 0.194;
  else if(abseta>2.4)			theAeff = 0.261;
  return theAeff;
}
   
  
  
  
  
 private:

  //! Pile-Up weights
  std::vector<double> PUWeights;
      
  //! Scale factors
  TH2F * scaleFactEl;
  TH2F * scaleFactMu;
  
  TGraphAsymmErrors * scaleFactMu_ID_abseta_inf0p9;
  TGraphAsymmErrors * scaleFactMu_ID_abseta_0p9_1p2;
  TGraphAsymmErrors * scaleFactMu_ID_abseta_1p2_2p1;
  TGraphAsymmErrors * scaleFactMu_ID_abseta_2p1_2p4;
  
  
  
  TGraphAsymmErrors * scaleFactMu_Iso12_abseta_inf0p9;
  TGraphAsymmErrors * scaleFactMu_Iso12_abseta_0p9_1p2;
  TGraphAsymmErrors * scaleFactMu_Iso12_abseta_1p2_2p1;
  TGraphAsymmErrors * scaleFactMu_Iso12_abseta_2p1_2p4;
  
  
  TGraphAsymmErrors * scaleFactMu_Iso20_abseta_inf0p9;
  TGraphAsymmErrors * scaleFactMu_Iso20_abseta_0p9_1p2;
  TGraphAsymmErrors * scaleFactMu_Iso20_abseta_1p2_2p1;
  TGraphAsymmErrors * scaleFactMu_Iso20_abseta_2p1_2p4;
  
  
  
  //! BTag variables
  int flag_btagweight_;
  SFBweight sfb_;
  int methodb_;
  int systb_;

  //! JES
  TH2F* histo_jesunc_;
  
      struct ltstr
      {
	bool operator()(const char* s1, const char* s2) const
	{
	  return strcmp(s1, s2) < 0;
	}
      };
      
      
  mutable map<const char*, JetCorrectionUncertainty*, ltstr> vsrc;
  char * jecSource;   


   



};

#endif
