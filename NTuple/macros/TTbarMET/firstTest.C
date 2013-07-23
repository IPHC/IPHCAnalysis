#include <iomanip>
#include <iostream>
#include <time.h>
#include "NTFormat/interface/NTEvent.h"

//NTupleAnalysis classes
#include "Selection/interface/SSDiLeptonSelection.h"
#include "Selection/interface/SelectionTable.h"
#include "Tools/interface/Dataset.h"
#include "Tools/interface/AnalysisEnvironmentLoader.h"
#include "Plots/interface/SSDiLepAnaHistoManager.h"


#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

using namespace IPHCTree;
using namespace std;

double pi=acos(-1.);

struct SortByPt
{
  bool operator()( TLorentzVector j1, TLorentzVector j2 ) const
  {
    return j1.Pt() > j2.Pt() ;
  }
};



int main (int argc, char *argv[])
{
 
  cout << "#########################" << endl;
  cout << "Beginning of the program" << endl;
  cout << "#########################" << endl;

  //////////////////////
  //Global variables
  //////////////////////
  vector < Dataset > datasets;
  SSDiLeptonSelection sel;
  float Luminosity = 0;
  float LumiError = 0.;
  // 0: MC - 1: Data - 2 Data & MC
  int DataType = 0;
  int verbosity = -1;

  int AnaStep = 6;//which defines the cuts that the events should pass to be considered as selected

  float Nobs_mumu = 0.;
  //////////////////////

  //////////////////////
  // Initialisation
  //////////////////////
  string xmlFileName;
  cout<<"argc "<<argc<<" "<<argv[0]<<endl;
  if (argc>1 ) xmlFileName = string(argv[1]);
  else xmlFileName = string ("../../config/TTbarMETAnalysis.xml");
  AnalysisEnvironmentLoader anaEL (xmlFileName);
  anaEL.LoadSamples (datasets);	// now the list of datasets written in the xml file is known
  anaEL.LoadSSDiLeptonSelection (sel);	// now the parameters for the selection are given to the selection
  anaEL.LoadGeneralInfo (DataType, Luminosity, LumiError, verbosity);
  int flagOriginal=sel.GetFlagb();
  int methodOriginal=sel.GetMethodb();
  int systOriginal= sel.GetSystb();
  std::cout << " For btag : flag " << flagOriginal << ", method " << methodOriginal << ", syst " << systOriginal << std::endl;
  IPHCTree::NTEvent * event = 0;
  //Selection table
  SelectionTable selTable_allChannels (sel.GetCutList (), datasets, string ("*"));
  SelectionTable selTable_ee (sel.GetCutList (), datasets, string ("ee"));
  SelectionTable selTable_emu (sel.GetCutList (), datasets, string ("emu"));
  SelectionTable selTable_mumu (sel.GetCutList (), datasets, string ("mumu"));
  
/*
  //Book keeping of standard histos
  bool doHistoManager = true;
  SSDiLepAnaHistoManager histoManager;
  if(doHistoManager){
  	histoManager.LoadDatasets (datasets);
  	histoManager.LoadSelectionSteps (sel.GetCutList ());
  	histoManager.LoadChannels (sel.GetChannelList ());
  	histoManager.CreateHistos ();
  }
*/

  TH1F* h_mult_lep_0  =  new TH1F("h_mult_lep_0",     "# of leptons (bf cuts)",   5,0.0,5.0);
  TH1F* h_mult_jet_0  =  new TH1F("h_mult_jet_0",     "# of jets (bf cuts)",     10,0.0,10.0);
  TH1F* h_mult_bjet_0 =  new TH1F("h_mult_bjet_0",    "# of b-jets (bf cuts)",    5,0.0,5.0);
  TH1F* h_pt_lep1_0   =  new TH1F("h_pt_lep1_0",      "PT of leading lepton (bf cuts)",60,0.0,300.0);
  TH1F* h_eta_lep1_0  =  new TH1F("h_eta_lep1_0",     "Eta of leading lepton (bf cuts)",66,-3.3,3.3);
  TH1F* h_pt_jet1_0   =  new TH1F("h_pt_jet1_0",      "PT of leading jet (bf cuts)",60,0.0,300.0);
  TH1F* h_pt_jet2_0   =  new TH1F("h_pt_jet2_0",      "PT of jet 2 (bf cuts)",60,0.0,300.0);
  TH1F* h_pt_jet3_0   =  new TH1F("h_pt_jet3_0",      "PT of jet 3 (bf cuts)",60,0.0,300.0);
  TH1F* h_pt_jet4_0   =  new TH1F("h_pt_jet4_0",      "PT of jet 4 (bf cuts)",60,0.0,300.0);
  TH1F* h_met_0       =  new TH1F("h_met_0",          "MET (bf cuts)",40,0.0,400.0);
  TH1F* h_mt_0        =  new TH1F("h_mt_0",           "MT(invisible+lepton) (bf cuts)",100,0.0,1000.0);
  TH1F* h_mt_0_2      =  new TH1F("h_mt_0_2",         "MT(invisible+lepton) (bf cuts)",20,0.0,1000.0);
  TH1F* h_ht_0        =  new TH1F("h_ht_0",           "HT (bf cuts)",100,0.0,500.0);
  TH1F* h_m3_0        =  new TH1F("h_m3_0",           "M3 (bf cuts)",100,0.0,1000.0);

  TH1F* h_mult_jet_1  =  new TH1F("h_mult_jet_1",     "# of jets (af cuts)",     10,0.0,10.0);
  TH1F* h_met_1       =  new TH1F("h_met_1",          "MET (af cuts)",40,0.0,400.0);
  TH1F* h_met_1_2     =  new TH1F("h_met_1_2",        "MET (af cuts)",10,0.0,400.0);
  TH1F* h_mt_1        =  new TH1F("h_mt_1",           "MT(invisible+lepton) (af cuts)",100,0.0,1000.0);
  TH1F* h_mt_1_2      =  new TH1F("h_mt_1_2",         "MT(invisible+lepton) (af cuts)",20,0.0,1000.0);
  TH1F* h_pt_lep1_1   =  new TH1F("h_pt_lep1_1",      "PT of leading lepton (af cuts)",60,0.0,300.0);
  TH1F* h_eta_lep1_1  =  new TH1F("h_eta_lep1_1",     "Eta of leading lepton (af cuts)",66,-3.3,3.3);
  TH1F* h_pt_jet1_1   =  new TH1F("h_pt_jet1_1",      "PT of leading jet (af cuts)",60,0.0,300.0);
  TH1F* h_pt_jet2_1   =  new TH1F("h_pt_jet2_1",      "PT of jet 2 (af cuts)",60,0.0,300.0);
  TH1F* h_pt_jet3_1   =  new TH1F("h_pt_jet3_1",      "PT of jet 3 (af cuts)",60,0.0,300.0);
  TH1F* h_pt_jet4_1   =  new TH1F("h_pt_jet4_1",      "PT of jet 4 (af cuts)",60,0.0,300.0);
  TH1F* h_ht_1        =  new TH1F("h_ht_1",           "HT (af cuts)",100,0.0,500.0);
  TH1F* h_m3_1        =  new TH1F("h_m3_1",           "M3 (af cuts)",100,0.0,1000.0);
  TH1F* h_mult_bjet_1 =  new TH1F("h_mult_bjet_1",    "# of b-jets (af cuts)",    5,0.0,5.0);
  TH1F* h_pt_bjet1_2  =  new TH1F("h_pt_bjet1_2",     "PT of leading b-jet (af cuts)",60,0.0,300.0);

  TH1F* h_tophad_pt_1      =  new TH1F("h_tophad_pt_1",     "PT of top1 (af cuts)",100,0.,1000.);
  TH1F* h_toplep_pt_1      =  new TH1F("h_toplep_pt_1",     "PT of top2 (af cuts)",100,0.,1000.);
  TH1F* h_wlep_pt_1        =  new TH1F("h_wlep_pt_1",       "PT of w2 (af cuts)",100,0.,1000.);
  TH1F* h_toplep_m_1       =  new TH1F("h_toplep_m_1",      "M of top2 (af cuts)",100,0.,1000.);
  TH1F* h_dphi_l_mis_1     =  new TH1F("h_dphi_l_mis_1",    "Dphi(l,MET) (af cuts)",30,-1.*pi,pi);
  TH1F* h_dphi_l_jet4th_1  =  new TH1F("h_dphi_l_jet4th_1", "Dphi(l,jet4th) (af cuts)",30,-1.*pi,pi);
  TH1F* h_dphi_2top_1      =  new TH1F("h_dphi_2top_1",     "Dphi(top1,top2) (af cuts)",30,-1.*pi,pi);
  TH1F* h_deta_l_tophad_1  =  new TH1F("h_deta_l_tophad_1", "Deta(l,top1) (af cuts)",50,0.,5.);
  TH1F* h_deta_l_jet4th_1  =  new TH1F("h_deta_l_jet4th_1", "Deta(l,jet4th) (af cuts)",50,0.,5.);
  TH1F* h_tophad_pt_2      =  new TH1F("h_tophad_pt_2",     "PT of top1 (af cuts)",100,0.,1000.);
  TH1F* h_toplep_pt_2      =  new TH1F("h_toplep_pt_2",     "PT of top2 (af cuts)",100,0.,1000.);
  TH1F* h_wlep_pt_2        =  new TH1F("h_wlep_pt_2",       "PT of w2 (af cuts)",100,0.,1000.);
  TH1F* h_toplep_m_2       =  new TH1F("h_toplep_m_2",      "M of top2 (af cuts)",100,0.,1000.);
  TH1F* h_dphi_l_mis_2     =  new TH1F("h_dphi_l_mis_2",    "Dphi(l,MET) (af cuts)",30,-1.*pi,pi);
  TH1F* h_dphi_l_jet4th_2  =  new TH1F("h_dphi_l_jet4th_2", "Dphi(l,jet4th) (af cuts)",30,-1.*pi,pi);
  TH1F* h_dphi_2top_2      =  new TH1F("h_dphi_2top_2",     "Dphi(top1,top2) (af cuts)",30,-1.*pi,pi);
  TH1F* h_deta_l_tophad_2  =  new TH1F("h_deta_l_tophad_2", "Deta(l,top1) (af cuts)",50,0.,5.);
  TH1F* h_deta_l_jet4th_2  =  new TH1F("h_deta_l_jet4th_2", "Deta(l,jet4th) (af cuts)",50,0.,5.);
  
  TH2F* h_top1_top2_1       =  new TH2F("h_top1_top2_1",      "M of top2 vs M of top1 (af cuts)",100,0.,1000.,100,0.,1000.);
  TH2F* h_top1_top2_2       =  new TH2F("h_top1_top2_2",      "M of top2 vs M of top1 (af cuts)",100,0.,1000.,100,0.,1000.);


  
  //////////////////////


  cout << "The verbosity mode is " << verbosity << endl;
  cout << "The luminosity is equal to " << Luminosity << endl;
  cout << "The DataType is ";
  switch (DataType) {
  case 0:
    cout << "MC" << endl;
    break;
  case 1:
    cout << "Data" << endl;
    break;
  case 2:
    cout << "Data & MC" << endl;
    break;
  default:
    cout << " unknown" << endl;
    break;
  }
  //////////////////////



  //////////////////////
  //LOOP OVER THE DATASETS
  //////////////////////
  if (verbosity > 0) {
    cout << "#########################" << endl;
    cout << " Loop over the datasets  " << endl;
    cout << "#########################" << endl;
  }

  for (unsigned int d = 0; d < datasets.size (); d++) 
  {

    if(verbosity>2) cout<<"TTbarMET> Dataset: "<<datasets[d].Name()<<endl;
    datasets[d].eventTree ()->SetBranchAddress ("NTEvent", &event);

    unsigned int nEvents = static_cast<unsigned int>(datasets[d].eventTree ()->GetEntries ());
    cout << "TTbarMET> NEvents = " << nEvents << endl;
    cout << "TTbarMET> NEvents to run over = "<<datasets[d].NofEvtsToRunOver()<<endl;


    //LOOP OVER THE EVENTS
    for (unsigned int ievt = 0; ievt < datasets[d].NofEvtsToRunOver(); ievt++)
    {
      float weight = 1.;
      if(datasets[d].isData() == false) weight = datasets[d].NormFactor()*Luminosity; //if Data , weight = 1
      //cout<<"weight "<<weight<<" "<<datasets[d].isData()<<endl;
      datasets[d].eventTree()->GetEntry(ievt);
      IPHCTree::NTTransient::InitializeAfterReading(event);

      if (ievt==0) {
	      cout << "TTbarMET> weight : " << weight << endl;
      }
      if (verbosity > 3)
      {
        std::cout << "event " << ievt 
                  <<" - event number=" << event->general.eventNb
                  <<" - run number="   << event->general.runNb
                  << std::endl;
      }
      if (ievt % 5000 == 0)
	      cout << "TTbarMET> Progress bar : " << ievt << endl;

      //Load event for the selection
      sel.LoadEvent(event);

      //Collection of selected objects
      vector < NTElectron > selElectrons = sel.GetSelectedElectrons ();
      vector < NTMuon > selMuons = sel.GetSelectedMuons ();
      vector < NTJet > selJets = sel.GetSelectedJets (selMuons,selElectrons);
      NTMET themet = sel.GetSelectedMET ();


      // DEFINE NEW OBJECTS
      TLorentzVector resetvector(0.,0.,0.,0.);
      // 1. Leptons
      vector <TLorentzVector> selLepton;
      for (int i=0; i<selElectrons.size(); i++) {
          selLepton.push_back(selElectrons[i].p4);
      }
      for (int i=0; i<selMuons.size(); i++) {
          selLepton.push_back(selMuons[i].p4);
      }
      sort(selLepton.begin(),selLepton.end(),SortByPt());


      // 2. Total hadronic activity
      TLorentzVector all_hadronic;
      for (UInt_t ind=0;ind<selJets.size();ind++) {
        all_hadronic+=selJets[ind].p4;
      }

      // 3. The leptonic W
      TLorentzVector w_leptonic;
      float misset = themet.met();
      float mt_e_mis=-1.;
//      if (ievt<5) std::cout << "MET " << misset << " " << themet.p2.Mod() << " " << themet.p2MuonCorr.Mod() << " " << themet.correction << std::endl;
      TLorentzVector mismis(themet.p2.Px(),themet.p2.Py(),0.,themet.p2.Mod());
      if (misset>0 && selLepton.size()>0) {
         w_leptonic = selLepton[0] + mismis;
         mt_e_mis = sqrt( 2.* selLepton[0].P() * misset *(1. - cos( selLepton[0].Phi() - mismis.Phi()) ));

      }


      //3. The hadronic Top (from M3)
      //4. The 4th jet (after M3)
      TLorentzVector top_hadronic;
      TLorentzVector tmp_b;
      float test_pt=0;
      float test_pt_2=0;

      bool secure_pt_jet_kin= false;
      if (selJets.size()>=3) {
       for (UInt_t ind=0;ind<selJets.size();ind++) {
        for (UInt_t ind1=ind+1; ind1<selJets.size();ind1++) {
         for (UInt_t ind2=ind1+1; ind2<selJets.size();ind2++) {
           TLorentzVector combi_test;
           combi_test=selJets[ind].p4;
           combi_test+=selJets[ind1].p4;
           combi_test+=selJets[ind2].p4;
           if (combi_test.Pt()>test_pt) {
             test_pt=combi_test.Pt();
             top_hadronic=combi_test;
             if (selJets.size()>3) {
              for (UInt_t ind3=0; ind3<selJets.size();ind3++) {
               if (ind3!=ind && ind3!=ind1 && ind3!=ind2 && selJets[ind3].p4.Pt()>test_pt_2) {
                   test_pt_2=selJets[ind3].p4.Pt();
                   tmp_b=selJets[ind3].p4;
                   // mini-test
                   if (selJets[ind].p4.Pt()>25 && selJets[ind1].p4.Pt()>25 && selJets[ind2].p4.Pt()>25 && selJets[ind3].p4.Pt()>25) {
                      secure_pt_jet_kin=true;
                   }
                   else {
                      secure_pt_jet_kin=false;
                   }
          
               }
              }
             }
           }
         }
        }
       }
      }

      //5. The leptonic Top
      TLorentzVector top_leptonic;
      if (misset>0 && selLepton.size()>0 && selJets.size()>3) top_leptonic = w_leptonic + tmp_b;

      //6. The b-jets avant cut en Pt
      vector < NTJet > selBtagJets;
      vector < float >btagDiscri;
      //cas CSV ou JP : recommandation de BTV  
      // CSVL wp = 0.244 ; CSVM wp = 0.679 ; CSVT wp = 0.898 combinedSecondaryVertexBJetTags
      // JPL  wp = 0.275 ; JPM  wp = 0.545 ; JPT  wp = 0.790 jetProbabilityBJetTags
      string btag_algo = "combinedSecondaryVertexBJetTags";
      float btag_DiscriCut=0.244;
      for (UInt_t j = 0; j < selJets.size (); j++) {
//        if (ievt<5) cout << "jet " << j << " algo CSV " << selJets[j].bTag[btag_algo] << " >= ? " << btag_DiscriCut << endl;
        if (selJets[j].bTag[btag_algo] >= btag_DiscriCut) {
            selBtagJets.push_back (selJets[j]);
            btagDiscri.push_back (selJets[j].bTag[btag_algo]);
        }
      }


      // FILL HISTO BEFORE CUTS
      h_mult_lep_0->Fill(selLepton.size(),weight);
      h_mult_jet_0->Fill(selJets.size(),weight);
      h_mult_bjet_0->Fill(selBtagJets.size());
      if (ievt<5) cout << "Bjets " <<  selBtagJets.size() << endl;

      if (selLepton.size()>0) {
        h_pt_lep1_0->Fill(selLepton[0].Pt(),weight);
        h_eta_lep1_0->Fill(selLepton[0].Eta(),weight);
      }

      if (selJets.size()>0) {
        h_pt_jet1_0->Fill(selJets[0].p4.Pt(),weight);
        if (selJets.size()>1) {
          h_pt_jet2_0->Fill(selJets[1].p4.Pt(),weight);
          if (selJets.size()>2) {
            h_pt_jet3_0->Fill(selJets[2].p4.Pt(),weight);
            if (selJets.size()>3) {
              h_pt_jet4_0->Fill(selJets[3].p4.Pt(),weight);
            }
          }
        }
      }  

      h_met_0->Fill(misset,weight);
      h_ht_0->Fill(all_hadronic.Pt(),weight);
      h_mt_0->Fill(mt_e_mis,weight);
      h_mt_0_2->Fill(mt_e_mis,weight);
      h_m3_0->Fill( top_hadronic.M() ,weight); 


      // APPLY SELECTION CUTS

      bool pass_cut = true;
      if (selLepton.size()==0) pass_cut = false ;
      if (selLepton.size()>0 && selLepton[0].Pt()<25.0) pass_cut = false;
      vector < NTJet > selJetsPt25;
      for (UInt_t j = 0; j < selJets.size (); j++) {
        if (selJets[j].p4.Pt()>=25.) {
            selJetsPt25.push_back (selJets[j]);
        }
      }
      if (selJetsPt25.size()<=3) pass_cut = false ;
      if (misset<100.) pass_cut = false; 
      vector < NTJet > selBtagJetsPt25;
      for (UInt_t j = 0; j < selBtagJets.size (); j++) {
        if (selBtagJets[j].p4.Pt()>=25.) {
            selBtagJetsPt25.push_back (selBtagJets[j]);
        }
      }

      if (pass_cut) {

          // 2. Total hadronic activity : recalcul
          all_hadronic = resetvector;
          for (UInt_t ind=0;ind<selJetsPt25.size();ind++) {
           all_hadronic+=selJets[ind].p4;
          }
          // question : est-ce que ces 4 jets selectionnes pour m3 et tmp_b ont bien un Pt>25 GeV???
          if (!secure_pt_jet_kin) {
            // recalcul des info kin
            //3. The hadronic Top (from M3) : recalcul
            //4. The 4th jet (after M3) : recalcul
            top_hadronic = resetvector;
            tmp_b = resetvector;
            test_pt=0;
            test_pt_2=0;
            for (UInt_t ind=0;ind<selJetsPt25.size();ind++) {
             for (UInt_t ind1=ind+1; ind1<selJetsPt25.size();ind1++) {
              for (UInt_t ind2=ind1+1; ind2<selJetsPt25.size();ind2++) {
                TLorentzVector combi_test;
                combi_test=selJetsPt25[ind].p4;
                combi_test+=selJetsPt25[ind1].p4;
                combi_test+=selJetsPt25[ind2].p4;
                if (combi_test.Pt()>test_pt) {
                  test_pt=combi_test.Pt();
                  top_hadronic=combi_test;
                   for (UInt_t ind3=0; ind3<selJetsPt25.size();ind3++) {
                    if (ind3!=ind && ind3!=ind1 && ind3!=ind2 && selJetsPt25[ind3].p4.Pt()>test_pt_2) {
                        test_pt_2=selJetsPt25[ind3].p4.Pt();
                        tmp_b=selJetsPt25[ind3].p4;
                    }
                   }
                  }
                }
               }
              }
              //5. The leptonic Top
              top_leptonic = w_leptonic + tmp_b;
          }

          // FILL HISTO AFTER CUTS

          h_mult_jet_1->Fill(selJetsPt25.size(),weight);
          h_mult_bjet_1->Fill(selBtagJetsPt25.size());

          h_pt_lep1_1->Fill(selLepton[0].Pt(),weight);
          h_eta_lep1_1->Fill(selLepton[0].Eta(),weight);

          h_pt_jet1_1->Fill(selJetsPt25[0].p4.Pt(),weight);
          h_pt_jet2_1->Fill(selJetsPt25[1].p4.Pt(),weight);
          h_pt_jet3_1->Fill(selJetsPt25[2].p4.Pt(),weight);
          h_pt_jet4_1->Fill(selJetsPt25[3].p4.Pt(),weight);

          h_met_1->Fill(misset,weight);
          h_met_1_2->Fill(misset,weight);

          h_mt_1->Fill( mt_e_mis ,weight);
          h_mt_1_2->Fill( mt_e_mis ,weight);

          h_ht_1->Fill(all_hadronic.Pt(),weight);

          h_m3_1->Fill( top_hadronic.M() ,weight);

          float delta_phi_e_miss= selLepton[0].DeltaPhi(mismis);
          float delta_phi_ejet4th= selLepton[0].DeltaPhi(tmp_b);
          float delta_phi_2top= top_hadronic.DeltaPhi(top_leptonic);
          float deta_e_tophad= selLepton[0].Eta() - top_hadronic.Eta();
          if (deta_e_tophad<0.) deta_e_tophad*=-1.;
          float deta_e_jet4th= selLepton[0].Eta() - tmp_b.Eta();
          if (deta_e_jet4th<0.) deta_e_jet4th*=-1.;

          h_tophad_pt_1->Fill( top_hadronic.Pt() ,weight);
          h_toplep_pt_1->Fill( top_leptonic.Pt() ,weight);
          h_wlep_pt_1->Fill( w_leptonic.Pt() ,weight);
          h_toplep_m_1->Fill( top_leptonic.M() ,weight);
          h_dphi_l_mis_1->Fill( delta_phi_e_miss ,weight);
          h_dphi_l_jet4th_1->Fill( delta_phi_ejet4th ,weight);
          h_dphi_2top_1->Fill( delta_phi_2top ,weight);
          h_deta_l_tophad_1->Fill (deta_e_tophad ,weight);
          h_deta_l_jet4th_1->Fill (deta_e_jet4th ,weight);
          h_top1_top2_1->Fill( top_hadronic.M(), top_leptonic.M() ,weight);


          if (selBtagJetsPt25.size()>0) {
           h_pt_bjet1_2->Fill(selBtagJetsPt25[0].p4.Pt());
           h_tophad_pt_2->Fill( top_hadronic.Pt() ,weight);
           h_toplep_pt_2->Fill( top_leptonic.Pt() ,weight);
           h_wlep_pt_2->Fill( w_leptonic.Pt() ,weight);
           h_toplep_m_2->Fill( top_leptonic.M() ,weight);
           h_dphi_l_mis_2->Fill( delta_phi_e_miss ,weight);
           h_dphi_l_jet4th_2->Fill( delta_phi_ejet4th ,weight);
           h_dphi_2top_2->Fill( delta_phi_2top ,weight);
           h_deta_l_tophad_2->Fill (deta_e_tophad ,weight);
           h_deta_l_jet4th_2->Fill (deta_e_jet4th ,weight);
           h_top1_top2_2->Fill( top_hadronic.M(), top_leptonic.M() ,weight);

          }


      }  // end of "after cuts"
      


    }				// end of loop over evts


  }				// end of loop over the datasets 
  cout << "#########################" << endl;
  cout << " Loop over the datasets  " << endl;
  cout << "#########################" << endl;



  ////////////////////////////
  //  Computation after loops
  ////////////////////////////

  ////////////////////////////
  //  Histograms
  ////////////////////////////

  string ofilename= string("CrossSection")+string(".root");
  TFile* fout  = new TFile(ofilename.c_str(),"RECREATE");


    h_mult_lep_0->Write();
    h_mult_jet_0->Write();
    h_mult_bjet_0->Write();
    h_pt_lep1_0->Write();
    h_eta_lep1_0->Write();
    h_pt_jet1_0->Write();
    h_pt_jet2_0->Write();
    h_pt_jet3_0->Write();
    h_pt_jet4_0->Write();
    h_met_0->Write();
    h_mt_0->Write();
    h_mt_0_2->Write();
    h_ht_0->Write();
    h_m3_0->Write();

    
    h_mult_jet_1->Write();
    h_met_1->Write();
    h_met_1_2->Write();    
    h_mt_1->Write();       
    h_mt_1_2->Write();     
    h_pt_lep1_1->Write();
    h_eta_lep1_1->Write(); 
    h_pt_jet1_1->Write();  
    h_pt_jet2_1->Write();
    h_pt_jet3_1->Write(); 
    h_pt_jet4_1->Write(); 
    h_ht_1->Write();       
    h_m3_1->Write();     
    h_mult_bjet_1->Write();
    h_pt_bjet1_2->Write();
    
    h_tophad_pt_1->Write();
    h_toplep_pt_1->Write();
    h_wlep_pt_1->Write();
    h_toplep_m_1->Write();
    h_dphi_l_mis_1->Write();
    h_dphi_l_jet4th_1->Write();
    h_dphi_2top_1->Write();
    h_deta_l_tophad_1->Write();
    h_deta_l_jet4th_1->Write();
    h_tophad_pt_2->Write();
    h_toplep_pt_2->Write();
    h_wlep_pt_2->Write();
    h_toplep_m_2->Write();
    h_dphi_l_mis_2->Write();
    h_dphi_l_jet4th_2->Write();
    h_dphi_2top_2->Write();
    h_deta_l_tophad_2->Write();
    h_deta_l_jet4th_2->Write();

    h_top1_top2_1->Write();
    h_top1_top2_2->Write();

  fout->Close();

  delete fout;


  if (verbosity > 0) {
    cout << "#########################" << endl;
    cout << "    End of the program   " << endl;
    cout << "#########################" << endl;
  }

  return (0);
}


