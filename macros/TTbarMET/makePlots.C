{

 gROOT->ProcessLine(".L PlotStack.C");
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);


 bool setlogy=0;
 bool ratio=1;
 bool norma=0;
// PlotStack(TString categorie, TString plotname, TString namechan, TString selection, bool setlogy, bool ratio){
  PlotStack("NofEvents","NofEvents","all","All",setlogy,norma,ratio);
  PlotStack("KinPlots","MET","all","NJets cut",setlogy,norma,ratio);

  PlotStack("KinPlots","HT","all","NJets cut",setlogy,norma,ratio);
  PlotStack("KinPlots","M3","all","NJets cut",setlogy,norma,ratio);
  PlotStack("KinPlots","Mtl","all","NJets cut",setlogy,norma,ratio);
  PlotStack("KinPlots","MTw","all","NJets cut",setlogy,norma,ratio);
  PlotStack("KinPlots","PT3","all","NJets cut",setlogy,norma,ratio);
  PlotStack("KinPlots","PTtl","all","NJets cut",setlogy,norma,ratio);
  PlotStack("KinPlots","PTw","all","NJets cut",setlogy,norma,ratio);
  PlotStack("KinPlots","Dphi_l_met","all","NJets cut",setlogy,norma,ratio);
  PlotStack("KinPlots","Dphi_l_4thjet","all","NJets cut",setlogy,norma,ratio);
  PlotStack("KinPlots","Dphi_2tops","all","NJets cut",setlogy,norma,ratio);
  PlotStack("KinPlots","Deta_l_thad","all","NJets cut",setlogy,norma,ratio);
  PlotStack("KinPlots","Deta_l_4thjet","all","NJets cut",setlogy,norma,ratio);

  PlotStack("Electrons","e_Pt","all","NJets cut",setlogy,norma,ratio);
  PlotStack("Electrons","e_Eta","all","NJets cut",setlogy,norma,ratio);

  PlotStack("Muons","mu_Pt","all","NJets cut",setlogy,norma,ratio);
  PlotStack("Muons","mu_Eta","all","NJets cut",setlogy,norma,ratio);

  PlotStack("Jets","j_Pt1","all","NJets cut",setlogy,norma,ratio);
  PlotStack("Jets","j_Eta1","all","NJets cut",setlogy,norma,ratio);
  PlotStack("Jets","j_Pt2","all","NJets cut",setlogy,norma,ratio);
  PlotStack("Jets","j_Eta2","all","NJets cut",setlogy,norma,ratio);
  PlotStack("Jets","j_Pt3","all","NJets cut",setlogy,norma,ratio);
  PlotStack("Jets","j_Eta3","all","NJets cut",setlogy,norma,ratio);
  PlotStack("Jets","j_Pt4","all","NJets cut",setlogy,norma,ratio);
  PlotStack("Jets","j_Eta4","all","NJets cut",setlogy,norma,ratio);




  PlotStack("KinPlots","MET","all","NbtagJets cut1",setlogy,norma,ratio);

  PlotStack("KinPlots","HT","all","NbtagJets cut1",setlogy,norma,ratio);
  PlotStack("KinPlots","M3","all","NbtagJets cut1",setlogy,norma,ratio);
  PlotStack("KinPlots","Mtl","all","NbtagJets cut1",setlogy,norma,ratio);
  PlotStack("KinPlots","MTw","all","NbtagJets cut1",setlogy,norma,ratio);
  PlotStack("KinPlots","PT3","all","NbtagJets cut1",setlogy,norma,ratio);
  PlotStack("KinPlots","PTtl","all","NbtagJets cut1",setlogy,norma,ratio);
  PlotStack("KinPlots","PTw","all","NbtagJets cut1",setlogy,norma,ratio);
  PlotStack("KinPlots","Dphi_l_met","all","NbtagJets cut1",setlogy,norma,ratio);
  PlotStack("KinPlots","Dphi_l_4thjet","all","NbtagJets cut1",setlogy,norma,ratio);
  PlotStack("KinPlots","Dphi_2tops","all","NbtagJets cut1",setlogy,norma,ratio);
  PlotStack("KinPlots","Deta_l_thad","all","NbtagJets cut1",setlogy,norma,ratio);
  PlotStack("KinPlots","Deta_l_4thjet","all","NbtagJets cut1",setlogy,norma,ratio);

/*
 for (int j=0; j<3; j++) {
 // loop over channels
  TString jchan;
  if (j==0) jchan="e";
  else if (j==1) jchan="mu";
  else if (j==2) jchan="all";

  for (int i=0; i<2; i++) {
  //  loop for SetLogy option 

//   for (int k=0; k<2; k++) {
    int k=1;
    // loop for ratio plot


    // after lepton pair cut
    PlotStack("Njets_",jchan,"", i, k);
    PlotStack("NBjets_",jchan,"", i, k);
    PlotStack("Inv",jchan,"MassPair", i, k);
    PlotStack("Met_",jchan,"", i, k);

    
    
//   } // end loop k

  } // end loop i
 } // end loop jchan
*/

}
