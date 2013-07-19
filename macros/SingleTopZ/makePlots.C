{

 gROOT->ProcessLine(".L PlotStack.C+");
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);


 for (int j=0; j<5; j++) {
 // loop over channels
  TString jchan;
  if (j==0) jchan="mumumu";
  else if (j==1) jchan="mumue";
  else if (j==2) jchan="eemu";
  else if (j==3) jchan="eee";
  else if (j==4) jchan="all";

  for (int i=0; i<2; i++) {
  //  loop for SetLogy option 

   for (int k=0; k<2; k++) {
    // loop for ratio plot

    //PlotStack("NVtx_",jchan,"_aftertrigsel", i, k);
    //PlotStack("NVtx_",jchan,"_afterleptsel", i, k);

    // plot met after btag
    //PlotStack("Mt_",jchan,"_afterbjetsel", i, k);
    /*PlotStack("Mt_",jchan,"_afterbjetveto", i, k);
    PlotStack("NJet_",jchan,"_afterZsel", i, k);
    
    PlotStack("NBJet_",jchan,"_afterZsel", i, k);
    PlotStack("NBJet_",jchan,"_afterjetsel", i, k);
    
    PlotStack("InvM_ll_",jchan,"_afterleptsel_highSumPt", i, k);
    PlotStack("InvM_ll_",jchan,"_afterleptsel_lowbin", i, k);
    PlotStack("InvM_ll_",jchan,"_afterjetsel", i, k);
    PlotStack("InvM_ll_",jchan,"_afterleptsel", i, k);*/
    //PlotStack("InvM_ll_",jchan,"_afterbjetsel", i, k);
    
    /*PlotStack("LeptPt_",jchan,"_afterleptsel", i, k);
    PlotStack("LeptPt_",jchan,"_afterjetsel", i, k);
    PlotStack("LeptPt_",jchan,"_afterbjetsel", i, k);*/
    
    //PlotStack("LeptPt_",jchan,"_afterbjetsel", i, k);
    
    PlotStack("NBJet_",jchan,"_afterjetsel", i, k);
    
    
    //PlotStack("LeptZPt_",jchan,"_afterleptsel", i, k);
    //PlotStack("LeptZPt_",jchan,"_afterjetsel", i, k);
    //PlotStack("LeptZPt_",jchan,"_afterbjetsel", i, k);
    
    //PlotStack("LeptWPt_",jchan,"_afterleptsel", i, k);
    //PlotStack("LeptWPt_",jchan,"_afterjetsel", i, k);
    //PlotStack("LeptWPt_",jchan,"_afterbjetsel", i, k);
    //PlotStack("LeptWPt_",jchan,"_afterbjetveto", i, k);
    
    //PlotStack("MET_",jchan,"_afterleptsel_mWT110", i, k);
    //PlotStack("MET_",jchan,"_afterleptsel", i, k);
    //PlotStack("MET_",jchan,"_afterjetsel", i, k);
    //PlotStack("MET_",jchan,"_afterbjetsel", i, k);
    
    //PlotStack("HT_",jchan,"_afterleptsel", i, k);
    //PlotStack("HT_",jchan,"_afterjetsel", i, k);
   // PlotStack("HT_",jchan,"_afterbjetsel", i, k);
    //PlotStack("HT_",jchan,"_afterbjetveto", i, k);
    
    
    //PlotStack("JetPt_",jchan,"_afterleptsel", i, k);
    //PlotStack("JetPt_",jchan,"_afterjetsel", i, k);
    //PlotStack("JetPt_",jchan,"_afterbjetsel", i, k);
    //PlotStack("JetPt_",jchan,"_afterbjetveto", i, k);
    
    
    //PlotStack("RecoTopMass_",jchan,"_afterbjetsel", i, k);
    //PlotStack("RecoTopMass_",jchan,"_afterbjetveto", i, k);
    //PlotStack("RecoPtZ_",jchan,"_afterbjetsel", i, k);
    //PlotStack("RecoPtZ_",jchan,"_afterbjetveto", i, k);
    //PlotStack("Asym_",jchan,"_afterbjetsel", i, k);
    //PlotStack("deltaPhilb_",jchan,"_afterbjetsel", i, k);
    //PlotStack("deltaPhilj_",jchan,"_afterbjetveto", i, k);
    
    //PlotStack("RecoPtZ_",jchan,"_afterbjetveto", i, k);
    
    //PlotStack("RecoPtZ_",jchan,"_afterleptsel", i, k);
     
     
    /* 
    PlotStack("mWT_",jchan,"_afterleptsel",i,k);
    PlotStack("mWT_",jchan,"_afterjetsel",i,k);*/
    //PlotStack("mWT_",jchan,"_afterbjetsel",i,k);
    
    //PlotStack("NJet_",jchan,"_afterbsel",i,k);
    //PlotStack("NLept_",jchan,"_afterbsel",i,k);
    
    
    //PlotStack("deltaRLeptJet_",jchan,"_afterleptsel_mWT110",i,k);
    //PlotStack("deltaRLeptMet_",jchan,"_afterleptsel_mWT110",i,k);
   
    //PlotStack("Nvtx_",jchan,"_afterleptsel",i,k);
    
    //PlotStack("DijetInvM_", jchan,"_afterleptsel_inZpeak",i,k);
    
    //PlotStacl("", jchan, "afterjetsel", i, k);
    
    
   } // end loop k
  } // end loop i
 } // end loop jchan
 }
