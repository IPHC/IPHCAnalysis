{

  gROOT->ProcessLine(".L Compare.C");
 
  char* rootfile="backup_outputProof08-10-12_12-05-06/proof.root"; 
  bool log=false;
  bool print=true;
  bool ratio=true;

  char* channel[4] = {"mumumu", "mumue", "eemu", "eee"};


  for(int i=0; i<4; i++)
  {
  cout<<"Channel "<<channel[i]<<endl;

  Draw2Histos(Form("RecoPtZ_%s_afterbjetsel_FCNCkut", channel[i]), Form("RecoPtZ_%s_afterbjetveto_FCNCkut", channel[i]),
              rootfile, "1 b", "no b", log, Form("RecoPtZ_bjetveto_%s_FCNCkut.gif", channel[i]), print, ratio);


  Draw2Histos(Form("mWT_%s_afterbjetsel_WZ", channel[i]), Form("mWT_%s_afterbjetveto_WZ", channel[i]), 
              rootfile, "1 b", "no b", log, Form("mWT_bjetveto_%s_WZ.gif", channel[i]), print, ratio);

  Draw2Histos(Form("RecoPtZ_%s_afterbjetsel_WZ", channel[i]), Form("RecoPtZ_%s_afterbjetveto_WZ", channel[i]),
              rootfile, "1 b", "no b", log, Form("RecoPtZ_bjetveto_%s_WZ.gif", channel[i]), print, ratio);
 
  Draw2Histos(Form("HT_%s_afterbjetsel_WZ", channel[i]), Form("HT_%s_afterbjetveto_WZ", channel[i]),
              rootfile, "1 b", "no b", log, Form("HT_bjetveto_%s_WZ.gif", channel[i]), print, ratio);

  Draw2Histos(Form("deltaPhilb_%s_afterbjetsel_WZ", channel[i]), Form("deltaPhilj_%s_afterbjetveto_WZ", channel[i]),
              rootfile, "1 b", "no b", log, Form("deltaPhilb_bjetveto_%s_WZ.gif", channel[i]), print, ratio);

  Draw2Histos(Form("Mt_%s_afterbjetsel_WZ", channel[i]), Form("Mt_%s_afterbjetveto_WZ", channel[i]),
              rootfile, "1 b", "no b", log, Form("Mt_bjetveto_%s_WZ.gif", channel[i]), print, ratio);

  Draw2Histos(Form("RecoTopMass_%s_afterbjetsel_WZ", channel[i]), Form("RecoTopMass_%s_afterbjetveto_WZ", channel[i]),
              rootfile, "1 b", "no b", log, Form("RecoTopMass_bjetveto_%s_WZ.gif", channel[i]), print, ratio);

  Draw2Histos(Form("LeptWPt_%s_afterbjetsel_WZ", channel[i]), Form("LeptWPt_%s_afterbjetveto_WZ", channel[i]),
              rootfile, "1 b", "no b", log, Form("LeptWPt_bjetveto_%s_WZ.gif", channel[i]), print, ratio);

  Draw2Histos(Form("JetPt_%s_afterbjetsel_WZ", channel[i]), Form("JetPt_%s_afterbjetveto_WZ", channel[i]),
              rootfile, "1 b", "no b", log, Form("JetPt_bjetveto_%s_WZ.gif", channel[i]), print, ratio);
  }

}
