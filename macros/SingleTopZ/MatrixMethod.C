void ComputeEstimation(double Nt, double Nl, double &NtS, double &NtZ, double Eff, double Fake)
{
  NtS = Eff*( Nt - Fake*Nl)/(Eff-Fake);
  NtZ = Fake*(-Nt + Eff*Nl)/(Eff-Fake);
}

void GenerateStatError(double Nt, double Nl, double &NtS_err, double &NtZ_err, double Eff, double Fake)
{
  double NtS, NtZ;
  TH1F* htS = new TH1F("htS", "htS", 11000, -100, 1000);
  TH1F* htZ = new TH1F("htZ", "htZ", 11000, -100, 1000);
  
  for(int itest=0; itest<1000; itest++)
  {
    double N_rand=gRandom->Poisson(Nl-Nt);
    double Nt_rand=gRandom->Poisson(Nt);
    double Nl_rand=N_rand+Nt_rand;
    
    ComputeEstimation(Nt_rand, Nl_rand, NtS, NtZ, Eff, Fake);
    htS->Fill(NtS);
    htZ->Fill(NtZ);
  }
  
  NtS_err = htS->GetRMS();
  NtZ_err = htZ->GetRMS();
  
  cout<<"S rms : "<<htS->GetRMS()<<endl;
  cout<<"Z rms : "<<htZ->GetRMS()<<endl;

}

void GenerateStatErrorMC(double Nt, double Nl, double Nt_err, double Nl_err, double &NtS_err, double &NtZ_err, double Eff, double Fake)
{
  double NtS, NtZ;
  TH1F* htS = new TH1F("htS", "htS", 11000, -100, 1000);
  TH1F* htZ = new TH1F("htZ", "htZ", 11000, -100, 1000);
  
  for(int itest=0; itest<1000; itest++)
  {
    double N_err = 0;
    if(Nt_err<=Nl_err) N_err = sqrt(Nl_err*Nl_err-Nt_err*Nt_err);
    else N_err = sqrt(Nl_err*Nl_err+Nt_err*Nt_err);
    double N_rand=gRandom->Gaus(Nl-Nt, N_err);
    double Nt_rand=gRandom->Gaus(Nt, Nt_err);
    double Nl_rand=N_rand+Nt_rand;
    
    ComputeEstimation(Nt_rand, Nl_rand, NtS, NtZ, Eff, Fake);
    htS->Fill(NtS);
    htZ->Fill(NtZ);
  }
  
  NtS_err = htS->GetRMS();
  NtZ_err = htZ->GetRMS();
  
  cout<<"S MC rms : "<<htS->GetRMS()<<endl;
  cout<<"Z MC rms : "<<htZ->GetRMS()<<endl;

}

double NtS_Stat_ErrorMC(double Nt, double Nl, double Nt_err, double Nl_err, double Eff, double Fake)
{
  double N_err = 0;
  if(Nt_err<=Nl_err) N_err = sqrt(Nl_err*Nl_err-Nt_err*Nt_err);
  else { cout<<"Warning : Nt_err: "<<Nt_err<<" > Nl_err "<<Nl_err<<endl; N_err = sqrt(Nl_err*Nl_err+Nt_err*Nt_err);}
  double Nt_err = TMath::Power(Eff*Fake/(Eff-Fake)*N_err, 2) + TMath::Power(Eff*(1-Fake)/(Eff-Fake)*Nt_err, 2);
  Nt_err=sqrt(Nt_err);
  return Nt_err;
}

double NtZ_Stat_ErrorMC(double Nt, double Nl, double Nt_err, double Nl_err, double Eff, double Fake)
{
  double N_err = 0;
  if(Nt_err<=Nl_err) N_err = sqrt(Nl_err*Nl_err-Nt_err*Nt_err);
  else N_err = sqrt(Nl_err*Nl_err+Nt_err*Nt_err);
  double Nt_err = TMath::Power(Fake*Eff/(Fake-Eff)*N_err, 2) + TMath::Power(Fake*(1-Eff)/(Fake-Eff)*Nt_err, 2);
  Nt_err=sqrt(Nt_err);
  return Nt_err;
}

double NtS_Stat_Error(double Nt, double Nl, double Eff, double Fake)
{
  return NtS_Stat_ErrorMC(Nt, Nl, sqrt(Nt), sqrt(Nl), Eff, Fake);
}

double NtZ_Stat_Error(double Nt, double Nl, double Eff, double Fake)
{
  return NtZ_Stat_ErrorMC(Nt, Nl, sqrt(Nt), sqrt(Nl), Eff, Fake);
}

double NtS_Syst_Error(double NtS, double Nl, double Eff, double Fake, double Eff_err, double Fake_err)
{
  double Nt_err = TMath::Power(NtS*Fake/Eff*(Eff-Fake)*Eff_err, 2) 
                   + TMath::Power((-Eff*Nl+NtS)/(Eff-Fake)*Fake_err, 2);
  Nt_err=sqrt(Nt_err);
  return Nt_err;
}

double NtZ_Syst_Error(double NtZ, double Nl, double Eff, double Fake, double Eff_err, double Fake_err)
{
  double Nt_err = TMath::Power(NtZ*Eff/Fake*(Fake-Eff)*Fake_err, 2) 
                   + TMath::Power((-Fake*Nl+NtZ)/(Fake-Eff)*Eff_err, 2);
  Nt_err=sqrt(Nt_err);
  return Nt_err;
}

//#################
//  MAIN FUNCTION
//#################

void MatrixMethod(bool isData=false){

  TFile* f = new TFile("MMethod.root", "update");

  string dataset="";
  string channels[] = {"mumumu", "mumue", "eemu", "eee"};


  // Loop over channels
  for(int iChannel=0; iChannel<4; iChannel++)
  {
    cout<<"-------------"<<endl; 
    cout<<"Channel "<<channels[iChannel]<<endl; 
    cout<<"-------------"<<endl; 
    
    TH1F* htight=0;
    TH1F* hloose=0;
    TH1F* htight_err=0;
    TH1F* hloose_err=0;
   
    // Get cut flow histos
    if(isData)
    {
      if(iChannel==0) dataset="DataMu";
      if(iChannel==1) dataset="DataMuEG";
      if(iChannel==2) dataset="DataMuEG";
      if(iChannel==3) dataset="DataEG";    
      htight = (TH1F*) f->Get(string("CutFlow_"+channels[iChannel]+"_tight_"+dataset).c_str());
      hloose = (TH1F*) f->Get(string("CutFlow_"+channels[iChannel]+"_loose_"+dataset).c_str());
      htight_err = (TH1F*) f->Get(string("ErrCutFlow_"+channels[iChannel]+"_tight_"+dataset).c_str());
      hloose_err = (TH1F*) f->Get(string("ErrCutFlow_"+channels[iChannel]+"_loose_"+dataset).c_str());
    }
    else
    {
      dataset="MC";
      TIter nextkey( gDirectory->GetListOfKeys() );
      TKey *key;
      while ( (key = (TKey*)nextkey()))
      {
	TObject *obj = key->ReadObj();
	if(! obj->IsA()->InheritsFrom( TH1F::Class() )) continue;
	
	TString name(obj->GetName());
	if(!name.Contains("CutFlow")) continue;
	if(name.Contains("Data") || name.Contains("FCNC")) continue;
	
	TH1F* h = (TH1F*) obj;
	
	if(!name.Contains(channels[iChannel])) continue;
	//cout<<name<<endl;
		
	if(name.Contains("CutFlow") && !name.Contains("ErrCutFlow"))
	{
	  if(name.Contains("tight"))
	  {
	    cout<<name<<endl;
	    if(!htight) htight = (TH1F*) h->Clone();
	    else htight->Add(h);
	    //htight->Print("all");
	  }
	  if(name.Contains("loose"))
	  {
	    if(!hloose) hloose = (TH1F*) h->Clone();
	    else hloose->Add(h);
	  }
	}
	
	if(name.Contains("ErrCutFlow"))
	{
	  if(name.Contains("tight"))
	  {
	    if(!htight_err) htight_err = (TH1F*) h->Clone();
	    else htight_err->Add(h);
	  }
	  if(name.Contains("loose"))
	  {
	    if(!hloose_err) hloose_err = (TH1F*) h->Clone();
	    else hloose_err->Add(h);
	  }
	}
	
      }
    
    }
  
    if(!htight || !hloose) {cout<<"Histos not found."<<endl; return;}
    if(!htight_err || !hloose_err) {cout<<"Error histos not found."<<endl; return;}
    
    htight->Print("all");


    // Prepare estimations histos
    TH1F* hestimation_S = (TH1F*) htight->Clone();
    hestimation_S->SetName(string("Estimation_"+channels[iChannel]+"_S_"+dataset).c_str());
    hestimation_S->Reset();
    TH1F* hestimation_Z = (TH1F*) htight->Clone();
    hestimation_Z->SetName(string("Estimation_"+channels[iChannel]+"_Z_"+dataset).c_str());
    hestimation_Z->Reset();

    double Nt, Nl, NtS, NtZ, Nt_err, Nl_err, NtS_err_eff, NtZ_err_eff, NtS_err_stat, NtZ_err_stat;
    double Eff=0.95;
    double Eff_err=0.05;
    double Fake=0.05;
    double Fake_err=0.02;
    
    //test
    //ComputeEstimation(105, 300, NtS, NtZ, Eff, Fake);
    //cout<<" Estimation Signal : "<<NtS<<endl;
    //cout<<" Estimation Z+jets : "<<NtZ<<endl;
    
    // Set efficiencies and fake rate
    
    double EffMu, EffMu_err, FakeMu, FakeMu_err; 
    double EffEl, EffEl_err, FakeEl, FakeEl_err;

    //--------------------------------
    // DATA efficiencies and Fake rate
    //--------------------------------

    if(isData)
    {
      // Estimation from Data Z and QCD
      EffMu = 0.976603;
      EffMu_err = 0.000101254;
      FakeMu = 0.0796632; 
      FakeMu_err = 0.0456416; 
      EffEl = 0.970103;
      EffEl_err = 0.000129266;
      FakeEl = 0.378883;
      FakeEl_err = 0.0227932;
    }

    //------------------------------
    // MC efficiencies and Fake rate
    //------------------------------
    if(!isData)
    {
      // Estimation from Z and W with 1 jet.
      // pt > 10 GeV
     /* EffMu = 0.969844;
      EffMu_err = 0.000108481;
      FakeMu = 0.119625;
      FakeMu_err = 0.0153273;
      EffEl = 0.966015;
      EffEl_err = 0.000132852;
      FakeEl = 0.369239;
      FakeEl_err = 0.0164294;*/

      // pt > 20 GeV
      EffMu = 0.976603;
      EffMu_err = 0.000101254;
      FakeMu = 0.0796632; // utilise val pour 2 jet, car pour 1 jet elle est exageree : 0.156161
      FakeMu_err = 0.0456416; // 1jet : 0.0391486
      //EffEl = 0.970103;
      //EffEl_err = 0.000129266;
      //FakeEl = 0.378883;
      //FakeEl_err = 0.0227932;
      // good El Id
      EffEl = 0.971209;
      EffEl_err = 0.000129555;
      FakeEl = 0.352737;
      FakeEl_err = 0.0258883;

      // pt > 20 GeV , loose Iso = 999
      /*EffMu = 0.974195;
      EffMu_err = 0.000106013;
      FakeMu = 0.0181453; // utilise val pour 2 jet, car pour 1 jet elle est exageree : 0.0342193
      FakeMu_err = 0.01117; // 1jet : 0.00917115
      EffEl = 0.970195;
      EffEl_err = 0.000131647;
      FakeEl = 0.329943;
      FakeEl_err = 0.0246272;*/

      // pt > 30 GeV (only for W lepton)
     /* EffMu = 0.984503;
      EffMu_err = 9.60599e-05;
      FakeMu = 0.144276; 
      FakeMu_err = 0.0743225; // low stat. Wjets seems to have mainly low pt fakes
      EffEl = 0.977502;
      EffEl_err = 0.000128996;
      FakeEl = 0.382805;
      FakeEl_err = 0.0340156;*/

     /* EffEl = 0.854;
      EffEl_err = 0.001;
      FakeEl = 0.157;
      FakeEl_err = 0.01;*/ //Test last bin eee
    }
    
    
    if(iChannel==0 || iChannel==2) { Eff=EffMu; Eff_err=EffMu_err; Fake=FakeMu; Fake_err=FakeMu_err; }
    if(iChannel==1 || iChannel==3) { Eff=EffEl; Eff_err=EffEl_err; Fake=FakeEl; Fake_err=FakeEl_err; }

    // Compute estimations
    // Loop on all steps of cut flow
    for(int iCut=0; iCut<8; iCut++) //htight->GetNbinsX()
    {
      cout<<"iCut "<<iCut<<endl; 
      Nt = htight->GetBinContent(iCut);
      Nl = hloose->GetBinContent(iCut);
      if(!htight_err) Nt_err=0;
      else Nt_err = sqrt(htight_err->GetBinContent(iCut));
      if(!hloose_err) Nl_err=0;
      else Nl_err = sqrt(hloose_err->GetBinContent(iCut));

      ComputeEstimation(Nt, Nl, NtS, NtZ, Eff, Fake);
      //GenerateStatError(Nt, Nl, NtS_err, NtZ_err, Eff, Fake);
      //GenerateStatErrorMC(Nt, Nl, Nt_err, Nl_err, NtS_err, NtZ_err, Eff, Fake);
      
      // Calcul des erreurs par propagation
      // Propagate stat error on Nt and Nl
      if(isData) NtS_err_stat = NtS_Stat_Error(Nt, Nl, Eff, Fake);
      else NtS_err_stat = NtS_Stat_ErrorMC(Nt, Nl, Nt_err, Nl_err, Eff, Fake);
      if(isData) NtZ_err_stat = NtZ_Stat_Error(Nt, Nl, Eff, Fake);
      else NtZ_err_stat = NtZ_Stat_ErrorMC(Nt, Nl, Nt_err, Nl_err, Eff, Fake);
      
      // Propagate stat error on Eff and Fake
      NtS_err_eff = NtS_Syst_Error(NtS, Nl, Eff, Fake, Eff_err, Fake_err);
      NtZ_err_eff = NtZ_Syst_Error(NtZ, Nl, Eff, Fake, Eff_err, Fake_err);

      cout<<" Estimation Signal : "<<NtS<<" +/- "<<NtS_err_stat<<" +/- "<<NtS_err_eff<<endl;
      cout<<" Estimation Z+jets : "<<NtZ<<" +/- "<<NtZ_err_stat<<" +/- "<<NtZ_err_eff<<endl;

      hestimation_S->SetBinContent(iCut, NtS); 
      hestimation_Z->SetBinContent(iCut, NtZ); 

      // Stat errors on N already plotted on histo so take only errors from Eff and Fake
      hestimation_S->SetBinError(iCut, NtS_err_eff); 
      hestimation_Z->SetBinError(iCut, NtZ_err_eff); 
    
    }

    hestimation_S->Write();
    hestimation_Z->Write();

  }
  

}
