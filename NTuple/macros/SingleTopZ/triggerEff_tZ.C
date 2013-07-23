void ComputeMCEff(float &eff, float &err, TH1F* hpass, TH1F* htot)
{
    float pass=hpass->Integral();
    float total=htot->Integral();
     eff = pass/total;
    float tt=htot->GetEntries(); // sans poids!
     err = sqrt(tt*eff*(1.-eff))/tt;
    if(eff==1) err = sqrt(tt*eff*(htot->GetEntries()-1)/htot->GetEntries()*(1.-eff*(htot->GetEntries()-1)/htot->GetEntries()))/tt;    
}

void MC_error(float &eff_MC, float &err_MC, float eff_Zjets, float err_Zjets, float eff_WZ, float err_WZ,
 float a, float err_a, float b, float err_b)
{

  eff_MC = ( eff_Zjets*a + eff_WZ*b )/(a+b);
  /*err_MC = pow(a*err_Zjets,2) + pow(b*(eff_Zjets-eff_WZ)/((a+b)*(a+b))*err_a,2) + 
           pow(b*err_WZ,2) + pow((a*(eff_WZ-eff_Zjets))/((a+b)*(a+b))*err_b,2);
  err_MC = sqrt(err_MC/pow(a+b,2));
  */
  // a+b = nb evts in data , a and b correlated -> use c=a+b and and b for error calculation
  // then replace err_c² = err_a²+err_b²
  err_MC = pow(a*err_Zjets,2) + pow(b*err_WZ,2) + pow((eff_Zjets-eff_WZ)*err_a,2);
  err_MC = sqrt(err_MC/pow(a+b,2));
  
}

void MC_error(float &eff_MC, float &err_MC, float eff_Zjets, float err_Zjets, float eff_WZ, float err_WZ,
 float a, float err_a, float a_sf, float err_a_sf, float c, float err_c)
{
  // c = a*a_sf+b*b_sf
  // a : DY
  // c : all data
  eff_MC = ( eff_Zjets*a*a_sf + eff_WZ*(c-a*a_sf) )/c;
  err_MC = pow(a*a_sf*err_Zjets,2) + pow((c-a*a_sf)*err_WZ,2) + pow(a_sf*(eff_Zjets-eff_WZ)*err_a,2)
           + pow(a*(eff_Zjets-eff_WZ)*err_a_sf,2) + pow(a*a_sf*(-eff_Zjets+eff_WZ)/c*err_c,2);
  err_MC = sqrt(err_MC/pow(c,2));

}

// pour eee + eemu
void MC_error(float &eff_MC, float &err_MC, float eff_Zjets, float err_Zjets, float eff_WZ, float err_WZ,
 float a1, float err_a1, float a1_sf, float err_a1_sf, float c1, float err_c1,
 float a2, float err_a2, float a2_sf, float err_a2_sf, float c2, float err_c2)
{
  // c = a*a_sf+b*b_sf
  eff_MC = ( eff_Zjets*a1*a1_sf + eff_WZ*(c1-a1*a1_sf) + eff_Zjets*a2*a2_sf + eff_WZ*(c2-a2*a2_sf) )/(c1 + c2);
  err_MC = pow((a1*a1_sf+a2*a2_sf)*err_Zjets,2) + pow((c1-a1*a1_sf+c2-a2*a2_sf)*err_WZ,2)
           + pow(a1_sf*(eff_Zjets-eff_WZ)*err_a1,2) + pow(a2_sf*(eff_Zjets-eff_WZ)*err_a2,2)
           + pow(a1*(eff_Zjets-eff_WZ)*err_a1_sf,2) + pow(a2*(eff_Zjets-eff_WZ)*err_a2_sf,2)
	   + pow((a1*a1_sf+a2*a2_sf)*(-eff_Zjets+eff_WZ)/(c1+c2)*err_c1,2)
	   + pow((a1*a1_sf+a2*a2_sf)*(-eff_Zjets+eff_WZ)/(c1+c2)*err_c2,2);
  err_MC = sqrt(err_MC/pow(c1+c2,2));

}


void Correl(double &correl, double &correl_err, TFile* f, string dataset, string chan, bool btag=false)
{
  f->cd();
  TString histo_correl_FCNCkut_name = "correlationMET_"+chan+"_"+dataset;
  TString histo_correl_FCNCkut_btag_name = "correlationBTAGMET_"+chan+"_"+dataset;
  TH2F* correl_FCNCkut;
  if(!btag) correl_FCNCkut = (TH2F*)gROOT->FindObject(histo_correl_FCNCkut_name);
  else correl_FCNCkut = (TH2F*)gROOT->FindObject(histo_correl_FCNCkut_btag_name);

  correl = 0;
    if(correl_FCNCkut->GetBinContent(2,2)!=0 && correl_FCNCkut->GetBinContent(1,1)!=0)
      correl = correl_FCNCkut->GetBinContent(2,2) / (correl_FCNCkut->GetBinContent(2,1) + correl_FCNCkut->GetBinContent(2,2))
               * (correl_FCNCkut->GetBinContent(1,1) + correl_FCNCkut->GetBinContent(1,2)) / correl_FCNCkut->GetBinContent(1,2);
      
  correl_err = 0;
    if(correl_FCNCkut->GetBinContent(2,2)!=0 && correl_FCNCkut->GetBinContent(1,1)!=0 &&
         correl_FCNCkut->GetBinContent(2,1)!=0 && correl_FCNCkut->GetBinContent(1,2)!=0)
  {
    correl_err = TMath::Power(correl_FCNCkut->GetBinError(1,1)/
                               (correl_FCNCkut->GetBinContent(1,1)+correl_FCNCkut->GetBinContent(1,2)), 2)
                + TMath::Power(correl_FCNCkut->GetBinError(1,2)*correl_FCNCkut->GetBinContent(1,1)/
		               correl_FCNCkut->GetBinContent(1,2)/(correl_FCNCkut->GetBinContent(1,1)+correl_FCNCkut->GetBinContent(1,2)), 2)
                + TMath::Power(correl_FCNCkut->GetBinError(2,1)/
                               (correl_FCNCkut->GetBinContent(2,1)+correl_FCNCkut->GetBinContent(2,2)), 2)
                + TMath::Power(correl_FCNCkut->GetBinError(2,2)*correl_FCNCkut->GetBinContent(2,1)/
		               correl_FCNCkut->GetBinContent(2,1)/(correl_FCNCkut->GetBinContent(2,2)+correl_FCNCkut->GetBinContent(2,1)), 2);
    correl_err = sqrt(correl_err);
    correl_err *= correl;
  }
      
}

triggerEff(string chan){

  TFile *f1 = new TFile("TriggerPlots_tZ/TriggerPlots_Data.root"); 
  if (chan=="eee") f1 = new TFile("TriggerPlots_tZ/TriggerPlots_Data_3lept_excl.root"); // eee + eemu triggered with ee
  f1->cd();
  
  TString name;
  if (chan=="mumumu") {
     name="Muon";
  }
  else if (chan=="eee") {
     name="Elec";
  }
  else if (chan=="mumue") {
     name="DiMuonElec";
  }
  else if (chan=="eemu") {
     name="DiElecMuon";
  }
  string name_2lept=name;
  if (chan=="mumue") name_2lept="DiElecMuon"; // trig emu for eemu and mumue, only eemu filled
  
  TString histo_Data_name_pt = "HTriggerEff_"+name+"Sel_pT_MET1";
  TH1F* selElectron_pt_data            = (TH1F*)gROOT->FindObject(histo_Data_name_pt);
  TString histo_Data_name_pttrig = "HTriggerEff_"+name+"Sel_selTrig_pT_MET1";
  TH1F* selElectron_ElTrigger_pt_data  = (TH1F*)gROOT->FindObject(histo_Data_name_pttrig);
  								    
  TString histo_Data_name_pt2 = "HTriggerEff_"+name+"Sel_pT_MET2";
  TH1F* selElectron_pt_data2            = (TH1F*)gROOT->FindObject(histo_Data_name_pt2);
  TString histo_Data_name_pttrig2 = "HTriggerEff_"+name+"Sel_selTrig_pT_MET2";
  TH1F* selElectron_ElTrigger_pt_data2  = (TH1F*)gROOT->FindObject(histo_Data_name_pttrig2);

  TString histo_Data_name_pt3 = "HTriggerEff_"+name+"Sel_pT_MET3";
  TH1F* selElectron_pt_data3            = (TH1F*)gROOT->FindObject(histo_Data_name_pt3);
  TString histo_Data_name_pttrig3 = "HTriggerEff_"+name+"Sel_selTrig_pT_MET3";
  TH1F* selElectron_ElTrigger_pt_data3  = (TH1F*)gROOT->FindObject(histo_Data_name_pttrig3);


  // 2>= leptons -> mostly Zjets 
  TFile *f1_Zjets = new TFile("TriggerPlots_tZ/TriggerPlots_Data_2lept_incl.root"); 
  f1_Zjets->cd();

  TString histo_Data_name_2lept_pt = "HTriggerEff_"+name_2lept+"Sel_pT_MET1";
  TH1F* selElectron_pt_data_Zjets            = (TH1F*)gROOT->FindObject(histo_Data_name_2lept_pt);
  TString histo_Data_name_2lept_pttrig = "HTriggerEff_"+name_2lept+"Sel_selTrig_pT_MET1";
  TH1F* selElectron_ElTrigger_pt_data_Zjets  = (TH1F*)gROOT->FindObject(histo_Data_name_2lept_pttrig);
  								    
  TString histo_Data_name_2lept_pt2 = "HTriggerEff_"+name_2lept+"Sel_pT_MET2";
  TH1F* selElectron_pt_data2_Zjets            = (TH1F*)gROOT->FindObject(histo_Data_name_2lept_pt2);
  TString histo_Data_name_2lept_pttrig2 = "HTriggerEff_"+name_2lept+"Sel_selTrig_pT_MET2";
  TH1F* selElectron_ElTrigger_pt_data2_Zjets  = (TH1F*)gROOT->FindObject(histo_Data_name_2lept_pttrig2);

  TString histo_Data_name_2lept_pt3 = "HTriggerEff_"+name_2lept+"Sel_pT_MET3";
  TH1F* selElectron_pt_data3_Zjets            = (TH1F*)gROOT->FindObject(histo_Data_name_2lept_pt3);
  TString histo_Data_name_2lept_pttrig3 = "HTriggerEff_"+name_2lept+"Sel_selTrig_pT_MET3";
  TH1F* selElectron_ElTrigger_pt_data3_Zjets  = (TH1F*)gROOT->FindObject(histo_Data_name_2lept_pttrig3);


  // because not enough statistic for Zjets MC -> use dilepton selection
  TFile *f2 = new TFile("TriggerPlots_tZ/TriggerPlots_MC_2lept_incl.root"); 
  f2->cd();

  TString histo_Zjets_name_pt = "HTriggerEff_"+name_2lept+"Sel_pT_Zjets";
  TH1F* selElectron_pt_Zjets            = (TH1F*)gROOT->FindObject(histo_Zjets_name_pt);
  TString histo_Zjets_name_pttrig = "HTriggerEff_"+name_2lept+"Sel_selTrig_pT_Zjets";
  TH1F* selElectron_ElTrigger_pt_Zjets  = (TH1F*)gROOT->FindObject(histo_Zjets_name_pttrig);
  								    
  TFile *f3 = new TFile("TriggerPlots_tZ/TriggerPlots_MC.root"); 
  f3->cd();

  TString histo_WZ_name_pt = "HTriggerEff_"+name+"Sel_pT_WZ";
  TH1F* selElectron_pt_WZ            = (TH1F*)gROOT->FindObject(histo_WZ_name_pt);
  TString histo_WZ_name_pttrig = "HTriggerEff_"+name+"Sel_selTrig_pT_WZ";
  TH1F* selElectron_ElTrigger_pt_WZ  = (TH1F*)gROOT->FindObject(histo_WZ_name_pttrig);

  TString histo_FCNCkut_name_pt = "HTriggerEff_"+name+"Sel_pT_FCNCkut";
  TH1F* selElectron_pt_FCNCkut            = (TH1F*)gROOT->FindObject(histo_FCNCkut_name_pt);
  TString histo_FCNCkut_name_pttrig = "HTriggerEff_"+name+"Sel_selTrig_pT_FCNCkut";
  TH1F* selElectron_ElTrigger_pt_FCNCkut  = (TH1F*)gROOT->FindObject(histo_FCNCkut_name_pttrig);
   
  								    
    // compute the inclusive trigger efficiencies and print them

    // merge the 3 periods because statistical error too large to observe differences
    float pass_data=0;
    if(selElectron_ElTrigger_pt_data) pass_data=selElectron_ElTrigger_pt_data->GetEntries()+selElectron_ElTrigger_pt_data2->GetEntries()+selElectron_ElTrigger_pt_data3->GetEntries();
    float total_data=0;
    if(selElectron_pt_data) total_data=selElectron_pt_data->GetEntries()+selElectron_pt_data2->GetEntries()+selElectron_pt_data3->GetEntries();

    float pass_data_Zjets=0;
    if(selElectron_ElTrigger_pt_data_Zjets) pass_data_Zjets=selElectron_ElTrigger_pt_data_Zjets->GetEntries()+selElectron_ElTrigger_pt_data2_Zjets->GetEntries()+selElectron_ElTrigger_pt_data3_Zjets->GetEntries();
    float total_data_Zjets=0;
    if(selElectron_pt_data_Zjets) total_data_Zjets=selElectron_pt_data_Zjets->GetEntries()+selElectron_pt_data2_Zjets->GetEntries()+selElectron_pt_data3_Zjets->GetEntries();


    float eff_data = pass_data/total_data;
    float err_data = sqrt(total_data*eff_data*(1.-eff_data))/total_data;
    if(eff_data==1) err_data = sqrt((pass_data-1)*(1.-(pass_data-1)/total_data))/total_data;
    cout << " pass_data " << pass_data << " total_data " << total_data << " eff " << eff_data << 
            " +/- " << err_data << endl; 
    float err_syst_data = 0;
    if(chan=="eee") err_syst_data = 0.02;
	    
    float eff_data_Zjets = pass_data_Zjets/total_data_Zjets;
    float err_data_Zjets = sqrt(total_data_Zjets*eff_data_Zjets*(1.-eff_data_Zjets))/total_data_Zjets;
    if(eff_data_Zjets==1) err_data_Zjets = sqrt((pass_data_Zjets-1)*(1.-(pass_data_Zjets-1)/total_data_Zjets))/total_data_Zjets;
    cout << " pass_data_Zjets " << pass_data_Zjets << " total_data_Zjets " << total_data_Zjets << " eff " << eff_data_Zjets << 
            " +/- " << err_data_Zjets << endl << endl; 


    float eff_Zjets, err_Zjets;
    ComputeMCEff(eff_Zjets, err_Zjets, selElectron_ElTrigger_pt_Zjets, selElectron_pt_Zjets);
    cout << " Eff " << "Zjets" << " " <<  eff_Zjets << " +/- " << err_Zjets << "     ("<<selElectron_ElTrigger_pt_Zjets->GetEntries()<<" evt MC pass)"<< endl;
        
    float eff_WZ, err_WZ;
    ComputeMCEff(eff_WZ, err_WZ, selElectron_ElTrigger_pt_WZ, selElectron_pt_WZ);
    cout << " Eff " << "WZ" << " " <<  eff_WZ << " +/- " << err_WZ << endl;
        
    float eff_FCNCkut, err_FCNCkut;
    ComputeMCEff(eff_FCNCkut, err_FCNCkut, selElectron_ElTrigger_pt_FCNCkut, selElectron_pt_FCNCkut);
    cout << " Eff " << "FCNCkut" << " " <<  eff_FCNCkut << " +/- " << err_FCNCkut << endl;

    // WZ and DY combination taken into account ratio after lepton selection
    // DY : 95.2+/-0.8 
    // WZ : 802.4+/-2.2 SF 0.60+/-0.20
    // SF mumumu 3.76+/-1.94 SF 0.66
    // SF mumue  1.52+/-0.25 SF 0.68
    // SF eemu   3.66+/-2.63 SF 0.60
    // SF eee    2.24+/-0.31 SF 0.71

    // Correct SF +/- stat. uncert.
    // mumumu 3.68+/-0.35
    // mumue  1.53+/-0.06
    // eemu   3.52+/-0.41
    // eee    2.15+/-0.06

    float eff_MC, err_MC, err_syst_MC;
    
    /*if (chan=="mumumu")
      MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 9.2*3.76, 9.2*1.94, 245.7*0.66, 245.7*0.2);
      //MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 9.2, 3.6, 3.76, 1.94, 196, 14); // quasiment pareil
    if (chan=="mumue")
      MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 145.9*1.52, 145.9*0.25, 216.1*0.68, 216.1*0.2);
    if (chan=="eemu")
      MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 3.9*3.66, 3.9*2.63, 180.3*0.6, 180.3*0.2);
    if (chan=="eee")
      MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 91.3*2.24, 91.3*0.31, 160.4*0.71, 160.4*0.2);
    */
    
    // Fully separate stat and syst errors
    // stat (data+MC)
    if (chan=="mumumu")
      //MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 9.2, 3.6, 3.76, 1.91, 196, 14);
      MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 9.2, 3.6, 3.69, 1.87, 196, 14);
    if (chan=="mumue")
      //MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 145.9, 9.9, 1.52, 0.24, 356, 18.9);
      MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 145.9, 9.9, 1.53, 0.24, 356, 18.9);
    if (chan=="eemu")
      //MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 3.9, 1.5, 3.66, 2.60, 133, 11.5);
      MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 3.9, 1.5, 3.59, 2.63, 133, 11.5);
    if (chan=="eee")
    {
      //MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 91.3, 7.8, 2.24, 0.30, 250, 15.8);
      // Can compute eff_MC for eee + eemu, but we are not sure that the proportion of WZ and DY is the same when eemu is 
      // triggered by ee or emu trigger. Anyway get a very small difference.
      //MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 3.9, 1.5, 3.66, 2.60, 133, 11.5, 91.3, 7.8, 2.24, 0.30, 250, 15.8);
      MC_error(eff_MC, err_MC, eff_Zjets, err_Zjets, eff_WZ, err_WZ, 91.3, 7.8, 2.15, 0.29, 250, 15.8);
    }
    // syst
    if (chan=="mumumu")
      //MC_error(eff_MC, err_syst_MC, eff_Zjets, (1-eff_Zjets), eff_WZ, 0, 9.2, 0, 3.76, 0.35, 196, 0);
      MC_error(eff_MC, err_syst_MC, eff_Zjets, (1-eff_Zjets), eff_WZ, 0, 9.2, 0, 3.69, 1.23, 196, 0);
    if (chan=="mumue")
      //MC_error(eff_MC, err_syst_MC, eff_Zjets, (1-eff_Zjets), eff_WZ, 0, 145.9, 0, 1.52, 0.06, 356, 0);
      MC_error(eff_MC, err_syst_MC, eff_Zjets, (1-eff_Zjets), eff_WZ, 0, 145.9, 0, 1.53, 0.35, 356, 0);
    if (chan=="eemu")
      //MC_error(eff_MC, err_syst_MC, eff_Zjets, (1-eff_Zjets), eff_WZ, 0, 3.9, 0, 3.66, 0.42, 133, 0);
      MC_error(eff_MC, err_syst_MC, eff_Zjets, (1-eff_Zjets), eff_WZ, 0, 3.9, 0, 3.59, 1.04, 133, 0);
    if (chan=="eee")
    {
      //MC_error(eff_MC, err_syst_MC, eff_Zjets, (1-eff_Zjets), eff_WZ, 0, 91.3, 0, 2.24, 0.06, 250, 0);
      //MC_error(eff_MC, err_syst_MC, eff_Zjets, (1-eff_Zjets), eff_WZ, 0, 3.9, 0, 3.66, 0.42, 133, 0, 91.3, 0, 2.24, 0.06, 250, 0);
      MC_error(eff_MC, err_syst_MC, eff_Zjets, (1-eff_Zjets), eff_WZ, 0, 91.3, 0, 2.15, 0.27, 250, 0);
    }
    cout << " Eff " << "MC" << " " <<  eff_MC << " +/- " << err_MC << " +/- " << err_syst_MC << " +/- " << 0.01 << endl;
    // With 1% for correlations
    cout << " Eff " << "MC" << " " <<  eff_MC << " +/- " << err_MC << " +/- " << sqrt(err_syst_MC*err_syst_MC + 0.01*0.01)<< endl;
    cout << " Eff " << "MC" << " " <<  eff_MC << " +/- " << sqrt(err_MC*err_MC + err_syst_MC*err_syst_MC + 0.01*0.01)<< endl;
    
    
    // Last implementation
    float SF_data_MC= eff_data/eff_MC;
    float err_SF = pow( pow(err_data/eff_data, 2.) + pow(err_MC/eff_MC,2) , 0.5 )*SF_data_MC;
    float err_syst_SF = pow( pow(err_syst_data/eff_data, 2.) + pow(err_syst_MC/eff_MC,2) , 0.5 )*SF_data_MC;
    float err_correl_SF = pow( pow(0.01/eff_data, 2.) + pow(0.01/eff_MC,2) , 0.5 )*SF_data_MC;
    
    float corr_signal = eff_FCNCkut/eff_MC;
    float err_corr_signal = pow( pow(err_FCNCkut/eff_FCNCkut, 2.) + pow(err_MC/eff_MC,2) , 0.5 )*corr_signal;
    
      
    cout << endl;
    cout << " Eff Data     " << eff_data << " +/- " << err_data << " +/- " << err_syst_data << " +/- " << 0.01 << endl;
    // With 1% for correlations
    cout << " Eff Data     " << " " << eff_data << " +/- " << err_data << " +/- " << sqrt(err_syst_data*err_syst_data + 0.01*0.01)<< endl;
    cout << " Eff Data     " << " " << eff_data << " +/- " << sqrt(err_data*err_data + err_syst_data*err_syst_data + 0.01*0.01)<< endl;
    cout << " Eff Data2lep " << eff_data_Zjets << " +/- " << err_data_Zjets << endl;
    cout << endl;
    cout << " SF           " << SF_data_MC << " +/- " << err_SF << " +/- " << err_syst_SF << " +/- " << err_correl_SF<<endl;
    // With 1% for correlations
    cout << " SF           " << SF_data_MC << " +/- " << err_SF << " +/- " << sqrt(err_syst_SF*err_syst_SF + err_correl_SF*err_correl_SF) <<endl;
    cout << " SF           " << SF_data_MC << " +/- " << sqrt(err_SF*err_SF + err_syst_SF*err_syst_SF + err_correl_SF*err_correl_SF) <<endl;
    cout << " ratio Eff_s/Eff_b  " << corr_signal << " +/- " << err_corr_signal << endl;
    cout << endl;

    float SF_data_MC_Zjets= eff_data_Zjets/eff_Zjets;
    float err_SF_Zjets = pow( pow(err_data_Zjets/eff_data_Zjets, 2.) + pow(err_Zjets/eff_Zjets,2) , 0.5 )*SF_data_MC_Zjets;
    cout << " SF Zjets     " << SF_data_MC_Zjets << " +/- " << err_SF_Zjets <<endl;
    cout << endl;

    double correl, correl_err;
    Correl(correl, correl_err, f3, "FCNCkut", chan, false);
    cout << " correl FCNC  " << correl << " +/- "<< correl_err << endl;
    Correl(correl, correl_err, f3, "FCNCkut", chan, true);
    cout << " correl FCNC btag " << correl << " +/- "<< correl_err << endl;
    Correl(correl, correl_err, f3, "WZ", chan, false);
    cout << " correl WZ  " << correl << " +/- "<< correl_err << endl;
    Correl(correl, correl_err, f3, "WZ", chan, true);
    cout << " correl WZ btag " << correl << " +/- "<< correl_err << endl;
    Correl(correl, correl_err, f2, "Zjets", chan, false);
    cout << " correl Z  " << correl << " +/- "<< correl_err << endl;
    Correl(correl, correl_err, f2, "Zjets", chan, true);
    cout << " correl Z btag " << correl << " +/- "<< correl_err << endl;

    cout << endl;
}
