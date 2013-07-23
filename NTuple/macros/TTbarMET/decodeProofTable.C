{


TFile *_file0 = TFile::Open("backup_outputProof13-07-12_fullmc/proof.root");
TFile *_file1 = TFile::Open("backup_outputProof16-07-12_12_data/proof.root");

const int DataSetMax=7;
TString DataSetName[DataSetMax];
DataSetName[0]="signal"; 
DataSetName[1]="ttbar"; 
DataSetName[2]="w+jets"; 
DataSetName[3]="z+jets"; 
DataSetName[4]="diboson"; 
DataSetName[5]="single t"; 
DataSetName[6]="Data"; 

const int chan=3; // e; mu; merged

const int cut=7;
vector<string> CutName;
CutName.push_back ("Filtre $\\ge$1lep+$\\ge$3jets");
CutName.push_back ("Trigger");
CutName.push_back ("1 Lepton");
CutName.push_back ("Veto on other leptons");
CutName.push_back ("NJets $\\ge$ 4");
CutName.push_back ("NbtagJets $\\ge$ 1");
CutName.push_back ("NbtagJets $\\ge$ 2");
if (CutName.size()!=cut) {cout << " CutName.size()!=cut " << endl; return 0;}

float table[DataSetMax][chan][cut];
float error[DataSetMax][chan][cut];
for (int i1=0; i1<DataSetMax; i1++) {
 for (int i2=0; i2<chan; i2++) {
  for (int i3=0; i3<cut; i3++) {
     table[i1][i2][i3]=0.;
     error[i1][i2][i3]=0.;
  }
 }
}


int index=-1;
const int DataSetNumber=19;

TH1F* hvide = new TH1F("hvide","hvide",cut,0,cut);
for (int i=0; i<DataSetNumber; i++) {
TString NameLecture;
 if (i==0)          {NameLecture="signal";         index=0;   }
 else if (i==1)     {NameLecture="ttbar";          index=1;   }
 else if (i==2)     {NameLecture="W1Jet";          index=2;   }
 else if (i==3)     {NameLecture="W2Jet";          index=2;   }
 else if (i==4)     {NameLecture="W3Jet";          index=2;   }
 else if (i==5)     {NameLecture="W4Jet";          index=2;   }
 else if (i==6)     {NameLecture="DY1";            index=3;   }
 else if (i==7)     {NameLecture="DY2";            index=3;   }
 else if (i==8)     {NameLecture="WW";             index=4;   }
 else if (i==9)     {NameLecture="WZ";             index=4;   }
 else if (i==10)    {NameLecture="ZZ";             index=4;   }
 else if (i==11)    {NameLecture="singleTop1";     index=5;   }
 else if (i==12)    {NameLecture="singleTop2";     index=5;   }
 else if (i==13)    {NameLecture="singleTop3";     index=5;   }
 else if (i==14)    {NameLecture="singleTop4";     index=5;   }
 else if (i==15)    {NameLecture="singleTop5";     index=5;   }
 else if (i==16)    {NameLecture="singleTop6";     index=5;   }
 else if (i==17)    {NameLecture="DataMu";     index=6;   }
 else if (i==18)    {NameLecture="DataEG";     index=6;   }
 else { cout << " DataSetNumber non attribue! " << endl; return 0; }
 TString ha1_name="TTbarMetAnaPlots/NofEvents/e/All/"+NameLecture+"/NofEvents";
 TString ha2_name="TTbarMetAnaPlots/NofEvents/mu/All/"+NameLecture+"/NofEvents";
 TString ha3_name="TTbarMetAnaPlots/NofEvents/merged/All/"+NameLecture+"/NofEvents";

 cout << i << " " << NameLecture << " --> " << DataSetName[index] << endl;
 // acces aux histo
 TH1F* ha1;
 TH1F* ha2;
 TH1F* ha3;
 if (NameLecture=="DataMu") {
  ha1 = (TH1F*)hvide->Clone();
  ha2 = (TH1F*)_file1->Get(ha2_name);
  ha3 = (TH1F*)_file1->Get(ha2_name);
 }
 else if (NameLecture=="DataEG") {
  ha1 = (TH1F*)_file1->Get(ha1_name);
  ha2 = (TH1F*)hvide->Clone();
  ha3 = (TH1F*)_file1->Get(ha1_name);
 } 
 else {
  ha1 = (TH1F*)_file0->Get(ha1_name);
  ha2 = (TH1F*)_file0->Get(ha2_name);
  ha3 = (TH1F*)_file0->Get(ha3_name);
 }

 int Ncut=ha1.GetNbinsX();
 if (Ncut>cut) { cout << "probleme avec le binning des histogrames : Ncut>cut! " << endl; return 0; }

 // remplissage des info dans des tables
 for (int j=0; j<Ncut; j++) {
   table[index][0][j]+= ha1.GetBinContent(j+1); // e
   table[index][1][j]+= ha2.GetBinContent(j+1); // mu
   table[index][2][j]+= ha3.GetBinContent(j+1); // lepton
   error[index][0][j]+= ha1.GetBinError(j+1)*ha1.GetBinError(j+1); // e
   error[index][1][j]+= ha2.GetBinError(j+1)*ha2.GetBinError(j+1); // mu
   error[index][2][j]+= ha3.GetBinError(j+1)*ha3.GetBinError(j+1); // lepton
/*
   cout << index << "," << j 
   << " table " << table[index][0][j] <<  " " << table[index][1][j] << " " << table[index][2][j] 
   << " error " << error[index][0][j] <<  " " << error[index][1][j] << " " << error[index][2][j] << endl;
*/
 }
}

// erreurs au bon format
for (int i=0; i<DataSetMax; i++) {
 for (int j=0; j<cut; j++) {
   error[i][0][j]=sqrt(error[i][0][j]);
   error[i][1][j]=sqrt(error[i][1][j]);
   error[i][2][j]=sqrt(error[i][2][j]);
 }
}


    float tabtotmc[chan][cut];
    float errtotmc[chan][cut];
    float tabtotbg[chan][cut];
    float errtotbg[chan][cut];
    //**************************
    //Compute total mc & bg
    for(int k1=0; k1<chan; ++k1) {
      for(int k2=0; k2<CutName.size(); ++k2) {
        for(int k0=0; k0<DataSetMax-1; ++k0) {  // pour le moment k0=DataSetMax ce sont les donnees
            if (k0>0) {  // pour le moment k0=0 c'est le signal
              tabtotbg[k1][k2] += table[k0][k1][k2];
              errtotbg[k1][k2] += error[k0][k1][k2]*error[k0][k1][k2];
            }
            else {
              tabtotbg[k1][k2] =0.;
              errtotbg[k1][k2] =0.;
            }
            tabtotmc[k1][k2] += table[k0][k1][k2];
            errtotmc[k1][k2] += error[k0][k1][k2]*error[k0][k1][k2];
         }  
         errtotbg[k1][k2] = sqrt(errtotbg[k1][k2]);
         errtotmc[k1][k2] = sqrt(errtotmc[k1][k2]);
      }
    }
    //*********************************


  // print 
  for (int i=0; i<DataSetMax; i++) {
   cout << DataSetName[i] <<   "  e+mu          e           mu " <<  endl;
   for (int j=0; j<cut; j++) {
    cout << " cut "<< j <<  "   " << table[i][2][j] << "      " << table[i][0][j] << "       " << table[i][1][j]<< endl;
   }
  }



  cout << " ECRITURE DE LA TABLE LATEX " << endl;
// ecriture de la table 
    string ofilenametex = "table.tex";
    ofstream ofile(ofilenametex.c_str());
        
    ofile<<"\\documentclass[8pt]{article}"<<endl;
    ofile<<"\\usepackage{lscape}"<<endl;
    ofile<<"\\begin{document}"<<endl;


    ofile.setf(ios::fixed);
    ofile.precision(1);
    


    int cutstart=0;
    for (int IChannel=0; IChannel<chan; IChannel++) {
      
    // Summary tables
//     ofile << "\\clearpage" << endl;
     ofile << "\\begin{landscape}" << endl;
     ofile << "\\begin{table}[p]" << endl;
      
     ofile << "\\begin{tabular}{|l|c|c|c|c|c|c|}" << endl;
     ofile << "\\hline" << endl;
     ofile << "\\hline" << endl; 
     ofile << "Cut & Data & Total MC & Signal  & Total Background & S/B  &  Data/SM MC\\\\" << endl;
     ofile << "\\hline" << endl;
     if (IChannel==chan-1) cutstart=2;
     for(int ic=cutstart; ic<CutName.size(); ++ic) {
          // en attendant d'avoir toutes les infos; [DataSetNumber]--> 0.
          ofile.precision(1);
          ofile << CutName[ic]
                << " & " <<  table[DataSetMax-1][IChannel][ic] << " $\\pm$ "<<  error[DataSetMax-1][IChannel][ic] ;
//          ofile.precision(3);
          ofile.precision(1);
          ofile << " & " <<  tabtotmc[IChannel][ic] << " $\\pm$ "<<  errtotmc[IChannel][ic]
                << " & " <<  table[0][IChannel][ic] << " $\\pm$ "<<  error[0][IChannel][ic]
                << " & " <<  tabtotbg[IChannel][ic] << " $\\pm$ "<<  errtotbg[IChannel][ic] ;
          ofile.precision(3);
          if (tabtotbg[IChannel][ic]>0.) ofile << " & " <<   table[0][IChannel][ic]/(1.*tabtotbg[IChannel][ic]) ;
          else ofile << " &  --- " ;
          if (tabtotbg[IChannel][ic]>0.) ofile << " & " <<   table[DataSetMax-1][IChannel][ic]/(1.*tabtotbg[IChannel][ic]) ;
          else ofile << " &  --- " ;
          ofile  << " \\\\" << endl;
      }
      ofile << "\\hline" << endl;
      ofile << "\\hline" << endl;
      ofile << "\\end{tabular}" << endl;


      if (DataSetNumber>1) {
       ofile << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;
       ofile << "\\hline" << endl;
       ofile << "\\hline" << endl;
       ofile << "Cut " ;
       int idmax=DataSetMax;
       if (DataSetMax>6) idmax=6;
       for(int id=1; id<idmax; ++id) {
          ofile << " & " << DataSetName[id];
       }
       if (DataSetNumber<6) {
        for(int id=idmax; id<6; ++id) {
          ofile << " & " ;
        }
       }
       ofile  << " \\\\" << endl;
       ofile << "\\hline" << endl;
       for(int ic=cutstart; ic<CutName.size(); ++ic) {
          // en attendant d'avoir toutes les infos; [DataSetNumber]--> 0.
          ofile.precision(1);
          ofile << CutName[ic];
          for(int id=1; id<idmax; ++id) { 
               ofile << " & " <<  table[id][IChannel][ic] << " $\\pm$ "<<  error[id][IChannel][ic] ;
          }
          if (DataSetNumber<6) {
            for(int id=idmax; id<6; ++id) {
              ofile << " & " ;
            }
          }
          //ofile.precision(3);
          ofile.precision(1);
          ofile  << " \\\\" << endl;
       }
       ofile << "\\hline " << endl;
       ofile << "\\hline" << endl;
       ofile << "\\end{tabular}" << endl;
      }
      if (IChannel==0) {
      ofile << "\\caption{Number of events - channel e}" << endl;
      ofile << "\\label{tab:SelectionTable_e}" << endl;
      }
      else if (IChannel==1) {
      ofile << "\\caption{Number of events - channel $\\mu$}" << endl;
      ofile << "\\label{tab:SelectionTable_mu}" << endl;
      }
      else if (IChannel==2) {
      ofile << "\\caption{Number of events - channel lepton} " << endl;
      ofile << "\\label{tab:SelectionTable_lep}" << endl;
      }
      ofile << "\\end{table}" << endl;
      ofile << "\\end{landscape}" << endl;
    } // end loop IChannel
    ofile<<"\\end{document}"<<endl;
    string prodpdf = string("pdflatex ")+ofilenametex;
    system(prodpdf.c_str());

}
