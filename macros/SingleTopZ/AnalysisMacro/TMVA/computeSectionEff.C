



double computeSectionEff(int type, double param){

  
  double sectionEff = 0;
  
  if(type == 0){
    sectionEff = 3.077*pow(10, -7) - 0.0001102*param+ param*param*8.562;
  }
  
  
  if(type == 1){
    sectionEff = 2.54*pow(10, -7 ) -4.823*pow(10, -5) *param+ param*param*2.733;
  }
  
  if(type == 2){
    sectionEff = -5.367*pow(10, -8) 9.449*param*pow(10, -6)+ param*param*1.045;
  }
  
  cout << "xs = " << sectionEff << endl;
  
  return sectionEff;
}



void computeSectionEff(double param){

  double xs_init =  0.0856093;
  
  double xs_new  = computeSectionEff(0, param);
  
  double fact = xs_new/xs_init;
   
  cout << "fact " << fact << endl;
  cout << "norm " << 170*fact << endl;


}

