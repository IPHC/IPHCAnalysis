#include "MatrixMethod_Results_Estimations.C"
//#include "MatrixMethod_Results_EMu.C"
#include "MatrixMethod_Results_Efficiencies.C"
//#include "MatrixMethod_Results_PullTest_EE.C"
//#include "MatrixMethod_Results_PullTest_MuMu.C"
//#include "MatrixMethod_Results_PullTest_EMu.C"
//#include "MatrixMethod_Results_EpsilonsTest_EE.C"
#include "MatrixMethod_Results_EpsilonsTest_MuMu.C"


void MatrixMethod_Results(){
/*
MatrixMethod_Results_PullTest_EE("Tight","Signal");
MatrixMethod_Results_PullTest_EE("Tight","Wlike");
MatrixMethod_Results_PullTest_EE("Tight","QCDlike");

MatrixMethod_Results_PullTest_EE("Medium","Signal");
MatrixMethod_Results_PullTest_EE("Medium","Wlike");
MatrixMethod_Results_PullTest_EE("Mediul","QCDlike");

MatrixMethod_Results_PullTest_EE("Loose","Signal");
MatrixMethod_Results_PullTest_EE("Loose","Wlike");
MatrixMethod_Results_PullTest_EE("Loose","QCDlike");


MatrixMethod_Results_PullTest_MuMu("Tight","Signal");
MatrixMethod_Results_PullTest_MuMu("Tight","Wlike");
MatrixMethod_Results_PullTest_MuMu("Tight","QCDlike");

MatrixMethod_Results_PullTest_MuMu("Medium","Signal");
MatrixMethod_Results_PullTest_MuMu("Medium","Wlike");
MatrixMethod_Results_PullTest_MuMu("Mediul","QCDlike");

MatrixMethod_Results_PullTest_MuMu("Loose","Signal");
MatrixMethod_Results_PullTest_MuMu("Loose","Wlike");
MatrixMethod_Results_PullTest_MuMu("Loose","QCDlike");


MatrixMethod_Results_PullTest_EMu("TT","Signal");
MatrixMethod_Results_PullTest_EMu("TT","Wlike");
MatrixMethod_Results_PullTest_EMu("TT","QCDlike");

MatrixMethod_Results_PullTest_EMu("TL","Signal");
MatrixMethod_Results_PullTest_EMu("TL","Wlike");
MatrixMethod_Results_PullTest_EMu("TL","QCDlike");

MatrixMethod_Results_PullTest_EMu("LT","Signal");
MatrixMethod_Results_PullTest_EMu("LT","Wlike");
MatrixMethod_Results_PullTest_EMu("LT","QCDlike");

MatrixMethod_Results_PullTest_EMu("LL","Signal");
MatrixMethod_Results_PullTest_EMu("LL","Wlike");
MatrixMethod_Results_PullTest_EMu("LL","QCDlike");



MatrixMethod_Results_EpsilonsTest_EE("Tight","Signal");
MatrixMethod_Results_EpsilonsTest_EE("Tight","Wlike");
MatrixMethod_Results_EpsilonsTest_EE("Tight","QCDlike");

MatrixMethod_Results_EpsilonsTest_EE("Medium","Signal");
MatrixMethod_Results_EpsilonsTest_EE("Medium","Wlike");
MatrixMethod_Results_EpsilonsTest_EE("Mediul","QCDlike");

MatrixMethod_Results_EpsilonsTest_EE("Loose","Signal");
MatrixMethod_Results_EpsilonsTest_EE("Loose","Wlike");
MatrixMethod_Results_EpsilonsTest_EE("Loose","QCDlike");
*/

/*MatrixMethod_Results_EpsilonsTest_MuMu("Tight","Signal");
MatrixMethod_Results_EpsilonsTest_MuMu("Tight","Wlike");
MatrixMethod_Results_EpsilonsTest_MuMu("Tight","QCDlike");

MatrixMethod_Results_EpsilonsTest_MuMu("Medium","Signal");
MatrixMethod_Results_EpsilonsTest_MuMu("Medium","Wlike");
MatrixMethod_Results_EpsilonsTest_MuMu("Mediul","QCDlike");

MatrixMethod_Results_EpsilonsTest_MuMu("Loose","Signal");
MatrixMethod_Results_EpsilonsTest_MuMu("Loose","Wlike");
MatrixMethod_Results_EpsilonsTest_MuMu("Loose","QCDlike");
*/


//MatrixMethod_Results_Efficiencies();


//string ntupleType[2] = {"_MC","_DATA"};
string ntupleType[1] = {"_MC"};
//string fast[2] = {"","_Fast"};
string fast[1] = {"_Fast"};
//string other[2] = {"","_other"};
string other[1] = {""};

for(unsigned int i = 0; i < 1; i++){
for(unsigned int j = 0; j < 1; j++){
for(unsigned int k = 0; k < 1; k++){


//MatrixMethod_Results_Estimations("Tight","Signal", ntupleType[i], fast[j], other[k], "EE");
//MatrixMethod_Results_Estimations("Tight","Wlike", ntupleType[i], fast[j], other[k], "EE");
//MatrixMethod_Results_Estimations("Tight","QCDlike", ntupleType[i], fast[j], other[k], "EE");
/*
MatrixMethod_Results_Estimations("Medium","Signal", ntupleType[i], fast[j], other[k], "EE");
MatrixMethod_Results_Estimations("Medium","Wlike", ntupleType[i], fast[j], other[k], "EE");
MatrixMethod_Results_Estimations"(Medium","QCDlike", ntupleType[i], fast[j], other[k], "EE");

MatrixMethod_Results_Estimations("Loose","Signal", ntupleType[i], fast[j], other[k], "EE");
MatrixMethod_Results_Estimations("Loose","Wlike", ntupleType[i], fast[j], other[k], "EE");
MatrixMethod_Results_Estimations("Loose","QCDlike", ntupleType[i], fast[j], other[k], "EE");
*/


MatrixMethod_Results_Estimations("Tight","Signal", ntupleType[i], fast[j], other[k], "MuMu");
MatrixMethod_Results_Estimations("Tight","Wlike", ntupleType[i], fast[j], other[k], "MuMu");
MatrixMethod_Results_Estimations("Tight","QCDlike", ntupleType[i], fast[j], other[k], "MuMu");

MatrixMethod_Results_Estimations("Medium","Signal", ntupleType[i], fast[j], other[k], "MuMu");
MatrixMethod_Results_Estimations("Medium","Wlike", ntupleType[i], fast[j], other[k], "MuMu");
MatrixMethod_Results_Estimations("Medium","QCDlike", ntupleType[i], fast[j], other[k], "MuMu");

MatrixMethod_Results_Estimations("Loose","Signal", ntupleType[i], fast[j], other[k], "MuMu");
MatrixMethod_Results_Estimations("Loose","Wlike", ntupleType[i], fast[j], other[k], "MuMu");
MatrixMethod_Results_Estimations("Loose","QCDlike", ntupleType[i], fast[j], other[k], "MuMu");



//MatrixMethod_Results_EMu("TT","Signal", ntupleType[i], fast[j], other[k]);
//MatrixMethod_Results_EMu("TT","Wlike", ntupleType[i], fast[j], other[k]);
//MatrixMethod_Results_EMu("TT","QCDlike", ntupleType[i], fast[j], other[k]);
/*
MatrixMethod_Results_EMu("TL","Signal", ntupleType[i], fast[j], other[k]);
MatrixMethod_Results_EMu("TL","Wlike", ntupleType[i], fast[j], other[k]);
MatrixMethod_Results_EMu("TL","QCDlike", ntupleType[i], fast[j], other[k]);

MatrixMethod_Results_EMu("LT","Signal", ntupleType[i], fast[j], other[k]);
MatrixMethod_Results_EMu("LT","Wlike", ntupleType[i], fast[j], other[k]);
MatrixMethod_Results_EMu("LT","QCDlike", ntupleType[i], fast[j], other[k]);

MatrixMethod_Results_EMu("LL","Signal", ntupleType[i], fast[j], other[k]);
MatrixMethod_Results_EMu("LL","Wlike", ntupleType[i], fast[j], other[k]);
MatrixMethod_Results_EMu("LL","QCDlike", ntupleType[i], fast[j], other[k]);
*/


}
}
}


}
