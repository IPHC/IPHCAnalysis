#ifndef BtagHardcodedConditions_h
#define BtagHardcodedConditions_h



#include <iostream>
#include <vector>
#include <algorithm>



class BtagHardcodedConditions{

 public:
    
  BtagHardcodedConditions();
  ~BtagHardcodedConditions();

  /**
   *   Returns the discriminant for a particular algo/OP
   */
  float getDiscriminant(const std::string & op);
  /**
   *   Returns the name of the b-tag algo used in PAT
   */
  std::string getAlgoName(const std::string & op);

  /**
   *   Returns the algorithm tag for a particular algo/OP, e.g. CSVM -> CSV
   */
  inline std::string getAlgoTag(const std::string & op){
    return op.substr(0,op.length()-1);
  }
  /**
   *   Returns the letter of the algo/OP, e.g. CSVM -> M
   */
  inline char getOPTag(const std::string & op) {
    return op[op.length()-1];
  }

  double GetBtagEfficiency(double pt, double eta, std::string tagger="CSVM");
  double GetBtagScaleFactor(double pt, double eta, std::string tagger="CSVM");
  double GetBtagSFUncertUp(double pt, double eta, std::string tagger="CSVM");
  double GetBtagSFUncertDown(double pt, double eta, std::string tagger="CSVM");

  double GetMistagRate(double pt, double eta, std::string tagger="CSVM");
  double GetMistagScaleFactor(double pt, double eta, std::string tagger="CSVM");
  double GetMistagSFUncertUp(double pt, double eta, std::string tagger="CSVM");
  double GetMistagSFUncertDown(double pt, double eta, std::string tagger="CSVM");
    
private:
  double GetBtagSFUncertainty2011(double pt, double eta, std::string tagger="CSVM");
  double GetMistagSF2011(double pt, double eta, std::string tagger,
	std::string meanminmax);
  double GetMistagCorr2012(double pt, double eta, std::string tagger);
  inline void fillArray(float* a, float* b, int n) {
    for (int i=0;i<n;++i) a[i] = b[i];
  }

  float SFb_TCHPT_error[14], SFb_CSVL_error[14], SFb_CSVM_error[14], SFb_CSVT_error[14], SFb_JPL_error[14], SFb_JPM_error[14], SFb_JPT_error[14];
  float ptmin, ptmax;
  typedef std::vector< float > FVec;
  typedef std::vector< float >::iterator FVecI;
  FVec ptRange;
  inline int findBin(float pt){
    return (std::upper_bound(ptRange.begin(), ptRange.end(), pt)-ptRange.begin())-1;
  }

};


#endif
	 
