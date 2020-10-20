#ifndef __MVAnalysis__h
#define __MVAnalysis__h

#include <fstream>
#include <string>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

typedef struct
{
  //float Kminpt_;
  //  float Kppt_;
  //float Kmaxtrkdcasigbs_;
  //  float Kptrkdcasigbs_;
  float KmtrkMinIPSig_;  
  float KptrkMinIPSig_;  
  /* float Max_Kpt_; */
  /* float Max_Mpt_; */
  /* float Max_trk_; */
  /* //float Max_K_MinIP_; */
  float BsIsot_;
  float Bsdcasigbs_;
  float Bcosalphabs2d_;
  float Bvtxcl_;
  float Blxysig_;
  float Phimass_;
  /* float bmass_; */
  /* float mumumass_; */
  /* float mumumasserr_; */
  /* float kmpt_; */
  /* float kppt_; */
  /* float kmtkdca_; */
  /* float kptkdca_; */
  

}  InputVariables;


class MVAnalysis {

 public:

  MVAnalysis(const std::string& mva_algo);//, const std::string& xmlfile);
  virtual ~MVAnalysis() {}
  double evaluate(const std::string& tag, const InputVariables& varList);

  InputVariables varList_;
  std::unique_ptr<TMVA::Reader> reader_;
};
#endif
