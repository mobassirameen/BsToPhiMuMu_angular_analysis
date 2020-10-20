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
  
  float Max_Kpt_;
  float Max_MuMinIPsig_ ;
  float Max_MinIPsig_ ;
  float Max_DCA_ ;
  //float Max_Mpt_;
  //float Max_trk_;
  //float Max_K_MinIP_;
  float Bcosalphabs2d_;
  float Blxysig_;
  float Bvtxcl_;
  float Bpt_ ;
  float Bsdcasigbs_ ; 
  //float Blxysig_;
  float Phimass_ ;
  float BsIso_ ;
  float K_Iso_ ;
  float bmass_;
  float mumumass_;
  float mumumasserr_;
  

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
