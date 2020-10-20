#include <iostream>
#include <memory>
#include "MVAnalysis.h"

using std::string;
using std::cout;
using std::endl;
// template<typename T, typename... Args>
// std::unique_ptr<T> make_unique(Args&&... args) {
//   return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
// }

MVAnalysis::MVAnalysis(const string& mva_algo)//, const string& xmlfile)
{
  //reader_ = std::make_unique<TMVA::Reader>("!Color:!Silent");
  reader_ = std::unique_ptr<TMVA::Reader>(new TMVA::Reader("!Color:!Silent"));
  reader_->AddVariable( "Max_Kpt := TMath::Max(Kmpt,Kppt)", &varList_.Max_Kpt_);
  reader_->AddVariable( "Max_MuMinIPsig := TMath::Max(MumMinIP/MumMinIPE, MupMinIP/MupMinIPE)", &varList_.Max_MuMinIPsig_ );
  reader_->AddVariable( "Max_MinIPsig := TMath::Max(KmtrkMinIP/KmtrkMinIPE, KptrkMinIP/KptrkMinIPE)", &varList_.Max_MinIPsig_ );

  reader_->AddVariable( "Max_DCA := TMath::Max(Kmtrkdcasigbs,Kptrkdcasigbs)",&varList_.Max_DCA_);
  reader_->AddVariable( "Bcosalphabs2d", &varList_.Bcosalphabs2d_);

  reader_->AddVariable( "Blxysig", &varList_.Blxysig_);
  reader_->AddVariable( "Bvtxcl", &varList_.Bvtxcl_);
  reader_->AddVariable( "Bpt", &varList_.Bpt_);
  reader_->AddVariable( "Bsdcasigbs", &varList_.Bsdcasigbs_);
  reader_->AddVariable( "Phimass", &varList_.Phimass_);
  reader_->AddVariable( "BsIso", &varList_.BsIso_);
  reader_->AddVariable( "K_Iso := TMath::Max(kmtrkIso, kptrkIso)", &varList_.K_Iso_);
  
  
  reader_->AddSpectator( "Bmass", &varList_.bmass_);
  reader_->AddSpectator( "Mumumass", &varList_.mumumass_);
  reader_->AddSpectator( "Mumumasserr",&varList_.mumumasserr_);

  // reader_->AddVariable( "Kminpt", &varList_.Kminpt_);
  // reader_->AddVariable( "Bcosalphabs2d", &varList_.Bcosalphabs2d_);
  // reader_->AddVariable( "Kmaxtrkdcasigbs",&varList_.Kmaxtrkdcasigbs_);
  // reader_->AddVariable( "Blxysig", &varList_.Blxysig_);
  // reader_->AddVariable( "Bvtxcl", &varList_.Bvtxcl_);
  //reader_->BookMVA(mva_algo.c_str(), "weights/TMVAClassification_BDT-s01.weights.xml");
  reader_->BookMVA(mva_algo.c_str(), "weights/TMVAClassification_BDTG.weights.xml");
  //reader_->BookMVA(mva_algo.c_str(), "weights/TMVAClassification_BDT.weights.xml");
  //  reader_->BookMVA(mva_algo.c_str(), "weights/TMVAClassification_BDT06.weights.xml");

}

double MVAnalysis::evaluate(const string& mva_algo, const InputVariables& varList) {
  memcpy(&varList_, &varList, sizeof(varList)); // use move syntax here                                                                                                      
  return reader_->EvaluateMVA(mva_algo.c_str());
}
