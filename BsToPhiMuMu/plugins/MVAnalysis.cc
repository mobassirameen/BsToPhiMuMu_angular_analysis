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

  // reader_->AddVariable( "Max_Kpt := TMath::Max(Kmpt,Kppt)", &varList_.Max_Kpt_);
  // reader_->AddVariable( "Max_DCA := TMath::Max(Kmtrkdcasigbs,Kptrkdcasigbs)",&varList_.Max_trk_);
  // reader_->AddVariable( "Max_Mpt := TMath::Max(Mumpt,Muppt)", &varList_.Max_Mpt_);
  //reader_->AddVariable( "Max_K_MinIP := TMath::Max(KmtrkMinIP,KptrkMinIP)", &varList_.Max_K_MinIP_);
  reader_->AddVariable( "KmtrkMinIPSig",&varList_.KmtrkMinIPSig_);
  reader_->AddVariable( "KptrkMinIPSig",&varList_.KptrkMinIPSig_);
  reader_->AddVariable( "Phimass", &varList_.Phimass_);
  reader_->AddVariable( "BsIso",&varList_.BsIsot_);
  reader_->AddVariable( "Bsdcasigbs",&varList_.Bsdcasigbs_);
  reader_->AddVariable( "Bcosalphabs2d", &varList_.Bcosalphabs2d_);
  reader_->AddVariable( "Blxysig", &varList_.Blxysig_);
  reader_->AddVariable( "Bvtxcl", &varList_.Bvtxcl_);

  // reader_->AddSpectator( "Bmass", &varList_.bmass_);
  // reader_->AddSpectator( "Mumumass", &varList_.mumumass_);
  // reader_->AddSpectator( "Mumumasserr",&varList_.mumumasserr_);
  // reader_->AddSpectator( "Kmpt", &varList_.kmpt_);      
  // reader_->AddSpectator( "Kppt", &varList_.kppt_); 
  // reader_->AddSpectator( "Kmtrkdcasigbs",&varList_.kmtkdca_);
  // reader_->AddSpectator( "Kptrkdcasigbs",&varList_.kptkdca_);

  // reader_->AddVariable( "Kminpt", &varList_.Kminpt_);
  // reader_->AddVariable( "Bcosalphabs2d", &varList_.Bcosalphabs2d_);
  // reader_->AddVariable( "Kmaxtrkdcasigbs",&varList_.Kmaxtrkdcasigbs_);
  // reader_->AddVariable( "Blxysig", &varList_.Blxysig_);
  // reader_->AddVariable( "Bvtxcl", &varList_.Bvtxcl_);
  //reader_->BookMVA(mva_algo.c_str(), "weights/TMVAClassification_BDT-s01.weights.xml");
  //reader_->BookMVA(mva_algo.c_str(), "weights/TMVAClassification_BDT-s03.weights.xml");
  //reader_->BookMVA(mva_algo.c_str(), "weights/TMVAClassification_BDT_2018_Iso.weights.xml");
  reader_->BookMVA(mva_algo.c_str(), "weights/TMVAClassification_BDT_aftercorrection_2018.weights.xml");
  //reader_->BookMVA(mva_algo.c_str(), "weights/TMVAClassification_BDT_aftercorrection_2017.weights.xml");
  //reader_->BookMVA(mva_algo.c_str(), "weights/TMVAClassification_BDT_aftercorrection.weights.xml");



}

double MVAnalysis::evaluate(const string& mva_algo, const InputVariables& varList) {
  memcpy(&varList_, &varList, sizeof(varList)); // use move syntax here                                                                                                      
  return reader_->EvaluateMVA(mva_algo.c_str());
}
