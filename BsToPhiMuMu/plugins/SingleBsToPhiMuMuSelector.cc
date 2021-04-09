
//----------------------------------------------
// SELECTOR CODE FOR BSTOPHIMUMU ANALYSIS 
// @author: N.Sahoo, NISER, BHUBANESWAR
//@author: D.K.Sahoo, IIT, BHUBANESWAR
// to do: add phi momentum info to the analyzer
// 2017-07-11: added phi angle info 
// 2018-03-02: (a) commented out the extra print statements, (b) added trigger efficiency print statements
// 2018-03-09: added more counters to dump more infos at each step  (for ex. how many candidates per event after selection ?)
//----------------------------------------------
// author: Chandiprasad Kar
//----------------------------------------------

#define SingleBsToPhiMuMuSelector_cxx

#include <iostream>
#include <sstream>
#include <map>
#include "SingleBsToPhiMuMuSelector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TProof.h>
#include <string.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "TMVA/Factory.h"
#include "TMVA/TMVAGui.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


//-------------------
// Global Constants
//-------------------
const double PHI_MASS  = 1.01946; // GeV
const double PHI_WIDTH = 0.0508; // GeV 
const double MUON_MASS = 0.10565837;
const double KAON_MASS = 0.493677;
const double PION_MASS = 0.13957061;
const double PROTON_MASS = 0.938272081;


//-------------------------
// user defined variables
//-------------------------
TDatime t_begin_ , t_now_ ;
int n_processed_, n_selected_; 
int n_triggers0, n_triggers1;
int n_passMuonID_bdtbkg_ ,n_total_bdtbkg_, n_passvtxcl_bdtbkg_ ;

int n_total_, n_passMuonID_, n_passSelCut_, n_passBestB_; 
int n_passBestB_bkg, n_passPresel_bdtbkg_;
TTree *tree_;
TTree *gentree_;

//bool bdt_bkg = true;
bool bdt_bkg = false;//


//--------------------------------
// Branch variables for new tree
//--------------------------------
int    Nb             = 0;
int    Npv            = 0;
double Mumumass       = 0;
double Mumumasserr    = 0;
double Phimass        = 0;
double Kmpt           = 0;
double Kppt           = 0;
double Kmeta          = 0;
double Kpeta           = 0;
double Kmphi           = 0;
double Kpphi           = 0;
double Kmtrkdcasigbs  = 0;
double Kptrkdcasigbs  = 0;
double Bsdcasigbs     = 0;
double Mumpt           = 0;
double Muppt           = 0;
double Mumeta           = 0;
double Mupeta           = 0;
double Mumphi           = 0;
double Mupphi           = 0;
double Mumdcasigbs  = 0;
double Mupdcasigbs  = 0;
double dimuvtxcl    = 0;
double dimulsig     = 0;
double dimucosalphabs = 0;
double dimulensig = 0;
double dimuDCA = 0;
double cosdimuon = 0;
double cosphimup = 0;
double cosphimum = 0;
double pmum_mass = 0;
double pmup_mass = 0;


double Bmass           = 0;
double BOmass          = 0;
double BOSmass         = 0;
double Kst0mass        = 0;
double Kstmass         = 0;

double Lmass           = 0;
double LBmass          = 0;
double Lsmass           = 0;
double LBsmass          = 0;

double Lhtmass         = 0;
double LBhtmass        = 0;
double Lpihtmass         = 0;
double LBpihtmass        = 0;

double Lhtptmass         = 0;
double LBhtptmass        = 0;

double Lltmass         = 0;
double LBltmass        = 0;
double Bmass_kmu1      = 0;
double Bmass_kmu2      = 0;
double Bmass_kk1_s     = 0;
double Bmass_kk2_s     = 0;
double Bmass_mumu1_s     = 0;
double Bmass_mumu2_s     = 0;

double Bmass_kkmumu    = 0;
double Bmass_kk_s      = 0;
double Bmass_mumu_s    = 0;

 
double Bpt            = 0;
double Beta           = 0;
double Bphi           = 0;
double Bvtxcl         = 0;
double Blxysig        = 0;
double Bcosalphabs    = 0;
double Bcosalphabs2d  = 0;
double Bctau          = 0;
double Phipt           = 0;
double Phiphi          = 0;
double Phieta          = 0;
double Bdt            = -99;
double Puw8           = -2;
double mumuangle      = -99.;
double kkangle        = -99.;
double mu      = -99.;

double MumMinIP = -99.;
double MupMinIP = -99.;
double MumMinIPE = -99.;
double MupMinIPE = -99.;
double KptrkMinIP = -99.;
double KptrkMinIPE = -99.;
double KmtrkMinIP = -99.;
double KmtrkMinIPE = -99.;

double MumMinIP2D = -99.;
double MupMinIP2D = -99.;
double MumMinIP2DE = -99.;
double MupMinIP2DE = -99.;
double KptrkMinIP2D = -99.;
double KptrkMinIP2DE = -99.;
double KmtrkMinIP2D = -99.;
double KmtrkMinIP2DE = -99.;
double mupIso_ = -99.;
double mumIso_ = -99.;
double kptrkIso_ = -99.;
double kmtrkIso_ = -99.;
double BsIso_   = -99.;

double Q2             = 0;
double dimupt         = 0;
double dimueta        = 0;
double CosThetaL      = 999;
double CosThetaK      = 999;
double Phi            = 999;

double dr0 = -99.0;
double dr1 = -99.0;
double dpt0 = -99.0;
double dpt1 = -99.0;

int  mtrkqual       = 99;
int  ptrkqual       = 99;
int    Triggers       = 0;
int    JpsiTriggers       = 0;
int    PsiPTriggers       = 0;
int    LMNTTriggers       = 0;

int    goodMuon       = 0;
int    goodPresel     = 0;

// Branches for Generator level information
double  genBpid      = 999;

double  genBPt       = 0;
double  genBEta      = 0;
double  genBPhi      = 0;
double  genBVtxX     = 0;
double  genBVtxY     = 0;
double  genBVtxZ     = 0;
double  genMupPt     = 0;
double  genMupEta    = 0;
double  genMupPhi    = 0;
double  genMumPt     = 0;
double  genMumEta    = 0;
double  genMumPhi    = 0;

double  genKpPt     = 0;
double  genKpEta    = 0;
double  genKpPhi    = 0;
double  genKmPt     = 0;
double  genKmEta    = 0;
double  genKmPhi    = 0;

double  gendimuPt    = 0;
double  gendimuEta   = 0;
double  gendimuPhi   = 0;

double  genPhipt    = 0;
double  genPhiphi   = 0;
double  genPhieta   = 0;

double  genQ2        = 0;
double  genCosThetaL = 999;
double  genCosThetaK = 999;
double  genPhi       = 999;


void ClearEvent()
{//{{{
  Nb             = 0;
  Npv            = 0;
  Mumumass       = 0;
  Mumumasserr    = 0;
  Phimass        = 0;
  Kmpt           = 0;
  Kppt           = 0;
  Kmeta           = 0;
  Kpeta           = 0;
  Kmphi           = 0;
  Kpphi           = 0;

  Kmtrkdcasigbs    = 0;
  Kptrkdcasigbs    = 0;
  Mumpt            = 0;
  Muppt            = 0;
  Mumeta           = 0;
  Mupeta           = 0;
  Mumphi           = 0;
  Mupphi           = 0;
  dimuvtxcl        = 0;
  dimulsig         = 0;
  dimucosalphabs   = 0;
  dimulensig = 0;
  dimuDCA =0;
  cosdimuon = 0;
  cosphimup = 0;
  cosphimum = 0;
  pmum_mass = 0;
  pmup_mass = 0;

  Mumdcasigbs      = 0;
  Mupdcasigbs      = 0;

  MumMinIP = -99.0;
  MupMinIP = -99.0;
  MumMinIPE = -99.0;
  MupMinIPE = -99.0;
  KptrkMinIP = -99.0;
  KptrkMinIPE = -99.0;
  KmtrkMinIP = -99.0;
  KmtrkMinIPE = -99.0;
  MumMinIP2D = -99.0;
  MupMinIP2D = -99.0;
  MumMinIP2DE = -99.0;
  MupMinIP2DE = -99.0;
  KptrkMinIP2D = -99.0;
  KptrkMinIP2DE = -99.0;
  KmtrkMinIP2D = -99.0;
  KmtrkMinIP2DE = -99.0;
  Bsdcasigbs     = 0;
  
  mupIso_ = -99.;
  mumIso_ = -99.;
  kptrkIso_ = -99.;
  kmtrkIso_ = -99.;
  BsIso_    = -99.;

  Bmass           = 0;
  BOmass          = 0;
  BOSmass         = 0;
  Kst0mass        = 0;
  Kstmass         = 0;
  LBmass          = 0;
  Lmass           = 0;
  LBsmass          = 0;
  Lsmass           = 0;
  Lhtmass         = 0;
  LBhtmass        = 0;
  Lpihtmass         = 0;
  LBpihtmass        = 0;
  Lhtptmass         = 0;
  LBhtptmass        = 0;
  Lltmass         = 0;
  LBltmass        = 0;
  Bmass_kmu1       = 0;
  Bmass_kmu2       = 0;
  Bmass_kkmumu     = 0;
  Bmass_kk_s       = 0;
  Bmass_mumu_s     = 0;
  Bmass_kk1_s      = 0;
  Bmass_mumu1_s      = 0;
  Bmass_kk2_s      = 0;
  Bmass_mumu2_s      = 0;

  Bpt            = 0;
  Beta           = 0;
  Bphi           = 0;
  Bvtxcl         = 0;
  Blxysig        = 0;
  Bcosalphabs    = 0;
  Bcosalphabs2d  = 0;
  Bctau          = 0;
  Puw8           = 0;
  Phipt          = 0;
  Phiphi         = 0;
  Phieta         = 0;
  Bdt            = -99.;
  Q2             = 0;
  dimupt         = 0;
  dimueta        = 0;
  CosThetaL      = 999;
  CosThetaK      = 999;
  Phi            = 999;

  dr0            = -99.0;
  dr1            = -99.0;
  dpt0            = -99.0;
  dpt1            = -99.0;
  mtrkqual       = 99;
  ptrkqual       = 99;
  Triggers       = 0;
  JpsiTriggers       = 0;
  PsiPTriggers       = 0;
  LMNTTriggers       = 0;

  goodMuon       = 0;
  goodPresel     = 0;

  //mc
  genBpid        = 999;
  genBPt         = 0;
  genBEta        = 0;
  genBPhi        = 0;
  genBVtxX       = 0;
  genBVtxY       = 0;
  genBVtxZ       = 0;

  genMupPt       = 0;
  genMupEta      = 0;
  genMupPhi      = 0;
  genMumPt       = 0;
  genMumEta      = 0;
  genMumPhi      = 0;

  genPhipt       = 0;
  genPhiphi      = 0;
  genPhieta      = 0;


  genKpPt     = 0;
  genKpEta    = 0;
  genKpPhi    = 0;
  genKmPt     = 0;
  genKmEta    = 0;
  genKmPhi    = 0;    

  gendimuPt      = 0;
  gendimuEta     = 0;
  gendimuPhi     = 0;

  genQ2          = 0;
  genCosThetaL   = 999;
  genCosThetaK   = 999;
  genPhi         = 999;

}//}}}

void str_replace(std::string& str, const std::string& oldStr, const std::string& newStr)
{//{{{
  size_t pos = 0;
  while((pos = str.find(oldStr, pos)) != std::string::npos)
    {
      str.replace(pos, oldStr.length(), newStr);
      pos += newStr.length();
    }
}//}}}

string get_option_value(string option, string name)
{//{{{
  vector<string> args;
  istringstream f(option);
  string s;    
  while (getline(f, s, ';')) {
    args.push_back(s);
  }
  
  string value; 
  for(vector<string>::iterator it = args.begin(); it != args.end(); ++it) {
    value = *it; 
    unsigned found = value.find(name);
    if (found == 0) {
      str_replace(value, name+"=", ""); 
      break; 
    }
  }
  return value; 
}//}}}

void SingleBsToPhiMuMuSelector::Begin(TTree * /*tree*/)
{//{{{
  t_begin_.Set(); 
  printf("\n ---------- Begin Job ---------- \n");
  t_begin_.Print();

  n_processed_ = 0;
  n_selected_  = 0;

  n_triggers0 = 0;
  n_triggers1 = 0;


  n_total_        = 0; 
  n_passMuonID_   = 0; 
  n_passSelCut_   = 0; 
  n_passBestB_    = 0;
  n_passBestB_bkg = 0;
  n_total_bdtbkg_ = 0;
  n_passMuonID_bdtbkg_ = 0;
  n_passvtxcl_bdtbkg_  = 0;
  n_passPresel_bdtbkg_ = 0;
}//}}}



void SingleBsToPhiMuMuSelector::SlaveBegin(TTree * /*tree*/)
{

  string option = GetOption();
  string cut = get_option_value(option, "cut");
  if(cut== "cut_bdt"){
    mvAna_   = std::unique_ptr<MVAnalysis>(new MVAnalysis(mvAlgo_));
  }
  tree_ = new TTree("tree", "tree"); 
  
  tree_->Branch("Npv"        , &Npv           , "Npv/I");
  tree_->Branch("Nb"        , &Nb           , "Nb/I");
  tree_->Branch("Mumumass"      , &Mumumass      , "Mumumass/D");
  tree_->Branch("Mumumasserr"   , &Mumumasserr   , "Mumumasserr/D");
  tree_->Branch("Phimass"       , &Phimass       , "Phimass/D");
  tree_->Branch("Kmpt"          , &Kmpt          , "Kmpt/D");
  tree_->Branch("Kppt"          , &Kppt          , "Kppt/D");
  tree_->Branch("Kmeta"          , &Kmeta          , "Kmeta/D");
  tree_->Branch("Kpeta"          , &Kpeta          , "Kpeta/D");
  tree_->Branch("Kmphi"          , &Kmphi          , "Kmphi/D");
  tree_->Branch("Kpphi"          , &Kpphi          , "Kpphi/D");

  tree_->Branch("Kmtrkdcasigbs" , &Kmtrkdcasigbs , "Kmtrkdcasigbs/D");
  tree_->Branch("Kptrkdcasigbs" , &Kptrkdcasigbs , "Kptrkdcasigbs/D");

  tree_->Branch("Mumpt"          , &Mumpt          , "Mumpt/D");
  tree_->Branch("Muppt"          , &Muppt          , "Muppt/D");
  tree_->Branch("Mumeta"          , &Mumeta          , "Mumeta/D");
  tree_->Branch("Mupeta"          , &Mupeta          , "Mupeta/D");
  tree_->Branch("Mumphi"          , &Mumphi          , "Mumphi/D");
  tree_->Branch("Mupphi"          , &Mupphi          , "Mupphi/D");
  tree_->Branch("Mumdcasigbs" , &Mumdcasigbs , "Mumdcasigbs/D");
  tree_->Branch("Mupdcasigbs" , &Mupdcasigbs , "Mupdcasigbs/D");
  tree_->Branch("Bsdcasigbs", &Bsdcasigbs, "Bsdcasigbs/D");
  tree_->Branch("MumMinIP",&MumMinIP, "MumMinIP/D" );
  tree_->Branch("MupMinIP" ,&MupMinIP, "MupMinIP/D" );
  tree_->Branch("MumMinIPE",&MumMinIPE, "MumMinIPE/D" );
  tree_->Branch("MupMinIPE" ,&MupMinIPE, "MupMinIPE/D" );
  tree_->Branch("KptrkMinIP", &KptrkMinIP, "KptrkMinIP/D" );
  tree_->Branch("KptrkMinIPE", &KptrkMinIPE, "KptrkMinIPE/D");
  tree_->Branch("KmtrkMinIP", &KmtrkMinIP, "KmtrkMinIP/D");
  tree_->Branch("KmtrkMinIPE" ,&KmtrkMinIPE, "KmtrkMinIPE/D");

  tree_->Branch("MumMinIP2D",&MumMinIP2D, "MumMinIP2D/D" );
  tree_->Branch("MupMinIP2D" ,&MupMinIP2D, "MupMinIP2D/D" );
  tree_->Branch("MumMinIP2DE",&MumMinIP2DE, "MumMinIP2DE/D" );
  tree_->Branch("MupMinIP2DE" ,&MupMinIP2DE, "MupMinIP2DE/D" );
  tree_->Branch("KptrkMinIP2D", &KptrkMinIP2D, "KptrkMinIP2D/D" );
  tree_->Branch("KptrkMinIP2DE", &KptrkMinIP2DE, "KptrkMinIP2DE/D");
  tree_->Branch("KmtrkMinIP2D", &KmtrkMinIP2D, "KmtrkMinIP2D/D");
  tree_->Branch("KmtrkMinIP2DE" ,&KmtrkMinIP2DE, "KmtrkMinIP2DE/D");

  tree_->Branch("mupIso", &mupIso_, "mupIso/D");
  tree_->Branch("mumIso", &mumIso_, "mumIso/D");
  tree_->Branch("BsIso" , &BsIso_, "BsIso/D");
  tree_->Branch("kptrkIso", &kptrkIso_, "kptrkIso/D");
  tree_->Branch("kmtrkIso", &kmtrkIso_, "kmtrkIso/D");
  
  tree_->Branch("Bmass", &Bmass, "Bmass/D");
  tree_->Branch("BOmass", &BOmass, "BOmass/D");
  tree_->Branch("BOSmass",&BOSmass, "BOSmass/D");
  tree_->Branch("Kst0mass", &Kst0mass, "Kst0mass/D");
  tree_->Branch("Kstmass", &Kstmass, "Kstmass/D");
  tree_->Branch("LBmass",&LBmass, "LBmass/D");
  tree_->Branch("Lmass", &Lmass, "Lmass/D");
  tree_->Branch("Lsmass", &Lsmass, "Lsmass/D");
  tree_->Branch("LBsmass",&LBsmass, "LBsmass/D");
  tree_->Branch("Lhtmass", &Lhtmass, "Lhtmass/D");
  tree_->Branch("LBhtmass", &LBhtmass, "LBhtmass/D");
  tree_->Branch("Lltmass", &Lltmass, "Lltmass/D");
  tree_->Branch("LBltmass", &LBltmass, "LBltmass/D");
  tree_->Branch("LBhtptmass", &LBhtptmass, "LBhtptmass/D");
  tree_->Branch("Lhtptmass", &Lhtptmass, "Lhtptmass/D");
  tree_->Branch("Lpihtmass", &Lpihtmass, "Lpihtmass/D");
  tree_->Branch("LBpihtmass", &LBpihtmass, "LBpihtmass/D");

  tree_->Branch("Bmass_kmu1", &Bmass_kmu1, "Bmass_kmu1/D");
  tree_->Branch("Bmass_kmu2", &Bmass_kmu2, "Bmass_kmu2/D");
  tree_->Branch("Bmass_kkmumu", &Bmass_kkmumu, "Bmass_kkmumu/D");

  tree_->Branch("Bmass_kk1_s", &Bmass_kk1_s, "Bmass_kk1_s/D");
  tree_->Branch("Bmass_mumu1_s", &Bmass_mumu1_s, "Bmass_mumu1_s/D");
  tree_->Branch("Bmass_kk2_s", &Bmass_kk2_s, "Bmass_kk2_s/D");
  tree_->Branch("Bmass_mumu2_s", &Bmass_mumu2_s, "Bmass_mumu2_s/D");

  tree_->Branch("Bmass_kk_s", &Bmass_kk_s, "Bmass_kk_s/D");
  tree_->Branch("Bmass_mumu_s", &Bmass_mumu_s, "Bmass_mumu_s/D");

  tree_->Branch("Bpt"           , &Bpt           , "Bpt/D");
  tree_->Branch("Beta"          , &Beta          , "Beta/D");
  tree_->Branch("Bphi"          , &Bphi          , "Bphi/D");
  tree_->Branch("Phipt"          , &Phipt          , "Phipt/D");
  tree_->Branch("Phieta"          , &Phieta          , "Phieta/D");
  tree_->Branch("Phiphi"          , &Phiphi          , "Phiphi/D");
  tree_->Branch("Bdt"           , &Bdt           , "Bdt/D");

  tree_->Branch("Bvtxcl"        , &Bvtxcl        , "Bvtxcl/D");
  tree_->Branch("Blxysig"       , &Blxysig       , "Blxysig/D");
  tree_->Branch("Bcosalphabs"   , &Bcosalphabs   , "Bcosalphabs/D");
  tree_->Branch("Bcosalphabs2d" , &Bcosalphabs2d , "Bcosalphabs2d/D");
  tree_->Branch("Bctau"         , &Bctau         , "Bctau/D");
  
  tree_->Branch("Q2"            , &Q2            , "Q2/D");
  tree_->Branch("dimupt"        , &dimupt        , "dimupt/D");
  tree_->Branch("dimueta"       , &dimueta       , "dimueta/D");
  tree_->Branch("dimuvtxcl"     , &dimuvtxcl     , "dimuvtxcl/D" );
  tree_->Branch("dimulsig"      , &dimulsig      , "dimulsig/D");
  tree_->Branch("dimucosalphabs", &dimucosalphabs, "dimucosalphabs/D");
  tree_->Branch("dimulensig", &dimulensig, "dimulensig/D");
  tree_->Branch("dimuDCA", &dimuDCA, "dimuDCA/D");
  tree_->Branch("cosdimuon", &cosdimuon, "cosdimuon/D");
  tree_->Branch("cosphimup", &cosphimup, "cosphimup/D");
  tree_->Branch("cosphimum", &cosphimum, "cosphimum/D");
  tree_->Branch("pmum_mass", &pmum_mass, "pmum_mass/D");
  tree_->Branch("pmup_mass", &pmup_mass, "pmup_mass/D");
  
  tree_->Branch("CosThetaL"     , &CosThetaL     , "CosThetaL/D");
  tree_->Branch("CosThetaK"     , &CosThetaK     , "CosThetaK/D");
  tree_->Branch("Phi"           , &Phi           , "Phi/D");
  tree_->Branch("ptrkqual"      , &ptrkqual      , "ptrkqual/I");
  tree_->Branch("mtrkqual"      , &mtrkqual      , "mtrkqual/I");
  tree_->Branch("dr0"           , &dr0           , "dr0/D");
  tree_->Branch("dr1"           , &dr1           , "dr1/D");
  tree_->Branch("dpt0"           , &dpt0           , "dpt0/D");
  tree_->Branch("dpt1"           , &dpt1           , "dpt1/D");
  tree_->Branch("Triggers"      , &Triggers      , "Triggers/I");
  tree_->Branch("JpsiTriggers"      , &JpsiTriggers      , "JpsiTriggers/I");
  tree_->Branch("PsiPTriggers"      , &PsiPTriggers      , "PsiPTriggers/I");
  tree_->Branch("LMNTTriggers"      , &LMNTTriggers      , "LMNTTriggers/I");

  tree_->Branch("goodMuon"      , &goodMuon      , "goodMuon/I");
  tree_->Branch("goddPresel"    , &goodPresel    , "goodPresel/I");

  gentree_ = new TTree("gentree", "gentree"); 
  string datatype = get_option_value(option, "datatype");
  std::map<string,int> maptype;
  maptype.insert(std::pair<string,int>("data",1));
  maptype.insert(std::pair<string,int>("mc.lite",2));
  maptype.insert(std::pair<string,int>("mc",999));
  maptype.insert(std::pair<string,int>("mc.nogen",998));
  maptype.insert(std::pair<string,int>("mc.hlt",997));
  switch (maptype[datatype]) {
  case 1:
    break;
  case 2:
    tree_->Branch("Puw8"       ,&Puw8, "Puw8/D");  
    tree_->Branch("genBpid"      , &genBpid      , "genBpid/D");
    // tree_->Branch("genBPhi"      , &genBPhi      , "genBPhi/D");
    // tree_->Branch("genMupPt"     , &genMupPt     , "genMupPt/D");
    // tree_->Branch("genMupEta"    , &genMupEta    , "genMupEta/D");
    // tree_->Branch("genMupPhi"    , &genMupPhi    , "genMupPhi/D");
    // tree_->Branch("genMumPt"     , &genMumPt     , "genMumPt/D");
    // tree_->Branch("genMumEta"    , &genMumEta    , "genMumEta/D");
    // tree_->Branch("genMumPhi"    , &genMumPhi    , "genMumPhi/D");
    
    // tree_->Branch("gendimuPt"    , &gendimuPt    , "gendimuPt/D");
    // tree_->Branch("gendimuEta"   , &gendimuEta   , "gendimuEta/D");
    // tree_->Branch("gendimuPhi"   , &gendimuPhi   , "gendimuPhi/D");
    // tree_->Branch("genPhipt"     , &genPhipt     , "genPhipt/D");
    // tree_->Branch("genPhieta"    , &genPhieta    , "genPhieta/D");
    // tree_->Branch("genPhiphi"    , &genPhiphi    , "genPhiphi/D");
    // tree_->Branch("genKpPt"     , &genKpPt     , "genKpPt/D");
    // tree_->Branch("genKpEta"    , &genKpEta    , "genKpEta/D");
    // tree_->Branch("genKmPt"     , &genKmPt     , "genKmPt/D");
    // tree_->Branch("genKmEta"    , &genKmEta    , "genKmEta/D");
    
    
    // tree_->Branch("genQ2"        , &genQ2        , "genQ2/D");
    // tree_->Branch("genCosThetaL" , &genCosThetaL , "genCosThetaL/D");
    // tree_->Branch("genCosThetaK" , &genCosThetaK , "genCosThetaK/D");



    gentree_->Branch("genBpid"      , &genBpid      , "genBpid/D");

    gentree_->Branch("genBPhi"      , &genBPhi      , "genBPhi/D");
    gentree_->Branch("genMupPt"     , &genMupPt     , "genMupPt/D");
    gentree_->Branch("genMupEta"    , &genMupEta    , "genMupEta/D");
    gentree_->Branch("genMupPhi"    , &genMupPhi    , "genMupPhi/D");
    gentree_->Branch("genMumPt"     , &genMumPt     , "genMumPt/D");
    gentree_->Branch("genMumEta"    , &genMumEta    , "genMumEta/D");
    gentree_->Branch("genMumPhi"    , &genMumPhi    , "genMumPhi/D");
    
    gentree_->Branch("gendimuPt"    , &gendimuPt    , "gendimuPt/D");
    gentree_->Branch("gendimuEta"   , &gendimuEta   , "gendimuEta/D");
    gentree_->Branch("gendimuPhi"   , &gendimuPhi   , "gendimuPhi/D");
    gentree_->Branch("genPhipt"     , &genPhipt     , "genPhipt/D");
    gentree_->Branch("genPhieta"    , &genPhieta    , "genPhieta/D");
    gentree_->Branch("genPhiphi"    , &genPhiphi    , "genPhiphi/D");
    gentree_->Branch("genKpPt"     , &genKpPt     , "genKpPt/D");
    gentree_->Branch("genKpEta"    , &genKpEta    , "genKpEta/D");
    gentree_->Branch("genKmPt"     , &genKmPt     , "genKmPt/D");
    gentree_->Branch("genKmEta"    , &genKmEta    , "genKmEta/D");

    
    gentree_->Branch("genQ2"        , &genQ2        , "genQ2/D");
    gentree_->Branch("genCosThetaL" , &genCosThetaL , "genCosThetaL/D");
    gentree_->Branch("genCosThetaK" , &genCosThetaK , "genCosThetaK/D");
    break;
    
  case 999:
    tree_->Branch("genBpid"      , &genBpid      , "genBpid/D");
    tree_->Branch("genBPt"       , &genBPt       , "genBPt/D");
    tree_->Branch("genBEta"      , &genBEta      , "genBEta/D");
    tree_->Branch("genBPhi"      , &genBPhi      , "genBPhi/D");
    tree_->Branch("genBVtxX"     , &genBVtxX     , "genBVtxX/D");
    tree_->Branch("genBVtxY"     , &genBVtxY     , "genBVtxY/D");
    tree_->Branch("genBVtxZ"     , &genBVtxZ     , "genBVtxZ/D");
    tree_->Branch("genMupPt"     , &genMupPt     , "genMupPt/D");
    tree_->Branch("genMupEta"    , &genMupEta    , "genMupEta/D");
    tree_->Branch("genMupPhi"    , &genMupPhi    , "genMupPhi/D");
    tree_->Branch("genMumPt"     , &genMumPt     , "genMumPt/D");
    tree_->Branch("genMumEta"    , &genMumEta    , "genMumEta/D");
    tree_->Branch("genMumPhi"    , &genMumPhi    , "genMumPhi/D");
    
    tree_->Branch("gendimuPt"    , &gendimuPt    , "gendimuPt/D");
    tree_->Branch("gendimuEta"   , &gendimuEta   , "gendimuEta/D");
    tree_->Branch("gendimuPhi"   , &gendimuPhi   , "gendimuPhi/D");
    tree_->Branch("genPhipt"     , &genPhipt     , "genPhipt/D");
    tree_->Branch("genPhieta"    , &genPhieta    , "genPhieta/D");
    tree_->Branch("genPhiphi"    , &genPhiphi    , "genPhiphi/D");
    
    tree_->Branch("genKpPt"     , &genKpPt     , "genKpPt/D");
    tree_->Branch("genKpEta"    , &genKpEta    , "genKpEta/D");
    tree_->Branch("genKpPhi"    , &genKpPhi    , "genKpPhi/D");
    tree_->Branch("genKmPt"     , &genKmPt     , "genKmPt/D");
    tree_->Branch("genKmEta"    , &genKmEta    , "genKmEta/D");
    tree_->Branch("genKmPhi"    , &genKmPhi    , "genKmPhi/D");
    tree_->Branch("genQ2"        , &genQ2        , "genQ2/D");
    tree_->Branch("genCosThetaL" , &genCosThetaL , "genCosThetaL/D");
    tree_->Branch("genCosThetaK" , &genCosThetaK , "genCosThetaK/D");
    tree_->Branch("genPhi"       , &genPhi       , "genPhi/D");

    gentree_->Branch("genBpid"      , &genBpid      , "genBpid/D");
    gentree_->Branch("genBPt"       , &genBPt       , "genBPt/D");
    gentree_->Branch("genBEta"      , &genBEta      , "genBEta/D");
    gentree_->Branch("genBPhi"      , &genBPhi      , "genBPhi/D");
    gentree_->Branch("genBVtxX"     , &genBVtxX     , "genBVtxX/D");
    gentree_->Branch("genBVtxY"     , &genBVtxY     , "genBVtxY/D");
    gentree_->Branch("genBVtxZ"     , &genBVtxZ     , "genBVtxZ/D");
    gentree_->Branch("genMupPt"     , &genMupPt     , "genMupPt/D");
    gentree_->Branch("genMupEta"    , &genMupEta    , "genMupEta/D");
    gentree_->Branch("genMupPhi"    , &genMupPhi    , "genMupPhi/D");
    gentree_->Branch("genMumPt"     , &genMumPt     , "genMumPt/D");
    gentree_->Branch("genMumEta"    , &genMumEta    , "genMumEta/D");
    gentree_->Branch("genMumPhi"    , &genMumPhi    , "genMumPhi/D");
    
    gentree_->Branch("gendimuPt"    , &gendimuPt    , "gendimuPt/D");
    gentree_->Branch("gendimuEta"   , &gendimuEta   , "gendimuEta/D");
    gentree_->Branch("gendimuPhi"   , &gendimuPhi   , "gendimuPhi/D");
    gentree_->Branch("genPhipt"     , &genPhipt     , "genPhipt/D");
    gentree_->Branch("genPhieta"    , &genPhieta    , "genPhieta/D");
    gentree_->Branch("genPhiphi"    , &genPhiphi    , "genPhiphi/D");
    
    gentree_->Branch("genKpPt"     , &genKpPt     , "genKpPt/D");
    gentree_->Branch("genKpEta"    , &genKpEta    , "genKpEta/D");
    gentree_->Branch("genKpPhi"    , &genKpPhi    , "genKpPhi/D");
    gentree_->Branch("genKmPt"     , &genKmPt     , "genKmPt/D");
    gentree_->Branch("genKmEta"    , &genKmEta    , "genKmEta/D");
    gentree_->Branch("genKmPhi"    , &genKmPhi    , "genKmPhi/D");
    gentree_->Branch("genQ2"        , &genQ2        , "genQ2/D");
    gentree_->Branch("genCosThetaL" , &genCosThetaL , "genCosThetaL/D");
    gentree_->Branch("genCosThetaK" , &genCosThetaK , "genCosThetaK/D");
    gentree_->Branch("genPhi"       , &genPhi       , "genPhi/D");
    break;
    
  case 998:
    break;
    
  case 997:
    break;
    
  default:
    printf("No compatible datatype found. Please check use following types...\n\t\t[");
    for (std::map<string,int>::iterator iType = maptype.begin(); iType != maptype.end(); iType++){
      if (iType->second != 0) printf("%s,",iType->first.c_str());
    }
    printf("]\n");
    break;
  }
  
  fOutput->AddAll(gDirectory->GetList()); 
  
}

Bool_t SingleBsToPhiMuMuSelector::Process(Long64_t entry)
{

  ClearEvent();

  string option = GetOption();
  string datatype = get_option_value(option, "datatype"); 
  string cut = get_option_value(option, "cut"); 
  string spec_data = get_option_value(option, "spec_data");
  GetEntry(entry); 
  if(entry%100000==0)std::cout << "Processing "<<entry<<" entries"<<" Run "<<run<<" event "<<event<< std::endl;
  n_processed_ += 1; 
  Nb = nb; 
  Npv = nprivtx;
  Puw8 = fpuw8;

  //spec_data

  if (datatype != "data") SaveGen();
  //cout<<"Event number: "<<event<<"\t number of B "<<nb<<endl;

  if(cut == "cut0"){
    if ( (datatype == "data" || datatype=="mc.lite")){
      n_selected_ += 1;
      for (int im = 0; im< nb; im++) {	
  	if ( ! HasGoodDimuon(im) ) continue;
  	n_passMuonID_bdtbkg_++;
  	if ( ! HasGoodPreselection(im) ) continue;
  	n_passPresel_bdtbkg_++;
  	//cout<<"BestB vertex "<<im<<"\t"<<bvtxcl->at(im)<<endl;
  	//if(datatype == "data" || datatype=="mc.lite"){
	if(datatype == "data" || istruebs->at(im)){
  	  n_total_bdtbkg_++;
  	  SaveEvent(im, spec_data);
  	  tree_->Fill();  
  	}
      }

    }
  }
  else if(cut == "cut_bdt" || cut == "cutopt"){
    if (datatype != "data") SaveGen();
    if (datatype != "data") gentree_->Fill();    
    int i = SelectB(cut);
    if ( i != -1 ) n_passBestB_++;
    //if ( i != -1 && (datatype == "data" || datatype=="mc.lite")){//for background
    if ( i != -1 && (datatype == "data" || istruebs->at(i))){
      n_selected_ += 1;      
      SaveEvent(i, spec_data);
      //SaveRecoGen(i);
      tree_->Fill();      
    }
  }
  
  if(cut=="genonly")gentree_->Fill();
  return kTRUE;

}

void SingleBsToPhiMuMuSelector::SlaveTerminate()
{

}

void SingleBsToPhiMuMuSelector::Terminate()
{

  string option = GetOption();
  TString outfile = get_option_value(option, "ofile"); 
  //string outfile = get_option_value(option, "ofile"); 
  //printf("option=%s\n",option.c_str());
  //printf("outfile=%s",outfile.Data());
    
  TFile file(outfile.Data(), "recreate"); 
  //TFile file(outfile.c_str(), "recreate"); 
  fOutput->Write();

  t_now_.Set(); 
  printf(" \n ---------- End Job ---------- \n" ) ;
  t_now_.Print();  
  printf(" processed: %i \n selected: %i \n \
           duration: %i sec \n rate: %g evts/sec\n",
	 n_processed_, n_selected_, 
	 t_now_.Convert() - t_begin_.Convert(), 
	 float(n_processed_)/(t_now_.Convert()-t_begin_.Convert()) );

  printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  
  printf("total events = %i\n", n_total_);
  printf("#evts passsing muonID = %i\n", n_passMuonID_);
  printf("muonID efficiency = %i/%i = %.2f  \n", n_passMuonID_, n_total_, 100*float(n_passMuonID_)/float(n_total_));
  //printf("#evts passing selection cuts = %i\n", n_passSelCut_);
  //printf("#cands/evt = %i/%i = %.2f \n", n_passSelCut_, n_total_, 100*float(n_passSelCut_)/float(n_total_));
  printf("#evts passing best B sel = %i \n", n_passBestB_);
  printf("sel. efficiency = %i/%i = %.2f  \n", n_selected_, n_processed_, 100*float(n_selected_)/float(n_processed_));

  printf("++++++++++++++++++++++++++++++++++\n");
  printf("Number of event passed goodmuon and bestvtxcl:%i \n ", n_passBestB_bkg);
  printf("Further check total event = %i\n",n_total_bdtbkg_);
  printf("Further check total event after muonid = %i  and Preselcetion = %i \n",n_passMuonID_bdtbkg_, n_passPresel_bdtbkg_);

}

    
int SingleBsToPhiMuMuSelector::SelectB(string cut)
{//{{{

  int best_idx = -1; 
  double best_bvtxcl = 0.0; 
  double best_bdt=-99.0;

  if (cut == "cut_bdt") {
    for (int i = 0; i < nb; i++) {

      n_total_++;
      if ( HasGoodDimuon(i) ) goodMuon = 1;
      if ( ! HasGoodDimuon(i) ) continue; 
      n_passMuonID_++;
      
      if ( HasGoodPreselection(i) ) goodPresel = 1;
      if ( ! HasGoodPreselection(i) ) continue;      
      //cout<<"BestB vertex "<<i<<"\t"<<bvtxcl->at(i)<<endl;
      InputVariables varList;
      
      double Kmpt_t = float(sqrt( ((kmpx->at(i))*(kmpx->at(i))) + ((kmpy->at(i))*(kmpy->at(i))) ));
      double Kppt_t = float(sqrt( ((kppx->at(i))*(kppx->at(i))) + ((kppy->at(i))*(kppy->at(i))) ));
      //double Mumpt_t = float(sqrt( ((mumpx->at(i))*(mumpx->at(i))) + ((mumpy->at(i))*(mumpy->at(i))) ));
      //double Muppt_t = float(sqrt( ((muppx->at(i))*(muppx->at(i))) + ((muppy->at(i))*(muppy->at(i))) ));
      double Kmtrkdcasigbs_t = float((kmtrkdcabs->at(i)/kmtrkdcabserr->at(i)));
      double Kptrkdcasigbs_t = float((kptrkdcabs->at(i)/kptrkdcabserr->at(i)));

      double MumMinIP_t = mumMinIP->at(i);
      double MupMinIP_t = mupMinIP->at(i);
      double MumMinIPE_t = mumMinIPE->at(i);
      double MupMinIPE_t = mupMinIPE->at(i);

      double KmtrkMinIP_t = kmtrkMinIP->at(i);
      double KmtrkMinIPE_t = kmtrkMinIPE->at(i);
      double KptrkMinIPE_t = kptrkMinIPE->at(i);
      double KptrkMinIP_t = kptrkMinIP->at(i);

      //double KptrkMinIP_t = float(KptrkMinIP, KmtrkMinIP);
      float sumBmpt_ =0;
      for(unsigned int kl=0; kl<bmassIsodR->at(i).size();kl++){
	if(bmassIsodR->at(i).at(kl)<0.7){
	  sumBmpt_ += bmassIsoPt->at(i).at(kl);
	}
      }
      TLorentzVector B_4vec_m;
      B_4vec_m.SetXYZM(bpx->at(i),bpy->at(i),bpz->at(i),bmass->at(i));
      double   Bpt_m = B_4vec_m.Pt();
      double BsIso_m = Bpt_m/(sumBmpt_ + Bpt_m);

      float sumkppt_ =0;
  for(unsigned int kl=0; kl<kptrkIsodR->at(i).size();kl++){
    if(kptrkIsodR->at(i).at(kl)<0.5){
      sumkppt_ += kptrkIsoPt->at(i).at(kl);
    }
  }
     double KptrkIso_ = rawTrkppt->at(i)/(sumkppt_ + rawTrkppt->at(i));

     float sumkmpt_ =0;
  for(unsigned int kl=0; kl<kmtrkIsodR->at(i).size();kl++){
    if(kmtrkIsodR->at(i).at(kl)<0.5){
      sumkmpt_ += kmtrkIsoPt->at(i).at(kl);
    }
  }
     double KmtrkIso_ = rawTrkmpt->at(i)/(sumkmpt_ + rawTrkmpt->at(i));

      varList.Max_Kpt_ = TMath::Max(Kmpt_t, Kppt_t);
      varList.Max_MuMinIPsig_ = TMath::Max(MumMinIP_t/MumMinIPE_t, MupMinIP_t/MupMinIPE_t);
      varList.Max_MinIPsig_  = TMath::Max(KmtrkMinIP_t/KmtrkMinIPE_t, KptrkMinIP_t/KptrkMinIPE_t);
      varList.Max_DCA_ = TMath::Max(Kmtrkdcasigbs_t, Kptrkdcasigbs_t);
      varList.Bcosalphabs2d_ = float(bcosalphabs2d->at(i));
      varList.Blxysig_  = float((blsbs->at(i)/blsbserr->at(i)));
      varList.Bvtxcl_ = float(bvtxcl->at(i));
      varList.Bpt_  = Bpt_m;
      varList.Bsdcasigbs_= fabs( bdcabs->at(i)/bdcabserr->at(i) );
   //   varList.Phimass_ = phimass->at(i);
      varList.BsIso_= BsIso_m;
      varList.K_Iso_         	 = TMath::Max(KmtrkIso_, KptrkIso_);

      //spectator variable
      varList.Bmass_         = bmass->at(i);
      varList.Mumumass_      = mumumass->at(i);
      varList.Mumumasserr_   = mumumasserr->at(i);
      
      // varList.Max_Kpt_ = TMath::Max(Kmpt_t, Kppt_t); 
      // varList.Max_trk_ = TMath::Max(Kmtrkdcasigbs_t, Kptrkdcasigbs_t);  
      // varList.Max_Mpt_ = TMath::Max(Mumpt_t, Muppt_t); 
      //varList.Max_K_MinIP_ = TMath::Max(KmtrkMinIP, KptrkMinIP);
      // varList.bmass_ = bmass->at(i);
      // varList.mumumass_ = mumumass->at(i);
      // varList.mumumasserr_ = mumumasserr->at(i);
      // varList.kmpt_ = Kmpt_t;
      // varList.kppt_ = Kppt_t;
      // varList.kmtkdca_ = Kmtrkdcasigbs_t;
      // varList.kptkdca_ = Kptrkdcasigbs_t;
      if(kmtrkMinIP->at(i)/kmtrkMinIPE->at(i) >0  && kptrkMinIP->at(i)/kptrkMinIPE->at(i) >0  && BsIso_m>0 && fabs( bdcabs->at(i)/bdcabserr->at(i))>0){
	  double mvareader = -99.0;
	  mvareader=mvAna_->evaluate(mvAlgo_, varList);
	  //cout<<"BDT: "<<mvareader<<endl;
	  if (mvareader > best_bdt) {
	    best_bdt = mvareader; 
	    //cout<<"best_bdt inside "<<best_bdt<<endl;
	    Bdt=best_bdt;
	    n_passBestB_++;
	    n_passBestB_bkg++;
	    best_idx = i; 
	  }
	}
	}
    }else if (cut == "cutopt") {
      
      for (int i = 0; i< nb; i++) {
      n_total_++;

      if ( ! HasGoodDimuon(i) ) continue;
      n_passMuonID_++;

      if ( HasGoodPreselection(i) ) goodPresel = 1;
      if ( ! HasGoodPreselection(i) ) continue;

      if (bvtxcl->at(i) > best_bvtxcl) {
	best_bvtxcl = bvtxcl->at(i);
	n_passBestB_++;
	n_passBestB_bkg++;
        best_idx = i;
      }
    }

  }else if (cut == "nocut") {
    for (int i = 0; i < nb; i++) {
      if (bvtxcl->at(i) > best_bvtxcl) {
	best_bvtxcl = bvtxcl->at(i);
	n_passBestB_++; 
	best_idx = i; 
      }
    }
  }else if (cut == "genonly") {
    best_idx = -1;
  }else{
    printf("WARNING: Unknown cut, apply 'genonly' by default.\n");
    best_idx = -1;
  }

  return best_idx;
}//}}}

bool SingleBsToPhiMuMuSelector::HasGoodDimuon(int i)
{//{{{
 
  if ( // new soft muon id
      mumisgoodmuon->at(i)
      && mupisgoodmuon->at(i) 
      && mumntrklayers->at(i) > 5  
      && mupntrklayers->at(i) > 5   
      && mumnpixlayers->at(i) > 0  
      && mupnpixlayers->at(i) > 0   

      && mumtrkqual->at(i)==1
      && muptrkqual->at(i)==1
      && fabs(mumdxyvtx->at(i)) < 0.3 
      && fabs(mupdxyvtx->at(i)) < 0.3 
      && fabs(mumdzvtx->at(i)) < 20    
      && fabs(mupdzvtx->at(i)) < 20   
      
      //&& mumloosemuon->at(i)
      //&& muploosemuon->at(i)

       ) return true; 
  return false; 
}//}}}


bool SingleBsToPhiMuMuSelector::HasGoodPreselection(int i)
{//{{{                                                                                                                                                                         
  double kmpt = sqrt( (kmpx->at(i))*(kmpx->at(i)) + (kmpy->at(i))*(kmpy->at(i)) );
  double kppt = sqrt( (kppx->at(i))*(kppx->at(i)) + (kppy->at(i))*(kppy->at(i)) );

  if ( bcosalphabs2d->at(i)>0.9 && bvtxcl->at(i)>0.01 && (fabs(kmtrkdcabs->at(i)/kmtrkdcabserr->at(i))>0.8) && (fabs(kptrkdcabs->at(i)/kptrkdcabserr->at(i))>0.8) && kmpt>0.8 && kppt>0.8) return true;

  return false;
}//}}}        

void SingleBsToPhiMuMuSelector::SaveEvent(int i, string spec)
{//{{{

  TLorentzVector B_4vec, Phi_4vec, Mup_4vec, Mum_4vec, Km_4vec, Kp_4vec, buff1, buff2, buff3;
  B_4vec.SetXYZM(bpx->at(i),bpy->at(i),bpz->at(i),bmass->at(i));
  Phi_4vec.SetXYZM(kmpx->at(i)+kppx->at(i),kmpy->at(i)+kppy->at(i),kmpz->at(i)+kppz->at(i),phimass->at(i));
  Mup_4vec.SetXYZM(muppx->at(i),muppy->at(i),muppz->at(i),MUON_MASS);
  Mum_4vec.SetXYZM(mumpx->at(i),mumpy->at(i),mumpz->at(i),MUON_MASS);
  Km_4vec.SetXYZM(kmpx->at(i),kmpy->at(i),kmpz->at(i),KAON_MASS);
  Kp_4vec.SetXYZM(kppx->at(i),kppy->at(i),kppz->at(i),KAON_MASS);

  Bmass = bmass->at(i); 
  Bvtxcl = bvtxcl->at(i); 
  Blxysig = (blsbs->at(i)/blsbserr->at(i)); 
  Bcosalphabs = bcosalphabs->at(i); 
  Bcosalphabs2d = bcosalphabs2d->at(i);
  Bctau = bctau->at(i); 

  Phipt  = Phi_4vec.Pt();
  Phieta = Phi_4vec.Eta();
  Phiphi = Phi_4vec.Phi();

  Mumpt = Mum_4vec.Pt();
  Mumeta = Mum_4vec.Eta();
  Mumphi = Mum_4vec.Phi();
  Muppt = Mup_4vec.Pt();
  Mupeta = Mup_4vec.Eta();
  Mupphi = Mup_4vec.Phi();

  Kmpt = sqrt( (kmpx->at(i))*(kmpx->at(i)) + (kmpy->at(i))*(kmpy->at(i)) );
  Kppt = sqrt( (kppx->at(i))*(kppx->at(i)) + (kppy->at(i))*(kppy->at(i)) );
  Kmeta = Km_4vec.Eta();
  Kmphi = Km_4vec.Phi();
  Kppt = Kp_4vec.Pt();
  Kpeta = Kp_4vec.Eta();
  Kpphi = Kp_4vec.Phi();
  TLorentzVector B0_4vec, B0s_4vec,Mup_4vec1, Mum_4vec1, Km_4vec1,  Kp_4vec1, pim_4vec1, pip_4vec1;
  TLorentzVector LB_4vec, Km_4vecl,  pp_4vecl, pm_4vecl, LBs_4vec, Kp_4vecl;

  Mup_4vec1.SetXYZM(muppx->at(i),muppy->at(i),muppz->at(i),MUON_MASS);
  Mum_4vec1.SetXYZM(mumpx->at(i),mumpy->at(i),mumpz->at(i),MUON_MASS);
  pim_4vec1.SetXYZM(kmpx->at(i),kmpy->at(i),kmpz->at(i),PION_MASS);
  Kp_4vec1.SetXYZM(kppx->at(i),kppy->at(i),kppz->at(i),KAON_MASS);

  Km_4vec1.SetXYZM(kmpx->at(i),kmpy->at(i),kmpz->at(i),KAON_MASS);
  pip_4vec1.SetXYZM(kppx->at(i),kppy->at(i),kppz->at(i),PION_MASS);

  Km_4vecl.SetXYZM(kmpx->at(i),kmpy->at(i),kmpz->at(i),KAON_MASS);
  Kp_4vecl.SetXYZM(kppx->at(i),kppy->at(i),kppz->at(i),KAON_MASS);
  pp_4vecl.SetXYZM(kppx->at(i),kppy->at(i),kppz->at(i),PROTON_MASS);
  pm_4vecl.SetXYZM(kmpx->at(i),kmpy->at(i),kmpz->at(i),PROTON_MASS);
  
  B0_4vec = Mup_4vec1 + Mum_4vec1 + pim_4vec1 + Kp_4vec1;
  B0s_4vec = Mup_4vec1 + Mum_4vec1 + pip_4vec1 + Km_4vec1;
  LB_4vec = Mup_4vec1 + Mum_4vec1 + pp_4vecl + Km_4vecl;
  LBs_4vec = Mup_4vec1 + Mum_4vec1 + pm_4vecl + Kp_4vecl;

  TLorentzVector Bs1_4vec, Bs2_4vec,Mup_S_4vec1, Mum_S_4vec1, Km_S_4vec1,  Kp_S_4vec1;
  Mup_S_4vec1.SetXYZM(muppx->at(i),muppy->at(i),muppz->at(i),KAON_MASS);
  Mum_S_4vec1.SetXYZM(mumpx->at(i),mumpy->at(i),mumpz->at(i),KAON_MASS);
  Km_S_4vec1.SetXYZM(kmpx->at(i),kmpy->at(i),kmpz->at(i),MUON_MASS);
  Kp_S_4vec1.SetXYZM(kppx->at(i),kppy->at(i),kppz->at(i),MUON_MASS);

  Bs1_4vec = Mup_4vec1+Mum_S_4vec1+Kp_4vec1+Km_S_4vec1;
  Bs2_4vec = Mup_S_4vec1+Mum_4vec1+Kp_S_4vec1+Km_4vec1;

  TLorentzVector Bs3_4vec= Mup_S_4vec1 +Mum_S_4vec1+Kp_S_4vec1+Km_S_4vec1;
  Bmass_kmu1 = Bs1_4vec.M();
  Bmass_kmu2 = Bs2_4vec.M();
  Bmass_kkmumu = Bs3_4vec.M();
  Bmass_mumu_s = (Mup_S_4vec1 +Mum_S_4vec1).M();
  Bmass_kk_s = (Kp_S_4vec1+Km_S_4vec1).M();

  Bmass_kk1_s = (Kp_4vec1+Mum_S_4vec1).M();
  Bmass_kk2_s = (Mup_S_4vec1+Km_4vec1).M();
  Bmass_mumu1_s = (Kp_S_4vec1+Mum_4vec1).M();
  Bmass_mumu2_s = (Km_S_4vec1+Mup_4vec1).M();
  
  TLorentzVector L_4vec, Kst0_4vecl, Kst_4vecl , Ls_4vec;

  L_4vec = pp_4vecl + Km_4vecl;
  Ls_4vec = pm_4vecl + Kp_4vecl;
  Kst0_4vecl =  pip_4vec1 + Km_4vec1;
  Kst_4vecl = pim_4vec1 + Kp_4vec1;
  Lmass = L_4vec.M();
  Lsmass = Ls_4vec.M();
  
  Kst0mass = Kst0_4vecl.M();
  Kstmass = Kst_4vecl.M();
  BOmass = B0_4vec.M();
  BOSmass = B0s_4vec.M();
  LBmass = LB_4vec.M();
  LBsmass = LBs_4vec.M();
  
  TLorentzVector pht_4vec, plt_4vec;
  TLorentzVector kht_4vec, klt_4vec;
  TLorentzVector piht_4vec;

  if(Kp_4vec.P()>Km_4vec.P()){
    pht_4vec.SetXYZM(kppx->at(i),kppy->at(i),kppz->at(i),PROTON_MASS);
    kht_4vec.SetXYZM(kmpx->at(i),kmpy->at(i),kmpz->at(i),KAON_MASS);
    piht_4vec.SetXYZM(kmpx->at(i),kmpy->at(i),kmpz->at(i),PION_MASS);
    plt_4vec.SetXYZM(kppx->at(i),kppy->at(i),kppz->at(i),KAON_MASS);
    klt_4vec.SetXYZM(kmpx->at(i),kmpy->at(i),kmpz->at(i),PROTON_MASS);
        
  } else {
    pht_4vec.SetXYZM(kmpx->at(i),kmpy->at(i),kmpz->at(i),PROTON_MASS);
    kht_4vec.SetXYZM(kppx->at(i),kppy->at(i),kppz->at(i),KAON_MASS);
    piht_4vec.SetXYZM(kppx->at(i),kppy->at(i),kppz->at(i),PION_MASS);
    plt_4vec.SetXYZM(kmpx->at(i),kmpy->at(i),kmpz->at(i),KAON_MASS);
    klt_4vec.SetXYZM(kppx->at(i),kppy->at(i),kppz->at(i),PROTON_MASS);
  }
  TLorentzVector Lht_4vec, Llt_4vec,Lpiht_4vec;
  TLorentzVector LBht_4vec, LBlt_4vec, LBpiht_4vec;
  Lht_4vec = pht_4vec+kht_4vec;
  Lhtmass = Lht_4vec.M();
  LBht_4vec = Mup_4vec1 + Mum_4vec1 + pht_4vec + kht_4vec;
  LBhtmass =LBht_4vec.M();

  Lpiht_4vec = pht_4vec+piht_4vec;
  Lpihtmass = Lpiht_4vec.M();
  LBpiht_4vec = Mup_4vec1 + Mum_4vec1 + pht_4vec + piht_4vec;
  LBpihtmass =LBpiht_4vec.M();


  LBlt_4vec = Mup_4vec1 + Mum_4vec1 + plt_4vec + klt_4vec;
  LBltmass =LBlt_4vec.M();

  Llt_4vec = plt_4vec+klt_4vec;
  Lltmass = Llt_4vec.M();
  //pt
  TLorentzVector phtpt_4vec,khtpt_4vec;
  if(Kppt > Kmpt){
    phtpt_4vec.SetXYZM(kppx->at(i),kppy->at(i),kppz->at(i),PROTON_MASS);
    khtpt_4vec.SetXYZM(kmpx->at(i),kmpy->at(i),kmpz->at(i),KAON_MASS);

  } else {
    phtpt_4vec.SetXYZM(kmpx->at(i),kmpy->at(i),kmpz->at(i),PROTON_MASS);
    khtpt_4vec.SetXYZM(kppx->at(i),kppy->at(i),kppz->at(i),KAON_MASS);

  }
  
  TLorentzVector Lhtpt_4vec;
  TLorentzVector LBhtpt_4vec;
  Lhtpt_4vec = phtpt_4vec+khtpt_4vec;
  Lhtptmass = Lhtpt_4vec.M();
  LBhtpt_4vec = Mup_4vec1 + Mum_4vec1 + phtpt_4vec + khtpt_4vec;
  LBhtptmass =LBhtpt_4vec.M();

  //std::cout<<"Bmass 4vector "<<B_4vec.M()<<" fromntuple "<<Bmass<<" caculating here "<<B0_4vec.M()<<" b0s "<<B0s_4vec.M()<<" lambdaB "<<LB_4vec.M()<<std::endl;  
  
  Bpt  = B_4vec.Pt(); 
  Beta = B_4vec.Eta();
  Bphi = B_4vec.Phi();

  MumMinIP = mumMinIP->at(i);
  MupMinIP = mupMinIP->at(i);
  MumMinIPE = mumMinIPE->at(i);
  MupMinIPE = mupMinIPE->at(i);
  KptrkMinIP = kptrkMinIP->at(i);
  KptrkMinIPE = kptrkMinIPE->at(i);
  KmtrkMinIP = kmtrkMinIP->at(i);
  KmtrkMinIPE = kmtrkMinIPE->at(i);

  MumMinIP2D = mumMinIP2D->at(i);
  MupMinIP2D = mupMinIP2D->at(i);
  MumMinIP2DE = mumMinIP2DE->at(i);
  MupMinIP2DE = mupMinIP2DE->at(i);
  KptrkMinIP2D = kptrkMinIP2D->at(i);
  KptrkMinIP2DE = kptrkMinIP2DE->at(i);
  KmtrkMinIP2D = kmtrkMinIP2D->at(i);
  KmtrkMinIP2DE = kmtrkMinIP2DE->at(i);


  Mumumass    = mumumass->at(i); 
  Mumumasserr = mumumasserr->at(i); 
  Phimass     = phimass->at(i); 

  Kmtrkdcasigbs = fabs( kmtrkdcabs->at(i)/kmtrkdcabserr->at(i) ); 
  Kptrkdcasigbs = fabs( kptrkdcabs->at(i)/kptrkdcabserr->at(i) ); 

  Mumdcasigbs = fabs( mumdcabs->at(i)/mumdcabserr->at(i) ); 
  Mupdcasigbs = fabs( mupdcabs->at(i)/mupdcabserr->at(i) ); 
  Bsdcasigbs  = fabs( bdcabs->at(i)/bdcabserr->at(i) );
  Q2 = pow(mumumass->at(i),2);
  dimuvtxcl = mumuvtxcl->at(i);
  dimulsig  = (mumulsbs->at(i))/(mumulsbserr->at(i));
  dimucosalphabs = mumucosalphabs->at(i);
  dimulensig = (mumulsbs->at(i))/(mumulsbserr->at(i));
  dimuDCA = mumudca->at(i);

  cosdimuon = (Mup_4vec.Vect().Dot(Mum_4vec.Vect()))/(Mup_4vec.Vect().Mag()*Mum_4vec.Vect().Mag());
  cosphimup = (Phi_4vec.Vect().Dot(Mup_4vec.Vect()))/(Phi_4vec.Vect().Mag()*Mup_4vec.Vect().Mag());
  cosphimum = (Phi_4vec.Vect().Dot(Mum_4vec.Vect()))/(Phi_4vec.Vect().Mag()*Mum_4vec.Vect().Mag());
  pmum_mass = (Phi_4vec+Mum_4vec).M();
  pmup_mass = (Phi_4vec+Mup_4vec).M();

  buff1 = B_4vec;
  buff2 = Mup_4vec+Mum_4vec;

  dimupt  = buff2.Pt();
  dimueta = buff2.Eta();

  buff1.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());

  // review this part again ??
  //  if ( Bchg > 0){
  buff3 = Mum_4vec; 
    //  }else{
    //buff3 = Mup_4vec;
    //  }
  buff3.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
  CosThetaL = buff1.Vect().Dot(buff3.Vect())/buff1.Vect().Mag()/buff3.Vect().Mag();
    
  buff1 = B_4vec;
  buff2 = Phi_4vec;
  buff1.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
  buff3 = Km_4vec; // double-check 
  buff3.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
  CosThetaK = buff1.Vect().Dot(buff3.Vect())/buff1.Vect().Mag()/buff3.Vect().Mag();

  // add phi angle
  TVector3 boostB = B_4vec.BoostVector();
  Mum_4vec.Boost(-boostB);
  Mup_4vec.Boost(-boostB);
  Km_4vec.Boost(-boostB);
  Kp_4vec.Boost(-boostB);
  TVector3 MuMuPlane = Mum_4vec.Vect().Cross(Mup_4vec.Vect());   /// cross product between mu- and mu+ vectors                                       
  /////cout << "cross product mag. of muons (at RECO level) = " << MuMuPlane.Mag() << endl;   /// print statement for x-check
  TVector3 PhiPlane = Kp_4vec.Vect().Cross(Km_4vec.Vect());
  /////cout << "cross product mag. of kaons (at RECO level) = " << PhiPlane.Mag() << endl;    /// print statement for x-check
  if (MuMuPlane.Cross(PhiPlane).Dot(-B_4vec.Vect()) > 0.0)
    Phi = MuMuPlane.Angle(PhiPlane);
  else
    Phi = -MuMuPlane.Angle(PhiPlane);
  mtrkqual   = kmtrkqual->at(i);
  ptrkqual   = kptrkqual->at(i);

  // no implementation
  // JpsiTriggers = tri_JpsiTk->at(i);
  // PsiPTriggers = tri_PsipTk->at(i);
  // LMNTTriggers = tri_LMNTk->at(i);

  // implement
  Int_t JpsiTrig_t=0.;
  Int_t PsiPTrig_t=0.;
  Int_t LMNTTrig_t=0.;
    
  JpsiTrig_t = tri_JpsiTk->at(i);
  PsiPTrig_t = tri_PsipTk->at(i);
  LMNTTrig_t = tri_LMNTk->at(i);


  if(spec == "Charmonium"){
    LMNTTrig_t =0;

  }
  if(spec == "DoubleMuonLowMass"){
    JpsiTrig_t = 0;
    PsiPTrig_t = 0;
  }
  double mumu_t = mumumass->at(i);
  //This block is used to test the triggers
  //method 1
  // int nDummy = (int)(mumu_t *1000.0);
  // mumu_t = ((double)nDummy) / 1000.0;
  mumu_t = round(mumu_t*1000)/1000; //rounding up upto3rd decimal

  if(JpsiTrig_t){
    //if(mumumass->at(i)>=2.9 && mumumass->at(i)<=3.305){
    if(mumu_t >= 2.899 && mumu_t <= 3.301){
      JpsiTriggers = 1;
      LMNTTriggers = 0;
      PsiPTriggers = 0;
    }
    else {
      JpsiTriggers = 0;    
    }
  } else {
    JpsiTriggers = 0;
  }
  if(PsiPTrig_t){
    if(mumu_t >= 3.299 && mumu_t <= 4.050){
      //if(mumumass->at(i)>=3.295 && mumumass->at(i)<=4.04){
      PsiPTriggers = 1;
      JpsiTriggers = 0;
      LMNTTriggers = 0;
    } else {
      PsiPTriggers = 0;
    }
  } else {
    PsiPTriggers = 0;
  }
  
  if(LMNTTrig_t){
    if(mumu_t <= 2.901 || mumu_t >=3.999){
      //if(mumumass->at(i)<=2.9 || mumumass->at(i)>=4.03){
      LMNTTriggers = 1;
      JpsiTriggers = 0;
      PsiPTriggers = 0;      
    } else {
      LMNTTriggers = 0;
    }
  } else {
    LMNTTriggers = 0;
  }

  // Method 2
  // mumu_t = round(mumu_t*100)/100; //rounding up upto 2nd decimal

  // if(JpsiTrig_t){
  //   if(mumu_t >= 2.89 && mumu_t <= 3.31){
  //     JpsiTriggers = 1;
  //     LMNTTriggers = 0;
  //     PsiPTriggers = 0;
  //   }
  //   else {
  //     JpsiTriggers = 0;    
  //   }
  // } else {
  //   JpsiTriggers = 0;
  // }
  // if(PsiPTrig_t){
  //   if(mumu_t >= 3.29 && mumu_t <= 4.05){
  //     PsiPTriggers = 1;
  //     JpsiTriggers = 0;
  //     LMNTTriggers = 0;
  //   } else {
  //     PsiPTriggers = 0;
  //   }
  // } else {
  //   PsiPTriggers = 0;
  // }
  
  // if(LMNTTrig_t){
  //   if(mumu_t <= 2.91 || mumu_t >=3.99){
  //     LMNTTriggers = 1;
  //     JpsiTriggers = 0;
  //     PsiPTriggers = 0;      
  //   } else {
  //     LMNTTriggers = 0;
  //   }
  // } else {
  //   LMNTTriggers = 0;
  // }
  // //---------------------------------
  // if(tri_PsipTk->at(i) == 1 && tri_JpsiTk->at(i) == 1 && (mumu_t>3.28 && mumu_t<3.31))
  //   cout<<" Run "<<run<<" event "<<event<< " mumu "<<mumu_t<<" jpsi "<<JpsiTriggers<<" psip "<<PsiPTriggers<<" lmnt "<<LMNTTriggers<<endl;
  // if(tri_JpsiTk->at(i) == 1 && tri_LMNTk->at(i) == 1 && (mumu_t>2.88 && mumu_t<2.92))
  //   cout<<" Run "<<run<<" event "<<event<< " mumu "<<mumu_t<<" jpsi "<<JpsiTriggers<<" psip "<<PsiPTriggers<<" lmnt "<<LMNTTriggers<<endl;

  // if((mumu_t>2.88 && mumu_t<2.92))
  //   cout<<spec <<" Run "<<run<<" event "<<event<< " mumu "<<mumu_t<<" jpsi "<<JpsiTriggers<<" psip "<<PsiPTriggers<<" lmnt "<<LMNTTriggers<<" orgjpsi "<<tri_JpsiTk->at(i)<<" orgpsi "<<tri_PsipTk->at(i)<<" orglmnt "<<tri_LMNTk->at(i)<<endl;

  if(spec == "DoubleMuonLowMass" && tri_LMNTk->at(i) == 1 && tri_PsipTk->at(i) == 1 &&(mumu_t>3.99 && mumu_t<4.05)){
    LMNTTriggers = 0;
  }
  //if(spec == "
  dr0 = tri_dr0->at(i);
  dr1 = tri_dr1->at(i);
  dpt0 = tri_dpt0->at(i);
  dpt1 = tri_dpt1->at(i);
  
  double sumppt =0;
  for(unsigned int kl=0; kl<mupIsodR->at(i).size();kl++){
    if(mupIsodR->at(i).at(kl)<0.5){
      sumppt += mupIsoPt->at(i).at(kl);
    }
  }
  mupIso_ = rawmuppt->at(i)/(sumppt + rawmuppt->at(i));
  float summpt =0;
  for(unsigned int kl=0; kl<mumIsodR->at(i).size();kl++){
    if(mumIsodR->at(i).at(kl)<0.5){
      summpt += mumIsoPt->at(i).at(kl);
    }
  }
  mumIso_ = rawmumpt->at(i)/(summpt + rawmumpt->at(i));
  float sumkppt =0;
  for(unsigned int kl=0; kl<kptrkIsodR->at(i).size();kl++){
    if(kptrkIsodR->at(i).at(kl)<0.5){
      sumkppt += kptrkIsoPt->at(i).at(kl);
    }
  }
  kptrkIso_ = rawTrkppt->at(i)/(sumkppt + rawTrkppt->at(i));

  float sumkmpt =0;
  for(unsigned int kl=0; kl<kmtrkIsodR->at(i).size();kl++){
    if(kmtrkIsodR->at(i).at(kl)<0.5){
      sumkmpt += kmtrkIsoPt->at(i).at(kl);
    }
  }
  kmtrkIso_ = rawTrkmpt->at(i)/(sumkmpt + rawTrkmpt->at(i));
  float sumBmpt =0;
  for(unsigned int kl=0; kl<bmassIsodR->at(i).size();kl++){
    if(bmassIsodR->at(i).at(kl)<0.7){
      sumBmpt += bmassIsoPt->at(i).at(kl);
    }
  }
  BsIso_ = Bpt/(sumBmpt + Bpt);



}//}}}

// void SingleBsToPhiMuMuSelector::SaveRecoGen(int i)
// {//{{{

//   TLorentzVector genB_4vec, genPhi_4vec, genMup_4vec, genMum_4vec, genKm_4vec, genKp_4vec, buff1, buff2, buff3;
//   genB_4vec.SetXYZM(genbpx->at(i),genbpy->at(i),genbpz->at(i),5.3668);
//   genPhi_4vec.SetXYZM(genphipx->at(i),genphipy->at(i),genphipz->at(i),PHI_MASS);
//   genMup_4vec.SetXYZM(genmuppx->at(i),genmuppy->at(i),genmuppz->at(i),MUON_MASS);
//   genMum_4vec.SetXYZM(genmumpx->at(i),genmumpy->at(i),genmumpz->at(i),MUON_MASS);
//   genKm_4vec.SetXYZM(genkmpx->at(i),genkmpy->at(i),genkmpz->at(i),KAON_MASS);
//   genKp_4vec.SetXYZM(genkppx->at(i),genkppy->at(i),genkppz->at(i),KAON_MASS);
  

//   genBpid      = genbpid->at(i);
//   genBPt       = genB_4vec.Pt();
//   genBEta      = genB_4vec.Eta();
//   genBPhi      = genB_4vec.Phi();
//   genBVtxX     = 0;//Should be at PV?
//   genBVtxY     = 0;
//   genBVtxZ     = 0;

//   genPhipt  = genPhi_4vec.Pt();
//   genPhieta = genPhi_4vec.Eta();
//   genPhiphi = genPhi_4vec.Phi();

//   genKpPt       = genKp_4vec.Pt();
//   genKpEta      = genKp_4vec.Eta();
//   genKpPhi      = genKp_4vec.Phi();
//   genKmPt       = genKm_4vec.Pt();
//   genKmEta      = genKm_4vec.Eta();
//   genKmPhi      = genKm_4vec.Phi();

//   genMupPt     = genMup_4vec.Pt();
//   genMupEta    = genMup_4vec.Eta();
//   genMupPhi    = genMup_4vec.Phi();
//   genMumPt     = genMum_4vec.Pt();
//   genMumEta    = genMum_4vec.Eta();
//   genMumPhi    = genMum_4vec.Phi();

//   genQ2        = (genMup_4vec+genMum_4vec).Mag2();
    
//   buff1        = genB_4vec;
//   buff2        = genMup_4vec+genMum_4vec;

//   gendimuPt    = buff2.Pt();
//   gendimuEta   = buff2.Eta();
//   gendimuPhi   = buff2.Phi();

//   buff1.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
//   //  if (genBChg > 0){
//     buff3 = genMum_4vec;//Take mu- to avoid extra minus sign.
//     //}else{
//     //buff3 = genMup_4vec;
//     //}
//   buff3.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
//   genCosThetaL = buff1.Vect().Dot(buff3.Vect())/buff1.Vect().Mag()/buff3.Vect().Mag();
    
//   buff1 = genB_4vec;
//   buff2 = genPhi_4vec;
//   buff1.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
//   buff3 = genKm_4vec; // double check
//   buff3.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
//   genCosThetaK = buff1.Vect().Dot(buff3.Vect())/buff1.Vect().Mag()/buff3.Vect().Mag();
 
//   // add phi angle
//   TVector3 boostB = genB_4vec.BoostVector();
//   genMum_4vec.Boost(-boostB);
//   genMup_4vec.Boost(-boostB);
//   genKm_4vec.Boost(-boostB);
//   genKp_4vec.Boost(-boostB);
//   TVector3 MuMuPlane = genMum_4vec.Vect().Cross(genMup_4vec.Vect());   /// cross product between mu- and mu+ vectors
//   /////cout << "cross product mag. of muons (at GEN level) = " << MuMuPlane.Mag() << endl;  ////  print statement for x-check
//   TVector3 PhiPlane = genKp_4vec.Vect().Cross(genKm_4vec.Vect());
//   /////cout << "cross product mag. of kaons (at GEN level) = " << PhiPlane.Mag() << endl;   ////  print statement for x-check
//   if (MuMuPlane.Cross(PhiPlane).Dot(-genB_4vec.Vect()) > 0.0) 
//     genPhi = MuMuPlane.Angle(PhiPlane);
//   else                                                        
//     genPhi = -MuMuPlane.Angle(PhiPlane);

// }//}}}
void SingleBsToPhiMuMuSelector::SaveGen()
{//{{{

  TLorentzVector genB_4vec, genPhi_4vec, genMup_4vec, genMum_4vec, genKm_4vec, genKp_4vec, buff1, buff2, buff3;
  genB_4vec.SetXYZM(genbpx,genbpy,genbpz,5.3668);
  
  genMup_4vec.SetXYZM(genmuppx,genmuppy,genmuppz,MUON_MASS);
  genMum_4vec.SetXYZM(genmumpx,genmumpy,genmumpz,MUON_MASS);
  genKm_4vec.SetXYZM(genkmpx,genkmpy,genkmpz,KAON_MASS);
  genKp_4vec.SetXYZM(genkppx,genkppy,genkppz,KAON_MASS);
  
  Double_t genphipx_t = genkmpx + genkppx;
  Double_t genphipy_t = genkmpy + genkppy;
  Double_t genphipz_t = genkmpz + genkppz;

  genPhi_4vec.SetXYZM(genphipx_t,genphipy_t,genphipz_t,PHI_MASS);  


  genBpid      = genbpid;
  genBPt       = genB_4vec.Pt();
  genBEta      = genB_4vec.Eta();
  genBPhi      = genB_4vec.Phi();
  genBVtxX     = 0;//Should be at PV?
  genBVtxY     = 0;
  genBVtxZ     = 0;

  genPhipt  = genPhi_4vec.Pt();
  genPhieta = genPhi_4vec.Eta();
  genPhiphi = genPhi_4vec.Phi();

  genKpPt       = genKp_4vec.Pt();
  genKpEta      = genKp_4vec.Eta();
  genKpPhi      = genKp_4vec.Phi();
  genKmPt       = genKm_4vec.Pt();
  genKmEta      = genKm_4vec.Eta();
  genKmPhi      = genKm_4vec.Phi();

  genMupPt     = genMup_4vec.Pt();
  genMupEta    = genMup_4vec.Eta();
  genMupPhi    = genMup_4vec.Phi();
  genMumPt     = genMum_4vec.Pt();
  genMumEta    = genMum_4vec.Eta();
  genMumPhi    = genMum_4vec.Phi();

  genQ2        = (genMup_4vec+genMum_4vec).Mag2();
    
  buff1        = genB_4vec;
  buff2        = genMup_4vec+genMum_4vec;

  gendimuPt    = buff2.Pt();
  gendimuEta   = buff2.Eta();
  gendimuPhi   = buff2.Phi();

  buff1.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
  //  if (genBChg > 0){
    buff3 = genMum_4vec;//Take mu- to avoid extra minus sign.
    //}else{
    //buff3 = genMup_4vec;
    //}
  buff3.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
  genCosThetaL = buff1.Vect().Dot(buff3.Vect())/buff1.Vect().Mag()/buff3.Vect().Mag();
    
  buff1 = genB_4vec;
  buff2 = genPhi_4vec;
  buff1.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
  buff3 = genKm_4vec; // double check
  buff3.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
  genCosThetaK = buff1.Vect().Dot(buff3.Vect())/buff1.Vect().Mag()/buff3.Vect().Mag();
 
  // add phi angle
  TVector3 boostB = genB_4vec.BoostVector();
  genMum_4vec.Boost(-boostB);
  genMup_4vec.Boost(-boostB);
  genKm_4vec.Boost(-boostB);
  genKp_4vec.Boost(-boostB);
  TVector3 MuMuPlane = genMum_4vec.Vect().Cross(genMup_4vec.Vect());   /// cross product between mu- and mu+ vectors
  /////cout << "cross product mag. of muons (at GEN level) = " << MuMuPlane.Mag() << endl;  ////  print statement for x-check
  TVector3 PhiPlane = genKp_4vec.Vect().Cross(genKm_4vec.Vect());
  /////cout << "cross product mag. of kaons (at GEN level) = " << PhiPlane.Mag() << endl;   ////  print statement for x-check
  if (MuMuPlane.Cross(PhiPlane).Dot(-genB_4vec.Vect()) > 0.0) 
    genPhi = MuMuPlane.Angle(PhiPlane);
  else                                                        
    genPhi = -MuMuPlane.Angle(PhiPlane);

}//}}}


#ifndef __CINT__ 
#include <algorithm>

char* get_option(char ** begin, char ** end, const std::string & option)
{//{{{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)  return *itr;
  return 0;
}//}}}

bool option_exists(char** begin, char** end, const std::string& option)
{//{{{
  return std::find(begin, end, option) != end;
}//}}}

void print_usage()
{//{{{
  cerr << "Usage: SingleBsToPhiMuMuSelector datatype cut infile outfile [-n] [-s] [-j] [-h]\n"
       << "  datatype: data, mc, mc.lite, mc.hlt\n"
       << "  cut     : cut0, nocut, genonly.\n"
       << "Options: \n" 
       << "  -h \t\tPrint thin_passBestB_bkgs info.\n"
       << "  -n \t\tNumber of entries.\n" 
       << "  -s \t\tStarting run number.\n"
       << "  -j \t\tNumber of workers.\n" 
       << endl; 
}//}}}

int main(int argc, char** argv) {
  if ( (argc < 3) or option_exists(argv, argv+argc, "-h") ){
    print_usage() ;  
    return -1; 
  }

  TString datatype = argv[1]; 
  TString spec_data    = argv[2]; 
  TString cut      = argv[3]; 
  TString infile   = argv[4]; 
  TString outfile  = argv[5]; 

  Printf("datatype: '%s'", datatype.Data());
  Printf("On which data: '%s'", spec_data.Data());
  Printf("cut: '%s'", cut.Data());
  Printf("input file: '%s'", infile.Data());
  Printf("output file: '%s'", outfile.Data());

  ///TChain *ch = new TChain("ntuple/tree"); 
  TChain *ch = new TChain("tree"); 
  ch->Add(infile.Data()); 

  char *j = get_option(argv, argv+argc, "-j");
  if (j) {
    TProof::Open(Form("workers=%s", j));
    ch->SetProof(); 
  }

  Long64_t nentries = 1000000000; 
  char * n = get_option(argv, argv+argc, "-n");  
  if (n){
    nentries = atoi(n);
  }
    
  int     iStart = 0;
  char *s = get_option(argv, argv+argc, "-s");
  printf("Number of entries is %lld.\n",ch->GetEntries());
  if (s) {
    iStart = atoi(s);
    if (iStart > ch->GetEntries()){
      printf("ERROR: Number of entries is %lld.\n",ch->GetEntries());
      return -1;
    }
  }

  TString option; 
  option.Form("datatype=%s;spec_data=%s;cut=%s;ofile=files/sel_%s_%s_%s_s%d.root", datatype.Data(), spec_data.Data(), cut.Data(), outfile.Data(), datatype.Data(), cut.Data(), iStart); 
    
  // // It's not allowed to run with fat trees!
  // if (datatype.Data() == "mc" && (!(s) || !(n))){
  //   printf("WARNING: You must specify #entries(-n) and start run(-s) for datatype '%s'.\n",datatype.Data());
  //   return -1;
  // }
  std::cout<<"entries "<<nentries<<" start "<<iStart<<std::endl;
  ch->Process("SingleBsToPhiMuMuSelector.cc+", option, nentries, iStart); 

  gSystem->Exit(0);

  return 0 ;
}

#endif
