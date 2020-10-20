//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 28 16:18:48 2017 by ROOT version 5.34/30
// from TTree tree/Bs->J/psi kaskey menos ntuple
// found on file: BsToPhiMuMu_2016B_10k.root
//////////////////////////////////////////////////////////

#ifndef SingleBsToPhiMuMuSelector_h
#define SingleBsToPhiMuMuSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TSelector.h>
#include <TTreePlayer.h>
// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"
#include "TVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"


#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "MVAnalysis.h"

using namespace std; 

// Fixed size dimensions of array or collections stored in the TTree if any.

class SingleBsToPhiMuMuSelector : public TSelector {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
  // Declaration of leaf types
  UInt_t          run;
  UInt_t          event;
  UInt_t          lumiblock;
  UInt_t          nprivtx;
  Double_t        fpuw8;
  vector<string>  *triggernames;
  vector<int>     *triggerprescales;
  vector<int>   *tri_JpsiTk;
  vector<int>   *tri_PsipTk;
  vector<int>   *tri_LMNTk;
  vector<double>   *tri_dr0;
  vector<double>   *tri_dr1;
  vector<double>   *tri_dpt0;
  vector<double>   *tri_dpt1;
  vector<double>  *mumdcabs;
  vector<double>  *mumdcabserr;
  vector<double>  *mumpx;
  vector<double>  *mumpy;
  vector<double>  *mumpz;
  vector<double>  *mupdcabs;
  vector<double>  *mupdcabserr;
  vector<double>  *muppx;
  vector<double>  *muppy;
  vector<double>  *muppz;
  vector<double>  *mumutrkr;
  vector<double>  *mumutrkz;
  vector<double>  *mumudca;
  vector<double>  *mumuvtxcl;
  vector<double>  *mumulsbs;
  vector<double>  *mumulsbserr;
  vector<double>  *mumucosalphabs;
  vector<double>  *mumucosalphabserr;
  vector<double>  *mumumass;
  vector<double>  *mumumasserr;
  
  vector<double>  *mumMinIP;
  vector<double>  *mupMinIP;
  vector<double>  *mumMinIPE;
  vector<double>  *mupMinIPE;
  vector<double>  *kptrkMinIP;
  vector<double>  *kptrkMinIPE;
  vector<double>  *kmtrkMinIP;
  vector<double>  *kmtrkMinIPE;
  vector<double>  *bdcabs;
  vector<double>  *bdcabserr;
  vector<double>  *mumMinIP2D;
  vector<double>  *mupMinIP2D;
  vector<double>  *mumMinIP2DE;
  vector<double>  *mupMinIP2DE;
  vector<double>  *kptrkMinIP2D;
  vector<double>  *kptrkMinIP2DE;
  vector<double>  *kmtrkMinIP2D;
  vector<double>  *kmtrkMinIP2DE;

  vector<vector<double> > * mupIsoPt ;
  vector<vector<double> > * mumIsoPt ;
  vector<vector<double> > * mupIsodR ;
  vector<vector<double> > * mumIsodR ;
  vector<vector<double> > * kptrkIsoPt ;
  vector<vector<double> > * kmtrkIsoPt ;
  vector<vector<double> > * kptrkIsodR ;
  vector<vector<double> > * kmtrkIsodR ;
  vector<vector<double> > * bmassIsoPt ;
  vector<vector<double> > * bmassIsodR ;

  vector<double >  *rawmuppt;
  vector<double >  *rawmumpt;
  vector<double >  *rawTrkppt;
  vector<double >  *rawTrkmpt;
  
  vector<bool>    *mumisgoodmuon;
  vector<bool>    *mupisgoodmuon;
  //vector<bool>    *mumloosemuon;
  // vector<bool>    *muploosemuon;
  vector<int>     *mumnpixhits;
  vector<int>     *mupnpixhits;
  vector<int>     *mumnpixlayers;
  vector<int>     *mupnpixlayers;
  vector<int>     *mumntrkhits;
  vector<int>     *mupntrkhits;
  vector<int>     *mumntrklayers;
  vector<int>     *mupntrklayers;
  vector<double>  *mumnormchi2;
  vector<double>  *mupnormchi2;
  vector<int>     *mumtrkqual;
  vector<int>     *muptrkqual;
  vector<double>  *mumdxyvtx;
  vector<double>  *mupdxyvtx;
  vector<double>  *mumdzvtx;
  vector<double>  *mupdzvtx;
  vector<string>  *mumtriglastfilter;
  vector<string>  *muptriglastfilter;
  vector<double>  *mumpt;
  vector<double>  *muppt;
  vector<double>  *mumeta;
  vector<double>  *mupeta;
  vector<double>  *kptrkdcabs;
  vector<double>  *kptrkdcabserr;
  vector<double>  *kmtrkdcabs;
  vector<double>  *kmtrkdcabserr;
  vector<int>     *kpchg;
  vector<int>     *kptrkqual;
  vector<int>     *kmtrkqual;
  
  vector<double>  *kppx;
  vector<double>  *kppy;
  vector<double>  *kppz;
  vector<int>     *kmchg;
  vector<double>  *kmpx;
  vector<double>  *kmpy;
  vector<double>  *kmpz;
  vector<double>  *phimass;
  Int_t           nb;
  vector<double>  *bpx;
  vector<double>  *bpxerr;
  vector<double>  *bpy;
  vector<double>  *bpyerr;
  vector<double>  *bpz;
  vector<double>  *bpzerr;
  vector<double>  *bmass;
  vector<double>  *bmasserr;
  vector<double>  *bvtxcl;
  vector<double>  *bvtxx;
  vector<double>  *bvtxxerr;
  vector<double>  *bvtxy;
  vector<double>  *bvtxyerr;
  vector<double>  *bvtxz;
  vector<double>  *bvtxzerr;
  vector<double>  *bcosalphabs;
  vector<double>  *bcosalphabserr;
  vector<double>  *bcosalphabs2d;
  vector<double>  *bcosalphabs2derr;
  vector<double>  *blsbs;
  vector<double>  *blsbserr;
  vector<double>  *bctau;
  vector<double>  *bctauerr;

  
  Double_t        genbpid;
  Double_t        genbpx;
  Double_t        genbpy;
  Double_t        genbpz;
  Double_t        genphipx;
  Double_t        genphipy;
  Double_t        genphipz;
  Double_t        genphivtxx;
  Double_t        genphivtxy;
  Double_t        genphivtxz;
  Int_t           genkpchg;
  Double_t        genkppx;
  Double_t        genkppy;
  Double_t        genkppz;
  Int_t           genkmchg;
  Double_t        genkmpx;
  Double_t        genkmpy;
  Double_t        genkmpz;
  Double_t        genmumpx;
  Double_t        genmumpy;
  Double_t        genmumpz;
  Double_t        genmuppx;
  Double_t        genmuppy;
  Double_t        genmuppz;
  string          *decname;
  vector<bool>    *istruemum;
  vector<bool>    *istruemup;
  vector<bool>    *istruekp;
  vector<bool>    *istruekm;
  vector<bool>    *istruebs;
  
   // List of branches
  TBranch        *b_run;   //!
  TBranch        *b_event;   //!
  TBranch        *b_lumiblock;   //!
  TBranch        *b_nprivtx;   //!
  TBranch        *b_triggernames;   //!
  TBranch        *b_triggerprescales;   //!

  TBranch        *b_tri_dr0;
  TBranch        *b_tri_dr1;
  TBranch        *b_tri_dpt0;
  TBranch        *b_tri_dpt1;
  
  TBranch        *b_tri_JpsiTk;
  TBranch        *b_tri_PsipTk;
  TBranch        *b_tri_LMNTk;
  
  TBranch        *b_mumdcabs;   //!
  TBranch        *b_mumdcabserr;   //!
  TBranch        *b_mumpx;   //!
  TBranch        *b_mumpy;   //!
  TBranch        *b_mumpz;   //!
  TBranch        *b_mupdcabs;   //!
  TBranch        *b_mupdcabserr;   //!
  TBranch        *b_muppx;   //!
  TBranch        *b_muppy;   //!
  TBranch        *b_muppz;   //!
  TBranch        *b_mumutrkr;   //!
  TBranch        *b_mumutrkz;   //!
  TBranch        *b_mumudca;   //!
  TBranch        *b_mumuvtxcl;   //!
  TBranch        *b_mumulsbs;   //!
  TBranch        *b_mumulsbserr;   //!
  TBranch        *b_mumucosalphabs;   //!
  TBranch        *b_mumucosalphabserr;   //!
  TBranch        *b_mumumass;   //!
  TBranch        *b_mumumasserr;   //!
  TBranch        *b_mumisgoodmuon;   //!
  TBranch        *b_mupisgoodmuon;   //!
  //TBranch        *b_mumloosemuon;   //!
  //TBranch        *b_muploosemuon;   //!

  TBranch        *b_mumMinIP;
  TBranch        *b_mupMinIP;
  TBranch        *b_mumMinIPE;
  TBranch        *b_mupMinIPE;
  TBranch        *b_kptrkMinIP;
  TBranch        *b_kptrkMinIPE;
  TBranch        *b_kmtrkMinIP;
  TBranch        *b_kmtrkMinIPE;
  TBranch        *b_bdcabs;
  TBranch        *b_bdcabserr;
  TBranch        *b_mumMinIP2D;
  TBranch        *b_mupMinIP2D;
  TBranch        *b_mumMinIP2DE;
  TBranch        *b_mupMinIP2DE;
  TBranch        *b_kptrkMinIP2D;
  TBranch        *b_kptrkMinIP2DE;
  TBranch        *b_kmtrkMinIP2D;
  TBranch        *b_kmtrkMinIP2DE;
  TBranch        * b_mupIsoPt ;
  TBranch        * b_mumIsoPt ;
  TBranch        * b_mupIsodR ;
  TBranch        * b_mumIsodR ;
  TBranch        * b_kptrkIsoPt ;
  TBranch        * b_kmtrkIsoPt ;
  TBranch        * b_kptrkIsodR ;
  TBranch        * b_kmtrkIsodR ;
  TBranch        * b_bmassIsodR ;
  TBranch        * b_bmassIsoPt ;

  TBranch        * b_rawmuppt;
  TBranch        * b_rawmumpt;
  TBranch        * b_rawTrkppt;
  TBranch        * b_rawTrkmpt;
  
  TBranch        *b_mumnpixhits;   //!
  TBranch        *b_mupnpixhits;   //!
  TBranch        *b_mumnpixlayers;   //!
  TBranch        *b_mupnpixlayers;   //!
  TBranch        *b_mumntrkhits;   //!
  TBranch        *b_mupntrkhits;   //!
  TBranch        *b_mumntrklayers;   //!
  TBranch        *b_mupntrklayers;   //!
  TBranch        *b_mumnormchi2;   //!
  TBranch        *b_mupnormchi2;   //!
  TBranch        *b_mumtrkqual;   //!
  TBranch        *b_muptrkqual;   //!
  TBranch        *b_mumdxyvtx;   //!
  TBranch        *b_mupdxyvtx;   //!
  TBranch        *b_mumdzvtx;   //!
  TBranch        *b_mupdzvtx;   //!
  TBranch        *b_mumtriglastfilter;   //!
  TBranch        *b_muptriglastfilter;   //!
  TBranch        *b_mumpt;   //!
  TBranch        *b_muppt;   //!
  TBranch        *b_mumeta;   //!
  TBranch        *b_mupeta;   //!
  TBranch        *b_kptrkdcabs;   //!
  TBranch        *b_kptrkdcabserr;   //!
  TBranch        *b_kmtrkdcabs;   //!
  TBranch        *b_kmtrkdcabserr;   //!
  TBranch        *b_kpchg;   //!
  TBranch        *b_kptrkqual;   //!
  TBranch        *b_kmtrkqual;   //!
  TBranch        *b_kppx;   //!
  TBranch        *b_kppy;   //!
  TBranch        *b_kppz;   //!
  TBranch        *b_kmchg;   //!
  TBranch        *b_kmpx;   //!
  TBranch        *b_kmpy;   //!
  TBranch        *b_kmpz;   //!
  TBranch        *b_phimass;   //!
  TBranch        *b_nb;   //!
  TBranch        *b_bpx;   //!
  TBranch        *b_bpxerr;   //!
  TBranch        *b_bpy;   //!
  TBranch        *b_bpyerr;   //!
  TBranch        *b_bpz;   //!
  TBranch        *b_bpzerr;   //!
  TBranch        *b_bmass;   //!
  TBranch        *b_bmasserr;   //!
  TBranch        *b_bvtxcl;   //!
  TBranch        *b_bvtxx;   //!
  TBranch        *b_bvtxxerr;   //!
  TBranch        *b_bvtxy;   //!
  TBranch        *b_bvtxyerr;   //!
  TBranch        *b_bvtxz;   //!
  TBranch        *b_bvtxzerr;   //!
  TBranch        *b_bcosalphabs;   //!
  TBranch        *b_bcosalphabserr;   //!
  TBranch        *b_bcosalphabs2d;   //!
  TBranch        *b_bcosalphabs2derr;   //!
  TBranch        *b_blsbs;   //!
  TBranch        *b_blsbserr;   //!
  TBranch        *b_bctau;   //!
  TBranch        *b_bctauerr;   //!
  TBranch        *b_fpuw8;   //!
  
  TBranch        *b_genbpid;   //!
  TBranch        *b_genbpx;   //!
  TBranch        *b_genbpy;   //!
  TBranch        *b_genbpz;   //!
  TBranch        *b_genphipx;   //!
  TBranch        *b_genphipy;   //!
  TBranch        *b_genphipz;   //!
  TBranch        *b_genphivtxx;   //!
  TBranch        *b_genphivtxy;   //!
  TBranch        *b_genphivtxz;   //!
  TBranch        *b_genkpchg;   //!
  TBranch        *b_genkppx;   //!
  TBranch        *b_genkppy;   //!
  TBranch        *b_genkppz;   //!
  TBranch        *b_genkmchg;   //!
  TBranch        *b_genkmpx;   //!
  TBranch        *b_genkmpy;   //!
  TBranch        *b_genkmpz;   //!
  TBranch        *b_genmumpx;   //!
  TBranch        *b_genmumpy;   //!
  TBranch        *b_genmumpz;   //!
  TBranch        *b_genmuppx;   //!
  TBranch        *b_genmuppy;   //!
  TBranch        *b_genmuppz;   //!
  TBranch        *b_decname;   //!
  TBranch        *b_istruemum;   //!
  TBranch        *b_istruemup;   //!
  TBranch        *b_istruekp;   //!
  TBranch        *b_istruekm;   //!
  TBranch        *b_istruebs; //!
  

   //singlebstophimumuselector(TTree * /*tree*/ =0) : fChain(0) { }
   SingleBsToPhiMuMuSelector(TTree * /*tree*/ =0) { }
   virtual ~SingleBsToPhiMuMuSelector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   //ClassDef(SingleBsToPhiMuMuSelector,0);
   int SelectB(string); 
   bool HasGoodDimuon(int); 
   bool HasGoodPreselection(int);
   void SaveEvent(int, string); 
   void SaveGen();
   //string get_option_value(string, string);
 private:

   bool createMVATree_ {false};
   bool readMVA_ {false};

   std::string mvaXMLFile_ {""};
   std::string mvAlgo_ {"BDT"};
   std::unique_ptr<MVAnalysis> mvAna_ {nullptr};


};

#endif

#ifdef SingleBsToPhiMuMuSelector_cxx
void SingleBsToPhiMuMuSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   triggernames = 0;
   triggerprescales = 0;
   tri_JpsiTk =0;
   tri_PsipTk =0;
   tri_LMNTk =0;
   tri_dr0 = 0;
   tri_dr1 = 0;
   tri_dpt0 = 0;
   tri_dpt1 = 0;


   mumdcabs = 0;
   mumdcabserr = 0;
   mumpx = 0;
   mumpy = 0;
   mumpz = 0;
   mupdcabs = 0;
   mupdcabserr = 0;
   muppx = 0;
   muppy = 0;
   muppz = 0;
   mumutrkr = 0;
   mumutrkz = 0;
   mumudca = 0;
   mumuvtxcl = 0;
   mumulsbs = 0;
   mumulsbserr = 0;
   mumucosalphabs = 0;
   mumucosalphabserr = 0;
   mumumass = 0;
   mumumasserr = 0;
   mumisgoodmuon = 0;
   mupisgoodmuon = 0;
   //mumloosemuon = 0;
   //muploosemuon = 0;
   mumMinIP = 0;
   mupMinIP = 0;
   mumMinIPE = 0;
   mupMinIPE = 0;
   kptrkMinIP = 0;
   kptrkMinIPE = 0;
   kmtrkMinIP = 0;
   kmtrkMinIPE = 0;
   bdcabs = 0;
   bdcabserr = 0;

   mumMinIP2D = 0;
   mupMinIP2D = 0;
   mumMinIP2DE = 0;
   mupMinIP2DE = 0;
   kptrkMinIP2D = 0;
   kptrkMinIP2DE = 0;
   kmtrkMinIP2D = 0;
   kmtrkMinIP2DE = 0;

   mupIsoPt = NULL;
   mumIsoPt = NULL;
   mupIsodR = NULL;
   mumIsodR = NULL;
   kptrkIsoPt = NULL;
   kmtrkIsoPt = NULL;
   kptrkIsodR = NULL;
   kmtrkIsodR = NULL;
   bmassIsoPt = NULL;
   bmassIsodR = NULL;

   rawmuppt = 0;
   rawmumpt = 0;
   rawTrkppt = 0;
   rawTrkmpt = 0;
   
   mumnpixhits = 0;
   mupnpixhits = 0;
   mumnpixlayers = 0;
   mupnpixlayers = 0;
   mumntrkhits = 0;
   mupntrkhits = 0;
   mumntrklayers = 0;
   mupntrklayers = 0;
   mumnormchi2 = 0;
   mupnormchi2 = 0;
   mumtrkqual = 0;
   muptrkqual = 0;
   mumdxyvtx = 0;
   mupdxyvtx = 0;
   mumdzvtx = 0;
   mupdzvtx = 0;
   mumtriglastfilter = 0;
   muptriglastfilter = 0;
   mumpt = 0;
   muppt = 0;
   mumeta = 0;
   mupeta = 0;
   kptrkdcabs = 0;
   kptrkdcabserr = 0;
   kmtrkdcabs = 0;
   kmtrkdcabserr = 0;
   kpchg = 0;
   kptrkqual = 0;
   kmtrkqual = 0;
   kppx = 0;
   kppy = 0;
   kppz = 0;
   kmchg = 0;
   kmpx = 0;
   kmpy = 0;
   kmpz = 0;
   phimass = 0;
   bpx = 0;
   bpxerr = 0;
   bpy = 0;
   bpyerr = 0;
   bpz = 0;
   bpzerr = 0;
   bmass = 0;
   bmasserr = 0;
   bvtxcl = 0;
   bvtxx = 0;
   bvtxxerr = 0;
   bvtxy = 0;
   bvtxyerr = 0;
   bvtxz = 0;
   bvtxzerr = 0;
   bcosalphabs = 0;
   bcosalphabserr = 0;
   bcosalphabs2d = 0;
   bcosalphabs2derr = 0;
   blsbs = 0;
   blsbserr = 0;
   bctau = 0;
   bctauerr = 0;


   decname = 0;
   istruemum = 0;
   istruemup = 0;
   istruekp = 0;
   istruekm = 0;
   istruebs = 0;
   
   // Set branch addr
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("nprivtx", &nprivtx, &b_nprivtx);
   fChain->SetBranchAddress("triggernames", &triggernames, &b_triggernames);
   fChain->SetBranchAddress("triggerprescales", &triggerprescales, &b_triggerprescales);
   fChain->SetBranchAddress("tri_JpsiTk",&tri_JpsiTk, &b_tri_JpsiTk);
   fChain->SetBranchAddress("tri_PsipTk",&tri_PsipTk, &b_tri_PsipTk);
   fChain->SetBranchAddress("tri_LMNTk",&tri_LMNTk, &b_tri_LMNTk);
   fChain->SetBranchAddress("tri_dr0",&tri_dr0, &b_tri_dr0);
   fChain->SetBranchAddress("tri_dr1",&tri_dr1, &b_tri_dr1);
   fChain->SetBranchAddress("tri_dpt0",&tri_dpt0, &b_tri_dpt0);
   fChain->SetBranchAddress("tri_dpt1",&tri_dpt1, &b_tri_dpt1);

   fChain->SetBranchAddress("mumdcabs", &mumdcabs, &b_mumdcabs);
   fChain->SetBranchAddress("mumdcabserr", &mumdcabserr, &b_mumdcabserr);
   fChain->SetBranchAddress("mumpx", &mumpx, &b_mumpx);
   fChain->SetBranchAddress("mumpy", &mumpy, &b_mumpy);
   fChain->SetBranchAddress("mumpz", &mumpz, &b_mumpz);
   fChain->SetBranchAddress("mupdcabs", &mupdcabs, &b_mupdcabs);
   fChain->SetBranchAddress("mupdcabserr", &mupdcabserr, &b_mupdcabserr);
   fChain->SetBranchAddress("muppx", &muppx, &b_muppx);
   fChain->SetBranchAddress("muppy", &muppy, &b_muppy);
   fChain->SetBranchAddress("muppz", &muppz, &b_muppz);
   fChain->SetBranchAddress("mumutrkr", &mumutrkr, &b_mumutrkr);
   fChain->SetBranchAddress("mumutrkz", &mumutrkz, &b_mumutrkz);
   fChain->SetBranchAddress("mumudca", &mumudca, &b_mumudca);
   fChain->SetBranchAddress("mumuvtxcl", &mumuvtxcl, &b_mumuvtxcl);
   fChain->SetBranchAddress("mumulsbs", &mumulsbs, &b_mumulsbs);
   fChain->SetBranchAddress("mumulsbserr", &mumulsbserr, &b_mumulsbserr);
   fChain->SetBranchAddress("mumucosalphabs", &mumucosalphabs, &b_mumucosalphabs);
   fChain->SetBranchAddress("mumucosalphabserr", &mumucosalphabserr, &b_mumucosalphabserr);
   fChain->SetBranchAddress("mumumass", &mumumass, &b_mumumass);
   fChain->SetBranchAddress("mumumasserr", &mumumasserr, &b_mumumasserr);
   fChain->SetBranchAddress("mumisgoodmuon", &mumisgoodmuon, &b_mumisgoodmuon);
   fChain->SetBranchAddress("mupisgoodmuon", &mupisgoodmuon, &b_mupisgoodmuon);
   //fChain->SetBranchAddress("mumloosemuon", &mumloosemuon, &b_mumloosemuon);
   //fChain->SetBranchAddress("muploosemuon", &muploosemuon, &b_muploosemuon);
   fChain->SetBranchAddress("mumMinIP",&mumMinIP,&b_mumMinIP);
   fChain->SetBranchAddress("mupMinIP",&mupMinIP, &b_mupMinIP);
   fChain->SetBranchAddress("mumMinIPE", &mumMinIPE, &b_mumMinIPE);
   fChain->SetBranchAddress("mupMinIPE", &mupMinIPE, &b_mupMinIPE);
   fChain->SetBranchAddress("kptrkMinIP", &kptrkMinIP, &b_kptrkMinIP);
   fChain->SetBranchAddress("kptrkMinIPE", &kptrkMinIPE, &b_kptrkMinIPE);
   fChain->SetBranchAddress("kmtrkMinIP", &kmtrkMinIP, &b_kmtrkMinIP);
   fChain->SetBranchAddress("kmtrkMinIPE", &kmtrkMinIPE, &b_kmtrkMinIPE);
   fChain->SetBranchAddress("bdcabs", &bdcabs, &b_bdcabs);
   fChain->SetBranchAddress("bdcabserr", &bdcabserr, &b_bdcabserr);
   fChain->SetBranchAddress("mumMinIP2D", &mumMinIP2D, &b_mumMinIP2D);
   fChain->SetBranchAddress("mupMinIP2D", &mupMinIP2D, &b_mupMinIP2D);
   fChain->SetBranchAddress("mumMinIP2DE", &mumMinIP2DE, &b_mumMinIP2DE);
   fChain->SetBranchAddress("mupMinIP2DE", &mupMinIP2DE, &b_mupMinIP2DE);
   fChain->SetBranchAddress("kptrkMinIP2D", &kptrkMinIP2D, &b_kptrkMinIP2D);
   fChain->SetBranchAddress("kptrkMinIP2DE", &kptrkMinIP2DE, &b_kptrkMinIP2DE);
   fChain->SetBranchAddress("kmtrkMinIP2D", &kmtrkMinIP2D, &b_kmtrkMinIP2D);
   fChain->SetBranchAddress("kmtrkMinIP2DE", &kmtrkMinIP2DE, &b_kmtrkMinIP2DE);

   fChain->SetBranchAddress("mupIsoPt",&mupIsoPt, &b_mupIsoPt);
   fChain->SetBranchAddress("mumIsoPt",&mumIsoPt, &b_mumIsoPt);
   fChain->SetBranchAddress("mupIsodR",&mupIsodR, &b_mupIsodR);
   fChain->SetBranchAddress("mumIsodR",&mumIsodR, &b_mumIsodR);

   fChain->SetBranchAddress("kptrkIsoPt",&kptrkIsoPt, &b_kptrkIsoPt);
   fChain->SetBranchAddress("kmtrkIsoPt",&kmtrkIsoPt, &b_kmtrkIsoPt);
   fChain->SetBranchAddress("kptrkIsodR",&kptrkIsodR, &b_kptrkIsodR);
   fChain->SetBranchAddress("kmtrkIsodR",&kmtrkIsodR, &b_kmtrkIsodR);
   fChain->SetBranchAddress("bmassIsoPt",&bmassIsoPt, &b_bmassIsoPt);
   fChain->SetBranchAddress("bmassIsodR",&bmassIsodR, &b_bmassIsodR);
   
   fChain->SetBranchAddress("rawmuppt",&rawmuppt,&b_rawmuppt);
   fChain->SetBranchAddress("rawmumpt",&rawmumpt,&b_rawmumpt);
   fChain->SetBranchAddress("rawTrkppt",&rawTrkppt,&b_rawTrkppt);
   fChain->SetBranchAddress("rawTrkmpt",&rawTrkmpt,&b_rawTrkmpt);

   fChain->SetBranchAddress("mumnpixhits", &mumnpixhits, &b_mumnpixhits);
   fChain->SetBranchAddress("mupnpixhits", &mupnpixhits, &b_mupnpixhits);
   fChain->SetBranchAddress("mumnpixlayers", &mumnpixlayers, &b_mumnpixlayers);
   fChain->SetBranchAddress("mupnpixlayers", &mupnpixlayers, &b_mupnpixlayers);
   fChain->SetBranchAddress("mumntrkhits", &mumntrkhits, &b_mumntrkhits);
   fChain->SetBranchAddress("mupntrkhits", &mupntrkhits, &b_mupntrkhits);
   fChain->SetBranchAddress("mumntrklayers", &mumntrklayers, &b_mumntrklayers);
   fChain->SetBranchAddress("mupntrklayers", &mupntrklayers, &b_mupntrklayers);
   fChain->SetBranchAddress("mumnormchi2", &mumnormchi2, &b_mumnormchi2);
   fChain->SetBranchAddress("mupnormchi2", &mupnormchi2, &b_mupnormchi2);
   fChain->SetBranchAddress("mumtrkqual", &mumtrkqual, &b_mumtrkqual);
   fChain->SetBranchAddress("muptrkqual", &muptrkqual, &b_muptrkqual);
   fChain->SetBranchAddress("mumdxyvtx", &mumdxyvtx, &b_mumdxyvtx);
   fChain->SetBranchAddress("mupdxyvtx", &mupdxyvtx, &b_mupdxyvtx);
   fChain->SetBranchAddress("mumdzvtx", &mumdzvtx, &b_mumdzvtx);
   fChain->SetBranchAddress("mupdzvtx", &mupdzvtx, &b_mupdzvtx);
   fChain->SetBranchAddress("mumtriglastfilter", &mumtriglastfilter, &b_mumtriglastfilter);
   fChain->SetBranchAddress("muptriglastfilter", &muptriglastfilter, &b_muptriglastfilter);
   fChain->SetBranchAddress("mumpt", &mumpt, &b_mumpt);
   fChain->SetBranchAddress("muppt", &muppt, &b_muppt);
   fChain->SetBranchAddress("mumeta", &mumeta, &b_mumeta);
   fChain->SetBranchAddress("mupeta", &mupeta, &b_mupeta);
   fChain->SetBranchAddress("kptrkdcabs", &kptrkdcabs, &b_kptrkdcabs);
   fChain->SetBranchAddress("kptrkdcabserr", &kptrkdcabserr, &b_kptrkdcabserr);
   fChain->SetBranchAddress("kmtrkdcabs", &kmtrkdcabs, &b_kmtrkdcabs);
   fChain->SetBranchAddress("kmtrkdcabserr", &kmtrkdcabserr, &b_kmtrkdcabserr);
   fChain->SetBranchAddress("kpchg", &kpchg, &b_kpchg);
   fChain->SetBranchAddress("kptrkqual", &kptrkqual, &b_kptrkqual);
   fChain->SetBranchAddress("kmtrkqual", &kmtrkqual, &b_kmtrkqual);
   fChain->SetBranchAddress("kppx", &kppx, &b_kppx);
   fChain->SetBranchAddress("kppy", &kppy, &b_kppy);
   fChain->SetBranchAddress("kppz", &kppz, &b_kppz);
   fChain->SetBranchAddress("kmchg", &kmchg, &b_kmchg);
   fChain->SetBranchAddress("kmpx", &kmpx, &b_kmpx);
   fChain->SetBranchAddress("kmpy", &kmpy, &b_kmpy);
   fChain->SetBranchAddress("kmpz", &kmpz, &b_kmpz);
   fChain->SetBranchAddress("phimass", &phimass, &b_phimass);
   fChain->SetBranchAddress("nb", &nb, &b_nb);
   fChain->SetBranchAddress("bpx", &bpx, &b_bpx);
   fChain->SetBranchAddress("bpxerr", &bpxerr, &b_bpxerr);
   fChain->SetBranchAddress("bpy", &bpy, &b_bpy);
   fChain->SetBranchAddress("bpyerr", &bpyerr, &b_bpyerr);
   fChain->SetBranchAddress("bpz", &bpz, &b_bpz);
   fChain->SetBranchAddress("bpzerr", &bpzerr, &b_bpzerr);
   fChain->SetBranchAddress("bmass", &bmass, &b_bmass);
   fChain->SetBranchAddress("bmasserr", &bmasserr, &b_bmasserr);
   fChain->SetBranchAddress("bvtxcl", &bvtxcl, &b_bvtxcl);
   fChain->SetBranchAddress("bvtxx", &bvtxx, &b_bvtxx);
   fChain->SetBranchAddress("bvtxxerr", &bvtxxerr, &b_bvtxxerr);
   fChain->SetBranchAddress("bvtxy", &bvtxy, &b_bvtxy);
   fChain->SetBranchAddress("bvtxyerr", &bvtxyerr, &b_bvtxyerr);
   fChain->SetBranchAddress("bvtxz", &bvtxz, &b_bvtxz);
   fChain->SetBranchAddress("bvtxzerr", &bvtxzerr, &b_bvtxzerr);
   fChain->SetBranchAddress("bcosalphabs", &bcosalphabs, &b_bcosalphabs);
   fChain->SetBranchAddress("bcosalphabserr", &bcosalphabserr, &b_bcosalphabserr);
   fChain->SetBranchAddress("bcosalphabs2d", &bcosalphabs2d, &b_bcosalphabs2d);
   fChain->SetBranchAddress("bcosalphabs2derr", &bcosalphabs2derr, &b_bcosalphabs2derr);
   fChain->SetBranchAddress("blsbs", &blsbs, &b_blsbs);
   fChain->SetBranchAddress("blsbserr", &blsbserr, &b_blsbserr);
   fChain->SetBranchAddress("bctau", &bctau, &b_bctau);
   fChain->SetBranchAddress("bctauerr", &bctauerr, &b_bctauerr);
   fChain->SetBranchAddress("fpuw8", &fpuw8, &b_fpuw8);

   fChain->SetBranchAddress("genbpid", &genbpid, &b_genbpid);
   fChain->SetBranchAddress("genbpx", &genbpx, &b_genbpx);
   fChain->SetBranchAddress("genbpy", &genbpy, &b_genbpy);
   fChain->SetBranchAddress("genbpz", &genbpz, &b_genbpz);
   fChain->SetBranchAddress("genphipx", &genphipx, &b_genphipx);
   fChain->SetBranchAddress("genphipy", &genphipy, &b_genphipy);
   fChain->SetBranchAddress("genphipz", &genphipz, &b_genphipz);
   fChain->SetBranchAddress("genphivtxx", &genphivtxx, &b_genphivtxx);
   fChain->SetBranchAddress("genphivtxy", &genphivtxy, &b_genphivtxy);
   fChain->SetBranchAddress("genphivtxz", &genphivtxz, &b_genphivtxz);
   fChain->SetBranchAddress("genkpchg", &genkpchg, &b_genkpchg);
   fChain->SetBranchAddress("genkppx", &genkppx, &b_genkppx);
   fChain->SetBranchAddress("genkppy", &genkppy, &b_genkppy);
   fChain->SetBranchAddress("genkppz", &genkppz, &b_genkppz);
   fChain->SetBranchAddress("genkmchg", &genkmchg, &b_genkmchg);
   fChain->SetBranchAddress("genkmpx", &genkmpx, &b_genkmpx);
   fChain->SetBranchAddress("genkmpy", &genkmpy, &b_genkmpy);
   fChain->SetBranchAddress("genkmpz", &genkmpz, &b_genkmpz);
   fChain->SetBranchAddress("genmumpx", &genmumpx, &b_genmumpx);
   fChain->SetBranchAddress("genmumpy", &genmumpy, &b_genmumpy);
   fChain->SetBranchAddress("genmumpz", &genmumpz, &b_genmumpz);
   fChain->SetBranchAddress("genmuppx", &genmuppx, &b_genmuppx);
   fChain->SetBranchAddress("genmuppy", &genmuppy, &b_genmuppy);
   fChain->SetBranchAddress("genmuppz", &genmuppz, &b_genmuppz);
   fChain->SetBranchAddress("decname", &decname, &b_decname);
   fChain->SetBranchAddress("istruemum", &istruemum, &b_istruemum);
   fChain->SetBranchAddress("istruemup", &istruemup, &b_istruemup);
   fChain->SetBranchAddress("istruekp", &istruekp, &b_istruekp);
   fChain->SetBranchAddress("istruekm", &istruekm, &b_istruekm);
   fChain->SetBranchAddress("istruebs", &istruebs, &b_istruebs);
}

Bool_t SingleBsToPhiMuMuSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef SingleBsToPhiMuMuSelector_cxx
