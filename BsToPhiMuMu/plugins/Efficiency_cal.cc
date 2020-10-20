#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>

using namespace std;
bool MC=true;
double dEff(int in, int iN) {
  double n = (double)in;
  double N = (double)iN;
  return sqrt(((n+1)*(N-n+1))/((N+3)*(N+2)*(N+2)));
}
double dEfferr(double in, double iN, double inErr, double iNErr) {
  return sqrt((pow(iN*inErr,2)+pow(in*iNErr,2)));
}

double rEfferr(double in, double iN, double inErr, double iNErr) {
  return (in/iN)* sqrt((pow((inErr/in),2)+pow((iNErr/iN),2)));
}


void Efficiency(){
  double BfJpsi = 0.05961; 
  double BfJpsiErr = 0.00033; 

  double Bfpsi = 0.0079;
  double BfpsiErr = 0.0006; 

  double ypsi = 6625.3;
  double yJpsi = 81823.5;

  double ypsiErr = sqrt(6625.3);
  double yJpsiErr = sqrt(81823.5);

  double AccEff[3] = {0.0459, 0.0583, 0.0564};//Jpsi, psi, phimm
  double AccEffErr[3] = {0.0000859, 0.0000883, 0.0000854};//Jpsi, psi, phimm

  double recEff[3] = {0.02234, 0.0231, 0.0178};
  double recEffErr[3] = {0.00002234, 0.000023, 0.000027};

  double totEff[3] = {AccEff[0]*recEff[0], AccEff[1]*recEff[1], AccEff[2]*recEff[2]};//Jpsi, psi, phimm
  double totEffErr[3] = {dEfferr(recEff[0], AccEff[0], recEffErr[0], AccEffErr[0]), dEfferr(recEff[1], AccEff[1], recEffErr[1], AccEffErr[1]), dEfferr(recEff[2], AccEff[2], recEffErr[2], AccEffErr[2])};

  cout<<"totEff : "<<totEff[0]<<" totEff1 "<<totEff[1]<<" totEff2 "<<totEff[2]<<endl;

  double RelBFpsi_Jpsi = (ypsi/yJpsi) * (BfJpsi/Bfpsi) *(totEff[0]/totEff[1]);
  cout<<"Relative Branching fraction Psi/Jpsi "<<RelBFpsi_Jpsi<<endl;

  double  yieldErr = rEfferr(ypsi, yJpsi, ypsiErr, yJpsiErr);
  double  BFErr = rEfferr(BfJpsi, Bfpsi, BfJpsiErr, BfpsiErr);
  double  EffErr = rEfferr(totEff[0], totEff[1], totEffErr[0], totEffErr[1]);

  double RelErr = RelBFpsi_Jpsi* sqrt(pow(yieldErr/(ypsi/yJpsi),2)+ pow(BFErr/(BfJpsi/Bfpsi),2)+ pow(EffErr/(totEff[0]/totEff[1]),2));
  cout<<"Relative Branching Fraction Error: "<<RelErr<<endl;
}

void convert_file(int isel){  
  Double_t        Mumumass;
  Double_t        Mumumasserr;
  Double_t        Phimass;
  Double_t        Kmpt;
  Double_t        Kppt;
  Double_t        Kmeta;
  Double_t        Kpeta;
  Double_t        Kmphi;
  Double_t        Kpphi;
  Double_t        Kmtrkdcasigbs;
  Double_t        Kptrkdcasigbs;
  Double_t        Mumpt;
  Double_t        Muppt;
  Double_t        Mumeta;
  Double_t        Mupeta;
  Double_t        Mumphi;
  Double_t        Mupphi;
  Double_t        Mumdcasigbs;
  Double_t        Mupdcasigbs;
  Int_t           Npv;

  Double_t        Bmass;
  Double_t        Bpt;
  Double_t        Beta;
  Double_t        Bphi;
  Double_t        Bvtxcl;
  Double_t        Blxysig;
  Double_t        Bcosalphabs;
  Double_t        Bcosalphabs2d;
  Double_t        Bctau;
  Double_t        Q2;
  Double_t        dimupt;
  Double_t        dimueta;
  Double_t        CosThetaL;
  Double_t        CosThetaK;
  Double_t        Phi;
  Double_t        Phipt;
  Double_t        Phieta;
  Double_t        Phiphi;
  Double_t        Bdt;
  Int_t           Triggers;

  Double_t        genBpid;
  Double_t        genBPhi;
  Double_t        genMupPt;
  Double_t        genMupEta;
  Double_t        genMupPhi;
  Double_t        genMumPt;
  Double_t        genMumEta;
  Double_t        genMumPhi;
  Double_t        gendimuPt;
  Double_t        gendimuEta;
  Double_t        gendimuPhi;
  Double_t        genQ2;
  Double_t        genCosThetaL;
  Double_t        genCosThetaK;

  Double_t        genBpid_g;
  Double_t        genBPhi_g;
  Double_t        genMupPt_g;
  Double_t        genMupEta_g;
  Double_t        genMupPhi_g;
  Double_t        genMumPt_g;
  Double_t        genMumEta_g;
  Double_t        genMumPhi_g;
  Double_t        gendimuPt_g;
  Double_t        gendimuEta_g;
  Double_t        gendimuPhi_g;
  Double_t        genQ2_g;
  Double_t        genCosThetaL_g;
  Double_t        genCosThetaK_g;
  Int_t           goodMuon;

  TChain *par_tree= new TChain("tree"); 
  TChain *gentree= new TChain("gentree"); 
  if(isel ==0){
    par_tree->Add("root://se01.indiacms.res.in//store/user/ckar/Analysis2016/sel_BsToJpsiPhi_2016MC_Official_Presel_mc.lite_cut_bdt.root/tree");
    gentree->Add("root://se01.indiacms.res.in//store/user/ckar/Analysis2016/sel_BsToJpsiPhi_2016MC_Official_Presel_mc.lite_cut_bdt.root/gentree");
  }else if(isel==1){
    par_tree->Add("~/t3store3/sel_BsToPsiPhi_2016MC_Official_Presel_mc.lite_cut_bdt.root/tree");
    gentree->Add("~/t3store3/sel_BsToPsiPhi_2016MC_Official_Presel_mc.lite_cut_bdt.root/gentree");
  }else if(isel==2){
    par_tree->Add("~/t3store3/sel_BsToPhiMuMu_OfficialMC_signal_2016_Presel_mc.lite_cut_bdt_genadded.root/tree");
    par_tree->Add("~/t3store3/sel_BsToPhiMuMu_OfficialMC_signal_2016_Jpsi_Psi_mc.lite_cut_bdt_addgentree.root/tree");
    
    gentree->Add("~/t3store3/sel_BsToPhiMuMu_OfficialMC_signal_2016_Presel_mc.lite_cut_bdt_genadded.root/gentree");
    gentree->Add("~/t3store3/sel_BsToPhiMuMu_OfficialMC_signal_2016_Jpsi_Psi_mc.lite_cut_bdt_addgentree.root/gentree");
  }
  
  par_tree->SetBranchAddress("Npv", &Npv);
  par_tree->SetBranchAddress("Mumumass", &Mumumass);
  par_tree->SetBranchAddress("Mumumasserr", &Mumumasserr);
  par_tree->SetBranchAddress("Phimass", &Phimass);
  par_tree->SetBranchAddress("Kmpt", &Kmpt);
  par_tree->SetBranchAddress("Kppt", &Kppt);
  par_tree->SetBranchAddress("Kmphi", &Kmphi);
  par_tree->SetBranchAddress("Kpphi", &Kpphi);
  par_tree->SetBranchAddress("Kmeta", &Kmeta);
  par_tree->SetBranchAddress("Kpeta", &Kpeta);
  par_tree->SetBranchAddress("Kmtrkdcasigbs", &Kmtrkdcasigbs);
  par_tree->SetBranchAddress("Kptrkdcasigbs", &Kptrkdcasigbs);
  par_tree->SetBranchAddress("Mumpt", &Mumpt);
  par_tree->SetBranchAddress("Muppt", &Muppt);
  par_tree->SetBranchAddress("Mumphi", &Mumphi);
  par_tree->SetBranchAddress("Mupphi", &Mupphi);
  par_tree->SetBranchAddress("Mumeta", &Mumeta);
  par_tree->SetBranchAddress("Mupeta", &Mupeta);
  par_tree->SetBranchAddress("Mumdcasigbs", &Mumdcasigbs);
  par_tree->SetBranchAddress("Mupdcasigbs", &Mupdcasigbs);

  par_tree->SetBranchAddress("Bmass", &Bmass);
  par_tree->SetBranchAddress("Bpt", &Bpt);
  par_tree->SetBranchAddress("Beta", &Beta);
  par_tree->SetBranchAddress("Bphi", &Bphi);
  par_tree->SetBranchAddress("Bvtxcl", &Bvtxcl);
  par_tree->SetBranchAddress("Blxysig", &Blxysig);
  par_tree->SetBranchAddress("Bcosalphabs", &Bcosalphabs);
  par_tree->SetBranchAddress("Bcosalphabs2d", &Bcosalphabs2d);
  par_tree->SetBranchAddress("Bctau", &Bctau);
  par_tree->SetBranchAddress("Q2", &Q2);
  par_tree->SetBranchAddress("dimupt", &dimupt);
  par_tree->SetBranchAddress("dimueta", &dimueta);
  par_tree->SetBranchAddress("CosThetaL", &CosThetaL);
  par_tree->SetBranchAddress("CosThetaK", &CosThetaK);
  par_tree->SetBranchAddress("Phi", &Phi);
  par_tree->SetBranchAddress("Phipt", &Phipt);
  par_tree->SetBranchAddress("Phieta", &Phieta);
  par_tree->SetBranchAddress("Phiphi", &Phiphi);
  par_tree->SetBranchAddress("Bdt", &Bdt);
  par_tree->SetBranchAddress("Triggers", &Triggers);
  if(isel==0 ||isel==1){
    par_tree->SetBranchAddress("goodMuon", &goodMuon);
  }
  
  par_tree->SetBranchAddress("genBpid", &genBpid);
  par_tree->SetBranchAddress("genBPhi", &genBPhi);
  par_tree->SetBranchAddress("genMupPt", &genMupPt);
  par_tree->SetBranchAddress("genMupEta", &genMupEta);
  par_tree->SetBranchAddress("genMupPhi", &genMupPhi);
  par_tree->SetBranchAddress("genMumPt", &genMumPt);
  par_tree->SetBranchAddress("genMumEta", &genMumEta);
  par_tree->SetBranchAddress("genMumPhi", &genMumPhi);
  par_tree->SetBranchAddress("gendimuPt", &gendimuPt);
  par_tree->SetBranchAddress("gendimuEta", &gendimuEta);
  par_tree->SetBranchAddress("gendimuPhi", &gendimuPhi);
  par_tree->SetBranchAddress("genQ2", &genQ2);
  par_tree->SetBranchAddress("genCosThetaL", &genCosThetaL);
  par_tree->SetBranchAddress("genCosThetaK", &genCosThetaK);

  // gentree->SetBranchAddress("genBpid", &genBpid_g);
  // gentree->SetBranchAddress("genBPhi", &genBPhi_g);
  // gentree->SetBranchAddress("genMupPt", &genMupPt_g);
  // gentree->SetBranchAddress("genMupEta", &genMupEta_g);
  // gentree->SetBranchAddress("genMupPhi", &genMupPhi_g);
  // gentree->SetBranchAddress("genMumPt", &genMumPt_g);
  // gentree->SetBranchAddress("genMumEta", &genMumEta_g);
  // gentree->SetBranchAddress("genMumPhi", &genMumPhi_g);
  // gentree->SetBranchAddress("gendimuPt", &gendimuPt_g);
  // gentree->SetBranchAddress("gendimuEta", &gendimuEta_g);
  // gentree->SetBranchAddress("gendimuPhi", &gendimuPhi_g);
  // gentree->SetBranchAddress("genQ2", &genQ2_g);
  // gentree->SetBranchAddress("genCosThetaL", &genCosThetaL_g);
  // gentree->SetBranchAddress("genCosThetaK", &genCosThetaK_g);
  

  
  int max=0 , max2=0;
  for(Long64_t lo_p = 0; lo_p <par_tree->GetEntries(); ++lo_p) {
    par_tree->GetEntry(lo_p);
    max++;
    if(lo_p%1000000==0) cout << "Processed " << lo_p << " out of " << par_tree->GetEntries() << " events" << endl;
    bool masscut =false;
    if(isel ==0){
      if((Mumumass < 3.096916+3.5*Mumumasserr) && (Mumumass > 3.096916-5.5*Mumumasserr) && goodMuon ==1) masscut=true;
    }else if(isel==1){
      if( (Mumumass < 3.6861+3.5*Mumumasserr) && (Mumumass > 3.6861-3.5*Mumumasserr) && goodMuon ==1) masscut=true;
    }else if(isel==2){
      if(((Mumumass > 3.096916+3.5*Mumumasserr) || (Mumumass < 3.096916-5.5*Mumumasserr)) && ((Mumumass > 3.686109+3.5*Mumumasserr) || (Mumumass < 3.686109-3.5*Mumumasserr))) masscut =true;//&& fabs(Bmass-Mumumass-2.270)>0.170 && fabs(Bmass-Mumumass-1.680)>0.080 
    }    
    if(Triggers>0 && Bdt>0.06 && masscut && Phimass>1.01 && Phimass<1.03 && Bmass>5.1 && Bmass<5.6){
      if(genMumPt>4.0 && genMupPt >4.0 && TMath::Abs(genMumEta)<2.2 && TMath::Abs(genMupEta)<2.2){
	max2++;
      }
    }
  }
  int gmax=0 , gmax2=0;
  // for(Long64_t lo_p = 0; lo_p <gentree->GetEntries(); ++lo_p) {
  //   gentree->GetEntry(lo_p);
  //   gmax++;
  //   if(lo_p%1000000==0) cout << "Processed " << lo_p << " out of " << gentree->GetEntries() << " events" << endl;
  //   if(genMumPt_g>4.0 && genMupPt_g >4.0 && TMath::Abs(genMumEta_g)<2.2 && TMath::Abs(genMupEta_g)<2.2){
  //     gmax2++;
  //   }
  // }
  gmax = gentree->GetEntries();
  gmax2 = gentree->GetEntries("genMumPt>2.5 && genMupPt >2.5 && TMath::Abs(genMumEta)<2.2 && TMath::Abs(genMupEta)<2.2"); 
  cout<<"max entry in tree loop "<<max<<"\t after cut "<<max2<<"\t eff "<<(double) max2/max<<endl;
  cout<<"max entry in gentree loop "<<gmax<<"\t after muoncut "<<gmax2<<"\t eff "<<(double) gmax2/gmax<<endl;

  double reff= (double) max2/max;
  double meff= (double) max/gmax2;
  double aeff;
  if(isel==2){
    aeff = (double) gmax2*2/gmax;
  }else
    aeff = (double) gmax2/gmax;
  double teff= reff*aeff*meff;

  cout<<"Total Efficiency: "<<(double)max2/gmax<<" check: reco "<<reff<<" acc "<<aeff<<" total "<<teff<<" Final "<<(double)max2/gmax2<<endl;
  


}
void Efficiency_cal(){
  cout<<"PsiPhi "<<endl;
  //convert_file(1);
  cout<<"Phimumu "<<endl;
  //convert_file(2);
  cout<<"JpsiPhi "<<endl;
  //convert_file(0);
  Efficiency();
}

