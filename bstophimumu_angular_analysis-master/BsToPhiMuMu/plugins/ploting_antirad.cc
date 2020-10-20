#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TEfficiency.h>
#include <TH1F.h>
#include <TH2.h>
#include <TCanvas.h>
#include <iostream>
#include <TString.h>
#include <TBox.h>

using namespace std;

const int nQ2Ranges = 9;

char genQ2range[nQ2Ranges][128] = {"genQ2 <4.30 && genQ2 > 1.00",//Bin 0
                                  "genQ2 < 8.68 && genQ2 > 4.30",//Bin 1
				  "genQ2 <10.09 && genQ2 > 8.68",//Jpsi (Bin 2)
                                  "genQ2 <12.86 && genQ2 >10.09",//Bin 3     
                                  "genQ2 <14.18 && genQ2 >12.86",//Psi2S (Bin 4)
                                  "genQ2 <16.00 && genQ2 >14.18",//Bin 5        
                                  "genQ2 <19.00 && genQ2 >16.00",//Bin 6       
                                  "genQ2 < 6.00 && genQ2 > 1.00",// Bin 7 Summary Bin 1
                                  "genQ2 <19.00 && genQ2 >1.00"};//Bin 8 Summary Bin 2 
double q2rangedn[nQ2Ranges] = {1.00 , 4.30 , 8.68 , 10.09 , 12.86 , 14.18 , 16.00 , 1.00 , 1.00};
double q2rangeup[nQ2Ranges] = {4.30 , 8.68 , 10.09 , 12.86 , 14.18 , 16.00 , 19.00 , 6.00 , 19.00};

Double_t        Mumumass;
Double_t        Mumumasserr;
Double_t        Phimass;
Double_t        Kmpt;
Double_t        Kppt;
Double_t        Kmtrkdcasigbs;
Double_t        Kptrkdcasigbs;
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
Double_t        Bdt;
Double_t        dimupt;
Double_t        dimueta;
Double_t        CosThetaL;
Double_t        CosThetaK;
Double_t        Phi;
Int_t           JpsiTriggers;
Int_t           PsiPTriggers;
Int_t           LMNTTriggers;

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

void angular_distr(TString filename= "Modified_sel_BsToPhiMuMu_reco.root", TString filename1= "Modified_sel_BsToPhiMuMu_gen.root"){

  TFile smalltree_f(filename.Data());
  //  TFile smalltree_f("Modified_sel_BsToPhiMuMu_OfficialMC_signal_2016_mc_lite_cut0.root","READ");
  //  TFile smalltree_f("Modified_sel_BsToPhiMuMu_OfficialMC_signal_2016_mc.lite_cutopt.root","READ");
  TTree* red_tree = (TTree*)smalltree_f.Get("tree");

  red_tree->SetBranchAddress("Mumumass", &Mumumass);
  red_tree->SetBranchAddress("Mumumasserr", &Mumumasserr);
  red_tree->SetBranchAddress("Phimass", &Phimass);
  red_tree->SetBranchAddress("Kmpt", &Kmpt);
  red_tree->SetBranchAddress("Kppt", &Kppt);
  red_tree->SetBranchAddress("Kmtrkdcasigbs", &Kmtrkdcasigbs);
  red_tree->SetBranchAddress("Kptrkdcasigbs", &Kptrkdcasigbs);
  red_tree->SetBranchAddress("Bmass", &Bmass);
  red_tree->SetBranchAddress("Bpt", &Bpt);
  red_tree->SetBranchAddress("Beta", &Beta);
  red_tree->SetBranchAddress("Bphi", &Bphi);
  red_tree->SetBranchAddress("Bvtxcl", &Bvtxcl);
  red_tree->SetBranchAddress("Blxysig", &Blxysig);
  red_tree->SetBranchAddress("Bcosalphabs", &Bcosalphabs);
  red_tree->SetBranchAddress("Bcosalphabs2d", &Bcosalphabs2d);
  red_tree->SetBranchAddress("Bctau", &Bctau);
  red_tree->SetBranchAddress("Q2", &Q2);
  red_tree->SetBranchAddress("dimupt", &dimupt);
  red_tree->SetBranchAddress("dimueta", &dimueta);
  red_tree->SetBranchAddress("CosThetaL", &CosThetaL);
  red_tree->SetBranchAddress("CosThetaK", &CosThetaK);
  red_tree->SetBranchAddress("Bdt", &Bdt);
  red_tree->SetBranchAddress("Phi", &Phi);
  red_tree->SetBranchAddress("JpsiTriggers", &JpsiTriggers);
  red_tree->SetBranchAddress("PsiPTriggers", &PsiPTriggers);
  red_tree->SetBranchAddress("LMNTTriggers", &LMNTTriggers);
  red_tree->SetBranchAddress("genBpid", &genBpid);
  red_tree->SetBranchAddress("genBPhi", &genBPhi);
  red_tree->SetBranchAddress("genMupPt", &genMupPt);
  red_tree->SetBranchAddress("genMupEta", &genMupEta);
  red_tree->SetBranchAddress("genMupPhi", &genMupPhi);
  red_tree->SetBranchAddress("genMumPt", &genMumPt);
  red_tree->SetBranchAddress("genMumEta", &genMumEta);
  red_tree->SetBranchAddress("genMumPhi", &genMumPhi);
  red_tree->SetBranchAddress("gendimuPt", &gendimuPt);
  red_tree->SetBranchAddress("gendimuEta", &gendimuEta);
  red_tree->SetBranchAddress("gendimuPhi", &gendimuPhi);
  red_tree->SetBranchAddress("genQ2", &genQ2);
  red_tree->SetBranchAddress("genCosThetaL", &genCosThetaL);
  red_tree->SetBranchAddress("genCosThetaK", &genCosThetaK);
  


  TH2D* mass_correl_MC=new TH2D("mass_correl_MC",";MuMumass;B mass ",50,2.8,5,50,4.8,5.9);
  TH2D* mass_correl_MC_rej=new TH2D("mass_correl_MC_rej","MC Rejected ;MuMumass;B mass ",50,2.8,5,50,4.8,5.9);

  TH2D* mass_correl_Data_No=new TH2D("mass_correl_Data_No","No cut;m_{#mu^{+}#mu^{-}}[GeV];B_{s} Mass[GeV] ",50,1.0,4.9,50,4.8,5.9);
  TH2D* mass_correl_Data=new TH2D("mass_correl_Data","+Anti Radiation;m_{#mu^{+}#mu^{-}}[GeV];B_{s} Mass[GeV] ",50,1.0,4.9,50,4.8,5.9);
  TH2D* mass_correl_Data_rej=new TH2D("mass_correl_Data_rej","J/#psi & #psi^{'} rejection ;m_{#mu^{+}#mu^{-}}[GeV]; B_{s} Mass[GeV] ",50,1.0,4.9,50,4.8,5.9);

  TH1D* mass_correl_1dData_No=new TH1D("mass_correl_1dData_No","No cut;B_{s} Mass[GeV] ",50,4.8,5.9);
  TH1D* mass_correl_1dData=new TH1D("mass_correl_1dData","+Anti Radiation;B_{s} Mass[GeV] ",50,4.8,5.9);
  TH1D* mass_correl_1dData_rej=new TH1D("mass_correl_1dData_rej","J/#psi & #psi^{'} rejection; B_{s} Mass[GeV] ",50,4.8,5.9);
  

  TH1F* hmass= new TH1F("hmass",";m_{#mu^{+}#mu^{-}}[GeV];Events",100,.5,5.0);
  int max=0 , max2=0, max3=0, nacc_red=0, nrec=0;
  int nocut=0;
  TCut No_cut="((JpsiTriggers > 0 || PsiPTriggers>0 || LMNTTriggers>0 )&& Bmass>4.8 && Bmass<5.9 &&  Phimass>1.01 && Phimass<1.03)";
  TCut rej_cut="(((Mumumass > 3.096916+3.5*Mumumasserr) || (Mumumass < 3.096916-5.5*Mumumasserr))&& ((Mumumass > 3.686109+3.5*Mumumasserr) || (Mumumass < 3.686109-3.5*Mumumasserr)))";
  TCut anti_cut="(fabs(Bmass-Mumumass-2.270)>0.180 && fabs(Bmass-Mumumass-1.680)>0.08)";

  TCut c0 = No_cut && rej_cut;
  TCut c1 = c0 && anti_cut;

  red_tree->Draw("Bmass:Mumumass>>mass_correl_Data_No", No_cut);
  red_tree->Draw("Bmass:Mumumass>>mass_correl_Data_rej", c0);
  red_tree->Draw("Bmass:Mumumass>>mass_correl_Data", c1);

  red_tree->Draw("Bmass>>mass_correl_1dData_No", No_cut);
  red_tree->Draw("Bmass>>mass_correl_1dData_rej", c0);
  red_tree->Draw("Bmass>>mass_correl_1dData", c1);
  
  TBox *box1 = new TBox(2.9461840,4.8, 3.1764760,5.9);
  box1->SetFillColor(kRed);
  box1->SetFillStyle(3357);
  //box1->SetLineColor(kRed);
  TBox *box2 = new TBox(3.5860842,4.8,3.7656341,5.9);
  box2->SetFillColor(kRed);
  box2->SetFillStyle(3357);
  //box2->SetLineColor(kRed);

  TBox *box3 = new TBox();//1.0,5.28,4.9,5.46);
  //TBox *box3 = new TBox(1.0,5.28,4.9,5.46);
  box3->SetFillColor(kRed);
  box3->SetFillStyle(3357);
  //box3->SetLineColor(kRed);

  mass_correl_Data_rej->SetMaximum(30000);
  mass_correl_Data->SetMaximum(30000);
  mass_correl_Data_No->SetMaximum(30000);
  TCanvas* canv_kin_rej1=new TCanvas("canv_kin_rej1","",1800,1200);
  canv_kin_rej1->Divide(3,2);
  canv_kin_rej1->cd(1); mass_correl_Data_No->Draw(); box1->Draw(); box2->Draw();box3->DrawBox(1.0,5.28,4.8,5.46);
  canv_kin_rej1->cd(2); mass_correl_Data_rej->Draw();box1->Draw(); box2->Draw();box3->DrawBox(1.0,5.28,4.8,5.46);
  canv_kin_rej1->cd(3); mass_correl_Data->Draw(); box1->Draw(); box2->Draw();box3->DrawBox(1.0,5.28,4.8,5.46);
  canv_kin_rej1->cd(4); mass_correl_1dData_No->Draw();
  canv_kin_rej1->cd(5); mass_correl_1dData_rej->Draw();
  canv_kin_rej1->cd(6); mass_correl_1dData->Draw();
  canv_kin_rej1->SaveAs("plot/histogram_antiradition2.pdf");

}

void ploting_antirad(){

  
  angular_distr("sel_Combine_2016_Mini_Presel_data_cut0.root","sel_BsToPhiMuMu_OfficialMC_signal_2016Mini_combine_Presel_mc.lite_cut0.root");

}
