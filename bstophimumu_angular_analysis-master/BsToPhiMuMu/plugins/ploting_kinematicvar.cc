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
  

  //kinematic varriables
  TH1F* h_bdt[3];
  TH1F* h_kmpt[3];
  TH1F* h_kppt[3];
  TH1F* h_kminpt[3];
  TH1F* h_kmaxpt[3];
  TH1F* h_mtrkdcasig[3];
  TH1F* h_ptrkdcasig[3];
  TH1F* h_mintrkdcasig[3];
  TH1F* h_maxtrkdcasig[3];
  TH1F* h_Bvtxcl[3];
  TH1F* h_Blxysig[3];
  TH1F* h_Bcosalphabs2d[3];
  for(int i=0;i<3;i++){
    TString d="Data_SB"; if (i==1)d="Data_SG";if (i==2)d="MC";
    h_bdt[i]=new TH1F(Form("h_bdt_%s",d.Data()),"BDT ",50,0.,1.);
    h_kmpt[i]=new TH1F(Form("h_kmpt_%s",d.Data()),"K - pt;K^{-} p_{T};Events ",30,0.,15.);
    h_kppt[i]=new TH1F(Form("h_kppt_%s",d.Data()),"K + pt;K^{+} p_{T};Events ",30,0.,15.);
    h_kminpt[i]=new TH1F(Form("h_kminpt_%s",d.Data()),"KMin pt;K_{min} p_{T};Events ",30,0.,15.);
    h_kmaxpt[i]=new TH1F(Form("h_kmaxpt_%s",d.Data()),"KMax pt;K_{max} p_{T};Events ",30,0.,15.);
    h_mtrkdcasig[i]=new TH1F(Form("h_mtrkdcasig_%s",d.Data()),"mTrack dca significance; [-]trkdcasig",40,0.,40.);
    h_ptrkdcasig[i]=new TH1F(Form("h_ptrkdcasig_%s",d.Data()),"pTrack dca significance;[+]trkdccasig ",40,0.,40.);  
    h_mintrkdcasig[i]=new TH1F(Form("h_mintrkdcasig_%s",d.Data()),"minTrack dca significance;[+]trkdccasig ",40,0.,40.);  
    h_maxtrkdcasig[i]=new TH1F(Form("h_maxtrkdcasig_%s",d.Data()),"maxTrack dca significance;[+]trkdccasig ",40,0.,40.);  
    h_Bvtxcl[i]=new TH1F(Form("h_Bvtxcl_%s",d.Data()),"B vtx cl; vtxcl ", 50,0,0.9);
    h_Blxysig[i]=new TH1F(Form("h_Blxysig_%s",d.Data()),"B lx sig;L_{xy}sig ", 500,0,50);
    h_Bcosalphabs2d[i]=new TH1F(Form("h_Bcosalphabs2d_%s",d.Data()),"B cos alpha 2d;cos(#alpha_{2D}) ", 40,0.9,1.);
    h_bdt[i]->SetStats(0);
    h_kmpt[i]->SetStats(0);
    h_kppt[i]->SetStats(0);
    h_kminpt[i]->SetStats(0);
    h_kmaxpt[i]->SetStats(0);
    h_mtrkdcasig[i]->SetStats(0);
    h_ptrkdcasig[i]->SetStats(0);
    h_mintrkdcasig[i]->SetStats(0);
    h_maxtrkdcasig[i]->SetStats(0);
    h_Bvtxcl[i]->SetStats(0);
    h_Blxysig[i]->SetStats(0);
    h_Bcosalphabs2d[i]->SetStats(0);

  }

  TH2D* mass_correl_MC=new TH2D("mass_correl_MC",";MuMumass;B mass ",50,2.8,5,50,4.9,5.9);
  TH2D* mass_correl_MC_rej=new TH2D("mass_correl_MC_rej","MC Signal region ;MuMumass;B mass ",50,2.8,5,50,4.9,5.9);

  TH2D* mass_correl_Data_No=new TH2D("mass_correl_Data_No","Data;MuMumass;B mass ",50,2.8,4.2,50,4.9,5.9);
  TH2D* mass_correl_Data=new TH2D("mass_correl_Data","Data;MuMumass;B mass ",50,2.8,4.2,50,4.9,5.9);
  TH2D* mass_correl_Data_rej=new TH2D("mass_correl_Data_rej","Data Sideband sideband;MuMumass;B mass ",50,2.8,4.2,50,4.9,5.9);

  int max=0 , max2=0, max3=0, nacc_red=0, nrec=0;
  int nocut=0;
  for(Long64_t lo_p = 0; lo_p <red_tree->GetEntries(); ++lo_p) {
    red_tree->GetEntry(lo_p);
    max++;
    //    if (genQ2 > q2rangeup[2] && genQ2 < q2rangedn[2])continue;
    //if (genQ2 > q2rangeup[4] && genQ2 < q2rangedn[4])continue;
    if( (JpsiTriggers > 0 ||PsiPTriggers>0 ||LMNTTriggers>0) && Bmass>4.9 && Bmass<5.9 && Phimass>1.01 && Phimass<1.03){mass_correl_Data_No->Fill(Mumumass ,Bmass); nocut++;}
    if((JpsiTriggers > 0 ||PsiPTriggers>0 ||LMNTTriggers>0) && Bmass>4.9 && Bmass<5.9 && ((Mumumass > 3.096916+3.5*Mumumasserr) || (Mumumass < 3.096916-5.5*Mumumasserr)) && ((Mumumass > 3.686109+3.5*Mumumasserr) || (Mumumass < 3.686109-3.5*Mumumasserr))&& Phimass>1.01 && Phimass<1.03){      
      max2++;
      mass_correl_Data->Fill(Mumumass ,Bmass);
      if(fabs(Bmass-Mumumass-2.270)>0.180 && fabs(Bmass-Mumumass-1.680)>0.080){ 
	//mass_correl_Data_rej->Fill(Mumumass ,Bmass);
	max3++;
	if((Bmass>5.14 && Bmass<5.23) || (Bmass>5.51 && Bmass<5.6)){
	  mass_correl_Data_rej->Fill(Mumumass ,Bmass);
	  nacc_red++;
	  h_bdt[0]->Fill(Bdt);
	  h_kmpt[0]->Fill(Kmpt);
	  h_kppt[0]->Fill(Kppt);
	  h_kminpt[0]->Fill(TMath::Min(Kmpt, Kppt));
          h_kmaxpt[0]->Fill(TMath::Max(Kmpt, Kppt));
	  h_mtrkdcasig[0]->Fill(Kmtrkdcasigbs);
	  h_ptrkdcasig[0]->Fill(Kptrkdcasigbs);
	  h_mintrkdcasig[0]->Fill(TMath::Min(Kmtrkdcasigbs, Kptrkdcasigbs));
          h_maxtrkdcasig[0]->Fill(TMath::Max(Kmtrkdcasigbs, Kptrkdcasigbs));
	  h_Bvtxcl[0]->Fill(Bvtxcl);
	  h_Blxysig[0]->Fill(Blxysig);
	  h_Bcosalphabs2d[0]->Fill(Bcosalphabs2d);
	  //	  if(Mumumass>2.946 && Mumumass<3.176)cout<<"Mumumass "<<Mumumass<<"\terror"<<Mumumasserr<<endl;
	}
	else if( Bmass>5.28 && Bmass<5.46 ){
	  nrec++;

	  h_bdt[1]->Fill(Bdt);
	  h_kmpt[1]->Fill(Kmpt);
	  h_kppt[1]->Fill(Kppt);
	  h_kminpt[1]->Fill(TMath::Min(Kmpt, Kppt));
	  h_kmaxpt[1]->Fill(TMath::Max(Kmpt, Kppt));
	  h_mtrkdcasig[1]->Fill(Kmtrkdcasigbs);
	  h_ptrkdcasig[1]->Fill(Kptrkdcasigbs);
	  h_mintrkdcasig[1]->Fill(TMath::Min(Kmtrkdcasigbs, Kptrkdcasigbs));
	  h_maxtrkdcasig[1]->Fill(TMath::Max(Kmtrkdcasigbs, Kptrkdcasigbs));
	  h_Bvtxcl[1]->Fill(Bvtxcl);
	  h_Blxysig[1]->Fill(Blxysig);
	  h_Bcosalphabs2d[1]->Fill(Bcosalphabs2d);
	}
      }
    }
  }
    
  
  cout<<"max entry in loop "<<max<<"\t signal "<<nrec<<"\t SB "<<nacc_red<<endl;
  cout<<"No cut "<<nocut<<"\t after rejection "<<max2<<"\t after anti "<<max3<<endl;
  
  TCanvas* canv_kin=new TCanvas("canv_kin","",1200,800);
  canv_kin->Divide(3,3);
  canv_kin->cd(1); canv_kin->cd(1)->SetLogy(true); h_bdt[0]->Draw();h_bdt[1]->Draw("e0same");
  canv_kin->cd(2); canv_kin->cd(2)->SetLogy(false); h_kmpt[0]->Draw();
  canv_kin->cd(3); canv_kin->cd(3)->SetLogy(false); h_kppt[0]->Draw();
  canv_kin->cd(4); canv_kin->cd(4)->SetLogy(false); h_mtrkdcasig[0]->Draw();
  canv_kin->cd(5); canv_kin->cd(5)->SetLogy(false); h_ptrkdcasig[0]->Draw();
  canv_kin->cd(6); canv_kin->cd(6)->SetLogy(false); h_Bvtxcl[0]->Draw();
  canv_kin->cd(7); canv_kin->cd(7)->SetLogy(false); h_Blxysig[0]->Draw();
  canv_kin->cd(8); canv_kin->cd(8)->SetLogy(false); h_Bcosalphabs2d[0]->Draw();
  canv_kin->cd(9); canv_kin->cd(9)->SetLogy(false); mass_correl_Data_rej->Draw();
  canv_kin->SaveAs("plot/histogram_Kinematics_Data_rej_Presel.pdf");
  
  //cout<<" Integral signal region "<<h_bdtData[1]->Integral()<<"\t bkg "<<h_bdtData[0]->Integral()<<endl;
  //h_bdtData[1]->Add(h_bdtData[0],-1);
  //h_bdtData[1]->Scale(1/h_bdtData[1]->Integral());
  //h_bdtData[0]->Scale(1/h_bdtData[0]->Integral());
  //  h_bdtData->Add(h_bdt[0],-1);

  //return;
  
  TFile gentree_f(filename1.Data());
  TTree* gentree = (TTree*)gentree_f.Get("tree");
  Double_t        Mumumass_;
  Double_t        Mumumasserr_;
  Double_t        Phimass_;
  Double_t        Kmpt_;
  Double_t        Kppt_;
  Double_t        Kmtrkdcasigbs_;
  Double_t        Kptrkdcasigbs_;
  Double_t        Bmass_;
  Double_t        Bpt_;
  Double_t        Beta_;
  Double_t        Bphi_;
  Double_t        Bvtxcl_;
  Double_t        Blxysig_;
  Double_t        Bcosalphabs_;
  Double_t        Bcosalphabs2d_;
  Double_t        Bctau_;
  Double_t        Q2_;
  Double_t        dimupt_;
  Double_t        dimueta_;
  Double_t        CosThetaL_;
  Double_t        CosThetaK_;
  Double_t        Phi_;
  Double_t        Bdt_;

  Int_t           JpsiTriggers_;
  Int_t           PsiPTriggers_;
  Int_t           LMNTTriggers_;
  Double_t        genBpid_;
  Double_t        genBPhi_;
  Double_t        genMupPt_;
  Double_t        genMupEta_;
  Double_t        genMupPhi_;
  Double_t        genMumPt_;
  Double_t        genMumEta_;
  Double_t        genMumPhi_;
  Double_t        gendimuPt_;
  Double_t        gendimuEta_;
  Double_t        gendimuPhi_;
  Double_t        genQ2_;
  Double_t        genCosThetaL_;
  Double_t        genCosThetaK_;
  gentree->SetBranchAddress("Mumumass", &Mumumass_);
  gentree->SetBranchAddress("Mumumasserr", &Mumumasserr_);
  gentree->SetBranchAddress("Phimass", &Phimass_);
  gentree->SetBranchAddress("Kmpt", &Kmpt_);
  gentree->SetBranchAddress("Kppt", &Kppt_);
  gentree->SetBranchAddress("Kmtrkdcasigbs", &Kmtrkdcasigbs_);
  gentree->SetBranchAddress("Kptrkdcasigbs", &Kptrkdcasigbs_);
  gentree->SetBranchAddress("Bmass", &Bmass_);
  gentree->SetBranchAddress("Bpt", &Bpt_);
  gentree->SetBranchAddress("Beta", &Beta_);
  gentree->SetBranchAddress("Bphi", &Bphi_);
  gentree->SetBranchAddress("Bvtxcl", &Bvtxcl_);
  gentree->SetBranchAddress("Blxysig", &Blxysig_);
  gentree->SetBranchAddress("Bcosalphabs", &Bcosalphabs_);
  gentree->SetBranchAddress("Bcosalphabs2d", &Bcosalphabs2d_);
  gentree->SetBranchAddress("Bctau", &Bctau_);
  gentree->SetBranchAddress("Q2", &Q2_);
  gentree->SetBranchAddress("dimupt", &dimupt_);
  gentree->SetBranchAddress("dimueta", &dimueta_);
  gentree->SetBranchAddress("CosThetaL", &CosThetaL_);
  gentree->SetBranchAddress("CosThetaK", &CosThetaK_);
  gentree->SetBranchAddress("Phi", &Phi_);
  gentree->SetBranchAddress("JpsiTriggers", &JpsiTriggers_);
  gentree->SetBranchAddress("PsiPTriggers", &PsiPTriggers_);
  gentree->SetBranchAddress("LMNTTriggers", &LMNTTriggers_);
  
  gentree->SetBranchAddress("Bdt", &Bdt_);
  gentree->SetBranchAddress("genBpid", &genBpid_);
  gentree->SetBranchAddress("genBPhi", &genBPhi_);
  gentree->SetBranchAddress("genMupPt", &genMupPt_);
  gentree->SetBranchAddress("genMupPt", &genMupPt_);
  gentree->SetBranchAddress("genMupEta", &genMupEta_);
  gentree->SetBranchAddress("genMupPhi", &genMupPhi_);
  gentree->SetBranchAddress("genMumPt", &genMumPt_);
  gentree->SetBranchAddress("genMumEta", &genMumEta_);
  gentree->SetBranchAddress("genMumPhi", &genMumPhi_);
  gentree->SetBranchAddress("gendimuPt", &gendimuPt_);
  gentree->SetBranchAddress("gendimuEta", &gendimuEta_);
  gentree->SetBranchAddress("gendimuPhi", &gendimuPhi_);
  gentree->SetBranchAddress("genQ2", &genQ2_);
  gentree->SetBranchAddress("genCosThetaL", &genCosThetaL_);
  gentree->SetBranchAddress("genCosThetaK", &genCosThetaK_);



  int max_gen=0 , nred=0;
  for(Long64_t lo_p = 0; lo_p <gentree->GetEntries(); lo_p++) {
    gentree->GetEntry(lo_p);
    max_gen++;
    //if (genQ2_ > q2rangeup[2] && genQ2_ < q2rangedn[2])continue;
    //if (genQ2_ > q2rangeup[4] && genQ2_ < q2rangedn[4])continue;
    //if( fabs(genMumEta_)<2.5 && fabs(genMupEta_)< 2.5 && genMumPt_ > 2.5 && genMupPt_ >2.5){          
    //if( Triggers_>0 && Bmass_>4.7 &&((Mumumass_ > 3.096916+3.5*Mumumasserr_|| Mumumass_ < 3.096916-5.5*Mumumasserr_) && (Mumumass_ > 3.686109+3.5*Mumumasserr_ || Mumumass_ < 3.686109-3.5*Mumumasserr_))){					     
    if( (JpsiTriggers_ > 0 || LMNTTriggers_>0 ||PsiPTriggers_>0)&& Bmass_ !=0 && ((Mumumass_ > 3.096916+3.5*Mumumasserr_) || (Mumumass_ < 3.096916-5.5*Mumumasserr_)) && ((Mumumass_ > 3.686109+3.5*Mumumasserr_) || (Mumumass_ < 3.686109-3.5*Mumumasserr_)) && Phimass_>1.01 && Phimass_<1.03){      
      if(Bmass_>5.28 && Bmass_<5.46){	
	  mass_correl_MC->Fill(Mumumass_ ,Bmass_);
	  if(fabs(Bmass_-Mumumass_-2.270)>0.180 && fabs(Bmass_-Mumumass_-1.680)>0.080){mass_correl_MC_rej->Fill(Mumumass_ ,Bmass_);
	    nred++;
	    h_bdt[2]->Fill(Bdt_);
	    h_kmpt[2]->Fill(Kmpt_);
	    h_kppt[2]->Fill(Kppt_);
	    h_kminpt[2]->Fill(TMath::Min(Kmpt_, Kppt_));
	    h_kmaxpt[2]->Fill(TMath::Max(Kmpt_, Kppt_));
   
	    h_mtrkdcasig[2]->Fill(Kmtrkdcasigbs_);
	    h_ptrkdcasig[2]->Fill(Kptrkdcasigbs_);
	    h_mintrkdcasig[2]->Fill(TMath::Min(Kmtrkdcasigbs_, Kptrkdcasigbs_));
	    h_maxtrkdcasig[2]->Fill(TMath::Max(Kmtrkdcasigbs_, Kptrkdcasigbs_));
	    h_Bvtxcl[2]->Fill(Bvtxcl_);
	    h_Blxysig[2]->Fill(Blxysig_);
	    h_Bcosalphabs2d[2]->Fill(Bcosalphabs2d_);
	  }
	}
      }
    }
    

  //  h_bdtData->Scale(h_bdt[1]->Integral()/h_bdtData->Integral());
  cout<<"max entry in loop "<<max_gen<<"\t after cut "<<nred<<endl;
  TCanvas* canv_kin_MC=new TCanvas("canv_kin_MC","",1800,1800);
  canv_kin_MC->Divide(3,3);
  canv_kin_MC->cd(1); canv_kin_MC->cd(1)->SetLogy(false); h_bdt[2]->SetLineColor(kBlue); h_bdt[2]->Draw(); h_bdt[1]->Draw("e0same");
  canv_kin_MC->cd(2); canv_kin_MC->cd(2)->SetLogy(false); h_kmpt[2]->SetLineColor(kBlue);h_kmpt[2]->Draw();
  canv_kin_MC->cd(3); canv_kin_MC->cd(3)->SetLogy(false); h_kppt[2]->SetLineColor(kBlue); h_kppt[2]->Draw();
  canv_kin_MC->cd(4); canv_kin_MC->cd(4)->SetLogy(false); h_mtrkdcasig[2]->SetLineColor(kBlue); h_mtrkdcasig[2]->Draw();
  canv_kin_MC->cd(5); canv_kin_MC->cd(5)->SetLogy(false); h_ptrkdcasig[2]->SetLineColor(kBlue); h_ptrkdcasig[2]->Draw();
  canv_kin_MC->cd(6); canv_kin_MC->cd(6)->SetLogy(false); h_Bvtxcl[2]->SetLineColor(kBlue); h_Bvtxcl[2]->Draw();
  canv_kin_MC->cd(7); canv_kin_MC->cd(7)->SetLogy(false); h_Blxysig[2]->SetLineColor(kBlue); h_Blxysig[2]->Draw();
  canv_kin_MC->cd(8); canv_kin_MC->cd(8)->SetLogy(false); h_Bcosalphabs2d[2]->SetLineColor(kBlue); h_Bcosalphabs2d[2]->Draw();
  canv_kin_MC->cd(9); canv_kin_MC->cd(9)->SetLogy(false); mass_correl_MC->Draw();
  canv_kin_MC->SaveAs("plot/histogram_Kinematics_MC_norej_Presel.pdf");

  //h_bdt[0]->Scale(1/h_bdt[0]->Integral());
  h_kmpt[0]->Scale(1/h_kmpt[0]->Integral());
  h_kppt[0]->Scale(1/h_kppt[0]->Integral());
  h_kminpt[0]->Scale(1/h_kminpt[0]->Integral());
  h_kmaxpt[0]->Scale(1/h_kmaxpt[0]->Integral());
  h_mtrkdcasig[0]->Scale(1/ h_mtrkdcasig[0]->Integral());
  h_ptrkdcasig[0]->Scale(1/ h_ptrkdcasig[0]->Integral());
  h_mintrkdcasig[0]->Scale(1/ h_mintrkdcasig[0]->Integral());
  h_maxtrkdcasig[0]->Scale(1/ h_maxtrkdcasig[0]->Integral());
  h_Bvtxcl[0]->Scale(1/h_Bvtxcl[0]->Integral());
  h_Blxysig[0]->Scale(1/h_Blxysig[0]->Integral());
  h_Bcosalphabs2d[0]->Scale(1/h_Bcosalphabs2d[0]->Integral());


  //  h_bdt[1]->Scale(1/h_bdt[1]->Integral());
  h_kmpt[1]->Scale(1/h_kmpt[1]->Integral());
  h_kppt[1]->Scale(1/h_kppt[1]->Integral());
  h_kminpt[1]->Scale(1/h_kminpt[1]->Integral());
  h_kmaxpt[1]->Scale(1/h_kmaxpt[1]->Integral());
  h_mtrkdcasig[1]->Scale(1/ h_mtrkdcasig[1]->Integral());
  h_ptrkdcasig[1]->Scale(1/ h_ptrkdcasig[1]->Integral());
  h_mintrkdcasig[1]->Scale(1/ h_mintrkdcasig[1]->Integral());
  h_maxtrkdcasig[1]->Scale(1/ h_maxtrkdcasig[1]->Integral());
  h_Bvtxcl[1]->Scale(1/h_Bvtxcl[1]->Integral());
  h_Blxysig[1]->Scale(1/h_Blxysig[1]->Integral());
  h_Bcosalphabs2d[1]->Scale(1/h_Bcosalphabs2d[1]->Integral());


  //h_bdt[2]->Scale(1/h_bdt[2]->Integral());
  h_kmpt[2]->Scale(1/h_kmpt[2]->Integral());
  h_kppt[2]->Scale(1/h_kppt[2]->Integral());
  h_kminpt[2]->Scale(1/h_kminpt[2]->Integral());
  h_kmaxpt[2]->Scale(1/h_kmaxpt[2]->Integral());
  h_mtrkdcasig[2]->Scale(1/ h_mtrkdcasig[2]->Integral());
  h_ptrkdcasig[2]->Scale(1/ h_ptrkdcasig[2]->Integral());
  h_mintrkdcasig[2]->Scale(1/ h_mintrkdcasig[2]->Integral());
  h_maxtrkdcasig[2]->Scale(1/ h_maxtrkdcasig[2]->Integral());
  h_Bvtxcl[2]->Scale(1/h_Bvtxcl[2]->Integral());
  h_Blxysig[2]->Scale(1/h_Blxysig[2]->Integral());
  h_Bcosalphabs2d[2]->Scale(1/h_Bcosalphabs2d[2]->Integral());



  h_bdt[0]->SetFillStyle(3365);h_bdt[0]->SetFillColor(kRed);h_bdt[0]->SetLineColor(kRed);
  h_kmpt[0]->SetFillStyle(3365);h_kmpt[0]->SetFillColor(kRed);h_kmpt[0]->SetLineColor(kRed);
  h_kppt[0]->SetFillStyle(3365);  h_kppt[0]->SetFillColor(kRed);h_kppt[0]->SetLineColor(kRed);
  h_kminpt[0]->SetFillStyle(3365);  h_kminpt[0]->SetFillColor(kRed);h_kminpt[0]->SetLineColor(kRed);
  h_kmaxpt[0]->SetFillStyle(3365);  h_kmaxpt[0]->SetFillColor(kRed);h_kmaxpt[0]->SetLineColor(kRed);
  h_mtrkdcasig[0]->SetFillStyle(3365);  h_mtrkdcasig[0]->SetFillColor(kRed);h_mtrkdcasig[0]->SetLineColor(kRed);
  h_ptrkdcasig[0]->SetFillStyle(3365);  h_ptrkdcasig[0]->SetFillColor(kRed);h_ptrkdcasig[0]->SetLineColor(kRed);
  h_mintrkdcasig[0]->SetFillStyle(3365);  h_mintrkdcasig[0]->SetFillColor(kRed);h_mintrkdcasig[0]->SetLineColor(kRed);
  h_maxtrkdcasig[0]->SetFillStyle(3365);  h_maxtrkdcasig[0]->SetFillColor(kRed);h_maxtrkdcasig[0]->SetLineColor(kRed);
  h_Bvtxcl[0]->SetFillStyle(3365);  h_Bvtxcl[0]->SetFillColor(kRed);h_Bvtxcl[0]->SetLineColor(kRed);
  h_Blxysig[0]->SetFillStyle(3365);  h_Blxysig[0]->SetFillColor(kRed);h_Blxysig[0]->SetLineColor(kRed);
  h_Bcosalphabs2d[0]->SetFillStyle(3365);  h_Bcosalphabs2d[0]->SetFillColor(kRed);h_Bcosalphabs2d[0]->SetLineColor(kRed);


  h_bdt[1]->SetFillStyle(3354);h_bdt[1]->SetFillColor(kBlue);
  h_kmpt[1]->SetFillStyle(3354);h_kmpt[1]->SetFillColor(kBlue);
  h_kppt[1]->SetFillStyle(3354);  h_kppt[1]->SetFillColor(kBlue);
  h_kminpt[1]->SetFillStyle(3354);  h_kminpt[1]->SetFillColor(kBlue);
  h_kmaxpt[1]->SetFillStyle(3354);  h_kmaxpt[1]->SetFillColor(kBlue);
  h_mtrkdcasig[1]->SetFillStyle(3354);  h_mtrkdcasig[1]->SetFillColor(kBlue);
  h_ptrkdcasig[1]->SetFillStyle(3354);  h_ptrkdcasig[1]->SetFillColor(kBlue);
  h_mintrkdcasig[1]->SetFillStyle(3354);  h_mintrkdcasig[1]->SetFillColor(kBlue);
  h_maxtrkdcasig[1]->SetFillStyle(3354);  h_maxtrkdcasig[1]->SetFillColor(kBlue);
  h_Bvtxcl[1]->SetFillStyle(3354);  h_Bvtxcl[1]->SetFillColor(kBlue);
  h_Blxysig[1]->SetFillStyle(3354);  h_Blxysig[1]->SetFillColor(kBlue);
  h_Bcosalphabs2d[1]->SetFillStyle(3354);  h_Bcosalphabs2d[1]->SetFillColor(kBlue);


  h_bdt[2]->SetFillStyle(3354);h_bdt[2]->SetFillColor(kBlue);
  h_kmpt[2]->SetFillStyle(3354);h_kmpt[2]->SetFillColor(kBlue);
  h_kppt[2]->SetFillStyle(3354);  h_kppt[2]->SetFillColor(kBlue);
  h_kminpt[2]->SetFillStyle(3354);  h_kminpt[2]->SetFillColor(kBlue);
  h_kmaxpt[2]->SetFillStyle(3354);  h_kmaxpt[2]->SetFillColor(kBlue);
  h_mtrkdcasig[2]->SetFillStyle(3354);  h_mtrkdcasig[2]->SetFillColor(kBlue);
  h_ptrkdcasig[2]->SetFillStyle(3354);  h_ptrkdcasig[2]->SetFillColor(kBlue);
  h_mintrkdcasig[2]->SetFillStyle(3354);  h_mintrkdcasig[2]->SetFillColor(kBlue);
  h_maxtrkdcasig[2]->SetFillStyle(3354);  h_maxtrkdcasig[2]->SetFillColor(kBlue);
  h_Bvtxcl[2]->SetFillStyle(3354);  h_Bvtxcl[2]->SetFillColor(kBlue);
  h_Blxysig[2]->SetFillStyle(3354);  h_Blxysig[2]->SetFillColor(kBlue);
  h_Bcosalphabs2d[2]->SetFillStyle(3354);  h_Bcosalphabs2d[2]->SetFillColor(kBlue);



  TLegend* leg3 = new TLegend(0.7,0.7,0.89,0.89);
  leg3->AddEntry(h_kmpt[2],"MC","f");
  leg3->AddEntry(h_kmpt[0],"Sideband data","f"); 

  TCanvas* canv_kin_Over=new TCanvas("canv_kin_Over","",2500,2500);
  canv_kin_Over->Divide(4,4);
  canv_kin_Over->cd(1); canv_kin_Over->cd(1)->SetLogy(false); mass_correl_Data_rej->Draw();
  canv_kin_Over->cd(2); canv_kin_Over->cd(2)->SetLogy(false); mass_correl_MC_rej->Draw();
  canv_kin_Over->cd(3); canv_kin_Over->cd(3)->SetLogy(true); h_kmpt[2]->Draw("hist");h_kmpt[0]->Draw("samehist");  leg3->Draw();
  canv_kin_Over->cd(4); canv_kin_Over->cd(4)->SetLogy(true); h_kppt[2]->Draw("hist"); h_kppt[0]->Draw("samehist"); leg3->Draw();
  canv_kin_Over->cd(5); canv_kin_Over->cd(5)->SetLogy(true); h_mtrkdcasig[2]->Draw("hist"); h_mtrkdcasig[0]->Draw("samehist");leg3->Draw();
  canv_kin_Over->cd(6); canv_kin_Over->cd(6)->SetLogy(true); h_ptrkdcasig[2]->Draw("hist"); h_ptrkdcasig[0]->Draw("samehist");leg3->Draw();
  canv_kin_Over->cd(7); canv_kin_Over->cd(7)->SetLogy(true); h_Bvtxcl[2]->Draw("hist"); h_Bvtxcl[0]->Draw("samehist");leg3->Draw();
  canv_kin_Over->cd(8); canv_kin_Over->cd(8)->SetLogy(true); h_Blxysig[2]->Draw("hist"); h_Blxysig[0]->Draw("samehist");leg3->Draw();
  canv_kin_Over->cd(9); canv_kin_Over->cd(9)->SetLogy(true); h_Bcosalphabs2d[2]->Draw("hist"); h_Bcosalphabs2d[0]->Draw("samehist");leg3->Draw();
  canv_kin_Over->cd(10); canv_kin_Over->cd(10)->SetLogy(true); h_mintrkdcasig[2]->Draw("hist"); h_mintrkdcasig[0]->Draw("samehist");leg3->Draw();
  canv_kin_Over->cd(11); canv_kin_Over->cd(11)->SetLogy(true); h_maxtrkdcasig[2]->Draw("hist"); h_maxtrkdcasig[0]->Draw("samehist");leg3->Draw();
  canv_kin_Over->cd(12); canv_kin_Over->cd(12)->SetLogy(true); h_kminpt[2]->Draw("hist");h_kminpt[0]->Draw("samehist");  leg3->Draw();
  canv_kin_Over->cd(13); canv_kin_Over->cd(13)->SetLogy(true); h_kmaxpt[2]->Draw("hist");h_kmaxpt[0]->Draw("samehist");  leg3->Draw();
  //   canv_kin_Over->cd(3); canv_kin_Over->cd(3)->SetLogy(true); h_kmpt[0]->Draw("hist");h_kmpt[2]->Draw("samehist");  leg3->Draw();
  // canv_kin_Over->cd(4); canv_kin_Over->cd(4)->SetLogy(true); h_kppt[0]->Draw("hist"); h_kppt[2]->Draw("samehist"); leg3->Draw();
  // canv_kin_Over->cd(5); canv_kin_Over->cd(5)->SetLogy(true); h_mtrkdcasig[0]->Draw("hist"); h_mtrkdcasig[2]->Draw("samehist");leg3->Draw();
  // canv_kin_Over->cd(6); canv_kin_Over->cd(6)->SetLogy(true); h_ptrkdcasig[0]->Draw("hist"); h_ptrkdcasig[2]->Draw("samehist");leg3->Draw();
  // canv_kin_Over->cd(7); canv_kin_Over->cd(7)->SetLogy(true); h_Bvtxcl[0]->Draw("hist"); h_Bvtxcl[2]->Draw("samehist");leg3->Draw();
  // canv_kin_Over->cd(8); canv_kin_Over->cd(8)->SetLogy(true); h_Blxysig[0]->Draw("hist"); h_Blxysig[2]->Draw("samehist");leg3->Draw();
  // canv_kin_Over->cd(9); canv_kin_Over->cd(9)->SetLogy(true); h_Bcosalphabs2d[0]->Draw("hist"); h_Bcosalphabs2d[2]->Draw("samehist");leg3->Draw();
  // //canv_kin_Over->cd(9); canv_kin_Over->cd(9)->SetLogy(true); h_Bcosalphabs[0]->Draw("hist"); h_Bcosalphabs[1]->Draw("samehist");
  // canv_kin_Over->cd(10); canv_kin_Over->cd(10)->SetLogy(true); h_mintrkdcasig[0]->Draw("hist"); h_mintrkdcasig[2]->Draw("samehist");leg3->Draw();
  // canv_kin_Over->cd(11); canv_kin_Over->cd(11)->SetLogy(true); h_maxtrkdcasig[0]->Draw("hist"); h_maxtrkdcasig[2]->Draw("samehist");leg3->Draw();
  // canv_kin_Over->cd(12); canv_kin_Over->cd(12)->SetLogy(true); h_kminpt[0]->Draw("hist");h_kminpt[2]->Draw("samehist");  leg3->Draw();
  // canv_kin_Over->cd(13); canv_kin_Over->cd(13)->SetLogy(true); h_kmaxpt[0]->Draw("hist");h_kmaxpt[2]->Draw("samehist");  leg3->Draw();
  canv_kin_Over->SaveAs("plot/histogram_Kinematics_Over_Presel.pdf");


  //*/
}

void ploting_kinematicvar(){

  //angular_distr("Modified_sel_BsToPhiMuMu_2016_Rereco07Aug17_Presel_data_cut0_s0.root","Modified_sel_BsToPhiMuMu_OfficialMC_signal_2016_mc.lite_cut0.root");
  angular_distr("sel_BsToPhiMuMu_OfficialMC_signal_2016Mini_combine_Presel_mc.lite_cut0.root","sel_Combine_2016_Mini_Presel_data_cut0.root");
}
