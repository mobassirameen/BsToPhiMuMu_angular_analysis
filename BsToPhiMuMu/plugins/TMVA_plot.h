#ifndef __TMVA_plot__h
#define __TMVA_plot__h


#include <TROOT.h>
#include "TMath.h"
#include <iostream>

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "TPaveText.h"
#include "TLegend.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>

#include <RooRandom.h>
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLatex.h"

#include "RooExponential.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TCut.h"
#include <TRandom3.h>
#include "RooRandom.h"
#include "fstream"



using namespace std;
using namespace RooFit;


double iF_gauss2(double *x, double *par) {
  // par[0] -> const                                                                                                                                                        
  // par[1] -> mean                                                                                                                                                         
  // par[2] -> sigma                                                                                                                                                        
  // par[3] -> fraction in second gaussian                                                                                                                                  
  // par[4] -> mean of second gaussian                                                                                                                                      
  // par[5] -> sigma of second gaussian                                                                                                                                     
  Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.);
  if (par[2] > 0) {
    arg1 = (x[0] - par[1]) / par[2];
    fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
  }
  if (par[5] > 0.) {
    arg2 = (x[0] - par[4]) / par[5];
    fitval2 =  par[3]*par[0]*TMath::Exp(-0.5*arg2*arg2);
  }
  Double_t fitval = fitval1 + fitval2;
  return fitval;
}
double iF_expo(double *x, double *par) {
  return par[0]*TMath::Exp(x[0]*par[1]);
}

double iF_Gauss(double *x, double *par) {
  // par[0] -> area                                                                                                                                                         
  // par[1] -> mean                                                                                                                                                         
  // par[2] -> sigma                                                                                                                                                        

  double sqrt2pi = 2.506628275;

  if (par[2] > 0.) {
    Double_t arg = (x[0] - par[1]) / par[2];
    Double_t fitval =  (par[0]/(sqrt2pi*par[2])) * TMath::Exp(-0.5*arg*arg);
    return fitval;
  }
  else {
    return -1.;
  }
}
double iF_gauss2c(double *x, double *par) {
  // constrained to have the same mean in the second gaussian                                                                                                               
  // par[0] -> const                                                                                                                                                        
  // par[1] -> mean                                                                                                                                                         
  // par[2] -> sigma                                                                                                                                                        
  // par[3] -> fraction in second gaussian                                                                                                                                  
  // par[4] -> sigma of second gaussian                                                                                                                                     
  Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.);
  if (par[2] > 0) {
    arg1 = (x[0] - par[1]) / par[2];
    fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
  }
  if (par[4] > 0.) {
    arg2 = (x[0] - par[1]) / par[4];
    fitval2 =  par[3]*par[0]*TMath::Exp(-0.5*arg2*arg2);
  }
  Double_t fitval = fitval1 + fitval2;
  return fitval;
}

double iF_expo_Gauss(double *x, double *par) {
  //   par[0] = normalization of gaussian                                                                                                                                   
  //   par[1] = mean of gaussian                                                                                                                                            
  //   par[2] = sigma of gaussian                                                                                                                                           
  //   par[3] = base of expo                                                                                                                                                
  //   par[4] = exponent of expo                                                                                                                                            
  return  (iF_expo(x, &par[3]) + iF_Gauss(x, &par[0]));
}

double iF_expo_gauss2(double *x, double *par) {
  // par[0] -> const                                                                                                                                                        
  // par[1] -> mean                                                                                                                                                         
  // par[2] -> sigma                                                                                                                                                        
  // par[3] -> fraction in second gaussian                                                                                                                                  
  // par[4] -> mean of second gaussian                                                                                                                                      
  // par[5] -> sigma of second gaussian                                                                                                                                     
  // par[6] = base                                                                                                                                                          
  // par[7] = exp                                                                                                                                                           
  return  (iF_expo(x, &par[5]) + iF_gauss2c(x, &par[0]));
}

void canv_withrat(TH1D* h1, TH1D* h2, TH1D* h3, TString var, TString year){

  h1->SetMaximum(1.5*h2->GetMaximum());
  h1->SetMinimum(0.);
  TCanvas canv_set("canv_Set","",1200,600);
  canv_set.Divide(2,1);
  canv_set.cd(1); h1->Draw();
  canv_set.cd(2); h2->Draw();
  canv_set.SaveAs(Form("plot/1D_Jpsi_plot_distribution_%s.pdf",var.Data()));
  canv_set.Close();


  TCanvas* lsbdt=new TCanvas("lsbdt","",600,600);
  lsbdt->cd();
  TPad *pad5 = new TPad("pad5","The pad with the histograms",0.05,0.4,0.95,0.93,21);
  TPad *pad6 = new TPad("pad6","The pad with the ratio",0.05,0.05,0.95,0.4,21);
  pad5->Draw();
  pad6->Draw();
  pad5->cd();
  pad5->SetBottomMargin(0);
  pad5->SetFillColor(10);
  pad5->SetGridx(1);
  pad5->SetGridy(1);
  pad5->SetTicky(1);
  pad5->SetTickx(1);
  //gPad->SetLogy();                                                                                                                                                      
  //  lsbplot->Draw();                                                                                                                                                     
  h1->SetTitle(" ");
  h1->GetYaxis()->SetTitle("Number of Events");
  h1->Draw("hist");
  h2->SetMarkerStyle(20);
  h2->Draw("same");
  TString title4 = " 35.9 fb^{-1}( 13 TeV)";
  TString title = "#font[61]{CMS}";
  TString title2= "#font[52]{Preliminary}";
  TLatex tex;
  tex.SetTextFont(42);
  tex.SetTextSize(0.035);
  tex.SetTextAlign(11);
  tex.SetNDC();
  tex.SetTextSize(0.035);
  tex.DrawLatex(0.14,0.94,title);
  tex.SetTextSize(0.028);
  tex.DrawLatex(0.22,0.94,title2);
  tex.DrawLatex(0.74,0.94,title4);
  tex.DrawLatex(0.54,0.94,year);
  TLegend *legs = new TLegend(0.60,0.70,0.91,0.91);
  legs->SetTextFont(43);
  legs->SetTextSize(16);
  legs->SetFillColor(0);
  legs->SetBorderSize(0);
  legs->SetFillStyle(0);
  legs->AddEntry(h2, "Data", "lep");
  legs->AddEntry(h1, "MC", "f");
  legs->Draw();

  pad6->cd();
  pad6->SetBottomMargin(0.2);
  pad6->SetTopMargin(0);
  pad6->SetTickx();
  pad6->SetFillColor(0);
  pad6->SetGridx(0);
  pad6->SetGridy(0);
  pad6->SetTitle(" ");
  //  double lo_edge=h3->getBinLowedge                                                                                                                                      
  double fLo   = h3->GetXaxis()->GetXmin();
  double fHi   = h3->GetXaxis()->GetXmax();

  cout<<" low "<<fLo<<"\t high "<<fHi<<endl;

  TBox *box1 = new TBox(fLo,0.8, fHi,1.2);
  box1->SetFillColor(kGreen);
  box1->SetFillStyle(3335);
  box1->SetLineColor(kGreen);
  h3->SetTitle(" ");
  h3->GetYaxis()->SetTitle("Data/MC");
  h3->GetYaxis()->CenterTitle();
  h3->GetXaxis()->CenterTitle();
  //  h3->SetTitleSize(0.1,"x");                                                                                                                                            
  //  h3->GetXaxis()->SetTitle(Form("%s",var.Data()));                                                                                                                      
  h3->GetYaxis()->SetTitleSize(0.043);
  h3->GetXaxis()->SetTitleSize(0.09);
  h3->SetMinimum(0.5);
  h3->SetMaximum(1.5);
  h3->SetMarkerStyle(20);
  h3->Draw();
  box1->Draw();
  lsbdt->SaveAs(Form("plot/Plot_bdt_comp_%s_Jpsi_%s.pdf",var.Data(), year.Data()));

  lsbdt->Close();
  lsbdt->Clear();

}

double Low_mass=5.05;
double High_mass=5.7;

double signal_lo=5.28;
double signal_hi=5.46;

double Lo_SB_lo=5.14;
double Lo_SB_hi=5.22;

double Hi_SB_lo=5.51;
double Hi_SB_hi=5.59;

#endif
