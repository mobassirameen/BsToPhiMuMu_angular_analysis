#include <iostream>
#include <TH1D.h>
#include <TH2D.h>


using namespace std;

void bdtplot() {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  TCanvas* c0=new TCanvas("c0","",700,700);  
  //TPad *pad5 = new TPad("pad5","The pad with the histograms",0.05,0.4,0.95,0.93,21);
  //pad5->cd();
  gPad->SetBottomMargin(0.2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetTopMargin(0.15);

  //TFile smalltree_f("TMVA_withweight_changes_Bmm.root");
  TFile smalltree_f("TMVA_withweight_Defaultsetting.root","READ");

  TTree* test_tree = (TTree*)smalltree_f.Get("TestTree");

  TH2 *h2 = (TH2*)smalltree_f.Get("dataset/CorrelationMatrixS");
  if (!h2) h2 = (TH2*)smalltree_f.Get("CorrelationMatrixS");
  h2->SetLabelSize(0.05, "x");
  h2->SetLabelSize(0.05, "y");
  h2->GetXaxis()->LabelsOption("v");

  h2->Draw("colztext");
  c0->SaveAs("CorrelationMatrixS.pdf");

  h2 = (TH2*)smalltree_f.Get("dataset/CorrelationMatrixB");
  if (!h2) h2 = (TH2*)smalltree_f.Get("CorrelationMatrixB");
  h2->SetLabelSize(0.05, "x");
  h2->SetLabelSize(0.05, "y");
  h2->GetXaxis()->LabelsOption("v");
  h2->Draw("colztext");
  c0->SaveAs("CorrelationMatrixB.pdf");


  //TDirectory* BDT_method = (TDirectory*)smalltree_f.Get("Method_BDT");
  TH1* sig = dynamic_cast<TH1*>(smalltree_f.Get( "dataset/Method_BDT/BDT/MVA_BDT_S" ));
  TH1* bgd = dynamic_cast<TH1*>(smalltree_f.Get( "dataset/Method_BDT/BDT/MVA_BDT_B" ));
  if (sig->GetSumw2N() == 0) sig->Sumw2();
  if (bgd && bgd->GetSumw2N() == 0) bgd->Sumw2();

  if(sig->GetSumOfWeights()!=0) {
    Float_t dx = (sig->GetXaxis()->GetXmax() - sig->GetXaxis()->GetXmin())/sig->GetNbinsX();
    sig->Scale( 1.0/sig->GetSumOfWeights()/dx );
  }
  if (bgd != 0 && bgd->GetSumOfWeights()!=0) {
    Float_t dx = (bgd->GetXaxis()->GetXmax() - bgd->GetXaxis()->GetXmin())/bgd->GetNbinsX();
    bgd->Scale( 1.0/bgd->GetSumOfWeights()/dx );
  }

  //cout << " containing " << hname << "_S/_B" << endl;
  // chop off useless stuff
  sig->SetTitle( Form("TMVA overtraining check for classifier: "));
  sig->SetLineColor(kBlue);
  bgd->SetLineColor(kRed);
  static Int_t c_SignalLine     = TColor::GetColor( "#0000ee" );
  static Int_t c_SignalFill     = TColor::GetColor( "#7d99d1" );
  static Int_t c_BackgroundLine = TColor::GetColor( "#ff0000" );
  static Int_t c_BackgroundFill = TColor::GetColor( "#ff0000" );

  Int_t FillColor__B = c_BackgroundFill;
  Int_t FillStyle__B = 3365; //3554;
  Int_t LineColor__B = c_BackgroundLine;
  Int_t LineWidth__B = 1;

  Int_t FillColor__S = c_SignalFill;
  Int_t FillStyle__S = 3356; //1001;
  Int_t LineColor__S = c_SignalLine;
  Int_t LineWidth__S = 1;

  sig->SetMarkerColor(LineColor__S );
  sig->SetLineColor( LineColor__S );
  sig->SetLineWidth( LineWidth__S );
  sig->SetFillStyle( FillStyle__S );
  sig->SetFillColor( FillColor__S );
  bgd->SetMarkerColor(LineColor__B );
  bgd->SetLineColor( LineColor__B );
  bgd->SetLineWidth( LineWidth__B );
  bgd->SetFillStyle( FillStyle__B );
  bgd->SetFillColor( FillColor__B );
  TCanvas* c=new TCanvas("c","",700,700); 
 
  // create new canvas
  Float_t nrms = 10;
  cout << "--- Mean and RMS (S): " << sig->GetMean() << ", " << sig->GetRMS() << endl;
  cout << "--- Mean and RMS (B): " << bgd->GetMean() << ", " << bgd->GetRMS() << endl;
  Float_t xmin = TMath::Max( TMath::Min(sig->GetMean() - nrms*sig->GetRMS(),
					bgd->GetMean() - nrms*bgd->GetRMS() ),
			     sig->GetXaxis()->GetXmin() );
  Float_t xmax = TMath::Min( TMath::Max(sig->GetMean() + nrms*sig->GetRMS(),
					bgd->GetMean() + nrms*bgd->GetRMS() ),
			     sig->GetXaxis()->GetXmax() );
  Float_t ymin = 0;
  Float_t maxMult = 1.5;
  Float_t ymax = TMath::Max( sig->GetMaximum(), bgd->GetMaximum() )*maxMult;
  //change the bin
  xmin = -0.5;
  xmax = 0.5;

  Int_t nb = 500;
  TH2F* frame = new TH2F("frame", sig->GetTitle(),
			  nb, xmin, xmax, nb, ymin, ymax );
  frame->GetXaxis()->SetTitle(  "BDT response"  );
  frame->GetYaxis()->SetTitle("(1/N) dN^{ }/^{ }dx");
  //SetFrameStyle( frame );

  // eventually: draw the frame
  frame->Draw();

  c->GetPad(0)->SetLeftMargin( 0.105 );
  frame->GetYaxis()->SetTitleOffset( 1.2 );

  // Draw legend
  TLegend *legend= new TLegend( c->GetLeftMargin(), 1 - c->GetTopMargin() - 0.12,
				c->GetLeftMargin() +  0.40, 1 - c->GetTopMargin() );
  legend->SetFillStyle( 1 );
  legend->AddEntry(sig,TString("Signal")     + " (test sample)", "F");
  legend->AddEntry(bgd,TString("Background") + " (test sample)", "F");
  legend->SetBorderSize(1);
  legend->SetMargin(0.2);
  legend->Draw("same");

  // overlay signal and background histograms
  sig->Draw("samehist");
  bgd->Draw("samehist");

  TH1* sigOv = 0;
  TH1* bgdOv = 0;


  sigOv = dynamic_cast<TH1*>(smalltree_f.Get( "dataset/Method_BDT/BDT/MVA_BDT_Train_S" ));
  bgdOv = dynamic_cast<TH1*>(smalltree_f.Get( "dataset/Method_BDT/BDT/MVA_BDT_Train_B" ));


  if (sigOv == 0 || bgdOv == 0) {
    cout << "+++ Problem in \"mvas.C\": overtraining check histograms do not exist" << endl;
  }
  else {
    cout << "--- Found comparison histograms for overtraining check" << endl;

    TLegend *legend2= new TLegend( 1 - c->GetRightMargin() - 0.42, 1 - c->GetTopMargin() - 0.12,
				   1 - c->GetRightMargin(), 1 - c->GetTopMargin() );
    legend2->SetFillStyle( 1 );
    legend2->SetBorderSize(1);
    legend2->AddEntry(sigOv,"Signal (training sample)","P");
    legend2->AddEntry(bgdOv,"Background (training sample)","P");
    legend2->SetMargin( 0.1 );
    legend2->Draw("same");
  }
  // normalise both signal and background
  
  if (sigOv->GetSumw2N() == 0) sigOv->Sumw2();
  if (bgdOv && bgdOv->GetSumw2N() == 0) bgdOv->Sumw2();

  if(sigOv->GetSumOfWeights()!=0) {
    Float_t dx = (sigOv->GetXaxis()->GetXmax() - sigOv->GetXaxis()->GetXmin())/sigOv->GetNbinsX();
    sigOv->Scale( 1.0/sigOv->GetSumOfWeights()/dx );
  }
  if (bgdOv != 0 && bgdOv->GetSumOfWeights()!=0) {
    Float_t dx = (bgdOv->GetXaxis()->GetXmax() - bgdOv->GetXaxis()->GetXmin())/bgdOv->GetNbinsX();
    bgdOv->Scale( 1.0/bgdOv->GetSumOfWeights()/dx );
  }

  Int_t col = sig->GetLineColor();
  sigOv->SetMarkerColor( col );
  sigOv->SetMarkerSize( 0.7 );
  sigOv->SetMarkerStyle( 20 );
  sigOv->SetLineWidth( 1 );
  sigOv->SetLineColor( col );
  sigOv->Draw("e1same");

  col = bgd->GetLineColor();
  bgdOv->SetMarkerColor( col );
  bgdOv->SetMarkerSize( 0.7 );
  bgdOv->SetMarkerStyle( 20 );
  bgdOv->SetLineWidth( 1 );
  bgdOv->SetLineColor( col );
  bgdOv->Draw("e1same");
  Double_t kolS = sig->KolmogorovTest( sigOv );
  Double_t kolB = bgd->KolmogorovTest( bgdOv );
  cout << "--- Goodness of signal (background) consistency: " << kolS << " (" << kolB << ")" << endl;
  
  //gPad->SetLogy();
  c->SaveAs("overtraining.pdf");
}

