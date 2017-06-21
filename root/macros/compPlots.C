// Macro to make plots comparing weighted MC events and data
// Histograms were produced via. mc_reweight.C
// Author: Eric Pooser, pooser@jlab.org

// Input ROOT file created via mc_reweight.C
TFile *inFile;
// Input ROOT file directories
TDirectory *mcWgtDir, *dataDir;
// Weighted 1D MC histos
TH1F *h_xFocalMCWgt, *h_xpFocalMCWgt, *h_yFocalMCWgt, *h_ypFocalMCWgt;
TH1F *h_yTarMCWgt, *h_xpTarMCWgt, *h_ypTarMCWgt, *h_deltaMCWgt;
TH1F *h_thetaMCWgt, *h_q2MCWgt, *h_w2MCWgt;
// Weighted 2D MC histos
TH2F *h2_xVxpFocalMCWgt, *h2_xVyFocalMCWgt, *h2_xVypFocalMCWgt;
TH2F *h2_xpVyFocalMCWgt, *h2_xpVypFocalMCWgt, *h2_yVypFocalMCWgt;
// Data 1D histos
TH1F *h_xFocalData, *h_xpFocalData, *h_yFocalData, *h_ypFocalData;
TH1F *h_yTarData, *h_xpTarData, *h_ypTarData, *h_deltaData;
TH1F *h_thetaData, *h_q2Data, *h_w2Data;
TH1F *h_yTarDataClone, *h_xpTarDataClone, *h_ypTarDataClone, *h_deltaDataClone;
TH1F *h_thetaDataClone, *h_q2DataClone, *h_w2DataClone;
// Data 2D histos
TH2F *h2_xVxpFocalData, *h2_xVyFocalData, *h2_xVypFocalData;
TH2F *h2_xpVyFocalData, *h2_xpVypFocalData, *h2_yVypFocalData;
// Comparison canvas'
TCanvas *c_tarComp, *c_kinComp, *c_focalComp;
// Legends
TLegend *l_yTarComp, *l_xpTarComp, *l_ypTarComp, *l_deltaComp;
TLegend *l_thetaComp, *l_q2Comp, *l_w2Comp;
// Lines
TLine *ln_yTarComp, *ln_xpTarComp, *ln_ypTarComp, *ln_deltaComp;
TLine *ln_thetaComp, *ln_q2Comp, *ln_w2Comp;

// Define constants
// Canvas size paramters
static const Double_t canWidth   = 1600.0;
static const Double_t canHeight  = 800.0;
// X-axis limits
static const Double_t yTarXMin   = -2.0;
static const Double_t yTarXMax   = 2.0;
static const Double_t xpTarXMin  = -50.0;
static const Double_t xpTarXMax  = 50.0;
static const Double_t ypTarXMin  = -30.0;
static const Double_t ypTarXMax  = 30.0;
static const Double_t deltaXMin  = -10.0;
static const Double_t deltaXMax  = 20.0;
static const Double_t thetaXMin  = 13.5;
static const Double_t thetaXMax  = 16.5;
static const Double_t q2XMin     = 0.95;
static const Double_t q2XMax     = 1.9;
static const Double_t w2XMin     = 4.2;
static const Double_t w2XMax     = 6.9;
// Y-axis limits
static const Double_t yTarYMin   = 0.0;
static const Double_t yTarYMax   = 2300.0;
static const Double_t xpTarYMin  = 0.0;
static const Double_t xpTarYMax  = 450.0;
static const Double_t ypTarYMin  = 0.0;
static const Double_t ypTarYMax  = 900.0;
static const Double_t deltaYMin  = 0.0;
static const Double_t deltaYMax  = 900.0;
static const Double_t thetaYMin  = 0.0;
static const Double_t thetaYMax  = 800.0;
static const Double_t q2YMin     = 0.0;
static const Double_t q2YMax     = 700.0;
static const Double_t w2YMin     = 0.0;
static const Double_t w2YMax     = 450.0;
// Ratio limits
static const Double_t idealRatio = 1.0;
static const Double_t ratioMin   = 0.0;
static const Double_t ratioMax   = 2.0;
// Scale factor for MC events to match data
static const Double_t scaleFactor = 1.5;


void compPlots() {

  // Global ROOT settings
  gStyle->SetOptStat(0);
  
  // Obtain the comparison ROOT file
  inFile = new TFile("../output/run_488_comp.root", "READ");
  // Obtain the directories which contain the histograms of interest
  mcWgtDir = dynamic_cast <TDirectory*> (inFile->FindObjectAny("mcWgtDir"));
  dataDir  = dynamic_cast <TDirectory*> (inFile->FindObjectAny("dataDir"));
  // Obtain the 1D histograms of interest
  h_yTarMCWgt  = dynamic_cast <TH1F*> (mcWgtDir->FindObjectAny("h_yTarMCWgt"));
  h_xpTarMCWgt = dynamic_cast <TH1F*> (mcWgtDir->FindObjectAny("h_xpTarMCWgt"));
  h_ypTarMCWgt = dynamic_cast <TH1F*> (mcWgtDir->FindObjectAny("h_ypTarMCWgt"));
  h_deltaMCWgt = dynamic_cast <TH1F*> (mcWgtDir->FindObjectAny("h_deltaMCWgt"));
  h_thetaMCWgt = dynamic_cast <TH1F*> (mcWgtDir->FindObjectAny("h_thetaMCWgt"));
  h_q2MCWgt    = dynamic_cast <TH1F*> (mcWgtDir->FindObjectAny("h_q2MCWgt"));
  h_w2MCWgt    = dynamic_cast <TH1F*> (mcWgtDir->FindObjectAny("h_w2MCWgt"));
  h_yTarData   = dynamic_cast <TH1F*> (dataDir->FindObjectAny("h_yTarData"));
  h_xpTarData  = dynamic_cast <TH1F*> (dataDir->FindObjectAny("h_xpTarData"));
  h_ypTarData  = dynamic_cast <TH1F*> (dataDir->FindObjectAny("h_ypTarData"));
  h_deltaData  = dynamic_cast <TH1F*> (dataDir->FindObjectAny("h_deltaData"));
  h_thetaData  = dynamic_cast <TH1F*> (dataDir->FindObjectAny("h_thetaData"));
  h_q2Data     = dynamic_cast <TH1F*> (dataDir->FindObjectAny("h_q2Data"));
  h_w2Data     = dynamic_cast <TH1F*> (dataDir->FindObjectAny("h_w2Data"));
  // Obtain the 2D histos of interest
  h2_xVxpFocalMCWgt  = dynamic_cast <TH2F*> (mcWgtDir->FindObjectAny("h2_xVxpFocalMCWgt"));
  h2_xVyFocalMCWgt   = dynamic_cast <TH2F*> (mcWgtDir->FindObjectAny("h2_xVyFocalMCWgt"));
  h2_xVypFocalMCWgt  = dynamic_cast <TH2F*> (mcWgtDir->FindObjectAny("h2_xVypFocalMCWgt"));
  h2_xpVyFocalMCWgt  = dynamic_cast <TH2F*> (mcWgtDir->FindObjectAny("h2_xpVyFocalMCWgt"));
  h2_xpVypFocalMCWgt = dynamic_cast <TH2F*> (mcWgtDir->FindObjectAny("h2_xpVypFocalMCWgt"));
  h2_yVypFocalMCWgt  = dynamic_cast <TH2F*> (mcWgtDir->FindObjectAny("h2_yVypFocalMCWgt"));
  h2_xVxpFocalData   = dynamic_cast <TH2F*> (dataDir->FindObjectAny("h2_xVxpFocalData"));
  h2_xVyFocalData    = dynamic_cast <TH2F*> (dataDir->FindObjectAny("h2_xVyFocalData"));
  h2_xVypFocalData   = dynamic_cast <TH2F*> (dataDir->FindObjectAny("h2_xVypFocalData"));
  h2_xpVyFocalData   = dynamic_cast <TH2F*> (dataDir->FindObjectAny("h2_xpVyFocalData"));
  h2_xpVypFocalData  = dynamic_cast <TH2F*> (dataDir->FindObjectAny("h2_xpVypFocalData"));
  h2_yVypFocalData   = dynamic_cast <TH2F*> (dataDir->FindObjectAny("h2_yVypFocalData"));
  // Clone histos for calculating ratios
  h_yTarDataClone  = dynamic_cast <TH1F*> (h_yTarData->Clone());
  h_xpTarDataClone = dynamic_cast <TH1F*> (h_xpTarData->Clone());
  h_ypTarDataClone = dynamic_cast <TH1F*> (h_ypTarData->Clone());
  h_deltaDataClone = dynamic_cast <TH1F*> (h_deltaData->Clone());
  h_thetaDataClone = dynamic_cast <TH1F*> (h_thetaData->Clone());
  h_q2DataClone    = dynamic_cast <TH1F*> (h_q2Data->Clone());
  h_w2DataClone    = dynamic_cast <TH1F*> (h_w2Data->Clone());

  c_tarComp = new TCanvas("c_tarComp", "Comparison of Target Variables", canWidth, canHeight); 
  c_tarComp->Divide(3, 2);
  c_tarComp->cd(1);
  h_xpTarMCWgt->Scale(scaleFactor);
  h_xpTarMCWgt->SetOption("HIST");
  h_xpTarMCWgt->SetFillColor(8);
  h_xpTarMCWgt->SetLineColor(8);
  h_xpTarMCWgt->SetFillStyle(3001);
  h_xpTarMCWgt->GetXaxis()->SetRangeUser(xpTarXMin, xpTarXMax);
  h_xpTarMCWgt->GetYaxis()->SetRangeUser(xpTarYMin, xpTarYMax);
  h_xpTarMCWgt->SetTitle("Weighted MC to Data Comparison: X'_{tar}");
  h_xpTarMCWgt->Draw();
  h_xpTarData->SetMarkerColor(4);
  h_xpTarData->SetMarkerStyle(22);
  h_xpTarData->Draw("E, SAME");
  l_xpTarComp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_xpTarComp->AddEntry(h_xpTarMCWgt, "Weighted MC", "F");
  l_xpTarComp->AddEntry(h_xpTarData, "Data", "EP");
  l_xpTarComp->Draw();

  c_tarComp->cd(2);
  h_ypTarMCWgt->Scale(scaleFactor);
  h_ypTarMCWgt->SetOption("HIST");
  h_ypTarMCWgt->SetFillColor(8);
  h_ypTarMCWgt->SetLineColor(8);
  h_ypTarMCWgt->SetFillStyle(3001);
  h_ypTarMCWgt->GetXaxis()->SetRangeUser(ypTarXMin, ypTarXMax);
  h_ypTarMCWgt->GetYaxis()->SetRangeUser(ypTarYMin, ypTarYMax);
  h_ypTarMCWgt->SetTitle("Weighted MC to Data Comparison: Y'_{tar}");
  h_ypTarMCWgt->Draw();
  h_ypTarData->SetMarkerColor(4);
  h_ypTarData->SetMarkerStyle(22);
  h_ypTarData->Draw("E, SAME");
  l_ypTarComp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_ypTarComp->AddEntry(h_ypTarMCWgt, "Weighted MC", "F");
  l_ypTarComp->AddEntry(h_ypTarData, "Data", "EP");
  l_ypTarComp->Draw();

  c_tarComp->cd(3);
  h_deltaMCWgt->Scale(scaleFactor);
  h_deltaMCWgt->SetOption("HIST");
  h_deltaMCWgt->SetFillColor(8);
  h_deltaMCWgt->SetLineColor(8);
  h_deltaMCWgt->SetFillStyle(3001);
  h_deltaMCWgt->GetXaxis()->SetRangeUser(deltaXMin, deltaXMax);
  h_deltaMCWgt->GetYaxis()->SetRangeUser(deltaYMin, deltaYMax);
  h_deltaMCWgt->SetTitle("Weighted MC to Data Comparison: #delta");
  h_deltaMCWgt->Draw();
  h_deltaData->SetMarkerColor(4);
  h_deltaData->SetMarkerStyle(22);
  h_deltaData->Draw("E, SAME");
  l_deltaComp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_deltaComp->AddEntry(h_deltaMCWgt, "Weighted MC", "F");
  l_deltaComp->AddEntry(h_deltaData, "Data", "EP");
  l_deltaComp->Draw();

  c_tarComp->cd(4);
  h_xpTarDataClone->Divide(h_xpTarMCWgt);
  h_xpTarDataClone->SetMarkerStyle(22);
  h_xpTarDataClone->SetMarkerSize(1.25);
  h_xpTarDataClone->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_xpTarDataClone->GetXaxis()->SetRangeUser(xpTarXMin, xpTarXMax);
  h_xpTarDataClone->GetYaxis()->SetTitle("X'_{tar}^{data} / X'_{tar}^{MC}");
  h_xpTarDataClone->SetTitle("Ratio of X'_{tar}^{data} to X'_{tar}^{MC}");
  h_xpTarDataClone->Draw();
  ln_xpTarComp = new TLine(xpTarXMin, idealRatio, xpTarXMax, idealRatio);
  ln_xpTarComp->SetLineStyle(9);
  ln_xpTarComp->SetLineWidth(2);
  ln_xpTarComp->SetLineColor(38);
  ln_xpTarComp->Draw();

  c_tarComp->cd(5);
  h_ypTarDataClone->Divide(h_ypTarMCWgt);
  h_ypTarDataClone->SetMarkerStyle(22);
  h_ypTarDataClone->SetMarkerSize(1.25);
  h_ypTarDataClone->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_ypTarDataClone->GetXaxis()->SetRangeUser(ypTarXMin, ypTarXMax);
  h_ypTarDataClone->GetYaxis()->SetTitle("Y'_{tar}^{data} / Y'_{tar}^{MC}");
  h_ypTarDataClone->SetTitle("Ratio of Y'_{tar}^{data} to Y'_{tar}^{MC}");
  h_ypTarDataClone->Draw();
  ln_ypTarComp = new TLine(ypTarXMin, idealRatio, ypTarXMax, idealRatio);
  ln_ypTarComp->SetLineStyle(9);
  ln_ypTarComp->SetLineWidth(2);
  ln_ypTarComp->SetLineColor(38);
  ln_ypTarComp->Draw();

  c_tarComp->cd(6);
  h_deltaDataClone->Divide(h_deltaMCWgt);
  h_deltaDataClone->SetMarkerStyle(22);
  h_deltaDataClone->SetMarkerSize(1.25);
  h_deltaDataClone->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_deltaDataClone->GetXaxis()->SetRangeUser(deltaXMin, deltaXMax);
  h_deltaDataClone->GetYaxis()->SetTitle("#delta_{data} / #delta_{MC}");
  h_deltaDataClone->SetTitle("Ratio of #delta_{data} to #delta_{MC}");
  h_deltaDataClone->Draw();
  ln_deltaComp = new TLine(deltaXMin, idealRatio, deltaXMax, idealRatio);
  ln_deltaComp->SetLineStyle(9);
  ln_deltaComp->SetLineWidth(2);
  ln_deltaComp->SetLineColor(38);
  ln_deltaComp->Draw();

  c_kinComp = new TCanvas("c_kinComp", "Comparison of Kinematic Variables", canWidth, canHeight); 
  c_kinComp->Divide(3, 2);
  c_kinComp->cd(1);
  h_thetaMCWgt->Scale(scaleFactor);
  h_thetaMCWgt->SetOption("HIST");
  h_thetaMCWgt->SetFillColor(8);
  h_thetaMCWgt->SetLineColor(8);
  h_thetaMCWgt->SetFillStyle(3001);
  h_thetaMCWgt->GetXaxis()->SetRangeUser(thetaXMin, thetaXMax);
  h_thetaMCWgt->GetYaxis()->SetRangeUser(thetaYMin, thetaYMax);
  h_thetaMCWgt->SetTitle("Weighted MC to Data Comparison: #theta");
  h_thetaMCWgt->Draw();
  h_thetaData->SetMarkerColor(4);
  h_thetaData->SetMarkerStyle(22);
  h_thetaData->Draw("E, SAME");
  l_thetaComp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_thetaComp->AddEntry(h_thetaMCWgt, "Weighted MC", "F");
  l_thetaComp->AddEntry(h_thetaData, "Data", "EP");
  l_thetaComp->Draw();

  c_kinComp->cd(2);
  h_q2MCWgt->Scale(scaleFactor);
  h_q2MCWgt->SetOption("HIST");
  h_q2MCWgt->SetFillColor(8);
  h_q2MCWgt->SetLineColor(8);
  h_q2MCWgt->SetFillStyle(3001);
  h_q2MCWgt->GetXaxis()->SetRangeUser(q2XMin, q2XMax);
  h_q2MCWgt->GetYaxis()->SetRangeUser(q2YMin, q2YMax);
  h_q2MCWgt->SetTitle("Weighted MC to Data Comparison: Q^{2}");
  h_q2MCWgt->Draw();
  h_q2Data->SetMarkerColor(4);
  h_q2Data->SetMarkerStyle(22);
  h_q2Data->Draw("E, SAME");
  l_q2Comp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_q2Comp->AddEntry(h_q2MCWgt, "Weighted MC", "F");
  l_q2Comp->AddEntry(h_q2Data, "Data", "EP");
  l_q2Comp->Draw();

  c_kinComp->cd(3);
  h_w2MCWgt->Scale(scaleFactor);
  h_w2MCWgt->SetOption("HIST");
  h_w2MCWgt->SetFillColor(8);
  h_w2MCWgt->SetLineColor(8);
  h_w2MCWgt->SetFillStyle(3001);
  h_w2MCWgt->GetXaxis()->SetRangeUser(w2XMin, w2XMax);
  h_w2MCWgt->GetYaxis()->SetRangeUser(w2YMin, w2YMax);
  h_w2MCWgt->SetTitle("Weighted MC to Data Comparison: W^{2}");
  h_w2MCWgt->Draw();
  h_w2Data->SetMarkerColor(4);
  h_w2Data->SetMarkerStyle(22);
  h_w2Data->Draw("E, SAME");
  l_w2Comp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_w2Comp->AddEntry(h_w2MCWgt, "Weighted MC", "F");
  l_w2Comp->AddEntry(h_w2Data, "Data", "EP");
  l_w2Comp->Draw();

  c_kinComp->cd(4);
  h_thetaDataClone->Divide(h_thetaMCWgt);
  h_thetaDataClone->SetMarkerStyle(22);
  h_thetaDataClone->SetMarkerSize(1.25);
  h_thetaDataClone->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_thetaDataClone->GetXaxis()->SetRangeUser(thetaXMin, thetaXMax);
  h_thetaDataClone->GetYaxis()->SetTitle("#theta_{data} / #theta_{MC}");
  h_thetaDataClone->SetTitle("Ratio of #theta_{data} to #theta_{MC}");
  h_thetaDataClone->Draw();
  ln_thetaComp = new TLine(thetaXMin, idealRatio, thetaXMax, idealRatio);
  ln_thetaComp->SetLineStyle(9);
  ln_thetaComp->SetLineWidth(2);
  ln_thetaComp->SetLineColor(38);
  ln_thetaComp->Draw();

  c_kinComp->cd(5);
  h_q2DataClone->Divide(h_q2MCWgt);
  h_q2DataClone->SetMarkerStyle(22);
  h_q2DataClone->SetMarkerSize(1.25);
  h_q2DataClone->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_q2DataClone->GetXaxis()->SetRangeUser(q2XMin, q2XMax);
  h_q2DataClone->GetYaxis()->SetTitle("Q^{2}_{data} / Q^{2}_{MC}");
  h_q2DataClone->SetTitle("Ratio of Q^{2}_{data} to Q^{2}_{MC}");
  h_q2DataClone->Draw();
  ln_q2Comp = new TLine(q2XMin, idealRatio, q2XMax, idealRatio);
  ln_q2Comp->SetLineStyle(9);
  ln_q2Comp->SetLineWidth(2);
  ln_q2Comp->SetLineColor(38);
  ln_q2Comp->Draw();

  c_kinComp->cd(6);
  h_w2DataClone->Divide(h_w2MCWgt);
  h_w2DataClone->SetMarkerStyle(22);
  h_w2DataClone->SetMarkerSize(1.25);
  h_w2DataClone->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_w2DataClone->GetXaxis()->SetRangeUser(w2XMin, w2XMax);
  h_w2DataClone->GetYaxis()->SetTitle("W^{2}_{data} / W^{2}_{MC}");
  h_w2DataClone->SetTitle("Ratio of W^{2}_{data} to W^{2}_{MC}");
  h_w2DataClone->Draw();
  ln_w2Comp = new TLine(w2XMin, idealRatio, w2XMax, idealRatio);
  ln_w2Comp->SetLineStyle(9);
  ln_w2Comp->SetLineWidth(2);
  ln_w2Comp->SetLineColor(38);
  ln_w2Comp->Draw();

  c_focalComp = new TCanvas("c_focalComp", "Comparison of Focal Plane Quantites", canWidth, canHeight); 
  c_focalComp->Divide(4, 3);
  c_focalComp->cd(1); gPad->SetLogz();
  h2_xVxpFocalMCWgt->Draw("COLZ");
  c_focalComp->cd(2); gPad->SetLogz();
  h2_xVxpFocalData->Draw("COLZ");
  c_focalComp->cd(3); gPad->SetLogz();
  h2_xVyFocalMCWgt->Draw("COLZ");
  c_focalComp->cd(4); gPad->SetLogz();
  h2_xVyFocalData->Draw("COLZ");
  c_focalComp->cd(5); gPad->SetLogz();
  h2_xVypFocalMCWgt->Draw("COLZ");
  c_focalComp->cd(6); gPad->SetLogz();
  h2_xVypFocalData->Draw("COLZ");
  c_focalComp->cd(7); gPad->SetLogz();
  h2_xpVyFocalMCWgt->Draw("COLZ");
  c_focalComp->cd(8); gPad->SetLogz();
  h2_xpVyFocalData->Draw("COLZ");
  c_focalComp->cd(9); gPad->SetLogz();
  h2_xpVypFocalMCWgt->Draw("COLZ");
  c_focalComp->cd(10); gPad->SetLogz();
  h2_xpVypFocalData->Draw("COLZ");
  c_focalComp->cd(11); gPad->SetLogz();
  h2_yVypFocalMCWgt->Draw("COLZ");
  c_focalComp->cd(12); gPad->SetLogz();
  h2_yVypFocalData->Draw("COLZ");
}
