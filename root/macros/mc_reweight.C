#include <TTree.h>
#include <TBranch.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TF1.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TVectorD.h>
#include <TNtuple.h>

#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;

// Declare data files
TFile *mcRawFile, *mcWgtFile, *dataInFile, *compFile;

// Declare trees
TTree *mcRawTree, *mcWgtTree, *dataTree;

// Declare text data files
ifstream reconDataFile, inputDataFile, targetDataFile;

// Declare variables
// Input text file data
Float_t target, E0, P, theta, psFactors, bcm1, bcm2, q1, q2, clt, elt, trackEff, tofEff, trigEff;
Float_t numEvents, phiRange, thetaRange, dPMax, dPmin, chargeSymmFlag, radCorrFlag;
Float_t tarID, tarDensity, tarThickness, tarAtomicMass, tarAtomicNum, tarArealDensity;
// Number of entries and scaling
Int_t    nEntriesMCRaw, nEntriesMCWgt, nEntriesData;
// Leaf variables
Double_t phgcerNpeSum, pngcerNpeSum, paeroNpeSum; 
Double_t numTracks, goldTrackBeta, trackEnergyNorm;
Float_t  xFocalMCRaw, xpFocalMCRaw, yFocalMCRaw, ypFocalMCRaw;
Float_t  yTarMCRaw, xpTarMCRaw, ypTarMCRaw, deltaMCRaw;
Float_t  failID, bornCS, intRadCorr;
Double_t scaleFactor;
Float_t  xFocalMCWgt, xpFocalMCWgt, yFocalMCWgt, ypFocalMCWgt;
Float_t  yTarMCWgt, xpTarMCWgt, ypTarMCWgt, deltaMCWgt;
Double_t xFocalData, xpFocalData, yFocalData, ypFocalData;
Double_t yTarData, xpTarData, ypTarData, deltaData;

// Declare histos
// MC histos
TH1F *h_xFocalMCRaw, *h_xpFocalMCRaw, *h_yFocalMCRaw, *h_ypFocalMCRaw;
TH1F *h_yTarMCRaw, *h_xpTarMCRaw, *h_ypTarMCRaw, *h_deltaMCRaw;
TH2F *h2_xVxpFocalMCRaw, *h2_xVyFocalMCRaw, *h2_xVypFocalMCRaw;
TH2F *h2_xpVyFocalMCRaw, *h2_xpVypFocalMCRaw, *h2_yVypFocalMCRaw;
TH2F *h2_yVxpTarMCRaw, *h2_yVypTarMCRaw, *h2_xpVypTarMCRaw;

TH1F *h_xFocalMCWgt, *h_xpFocalMCWgt, *h_yFocalMCWgt, *h_ypFocalMCWgt;
TH1F *h_yTarMCWgt, *h_xpTarMCWgt, *h_ypTarMCWgt, *h_deltaMCWgt;
TH2F *h2_xVxpFocalMCWgt, *h2_xVyFocalMCWgt, *h2_xVypFocalMCWgt;
TH2F *h2_xpVyFocalMCWgt, *h2_xpVypFocalMCWgt, *h2_yVypFocalMCWgt;
TH2F *h2_yVxpTarMCWgt, *h2_yVypTarMCWgt, *h2_xpVypTarMCWgt;  
// Data histos
TH1F *h_xFocalData, *h_xpFocalData, *h_yFocalData, *h_ypFocalData;
TH1F *h_yTarData, *h_xpTarData, *h_ypTarData, *h_deltaData;
TH2F *h2_xVxpFocalData, *h2_xVyFocalData, *h2_xVypFocalData;
TH2F *h2_xpVyFocalData, *h2_xpVypFocalData, *h2_yVypFocalData;
TH2F *h2_yVxpTarData, *h2_yVypTarData, *h2_xpVypTarData;

// Declare constants
static const Double_t protonMass = 0.938272;  // GeV

static const Double_t mcScaleFactor = 1.75554138E-02;
static const Double_t bsFactor = 14.0;

// Define functions to calculate kinematic variables
Double_t calc_nu(Double_t beamEnergy, Double_t specMom);
Double_t calc_q2(Double_t beamEnergy, Double_t specMom, Double_t specAngle);
Double_t calc_w2(Double_t beamEnergy, Double_t specMom, Double_t specAngle);

void mc_reweight() {

  mcRawFile  = new TFile("../../input/monte-carlo/kpp_shms_488.root");
  mcWgtFile  = new TFile("../../output/mc-ntuples/mc488.root");
  dataInFile = new TFile("../../input/shms-data/shms_replay_488_100000.root");
  
  // Obtain the ROOT trees
  mcRawTree = dynamic_cast <TTree*> (mcRawFile->Get("h1411"));
  mcWgtTree = dynamic_cast <TTree*> (mcWgtFile->Get("h9040"));
  dataTree  = dynamic_cast <TTree*> (dataInFile->Get("T"));
  // Open the data text files
  reconDataFile.open("../../input/recon-mc/kpp_shms_488_recon.in");
  inputDataFile.open("../../input/recon-mc/kpp_shms_488_input_root.dat");
  targetDataFile.open("../../input/target/kpp_shms_488_target.dat");
  // Obtain data from the input files
  reconDataFile >> target >> E0 >> P >> theta >> psFactors >> bcm1 >> bcm2 
		>> q1 >> q2 >> clt >> elt >> trackEff >> tofEff >> trigEff;
  inputDataFile >> numEvents >> phiRange >> thetaRange >> dPMax >> dPmin 
		>> chargeSymmFlag >> radCorrFlag;
  targetDataFile >> tarID >> tarDensity >> tarThickness >> tarAtomicMass 
		 >> tarAtomicNum >> tarArealDensity;
  // Close the data text files
  reconDataFile.close(); inputDataFile.close(), targetDataFile.close();
  // cout << target << "\t" << E0 << "\t" << P  << "\t" << theta  << endl;
  // cout << numEvents << "\t" << phiRange << "\t" << thetaRange << endl;
  // cout << tarID << "\t" << tarDensity << "\t" << tarThickness << endl;

  // Acquire the number of entries for each tree, and calculate scale factor
  nEntriesMCRaw = mcRawTree->GetEntries();
  nEntriesMCWgt = mcWgtTree->GetEntries();
  nEntriesData  = dataTree->GetEntries();

  // Acquire leafs of interest
  // MC leafs
  mcRawTree->SetBranchAddress("hsxfp",   &xFocalMCRaw);
  mcRawTree->SetBranchAddress("hsxpfp",  &xpFocalMCRaw);
  mcRawTree->SetBranchAddress("hsyfp",   &yFocalMCRaw);
  mcRawTree->SetBranchAddress("hsypfp",  &ypFocalMCRaw);
  mcRawTree->SetBranchAddress("hsytar",  &yTarMCRaw);
  mcRawTree->SetBranchAddress("hsxptar", &xpTarMCRaw);
  mcRawTree->SetBranchAddress("hsyptar", &ypTarMCRaw);
  mcRawTree->SetBranchAddress("hsdelta", &deltaMCRaw);

  mcWgtTree->SetBranchAddress("fail_id",  &failID);
  mcWgtTree->SetBranchAddress("xfoc",   &xFocalMCWgt);
  mcWgtTree->SetBranchAddress("dxdz",  &xpFocalMCWgt);
  mcWgtTree->SetBranchAddress("yfoc",   &yFocalMCWgt);
  mcWgtTree->SetBranchAddress("dydz",  &ypFocalMCWgt);
  mcWgtTree->SetBranchAddress("yrec",  &yTarMCWgt);
  mcWgtTree->SetBranchAddress("xprec", &xpTarMCWgt);
  mcWgtTree->SetBranchAddress("yprec", &ypTarMCWgt);
  mcWgtTree->SetBranchAddress("dppr", &deltaMCWgt);

  mcWgtTree->SetBranchAddress("born", &bornCS);
  mcWgtTree->SetBranchAddress("rci", &intRadCorr);
  // Data leafs
  // Cherenkov quantities
  dataTree->SetBranchAddress("P.hgcer.npeSum", &phgcerNpeSum);
  dataTree->SetBranchAddress("P.ngcer.npeSum", &pngcerNpeSum);
  dataTree->SetBranchAddress("P.aero.npeSum",  &paeroNpeSum);
  // Track quantities
  dataTree->SetBranchAddress("P.tr.n",           &numTracks);
  dataTree->SetBranchAddress("P.gtr.beta",       &goldTrackBeta);
  dataTree->SetBranchAddress("P.cal.etracknorm", &trackEnergyNorm);
  // Focal plane quantities (golden track)
  dataTree->SetBranchAddress("P.dc.x_fp",  &xFocalData);
  dataTree->SetBranchAddress("P.dc.xp_fp", &xpFocalData);
  dataTree->SetBranchAddress("P.dc.y_fp",  &yFocalData);
  dataTree->SetBranchAddress("P.dc.yp_fp", &ypFocalData);

  // Target quantities (golden track)
  dataTree->SetBranchAddress("P.gtr.y",  &yTarData);
  dataTree->SetBranchAddress("P.gtr.th", &xpTarData);
  dataTree->SetBranchAddress("P.gtr.ph", &ypTarData);
  dataTree->SetBranchAddress("P.gtr.dp", &deltaData);

  // Book histos
  // MC histos
  // 1D histos
  h_xFocalMCRaw  = new TH1F("h_xFocalMCRaw",  "Monte-Carlo Raw: X_{fp}; X_{fp} (cm); Number of Entries / 5 mm", 160, -40, 40);
  h_xpFocalMCRaw = new TH1F("h_xpFocalMCRaw", "Monte-Carlo Raw: X'_{fp}; X'_{fp}; Number of Entries / 0.001", 200, -0.1, 0.1);
  h_yFocalMCRaw  = new TH1F("h_yFocalMCRaw",  "Monte-Carlo Raw: Y_{fp}; Y_{fp} (cm); Number of Entries / 5 mm", 160, -40, 40);
  h_ypFocalMCRaw = new TH1F("h_ypFocalMCRaw", "Monte-Carlo Raw: Y'_{fp}; Y'_{fp}; Number of Entries / 0.001", 200, -0.1, 0.1);
  h_yTarMCRaw    = new TH1F("h_yTarMCRaw",    "Monte-Carlo Raw: Y_{tar}; Y_{tar} (cm); Number of Entries / 1 mm", 200, -10, 10);
  h_xpTarMCRaw   = new TH1F("h_xpTarMCRaw",   "Monte-Carlo Raw: X'_{tar}; X'_{tar}; Number of Entries / 0.001", 200, -0.1, 0.1);
  h_ypTarMCRaw   = new TH1F("h_ypTarMCRaw",   "Monte-Carlo Raw: Y'_{tar}; Y'_{tar}; Number of Entries / 0.001", 200, -0.1, 0.1);

  h_deltaMCRaw   = new TH1F("h_deltaMCRaw",   "Monte-Carlo Raw: #delta P; #delta P; Number of Entries", 160, -40, 40);
  // 2D histos
  h2_xVxpFocalMCRaw  = new TH2F("h2_xVxpFocalMCRaw",  "Monte-Carlo Raw: X_{fp} vs. X'_{fp}; X'_{fp} / 0.001; X_{fp} (cm) / 5 mm", 200, -0.1, 0.1, 160, -40, 40);
  h2_xVyFocalMCRaw   = new TH2F("h2_xVyFocalMCRaw",   "Monte-Carlo Raw: X_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X_{fp} (cm) / 5 mm", 160, -40, 40, 160, -40, 40);
  h2_xVypFocalMCRaw  = new TH2F("h2_xVypFocalMCRaw",  "Monte-Carlo Raw: X_{fp} vs. Y'_{fp}; Y'_{fp} / 0.001; X_{fp} (cm) / 5 mm", 200, -0.1, 0.1, 160, -40, 40);
  h2_xpVyFocalMCRaw  = new TH2F("h2_xpVyFocalMCRaw",  "Monte-Carlo Raw: X'_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X'_{fp} / 0.001", 160, -40, 40, 200, -0.1, 0.1);
  h2_xpVypFocalMCRaw = new TH2F("h2_xpVypFocalMCRaw", "Monte-Carlo Raw: X'_{fp} vs. Y'_{fp}; Y'_{fp} / 0.001; X'_{fp} / 0.001", 200, -0.1, 0.1, 200, -0.1, 0.1);
  h2_yVypFocalMCRaw  = new TH2F("h2_yVypFocalMCRaw",  "Monte-Carlo Raw: Y_{fp} vs. Y'_{fp}; Y'_{fp} / 0.001; Y_{fp} (cm) / 5 mm", 200, -0.1, 0.1, 160, -40, 40);
  h2_yVxpTarMCRaw    = new TH2F("h2_yVxpTarMCRaw",    "Monte-Carlo Raw: Y_{tar} vs. X'_{tar}; X'_{tar} / 0.001; Y_{tar} / 1 mm", 200, -0.1, 0.1, 200, -10, 10);
  h2_yVypTarMCRaw    = new TH2F("h2_yVypTarMCRaw",    "Monte-Carlo Raw: Y_{tar} vs. Y'_{tar}; Y'_{tar} / 0.001; Y_{tar} / 1 mm", 200, -0.1, 0.1, 200, -10, 10);
  h2_xpVypTarMCRaw   = new TH2F("h2_xpVypTarMCRaw",   "Monte-Carlo Raw: X'_{tar} vs. Y'_{tar}; Y'_{tar} / 0.001; X'_{tar} / 0.001", 200, -0.1, 0.1, 200, -0.1, 0.1);

  // MC histos
  // 1D histos
  h_xFocalMCWgt  = new TH1F("h_xFocalMCWgt",  "Monte-Carlo Weighted: X_{fp}; X_{fp} (cm); Number of Entries / 5 mm", 160, -40, 40);
  h_xpFocalMCWgt = new TH1F("h_xpFocalMCWgt", "Monte-Carlo Weighted: X'_{fp}; X'_{fp}; Number of Entries / 0.001", 200, -0.1, 0.1);
  h_yFocalMCWgt  = new TH1F("h_yFocalMCWgt",  "Monte-Carlo Weighted: Y_{fp}; Y_{fp} (cm); Number of Entries / 5 mm", 160, -40, 40);
  h_ypFocalMCWgt = new TH1F("h_ypFocalMCWgt", "Monte-Carlo Weighted: Y'_{fp}; Y'_{fp}; Number of Entries / 0.001", 200, -0.1, 0.1);
  h_yTarMCWgt    = new TH1F("h_yTarMCWgt",    "Monte-Carlo Weighted: Y_{tar}; Y_{tar} (cm); Number of Entries / 1 mm", 200, -10, 10);
  h_xpTarMCWgt   = new TH1F("h_xpTarMCWgt",   "Monte-Carlo Weighted: X'_{tar}; X'_{tar}; Number of Entries / 0.001", 200, -0.1, 0.1);
  h_ypTarMCWgt   = new TH1F("h_ypTarMCWgt",   "Monte-Carlo Weighted: Y'_{tar}; Y'_{tar}; Number of Entries / 0.001", 200, -0.1, 0.1);
  h_deltaMCWgt   = new TH1F("h_deltaMCWgt",   "Monte-Carlo Weighted: #delta P; #delta P; Number of Entries", 160, -40, 40);
  // 2D histos
  h2_xVxpFocalMCWgt  = new TH2F("h2_xVxpFocalMCWgt",  "Monte-Carlo Weighted: X_{fp} vs. X'_{fp}; X'_{fp} / 0.001; X_{fp} (cm) / 5 mm", 200, -0.1, 0.1, 160, -40, 40);
  h2_xVyFocalMCWgt   = new TH2F("h2_xVyFocalMCWgt",   "Monte-Carlo Weighted: X_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X_{fp} (cm) / 5 mm", 160, -40, 40, 160, -40, 40);
  h2_xVypFocalMCWgt  = new TH2F("h2_xVypFocalMCWgt",  "Monte-Carlo Weighted: X_{fp} vs. Y'_{fp}; Y'_{fp} / 0.001; X_{fp} (cm) / 5 mm", 200, -0.1, 0.1, 160, -40, 40);
  h2_xpVyFocalMCWgt  = new TH2F("h2_xpVyFocalMCWgt",  "Monte-Carlo Weighted: X'_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X'_{fp} / 0.001", 160, -40, 40, 200, -0.1, 0.1);
  h2_xpVypFocalMCWgt = new TH2F("h2_xpVypFocalMCWgt", "Monte-Carlo Weighted: X'_{fp} vs. Y'_{fp}; Y'_{fp} / 0.001; X'_{fp} / 0.001", 200, -0.1, 0.1, 200, -0.1, 0.1);
  h2_yVypFocalMCWgt  = new TH2F("h2_yVypFocalMCWgt",  "Monte-Carlo Weighted: Y_{fp} vs. Y'_{fp}; Y'_{fp} / 0.001; Y_{fp} (cm) / 5 mm", 200, -0.1, 0.1, 160, -40, 40);
  h2_yVxpTarMCWgt    = new TH2F("h2_yVxpTarMCWgt",    "Monte-Carlo Weighted: Y_{tar} vs. X'_{tar}; X'_{tar} / 0.001; Y_{tar} / 1 mm", 200, -0.1, 0.1, 200, -10, 10);
  h2_yVypTarMCWgt    = new TH2F("h2_yVypTarMCWgt",    "Monte-Carlo Weighted: Y_{tar} vs. Y'_{tar}; Y'_{tar} / 0.001; Y_{tar} / 1 mm", 200, -0.1, 0.1, 200, -10, 10);
  h2_xpVypTarMCWgt   = new TH2F("h2_xpVypTarMCWgt",   "Monte-Carlo Weighted: X'_{tar} vs. Y'_{tar}; Y'_{tar} / 0.001; X'_{tar} / 0.001", 200, -0.1, 0.1, 200, -0.1, 0.1);

  // Add in target permutations 

  // Data histos
  // 1D histos
  h_xFocalData  = new TH1F("h_xFocalData",  "Data: X_{fp}; X_{fp} (cm); Number of Entries / 5 mm", 160, -40, 40);
  h_xpFocalData = new TH1F("h_xpFocalData", "Data: X'_{fp}; X'_{fp} (cm); Number of Entries / 0.001", 200, -0.1, 0.1);
  h_yFocalData  = new TH1F("h_yFocalData",  "Data: Y_{fp}; Y_{fp} (cm); Number of Entries / 5 mm", 160, -40, 40);
  h_ypFocalData = new TH1F("h_ypFocalData", "Data: Y'_{fp}; Y'_{fp} (cm); Number of Entries / 0.001", 200, -0.1, 0.1);
  h_yTarData    = new TH1F("h_yTarData",    "Data: Y_{tar}; Y_{tar} (cm); Number of Entries / 1 mm", 200, -10, 10);
  h_xpTarData   = new TH1F("h_xpTarData",   "Data: X'_{tar}; X'_{tar} (cm); Number of Entries / 0.001", 200, -0.1, 0.1);
  h_ypTarData   = new TH1F("h_ypTarData",   "Data: Y'_{tar}; Y'_{tar} (cm); Number of Entries / 0.001", 200, -0.1, 0.1);
  h_deltaData   = new TH1F("h_deltaData",   "Data: #delta P; #delta P; Number of Entries", 160, -40, 40);
  // 2D histos
  h2_xVxpFocalData  = new TH2F("h2_xVxpFocalData",  "Data: X_{fp} vs. X'_{fp}; X'_{fp} / 0.001; X_{fp} (cm) / 5 mm", 200, -0.1, 0.1, 160, -40, 40);
  h2_xVyFocalData   = new TH2F("h2_xVyFocalData",   "Data: X_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X_{fp} (cm) / 5 mm", 160, -40, 40, 160, -40, 40);
  h2_xVypFocalData  = new TH2F("h2_xVypFocalData",  "Data: X_{fp} vs. Y'_{fp}; Y'_{fp} / 0.001; X_{fp} (cm) / 5 mm", 200, -0.1, 0.1, 160, -40, 40);
  h2_xpVyFocalData  = new TH2F("h2_xpVyFocalData",  "Data: X'_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X'_{fp} / 0.001", 160, -40, 40, 200, -0.1, 0.1);
  h2_xpVypFocalData = new TH2F("h2_xpVypFocalData", "Data: X'_{fp} vs. Y'_{fp}; Y'_{fp} / 0.001; X'_{fp} / 0.001", 200, -0.1, 0.1, 200, -0.1, 0.1);
  h2_yVypFocalData  = new TH2F("h2_yVypFocalData",  "Data: Y_{fp} vs. Y'_{fp}; Y'_{fp} / 0.001; Y_{fp} (cm) / 5 mm", 200, -0.1, 0.1, 160, -40, 40);
  h2_yVxpTarData    = new TH2F("h2_yVxpTarData",    "Data: Y_{tar} vs. X'_{tar}; X'_{tar} / 0.001; Y_{tar} / 1 mm", 200, -0.1, 0.1, 200, -10, 10);
  h2_yVypTarData    = new TH2F("h2_yVypTarData",    "Data: Y_{tar} vs. Y'_{tar}; Y'_{tar} / 0.001; Y_{tar} / 1 mm", 200, -0.1, 0.1, 200, -10, 10);
  h2_xpVypTarData   = new TH2F("h2_xpVypTarData",   "Data: X'_{tar} vs. Y'_{tar}; Y'_{tar} / 0.001; X'_{tar} / 0.001", 200, -0.1, 0.1, 200, -0.1, 0.1);

  UInt_t mcRawEventCntr = 0;
  // Loop over raw MC data
  for (UInt_t imcRaw = 0; imcRaw < nEntriesMCRaw; imcRaw++) {
    // Obtain the data entry
    mcRawTree->GetEntry(imcRaw);
    if (mcRawEventCntr % 10000 == 0 && mcRawEventCntr != 0) 
      cout << mcRawEventCntr << " Raw Monte-Carlo events have been processed..." << endl;
    mcRawEventCntr++;
    // Fill histos
    // 1D histos
    h_xFocalMCRaw->Fill(xFocalMCRaw);
    h_xpFocalMCRaw->Fill(xpFocalMCRaw);
    h_yFocalMCRaw->Fill(yFocalMCRaw);
    h_ypFocalMCRaw->Fill(ypFocalMCRaw);
    h_yTarMCRaw->Fill(yTarMCRaw);
    h_xpTarMCRaw->Fill(xpTarMCRaw);
    h_ypTarMCRaw->Fill(ypTarMCRaw);
    h_deltaMCRaw->Fill(deltaMCRaw);
    // 2D histos
    h2_xVxpFocalMCRaw->Fill(xpFocalMCRaw, xFocalMCRaw);
    h2_xVyFocalMCRaw->Fill(yFocalMCRaw, xFocalMCRaw);
    h2_xVypFocalMCRaw->Fill(ypFocalMCRaw, xFocalMCRaw);
    h2_xpVyFocalMCRaw->Fill(yFocalMCRaw, xpFocalMCRaw);
    h2_xpVypFocalMCRaw->Fill(ypFocalMCRaw, xpFocalMCRaw);
    h2_yVypFocalMCRaw->Fill(ypFocalMCRaw, yFocalMCRaw);
    h2_yVxpTarMCRaw->Fill(xpTarMCRaw, yTarMCRaw);
    h2_yVypTarMCRaw->Fill(ypTarMCRaw, yTarMCRaw);
    h2_xpVypTarMCRaw->Fill(ypTarMCRaw, xpTarMCRaw);

  }  // Raw MC loop

  UInt_t mcWgtEventCntr = 0;
  // Loop over raw MC data
  for (UInt_t imcWgt = 0; imcWgt < nEntriesMCWgt; imcWgt++) {
    // Obtain the data entry
    mcWgtTree->GetEntry(imcWgt);
    if (mcWgtEventCntr % 10000 == 0 && mcWgtEventCntr != 0) 
      cout << mcWgtEventCntr << " Weighted Monte-Carlo events have been processed..." << endl;
    mcWgtEventCntr++;
    Bool_t failIDCut = failID == 0.0;
    if (!failIDCut) continue;

    Bool_t deltaPCut = TMath::Abs(deltaMCWgt - 5.0) < 20.0;
    if (!deltaPCut) continue;

    scaleFactor = (bornCS * mcScaleFactor) / (bsFactor * intRadCorr);
    // Fill histos
    // 1D histos
    h_xFocalMCWgt->Fill(xFocalMCWgt, scaleFactor);
    h_xpFocalMCWgt->Fill(xpFocalMCWgt, scaleFactor);
    h_yFocalMCWgt->Fill(yFocalMCWgt, scaleFactor);
    h_ypFocalMCWgt->Fill(ypFocalMCWgt, scaleFactor);
    h_yTarMCWgt->Fill(yTarMCWgt, scaleFactor);
    h_xpTarMCWgt->Fill(xpTarMCWgt, scaleFactor);
    h_ypTarMCWgt->Fill(ypTarMCWgt, scaleFactor);
    h_deltaMCWgt->Fill(deltaMCWgt, scaleFactor);
    // 2D histos
    h2_xVxpFocalMCWgt->Fill(xpFocalMCWgt, xFocalMCWgt, scaleFactor);
    h2_xVyFocalMCWgt->Fill(yFocalMCWgt, xFocalMCWgt, scaleFactor);
    h2_xVypFocalMCWgt->Fill(ypFocalMCWgt, xFocalMCWgt, scaleFactor);
    h2_xpVyFocalMCWgt->Fill(yFocalMCWgt, xpFocalMCWgt, scaleFactor);
    h2_xpVypFocalMCWgt->Fill(ypFocalMCWgt, xpFocalMCWgt, scaleFactor);
    h2_yVypFocalMCWgt->Fill(ypFocalMCWgt, yFocalMCWgt, scaleFactor);
    h2_yVxpTarMCWgt->Fill(xpTarMCWgt, yTarMCWgt, scaleFactor);
    h2_yVypTarMCWgt->Fill(ypTarMCWgt, yTarMCWgt, scaleFactor);
    h2_xpVypTarMCWgt->Fill(ypTarMCWgt, xpTarMCWgt, scaleFactor);

  }  // Wgt MC loop

  UInt_t dataEventCntr = 0;
  // Loop over data
  for (UInt_t idata = 0; idata < nEntriesData; idata++) {
    // Obtain the data entry
    dataTree->GetEntry(idata);
    if (dataEventCntr % 10000 == 0 && dataEventCntr != 0) 
      cout << dataEventCntr << " Data events have been processed..." << endl;
    dataEventCntr++;

    // Number of photo-electron cuts
    Bool_t phgcerNpeSumCut = phgcerNpeSum >= 1.0;
    Bool_t pngcerNpeSumCut = pngcerNpeSum >= 2.0;
    Bool_t paeroNpeSumCut  = paeroNpeSum  >= 2.0;
    Bool_t npeCuts         = phgcerNpeSumCut && pngcerNpeSumCut && paeroNpeSum;
    //if (!npeCuts) continue;
    if (!phgcerNpeSumCut) continue;
    
    // Tracking cuts
    Bool_t numTracksCut       = numTracks > 0.0;
    Bool_t goldTrackBetaCut   = TMath::Abs(goldTrackBeta - 1.0) < 0.2;
    //Bool_t trackEnergyNormCut = TMath::Abs(trackEnergyNorm - 1.1) < 0.3;
    Bool_t trackCuts          = numTracksCut && goldTrackBetaCut;
    //Bool_t trackCuts          = numTracksCut && goldTrackBetaCut && trackEnergyNormCut;
    if (!trackCuts) continue;

    Bool_t deltaPCut = TMath::Abs(deltaData - 5.0) < 20.0;
    if (!deltaPCut) continue;

    // Fill histos
    // 1D histos
    h_xFocalData->Fill(xFocalData);
    h_xpFocalData->Fill(xpFocalData);
    h_yFocalData->Fill(yFocalData);
    h_ypFocalData->Fill(ypFocalData);
    h_yTarData->Fill(yTarData);
    h_xpTarData->Fill(xpTarData);
    h_ypTarData->Fill(ypTarData);
    h_deltaData->Fill(deltaData);
    // 2D histos
    h2_xVxpFocalData->Fill(xpFocalData, xFocalData);
    h2_xVyFocalData->Fill(yFocalData, xFocalData);
    h2_xVypFocalData->Fill(ypFocalData, xFocalData);
    h2_xpVyFocalData->Fill(yFocalData, xpFocalData);
    h2_xpVypFocalData->Fill(ypFocalData, xpFocalData);
    h2_yVypFocalData->Fill(ypFocalData, yFocalData);
    h2_yVxpTarData->Fill(xpTarData, yTarData);
    h2_yVypTarData->Fill(ypTarData, yTarData);
    h2_xpVypTarData->Fill(ypTarData, xpTarData);
  }  // Data loop

}

// User specific funtions
Double_t calc_nu(Double_t beamEnergy, Double_t specMom) {
  return beamEnergy - specMom;
}

Double_t calc_q2(Double_t beamEnergy, Double_t specMom, Double_t specAngle) {
  return 2.0*beamEnergy*specMom*(1.0 - TMath::Cos(specAngle*TMath::DegToRad()));
}

Double_t calc_w2(Double_t beamEnergy, Double_t specMom, Double_t specAngle) {
  return protonMass*protonMass + 2.0*protonMass*(beamEnergy - specMom) - 2.0*beamEnergy*specMom*(1.0 - TMath::Cos(specAngle*TMath::DegToRad()));
}
