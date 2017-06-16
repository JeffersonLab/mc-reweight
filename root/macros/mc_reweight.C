// Macro to create histograms of MC events and Data
// Author: Eric Pooser, pooser@jlab.org

// Declare data files and directories
TFile *mcRawFile, *mcWgtFile, *dataFile, *compFile;
TDirectory *mcRawDir, *mcWgtDir, *dataDir;
// Declare trees
TTree *mcRawTree, *mcWgtTree, *dataTree;
// Declare text data files
ifstream reconDataFile, inputDataFile, targetDataFile;
// Declare variables
// Input text file data
Float_t target, E0, P, theta, psFactors, bcm1, bcm2, q1, q2, clt, elt, trackEff, tofEff, trigEff;
Float_t numEvents, phiRange, thetaRange, dPMax, dPmin, chargeSymmFlag, radCorrFlag;
Float_t tarID, tarDensity, tarThickness, tarAtomicMass, tarAtomicNum, tarArealDensity;
// Number of entries
Int_t    nEntriesMCRaw, nEntriesMCWgt, nEntriesData;
// Leaf variables
Double_t phgcerNpeSum, pngcerNpeSum, paeroNpeSum; 
Double_t gtrBetaData, gtrMomData, trackEnergyNorm;
// Input MC variables
Float_t  xFocalMCRaw, xpFocalMCRaw, yFocalMCRaw, ypFocalMCRaw;
Float_t  yTarMCRaw, xpTarMCRaw, ypTarMCRaw, deltaMCRaw;
// Weighted MC variables
Float_t  xFocalMCWgt, xpFocalMCWgt, yFocalMCWgt, ypFocalMCWgt;
Float_t  yTarMCWgt, xpTarMCWgt, ypTarMCWgt, deltaMCWgt;
Float_t  thetaMCWgt, q2MCWgt, w2MCWgt;
Float_t  failID, bornCS, intRadCorr;
Double_t scaleFactor;
// Data variables
Double_t xFocalData, xpFocalData, yFocalData, ypFocalData, gtrOkData;
Double_t yTarData, xpTarData, ypTarData, deltaData;
Double_t thetaData, q2Data, w2Data;
// Declare histos
// Input MC histos
TH1F *h_xFocalMCRaw, *h_xpFocalMCRaw, *h_yFocalMCRaw, *h_ypFocalMCRaw;
TH1F *h_yTarMCRaw, *h_xpTarMCRaw, *h_ypTarMCRaw, *h_deltaMCRaw;
TH2F *h2_xVxpFocalMCRaw, *h2_xVyFocalMCRaw, *h2_xVypFocalMCRaw;
TH2F *h2_xpVyFocalMCRaw, *h2_xpVypFocalMCRaw, *h2_yVypFocalMCRaw;
TH2F *h2_yVxpTarMCRaw, *h2_yVypTarMCRaw, *h2_xpVypTarMCRaw;
// Weighted MC histos
TH1F *h_xFocalMCWgt, *h_xpFocalMCWgt, *h_yFocalMCWgt, *h_ypFocalMCWgt;
TH1F *h_yTarMCWgt, *h_xpTarMCWgt, *h_ypTarMCWgt, *h_deltaMCWgt;
TH1F *h_thetaMCWgt, *h_q2MCWgt, *h_w2MCWgt;
TH2F *h2_xVxpFocalMCWgt, *h2_xVyFocalMCWgt, *h2_xVypFocalMCWgt;
TH2F *h2_xpVyFocalMCWgt, *h2_xpVypFocalMCWgt, *h2_yVypFocalMCWgt;
TH2F *h2_yVxpTarMCWgt, *h2_yVypTarMCWgt, *h2_xpVypTarMCWgt;  
// Data histos
TH1F *h_xFocalData, *h_xpFocalData, *h_yFocalData, *h_ypFocalData;
TH1F *h_yTarData, *h_xpTarData, *h_ypTarData, *h_deltaData;
TH1F *h_thetaData, *h_q2Data, *h_w2Data;
TH2F *h2_xVxpFocalData, *h2_xVyFocalData, *h2_xVypFocalData;
TH2F *h2_xpVyFocalData, *h2_xpVypFocalData, *h2_yVypFocalData;
TH2F *h2_yVxpTarData, *h2_yVypTarData, *h2_xpVypTarData;
// Declare constants
static const Double_t protonMass    = 0.938272;  // GeV
static const Double_t mcScaleFactor = 1.75554138E-02;  // Output from recon_mc.f
static const Double_t rad2mrad      = 1000.0;
static const Double_t fudgeFactor   = 14.0;  // Arbitrary scaling factor
// Define functions to calculate kinematic variables
Double_t calc_q2(Double_t beamEnergy, Double_t scatMom, Double_t scatAngle);
Double_t calc_w2(Double_t beamEnergy, Double_t scatMom, Double_t scatAngle);

void mc_reweight() {
  TH1::SetDefaultSumw2(false);
  // Open the ROOT files
  mcRawFile = new TFile("../../input/monte-carlo/kpp_shms_488.root");
  mcWgtFile = new TFile("../../output/mc-ntuples/mc488.root");
  dataFile  = new TFile("../../input/shms-data/shms_replay_488_100000.root");
  compFile  = new TFile("../output/run_488_comp.root", "RECREATE");
  // Obtain the ROOT trees
  mcRawTree = dynamic_cast <TTree*> (mcRawFile->Get("h1411"));
  mcWgtTree = dynamic_cast <TTree*> (mcWgtFile->Get("h9040"));
  dataTree  = dynamic_cast <TTree*> (dataFile->Get("T"));
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
  // Input MC leafs
  mcRawTree->SetBranchAddress("hsxfp",   &xFocalMCRaw);
  mcRawTree->SetBranchAddress("hsxpfp",  &xpFocalMCRaw);
  mcRawTree->SetBranchAddress("hsyfp",   &yFocalMCRaw);
  mcRawTree->SetBranchAddress("hsypfp",  &ypFocalMCRaw);
  mcRawTree->SetBranchAddress("hsytar",  &yTarMCRaw);
  mcRawTree->SetBranchAddress("hsxptar", &xpTarMCRaw);
  mcRawTree->SetBranchAddress("hsyptar", &ypTarMCRaw);
  mcRawTree->SetBranchAddress("hsdelta", &deltaMCRaw);
  // Weighted MC leafs
  mcWgtTree->SetBranchAddress("fail_id", &failID);
  mcWgtTree->SetBranchAddress("xfoc",    &xFocalMCWgt);
  mcWgtTree->SetBranchAddress("dxdz",    &xpFocalMCWgt);
  mcWgtTree->SetBranchAddress("yfoc",    &yFocalMCWgt);
  mcWgtTree->SetBranchAddress("dydz",    &ypFocalMCWgt);
  mcWgtTree->SetBranchAddress("yrec",    &yTarMCWgt);
  mcWgtTree->SetBranchAddress("xprec",   &xpTarMCWgt);
  mcWgtTree->SetBranchAddress("yprec",   &ypTarMCWgt);
  mcWgtTree->SetBranchAddress("dppr",    &deltaMCWgt);
  mcWgtTree->SetBranchAddress("hstheta", &thetaMCWgt);
  mcWgtTree->SetBranchAddress("q2",      &q2MCWgt);
  mcWgtTree->SetBranchAddress("w2",      &w2MCWgt);
  mcWgtTree->SetBranchAddress("born",    &bornCS);
  mcWgtTree->SetBranchAddress("rci",     &intRadCorr);
  // Data leafs
  // Cherenkov quantities
  dataTree->SetBranchAddress("P.hgcer.npeSum", &phgcerNpeSum);
  dataTree->SetBranchAddress("P.ngcer.npeSum", &pngcerNpeSum);
  dataTree->SetBranchAddress("P.aero.npeSum",  &paeroNpeSum);
  // Track quantities
  dataTree->SetBranchAddress("P.gtr.beta",       &gtrBetaData);
  dataTree->SetBranchAddress("P.cal.etracknorm", &trackEnergyNorm);
  // Focal plane quantities
  dataTree->SetBranchAddress("P.dc.x_fp",  &xFocalData);
  dataTree->SetBranchAddress("P.dc.xp_fp", &xpFocalData);
  dataTree->SetBranchAddress("P.dc.y_fp",  &yFocalData);
  dataTree->SetBranchAddress("P.dc.yp_fp", &ypFocalData);
  // Target quantities (golden track)
  dataTree->SetBranchAddress("P.gtr.p",  &gtrMomData);
  dataTree->SetBranchAddress("P.gtr.y",  &yTarData);
  dataTree->SetBranchAddress("P.gtr.th", &xpTarData);
  dataTree->SetBranchAddress("P.gtr.ph", &ypTarData);
  dataTree->SetBranchAddress("P.gtr.dp", &deltaData);
  dataTree->SetBranchAddress("P.gtr.ok", &gtrOkData);
  
  // Create input MC directory and descend into it
  mcRawDir = dynamic_cast <TDirectory*> (compFile->Get("mcRawDir"));
  if(!mcRawDir) {mcRawDir = compFile->mkdir("mcRawDir"); mcRawDir->cd();}
  else compFile->cd("mcRawDir");
  // Book input 1D MC histos
  h_xFocalMCRaw  = new TH1F("h_xFocalMCRaw",  "Input Monte-Carlo: X_{fp}; X_{fp} (cm); Number of Entries / 5 mm",   160, -40, 40);
  h_xpFocalMCRaw = new TH1F("h_xpFocalMCRaw", "Input Monte-Carlo: X'_{fp}; X'_{fp}; Number of Entries / 2 mrad",    60, -60.0, 60.0);
  h_yFocalMCRaw  = new TH1F("h_yFocalMCRaw",  "Input Monte-Carlo: Y_{fp}; Y_{fp} (cm); Number of Entries / 5 mm",   160, -40, 40);
  h_ypFocalMCRaw = new TH1F("h_ypFocalMCRaw", "Input Monte-Carlo: Y'_{fp}; Y'_{fp}; Number of Entries / 2 mrad",    60, -60.0, 60.0);
  h_yTarMCRaw    = new TH1F("h_yTarMCRaw",    "Input Monte-Carlo: Y_{tar}; Y_{tar} (cm); Number of Entries / 1 mm", 100, -5, 5);
  h_xpTarMCRaw   = new TH1F("h_xpTarMCRaw",   "Input Monte-Carlo: X'_{tar}; X'_{tar}; Number of Entries / 2 mrad",  60, -60.0, 60.0);
  h_ypTarMCRaw   = new TH1F("h_ypTarMCRaw",   "Input Monte-Carlo: Y'_{tar}; Y'_{tar}; Number of Entries / 2 mrad",  60, -60.0, 60.0);
  h_deltaMCRaw   = new TH1F("h_deltaMCRaw",   "Input Monte-Carlo: #delta; #delta; Number of Entries",             80, -40, 40);
  // Book input 2D MC histos
  h2_xVxpFocalMCRaw  = new TH2F("h2_xVxpFocalMCRaw",  "Input Monte-Carlo: X_{fp} vs. X'_{fp}; X'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm",    100, -100.0, 100.0, 160, -40, 40);
  h2_xVyFocalMCRaw   = new TH2F("h2_xVyFocalMCRaw",   "Input Monte-Carlo: X_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X_{fp} (cm) / 5 mm",   160, -40, 40, 160, -40, 40);
  h2_xVypFocalMCRaw  = new TH2F("h2_xVypFocalMCRaw",  "Input Monte-Carlo: X_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm",    60, -60.0, 60.0, 160, -40, 40);
  h2_xpVyFocalMCRaw  = new TH2F("h2_xpVyFocalMCRaw",  "Input Monte-Carlo: X'_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X'_{fp} / 2 mrad",    160, -40, 40, 100, -100.0, 100.0);
  h2_xpVypFocalMCRaw = new TH2F("h2_xpVypFocalMCRaw", "Input Monte-Carlo: X'_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X'_{fp} / 2 mrad",     60, -60.0, 60.0, 100, -100.0, 100.0);
  h2_yVypFocalMCRaw  = new TH2F("h2_yVypFocalMCRaw",  "Input Monte-Carlo: Y_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; Y_{fp} (cm) / 5 mm",    60, -60.0, 60.0, 160, -40, 40);
  h2_yVxpTarMCRaw    = new TH2F("h2_yVxpTarMCRaw",    "Input Monte-Carlo: Y_{tar} vs. X'_{tar}; X'_{tar} / 2 mrad; Y_{tar} / 1 mm",     200, -100.0, 100.0, 100, -5, 5);
  h2_yVypTarMCRaw    = new TH2F("h2_yVypTarMCRaw",    "Input Monte-Carlo: Y_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; Y_{tar} / 1 mm",     200, -100.0, 100.0, 100, -5, 5);
  h2_xpVypTarMCRaw   = new TH2F("h2_xpVypTarMCRaw",   "Input Monte-Carlo: X'_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; X'_{tar} / 2 mrad", 60, -60.0, 60.0, 100, -100.0, 100.0);
  compFile->cd("../");

  // Create weighted MC directory and descend into it
  mcWgtDir = dynamic_cast <TDirectory*> (compFile->Get("mcWgtDir"));
  if(!mcWgtDir) {mcWgtDir = compFile->mkdir("mcWgtDir"); mcWgtDir->cd();}
  else compFile->cd("mcWgtDir");
  // Book weighted 1D MC histos
  h_xFocalMCWgt  = new TH1F("h_xFocalMCWgt",  "Weighted Monte-Carlo: X_{fp}; X_{fp} (cm); Number of Entries / 5 mm",             160, -40, 40);    h_xFocalMCWgt->SetBit(TH1::kIsNotW);
  h_xpFocalMCWgt = new TH1F("h_xpFocalMCWgt", "Weighted Monte-Carlo: X'_{fp}; X'_{fp}; Number of Entries / 2 mrad",              60, -60.0, 60.0); h_xpFocalMCWgt->SetBit(TH1::kIsNotW); 
  h_yFocalMCWgt  = new TH1F("h_yFocalMCWgt",  "Weighted Monte-Carlo: Y_{fp}; Y_{fp} (cm); Number of Entries / 5 mm",             160, -40, 40);    h_yFocalMCWgt->SetBit(TH1::kIsNotW);
  h_ypFocalMCWgt = new TH1F("h_ypFocalMCWgt", "Weighted Monte-Carlo: Y'_{fp}; Y'_{fp}; Number of Entries / 2 mrad",              60, -60.0, 60.0); h_ypFocalMCWgt->SetBit(TH1::kIsNotW);
  h_yTarMCWgt    = new TH1F("h_yTarMCWgt",    "Weighted Monte-Carlo: Y_{tar}; Y_{tar} (cm); Number of Entries / 1 mm",           100, -5, 5);      h_yTarMCWgt->SetBit(TH1::kIsNotW);
  h_xpTarMCWgt   = new TH1F("h_xpTarMCWgt",   "Weighted Monte-Carlo: X'_{tar}; X'_{tar}; Number of Entries / 2 mrad",            60, -60.0, 60.0); h_xpTarMCWgt->SetBit(TH1::kIsNotW); 
  h_ypTarMCWgt   = new TH1F("h_ypTarMCWgt",   "Weighted Monte-Carlo: Y'_{tar}; Y'_{tar}; Number of Entries / 2 mrad",            60, -60.0, 60.0); h_ypTarMCWgt->SetBit(TH1::kIsNotW);
  h_deltaMCWgt   = new TH1F("h_deltaMCWgt",   "Weighted Monte-Carlo: #delta; #delta; Number of Entries",                       80, -40, 40);     h_deltaMCWgt->SetBit(TH1::kIsNotW);
  h_thetaMCWgt   = new TH1F("h_thetaMCWgt",   "Weighted Monte-Carlo: #theta; #theta; Number of Entries / 0.01 deg",              100, 10.0, 20.0); h_thetaMCWgt->SetBit(TH1::kIsNotW);
  h_q2MCWgt      = new TH1F("h_q2MCWgt",      "Weighted Monte-Carlo: Q^{2}; Q^{2} (GeV^{2}); Number of Entries / 0.025 GeV^{2}", 120, 0.0, 3.0);   h_q2MCWgt->SetBit(TH1::kIsNotW);
  h_w2MCWgt      = new TH1F("h_w2MCWgt",      "Weighted Monte-Carlo: W^{2}; W^{2} (GeV^{2}); Number of Entries / 0.050 GeV^{2}", 180, 0.0, 9.0);   h_w2MCWgt->SetBit(TH1::kIsNotW);
  // Weighted 2D MC histos
  h2_xVxpFocalMCWgt  = new TH2F("h2_xVxpFocalMCWgt",  "Weighted Monte-Carlo: X_{fp} vs. X'_{fp}; X'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm",    100, -100.0, 100.0, 160, -40, 40);
  h2_xVyFocalMCWgt   = new TH2F("h2_xVyFocalMCWgt",   "Weighted Monte-Carlo: X_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X_{fp} (cm) / 5 mm",   160, -40, 40, 160, -40, 40);
  h2_xVypFocalMCWgt  = new TH2F("h2_xVypFocalMCWgt",  "Weighted Monte-Carlo: X_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm",    60, -60.0, 60.0, 160, -40, 40);
  h2_xpVyFocalMCWgt  = new TH2F("h2_xpVyFocalMCWgt",  "Weighted Monte-Carlo: X'_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X'_{fp} / 2 mrad",    160, -40, 40, 100, -100.0, 100.0);
  h2_xpVypFocalMCWgt = new TH2F("h2_xpVypFocalMCWgt", "Weighted Monte-Carlo: X'_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X'_{fp} / 2 mrad",     60, -60.0, 60.0, 100, -100.0, 100.0);
  h2_yVypFocalMCWgt  = new TH2F("h2_yVypFocalMCWgt",  "Weighted Monte-Carlo: Y_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; Y_{fp} (cm) / 5 mm",    60, -60.0, 60.0, 160, -40, 40);
  h2_yVxpTarMCWgt    = new TH2F("h2_yVxpTarMCWgt",    "Weighted Monte-Carlo: Y_{tar} vs. X'_{tar}; X'_{tar} / 2 mrad; Y_{tar} / 1 mm",     200, -100.0, 100.0, 100, -5, 5);
  h2_yVypTarMCWgt    = new TH2F("h2_yVypTarMCWgt",    "Weighted Monte-Carlo: Y_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; Y_{tar} / 1 mm",     200, -100.0, 100.0, 100, -5, 5);
  h2_xpVypTarMCWgt   = new TH2F("h2_xpVypTarMCWgt",   "Weighted Monte-Carlo: X'_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; X'_{tar} / 2 mrad", 60, -60.0, 60.0, 100, -100.0, 100.0);
  compFile->cd("../");

  // Create data directory and descend into it
  dataDir = dynamic_cast <TDirectory*> (compFile->Get("dataDir"));
  if(!dataDir) {dataDir = compFile->mkdir("dataDir"); dataDir->cd();}
  else compFile->cd("dataDir");
  // Data 1D histos
  h_xFocalData  = new TH1F("h_xFocalData",  "Data: X_{fp}; X_{fp} (cm); Number of Entries / 5 mm",             160, -40, 40);
  h_xpFocalData = new TH1F("h_xpFocalData", "Data: X'_{fp}; X'_{fp}; Number of Entries / 2 mrad",              60, -60.0, 60.0);
  h_yFocalData  = new TH1F("h_yFocalData",  "Data: Y_{fp}; Y_{fp} (cm); Number of Entries / 5 mm",             160, -40, 40);
  h_ypFocalData = new TH1F("h_ypFocalData", "Data: Y'_{fp}; Y'_{fp}; Number of Entries / 2 mrad",              60, -60.0, 60.0);
  h_yTarData    = new TH1F("h_yTarData",    "Data: Y_{tar}; Y_{tar} (cm); Number of Entries / 1 mm",           100, -5, 5);
  h_xpTarData   = new TH1F("h_xpTarData",   "Data: X'_{tar}; X'_{tar}; Number of Entries / 2 mrad",            60, -60.0, 60.0);
  h_ypTarData   = new TH1F("h_ypTarData",   "Data: Y'_{tar}; Y'_{tar}; Number of Entries / 2 mrad",            60, -60.0, 60.0);
  h_deltaData   = new TH1F("h_deltaData",   "Data: #delta; #delta; Number of Entries",                       80, -40, 40);
  h_thetaData   = new TH1F("h_thetaData",   "Data: #theta; #theta; Number of Entries / 0.01 deg",              100, 10.0, 20.0);
  h_q2Data      = new TH1F("h_q2Data",      "Data: Q^{2}; Q^{2} (GeV^{2}); Number of Entries / 0.025 GeV^{2}", 120, 0.0, 3.0);
  h_w2Data      = new TH1F("h_w2Data",      "Data: W^{2}; W^{2} (GeV^{2}); Number of Entries / 0.050 GeV^{2}", 180, 0.0, 9.0);
  // Data 2D histos
  h2_xVxpFocalData  = new TH2F("h2_xVxpFocalData",  "Data: X_{fp} vs. X'_{fp}; X'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm",    100, -100.0, 100.0, 160, -40, 40);
  h2_xVyFocalData   = new TH2F("h2_xVyFocalData",   "Data: X_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X_{fp} (cm) / 5 mm",   160, -40, 40, 160, -40, 40);
  h2_xVypFocalData  = new TH2F("h2_xVypFocalData",  "Data: X_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm",    60, -60.0, 60.0, 160, -40, 40);
  h2_xpVyFocalData  = new TH2F("h2_xpVyFocalData",  "Data: X'_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X'_{fp} / 2 mrad",    160, -40, 40, 100, -100.0, 100.0);
  h2_xpVypFocalData = new TH2F("h2_xpVypFocalData", "Data: X'_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X'_{fp} / 2 mrad",     60, -60.0, 60.0, 100, -100.0, 100.0);
  h2_yVypFocalData  = new TH2F("h2_yVypFocalData",  "Data: Y_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; Y_{fp} (cm) / 5 mm",    60, -60.0, 60.0, 160, -40, 40);
  h2_yVxpTarData    = new TH2F("h2_yVxpTarData",    "Data: Y_{tar} vs. X'_{tar}; X'_{tar} / 2 mrad; Y_{tar} / 1 mm",     200, -100.0, 100.0, 100, -5, 5);
  h2_yVypTarData    = new TH2F("h2_yVypTarData",    "Data: Y_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; Y_{tar} / 1 mm",     200, -100.0, 100.0, 100, -5, 5);
  h2_xpVypTarData   = new TH2F("h2_xpVypTarData",   "Data: X'_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; X'_{tar} / 2 mrad", 60, -60.0, 60.0, 100, -100.0, 100.0);
  compFile->cd("../");

  UInt_t mcRawEventCntr = 0;
  // Loop over raw MC data
  for (UInt_t imcRaw = 0; imcRaw < nEntriesMCRaw; imcRaw++) {
    // Obtain the data entry
    mcRawTree->GetEntry(imcRaw);
    if (mcRawEventCntr % 10000 == 0 && mcRawEventCntr != 0) 
      cout << mcRawEventCntr << " Input Monte-Carlo events have been processed..." << endl;
    mcRawEventCntr++;
    // Fill histos
    // 1D histos
    h_xFocalMCRaw->Fill(xFocalMCRaw);
    h_xpFocalMCRaw->Fill(xpFocalMCRaw*rad2mrad);
    h_yFocalMCRaw->Fill(yFocalMCRaw);
    h_ypFocalMCRaw->Fill(ypFocalMCRaw*rad2mrad);
    h_yTarMCRaw->Fill(yTarMCRaw);
    h_xpTarMCRaw->Fill(xpTarMCRaw*rad2mrad);
    h_ypTarMCRaw->Fill(ypTarMCRaw*rad2mrad);
    h_deltaMCRaw->Fill(deltaMCRaw);
    // 2D histos
    h2_xVxpFocalMCRaw->Fill(xpFocalMCRaw*rad2mrad, xFocalMCRaw);
    h2_xVyFocalMCRaw->Fill(yFocalMCRaw, xFocalMCRaw);
    h2_xVypFocalMCRaw->Fill(ypFocalMCRaw*rad2mrad, xFocalMCRaw);
    h2_xpVyFocalMCRaw->Fill(yFocalMCRaw, xpFocalMCRaw*rad2mrad);
    h2_xpVypFocalMCRaw->Fill(ypFocalMCRaw*rad2mrad, xpFocalMCRaw*rad2mrad);
    h2_yVypFocalMCRaw->Fill(ypFocalMCRaw*rad2mrad, yFocalMCRaw);
    h2_yVxpTarMCRaw->Fill(xpTarMCRaw*rad2mrad, yTarMCRaw);
    h2_yVypTarMCRaw->Fill(ypTarMCRaw*rad2mrad, yTarMCRaw);
    h2_xpVypTarMCRaw->Fill(ypTarMCRaw*rad2mrad, xpTarMCRaw*rad2mrad);
  }  // Raw MC loop

  UInt_t mcWgtEventCntr = 0;
  // Loop over weighted MC data
  for (UInt_t imcWgt = 0; imcWgt < nEntriesMCWgt; imcWgt++) {
    // Obtain the data entry
    mcWgtTree->GetEntry(imcWgt);
    if (mcWgtEventCntr % 10000 == 0 && mcWgtEventCntr != 0) 
      cout << mcWgtEventCntr << " Weighted Monte-Carlo events have been processed..." << endl;
    mcWgtEventCntr++;
    // Define and implement cuts
    Bool_t failIDCut = failID == 0.0;
    if (!failIDCut) continue;
    Bool_t deltaPCut = TMath::Abs(deltaMCWgt - 5.0) < 20.0;
    if (!deltaPCut) continue;
    // Caluclate the scale factor
    scaleFactor = (bornCS * mcScaleFactor) / (fudgeFactor * intRadCorr);
    // Fill histos
    // 1D histos
    h_xFocalMCWgt->Fill(xFocalMCWgt, scaleFactor);
    h_xpFocalMCWgt->Fill(xpFocalMCWgt*rad2mrad, scaleFactor);
    h_yFocalMCWgt->Fill(yFocalMCWgt, scaleFactor);
    h_ypFocalMCWgt->Fill(ypFocalMCWgt*rad2mrad, scaleFactor);
    h_yTarMCWgt->Fill(yTarMCWgt, scaleFactor);
    h_xpTarMCWgt->Fill(xpTarMCWgt*rad2mrad, scaleFactor);
    h_ypTarMCWgt->Fill(ypTarMCWgt*rad2mrad, scaleFactor);
    h_deltaMCWgt->Fill(deltaMCWgt, scaleFactor);
    h_thetaMCWgt->Fill(thetaMCWgt*TMath::RadToDeg(), scaleFactor);
    h_q2MCWgt->Fill(q2MCWgt, scaleFactor);
    h_w2MCWgt->Fill(w2MCWgt, scaleFactor);
    //h_thetaMCWgt->Fill(thetaMCWgt, scaleFactor);
    // 2D histos
    h2_xVxpFocalMCWgt->Fill(xpFocalMCWgt*rad2mrad, xFocalMCWgt, scaleFactor);
    h2_xVyFocalMCWgt->Fill(yFocalMCWgt, xFocalMCWgt, scaleFactor);
    h2_xVypFocalMCWgt->Fill(ypFocalMCWgt*rad2mrad, xFocalMCWgt, scaleFactor);
    h2_xpVyFocalMCWgt->Fill(yFocalMCWgt, xpFocalMCWgt*rad2mrad, scaleFactor);
    h2_xpVypFocalMCWgt->Fill(ypFocalMCWgt*rad2mrad, xpFocalMCWgt*rad2mrad, scaleFactor);
    h2_yVypFocalMCWgt->Fill(ypFocalMCWgt*rad2mrad, yFocalMCWgt, scaleFactor);
    h2_yVxpTarMCWgt->Fill(xpTarMCWgt*rad2mrad, yTarMCWgt, scaleFactor);
    h2_yVypTarMCWgt->Fill(ypTarMCWgt*rad2mrad, yTarMCWgt, scaleFactor);
    h2_xpVypTarMCWgt->Fill(ypTarMCWgt*rad2mrad, xpTarMCWgt*rad2mrad, scaleFactor);
  }  // Weighted MC loop
  
  UInt_t dataEventCntr = 0;
  // Loop over data
  for (UInt_t idata = 0; idata < nEntriesData; idata++) {
    // Obtain the data entry
    dataTree->GetEntry(idata);
    if (dataEventCntr % 10000 == 0 && dataEventCntr != 0) 
      cout << dataEventCntr << " Data events have been processed..." << endl;
    dataEventCntr++;
    // Define and implement cuts
    // PID cut in heavy gas Cherenkov
    Bool_t phgcerNpeSumCut = phgcerNpeSum >= 1.0;
    if (!phgcerNpeSumCut) continue;
    // Tracking cuts
    Bool_t gtrOkCut   = gtrOkData == 1.0;
    Bool_t gtrBetaCut = TMath::Abs(gtrBetaData - 1.0) < 0.2;
    Bool_t deltaPCut  = TMath::Abs(deltaData - 5.0) < 20.0;
    Bool_t trackCuts  = gtrOkData && gtrBetaCut && deltaPCut;
    //Bool_t trackEnergyNormCut = TMath::Abs(trackEnergyNorm - 1.1) < 0.3;
    //Bool_t trackCuts          = numTracksCut && goldTrackBetaCut && trackEnergyNormCut;
    if (!trackCuts) continue;

    // Fill histos
    // 1D histos
    h_xFocalData->Fill(xFocalData);
    h_xpFocalData->Fill(xpFocalData*rad2mrad);
    h_yFocalData->Fill(yFocalData);
    h_ypFocalData->Fill(ypFocalData*rad2mrad);
    h_yTarData->Fill(yTarData);
    h_xpTarData->Fill(xpTarData*rad2mrad);
    h_ypTarData->Fill(ypTarData*rad2mrad);
    h_deltaData->Fill(deltaData);
    h_thetaData->Fill((theta*TMath::DegToRad() + ypTarData)*TMath::RadToDeg());
    h_q2Data->Fill(calc_q2(E0, gtrMomData, theta*TMath::DegToRad() + ypTarData));
    h_w2Data->Fill(calc_w2(E0, gtrMomData, theta*TMath::DegToRad() + ypTarData));
    // 2D histos
    h2_xVxpFocalData->Fill(xpFocalData*rad2mrad, xFocalData);
    h2_xVyFocalData->Fill(yFocalData, xFocalData);
    h2_xVypFocalData->Fill(ypFocalData*rad2mrad, xFocalData);
    h2_xpVyFocalData->Fill(yFocalData, xpFocalData*rad2mrad);
    h2_xpVypFocalData->Fill(ypFocalData*rad2mrad, xpFocalData*rad2mrad);
    h2_yVypFocalData->Fill(ypFocalData*rad2mrad, yFocalData);
    h2_yVxpTarData->Fill(xpTarData*rad2mrad, yTarData);
    h2_yVypTarData->Fill(ypTarData*rad2mrad, yTarData);
    h2_xpVypTarData->Fill(ypTarData*rad2mrad, xpTarData*rad2mrad);
  }  // Data loop

  // Write comparison ROOT file
  compFile->Write();
  // Close comparison ROOT file
  compFile->Close();

}

// User specific funtions
Double_t calc_q2(Double_t beamEnergy, Double_t scatMom, Double_t scatAngle) {
  return 2.0*beamEnergy*scatMom*(1.0 - TMath::Cos(scatAngle));
}

Double_t calc_w2(Double_t beamEnergy, Double_t scatMom, Double_t scatAngle) {
  return protonMass*protonMass + 2.0*protonMass*(beamEnergy - scatMom) - 2.0*beamEnergy*scatMom*(1.0 - TMath::Cos(scatAngle));
}
