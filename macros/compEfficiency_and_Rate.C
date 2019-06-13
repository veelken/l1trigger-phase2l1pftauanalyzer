
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include <limits>

TFile* openFile(const std::string& inputFilePath, const std::string& inputFileName)
{
  std::string inputFileName_full = inputFilePath;
  if ( inputFileName_full.find_last_of("/") != (inputFileName_full.size() - 1) ) inputFileName_full.append("/");
  inputFileName_full.append(inputFileName);
  TFile* inputFile = new TFile(inputFileName_full.data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = '" << inputFileName_full << "' !!" << std::endl;
    assert(0);
  }
  return inputFile;
}

TH2* loadHistogram2d(TFile* inputFile, const std::string& histogramName)
{
  TH2* histogram = dynamic_cast<TH2*>(inputFile->Get(histogramName.data()));
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  return histogram;
}

void printEfficiency_and_Rate(const TH2* histogram2d_signal, const TH2* histogram2d_background, const std::string& header, 
			      double rate_target, double rate_accCorrFactor)
{
  assert(histogram2d_signal->GetNbinsX() == histogram2d_background->GetNbinsX());
  const TAxis* xAxis = histogram2d_signal->GetXaxis();
  int numBinsX = xAxis->GetNbins();
  assert(histogram2d_signal->GetNbinsY() == histogram2d_background->GetNbinsY());
  const TAxis* yAxis = histogram2d_signal->GetYaxis();
  int numBinsY = yAxis->GetNbins();

  double tau_ptThreshold_symmetric = -1.;
  double best_rate_symmetric = -1.;
  double best_efficiency_symmetric = -1.;
  double leadingTau_ptThreshold_asymmetric = -1.;
  double subleadingTau_ptThreshold_asymmetric = -1.;
  double best_rate_asymmetric = -1.;
  double best_efficiency_asymmetric = -1.;

  for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX ) {
    double binCenterX = xAxis->GetBinCenter(idxBinX);
    for ( int idxBinY = 1; idxBinY <= numBinsY; ++idxBinY ) {
      double binCenterY = yAxis->GetBinCenter(idxBinY);
      
      double rate = histogram2d_background->GetBinContent(idxBinX, idxBinY);
      rate *= rate_accCorrFactor; // extrapolation from |eta| < 1.4 to full HL-LHC tracking acceptance (squared for ditau trigger)
      rate *= 2.8e+7;             // bunch-crossing frequency of colliding bunches = 28 MHz
      double efficiency = histogram2d_signal->GetBinContent(idxBinX, idxBinY);
      
      const double epsilon = 1.e-3;
      if ( rate < rate_target && efficiency > best_efficiency_symmetric && TMath::Abs(binCenterX - binCenterY) < epsilon ) {
	tau_ptThreshold_symmetric = binCenterX;
	best_rate_symmetric = rate;
	best_efficiency_symmetric = efficiency;
      }
      
      if ( rate < rate_target && efficiency > best_efficiency_asymmetric ) {
	leadingTau_ptThreshold_asymmetric = binCenterX;
	subleadingTau_ptThreshold_asymmetric = binCenterY;
	best_rate_asymmetric = rate;
	best_efficiency_asymmetric = efficiency;
      }
    }	
  }
  
  std::cout << header << std::endl;
  std::cout << "results for symmetric pT thresholds " 
	    << "(tau pT > " << tau_ptThreshold_symmetric << " GeV)" << std::endl;
  std::cout << " rate = " << best_rate_symmetric << std::endl;
  std::cout << " efficiency = " << best_efficiency_symmetric << std::endl;
  std::cout << "results for asymmetric pT thresholds " 
	    << "(leading tau pT > " << leadingTau_ptThreshold_asymmetric << "," 
	    << " subleading tau pT > " << subleadingTau_ptThreshold_asymmetric << " GeV)" << std::endl;
  std::cout << " rate = " << best_rate_asymmetric << std::endl;
  std::cout << " efficiency = " << best_efficiency_asymmetric << std::endl;	
}

void compEfficiency_and_Rate()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/L1Trigger/TallinnL1PFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName_signal = "TallinnL1PFTauAnalyzer_signal_2019Jun12.root";
  TFile* inputFile_signal = openFile(inputFilePath, inputFileName_signal);
  std::string inputFileName_background = "TallinnL1PFTauAnalyzer_background_2019Jun12.root";
  TFile* inputFile_background = openFile(inputFilePath, inputFileName_background);

  std::vector<std::string> pfAlgos;
  pfAlgos.push_back("WithStripsAndPreselectionPF");
  pfAlgos.push_back("WithStripsWithoutPreselectionPF");
  pfAlgos.push_back("WithoutStripsWithPreselectionPF");
  pfAlgos.push_back("WithoutStripsAndPreselectionPF");
  pfAlgos.push_back("WithStripsAndPreselectionPuppi");
  pfAlgos.push_back("WithStripsWithoutPreselectionPuppi");
  pfAlgos.push_back("WithoutStripsWithPreselectionPuppi");
  pfAlgos.push_back("WithoutStripsAndPreselectionPuppi");

  std::vector<std::string> absEtaRanges;
  //absEtaRanges.push_back("absEtaLt1p00");
  absEtaRanges.push_back("absEtaLt1p40");

  std::vector<std::string> isolationWPs;
  isolationWPs.push_back("relChargedIsoLt0p40");
  isolationWPs.push_back("relChargedIsoLt0p20");
  isolationWPs.push_back("relChargedIsoLt0p10");
  isolationWPs.push_back("relChargedIsoLt0p05");

  double rate_target        = 12.e+3; // 12 kHz
  double rate_accCorrFactor =  4.;    // extrapolation from |eta| < 1.4 to full HL-LHC tracking acceptance (squared for ditau trigger)

  std::string dqmDirectory = "DQMData/TallinnL1PFTauPairAnalyzer";

  std::vector<std::string> labelTextLines;

  int colors[6] = { 1, 2, 8, 4, 6, 7 };
  int lineStyles[6] = { 1, 1, 1, 1, 1, 1 };

  // TallinnL1PFTaus 
  for ( std::vector<std::string>::const_iterator pfAlgo = pfAlgos.begin();
	pfAlgo != pfAlgos.end(); ++pfAlgo ) {
    for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	  absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
      for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs.begin();
	    isolationWP != isolationWPs.end(); ++isolationWP ) {
	std::string histogram2dName_signal = Form("%s%s_wrtOfflineTaus/efficiency_or_rate_%s_%s", 
          dqmDirectory.data(), pfAlgo->data(), absEtaRange->data(), isolationWP->data());
        TH2* histogram2d_signal = loadHistogram2d(inputFile_signal, histogram2dName_signal);
	std::string histogram2dName_background = Form("%s%s/efficiency_or_rate_%s_%s", 
          dqmDirectory.data(), pfAlgo->data(), absEtaRange->data(), isolationWP->data());
	TH2* histogram2d_background = loadHistogram2d(inputFile_background, histogram2dName_background);

	std::string header = Form("srcL1PFTaus = TallinnL1PFTauProducer%s, absEtaRange = %s, isolationWP = %s:", pfAlgo->data(), absEtaRange->data(), isolationWP->data());
	printEfficiency_and_Rate(histogram2d_signal, histogram2d_background, header, rate_target, rate_accCorrFactor);
      }
    }
  }

  // Isobel's L1PFTaus 
  std::vector<std::string> isolationWPs_isobel;
  isolationWPs_isobel.push_back("vLooseIso");
  isolationWPs_isobel.push_back("LooseIso");
  isolationWPs_isobel.push_back("MediumIso");
  isolationWPs_isobel.push_back("TightIso");

  std::string dqmDirectory_isobel = "DQMData/L1PFTauPairAnalyzer";

  for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs_isobel.begin();
	isolationWP != isolationWPs_isobel.end(); ++isolationWP ) {
    std::string histogram2dName_signal = Form("%sPF_wrtOfflineTaus/efficiency_or_rate_absEtaLt1p40_%s", 
      dqmDirectory_isobel.data(), isolationWP->data());
    TH2* histogram2d_signal = loadHistogram2d(inputFile_signal, histogram2dName_signal);
    std::string histogram2dName_background = Form("%sPF/efficiency_or_rate_absEtaLt1p40_%s", 
      dqmDirectory_isobel.data(), isolationWP->data());
    TH2* histogram2d_background = loadHistogram2d(inputFile_background, histogram2dName_background);

    std::string header = Form("srcL1PFTaus = L1PFTauProducer:L1PFTaus, isolationWP = %s:", isolationWP->data());
    printEfficiency_and_Rate(histogram2d_signal, histogram2d_background, header, rate_target, rate_accCorrFactor);
  }

  delete inputFile_signal;
  delete inputFile_background;
}

