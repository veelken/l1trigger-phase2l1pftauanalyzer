
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
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

TH2* makeEfficiencyHistogram2d(TH2* histogram_numerator, TH2* histogram_denominator)
{  
  assert(histogram_numerator->GetNbinsX() == histogram_denominator->GetNbinsX());
  int numBinsX = histogram_numerator->GetNbinsX();
  assert(histogram_numerator->GetNbinsY() == histogram_denominator->GetNbinsY());
  int numBinsY = histogram_numerator->GetNbinsY();
  TString histogramName_efficiency = TString(histogram_numerator->GetName()).ReplaceAll("_numerator", "");
  TH2* histogram_efficiency = (TH2*)histogram_numerator->Clone(histogramName_efficiency.Data());
  histogram_efficiency->Reset();
  for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX ) { 
    for ( int idxBinY = 1; idxBinY <= numBinsY; ++idxBinY ) {
      double binContent_numerator = histogram_numerator->GetBinContent(idxBinX, idxBinY);
      double binError_numerator = histogram_numerator->GetBinError(idxBinX, idxBinY);
      double binContent_denominator = histogram_denominator->GetBinContent(idxBinX, idxBinY);
      double binError_denominator = histogram_denominator->GetBinError(idxBinX, idxBinY);
      std::cout << " bin #" << idxBinX << ": numerator = " << binContent_numerator << ", denominator = " << binContent_denominator << std::endl;
      if( binContent_numerator > binContent_denominator ) {
	std::cerr << "Error in <makeEfficiencyGraph>: numerator = " << binContent_numerator << " exceeds denominator = " << binContent_denominator 
		  << " @ x = " << histogram_denominator->GetBinCenter(idxBinX) << " !!" << std::endl;
	assert(0);
      }
      TH1* histogram_numerator_dummy = new TH1D("histogram_numerator_dummy", "histogram_numerator_dummy", 1, -0.5, +0.5);
      histogram_numerator_dummy->SetBinContent(1, binContent_numerator);
      histogram_numerator_dummy->SetBinError(1, binError_numerator);
      TH1* histogram_denominator_dummy = new TH1D("histogram_denominator_dummy", "histogram_denominator_dummy", 1, -0.5, +0.5);
      histogram_denominator_dummy->SetBinContent(1, binContent_denominator);
      histogram_denominator_dummy->SetBinError(1, binError_denominator);
      TGraphAsymmErrors* graph_efficiency_dummy = new TGraphAsymmErrors(histogram_numerator_dummy, histogram_denominator_dummy, "w");
      double dummy, efficiency;
      graph_efficiency_dummy->GetPoint(0, dummy, efficiency);
      double efficiencyErrUp = graph_efficiency_dummy->GetErrorYhigh(0);
      double efficiencyErrDown = graph_efficiency_dummy->GetErrorYlow(0);
      double efficiencyErr = TMath::Sqrt(0.5*(efficiencyErrUp*efficiencyErrUp + efficiencyErrDown*efficiencyErrDown));
      histogram_efficiency->SetBinContent(idxBinX, idxBinY, efficiency);
      histogram_efficiency->SetBinError(idxBinX, idxBinY, efficiencyErr);
      delete histogram_numerator_dummy;
      delete histogram_denominator_dummy;
      delete graph_efficiency_dummy;
    }
  }
  return histogram_efficiency;
}

void showHistogram2d(double canvasSizeX, double canvasSizeY,
		     TH2* histogram,
		     std::vector<std::string>& labelTextLines, double labelTextSize,
		     double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		     const std::string& xAxisTitle, double xAxisOffset,
		     const std::string& yAxisTitle, double yAxisOffset,
		     const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  
  canvas->SetTopMargin(0.05);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.14);
  canvas->SetRightMargin(0.05);

  histogram->SetTitle("");
  histogram->SetStats(false);

  TAxis* xAxis = histogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  

  TAxis* yAxis = histogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  //histogram->Draw("COLZ");
  histogram->SetMarkerSize(0.5);
  histogram->Draw("TEXT");

  TPaveText* label = 0;
  if ( labelTextLines.size() > 0 ) {
    label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "brNDC");
    for ( std::vector<std::string>::const_iterator labelTextLine = labelTextLines.begin();
          labelTextLine != labelTextLines.end(); ++labelTextLine ) {
      label->AddText(labelTextLine->data());
    }
    label->SetFillColor(10);
    label->SetBorderSize(0);
    label->SetTextColor(1);
    label->SetTextAlign(12);
    label->SetTextSize(labelTextSize);
    label->Draw();
  }

  histogram->Draw("axissame");

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete label;
  delete canvas;  
}

void makeTrackingEfficiencyPlots2d()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/L1Trigger/TallinnL1PFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName = "TallinnL1PFTauAnalyzer_signal_ggH_2019Jul05.root";
  std::string inputFileName_full = inputFilePath;
  if ( inputFileName_full.find_last_of("/") != (inputFileName_full.size() - 1) ) inputFileName_full.append("/");
  inputFileName_full.append(inputFileName);
  TFile* inputFile = new TFile(inputFileName_full.data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = '" << inputFileName_full << "' !!" << std::endl;
    assert(0);
  }

  std::vector<std::string> recTrack_types;
  recTrack_types.push_back("l1Track");
  recTrack_types.push_back("l1PFCandTrack");
  recTrack_types.push_back("offlineTrack");
  recTrack_types.push_back("offlinePFCandTrack");

  std::vector<std::string> recTrack_options;
  recTrack_options.push_back("woQualityCuts");
  recTrack_options.push_back("wQualityCuts");

  std::vector<std::string> vtxModes;
  vtxModes.push_back("GenVertex");
  vtxModes.push_back("RecVertex");

  std::vector<std::string> observables;
  observables.push_back("pt_vs_absEta");

  std::vector<std::string> decayModes;
  decayModes.push_back("oneProng0Pi0");
  decayModes.push_back("oneProng1Pi0");
  decayModes.push_back("oneProng2Pi0");
  decayModes.push_back("threeProng0Pi0");
  decayModes.push_back("threeProng1Pi0");
  decayModes.push_back("all");

  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["pt_vs_absEta"] = "|#eta|";

  std::map<std::string, std::string> yAxisTitles; // key = observable
  yAxisTitles["pt_vs_absEta"] = "p_{T} [GeV]";

  std::string dqmDirectory = "DQMData/L1TrackAnalyzer";

  for ( std::vector<std::string>::const_iterator vtxMode = vtxModes.begin();
	vtxMode != vtxModes.end(); ++vtxMode ) {
    for ( std::vector<std::string>::const_iterator recTrack_type = recTrack_types.begin();
	  recTrack_type != recTrack_types.end(); ++recTrack_type ) {
      for ( std::vector<std::string>::const_iterator recTrack_option = recTrack_options.begin();
	    recTrack_option != recTrack_options.end(); ++recTrack_option ) {
        for ( std::vector<std::string>::const_iterator observable = observables.begin();
	      observable != observables.end(); ++observable ) {
  	  for ( std::vector<std::string>::const_iterator decayMode = decayModes.begin();
	        decayMode != decayModes.end(); ++decayMode ) {
	    std::string recTrack_type_capitalized = *recTrack_type;
	    recTrack_type_capitalized[0] = toupper(recTrack_type_capitalized[0]);
	    std::string decayMode_capitalized = *decayMode;
	    decayMode_capitalized[0] = toupper(decayMode_capitalized[0]);
	    std::string dqmDirectory_full = Form("%sWrt%s/%s/gen%sTau/%s_%s", 
              dqmDirectory.data(), vtxMode->data(), "absEtaLt2p40", decayMode_capitalized.data(), recTrack_type->data(), recTrack_option->data());	    
	    std::string histogramName_numerator = Form("%s/eff%s_%s_vs_%s_numerator_%s_pt1p00to1000p00", 
	      dqmDirectory_full.data(), recTrack_type_capitalized.data(), recTrack_option->data(), observable->data(), "absEtaLt2p40");
	    if ( (*decayMode) != "all" ) histogramName_numerator.append(Form("_gen%sTau", decayMode_capitalized.data()));
	    TH2* histogram_numerator = loadHistogram2d(inputFile, histogramName_numerator);
	    std::string histogramName_denominator = Form("%s/eff%s_%s_vs_%s_denominator_%s_pt1p00to1000p00", 
	      dqmDirectory_full.data(), recTrack_type_capitalized.data(), recTrack_option->data(), observable->data(), "absEtaLt2p40");
	    if ( (*decayMode) != "all" ) histogramName_denominator.append(Form("_gen%sTau", decayMode_capitalized.data()));
	    TH2* histogram_denominator = loadHistogram2d(inputFile, histogramName_denominator);
	    TH2* histogram_efficiency = makeEfficiencyHistogram2d(histogram_numerator, histogram_denominator);

	    std::vector<std::string> labelTextLines;
	    std::string outputFileName = Form("makeTrackingEfficiencyPlots2d_%s_%s_%s_wrt%s_%s.png", 
	    observable->data(), recTrack_type->data(), recTrack_option->data(), vtxMode->data(), decayMode->data());
	    showHistogram2d(1150, 1150,
		  	    histogram_efficiency,
			    labelTextLines, 0.050,
			    0.63, 0.65, 0.26, 0.07, 
			    xAxisTitles[*observable], 1.2, 
			    yAxisTitles[*observable], 1.4, 
			    outputFileName);
	  }
	}
      }
    }
  }

  delete inputFile;
}

