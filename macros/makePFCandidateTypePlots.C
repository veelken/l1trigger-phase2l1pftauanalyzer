
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
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

TH1* loadHistogram(TFile* inputFile, const std::string& histogramName)
{
  TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  return histogram;
}

void multiplyByBinWidth(TH1* histogram)
{
  if ( !histogram ) return;
  TAxis* xAxis = histogram->GetXaxis();
  int numBins = xAxis->GetNbins();
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin ) {
    double binContent = histogram->GetBinContent(idxBin);
    double binError = histogram->GetBinError(idxBin);
    double binWidth = xAxis->GetBinWidth(idxBin);
    histogram->SetBinContent(idxBin, binContent*binWidth);
    histogram->SetBinError(idxBin, binError*binWidth);
  }
}

void divideByBinWidth(TH1* histogram)
{
  if ( !histogram ) return;
  TAxis* xAxis = histogram->GetXaxis();
  int numBins = xAxis->GetNbins();
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin ) {
    double binContent = histogram->GetBinContent(idxBin);
    double binError = histogram->GetBinError(idxBin);
    double binWidth = xAxis->GetBinWidth(idxBin);
    histogram->SetBinContent(idxBin, binContent/binWidth);
    histogram->SetBinError(idxBin, binError/binWidth);
  }
}

void showHistograms_stacked(double canvasSizeX, double canvasSizeY,
			    const std::string& pfCandType1, TH1* histogram1,
			    const std::string& pfCandType2, TH1* histogram2,
			    const std::string& pfCandType3, TH1* histogram3,
			    const std::string& pfCandType4, TH1* histogram4,
			    const std::string& pfCandType5, TH1* histogram5,
			    const std::string& pfCandType6, TH1* histogram6,
			    std::map<std::string, std::string>& legendEntries,
			    bool doNormalize,
			    std::map<std::string, int>& colors, std::map<std::string, int>& fillStyles, 
			    double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
			    std::vector<std::string>& labelTextLines, double labelTextSize,
			    double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
			    bool useLogScaleX, double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
			    bool useLogScaleY, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
			    const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  
  canvas->SetTopMargin(0.05);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.14);
  canvas->SetRightMargin(0.05);

  canvas->SetLogx(useLogScaleX);
  canvas->SetLogy(useLogScaleY);
  
  canvas->SetGridx(1);
  canvas->SetGridy(1);

  if ( !histogram1 ) {
    std::cerr << "<showHistograms>: histogram1 = NULL --> skipping !!" << std::endl;
    return;
  }

  histogram1->SetTitle("");
  histogram1->SetStats(false);
  histogram1->SetMinimum(yMin);
  histogram1->SetMaximum(yMax);

  TAxis* xAxis = histogram1->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  
  if ( xMax > xMin ) {
    std::cout << "limiting x-axis range to " << xMin << ".." << xMax << std::endl;
    xAxis->SetRangeUser(xMin, xMax);
  }

  TAxis* yAxis = histogram1->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  histogram1->Draw("axis");

  TH1* histogramSum = nullptr;
  if ( doNormalize ) {
    histogramSum = (TH1*)histogram1->Clone("histogramSum");
    if ( histogram2 ) histogramSum->Add(histogram2);
    if ( histogram3 ) histogramSum->Add(histogram3);
    if ( histogram4 ) histogramSum->Add(histogram4);
    if ( histogram5 ) histogramSum->Add(histogram5);
    if ( histogram6 ) histogramSum->Add(histogram6);
  }

  THStack* histogramStack = new THStack("histogramStack", "histogramStack");

  if ( doNormalize ) {
    histogram1->Divide(histogramSum);
  }
  histogram1->SetFillColor(colors[pfCandType1]);
  histogram1->SetFillStyle(fillStyles[pfCandType1]);
  histogram1->SetLineColor(1);
  histogram1->SetLineWidth(1);
  histogram1->SetLineStyle(1);
  histogramStack->Add(histogram1);

  if ( histogram2 ) {
    if ( doNormalize ) {
      histogram2->Divide(histogramSum);
    }
    histogram2->SetFillColor(colors[pfCandType2]);
    histogram2->SetFillStyle(fillStyles[pfCandType2]);
    histogram2->SetLineColor(1);
    histogram2->SetLineWidth(1);
    histogram2->SetLineStyle(1);
    histogramStack->Add(histogram2);
  }

  if ( histogram3 ) {
    if ( doNormalize ) {
      histogram3->Divide(histogramSum);
    }
    histogram3->SetFillColor(colors[pfCandType3]);
    histogram3->SetFillStyle(fillStyles[pfCandType3]);
    histogram3->SetLineColor(1);
    histogram3->SetLineWidth(1);
    histogram3->SetLineStyle(1);
    histogramStack->Add(histogram1);
  }

  if ( histogram4 ) {
    if ( doNormalize ) {
      histogram4->Divide(histogramSum);
    }
    histogram4->SetFillColor(colors[pfCandType4]);
    histogram4->SetFillStyle(fillStyles[pfCandType4]);
    histogram4->SetLineColor(1);
    histogram4->SetLineWidth(1);
    histogram4->SetLineStyle(1);
    histogramStack->Add(histogram4);
  }

  if ( histogram5 ) {
    if ( doNormalize ) {
      histogram5->Divide(histogramSum);
    }
    histogram5->SetFillColor(colors[pfCandType5]);
    histogram5->SetFillStyle(fillStyles[pfCandType5]);
    histogram5->SetLineColor(1);
    histogram5->SetLineWidth(1);
    histogram5->SetLineStyle(1);
    histogramStack->Add(histogram5);
  }

  if ( histogram6 ) {
    if ( doNormalize ) {
      histogram6->Divide(histogramSum);
    }
    histogram6->SetFillColor(colors[pfCandType6]);
    histogram6->SetFillStyle(fillStyles[pfCandType6]);
    histogram6->SetLineColor(1);
    histogram6->SetLineWidth(1);
    histogram6->SetLineStyle(1);
    histogramStack->Add(histogram6);
  }

  histogramStack->Draw("same");

  TLegend* legend = 0;
  if ( legendEntries[pfCandType1] != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    legend->AddEntry(histogram1, legendEntries[pfCandType1].data(), "f");
    if ( histogram2 ) legend->AddEntry(histogram2, legendEntries[pfCandType2].data(), "f");
    if ( histogram3 ) legend->AddEntry(histogram3, legendEntries[pfCandType3].data(), "f");
    if ( histogram4 ) legend->AddEntry(histogram4, legendEntries[pfCandType4].data(), "f");
    if ( histogram5 ) legend->AddEntry(histogram5, legendEntries[pfCandType5].data(), "f");
    if ( histogram6 ) legend->AddEntry(histogram6, legendEntries[pfCandType6].data(), "f");
    legend->Draw();
  }

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

  histogram1->Draw("axissame");

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete histogramSum;
  delete histogramStack;
  delete label;
  delete legend;
  delete canvas;  
}

double compIntegral_within_Range(const TH1* histogram, double xMin, double xMax)
{
  const TAxis* xAxis = histogram->GetXaxis();
  double integral = 0.;
  int numBinsX = xAxis->GetNbins();
  for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX ) { 
    double binCenter = xAxis->GetBinCenter(idxBinX);
    double binContent = histogram->GetBinContent(idxBinX);
    if ( binCenter >= xMin && binCenter < xMax ) {
      integral += binContent;
    }
  }
  return integral;
}

void makePFCandidateTypePlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/L1Trigger/TallinnL1PFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName = "L1PFCandidateTypeAnalyzer_signal_2019May28.root";
  TFile* inputFile = openFile(inputFilePath, inputFileName);

  std::vector<std::string> pfAlgos;
  pfAlgos.push_back("PF");
  pfAlgos.push_back("Puppi");

  std::vector<std::string> observables;
  observables.push_back("eta");
  observables.push_back("pt");

  std::map<std::string, int> rebin; // key = observable
  rebin["eta"]                         = 1;
  rebin["pt"]                          = 1;

  std::map<std::string, bool> useLogScaleX; // key = observable
  useLogScaleX["eta"]                  = false;
  useLogScaleX["pt"]                   = true;

  std::map<std::string, double> xMin; // key = observable
  xMin["pt"]                           =    0.;
  xMin["eta"]                          =   -3.0;

  std::map<std::string, double> xMax; // key = observable
  xMax["pt"]                           = 1000.;
  xMax["eta"]                          =   +3.0;
  
  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["pt"]                    = "p_{T} [GeV]";
  xAxisTitles["eta"]                   = "#eta";

  std::string dqmDirectory = "DQMData/L1PFCandidateTypeAnalyzer";
  
  std::vector<std::string> pfCandTypes;
  pfCandTypes.push_back("chargedHadronPileup");
  pfCandTypes.push_back("chargedHadron");
  pfCandTypes.push_back("photon");
  pfCandTypes.push_back("neutralHadron");
  pfCandTypes.push_back("electron");
  pfCandTypes.push_back("muon");

  std::map<std::string, int> colors; // key = pfCandType
  colors["chargedHadronPileup"]        = kRed - 4;
  colors["chargedHadron"]              = kRed;
  colors["photon"]                     = 9;
  colors["neutralHadron"]              = 8;
  colors["electron"]                   = kCyan;
  colors["muon"]                       = kCyan + 2;

  std::map<std::string, int> fillStyles; // key = pfCandType
  fillStyles["chargedHadronPileup"]    = 3001;
  fillStyles["chargedHadron"]          = 3001;
  fillStyles["photon"]                 = 3001;
  fillStyles["neutralHadron"]          = 3001;
  fillStyles["electron"]               = 3001;
  fillStyles["muon"]                   = 3001;

  std::map<std::string, std::string> legendEntries; // key = pfCandType
  legendEntries["chargedHadronPileup"] = "Charged PU hadrons";
  legendEntries["chargedHadron"]       = "Charged hadrons";
  legendEntries["photon"]              = "Photons";
  legendEntries["neutralHadron"]       = "Neutral hadrons";
  legendEntries["electron"]            = "Electrons";
  legendEntries["muon"]                = "Muons";

  std::vector<std::string> labelTextLines;

  for ( std::vector<std::string>::const_iterator pfAlgo = pfAlgos.begin();
	pfAlgo != pfAlgos.end(); ++pfAlgo ) {
    for ( std::vector<std::string>::const_iterator observable = observables.begin();
	  observable != observables.end(); ++observable ) {
      std::map<std::string, TH1*> histograms_ptFraction_rebinned; // key = pfCandType
      for ( std::vector<std::string>::const_iterator pfCandType = pfCandTypes.begin();
	    pfCandType != pfCandTypes.end(); ++pfCandType ) {
	std::string histogramName_ptFraction = Form("%s%s/%sPtFraction_vs_%s", dqmDirectory.data(), pfAlgo->data(), pfCandType->data(), observable->data());
	TH1* histogram_ptFraction = loadHistogram(inputFile, histogramName_ptFraction);
	TH1* histogram_ptFraction_rebinned = (TH1*)histogram_ptFraction->Clone(Form("%s_rebinned", histogram_ptFraction->GetName()));
	if ( rebin[*observable] > 1 ) {
	  multiplyByBinWidth(histogram_ptFraction_rebinned);
	  histogram_ptFraction_rebinned = histogram_ptFraction->Rebin(rebin[*observable]);
	  divideByBinWidth(histogram_ptFraction_rebinned);
	}
      }

      std::string outputFileName_unnormalized = Form("makePFCandidateTypePlots_%s_%s_unnormalized.png", pfAlgo->data(), observable->data());
      showHistograms_stacked(1150, 1150,
			     "chargedHadronPileup", histograms_ptFraction_rebinned["chargedHadronPileup"],
			     "chargedHadron",       histograms_ptFraction_rebinned["chargedHadron"], 
			     "photon",              histograms_ptFraction_rebinned["photon"], 
			     "neutralHadron",       histograms_ptFraction_rebinned["neutralHadron"],
			     "electron",            histograms_ptFraction_rebinned["electron"], 
			     "muon",                histograms_ptFraction_rebinned["muon"], 
			     legendEntries,
			     false,
			     colors, fillStyles, 
			     0.040, 0.17, 0.17, 0.23, 0.27, 
			     labelTextLines, 0.040,
			     0.70, 0.21, 0.23, 0.06, 
			     useLogScaleX[*observable], xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
			     true, 1.e-1, 1.99e+2, "#Sigma p_{T} [GeV]", 1.4, 
			     outputFileName_unnormalized);

      std::map<std::string, TH1*> histograms_energyFraction_rebinned; // key = pfCandType
      for ( std::vector<std::string>::const_iterator pfCandType = pfCandTypes.begin();
	    pfCandType != pfCandTypes.end(); ++pfCandType ) {
	std::string histogramName_energyFraction = Form("%s%s/%sPtFraction_vs_%s", dqmDirectory.data(), pfAlgo->data(), pfCandType->data(), observable->data());
	TH1* histogram_energyFraction = loadHistogram(inputFile, histogramName_energyFraction);
	TH1* histogram_energyFraction_rebinned = (TH1*)histogram_energyFraction->Clone(Form("%s_rebinned", histogram_energyFraction->GetName()));
	if ( rebin[*observable] > 1 ) {
	  multiplyByBinWidth(histogram_energyFraction_rebinned);
	  histogram_energyFraction_rebinned = histogram_energyFraction->Rebin(rebin[*observable]);
	  divideByBinWidth(histogram_energyFraction_rebinned);
	}
      }

      std::string outputFileName_normalized = Form("makePFCandidateTypePlots_%s_%s_unnormalized.png", pfAlgo->data(), observable->data());
      showHistograms_stacked(1150, 1150,
			     "chargedHadronPileup", histograms_energyFraction_rebinned["chargedHadronPileup"], 
			     "chargedHadron",       histograms_energyFraction_rebinned["chargedHadron"],
			     "photon",              histograms_energyFraction_rebinned["photon"],       
			     "neutralHadron",       histograms_energyFraction_rebinned["neutralHadron"],
			     "electron",            histograms_energyFraction_rebinned["electron"],
			     "muon",                histograms_energyFraction_rebinned["muon"],
			     legendEntries,
			     true,
			     colors, fillStyles, 
			     0.040, 0.17, 0.17, 0.23, 0.27, 
			     labelTextLines, 0.040,
			     0.70, 0.21, 0.23, 0.06, 
			     useLogScaleX[*observable], xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
			     false, 0., 1., "Energy fraction", 1.4, 
			     outputFileName_normalized);
    }
  }

  delete inputFile;
}
