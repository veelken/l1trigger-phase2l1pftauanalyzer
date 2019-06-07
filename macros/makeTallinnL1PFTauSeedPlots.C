
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
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

double compIntegral(const TH1* histogram, bool includeUnderflowBin = false, bool includeOverflowBin = false)
{
  int firstBin = 1;
  if ( includeUnderflowBin ) firstBin -= 1;
  int lastBin = histogram->GetNbinsX();
  if ( includeOverflowBin ) lastBin += 1;
  double integral = 0.;
  for ( int idxBin = firstBin; idxBin <= lastBin; ++idxBin ) {
    integral += histogram->GetBinContent(idxBin);
  }
  return integral;
}

TH1* loadHistogram(TFile* inputFile, const std::string& histogramName, bool doNormalize)
{
  TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  if ( doNormalize ) {
    double integral = compIntegral(histogram);
    if ( integral > 0. ) {
      histogram->Scale(1./integral);
    }
  }
  return histogram;
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

void showHistograms(double canvasSizeX, double canvasSizeY,
                    TH1* histogram1, const std::string& legendEntry1,
                    TH1* histogram2, const std::string& legendEntry2,
                    TH1* histogram3, const std::string& legendEntry3,
                    TH1* histogram4, const std::string& legendEntry4,
                    TH1* histogram5, const std::string& legendEntry5,
                    TH1* histogram6, const std::string& legendEntry6,
                    int colors[], int lineStyles[], 
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

  histogram1->SetLineColor(colors[0]);
  histogram1->SetLineWidth(2);
  histogram1->SetLineStyle(lineStyles[0]);
  histogram1->Draw("hist");

  if ( histogram2 ) {
    histogram2->SetLineColor(colors[1]);
    histogram2->SetLineWidth(2);
    histogram2->SetLineStyle(lineStyles[1]);
    histogram2->Draw("histsame");
  }

  if ( histogram3 ) {
    histogram3->SetLineColor(colors[2]);
    histogram3->SetLineWidth(2);
    histogram3->SetLineStyle(lineStyles[2]);
    histogram3->Draw("histsame");
  }

  if ( histogram4 ) {
    histogram4->SetLineColor(colors[3]);
    histogram4->SetLineWidth(2);
    histogram4->SetLineStyle(lineStyles[3]);
    histogram4->Draw("histsame");
  }

  if ( histogram5 ) {
    histogram5->SetLineColor(colors[4]);
    histogram5->SetLineWidth(2);
    histogram5->SetLineStyle(lineStyles[4]);
    histogram5->Draw("histsame");
  }

  if ( histogram6 ) {
    histogram6->SetLineColor(colors[5]);
    histogram6->SetLineWidth(2);
    histogram6->SetLineStyle(lineStyles[5]);
    histogram6->Draw("histsame");
  }

  TLegend* legend = 0;
  if ( legendEntry1 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    legend->AddEntry(histogram1, legendEntry1.data(), "l");
    if ( histogram2 ) legend->AddEntry(histogram2, legendEntry2.data(), "l");
    if ( histogram3 ) legend->AddEntry(histogram3, legendEntry3.data(), "l");
    if ( histogram4 ) legend->AddEntry(histogram4, legendEntry4.data(), "l");
    if ( histogram5 ) legend->AddEntry(histogram5, legendEntry5.data(), "l");
    if ( histogram6 ) legend->AddEntry(histogram6, legendEntry6.data(), "l");
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

void makeTallinnL1PFTauSeedPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/L1Trigger/TallinnL1PFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName = "TallinnL1PFTauSeedAnalyzer_background_2019Jun03.root";
  TFile* inputFile = openFile(inputFilePath, inputFileName);

  std::vector<std::string> observables;
  observables.push_back("seedChargedPFCand_pt");
  observables.push_back("seedChargedPFCand_eta");
  observables.push_back("seedChargedPFCand_phi");
  observables.push_back("seedChargedPFCand_dz");
  observables.push_back("numSeedChargedPFCands");
  observables.push_back("seedPFJet_pt");
  observables.push_back("seedPFJet_eta");
  observables.push_back("seedPFJet_phi");
  observables.push_back("numSeedPFJets");

  std::map<std::string, int> rebin; // key   = observable
  rebin["seedChargedPFCand_pt"]         = 1;
  rebin["seedChargedPFCand_eta"]        = 1;
  rebin["seedChargedPFCand_phi"]        = 1;
  rebin["seedChargedPFCand_dz"]         = 1;
  rebin["numSeedChargedPFCands"]        = 1;
  rebin["seedPFJet_pt"]                 = 1;
  rebin["seedPFJet_eta"]                = 1;
  rebin["seedPFJet_phi"]                = 1;
  rebin["numSeedPFJets"]                = 1;
  rebin[""]                             = 1;

  std::map<std::string, bool> useLogScaleX; // key = observable
  useLogScaleX["seedChargedPFCand_pt"]  = false;
  useLogScaleX["seedChargedPFCand_eta"] = false;
  useLogScaleX["seedChargedPFCand_phi"] = false;
  useLogScaleX["seedChargedPFCand_dz"]  = false;
  useLogScaleX["numSeedChargedPFCands"] = false;
  useLogScaleX["seedPFJet_pt"]          = false;
  useLogScaleX["seedPFJet_eta"]         = false;
  useLogScaleX["seedPFJet_phi"]         = false;
  useLogScaleX["numSeedPFJets"]         = false;

  std::map<std::string, double> xMin; // key = observable
  xMin["seedChargedPFCand_pt"]          =   0.;
  xMin["seedChargedPFCand_eta"]         =  -3.0;
  xMin["seedChargedPFCand_phi"]         = -TMath::Pi();
  xMin["seedChargedPFCand_dz"]          =  -5.0;
  xMin["numSeedChargedPFCands"]         =  -0.5;
  xMin["seedPFJet_pt"]                  =   0.;
  xMin["seedPFJet_eta"]                 =  -3.0;
  xMin["seedPFJet_phi"]                 = -TMath::Pi();
  xMin["numSeedPFJets"]                 =  -0.5;

  std::map<std::string, double> xMax; // key = observable
  xMax["seedChargedPFCand_pt"]          =  50.;
  xMax["seedChargedPFCand_eta"]         =  +3.0;
  xMax["seedChargedPFCand_phi"]         = -TMath::Pi();
  xMax["seedChargedPFCand_dz"]          =  +5.0;
  xMax["numSeedChargedPFCands"]         =  24.5;
  xMax["seedPFJet_pt"]                  = 150.;
  xMax["seedPFJet_eta"]                 =  +3.0;
  xMax["seedPFJet_phi"]                 = -TMath::Pi();
  xMax["numSeedPFJets"]                 =  14.5;
  
  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["seedChargedPFCand_pt"]   = "seed charged PFCand p_{T} [GeV]";
  xAxisTitles["seedChargedPFCand_eta"]  = "seed charged PFCand #eta";
  xAxisTitles["seedChargedPFCand_phi"]  = "seed charged PFCand #phi";
  xAxisTitles["seedChargedPFCand_dz"]   = "seed charged PFCand #Delta_{z} [cm]";
  xAxisTitles["numSeedChargedPFCands"]  = "Number of seed charged PFCands";
  xAxisTitles["seedPFJet_pt"]           = "seed PFJet p_{T} [GeV]";
  xAxisTitles["seedPFJet_eta"]          = "seed PFJet #eta";
  xAxisTitles["seedPFJet_phi"]          = "seed PFJet #phi";
  xAxisTitles["numSeedPFJets"]          = "Number of seed PFJets";

  std::map<std::string, std::string> yAxisTitles; // key = observable
  yAxisTitles["seedChargedPFCand_pt"]   = "#frac{dN}{dp_{T}} [1/GeV]";
  yAxisTitles["seedChargedPFCand_eta"]  = "#frac{dN}{d#eta}";
  yAxisTitles["seedChargedPFCand_phi"]  = "#frac{dN}{d#phi}";
  yAxisTitles["seedChargedPFCand_dz"]   = "#frac{dN}{d#Delta_{z}} [1/cm]";
  yAxisTitles["numSeedChargedPFCands"]  = "Events";
  yAxisTitles["seedPFJet_pt"]           = "#frac{dN}{dp_{T}} [1/GeV]";
  yAxisTitles["seedPFJet_eta"]          = "#frac{dN}{d#eta}";
  yAxisTitles["seedPFJet_phi"]          = "#frac{dN}{d#phi}";
  yAxisTitles["numSeedPFJets"]          = "Events";

  std::string dqmDirectory = "DQMData/TallinnL1PFTauSeedAnalyzer";
  
  int colors[6] = { 1, 2, 8, 4, 6, 7 };
  int lineStyles[6] = { 1, 1, 1, 1, 1, 1 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  std::vector<std::string> labelTextLines;

  for ( std::vector<std::string>::const_iterator observable = observables.begin();
	observable != observables.end(); ++observable ) {
    std::string histogramName = Form("%s/%s", dqmDirectory.data(), observable->data()); 
    TH1* histogram = loadHistogram(inputFile, histogramName, true);
    TH1* histogram_rebinned = (TH1*)histogram->Clone(Form("%s_rebinned", histogram->GetName()));
    if ( rebin[*observable] > 1 ) histogram_rebinned->Rebin(rebin[*observable]);
    divideByBinWidth(histogram_rebinned);

    std::string outputFileName = Form("makeTallinnL1PFTauSeedPlots_%s.png", observable->data());
    showHistograms(1150, 850,
		   histogram_rebinned, "",
		   0, "",
		   0, "",
		   0, "",
		   0, "",
		   0, "",
		   colors, lineStyles, 
		   0.050, 0.77, 0.74, 0.16, 0.18, 
		   labelTextLines, 0.050,
		   0.70, 0.62, 0.23, 0.06, 
		   useLogScaleX[*observable], xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.3,
		   true, 1.e-5, 1.99, yAxisTitles[*observable], 1.4, 
		   outputFileName);

    delete histogram_rebinned;
  }

  delete inputFile;
}

