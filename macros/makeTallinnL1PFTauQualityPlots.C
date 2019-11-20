
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
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

TH1* loadHistogram(TFile* inputFile, const std::string& histogramName)
{
  TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  double integral = compIntegral(histogram);
  if ( integral > 0. ) {
    histogram->Scale(1./integral);
  }
  return histogram;
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
                    double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
                    bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
                    const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  
  canvas->SetTopMargin(0.05);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.14);
  canvas->SetRightMargin(0.05);

  canvas->SetLogy(useLogScale);
  
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

void makeTallinnL1PFTauQualityPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/L1Trigger/TallinnL1PFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName_signal = "TallinnL1PFTauQualityAnalyzer_signal_qqH_2019Nov20.root";
  TFile* inputFile_signal = openFile(inputFilePath, inputFileName_signal);
  std::string inputFileName_background = "TallinnL1PFTauQualityAnalyzer_background_2019Nov20.root";
  TFile* inputFile_background = openFile(inputFilePath, inputFileName_background);

  std::vector<std::string> pfAlgos;
  //pfAlgos.push_back("WithStripsAndPreselectionPF");
  pfAlgos.push_back("WithStripsWithoutPreselectionPF");
  //pfAlgos.push_back("WithoutStripsWithPreselectionPF");
  pfAlgos.push_back("WithoutStripsAndPreselectionPF");
  
  std::vector<std::string> absEtaRanges;
  absEtaRanges.push_back("absEtaLt1p40");
  absEtaRanges.push_back("absEta1p40to2p40");

  std::vector<std::string> ptRanges;
  ptRanges.push_back("pt20to25");
  ptRanges.push_back("pt25to30");
  ptRanges.push_back("pt30to35");
  ptRanges.push_back("pt35to40");

  std::vector<std::string> isolationWPs;
  isolationWPs.push_back("relChargedIsoLt0p20");

  std::vector<std::string> observables;
  observables.push_back("leadTrack_pt");
  observables.push_back("leadTrack_chi2");
  observables.push_back("leadTrack_numStubs");
  observables.push_back("signalTrack_pt");
  observables.push_back("signalTrack_chi2");
  observables.push_back("signalTrack_numStubs");
  observables.push_back("isolationTrack_pt");
  observables.push_back("isolationTrack_chi2");
  observables.push_back("isolationTrack_numStubs");
  observables.push_back("mTracks");
  observables.push_back("mTracks_plus_Strip");
  observables.push_back("strip_pt");
  observables.push_back("strip_numPhotons");
  observables.push_back("strip_numElectrons");
  observables.push_back("strip_numPhotons_plus_Electrons");

  std::map<std::string, int> rebin;   // key = observable
  rebin["leadTrack_pt"]                          =   1;
  rebin["leadTrack_chi2"]                        =   5;
  rebin["leadTrack_numStubs"]                    =   1;
  rebin["signalTrack_pt"]                        =   1;
  rebin["signalTrack_chi2"]                      =   5;
  rebin["signalTrack_numStubs"]                  =   1;
  rebin["isolationTrack_pt"]                     =   1;
  rebin["isolationTrack_chi2"]                   =   5;
  rebin["isolationTrack_numStubs"]               =   1;
  rebin["mTracks"]                               =   1;
  rebin["mTracks_plus_Strip"]                    =   1;
  rebin["strip_pt"]                              =   1;
  rebin["strip_numPhotons"]                      =   1;
  rebin["strip_numElectrons"]                    =   1;
  rebin["strip_numPhotons_plus_Electrons"]       =   1;

  std::map<std::string, double> xMin; // key = observable
  xMin["leadTrack_pt"]                           =   0.;
  xMin["leadTrack_chi2"]                         =   0.;
  xMin["leadTrack_numStubs"]                     =   2.5;
  xMin["signalTrack_pt"]                         =   0.;
  xMin["signalTrack_chi2"]                       =   0.;
  xMin["signalTrack_numStubs"]                   =   2.5;
  xMin["isolationTrack_pt"]                      =   0.;
  xMin["isolationTrack_chi2"]                    =   0.;
  xMin["isolationTrack_numStubs"]                =   2.5;
  xMin["mTracks"]                                =   0.;
  xMin["mTracks_plus_Strip"]                     =   0.;
  xMin["strip_pt"]                               =   0.;
  xMin["strip_numPhotons"]                       =   0.;
  xMin["strip_numElectrons"]                     =   0.;
  xMin["strip_numPhotons_plus_Electrons"]        =   0.;

  std::map<std::string, double> xMax; // key = observable
  xMax["leadTrack_pt"]                           =  50.;
  xMax["leadTrack_chi2"]                         = 100.;
  xMax["leadTrack_numStubs"]                     =  10.5;
  xMax["signalTrack_pt"]                         =  50.;
  xMax["signalTrack_chi2"]                       = 100.;
  xMax["signalTrack_numStubs"]                   =  10.5;
  xMax["isolationTrack_pt"]                      =  10.;
  xMax["isolationTrack_chi2"]                    = 100.;
  xMax["isolationTrack_numStubs"]                =  10.5;
  xMax["mTracks"]                                =  10.;
  xMax["mTracks_plus_Strip"]                     =  10.;
  xMax["strip_pt"]                               =  50.;
  xMax["strip_numPhotons"]                       =  24.5;
  xMax["strip_numElectrons"]                     =  24.5;
  xMax["strip_numPhotons_plus_Electrons"]        =  24.5;
  
  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["leadTrack_pt"]                    = "p_{T} [GeV]";
  xAxisTitles["leadTrack_chi2"]                  = "#chi^{2}";
  xAxisTitles["leadTrack_numStubs"]              = "N_{stubs}";
  xAxisTitles["signalTrack_pt"]                  = "p_{T} [GeV]";
  xAxisTitles["signalTrack_chi2"]                = "#chi^{2}";
  xAxisTitles["signalTrack_numStubs"]            = "N_{stubs}";
  xAxisTitles["isolationTrack_pt"]               = "p_{T} [GeV]";
  xAxisTitles["isolationTrack_chi2"]             = "#chi^{2}";
  xAxisTitles["isolationTrack_numStubs"]         = "N_{stubs}";
  xAxisTitles["mTracks"]                         = "m [GeV]";
  xAxisTitles["mTracks_plus_Strip"]              = "m [GeV]";
  xAxisTitles["strip_pt"]                        = "p_{T} [GeV]";
  xAxisTitles["strip_numPhotons"]                = "N_{#gamma}";
  xAxisTitles["strip_numElectrons"]              = "N_{e}";
  xAxisTitles["strip_numPhotons_plus_Electrons"] = "N_{#gamma} + N_{e}";

  std::string dqmDirectory = "DQMData/TallinnL1PFTauQualityAnalyzer";
  
  int colors[6] = { 1, 2, 8, 4, 6, 7 };
  int lineStyles[6] = { 1, 1, 1, 1, 1, 1 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  // TallinnL1PFTaus 
  for ( std::vector<std::string>::const_iterator pfAlgo = pfAlgos.begin();
	pfAlgo != pfAlgos.end(); ++pfAlgo ) {
    for ( std::vector<std::string>::const_iterator ptRange = ptRanges.begin();
	  ptRange != ptRanges.end(); ++ptRange ) {
      for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	    absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
	for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs.begin();
	      isolationWP != isolationWPs.end(); ++isolationWP ) {	    
	  for ( std::vector<std::string>::const_iterator observable = observables.begin();
		observable != observables.end(); ++observable ) {      
	    std::string histogramName = Form("%s%s/%s/%s/%s_%s_%s_%s", 
              dqmDirectory.data(), pfAlgo->data(), ptRange->data(), absEtaRange->data(), observable->data(), ptRange->data(), absEtaRange->data(), isolationWP->data());
	    TH1* histogram_signal = loadHistogram(inputFile_signal, histogramName);
	    TH1* histogram_signal_rebinned = ( rebin[*observable] > 1 ) ? histogram_signal->Rebin(rebin[*observable]) : histogram_signal;
	    TH1* histogram_background = loadHistogram(inputFile_background, histogramName);
	    TH1* histogram_background_rebinned = ( rebin[*observable] > 1 ) ? histogram_background->Rebin(rebin[*observable]) : histogram_background;
	    std::vector<std::string> labelTextLines;
  	    std::string outputFileName = Form("makeTallinnL1PFTau%sQualityPlots_%s_%s_%s_%s.png", 
	      pfAlgo->data(), observable->data(), ptRange->data(), absEtaRange->data(), isolationWP->data());
	    showHistograms(1150, 850,
			   histogram_signal_rebinned, "Signal",
			   histogram_background_rebinned, "Background",
			   0, "",	     
			   0, "",
			   0, "",
			   0, "",
			   colors, lineStyles, 
			   0.050, 0.70, 0.74, 0.23, 0.18, 
			   labelTextLines, 0.050,
			   0.71, 0.76, 0.23, 0.15, 
			   xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
			   true, 1.e-4, 1.99, "Events", 1.4, 
			   outputFileName);
	  }
	}
      }
    }
  }

  delete inputFile_signal;
  delete inputFile_background;
}

