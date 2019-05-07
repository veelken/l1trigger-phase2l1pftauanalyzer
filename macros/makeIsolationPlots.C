
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
  canvas->SetBottomMargin(0.12);
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

TGraph* compGraph_efficiency(TH1* histogram)
{
  TAxis* xAxis = histogram->GetXaxis();
  int numPoints = xAxis->GetNbins() + 1;
  TGraph* graph_Efficiency = new TGraph(numPoints);
  graph_Efficiency->SetPoint(0, 0., 1.0);
  double integral = histogram->Integral(1, histogram->GetNbinsX());
  double sum = 0.;
  for ( int idxPoint = 0; idxPoint <= numPoints; ++idxPoint ) {
    double binCenter = 0.;
    if ( idxPoint >= 1 ) {      
      int idxBin = idxPoint;
      binCenter = xAxis->GetBinCenter(idxBin);
      double binContent = histogram->GetBinContent(idxBin);
      sum += binContent;
    }
    std::string histogramName = histogram->GetName();
    if ( histogramName.find("nllKinFit") != std::string::npos ) graph_Efficiency->SetPoint(idxPoint, binCenter, sum/integral);
    else graph_Efficiency->SetPoint(idxPoint, binCenter, 1.0 - (sum/integral));
  }
  return graph_Efficiency;
}

TGraph* compGraph_rocCurve(TH1* histogram_signal, TH1* histogram_background)
{
  TGraph* graph_signal = compGraph_efficiency(histogram_signal);
  TGraph* graph_background = compGraph_efficiency(histogram_background);
  int numPoints = graph_signal->GetN();
  assert(graph_background->GetN() == numPoints);
  TGraph* graph_ROCcurve = new TGraph(numPoints);
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x_signal, y_signal;
    graph_signal->GetPoint(idxPoint, x_signal, y_signal);
    //double y_background = graph_background->Eval(x_signal);
    double x_background, y_background;
    graph_background->GetPoint(idxPoint, x_background, y_background);
    assert(TMath::Abs(x_signal - x_background) < 1.e-4);
    double x_ROCcurve = y_signal;
    //double y_ROCcurve = 1.0 - y_background;
    double y_ROCcurve = y_background;
    graph_ROCcurve->SetPoint(idxPoint, x_ROCcurve, y_ROCcurve);
  }
  return graph_ROCcurve;
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraph* graph1, const std::string& legendEntry1,
		TGraph* graph2, const std::string& legendEntry2,
		TGraph* graph3, const std::string& legendEntry3,
		TGraph* graph4, const std::string& legendEntry4,
		TGraph* graph5, const std::string& legendEntry5,
		TGraph* graph6, const std::string& legendEntry6,
		int colors[], int markerStyles[], int lineStyles[], 
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
  canvas->SetBottomMargin(0.12);
  canvas->SetRightMargin(0.05);

  canvas->SetLogy(useLogScale);
  
  canvas->SetGridx(1);
  canvas->SetGridy(1);

  if ( !graph1 ) {
    std::cerr << "<showGraphs>: histogram1 = NULL --> skipping !!" << std::endl;
    return;
  }

  TH1* dummyHistogram = new TH1F("dummyHistogram", "dummyHistogram", 10, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  
  if ( xMax > xMin ) {
    std::cout << "limiting x-axis range to " << xMin << ".." << xMax << std::endl;
    xAxis->SetRangeUser(xMin, xMax);
  }

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  dummyHistogram->Draw("axis");

  graph1->SetMarkerColor(colors[0]);
  graph1->SetMarkerSize(2);
  graph1->SetMarkerStyle(markerStyles[0]);
  graph1->SetLineColor(colors[0]);
  graph1->SetLineWidth(2);
  graph1->SetLineStyle(lineStyles[0]);
  graph1->Draw("L");

  if ( graph2 ) {
    graph2->SetMarkerColor(colors[1]);
    graph2->SetMarkerSize(2);
    graph2->SetMarkerStyle(markerStyles[1]);
    graph2->SetLineColor(colors[1]);
    graph2->SetLineWidth(2);
    graph2->SetLineStyle(lineStyles[1]);
    graph2->Draw("L");
  }

  if ( graph3 ) {
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerSize(2);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->SetLineColor(colors[2]);
    graph3->SetLineWidth(2);
    graph3->SetLineStyle(lineStyles[2]);
    graph3->Draw("L");
  }

  if ( graph4 ) {
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerSize(2);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->SetLineColor(colors[3]);
    graph4->SetLineWidth(2);
    graph4->SetLineStyle(lineStyles[3]);
    graph4->Draw("L");
  }

  if ( graph5 ) {
    graph5->SetMarkerColor(colors[4]);
    graph5->SetMarkerSize(2);
    graph5->SetMarkerStyle(markerStyles[4]);
    graph5->SetLineColor(colors[4]);
    graph5->SetLineWidth(2);
    graph5->SetLineStyle(lineStyles[4]);
    graph5->Draw("L");
  }

  if ( graph6 ) {
    graph6->SetMarkerColor(colors[5]);
    graph6->SetMarkerSize(2);
    graph6->SetMarkerStyle(markerStyles[5]);
    graph6->SetLineColor(colors[5]);
    graph6->SetLineWidth(2);
    graph6->SetLineStyle(lineStyles[5]);
    graph6->Draw("L");
  }

  TLegend* legend = 0;
  if ( legendEntry1 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    legend->AddEntry(graph1, legendEntry1.data(), "l");
    if ( graph2 ) legend->AddEntry(graph2, legendEntry2.data(), "l");
    if ( graph3 ) legend->AddEntry(graph3, legendEntry3.data(), "l");
    if ( graph4 ) legend->AddEntry(graph4, legendEntry4.data(), "l");
    if ( graph5 ) legend->AddEntry(graph5, legendEntry5.data(), "l");
    if ( graph6 ) legend->AddEntry(graph6, legendEntry6.data(), "l");
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

  dummyHistogram->Draw("axissame");

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
  delete dummyHistogram;
}

std::vector<std::string> getLabelTextLines(const std::string& ptThreshold)
{
  std::vector<std::string> labelTextLines;
  if ( ptThreshold == "ptGt20" ) labelTextLines.push_back("p_{T} > 20 GeV");
  if ( ptThreshold == "ptGt30" ) labelTextLines.push_back("p_{T} > 30 GeV");
  if ( ptThreshold == "ptGt40" ) labelTextLines.push_back("p_{T} > 40 GeV");
  else assert(0);
  return labelTextLines;
}

void makeIsolationPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/L1Trigger/TallinnL1PFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName_signal = "TallinnL1PFTauAnalyzer_signal_2019May02.root";
  TFile* inputFile_signal = openFile(inputFilePath, inputFileName_signal);
  std::string inputFileName_background = "TallinnL1PFTauAnalyzer_background_2019May02.root";
  TFile* inputFile_background = openFile(inputFilePath, inputFileName_background);

  std::vector<std::string> pfAlgos;
  pfAlgos.push_back("PF");
  pfAlgos.push_back("Puppi");

  std::vector<std::string> absEtaRanges;
  absEtaRanges.push_back("absEtaLt1p00");
  absEtaRanges.push_back("absEtaLt1p40");

  std::vector<std::string> ptThresholds;
  ptThresholds.push_back("ptGt20");
  ptThresholds.push_back("ptGt30");
  ptThresholds.push_back("ptGt40");

  std::vector<std::string> observables;
  observables.push_back("absChargedIso");
  observables.push_back("absNeutralIso");
  observables.push_back("absCombinedIso");
  observables.push_back("relChargedIso");
  observables.push_back("relNeutralIso");
  observables.push_back("relCombinedIso");
  observables.push_back("pt");

  std::map<std::string, int> rebin; // key = observable
  rebin["pt"]             =   5;
  rebin["absChargedIso"]  =   5;
  rebin["absNeutralIso"]  =   5;
  rebin["absCombinedIso"] =   5;
  rebin["relChargedIso"]  =   2;
  rebin["relNeutralIso"]  =   2;
  rebin["relCombinedIso"] =   2;

  std::map<std::string, double> xMin; // key = observable
  xMin["pt"]             =   0.;
  xMin["absChargedIso"]  =   0.;
  xMin["absNeutralIso"]  =   0.;
  xMin["absCombinedIso"] =   0.;
  xMin["relChargedIso"]  =   0.;
  xMin["relNeutralIso"]  =   0.;
  xMin["relCombinedIso"] =   0.;

  std::map<std::string, double> xMax; // key = observable
  xMax["pt"]             = 100.;
  xMax["absChargedIso"]  =  25.;
  xMax["absNeutralIso"]  =  25.;
  xMax["absCombinedIso"] =  25.;
  xMax["relChargedIso"]  =   1.;
  xMax["relNeutralIso"]  =   1.;
  xMax["relCombinedIso"] =   1.;
  
  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["pt"]             = "L1 #tau_{h} p_{T} [GeV]";
  xAxisTitles["absChargedIso"]  = "L1 #tau I_{ch} [GeV]";
  xAxisTitles["absNeutralIso"]  = "L1 #tau I_{neu} [GeV]";
  xAxisTitles["absCombinedIso"] = "L1 #tau I_{cmb} [GeV]";
  xAxisTitles["relChargedIso"]  = "L1 #tau I_{ch} / p_{T}";
  xAxisTitles["relNeutralIso"]  = "L1 #tau I_{neu} / p_{T}";
  xAxisTitles["relCombinedIso"] = "L1 #tau I_{cmb} / p_{T}";

  std::map<std::string, std::string> legendEntries; // key = observable
  legendEntries["pt"]             = "p_{T}";
  legendEntries["absChargedIso"]  = "I_{ch}";
  legendEntries["absNeutralIso"]  = "I_{neu}";
  legendEntries["absCombinedIso"] = "I_{cmb}";
  legendEntries["relChargedIso"]  = "I_{ch}";
  legendEntries["relNeutralIso"]  = "I_{neu}";
  legendEntries["relCombinedIso"] = "I_{cmb}";

  std::string dqmDirectory = "DQMData/TallinnL1PFTauIsolationAnalyzer";
  
  int colors[6] = { 1, 2, 8, 4, 6, 7 };
  int lineStyles[6] = { 1, 1, 1, 1, 1, 1 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  for ( std::vector<std::string>::const_iterator pfAlgo = pfAlgos.begin();
	pfAlgo != pfAlgos.end(); ++pfAlgo ) {
    for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	  absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
      for ( std::vector<std::string>::const_iterator ptThreshold = ptThresholds.begin();
	  ptThreshold != ptThresholds.end(); ++ptThreshold ) {
        std::map<std::string, TH1*> histograms_signal;     // key = observable
        std::map<std::string, TH1*> histograms_background; // key = observable
	std::map<std::string, TGraph*> graphs_roc;         // key = observable
        for ( std::vector<std::string>::const_iterator observable = observables.begin();
	      observable != observables.end(); ++observable ) {
          std::string histogramName = Form("%s%s/%s_all_%s_%s", 
            dqmDirectory.data(), pfAlgo->data(), observable->data(), absEtaRange->data(), ptThreshold->data());
	  TH1* histogram_signal = loadHistogram(inputFile_signal, histogramName);
	  histograms_signal[*observable] = histogram_signal;
	  TH1* histogram_signal_rebinned = ( rebin[*observable] > 1 ) ? histogram_signal->Rebin(rebin[*observable]) : histogram_signal;

	  TH1* histogram_background = loadHistogram(inputFile_background, histogramName);
	  histograms_background[*observable] = histogram_background;
	  TH1* histogram_background_rebinned = ( rebin[*observable] > 1 ) ? histogram_background->Rebin(rebin[*observable]) : histogram_signal;

	  TGraph* graph_roc = compGraph_rocCurve(histogram_signal, histogram_background);
	  graphs_roc[*observable] = graph_roc;

	  std::vector<std::string> labelTextLines = getLabelTextLines(*ptThreshold);
          std::string outputFileName_distribution = Form("makeIsolationPlots_%s_%s_%s_%s.png", 
	    observable->data(), pfAlgo->data(), absEtaRange->data(), ptThreshold->data());
          showHistograms(1150, 850,
			 histogram_signal_rebinned,     "Signal",
			 histogram_background_rebinned, "Background",
			 0, "",
			 0, "",
			 0, "",
			 0, "",
			 colors, lineStyles, 
			 0.050, 0.66, 0.74, 0.23, 0.15, 
			 labelTextLines, 0.050,
			 0.70, 0.62, 0.23, 0.06, 
			 xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
			 false, 0., 1.09, "Events", 1.4, 
			 outputFileName_distribution);
	}

	std::vector<std::string> labelTextLines = getLabelTextLines(*ptThreshold);
	std::string outputFileName_roc = Form("makeIsolationPlots_rocCurves_%s_%s_%s.png", 
	  pfAlgo->data(), absEtaRange->data(), ptThreshold->data());
	showGraphs(1150, 1150,
		   graphs_roc["absChargedIso"],  legendEntries["absChargedIso"],
		   graphs_roc["absNeutralIso"],  legendEntries["absNeutralIso"],
		   graphs_roc["absCombinedIso"], legendEntries["absCombinedIso"],
		   graphs_roc["relChargedIso"],  legendEntries["relChargedIso"],
		   graphs_roc["relNeutralIso"],  legendEntries["relNeutralIso"],
		   graphs_roc["relCombinedIso"], legendEntries["relCombinedIso"],
		   colors, markerStyles, lineStyles, 
		   0.045, 0.18, 0.17, 0.23, 0.26, 
		   labelTextLines, 0.050,
		   0.63, 0.65, 0.26, 0.07, 
		   0., 1.09, "Signal Efficiency", 1.2, 
		   false, 0., 1.09, "Background Rejection", 1.4, 
		   outputFileName_roc);
      }
    }
  }

  delete inputFile_signal;
  delete inputFile_background;
}

