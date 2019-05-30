
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
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

void showHistograms(double canvasSizeX, double canvasSizeY,
                    TH1* histogram1, const std::string& legendEntry1,
                    TH1* histogram2, const std::string& legendEntry2,
                    TH1* histogram3, const std::string& legendEntry3,
                    TH1* histogram4, const std::string& legendEntry4,
                    TH1* histogram5, const std::string& legendEntry5,
                    TH1* histogram6, const std::string& legendEntry6,
                    int colors[], int lineStyles[], int lineWidths[], 
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
  histogram1->SetLineStyle(lineStyles[0]);
  histogram1->SetLineWidth(lineWidths[0]);
  histogram1->Draw("hist");

  if ( histogram2 ) {
    histogram2->SetLineColor(colors[1]);
    histogram2->SetLineStyle(lineStyles[1]);
    histogram2->SetLineWidth(lineWidths[1]);
    histogram2->Draw("histsame");
  }

  if ( histogram3 ) {
    histogram3->SetLineColor(colors[2]);
    histogram3->SetLineStyle(lineStyles[2]);
    histogram3->SetLineWidth(lineWidths[2]);
    histogram3->Draw("histsame");
  }

  if ( histogram4 ) {
    histogram4->SetLineColor(colors[3]);
    histogram4->SetLineStyle(lineStyles[3]);
    histogram4->SetLineWidth(lineWidths[3]);
    histogram4->Draw("histsame");
  }

  if ( histogram5 ) {
    histogram5->SetLineColor(colors[4]);
    histogram5->SetLineStyle(lineStyles[4]);
    histogram5->SetLineWidth(lineWidths[4]);
    histogram5->Draw("histsame");
  }

  if ( histogram6 ) {
    histogram6->SetLineColor(colors[5]);
    histogram6->SetLineStyle(lineStyles[5]);
    histogram6->SetLineWidth(lineWidths[5]);
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

enum { kGtCut, kLtCut };

TGraph* compGraph_efficiency(TH1* histogram, int mode)
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
    if      ( mode == kGtCut ) graph_Efficiency->SetPoint(idxPoint, binCenter, sum/integral);
    else if ( mode == kLtCut ) graph_Efficiency->SetPoint(idxPoint, binCenter, 1.0 - (sum/integral));
    else assert(0);
  }
  return graph_Efficiency;
}

TGraph* compGraph_rocCurve(TH1* histogram_signal, TH1* histogram_background, int mode)
{
  TGraph* graph_signal = compGraph_efficiency(histogram_signal, mode);
  TGraph* graph_background = compGraph_efficiency(histogram_background, mode);
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
		int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], 
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
  graph1->SetMarkerStyle(markerStyles[0]);
  graph1->SetMarkerSize(markerSizes[0]);
  graph1->SetLineColor(colors[0]);
  graph1->SetLineStyle(lineStyles[0]);
  graph1->SetLineWidth(lineWidths[0]);
  graph1->Draw("L");

  if ( graph2 ) {
    graph2->SetMarkerColor(colors[1]);
    graph2->SetMarkerStyle(markerStyles[1]);
    graph2->SetMarkerSize(markerSizes[1]);
    graph2->SetLineColor(colors[1]);
    graph2->SetLineStyle(lineStyles[1]);
    graph2->SetLineWidth(lineWidths[1]);
    graph2->Draw("L");
  }

  if ( graph3 ) {
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->SetMarkerSize(markerSizes[2]);
    graph3->SetLineColor(colors[2]);
    graph3->SetLineStyle(lineStyles[2]);
    graph3->SetLineWidth(lineWidths[2]);
    graph3->Draw("L");
  }

  if ( graph4 ) {
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->SetMarkerSize(markerSizes[3]);
    graph4->SetLineColor(colors[3]);
    graph4->SetLineStyle(lineStyles[3]);
    graph4->SetLineWidth(lineWidths[3]);
    graph4->Draw("L");
  }

  if ( graph5 ) {
    graph5->SetMarkerColor(colors[4]);
    graph5->SetMarkerStyle(markerStyles[4]);
    graph5->SetMarkerSize(markerSizes[4]);
    graph5->SetLineColor(colors[4]);
    graph5->SetLineStyle(lineStyles[4]);
    graph5->SetLineWidth(lineWidths[4]);
    graph5->Draw("L");
  }

  if ( graph6 ) {
    graph6->SetMarkerColor(colors[5]);
    graph6->SetMarkerStyle(markerStyles[5]);
    graph6->SetMarkerSize(markerSizes[5]);
    graph6->SetLineColor(colors[5]);
    graph6->SetLineStyle(lineStyles[5]);
    graph6->SetLineWidth(lineWidths[5]);
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

void showFitFunctions(double canvasSizeX, double canvasSizeY,
		      TProfile* profile1, TF1* fitFunction1, const std::string& legendEntry1,
		      TProfile* profile2, TF1* fitFunction2, const std::string& legendEntry2,
		      TProfile* profile3, TF1* fitFunction3, const std::string& legendEntry3,
		      TProfile* profile4, TF1* fitFunction4, const std::string& legendEntry4,
		      TProfile* profile5, TF1* fitFunction5, const std::string& legendEntry5,
		      TProfile* profile6, TF1* fitFunction6, const std::string& legendEntry6,
		      int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], 
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

  if ( !(profile1 && fitFunction1) ) {
    std::cerr << "<showHistograms>: profile1 = NULL or fitFunction1 = NULL --> skipping !!" << std::endl;
    return;
  }

  profile1->SetTitle("");
  profile1->SetStats(false);
  profile1->SetMinimum(yMin);
  profile1->SetMaximum(yMax);

  TAxis* xAxis = profile1->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  
  if ( xMax > xMin ) {
    std::cout << "limiting x-axis range to " << xMin << ".." << xMax << std::endl;
    xAxis->SetRangeUser(xMin, xMax);
  }

  TAxis* yAxis = profile1->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  profile1->SetLineColor(colors[0]);
  profile1->SetLineStyle(1);
  profile1->SetLineWidth(lineWidths[0]);
  profile1->SetMarkerColor(colors[0]);
  profile1->SetMarkerStyle(markerStyles[0]);
  profile1->SetMarkerSize(markerSizes[0]);
  profile1->Draw("e1p");

  fitFunction1->SetLineColor(colors[0]);
  fitFunction1->SetLineStyle(lineStyles[0]);
  fitFunction1->SetLineWidth(lineWidths[0]);
  fitFunction1->Draw("Lsame");

  if ( profile2 && fitFunction2 ) {
    profile2->SetLineColor(colors[1]);
    profile2->SetLineStyle(1);
    profile2->SetLineWidth(lineWidths[1]);
    profile2->SetMarkerColor(colors[1]);
    profile2->SetMarkerStyle(markerStyles[1]);
    profile2->SetMarkerSize(markerSizes[1]);
    profile2->Draw("e1psame");

    fitFunction2->SetLineColor(colors[1]);
    fitFunction2->SetLineStyle(lineStyles[1]);
    fitFunction2->SetLineWidth(lineWidths[1]);
    fitFunction2->Draw("Lsame");
  }

  if ( profile3 && fitFunction3 ) {
    profile3->SetLineColor(colors[2]);
    profile3->SetLineStyle(1);
    profile3->SetLineWidth(lineWidths[2]);
    profile3->SetMarkerColor(colors[2]);
    profile3->SetMarkerStyle(markerStyles[2]);
    profile3->SetMarkerSize(markerSizes[2]);
    profile3->Draw("e1psame");

    fitFunction3->SetLineColor(colors[2]);
    fitFunction3->SetLineStyle(lineStyles[2]);
    fitFunction3->SetLineWidth(lineWidths[2]);
    fitFunction3->Draw("Lsame");
  }
  
  if ( profile4 && fitFunction4 ) {
    profile4->SetLineColor(colors[3]);
    profile4->SetLineStyle(1);
    profile4->SetLineWidth(lineWidths[3]);
    profile4->SetMarkerColor(colors[3]);
    profile4->SetMarkerStyle(markerStyles[3]);
    profile4->SetMarkerSize(markerSizes[3]);
    profile4->Draw("e1psame");

    fitFunction4->SetLineColor(colors[3]);
    fitFunction4->SetLineStyle(lineStyles[3]);
    fitFunction4->SetLineWidth(lineWidths[3]);
    fitFunction4->Draw("Lsame");
  }

  if ( profile5 && fitFunction5 ) {
    profile5->SetLineColor(colors[4]);
    profile5->SetLineStyle(1);
    profile5->SetLineWidth(lineWidths[4]);
    profile5->SetMarkerColor(colors[4]);
    profile5->SetMarkerStyle(markerStyles[4]);
    profile5->SetMarkerSize(markerSizes[4]);
    profile5->Draw("e1psame");

    fitFunction5->SetLineColor(colors[4]);
    fitFunction5->SetLineStyle(lineStyles[4]);
    fitFunction5->SetLineWidth(lineWidths[4]);
    fitFunction5->Draw("Lsame");
  }
  
  if ( profile6 && fitFunction6 ) {
    profile6->SetLineColor(colors[5]);
    profile6->SetLineStyle(1);
    profile6->SetLineWidth(lineWidths[5]);
    profile6->SetMarkerColor(colors[5]);
    profile6->SetMarkerStyle(markerStyles[5]);
    profile6->SetMarkerSize(markerSizes[5]);
    profile6->Draw("e1psame");

    fitFunction6->SetLineColor(colors[5]);
    fitFunction6->SetLineStyle(lineStyles[5]);
    fitFunction6->SetLineWidth(lineWidths[5]);
    fitFunction6->Draw("Lsame");
  }

  TLegend* legend = 0;
  if ( legendEntry1 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    legend->AddEntry(fitFunction1, legendEntry1.data(), "l");
    if ( profile2 && fitFunction2 ) legend->AddEntry(fitFunction2, legendEntry2.data(), "l");
    if ( profile3 && fitFunction3 ) legend->AddEntry(fitFunction3, legendEntry3.data(), "l");
    if ( profile4 && fitFunction4 ) legend->AddEntry(fitFunction4, legendEntry4.data(), "l");
    if ( profile5 && fitFunction5 ) legend->AddEntry(fitFunction5, legendEntry5.data(), "l");
    if ( profile6 && fitFunction6 ) legend->AddEntry(fitFunction6, legendEntry6.data(), "l");
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

  profile1->Draw("axissame");

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

double square(double x)
{
  return x*x;
}

double compFitRMS(const TH2* histogram, const TF1* fitFunction)
{
  const TAxis* xAxis = histogram->GetXaxis();
  int numBinsX = xAxis->GetNbins();
  const TAxis* yAxis = histogram->GetYaxis();
  int numBinsY = yAxis->GetNbins();
  double rms2 = 0.;
  double sum  = 0.;
  for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX )
  {
    double binCenterX   = xAxis->GetBinCenter(idxBinX);
    double fitFunctionX = binCenterX;
    double fitFunctionY = fitFunction->Eval(fitFunctionX);
    double binX_rms2 = 0.;
    double binX_sum  = 0.;
    for ( int idxBinY = 1; idxBinY <= numBinsY; ++idxBinY )
    {
      double binCenterY = yAxis->GetBinCenter(idxBinY);
      double binContent = histogram->GetBinContent(idxBinX, idxBinY);
      binX_rms2 += binContent*square(binCenterY - fitFunctionY);
      binX_sum  += binContent;
    }
    //std::cout << "binX #" << idxBinX << " (x = " << binCenterX << "): y(fitted) = " << fitFunctionY << ", rms^s = ";
    //if ( binX_sum > 0. ) std::cout << binX_rms2/binX_sum;
    //else std::cout << "N/A";
    //std::cout << std::endl;
    rms2 += binX_rms2;
    sum  += binX_sum;
  }
  if ( sum > 0. ) 
  {
    rms2 /= sum;
  }
  double rms = TMath::Sqrt(rms2);
  //std::cout << "--> returning rms = " << rms << std::endl;
  return rms;
}

void makeIsolationPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/L1Trigger/TallinnL1PFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName_signal = "TallinnL1PFTauAnalyzer_signal_2019May27v2.root";
  TFile* inputFile_signal = openFile(inputFilePath, inputFileName_signal);
  std::string inputFileName_background = "TallinnL1PFTauAnalyzer_background_2019May27v2.root";
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

  std::vector<std::string> binning_absEta;
  binning_absEta.push_back("absEta0p00to0p30");
  binning_absEta.push_back("absEta0p30to0p60");
  binning_absEta.push_back("absEta0p60to0p90");
  binning_absEta.push_back("absEta0p90to1p20");
  binning_absEta.push_back("absEta1p20to1p50");
  
  std::vector<std::string> binning_pt;
  binning_pt.push_back("pt20to30");
  binning_pt.push_back("pt30to40");
  binning_pt.push_back("pt40to60");
  binning_pt.push_back("pt60to100");

  std::vector<std::string> observables;
  observables.push_back("absChargedIso");
  observables.push_back("absNeutralIso");
  observables.push_back("absNeutralIso_wDeltaBetaCorr");
  observables.push_back("absNeutralIso_wRhoCorr");
  observables.push_back("absCombinedIso");
  observables.push_back("absCombinedIso_wDeltaBetaCorr");
  observables.push_back("absCombinedIso_wRhoCorr");
  observables.push_back("relChargedIso");
  observables.push_back("relNeutralIso");
  observables.push_back("relNeutralIso_wDeltaBetaCorr");
  observables.push_back("relNeutralIso_wRhoCorr"); 
  observables.push_back("relCombinedIso");
  observables.push_back("relCombinedIso_wDeltaBetaCorr");
  observables.push_back("relCombinedIso_wRhoCorr"); 
  observables.push_back("pt");

  std::map<std::string, int> mode; // key = observable
  mode["pt"]                                     = kGtCut;
  mode["absChargedIso"]                          = kLtCut;
  mode["absNeutralIso"]                          = kLtCut;
  mode["absNeutralIso_wDeltaBetaCorr"]           = kLtCut;
  mode["absNeutralIso_wRhoCorr"]                 = kLtCut;
  mode["absCombinedIso"]                         = kLtCut;
  mode["absCombinedIso_wDeltaBetaCorr"]          = kLtCut;
  mode["absCombinedIso_wRhoCorr"]                = kLtCut;
  mode["relChargedIso"]                          = kLtCut;
  mode["relNeutralIso"]                          = kLtCut;
  mode["relNeutralIso_wDeltaBetaCorr"]           = kLtCut;
  mode["relNeutralIso_wRhoCorr"]                 = kLtCut;
  mode["relCombinedIso"]                         = kLtCut;
  mode["relCombinedIso_wDeltaBetaCorr"]          = kLtCut;
  mode["relCombinedIso_wRhoCorr"]                = kLtCut;
  mode["sumChargedIsoPileup"]                    = kLtCut;

  std::map<std::string, int> rebin; // key = observable
  rebin["pt"]                                    =   5;
  rebin["absChargedIso"]                         =   5;
  rebin["absNeutralIso"]                         =   5;
  rebin["absNeutralIso_wDeltaBetaCorr"]          =   5;
  rebin["absNeutralIso_wRhoCorr"]                =   5;
  rebin["absCombinedIso"]                        =   5;
  rebin["absCombinedIso_wDeltaBetaCorr"]         =   5;
  rebin["absCombinedIso_wRhoCorr"]               =   5;
  rebin["relChargedIso"]                         =   2;
  rebin["relNeutralIso"]                         =   2;
  rebin["relNeutralIso_wDeltaBetaCorr"]          =   2;
  rebin["relNeutralIso_wRhoCorr"]                =   2;
  rebin["relCombinedIso"]                        =   2;
  rebin["relCombinedIso_wDeltaBetaCorr"]         =   2;
  rebin["relCombinedIso_wRhoCorr"]               =   2;
  rebin["sumChargedIsoPileup"]                   =   5;

  std::map<std::string, double> xMin; // key = observable
  xMin["pt"]                                     =   0.;
  xMin["absChargedIso"]                          =   0.;
  xMin["absNeutralIso"]                          =   0.;
  xMin["absNeutralIso_wDeltaBetaCorr"]           =   0.;
  xMin["absNeutralIso_wRhoCorr"]                 =   0.;
  xMin["absCombinedIso"]                         =   0.;
  xMin["absCombinedIso_wDeltaBetaCorr"]          =   0.;
  xMin["absCombinedIso_wRhoCorr"]                =   0.;
  xMin["relChargedIso"]                          =   0.;
  xMin["relNeutralIso"]                          =   0.;
  xMin["relNeutralIso_wDeltaBetaCorr"]           =   0.;  
  xMin["relNeutralIso_wRhoCorr"]                 =   0.; 
  xMin["relCombinedIso"]                         =   0.;
  xMin["relCombinedIso_wDeltaBetaCorr"]          =   0.;  
  xMin["relCombinedIso_wRhoCorr"]                =   0.; 
  xMin["sumChargedIsoPileup"]                    =   0.; 

  std::map<std::string, double> xMax; // key = observable
  xMax["pt"]                                     = 100.;
  xMax["absChargedIso"]                          =  25.;
  xMax["absNeutralIso"]                          =  25.;
  xMax["absNeutralIso_wDeltaBetaCorr"]           =  25.;
  xMax["absNeutralIso_wRhoCorr"]                 =  25.;
  xMax["absCombinedIso"]                         =  25.;
  xMax["absCombinedIso_wDeltaBetaCorr"]          =  25.;
  xMax["absCombinedIso_wRhoCorr"]                =  25.;
  xMax["relChargedIso"]                          =   1.;
  xMax["relNeutralIso"]                          =   1.;
  xMax["relCombinedIso"]                         =   1.;
  xMax["relCombinedIso_wDeltaBetaCorr"]          =   1.;
  xMax["relCombinedIso_wRhoCorr"]                =   1.;
  xMax["sumChargedIsoPileup"]                    =  25.;
  
  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["pt"]                              = "L1 #tau_{h} p_{T} [GeV]";
  xAxisTitles["absChargedIso"]                   = "L1 #tau I_{ch} [GeV]";
  xAxisTitles["absNeutralIso"]                   = "L1 #tau I_{neu} [GeV]";
  xAxisTitles["absNeutralIso_wDeltaBetaCorr"]    = "#Delta#beta-corrected L1 #tau I_{neu} [GeV]";
  xAxisTitles["absNeutralIso_wRhoCorr"]          = "#rho-corrected L1 #tau I_{neu} [GeV]";
  xAxisTitles["absCombinedIso"]                  = "L1 #tau I_{cmb} [GeV]";
  xAxisTitles["absCombinedIso_wDeltaBetaCorr"]   = "#Delta#beta-corrected L1 #tau I_{cmb} [GeV]";
  xAxisTitles["absCombinedIso_wRhoCorr"]         = "#rho-corrected L1 #tau I_{cmb} [GeV]";
  xAxisTitles["relChargedIso"]                   = "L1 #tau I_{ch} / p_{T}";
  xAxisTitles["relNeutralIso"]                   = "L1 #tau I_{neu} / p_{T}";
  xAxisTitles["relNeutralIso_wDeltaBetaCorr"]    = "#Delta#beta-corrected L1 #tau I_{neu} / p_{T}";
  xAxisTitles["relNeutralIso_wRhoCorr"]          = "#rho-corrected L1 #tau I_{neu} / p_{T}";
  xAxisTitles["relCombinedIso"]                  = "L1 #tau I_{cmb} / p_{T}";
  xAxisTitles["relCombinedIso_wDeltaBetaCorr"]   = "#Delta#beta-corrected L1 #tau I_{cmb} / p_{T}";
  xAxisTitles["relCombinedIso_wRhoCorr"]         = "#rho-corrected L1 #tau I_{cmb} / p_{T}";

  std::map<std::string, std::string> legendEntries; // key = observable
  legendEntries["pt"]                            = "p_{T}";
  legendEntries["absChargedIso"]                 = "I_{ch}";
  legendEntries["absNeutralIso"]                 = "I_{neu}";
  legendEntries["absNeutralIso_wDeltaBetaCorr"]  = "I_{neu} (#Delta#beta-corr.)";
  legendEntries["absNeutralIso_wRhoCorr"]        = "I_{neu} (#rho-corr.)";
  legendEntries["absCombinedIso"]                = "I_{cmb}";
  legendEntries["absCombinedIso_wDeltaBetaCorr"] = "I_{cmb} (#Delta#beta-corr.)";
  legendEntries["absCombinedIso_wRhoCorr"]       = "I_{cmb} (#rho-corr.)";
  legendEntries["relChargedIso"]                 = "I_{ch}";
  legendEntries["relNeutralIso"]                 = "I_{neu}";
  legendEntries["relNeutralIso_wDeltaBetaCorr"]  = "I_{neu}(#Delta#beta-corr.)";
  legendEntries["relNeutralIso_wRhoCorr"]        = "I_{neu} (#rho-corr.)";
  legendEntries["relCombinedIso"]                = "I_{cmb}";
  legendEntries["relCombinedIso_wDeltaBetaCorr"] = "I_{cmb}(#Delta#beta-corr.)";
  legendEntries["relCombinedIso_wRhoCorr"]       = "I_{cmb} (#rho-corr.)";

  std::string dqmDirectory = "DQMData/TallinnL1PFTauIsolationAnalyzerWithStripsWithoutPreselection";
  
  int colors[6]                 = {  1,  2,  8,  4,  6,  7 };
  int lineStyles[6]             = {  1,  1,  1,  1,  1,  1 };
  int lineWidths_histogram[6]   = {  2,  2,  2,  2,  2,  2 };
  int lineWidths_graph[6]       = {  3,  3,  3,  3,  3,  3 };
  int lineWidths_fitFunction[6] = {  2,  2,  2,  2,  2,  2 };
  int markerStyles[6]           = { 22, 32, 20, 24, 21, 25 };
  int markerSizes[6]            = {  2,  2,  2,  2,  2,  2 };

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
	  TH1* histogram_signal_rebinned = ( rebin[*observable] > 1 ) ? 
	    histogram_signal->Rebin(rebin[*observable]) : histogram_signal;

	  TH1* histogram_background = loadHistogram(inputFile_background, histogramName);
	  histograms_background[*observable] = histogram_background;
	  TH1* histogram_background_rebinned = ( rebin[*observable] > 1 ) ? 
	    histogram_background->Rebin(rebin[*observable]) : histogram_signal;

	  TGraph* graph_roc = compGraph_rocCurve(histogram_signal, histogram_background, mode[*observable]);
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
			 colors, lineStyles, lineWidths_histogram,
			 0.050, 0.68, 0.74, 0.23, 0.17, 
			 labelTextLines, 0.050,
			 0.70, 0.62, 0.23, 0.06, 
			 xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
			 true, 1.e-3, 1.99e0, "Events", 1.4, 
			 outputFileName_distribution);
	}

	std::vector<std::string> labelTextLines = getLabelTextLines(*ptThreshold);
	std::string outputFileName_roc_absIso = Form("makeIsolationPlots_rocCurves_absIso_%s_%s_%s.png", 
	  pfAlgo->data(), absEtaRange->data(), ptThreshold->data());

	showGraphs(1150, 1150,
		   graphs_roc["absChargedIso"],                 legendEntries["absChargedIso"],
		   graphs_roc["absNeutralIso"],                 legendEntries["absNeutralIso"],
		   graphs_roc["absCombinedIso"],                legendEntries["absCombinedIso"],
		   graphs_roc["absCombinedIso_wDeltaBetaCorr"], legendEntries["absCombinedIso_wDeltaBetaCorr"],
		   graphs_roc["absCombinedIso_wRhoCorr"],       legendEntries["absCombinedIso_wRhoCorr"],
		   graphs_roc["pt"],                            legendEntries["pt"],
		   colors, markerStyles, markerSizes, lineStyles, lineWidths_graph, 
		   0.040, 0.67, 0.17, 0.23, 0.27, 
		   labelTextLines, 0.045,
		   0.18, 0.86, 0.26, 0.05, 
		   0., 1.09, "Signal Efficiency", 1.2, 
		   false, 0., 1.09, "Background Rejection", 1.4, 
		   outputFileName_roc_absIso);
	std::string outputFileName_roc_relIso = Form("makeIsolationPlots_rocCurves_relIso_%s_%s_%s.png", 
	  pfAlgo->data(), absEtaRange->data(), ptThreshold->data());
	showGraphs(1150, 1150,
		   graphs_roc["relChargedIso"],                 legendEntries["relChargedIso"],
		   graphs_roc["relNeutralIso"],                 legendEntries["relNeutralIso"],
		   graphs_roc["relCombinedIso"],                legendEntries["relCombinedIso"],
		   graphs_roc["relCombinedIso_wDeltaBetaCorr"], legendEntries["relCombinedIso_wDeltaBetaCorr"],
		   graphs_roc["relCombinedIso_wRhoCorr"],       legendEntries["relCombinedIso_wRhoCorr"],
		   graphs_roc["pt"],                            legendEntries["pt"],
		   colors, markerStyles, markerSizes, lineStyles, lineWidths_graph, 
		   0.040, 0.67, 0.17, 0.23, 0.27, 
		   labelTextLines, 0.045,
		   0.18, 0.86, 0.26, 0.05, 
		   0., 1.09, "Signal Efficiency", 1.2, 
		   false, 0., 1.09, "Background Rejection", 1.4, 
		   outputFileName_roc_relIso);	
      }
    }
  }

  for ( std::vector<std::string>::const_iterator pfAlgo = pfAlgos.begin();
	pfAlgo != pfAlgos.end(); ++pfAlgo ) {
    for ( std::vector<std::string>::const_iterator bin_absEta = binning_absEta.begin();
	  bin_absEta != binning_absEta.end(); ++bin_absEta ) {
      for ( std::vector<std::string>::const_iterator bin_pt = binning_pt.begin();
	    bin_pt != binning_pt.end(); ++bin_pt ) {
	std::string histogramName_absNeutralIso = Form("%s%s/absNeutralIso_all_%s_%s", 
	  dqmDirectory.data(), pfAlgo->data(), bin_absEta->data(), bin_pt->data());
	TH1* histogram_absNeutralIso = loadHistogram(inputFile_signal, histogramName_absNeutralIso);
	TH1* histogram_absNeutralIso_rebinned = ( rebin["absNeutralIso"] > 1 ) ? 
	  histogram_absNeutralIso->Rebin(rebin["absNeutralIso"]) : histogram_absNeutralIso;

	std::string histogramName_sumChargedIsoPileup = Form("%s%s/sumChargedIsoPileup_all_%s_%s", 
	  dqmDirectory.data(), pfAlgo->data(), bin_absEta->data(), bin_pt->data());
	TH1* histogram_sumChargedIsoPileup = loadHistogram(inputFile_signal, histogramName_sumChargedIsoPileup);
	TH1* histogram_sumChargedIsoPileup_rebinned = ( rebin["sumChargedIsoPileup"] > 1 ) ? 
	  histogram_sumChargedIsoPileup->Rebin(rebin["sumChargedIsoPileup"]) : histogram_sumChargedIsoPileup;

	std::vector<std::string> labelTextLines;
	labelTextLines.push_back(Form("%s & %s", bin_pt->data(), bin_absEta->data()));
	labelTextLines.push_back(Form(" I_{neu}: mean = %1.2f, rms = %1.2f", histogram_absNeutralIso->GetMean(), histogram_absNeutralIso->GetRMS()));
	labelTextLines.push_back(Form(" I_{ch}^{PU}: mean = %1.2f, rms = %1.2f", histogram_sumChargedIsoPileup->GetMean(), histogram_sumChargedIsoPileup->GetRMS()));
	std::string outputFileName = Form("makeIsolationPlots_absNeutralIso_vs_sumChargedIsoPileup_%s_%s_%s.png", 
	  pfAlgo->data(), bin_absEta->data(), bin_pt->data());
	showHistograms(1150, 850,
		       histogram_absNeutralIso_rebinned,       "I_{neu}",
		       histogram_sumChargedIsoPileup_rebinned, "I_{ch}^{PU}",
		       0, "",
		       0, "",
		       0, "",
		       0, "",
		       colors, lineStyles, lineWidths_histogram, 
		       0.050, 0.80, 0.74, 0.08, 0.17, 
		       labelTextLines, 0.040,
		       0.18, 0.75, 0.33, 0.17, 
		       0., 25., "L1 #tau Isolation [GeV]", 1.2, 
		       true, 1.e-3, 1.99e0, "Events", 1.4, 
		       outputFileName);
      }
    }
  }

  std::cout << "correlations:" << std::endl;
  for ( std::vector<std::string>::const_iterator pfAlgo = pfAlgos.begin();
	pfAlgo != pfAlgos.end(); ++pfAlgo ) {
    for ( std::vector<std::string>::const_iterator bin_absEta = binning_absEta.begin();
	  bin_absEta != binning_absEta.end(); ++bin_absEta ) {
      for ( std::vector<std::string>::const_iterator bin_pt = binning_pt.begin();
	    bin_pt != binning_pt.end(); ++bin_pt ) {
	std::string histogramName_sumNeutralIso_vs_sumChargedIsoPileup = Form("%s%s/sumNeutralIso_vs_sumChargedIsoPileup_all_%s_%s", 
	  dqmDirectory.data(), pfAlgo->data(), bin_absEta->data(), bin_pt->data());
	TH2* histogram_sumNeutralIso_vs_sumChargedIsoPileup = loadHistogram2d(inputFile_signal, histogramName_sumNeutralIso_vs_sumChargedIsoPileup);
	double deltaBeta_correlation = histogram_sumNeutralIso_vs_sumChargedIsoPileup->GetCorrelationFactor();
	
	std::string histogramName_sumNeutralIso_vs_rhoCorr = Form("%s%s/sumNeutralIso_vs_rhoCorr_all_%s_%s", 
	  dqmDirectory.data(), pfAlgo->data(), bin_absEta->data(), bin_pt->data());
	TH2* histogram_sumNeutralIso_vs_rhoCorr = loadHistogram2d(inputFile_signal, histogramName_sumNeutralIso_vs_rhoCorr);
	double rho_correlation = histogram_sumNeutralIso_vs_rhoCorr->GetCorrelationFactor();

	std::cout << (*bin_pt) << " & " << bin_absEta->data() << ":" 
		  << " deltaBeta = " << deltaBeta_correlation << ","	  
		  << " rho = " << rho_correlation << std::endl;
      }
    }
  }

  std::cout << "linear fits:" << std::endl;
  for ( std::vector<std::string>::const_iterator pfAlgo = pfAlgos.begin();
	pfAlgo != pfAlgos.end(); ++pfAlgo ) {
    for ( std::vector<std::string>::const_iterator bin_absEta = binning_absEta.begin();
	  bin_absEta != binning_absEta.end(); ++bin_absEta ) {
      for ( std::vector<std::string>::const_iterator bin_pt = binning_pt.begin();
	    bin_pt != binning_pt.end(); ++bin_pt ) {
	std::string histogramName_sumNeutralIso_vs_sumChargedIsoPileup = Form("%s%s/sumNeutralIso_vs_sumChargedIsoPileup_all_%s_%s", 
	  dqmDirectory.data(), pfAlgo->data(), bin_absEta->data(), bin_pt->data());
	TH2* histogram_sumNeutralIso_vs_sumChargedIsoPileup = loadHistogram2d(inputFile_signal, histogramName_sumNeutralIso_vs_sumChargedIsoPileup);
	std::string profileName_sumNeutralIso_vs_sumChargedIsoPileup = Form("sumNeutralIso_vs_sumChargedIsoPileup_all_%s_%s_profile_%s", 
	  bin_absEta->data(), bin_pt->data(), pfAlgo->data());
	TProfile* profile_sumNeutralIso_vs_sumChargedIsoPileup = histogram_sumNeutralIso_vs_sumChargedIsoPileup->ProfileX(
	  histogramName_sumNeutralIso_vs_sumChargedIsoPileup.data(), 1, histogram_sumNeutralIso_vs_sumChargedIsoPileup->GetNbinsY());
        std::string fitFunctionName_sumNeutralIso_vs_sumChargedIsoPileup = Form("fitFunction_sumNeutralIso_vs_sumChargedIsoPileup_all_%s_%s_%s", 
          bin_absEta->data(), bin_pt->data(), pfAlgo->data());
	TF1* fitFunction_sumNeutralIso_vs_sumChargedIsoPileup = new TF1(fitFunctionName_sumNeutralIso_vs_sumChargedIsoPileup.data(), "[0] + [1]*x", 0., 25.);
	fitFunction_sumNeutralIso_vs_sumChargedIsoPileup->SetParameter(0, 0.);
	fitFunction_sumNeutralIso_vs_sumChargedIsoPileup->SetParameter(1, 0.5);
	profile_sumNeutralIso_vs_sumChargedIsoPileup->Fit(fitFunction_sumNeutralIso_vs_sumChargedIsoPileup);
	double deltaBeta_slope  = fitFunction_sumNeutralIso_vs_sumChargedIsoPileup->GetParameter(1);
	double deltaBeta_offset = fitFunction_sumNeutralIso_vs_sumChargedIsoPileup->GetParameter(0);
	double deltaBeta_rms    = compFitRMS(histogram_sumNeutralIso_vs_sumChargedIsoPileup, fitFunction_sumNeutralIso_vs_sumChargedIsoPileup);

	std::string histogramName_sumNeutralIso_vs_rhoCorr = Form("%s%s/sumNeutralIso_vs_rhoCorr_all_%s_%s", 
	  dqmDirectory.data(), pfAlgo->data(), bin_absEta->data(), bin_pt->data());
	TH2* histogram_sumNeutralIso_vs_rhoCorr = loadHistogram2d(inputFile_signal, histogramName_sumNeutralIso_vs_rhoCorr);
        std::string profileName_sumNeutralIso_vs_rhoCorr = Form("sumNeutralIso_vs_rhoCorr_all_%s_%s_profile_%s", 
          bin_absEta->data(), bin_pt->data(), pfAlgo->data());
	TProfile* profile_sumNeutralIso_vs_rhoCorr = histogram_sumNeutralIso_vs_rhoCorr->ProfileX(
	  histogramName_sumNeutralIso_vs_rhoCorr.data(), 1, histogram_sumNeutralIso_vs_rhoCorr->GetNbinsY());
	std::string fitFunctionName_sumNeutralIso_vs_rhoCorr = Form("fitFunction_sumNeutralIso_vs_rhoCorr_all_%s_%s_%s", 
          bin_absEta->data(), bin_pt->data(), pfAlgo->data());
	TF1* fitFunction_sumNeutralIso_vs_rhoCorr = new TF1(fitFunctionName_sumNeutralIso_vs_rhoCorr.data(), "[0] + [1]*x", 0., 25.);
	fitFunction_sumNeutralIso_vs_rhoCorr->SetParameter(0, 0.);
	fitFunction_sumNeutralIso_vs_rhoCorr->SetParameter(1, 1.0);
	profile_sumNeutralIso_vs_rhoCorr->Fit(fitFunction_sumNeutralIso_vs_rhoCorr);
	double rho_slope  = fitFunction_sumNeutralIso_vs_rhoCorr->GetParameter(1);
	double rho_offset = fitFunction_sumNeutralIso_vs_rhoCorr->GetParameter(0);
	double rho_rms    = compFitRMS(histogram_sumNeutralIso_vs_rhoCorr, fitFunction_sumNeutralIso_vs_rhoCorr);
	
	std::vector<std::string> labelTextLines;
	labelTextLines.push_back(Form("%s & %s", bin_pt->data(), bin_absEta->data()));
	labelTextLines.push_back(Form(" #Delta#beta: slope = %1.2f, offset = %1.2f (rms of fit = %1.2f)", deltaBeta_slope, deltaBeta_offset, deltaBeta_rms));
	labelTextLines.push_back(Form(" #rho: slope = %1.2f, offset = %1.2f (rms of fit = %1.2f)", rho_slope, rho_offset, rho_rms));
	std::string outputFileName = Form("makeIsolationPlots_pileupCorrFits_%s_%s_%s.png", 
	  pfAlgo->data(), bin_absEta->data(), bin_pt->data());
	showFitFunctions(1150, 850,
			 profile_sumNeutralIso_vs_sumChargedIsoPileup, fitFunction_sumNeutralIso_vs_sumChargedIsoPileup, "#Delta#beta",
			 profile_sumNeutralIso_vs_rhoCorr,             fitFunction_sumNeutralIso_vs_rhoCorr,             "#rho",
			 0, 0, "",
			 0, 0, "",
			 0, 0, "",
			 0, 0, "",
			 colors, markerStyles, markerSizes, lineStyles, lineWidths_fitFunction, 
			 0.050, 0.80, 0.74, 0.08, 0.17, 
			 labelTextLines, 0.030,
			 0.18, 0.75, 0.43, 0.17, 
			 0., 25., "I_{ch}^{PU} or #rho^{corr} [GeV]", 1.2, 
			 false, 0., 5.e+1, "I_{neu}", 1.4, 
			 outputFileName);

	std::cout << (*bin_pt) << " & " << bin_absEta->data() << std::endl;
	std::cout << " deltaBeta: slope = " << deltaBeta_slope << ", offset = " << deltaBeta_offset << " (rms of fit = " << deltaBeta_rms << ")" << std::endl;
	std::cout << " rho: slope = " << rho_slope << ", offset = " << rho_offset << " (rms of fit = " << rho_rms << ")" << std::endl;

	delete fitFunction_sumNeutralIso_vs_sumChargedIsoPileup;
	delete profile_sumNeutralIso_vs_sumChargedIsoPileup;
        delete fitFunction_sumNeutralIso_vs_rhoCorr;
	delete profile_sumNeutralIso_vs_rhoCorr;
      }
    }
  }

  delete inputFile_signal;
  delete inputFile_background;
}

