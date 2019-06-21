
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
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

TGraph* makeGraph(const std::vector<double>& graph_pointsX, const std::vector<double>& graph_pointsY)
{
  assert(graph_pointsX.size() == graph_pointsY.size());
  size_t numPoints = graph_pointsX.size();
  TGraph* graph = new TGraph(numPoints);
  for ( size_t idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    graph->SetPoint(idxPoint, graph_pointsX[idxPoint]*1.e-3, graph_pointsY[idxPoint]); // CV: convert rate to units of kHz
  }
  return graph;
}

std::pair<TGraph*, TGraph*> makeEfficiency_and_RateGraphs(const TH2* histogram2d_signal, const TH2* histogram2d_background, 
							  double rate_min, double rate_max, double rate_step, double rate_accCorrFactor)
{
  assert(histogram2d_signal->GetNbinsX() == histogram2d_background->GetNbinsX());
  const TAxis* xAxis = histogram2d_signal->GetXaxis();
  int numBinsX = xAxis->GetNbins();
  assert(histogram2d_signal->GetNbinsY() == histogram2d_background->GetNbinsY());
  const TAxis* yAxis = histogram2d_signal->GetYaxis();
  int numBinsY = yAxis->GetNbins();

  std::vector<double> graph_symmetric_pointsX; 
  std::vector<double> graph_symmetric_pointsY; 
  std::vector<double> graph_asymmetric_pointsX; 
  std::vector<double> graph_asymmetric_pointsY; 
  
  for ( double rate_target = rate_min; rate_target <= rate_max; rate_target += rate_step ) {

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

    graph_symmetric_pointsX.push_back(best_rate_symmetric);
    graph_symmetric_pointsY.push_back(best_efficiency_symmetric);
    graph_asymmetric_pointsX.push_back(best_rate_asymmetric);
    graph_asymmetric_pointsY.push_back(best_efficiency_asymmetric);	  
  }

  TGraph* graph_symmetric = makeGraph(graph_symmetric_pointsX, graph_symmetric_pointsY);
  TGraph* graph_asymmetric = makeGraph(graph_asymmetric_pointsX, graph_asymmetric_pointsY);
  return std::pair<TGraph*, TGraph*>(graph_symmetric, graph_asymmetric);
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraph* graph1, const std::string& legendEntry1,
		TGraph* graph2, const std::string& legendEntry2,
		TGraph* graph3, const std::string& legendEntry3,
		TGraph* graph4, const std::string& legendEntry4,
		TGraph* graph5, const std::string& legendEntry5,
		TGraph* graph6, const std::string& legendEntry6,
		TGraph* graph7, const std::string& legendEntry7,
		TGraph* graph8, const std::string& legendEntry8,
		bool addFitFunctions,
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
    std::cerr << "<showGraphs>: graph1 = NULL --> skipping !!" << std::endl;
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
  graph1->Draw("P");

  if ( graph2 ) {
    graph2->SetMarkerColor(colors[1]);
    graph2->SetMarkerSize(2);
    graph2->SetMarkerStyle(markerStyles[1]);
    graph2->SetLineColor(colors[1]);
    graph2->SetLineWidth(2);
    graph2->SetLineStyle(lineStyles[1]);
    graph2->Draw("P");
  }

  if ( graph3 ) {
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerSize(2);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->SetLineColor(colors[2]);
    graph3->SetLineWidth(2);
    graph3->SetLineStyle(lineStyles[2]);
    graph3->Draw("P");
  }

  if ( graph4 ) {
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerSize(2);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->SetLineColor(colors[3]);
    graph4->SetLineWidth(2);
    graph4->SetLineStyle(lineStyles[3]);
    graph4->Draw("P");
  }

  if ( graph5 ) {
    graph5->SetMarkerColor(colors[4]);
    graph5->SetMarkerSize(2);
    graph5->SetMarkerStyle(markerStyles[4]);
    graph5->SetLineColor(colors[4]);
    graph5->SetLineWidth(2);
    graph5->SetLineStyle(lineStyles[4]);
    graph5->Draw("P");
  }

  if ( graph6 ) {
    graph6->SetMarkerColor(colors[5]);
    graph6->SetMarkerSize(2);
    graph6->SetMarkerStyle(markerStyles[5]);
    graph6->SetLineColor(colors[5]);
    graph6->SetLineWidth(2);
    graph6->SetLineStyle(lineStyles[5]);
    graph6->Draw("P");
  }

  if ( graph7 ) {
    graph7->SetMarkerColor(colors[6]);
    graph7->SetMarkerSize(2);
    graph7->SetMarkerStyle(markerStyles[6]);
    graph7->SetLineColor(colors[6]);
    graph7->SetLineWidth(2);
    graph7->SetLineStyle(lineStyles[6]);
    graph7->Draw("P");
  }

  if ( graph8 ) {
    graph8->SetMarkerColor(colors[7]);
    graph8->SetMarkerSize(2);
    graph8->SetMarkerStyle(markerStyles[7]);
    graph8->SetLineColor(colors[7]);
    graph8->SetLineWidth(2);
    graph8->SetLineStyle(lineStyles[7]);
    graph8->Draw("P");
  }

  TLegend* legend = 0;
  if ( legendEntry1 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    legend->SetNColumns(2);
    legend->AddEntry(graph1, legendEntry1.data(), "p");
    if ( graph5 && legendEntry5 != "" ) legend->AddEntry(graph5, legendEntry5.data(), "p");
    if ( graph2                       ) legend->AddEntry(graph2, legendEntry2.data(), "p");
    if ( graph6 && legendEntry6 != "" ) legend->AddEntry(graph6, legendEntry6.data(), "p");
    if ( graph3                       ) legend->AddEntry(graph3, legendEntry3.data(), "p");
    if ( graph7 && legendEntry7 != "" ) legend->AddEntry(graph7, legendEntry7.data(), "p");
    if ( graph4                       ) legend->AddEntry(graph4, legendEntry4.data(), "p");
    if ( graph8 && legendEntry8 != "" ) legend->AddEntry(graph8, legendEntry8.data(), "p");
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
  
  delete dummyHistogram;
  delete label;
  delete legend;
  delete canvas;  
}

void makeEfficiency_vs_RatePlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/L1Trigger/TallinnL1PFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName_signal = "TallinnL1PFTauAnalyzer_signal_qqH_2019Jun21.root";
  TFile* inputFile_signal = openFile(inputFilePath, inputFileName_signal);
  std::string inputFileName_background = "TallinnL1PFTauAnalyzer_background_2019Jun21.root";
  TFile* inputFile_background = openFile(inputFilePath, inputFileName_background);

  std::vector<std::string> pfAlgos;
  pfAlgos.push_back("WithStripsWithoutPreselectionPF");
  pfAlgos.push_back("WithoutStripsAndPreselectionPF");

  std::vector<std::string> absEtaRanges;
  absEtaRanges.push_back("absEtaLt1p40");
  absEtaRanges.push_back("absEta1p40to2p17");
  absEtaRanges.push_back("absEta1p40to2p40");
  absEtaRanges.push_back("absEtaLt2p17");
  absEtaRanges.push_back("absEtaLt2p40");

  std::vector<std::string> isolationWPs;
  isolationWPs.push_back("relChargedIsoLt0p40");
  isolationWPs.push_back("relChargedIsoLt0p20");
  isolationWPs.push_back("relChargedIsoLt0p10");
  isolationWPs.push_back("relChargedIsoLt0p05");

  std::map<std::string, std::string> legendEntries_vs_isolationWPs; // key = isolationWP
  legendEntries_vs_isolationWPs["relChargedIsoLt0p40"] = "I_{ch} < 0.40*p_{T}";
  legendEntries_vs_isolationWPs["relChargedIsoLt0p20"] = "I_{ch} < 0.20*p_{T}";
  legendEntries_vs_isolationWPs["relChargedIsoLt0p10"] = "I_{ch} < 0.10*p_{T}";
  legendEntries_vs_isolationWPs["relChargedIsoLt0p05"] = "I_{ch} < 0.05*p_{T}";

  double rate_min           =  1.e+3; //  1 kHz
  double rate_max           = 20.e+3; // 20 kHz
  double rate_step          =  1.e+3; //  1 kHz
  //double rate_accCorrFactor =  4.;    // extrapolation from |eta| < 1.4 to full HL-LHC tracking acceptance (squared for ditau trigger)
  double rate_accCorrFactor =  1.;    // no extrapolation

  std::string dqmDirectory = "DQMData/TallinnL1PFTauPairAnalyzer";

  int colors[8] = { 8, 4, 1, 2, 8, 4, 1, 2 };
  int lineStyles[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
  int markerStyles[8] = { 20, 22, 21, 23, 24, 26, 25, 32 };

  // TallinnL1PFTaus 
  std::map<std::string, std::map<std::string, std::map<std::string, TGraph*>>> graphs_symmetric;  // keys = pfAlgo, absEtaRange, isolationWP
  std::map<std::string, std::map<std::string, std::map<std::string, TGraph*>>> graphs_asymmetric; // keys = pfAlgo, absEtaRange, isolationWP
  
  for ( std::vector<std::string>::const_iterator pfAlgo = pfAlgos.begin();
	pfAlgo != pfAlgos.end(); ++pfAlgo ) {
    for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	  absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
      for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs.begin();
	    isolationWP != isolationWPs.end(); ++isolationWP ) {
	//std::string histogram2dName_signal = Form("%s%s_wrtOfflineTaus/efficiency_or_rate_%s_%s", 
        //  dqmDirectory.data(), pfAlgo->data(), absEtaRange->data(), isolationWP->data());
	std::string histogram2dName_signal = Form("%s%s_wrtGenHadTaus/efficiency_or_rate_%s_%s", 
          dqmDirectory.data(), pfAlgo->data(), absEtaRange->data(), isolationWP->data());
        TH2* histogram2d_signal = loadHistogram2d(inputFile_signal, histogram2dName_signal);
	std::string histogram2dName_background = Form("%s%s/efficiency_or_rate_%s_%s", 
          dqmDirectory.data(), pfAlgo->data(), absEtaRange->data(), isolationWP->data());
	TH2* histogram2d_background = loadHistogram2d(inputFile_background, histogram2dName_background);

	std::pair<TGraph*, TGraph*> graphs = makeEfficiency_and_RateGraphs(histogram2d_signal, histogram2d_background, rate_min, rate_max, rate_step, rate_accCorrFactor);
	TGraph* graph_symmetric = graphs.first;
	graphs_symmetric[*pfAlgo][*absEtaRange][*isolationWP] = graph_symmetric;
	TGraph* graph_asymmetric = graphs.second;
	graphs_asymmetric[*pfAlgo][*absEtaRange][*isolationWP] = graph_asymmetric;
      }

      std::vector<std::string> labelTextLines1;
      labelTextLines1.push_back("Solid (open) symbols: Asymmetric (symmetric) p_{T} thresholds");
      std::string outputFileName1 = Form("makeEfficiency_vs_RatePlots_HPSatL1_%s_%s.png", 
	pfAlgo->data(), absEtaRange->data());
      showGraphs(1150, 850,
		 graphs_asymmetric[*pfAlgo][*absEtaRange]["relChargedIsoLt0p40"], legendEntries_vs_isolationWPs["relChargedIsoLt0p40"],
		 graphs_asymmetric[*pfAlgo][*absEtaRange]["relChargedIsoLt0p20"], legendEntries_vs_isolationWPs["relChargedIsoLt0p20"],
		 graphs_asymmetric[*pfAlgo][*absEtaRange]["relChargedIsoLt0p10"], legendEntries_vs_isolationWPs["relChargedIsoLt0p10"],
		 graphs_asymmetric[*pfAlgo][*absEtaRange]["relChargedIsoLt0p05"], legendEntries_vs_isolationWPs["relChargedIsoLt0p05"],
		 graphs_symmetric[*pfAlgo][*absEtaRange]["relChargedIsoLt0p40"],  "",
		 graphs_symmetric[*pfAlgo][*absEtaRange]["relChargedIsoLt0p20"],  "",
		 graphs_symmetric[*pfAlgo][*absEtaRange]["relChargedIsoLt0p10"],  "",
		 graphs_symmetric[*pfAlgo][*absEtaRange]["relChargedIsoLt0p05"],  "",		 
		 false,
		 colors, markerStyles, lineStyles, 
		 0.045, 0.23, 0.785, 0.70, 0.135, 
		 labelTextLines1, 0.038,
		 0.145, 0.71, 0.785, 0.042, 
		 rate_min*1.e-3, rate_max*1.e-3, "Rate [kHz]", 1.2, // CV: convert rates to units of kHz
		 false, 0., 1.09, "Efficiency", 1.4, 
		 outputFileName1);
    }
  }

  std::vector<std::string> labelTextLines2;
  labelTextLines2.push_back("Solid (open) symbols: Without (with) strips");
  showGraphs(1150, 850,
	     graphs_symmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p40"],   legendEntries_vs_isolationWPs["relChargedIsoLt0p40"],
	     graphs_symmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p20"],   legendEntries_vs_isolationWPs["relChargedIsoLt0p20"],
	     graphs_symmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p10"],   legendEntries_vs_isolationWPs["relChargedIsoLt0p10"],
	     graphs_symmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p05"],   legendEntries_vs_isolationWPs["relChargedIsoLt0p05"],
	     graphs_symmetric["WithStripsWithoutPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p40"],  "",
	     graphs_symmetric["WithStripsWithoutPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p20"],  "",
	     graphs_symmetric["WithStripsWithoutPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p10"],  "",
	     graphs_symmetric["WithStripsWithoutPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p05"],  "",
	     false,
	     colors, markerStyles, lineStyles, 
	     0.045, 0.23, 0.785, 0.70, 0.135, 
	     labelTextLines2, 0.040,
	     0.23, 0.71, 0.70, 0.045, 
	     rate_min*1.e-3, rate_max*1.e-3, "Rate [kHz]", 1.2, // CV: convert rates to units of kHz
	     false, 0., 1.09, "Efficiency", 1.4, 
	     "makeEfficiency_vs_RatePlots_HPSatL1_symmetricPtThresholds.png");
  showGraphs(1150, 850,
	     graphs_asymmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p40"],  legendEntries_vs_isolationWPs["relChargedIsoLt0p40"],
	     graphs_asymmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p20"],  legendEntries_vs_isolationWPs["relChargedIsoLt0p20"],
	     graphs_asymmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p10"],  legendEntries_vs_isolationWPs["relChargedIsoLt0p10"],
	     graphs_asymmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p05"],  legendEntries_vs_isolationWPs["relChargedIsoLt0p05"],
	     graphs_asymmetric["WithStripsWithoutPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p40"], "",
	     graphs_asymmetric["WithStripsWithoutPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p20"], "",
	     graphs_asymmetric["WithStripsWithoutPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p10"], "",
	     graphs_asymmetric["WithStripsWithoutPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p05"], "",
	     false,
	     colors, markerStyles, lineStyles, 
	     0.045, 0.23, 0.785, 0.70, 0.135, 
	     labelTextLines2, 0.040,
	     0.23, 0.71, 0.70, 0.045, 
	     rate_min*1.e-3, rate_max*1.e-3, "Rate [kHz]", 1.2, // CV: convert rates to units of kHz
	     false, 0., 1.09, "Efficiency", 1.4, 
	     "makeEfficiency_vs_RatePlots_HPSatL1_asymmetricPtThresholds.png");

  // Isobel's L1PFTaus 
  std::vector<std::string> isolationWPs_isobel;
  isolationWPs_isobel.push_back("vLooseIso");
  isolationWPs_isobel.push_back("LooseIso");
  isolationWPs_isobel.push_back("MediumIso");
  isolationWPs_isobel.push_back("TightIso");

  std::map<std::string, std::string> legendEntries_vs_isolationWPs_isobel; // key = isolationWP
  legendEntries_vs_isolationWPs_isobel["vLooseIso"] = "very Loose";
  legendEntries_vs_isolationWPs_isobel["LooseIso"]  = "Loose";
  legendEntries_vs_isolationWPs_isobel["MediumIso"] = "Medium";
  legendEntries_vs_isolationWPs_isobel["TightIso"]  = "Tight";

  std::string dqmDirectory_isobel = "DQMData/L1PFTauPairAnalyzer";

  std::map<std::string, TGraph*> graphs_isobel_symmetric;  // keys = isolationWP
  std::map<std::string, TGraph*> graphs_isobel_asymmetric; // keys = isolationWP

  for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs_isobel.begin();
	isolationWP != isolationWPs_isobel.end(); ++isolationWP ) {
    //std::string histogram2dName_signal = Form("%sPF_wrtOfflineTaus/efficiency_or_rate_absEtaLt1p40_%s", 
    //  dqmDirectory_isobel.data(), isolationWP->data());
    std::string histogram2dName_signal = Form("%sPF_wrtGenHadTaus/efficiency_or_rate_absEtaLt1p40_%s", 
      dqmDirectory_isobel.data(), isolationWP->data());
    TH2* histogram2d_signal = loadHistogram2d(inputFile_signal, histogram2dName_signal);
    std::string histogram2dName_background = Form("%sPF/efficiency_or_rate_absEtaLt1p40_%s", 
      dqmDirectory_isobel.data(), isolationWP->data());
    TH2* histogram2d_background = loadHistogram2d(inputFile_background, histogram2dName_background);

    std::pair<TGraph*, TGraph*> graphs = makeEfficiency_and_RateGraphs(histogram2d_signal, histogram2d_background, rate_min, rate_max, rate_step, rate_accCorrFactor);
    TGraph* graph_symmetric = graphs.first;
    graphs_isobel_symmetric[*isolationWP] = graph_symmetric;
    TGraph* graph_asymmetric = graphs.second;
    graphs_isobel_asymmetric[*isolationWP] = graph_asymmetric;
  }

  std::vector<std::string> labelTextLines3;
  labelTextLines3.push_back("Solid (open) symbols: Asymmetric (symmetric) p_{T} thresholds");
  showGraphs(1150, 850,
	     graphs_isobel_asymmetric["vLooseIso"], legendEntries_vs_isolationWPs_isobel["vLooseIso"],
	     graphs_isobel_asymmetric["LooseIso"],  legendEntries_vs_isolationWPs_isobel["LooseIso"],
	     graphs_isobel_asymmetric["MediumIso"], legendEntries_vs_isolationWPs_isobel["MediumIso"],
	     graphs_isobel_asymmetric["TightIso"],  legendEntries_vs_isolationWPs_isobel["TightIso"],
	     graphs_isobel_symmetric["vLooseIso"],  "",
	     graphs_isobel_symmetric["LooseIso"],   "",
	     graphs_isobel_symmetric["MediumIso"],  "",
	     graphs_isobel_symmetric["TightIso"],   "",		 
	     false,
	     colors, markerStyles, lineStyles, 
	     0.045, 0.23, 0.785, 0.70, 0.135, 
	     labelTextLines3, 0.038,
	     0.145, 0.71, 0.785, 0.042, 
	     rate_min*1.e-3, rate_max*1.e-3, "Rate [kHz]", 1.2, // CV: convert rates to units of kHz
	     false, 0., 1.09, "Efficiency", 1.4, 
	     "makeEfficiency_vs_RatePlots_L1PFTau.png");
  
  std::vector<std::string> labelTextLines4;
  labelTextLines4.push_back("Solid (open) symbols: HPS@L1 (L1PFTau) algorithm");
  showGraphs(1150, 850,
	     graphs_symmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p40"],  legendEntries_vs_isolationWPs["relChargedIsoLt0p40"],
	     graphs_symmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p20"],  legendEntries_vs_isolationWPs["relChargedIsoLt0p20"],
	     graphs_symmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p10"],  legendEntries_vs_isolationWPs["relChargedIsoLt0p10"],
	     graphs_symmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p05"],  legendEntries_vs_isolationWPs["relChargedIsoLt0p05"],
	     graphs_isobel_symmetric["vLooseIso"],                                                       legendEntries_vs_isolationWPs_isobel["vLooseIso"],
	     graphs_isobel_symmetric["LooseIso"],                                                        legendEntries_vs_isolationWPs_isobel["LooseIso"],
	     graphs_isobel_symmetric["MediumIso"],                                                       legendEntries_vs_isolationWPs_isobel["MediumIso"],
	     graphs_isobel_symmetric["TightIso"],                                                        legendEntries_vs_isolationWPs_isobel["TightIso"],
	     false,
	     colors, markerStyles, lineStyles, 
	     0.045, 0.23, 0.685, 0.70, 0.235, 
	     labelTextLines4, 0.040,
	     0.19, 0.61, 0.70, 0.045, 
	     rate_min*1.e-3, rate_max*1.e-3, "Rate [kHz]", 1.2, // CV: convert rates to units of kHz
	     false, 0., 1.09, "Efficiency", 1.4, 
	     "makeEfficiency_vs_RatePlots_HPSatL1_vs_L1PFTau_symmetricPtThresholds.png");
  showGraphs(1150, 850,
	     graphs_asymmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p40"], legendEntries_vs_isolationWPs["relChargedIsoLt0p40"],
	     graphs_asymmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p20"], legendEntries_vs_isolationWPs["relChargedIsoLt0p20"],
	     graphs_asymmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p10"], legendEntries_vs_isolationWPs["relChargedIsoLt0p10"],
	     graphs_asymmetric["WithoutStripsAndPreselectionPF"]["absEtaLt1p40"]["relChargedIsoLt0p05"], legendEntries_vs_isolationWPs["relChargedIsoLt0p05"],
	     graphs_isobel_asymmetric["vLooseIso"],                                                      legendEntries_vs_isolationWPs_isobel["vLooseIso"],
	     graphs_isobel_asymmetric["LooseIso"],                                                       legendEntries_vs_isolationWPs_isobel["LooseIso"],
	     graphs_isobel_asymmetric["MediumIso"],                                                      legendEntries_vs_isolationWPs_isobel["MediumIso"],
	     graphs_isobel_asymmetric["TightIso"],                                                       legendEntries_vs_isolationWPs_isobel["TightIso"],
	     false,
	     colors, markerStyles, lineStyles, 
	     0.045, 0.23, 0.685, 0.70, 0.235, 
	     labelTextLines4, 0.040,
	     0.19, 0.61, 0.70, 0.045, 
	     rate_min*1.e-3, rate_max*1.e-3, "Rate [kHz]", 1.2, // CV: convert rates to units of kHz
	     false, 0., 1.09, "Efficiency", 1.4, 
	     "makeEfficiency_vs_RatePlots_HPSatL1_vs_L1PFTau_asymmetricPtThresholds.png");

  delete inputFile_signal;
  delete inputFile_background;
}

