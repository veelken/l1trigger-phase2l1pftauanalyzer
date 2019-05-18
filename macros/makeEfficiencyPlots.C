
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
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

TGraph* makeEfficiencyGraph(TH1* histogram_numerator, TH1* histogram_denominator)
{
  assert(histogram_numerator->GetNbinsX() == histogram_denominator->GetNbinsX());
  int numBinsX = histogram_numerator->GetNbinsX();
  for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX ) { 
    double binContent_numerator = histogram_numerator->GetBinContent(idxBinX);
    double binContent_denominator = histogram_denominator->GetBinContent(idxBinX);
    if( binContent_numerator > binContent_denominator ) {
      std::cerr << "Error in <makeEfficiencyGraph>: numerator = " << binContent_numerator << " exceeds denominator = " << binContent_denominator 
		<< " @ x = " << histogram_denominator->GetBinCenter(idxBinX) << " !!" << std::endl;
      assert(0);
    }
  }
  TString graphName_efficiency = TString(histogram_numerator->GetName()).ReplaceAll("_numerator", "");
  TGraphAsymmErrors* graph_efficiency = new TGraphAsymmErrors(numBinsX);
  graph_efficiency->Divide(histogram_numerator, histogram_denominator, "w");
  //TH1* histogram2_numerator = new TH1F("numerator", "numerator", 1, -0.5, +0.5);
  //histogram2_numerator->Fill(0);
  //TH1* histogram2_denominator = new TH1F("denominator", "denominator", 1, -0.5, +0.5);
  //histogram2_denominator->Fill(0);
  //histogram2_denominator->Fill(0);
  //graph_efficiency->Divide(histogram2_numerator, histogram2_denominator, "w");
  //double x, y;
  //graph_efficiency->GetPoint(0, x, y);
  //std::cout << "efficiency (x = " << x << "): " << y << std::endl;
  //assert(0);
  graph_efficiency->SetName(graphName_efficiency.Data());
  return graph_efficiency;
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

  TLegend* legend = 0;
  if ( legendEntry1 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    legend->AddEntry(graph1, legendEntry1.data(), "p");
    if ( graph2 ) legend->AddEntry(graph2, legendEntry2.data(), "p");
    if ( graph3 ) legend->AddEntry(graph3, legendEntry3.data(), "p");
    if ( graph4 ) legend->AddEntry(graph4, legendEntry4.data(), "p");
    if ( graph5 ) legend->AddEntry(graph5, legendEntry5.data(), "p");
    if ( graph6 ) legend->AddEntry(graph6, legendEntry6.data(), "p");
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

void makeEfficiencyPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/L1Trigger/TallinnL1PFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName = "TallinnL1PFTauAnalyzer_signal_2019May14.root";
  std::string inputFileName_full = inputFilePath;
  if ( inputFileName_full.find_last_of("/") != (inputFileName_full.size() - 1) ) inputFileName_full.append("/");
  inputFileName_full.append(inputFileName);
  TFile* inputFile = new TFile(inputFileName_full.data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = '" << inputFileName_full << "' !!" << std::endl;
    assert(0);
  }

  std::vector<std::string> pfAlgos;
  pfAlgos.push_back("WithStripsAndPreselectionPF");
  pfAlgos.push_back("WithStripsWithoutPreselectionPF");
  pfAlgos.push_back("WithoutStripsWithPreselectionPF");
  pfAlgos.push_back("WithoutStripsAndPreselectionPF");
  pfAlgos.push_back("WithStripsAndPreselectionPuppi");
  pfAlgos.push_back("WithStripsWithoutPreselectionPuppi");
  pfAlgos.push_back("WithoutStripsWithPreselectionPuppi");
  pfAlgos.push_back("WithoutStripsAndPreselectionPuppi");

  std::vector<std::string> observables;
  observables.push_back("pt");
  observables.push_back("eta");

  std::vector<std::string> absEtaRanges;
  absEtaRanges.push_back("absEtaLt1p00");
  absEtaRanges.push_back("absEtaLt1p40");

  std::vector<std::string> decayModes;
  decayModes.push_back("oneProng0Pi0");
  decayModes.push_back("oneProng1Pi0");
  decayModes.push_back("oneProng2Pi0");
  decayModes.push_back("threeProng0Pi0");
  decayModes.push_back("threeProng1Pi0");
  decayModes.push_back("all");

  std::vector<std::string> ptThresholds;
  ptThresholds.push_back("ptGt20");
  ptThresholds.push_back("ptGt30");
  ptThresholds.push_back("ptGt40");

  std::vector<std::string> isolationWPs;
  isolationWPs.push_back("relChargedIsoLt0p40");
  isolationWPs.push_back("relChargedIsoLt0p20");
  isolationWPs.push_back("relChargedIsoLt0p10");
  isolationWPs.push_back("relChargedIsoLt0p05");

  std::map<std::string, double> xMin; // key = observable
  xMin["pt"]  =   0.;
  xMin["eta"] =  -3.0;

  std::map<std::string, double> xMax; // key = observable
  xMax["pt"]  = 100.;
  xMax["eta"] =  +3.0;

  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["pt"] = "True #tau_{h} p_{T} [GeV]";
  xAxisTitles["eta"] = "True #tau_{h} #eta";

  std::map<std::string, std::string> legendEntries_vs_isolationWPs; // key = isolationWP
  legendEntries_vs_isolationWPs["relChargedIsoLt0p40"] = "I_{ch} < 0.40*p_{T}";
  legendEntries_vs_isolationWPs["relChargedIsoLt0p20"] = "I_{ch} < 0.20*p_{T}";
  legendEntries_vs_isolationWPs["relChargedIsoLt0p10"] = "I_{ch} < 0.10*p_{T}";
  legendEntries_vs_isolationWPs["relChargedIsoLt0p05"] = "I_{ch} < 0.05*p_{T}";

  std::map<std::string, std::string> legendEntries_vs_decayModes; // key = decayMode
  legendEntries_vs_decayModes["oneProng0Pi0"]   = "h^{#pm}";
  legendEntries_vs_decayModes["oneProng1Pi0"]   = "h^{#pm}#pi^{0}";
  legendEntries_vs_decayModes["oneProng2Pi0"]   = "h^{#pm}#pi^{0}#pi^{0}";
  legendEntries_vs_decayModes["threeProng0Pi0"] = "h^{#pm}h^{#mp}h^{#pm}";
  legendEntries_vs_decayModes["threeProng1Pi0"] = "h^{#pm}h^{#mp}h^{#pm}#pi^{0}";
  legendEntries_vs_decayModes["all"]            = "all";

  std::string dqmDirectory = "DQMData/TallinnL1PFTauAnalyzerSignal";

  int colors[6] = { 1, 2, 8, 4, 6, 7 };
  int lineStyles[6] = { 1, 1, 1, 1, 1, 1 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  for ( std::vector<std::string>::const_iterator pfAlgo = pfAlgos.begin();
	pfAlgo != pfAlgos.end(); ++pfAlgo ) {
    for ( std::vector<std::string>::const_iterator observable = observables.begin();
	observable != observables.end(); ++observable ) {
      for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	    absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
        for ( std::vector<std::string>::const_iterator ptThreshold = ptThresholds.begin();
	      ptThreshold != ptThresholds.end(); ++ptThreshold ) {      
          std::map<std::string, TGraph*> graphs_efficiency_vs_isolationWPs; // key = isolationWP
          for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs.begin();
	        isolationWP != isolationWPs.end(); ++isolationWP ) {
            std::string histogramName_numerator = Form("%s%s/effL1PFTau_vs_%s_numerator_all_%s_%s_%s", 
              dqmDirectory.data(), pfAlgo->data(), observable->data(), absEtaRange->data(), ptThreshold->data(), isolationWP->data());
            TH1* histogram_numerator = loadHistogram(inputFile, histogramName_numerator);
	    std::string histogramName_denominator = Form("%s%s/effL1PFTau_vs_%s_denominator_all_%s_%s_%s", 
              dqmDirectory.data(), pfAlgo->data(), observable->data(), absEtaRange->data(), ptThreshold->data(), isolationWP->data());
            TH1* histogram_denominator = loadHistogram(inputFile, histogramName_denominator);
	    TGraph* graph_efficiency = makeEfficiencyGraph(histogram_numerator, histogram_denominator);
	    graphs_efficiency_vs_isolationWPs[*isolationWP] = graph_efficiency;
          }

	  std::vector<std::string> labelTextLines = getLabelTextLines(*ptThreshold);
          std::string outputFileName_efficiency_vs_isolationWPs = Form("makeEfficiencyPlots_%s_vs_%s_%s_%s.png", 
            pfAlgo->data(), observable->data(), absEtaRange->data(), ptThreshold->data());
	  showGraphs(1150, 850,
		     graphs_efficiency_vs_isolationWPs["relChargedIsoLt0p40"], legendEntries_vs_isolationWPs["relChargedIsoLt0p40"],
		     graphs_efficiency_vs_isolationWPs["relChargedIsoLt0p20"], legendEntries_vs_isolationWPs["relChargedIsoLt0p20"],
		     graphs_efficiency_vs_isolationWPs["relChargedIsoLt0p10"], legendEntries_vs_isolationWPs["relChargedIsoLt0p10"],
		     graphs_efficiency_vs_isolationWPs["relChargedIsoLt0p05"], legendEntries_vs_isolationWPs["relChargedIsoLt0p05"],
		     0, "",
		     0, "",
		     colors, markerStyles, lineStyles, 
		     0.045, 0.18, 0.17, 0.23, 0.26, 
		     labelTextLines, 0.050,
		     0.63, 0.65, 0.26, 0.07, 
		     xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
		     false, 0., 1.09, "Efficiency", 1.4, 
		     outputFileName_efficiency_vs_isolationWPs);
	
  	  for ( std::map<std::string, TGraph*>::const_iterator it = graphs_efficiency_vs_isolationWPs.begin(); it != graphs_efficiency_vs_isolationWPs.end(); ++it )
          {
  	    delete it->second;
          }

  	  for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs.begin();
	        isolationWP != isolationWPs.end(); ++isolationWP ) {  
	    std::map<std::string, TGraph*> graphs_efficiency_vs_decayModes; // key = decayMode
	    for ( std::vector<std::string>::const_iterator decayMode = decayModes.begin();
		  decayMode != decayModes.end(); ++decayMode ) {
  	      std::string histogramName_numerator = Form("%s%s/effL1PFTau_vs_%s_numerator_%s_%s_%s_%s", 
	        dqmDirectory.data(), pfAlgo->data(), observable->data(), decayMode->data(), absEtaRange->data(), ptThreshold->data(), isolationWP->data());
              TH1* histogram_numerator = loadHistogram(inputFile, histogramName_numerator);
	      std::string histogramName_denominator = Form("%s%s/effL1PFTau_vs_%s_denominator_%s_%s_%s_%s", 
                dqmDirectory.data(), pfAlgo->data(), observable->data(), decayMode->data(), absEtaRange->data(), ptThreshold->data(), isolationWP->data());
              TH1* histogram_denominator = loadHistogram(inputFile, histogramName_denominator);
	      TGraph* graph_efficiency = makeEfficiencyGraph(histogram_numerator, histogram_denominator);
	      graphs_efficiency_vs_decayModes[*decayMode] = graph_efficiency;
            }
	  
	    std::vector<std::string> labelTextLines = getLabelTextLines(*ptThreshold);
	    std::string outputFileName_efficiency_vs_decayModes = Form("makeEfficiencyPlots_%s_vs_%s_%s_%s_%s.png", 
              pfAlgo->data(), observable->data(), absEtaRange->data(), ptThreshold->data(), isolationWP->data());
	    showGraphs(1150, 850,
		       graphs_efficiency_vs_decayModes["oneProng0Pi0"], legendEntries_vs_decayModes["oneProng0Pi0"],
		       graphs_efficiency_vs_decayModes["oneProng1Pi0"], legendEntries_vs_decayModes["oneProng1Pi0"],
		       graphs_efficiency_vs_decayModes["oneProng2Pi0"], legendEntries_vs_decayModes["oneProng2Pi0"],
		       graphs_efficiency_vs_decayModes["threeProng0Pi0"], legendEntries_vs_decayModes["threeProng0Pi0"],
		       graphs_efficiency_vs_decayModes["threeProng1Pi0"], legendEntries_vs_decayModes["threeProng1Pi0"],
		       graphs_efficiency_vs_decayModes["threeProng1Pi0"], legendEntries_vs_decayModes["threeProng1Pi0"],
		       colors, markerStyles, lineStyles, 
		       0.045, 0.18, 0.17, 0.23, 0.26, 
		       labelTextLines, 0.050,
		       0.63, 0.65, 0.26, 0.07, 
		       xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
		       false, 0., 1.09, "Efficiency", 1.4, 
		       outputFileName_efficiency_vs_decayModes);
	
  	    for ( std::map<std::string, TGraph*>::const_iterator it = graphs_efficiency_vs_decayModes.begin(); it != graphs_efficiency_vs_decayModes.end(); ++it )
	    {
	      delete it->second;
            }
	  }
	}
      }
    }
  }

  delete inputFile;
}

