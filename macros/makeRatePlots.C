
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

double compRate(TH2* histogram2d, int idxBinX, int minNumTaus)
{
  const TAxis* yAxis = histogram2d->GetYaxis();
  int numBinsY = yAxis->GetNbins();
  double integral = 0.;
  double integral_failed = 0.;
  for ( int idxBinY = 1; idxBinY <= numBinsY; ++idxBinY ) { 
    double binContent = histogram2d->GetBinContent(idxBinX, idxBinY);
    integral += binContent;
    double numTaus = yAxis->GetBinCenter(idxBinY);
    if ( numTaus < minNumTaus ) integral_failed += binContent;
  }
  double integral_passed = integral - integral_failed;
  double rate = 2.8e+7; // bunch-crossing frequency of colliding bunches = 28 MHz
  rate *= (integral_passed/integral);
  return rate;
}

TH1* makeRateHistogram(TH2* histogram2d, int minNumTaus)
{
  std::string histogramName_rate = Form("%s_rate_minNumTausEq%i", histogram2d->GetName(), minNumTaus);
  TH1* histogram_rate = histogram2d->ProjectionX(histogramName_rate.data());
  histogram_rate->Reset();
  const TAxis* xAxis = histogram2d->GetXaxis();
  int numBinsX = xAxis->GetNbins();
  for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX ) { 
    double rate = compRate(histogram2d, idxBinX, minNumTaus);
    histogram_rate->SetBinContent(idxBinX, rate);
    histogram_rate->SetBinError(idxBinX, 0.); // CV: error computation not implemented (needed) yet
  }
  return histogram_rate;
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

void makeRatePlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/L1Trigger/TallinnL1PFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName = "TallinnL1PFTauAnalyzer_background_2019May31.root";
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

  std::vector<std::string> absEtaRanges;
  absEtaRanges.push_back("absEtaLt1p40");
  absEtaRanges.push_back("absEta1p40to2p17");
  absEtaRanges.push_back("absEtaLt2p17");
  absEtaRanges.push_back("absEtaLt2p40");

  std::vector<std::string> isolationWPs;
  isolationWPs.push_back("relChargedIsoLt0p40");
  isolationWPs.push_back("relChargedIsoLt0p20");
  isolationWPs.push_back("relChargedIsoLt0p10");
  isolationWPs.push_back("relChargedIsoLt0p05");

  std::map<std::string, std::string> legendEntries; // key = isolationWP
  legendEntries["relChargedIsoLt0p40"] = "I_{ch} < 0.40*p_{T}";
  legendEntries["relChargedIsoLt0p20"] = "I_{ch} < 0.20*p_{T}";
  legendEntries["relChargedIsoLt0p10"] = "I_{ch} < 0.10*p_{T}";
  legendEntries["relChargedIsoLt0p05"] = "I_{ch} < 0.05*p_{T}";

  std::string dqmDirectory = "DQMData/TallinnL1PFTauAnalyzerBackground";

  std::vector<std::string> labelTextLines;

  int colors[6] = { 1, 2, 8, 4, 6, 7 };
  int lineStyles[6] = { 1, 1, 1, 1, 1, 1 };

  // TallinnL1PFTaus 
  typedef std::map<std::string, TH1*>              string_to_TH1Map1;
  typedef std::map<std::string, string_to_TH1Map1> string_to_TH1Map2;
  typedef std::map<std::string, string_to_TH1Map2> string_to_TH1Map3;
  string_to_TH1Map3 histograms_rateSingleTau; // key = pfAlgo, absEtaRange, isolationWP
  string_to_TH1Map3 histograms_rateDoubleTau; // key = pfAlgo, absEtaRange, isolationWP

  for ( std::vector<std::string>::const_iterator pfAlgo = pfAlgos.begin();
	pfAlgo != pfAlgos.end(); ++pfAlgo ) {
    for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	  absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
      for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs.begin();
	    isolationWP != isolationWPs.end(); ++isolationWP ) {
        std::string histogram2dName = Form("%s%s/numL1PFTaus_vs_ptThreshold_%s_%s", 
          dqmDirectory.data(), pfAlgo->data(), absEtaRange->data(), isolationWP->data());
        TH2* histogram2d = loadHistogram2d(inputFile, histogram2dName);

        TH1* histogram_rateSingleTau = makeRateHistogram(histogram2d, 1);
        histograms_rateSingleTau[*pfAlgo][*absEtaRange][*isolationWP] = histogram_rateSingleTau;
        TH1* histogram_rateDoubleTau = makeRateHistogram(histogram2d, 2);
        histograms_rateDoubleTau[*pfAlgo][*absEtaRange][*isolationWP] = histogram_rateDoubleTau;
      }
      
      string_to_TH1Map1 histograms1 = histograms_rateSingleTau[*pfAlgo][*absEtaRange];
      std::string outputFileName1 = Form("makeRatePlots_HPSatL1_SingleTau_%s_%s.png", 
        pfAlgo->data(), absEtaRange->data());
      showHistograms(1150, 1150,
		     histograms1["relChargedIsoLt0p40"], legendEntries["relChargedIsoLt0p40"],
		     histograms1["relChargedIsoLt0p20"], legendEntries["relChargedIsoLt0p20"],
		     histograms1["relChargedIsoLt0p10"], legendEntries["relChargedIsoLt0p10"],
		     histograms1["relChargedIsoLt0p05"], legendEntries["relChargedIsoLt0p05"],
		     0, "",
		     0, "",
		     colors, lineStyles, 
		     0.045, 0.65, 0.70, 0.23, 0.21,
		     labelTextLines, 0.050,
		     0.63, 0.65, 0.26, 0.07, 
		     -1., -1., "L1 #tau p_{T} Threshold [GeV]", 1.2, 
		     true, 1.e+3, 1.e+8, "Single #tau Trigger Rate [Hz]", 1.4, 
		     outputFileName1);
      
      string_to_TH1Map1 histograms2 = histograms_rateDoubleTau[*pfAlgo][*absEtaRange];
      std::string outputFileName_rateDoubleTau = Form("makeRatePlots_HPSatL1_DoubleTau_%s_%s.png", 
        pfAlgo->data(), absEtaRange->data());
      showHistograms(1150, 1150,
		     histograms2["relChargedIsoLt0p40"], legendEntries["relChargedIsoLt0p40"],
		     histograms2["relChargedIsoLt0p20"], legendEntries["relChargedIsoLt0p20"],
		     histograms2["relChargedIsoLt0p10"], legendEntries["relChargedIsoLt0p10"],
		     histograms2["relChargedIsoLt0p05"], legendEntries["relChargedIsoLt0p05"],
		     0, "",
		     0, "",
		     colors, lineStyles, 
		     0.045, 0.65, 0.70, 0.23, 0.21,
		     labelTextLines, 0.050,
		     0.63, 0.65, 0.26, 0.07, 
		     -1., -1., "L1 #tau p_{T} Threshold [GeV]", 1.2, 
		     true, 1.e+3, 1.e+8, "Double #tau Trigger Rate [Hz]", 1.4, 
		     outputFileName_rateDoubleTau);
    }
  }

  // Isobel's L1PFTaus 
  std::vector<std::string> isolationWPs_isobel;
  isolationWPs_isobel.push_back("vLooseIso");
  isolationWPs_isobel.push_back("LooseIso");
  isolationWPs_isobel.push_back("MediumIso");
  isolationWPs_isobel.push_back("TightIso");

  std::map<std::string, std::string> legendEntries_isobel; // key = isolationWP
  legendEntries_isobel["vLooseIso"] = "very Loose";
  legendEntries_isobel["LooseIso"]  = "Loose";
  legendEntries_isobel["MediumIso"] = "Medium";
  legendEntries_isobel["TightIso"]  = "Tight";

  std::string dqmDirectory_isobel = "DQMData/L1PFTauAnalyzerBackground";

  string_to_TH1Map2 histograms_rateSingleTau_isobel; // key = absEtaRange, isolationWP
  string_to_TH1Map2 histograms_rateDoubleTau_isobel; // key = absEtaRange, isolationWP

  for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
    for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs_isobel.begin();
	  isolationWP != isolationWPs_isobel.end(); ++isolationWP ) {
      std::string histogram2dName = Form("%s/numL1PFTaus_vs_ptThreshold_%s_%s", 
	dqmDirectory_isobel.data(), absEtaRange->data(), isolationWP->data());
      TH2* histogram2d = loadHistogram2d(inputFile, histogram2dName);

      TH1* histogram_rateSingleTau = makeRateHistogram(histogram2d, 1);
      histograms_rateSingleTau_isobel[*absEtaRange][*isolationWP] = histogram_rateSingleTau;
      TH1* histogram_rateDoubleTau = makeRateHistogram(histogram2d, 2);
      histograms_rateDoubleTau_isobel[*absEtaRange][*isolationWP] = histogram_rateDoubleTau;
    }
      
    string_to_TH1Map1 histograms3 = histograms_rateSingleTau_isobel[*absEtaRange];
    std::string outputFileName3 = Form("makeRatePlots_L1PFTau_SingleTau_%s.png", 
      absEtaRange->data());
    showHistograms(1150, 1150,
		   histograms3["vLooseIso"], legendEntries_isobel["vLooseIso"],
		   histograms3["LooseIso"],  legendEntries_isobel["LooseIso"],
		   histograms3["MediumIso"], legendEntries_isobel["MediumIso"],
		   histograms3["TightIso"],  legendEntries_isobel["TightIso"],
		   0, "",
		   0, "",
		   colors, lineStyles, 
		   0.045, 0.65, 0.70, 0.23, 0.21,
		   labelTextLines, 0.050,
		   0.63, 0.65, 0.26, 0.07, 
		   -1., -1., "L1 #tau p_{T} Threshold [GeV]", 1.2, 
		   true, 1.e+3, 1.e+8, "Single #tau Trigger Rate [Hz]", 1.4, 
		   outputFileName3);
      
    string_to_TH1Map1 histograms4 = histograms_rateDoubleTau_isobel[*absEtaRange];							  
    std::string outputFileName4 = Form("makeRatePlots_L1PFTau_DoubleTau_%s.png", 
      absEtaRange->data());
    showHistograms(1150, 1150,
		   histograms4["vLooseIso"], legendEntries_isobel["vLooseIso"],
		   histograms4["LooseIso"],  legendEntries_isobel["LooseIso"],
		   histograms4["MediumIso"], legendEntries_isobel["MediumIso"],
		   histograms4["TightIso"],  legendEntries_isobel["TightIso"],
		   0, "",
		   0, "",
		   colors, lineStyles, 
		   0.045, 0.65, 0.70, 0.23, 0.21,
		   labelTextLines, 0.050,
		   0.63, 0.65, 0.26, 0.07, 
		   -1., -1., "L1 #tau p_{T} Threshold [GeV]", 1.2, 
		   true, 1.e+3, 1.e+8, "Double #tau Trigger Rate [Hz]", 1.4, 
		   outputFileName4);
  }

  for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
    assert(isolationWPs.size() == isolationWPs_isobel.size());
    for ( int idxIsolationWP = 0; idxIsolationWP < 4; ++idxIsolationWP ) {
      const std::string& isolationWP = isolationWPs[idxIsolationWP];
      TH1* histogram_rateSingleTau = histograms_rateSingleTau["WithoutStripsAndPreselectionPF"][*absEtaRange][isolationWP];
      TH1* histogram_rateDoubleTau = histograms_rateDoubleTau["WithoutStripsAndPreselectionPF"][*absEtaRange][isolationWP];
      const std::string& isolationWP_isobel = isolationWPs_isobel[idxIsolationWP];
      TH1* histogram_rateSingleTau_isobel = histograms_rateSingleTau_isobel[*absEtaRange][isolationWP_isobel];
      TH1* histogram_rateDoubleTau_isobel = histograms_rateDoubleTau_isobel[*absEtaRange][isolationWP_isobel];
      
      std::string outputFileName5 = Form("makeRatePlots_HPSatL1_vs_L1PFTau_SingleTau_%s.png", 
        absEtaRange->data());
      showHistograms(1150, 1150,
		     histogram_rateSingleTau,        "HPS@L1 (Tallinn)",
		     histogram_rateSingleTau_isobel, "L1PFTau",
		     0, "",
		     0, "",
		     0, "",
		     0, "",
		     colors, lineStyles, 
		     0.045, 0.65, 0.70, 0.23, 0.21,
		     labelTextLines, 0.050,
		     0.63, 0.65, 0.26, 0.07, 
		     -1., -1., "L1 #tau p_{T} Threshold [GeV]", 1.2, 
		     true, 1.e+3, 1.e+8, "Single #tau Trigger Rate [Hz]", 1.4, 
		     outputFileName5);
      
      std::string outputFileName6 = Form("makeRatePlots_HPSatL1_vs_L1PFTau_DoubleTau_%s.png", 
        absEtaRange->data());
      showHistograms(1150, 1150,
		     histogram_rateDoubleTau,        "HPS@L1 (Tallinn)",
		     histogram_rateDoubleTau_isobel, "L1PFTau",
		     0, "",
		     0, "",
		     0, "",
		     0, "",
		     colors, lineStyles, 
		     0.045, 0.65, 0.70, 0.23, 0.21,
		     labelTextLines, 0.050,
		     0.63, 0.65, 0.26, 0.07, 
		     -1., -1., "L1 #tau p_{T} Threshold [GeV]", 1.2, 
		     true, 1.e+3, 1.e+8, "Double #tau Trigger Rate [Hz]", 1.4, 
		     outputFileName6);
    }
  }

  delete inputFile;
}

