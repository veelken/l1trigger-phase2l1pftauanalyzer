#ifndef L1Trigger_TallinnL1PFTauAnalyzer_L1PFTauAnalyzerBackground_h
#define L1Trigger_TallinnL1PFTauAnalyzer_L1PFTauAnalyzerBackground_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/L1Trigger/interface/L1PFTau.h"                          // l1t::L1PFTau, l1t::L1PFTauCollection
#include "L1Trigger/TallinnL1PFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverFlow()

#include <TH1.h>     // TH1
#include <TH2.h>     // TH2
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs(), TMath::Nint()

#include <vector>    // std::vector
#include <string>    // std::string
#include <algorithm> // std::sort

namespace
{
  bool
  isHigherPt(const l1t::L1PFTau* l1PFTau1,
	     const l1t::L1PFTau* l1PFTau2)
  {
    return l1PFTau1->pt() > l1PFTau2->pt();
  }
}

class L1PFTauAnalyzerBackground : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit L1PFTauAnalyzerBackground(const edm::ParameterSet&);
    
  // destructor
  ~L1PFTauAnalyzerBackground();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcL1PFTaus_;
  edm::EDGetTokenT<l1t::L1PFTauCollection> tokenL1PFTaus_;

  std::string dqmDirectory_;

  struct ratePlotEntryType
  {
    ratePlotEntryType(double min_absEta, double max_absEta, const std::string& isolation)
      : me_numL1PFTaus_vs_ptThreshold_(nullptr)
      , histogram_numL1PFTaus_vs_ptThreshold_(nullptr)
      , me_numL1PFTausPtGt20_(nullptr)
      , histogram_numL1PFTausPtGt20_(nullptr)
      , me_numL1PFTausPtGt25_(nullptr)
      , histogram_numL1PFTausPtGt25_(nullptr)
      , me_numL1PFTausPtGt30_(nullptr)
      , histogram_numL1PFTausPtGt30_(nullptr)
      , me_numL1PFTausPtGt35_(nullptr)
      , histogram_numL1PFTausPtGt35_(nullptr)
      , me_numL1PFTausPtGt40_(nullptr)
      , histogram_numL1PFTausPtGt40_(nullptr) 
      , min_absEta_(min_absEta)
      , max_absEta_(max_absEta)
      , isolation_string_(isolation)
    {
      if      ( isolation_string_ == ""       ) isolation_ = kNone;
      else if ( isolation_string_ == "vLoose" ) isolation_ = kVLooseIso;
      else if ( isolation_string_ == "Loose"  ) isolation_ = kLooseIso;
      else if ( isolation_string_ == "Medium" ) isolation_ = kMediumIso;
      else if ( isolation_string_ == "Tight"  ) isolation_ = kTightIso; 
      else throw cms::Exception("L1PFTauAnalyzerBackground") 
	<< " Invalid Configuration parameter 'isolation' = " << isolation_string_ << " !!\n";
    }
    ~ratePlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix;
      if      ( min_absEta_ >= 0. && max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEta%1.2fto%1.2f", min_absEta_, max_absEta_));
      else if ( min_absEta_ >= 0.                     ) histogramName_suffix.Append(Form("_absEtaGt%1.2f", min_absEta_));
      else if (                      max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f", max_absEta_));
      if ( isolation_string_ != "" ) histogramName_suffix.Append(Form("_%sIso", isolation_string_.data()));
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_numL1PFTaus_vs_ptThreshold = Form("numL1PFTaus_vs_ptThreshold%s", histogramName_suffix.Data());
      me_numL1PFTaus_vs_ptThreshold_ = dqmStore.book2D(histogramName_numL1PFTaus_vs_ptThreshold.Data(), histogramName_numL1PFTaus_vs_ptThreshold.Data(), 101, -0.5, 100.5, 11, -0.5, +10.5);
      histogram_numL1PFTaus_vs_ptThreshold_ = dynamic_cast<TH2*>(me_numL1PFTaus_vs_ptThreshold_->getTH1());
      assert(histogram_numL1PFTaus_vs_ptThreshold_);

      TString histogramName_numL1PFTausPtGt20 = Form("numL1PFTausPtGt20%s", histogramName_suffix.Data());
      me_numL1PFTausPtGt20_ = dqmStore.book1D(histogramName_numL1PFTausPtGt20.Data(), histogramName_numL1PFTausPtGt20.Data(), 21, -0.5, +20.5);
      histogram_numL1PFTausPtGt20_ = me_numL1PFTausPtGt20_->getTH1();
      assert(histogram_numL1PFTausPtGt20_);
      TString histogramName_numL1PFTausPtGt25 = Form("numL1PFTausPtGt25%s", histogramName_suffix.Data());
      me_numL1PFTausPtGt25_ = dqmStore.book1D(histogramName_numL1PFTausPtGt25.Data(), histogramName_numL1PFTausPtGt25.Data(), 21, -0.5, +20.5);
      histogram_numL1PFTausPtGt25_ = me_numL1PFTausPtGt25_->getTH1();
      assert(histogram_numL1PFTausPtGt25_);
      TString histogramName_numL1PFTausPtGt30 = Form("numL1PFTausPtGt30%s", histogramName_suffix.Data());
      me_numL1PFTausPtGt30_ = dqmStore.book1D(histogramName_numL1PFTausPtGt30.Data(), histogramName_numL1PFTausPtGt30.Data(), 21, -0.5, +20.5);
      histogram_numL1PFTausPtGt30_ = me_numL1PFTausPtGt30_->getTH1();
      assert(histogram_numL1PFTausPtGt30_);
      TString histogramName_numL1PFTausPtGt35 = Form("numL1PFTausPtGt35%s", histogramName_suffix.Data());
      me_numL1PFTausPtGt35_ = dqmStore.book1D(histogramName_numL1PFTausPtGt35.Data(), histogramName_numL1PFTausPtGt35.Data(), 21, -0.5, +20.5);
      histogram_numL1PFTausPtGt35_ = me_numL1PFTausPtGt35_->getTH1();
      assert(histogram_numL1PFTausPtGt35_);
      TString histogramName_numL1PFTausPtGt40 = Form("numL1PFTausPtGt40%s", histogramName_suffix.Data());
      me_numL1PFTausPtGt40_ = dqmStore.book1D(histogramName_numL1PFTausPtGt40.Data(), histogramName_numL1PFTausPtGt40.Data(), 21, -0.5, +20.5);
      histogram_numL1PFTausPtGt40_ = me_numL1PFTausPtGt40_->getTH1();
      assert(histogram_numL1PFTausPtGt40_);
    }
    void fillHistograms(const l1t::L1PFTauCollection& l1PFTaus, double evtWeight)
    {
      std::vector<const l1t::L1PFTau*> l1PFTaus_passingAbsEta;
      for ( l1t::L1PFTauCollection::const_iterator l1PFTau = l1PFTaus.begin(); 
	    l1PFTau != l1PFTaus.end(); ++l1PFTau ) {
	double l1PFTau_absEta = TMath::Abs(l1PFTau->eta());
	if ( (min_absEta_        < 0. || l1PFTau_absEta             >=  min_absEta_                      ) &&
	     (max_absEta_        < 0. || l1PFTau_absEta             <=  max_absEta_                      ) &&
	     ( isolation_ == kNone                                                                       || 
	      (isolation_ == kVLooseIso && l1PFTau->passVLooseIso())                                     ||
	      (isolation_ == kLooseIso  && l1PFTau->passLooseIso())                                      ||
	      (isolation_ == kMediumIso && l1PFTau->passMediumIso())                                     ||
	      (isolation_ == kTightIso  && l1PFTau->passTightIso())                                      ) )
	{
	  l1PFTaus_passingAbsEta.push_back(&(*l1PFTau));
	}
      }

      // CV: sort L1PFTaus passing abs(eta) cut by decreasing pT
      std::sort(l1PFTaus_passingAbsEta.begin(), l1PFTaus_passingAbsEta.end(), isHigherPt);

      int numL1PFTausPtGt20 = 0;
      int numL1PFTausPtGt25 = 0;
      int numL1PFTausPtGt30 = 0;
      int numL1PFTausPtGt35 = 0;
      int numL1PFTausPtGt40 = 0;
      for ( std::vector<const l1t::L1PFTau*>::const_iterator l1PFTau = l1PFTaus_passingAbsEta.begin();
	    l1PFTau != l1PFTaus_passingAbsEta.end(); ++l1PFTau ) {
	if ( (*l1PFTau)->pt() > 20. ) ++numL1PFTausPtGt20;
	if ( (*l1PFTau)->pt() > 25. ) ++numL1PFTausPtGt25;
	if ( (*l1PFTau)->pt() > 30. ) ++numL1PFTausPtGt30;
	if ( (*l1PFTau)->pt() > 35. ) ++numL1PFTausPtGt35;
	if ( (*l1PFTau)->pt() > 40. ) ++numL1PFTausPtGt40;	
      }
      fillWithOverFlow(histogram_numL1PFTausPtGt20_, numL1PFTausPtGt20, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt25_, numL1PFTausPtGt25, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt30_, numL1PFTausPtGt30, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt35_, numL1PFTausPtGt35, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt40_, numL1PFTausPtGt40, evtWeight);

      TAxis* xAxis = histogram_numL1PFTaus_vs_ptThreshold_->GetXaxis();
      TAxis* yAxis = histogram_numL1PFTaus_vs_ptThreshold_->GetYaxis();
      int numBinsX = xAxis->GetNbins();
      bool numL1PFTaus_passingPt_isZero = false;
      for ( int idxBin = 1; idxBin <= numBinsX; ++idxBin )
      {
	double ptThreshold = xAxis->GetBinCenter(idxBin);

	int numL1PFTaus_passingPt = 0;
        if ( !numL1PFTaus_passingPt_isZero ) 
        {
	  for ( std::vector<const l1t::L1PFTau*>::const_iterator l1PFTau = l1PFTaus_passingAbsEta.begin(); 
		l1PFTau != l1PFTaus_passingAbsEta.end(); ++l1PFTau ) {
	    if ( (*l1PFTau)->pt() > ptThreshold )
	    {
	      ++numL1PFTaus_passingPt;
	    }
	  }
	}

	if ( numL1PFTaus_passingPt > yAxis->GetXmax() ) numL1PFTaus_passingPt = TMath::Nint(yAxis->GetXmax() - 0.5);
	histogram_numL1PFTaus_vs_ptThreshold_->Fill(ptThreshold, numL1PFTaus_passingPt, evtWeight);

	if ( numL1PFTaus_passingPt == 0 ) 
	{
	  numL1PFTaus_passingPt_isZero = true;
	}
      }
    }
    MonitorElement* me_numL1PFTaus_vs_ptThreshold_;
    TH2* histogram_numL1PFTaus_vs_ptThreshold_;
    MonitorElement* me_numL1PFTausPtGt20_;
    TH1* histogram_numL1PFTausPtGt20_;
    MonitorElement* me_numL1PFTausPtGt25_;
    TH1* histogram_numL1PFTausPtGt25_;
    MonitorElement* me_numL1PFTausPtGt30_;
    TH1* histogram_numL1PFTausPtGt30_;
    MonitorElement* me_numL1PFTausPtGt35_;
    TH1* histogram_numL1PFTausPtGt35_;
    MonitorElement* me_numL1PFTausPtGt40_;
    TH1* histogram_numL1PFTausPtGt40_;    
    double min_absEta_;    
    double max_absEta_;    
    std::string isolation_string_;
    enum { kNone, kVLooseIso, kLooseIso, kMediumIso, kTightIso };
    int isolation_;
  };
  std::vector<ratePlotEntryType*> ratePlots_;
};

#endif   
