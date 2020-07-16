#ifndef L1Trigger_TallinnL1PFTauAnalyzer_L1HPSPFTauAnalyzerBackground_h
#define L1Trigger_TallinnL1PFTauAnalyzer_L1HPSPFTauAnalyzerBackground_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/Phase2L1Taus/interface/L1HPSPFTau.h"                    // l1t::L1HPSPFTau
#include "DataFormats/Phase2L1Taus/interface/L1HPSPFTauFwd.h"                 // l1t::L1HPSPFTauCollection
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
  isHigherPt(const l1t::L1HPSPFTau* l1PFTau1,
	     const l1t::L1HPSPFTau* l1PFTau2)
  {
    return l1PFTau1->pt() > l1PFTau2->pt();
  }
}

using namespace dqm::implementation;

class L1HPSPFTauAnalyzerBackground : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit L1HPSPFTauAnalyzerBackground(const edm::ParameterSet&);
    
  // destructor
  ~L1HPSPFTauAnalyzerBackground();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcL1PFTaus_;
  edm::EDGetTokenT<l1t::L1HPSPFTauCollection> tokenL1PFTaus_;

  std::string dqmDirectory_;

  struct ratePlotEntryType
  {
    ratePlotEntryType(double min_absEta, double max_absEta, double max_relChargedIso, double max_absChargedIso, double max_dz)
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
      , me_numL1PFTausPtGt45_(nullptr)
      , histogram_numL1PFTausPtGt45_(nullptr)
      , me_numL1PFTausPtGt50_(nullptr)
      , histogram_numL1PFTausPtGt50_(nullptr)
      , min_absEta_(min_absEta)
      , max_absEta_(max_absEta)
      , max_relChargedIso_(max_relChargedIso)
      , max_absChargedIso_(max_absChargedIso)
      , max_dz_(max_dz)
    {}
    ~ratePlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix;
      if      ( min_absEta_ >= 0. && max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEta%1.2fto%1.2f", min_absEta_, max_absEta_));
      else if ( min_absEta_ >= 0.                     ) histogramName_suffix.Append(Form("_absEtaGt%1.2f", min_absEta_));
      else if (                      max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f", max_absEta_));
      if ( max_relChargedIso_ > 0. ) histogramName_suffix.Append(Form("_relChargedIsoLt%1.2f", max_relChargedIso_));
      if ( max_absChargedIso_ > 0. ) histogramName_suffix.Append(Form("_absChargedIsoLt%1.2f", max_absChargedIso_));
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
      TString histogramName_numL1PFTausPtGt45 = Form("numL1PFTausPtGt45%s", histogramName_suffix.Data());
      me_numL1PFTausPtGt45_ = dqmStore.book1D(histogramName_numL1PFTausPtGt45.Data(), histogramName_numL1PFTausPtGt45.Data(), 21, -0.5, +20.5);
      histogram_numL1PFTausPtGt45_ = me_numL1PFTausPtGt45_->getTH1();
      assert(histogram_numL1PFTausPtGt45_);
      TString histogramName_numL1PFTausPtGt50 = Form("numL1PFTausPtGt50%s", histogramName_suffix.Data());
      me_numL1PFTausPtGt50_ = dqmStore.book1D(histogramName_numL1PFTausPtGt50.Data(), histogramName_numL1PFTausPtGt50.Data(), 21, -0.5, +20.5);
      histogram_numL1PFTausPtGt50_ = me_numL1PFTausPtGt50_->getTH1();
      assert(histogram_numL1PFTausPtGt50_);
    }
    void fillHistograms(const l1t::L1HPSPFTauCollection& l1PFTaus, double evtWeight)
    {
      std::vector<const l1t::L1HPSPFTau*> l1PFTaus_passingAbsEta;
      for ( l1t::L1HPSPFTauCollection::const_iterator l1PFTau = l1PFTaus.begin(); 
	    l1PFTau != l1PFTaus.end(); ++l1PFTau ) {
	double l1PFTau_absEta = TMath::Abs(l1PFTau->eta());
	if ( (min_absEta_        < 0. || l1PFTau_absEta             >=  min_absEta_                      ) &&
	     (max_absEta_        < 0. || l1PFTau_absEta             <=  max_absEta_                      ) &&
	     (max_relChargedIso_ < 0. || l1PFTau->sumChargedIso()   <= (max_relChargedIso_*l1PFTau->pt())) &&
	     (max_absChargedIso_ < 0. || l1PFTau->sumChargedIso()   <=  max_absChargedIso_               ) )
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
      int numL1PFTausPtGt45 = 0;
      int numL1PFTausPtGt50 = 0;
      for ( std::vector<const l1t::L1HPSPFTau*>::const_iterator l1PFTau = l1PFTaus_passingAbsEta.begin();
	    l1PFTau != l1PFTaus_passingAbsEta.end(); ++l1PFTau ) {
	if ( (*l1PFTau)->pt() > 20. ) ++numL1PFTausPtGt20;
	if ( (*l1PFTau)->pt() > 25. ) ++numL1PFTausPtGt25;
	if ( (*l1PFTau)->pt() > 30. ) ++numL1PFTausPtGt30;
	if ( (*l1PFTau)->pt() > 35. ) ++numL1PFTausPtGt35;
	if ( (*l1PFTau)->pt() > 40. ) ++numL1PFTausPtGt40;
	if ( (*l1PFTau)->pt() > 45. ) ++numL1PFTausPtGt45;
	if ( (*l1PFTau)->pt() > 50. ) ++numL1PFTausPtGt50;
      }
      fillWithOverFlow(histogram_numL1PFTausPtGt20_, numL1PFTausPtGt20, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt25_, numL1PFTausPtGt25, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt30_, numL1PFTausPtGt30, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt35_, numL1PFTausPtGt35, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt40_, numL1PFTausPtGt40, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt45_, numL1PFTausPtGt45, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt50_, numL1PFTausPtGt50, evtWeight);

      TAxis* xAxis = histogram_numL1PFTaus_vs_ptThreshold_->GetXaxis();
      TAxis* yAxis = histogram_numL1PFTaus_vs_ptThreshold_->GetYaxis();
      int numBinsX = xAxis->GetNbins();
      bool max_numL1PFTaus_passingPt_isZero = false;
      for ( int idxBin = 1; idxBin <= numBinsX; ++idxBin )
      {
	double ptThreshold = xAxis->GetBinCenter(idxBin);

	int max_numL1PFTaus_passingPt = 0;
        if ( !max_numL1PFTaus_passingPt_isZero ) 
        {
	  for ( std::vector<const l1t::L1HPSPFTau*>::const_iterator l1PFTau_zVtxRef = l1PFTaus_passingAbsEta.begin(); 
	        l1PFTau_zVtxRef != l1PFTaus_passingAbsEta.end(); ++l1PFTau_zVtxRef ) {
	    int numL1PFTaus_passingPt = 0;
	    for ( std::vector<const l1t::L1HPSPFTau*>::const_iterator l1PFTau_toMatch = l1PFTaus_passingAbsEta.begin(); 
		  l1PFTau_toMatch != l1PFTaus_passingAbsEta.end(); ++l1PFTau_toMatch ) {
	      if ( (*l1PFTau_toMatch)->pt() > ptThreshold )
	      {
	        if ( (*l1PFTau_zVtxRef)->leadChargedPFCand().isNonnull() && (*l1PFTau_zVtxRef)->leadChargedPFCand()->pfTrack().isNonnull() &&
	             (*l1PFTau_toMatch)->leadChargedPFCand().isNonnull() && (*l1PFTau_toMatch)->leadChargedPFCand()->pfTrack().isNonnull() )
  	        {
		  double dz = TMath::Abs((*l1PFTau_zVtxRef)->leadChargedPFCand()->pfTrack()->vertex().z() - (*l1PFTau_toMatch)->leadChargedPFCand()->pfTrack()->vertex().z());
		  if ( dz < max_dz_ ) 
	          {
		    ++numL1PFTaus_passingPt;
		  }
	        }
	      }
	    }
	    if ( numL1PFTaus_passingPt > max_numL1PFTaus_passingPt ) 
	    {
              max_numL1PFTaus_passingPt = numL1PFTaus_passingPt;
            }
	  }
	}

	if ( max_numL1PFTaus_passingPt > yAxis->GetXmax() ) max_numL1PFTaus_passingPt = TMath::Nint(yAxis->GetXmax() - 0.5);
	histogram_numL1PFTaus_vs_ptThreshold_->Fill(ptThreshold, max_numL1PFTaus_passingPt, evtWeight);

	if ( max_numL1PFTaus_passingPt == 0 ) 
	{
	  max_numL1PFTaus_passingPt_isZero = true;
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
    MonitorElement* me_numL1PFTausPtGt45_;
    TH1* histogram_numL1PFTausPtGt45_;
    MonitorElement* me_numL1PFTausPtGt50_;
    TH1* histogram_numL1PFTausPtGt50_;  
    double min_absEta_;    
    double max_absEta_;    
    double max_relChargedIso_;
    double max_absChargedIso_;
    double max_dz_;
  };
  std::vector<ratePlotEntryType*> ratePlots_;
};

#endif   
