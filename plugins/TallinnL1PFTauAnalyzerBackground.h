#ifndef L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauAnalyzerBackground_h
#define L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauAnalyzerBackground_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTau.h"             // l1t::TallinnL1PFTau
#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTauFwd.h"          // l1t::TallinnL1PFTauCollection
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
  isHigherPt(const l1t::TallinnL1PFTau* l1PFTau1,
	     const l1t::TallinnL1PFTau* l1PFTau2)
  {
    return l1PFTau1->pt() > l1PFTau2->pt();
  }
}

class TallinnL1PFTauAnalyzerBackground : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit TallinnL1PFTauAnalyzerBackground(const edm::ParameterSet&);
    
  // destructor
  ~TallinnL1PFTauAnalyzerBackground();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcTallinnL1PFTaus_;
  edm::EDGetTokenT<l1t::TallinnL1PFTauCollection> tokenTallinnL1PFTaus_;

  std::string dqmDirectory_;

  struct ratePlotEntryType
  {
    ratePlotEntryType(double max_absEta, double max_relChargedIso, double max_absChargedIso, double max_dz)
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
      if ( max_absEta_        > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f",        max_absEta_));
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
    }
    void fillHistograms(const l1t::TallinnL1PFTauCollection& l1PFTaus, double evtWeight)
    {
      std::vector<const l1t::TallinnL1PFTau*> l1PFTaus_passingAbsEta;
      int idx = 0;
      for ( l1t::TallinnL1PFTauCollection::const_iterator l1PFTau = l1PFTaus.begin(); l1PFTau != l1PFTaus.end(); ++l1PFTau ) 
      {
	std::cout << "TallinL1PFTau #" << idx << ":" << (*l1PFTau);
	++idx;
	if ( (max_absEta_        < 0. || TMath::Abs(l1PFTau->eta()) <=  max_absEta_                      ) &&
	     (max_relChargedIso_ < 0. || l1PFTau->sumChargedIso()   <= (max_relChargedIso_*l1PFTau->pt())) &&
	     (max_absChargedIso_ < 0. || l1PFTau->sumChargedIso()   <=  max_absChargedIso_               ) )
	{
	  std::cout << " passinggAbsEta cuts." << std::endl;
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
      for ( auto l1PFTau : l1PFTaus_passingAbsEta ) 
      {
	if ( l1PFTau->pt() > 20. ) ++numL1PFTausPtGt20;
	if ( l1PFTau->pt() > 25. ) ++numL1PFTausPtGt25;
	if ( l1PFTau->pt() > 30. ) ++numL1PFTausPtGt30;
	if ( l1PFTau->pt() > 35. ) ++numL1PFTausPtGt35;
	if ( l1PFTau->pt() > 40. ) ++numL1PFTausPtGt40;	
      }
      fillWithOverFlow(histogram_numL1PFTausPtGt20_, numL1PFTausPtGt20, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt25_, numL1PFTausPtGt25, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt30_, numL1PFTausPtGt30, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt35_, numL1PFTausPtGt35, evtWeight);
      fillWithOverFlow(histogram_numL1PFTausPtGt40_, numL1PFTausPtGt40, evtWeight);

      TAxis* xAxis = histogram_numL1PFTaus_vs_ptThreshold_->GetXaxis();
      TAxis* yAxis = histogram_numL1PFTaus_vs_ptThreshold_->GetYaxis();
      int numBinsX = xAxis->GetNbins();
      bool max_numL1PFTaus_passingPt_isZero = false;
      for ( int idxBin = 1; idxBin <= numBinsX; ++idxBin )
      {
	double ptThreshold = xAxis->GetBinCenter(idxBin);
	std::cout << "ptThreshold = " << ptThreshold << std::endl;
	int max_numL1PFTaus_passingPt = 0;
        if ( !max_numL1PFTaus_passingPt_isZero ) {
	int idx_zVtxRef = 0;
	for ( std::vector<const l1t::TallinnL1PFTau*>::const_iterator l1PFTau_zVtxRef = l1PFTaus_passingAbsEta.begin(); 
	      l1PFTau_zVtxRef != l1PFTaus_passingAbsEta.end(); ++l1PFTau_zVtxRef ) {
	  std::cout << "zVtxRef TallinL1PFTau #" << idx_zVtxRef << ":" << std::endl;
	  int numL1PFTaus_passingPt = 0;
	  int idx_test = 0;
	  for ( std::vector<const l1t::TallinnL1PFTau*>::const_iterator l1PFTau_test = l1PFTaus_passingAbsEta.begin(); 
		l1PFTau_test != l1PFTaus_passingAbsEta.end(); ++l1PFTau_test ) {
	    if ( (*l1PFTau_test)->pt() > ptThreshold )
	    {
	      if ( (*l1PFTau_zVtxRef)->leadChargedPFCand().isNonnull() && (*l1PFTau_zVtxRef)->leadChargedPFCand()->pfTrack().isNonnull() &&
	           (*l1PFTau_test)->leadChargedPFCand().isNonnull() && (*l1PFTau_test)->leadChargedPFCand()->pfTrack().isNonnull() )
	      {
		double dz = TMath::Abs((*l1PFTau_zVtxRef)->leadChargedPFCand()->pfTrack()->vertex().z() - (*l1PFTau_test)->leadChargedPFCand()->pfTrack()->vertex().z());
		std::cout << "idx_test = " << idx_test << ": dz = " << dz << std::endl;
		if ( dz < max_dz_ ) 
	        {
		  std::cout << " matching test TallinL1PFTau #" << idx_test << std::endl;
		  ++numL1PFTaus_passingPt;
		}
	      }
	    }
	    ++idx_test;
	  }
	  std::cout << "--> numL1PFTaus_passingPt = " << numL1PFTaus_passingPt << std::endl;
	  if ( numL1PFTaus_passingPt > max_numL1PFTaus_passingPt ) 
	  {
            max_numL1PFTaus_passingPt = numL1PFTaus_passingPt;
          }
	  ++idx_zVtxRef;
	}
	}
	std::cout << "ptThreshold = " << ptThreshold << ": max_numL1PFTaus_passingPt = " << max_numL1PFTaus_passingPt << std::endl;
	if ( (ptThreshold < 20. && max_numL1PFTaus_passingPt < numL1PFTausPtGt20) ||
	     (ptThreshold < 25. && max_numL1PFTaus_passingPt < numL1PFTausPtGt25) ||
	     (ptThreshold < 30. && max_numL1PFTaus_passingPt < numL1PFTausPtGt30) ||
	     (ptThreshold < 35. && max_numL1PFTaus_passingPt < numL1PFTausPtGt35) ||
	     (ptThreshold < 40. && max_numL1PFTaus_passingPt < numL1PFTausPtGt40) )
	{
	  std::cout << "Internal logic error in <ratePlotEntryType::fillHistograms>:" << std::endl;
	  std::cout << " max_numL1PFTaus_passingPt = " << max_numL1PFTaus_passingPt << " @ ptThreshold = " << ptThreshold << ","
		    << " while numL1PFTausPtGt20/25/30/35/40 = " << numL1PFTausPtGt20 << "/" << numL1PFTausPtGt25 
		    << "/" << numL1PFTausPtGt30 << "/" << numL1PFTausPtGt35 << "/" << numL1PFTausPtGt40 << " --> CHECK !!" << std::endl;
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
    double max_absEta_;    
    double max_relChargedIso_;
    double max_absChargedIso_;
    double max_dz_;
  };
  std::vector<ratePlotEntryType*> ratePlots_;
};

#endif   
