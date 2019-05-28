#ifndef L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauAnalyzerBackground_h
#define L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauAnalyzerBackground_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTau.h"    // l1t::TallinnL1PFTau
#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTauFwd.h" // l1t::TallinnL1PFTauCollection

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
      : me_(nullptr)
      , histogram_(nullptr)
      , max_absEta_(max_absEta)
      , max_relChargedIso_(max_relChargedIso)
      , max_absChargedIso_(max_absChargedIso)
      , max_dz_(max_dz)
    {}
    ~ratePlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName = "numL1PFTaus_vs_L1PFTauPt";
      if ( max_absEta_        > 0. ) histogramName.Append(Form("_absEtaLt%1.2f",        max_absEta_));
      if ( max_relChargedIso_ > 0. ) histogramName.Append(Form("_relChargedIsoLt%1.2f", max_relChargedIso_));
      if ( max_absChargedIso_ > 0. ) histogramName.Append(Form("_absChargedIsoLt%1.2f", max_absChargedIso_));
      histogramName = histogramName.ReplaceAll(".", "p");
      me_ = dqmStore.book2D(histogramName.Data(), histogramName.Data(), 101, -0.5, 100.5, 11, -0.5, +10.5);
      histogram_ = dynamic_cast<TH2*>(me_->getTH1());
      assert(histogram_);
    }
    void fillHistograms(const l1t::TallinnL1PFTauCollection& l1PFTaus, double evtWeight)
    {
      std::vector<const l1t::TallinnL1PFTau*> l1PFTaus_passingAbsEta;
      for ( auto l1PFTau : l1PFTaus ) 
      {
	if ( (max_absEta_        < 0. || TMath::Abs(l1PFTau.eta()) <=  max_absEta_                     ) &&
	     (max_relChargedIso_ < 0. || l1PFTau.sumChargedIso()   <= (max_relChargedIso_*l1PFTau.pt())) &&
	     (max_absChargedIso_ < 0. || l1PFTau.sumChargedIso()   <=  max_absChargedIso_              ) )
	{
	  l1PFTaus_passingAbsEta.push_back(&l1PFTau);
	}
      }
      // CV: sort L1PFTaus passing abs(eta) cut by decreasing pT
      std::sort(l1PFTaus_passingAbsEta.begin(), l1PFTaus_passingAbsEta.end(), isHigherPt);

      TAxis* xAxis = histogram_->GetXaxis();
      TAxis* yAxis = histogram_->GetYaxis();
      int numBinsX = xAxis->GetNbins();
      bool numL1PFTaus_isZero = false;
      for ( int idxBin = 1; idxBin <= numBinsX; ++idxBin )
      {
	double ptThreshold = xAxis->GetBinCenter(idxBin);
	int max_numL1PFTaus_passingPt = 0;
	if ( !numL1PFTaus_isZero ) 
	{
	  for ( auto l1PFTau1 : l1PFTaus_passingAbsEta ) 
          {
	    int numL1PFTaus_passingPt = 0;
	    for ( auto l1PFTau2 : l1PFTaus_passingAbsEta )
	    {
	      if ( l1PFTau2->pt() > ptThreshold )
	      {
	        if ( l1PFTau1->leadChargedPFCand().isNonnull() && l1PFTau1->leadChargedPFCand()->pfTrack().isNonnull() &&
		     l1PFTau2->leadChargedPFCand().isNonnull() && l1PFTau2->leadChargedPFCand()->pfTrack().isNonnull() )
		{
		  double dz = TMath::Abs(l1PFTau1->leadChargedPFCand()->pfTrack()->vertex().z() - l1PFTau2->leadChargedPFCand()->pfTrack()->vertex().z());
		  if ( dz < max_dz_ ) ++numL1PFTaus_passingPt;
		}
	      }
	    }
	    if ( numL1PFTaus_passingPt > max_numL1PFTaus_passingPt ) 
	    {
	      max_numL1PFTaus_passingPt = numL1PFTaus_passingPt;
	    }
	  }
	  if ( max_numL1PFTaus_passingPt == 0 ) numL1PFTaus_isZero = true;
	}
	if ( max_numL1PFTaus_passingPt > yAxis->GetXmax() ) max_numL1PFTaus_passingPt = TMath::Nint(yAxis->GetXmax() - 0.5);
	histogram_->Fill(ptThreshold, max_numL1PFTaus_passingPt, evtWeight);
      }
    }
    MonitorElement* me_;
    TH2* histogram_;
    double max_absEta_;    
    double max_relChargedIso_;
    double max_absChargedIso_;
    double max_dz_;
  };
  std::vector<ratePlotEntryType*> ratePlots_;
};

#endif   
