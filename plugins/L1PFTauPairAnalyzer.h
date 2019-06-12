#ifndef L1Trigger_TallinnL1PFTauAnalyzer_L1PFTauPairAnalyzer_h
#define L1Trigger_TallinnL1PFTauAnalyzer_L1PFTauPairAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/L1Trigger/interface/L1PFTau.h"                          // l1t::L1PFTau, l1t::L1PFTauCollection
#include "DataFormats/Candidate/interface/CandidateFwd.h"                     // reco::CandidateView
#include "DataFormats/Candidate/interface/Candidate.h"                        // reco::Candidate
#include "L1Trigger/TallinnL1PFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverFlow()

#include <TH1.h>     // TH1
#include <TH2.h>     // TH2
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs(), TMath::Nint()

#include <vector>    // std::vector
#include <string>    // std::string
#include <algorithm> // std::sort

namespace l1t
{
  class L1PFTauPair
  {
   public:
    L1PFTauPair(const l1t::L1PFTau* leadL1PFTau, const l1t::L1PFTau* subleadL1PFTau)
      : leadL1PFTau_(leadL1PFTau)
      , subleadL1PFTau_(subleadL1PFTau)
    {}
    ~L1PFTauPair() 
    {}
    const l1t::L1PFTau* leadL1PFTau()    const { return leadL1PFTau_;    }
    const l1t::L1PFTau* subleadL1PFTau() const { return subleadL1PFTau_; }
   private:
    const l1t::L1PFTau* leadL1PFTau_;
    const l1t::L1PFTau* subleadL1PFTau_;
  };

  typedef std::vector<L1PFTauPair> L1PFTauPairCollection;
}

class L1PFTauPairAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit L1PFTauPairAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~L1PFTauPairAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcL1PFTaus_;
  edm::EDGetTokenT<l1t::L1PFTauCollection> tokenL1PFTaus_;
  edm::InputTag srcRefTaus_;
  edm::EDGetTokenT<reco::CandidateView> tokenRefTaus_;

  double min_refTau_pt_;
  double max_refTau_eta_;
  double dRmatch_;

  std::string dqmDirectory_;

  struct efficiency_or_ratePlotEntryType
  {
    efficiency_or_ratePlotEntryType(double max_absEta, const std::string& isolation)
      : me_efficiency_or_rate_(nullptr)
      , histogram_efficiency_or_rate_(nullptr)
      , me_denominator_(nullptr)
      , histogram_denominator_(nullptr)
      , max_absEta_(max_absEta)
      , isolation_string_(isolation)
    {
      if      ( isolation_string_ == ""       ) isolation_ = kNone;
      else if ( isolation_string_ == "vLoose" ) isolation_ = kVLooseIso;
      else if ( isolation_string_ == "Loose"  ) isolation_ = kLooseIso;
      else if ( isolation_string_ == "Medium" ) isolation_ = kMediumIso;
      else if ( isolation_string_ == "Tight"  ) isolation_ = kTightIso; 
      else throw cms::Exception("L1PFTauPairAnalyzer") 
	<< " Invalid Configuration parameter 'isolation' = " << isolation_string_ << " !!\n";;
    }
    ~efficiency_or_ratePlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix;
      if ( max_absEta_        > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f", max_absEta_));
      if ( isolation_string_ != "" ) histogramName_suffix.Append(Form("_%sIso",         isolation_string_.data()));
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_efficiency_or_rate = Form("efficiency_or_rate%s", histogramName_suffix.Data());
      me_efficiency_or_rate_ = dqmStore.book2D(histogramName_efficiency_or_rate.Data(), histogramName_efficiency_or_rate.Data(), 101, -0.5, 100.5, 101, -0.5, 100.5);
      histogram_efficiency_or_rate_ = dynamic_cast<TH2*>(me_efficiency_or_rate_->getTH1());
      assert(histogram_efficiency_or_rate_);

      TString histogramName_denominator = Form("denominator%s", histogramName_suffix.Data());
      me_denominator_ = dqmStore.book1D(histogramName_denominator.Data(), histogramName_denominator.Data(), 1, -0.5, +0.5);
      histogram_denominator_ = me_denominator_->getTH1();
      assert(histogram_denominator_);
    }
    void fillHistograms(const l1t::L1PFTauPairCollection& l1PFTauPairs, double evtWeight)
    {
      std::vector<const l1t::L1PFTauPair*> l1PFTauPairs_passingAbsEta;
      for ( l1t::L1PFTauPairCollection::const_iterator l1PFTauPair = l1PFTauPairs.begin(); 
	    l1PFTauPair != l1PFTauPairs.end(); ++l1PFTauPair ) 
      {
	const l1t::L1PFTau* leadL1PFTau    = l1PFTauPair->leadL1PFTau();
	const l1t::L1PFTau* subleadL1PFTau = l1PFTauPair->subleadL1PFTau();
	if ( ( max_absEta_ < 0.                                                                             || 
	      (TMath::Abs(leadL1PFTau->eta())    <=  max_absEta_                                          && 
	       TMath::Abs(subleadL1PFTau->eta()) <=  max_absEta_                                          ) ) &&
	     ( isolation_ == kNone                                                                          || 
	      (isolation_ == kVLooseIso && leadL1PFTau->passVLooseIso() && subleadL1PFTau->passVLooseIso()) ||
	      (isolation_ == kLooseIso  && leadL1PFTau->passLooseIso()  && subleadL1PFTau->passLooseIso())  ||
	      (isolation_ == kMediumIso && leadL1PFTau->passMediumIso() && subleadL1PFTau->passMediumIso()) ||
	      (isolation_ == kTightIso  && leadL1PFTau->passTightIso()  && subleadL1PFTau->passTightIso())  ) )
	{
	  l1PFTauPairs_passingAbsEta.push_back(&(*l1PFTauPair));
	}
      }

      TAxis* xAxis = histogram_efficiency_or_rate_->GetXaxis();
      int numBinsX = xAxis->GetNbins();
      TAxis* yAxis = histogram_efficiency_or_rate_->GetYaxis();
      int numBinsY = yAxis->GetNbins();
      bool max_numL1PFTauPairs_passingPt_isZero = false;
      for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX )
      {
	double leadL1PFTau_ptThreshold = xAxis->GetBinCenter(idxBinX);
	int max_numL1PFTauPairs_passingPt = 0;
	if ( !max_numL1PFTauPairs_passingPt_isZero ) 
        {
	  for ( int idxBinY = 1; idxBinY <= numBinsY; ++idxBinY )
          {
	    double subleadL1PFTau_ptThreshold = yAxis->GetBinCenter(idxBinY);
	    int numL1PFTauPairs_passingPt = 0;
	    for ( std::vector<const l1t::L1PFTauPair*>::const_iterator l1PFTauPair = l1PFTauPairs_passingAbsEta.begin(); 
		  l1PFTauPair != l1PFTauPairs_passingAbsEta.end(); ++l1PFTauPair ) {
	      if ( (*l1PFTauPair)->leadL1PFTau()->pt()    > leadL1PFTau_ptThreshold    &&
		   (*l1PFTauPair)->subleadL1PFTau()->pt() > subleadL1PFTau_ptThreshold )
	      {
	        ++numL1PFTauPairs_passingPt;
	      }
	    }
	    if ( numL1PFTauPairs_passingPt >= 1 )
	    {
	      histogram_efficiency_or_rate_->Fill(leadL1PFTau_ptThreshold, subleadL1PFTau_ptThreshold, evtWeight);
	    }
	    if ( numL1PFTauPairs_passingPt > max_numL1PFTauPairs_passingPt ) 
	    {
	      max_numL1PFTauPairs_passingPt = numL1PFTauPairs_passingPt;
	    }
  	  }
	  if ( max_numL1PFTauPairs_passingPt == 0 ) 
	  {
	    max_numL1PFTauPairs_passingPt_isZero = true;
	  }
	}
      }

      histogram_denominator_->Fill(0., evtWeight);
    }
    void normalizeHistograms()
    {
      if ( histogram_denominator_->Integral() > 0. ) 
      {
	histogram_efficiency_or_rate_->Scale(1./histogram_denominator_->Integral());
      }
    }
    MonitorElement* me_efficiency_or_rate_;
    TH2* histogram_efficiency_or_rate_;
    MonitorElement* me_denominator_;
    TH1* histogram_denominator_;
    double max_absEta_;    
    std::string isolation_string_;
    enum { kNone, kVLooseIso, kLooseIso, kMediumIso, kTightIso };
    int isolation_;
  };
  std::vector<efficiency_or_ratePlotEntryType*> efficiency_or_ratePlots_;
};

#endif   
