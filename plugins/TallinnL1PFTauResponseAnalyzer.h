#ifndef L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauResponseAnalyzer_h
#define L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauResponseAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTau.h"    // l1t::TallinnL1PFTau
#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTauFwd.h" // l1t::TallinnL1PFTauCollection
#include "DataFormats/Math/interface/deltaR.h"                       // reco::deltaR
#include "DataFormats/JetReco/interface/GenJet.h"                    // reco::GenJet
#include "DataFormats/JetReco/interface/GenJetCollection.h"          // reco::GenJetCollection
#include "DataFormats/PatCandidates/interface/Tau.h"                 // pat::Tau, pat::TauCollection
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"              // JetMCTagUtils::genTauDecayMode()
#include "DataFormats/TauReco/interface/PFTau.h"                     // reco::PFTau::kOneProng0PiZero, reco::PFTau::kOneProng1PiZero,...

#include <TH1.h>     // TH1
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs(), TMath::Pi()

#include <vector>    // std::vector
#include <string>    // std::string

class TallinnL1PFTauResponseAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit TallinnL1PFTauResponseAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~TallinnL1PFTauResponseAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcL1PFTaus_;
  edm::EDGetTokenT<l1t::TallinnL1PFTauCollection> tokenL1PFTaus_;
  edm::InputTag srcRefTaus_;
  enum { kGen, kOffline };
  int typeRefTaus_;
  edm::EDGetTokenT<reco::GenJetCollection> tokenRefTaus_gen_;
  edm::EDGetTokenT<pat::TauCollection> tokenRefTaus_offline_;

  std::string dqmDirectory_;

  struct responsePlotEntryType
  {
    responsePlotEntryType(double min_pt, double max_absEta, const std::string& decayMode, double max_relChargedIso, double max_absChargedIso)
      : me_response_(nullptr)
      , histogram_response_(nullptr)
      , min_pt_(min_pt)
      , max_absEta_(max_absEta)
      , decayMode_(decayMode)
      , max_relChargedIso_(max_relChargedIso)
      , max_absChargedIso_(max_absChargedIso)
      , dRmatch_(0.3)
    {}
    ~responsePlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix = decayMode_.data();      
      if ( max_absEta_        > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f",        max_absEta_));
      if ( min_pt_            > 0. ) histogramName_suffix.Append(Form("_ptGt%1.0f",            min_pt_));
      if ( max_relChargedIso_ > 0. ) histogramName_suffix.Append(Form("_relChargedIsoLt%1.2f", max_relChargedIso_));
      if ( max_absChargedIso_ > 0. ) histogramName_suffix.Append(Form("_absChargedIsoLt%1.2f", max_absChargedIso_));
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_response = Form("response_%s", histogramName_suffix.Data());
      me_response_ = dqmStore.book1D(histogramName_response.Data(), histogramName_response.Data(), 200, 0., 2.);
      histogram_response_ = me_response_->getTH1();
      assert(histogram_response_);
    }
    void fillHistograms(const l1t::TallinnL1PFTauCollection& l1PFTaus, const reco::GenJetCollection& refTaus, double evtWeight)
    {
      for ( auto refTau : refTaus )
      {
	if ( refTau.pt() > min_pt_ && TMath::Abs(refTau.eta()) < max_absEta_ ) 
	{
  	  std::string refTau_decayMode = JetMCTagUtils::genTauDecayMode(refTau);
	  if ( decayMode_ != "all" && refTau_decayMode != decayMode_ ) continue;

	  for ( auto l1PFTau : l1PFTaus )
	  {
  	    if ( (                           l1PFTau.leadChargedPFCand().isNonnull()                     && 
		  			     l1PFTau.leadChargedPFCand()->pfTrack().isNonnull()          ) &&
               (max_relChargedIso_ < 0. ||   l1PFTau.sumChargedIso() <= (max_relChargedIso_*l1PFTau.pt())) &&
	       (max_absChargedIso_ < 0. ||   l1PFTau.sumChargedIso() <=  max_absChargedIso_              ) )
	    {
	      double dR = reco::deltaR(refTau.eta(), refTau.phi(), l1PFTau.eta(), l1PFTau.phi());
	      if ( dR < dRmatch_ ) 
	      {
  	        histogram_response_->Fill(l1PFTau.pt()/TMath::Max(1., refTau.pt()), evtWeight);
	      }
	    }
	  }
	}
      }
    }
    void fillHistograms(const l1t::TallinnL1PFTauCollection& l1PFTaus, const pat::TauCollection& refTaus, double evtWeight)
    {
      for ( auto refTau : refTaus )
      {
	if ( refTau.pt() > min_pt_ && TMath::Abs(refTau.eta()) < max_absEta_ ) 
	{
	  if ( (decayMode_ == "oneProng0Pi0"   && refTau.decayMode() != reco::PFTau::kOneProng0PiZero)   ||
	       (decayMode_ == "oneProng1Pi0"   && refTau.decayMode() != reco::PFTau::kOneProng1PiZero)   ||
	       (decayMode_ == "oneProng1Pi0"   && refTau.decayMode() != reco::PFTau::kOneProng2PiZero)   ||
	       (decayMode_ == "threeProng0Pi0" && refTau.decayMode() != reco::PFTau::kThreeProng0PiZero) ||
	       (decayMode_ == "threeProng1Pi0" && refTau.decayMode() != reco::PFTau::kThreeProng1PiZero) ) continue;	       

	  for ( auto l1PFTau : l1PFTaus )
	  {
	    if ( (                           l1PFTau.leadChargedPFCand().isNonnull()                     && 
		  			     l1PFTau.leadChargedPFCand()->pfTrack().isNonnull()          ) &&
               (max_relChargedIso_ < 0. ||   l1PFTau.sumChargedIso() <= (max_relChargedIso_*l1PFTau.pt())) &&
	       (max_absChargedIso_ < 0. ||   l1PFTau.sumChargedIso() <=  max_absChargedIso_              ) )
	    {
	      double dR = reco::deltaR(refTau.eta(), refTau.phi(), l1PFTau.eta(), l1PFTau.phi());
	      if ( dR < dRmatch_ ) 
	      {
  	        histogram_response_->Fill(l1PFTau.pt()/TMath::Max(1., refTau.pt()), evtWeight);
	      }
	    }
	  }
	}
      }
    }
    MonitorElement* me_response_;
    TH1* histogram_response_;
    // cuts applied to offline and generator-level taus
    double min_pt_;
    double max_absEta_;
    std::string decayMode_;
    // cuts applied to L1 trigger taus
    double max_relChargedIso_;
    double max_absChargedIso_;
    // matching between L1 trigger taus and either offline or generator-level taus
    double dRmatch_; 
  };
  std::vector<responsePlotEntryType*> responsePlots_;
};

#endif   
