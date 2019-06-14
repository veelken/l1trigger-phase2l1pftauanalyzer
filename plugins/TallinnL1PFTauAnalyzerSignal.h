#ifndef L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauAnalyzerSignal_h
#define L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauAnalyzerSignal_h

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

class TallinnL1PFTauAnalyzerSignal : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit TallinnL1PFTauAnalyzerSignal(const edm::ParameterSet&);
    
  // destructor
  ~TallinnL1PFTauAnalyzerSignal();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcNumerator_;
  edm::EDGetTokenT<l1t::TallinnL1PFTauCollection> tokenNumerator_;
  edm::InputTag srcDenominator_;
  enum { kGen, kOffline };
  int typeDenominator_;
  edm::EDGetTokenT<reco::GenJetCollection> tokenDenominator_gen_;
  edm::EDGetTokenT<pat::TauCollection> tokenDenominator_offline_;

  std::string dqmDirectory_;

  struct efficiencyPlotEntryType
  {
    efficiencyPlotEntryType(double min_pt, double max_pt, double min_absEta, double max_absEta, const std::string& decayMode, 
			    double ptThreshold, double max_relChargedIso, double max_absChargedIso)
      : me_pt_numerator_(nullptr)
      , histogram_pt_numerator_(nullptr)
      , me_pt_denominator_(nullptr)
      , histogram_pt_denominator_(nullptr)
      , me_eta_numerator_(nullptr)
      , histogram_eta_numerator_(nullptr)
      , me_eta_denominator_(nullptr)
      , histogram_eta_denominator_(nullptr)
      , me_phi_numerator_(nullptr)
      , histogram_phi_numerator_(nullptr)
      , me_phi_denominator_(nullptr)
      , histogram_phi_denominator_(nullptr)
      , min_pt_(min_pt)
      , max_pt_(max_pt)
      , min_absEta_(min_absEta)
      , max_absEta_(max_absEta)
      , decayMode_(decayMode)
      , ptThreshold_(ptThreshold)
      , max_relChargedIso_(max_relChargedIso)
      , max_absChargedIso_(max_absChargedIso)
      , dRmatch_(0.3)
    {}
    ~efficiencyPlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix = decayMode_.data(); 
      if      ( min_absEta_ >= 0. && max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEta%1.2fto%1.2f", min_absEta_, max_absEta_));
      else if ( min_absEta_ >= 0.                     ) histogramName_suffix.Append(Form("_absEtaGt%1.2f", min_absEta_));
      else if (                      max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f", max_absEta_));
      if ( max_relChargedIso_ > 0. ) histogramName_suffix.Append(Form("_relChargedIsoLt%1.2f", max_relChargedIso_));
      if ( max_absChargedIso_ > 0. ) histogramName_suffix.Append(Form("_absChargedIsoLt%1.2f", max_absChargedIso_));
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_pt_numerator = Form("effL1PFTau_vs_pt_numerator_%s", histogramName_suffix.Data());
      me_pt_numerator_ = dqmStore.book1D(histogramName_pt_numerator.Data(), histogramName_pt_numerator.Data(), 16, 20., 100.);
      histogram_pt_numerator_ = me_pt_numerator_->getTH1();
      assert(histogram_pt_numerator_);
      TString histogramName_pt_denominator = Form("effL1PFTau_vs_pt_denominator_%s", histogramName_suffix.Data());
      me_pt_denominator_ = dqmStore.book1D(histogramName_pt_denominator.Data(), histogramName_pt_denominator.Data(), 16, 20., 100.);
      histogram_pt_denominator_ = me_pt_denominator_->getTH1();
      assert(histogram_pt_denominator_);

      TString histogramName_eta_numerator = Form("effL1PFTau_vs_eta_numerator_%s", histogramName_suffix.Data());
      me_eta_numerator_ = dqmStore.book1D(histogramName_eta_numerator.Data(), histogramName_eta_numerator.Data(), 30, -3., +3.);
      histogram_eta_numerator_ = me_eta_numerator_->getTH1();
      assert(histogram_eta_numerator_);
      TString histogramName_eta_denominator = Form("effL1PFTau_vs_eta_denominator_%s", histogramName_suffix.Data());
      me_eta_denominator_ = dqmStore.book1D(histogramName_eta_denominator.Data(), histogramName_eta_denominator.Data(), 30, -3., +3.);
      histogram_eta_denominator_ = me_eta_denominator_->getTH1();
      assert(histogram_eta_denominator_);

      TString histogramName_phi_numerator = Form("effL1PFTau_vs_phi_numerator_%s", histogramName_suffix.Data());
      me_phi_numerator_ = dqmStore.book1D(histogramName_phi_numerator.Data(), histogramName_phi_numerator.Data(), 18, -TMath::Pi(), +TMath::Pi());
      histogram_phi_numerator_ = me_phi_numerator_->getTH1();
      assert(histogram_phi_numerator_);
      TString histogramName_phi_denominator = Form("effL1PFTau_vs_phi_denominator_%s", histogramName_suffix.Data());
      me_phi_denominator_ = dqmStore.book1D(histogramName_phi_denominator.Data(), histogramName_phi_denominator.Data(), 18, -TMath::Pi(), +TMath::Pi());
      histogram_phi_denominator_ = me_phi_denominator_->getTH1();
      assert(histogram_phi_denominator_);
    }
    void fillHistograms(const l1t::TallinnL1PFTauCollection& numeratorTaus, const reco::GenJetCollection& denominatorTaus, double evtWeight)
    {
      for ( auto denominatorTau : denominatorTaus )
      {
	std::string denominatorTau_decayMode = JetMCTagUtils::genTauDecayMode(denominatorTau);
	if ( decayMode_ != "all" && denominatorTau_decayMode != decayMode_ ) continue;

	bool isMatched = false;
	for ( auto numeratorTau : numeratorTaus )
	{
  	  if ( (ptThreshold_       < 0. || numeratorTau.pt()            >=  ptThreshold_                         ) &&
	       (                           numeratorTau.leadChargedPFCand().isNonnull()                          && 
					   numeratorTau.leadChargedPFCand()->pfTrack().isNonnull()               ) &&
               (max_relChargedIso_ < 0. || numeratorTau.sumChargedIso() <= (max_relChargedIso_*numeratorTau.pt())) &&
	       (max_absChargedIso_ < 0. || numeratorTau.sumChargedIso() <=  max_absChargedIso_                   ) )
	  {
	    double dR = reco::deltaR(denominatorTau.eta(), denominatorTau.phi(), numeratorTau.eta(), numeratorTau.phi());
	    if ( dR < dRmatch_ ) isMatched = true;
	  }
	}

	double denominatorTau_absEta = TMath::Abs(denominatorTau.eta());
	if ( denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ ) 
	{
	  histogram_pt_denominator_->Fill(denominatorTau.pt(), evtWeight);
	  if ( isMatched ) histogram_pt_numerator_->Fill(denominatorTau.pt(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_ && denominatorTau.pt() < max_pt_ ) 
	{
	  histogram_eta_denominator_->Fill(denominatorTau.eta(), evtWeight);
	  if ( isMatched ) histogram_eta_numerator_->Fill(denominatorTau.eta(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_ && denominatorTau.pt() < max_pt_ && denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ )
	{
	  histogram_phi_denominator_->Fill(denominatorTau.phi(), evtWeight);
	  if ( isMatched ) histogram_phi_numerator_->Fill(denominatorTau.phi(), evtWeight);
	}
      }
    }
    void fillHistograms(const l1t::TallinnL1PFTauCollection& numeratorTaus, const pat::TauCollection& denominatorTaus, double evtWeight)
    {
      for ( auto denominatorTau : denominatorTaus )
      {
	if ( (decayMode_ == "oneProng0Pi0"   && denominatorTau.decayMode() != reco::PFTau::kOneProng0PiZero)   ||
	     (decayMode_ == "oneProng1Pi0"   && denominatorTau.decayMode() != reco::PFTau::kOneProng1PiZero)   ||
	     (decayMode_ == "oneProng1Pi0"   && denominatorTau.decayMode() != reco::PFTau::kOneProng2PiZero)   ||
	     (decayMode_ == "threeProng0Pi0" && denominatorTau.decayMode() != reco::PFTau::kThreeProng0PiZero) ||
	     (decayMode_ == "threeProng1Pi0" && denominatorTau.decayMode() != reco::PFTau::kThreeProng1PiZero) ) continue;	       

	bool isMatched = false;
	for ( auto numeratorTau : numeratorTaus )
	{
  	  if ( (ptThreshold_       < 0. || numeratorTau.pt()            >=  ptThreshold_                         ) &&
	       (                           numeratorTau.leadChargedPFCand().isNonnull()                          && 
					   numeratorTau.leadChargedPFCand()->pfTrack().isNonnull()               ) &&
               (max_relChargedIso_ < 0. || numeratorTau.sumChargedIso() <= (max_relChargedIso_*numeratorTau.pt())) &&
	       (max_absChargedIso_ < 0. || numeratorTau.sumChargedIso() <=  max_absChargedIso_                   ) )
	  {
	    double dR = reco::deltaR(denominatorTau.eta(), denominatorTau.phi(), numeratorTau.eta(), numeratorTau.phi());
	    if ( dR < dRmatch_ ) isMatched = true;
	  }
	}

	double denominatorTau_absEta = TMath::Abs(denominatorTau.eta());
	if ( denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ ) 
	{
	  histogram_pt_denominator_->Fill(denominatorTau.pt(), evtWeight);
	  if ( isMatched ) histogram_pt_numerator_->Fill(denominatorTau.pt(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_ && denominatorTau.pt() < max_pt_ ) 
	{
	  histogram_eta_denominator_->Fill(denominatorTau.eta(), evtWeight);
	  if ( isMatched ) histogram_eta_numerator_->Fill(denominatorTau.eta(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_ && denominatorTau.pt() < max_pt_ && denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ )
	{
	  histogram_phi_denominator_->Fill(denominatorTau.phi(), evtWeight);
	  if ( isMatched ) histogram_phi_numerator_->Fill(denominatorTau.phi(), evtWeight);
	}
      }
    }
    MonitorElement* me_pt_numerator_;
    TH1* histogram_pt_numerator_;
    MonitorElement* me_pt_denominator_;
    TH1* histogram_pt_denominator_;
    MonitorElement* me_eta_numerator_;
    TH1* histogram_eta_numerator_;
    MonitorElement* me_eta_denominator_;
    TH1* histogram_eta_denominator_;
    MonitorElement* me_phi_numerator_;
    TH1* histogram_phi_numerator_;
    MonitorElement* me_phi_denominator_;
    TH1* histogram_phi_denominator_;
    // cuts applied to offline and generator-level taus in denominator 
    double min_pt_;
    double max_pt_;
    double min_absEta_;
    double max_absEta_;
    std::string decayMode_;
    // cuts applied to L1 trigger taus in numerator 
    double ptThreshold_;
    double max_relChargedIso_;
    double max_absChargedIso_;
    // matching between offline or generator-level taus in denominator and L1 trigger taus in numerator
    double dRmatch_; 
  };
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_;
};

#endif   
