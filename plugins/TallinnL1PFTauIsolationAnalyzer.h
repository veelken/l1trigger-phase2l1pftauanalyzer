#ifndef L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauIsolationAnalyzer_h
#define L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauIsolationAnalyzer_h

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
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"              // JetMCTagUtils::genTauDecayMode()

#include <TH1.h>
#include <TString.h> // TString, Form()
#include <TMath.h> // TMath::Abs(), TMath::Pi()

#include <vector>
#include <string>

class TallinnL1PFTauIsolationAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit TallinnL1PFTauIsolationAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~TallinnL1PFTauIsolationAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<l1t::TallinnL1PFTauCollection> token_;
  edm::InputTag srcGenTaus_;
  edm::EDGetTokenT<reco::GenJetCollection> tokenGenTaus_;

  std::string dqmDirectory_;

  struct isolationPlotEntryType
  {
    isolationPlotEntryType(double ptThreshold, double max_absEta, const std::string& decayMode)
      : me_absChargedIso_(nullptr)
      , histogram_absChargedIso_(nullptr)
      , me_absNeutralIso_(nullptr)
      , histogram_absNeutralIso_(nullptr)
      , me_absCombinedIso_(nullptr)
      , histogram_absCombinedIso_(nullptr)
      , me_relChargedIso_(nullptr)
      , histogram_relChargedIso_(nullptr)
      , me_relNeutralIso_(nullptr)
      , histogram_relNeutralIso_(nullptr)
      , me_relCombinedIso_(nullptr)
      , histogram_relCombinedIso_(nullptr)
      , ptThreshold_(ptThreshold)
      , max_absEta_(max_absEta)
      , decayMode_(decayMode)
      , dRmatch_(0.3)
    {}
    ~isolationPlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix = decayMode_.data();      
      if ( max_absEta_        > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f", max_absEta_));
      if ( ptThreshold_ > 0.       ) histogramName_suffix.Append(Form("_ptGt%1.0f",     ptThreshold_));
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_absChargedIso = Form("absChargedIso_%s", histogramName_suffix.Data());
      me_absChargedIso_ = dqmStore.book1D(histogramName_absChargedIso.Data(), histogramName_absChargedIso.Data(), 250, 0., 25.);
      histogram_absChargedIso_ = me_absChargedIso_->getTH1();
      assert(histogram_absChargedIso_);

      TString histogramName_absNeutralIso = Form("absNeutralIso_%s", histogramName_suffix.Data());
      me_absNeutralIso_ = dqmStore.book1D(histogramName_absNeutralIso.Data(), histogramName_absNeutralIso.Data(), 250, 0., 25.);
      histogram_absNeutralIso_ = me_absNeutralIso_->getTH1();
      assert(histogram_absNeutralIso_);

      TString histogramName_absCombinedIso = Form("absCombinedIso_%s", histogramName_suffix.Data());
      me_absCombinedIso_ = dqmStore.book1D(histogramName_absCombinedIso.Data(), histogramName_absCombinedIso.Data(), 250, 0., 25.);
      histogram_absCombinedIso_ = me_absCombinedIso_->getTH1();
      assert(histogram_absCombinedIso_);

      TString histogramName_relChargedIso = Form("relChargedIso_%s", histogramName_suffix.Data());
      me_relChargedIso_ = dqmStore.book1D(histogramName_relChargedIso.Data(), histogramName_relChargedIso.Data(), 100, 0., 1.);
      histogram_relChargedIso_ = me_relChargedIso_->getTH1();
      assert(histogram_relChargedIso_);

      TString histogramName_relNeutralIso = Form("relNeutralIso_%s", histogramName_suffix.Data());
      me_relNeutralIso_ = dqmStore.book1D(histogramName_relNeutralIso.Data(), histogramName_relNeutralIso.Data(), 100, 0., 1.);
      histogram_relNeutralIso_ = me_relNeutralIso_->getTH1();
      assert(histogram_relNeutralIso_);

      TString histogramName_relCombinedIso = Form("relCombinedIso_%s", histogramName_suffix.Data());
      me_relCombinedIso_ = dqmStore.book1D(histogramName_relCombinedIso.Data(), histogramName_relCombinedIso.Data(), 100, 0., 1.);
      histogram_relCombinedIso_ = me_relCombinedIso_->getTH1();
      assert(histogram_relCombinedIso_);

      TString histogramName_pt = Form("pt_%s", histogramName_suffix.Data());
      me_pt_ = dqmStore.book1D(histogramName_pt.Data(), histogramName_pt.Data(), 250, 0., 250.);
      histogram_pt_ = me_pt_->getTH1();
      assert(histogram_pt_);
    }
    void fillHistograms(const l1t::TallinnL1PFTau& l1Tau, double evtWeight)
    {
      if ( l1Tau.pt() > ptThreshold_ && TMath::Abs(l1Tau.eta()) < max_absEta_ ) 
      {
	double sumChargedIso = l1Tau.sumChargedIso();
	histogram_absChargedIso_->Fill(sumChargedIso, evtWeight);
	histogram_relChargedIso_->Fill(sumChargedIso/TMath::Max(1., l1Tau.pt()), evtWeight);
	double sumNeutralIso = l1Tau.sumNeutralIso();
	histogram_absNeutralIso_->Fill(sumNeutralIso, evtWeight);
	histogram_relNeutralIso_->Fill(sumNeutralIso/TMath::Max(1., l1Tau.pt()), evtWeight);
	double sumCombinedIso = sumChargedIso + sumNeutralIso;
	histogram_absCombinedIso_->Fill(sumCombinedIso, evtWeight);
	histogram_relCombinedIso_->Fill(sumCombinedIso/TMath::Max(1., l1Tau.pt()), evtWeight);
	histogram_pt_->Fill(l1Tau.pt(), evtWeight);
      }
    }
    void fillHistograms_woGenMatching(const l1t::TallinnL1PFTauCollection& l1Taus, double evtWeight)
    {
      for ( auto l1Tau : l1Taus )
      {
	fillHistograms(l1Tau, evtWeight);
      }
    }
    void fillHistograms_wGenMatching(const l1t::TallinnL1PFTauCollection& l1Taus, const reco::GenJetCollection& genTaus, double evtWeight)
    {
      for ( auto l1Tau : l1Taus )
      {
	std::string denominatorTau_decayMode;
	bool isMatched = false;
	for ( auto genTau : genTaus )
        {
	  double dR = reco::deltaR(genTau.eta(), genTau.phi(), l1Tau.eta(), l1Tau.phi());
	  if ( dR < dRmatch_ ) 
	  { 
	    std::string denominatorTau_decayMode = JetMCTagUtils::genTauDecayMode(genTau);
	    isMatched = true;
	  }
	}
	if ( isMatched && (decayMode_ == "all" || denominatorTau_decayMode == decayMode_) )
	{
	  fillHistograms(l1Tau, evtWeight);
	}
      }
    }
    MonitorElement* me_absChargedIso_;
    TH1* histogram_absChargedIso_;
    MonitorElement* me_absNeutralIso_;
    TH1* histogram_absNeutralIso_;
    MonitorElement* me_absCombinedIso_;
    TH1* histogram_absCombinedIso_;
    MonitorElement* me_relChargedIso_;
    TH1* histogram_relChargedIso_;
    MonitorElement* me_relNeutralIso_;
    TH1* histogram_relNeutralIso_;
    MonitorElement* me_relCombinedIso_;
    TH1* histogram_relCombinedIso_;
    MonitorElement* me_pt_;
    TH1* histogram_pt_;
    double ptThreshold_;
    double max_absEta_;
    std::string decayMode_;
    double dRmatch_; 
  };
  std::vector<isolationPlotEntryType*> isolationPlots_;
};

#endif   
