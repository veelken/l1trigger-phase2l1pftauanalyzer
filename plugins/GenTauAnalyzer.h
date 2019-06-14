#ifndef L1Trigger_TallinnL1PFTauAnalyzer_GenTauAnalyzer_h
#define L1Trigger_TallinnL1PFTauAnalyzer_GenTauAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/JetReco/interface/GenJet.h"           // reco::GenJet
#include "DataFormats/JetReco/interface/GenJetCollection.h" // reco::GenJetCollection
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"     // JetMCTagUtils::genTauDecayMode()

#include <TH1.h>     // TH1
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs(), TMath::Pi()

#include <vector>    // std::vector
#include <string>    // std::string

class GenTauAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit GenTauAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~GenTauAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<reco::GenJetCollection> token_;

  double min_pt_;
  double max_pt_;
  double min_absEta_;
  double max_absEta_;

  std::string dqmDirectory_;

  struct plotEntryType
  {
    plotEntryType(double min_pt, double max_pt, double min_absEta, double max_absEta, const std::string& decayMode, int idx)
      : me_pt_(nullptr)
      , histogram_pt_(nullptr)
      , me_eta_(nullptr)
      , histogram_eta_(nullptr)
      , me_phi_(nullptr)
      , histogram_phi_(nullptr)
      , min_pt_(min_pt)
      , max_pt_(max_pt)
      , min_absEta_(min_absEta)
      , max_absEta_(max_absEta)
      , decayMode_(decayMode)
      , idx_(idx)
    {}
    ~plotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_prefix;
      if      ( idx_ == 0 ) histogramName_prefix = "leadingTau";
      else if ( idx_ == 1 ) histogramName_prefix = "subleadingTau";
      else throw cms::Exception("GenTauAnalyzer") 
	<< " Invalid Configuration parameter 'idx' = " << idx_ << " !!\n";;
      TString histogramName_suffix = decayMode_.data();      
      if      ( min_absEta_ >= 0. && max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEta%1.2fto%1.2f", min_absEta_, max_absEta_));
      else if ( min_absEta_ >= 0.                     ) histogramName_suffix.Append(Form("_absEtaGt%1.2f", min_absEta_));
      else if (                      max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f", max_absEta_));
      if      ( min_pt_     >= 0. && max_pt_     > 0. ) histogramName_suffix.Append(Form("_pt%1.2fto%1.2f", min_pt_, max_pt_));
      else if ( min_pt_     >= 0.                     ) histogramName_suffix.Append(Form("_ptGt%1.2f", min_pt_));
      else if (                      max_pt_     > 0. ) histogramName_suffix.Append(Form("_ptLt%1.2f", max_pt_));
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_pt = Form("%s_pt_%s", histogramName_prefix.Data(), histogramName_suffix.Data());
      me_pt_ = dqmStore.book1D(histogramName_pt.Data(), histogramName_pt.Data(), 40, 0., 100.);
      histogram_pt_ = me_pt_->getTH1();
      assert(histogram_pt_);

      TString histogramName_eta = Form("%s_eta_%s", histogramName_prefix.Data(), histogramName_suffix.Data());
      me_eta_ = dqmStore.book1D(histogramName_eta.Data(), histogramName_eta.Data(), 30, -3., +3.);
      histogram_eta_ = me_eta_->getTH1();
      assert(histogram_eta_);
      
      TString histogramName_phi = Form("%s_phi_%s", histogramName_prefix.Data(), histogramName_suffix.Data());
      me_phi_ = dqmStore.book1D(histogramName_phi.Data(), histogramName_phi.Data(), 18, -TMath::Pi(), +TMath::Pi());
      histogram_phi_ = me_phi_->getTH1();
      assert(histogram_phi_);
    }
    void fillHistograms(const reco::GenJetCollection& genTaus, double evtWeight)
    {
      if ( (int)genTaus.size() > idx_ ) {
	const reco::GenJet& genTau = genTaus.at(idx_);

	std::string genTau_decayMode = JetMCTagUtils::genTauDecayMode(genTau);
	if ( decayMode_ != "all" && genTau_decayMode != decayMode_ ) return;

	double genTau_absEta = TMath::Abs(genTau.eta());
	if ( genTau_absEta > min_absEta_ && genTau_absEta < max_absEta_ ) 
	{
	  histogram_pt_->Fill(genTau.pt(), evtWeight);
	}
	if ( genTau.pt() > min_pt_ && genTau.pt() < max_pt_ )
	{
	  histogram_eta_->Fill(genTau.eta(), evtWeight);
	}
	if ( genTau.pt() > min_pt_ && genTau.pt() < max_pt_ && genTau_absEta > min_absEta_ && genTau_absEta < max_absEta_ ) 
	{
	  histogram_phi_->Fill(genTau.phi(), evtWeight);
	}
      }
    }
    MonitorElement* me_pt_;
    TH1* histogram_pt_;
    MonitorElement* me_eta_;
    TH1* histogram_eta_;
    MonitorElement* me_phi_;
    TH1* histogram_phi_;
    double min_pt_;
    double max_pt_;
    double min_absEta_;
    double max_absEta_;
    std::string decayMode_;
    int idx_;
  };
  std::vector<plotEntryType*> plots_;
};

#endif   
