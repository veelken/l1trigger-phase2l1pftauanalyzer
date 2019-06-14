#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/GenTauAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

GenTauAnalyzer::GenTauAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::GenJetCollection>(src_);

  min_pt_ = cfg.getParameter<double>("min_pt");
  max_pt_ = cfg.getParameter<double>("max_pt");
  min_absEta_ = cfg.getParameter<double>("min_absEta");
  max_absEta_ = cfg.getParameter<double>("max_absEta");

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

GenTauAnalyzer::~GenTauAnalyzer()
{
  for ( auto plot : plots_ ) 
  {
    delete plot;
  }
}

void GenTauAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("GenTauAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  std::vector<std::string> decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  for ( auto decayMode : decayModes )
  {    
    plots_.push_back(new plotEntryType(min_pt_, max_pt_, min_absEta_, max_absEta_, decayMode, 0)); // leadingTau
    plots_.push_back(new plotEntryType(min_pt_, max_pt_, min_absEta_, max_absEta_, decayMode, 1)); // subleadingTau
  }

  for ( auto plot : plots_ ) 
  {
    plot->bookHistograms(dqmStore);
  }
}

namespace
{
  bool
  isHigherPt(const reco::GenJet& genTau1,
	     const reco::GenJet& genTau2)
  {
    return genTau1.pt() > genTau2.pt();
  }
}

void GenTauAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::GenJetCollection> genTaus;
  evt.getByToken(token_, genTaus);
  if ( !(genTaus->size() >= 2) ) return;
  
  // sort genTau collection by decreasing pT
  reco::GenJetCollection genTaus_sorted = *genTaus;
  std::sort(genTaus_sorted.begin(), genTaus_sorted.end(), isHigherPt);

  const double evtWeight = 1.;

  for ( auto plot : plots_ ) 
  {    
    plot->fillHistograms(genTaus_sorted, evtWeight);
  }
}

void GenTauAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GenTauAnalyzer);
