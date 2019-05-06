#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/TallinnL1PFTauIsolationAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

TallinnL1PFTauIsolationAnalyzer::TallinnL1PFTauIsolationAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<l1t::TallinnL1PFTauCollection>(src_ );
  srcGenTaus_ = cfg.getParameter<edm::InputTag>("srcGenTaus");
  tokenGenTaus_ = consumes<reco::GenJetCollection>(srcGenTaus_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

TallinnL1PFTauIsolationAnalyzer::~TallinnL1PFTauIsolationAnalyzer()
{
  for ( auto isolationPlot : isolationPlots_ ) 
  {
    delete isolationPlot;
  }
}

void TallinnL1PFTauIsolationAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("TallinnL1PFTauIsolationAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  std::vector<std::string> decayModes;
  if ( srcGenTaus_.label() != "" ) 
  {
    decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  } 
  else
  {
    decayModes = { "all" };
  }
  std::vector<double> ptThresholds = { 20., 30., 40. };
  for ( auto decayMode : decayModes )
  {
    for ( auto ptThreshold : ptThresholds )
    {
      isolationPlots_.push_back(new isolationPlotEntryType(ptThreshold, 1.0, decayMode)); 
      isolationPlots_.push_back(new isolationPlotEntryType(ptThreshold, 1.4, decayMode));
    }
  }

  for ( auto isolationPlot : isolationPlots_ ) 
  {
    isolationPlot->bookHistograms(dqmStore);
  }
}

void TallinnL1PFTauIsolationAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<l1t::TallinnL1PFTauCollection> l1Taus;
  evt.getByToken(token_, l1Taus);
  
  edm::Handle<reco::GenJetCollection> genTaus;
  if ( srcGenTaus_.label() != "" ) 
  {
    evt.getByToken(tokenGenTaus_, genTaus);
  }

  const double evtWeight = 1.;

  for ( auto isolationPlot : isolationPlots_ ) 
  {    
    if ( srcGenTaus_.label() != "" ) 
    {
      isolationPlot->fillHistograms_wGenMatching(*l1Taus, *genTaus, evtWeight);
    }
    else
    {
      isolationPlot->fillHistograms_woGenMatching(*l1Taus, evtWeight);
    }  
  }
}

void TallinnL1PFTauIsolationAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TallinnL1PFTauIsolationAnalyzer);
