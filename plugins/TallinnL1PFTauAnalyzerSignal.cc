#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/TallinnL1PFTauAnalyzerSignal.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

TallinnL1PFTauAnalyzerSignal::TallinnL1PFTauAnalyzerSignal(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  srcNumerator_ = cfg.getParameter<edm::InputTag>("srcNumerator");
  tokenNumerator_ = consumes<l1t::TallinnL1PFTauCollection>(srcNumerator_ );
  srcDenominator_ = cfg.getParameter<edm::InputTag>("srcDenominator");
  tokenDenominator_ = consumes<reco::GenJetCollection>(srcDenominator_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

TallinnL1PFTauAnalyzerSignal::~TallinnL1PFTauAnalyzerSignal()
{
  for ( auto efficiencyPlot : efficiencyPlots_ ) 
  {
    delete efficiencyPlot;
  }
}

void TallinnL1PFTauAnalyzerSignal::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("TallinnL1PFTauAnalyzerSignal") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  std::vector<std::string> decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  std::vector<double> ptThresholds = { 20., 30., 40. };
  for ( auto decayMode : decayModes )
  {
    for ( auto ptThreshold : ptThresholds )
    {
      efficiencyPlots_.push_back(new efficiencyPlotEntryType(45., 1.0, decayMode, ptThreshold, 0.40, -1.)); // vLoose
      efficiencyPlots_.push_back(new efficiencyPlotEntryType(45., 1.0, decayMode, ptThreshold, 0.20, -1.)); // Loose
      efficiencyPlots_.push_back(new efficiencyPlotEntryType(45., 1.0, decayMode, ptThreshold, 0.10, -1.)); // Medium
      efficiencyPlots_.push_back(new efficiencyPlotEntryType(45., 1.0, decayMode, ptThreshold, 0.05, -1.)); // Tight
    
      efficiencyPlots_.push_back(new efficiencyPlotEntryType(45., 1.4, decayMode, ptThreshold, 0.40, -1.)); // vLoose
      efficiencyPlots_.push_back(new efficiencyPlotEntryType(45., 1.4, decayMode, ptThreshold, 0.20, -1.)); // Loose
      efficiencyPlots_.push_back(new efficiencyPlotEntryType(45., 1.4, decayMode, ptThreshold, 0.10, -1.)); // Medium
      efficiencyPlots_.push_back(new efficiencyPlotEntryType(45., 1.4, decayMode, ptThreshold, 0.05, -1.)); // Tight
    }
  }

  for ( auto efficiencyPlot : efficiencyPlots_ ) 
  {
    efficiencyPlot->bookHistograms(dqmStore);
  }
}

void TallinnL1PFTauAnalyzerSignal::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<l1t::TallinnL1PFTauCollection> numeratorTaus;
  evt.getByToken(tokenNumerator_, numeratorTaus);
  
  edm::Handle<reco::GenJetCollection> denominatorTaus;
  evt.getByToken(tokenDenominator_, denominatorTaus);
  
  const double evtWeight = 1.;

  for ( auto efficiencyPlot : efficiencyPlots_ ) 
  {    
    efficiencyPlot->fillHistograms(*numeratorTaus, *denominatorTaus, evtWeight);
  }
}

void TallinnL1PFTauAnalyzerSignal::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TallinnL1PFTauAnalyzerSignal);
