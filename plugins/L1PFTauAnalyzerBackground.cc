#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/L1PFTauAnalyzerBackground.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

L1PFTauAnalyzerBackground::L1PFTauAnalyzerBackground(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  srcL1PFTaus_ = cfg.getParameter<edm::InputTag>("srcL1PFTaus");
  tokenL1PFTaus_ = consumes<l1t::L1PFTauCollection>(srcL1PFTaus_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

L1PFTauAnalyzerBackground::~L1PFTauAnalyzerBackground()
{
  for ( auto ratePlot : ratePlots_ ) 
  {
    delete ratePlot;
  }
}

void L1PFTauAnalyzerBackground::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("L1PFTauAnalyzerBackground") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   "vLoose")); // vLoose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   "Loose"));  // Loose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   "Medium")); // Medium
  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   "Tight"));  // Tight
  
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, "vLoose")); // vLoose
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, "Loose"));  // Loose
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, "Medium")); // Medium
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, "Tight"));  // Tight

  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   "vLoose")); // vLoose
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   "Loose"));  // Loose
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   "Medium")); // Medium
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   "Tight"));  // Tight

  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, "vLoose")); // vLoose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, "Loose"));  // Loose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, "Medium")); // Medium
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, "Tight"));  // Tight

  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   "vLoose")); // vLoose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   "Loose"));  // Loose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   "Medium")); // Medium
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   "Tight"));  // Tight

  for ( auto ratePlot : ratePlots_ ) 
  {
    ratePlot->bookHistograms(dqmStore);
  }
}

void L1PFTauAnalyzerBackground::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<l1t::L1PFTauCollection> l1PFTaus;
  evt.getByToken(tokenL1PFTaus_, l1PFTaus);
  
  const double evtWeight = 1.;

  for ( auto ratePlot : ratePlots_ ) 
  {
    ratePlot->fillHistograms(*l1PFTaus, evtWeight);
  }
}

void L1PFTauAnalyzerBackground::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(L1PFTauAnalyzerBackground);
