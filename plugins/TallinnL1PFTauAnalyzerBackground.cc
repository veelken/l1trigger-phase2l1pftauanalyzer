#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/TallinnL1PFTauAnalyzerBackground.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

TallinnL1PFTauAnalyzerBackground::TallinnL1PFTauAnalyzerBackground(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  srcL1PFTaus_ = cfg.getParameter<edm::InputTag>("srcL1PFTaus");
  tokenL1PFTaus_ = consumes<l1t::TallinnL1PFTauCollection>(srcL1PFTaus_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

TallinnL1PFTauAnalyzerBackground::~TallinnL1PFTauAnalyzerBackground()
{
  for ( auto ratePlot : ratePlots_ ) 
  {
    delete ratePlot;
  }
}

void TallinnL1PFTauAnalyzerBackground::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("TallinnL1PFTauAnalyzerBackground") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   0.40, -1., 0.4)); // vLoose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   0.20, -1., 0.4)); // Loose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   0.10, -1., 0.4)); // Medium
  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   0.05, -1., 0.4)); // Tight
  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   0.02, -1., 0.4)); // vTight
  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   0.01, -1., 0.4)); // vvTight
  
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, 0.40, -1., 0.4)); // vLoose
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, 0.20, -1., 0.4)); // Loose
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, 0.10, -1., 0.4)); // Medium
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, 0.05, -1., 0.4)); // Tight
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, 0.02, -1., 0.4)); // vTight
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, 0.01, -1., 0.4)); // vvTight

  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   0.40, -1., 0.4)); // vLoose
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   0.20, -1., 0.4)); // Loose
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   0.10, -1., 0.4)); // Medium
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   0.05, -1., 0.4)); // Tight
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   0.02, -1., 0.4)); // vTight
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   0.01, -1., 0.4)); // vvTight

  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, 0.40, -1., 0.4)); // vLoose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, 0.20, -1., 0.4)); // Loose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, 0.10, -1., 0.4)); // Medium
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, 0.05, -1., 0.4)); // Tight
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, 0.02, -1., 0.4)); // vTight
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, 0.01, -1., 0.4)); // vvTight

  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   0.40, -1., 0.4)); // vLoose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   0.20, -1., 0.4)); // Loose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   0.10, -1., 0.4)); // Medium
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   0.05, -1., 0.4)); // Tight
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   0.02, -1., 0.4)); // vTight
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   0.01, -1., 0.4)); // vvTight

  for ( auto ratePlot : ratePlots_ ) 
  {
    ratePlot->bookHistograms(dqmStore);
  }
}

void TallinnL1PFTauAnalyzerBackground::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<l1t::TallinnL1PFTauCollection> l1PFTaus;
  evt.getByToken(tokenL1PFTaus_, l1PFTaus);
  
  const double evtWeight = 1.;

  for ( auto ratePlot : ratePlots_ ) 
  {
    ratePlot->fillHistograms(*l1PFTaus, evtWeight);
  }
}

void TallinnL1PFTauAnalyzerBackground::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TallinnL1PFTauAnalyzerBackground);
