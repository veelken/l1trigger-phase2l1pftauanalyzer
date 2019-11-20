#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/TallinnL1PFTauQualityAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

TallinnL1PFTauQualityAnalyzer::TallinnL1PFTauQualityAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<l1t::TallinnL1PFTauCollection>(src_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

TallinnL1PFTauQualityAnalyzer::~TallinnL1PFTauQualityAnalyzer()
{
  for ( auto plot : plots_ ) 
  {
    delete plot;
  }
}

void TallinnL1PFTauQualityAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("TallinnL1PFTauQualityAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  
  std::vector<double> min_ptValues = { 20., 25., 30., 35., 40., 50., 60.,  80., 20. };
  std::vector<double> max_ptValues = { 25., 30., 35., 40., 50., 60., 80., 100., -1. };
  assert(min_ptValues.size() == max_ptValues.size());
  std::vector<double> min_absEtaValues = { -1.,   1.4,   1.4, -1.,    -1.  };
  std::vector<double> max_absEtaValues = {  1.4,  2.172, 2.4,  2.172,  2.4 };
  assert(min_absEtaValues.size() == max_absEtaValues.size());
  size_t numPtRanges = min_ptValues.size();
  for ( size_t idxPtRange = 0; idxPtRange < numPtRanges; ++idxPtRange )
  {
    double min_pt = min_ptValues[idxPtRange];
    double max_pt = max_ptValues[idxPtRange];
    size_t numAbsEtaRanges = min_absEtaValues.size();
    for ( size_t idxAbsEtaRange = 0; idxAbsEtaRange < numAbsEtaRanges; ++idxAbsEtaRange )
    {
      double min_absEta = min_absEtaValues[idxAbsEtaRange];
      double max_absEta = max_absEtaValues[idxAbsEtaRange];
    
      TString dqmDirectory = dqmDirectory_.data();
      if ( min_pt >= 0. && max_pt > 0. ) 
      { 
	dqmDirectory.Append(Form("/pt%1.2fto%1.2f", min_pt, max_pt));
      }
      else if ( min_pt >= 0. ) 
      {
	dqmDirectory.Append(Form("/ptGt%1.2f", min_pt));
      }
      else if ( max_pt > 0. ) 
      {
	dqmDirectory.Append(Form("/ptLt%1.2f", max_pt));
      }
      if ( min_absEta >= 0. && max_absEta > 0. ) 
      { 
	dqmDirectory.Append(Form("/absEta%1.2fto%1.2f", min_absEta, max_absEta));
      }
      else if ( min_absEta >= 0. ) 
      {
	dqmDirectory.Append(Form("/absEtaGt%1.2f", min_absEta));
      }
      else if ( max_absEta > 0. ) 
      {
	dqmDirectory.Append(Form("/absEtaLt%1.2f", max_absEta));
      }
      dqmDirectory = dqmDirectory.ReplaceAll(".", "p");

      dqmStore.setCurrentFolder(dqmDirectory.Data());
      plotEntryType* plot = new plotEntryType(min_pt, max_pt, min_absEta, max_absEta);
      plot->bookHistograms(dqmStore);
      plots_.push_back(plot);
    }
  }
}

void TallinnL1PFTauQualityAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<l1t::TallinnL1PFTauCollection> l1PFTaus;
  evt.getByToken(token_, l1PFTaus);
  
  const double evtWeight = 1.;

  for ( auto plot : plots_ ) 
  {    
    plot->fillHistograms(*l1PFTaus, evtWeight);
  }
}

void TallinnL1PFTauQualityAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TallinnL1PFTauQualityAnalyzer);
