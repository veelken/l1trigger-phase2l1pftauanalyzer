#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/L1HPSPFTauQualityAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

L1HPSPFTauQualityAnalyzer::L1HPSPFTauQualityAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<l1t::L1HPSPFTauCollection>(src_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

L1HPSPFTauQualityAnalyzer::~L1HPSPFTauQualityAnalyzer()
{
  for ( auto plot : plots_ ) 
  {
    delete plot;
  }
}

void L1HPSPFTauQualityAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("L1HPSPFTauQualityAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());
  
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
	dqmDirectory.Append(Form("/pt%1.0fto%1.0f", min_pt, max_pt));
      }
      else if ( min_pt >= 0. ) 
      {
	dqmDirectory.Append(Form("/ptGt%1.0f", min_pt));
      }
      else if ( max_pt > 0. ) 
      {
	dqmDirectory.Append(Form("/ptLt%1.0f", max_pt));
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
      plotEntryType* plots_vLoose  = new plotEntryType(min_pt, max_pt, min_absEta, max_absEta, 0.40, -1.); // vLoose
      plots_vLoose->bookHistograms(dqmStore);
      plots_.push_back(plots_vLoose);
      plotEntryType* plots_Loose   = new plotEntryType(min_pt, max_pt, min_absEta, max_absEta, 0.20, -1.); // Loose
      plots_Loose->bookHistograms(dqmStore);
      plots_.push_back(plots_Loose);
      plotEntryType* plots_Medium  = new plotEntryType(min_pt, max_pt, min_absEta, max_absEta, 0.10, -1.); // Medium
      plots_Medium->bookHistograms(dqmStore);
      plots_.push_back(plots_Medium);
      plotEntryType* plots_Tight   = new plotEntryType(min_pt, max_pt, min_absEta, max_absEta, 0.05, -1.); // Tight
      plots_Tight->bookHistograms(dqmStore);
      plots_.push_back(plots_Tight);
      plotEntryType* plots_vTight  = new plotEntryType(min_pt, max_pt, min_absEta, max_absEta, 0.02, -1.); // vTight
      plots_vTight->bookHistograms(dqmStore);
      plots_.push_back(plots_vTight);
      plotEntryType* plots_vvTight = new plotEntryType(min_pt, max_pt, min_absEta, max_absEta, 0.01, -1.); // vvTight
      plots_vvTight->bookHistograms(dqmStore);
      plots_.push_back(plots_vvTight);
    }
  }
}

void L1HPSPFTauQualityAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<l1t::L1HPSPFTauCollection> l1PFTaus;
  evt.getByToken(token_, l1PFTaus);
  
  const double evtWeight = 1.;

  for ( auto plot : plots_ ) 
  {    
    plot->fillHistograms(*l1PFTaus, evtWeight);
  }
}

void L1HPSPFTauQualityAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(L1HPSPFTauQualityAnalyzer);
