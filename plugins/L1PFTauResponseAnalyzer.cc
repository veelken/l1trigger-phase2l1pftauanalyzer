#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/L1PFTauResponseAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

L1PFTauResponseAnalyzer::L1PFTauResponseAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  srcL1PFTaus_ = cfg.getParameter<edm::InputTag>("srcL1PFTaus");
  tokenL1PFTaus_ = consumes<l1t::L1PFTauCollection>(srcL1PFTaus_);
  srcRefTaus_ = cfg.getParameter<edm::InputTag>("srcRefTaus");
  std::string typeRefTaus_string = cfg.getParameter<std::string>("typeRefTaus");
  if ( typeRefTaus_string == "gen"     ) 
  {
    typeRefTaus_ = kGen;
    tokenRefTaus_gen_ = consumes<reco::GenJetCollection>(srcRefTaus_);
  }
  else if ( typeRefTaus_string == "offline" ) 
  {
    typeRefTaus_ = kOffline;
    tokenRefTaus_offline_ = consumes<pat::TauCollection>(srcRefTaus_);
  }
  else
  {
    throw cms::Exception("L1PFTauResponseAnalyzer") 
      << " Invalid Configuration parameter 'typeRefTaus' = " << typeRefTaus_string << " !!\n";;
  }

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

L1PFTauResponseAnalyzer::~L1PFTauResponseAnalyzer()
{
  for ( auto responsePlot : responsePlots_ ) 
  {
    delete responsePlot;
  }
}

void L1PFTauResponseAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("L1PFTauResponseAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  std::vector<std::string> decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  std::vector<double> min_pts = { 20., 25., 30., 35., 40. };
  for ( auto decayMode : decayModes )
  {
    for ( auto min_pt : min_pts )
    {
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.0, decayMode, "vLoose")); // vLoose
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.0, decayMode, "Loose"));  // Loose
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.0, decayMode, "Medium")); // Medium
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.0, decayMode, "Tight"));  // Tight
    
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.4, decayMode, "vLoose")); // vLoose
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.4, decayMode, "Loose"));  // Loose
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.4, decayMode, "Medium")); // Medium
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.4, decayMode, "Tight"));  // Tight
    }
  }

  for ( auto responsePlot : responsePlots_ ) 
  {
    responsePlot->bookHistograms(dqmStore);
  }
}

void L1PFTauResponseAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<l1t::L1PFTauCollection> l1PFTaus;
  evt.getByToken(tokenL1PFTaus_, l1PFTaus);
  
  const double evtWeight = 1.;

  if ( typeRefTaus_ == kGen )
  {
    edm::Handle<reco::GenJetCollection> refTaus_gen;
    evt.getByToken(tokenRefTaus_gen_, refTaus_gen);

    for ( auto responsePlot : responsePlots_ ) 
    {    
      responsePlot->fillHistograms(*l1PFTaus, *refTaus_gen, evtWeight);
    }
  }

  if ( typeRefTaus_ == kOffline )
  {
    edm::Handle<pat::TauCollection> refTaus_offline;
    evt.getByToken(tokenRefTaus_offline_, refTaus_offline);

    for ( auto responsePlot : responsePlots_ ) 
    {    
      responsePlot->fillHistograms(*l1PFTaus, *refTaus_offline, evtWeight);
    }
  }
}

void L1PFTauResponseAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(L1PFTauResponseAnalyzer);
