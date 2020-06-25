#ifndef L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauSelector_h
#define L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauSelector_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTau.h"    // l1t::TallinnL1PFTau
#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTauFwd.h" // l1t::TallinnL1PFTauCollection

#include <string>
#include <vector>

class TallinnL1PFTauSelector : public edm::EDProducer 
{
 public:
  explicit TallinnL1PFTauSelector(const edm::ParameterSet& cfg);
  ~TallinnL1PFTauSelector();

 private:
  void produce(edm::Event& evt, const edm::EventSetup& es);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<l1t::TallinnL1PFTauCollection> token_;

  double min_pt_;
  double max_pt_;
  double min_absEta_;
  double max_absEta_;
  double min_leadTrack_pt_;
  double max_leadTrack_pt_;
  double max_relChargedIso_;
  double max_absChargedIso_;
};

#endif
