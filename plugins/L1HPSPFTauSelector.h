#ifndef L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauSelector_h
#define L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauSelector_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/Phase2L1Taus/interface/L1HPSPFTau.h"    // l1t::L1HPSPFTau
#include "DataFormats/Phase2L1Taus/interface/L1HPSPFTauFwd.h" // l1t::L1HPSPFTauCollection

#include <string>
#include <vector>

class L1HPSPFTauSelector : public edm::EDProducer 
{
 public:
  explicit L1HPSPFTauSelector(const edm::ParameterSet& cfg);
  ~L1HPSPFTauSelector();

 private:
  void produce(edm::Event& evt, const edm::EventSetup& es);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<l1t::L1HPSPFTauCollection> token_;

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
