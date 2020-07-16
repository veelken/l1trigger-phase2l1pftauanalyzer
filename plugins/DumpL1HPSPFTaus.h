#ifndef L1Trigger_TallinnL1PFTauAnalyzer_DumpL1HPSPFTaus_h
#define L1Trigger_TallinnL1PFTauAnalyzer_DumpL1HPSPFTaus_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Phase2L1Taus/interface/L1HPSPFTau.h"    // l1t::L1HPSPFTau
#include "DataFormats/Phase2L1Taus/interface/L1HPSPFTauFwd.h" // l1t::L1HPSPFTauCollection

#include <vector>
#include <string>

class DumpL1HPSPFTaus : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpL1HPSPFTaus(const edm::ParameterSet&);
    
  // destructor
  ~DumpL1HPSPFTaus();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<l1t::L1HPSPFTauCollection> token_;

  bool debug_;
};

#endif   
