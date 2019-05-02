#ifndef L1Trigger_TallinnL1PFTauAnalyzer_DumpTallinL1PFTaus_h
#define L1Trigger_TallinnL1PFTauAnalyzer_DumpTallinL1PFTaus_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTau.h"         // l1t::TallinnL1PFTau
#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTauFwd.h"      // l1t::TallinnL1PFTauCollection

#include <vector>
#include <string>

class DumpTallinL1PFTaus : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpTallinL1PFTaus(const edm::ParameterSet&);
    
  // destructor
  ~DumpTallinL1PFTaus();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<l1t::TallinnL1PFTauCollection> token_;
};

#endif   
