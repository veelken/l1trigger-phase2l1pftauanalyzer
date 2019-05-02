#ifndef L1Trigger_TallinnL1PFTauAnalyzer_DumpL1PFCandidates_h
#define L1Trigger_TallinnL1PFTauAnalyzer_DumpL1PFCandidates_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidateFwd.h"

#include <vector>
#include <string>

class DumpL1PFCandidates : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpL1PFCandidates(const edm::ParameterSet&);
    
  // destructor
  ~DumpL1PFCandidates();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<l1t::PFCandidateCollection> token_;

  double minPt_;
  bool printPuppiWeights_;
};

#endif   
