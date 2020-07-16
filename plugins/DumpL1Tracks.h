#ifndef L1Trigger_TallinnL1PFTauAnalyzer_DumpL1Tracks_h
#define L1Trigger_TallinnL1PFTauAnalyzer_DumpL1Tracks_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/L1TParticleFlow/interface/PFTrack.h" // l1t::PFTrack::L1TTTrackType

#include <vector>
#include <string>

class DumpL1Tracks : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpL1Tracks(const edm::ParameterSet&);
    
  // destructor
  ~DumpL1Tracks();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<std::vector<l1t::PFTrack::L1TTTrackType>> token_;
};

#endif   
