#ifndef L1Trigger_TallinnL1PFTauAnalyzer_L1TrackAntiOverlapSelector_h
#define L1Trigger_TallinnL1PFTauAnalyzer_L1TrackAntiOverlapSelector_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h" 
#include "DataFormats/Candidate/interface/Candidate.h" 
#include "DataFormats/Phase2L1ParticleFlow/interface/PFTrack.h" // l1t::PFTrack::L1TTTrackType

#include "DataFormats/Math/interface/deltaR.h"

#include <vector>

class L1TrackAntiOverlapSelector : public edm::EDProducer 
{
 public:
  explicit L1TrackAntiOverlapSelector(const edm::ParameterSet& cfg);
  ~L1TrackAntiOverlapSelector();

 private:
  void produce(edm::Event& evt, const edm::EventSetup& es);

  typedef l1t::PFTrack::L1TTTrackType L1Track;
  typedef std::vector<L1Track> L1TrackCollection;

  edm::InputTag src_;
  edm::EDGetTokenT<L1TrackCollection> token_;

  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcNotToBeFiltered_;
  typedef std::vector<edm::EDGetTokenT<reco::CandidateView>> vCandidateViewToken;
  vCandidateViewToken tokensNotToBeFiltered_;

  //when invert is TRUE the selector looks for overlapping objects
  bool invert_;

  double dRmin_;
};

#endif
