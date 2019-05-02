#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/DumpL1Tracks.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpL1Tracks::DumpL1Tracks(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<std::vector<l1t::PFTrack::L1TTTrackType>>(src_);
}

DumpL1Tracks::~DumpL1Tracks()
{}

void DumpL1Tracks::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpL1Tracks::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<std::vector<l1t::PFTrack::L1TTTrackType>> l1Tracks;
  evt.getByToken(token_, l1Tracks);
  
  size_t numTracks = l1Tracks->size();
  for ( size_t idxTrack = 0; idxTrack < numTracks; ++idxTrack ) 
  {
    const l1t::PFTrack::L1TTTrackType& l1Track = l1Tracks->at(idxTrack);
    const unsigned nParam = 4;
    double l1Track_pt  = l1Track.getMomentum(nParam).perp();
    double l1Track_eta = l1Track.getMomentum(nParam).eta();
    double l1Track_phi = l1Track.getMomentum(nParam).phi();
    std::cout << "L1Track #" << idxTrack << ":" 
	      << " pT = " << l1Track_pt << ", eta = " << l1Track_eta << ", phi = " << l1Track_phi << std::endl;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpL1Tracks);





