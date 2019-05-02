#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/DumpL1Vertices.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpL1Vertices::DumpL1Vertices(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<l1t::VertexCollection>(src_);
}

DumpL1Vertices::~DumpL1Vertices()
{}

namespace
{
  double compSumTrackPt2(const std::vector<edm::Ptr<l1t::Vertex::Track_t>>& l1Tracks)
  {
    double sumTrackPt2 = 0.;    
    for ( size_t idxTrack = 0; idxTrack < l1Tracks.size(); ++idxTrack ) 
    {
      const l1t::Vertex::Track_t& l1Track = *l1Tracks.at(idxTrack);      
      const unsigned nParam = 4;
      double l1Track_pt = l1Track.getMomentum(nParam).perp();
      sumTrackPt2 += (l1Track_pt*l1Track_pt);
    }
    return sumTrackPt2;
  }
}

void DumpL1Vertices::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpL1Vertices::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<l1t::VertexCollection> l1Vertices;
  evt.getByToken(token_, l1Vertices);
  
  size_t numVertices = l1Vertices->size();
  for ( size_t idxVertex = 0; idxVertex < numVertices; ++idxVertex ) 
  {
    const l1t::Vertex& l1Vertex = l1Vertices->at(idxVertex);
    std::cout << "L1Vertex #" << idxVertex << ":" 
	      << " z0 = " << l1Vertex.z0() << " (sum track pT^2 = " << compSumTrackPt2(l1Vertex.tracks()) << ")" << std::endl;
    const std::vector<edm::Ptr<l1t::Vertex::Track_t>>& l1Tracks = l1Vertex.tracks();
    size_t numTracks = l1Tracks.size();
    for ( size_t idxTrack = 0; idxTrack < numTracks; ++idxTrack ) 
    {
      const l1t::Vertex::Track_t& l1Track = *l1Tracks.at(idxTrack);      
      const unsigned nParam = 4;
      double l1Track_pt  = l1Track.getMomentum(nParam).perp();
      double l1Track_eta = l1Track.getMomentum(nParam).eta();
      double l1Track_phi = l1Track.getMomentum(nParam).phi();
      std::cout << " L1Track #" << idxTrack << ":" 
		<< " pT = " << l1Track_pt << ", eta = " << l1Track_eta << ", phi = " << l1Track_phi << std::endl;
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpL1Vertices);





