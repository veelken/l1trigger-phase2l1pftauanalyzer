#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/DumpL1Vertices.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpL1Vertices::DumpL1Vertices(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<l1t::TkPrimaryVertexCollection>(src_);
}

DumpL1Vertices::~DumpL1Vertices()
{}

void DumpL1Vertices::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpL1Vertices::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<l1t::TkPrimaryVertexCollection> l1Vertices;
  evt.getByToken(token_, l1Vertices);
  
  size_t numVertices = l1Vertices->size();
  for ( size_t idxVertex = 0; idxVertex < numVertices; ++idxVertex ) 
  {
    const l1t::TkPrimaryVertex& l1Vertex = l1Vertices->at(idxVertex);
    std::cout << "L1Vertex #" << idxVertex << ":" 
	      << " z0 = " << l1Vertex.zvertex() << " (sum track pT^2 = " << l1Vertex.sum() << ")" << std::endl;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpL1Vertices);





