#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/DumpGenParticles.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>
#include <algorithm> // std::sort

DumpGenParticles::DumpGenParticles(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::GenParticleCollection>(src_);
}

DumpGenParticles::~DumpGenParticles()
{}

namespace
{
  bool
  isHigherPt(const reco::GenParticle& genParticle1,
	     const reco::GenParticle& genParticle2)
  {
    return genParticle1.pt() > genParticle2.pt();
  }
}

void DumpGenParticles::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpGenParticles::analyze> (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByLabel(src_, genParticles);

  // sort genTau collection by decreasing pT
  reco::GenParticleCollection genParticles_sorted = *genParticles;
  std::sort(genParticles_sorted.begin(), genParticles_sorted.end(), isHigherPt);
  
  size_t numParticles = genParticles_sorted.size();
  for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) 
  {
    const reco::GenParticle& genParticle = genParticles_sorted.at(idxParticle);
    std::cout << "genParticle #" << idxParticle  << ":" 
	      << " pT = " << genParticle.pt() << ", eta = " << genParticle.eta() << ", phi = " << genParticle.phi() 
	      << " pdgId = " << genParticle.pdgId() << std::endl;
    std::cout << " vertex: x = " << genParticle.vertex().x() << ", y = " << genParticle.vertex().y() << ", z = " << genParticle.vertex().z() << std::endl;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpGenParticles);





