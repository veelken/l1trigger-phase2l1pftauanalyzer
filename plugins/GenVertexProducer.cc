#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/GenVertexProducer.h"

#include "FWCore/Utilities/interface/InputTag.h"

GenVertexProducer::GenVertexProducer(const edm::ParameterSet& cfg) 
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::GenParticleCollection>(src_);

  if ( cfg.exists("pdgIds") ) 
  {
    pdgIds_ = cfg.getParameter<vint>("pdgIds");
  }

  produces<float>("z0");
}

GenVertexProducer::~GenVertexProducer()
{}

void GenVertexProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  std::unique_ptr<float> genVertex_z0(new float());

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByLabel(src_, genParticles);

  double max_genParticle_pt = -1.;
  for ( auto genParticle : *genParticles ) 
  {
    if ( pdgIds_.size() > 0 ) 
    {
      bool isSelected = false;
      for ( auto pdgId : pdgIds_ ) 
      {
	if ( genParticle.pdgId() == pdgId ) isSelected = true;
      }
      if ( !isSelected ) continue;
    }

    if ( genParticle.pt() > max_genParticle_pt ) 
    {
      *genVertex_z0 = genParticle.vertex().z();
    }
  }

  evt.put(std::move(genVertex_z0), "z0");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenVertexProducer);
