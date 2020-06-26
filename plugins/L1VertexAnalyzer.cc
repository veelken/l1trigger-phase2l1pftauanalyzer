#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/L1VertexAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

L1VertexAnalyzer::L1VertexAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_genVertex_z_ = cfg.getParameter<edm::InputTag>("srcGenVertex_z");
  token_genVertex_z_ = consumes<float>(src_genVertex_z_);

  src_l1Vertices_ = cfg.getParameter<edm::InputTag>("srcL1Vertices");
  token_l1Vertices_ = consumes<l1t::VertexCollection>(src_l1Vertices_);

  src_l1PFVertex_z_ = cfg.getParameter<edm::InputTag>("srcL1PFVertex_z");
  token_l1PFVertex_z_ = consumes<float>(src_l1PFVertex_z_);

  src_offlineVertices_ = cfg.getParameter<edm::InputTag>("srcOfflineVertices");
  token_offlineVertices_ = consumes<reco::VertexCollection>(src_offlineVertices_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

L1VertexAnalyzer::~L1VertexAnalyzer()
{
  delete l1VertexPlots_;

  delete l1PFVertexPlots_;

  delete offlineVertexPlots_;
}

void L1VertexAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("JetToTauFakeRateAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());
  dqmStore.setCurrentFolder(Form("%s/%s", dqmDirectory_.data(), "L1Vertex"));
  l1VertexPlots_ = new vertexPlotEntryType();
  l1VertexPlots_->bookHistograms(dqmStore); 
  dqmStore.setCurrentFolder(Form("%s/%s", dqmDirectory_.data(), "L1PFVertex"));
  l1PFVertexPlots_ = new vertexPlotEntryType();
  l1PFVertexPlots_->bookHistograms(dqmStore);
  dqmStore.setCurrentFolder(Form("%s/%s", dqmDirectory_.data(), "offlineVertex"));
  offlineVertexPlots_ = new vertexPlotEntryType();
  offlineVertexPlots_->bookHistograms(dqmStore);
}

namespace
{
  std::vector<double> get_recVertex_z(const l1t::VertexCollection& recVertices)
  {
    std::vector<double> recVertices_z;
    for ( auto recVertex : recVertices )
    {
      recVertices_z.push_back(recVertex.z0());
    }
    return recVertices_z;
  }

  std::vector<double> get_recVertex_z(const reco::VertexCollection& recVertices)
  {
    std::vector<double> recVertices_z;
    for ( auto recVertex : recVertices )
    {
      recVertices_z.push_back(recVertex.position().z());
    }
    return recVertices_z;
  }
}

void L1VertexAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<float> genVertex_z;
  evt.getByToken(token_genVertex_z_, genVertex_z);

  const double evtWeight = 1.;

  edm::Handle<l1t::VertexCollection> l1Vertices;
  evt.getByToken(token_l1Vertices_, l1Vertices);
  std::vector<double> l1Vertices_z = get_recVertex_z(*l1Vertices);
  l1VertexPlots_->fillHistograms(*genVertex_z, l1Vertices_z, evtWeight);

  edm::Handle<float> l1PFVertex_z;
  evt.getByToken(token_l1PFVertex_z_, l1PFVertex_z);
  std::vector<double> l1PFVertices_z = { *l1PFVertex_z };
  l1PFVertexPlots_->fillHistograms(*genVertex_z, l1PFVertices_z, evtWeight);

  edm::Handle<reco::VertexCollection> offlineVertices;
  evt.getByToken(token_offlineVertices_, offlineVertices);
  std::vector<double> offlineVertices_z = get_recVertex_z(*offlineVertices);
  offlineVertexPlots_->fillHistograms(*genVertex_z, offlineVertices_z, evtWeight);
}

void L1VertexAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(L1VertexAnalyzer);
