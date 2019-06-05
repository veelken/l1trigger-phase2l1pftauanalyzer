#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/TallinnL1PFTauSeedAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "L1Trigger/TallinnL1PFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverFlow()

#include <TMath.h> // TMath::Abs(), TMath::Pi()

#include <iostream>
#include <iomanip>

TallinnL1PFTauSeedAnalyzer::TallinnL1PFTauSeedAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , min_seedChargedPFCand_pt_(cfg.getParameter<double>("min_seedChargedPFCand_pt"))
  , max_seedChargedPFCand_eta_(cfg.getParameter<double>("max_seedChargedPFCand_eta"))
  , max_seedChargedPFCand_dz_(cfg.getParameter<double>("max_seedChargedPFCand_dz"))
  , min_seedPFJet_pt_(cfg.getParameter<double>("min_seedPFJet_pt"))
  , max_seedPFJet_eta_(cfg.getParameter<double>("max_seedPFJet_eta"))
  , me_seedChargedPFCand_pt_(nullptr)
  , histogram_seedChargedPFCand_pt_(nullptr)
  , me_seedChargedPFCand_eta_(nullptr)
  , histogram_seedChargedPFCand_eta_(nullptr)
  , me_seedChargedPFCand_phi_(nullptr)
  , histogram_seedChargedPFCand_phi_(nullptr)
  , me_seedChargedPFCand_dz_(nullptr)
  , histogram_seedChargedPFCand_dz_(nullptr)
  , me_numSeedChargedPFCands_(nullptr)
  , histogram_numSeedChargedPFCands_(nullptr)
{
  srcL1PFCands_ = cfg.getParameter<edm::InputTag>("srcL1PFCands");
  tokenL1PFCands_ = consumes<l1t::PFCandidateCollection>(srcL1PFCands_);
  srcL1PFJets_ = cfg.getParameter<edm::InputTag>("srcL1PFJets");
  tokenL1PFJets_ = consumes<l1t::PFJetCollection>(srcL1PFJets_);
  srcL1Vertices_ = cfg.getParameter<edm::InputTag>("srcL1Vertices");
  if ( srcL1Vertices_.label() != "" ) 
  {
    tokenL1Vertices_ = consumes<l1t::VertexCollection>(srcL1Vertices_);
  }

  edm::ParameterSet cfg_signalQualityCuts = cfg.getParameter<edm::ParameterSet>("signalQualityCuts");
  signalQualityCuts_dzCut_disabled_ = readL1PFTauQualityCuts(cfg_signalQualityCuts, "disabled");

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

TallinnL1PFTauSeedAnalyzer::~TallinnL1PFTauSeedAnalyzer()
{}

void TallinnL1PFTauSeedAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("TallinnL1PFTauSeedAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  const int numBins_pt = 24;
  float binning_pt[numBins_pt + 1] = { 
    0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12.5, 15., 17.5, 20., 25., 30., 35., 40., 50., 60., 80., 100., 150., 250.
  };

  me_seedChargedPFCand_pt_ = dqmStore.book1D("seedChargedPFCand_pt", "seedChargedPFCand_pt", numBins_pt, binning_pt);
  histogram_seedChargedPFCand_pt_ = me_seedChargedPFCand_pt_->getTH1();
  assert(histogram_seedChargedPFCand_pt_);
  me_seedChargedPFCand_eta_ = dqmStore.book1D("seedChargedPFCand_eta", "seedChargedPFCand_eta", 30, -3.0, +3.0);
  histogram_seedChargedPFCand_eta_ = me_seedChargedPFCand_eta_->getTH1();
  assert(histogram_seedChargedPFCand_eta_);
  me_seedChargedPFCand_phi_ = dqmStore.book1D("seedChargedPFCand_phi", "seedChargedPFCand_phi", 18, -TMath::Pi(), +TMath::Pi());
  histogram_seedChargedPFCand_phi_ = me_seedChargedPFCand_phi_->getTH1();
  assert(histogram_seedChargedPFCand_phi_);
  me_seedChargedPFCand_dz_ = dqmStore.book1D("seedChargedPFCand_dz", "seedChargedPFCand_dz", 1000, -50., +50.);
  histogram_seedChargedPFCand_dz_ = me_seedChargedPFCand_dz_->getTH1();
  assert(histogram_seedChargedPFCand_dz_);
  me_numSeedChargedPFCands_ = dqmStore.book1D("numSeedChargedPFCands", "numSeedChargedPFCands", 100, -0.5, 99.5);
  histogram_numSeedChargedPFCands_ = me_numSeedChargedPFCands_->getTH1();
  assert(histogram_numSeedChargedPFCands_);

  me_seedPFJet_pt_ = dqmStore.book1D("seedPFJet_pt", "seedPFJet_pt", numBins_pt, binning_pt);
  histogram_seedPFJet_pt_ = me_seedPFJet_pt_->getTH1();
  assert(histogram_seedPFJet_pt_);
  me_seedPFJet_eta_ = dqmStore.book1D("seedPFJet_eta", "seedPFJet_eta", 30, -3.0, +3.0);
  histogram_seedPFJet_eta_ = me_seedPFJet_eta_->getTH1();
  assert(histogram_seedPFJet_eta_);
  me_seedPFJet_phi_ = dqmStore.book1D("seedPFJet_phi", "seedPFJet_phi", 18, -TMath::Pi(), +TMath::Pi());
  histogram_seedPFJet_phi_ = me_seedPFJet_phi_->getTH1();
  me_numSeedPFJets_ = dqmStore.book1D("numSeedPFJets", "numSeedPFJets", 100, -0.5, 99.5);
  histogram_numSeedPFJets_ = me_numSeedPFJets_->getTH1();
  assert(histogram_numSeedPFJets_);
}

void TallinnL1PFTauSeedAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<l1t::PFCandidateCollection> l1PFCands;
  evt.getByToken(tokenL1PFCands_, l1PFCands);

  l1t::VertexRef primaryVertex;
  float primaryVertex_z = 0.;
  if ( srcL1Vertices_.label() != "" ) 
  {
    edm::Handle<l1t::VertexCollection> vertices;
    evt.getByToken(tokenL1Vertices_, vertices);
    if ( vertices->size() > 0 ) 
    {
      primaryVertex = l1t::VertexRef(vertices, 0);
      primaryVertex_z = primaryVertex->z0();
    }
  }

  edm::Handle<l1t::PFJetCollection> l1PFJets;
  evt.getByToken(tokenL1PFJets_, l1PFJets);

  const double evtWeight = 1.;

  int numSeedChargedPFCands = 0;
  for ( l1t::PFCandidateCollection::const_iterator l1PFCand = l1PFCands->begin(); l1PFCand != l1PFCands->end(); ++l1PFCand )
  {
    if ( isSelected(signalQualityCuts_dzCut_disabled_, *l1PFCand, primaryVertex_z) )
    {
      if ( l1PFCand->charge() != 0 && l1PFCand->pt() > min_seedChargedPFCand_pt_ && TMath::Abs(l1PFCand->eta()) < max_seedChargedPFCand_eta_ )
      {
	fillWithOverFlow(histogram_seedChargedPFCand_pt_,  l1PFCand->pt(),  evtWeight);
	fillWithOverFlow(histogram_seedChargedPFCand_eta_, l1PFCand->eta(), evtWeight);
	fillWithOverFlow(histogram_seedChargedPFCand_phi_, l1PFCand->phi(), evtWeight);
        if ( primaryVertex.get() ) 
        {
          l1t::PFTrackRef l1PFTrack = l1PFCand->pfTrack();
          double dz = l1PFTrack->vertex().z() - primaryVertex_z;
	  fillWithOverFlow(histogram_seedChargedPFCand_dz_, dz, evtWeight);
	  if ( TMath::Abs(dz) < max_seedChargedPFCand_dz_ ) 
          {
  	    ++numSeedChargedPFCands;
          }
	}
      }
    }
  }
  fillWithOverFlow(histogram_numSeedChargedPFCands_, numSeedChargedPFCands, evtWeight);

  int numSeedPFJets = 0;
  for ( l1t::PFJetCollection::const_iterator l1PFJet = l1PFJets->begin(); l1PFJet != l1PFJets->end(); ++l1PFJet )
  {
    if ( l1PFJet->pt() > min_seedPFJet_pt_ && std::fabs(l1PFJet->eta()) < max_seedPFJet_eta_ )
    {
      fillWithOverFlow(histogram_seedPFJet_pt_,  l1PFJet->pt(),  evtWeight);
      fillWithOverFlow(histogram_seedPFJet_eta_, l1PFJet->eta(), evtWeight);
      fillWithOverFlow(histogram_seedPFJet_phi_, l1PFJet->phi(), evtWeight);
      ++numSeedPFJets;
    }
  }
  fillWithOverFlow(histogram_numSeedPFJets_, numSeedPFJets, evtWeight);
}

void TallinnL1PFTauSeedAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TallinnL1PFTauSeedAnalyzer);
