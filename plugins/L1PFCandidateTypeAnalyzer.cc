#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/L1PFCandidateTypeAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

L1PFCandidateTypeAnalyzer::L1PFCandidateTypeAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , pfChargedHadronPlots_(nullptr)
  , pfChargedHadronPileupPlots_(nullptr)
  , pfElectronPlots_(nullptr)
  , pfNeutralHadronPlots_(nullptr)
  , pfPhotonPlots_(nullptr)
  , pfMuonPlots_(nullptr)
{
  src_l1PFCands_ = cfg.getParameter<edm::InputTag>("srcL1PFCands");
  token_l1PFCands_ = consumes<l1t::PFCandidateCollection>(src_l1PFCands_);

  srcL1Vertices_ = cfg.getParameter<edm::InputTag>("srcL1Vertices");
  if ( srcL1Vertices_.label() != "" ) 
  {
    tokenL1Vertices_ = consumes<l1t::VertexCollection>(srcL1Vertices_);
  }

  srcPileupSummaryInfo_ = cfg.getParameter<edm::InputTag>("srcPileupSummaryInfo");
  if ( srcPileupSummaryInfo_.label() != "" ) 
  {
    tokenPileupSummaryInfo_ = consumes<PileupSummaryInfoCollection>(srcPileupSummaryInfo_);
  }

  edm::ParameterSet cfg_isolationQualityCuts = cfg.getParameter<edm::ParameterSet>("isolationQualityCuts");
  isolationQualityCuts_dzCut_disabled_        = readL1PFTauQualityCuts(cfg_isolationQualityCuts, "disabled",        false);
  isolationQualityCuts_dzCut_enabled_primary_ = readL1PFTauQualityCuts(cfg_isolationQualityCuts, "enabled_primary", false);
  isolationQualityCuts_dzCut_enabled_pileup_  = readL1PFTauQualityCuts(cfg_isolationQualityCuts, "enabled_pileup",  false);

  pfChargedHadronPlots_       = new pfCandTypePlotEntryType(-1., 1.4, "chargedHadron");
  pfChargedHadronPileupPlots_ = new pfCandTypePlotEntryType(-1., 1.4, "chargedHadronPileup");
  pfElectronPlots_            = new pfCandTypePlotEntryType(-1., 1.4, "electron");
  pfNeutralHadronPlots_       = new pfCandTypePlotEntryType(-1., 1.4, "neutralHadron");
  pfPhotonPlots_              = new pfCandTypePlotEntryType(-1., 1.4, "photon");
  pfMuonPlots_                = new pfCandTypePlotEntryType(-1., 1.4, "muon");

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

L1PFCandidateTypeAnalyzer::~L1PFCandidateTypeAnalyzer()
{}

void L1PFCandidateTypeAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("L1PFCandidateTypeAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());

  dqmStore.setCurrentFolder(dqmDirectory_);

  pfChargedHadronPlots_->bookHistograms(dqmStore);
  pfChargedHadronPileupPlots_->bookHistograms(dqmStore);
  pfElectronPlots_->bookHistograms(dqmStore);
  pfNeutralHadronPlots_->bookHistograms(dqmStore);
  pfPhotonPlots_->bookHistograms(dqmStore);
  pfMuonPlots_->bookHistograms(dqmStore);

  me_EventCounter_ = dqmStore.book1D("EventCounter", "EventCounter", 1, -0.5, +0.5);
  histogram_EventCounter_ = me_EventCounter_->getTH1();
  assert(histogram_EventCounter_);
}
    
void L1PFCandidateTypeAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<L1PFCandidateTypeAnalyzer::analyze>:" << std::endl;

  edm::Handle<l1t::PFCandidateCollection> l1PFCands;
  evt.getByToken(token_l1PFCands_, l1PFCands);

  float primaryVertex_z = 0.;
  if ( srcL1Vertices_.label() != "" ) 
  {
    edm::Handle<l1t::VertexCollection> vertices;
    evt.getByToken(tokenL1Vertices_, vertices);
    if ( vertices->size() > 0 ) 
    {
      edm::Ref<l1t::VertexCollection> primaryVertex = l1t::VertexRef(vertices, 0);
      primaryVertex_z = primaryVertex->z0();
    }
  }
  //std::cout << "primaryVertex_z = " << primaryVertex_z << std::endl;

  int numPileup = -1;
  if ( srcPileupSummaryInfo_.label() != "" ) 
  {
    edm::Handle<PileupSummaryInfoCollection> pileupSummaryInfos;
    evt.getByToken(tokenPileupSummaryInfo_, pileupSummaryInfos);
    for ( PileupSummaryInfoCollection::const_iterator pileupSummaryInfo = pileupSummaryInfos->begin(); pileupSummaryInfo != pileupSummaryInfos->end(); ++pileupSummaryInfo ) 
    {
      if ( pileupSummaryInfo->getBunchCrossing() == 0 ) 
      {
	numPileup = pileupSummaryInfo->getPU_NumInteractions();
      }
    }
  }

  const double evtWeight = 1.;
  
  //int idxL1PFCand = 0;
  for ( l1t::PFCandidateCollection::const_iterator l1PFCand = l1PFCands->begin(); l1PFCand != l1PFCands->end(); ++l1PFCand )
  {
    //std::string type_string;
    //if      ( l1PFCand->id() == l1t::PFCandidate::ChargedHadron ) type_string = "PFChargedHadron";
    //else if ( l1PFCand->id() == l1t::PFCandidate::Electron      ) type_string = "PFElectron";
    //else if ( l1PFCand->id() == l1t::PFCandidate::NeutralHadron ) type_string = "PFNeutralHadron";
    //else if ( l1PFCand->id() == l1t::PFCandidate::Photon        ) type_string = "PFPhoton";
    //else if ( l1PFCand->id() == l1t::PFCandidate::Muon          ) type_string = "PFMuon";
    //else                                                          type_string = "N/A";
    //std::cout << "L1PFCandidate #" << idxL1PFCand << " (type = " << type_string << "):" 
    //	        << " pT = " << l1PFCand->pt() << ", eta = " << l1PFCand->eta() << ", phi = " << l1PFCand->phi() << ","
    //	        << " isSelected (dZcut disabled) = " << isSelected(isolationQualityCuts_dzCut_disabled_, *l1PFCand, primaryVertex_z) << std::endl;
    //++idxL1PFCand;
    if ( l1PFCand->id() == l1t::PFCandidate::ChargedHadron && isSelected(isolationQualityCuts_dzCut_enabled_primary_, *l1PFCand, primaryVertex_z) )
    {
      pfChargedHadronPlots_->fillHistograms(*l1PFCand, numPileup, evtWeight);
    }
    if ( l1PFCand->id() == l1t::PFCandidate::ChargedHadron && isSelected(isolationQualityCuts_dzCut_enabled_pileup_, *l1PFCand, primaryVertex_z) )
    {
      pfChargedHadronPileupPlots_->fillHistograms(*l1PFCand, numPileup, evtWeight);
    }
    if ( l1PFCand->id() == l1t::PFCandidate::Electron && isSelected(isolationQualityCuts_dzCut_disabled_, *l1PFCand, primaryVertex_z) )
    {
      pfElectronPlots_->fillHistograms(*l1PFCand, numPileup, evtWeight);
    }
    if ( l1PFCand->id() == l1t::PFCandidate::NeutralHadron && isSelected(isolationQualityCuts_dzCut_disabled_, *l1PFCand, primaryVertex_z) )
    {
      pfNeutralHadronPlots_->fillHistograms(*l1PFCand, numPileup, evtWeight);
    }
    if ( l1PFCand->id() == l1t::PFCandidate::Photon && isSelected(isolationQualityCuts_dzCut_disabled_, *l1PFCand, primaryVertex_z) )
    {
      pfPhotonPlots_->fillHistograms(*l1PFCand, numPileup, evtWeight);
    }
    if ( l1PFCand->id() == l1t::PFCandidate::Muon && isSelected(isolationQualityCuts_dzCut_disabled_, *l1PFCand, primaryVertex_z) )
    {
      pfMuonPlots_->fillHistograms(*l1PFCand, numPileup, evtWeight);
    }
  }

  histogram_EventCounter_->Fill(0., evtWeight);
}

void L1PFCandidateTypeAnalyzer::endJob()
{
  if ( histogram_EventCounter_->Integral() > 0. )
  {
    std::cout << "<L1PFCandidateTypeAnalyzer::endJob (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
    pfChargedHadronPlots_->normalizeHistograms(histogram_EventCounter_->Integral());
    pfChargedHadronPileupPlots_->normalizeHistograms(histogram_EventCounter_->Integral());
    pfElectronPlots_->normalizeHistograms(histogram_EventCounter_->Integral());
    pfNeutralHadronPlots_->normalizeHistograms(histogram_EventCounter_->Integral());
    std::cout << "scaling histogram = " << pfPhotonPlots_->histogram_ptFraction_vs_absEta_->GetName() << " by factor = " << (1./histogram_EventCounter_->Integral()) << std::endl;    
    pfPhotonPlots_->normalizeHistograms(histogram_EventCounter_->Integral());
    pfMuonPlots_->normalizeHistograms(histogram_EventCounter_->Integral());
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(L1PFCandidateTypeAnalyzer);


