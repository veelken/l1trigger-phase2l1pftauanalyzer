#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/PackedCandidateTypeAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h> // TMath::Abs()

PackedCandidateTypeAnalyzer::PackedCandidateTypeAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , applyPuppiWeights_(cfg.getParameter<bool>("applyPuppiWeights"))
  , pfChargedHadronPlots_(nullptr)
  , pfChargedHadronPileupPlots_(nullptr)
  , pfElectronPlots_(nullptr)
  , pfNeutralHadronPlots_(nullptr)
  , pfPhotonPlots_(nullptr)
  , pfMuonPlots_(nullptr)
{
  src_packedCands_ = cfg.getParameter<edm::InputTag>("srcPackedCands");
  token_packedCands_ = consumes<pat::PackedCandidateCollection>(src_packedCands_);

  srcOfflineVertices_ = cfg.getParameter<edm::InputTag>("srcOfflineVertices");
  if ( srcOfflineVertices_.label() != "" ) 
  {
    tokenOfflineVertices_ = consumes<reco::VertexCollection>(srcOfflineVertices_);
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

  pfChargedHadronPlots_       = new pfCandTypePlotEntryType("chargedHadron",       applyPuppiWeights_);
  pfChargedHadronPileupPlots_ = new pfCandTypePlotEntryType("chargedHadronPileup", applyPuppiWeights_);
  pfElectronPlots_            = new pfCandTypePlotEntryType("electron",            applyPuppiWeights_);
  pfNeutralHadronPlots_       = new pfCandTypePlotEntryType("neutralHadron",       applyPuppiWeights_);
  pfPhotonPlots_              = new pfCandTypePlotEntryType("photon",              applyPuppiWeights_);
  pfMuonPlots_                = new pfCandTypePlotEntryType("muon",                applyPuppiWeights_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

PackedCandidateTypeAnalyzer::~PackedCandidateTypeAnalyzer()
{}

void PackedCandidateTypeAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("PackedCandidateTypeAnalyzer") 
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

namespace
{
  const int pdgId_chargedHadron = 211;
  const int pdgId_electron      =  11;
  const int pdgId_neutralHadron = 130;
  const int pdgId_photon        =  22;
  const int pdgId_muon          =  13;

  bool passesQualityCut(const TallinnL1PFTauQualityCut& qualityCut, const pat::PackedCandidate& packedCand, float_t primaryVertex_z)
  {
    int packedCand_absPdgId = TMath::Abs(packedCand.pdgId());
    l1t::PFCandidate::Kind qualityCut_pfCandType = qualityCut.pfCandType();
    if ( (packedCand_absPdgId == pdgId_chargedHadron && qualityCut_pfCandType == l1t::PFCandidate::ChargedHadron) ||
	 (packedCand_absPdgId == pdgId_electron      && qualityCut_pfCandType == l1t::PFCandidate::Electron     ) ||
	 (packedCand_absPdgId == pdgId_neutralHadron && qualityCut_pfCandType == l1t::PFCandidate::NeutralHadron) ||
	 (packedCand_absPdgId == pdgId_photon        && qualityCut_pfCandType == l1t::PFCandidate::Photon       ) ||
	 (packedCand_absPdgId == pdgId_muon          && qualityCut_pfCandType == l1t::PFCandidate::Muon         ) )
    {  
      float_t qualityCut_min_pt = qualityCut.min_pt();
      if ( packedCand.pt() < qualityCut_min_pt ) 
      {
        return false;
      }

      int qualityCut_dzCut = qualityCut.dzCut();
      if ( packedCand.charge() != 0 )
      {
        if ( qualityCut_dzCut == TallinnL1PFTauQualityCut::kEnabled_primary || qualityCut_dzCut == TallinnL1PFTauQualityCut::kEnabled_pileup )
        {
	  double dz = std::fabs(packedCand.vertex().z() - primaryVertex_z);  
	  float_t qualityCut_max_dz = qualityCut.max_dz();
	  if ( qualityCut_dzCut == TallinnL1PFTauQualityCut::kEnabled_primary && dz >  qualityCut_max_dz ) return false;
	  if ( qualityCut_dzCut == TallinnL1PFTauQualityCut::kEnabled_pileup  && dz <= qualityCut_max_dz ) return false;
        }
      }
      else if ( qualityCut_dzCut == TallinnL1PFTauQualityCut::kEnabled_pileup )
      {
        return false; // CV: only consider charged PFCands as originating from pileup
      }
    }
    return true;
  }
    
  bool isSelected(const std::vector<TallinnL1PFTauQualityCut>& qualityCuts, const pat::PackedCandidate& packedCand, float_t primaryVertex_z)
  {
    for ( auto qualityCut : qualityCuts ) 
    {
      if ( !passesQualityCut(qualityCut, packedCand, primaryVertex_z) ) return false;
    }
    return true;
  }
}

void PackedCandidateTypeAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<pat::PackedCandidateCollection> packedCands;
  evt.getByToken(token_packedCands_, packedCands);

  float primaryVertex_z = 0.;
  if ( srcOfflineVertices_.label() != "" ) 
  {
    edm::Handle<reco::VertexCollection> vertices;
    evt.getByToken(tokenOfflineVertices_, vertices);
    if ( vertices->size() > 0 ) 
    {
      edm::Ref<reco::VertexCollection> primaryVertex = reco::VertexRef(vertices, 0);
      primaryVertex_z = primaryVertex->position().z();
    }
  }

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
  
  for ( pat::PackedCandidateCollection::const_iterator packedCand = packedCands->begin(); packedCand != packedCands->end(); ++packedCand )
  {    
    int packedCand_absPdgId = TMath::Abs(packedCand->pdgId());
    if ( packedCand_absPdgId == pdgId_chargedHadron && isSelected(isolationQualityCuts_dzCut_enabled_primary_, *packedCand, primaryVertex_z) )
    {
      pfChargedHadronPlots_->fillHistograms(*packedCand, numPileup, evtWeight);
    }
    if ( packedCand_absPdgId == pdgId_chargedHadron && isSelected(isolationQualityCuts_dzCut_enabled_pileup_, *packedCand, primaryVertex_z) )
    {
      pfChargedHadronPileupPlots_->fillHistograms(*packedCand, numPileup, evtWeight);
    }
    if ( packedCand_absPdgId == pdgId_electron && isSelected(isolationQualityCuts_dzCut_disabled_, *packedCand, primaryVertex_z) )
    {
      pfElectronPlots_->fillHistograms(*packedCand, numPileup, evtWeight);
    }
    if ( packedCand_absPdgId == pdgId_neutralHadron && isSelected(isolationQualityCuts_dzCut_disabled_, *packedCand, primaryVertex_z) )
    {
      pfNeutralHadronPlots_->fillHistograms(*packedCand, numPileup, evtWeight);
    }
    if ( packedCand_absPdgId == pdgId_photon && isSelected(isolationQualityCuts_dzCut_disabled_, *packedCand, primaryVertex_z) )
    {
      pfPhotonPlots_->fillHistograms(*packedCand, numPileup, evtWeight);
    }
    if ( packedCand_absPdgId == pdgId_muon && isSelected(isolationQualityCuts_dzCut_disabled_, *packedCand, primaryVertex_z) )
    {
      pfMuonPlots_->fillHistograms(*packedCand, numPileup, evtWeight);
    }
  }

  histogram_EventCounter_->Fill(0., evtWeight);
}

void PackedCandidateTypeAnalyzer::endJob()
{
  if ( histogram_EventCounter_->Integral() > 0. )
  {
    std::cout << "<PackedCandidateTypeAnalyzer::endJob (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
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

DEFINE_FWK_MODULE(PackedCandidateTypeAnalyzer);


