#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/RhoCorrAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

RhoCorrAnalyzer::RhoCorrAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_rho_ = cfg.getParameter<edm::InputTag>("srcRho");
  token_rho_ = consumes<float>(src_rho_);

  src_rhoNeutral_ = cfg.getParameter<edm::InputTag>("srcRhoNeutral");
  token_rhoNeutral_ = consumes<float>(src_rhoNeutral_);

  src_l1PFCands_ = cfg.getParameter<edm::InputTag>("srcL1PFCands");
  token_l1PFCands_ = consumes<l1t::PFCandidateCollection>(src_l1PFCands_);

  edm::ParameterSet cfg_isolationQualityCuts = cfg.getParameter<edm::ParameterSet>("isolationQualityCuts");
  isolationQualityCuts_ = readL1PFTauQualityCuts(cfg_isolationQualityCuts, "disabled", false);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

RhoCorrAnalyzer::~RhoCorrAnalyzer()
{}

void RhoCorrAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("RhoCorrAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());

  dqmStore.setCurrentFolder(dqmDirectory_);

  me_rho_ = dqmStore.book1D("rho", "rho", 200, 0., 100.);
  histogram_rho_ = me_rho_->getTH1();
  assert(histogram_rho_);

  me_rhoNeutral_ = dqmStore.book1D("rhoNeutral", "rhoNeutral", 200, 0., 100.);
  histogram_rhoNeutral_ = me_rhoNeutral_->getTH1();
  assert(histogram_rhoNeutral_);

  me_neutralPFCandPt_vs_absEta_ = dqmStore.book1D("neutralPFCandPt_vs_absEta", "neutralPFCandPt_vs_absEta", 30, 0., 3.0);
  histogram_neutralPFCandPt_vs_absEta_ = me_neutralPFCandPt_vs_absEta_->getTH1();
  assert(histogram_neutralPFCandPt_vs_absEta_);

  me_EventCounter_ = dqmStore.book1D("EventCounter", "EventCounter", 1, -0.5, +0.5);
  histogram_EventCounter_ = me_EventCounter_->getTH1();
  assert(histogram_EventCounter_);
}
    
void L1TrackAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<float> rho;
  evt.getByToken(token_rho_, rho);
  
  edm::Handle<float> rhoNeutral;
  evt.getByToken(token_rhoNeutral_, rhoNeutral);

  edm::Handle<l1t::PFCandidateCollection> l1PFCands;
  evt.getByToken(token_l1PFCands_, l1PFCands);

  const double evtWeight = 1.;

  histogram_rho_->Fill(*rho, evtWeight);
  
  histogram_rhoNeutral_->Fill(*rhoNeutral, evtWeight);
  
  for ( l1t::PFCandidateCollection::const_iterator l1PFCand = l1PFCands->begin(); l1PFCand != l1PFCands->end(); ++l1PFCand )
  {
    if ( l1PFCand->id() == l1t::PFCandidate::Photon && isSelected(isolationQualityCuts_, *l1PFCand, 0.) )
    {
      histogram_neutralPFCandPt_vs_absEta_->Fill(TMath::Abs(l1PFCand->eta()), l1PFCand->pt()*evtWeight);
    }
  }

  histogram_EventCounter_->Fill(0., evtWeight);
}

void L1TrackAnalyzer::endJob()
{
  if ( histogram_EventCounter_->Integral() > 0. )
  {
    std::cout << "<L1TrackAnalyzer::endJob (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
    std::cout << "scaling histogram = " << histogram_neutralPFCandPt_vs_absEta_->GetName() << " by factor = " << (1./histogram_EventCounter_->Integral()) << std::endl;
    histogram_neutralPFCandPt_vs_absEta_->Scale(1./histogram_EventCounter_->Integral());
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RhoCorrAnalyzer);


