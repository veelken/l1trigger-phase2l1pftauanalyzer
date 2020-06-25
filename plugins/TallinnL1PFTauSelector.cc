#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/TallinnL1PFTauSelector.h"

#include <cmath>

TallinnL1PFTauSelector::TallinnL1PFTauSelector(const edm::ParameterSet& cfg)
 : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<l1t::TallinnL1PFTauCollection>(src_);
 
  min_pt_ = cfg.getParameter<double>("min_pt");
  max_pt_ = cfg.getParameter<double>("max_pt");
  min_absEta_ = cfg.getParameter<double>("min_absEta");
  max_absEta_ = cfg.getParameter<double>("max_absEta");
  min_leadTrack_pt_ = cfg.getParameter<double>("min_leadTrack_pt");
  max_leadTrack_pt_ = cfg.getParameter<double>("max_leadTrack_pt");
  max_relChargedIso_ = cfg.getParameter<double>("max_relChargedIso");
  max_absChargedIso_ = cfg.getParameter<double>("max_absChargedIso");
}

TallinnL1PFTauSelector::~TallinnL1PFTauSelector()
{}

void 
TallinnL1PFTauSelector::produce(edm::Event& evt, const edm::EventSetup& es)
{
  std::unique_ptr<l1t::TallinnL1PFTauCollection> l1PFTaus_selected(new l1t::TallinnL1PFTauCollection());

  edm::Handle<l1t::TallinnL1PFTauCollection> l1PFTaus;
  evt.getByToken(token_, l1PFTaus);

  for ( l1t::TallinnL1PFTauCollection::const_iterator l1PFTau = l1PFTaus->begin();
        l1PFTau != l1PFTaus->end(); ++l1PFTau ) {
    if ( (min_pt_            < 0. || l1PFTau->pt()                      >=  min_pt_                          ) &&
         (max_pt_            < 0. || l1PFTau->pt()                      <=  max_pt_                          ) &&
         (min_absEta_        < 0. || std::fabs(l1PFTau->eta())          >=  min_absEta_                      ) &&
         (max_absEta_        < 0. || std::fabs(l1PFTau->eta())          <=  max_absEta_                      ) &&
         (                           l1PFTau->leadChargedPFCand().isNonnull()                                ) && 
         (min_leadTrack_pt_  < 0. || l1PFTau->leadChargedPFCand()->pt() >=  min_leadTrack_pt_                ) &&
         (max_leadTrack_pt_  < 0. || l1PFTau->leadChargedPFCand()->pt() <=  max_leadTrack_pt_                ) &&
         (max_relChargedIso_ < 0. || l1PFTau->sumChargedIso()           <= (max_relChargedIso_*l1PFTau->pt())) &&
         (max_absChargedIso_ < 0. || l1PFTau->sumChargedIso()           <=  max_absChargedIso_               ) ) 
    {
      l1PFTaus_selected->push_back(*l1PFTau);
    }
  }

  evt.put(std::move(l1PFTaus_selected));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TallinnL1PFTauSelector);

