#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/L1HPSPFTauSelector.h"

#include "FWCore/Utilities/interface/InputTag.h"

L1HPSPFTauSelector::L1HPSPFTauSelector(const edm::ParameterSet& cfg) 
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<l1t::L1HPSPFTauCollection>(src_);

  min_pt_ = cfg.getParameter<double>("min_pt");
  max_pt_ = cfg.getParameter<double>("max_pt");
  min_absEta_ = cfg.getParameter<double>("min_absEta");
  max_absEta_ = cfg.getParameter<double>("max_absEta");
  min_leadTrackPt_ = cfg.getParameter<double>("min_leadTrackPt");
  max_leadTrackPt_ = cfg.getParameter<double>("max_leadTrackPt");
  min_relChargedIso_ = cfg.getParameter<double>("min_relChargedIso");
  max_relChargedIso_ = cfg.getParameter<double>("max_relChargedIso");
  min_absChargedIso_ = cfg.getParameter<double>("min_absChargedIso");
  max_absChargedIso_ = cfg.getParameter<double>("max_absChargedIso");

  invert_ = cfg.getParameter<bool>("invert");

  produces<l1t::L1HPSPFTauCollection>("");
}

L1HPSPFTauSelector::~L1HPSPFTauSelector()
{}

void L1HPSPFTauSelector::produce(edm::Event& evt, const edm::EventSetup& es)
{
  std::unique_ptr<l1t::L1HPSPFTauCollection> pfTaus_selected(new l1t::L1HPSPFTauCollection());

  edm::Handle<l1t::L1HPSPFTauCollection> pfTaus;
  evt.getByToken(token_, pfTaus);
  
  size_t numPFTaus = pfTaus->size();
  for ( size_t idxPFTau = 0; idxPFTau < numPFTaus; ++idxPFTau ) 
  {  
    const l1t::L1HPSPFTau& pfTau = pfTaus->at(idxPFTau);
    double pfTau_absEta = std::fabs(pfTau.eta());
    double sumChargedIso = pfTau.sumChargedIso();
    bool isSelected = false;
    if ( (min_pt_            < 0. || pfTau.pt()                      >=  min_pt_                       ) &&
         (max_pt_            < 0. || pfTau.pt()                      <=  max_pt_                       ) &&
         (min_absEta_        < 0. || pfTau_absEta                    >=  min_absEta_                   ) &&
         (max_absEta_        < 0. || pfTau_absEta                    <=  max_absEta_                   ) &&
         (                           pfTau.leadChargedPFCand().isNonnull()                             ) &&
         (min_leadTrackPt_   < 0. || pfTau.leadChargedPFCand()->pt() >=  min_leadTrackPt_              ) &&
         (max_leadTrackPt_   < 0. || pfTau.leadChargedPFCand()->pt() <=  max_leadTrackPt_              ) &&
         (min_relChargedIso_ < 0. || sumChargedIso                   >= (min_relChargedIso_*pfTau.pt())) &&
         (max_relChargedIso_ < 0. || sumChargedIso                   <= (max_relChargedIso_*pfTau.pt())) &&
         (min_absChargedIso_ < 0. || sumChargedIso                   >=  min_absChargedIso_            ) &&
	 (max_absChargedIso_ < 0. || sumChargedIso                   <=  max_absChargedIso_            ) )
    {
      isSelected = true;
    }
    if ( (isSelected && !invert_) || (!isSelected && invert_) ) 
    {
      pfTaus_selected->push_back(pfTau);
    }
  }

  evt.put(std::move(pfTaus_selected));
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(L1HPSPFTauSelector);
