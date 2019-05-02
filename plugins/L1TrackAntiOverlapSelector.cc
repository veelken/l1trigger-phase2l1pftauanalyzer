#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/L1TrackAntiOverlapSelector.h"

L1TrackAntiOverlapSelector::L1TrackAntiOverlapSelector(const edm::ParameterSet& cfg) 
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<L1TrackCollection>(src_);
  srcNotToBeFiltered_ = cfg.getParameter<vInputTag>("srcNotToBeFiltered");
  for ( auto it : srcNotToBeFiltered_ )
  {
    tokensNotToBeFiltered_.push_back(consumes<reco::CandidateView>(it));
  }
  dRmin_ = cfg.getParameter<double>("dRmin");
  invert_ = ( cfg.exists("invert") ) ?
    cfg.getParameter<bool>("invert") : false;

  produces<L1TrackCollection>();
}

L1TrackAntiOverlapSelector::~L1TrackAntiOverlapSelector()
{}

void L1TrackAntiOverlapSelector::produce(edm::Event& evt, const edm::EventSetup& es)
{
  std::unique_ptr<L1TrackCollection> l1Tracks_selected(new L1TrackCollection());

  edm::Handle<L1TrackCollection> l1Tracks;
  evt.getByToken(token_, l1Tracks);

  std::vector<bool> isOverlap(l1Tracks->size());
    
  for ( auto it : tokensNotToBeFiltered_ )
  {
    edm::Handle<reco::CandidateView> particlesNotToBeFiltered;
    evt.getByToken(it, particlesNotToBeFiltered);
      
    for ( reco::CandidateView::const_iterator particleNotToBeFiltered = particlesNotToBeFiltered->begin(); 
	  particleNotToBeFiltered != particlesNotToBeFiltered->end();  ++particleNotToBeFiltered ) 
    {
      size_t numL1Tracks = l1Tracks->size();
      for ( size_t idxL1Track = 0; idxL1Track < numL1Tracks; ++idxL1Track )
      {
	const L1Track& l1Track = l1Tracks->at(idxL1Track);
	const unsigned nParam = 4;
	double l1Track_eta = l1Track.getMomentum(nParam).eta();
	double l1Track_phi = l1Track.getMomentum(nParam).phi();
	double dR = reco::deltaR(l1Track_eta, l1Track_phi, particleNotToBeFiltered->eta(), particleNotToBeFiltered->phi());	  
	if ( dR < dRmin_ ) isOverlap[idxL1Track] = true;
      }
    }
  }
    
  size_t numL1Tracks = l1Tracks->size();
  for ( size_t idxL1Track = 0; idxL1Track < numL1Tracks; ++idxL1Track )
  {
    const L1Track& l1Track = l1Tracks->at(idxL1Track);
    if ( (invert_ == false && !isOverlap[idxL1Track]) ||
	 (invert_ == true  &&  isOverlap[idxL1Track]) ) 
    {
      l1Tracks_selected->push_back(l1Track); 
    }
  }

  evt.put(std::move(l1Tracks_selected));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TrackAntiOverlapSelector);
