#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/L1PFTauPairAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"      // edm::Handle
#include "DataFormats/Math/interface/deltaR.h"        // reco::deltaR

#include "TMath.h"   // TMath::Abs()

#include <iostream>
#include <iomanip>
#include <algorithm> // std::sort

L1PFTauPairAnalyzer::L1PFTauPairAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , min_refTau_pt_(-1.)
  , max_refTau_eta_(-1.)
  , dRmatch_(0.3)
{
  srcL1PFTaus_ = cfg.getParameter<edm::InputTag>("srcL1PFTaus");
  tokenL1PFTaus_ = consumes<l1t::L1PFTauCollection>(srcL1PFTaus_);
  srcRefTaus_ = cfg.getParameter<edm::InputTag>("srcRefTaus");
  if ( srcRefTaus_.label() != "" )
  {
    tokenRefTaus_ = consumes<reco::CandidateView>(srcRefTaus_);
    min_refTau_pt_ = cfg.getParameter<double>("min_refTau_pt");
    max_refTau_eta_ = cfg.getParameter<double>("max_refTau_eta");
  }

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

L1PFTauPairAnalyzer::~L1PFTauPairAnalyzer()
{
  for ( auto efficiency_or_ratePlot : efficiency_or_ratePlots_ ) 
  {
    delete efficiency_or_ratePlot;
  }
}

void L1PFTauPairAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("L1PFTauPairAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(1.0, "vLoose")); // vLoose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(1.0, "Loose"));  // Loose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(1.0, "Medium")); // Medium
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(1.0, "Tight"));  // Tight
  
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(1.4, "vLoose")); // vLoose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(1.4, "Loose"));  // Loose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(1.4, "Medium")); // Medium
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(1.4, "Tight"));  // Tight

  for ( auto efficiency_or_ratePlot : efficiency_or_ratePlots_ ) 
  {
    efficiency_or_ratePlot->bookHistograms(dqmStore);
  }
}

namespace
{
  bool
  isHigherPt(const l1t::L1PFTau* l1PFTau1,
	     const l1t::L1PFTau* l1PFTau2)
  {
    return l1PFTau1->pt() > l1PFTau2->pt();
  }

  bool
  isGenMatched(const l1t::L1PFTau& l1PFTau, const std::vector<const reco::Candidate*>& refTaus, double dRmatch)
  {
    for ( std::vector<const reco::Candidate*>::const_iterator refTau = refTaus.begin();
	  refTau != refTaus.end(); ++refTau ) {
      double dR = reco::deltaR(l1PFTau.eta(), l1PFTau.phi(), (*refTau)->eta(), (*refTau)->phi());
      if ( dR < dRmatch ) return true;
    }
    return false;
  }
}

void L1PFTauPairAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<l1t::L1PFTauCollection> l1PFTaus;
  evt.getByToken(tokenL1PFTaus_, l1PFTaus);

  std::vector<const l1t::L1PFTau*> l1PFTaus_sorted;
  for ( l1t::L1PFTauCollection::const_iterator l1PFTau = l1PFTaus->begin();
	l1PFTau != l1PFTaus->end(); ++l1PFTau ) {
    l1PFTaus_sorted.push_back(&(*l1PFTau));
  }
  std::sort(l1PFTaus_sorted.begin(), l1PFTaus_sorted.end(), isHigherPt);

  l1t::L1PFTauPairCollection l1PFTauPairs;
  if ( srcRefTaus_.label() != "" ) 
  {
    edm::Handle<reco::CandidateView> refTaus;
    evt.getByToken(tokenRefTaus_, refTaus);

    std::vector<const reco::Candidate*> refTaus_passingAbsEtaAndPt;
    for ( reco::CandidateView::const_iterator refTau = refTaus->begin();
	  refTau != refTaus->end(); ++refTau ) {
      if ( refTau->pt() > min_refTau_pt_ && TMath::Abs(refTau->eta()) < max_refTau_eta_ )
      {
	refTaus_passingAbsEtaAndPt.push_back(&(*refTau));
      }
    }
    if ( !(refTaus_passingAbsEtaAndPt.size() >= 2) ) return;

    for ( std::vector<const l1t::L1PFTau*>::const_iterator leadL1PFTau = l1PFTaus_sorted.begin();
	  leadL1PFTau != l1PFTaus_sorted.end(); ++leadL1PFTau ) {
      if ( !isGenMatched(**leadL1PFTau, refTaus_passingAbsEtaAndPt, dRmatch_) ) continue;
      for ( std::vector<const l1t::L1PFTau*>::const_iterator subleadL1PFTau = leadL1PFTau + 1;
	    subleadL1PFTau != l1PFTaus_sorted.end(); ++subleadL1PFTau ) {
	if ( !isGenMatched(**subleadL1PFTau, refTaus_passingAbsEtaAndPt, dRmatch_) ) continue;
	l1PFTauPairs.push_back(l1t::L1PFTauPair(*leadL1PFTau, *subleadL1PFTau));
      }
    }
  } 
  else
  {
    for ( std::vector<const l1t::L1PFTau*>::const_iterator leadL1PFTau = l1PFTaus_sorted.begin();
	  leadL1PFTau != l1PFTaus_sorted.end(); ++leadL1PFTau ) {
      for ( std::vector<const l1t::L1PFTau*>::const_iterator subleadL1PFTau = leadL1PFTau + 1;
	    subleadL1PFTau != l1PFTaus_sorted.end(); ++subleadL1PFTau ) {
	l1PFTauPairs.push_back(l1t::L1PFTauPair(*leadL1PFTau, *subleadL1PFTau));
      }
    }
  }

  const double evtWeight = 1.;

  for ( auto efficiency_or_ratePlot : efficiency_or_ratePlots_ ) 
  {
    efficiency_or_ratePlot->fillHistograms(l1PFTauPairs, evtWeight);
  }
}

void L1PFTauPairAnalyzer::endJob()
{
  for ( auto efficiency_or_ratePlot : efficiency_or_ratePlots_ ) 
  {
    efficiency_or_ratePlot->normalizeHistograms();
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(L1PFTauPairAnalyzer);
