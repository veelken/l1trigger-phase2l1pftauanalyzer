#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/DumpL1PFCandidates.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpL1PFCandidates::DumpL1PFCandidates(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<l1t::PFCandidateCollection>(src_);
  minPt_ = ( cfg.exists("min_pt") ) ? cfg.getParameter<double>("min_pt") : -1.;
  printPuppiWeights_ = cfg.getParameter<bool>("printPuppiWeights");
}

DumpL1PFCandidates::~DumpL1PFCandidates()
{}

void DumpL1PFCandidates::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpL1PFCandidates::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<l1t::PFCandidateCollection> l1PFCands;
  evt.getByToken(token_, l1PFCands);
  
  size_t numL1PFCands = l1PFCands->size();
  for ( size_t idxL1PFCand = 0; idxL1PFCand < numL1PFCands; ++idxL1PFCand ) 
  {
    const l1t::PFCandidate& l1PFCand = l1PFCands->at(idxL1PFCand);
    if ( !(l1PFCand.pt() > minPt_) ) continue;
    std::string type_string;
    if      ( l1PFCand.id() == l1t::PFCandidate::ChargedHadron ) type_string = "PFChargedHadron";
    else if ( l1PFCand.id() == l1t::PFCandidate::Electron      ) type_string = "PFElectron";
    else if ( l1PFCand.id() == l1t::PFCandidate::NeutralHadron ) type_string = "PFNeutralHadron";
    else if ( l1PFCand.id() == l1t::PFCandidate::Photon        ) type_string = "PFPhoton";
    else if ( l1PFCand.id() == l1t::PFCandidate::Muon          ) type_string = "PFMuon";
    else                                                         type_string = "N/A";
    std::cout << "L1PFCandidate #" << idxL1PFCand << " (type = " << type_string << "):" 
	      << " pT = " << l1PFCand.pt() << ", eta = " << l1PFCand.eta() << ", phi = " << l1PFCand.phi();
    if ( printPuppiWeights_ )
    {
      std::cout << " (PUPPI weight = " << l1PFCand.puppiWeight() << ")";
    }
    std::cout << std::endl;
    if ( l1PFCand.pfTrack().isNonnull() )
    {
      std::cout << " PFTrack: pT = " << l1PFCand.pfTrack()->pt() << ", eta = " << l1PFCand.pfTrack()->eta() << ", phi = " << l1PFCand.pfTrack()->phi() << std::endl;
    }
    if ( l1PFCand.pfCluster().isNonnull() )
    {
      std::cout << " PFCluster: pT = " << l1PFCand.pfCluster()->pt() << ", eta = " << l1PFCand.pfCluster()->eta() << ", phi = " << l1PFCand.pfCluster()->phi() 
		<< " (isEM = " << l1PFCand.pfCluster()->isEM() << ", H/E = " << l1PFCand.pfCluster()->hOverE() << ")" << std::endl;
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpL1PFCandidates);





