#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/DumpTallinL1PFTaus.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpTallinL1PFTaus::DumpTallinL1PFTaus(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<l1t::TallinnL1PFTauCollection>(src_);
}

DumpTallinL1PFTaus::~DumpTallinL1PFTaus()
{}

void DumpTallinL1PFTaus::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpTallinL1PFTaus::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<l1t::TallinnL1PFTauCollection> taus;
  evt.getByToken(token_, taus);
  
  size_t numTaus = taus->size();
  for ( size_t idxTau = 0; idxTau < numTaus; ++idxTau ) 
  {
    const l1t::TallinnL1PFTau& tau = taus->at(idxTau);
    std::cout << "TallinnL1PFTau #" << idxTau << ": " << tau;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpTallinL1PFTaus);





