#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/DumpL1HPSPFTaus.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpL1HPSPFTaus::DumpL1HPSPFTaus(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , debug_(cfg.getUntrackedParameter<bool>("debug", false))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<l1t::L1HPSPFTauCollection>(src_);
}

DumpL1HPSPFTaus::~DumpL1HPSPFTaus()
{}

namespace
{
  bool isSelected(const l1t::L1HPSPFTau& tau, bool debug)
  {
    const double min_PFTau_pt              = 20.;
    const double max_PFTau_eta             = 2.1;
    const double min_leadChargedPFCand_pt  = 5.;
    const double max_leadChargedPFCand_eta = 2.1;
    const double max_leadChargedPFCand_dz  = 0.2;
    const double max_chargedIso            = 1.e+3;
    const double max_chargedRelIso         = 1.0;
    bool retVal = true;
    if ( !(tau.pt() > min_PFTau_pt) )
    {
      if ( debug )
      {
	std::cout << " FAILS min_PFTau_pt cut." << std::endl;
      }
      retVal = false;
    }
    if ( !(std::fabs(tau.eta()) < max_PFTau_eta) )
    {
      if ( debug )
      {
	std::cout << " FAILS max_PFTau_eta cut." << std::endl;
      }
      retVal = false;
    }
    if ( !(tau.leadChargedPFCand().isNonnull() && 
	   tau.leadChargedPFCand()->pt() > min_leadChargedPFCand_pt) )
    {
      if ( debug )
      {
	std::cout << " FAILS min_leadChargedPFCand_pt cut." << std::endl;
      }
      retVal = false;
    }
    if ( !(tau.leadChargedPFCand().isNonnull() && 
	   std::fabs(tau.leadChargedPFCand()->eta()) < max_leadChargedPFCand_eta) )
    {
      if ( debug )
      {
	std::cout << " FAILS max_leadChargedPFCand_eta cut." << std::endl;
      }
      retVal = false;
    }
    if ( tau.leadChargedPFCand().isNonnull() ) 
    {
      std::cout << "tau.leadChargedPFCand()->pfTrack()->vertex().z() = " << tau.leadChargedPFCand()->pfTrack()->vertex().z() << std::endl;
      if ( tau.primaryVertex().isNonnull() ) 
      {
        std::cout << "tau.primaryVertex()->zvertex() = " << tau.primaryVertex()->zvertex() << std::endl;
      } 
      else 
      {
        std::cout << "tau.primaryVertex()->zvertex() = N/A" << std::endl; 
      }
    }
    if ( !(tau.primaryVertex().isNonnull() && tau.leadChargedPFCand().isNonnull() && 
	   std::fabs(tau.leadChargedPFCand()->pfTrack()->vertex().z() - tau.primaryVertex()->zvertex()) < max_leadChargedPFCand_dz) )
    { 
      if ( debug )
      {
	std::cout << " FAILS max_leadChargedPFCand_dz cut." << std::endl;
      }
      retVal = false;
    }
    if ( !(tau.sumChargedIso() < max_chargedIso) )
    {
      if ( debug )
      {
	std::cout << " FAILS max_chargedIso cut." << std::endl;
      }
      retVal = false;
    }
    if ( !(tau.sumChargedIso() < max_chargedRelIso*tau.pt()) )
    {
      if ( debug )
      {
	std::cout << " FAILS max_chargedRelIso cut." << std::endl;
      }
      retVal = false;
    }
    if ( debug && retVal )
    {
      std::cout << " PASSES selection." << std::endl;
    }
    return retVal;
  }
}

void DumpL1HPSPFTaus::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpL1HPSPFTaus::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<l1t::L1HPSPFTauCollection> taus;
  evt.getByToken(token_, taus);
  
  size_t numTaus = taus->size();
  for ( size_t idxTau = 0; idxTau < numTaus; ++idxTau ) 
  {
    const l1t::L1HPSPFTau& tau = taus->at(idxTau);
    std::cout << "L1HPSPFTau #" << idxTau << ": " << tau;
    if ( debug_ ) 
    {
      if ( !isSelected(tau, debug_) )
      {
	std::cout << "--> CHECK !!" << std::endl;
      }
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpL1HPSPFTaus);





