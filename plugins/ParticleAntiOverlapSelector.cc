#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/ParticleAntiOverlapSelector.h"

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/Phase2L1Taus/interface/L1HPSPFTau.h"     // l1t::L1HPSPFTau
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h" // l1t::PFCandidate

typedef ObjectSelector<ParticleAntiOverlapSelector<l1t::L1HPSPFTau>> L1HPSPFTauAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<l1t::PFCandidate>> L1PFCandidateAntiOverlapSelector;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(L1HPSPFTauAntiOverlapSelector);
DEFINE_FWK_MODULE(L1PFCandidateAntiOverlapSelector);
