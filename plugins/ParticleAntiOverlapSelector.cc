#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/ParticleAntiOverlapSelector.h"

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTau.h"   // l1t::TallinnL1PFTau
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h" // l1t::PFCandidate

typedef ObjectSelector<ParticleAntiOverlapSelector<l1t::TallinnL1PFTau>> TallinnL1PFTauAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<l1t::PFCandidate>> L1PFCandidateAntiOverlapSelector;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TallinnL1PFTauAntiOverlapSelector);
DEFINE_FWK_MODULE(L1PFCandidateAntiOverlapSelector);
