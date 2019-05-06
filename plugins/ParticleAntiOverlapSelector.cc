#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/ParticleAntiOverlapSelector.h"

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Tau.h"                // pat::Tau
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"    // pat::PackedCandidate
#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTau.h"   // l1t::TallinnL1PFTau
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h" // l1t::PFCandidate
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"       // reco::GenParticle
#include "DataFormats/JetReco/interface/GenJet.h"                   // reco::GenJet

typedef ObjectSelector<ParticleAntiOverlapSelector<pat::Tau>> PATTauAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<pat::PackedCandidate>> PackedCandidateAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<l1t::TallinnL1PFTau>> TallinnL1PFTauAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<l1t::PFCandidate>> L1PFCandidateAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<reco::GenParticle>> GenParticleAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<reco::GenJet>> GenJetAntiOverlapSelector;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATTauAntiOverlapSelector);
DEFINE_FWK_MODULE(PackedCandidateAntiOverlapSelector);
DEFINE_FWK_MODULE(TallinnL1PFTauAntiOverlapSelector);
DEFINE_FWK_MODULE(L1PFCandidateAntiOverlapSelector);
DEFINE_FWK_MODULE(GenParticleAntiOverlapSelector);
DEFINE_FWK_MODULE(GenJetAntiOverlapSelector);
