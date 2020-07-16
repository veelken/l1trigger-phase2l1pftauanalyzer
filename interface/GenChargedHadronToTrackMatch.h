#ifndef L1Trigger_TallinnL1PFTauAnalyzer_GenChargedHadronToTrackMatch_h
#define L1Trigger_TallinnL1PFTauAnalyzer_GenChargedHadronToTrackMatch_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"   // reco::GenParticle
#include "DataFormats/TrackReco/interface/Track.h"              // reco::Track
#include "DataFormats/L1TParticleFlow/interface/PFTrack.h"      // l1t::PFTrack::L1TTTrackType

namespace l1t
{
  typedef PFTrack::L1TTTrackType Track;
  typedef std::vector<Track> TrackCollection;
}

class GenChargedHadronToTrackMatchBase
{
 public:
  GenChargedHadronToTrackMatchBase(const reco::Candidate* genChargedHadron);
  virtual ~GenChargedHadronToTrackMatchBase();

  const reco::Candidate* genChargedHadron() const;
  bool hasGenChargedHadron() const;
  double genChargedHadron_pt() const;
  double genChargedHadron_eta() const;
  double genChargedHadron_absEta() const;
  double genChargedHadron_phi() const;

  bool hasRecTrack() const;
  double recTrack_pt() const;
  double recTrack_eta() const;
  double recTrack_absEta() const;
  double recTrack_phi() const;

 protected:
  const reco::Candidate* genChargedHadron_;
  bool hasGenChargedHadron_;
  double genChargedHadron_pt_;
  double genChargedHadron_eta_;
  double genChargedHadron_absEta_;
  double genChargedHadron_phi_;
  
  bool hasRecTrack_;
  double recTrack_pt_;
  double recTrack_eta_;
  double recTrack_absEta_;
  double recTrack_phi_;
};

class GenChargedHadronToOfflineTrackMatch : public GenChargedHadronToTrackMatchBase
{
 public:
  GenChargedHadronToOfflineTrackMatch(const reco::Candidate* genChargedHadron, const reco::Track* recTrack);
  ~GenChargedHadronToOfflineTrackMatch();
  
  const reco::Track* recTrack() const;

  bool isOverlap(const GenChargedHadronToTrackMatchBase* other);

 private:
  const reco::Track* recTrack_;
};

bool isOverlap(const GenChargedHadronToOfflineTrackMatch& match1, const GenChargedHadronToOfflineTrackMatch& match2);

class GenChargedHadronToL1TrackMatch : public GenChargedHadronToTrackMatchBase
{
 public:
  GenChargedHadronToL1TrackMatch(const reco::Candidate* genChargedHadron, const l1t::Track* recTrack);
  ~GenChargedHadronToL1TrackMatch();
  
  const l1t::Track* recTrack() const;

 private:
  const l1t::Track* recTrack_;
};

bool isOverlap(const GenChargedHadronToL1TrackMatch& match1, const GenChargedHadronToL1TrackMatch& match2);

#endif
