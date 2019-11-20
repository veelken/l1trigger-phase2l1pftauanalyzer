#include "L1Trigger/TallinnL1PFTauAnalyzer/interface/GenChargedHadronToTrackMatch.h"

#include <TMath.h> // TMath::Abs()

//-------------------------------------------------------------------------------
// Implementation of GenChargedHadronToTrackMatchBase class
GenChargedHadronToTrackMatchBase::GenChargedHadronToTrackMatchBase(const reco::Candidate* genChargedHadron)
  : genChargedHadron_(genChargedHadron)
  , hasRecTrack_(false)
  , recTrack_pt_(0.)
  , recTrack_eta_(0.)
  , recTrack_absEta_(0.)
  , recTrack_phi_(0.)
{
  assert(genChargedHadron_);
  hasGenChargedHadron_  = true;
  genChargedHadron_pt_  = genChargedHadron_->pt();
  genChargedHadron_eta_ = genChargedHadron_->eta();
  genChargedHadron_absEta_ = TMath::Abs(genChargedHadron_eta_);
  genChargedHadron_phi_ = genChargedHadron_->phi();
}
   
GenChargedHadronToTrackMatchBase::~GenChargedHadronToTrackMatchBase()
{}

const reco::Candidate* GenChargedHadronToTrackMatchBase::genChargedHadron() const
{
  return genChargedHadron_;
}

bool GenChargedHadronToTrackMatchBase::hasGenChargedHadron() const
{
  return hasGenChargedHadron_;
}
 
double GenChargedHadronToTrackMatchBase::genChargedHadron_pt() const
{
  return genChargedHadron_pt_;
}
 
double GenChargedHadronToTrackMatchBase::genChargedHadron_eta() const
{
  return genChargedHadron_eta_;
}

double GenChargedHadronToTrackMatchBase::genChargedHadron_absEta() const
{
  return genChargedHadron_absEta_;
}
 
double GenChargedHadronToTrackMatchBase::genChargedHadron_phi() const
{
  return genChargedHadron_phi_;
}

bool GenChargedHadronToTrackMatchBase::hasRecTrack() const
{
  return hasRecTrack_;
}

double GenChargedHadronToTrackMatchBase::recTrack_pt() const
{
  return recTrack_pt_;
}

double GenChargedHadronToTrackMatchBase::recTrack_eta() const
{
  return recTrack_eta_;
}

double GenChargedHadronToTrackMatchBase::recTrack_absEta() const
{
  return recTrack_absEta_;
}
 
double GenChargedHadronToTrackMatchBase::recTrack_phi() const
{
  return recTrack_phi_;
}
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
// Implementation of GenChargedHadronToOfflineTrackMatch class
GenChargedHadronToOfflineTrackMatch::GenChargedHadronToOfflineTrackMatch(const reco::Candidate* genChargedHadron, const reco::Track* recTrack)
  : GenChargedHadronToTrackMatchBase(genChargedHadron)
  , recTrack_(recTrack)
{
  if ( recTrack_ ) 
  {
    hasRecTrack_  = true;
    recTrack_pt_  = recTrack_->pt();
    recTrack_eta_ = recTrack_->eta();
    recTrack_absEta_ = TMath::Abs(recTrack_eta_);
    recTrack_phi_ = recTrack_->phi();
  }
}

GenChargedHadronToOfflineTrackMatch::~GenChargedHadronToOfflineTrackMatch() 
{}
  
const reco::Track* GenChargedHadronToOfflineTrackMatch::recTrack() const
{
  return recTrack_;
}

bool isOverlap(const GenChargedHadronToOfflineTrackMatch& match1, const GenChargedHadronToOfflineTrackMatch& match2)
{
  if ( (match1.genChargedHadron() == match2.genChargedHadron()) || (match1.recTrack() != nullptr && match1.recTrack() == match2.recTrack()) ) return true;
  else return false;
}
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
// Implementation of GenChargedHadronToL1TrackMatch class
GenChargedHadronToL1TrackMatch::GenChargedHadronToL1TrackMatch(const reco::Candidate* genChargedHadron, const l1t::Track* recTrack)
  : GenChargedHadronToTrackMatchBase(genChargedHadron)
  , recTrack_(recTrack)
{
  if ( recTrack_ )
  {
    hasRecTrack_  = true;
    const unsigned nParam = 4;
    recTrack_pt_  = recTrack_->getMomentum(nParam).perp();
    recTrack_eta_ = recTrack_->getMomentum(nParam).eta();
    recTrack_absEta_ = TMath::Abs(recTrack_eta_);
    recTrack_phi_ = recTrack_->getMomentum(nParam).phi();
  }
}

GenChargedHadronToL1TrackMatch::~GenChargedHadronToL1TrackMatch() 
{}
  
const l1t::Track* GenChargedHadronToL1TrackMatch::recTrack() const
{
  return recTrack_;
}

bool isOverlap(const GenChargedHadronToL1TrackMatch& match1, const GenChargedHadronToL1TrackMatch& match2)
{
  if ( (match1.genChargedHadron() == match2.genChargedHadron()) || (match1.recTrack() != nullptr && match1.recTrack() == match2.recTrack()) ) return true;
  else return false;
}
//-------------------------------------------------------------------------------
