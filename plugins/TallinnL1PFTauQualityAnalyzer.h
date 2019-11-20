#ifndef L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauQualityAnalyzer_h
#define L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauQualityAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTau.h"             // l1t::TallinnL1PFTau
#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTauFwd.h"          // l1t::TallinnL1PFTauCollection
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"           // l1t::PFCandidate
#include "DataFormats/Phase2L1ParticleFlow/interface/PFTrack.h"               // l1t::PFTrack
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"                     // TTTrack
#include "L1Trigger/TallinnL1PFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverflow

#include <TH1.h>     // TH1
#include <TH2.h>     // TH2
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Sqrt(), TMath::Abs(), TMath::Pi()

#include <vector>    // std::vector
#include <string>    // std::string

reco::Candidate::LorentzVector getP4_with_massConstraint(const reco::Candidate::LorentzVector& p4, double mass)
{
  double px = p4.px();
  double py = p4.py();
  double pz = p4.pz();
  double energy = TMath::Sqrt(px*px + py*py + pz*pz + mass*mass);
  reco::Candidate::LorentzVector p4_with_massConstraint(px, py, pz, energy);
  return p4_with_massConstraint;
}

reco::Candidate::LorentzVector getChargedPionP4(const reco::Candidate::LorentzVector& p4)
{
  const double chargedPionMass = 0.1396; // GeV
  return getP4_with_massConstraint(p4, chargedPionMass);
}

reco::Candidate::LorentzVector getNeutralPionP4(const reco::Candidate::LorentzVector& p4)
{
  const double neutralPionMass = 0.1350; // GeV
  return getP4_with_massConstraint(p4, neutralPionMass);
}

class TallinnL1PFTauQualityAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit TallinnL1PFTauQualityAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~TallinnL1PFTauQualityAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<l1t::TallinnL1PFTauCollection> token_;

  std::string dqmDirectory_;

  struct plotEntryType
  {
    plotEntryType(double min_pt, double max_pt, double min_absEta, double max_absEta)
      : me_leadTrack_pt_(nullptr)
      , histogram_leadTrack_pt_(nullptr)
      , me_leadTrack_chi2_(nullptr)
      , histogram_leadTrack_chi2_(nullptr)
      , me_leadTrack_numStubs_(nullptr)
      , histogram_leadTrack_numStubs_(nullptr)
      , me_signalTrack_pt_(nullptr)
      , histogram_signalTrack_pt_(nullptr)
      , me_signalTrack_chi2_(nullptr)
      , histogram_signalTrack_chi2_(nullptr)
      , me_signalTrack_numStubs_(nullptr)
      , histogram_signalTrack_numStubs_(nullptr)
      , me_isolationTrack_pt_(nullptr)
      , histogram_isolationTrack_pt_(nullptr)
      , me_isolationTrack_chi2_(nullptr)
      , histogram_isolationTrack_chi2_(nullptr)
      , me_isolationTrack_numStubs_(nullptr)
      , histogram_isolationTrack_numStubs_(nullptr)
      , me_mTracks_(nullptr)
      , histogram_mTracks_(nullptr)
      , me_mTracks_plus_Strip_(nullptr)
      , histogram_mTracks_plus_Strip_(nullptr)
      , me_strip_pt_(nullptr)
      , histogram_strip_pt_(nullptr)
      , me_strip_numPhotons_(nullptr)
      , histogram_strip_numPhotons_(nullptr)
      , min_pt_(min_pt)
      , max_pt_(max_pt)
      , min_absEta_(min_absEta)
      , max_absEta_(max_absEta)
    {}
    ~plotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix;
      if      ( min_pt_     >= 0. && max_pt_     > 0. ) histogramName_suffix.Append(Form("_pt%1.2fto%1.2f", min_pt_, max_pt_));
      else if ( min_pt_     >= 0.                     ) histogramName_suffix.Append(Form("_ptGt%1.2f", min_pt_));
      else if (                      max_pt_     > 0. ) histogramName_suffix.Append(Form("_ptLt%1.2f", max_pt_));
      if      ( min_absEta_ >= 0. && max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEta%1.2fto%1.2f", min_absEta_, max_absEta_));
      else if ( min_absEta_ >= 0.                     ) histogramName_suffix.Append(Form("_absEtaGt%1.2f", min_absEta_));
      else if (                      max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f", max_absEta_));
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_leadTrack_pt = Form("leadTrack_pt_%s", histogramName_suffix.Data());
      me_leadTrack_pt_ = dqmStore.book1D(histogramName_leadTrack_pt.Data(), histogramName_leadTrack_pt.Data(), 100, 0., 100.);
      histogram_leadTrack_pt_ = me_leadTrack_pt_->getTH1();
      assert(histogram_leadTrack_pt_);
      TString histogramName_leadTrack_chi2 = Form("leadTrack_chi2_%s", histogramName_suffix.Data());
      me_leadTrack_chi2_ = dqmStore.book1D(histogramName_leadTrack_chi2.Data(), histogramName_leadTrack_chi2.Data(), 100, 0., 100.);
      histogram_leadTrack_chi2_ = me_leadTrack_chi2_->getTH1();
      assert(histogram_leadTrack_chi2_);
      TString histogramName_leadTrack_numStubs = Form("leadTrack_numStubs_%s", histogramName_suffix.Data());
      me_leadTrack_numStubs_ = dqmStore.book1D(histogramName_leadTrack_numStubs.Data(), histogramName_leadTrack_numStubs.Data(), 25, -0.5, 24.5);
      histogram_leadTrack_numStubs_ = me_leadTrack_numStubs_->getTH1();
      assert(histogram_leadTrack_numStubs_);

      TString histogramName_signalTrack_pt = Form("signalTrack_pt_%s", histogramName_suffix.Data());
      me_signalTrack_pt_ = dqmStore.book1D(histogramName_signalTrack_pt.Data(), histogramName_signalTrack_pt.Data(), 100, 0., 100.);
      histogram_signalTrack_pt_ = me_signalTrack_pt_->getTH1();
      assert(histogram_signalTrack_pt_);
      TString histogramName_signalTrack_chi2 = Form("signalTrack_chi2_%s", histogramName_suffix.Data());
      me_signalTrack_chi2_ = dqmStore.book1D(histogramName_signalTrack_chi2.Data(), histogramName_signalTrack_chi2.Data(), 100, 0., 100.);
      histogram_signalTrack_chi2_ = me_signalTrack_chi2_->getTH1();
      assert(histogram_signalTrack_chi2_);
      TString histogramName_signalTrack_numStubs = Form("signalTrack_numStubs_%s", histogramName_suffix.Data());
      me_signalTrack_numStubs_ = dqmStore.book1D(histogramName_signalTrack_numStubs.Data(), histogramName_signalTrack_numStubs.Data(), 25, -0.5, 24.5);
      histogram_signalTrack_numStubs_ = me_signalTrack_numStubs_->getTH1();
      assert(histogram_signalTrack_numStubs_);

      TString histogramName_isolationTrack_pt = Form("isolationTrack_pt_%s", histogramName_suffix.Data());
      me_isolationTrack_pt_ = dqmStore.book1D(histogramName_isolationTrack_pt.Data(), histogramName_isolationTrack_pt.Data(), 100, 0., 100.);
      histogram_isolationTrack_pt_ = me_isolationTrack_pt_->getTH1();
      assert(histogram_isolationTrack_pt_);
      TString histogramName_isolationTrack_chi2 = Form("isolationTrack_chi2_%s", histogramName_suffix.Data());
      me_isolationTrack_chi2_ = dqmStore.book1D(histogramName_isolationTrack_chi2.Data(), histogramName_isolationTrack_chi2.Data(), 100, 0., 100.);
      histogram_isolationTrack_chi2_ = me_isolationTrack_chi2_->getTH1();
      assert(histogram_isolationTrack_chi2_);
      TString histogramName_isolationTrack_numStubs = Form("isolationTrack_numStubs_%s", histogramName_suffix.Data());
      me_isolationTrack_numStubs_ = dqmStore.book1D(histogramName_isolationTrack_numStubs.Data(), histogramName_isolationTrack_numStubs.Data(), 25, -0.5, 24.5);
      histogram_isolationTrack_numStubs_ = me_isolationTrack_numStubs_->getTH1();
      assert(histogram_isolationTrack_numStubs_);

      TString histogramName_mTracks = Form("mTracks_%s", histogramName_suffix.Data());
      me_mTracks_ = dqmStore.book1D(histogramName_mTracks.Data(), histogramName_mTracks.Data(), 100, 0., 10.);
      histogram_mTracks_ = me_mTracks_->getTH1();
      assert(histogram_mTracks_);
      TString histogramName_mTracks_plus_Strip = Form("mTracks_plus_Strip_%s", histogramName_suffix.Data());
      me_mTracks_plus_Strip_ = dqmStore.book1D(histogramName_mTracks_plus_Strip.Data(), histogramName_mTracks_plus_Strip.Data(), 100, 0., 10.);
      histogram_mTracks_plus_Strip_ = me_mTracks_plus_Strip_->getTH1();
      assert(histogram_mTracks_plus_Strip_);

      TString histogramName_strip_pt = Form("strip_pt_%s", histogramName_suffix.Data());
      me_strip_pt_ = dqmStore.book1D(histogramName_strip_pt.Data(), histogramName_strip_pt.Data(), 100, 0., 100.);
      histogram_strip_pt_ = me_strip_pt_->getTH1();
      assert(histogram_strip_pt_);
      TString histogramName_strip_numPhotons = Form("strip_numPhotons_%s", histogramName_suffix.Data());
      me_strip_numPhotons_ = dqmStore.book1D(histogramName_strip_numPhotons.Data(), histogramName_strip_numPhotons.Data(), 25, -0.5, 24.5);
      histogram_strip_numPhotons_ = me_strip_numPhotons_->getTH1();
      assert(histogram_strip_numPhotons_);
    }
    void fillHistograms(const l1t::TallinnL1PFTauCollection& l1PFTaus, double evtWeight)
    {
      for ( l1t::TallinnL1PFTauCollection::const_iterator l1PFTau = l1PFTaus.begin(); 
	    l1PFTau != l1PFTaus.end(); ++l1PFTau ) {	
	double l1PFTau_absEta = TMath::Abs(l1PFTau->eta());
	if ( (min_pt_     < 0. || l1PFTau->pt()  > min_absEta_) &&
	     (max_pt_     < 0. || l1PFTau->pt()  < max_absEta_) &&
	     (min_absEta_ < 0. || l1PFTau_absEta > min_pt_)     && 
	     (max_absEta_ < 0. || l1PFTau_absEta < max_pt_)     )
	{
	  if ( l1PFTau->leadChargedPFCand().isNonnull() )
	  {
	    const l1t::PFCandidateRef& leadPFCand = l1PFTau->leadChargedPFCand();
	    if ( leadPFCand->pfTrack().isNonnull() && leadPFCand->pfTrack()->track().isNonnull() )
	    {
	      const l1t::PFTrack::L1TTTrackType& track = (*leadPFCand->pfTrack()->track());
	      fillWithOverFlow(histogram_leadTrack_pt_, track.getMomentum().perp(), evtWeight);
	      fillWithOverFlow(histogram_leadTrack_chi2_, track.getChi2(), evtWeight);
	      fillWithOverFlow(histogram_leadTrack_numStubs_, track.getStubRefs().size(), evtWeight);
	    }
	  }

	  const l1t::PFCandidateRefVector& signalPFCands = l1PFTau->signalAllL1PFCandidates();
	  for ( l1t::PFCandidateRefVector::const_iterator signalPFCand = signalPFCands.begin(); 
		signalPFCand != signalPFCands.end(); ++signalPFCand ) {
	    if ( (*signalPFCand)->pfTrack().isNonnull() && (*signalPFCand)->pfTrack()->track().isNonnull() ) 
	    {
	      const l1t::PFTrack::L1TTTrackType& track = *(*signalPFCand)->pfTrack()->track();
	      fillWithOverFlow(histogram_signalTrack_pt_, track.getMomentum().perp(), evtWeight);
	      fillWithOverFlow(histogram_signalTrack_chi2_, track.getChi2(), evtWeight);
	      fillWithOverFlow(histogram_signalTrack_numStubs_, track.getStubRefs().size(), evtWeight);
	    }
	  }

	  const l1t::PFCandidateRefVector& isolationPFCands = l1PFTau->isoAllL1PFCandidates();
	  for ( l1t::PFCandidateRefVector::const_iterator isolationPFCand = isolationPFCands.begin(); 
		isolationPFCand != isolationPFCands.end(); ++isolationPFCand ) {
	    if ( (*isolationPFCand)->pfTrack().isNonnull() && (*isolationPFCand)->pfTrack()->track().isNonnull() ) 
	    {
	      const l1t::PFTrack::L1TTTrackType& track = *(*isolationPFCand)->pfTrack()->track();
	      fillWithOverFlow(histogram_isolationTrack_pt_, track.getMomentum().perp(), evtWeight);
	      fillWithOverFlow(histogram_isolationTrack_chi2_, track.getChi2(), evtWeight);
	      fillWithOverFlow(histogram_isolationTrack_numStubs_, track.getStubRefs().size(), evtWeight);
	    }
	  }

	  reco::Candidate::LorentzVector tracksP4;
	  for ( l1t::PFCandidateRefVector::const_iterator signalPFCand = signalPFCands.begin(); 
		signalPFCand != signalPFCands.end(); ++signalPFCand ) {
	    if ( (*signalPFCand)->pfTrack().isNonnull() && (*signalPFCand)->pfTrack()->track().isNonnull() ) 
	    {
	      tracksP4 += getChargedPionP4((*signalPFCand)->p4());
	    }
	  }
	  fillWithOverFlow(histogram_mTracks_, tracksP4.mass(), evtWeight);
	  reco::Candidate::LorentzVector stripP4;
	  int strip_numPhotons   = 0;
	  int strip_numElectrons = 0;
	  const l1t::PFCandidateRefVector& stripPFCands = l1PFTau->stripAllL1PFCandidates();
	  for ( l1t::PFCandidateRefVector::const_iterator stripPFCand = stripPFCands.begin(); 
		stripPFCand != stripPFCands.end(); ++stripPFCand ) {
	    stripP4 += (*stripPFCand)->p4();
	    if      ( (*stripPFCand)->id() == l1t::PFCandidate::Photon   ) ++strip_numPhotons;
	    else if ( (*stripPFCand)->id() == l1t::PFCandidate::Electron ) ++strip_numElectrons;
	  }
	  stripP4 = getNeutralPionP4(stripP4);
	  reco::Candidate::LorentzVector tracks_plus_stripP4 = tracksP4 + stripP4;
	  fillWithOverFlow(histogram_mTracks_plus_Strip_, tracks_plus_stripP4.mass(), evtWeight);
	  fillWithOverFlow(histogram_strip_pt_, stripP4.pt(), evtWeight);
	  fillWithOverFlow(histogram_strip_numPhotons_, strip_numPhotons, evtWeight);
	  fillWithOverFlow(histogram_strip_numElectrons_, strip_numElectrons, evtWeight);
	}
      }
    }
    MonitorElement* me_leadTrack_pt_;
    TH1* histogram_leadTrack_pt_;
    MonitorElement* me_leadTrack_chi2_;
    TH1* histogram_leadTrack_chi2_;
    MonitorElement* me_leadTrack_numStubs_;
    TH1* histogram_leadTrack_numStubs_;
    MonitorElement* me_signalTrack_pt_;
    TH1* histogram_signalTrack_pt_;
    MonitorElement* me_signalTrack_chi2_;
    TH1* histogram_signalTrack_chi2_;
    MonitorElement* me_signalTrack_numStubs_;
    TH1* histogram_signalTrack_numStubs_;
    MonitorElement* me_isolationTrack_pt_;
    TH1* histogram_isolationTrack_pt_;
    MonitorElement* me_isolationTrack_chi2_;
    TH1* histogram_isolationTrack_chi2_;
    MonitorElement* me_isolationTrack_numStubs_;
    TH1* histogram_isolationTrack_numStubs_;
    MonitorElement* me_mTracks_;
    TH1* histogram_mTracks_;
    MonitorElement* me_mTracks_plus_Strip_;
    TH1* histogram_mTracks_plus_Strip_;
    MonitorElement* me_strip_pt_;
    TH1* histogram_strip_pt_;
    MonitorElement* me_strip_numPhotons_;
    TH1* histogram_strip_numPhotons_;
    MonitorElement* me_strip_numElectrons_;
    TH1* histogram_strip_numElectrons_;
    // cuts applied to offline and generator-level taus in denominator 
    double min_pt_;
    double max_pt_;
    double min_absEta_;
    double max_absEta_;
  };
  std::vector<plotEntryType*> plots_;
};

#endif   
