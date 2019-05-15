#ifndef L1Trigger_TallinnL1PFTauAnalyzer_L1VertexAnalyzer_h
#define L1Trigger_TallinnL1PFTauAnalyzer_L1VertexAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/L1TVertex/interface/Vertex.h"                                  // l1t::Vertex, l1t::VertexCollection
#include "DataFormats/Phase2L1ParticleFlow/interface/PFTrack.h"                      // l1t::PFTrack::L1TTTrackType
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"                  // l1t::PFCandidate
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidateFwd.h"               // l1t::PFCandidateCollection
#include "DataFormats/VertexReco/interface/Vertex.h"                                 // reco::Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"                              // reco::VertexCollection
#include "DataFormats/TrackReco/interface/Track.h"                                   // reco::Track
#include "DataFormats/TrackReco/interface/TrackFwd.h"                                // reco::TrackCollection
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"                 // reco::PFCandidate
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"              // reco::PFCandidateCollection  
#include "DataFormats/JetReco/interface/GenJet.h"                                    // reco::GenJet
#include "DataFormats/JetReco/interface/GenJetCollection.h"                          // reco::GenJetCollection
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"                        // reco::GenParticle
#include "DataFormats/Math/interface/deltaR.h"                                       // deltaR

#include "L1Trigger/TallinnL1PFTauAnalyzer/interface/GenChargedHadronToTrackMatch.h" // GenChargedHadronToOfflineTrackMatch, GenChargedHadronToL1TrackMatch

#include <TH1.h>
#include <TString.h> // TString, Form()
#include <TMath.h> // TMath::Abs()

#include <vector>
#include <string>

namespace 
{
  // auxiliary struct for generator-level charged hadrons produced in tau decays
  class GenChargedHadron_and_genTau_decayMode
  {
   public:
    GenChargedHadron_and_genTau_decayMode(const reco::GenParticle* genChargedHadron, const std::string& genTau_decayMode)
      : genChargedHadron_(genChargedHadron)
      , genTau_decayMode_(genTau_decayMode)
    {
      assert(genChargedHadron_);
    }
    ~GenChargedHadron_and_genTau_decayMode()
    {}

    const reco::GenParticle* genChargedHadron()     const { return genChargedHadron_;        }          
    double                   genChargedHadron_pt()  const { return genChargedHadron_->pt();  }
    double                   genChargedHadron_eta() const { return genChargedHadron_->eta(); }
    double                   genChargedHadron_phi() const { return genChargedHadron_->phi(); }
    const std::string&       genTau_decayMode()     const { return genTau_decayMode_;        }

   private:
    const reco::GenParticle* genChargedHadron_;
    std::string genTau_decayMode_;
  };

  // auxiliary struct for matches between generator-level charged hadrons produced in tau decays and tracks reconstructed offline or on L1 trigger level
  template <class T1, class T2>
  class GenChargedHadronToTrackMatch_and_genTau_decayMode
  {
   public:
    GenChargedHadronToTrackMatch_and_genTau_decayMode(const T2& genChargedHadronToTrackMatch, const std::string& genTau_decayMode, double dR)
      : genChargedHadronToTrackMatch_(genChargedHadronToTrackMatch)
      , genTau_decayMode_(genTau_decayMode)
      , dR_(dR)
    {}
    ~GenChargedHadronToTrackMatch_and_genTau_decayMode()
    {}

    const reco::Candidate* genChargedHadron()             const { return genChargedHadronToTrackMatch_.genChargedHadron();     }
    bool                   hasGenChargedHadron()          const { return genChargedHadronToTrackMatch_.hasGenChargedHadron();  }
    double                 genChargedHadron_pt()          const { return genChargedHadronToTrackMatch_.genChargedHadron_pt();  }
    double                 genChargedHadron_eta()         const { return genChargedHadronToTrackMatch_.genChargedHadron_eta(); }
    double                 genChargedHadron_phi()         const { return genChargedHadronToTrackMatch_.genChargedHadron_phi(); }
    const T1*              recTrack()                     const { return genChargedHadronToTrackMatch_.recTrack();             }
    bool                   hasRecTrack()                  const { return genChargedHadronToTrackMatch_.hasRecTrack();          }
    double                 recTrack_pt()                  const { return genChargedHadronToTrackMatch_.recTrack_pt();          }
    double                 recTrack_eta()                 const { return genChargedHadronToTrackMatch_.recTrack_eta();         }
    double                 recTrack_phi()                 const { return genChargedHadronToTrackMatch_.recTrack_phi();         }
    const std::string&     genTau_decayMode()             const { return genTau_decayMode_;                                    }
    double                 dR()                           const { return dR_;                                                  }
    const T2&              genChargedHadronToTrackMatch() const { return genChargedHadronToTrackMatch_;                        }
    
  private:
    T2 genChargedHadronToTrackMatch_;
    std::string genTau_decayMode_;
    double dR_; // used to sort (rank) matches
  };
  
  typedef GenChargedHadronToTrackMatch_and_genTau_decayMode<reco::Track, GenChargedHadronToOfflineTrackMatch> GenChargedHadronToOfflineTrackMatch_and_genTau_decayMode;
  typedef GenChargedHadronToTrackMatch_and_genTau_decayMode<l1t::Track, GenChargedHadronToL1TrackMatch> GenChargedHadronToL1TrackMatch_and_genTau_decayMode;
}

class L1TrackAnalyzer : public edm::EDAnalyzer 
{
 public:
  L1TrackAnalyzer(const edm::ParameterSet& cfg);
  ~L1TrackAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_genTaus_;
  edm::EDGetTokenT<reco::GenJetCollection> token_genTaus_;

  edm::InputTag src_offlineVertices_;
  edm::EDGetTokenT<reco::VertexCollection> token_offlineVertices_;

  edm::InputTag src_offlineTracks_;
  edm::EDGetTokenT<reco::TrackCollection> token_offlineTracks_;

  edm::InputTag src_offlinePFCands_;
  edm::EDGetTokenT<reco::PFCandidateCollection> token_offlinePFCands_;

  edm::InputTag src_l1Vertices_;
  edm::EDGetTokenT<l1t::VertexCollection> token_l1Vertices_;

  edm::InputTag src_l1Tracks_;
  edm::EDGetTokenT<l1t::TrackCollection> token_l1Tracks_;

  edm::InputTag src_l1PFVertex_z_;
  edm::EDGetTokenT<float> token_l1PFVertex_z_;

  edm::InputTag src_l1PFCands_;
  edm::EDGetTokenT<l1t::PFCandidateCollection> token_l1PFCands_;

  std::string dqmDirectory_;

  struct efficiencyPlotEntryType
  {
    efficiencyPlotEntryType(const std::string& recTrack_type, double genChargedHadron_min_pt, double genChargedHadron_max_absEta, const std::string& genTau_decayMode)
      : me_pt_numerator_(nullptr)
      , histogram_pt_numerator_(nullptr)
      , me_pt_denominator_(nullptr)
      , histogram_pt_denominator_(nullptr)
      , me_eta_numerator_(nullptr)
      , histogram_eta_numerator_(nullptr)
      , me_eta_denominator_(nullptr)
      , histogram_eta_denominator_(nullptr)
      , me_phi_numerator_(nullptr)
      , histogram_phi_numerator_(nullptr)
      , me_phi_denominator_(nullptr)
      , histogram_phi_denominator_(nullptr)
      , me_minDeltaR_numerator_(nullptr)
      , histogram_minDeltaR_numerator_(nullptr)
      , me_minDeltaR_denominator_(nullptr)
      , histogram_minDeltaR_denominator_(nullptr)
      , me_resolution_(nullptr)
      , histogram_resolution_(nullptr)
      , recTrack_type_(recTrack_type)
      , genChargedHadron_min_pt_(genChargedHadron_min_pt)
      , genChargedHadron_max_absEta_(genChargedHadron_max_absEta)
      , genTau_decayMode_(genTau_decayMode)
    {
      if ( recTrack_type_.length() > 0 ) 
      {	
	recTrack_type_capitalized_ = recTrack_type_;
	recTrack_type_capitalized_[0] = toupper(recTrack_type_capitalized_[0]);	
      }
      if ( genTau_decayMode_.length() > 0 ) 
      {	
	genTau_decayMode_capitalized_ = genTau_decayMode_;
	genTau_decayMode_capitalized_[0] = toupper(genTau_decayMode_capitalized_[0]);	
      }
    }
    ~efficiencyPlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix;
      if ( genChargedHadron_max_absEta_ > 0.     ) histogramName_suffix.Append(Form("_absEtaLt%1.2f", genChargedHadron_max_absEta_));
      if ( genChargedHadron_min_pt_     > 0.     ) histogramName_suffix.Append(Form("_ptGt%1.0f",     genChargedHadron_min_pt_));
      if ( genTau_decayMode_            != "all" ) histogramName_suffix.Append(Form("_gen%sTau",      genTau_decayMode_capitalized_.data()));      
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      const int numBins_pt = 19;
      float binning_pt[numBins_pt + 1] = { 
        0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12.5, 15., 17.5, 20., 25., 30., 40., 60., 100.
      };

      TString histogramName_pt_numerator = Form("eff%s_vs_pt_numerator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_pt_numerator_ = dqmStore.book1D(histogramName_pt_numerator.Data(), histogramName_pt_numerator.Data(), numBins_pt, binning_pt);
      histogram_pt_numerator_ = me_pt_numerator_->getTH1();
      assert(histogram_pt_numerator_);
      TString histogramName_pt_denominator = Form("eff%s_vs_pt_denominator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_pt_denominator_ = dqmStore.book1D(histogramName_pt_denominator.Data(), histogramName_pt_denominator.Data(), numBins_pt, binning_pt);
      histogram_pt_denominator_ = me_pt_denominator_->getTH1();
      assert(histogram_pt_denominator_);

      TString histogramName_eta_numerator = Form("eff%s_vs_eta_numerator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_eta_numerator_ = dqmStore.book1D(histogramName_eta_numerator.Data(), histogramName_eta_numerator.Data(), 30, -3., +3.);
      histogram_eta_numerator_ = me_eta_numerator_->getTH1();
      assert(histogram_eta_numerator_);
      TString histogramName_eta_denominator = Form("eff%s_vs_eta_denominator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_eta_denominator_ = dqmStore.book1D(histogramName_eta_denominator.Data(), histogramName_eta_denominator.Data(), 30, -3., +3.);
      histogram_eta_denominator_ = me_eta_denominator_->getTH1();
      assert(histogram_eta_denominator_);

      TString histogramName_phi_numerator = Form("eff%s_vs_phi_numerator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_phi_numerator_ = dqmStore.book1D(histogramName_phi_numerator.Data(), histogramName_phi_numerator.Data(), 18, -TMath::Pi(), +TMath::Pi());
      histogram_phi_numerator_ = me_phi_numerator_->getTH1();
      assert(histogram_phi_numerator_);
      TString histogramName_phi_denominator = Form("eff%s_vs_phi_denominator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_phi_denominator_ = dqmStore.book1D(histogramName_phi_denominator.Data(), histogramName_phi_denominator.Data(), 18, -TMath::Pi(), +TMath::Pi());
      histogram_phi_denominator_ = me_phi_denominator_->getTH1();
      assert(histogram_phi_denominator_);

      TString histogramName_minDeltaR_numerator = Form("eff%s_vs_minDeltaR_numerator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_minDeltaR_numerator_ = dqmStore.book1D(histogramName_minDeltaR_numerator.Data(), histogramName_minDeltaR_numerator.Data(), 100, 0., 0.100);
      histogram_minDeltaR_numerator_ = me_minDeltaR_numerator_->getTH1();
      assert(histogram_minDeltaR_numerator_);
      TString histogramName_minDeltaR_denominator = Form("eff%s_vs_minDeltaR_denominator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_minDeltaR_denominator_ = dqmStore.book1D(histogramName_minDeltaR_denominator.Data(), histogramName_minDeltaR_denominator.Data(), 100, 0., 0.100);
      histogram_minDeltaR_denominator_ = me_minDeltaR_denominator_->getTH1();
      assert(histogram_minDeltaR_denominator_);

      TString histogramName_resolution = Form("resolution%s%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_resolution_ = dqmStore.book1D(histogramName_resolution.Data(), histogramName_resolution.Data(), 100, -0.05, +0.05);  
      histogram_resolution_ = me_resolution_->getTH1();
      assert(histogram_resolution_);
    }
    template <class T>
    void fillHistograms(const std::vector<const T*>& cleanedGenChargedHadronToTrackMatches, double evtWeight)
    {
      double minDeltaR = 1.e+3;
      for ( auto match1 : cleanedGenChargedHadronToTrackMatches )
      {
        for ( auto match2 : cleanedGenChargedHadronToTrackMatches )
        {
	  double dR = deltaR(match1->genChargedHadron_eta(), match1->genChargedHadron_phi(), match2->genChargedHadron_eta(), match2->genChargedHadron_phi());
	  if ( dR < minDeltaR ) minDeltaR = dR;
        }
      }

      for ( auto match : cleanedGenChargedHadronToTrackMatches )
      {
	if ( genTau_decayMode_ != "all" && match->genTau_decayMode() != genTau_decayMode_ ) continue;

	if ( match->genChargedHadron_pt() > genChargedHadron_min_pt_ && TMath::Abs(match->genChargedHadron_eta()) < genChargedHadron_max_absEta_ )
	{
	  if ( match->hasRecTrack() )
	  {
	    histogram_pt_numerator_->Fill(match->genChargedHadron_pt(), evtWeight);
	    histogram_eta_numerator_->Fill(match->genChargedHadron_eta(), evtWeight);
	    histogram_phi_numerator_->Fill(match->genChargedHadron_phi(), evtWeight);
	    histogram_minDeltaR_numerator_->Fill(minDeltaR, evtWeight);
	  }

	  histogram_pt_denominator_->Fill(match->genChargedHadron_pt(), evtWeight);
	  histogram_eta_denominator_->Fill(match->genChargedHadron_eta(), evtWeight);
	  histogram_phi_denominator_->Fill(match->genChargedHadron_phi(), evtWeight);
	  histogram_minDeltaR_denominator_->Fill(minDeltaR, evtWeight);

	  if ( match->hasRecTrack() )
  	  {
	    histogram_resolution_->Fill((1./match->recTrack_pt()) - (1./match->genChargedHadron_pt()), evtWeight);
	  }
	}
      }
    }
    MonitorElement* me_pt_numerator_;
    TH1* histogram_pt_numerator_;
    MonitorElement* me_pt_denominator_;
    TH1* histogram_pt_denominator_;
    MonitorElement* me_eta_numerator_;
    TH1* histogram_eta_numerator_;
    MonitorElement* me_eta_denominator_;
    TH1* histogram_eta_denominator_;
    MonitorElement* me_phi_numerator_;
    TH1* histogram_phi_numerator_;
    MonitorElement* me_phi_denominator_;
    TH1* histogram_phi_denominator_;
    MonitorElement* me_minDeltaR_numerator_;
    TH1* histogram_minDeltaR_numerator_;
    MonitorElement* me_minDeltaR_denominator_;
    TH1* histogram_minDeltaR_denominator_;
    MonitorElement* me_resolution_;
    TH1* histogram_resolution_;
    std::string recTrack_type_;
    std::string recTrack_type_capitalized_;
    double genChargedHadron_min_pt_;
    double genChargedHadron_max_absEta_;
    std::string genTau_decayMode_;
    std::string genTau_decayMode_capitalized_;
  };
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_offlineTracks_woQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_offlineTracks_wQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_offlinePFCandTracks_woQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_offlinePFCandTracks_wQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_l1Tracks_woQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_l1Tracks_wQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_l1PFCandTracks_woQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_l1PFCandTracks_wQualityCuts_;
};

#endif   

