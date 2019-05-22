#ifndef L1Trigger_TallinnL1PFTauAnalyzer_RhoCorrAnalyzer_h
#define L1Trigger_TallinnL1PFTauAnalyzer_RhoCorrAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"    // l1t::PFCandidate
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidateFwd.h" // l1t::PFCandidateCollection

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

  template <class T>
  void printGenChargedHadronToTrackMatches(const std::vector<const T*>& genChargedHadronToTrackMatches)
  {
    size_t idx = 0;
    for ( auto match : genChargedHadronToTrackMatches )
    {
      std::cout << " genChargedHadron (#" << idx << "):" 
		<< " pT = " << match->genChargedHadron_pt() << ","
		<< " eta = " << match->genChargedHadron_eta() << ","
		<< " phi = " << match->genChargedHadron_phi() << ","
		<< " which is matched to recTrack:";
      if ( match->hasRecTrack() )
      {
	std::cout << " pT = " << match->recTrack_pt() << ","
		  << " eta = " << match->recTrack_eta() << ","
		  << " phi = " << match->recTrack_phi() 
		  << " (dR = " << match->dR() << ")";
      }
      else
      {
	std::cout << " N/A";
      }
      std::cout << std::endl;
      ++idx;
    }
  }

  template <class T>
  void printGenChargedHadronToTrackMatches(const std::vector<T>& genChargedHadronToTrackMatches)
  {
    std::vector<const T*> genChargedHadronToTrackMatches_ptr;
    for ( typename std::vector<T>::const_iterator match = genChargedHadronToTrackMatches.begin();
	  match != genChargedHadronToTrackMatches.end(); ++match )
    {
      genChargedHadronToTrackMatches_ptr.push_back(&(*match));
    }
    printGenChargedHadronToTrackMatches(genChargedHadronToTrackMatches_ptr);
  }
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

  edm::InputTag src_rho_;
  edm::EDGetTokenT<float> token_rho_;

  edm::InputTag src_rhoNeutral_;
  edm::EDGetTokenT<float> token_rhoNeutral_;

  edm::InputTag src_l1PFCands_;
  edm::EDGetTokenT<l1t::PFCandidateCollection> token_l1PFCands_;

  std::vector<TallinnL1PFTauQualityCut> isolationQualityCuts_;

  std::string dqmDirectory_;

  MonitorElement* me_rho_;
  TH1* histogram_rho_;
  
  MonitorElement* me_rhoNeutral_;
  TH1* histogram_rhoNeutral_;

  MonitorElement* me_neutralPFCandPt_vs_absEta_;
  TH1* histogram_neutralPFCandPt_vs_absEta_;

  MonitorElement* me_EventCounter_;
  TH1* histogram_EventCounter_;
};

#endif   

