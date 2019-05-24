#ifndef L1Trigger_TallinnL1PFTauAnalyzer_RhoCorrAnalyzer_h
#define L1Trigger_TallinnL1PFTauAnalyzer_RhoCorrAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "L1Trigger/TallinnL1PFTaus/interface/TallinnL1PFTauQualityCut.h" // TallinnL1PFTauQualityCut
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"       // l1t::PFCandidate
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidateFwd.h"    // l1t::PFCandidateCollection

#include <TH1.h>
#include <TString.h> // TString, Form()

#include <vector>
#include <string>

class RhoCorrAnalyzer : public edm::EDAnalyzer 
{
 public:
  RhoCorrAnalyzer(const edm::ParameterSet& cfg);
  ~RhoCorrAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_rho_;
  edm::EDGetTokenT<double> token_rho_;

  edm::InputTag src_rhoNeutral_;
  edm::EDGetTokenT<double> token_rhoNeutral_;

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

