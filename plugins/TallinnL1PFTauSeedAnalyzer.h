#ifndef L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauSeedAnalyzer_h
#define L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauSeedAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"       // l1t::PFCandidate
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidateFwd.h"    // l1t::PFCandidateCollection, l1t::PFCandidateRef
#include "DataFormats/Phase2L1ParticleFlow/interface/PFJet.h"             // l1t::PFJet
#include "DataFormats/Phase2L1ParticleFlow/interface/PFJetFwd.h"          // l1t::PFJetCollection, l1t::PFJetRef
#include "DataFormats/L1TVertex/interface/Vertex.h"                       // l1t::Vertex
#include "DataFormats/L1TVertex/interface/VertexFwd.h"                    // l1t::VertexCollection

#include "L1Trigger/TallinnL1PFTaus/interface/LocalFileInPath.h"              // LocalFileInPath
#include "L1Trigger/TallinnL1PFTaus/interface/TallinnL1PFTauQualityCut.h"     // TallinnL1PFTauQualityCut
#include "L1Trigger/TallinnL1PFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverFlow(), fillWithOverFlow2d()

#include <TFile.h>   // TFile
#include <TH1.h>     // TH1
#include <TH2.h>     // TH2
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs(), TMath::Pi()

#include <vector>
#include <string>

class TallinnL1PFTauSeedAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit TallinnL1PFTauSeedAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~TallinnL1PFTauSeedAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcL1PFCands_;
  edm::EDGetTokenT<l1t::PFCandidateCollection> tokenL1PFCands_;
  edm::InputTag srcL1PFJets_;
  edm::EDGetTokenT<l1t::PFJetCollection> tokenL1PFJets_;
  edm::InputTag srcL1Vertices_;
  edm::EDGetTokenT<l1t::VertexCollection> tokenL1Vertices_;

  std::vector<TallinnL1PFTauQualityCut> signalQualityCuts_dzCut_disabled_;

  std::string dqmDirectory_;

  double min_seedChargedPFCand_pt_;
  double max_seedChargedPFCand_eta_;
  double max_seedChargedPFCand_dz_;

  double min_seedPFJet_pt_;
  double max_seedPFJet_eta_;

  MonitorElement* me_seedChargedPFCand_pt_;
  TH1* histogram_seedChargedPFCand_pt_;
  MonitorElement* me_seedChargedPFCand_eta_;
  TH1* histogram_seedChargedPFCand_eta_;
  MonitorElement* me_seedChargedPFCand_phi_;
  TH1* histogram_seedChargedPFCand_phi_;
  MonitorElement* me_seedChargedPFCand_dz_;
  TH1* histogram_seedChargedPFCand_dz_;
  MonitorElement* me_numSeedChargedPFCands_;
  TH1* histogram_numSeedChargedPFCands_;

  MonitorElement* me_seedPFJet_pt_;
  TH1* histogram_seedPFJet_pt_;
  MonitorElement* me_seedPFJet_eta_;
  TH1* histogram_seedPFJet_eta_;
  MonitorElement* me_seedPFJet_phi_;
  TH1* histogram_seedPFJet_phi_;
  MonitorElement* me_numSeedPFJets_;
  TH1* histogram_numSeedPFJets_;
};

#endif   
