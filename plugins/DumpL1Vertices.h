#ifndef L1Trigger_TallinnL1PFTauAnalyzer_DumpL1Vertices_h
#define L1Trigger_TallinnL1PFTauAnalyzer_DumpL1Vertices_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/L1TCorrelator/interface/TkPrimaryVertex.h" // l1t::TkPrimaryVertex, l1t::TkPrimaryVertexCollection

#include <vector>
#include <string>

class DumpL1Vertices : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpL1Vertices(const edm::ParameterSet&);
    
  // destructor
  ~DumpL1Vertices();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<l1t::TkPrimaryVertexCollection> token_;
};

#endif   
