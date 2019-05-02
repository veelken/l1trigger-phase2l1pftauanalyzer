#ifndef L1Trigger_TallinnL1PFTauAnalyzer_DQMSimpleFileSaver_h
#define L1Trigger_TallinnL1PFTauAnalyzer_DQMSimpleFileSaver_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "L1Trigger/TallinnL1PFTauAnalyzer/interface/dqmAuxFunctions.h"

#include <string>
#include <vector>

class DQMSimpleFileSaver : public edm::EDAnalyzer
{
 public:
  explicit DQMSimpleFileSaver(const edm::ParameterSet&);
  virtual ~DQMSimpleFileSaver();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();  

private:
  std::vector<outputCommandEntry>* outputCommands_;

  std::string outputFileName_;

  int cfgError_;
};

#endif


