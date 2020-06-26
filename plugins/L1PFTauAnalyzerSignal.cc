#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/L1PFTauAnalyzerSignal.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

L1PFTauAnalyzerSignal::L1PFTauAnalyzerSignal(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  srcNumerator_ = cfg.getParameter<edm::InputTag>("srcNumerator");
  tokenNumerator_ = consumes<l1t::L1PFTauCollection>(srcNumerator_);
  srcDenominator_ = cfg.getParameter<edm::InputTag>("srcDenominator");
  std::string typeDenominator_string = cfg.getParameter<std::string>("typeDenominator");
  if ( typeDenominator_string == "gen" ) 
  {
    typeDenominator_ = kGen;
    tokenDenominator_gen_ = consumes<reco::GenJetCollection>(srcDenominator_);
  }
  else if ( typeDenominator_string == "offline" ) 
  {
    typeDenominator_ = kOffline;
    tokenDenominator_offline_ = consumes<pat::TauCollection>(srcDenominator_);
  }
  else
  {
    throw cms::Exception("L1PFTauAnalyzerSignal") 
      << " Invalid Configuration parameter 'typeDenominator' = " << typeDenominator_string << " !!\n";;
  }

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

L1PFTauAnalyzerSignal::~L1PFTauAnalyzerSignal()
{
  for ( auto efficiencyPlot : efficiencyPlots_ ) 
  {
    delete efficiencyPlot;
  }
}

void L1PFTauAnalyzerSignal::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("L1PFTauAnalyzerSignal") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());

  std::vector<std::string> decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  std::vector<double> min_absEtaValues = { -1.,   1.4,   1.4, -1.,    -1.  };
  std::vector<double> max_absEtaValues = {  1.4,  2.172, 2.4,  2.172,  2.4 };
  assert(min_absEtaValues.size() == max_absEtaValues.size());
  std::vector<double> ptThresholds = { 20., 25., 30., 35., 40., 45., 50. };
    size_t numAbsEtaRanges = min_absEtaValues.size();
  for ( size_t idxAbsEtaRange = 0; idxAbsEtaRange < numAbsEtaRanges; ++idxAbsEtaRange )
  {
    double min_absEta = min_absEtaValues[idxAbsEtaRange];
    double max_absEta = max_absEtaValues[idxAbsEtaRange];
    for ( auto decayMode : decayModes )
    {
      for ( auto ptThreshold : ptThresholds )
      {
	TString dqmDirectory = dqmDirectory_.data();
	if ( min_absEta >= 0. && max_absEta > 0. ) 
	{ 
	  dqmDirectory.Append(Form("/absEta%1.2fto%1.2f", min_absEta, max_absEta));
        }
        else if ( min_absEta >= 0. ) 
        {
          dqmDirectory.Append(Form("/absEtaGt%1.2f", min_absEta));
        }
        else if ( max_absEta > 0. ) 
        {
          dqmDirectory.Append(Form("/absEtaLt%1.2f", max_absEta));
        }
        std::string decayMode_capitalized = decayMode;
        decayMode_capitalized[0] = toupper(decayMode_capitalized[0]);	
        dqmDirectory.Append(Form("/gen%sTau", decayMode_capitalized.data()));
        dqmDirectory = dqmDirectory.ReplaceAll(".", "p");

	dqmStore.setCurrentFolder(dqmDirectory.Data());
        efficiencyPlotEntryType* efficiencyPlots_vLoose  = new efficiencyPlotEntryType(45., 1.e+3, min_absEta, max_absEta, decayMode, ptThreshold, "vLoose"); // vLoose
	efficiencyPlots_vLoose->bookHistograms(dqmStore);
	efficiencyPlots_.push_back(efficiencyPlots_vLoose);
	efficiencyPlotEntryType* efficiencyPlots_Loose   = new efficiencyPlotEntryType(45., 1.e+3, min_absEta, max_absEta, decayMode, ptThreshold, "Loose");  // Loose
	efficiencyPlots_Loose->bookHistograms(dqmStore);
	efficiencyPlots_.push_back(efficiencyPlots_Loose);
	efficiencyPlotEntryType* efficiencyPlots_Medium  = new efficiencyPlotEntryType(45., 1.e+3, min_absEta, max_absEta, decayMode, ptThreshold, "Medium"); // Medium
	efficiencyPlots_Medium->bookHistograms(dqmStore);
	efficiencyPlots_.push_back(efficiencyPlots_Medium);
	efficiencyPlotEntryType* efficiencyPlots_Tight   = new efficiencyPlotEntryType(45., 1.e+3, min_absEta, max_absEta, decayMode, ptThreshold, "Tight");  // Tight
	efficiencyPlots_Tight->bookHistograms(dqmStore);
	efficiencyPlots_.push_back(efficiencyPlots_Tight);
      }
    }
  }
}

void L1PFTauAnalyzerSignal::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<l1t::L1PFTauCollection> numeratorTaus;
  evt.getByToken(tokenNumerator_, numeratorTaus);
  
  const double evtWeight = 1.;

  if ( typeDenominator_ == kGen )
  {
    edm::Handle<reco::GenJetCollection> denominatorTaus_gen;
    evt.getByToken(tokenDenominator_gen_, denominatorTaus_gen);

    for ( auto efficiencyPlot : efficiencyPlots_ ) 
    {    
      efficiencyPlot->fillHistograms(*numeratorTaus, *denominatorTaus_gen, evtWeight);
    }
  }

  if ( typeDenominator_ == kOffline )
  {
    edm::Handle<pat::TauCollection> denominatorTaus_offline;
    evt.getByToken(tokenDenominator_offline_, denominatorTaus_offline);

    for ( auto efficiencyPlot : efficiencyPlots_ ) 
    {    
      efficiencyPlot->fillHistograms(*numeratorTaus, *denominatorTaus_offline, evtWeight);
    }
  }
}

void L1PFTauAnalyzerSignal::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(L1PFTauAnalyzerSignal);
