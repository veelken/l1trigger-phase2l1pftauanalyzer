#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/L1HPSPFTauIsolationAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/lutAuxFunctions.h" // openFile(), loadTH1()

#include <iostream>
#include <iomanip>

L1HPSPFTauIsolationAnalyzer::L1HPSPFTauIsolationAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , dRmatch_(cfg.getParameter<double>("dRmatch"))
  , inputFileName_rhoCorr_(cfg.getParameter<std::string>("inputFileName_rhoCorr"))
  , inputFile_rhoCorr_(nullptr)
  , histogramName_rhoCorr_(cfg.getParameter<std::string>("histogramName_rhoCorr"))
  , histogram_rhoCorr_(nullptr)
  , histogram_rhoCorr_yMax_(-1.)
{
  src_l1PFTaus_ = cfg.getParameter<edm::InputTag>("srcL1PFTaus");
  token_l1PFTaus_ = consumes<l1t::L1HPSPFTauCollection>(src_l1PFTaus_);
  src_genTaus_ = cfg.getParameter<edm::InputTag>("srcGenTaus");
  token_genTaus_ = consumes<reco::GenJetCollection>(src_genTaus_);
  src_rho_ = cfg.getParameter<edm::InputTag>("srcRho");
  token_rho_ = consumes<double>(src_rho_);

  if ( inputFileName_rhoCorr_ != "" && histogramName_rhoCorr_ != "" )
  {
    LocalFileInPath inputFileName_rhoCorr_fip(inputFileName_rhoCorr_);
    inputFile_rhoCorr_ = openFile(inputFileName_rhoCorr_fip);
    TH1* histogram_rhoCorr_temp = loadTH1(inputFile_rhoCorr_, histogramName_rhoCorr_);
    std::string histogramName_rhoCorr = Form("%s_cloned", histogram_rhoCorr_temp->GetName());
    histogram_rhoCorr_ = (TH1*)histogram_rhoCorr_temp->Clone(histogramName_rhoCorr.data());
    int numBins = histogram_rhoCorr_->GetNbinsX();
    for ( int idxBin = 1; idxBin <= numBins; ++idxBin )
    {
      double binContent = histogram_rhoCorr_->GetBinContent(idxBin);
      if ( binContent > histogram_rhoCorr_yMax_ ) 
      {
	histogram_rhoCorr_yMax_ = binContent;
      }
    }
    delete inputFile_rhoCorr_;
    inputFile_rhoCorr_ = nullptr;
  }

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

L1HPSPFTauIsolationAnalyzer::~L1HPSPFTauIsolationAnalyzer()
{
  delete histogram_rhoCorr_;

  for ( auto isolationPlot : isolationPlots_ ) 
  {
    delete isolationPlot;
  }
}

void L1HPSPFTauIsolationAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("L1HPSPFTauIsolationAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  std::vector<std::string> decayModes;
  if ( src_genTaus_.label() != "" ) 
  {
    decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  } 
  else
  {
    decayModes = { "all" };
  }
  std::vector<double> ptThresholds = { 20., 25., 30., 35., 40., 45., 50. };
  for ( auto decayMode : decayModes )
  {
    for ( auto ptThreshold : ptThresholds )
    {
      isolationPlots_.push_back(new isolationPlotEntryType(ptThreshold, -1., -1.,  1.4,   decayMode));
      isolationPlots_.push_back(new isolationPlotEntryType(ptThreshold, -1.,  1.4, 2.172, decayMode));
      isolationPlots_.push_back(new isolationPlotEntryType(ptThreshold, -1., -1.,  2.172, decayMode));
      isolationPlots_.push_back(new isolationPlotEntryType(ptThreshold, -1., -1.,  2.4,   decayMode));
    }
  }

  const int numBins_absEta = 10;
  float binning_absEta[numBins_absEta + 1] = { 
    0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0 
  };
  const int numBins_pt = 4;
  float binning_pt[numBins_pt + 1] = { 
    20., 30., 40., 60., 100.
  };
  for ( size_t idxBin_absEta = 0; idxBin_absEta < numBins_absEta; ++idxBin_absEta )
  {
    float min_absEta = binning_absEta[idxBin_absEta];
    float max_absEta = binning_absEta[idxBin_absEta + 1];
    for ( size_t idxBin_pt = 0; idxBin_pt < numBins_pt; ++idxBin_pt )
    {  
      float min_pt = binning_pt[idxBin_pt];
      float max_pt = binning_pt[idxBin_pt + 1];
      isolationPlots_.push_back(new isolationPlotEntryType(min_pt, max_pt, min_absEta, max_absEta, "all")); 
    }
  }

  for ( auto isolationPlot : isolationPlots_ ) 
  {
    isolationPlot->bookHistograms(dqmStore);
  }
}

void L1HPSPFTauIsolationAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<l1t::L1HPSPFTauCollection> l1PFTaus;
  evt.getByToken(token_l1PFTaus_, l1PFTaus);
  
  edm::Handle<reco::GenJetCollection> genTaus;
  if ( src_genTaus_.label() != "" ) 
  {
    evt.getByToken(token_genTaus_, genTaus);
  }

  edm::Handle<double> rho;
  if ( src_rho_.label() != "" ) 
  {
    evt.getByToken(token_rho_, rho);
  }
  
  const double evtWeight = 1.;
  
  for ( l1t::L1HPSPFTauCollection::const_iterator l1PFTau = l1PFTaus->begin(); l1PFTau != l1PFTaus->end();  ++l1PFTau )
  {
    std::string genTau_decayMode;
    double dRmin = 1.e+3;
    bool isMatched = false;
    if ( src_genTaus_.label() != "" ) 
    {
      for ( reco::GenJetCollection::const_iterator genTau = genTaus->begin(); genTau != genTaus->end(); ++genTau )
      {
	double dR = reco::deltaR(genTau->eta(), genTau->phi(), l1PFTau->eta(), l1PFTau->phi());
	if ( dR < dRmatch_ && dR < dRmin ) 
	{ 
	  genTau_decayMode = JetMCTagUtils::genTauDecayMode(*genTau);
	  dRmin = dR;
	  isMatched = true;
	}
      }
    }

    double rhoCorr = 0.;
    if ( src_rho_.label() != "" ) 
    {
      const double isolationConeArea = TMath::Pi()*0.4*0.4;
      rhoCorr = isolationConeArea*(*rho);
      if ( histogram_rhoCorr_ && histogram_rhoCorr_yMax_ > 0. ) 
      {
        int idxBin = histogram_rhoCorr_->FindBin(TMath::Abs(l1PFTau->eta()));
        if ( !(idxBin >= 1 && idxBin <= histogram_rhoCorr_->GetNbinsX()) )
        {
	  std::cerr << "Warning in <L1HPSPFTauIsolationAnalyzer::analyze>:" 
            	    << " Failed to compute rho correction for abs(eta) = " << l1PFTau->eta() << " --> skipping event !!" << std::endl;
          return;
        }
        rhoCorr *= histogram_rhoCorr_->GetBinContent(idxBin)/histogram_rhoCorr_yMax_;
      }
      //std::cout << "rho = " << (*rho) << ": rhoCorr = " << rhoCorr << std::endl;
    }
    
    for ( auto isolationPlot : isolationPlots_ ) 
    {    
      if ( src_genTaus_.label() != "" ) 
      {
        isolationPlot->fillHistograms_wGenMatching(*l1PFTau, rhoCorr, isMatched, genTau_decayMode, evtWeight);
      }
      else
      {
        isolationPlot->fillHistograms_woGenMatching(*l1PFTau, rhoCorr, evtWeight);
      }
    }  
  }
}

void L1HPSPFTauIsolationAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(L1HPSPFTauIsolationAnalyzer);
