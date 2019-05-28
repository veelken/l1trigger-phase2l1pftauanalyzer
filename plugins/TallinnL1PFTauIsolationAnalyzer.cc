#include "L1Trigger/TallinnL1PFTauAnalyzer/plugins/TallinnL1PFTauIsolationAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "L1Trigger/TallinnL1PFTaus/interface/lutAuxFunctions.h"

#include <iostream>
#include <iomanip>

TallinnL1PFTauIsolationAnalyzer::TallinnL1PFTauIsolationAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , dRmatch_(cfg.getParameter<double>("dRmatch"))
  , inputFileName_rhoCorr_(cfg.getParameter<std::string>("inputFileName_rhoCorr"))
  , inputFile_rhoCorr_(nullptr)
  , histogramName_rhoCorr_(cfg.getParameter<std::string>("histogramName_rhoCorr"))
  , histogram_rhoCorr_(nullptr)
  , histogram_rhoCorr_yMax_(-1.)
{
  src_l1Taus_ = cfg.getParameter<edm::InputTag>("srcL1Taus");
  token_l1Taus_ = consumes<l1t::TallinnL1PFTauCollection>(src_l1Taus_);
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

TallinnL1PFTauIsolationAnalyzer::~TallinnL1PFTauIsolationAnalyzer()
{
  delete histogram_rhoCorr_;

  for ( auto isolationPlot : isolationPlots_ ) 
  {
    delete isolationPlot;
  }
}

void TallinnL1PFTauIsolationAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("TallinnL1PFTauIsolationAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
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
  std::vector<double> ptThresholds = { 20., 30., 40. };
  for ( auto decayMode : decayModes )
  {
    for ( auto ptThreshold : ptThresholds )
    {
      isolationPlots_.push_back(new isolationPlotEntryType(ptThreshold, -1., -1., 1.0, decayMode)); 
      isolationPlots_.push_back(new isolationPlotEntryType(ptThreshold, -1., -1., 1.4, decayMode));
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

void TallinnL1PFTauIsolationAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<l1t::TallinnL1PFTauCollection> l1Taus;
  evt.getByToken(token_l1Taus_, l1Taus);
  
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
  
  for ( l1t::TallinnL1PFTauCollection::const_iterator l1Tau = l1Taus->begin(); l1Tau != l1Taus->end();  ++l1Tau )
  {
    std::string genTau_decayMode;
    double dRmin = 1.e+3;
    bool isMatched = false;
    if ( src_genTaus_.label() != "" ) 
    {
      for ( reco::GenJetCollection::const_iterator genTau = genTaus->begin(); genTau != genTaus->end(); ++genTau )
      {
	double dR = reco::deltaR(genTau->eta(), genTau->phi(), l1Tau->eta(), l1Tau->phi());
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
        int idxBin = histogram_rhoCorr_->FindBin(TMath::Abs(l1Tau->eta()));
        if ( !(idxBin >= 1 && idxBin <= histogram_rhoCorr_->GetNbinsX()) )
        {
	  std::cerr << "Warning in <TallinnL1PFTauIsolationAnalyzer::analyze>:" 
            	    << " Failed to compute rho correction for abs(eta) = " << l1Tau->eta() << " --> skipping event !!" << std::endl;
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
        isolationPlot->fillHistograms_wGenMatching(*l1Tau, rhoCorr, isMatched, genTau_decayMode, evtWeight);
      }
      else
      {
        isolationPlot->fillHistograms_woGenMatching(*l1Tau, rhoCorr, evtWeight);
      }
    }  
  }
}

void TallinnL1PFTauIsolationAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TallinnL1PFTauIsolationAnalyzer);
