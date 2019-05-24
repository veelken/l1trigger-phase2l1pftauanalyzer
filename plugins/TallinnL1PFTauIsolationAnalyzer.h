#ifndef L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauIsolationAnalyzer_h
#define L1Trigger_TallinnL1PFTauAnalyzer_TallinnL1PFTauIsolationAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTau.h"         // l1t::TallinnL1PFTau
#include "DataFormats/TallinnL1PFTaus/interface/TallinnL1PFTauFwd.h"      // l1t::TallinnL1PFTauCollection
#include "DataFormats/Math/interface/deltaR.h"                            // reco::deltaR
#include "DataFormats/JetReco/interface/GenJet.h"                         // reco::GenJet
#include "DataFormats/JetReco/interface/GenJetCollection.h"               // reco::GenJetCollection
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"                   // JetMCTagUtils::genTauDecayMode()

#include "L1Trigger/TallinnL1PFTaus/interface/LocalFileInPath.h"          // LocalFileInPath
#include "L1Trigger/TallinnL1PFTaus/interface/TallinnL1PFTauQualityCut.h" // TallinnL1PFTauQualityCut

#include <TFile.h>   // TFile
#include <TH1.h>     // TH1, TH2
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs()

#include <vector>
#include <string>

class TallinnL1PFTauIsolationAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit TallinnL1PFTauIsolationAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~TallinnL1PFTauIsolationAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_l1Taus_;
  edm::EDGetTokenT<l1t::TallinnL1PFTauCollection> token_l1Taus_;
  edm::InputTag src_genTaus_;
  edm::EDGetTokenT<reco::GenJetCollection> token_genTaus_;
  double dRmatch_;
  edm::InputTag src_rho_;
  edm::EDGetTokenT<double> token_rho_;

  std::string inputFileName_rhoCorr_;
  TFile* inputFile_rhoCorr_;
  std::string histogramName_rhoCorr_;
  TH1* histogram_rhoCorr_;

  std::string dqmDirectory_;

  struct isolationPlotEntryType
  {
    isolationPlotEntryType(double min_pt, double max_pt, double min_absEta, double max_absEta, const std::string& decayMode)
      : me_absChargedIso_(nullptr)
      , histogram_absChargedIso_(nullptr)
      , me_absNeutralIso_(nullptr)
      , histogram_absNeutralIso_(nullptr)
      , me_absCombinedIso_(nullptr)
      , histogram_absCombinedIso_(nullptr)
      , me_relChargedIso_(nullptr)
      , histogram_relChargedIso_(nullptr)
      , me_relNeutralIso_(nullptr)
      , histogram_relNeutralIso_(nullptr)
      , me_relCombinedIso_(nullptr)
      , histogram_relCombinedIso_(nullptr)
      , min_pt_(min_pt)
      , max_pt_(max_pt)
      , min_absEta_(min_absEta)
      , max_absEta_(max_absEta)
      , decayMode_(decayMode)
      , dRmatch_(0.3)
    {}
    ~isolationPlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix = decayMode_.data();      
      if      ( min_absEta_ >= 0. && max_absEta_ > min_absEta_ ) histogramName_suffix.Append(Form("_absEta%1.2fto%1.2f", min_absEta_, max_absEta_));
      else if ( min_absEta_ >= 0.                              ) histogramName_suffix.Append(Form("_absEtaGt%1.2f",      min_absEta_             ));
      else if (                      max_absEta_ >          0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f",                   max_absEta_));
      else throw cms::Exception("isolationPlotEntryType") 
	     << " Invalid Configuration parameters min_absEta = " << min_absEta_ << " and max_absEta = " << max_absEta_ << "!!\n";
      if      ( min_pt_     >  0. && max_pt_     > min_pt_     ) histogramName_suffix.Append(Form("_pt%1.2fto%1.2f",     min_pt_,     max_pt_    ));
      else if ( min_pt_     >  0.                              ) histogramName_suffix.Append(Form("_ptGt%1.2f",          min_pt_                 ));
      else if (                      max_pt_     >          0. ) histogramName_suffix.Append(Form("_ptLt%1.2f",                       max_pt_    ));
      else throw cms::Exception("isolationPlotEntryType") 
	     << " Invalid Configuration parameters min_pt = "     << min_pt_     << " and max_pt = "     << max_pt_     << "!!\n";      
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_absChargedIso = Form("absChargedIso_%s", histogramName_suffix.Data());
      me_absChargedIso_ = dqmStore.book1D(histogramName_absChargedIso.Data(), histogramName_absChargedIso.Data(), 250, 0., 25.);
      histogram_absChargedIso_ = me_absChargedIso_->getTH1();
      assert(histogram_absChargedIso_);

      TString histogramName_absNeutralIso = Form("absNeutralIso_%s", histogramName_suffix.Data());
      me_absNeutralIso_ = dqmStore.book1D(histogramName_absNeutralIso.Data(), histogramName_absNeutralIso.Data(), 250, 0., 25.);
      histogram_absNeutralIso_ = me_absNeutralIso_->getTH1();
      assert(histogram_absNeutralIso_);

      TString histogramName_absCombinedIso = Form("absCombinedIso_%s", histogramName_suffix.Data());
      me_absCombinedIso_ = dqmStore.book1D(histogramName_absCombinedIso.Data(), histogramName_absCombinedIso.Data(), 250, 0., 25.);
      histogram_absCombinedIso_ = me_absCombinedIso_->getTH1();
      assert(histogram_absCombinedIso_);

      TString histogramName_absCombinedIso_wDeltaBetaCorr = Form("absCombinedIso_wDeltaBetaCorr_%s", histogramName_suffix.Data());
      me_absCombinedIso_wDeltaBetaCorr_ = dqmStore.book1D(histogramName_absCombinedIso_wDeltaBetaCorr.Data(), histogramName_absCombinedIso_wDeltaBetaCorr.Data(), 250, 0., 25.);
      histogram_absCombinedIso_wDeltaBetaCorr_ = me_absCombinedIso_wDeltaBetaCorr_->getTH1();
      assert(histogram_absCombinedIso_wDeltaBetaCorr_);

      TString histogramName_absCombinedIso_wRhoCorr = Form("absCombinedIso_wRhoCorr_%s", histogramName_suffix.Data());
      me_absCombinedIso_wRhoCorr_ = dqmStore.book1D(histogramName_absCombinedIso_wRhoCorr.Data(), histogramName_absCombinedIso_wRhoCorr.Data(), 250, 0., 25.);
      histogram_absCombinedIso_wRhoCorr_ = me_absCombinedIso_wRhoCorr_->getTH1();
      assert(histogram_absCombinedIso_wRhoCorr_);

      TString histogramName_relChargedIso = Form("relChargedIso_%s", histogramName_suffix.Data());
      me_relChargedIso_ = dqmStore.book1D(histogramName_relChargedIso.Data(), histogramName_relChargedIso.Data(), 100, 0., 1.);
      histogram_relChargedIso_ = me_relChargedIso_->getTH1();
      assert(histogram_relChargedIso_);

      TString histogramName_relNeutralIso = Form("relNeutralIso_%s", histogramName_suffix.Data());
      me_relNeutralIso_ = dqmStore.book1D(histogramName_relNeutralIso.Data(), histogramName_relNeutralIso.Data(), 100, 0., 1.);
      histogram_relNeutralIso_ = me_relNeutralIso_->getTH1();
      assert(histogram_relNeutralIso_);

      TString histogramName_relCombinedIso = Form("relCombinedIso_%s", histogramName_suffix.Data());
      me_relCombinedIso_ = dqmStore.book1D(histogramName_relCombinedIso.Data(), histogramName_relCombinedIso.Data(), 100, 0., 1.);
      histogram_relCombinedIso_ = me_relCombinedIso_->getTH1();
      assert(histogram_relCombinedIso_);

      TString histogramName_relCombinedIso_wDeltaBetaCorr = Form("relCombinedIso_wDeltaBetaCorr_%s", histogramName_suffix.Data());
      me_relCombinedIso_wDeltaBetaCorr_ = dqmStore.book1D(histogramName_relCombinedIso_wDeltaBetaCorr.Data(), histogramName_relCombinedIso_wDeltaBetaCorr.Data(), 100, 0., 1.);
      histogram_relCombinedIso_wDeltaBetaCorr_ = me_relCombinedIso_wDeltaBetaCorr_->getTH1();
      assert(histogram_relCombinedIso_wDeltaBetaCorr_);
      
      TString histogramName_relCombinedIso_wRhoCorr = Form("relCombinedIso_wRhoCorr_%s", histogramName_suffix.Data());
      me_relCombinedIso_wRhoCorr_ = dqmStore.book1D(histogramName_relCombinedIso_wRhoCorr.Data(), histogramName_relCombinedIso_wRhoCorr.Data(), 100, 0., 1.);
      histogram_relCombinedIso_wRhoCorr_ = me_relCombinedIso_wRhoCorr_->getTH1();
      assert(histogram_relCombinedIso_wRhoCorr_);

      TString histogramName_sumChargedIsoPileup = Form("sumChargedIsoPileup_%s", histogramName_suffix.Data());
      me_sumChargedIsoPileup_ = dqmStore.book1D(histogramName_sumChargedIsoPileup.Data(), histogramName_sumChargedIsoPileup.Data(), 250, 0., 25.);
      histogram_sumChargedIsoPileup_ = me_sumChargedIsoPileup_->getTH1();
      assert(histogram_sumChargedIsoPileup_);

      TString histogramName_sumNeutralIso_vs_sumChargedIsoPileup = Form("sumNeutralIso_vs_sumChargedIsoPileup_%s", histogramName_suffix.Data());
      me_sumNeutralIso_vs_sumChargedIsoPileup_ = dqmStore.book2D(histogramName_sumNeutralIso_vs_sumChargedIsoPileup.Data(), histogramName_sumNeutralIso_vs_sumChargedIsoPileup.Data(), 50, 0., 25., 50, 0., 25.);
      histogram_sumNeutralIso_vs_sumChargedIsoPileup_ = dynamic_cast<TH2*>(me_sumNeutralIso_vs_sumChargedIsoPileup_->getTH1());
      assert(histogram_sumNeutralIso_vs_sumChargedIsoPileup_);

      TString histogramName_sumNeutralIso_vs_rhoCorr = Form("sumNeutralIso_vs_rhoCorr_%s", histogramName_suffix.Data());
      me_sumNeutralIso_vs_rhoCorr_ = dqmStore.book2D(histogramName_sumNeutralIso_vs_rhoCorr.Data(), histogramName_sumNeutralIso_vs_rhoCorr.Data(), 50, 0., 25., 50, 0., 25.);
      histogram_sumNeutralIso_vs_rhoCorr_ = dynamic_cast<TH2*>(me_sumNeutralIso_vs_rhoCorr_->getTH1());
      assert(histogram_sumNeutralIso_vs_rhoCorr_);

      TString histogramName_pt = Form("pt_%s", histogramName_suffix.Data());
      me_pt_ = dqmStore.book1D(histogramName_pt.Data(), histogramName_pt.Data(), 250, 0., 250.);
      histogram_pt_ = me_pt_->getTH1();
      assert(histogram_pt_);
    }
    void fillHistograms(const l1t::TallinnL1PFTau& l1Tau, double rhoCorr, double evtWeight)
    {
      if ( (l1Tau.pt()              > min_pt_     || min_pt_     <= 0.) && (l1Tau.pt()              < max_pt_     || max_pt_     <= 0.) &&
	   (TMath::Abs(l1Tau.eta()) > min_absEta_ || min_absEta_ <= 0.) && (TMath::Abs(l1Tau.eta()) < max_absEta_ || max_absEta_ <= 0.) )
      {
	double sumChargedIso = l1Tau.sumChargedIso();
	histogram_absChargedIso_->Fill(sumChargedIso, evtWeight);
	histogram_relChargedIso_->Fill(sumChargedIso/TMath::Max(1., l1Tau.pt()), evtWeight);

	double sumNeutralIso = l1Tau.sumNeutralIso();
	histogram_absNeutralIso_->Fill(sumNeutralIso, evtWeight);
	histogram_relNeutralIso_->Fill(sumNeutralIso/TMath::Max(1., l1Tau.pt()), evtWeight);

	double sumCombinedIso = sumChargedIso + sumNeutralIso;
	histogram_absCombinedIso_->Fill(sumCombinedIso, evtWeight);
	histogram_relCombinedIso_->Fill(sumCombinedIso/TMath::Max(1., l1Tau.pt()), evtWeight);

	double sumChargedIsoPileup = l1Tau.sumChargedIsoPileup();
	histogram_sumChargedIsoPileup_->Fill(sumChargedIsoPileup, evtWeight);

	double sumCombinedIso_wDeltaBetaCorr = sumChargedIso + TMath::Max(0., sumNeutralIso - 0.5*sumChargedIsoPileup);
	histogram_absCombinedIso_wDeltaBetaCorr_->Fill(sumCombinedIso_wDeltaBetaCorr, evtWeight);
	histogram_relCombinedIso_wDeltaBetaCorr_->Fill(sumCombinedIso_wDeltaBetaCorr/TMath::Max(1., l1Tau.pt()), evtWeight);

	double sumCombinedIso_wRhoCorr = sumChargedIso + TMath::Max(0., sumNeutralIso - l1Tau.rhoCorr());
	histogram_absCombinedIso_wRhoCorr_->Fill(sumCombinedIso_wRhoCorr, evtWeight);
	histogram_relCombinedIso_wRhoCorr_->Fill(sumCombinedIso_wRhoCorr/TMath::Max(1., l1Tau.pt()), evtWeight);

	histogram_sumNeutralIso_vs_sumChargedIsoPileup_->Fill(sumChargedIsoPileup, sumNeutralIso, evtWeight);
	histogram_sumNeutralIso_vs_rhoCorr_->Fill(rhoCorr, sumNeutralIso, evtWeight);

	histogram_pt_->Fill(l1Tau.pt(), evtWeight);
      }
    }
    void fillHistograms_woGenMatching(const l1t::TallinnL1PFTau& l1Tau, double rhoCorr, double evtWeight)
    {
      fillHistograms(l1Tau, rhoCorr, evtWeight);
    }
    void fillHistograms_wGenMatching(const l1t::TallinnL1PFTau& l1Tau, double rhoCorr, bool isMatched, const std::string& genTau_decayMode, double evtWeight)
    {
      if ( isMatched && (decayMode_ == "all" || genTau_decayMode == decayMode_) )
      {
	fillHistograms(l1Tau, rhoCorr, evtWeight);
      }
    }
    MonitorElement* me_absChargedIso_;
    TH1* histogram_absChargedIso_;
    MonitorElement* me_absNeutralIso_;
    TH1* histogram_absNeutralIso_;
    MonitorElement* me_absCombinedIso_;
    TH1* histogram_absCombinedIso_;
    MonitorElement* me_absCombinedIso_wDeltaBetaCorr_;
    TH1* histogram_absCombinedIso_wDeltaBetaCorr_;
    MonitorElement* me_absCombinedIso_wRhoCorr_;
    TH1* histogram_absCombinedIso_wRhoCorr_;
    MonitorElement* me_relChargedIso_;
    TH1* histogram_relChargedIso_;
    MonitorElement* me_relNeutralIso_;
    TH1* histogram_relNeutralIso_;
    MonitorElement* me_relCombinedIso_;
    TH1* histogram_relCombinedIso_;
    MonitorElement* me_sumChargedIsoPileup_;
    TH1* histogram_sumChargedIsoPileup_;
    MonitorElement* me_relCombinedIso_wDeltaBetaCorr_;
    TH1* histogram_relCombinedIso_wDeltaBetaCorr_;
    MonitorElement* me_relCombinedIso_wRhoCorr_;
    TH1* histogram_relCombinedIso_wRhoCorr_;
    MonitorElement* me_sumNeutralIso_vs_sumChargedIsoPileup_;
    TH2* histogram_sumNeutralIso_vs_sumChargedIsoPileup_;    
    MonitorElement* me_sumNeutralIso_vs_rhoCorr_;
    TH2* histogram_sumNeutralIso_vs_rhoCorr_; 
    MonitorElement* me_pt_;
    TH1* histogram_pt_;    
    double min_pt_;
    double max_pt_;
    double min_absEta_;
    double max_absEta_;
    std::string decayMode_;
    double dRmatch_; 
  };
  std::vector<isolationPlotEntryType*> isolationPlots_;
};

#endif   
