#ifndef L1Trigger_TallinnL1PFTauAnalyzer_L1PFCandidateTypeAnalyzer_h
#define L1Trigger_TallinnL1PFTauAnalyzer_L1PFCandidateTypeAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "L1Trigger/TallinnL1PFTaus/interface/TallinnL1PFTauQualityCut.h"     // TallinnL1PFTauQualityCut
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"           // l1t::PFCandidate
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidateFwd.h"        // l1t::PFCandidateCollection
#include "DataFormats/L1TVertex/interface/Vertex.h"                           // l1t::Vertex
#include "DataFormats/L1TVertex/interface/VertexFwd.h"                        // l1t::VertexCollection

#include "L1Trigger/TallinnL1PFTauAnalyzer/interface/histogramAuxFunctions.h" // divideByBinWidth

#include <TH1.h>     // TH1
#include <TString.h> // TString, Form()

#include <vector>    // std::vector
#include <string>    // std::string

class L1PFCandidateTypeAnalyzer : public edm::EDAnalyzer 
{
 public:
  L1PFCandidateTypeAnalyzer(const edm::ParameterSet& cfg);
  ~L1PFCandidateTypeAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_l1PFCands_;
  edm::EDGetTokenT<l1t::PFCandidateCollection> token_l1PFCands_;

  edm::InputTag srcL1Vertices_;
  edm::EDGetTokenT<l1t::VertexCollection> tokenL1Vertices_;

  std::vector<TallinnL1PFTauQualityCut> isolationQualityCuts_dzCut_disabled_;
  std::vector<TallinnL1PFTauQualityCut> isolationQualityCuts_dzCut_enabled_primary_;
  std::vector<TallinnL1PFTauQualityCut> isolationQualityCuts_dzCut_enabled_pileup_;

  std::string dqmDirectory_;

  struct pfCandTypePlotEntryType
  {
    pfCandTypePlotEntryType(const std::string& label)
      : label_(label)
      , me_energyFraction_vs_eta_(nullptr)
      , histogram_energyFraction_vs_eta_(nullptr)
      , me_ptFraction_vs_eta_(nullptr)
      , histogram_ptFraction_vs_eta_(nullptr)
      , me_energyFraction_vs_eta_fine_binning_(nullptr)
      , histogram_energyFraction_vs_eta_fine_binning_(nullptr)
      , me_ptFraction_vs_eta_fine_binning_(nullptr)
      , histogram_ptFraction_vs_eta_fine_binning_(nullptr)
      , me_energyFraction_vs_absEta_(nullptr)
      , histogram_energyFraction_vs_absEta_(nullptr)
      , me_ptFraction_vs_absEta_(nullptr)
      , histogram_ptFraction_vs_absEta_(nullptr)
      , me_energyFraction_vs_pt_(nullptr)
      , histogram_energyFraction_vs_pt_(nullptr)
      , me_ptFraction_vs_pt_(nullptr)
      , histogram_ptFraction_vs_pt_(nullptr)
    {}
    ~pfCandTypePlotEntryType()
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      std::string histogramName_energyFraction_vs_eta = Form("%sEnergyFraction_vs_eta", label_.data());
      me_energyFraction_vs_eta_ = dqmStore.book1D(
        histogramName_energyFraction_vs_eta.data(), histogramName_energyFraction_vs_eta.data(), 
	68, -34*0.087, +34*0.087); // 0.087 = size of one CaloTower (5 ECAL crystals) in eta
      histogram_energyFraction_vs_eta_ = me_energyFraction_vs_eta_->getTH1();
      assert(histogram_energyFraction_vs_eta_);

      std::string histogramName_ptFraction_vs_eta = Form("%sPtFraction_vs_eta", label_.data());
      me_ptFraction_vs_eta_ = dqmStore.book1D(
        histogramName_ptFraction_vs_eta.data(), histogramName_ptFraction_vs_eta.data(), 
	68, -34*0.087, +34*0.087); // 0.087 = size of one CaloTower (5 ECAL crystals) in eta
      histogram_ptFraction_vs_eta_ = me_ptFraction_vs_eta_->getTH1();
      assert(histogram_ptFraction_vs_eta_);

      std::string histogramName_energyFraction_vs_eta_fine_binning = Form("%sEnergyFraction_vs_eta_fine_binning", label_.data());
      me_energyFraction_vs_eta_fine_binning_ = dqmStore.book1D(
        histogramName_energyFraction_vs_eta_fine_binning.data(), histogramName_energyFraction_vs_eta_fine_binning.data(), 
	6000, -3.0, +3.0);
      histogram_energyFraction_vs_eta_fine_binning_ = me_energyFraction_vs_eta_fine_binning_->getTH1();
      assert(histogram_energyFraction_vs_eta_fine_binning_);

      std::string histogramName_ptFraction_vs_eta_fine_binning = Form("%sPtFraction_vs_eta_fine_binning", label_.data());
      me_ptFraction_vs_eta_fine_binning_ = dqmStore.book1D(
        histogramName_ptFraction_vs_eta_fine_binning.data(), histogramName_ptFraction_vs_eta_fine_binning.data(), 
	6000, -3.0, +3.0);
      histogram_ptFraction_vs_eta_fine_binning_ = me_ptFraction_vs_eta_fine_binning_->getTH1();
      assert(histogram_ptFraction_vs_eta_fine_binning_);

      std::string histogramName_energyFraction_vs_absEta = Form("%sEnergyFraction_vs_absEta", label_.data());
      me_energyFraction_vs_absEta_ = dqmStore.book1D(
        histogramName_energyFraction_vs_absEta.data(), histogramName_energyFraction_vs_absEta.data(), 
	34, 0., 34*0.087); // 0.087 = size of one CaloTower (5 ECAL crystals) in eta
      histogram_energyFraction_vs_absEta_ = me_energyFraction_vs_absEta_->getTH1();
      assert(histogram_energyFraction_vs_absEta_);

      std::string histogramName_ptFraction_vs_absEta = Form("%sPtFraction_vs_absEta", label_.data());
      me_ptFraction_vs_absEta_ = dqmStore.book1D(
        histogramName_ptFraction_vs_absEta.data(), histogramName_ptFraction_vs_absEta.data(), 
	34, 0., 34*0.087); // 0.087 = size of one CaloTower (5 ECAL crystals) in eta
      histogram_ptFraction_vs_absEta_ = me_ptFraction_vs_absEta_->getTH1();
      assert(histogram_ptFraction_vs_absEta_);

      const int numBins_pt = 15;
      float binning_pt[numBins_pt + 1] = { 
        0., 1., 2., 4., 6., 10., 15., 25., 40., 60., 100., 150., 250., 400., 600., 1000.
      };

      std::string histogramName_energyFraction_vs_pt = Form("%sEnergyFraction_vs_pt", label_.data());
      me_energyFraction_vs_pt_ = dqmStore.book1D(
        histogramName_energyFraction_vs_pt.data(), histogramName_energyFraction_vs_pt.data(), 
	numBins_pt, binning_pt);
      histogram_energyFraction_vs_pt_ = me_energyFraction_vs_pt_->getTH1();
      assert(histogram_energyFraction_vs_pt_);

      std::string histogramName_ptFraction_vs_pt = Form("%sPtFraction_vs_pt", label_.data());
      me_ptFraction_vs_pt_ = dqmStore.book1D(
        histogramName_ptFraction_vs_pt.data(), histogramName_ptFraction_vs_pt.data(), 
	numBins_pt, binning_pt);
      histogram_ptFraction_vs_pt_ = me_ptFraction_vs_pt_->getTH1();
      assert(histogram_ptFraction_vs_pt_);
    }
    void fillHistograms(const l1t::PFCandidate& l1PFCand, double evtWeight)
    {
      histogram_energyFraction_vs_eta_->Fill(l1PFCand.eta(), l1PFCand.energy()*evtWeight);
      histogram_ptFraction_vs_eta_->Fill(l1PFCand.eta(), l1PFCand.pt()*evtWeight);
      histogram_energyFraction_vs_eta_fine_binning_->Fill(l1PFCand.eta(), l1PFCand.energy()*evtWeight);
      histogram_ptFraction_vs_eta_fine_binning_->Fill(l1PFCand.eta(), l1PFCand.pt()*evtWeight);
      histogram_energyFraction_vs_absEta_->Fill(TMath::Abs(l1PFCand.eta()), l1PFCand.energy()*evtWeight);
      histogram_ptFraction_vs_absEta_->Fill(TMath::Abs(l1PFCand.eta()), l1PFCand.pt()*evtWeight);
      histogram_energyFraction_vs_pt_->Fill(l1PFCand.pt(), l1PFCand.energy()*evtWeight);
      histogram_ptFraction_vs_pt_->Fill(l1PFCand.pt(), l1PFCand.pt()*evtWeight);
    }
    void normalizeHistograms(double numEvents_processed)
    {
      if ( numEvents_processed > 0. ) 
      {
	divideByBinWidth(histogram_energyFraction_vs_eta_);
	histogram_energyFraction_vs_eta_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_ptFraction_vs_eta_);
	histogram_ptFraction_vs_eta_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_energyFraction_vs_eta_fine_binning_);
	histogram_energyFraction_vs_eta_fine_binning_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_ptFraction_vs_eta_fine_binning_);
	histogram_ptFraction_vs_eta_fine_binning_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_energyFraction_vs_absEta_);
	histogram_energyFraction_vs_absEta_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_ptFraction_vs_absEta_);
	histogram_ptFraction_vs_absEta_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_energyFraction_vs_pt_);
	histogram_energyFraction_vs_pt_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_ptFraction_vs_pt_);
	histogram_ptFraction_vs_pt_->Scale(1./numEvents_processed);
      }
    }
    std::string label_;
    MonitorElement* me_energyFraction_vs_eta_;
    TH1* histogram_energyFraction_vs_eta_;
    MonitorElement* me_ptFraction_vs_eta_;
    TH1* histogram_ptFraction_vs_eta_;
    MonitorElement* me_energyFraction_vs_eta_fine_binning_;
    TH1* histogram_energyFraction_vs_eta_fine_binning_;
    MonitorElement* me_ptFraction_vs_eta_fine_binning_;
    TH1* histogram_ptFraction_vs_eta_fine_binning_;
    MonitorElement* me_energyFraction_vs_absEta_;
    TH1* histogram_energyFraction_vs_absEta_;
    MonitorElement* me_ptFraction_vs_absEta_;
    TH1* histogram_ptFraction_vs_absEta_;
    MonitorElement* me_energyFraction_vs_pt_;
    TH1* histogram_energyFraction_vs_pt_;
    MonitorElement* me_ptFraction_vs_pt_;
    TH1* histogram_ptFraction_vs_pt_;
  };
  pfCandTypePlotEntryType* pfChargedHadronPlots_;
  pfCandTypePlotEntryType* pfChargedHadronPileupPlots_;
  pfCandTypePlotEntryType* pfElectronPlots_;
  pfCandTypePlotEntryType* pfNeutralHadronPlots_;
  pfCandTypePlotEntryType* pfPhotonPlots_;
  pfCandTypePlotEntryType* pfMuonPlots_;
  
  MonitorElement* me_EventCounter_;
  TH1* histogram_EventCounter_;
};

#endif   

