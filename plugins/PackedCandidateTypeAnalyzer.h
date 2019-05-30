#ifndef L1Trigger_TallinnL1PFTauAnalyzer_PackedCandidateTypeAnalyzer_h
#define L1Trigger_TallinnL1PFTauAnalyzer_PackedCandidateTypeAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "L1Trigger/TallinnL1PFTaus/interface/TallinnL1PFTauQualityCut.h"     // TallinnL1PFTauQualityCut
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"              // pat::PackedCandidate, pat::PackedCandidateCollection
#include "DataFormats/VertexReco/interface/Vertex.h"                          // reco::Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"                       // reco::VertexCollection
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"     // PileupSummaryInfo

#include "L1Trigger/TallinnL1PFTauAnalyzer/interface/histogramAuxFunctions.h" // divideByBinWidth

#include <TH1.h>     // TH1
#include <TString.h> // TString, Form()

#include <vector>    // std::vector
#include <string>    // std::string

class PackedCandidateTypeAnalyzer : public edm::EDAnalyzer 
{
 public:
  PackedCandidateTypeAnalyzer(const edm::ParameterSet& cfg);
  ~PackedCandidateTypeAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_packedCands_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> token_packedCands_;
  
  edm::InputTag srcOfflineVertices_;
  edm::EDGetTokenT<reco::VertexCollection> tokenOfflineVertices_;

  typedef std::vector<PileupSummaryInfo> PileupSummaryInfoCollection;
  edm::InputTag srcPileupSummaryInfo_;
  edm::EDGetTokenT<PileupSummaryInfoCollection> tokenPileupSummaryInfo_;

  std::vector<TallinnL1PFTauQualityCut> isolationQualityCuts_dzCut_disabled_;
  std::vector<TallinnL1PFTauQualityCut> isolationQualityCuts_dzCut_enabled_primary_;
  std::vector<TallinnL1PFTauQualityCut> isolationQualityCuts_dzCut_enabled_pileup_;

  bool applyPuppiWeights_;

  std::string dqmDirectory_;

  struct pfCandTypePlotEntryType
  {
    pfCandTypePlotEntryType(const std::string& label, bool applyPuppiWeights)
      : label_(label)
      , applyPuppiWeights_(applyPuppiWeights)
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

      std::string histogramName_energyFraction_vs_numPileup = Form("%sEnergyFraction_vs_numPileup", label_.data());
      me_energyFraction_vs_numPileup_ = dqmStore.book1D(
        histogramName_energyFraction_vs_numPileup.data(), histogramName_energyFraction_vs_numPileup.data(), 
	40, 0., 400.);
      histogram_energyFraction_vs_numPileup_ = me_energyFraction_vs_numPileup_->getTH1();
      assert(histogram_energyFraction_vs_numPileup_);

      std::string histogramName_ptFraction_vs_numPileup = Form("%sPtFraction_vs_numPileup", label_.data());
      me_ptFraction_vs_numPileup_ = dqmStore.book1D(
        histogramName_ptFraction_vs_numPileup.data(), histogramName_ptFraction_vs_numPileup.data(), 
	40, 0., 400.);
      histogram_ptFraction_vs_numPileup_ = me_ptFraction_vs_numPileup_->getTH1();
      assert(histogram_ptFraction_vs_numPileup_);
    }
    void fillHistograms(const pat::PackedCandidate& packedCand, int numPileup, double evtWeight)
    {
      double weight_energyFraction = packedCand.energy()*evtWeight;
      double weight_ptFraction = packedCand.pt()*evtWeight;
      if ( applyPuppiWeights_ ) 
      {
	double packedCand_puppiWeight = packedCand.puppiWeight();
	weight_energyFraction *= packedCand_puppiWeight;
	weight_ptFraction *= packedCand_puppiWeight;
      }
      histogram_energyFraction_vs_eta_->Fill(packedCand.eta(), weight_energyFraction);
      histogram_ptFraction_vs_eta_->Fill(packedCand.eta(), weight_ptFraction);
      histogram_energyFraction_vs_eta_fine_binning_->Fill(packedCand.eta(), weight_energyFraction);
      histogram_ptFraction_vs_eta_fine_binning_->Fill(packedCand.eta(), weight_ptFraction);
      histogram_energyFraction_vs_absEta_->Fill(TMath::Abs(packedCand.eta()), weight_energyFraction);
      histogram_ptFraction_vs_absEta_->Fill(TMath::Abs(packedCand.eta()), weight_ptFraction);
      histogram_energyFraction_vs_pt_->Fill(packedCand.pt(), weight_energyFraction);
      histogram_ptFraction_vs_pt_->Fill(packedCand.pt(), weight_ptFraction);
      if ( numPileup >= 0 ) 
      {
	histogram_energyFraction_vs_numPileup_->Fill(numPileup, weight_energyFraction);
	histogram_ptFraction_vs_numPileup_->Fill(numPileup, weight_ptFraction);
      }
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
	histogram_energyFraction_vs_numPileup_->Scale(1./numEvents_processed);
	histogram_ptFraction_vs_numPileup_->Scale(1./numEvents_processed);
      }
    }
    std::string label_;
    bool applyPuppiWeights_;
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
    MonitorElement* me_energyFraction_vs_numPileup_;
    TH1* histogram_energyFraction_vs_numPileup_;
    MonitorElement* me_ptFraction_vs_numPileup_;
    TH1* histogram_ptFraction_vs_numPileup_;
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

