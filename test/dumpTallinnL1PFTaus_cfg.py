
import FWCore.ParameterSet.Config as cms

process = cms.Process("dumpTallinnL1PFTaus")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/home/veelken/CMSSW_10_5_0_pre1/src/L1Trigger/TallinnL1PFTaus/test/NTuple_TallinnL1PFTauProducer.root'                        
        #'file:/home/veelken/CMSSW_10_5_0_pre1/src/L1Trigger/TallinnL1PFTaus/test/NTuple_TallinnL1PFTauProducer_2019May13.root'
        #'file:/home/veelken/CMSSW_10_5_0_pre1/src/L1Trigger/TallinnL1PFTaus/test/NTuple_TallinnL1PFTauProducer_DEBUG.root'
    )
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
#process.printGenParticleList = cms.EDAnalyzer("ParticleListDrawer",
#    src = cms.InputTag('prunedGenParticles'),
#    maxEventsToPrint = cms.untracked.int32(100)
#)
#process.analysisSequence += process.analysisSequence
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# CV: restrict analysis to ECAL barrel region (|eta| < 1.47) minus tau isolation cone size, i.e. to |eta| < 1.0,
#     as L1 particle-flow has not been implemented for ECAL endcap region yet
process.selectedGenHadTaus = cms.EDFilter("GenJetSelector",
  src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
  cut = cms.string('pt > 20. & abs(eta) < 1.0'),
  filter = cms.bool(False)
)
process.analysisSequence += process.selectedGenHadTaus

process.genTauLeptons = cms.EDFilter("GenParticleSelector",
  src = cms.InputTag('prunedGenParticles'),
  cut = cms.string("abs(pdgId) = 15"),
  stableOnly = cms.bool(False),
  filter = cms.bool(False)
)
process.analysisSequence += process.genTauLeptons

process.selectedGenTauLeptons = cms.EDFilter("GenParticleAntiOverlapSelector",
  src = cms.InputTag('genTauLeptons'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.5),
  invert = cms.bool(True),
  filter = cms.bool(True)                                                          
)
process.analysisSequence += process.selectedGenTauLeptons

process.genMatchedOfflinePFTaus = cms.EDFilter("PATTauAntiOverlapSelector",
  src = cms.InputTag('slimmedTaus'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(True)                                                          
)
process.analysisSequence += process.genMatchedOfflinePFTaus
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
process.dumpGenTauLeptons = cms.EDAnalyzer("DumpGenParticles",
  src = cms.InputTag('selectedGenTauLeptons')                           
)
process.analysisSequence += process.dumpGenTauLeptons

process.dumpGenHadTaus = cms.EDAnalyzer("DumpGenTaus",
  src = cms.InputTag('selectedGenHadTaus')                           
)
process.analysisSequence += process.dumpGenHadTaus

process.dumpOfflinePFTaus = cms.EDAnalyzer("DumpPATTaus",
  src = cms.InputTag('genMatchedOfflinePFTaus')                           
)
process.analysisSequence += process.dumpOfflinePFTaus
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
process.genMatchedTallinL1PFTausPF = cms.EDFilter("TallinnL1PFTauAntiOverlapSelector",
  src = cms.InputTag('TallinnL1PFTauProducerWithStripsAndPreselectionPF'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(False)                                                          
)
process.analysisSequence += process.genMatchedTallinL1PFTausPF

process.dumpTallinL1PFTausPF = cms.EDAnalyzer("DumpTallinL1PFTaus",
  src = cms.InputTag('genMatchedTallinL1PFTausPF'),
  debug = cms.untracked.bool(True)                                               
)
process.analysisSequence += process.dumpTallinL1PFTausPF

process.genMatchedTallinL1PFTausPuppi = process.genMatchedTallinL1PFTausPF.clone(
  src = cms.InputTag('TallinnL1PFTauProducerWithStripsAndPreselectionPuppi'),
  debug = cms.untracked.bool(True)
)    
process.analysisSequence += process.genMatchedTallinL1PFTausPuppi

process.dumpTallinL1PFTausPuppi = process.dumpTallinL1PFTausPF.clone(
  src = cms.InputTag('genMatchedTallinL1PFTausPuppi')                           
)
process.analysisSequence += process.dumpTallinL1PFTausPuppi

process.genMatchedL1PFCands = cms.EDFilter("L1PFCandidateAntiOverlapSelector",
  src = cms.InputTag('l1pfCandidates:PF'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.4),
  invert = cms.bool(True),
  filter = cms.bool(False)                                                          
)
process.analysisSequence += process.genMatchedL1PFCands

process.dumpL1PFCands = cms.EDAnalyzer("DumpL1PFCandidates",
  src = cms.InputTag('genMatchedL1PFCands'),
  min_pt = cms.double(1.),
  printPuppiWeights = cms.bool(True)                                       
)
process.analysisSequence += process.dumpL1PFCands

process.genMatchedL1PFCandsPuppi = process.genMatchedL1PFCands.clone(
  src = cms.InputTag('l1pfCandidates:Puppi')
)
process.analysisSequence += process.genMatchedL1PFCandsPuppi

process.dumpL1PFCandsPuppi = process.dumpL1PFCands.clone(
  src = cms.InputTag('genMatchedL1PFCandsPuppi'),
  printPuppiWeights = cms.bool(True)
)
process.analysisSequence += process.dumpL1PFCandsPuppi

process.genMatchedL1Tracks = cms.EDProducer("L1TrackAntiOverlapSelector",
  src = cms.InputTag('TTTracksFromTracklet:Level1TTTracks'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.4),
  invert = cms.bool(True)
)
process.analysisSequence += process.genMatchedL1Tracks

process.dumpL1Tracks = cms.EDAnalyzer("DumpL1Tracks",
  src = cms.InputTag('genMatchedL1Tracks')
)
process.analysisSequence += process.dumpL1Tracks

process.dumpL1Vertices = cms.EDAnalyzer("DumpL1Vertices",
  src = cms.InputTag('VertexProducer:l1vertextdr')
)
process.analysisSequence += process.dumpL1Vertices

process.dumpPFVertex = cms.EDAnalyzer("DumpFloat",
  src = cms.InputTag('l1pfProducerBarrel:z0')
)
process.analysisSequence += process.dumpPFVertex

process.genVertex = cms.EDProducer("GenVertexProducer",
  src = cms.InputTag('prunedGenParticles'),
  pdgIds = cms.vint32(-15, +15) # CV: use { -15, +15 } for signal, empty list for background                                    
) 
process.analysisSequence += process.genVertex

process.dumpGenVertex = cms.EDAnalyzer("DumpFloat",
  src = cms.InputTag('genVertex:z0')
)
process.analysisSequence += process.dumpGenVertex

process.dumpOfflineVertices = cms.EDAnalyzer("DumpRecoVertices",
  src = cms.InputTag('offlineSlimmedPrimaryVertices')
)
process.analysisSequence += process.dumpOfflineVertices

process.genMatchedOfflinePFCands = cms.EDFilter("PackedCandidateAntiOverlapSelector",
  src = cms.InputTag('packedPFCandidates'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.4),
  invert = cms.bool(True),
  filter = cms.bool(False)                                                          
)
process.analysisSequence += process.genMatchedOfflinePFCands

process.dumpOfflinePFCands = cms.EDAnalyzer("DumpPackedCandidates",
  src = cms.InputTag('genMatchedOfflinePFCands'),
  min_pt = cms.double(1.)                                                    
)
process.analysisSequence += process.dumpOfflinePFCands
#--------------------------------------------------------------------------------

process.p = cms.Path(process.analysisSequence)


