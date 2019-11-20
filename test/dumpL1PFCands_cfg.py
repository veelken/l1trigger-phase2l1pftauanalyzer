
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
        #'file:/home/sbhowmik/Phase2/Test/test1/CMSSW_10_5_0_pre1/src/L1Trigger/TallinnL1PFTaus/test/NTuple_TallinnL1PFTauProducer_Tallinn.root'
        #'file:/home/sbhowmik/Phase2/Test/test1/CMSSW_10_5_0_pre1/src/L1Trigger/TallinnL1PFTaus/test/NTuple_TallinnL1PFTauProducer_Tallinn_v2.root'
        'file:/home/sbhowmik/Phase2/Test/test1/CMSSW_10_5_0_pre1/src/L1Trigger/TallinnL1PFTaus/test/NTuple_TallinnL1PFTauProducer_Tallinn_v3.root'                           
        #'file:/home/sbhowmik/Phase2/Test/test1/CMSSW_10_5_0_pre1/src/L1Trigger/Phase2L1Taus/test/NTuple_TallinnL1PFTauProducer_Isobel.root'
    )
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
process.dumpGenParticles = cms.EDAnalyzer("DumpGenParticles",
  src = cms.InputTag('prunedGenParticles')                           
)
process.analysisSequence += process.dumpGenParticles

process.dumpL1PFCands = cms.EDAnalyzer("DumpL1PFCandidates",
  src = cms.InputTag('l1pfCandidates:PF'),
  min_pt = cms.double(1.),
  printPuppiWeights = cms.bool(True)                                       
)
process.analysisSequence += process.dumpL1PFCands

process.dumpOfflinePFCands = cms.EDAnalyzer("DumpPackedCandidates",
  src = cms.InputTag('packedPFCandidates'),
  min_pt = cms.double(1.)                                                    
)
process.analysisSequence += process.dumpOfflinePFCands
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
from L1Trigger.TallinnL1PFTaus.TallinnL1PFTauProducerPF_cff import TallinnL1PFTauProducerPF
process.analyzeL1PFCandidateTypePF = cms.EDAnalyzer("L1PFCandidateTypeAnalyzer",
  srcL1PFCands = cms.InputTag('l1pfCandidates:PF'),
  srcL1Vertices = cms.InputTag('VertexProducer:l1vertextdr'),
  srcPileupSummaryInfo = cms.InputTag(''),              
  isolationQualityCuts = TallinnL1PFTauProducerPF.isolationQualityCuts,             
  dqmDirectory = cms.string("L1PFCandidateTypeAnalyzerPF"),
)
process.analysisSequence += process.analyzeL1PFCandidateTypePF

process.analyzeOfflinePFCandidateTypePF = cms.EDAnalyzer("PackedCandidateTypeAnalyzer",
  srcPackedCands = cms.InputTag('packedPFCandidates'),
  srcOfflineVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
  srcPileupSummaryInfo = cms.InputTag(''),              
  isolationQualityCuts = TallinnL1PFTauProducerPF.isolationQualityCuts,
  applyPuppiWeights = cms.bool(False),
  dqmDirectory = cms.string("PackedCandidateTypeAnalyzerPF"),
)
process.analysisSequence += process.analyzeOfflinePFCandidateTypePF

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('dumpL1PFCands_Tallinn_v3.root')                                
    #outputFileName = cms.string('dumpL1PFCands_Isobel.root')
)
#--------------------------------------------------------------------------------

process.p = cms.Path(process.analysisSequence + process.savePlots)


