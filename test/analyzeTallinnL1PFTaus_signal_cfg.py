
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeTallinnL1PFTausSignal")

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
    )
)

#--------------------------------------------------------------------------------
# set input files

import os
import re

inputFilePath = '/hdfs/cms/store/user/sbhowmik/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack/PhaseIIMTDTDRAutumn18MiniAOD_20190505/190505_093730/0000/'
inputFile_regex = r"[a-zA-Z0-9_/:.-]*NTuple_TallinnL1PFTauProducer_[a-zA-Z0-9-_]+.root"

# check if name of inputFile matches regular expression
inputFileNames = []
files = [ "".join([ "file:", inputFilePath, file ]) for file in os.listdir(inputFilePath) ]
for file in files:
    inputFile_matcher = re.compile(inputFile_regex)
    if inputFile_matcher.match(file):
        inputFileNames.append(file)
print "inputFileNames = %s" % inputFileNames 

#process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
process.selectedGenHadTaus = cms.EDFilter("GenJetSelector",
  src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
  cut = cms.string('pt > 20. & abs(eta) < 1.0'),
  filter = cms.bool(False)
)
process.analysisSequence += process.selectedGenHadTaus

process.genMatchedOfflinePFTaus = cms.EDFilter("PATTauAntiOverlapSelector",
  src = cms.InputTag('slimmedTaus'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(True)                                                          
)
process.analysisSequence += process.genMatchedOfflinePFTaus

process.selectedOfflinePFTaus = cms.EDFilter("PATTauSelector",
  src = cms.InputTag('genMatchedOfflinePFTaus'),
  cut = cms.string("pt > 20. & abs(eta) < 1.0 & tauID('decayModeFinding') > 0.5 & tauID('byLooseCombinedIsolationDeltaBetaCorr3Hits') > 0.5")
)
process.analysisSequence += process.selectedOfflinePFTaus

process.selectedOfflinePFTauFilter = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('selectedOfflinePFTaus'),
  minNumber = cms.uint32(1)
)
process.analysisSequence += process.selectedOfflinePFTauFilter

process.offlineMatchedGenHadTaus = cms.EDFilter("GenJetAntiOverlapSelector",
  src = cms.InputTag('selectedGenHadTaus'),
  srcNotToBeFiltered = cms.VInputTag('selectedOfflinePFTaus'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(True)                                                          
)
process.analysisSequence += process.offlineMatchedGenHadTaus

process.genVertex = cms.EDProducer("GenVertexProducer",
  src = cms.InputTag('prunedGenParticles'),
  pdgIds = cms.vint32(-15, +15) # CV: use { -15, +15 } for signal, empty list for background                                    
) 
process.analysisSequence += process.genVertex

process.analyzeVertices = cms.EDAnalyzer("L1VertexAnalyzer",
  srcGenVertex_z = cms.InputTag('genVertex:z0'),
  srcL1Vertices = cms.InputTag('VertexProducer:l1vertextdr'),                        
  srcL1PFVertex_z = cms.InputTag('l1pfProducerBarrel:z0'),
  srcOfflineVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
  dqmDirectory = cms.string("L1VertexAnalyzer")
)
process.analysisSequence += process.analyzeVertices

process.analyzeTracks = cms.EDAnalyzer("L1TrackAnalyzer",
  srcGenTaus = cms.InputTag('selectedGenHadTaus'),
  srcOfflineVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
  #srcOfflineTracks = cms.InputTag('generalTracks'),
  srcOfflineTracks = cms.InputTag(''),                                        
  #srcOfflinePFCands = cms.InputTag('particleFlow'),
  srcOfflinePFCands = cms.InputTag(''), 
  srcL1Vertices = cms.InputTag('VertexProducer:l1vertextdr'),                           
  srcL1Tracks = cms.InputTag('TTTracksFromTracklet:Level1TTTracks'),
  srcL1PFVertex_z = cms.InputTag('l1pfProducerBarrel:z0'),
  #srcL1PFCands = cms.InputTag('l1pfCandidates:PF'),
  srcL1PFCands = cms.InputTag(''),                                      
  dqmDirectory = cms.string("L1TrackAnalyzer")
)
process.analysisSequence += process.analyzeTracks

process.analyzeTallinL1PFTausWithStripsPF = cms.EDAnalyzer("TallinnL1PFTauAnalyzerSignal",
  srcNumerator = cms.InputTag('TallinnL1PFTauProducerWithStripsPF'),
  srcDenominator = cms.InputTag('offlineMatchedGenHadTaus'),                                               
  dqmDirectory = cms.string("TallinnL1PFTauAnalyzerSignalWithStripsPF")
)
process.analysisSequence += process.analyzeTallinL1PFTausWithStripsPF

process.analyzeTallinL1PFTausWithoutStripsPF = process.analyzeTallinL1PFTausWithStripsPF.clone(
  srcNumerator = cms.InputTag('TallinnL1PFTauProducerWithoutStripsPF'),
  dqmDirectory = cms.string("TallinnL1PFTauAnalyzerSignalWithoutStripsPF")
)
process.analysisSequence += process.analyzeTallinL1PFTausWithoutStripsPF

process.analyzeTallinL1PFTauIsolationWithStripsPF = cms.EDAnalyzer("TallinnL1PFTauIsolationAnalyzer",
  src = cms.InputTag('TallinnL1PFTauProducerWithStripsPF'),
  srcGenTaus = cms.InputTag('offlineMatchedGenHadTaus'),                                               
  dqmDirectory = cms.string("TallinnL1PFTauIsolationAnalyzerWithStripsPF")
)
process.analysisSequence += process.analyzeTallinL1PFTauIsolationWithStripsPF

process.analyzeTallinL1PFTauIsolationWithoutStripsPF = process.analyzeTallinL1PFTauIsolationWithStripsPF.clone(
  src = cms.InputTag('TallinnL1PFTauProducerWithoutStripsPF'),
  dqmDirectory = cms.string("TallinnL1PFTauIsolationAnalyzerWithoutStripsPF")
)
process.analysisSequence += process.analyzeTallinL1PFTauIsolationWithoutStripsPF

process.analyzeTallinL1PFTausWithStripsPuppi = process.analyzeTallinL1PFTausWithStripsPF.clone(
  srcNumerator = cms.InputTag('TallinnL1PFTauProducerWithStripsPuppi'),
  dqmDirectory = cms.string("TallinnL1PFTauAnalyzerSignalWithStripsPuppi")
)
process.analysisSequence += process.analyzeTallinL1PFTausWithStripsPuppi

process.analyzeTallinL1PFTausWithoutStripsPuppi = process.analyzeTallinL1PFTausWithStripsPuppi.clone(
  srcNumerator = cms.InputTag('TallinnL1PFTauProducerWithoutStripsPuppi'),
  dqmDirectory = cms.string("TallinnL1PFTauAnalyzerSignalWithoutStripsPuppi")
)
process.analysisSequence += process.analyzeTallinL1PFTausWithoutStripsPuppi

process.analyzeTallinL1PFTauIsolationWithStripsPuppi = process.analyzeTallinL1PFTauIsolationWithStripsPF.clone(
  src = cms.InputTag('TallinnL1PFTauProducerWithStripsPuppi'),
  dqmDirectory = cms.string("TallinnL1PFTauIsolationAnalyzerWithStripsPuppi")
)
process.analysisSequence += process.analyzeTallinL1PFTauIsolationWithStripsPuppi

process.analyzeTallinL1PFTauIsolationWithoutStripsPuppi = process.analyzeTallinL1PFTauIsolationWithStripsPuppi.clone(
  src = cms.InputTag('TallinnL1PFTauProducerWithoutStripsPuppi'),
  dqmDirectory = cms.string("TallinnL1PFTauIsolationAnalyzerWithoutStripsPuppi")
)
process.analysisSequence += process.analyzeTallinL1PFTauIsolationWithoutStripsPuppi
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('TallinnL1PFTauAnalyzer_signal_2019May14.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

#process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True)
#)
