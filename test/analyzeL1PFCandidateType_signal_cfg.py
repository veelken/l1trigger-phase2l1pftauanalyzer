
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeL1PFCandidateType")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/home/veelken/Phase2HLT/CMSSW_11_1_0/src/HLTrigger/Phase2HLTPFTaus/test/step3_RAW2DIGI_RECO.root'
    )
)

inputFilePath = '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5_wOfflineVtx_wL1_2FM/'
inputFileNames = None

#--------------------------------------------------------------------------------
# set input files
if inputFilePath:
    from HLTrigger.TallinnHLTPFTauAnalyzer.tools.jobTools import getInputFileNames
    print("Searching for input files in path = '%s'" % inputFilePath)
    inputFileNames = getInputFileNames(inputFilePath)
    print("Found %i input files." % len(inputFileNames))
    process.source.fileNames = cms.untracked.vstring(inputFileNames)
else:
    print("Processing %i input files: %s" % (len(inputFileNames), inputFileNames))
    process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.analysisSequence = cms.Sequence()

from L1Trigger.Phase2L1Taus.L1HPSPFTauProducerPF_cfi import L1HPSPFTauProducerPF
process.analyzeL1PFCandidateTypePF = cms.EDAnalyzer("L1PFCandidateTypeAnalyzer",
  srcL1PFCands = cms.InputTag('l1pfCandidates:PF'),
  srcL1Vertices = cms.InputTag('L1TkPrimaryVertex'),
  #srcPileupSummaryInfo = cms.InputTag('slimmedAddPileupInfo'),
  srcPileupSummaryInfo = cms.InputTag(''),              
  isolationQualityCuts = L1HPSPFTauProducerPF.isolationQualityCuts,             
  dqmDirectory = cms.string("L1PFCandidateTypeAnalyzerPF"),
)
process.analysisSequence += process.analyzeL1PFCandidateTypePF

from L1Trigger.Phase2L1Taus.L1HPSPFTauProducerPuppi_cfi import L1HPSPFTauProducerPuppi
process.analyzeL1PFCandidateTypePuppi = cms.EDAnalyzer("L1PFCandidateTypeAnalyzer",
  srcL1PFCands = cms.InputTag('l1pfCandidates:Puppi'),
  srcL1Vertices = cms.InputTag('L1TkPrimaryVertex'),
  #srcPileupSummaryInfo = cms.InputTag('slimmedAddPileupInfo'),
  srcPileupSummaryInfo = cms.InputTag(''),                                                       
  isolationQualityCuts = L1HPSPFTauProducerPuppi.isolationQualityCuts,                                         
  dqmDirectory = cms.string("L1PFCandidateTypeAnalyzerPuppi"),
)
process.analysisSequence += process.analyzeL1PFCandidateTypePuppi

from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts
process.analyzeOfflinePFCandidateTypePF = cms.EDAnalyzer("PackedCandidateTypeAnalyzer",
  srcPackedCands = cms.InputTag('packedPFCandidates'),
  srcVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
  #srcPileupSummaryInfo = cms.InputTag('slimmedAddPileupInfo'),
  srcPileupSummaryInfo = cms.InputTag(''),              
  isolationQualityCuts = PFTauQualityCuts.isolationQualityCuts,
  applyPuppiWeights = cms.bool(False),
  dqmDirectory = cms.string("PackedCandidateTypeAnalyzerPF"),
)
process.analysisSequence += process.analyzeOfflinePFCandidateTypePF

process.analyzeOfflinePFCandidateTypePuppi = cms.EDAnalyzer("PackedCandidateTypeAnalyzer",
  srcPackedCands = cms.InputTag('packedPFCandidates'),
  srcVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
  #srcPileupSummaryInfo = cms.InputTag('slimmedAddPileupInfo'),
  srcPileupSummaryInfo = cms.InputTag(''),              
  isolationQualityCuts = PFTauQualityCuts.isolationQualityCuts,
  applyPuppiWeights = cms.bool(True),
  dqmDirectory = cms.string("PackedCandidateTypeAnalyzerPuppi"),
)
process.analysisSequence += process.analyzeOfflinePFCandidateTypePuppi

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('L1PFCandidateTypeAnalyzer_signal_2020Jul23.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
