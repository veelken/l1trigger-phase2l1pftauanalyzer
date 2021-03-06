
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeRho")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/home/veelken/CMSSW_10_5_0_pre1/src/L1Trigger/TallinnL1PFTaus/test/NTuple_TallinnL1PFTauProducer.root'
        #'file:/afs/cern.ch/user/v/veelken/cmssw/CMSSW_10_5_0_pre1/src/L1Trigger/TallinnL1PFTaus/test/NTuple_TallinnL1PFTauProducer.root'
    )
)

#--------------------------------------------------------------------------------
# set input files

import os
import re

inputFilePath = '/hdfs/cms/store/user/sbhowmik/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack/PhaseIIMTDTDRAutumn18MiniAOD_20190524/190524_111901/0000/'
inputFile_regex = r"[a-zA-Z0-9_/:.-]*NTuple_TallinnL1PFTauProducer_[a-zA-Z0-9-_]+.root"

# check if name of inputFile matches regular expression
inputFileNames = []
files = [ "".join([ "file:", inputFilePath, file ]) for file in os.listdir(inputFilePath) ]
for file in files:
    inputFile_matcher = re.compile(inputFile_regex)
    if inputFile_matcher.match(file):
        inputFileNames.append(file)
print "inputFileNames = %s" % inputFileNames 

process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')

process.analysisSequence = cms.Sequence()

from L1Trigger.TallinnL1PFTaus.TallinnL1PFTauProducerPF_cff import TallinnL1PFTauProducerPF
process.analyzeRhoPF = cms.EDAnalyzer("RhoAnalyzer",
  srcRho = cms.InputTag('kt6L1PFJetsPF:rho'),
  srcRhoNeutral = cms.InputTag('kt6L1PFJetsNeutralsPF:rho'),
  dqmDirectory = cms.string("RhoAnalyzerPF"),
)
process.analysisSequence += process.analyzeRhoPF

from L1Trigger.TallinnL1PFTaus.TallinnL1PFTauProducerPuppi_cff import TallinnL1PFTauProducerPuppi
process.analyzeRhoPuppi = cms.EDAnalyzer("RhoAnalyzer",
  srcRho = cms.InputTag('kt6L1PFJetsPuppi:rho'),
  srcRhoNeutral = cms.InputTag('kt6L1PFJetsNeutralsPuppi:rho'),
  dqmDirectory = cms.string("RhoAnalyzerPuppi"),
)
process.analysisSequence += process.analyzeRhoPuppi

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('RhoAnalyzer_signal_2019May28.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
