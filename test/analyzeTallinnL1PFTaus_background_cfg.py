
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeTallinnL1PFTausBackground")

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

inputFilePath = '/hdfs/cms/store/user/sbhowmik/NeutrinoGun_E_10GeV/PhaseIIMTDTDRAutumn18MiniAOD_20190505/190505_093529/0000/'
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
process.analyzeTallinL1PFTausWithStripsPF = cms.EDAnalyzer("TallinnL1PFTauAnalyzerBackground",
  srcTallinnL1PFTaus = cms.InputTag('TallinnL1PFTauProducerWithStripsPF'),
  dqmDirectory = cms.string("TallinnL1PFTauAnalyzerBackgroundWithStripsPF")
)
process.analysisSequence += process.analyzeTallinL1PFTausWithStripsPF

process.analyzeTallinL1PFTausWithoutStripsPF = process.analyzeTallinL1PFTausWithStripsPF.clone(
  srcTallinnL1PFTaus = cms.InputTag('TallinnL1PFTauProducerWithoutStripsPF'),
  dqmDirectory = cms.string("TallinnL1PFTauAnalyzerBackgroundWithoutStripsPF")
)
process.analysisSequence += process.analyzeTallinL1PFTausWithoutStripsPF

process.analyzeTallinL1PFTauIsolationWithStripsPF = cms.EDAnalyzer("TallinnL1PFTauIsolationAnalyzer",
  src = cms.InputTag('TallinnL1PFTauProducerWithStripsPF'),
  srcGenTaus = cms.InputTag(''),                                               
  dqmDirectory = cms.string("TallinnL1PFTauIsolationAnalyzerWithStripsPF")
)
process.analysisSequence += process.analyzeTallinL1PFTauIsolationWithStripsPF

process.analyzeTallinL1PFTauIsolationWithoutStripsPF = process.analyzeTallinL1PFTauIsolationWithStripsPF.clone(
  src = cms.InputTag('TallinnL1PFTauProducerWithoutStripsPF'),
  dqmDirectory = cms.string("TallinnL1PFTauIsolationAnalyzerWithoutStripsPF")
)
process.analysisSequence += process.analyzeTallinL1PFTauIsolationWithoutStripsPF

process.analyzeTallinL1PFTausWithStripsPuppi = process.analyzeTallinL1PFTausWithStripsPF.clone(
  srcTallinnL1PFTaus = cms.InputTag('TallinnL1PFTauProducerWithStripsPuppi'),
  dqmDirectory = cms.string("TallinnL1PFTauAnalyzerBackgroundWithStripsPuppi")
)
process.analysisSequence += process.analyzeTallinL1PFTausWithStripsPuppi

process.analyzeTallinL1PFTausWithoutStripsPuppi = process.analyzeTallinL1PFTausWithStripsPuppi.clone(
  srcTallinnL1PFTaus = cms.InputTag('TallinnL1PFTauProducerWithoutStripsPuppi'),
  dqmDirectory = cms.string("TallinnL1PFTauAnalyzerBackgroundWithoutStripsPuppi")    
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
    outputFileName = cms.string('TallinnL1PFTauAnalyzer_background_2019May14.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

#process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True)
#)
