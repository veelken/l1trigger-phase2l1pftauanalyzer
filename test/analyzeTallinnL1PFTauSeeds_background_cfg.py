
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeTallinnL1PFTauSeedsBackground")

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
    )
)

#--------------------------------------------------------------------------------
# set input files

import os
import re

#inputFilePaths = [ '/hdfs/cms/store/user/sbhowmik/NeutrinoGun_E_10GeV/PhaseIIMTDTDRAutumn18MiniAOD_20190505/190505_093529/0000/' ]
inputFilePaths = [
    '/hdfs/cms/store/user/sbhowmik/NeutrinoGun_E_10GeV/PhaseIIMTDTDRAutumn18MiniAOD_20190524/190524_111617/0000/',
    '/hdfs/cms/store/user/sbhowmik/NeutrinoGun_E_10GeV/PhaseIIMTDTDRAutumn18MiniAOD_20190524/190524_111617/0001/',
    '/hdfs/cms/store/user/sbhowmik/NeutrinoGun_E_10GeV/PhaseIIMTDTDRAutumn18MiniAOD_20190524/190524_111617/0002/',
    '/hdfs/cms/store/user/sbhowmik/NeutrinoGun_E_10GeV/PhaseIIMTDTDRAutumn18MiniAOD_20190524/190524_111617/0003/',
    '/hdfs/cms/store/user/sbhowmik/NeutrinoGun_E_10GeV/PhaseIIMTDTDRAutumn18MiniAOD_20190524/190524_111617/0004/',
]    
inputFile_regex = r"[a-zA-Z0-9_/:.-]*NTuple_TallinnL1PFTauProducer_[a-zA-Z0-9-_]+.root"

# check if name of inputFile matches regular expression
inputFileNames = []
for inputFilePath in inputFilePaths:
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

#--------------------------------------------------------------------------------
from L1Trigger.TallinnL1PFTaus.TallinnL1PFTauProducerPF_cff import TallinnL1PFTauProducerPF
process.analyzeTallinnL1PFTauSeeds = cms.EDAnalyzer("TallinnL1PFTauSeedAnalyzer",
  srcL1PFCands = cms.InputTag("l1pfCandidates:PF"),
  srcL1PFJets = cms.InputTag("ak4PFL1PFCorrected"),                      
  srcL1Vertices = cms.InputTag("VertexProducer:l1vertextdr"),
  min_seedChargedPFCand_pt = TallinnL1PFTauProducerPF.min_seedChargedPFCand_pt,
  max_seedChargedPFCand_eta = TallinnL1PFTauProducerPF.max_seedChargedPFCand_eta,
  max_seedChargedPFCand_dz = TallinnL1PFTauProducerPF.max_seedChargedPFCand_dz,
  min_seedPFJet_pt = TallinnL1PFTauProducerPF.min_seedPFJet_pt,
  max_seedPFJet_eta = TallinnL1PFTauProducerPF.max_seedPFJet_eta,
  signalQualityCuts = TallinnL1PFTauProducerPF.signalQualityCuts,             
  dqmDirectory = cms.string("TallinnL1PFTauSeedAnalyzer"),
)
process.analysisSequence += process.analyzeTallinnL1PFTauSeeds
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('TallinnL1PFTauSeedAnalyzer_background_2019Jun03.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
