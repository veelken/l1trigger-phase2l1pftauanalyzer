
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeTallinnL1PFTausBackground")

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

inputFilePaths = [
    '/hdfs/cms/store/user/sbhowmik/NeutrinoGun_E_10GeV/PhaseIIMTDTDRAutumn18MiniAOD_20190617/190618_085348/0000/',
    '/hdfs/cms/store/user/sbhowmik/NeutrinoGun_E_10GeV/PhaseIIMTDTDRAutumn18MiniAOD_20190617/190618_085348/0001/',
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
#process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
for useStrips in [ True, False ]:
    for applyPreselection in [ True, False ]:
        moduleNameBase = "TallinnL1PFTauProducer"
        moduleLabel = None
        if useStrips and applyPreselection:
            moduleLabel = "WithStripsAndPreselection"
        elif useStrips and not applyPreselection:
            moduleLabel = "WithStripsWithoutPreselection"
        elif not useStrips and applyPreselection:
            moduleLabel = "WithoutStripsWithPreselection"
        elif not useStrips and not applyPreselection:
            moduleLabel = "WithoutStripsAndPreselection"
        else:
            raise ValueError("Invalid Combination of 'useStrips' and 'applyPreselection' Configuration parameters !!")

        # TallinnL1PFTaus built from PFCandidates without PUPPI weights
        moduleNamePF_TallinnL1PFTauQualityAnalyzer = "analyzeTallinnL1PFTauQuality" + moduleLabel + "PF"
        modulePF_TallinnL1PFTauQualityAnalyzer = cms.EDAnalyzer("TallinnL1PFTauQualityAnalyzer",
          src = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
          dqmDirectory = cms.string("TallinnL1PFTauQualityAnalyzer" + moduleLabel + "PF")
        )
        setattr(process, moduleNamePF_TallinnL1PFTauQualityAnalyzer, modulePF_TallinnL1PFTauQualityAnalyzer)
        process.analysisSequence += getattr(process, moduleNamePF_TallinnL1PFTauQualityAnalyzer)
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('TallinnL1PFTauQualityAnalyzer_background_2019Nov20.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

#process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True)
#)
