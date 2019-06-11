
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

        moduleNamePF_TallinnL1PFTauAnalyzerBackground = "analyzeTallinnL1PFTaus" + moduleLabel + "PF"
        modulePF_TallinnL1PFTauAnalyzerBackground = cms.EDAnalyzer("TallinnL1PFTauAnalyzerBackground",
          srcL1PFTaus = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
          dqmDirectory = cms.string("TallinnL1PFTauAnalyzerBackground" + moduleLabel + "PF")
        )
        setattr(process, moduleNamePF_TallinnL1PFTauAnalyzerBackground, modulePF_TallinnL1PFTauAnalyzerBackground)
        process.analysisSequence += getattr(process, moduleNamePF_TallinnL1PFTauAnalyzerBackground)

        moduleNamePF_TallinnL1PFTauPairAnalyzer = "analyzeTallinnL1PFTauPairs" + moduleLabel + "PF"
        modulePF_TallinnL1PFTauPairAnalyzer = cms.EDAnalyzer("TallinnL1PFTauPairAnalyzer",
          srcL1PFTaus = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
          srcRefTaus = cms.InputTag(''),
          dqmDirectory = cms.string("TallinnL1PFTauPairAnalyzer" + moduleLabel + "PF")
        )
        setattr(process, moduleNamePF_TallinnL1PFTauPairAnalyzer, modulePF_TallinnL1PFTauPairAnalyzer)
        process.analysisSequence += getattr(process, moduleNamePF_TallinnL1PFTauPairAnalyzer)

        moduleNamePF_TallinnL1PFTauIsolationAnalyzer = "analyzeTallinL1PFTauIsolation" + moduleLabel + "PF"
        modulePF_TallinnL1PFTauIsolationAnalyzer = cms.EDAnalyzer("TallinnL1PFTauIsolationAnalyzer",
          srcL1PFTaus  = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
          srcGenTaus = cms.InputTag(''),
          dRmatch = cms.double(0.3),                                                            
          srcRho = cms.InputTag('kt6L1PFJetsPF:rho'),
          inputFileName_rhoCorr = cms.string("L1Trigger/TallinnL1PFTauAnalyzer/data/rhoCorr.root"),
          histogramName_rhoCorr = cms.string("DQMData/RhoCorrAnalyzerPF/neutralPFCandPt_vs_absEta"),                                             
          #inputFileName_rhoCorr = cms.string(""),
          #histogramName_rhoCorr = cms.string(""),                                                      
          dqmDirectory = cms.string("TallinnL1PFTauIsolationAnalyzer" + moduleLabel + "PF")
        )
        setattr(process, moduleNamePF_TallinnL1PFTauIsolationAnalyzer, modulePF_TallinnL1PFTauIsolationAnalyzer)
        process.analysisSequence += getattr(process, moduleNamePF_TallinnL1PFTauIsolationAnalyzer)
        
        moduleNamePuppi_TallinnL1PFTauAnalyzerBackground = "analyzeTallinnL1PFTaus" + moduleLabel + "Puppi"
        modulePuppi_TallinnL1PFTauAnalyzerBackground = cms.EDAnalyzer("TallinnL1PFTauAnalyzerBackground",
          srcL1PFTaus = cms.InputTag(moduleNameBase + moduleLabel + "Puppi"),
          dqmDirectory = cms.string("TallinnL1PFTauAnalyzerBackground" + moduleLabel + "Puppi")
        )
        setattr(process, moduleNamePuppi_TallinnL1PFTauAnalyzerBackground, modulePuppi_TallinnL1PFTauAnalyzerBackground)
        process.analysisSequence += getattr(process, moduleNamePuppi_TallinnL1PFTauAnalyzerBackground)

        moduleNamePuppi_TallinnL1PFTauPairAnalyzer = "analyzeTallinnL1PFTauPairs" + moduleLabel + "Puppi"
        modulePuppi_TallinnL1PFTauPairAnalyzer = cms.EDAnalyzer("TallinnL1PFTauPairAnalyzer",
          srcL1PFTaus = cms.InputTag(moduleNameBase + moduleLabel + "Puppi"),
          srcRefTaus = cms.InputTag(''),
          dqmDirectory = cms.string("TallinnL1PFTauPairAnalyzer" + moduleLabel + "Puppi")
        )
        setattr(process, moduleNamePuppi_TallinnL1PFTauPairAnalyzer, modulePuppi_TallinnL1PFTauPairAnalyzer)
        process.analysisSequence += getattr(process, moduleNamePuppi_TallinnL1PFTauPairAnalyzer)

        moduleNamePuppi_TallinnL1PFTauIsolationAnalyzer = "analyzeTallinL1PFTauIsolation" + moduleLabel + "Puppi"
        modulePuppi_TallinnL1PFTauIsolationAnalyzer = cms.EDAnalyzer("TallinnL1PFTauIsolationAnalyzer",
          srcL1PFTaus = cms.InputTag(moduleNameBase + moduleLabel + "Puppi"),
          srcGenTaus = cms.InputTag(''),
          dRmatch = cms.double(0.3),                                                
          srcRho = cms.InputTag('kt6L1PFJetsPuppi:rho'),
          inputFileName_rhoCorr = cms.string("L1Trigger/TallinnL1PFTauAnalyzer/data/rhoCorr.root"),
          histogramName_rhoCorr = cms.string("DQMData/RhoCorrAnalyzerPuppi/neutralPFCandPt_vs_absEta"),
          #inputFileName_rhoCorr = cms.string(""),
          #histogramName_rhoCorr = cms.string(""),                                                                          
          dqmDirectory = cms.string("TallinnL1PFTauIsolationAnalyzer" + moduleLabel + "Puppi")                                                           
        )
        setattr(process, moduleNamePuppi_TallinnL1PFTauIsolationAnalyzer, modulePuppi_TallinnL1PFTauIsolationAnalyzer)
        process.analysisSequence += getattr(process, moduleNamePuppi_TallinnL1PFTauIsolationAnalyzer)
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('TallinnL1PFTauAnalyzer_background_2019Jun07.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

#process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True)
#)
