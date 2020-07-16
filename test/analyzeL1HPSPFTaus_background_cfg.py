
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeTallinnL1PFTausBackground")

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

#--------------------------------------------------------------------------------
# set input files
##if inputFilePath:
##    from HLTTrigger.TallinnHLTPFTauAnalyzer.tools.jobTools import getInputFileNames
##    print("Searching for input files in path = '%s'" % inputFilePath)
##    inputFileNames = getInputFileNames(inputFilePath)
##    print("Found %i input files." % len(inputFileNames))
##else:
##    print("Processing %i input files: %s" % (len(inputFileNames), inputFileNames))
##process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
for useStrips in [ True, False ]:
    moduleNameBase = "L1HPSPFTauProducer"
    moduleLabel = None
    if useStrips:
        moduleLabel = "WithStrips"
    else:
        moduleLabel = "WithoutStrips"

    moduleNamePF_L1HPSPFTauAnalyzerBackground = "analyzeL1HPSPFTaus" + moduleLabel + "PF"
    modulePF_L1HPSPFTauAnalyzerBackground = cms.EDAnalyzer("L1HPSPFTauAnalyzerBackground",
      srcL1PFTaus = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
      dqmDirectory = cms.string("L1HPSPFTauAnalyzerBackground" + moduleLabel + "PF")
    )
    setattr(process, moduleNamePF_L1HPSPFTauAnalyzerBackground, modulePF_L1HPSPFTauAnalyzerBackground)
    process.analysisSequence += getattr(process, moduleNamePF_L1HPSPFTauAnalyzerBackground)

    moduleNamePF_L1HPSPFTauPairAnalyzer = "analyzeL1HPSPFTauPairs" + moduleLabel + "PF"
    modulePF_L1HPSPFTauPairAnalyzer = cms.EDAnalyzer("L1HPSPFTauPairAnalyzer",
      srcL1PFTaus = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
      srcRefTaus = cms.InputTag(''),
      dqmDirectory = cms.string("L1HPSPFTauPairAnalyzer" + moduleLabel + "PF")
    )
    setattr(process, moduleNamePF_L1HPSPFTauPairAnalyzer, modulePF_L1HPSPFTauPairAnalyzer)
    process.analysisSequence += getattr(process, moduleNamePF_L1HPSPFTauPairAnalyzer)

    moduleNamePF_L1HPSPFTauIsolationAnalyzer = "analyzeTallinL1PFTauIsolation" + moduleLabel + "PF"
    modulePF_L1HPSPFTauIsolationAnalyzer = cms.EDAnalyzer("L1HPSPFTauIsolationAnalyzer",
      srcL1PFTaus  = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
      srcGenTaus = cms.InputTag(''),
      dRmatch = cms.double(0.3),                                                            
      srcRho = cms.InputTag('kt6L1PFJetsPF:rho'),
      #inputFileName_rhoCorr = cms.string("L1Trigger/L1HPSPFTauAnalyzer/data/rhoCorr.root"),
      #histogramName_rhoCorr = cms.string("DQMData/RhoCorrAnalyzerPF/neutralPFCandPt_vs_absEta"),                                             
      inputFileName_rhoCorr = cms.string(""),
      histogramName_rhoCorr = cms.string(""),                                                      
      dqmDirectory = cms.string("L1HPSPFTauIsolationAnalyzer" + moduleLabel + "PF")
    )
    setattr(process, moduleNamePF_L1HPSPFTauIsolationAnalyzer, modulePF_L1HPSPFTauIsolationAnalyzer)
    process.analysisSequence += getattr(process, moduleNamePF_L1HPSPFTauIsolationAnalyzer)
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('analyzeL1HPSPFTaus_background_2020Jul16.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

#process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True)
#)
