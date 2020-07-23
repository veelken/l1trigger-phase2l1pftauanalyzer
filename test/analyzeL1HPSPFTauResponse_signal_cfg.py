
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeTallinnL1PFTausSignal")

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

#--------------------------------------------------------------------------------
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.analysisSequence += process.genParticles

process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.tauGenJets.GenParticles = cms.InputTag('genParticles')
process.analysisSequence += process.tauGenJets

process.load("PhysicsTools.JetMCAlgos.TauGenJetsDecayModeSelectorAllHadrons_cfi")
process.analysisSequence += process.tauGenJetsSelectorAllHadrons

process.selectedGenHadTaus = cms.EDFilter("GenJetSelector",
  src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
  cut = cms.string('pt > 20. & abs(eta) < 2.4'),
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
  cut = cms.string("pt > 20. & abs(eta) < 2.4 & tauID('decayModeFinding') > 0.5 & tauID('byLooseCombinedIsolationDeltaBetaCorr3Hits') > 0.5")
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

for useStrips in [ True, False ]:
    moduleNameBase = "L1HPSPFTauProducer"
    moduleLabel = None
    if useStrips:
        moduleLabel = "WithStrips"
    else:
        moduleLabel = "WithoutStrips"

    moduleNamePF_L1HPSPFTauResponseAnalyzerWrtGenHadTaus = "analyzeL1HPSPFTauResponse" + moduleLabel + "PFWrtGenHadTaus"
    modulePF_L1HPSPFTauResponseAnalyzerWrtGenHadTaus = cms.EDAnalyzer("L1HPSPFTauResponseAnalyzer",
      srcL1PFTaus = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
      srcRefTaus = cms.InputTag('offlineMatchedGenHadTaus'),
      typeRefTaus = cms.string("gen"),                                                                            
      dqmDirectory = cms.string("L1HPSPFTauResponseAnalyzer" + moduleLabel + "PF_wrtGenHadTaus")
    )
    setattr(process, moduleNamePF_L1HPSPFTauResponseAnalyzerWrtGenHadTaus, modulePF_L1HPSPFTauResponseAnalyzerWrtGenHadTaus)
    process.analysisSequence += getattr(process, moduleNamePF_L1HPSPFTauResponseAnalyzerWrtGenHadTaus)

    moduleNamePF_L1HPSPFTauResponseAnalyzerWrtOfflineTaus = "analyzeL1HPSPFTauResponse" + moduleLabel + "PFWrtOfflineTaus"
    modulePF_L1HPSPFTauResponseAnalyzerWrtOfflineTaus = cms.EDAnalyzer("L1HPSPFTauResponseAnalyzer",
      srcL1PFTaus = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
      srcRefTaus = cms.InputTag('selectedOfflinePFTaus'),
      typeRefTaus = cms.string("offline"),                                                                   
      dqmDirectory = cms.string("L1HPSPFTauResponseAnalyzer" + moduleLabel + "PF_wrtOfflineTaus")
    )
    setattr(process, moduleNamePF_L1HPSPFTauResponseAnalyzerWrtOfflineTaus, modulePF_L1HPSPFTauResponseAnalyzerWrtOfflineTaus)
    process.analysisSequence += getattr(process, moduleNamePF_L1HPSPFTauResponseAnalyzerWrtOfflineTaus)
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('analyzeL1HPSPFTauResponse_signal_2020Jul23.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
