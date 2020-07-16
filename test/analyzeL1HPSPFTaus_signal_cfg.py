
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

process.analyzeVertices = cms.EDAnalyzer("L1VertexAnalyzer",
  srcGenVertex_z = cms.InputTag('genVertex:z0'),
  srcL1Vertices = cms.InputTag('L1TkPrimaryVertex'),                        
  srcL1PFVertex_z = cms.InputTag('l1pfProducerBarrel:z0'),
  srcOfflineVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
  dqmDirectory = cms.string("L1VertexAnalyzer")
)
process.analysisSequence += process.analyzeVertices

process.analyzeTracksWrtRecVertex = cms.EDAnalyzer("L1TrackAnalyzer",
  srcGenTaus = cms.InputTag('selectedGenHadTaus'),
  vtxMode = cms.string("recVtx"),
  srcOfflineVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),                                       
  srcOfflineTracks = cms.InputTag('generalTracks'),                                                     
  srcOfflinePFCands = cms.InputTag('particleFlow'),
  srcL1Vertices = cms.InputTag('L1TkPrimaryVertex'),                                                  
  srcL1Tracks = cms.InputTag('TTTracksFromTracklet:Level1TTTracks'),
  srcL1PFVertex_z = cms.InputTag('l1pfProducerBarrel:z0'),                                                     
  srcL1PFCands = cms.InputTag('l1pfCandidates:PF'),
  dqmDirectory = cms.string("L1TrackAnalyzerWrtRecVertex"),
  debug = cms.bool(False)                                     
)
process.analysisSequence += process.analyzeTracksWrtRecVertex

process.genVertex = cms.EDProducer("GenVertexProducer",
  src = cms.InputTag('prunedGenParticles'),
  pdgIds = cms.vint32(-15, +15) # CV: use { -15, +15 } for signal, empty list for background                                    
) 
process.analysisSequence += process.genVertex

process.analyzeTracksWrtGenVertex = cms.EDAnalyzer("L1TrackAnalyzer",
  srcGenTaus = cms.InputTag('selectedGenHadTaus'),
  vtxMode = cms.string("genVtx"),
  srcGenVertex_position = cms.InputTag('genVertex:position'),                                                   
  srcOfflineTracks = cms.InputTag('generalTracks'),
  srcOfflinePFCands = cms.InputTag('particleFlow'),
  srcL1Tracks = cms.InputTag('TTTracksFromTracklet:Level1TTTracks'),
  srcL1PFCands = cms.InputTag('l1pfCandidates:PF'),                                                       
  dqmDirectory = cms.string("L1TrackAnalyzerWrtGenVertex"),
  debug = cms.bool(False)                                     
)
process.analysisSequence += process.analyzeTracksWrtGenVertex

for useStrips in [ True, False ]:
    moduleNameBase = "L1HPSPFTauProducer"
    moduleLabel = None
    if useStrips:
        moduleLabel = "WithStrips"
    else:
        moduleLabel = "WithoutStrips"

    moduleNamePF_L1HPSPFTauAnalyzerSignalWrtGenHadTaus = "analyzeL1HPSPFTaus" + moduleLabel + "PFWrtGenHadTaus"
    modulePF_L1HPSPFTauAnalyzerSignalWrtGenHadTaus = cms.EDAnalyzer("L1HPSPFTauAnalyzerSignal",
      srcNumerator = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
      srcDenominator = cms.InputTag('offlineMatchedGenHadTaus'),
      typeDenominator = cms.string("gen"),                                                                            
      dqmDirectory = cms.string("L1HPSPFTauAnalyzerSignal" + moduleLabel + "PF_wrtGenHadTaus")
    )
    setattr(process, moduleNamePF_L1HPSPFTauAnalyzerSignalWrtGenHadTaus, modulePF_L1HPSPFTauAnalyzerSignalWrtGenHadTaus)
    process.analysisSequence += getattr(process, moduleNamePF_L1HPSPFTauAnalyzerSignalWrtGenHadTaus)

    moduleNamePF_L1HPSPFTauPairAnalyzerWrtGenHadTaus = "analyzeL1HPSPFTauPairs" + moduleLabel + "PFWrtGenHadTaus"
    modulePF_L1HPSPFTauPairAnalyzerWrtGenHadTaus = cms.EDAnalyzer("L1HPSPFTauPairAnalyzer",
      srcL1PFTaus = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
      srcRefTaus = cms.InputTag('offlineMatchedGenHadTaus'),
      min_refTau_pt = cms.double(20.),
      max_refTau_pt = cms.double(1.e+3),                                                                
      min_refTau_absEta = cms.double(-1.),
      max_refTau_absEta = cms.double(2.4),                                                                
      dqmDirectory = cms.string("L1HPSPFTauPairAnalyzer" + moduleLabel + "PF_wrtGenHadTaus")
    )
    setattr(process, moduleNamePF_L1HPSPFTauPairAnalyzerWrtGenHadTaus, modulePF_L1HPSPFTauPairAnalyzerWrtGenHadTaus)
    process.analysisSequence += getattr(process, moduleNamePF_L1HPSPFTauPairAnalyzerWrtGenHadTaus)

    moduleNamePF_L1HPSPFTauAnalyzerSignalWrtOfflineTaus = "analyzeL1HPSPFTaus" + moduleLabel + "PFWrtOfflineTaus"
    modulePF_L1HPSPFTauAnalyzerSignalWrtOfflineTaus = cms.EDAnalyzer("L1HPSPFTauAnalyzerSignal",
      srcNumerator = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
      srcDenominator = cms.InputTag('selectedOfflinePFTaus'),
      typeDenominator = cms.string("offline"),                                                                   
      dqmDirectory = cms.string("L1HPSPFTauAnalyzerSignal" + moduleLabel + "PF_wrtOfflineTaus")
    )
    setattr(process, moduleNamePF_L1HPSPFTauAnalyzerSignalWrtOfflineTaus, modulePF_L1HPSPFTauAnalyzerSignalWrtOfflineTaus)
    process.analysisSequence += getattr(process, moduleNamePF_L1HPSPFTauAnalyzerSignalWrtOfflineTaus)

    moduleNamePF_L1HPSPFTauPairAnalyzerWrtOfflineTaus = "analyzeL1HPSPFTauPairs" + moduleLabel + "PFWrtOfflineTaus"
    modulePF_L1HPSPFTauPairAnalyzerWrtOfflineTaus = cms.EDAnalyzer("L1HPSPFTauPairAnalyzer",
      srcL1PFTaus = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
      srcRefTaus = cms.InputTag('selectedOfflinePFTaus'),
      min_refTau_pt = cms.double(20.),
      max_refTau_pt = cms.double(1.e+3),                                                                
      min_refTau_absEta = cms.double(-1.),
      max_refTau_absEta = cms.double(2.4),
      dqmDirectory = cms.string("L1HPSPFTauPairAnalyzer" + moduleLabel + "PF_wrtOfflineTaus")
    )
    setattr(process, moduleNamePF_L1HPSPFTauPairAnalyzerWrtOfflineTaus, modulePF_L1HPSPFTauPairAnalyzerWrtOfflineTaus)
    process.analysisSequence += getattr(process, moduleNamePF_L1HPSPFTauPairAnalyzerWrtOfflineTaus)

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
    outputFileName = cms.string('analyzeL1HPSPFTaus_signal_2020Jul16.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
