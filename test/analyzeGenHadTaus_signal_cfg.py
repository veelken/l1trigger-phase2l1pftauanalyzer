
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeTallinnL1PFTausSignal")

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

sample = "qqH" # SM VBF Higgs->tautau
#sample = "ggH" #  SM Higgs->tautau produced via gluon fusion 

#--------------------------------------------------------------------------------
# set input files

import os
import re

inputFilePath = None
if sample == "qqH":
    inputFilePath = '/hdfs/cms/store/user/sbhowmik/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack/PhaseIIMTDTDRAutumn18MiniAOD_20190617/190618_084235/0000/'
elif sample == "ggH":
    inputFilePath = '/hdfs/cms/store/user/sbhowmik/GluGluHToTauTau_M125_14TeV_powheg_pythia8/GluGluHToTauTau_PhaseIIMTDTDRAutumn18MiniAOD_20190617/190618_090007/0000/'
else:
    raise ValueError("Invalid sample = '%s' !!" % sample)
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
#process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
process.selectedGenHadTaus = cms.EDFilter("GenJetSelector",
  src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
  # CV: apply relaxed pT and eta cuts to see full pT and eta distributions
  cut = cms.string('pt > 1. & abs(eta) < 3.'),
  filter = cms.bool(False)
)
process.analysisSequence += process.selectedGenHadTaus

#process.selectedGenHadTauFilter = cms.EDFilter("CandViewCountFilter",
#  src = cms.InputTag('selectedGenHadTaus'),
#  minNumber = cms.uint32(2)
#)
#process.analysisSequence += process.selectedGenHadTauFilter

process.genMatchedOfflinePFTaus = cms.EDFilter("PATTauAntiOverlapSelector",
  src = cms.InputTag('slimmedTaus'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(False)                                                          
)
process.analysisSequence += process.genMatchedOfflinePFTaus

process.selectedOfflinePFTaus = cms.EDFilter("PATTauSelector",
  src = cms.InputTag('genMatchedOfflinePFTaus'),
  # CV: apply relaxed pT and eta cuts to see full pT and eta distributions
  cut = cms.string("pt > 1. & abs(eta) < 3. & tauID('decayModeFinding') > 0.5 & tauID('chargedIsoPtSum') < 1.5")
)
process.analysisSequence += process.selectedOfflinePFTaus

#process.selectedOfflinePFTauFilter = cms.EDFilter("CandViewCountFilter",
#  src = cms.InputTag('selectedOfflinePFTaus'),
#  minNumber = cms.uint32(2)
#)
#process.analysisSequence += process.selectedOfflinePFTauFilter

process.offlineMatchedGenHadTaus = cms.EDFilter("GenJetAntiOverlapSelector",
  src = cms.InputTag('selectedGenHadTaus'),
  srcNotToBeFiltered = cms.VInputTag('selectedOfflinePFTaus'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(False)                                                          
)
process.analysisSequence += process.offlineMatchedGenHadTaus

#process.offlineMatchedGenHadTauFilter = cms.EDFilter("CandViewCountFilter",
#  src = cms.InputTag('offlineMatchedGenHadTaus'),
#  minNumber = cms.uint32(2)
#)
#process.analysisSequence += process.offlineMatchedGenHadTauFilter

process.analyzeGenHadTaus = cms.EDAnalyzer("GenHadTauAnalyzer",
  src = cms.InputTag('selectedGenHadTaus'),
  min_pt = cms.double(20.),
  max_pt = cms.double(1.e+3),
  min_absEta = cms.double(-1.),                                          
  max_absEta = cms.double(2.4),
  dqmDirectory = cms.string("GenHadTauAnalyzer")
)
process.analysisSequence += process.analyzeGenHadTaus

process.analyzeGenHadTausOfflineMatched = cms.EDAnalyzer("GenHadTauAnalyzer",
  src = cms.InputTag('offlineMatchedGenHadTaus'),
  min_pt = cms.double(20.),
  max_pt = cms.double(1.e+3),
  min_absEta = cms.double(-1.),                                          
  max_absEta = cms.double(2.4),
  dqmDirectory = cms.string("GenHadTauAnalyzer_offlineMatched")
)
process.analysisSequence += process.analyzeGenHadTausOfflineMatched

process.prunedGenTaus = cms.EDFilter("GenParticleSelector",
  src = cms.InputTag('prunedGenParticles'),
  cut = cms.string("abs(pdgId) = 15"),                         
  stableOnly = cms.bool(False),
  filter = cms.bool(False)
)
process.analysisSequence += process.prunedGenTaus

process.offlineMatchedGenTaus = cms.EDFilter("GenParticleAntiOverlapSelector",
  src = cms.InputTag('prunedGenTaus'),
  srcNotToBeFiltered = cms.VInputTag('selectedOfflinePFTaus'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(True)                                                          
)
process.analysisSequence += process.offlineMatchedGenTaus

#process.offlineMatchedGenTauFilter = cms.EDFilter("CandViewCountFilter",
#  src = cms.InputTag('offlineMatchedGenTaus'),
#  minNumber = cms.uint32(2)
#)
#process.analysisSequence += process.offlineMatchedGenTauFilter

process.analyzeGenTaus = cms.EDAnalyzer("GenTauAnalyzer",
  src = cms.InputTag('prunedGenTaus'),
  min_pt = cms.double(20.),
  max_pt = cms.double(1.e+3),
  min_absEta = cms.double(-1.),                                          
  max_absEta = cms.double(2.4),
  dqmDirectory = cms.string("GenTauAnalyzer")
)
process.analysisSequence += process.analyzeGenTaus

process.analyzeGenTausOfflineMatched = cms.EDAnalyzer("GenTauAnalyzer",
  src = cms.InputTag('offlineMatchedGenTaus'),
  min_pt = cms.double(20.),
  max_pt = cms.double(1.e+3),
  min_absEta = cms.double(-1.),                                          
  max_absEta = cms.double(2.4),
  dqmDirectory = cms.string("GenTauAnalyzer_offlineMatched")
)
process.analysisSequence += process.analyzeGenTausOfflineMatched

process.analyzeOfflineTaus = cms.EDAnalyzer("PATTauAnalyzer",
  src = cms.InputTag('selectedOfflinePFTaus'),
  min_pt = cms.double(20.),
  max_pt = cms.double(1.e+3),
  min_absEta = cms.double(-1.),                                          
  max_absEta = cms.double(2.4),                                             
  dqmDirectory = cms.string("PATTauAnalyzer")
)
process.analysisSequence += process.analyzeOfflineTaus
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('GenHadTauAnalyzer_signal_%s_2019Oct16.root' % sample)
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
