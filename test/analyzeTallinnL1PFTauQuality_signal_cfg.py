
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
        moduleLabel_genMatched = "%sGenMatched" % moduleLabel
        module_genMatched = cms.EDFilter("TallinnL1PFTauAntiOverlapSelector",
          src = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
          srcNotToBeFiltered = cms.VInputTag('offlineMatchedGenHadTaus'),
          dRmin = cms.double(0.3),
          invert = cms.bool(True),
          filter = cms.bool(False)                                                          
        )
        setattr(process, moduleNameBase + moduleLabel_genMatched + "PF", module_genMatched)
        process.analysisSequence += getattr(process, moduleNameBase + moduleLabel_genMatched + "PF")
        
        moduleNamePF_TallinnL1PFTauQualityAnalyzer = "analyzeTallinnL1PFTauQuality" + moduleLabel + "PF"
        modulePF_TallinnL1PFTauQualityAnalyzer = cms.EDAnalyzer("TallinnL1PFTauQualityAnalyzer",
          src = cms.InputTag(moduleNameBase + moduleLabel_genMatched + "PF"),
          dqmDirectory = cms.string("TallinnL1PFTauQualityAnalyzer" + moduleLabel + "PF")
        )
        setattr(process, moduleNamePF_TallinnL1PFTauQualityAnalyzer, modulePF_TallinnL1PFTauQualityAnalyzer)
        process.analysisSequence += getattr(process, moduleNamePF_TallinnL1PFTauQualityAnalyzer)
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('TallinnL1PFTauQualityAnalyzer_signal_%s_2019Nov20.root' % sample)
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
