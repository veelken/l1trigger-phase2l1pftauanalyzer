
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
#
#inputFilePath = '/store/user/veelken/'
#inputFile_regex = r"[a-zA-Z0-9_/:.-]*trigPATtuple_[a-zA-Z0-9-_]+.root"
#
# check if name of inputFile matches regular expression
#inputFileNames = []
#files = [ "".join([ "file:", inputFilePath, file ]) for file in os.listdir(inputFilePath) ]
#for file in files:
#    #print "file = %s" % file
#    inputFile_matcher = re.compile(inputFile_regex)
#    if inputFile_matcher.match(file):
#        inputFileNames.append(file)
#print "inputFileNames = %s" % inputFileNames 
#
#process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
process.analyzeTallinL1PFTaus = cms.EDAnalyzer("TallinnL1PFTauAnalyzerBackground",
  srcTallinnL1PFTaus = cms.InputTag('TallinnL1PFTauProducer:L1PFTaus'),
  dqmDirectory = cms.string("TallinnL1PFTauAnalyzerBackground")
)
process.analysisSequence += process.analyzeTallinL1PFTaus
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('TallinnL1PFTauAnalyzer_background_2019May02.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

#process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True)
#)
