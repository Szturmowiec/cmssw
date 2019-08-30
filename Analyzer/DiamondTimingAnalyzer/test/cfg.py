import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from RecoCTPPS.TotemRPLocal.ctppsDiamondLocalReconstruction_cff import *
import string

process = cms.Process("TIMINGSTUDY")
options = VarParsing ('analysis')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)


options.register ('calibFile',
				  0,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "Calibration input file")

options.register ('geometryFile',
				  0,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "Geometry input file")

options.register ('validOOT',
				  0,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "valid OOT slice")

options.register ('debug',
				  0,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "debug option")

options.parseArguments()


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/p/pcwiklic/Project/Modular/CMSSW_11_0_0_pre4/src/produced.root'
        #'file:/eos/cms/store/group/dpg_ctpps/comm_ctpps/CalibrationDevel/10B9DB48-1EB4-E811-956C-FA163ECAADB8.root',
        )
        #secondaryFileNames= cms.untracked.vstring(
        #        'file:/afs/cern.ch/user/p/pcwiklic/Project/Modular/CMSSW_11_0_0_pre4/src/produced.root'
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# rechits production
process.load(options.geometryFile)
process.load('RecoCTPPS.TotemRPLocal.ctppsDiamondLocalReconstruction_cff')

process.DiamondTimingAnalyzer = cms.EDAnalyzer('DiamondTimingAnalyzer',
 tagPixelMuxMap = cms.InputTag("DiamondTimingProducer","PixelMuxmap"),
 tagSectorTBA = cms.InputTag("DiamondTimingProducer","SectorTBA"),
 #tagDigi = cms.InputTag("ctppsDiamondRawToDigi", "TimingDiamond"),
 tagRecHit = cms.InputTag("ctppsDiamondRecHits"),
 tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracks"),
 tagLocalTrack = cms.InputTag("ctppsDiamondLocalTracks"),
 tagCalibrationFile = cms.string(options.calibFile),
 tagValidOOT = cms.int32(options.validOOT),
 debug = cms.int32(options.debug),
)

process.DiamondTimingProducer = cms.EDProducer('DiamondTimingProducer',
    #tagDigi = cms.InputTag("ctppsDiamondRawToDigi", "TimingDiamond"),
    #tagLocalTrack = cms.InputTag("ctppsDiamondLocalTracks"),
    #tagRecHit = cms.InputTag("ctppsDiamondRecHits"),
    tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracks"),
    Ntracks_Lcuts = cms.vint32([-1,1,-1,1]),
    Ntracks_Ucuts = cms.vint32([-1,6,-1,6])
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string(options.outputFile)
)


process.path1 = cms.Path(
   ctppsDiamondLocalReconstruction
)

process.path2 = cms.Path(
   process.DiamondTimingProducer
)

process.path3 = cms.Path(
   process.DiamondTimingAnalyzer
)



process.schedule = cms.Schedule(
#process.path
process.path1,
#process.path2,
process.path3
)
