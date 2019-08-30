import FWCore.ParameterSet.Config as cms
from RecoCTPPS.TotemRPLocal.ctppsDiamondLocalReconstruction_cff import *

process = cms.Process("DIAMONDPRODUCER")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/group/dpg_ctpps/comm_ctpps/CalibrationDevel/10B9DB48-1EB4-E811-956C-FA163ECAADB8.root'
    )
)

process.DiamondTimingProducer = cms.EDProducer('DiamondTimingProducer',
    #tagDigi = cms.InputTag("ctppsDiamondRawToDigi", "TimingDiamond"),
    #tagRecHit = cms.InputTag("ctppsDiamondRecHits"),
    tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracks"),
    #tagLocalTrack = cms.InputTag("ctppsDiamondLocalTracks"),
    Ntracks_Lcuts = cms.vint32([-1,1,-1,1]),
    Ntracks_Ucuts = cms.vint32([-1,6,-1,6])
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('produced.root'),
    outputCommands = cms.untracked.vstring('drop *',
      "keep *_ctppsPixelLocalTracks_*_*",
      "keep *_ctppsDiamondRecHits_*_*",
      "keep *_ctppsDiamondLocalTracks_*_*",
      "keep *_ctppsDiamondRawToDigi_*_*",
      "keep *_DiamondTimingProducer_PixelMuxmap_*",
      "keep *_DiamondTimingProducer_SectorTBA_*")
)

process.p = cms.Path(process.DiamondTimingProducer)

process.e = cms.EndPath(process.out)
