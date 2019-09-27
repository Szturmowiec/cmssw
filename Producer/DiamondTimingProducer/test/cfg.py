import FWCore.ParameterSet.Config as cms
from RecoCTPPS.TotemRPLocal.ctppsDiamondLocalReconstruction_cff import *

from Configuration.StandardSequences.Eras import eras
process = cms.Process("DIAMONDPRODUCER",eras.Run2_2018)

    # import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'/store/data/Run2018D/ZeroBias/RAW/v1/000/324/021/00000/CED9C85E-7F4B-C545-BACC-9ECF0459EC04.root',
'/store/data/Run2018D/ZeroBias/RAW/v1/000/324/021/00000/47B83E15-4E53-A546-B34B-9708D5A2CDC6.root',
'/store/data/Run2018D/ZeroBias/RAW/v1/000/324/021/00000/2C1FD412-7B58-B449-922B-DA0A4CC66E02.root',
'/store/data/Run2018D/ZeroBias/RAW/v1/000/324/021/00000/8A77F857-0022-FE4B-BD3E-1382A082868F.root'
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/324/021/00000/7665DF69-4838-4447-A55E-C7714303E359.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/324/021/00000/619587B5-B7F1-414E-A9D1-BDE9DC1397E2.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/324/021/00000/22B32972-F03B-E048-896B-C9FDC035F850.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/324/021/00000/54D65CF5-D04B-6B46-A058-4CDBE28FAAD0.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/324/021/00000/2A637F3F-212E-EB40-898F-299E5345CEE7.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/324/021/00000/AF7CEED2-CAAE-B54B-9CA9-5F605CFCCB3E.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/324/021/00000/50A2EC60-3F07-6844-8C52-4ED8114C0420.root'
    #'/store/group/dpg_ctpps/comm_ctpps/CalibrationDevel/10B9DB48-1EB4-E811-956C-FA163ECAADB8.root'
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/FE86FF5E-99B1-E811-B512-FA163EBBBA60.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/FE4DDC18-87B1-E811-BA58-FA163E7F2CAB.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/FCC303F6-7CB1-E811-8D9D-FA163ECA48A6.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/FCA883E4-8DB1-E811-92C2-FA163E62F8E9.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/FC4674E1-84B1-E811-8CA5-FA163EB23546.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/FC2BC92C-96B1-E811-87AD-FA163EA65232.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/FAA38E94-7FB1-E811-BB21-FA163EBFA384.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/FA3C163E-7EB1-E811-914F-FA163E29065E.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/FA36B864-ACB1-E811-B311-FA163E20F796.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/FA27BEE3-B3B1-E811-A054-FA163E2415B4.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/FA154BCE-ADB1-E811-9978-FA163EA15274.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F848293E-ACB1-E811-8645-FA163E4200C7.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F69981E5-7BB1-E811-883F-FA163ED92461.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F4FD0C38-82B1-E811-ABDD-FA163E62F8E9.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F4A427BC-82B1-E811-8D9F-FA163E3F4095.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F49BE245-97B1-E811-9951-FA163E45F5DC.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F47358B1-83B1-E811-B368-FA163EAF1239.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F419DEB7-8DB1-E811-B17D-FA163E72BD86.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F4100073-8CB1-E811-A09B-FA163E49433D.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F291CDF2-97B1-E811-983F-02163E019EF4.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F0C1C6E8-8EB1-E811-8533-FA163E91FAAB.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F0BA6057-92B1-E811-A730-FA163ECA48A6.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F099AE93-8DB1-E811-BD9C-FA163E5F28A9.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F070336A-88B1-E811-8DA6-FA163E8A2D90.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F05B53A7-9AB1-E811-8433-FA163EA2499B.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/F015D799-92B1-E811-AD89-FA163E55495F.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/EE7577A3-ADB1-E811-A7AF-FA163E995805.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/EE21A2A4-94B1-E811-9189-FA163EA53599.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/EC929746-ACB1-E811-AAA6-FA163E97BB97.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/EAD705F7-7CB1-E811-9FF0-FA163EF96F9B.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/E8EA8135-ACB1-E811-82CB-FA163E9875C8.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/E8C0F542-87B1-E811-BC81-FA163E42743B.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/E8A607C1-84B1-E811-98D1-FA163E2B34AA.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/E6085C9D-94B1-E811-A22F-FA163E9681D4.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/E4CE5CF4-87B1-E811-953D-FA163EE65050.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/E297F714-9CB1-E811-A034-FA163E6A1015.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/E2319F62-8FB1-E811-B0BC-FA163E679CB5.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/E2135A53-ACB1-E811-9164-FA163E43FD94.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/E0D4ABF9-7CB1-E811-8B92-FA163E0303D4.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/DEC578CA-81B1-E811-A2A8-FA163EBCA265.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/DE2D2082-7EB1-E811-8BE7-02163E014945.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/DC8FBE28-89B1-E811-9F29-FA163EBB8D74.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/DAE9DC93-7FB1-E811-9DD0-02163E01A0A2.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/DAE69A64-99B1-E811-AD36-FA163E17C97B.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/DA8FD654-BAB1-E811-BE19-FA163E417864.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/D8E854F0-99B1-E811-9C5B-FA163EB1CCFA.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/D8A1372A-7EB1-E811-9D32-FA163EFAEDB1.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/D811C0FA-9BB1-E811-B167-A4BF01277792.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/D6F3089F-A6B1-E811-BC76-FA163EA8758B.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/D6DA6AEF-99B1-E811-967D-FA163E2E4B1A.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/D684CCD7-91B1-E811-ACBD-FA163E52EF36.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/D6732A08-9CB1-E811-A945-02163E01A14D.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/D263927D-AFB1-E811-A3F8-FA163E228A06.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/D21DA20D-8DB1-E811-A182-FA163E86F545.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/D0A75554-83B1-E811-AAE2-FA163E0683BD.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/D07B9883-92B1-E811-BDE3-02163E019E95.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/CECF0190-98B1-E811-AB91-FA163E25A55A.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/CE4874B9-88B1-E811-B271-FA163E20727F.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/CCF6012A-8EB1-E811-8E75-FA163E0303D4.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/CA388B2F-A6B1-E811-99D7-02163E019F0E.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/C8BBA6E8-82B1-E811-83DB-02163E010C3C.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/C6F04069-80B1-E811-AB17-FA163E116871.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/C6568727-8EB1-E811-AC5A-FA163E0CC6EF.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/C41F05D0-96B1-E811-BEFC-FA163E6DC2D0.root',
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/322/332/00000/C28DCBF4-8AB1-E811-A88F-FA163EE81CA0.root'
    )
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '110X_dataRun2_v3', '')

# raw-to-digi conversion
process.load("EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff")

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")
process.totemDAQMappingESSourceXML_TimingDiamond = cms.ESSource("TotemDAQMappingESSourceXML",
           verbosity = cms.untracked.uint32(0),
           subSystem = cms.untracked.string("TimingDiamond"),
           configuration = cms.VPSet(
                 cms.PSet(
                         validityRange=cms.EventRange("1:min -999999999:max"),
                         mappingFileNames=cms.vstring("CondFormats/CTPPSReadoutObjects/xml/mapping_timing_diamond_2021_fake.xml"),
                         maskFileNames=cms.vstring()
)
)
)
process.load('Geometry.VeryForwardGeometry.geometryRPFromDD_2021_cfi')
process.ctppsDiamondRecHits.applyCalibration = False
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

process.path = cms.Path(
    process.ctppsRawToDigi *
    process.recoCTPPS *
    process.DiamondTimingProducer
)

process.end_path = cms.EndPath(process.out)

process.schedule = cms.Schedule(
    process.path,
    process.end_path
)
