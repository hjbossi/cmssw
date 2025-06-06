import FWCore.ParameterSet.Config as cms

from L1Trigger.TrackFindingTracklet.Producer_cfi import TrackFindingTrackletProducer_params
from L1Trigger.TrackFindingTracklet.ChannelAssignment_cff import ChannelAssignment
from L1Trigger.TrackerTFP.TrackQuality_cff import *
from L1Trigger.TrackerTFP.LayerEncoding_cff import TrackTriggerLayerEncoding

l1tTTTracksFromTrackletEmulation = cms.EDProducer("L1FPGATrackProducer",
                                               TrackFindingTrackletProducer_params,
                                               TTStubSource = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
                                               InputTagTTDTC = cms.InputTag("ProducerDTC", "StubAccepted"),
                                               readMoreMcTruth = cms.bool(False),
                                               MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
                                               MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
                                               TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                               BeamSpotSource = cms.InputTag("offlineBeamSpot"),
                                               asciiFileName = cms.untracked.string(""),
                                               FailScenario = cms.untracked.int32(0),
                                               Extended = cms.bool(False),
                                               Reduced = cms.bool(False),
                                               Hnpar = cms.uint32(4),
                                               # These 3 files only used for extended or reduced mode.
                                               memoryModulesFile = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/memorymodules_hourglassExtendedAllCombined.dat'),
                                               processingModulesFile = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/processingmodules_hourglassExtendedAllCombined.dat'),
                                               wiresFile = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/wires_hourglassExtendedAllCombined.dat'),
                                               # Quality Flag and Quality params
                                               TrackQuality = cms.bool(True),
                                               Fakefit = cms.bool(False), # True causes Tracklet reco to output TTTracks before DR & KF
                                               StoreTrackBuilderOutput = cms.bool(False), # if True EDProducts for TrackBuilder tracks and stubs will be filled
                                               RemovalType = cms.string("merge"), # Duplicate track removal
                                               DoMultipleMatches = cms.bool(True) # Allow tracklet tracks multiple stubs per layer
    )

l1tTTTracksFromExtendedTrackletEmulation = l1tTTTracksFromTrackletEmulation.clone(
                                               Extended = cms.bool(True),
                                               Reduced = cms.bool(False),
                                               Hnpar = cms.uint32(5),
                                               # specifying where the TrackletEngineDisplaced(TED)/TripletEngine(TRE) tables are located
                                               tableTEDFile = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/table_TED/table_TED_D1PHIA1_D2PHIA1.txt'),
                                               tableTREFile = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/table_TRE/table_TRE_D1AD2A_1.txt'),
                                               # Quality Flag and Quality params
                                               TrackQuality = cms.bool(False)
    )
