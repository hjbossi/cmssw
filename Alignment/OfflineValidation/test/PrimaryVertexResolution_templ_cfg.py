#! /bin/env cmsRun

'''
cfg to produce pv resolution plots
here doing refit of tracks and vertices using latest alignment 
'''

import FWCore.ParameterSet.Config as cms
from fnmatch import fnmatch
import FWCore.ParameterSet.VarParsing as VarParsing
from pdb import set_trace

process = cms.Process("PrimaryVertexResolution")

###################################################################
# Set the process to run multi-threaded
###################################################################
process.options.numberOfThreads = 8

###################################################################
def best_match(rcd):
###################################################################
    '''
    find out where to best match the input conditions
    '''
    print(rcd)
    for pattern, string in connection_map:
        print(pattern, fnmatch(rcd, pattern))
        if fnmatch(rcd, pattern):
            return string

options = VarParsing.VarParsing("analysis")

options.register('lumi',
                 1.,
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.float,          # string, int, or float
                 "luminosity used")

options.register ('outputRootFile',
                  "pvresolution_YYY_KEY_YYY_XXX_RUN_XXX.root",
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "output root file")

options.register ('records',
                  [],
                  VarParsing.VarParsing.multiplicity.list, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "record:tag names to be used/changed from GT")

options.register ('external',
                  [],
                  VarParsing.VarParsing.multiplicity.list, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "record:fle.db picks the following record from this external file")

options.register ('TrackCollection',
                  'ALCARECOTkAlMinBias',
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "track collection to use")

options.register ('GlobalTag',
                  '110X_dataRun3_Prompt_v3',
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Global Tag to be used")

options.parseArguments()

print("TrackCollection   : ", options.TrackCollection)
print("conditionGT       : ", options.GlobalTag)
print("conditionOverwrite: ", options.records)
print("external conditions:", options.external)
print("outputFile        : ", options.outputRootFile)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr = cms.untracked.PSet(enable = cms.untracked.bool(True)) #False to silence errors
process.MessageLogger.cout = cms.untracked.PSet(INFO = cms.untracked.PSet(
        reportEvery = cms.untracked.int32(1000) # every 100th only
        #    limit = cms.untracked.int32(10)       # or limit to 10 printouts...
    ))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(150000) )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')

process.load('Configuration/StandardSequences/Services_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(XXX_FILES_XXX)
                            )

###################################################################
# Tell the program where to find the conditons
connection_map = [
    ('Tracker*', 'frontier://PromptProd/CMS_CONDITIONS'),
    ('SiPixel*', 'frontier://PromptProd/CMS_CONDITIONS'),
    ('SiStrip*', 'frontier://PromptProd/CMS_CONDITIONS'),
    ('Beam*', 'frontier://PromptProd/CMS_CONDITIONS'),
    ]

if options.external:
    connection_map.extend(
        (i.split(':')[0], 'sqlite_file:%s' % i.split(':')[1]) for i in options.external
        )

connection_map.sort(key=lambda x: -1*len(x[0]))

###################################################################
# creat the map for the GT toGet
records = []
if options.records:
    for record in options.records:
        rcd, tag = tuple(record.split(':'))
        records.append(
            cms.PSet(
                record = cms.string(rcd),
                tag    = cms.string(tag),
                connect = cms.string(best_match(rcd))
                )
            )

process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.GlobalTag.globaltag  = options.GlobalTag
#process.GlobalTag.DumpStat = cms.untracked.bool(True)
process.GlobalTag.toGet = cms.VPSet(*records)

process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
# remove the following lines if you run on RECO files
process.TrackRefitter.src =  options.TrackCollection
process.TrackRefitter.NavigationSchool = ''

####################################################################
# Refitting Sequence
####################################################################
process.load("RecoLocalTracker.SiPixelRecHits.SiPixelTemplateStoreESProducer_cfi")
process.seqTrackselRefit = cms.Sequence(process.offlineBeamSpot*
                                        process.TrackRefitter,
                                        cms.Task(process.SiPixelTemplateStoreESProducer))

####################################################################
## PV refit
####################################################################
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import offlinePrimaryVertices 
process.offlinePrimaryVerticesFromRefittedTrks  = offlinePrimaryVertices.clone()
process.offlinePrimaryVerticesFromRefittedTrks.TrackLabel                                       = cms.InputTag("TrackRefitter") 
process.offlinePrimaryVerticesFromRefittedTrks.vertexCollections.maxDistanceToBeam              = 1
process.offlinePrimaryVerticesFromRefittedTrks.TkFilterParameters.maxNormalizedChi2             = 20
process.offlinePrimaryVerticesFromRefittedTrks.TkFilterParameters.minSiliconLayersWithHits      = 5
process.offlinePrimaryVerticesFromRefittedTrks.TkFilterParameters.maxD0Significance             = 5.0
# as it was prior to https://github.com/cms-sw/cmssw/commit/c8462ae4313b6be3bbce36e45373aa6e87253c59
process.offlinePrimaryVerticesFromRefittedTrks.TkFilterParameters.maxD0Error                    = 1.0
process.offlinePrimaryVerticesFromRefittedTrks.TkFilterParameters.maxDzError                    = 1.0
process.offlinePrimaryVerticesFromRefittedTrks.TkFilterParameters.minPixelLayersWithHits        = 2   

process.PrimaryVertexResolution = cms.EDAnalyzer('SplitVertexResolution',
                                                 storeNtuple         = cms.bool(False),
                                                 intLumi             = cms.untracked.double(options.lumi),
                                                 vtxCollection       = cms.InputTag("offlinePrimaryVerticesFromRefittedTrks"),
                                                 trackCollection     = cms.InputTag("TrackRefitter"),		
                                                 minVertexNdf        = cms.untracked.double(10.),
                                                 minVertexMeanWeight = cms.untracked.double(0.5),
                                                 runControl = cms.untracked.bool(True),
                                                 runControlNumber = cms.untracked.vuint32(int(XXX_RUN_XXX)),
                                                 sumpTStartScale = cms.untracked.double(1.),
                                                 sumpTEndScale = cms.untracked.double(1000.),
                                                 nTrackBins = cms.untracked.double(60.),
                                                 nVtxBins = cms.untracked.double(40.)
                                                 )

###################################################################
# TFileService
###################################################################
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outputRootFile),	
                                   closeFileFast = cms.untracked.bool(False)
                                   )

###################################################################
# Beamspot compatibility check
###################################################################
from RecoVertex.BeamSpotProducer.beamSpotCompatibilityChecker_cfi import beamSpotCompatibilityChecker
process.BeamSpotChecker = beamSpotCompatibilityChecker.clone(
    bsFromFile = "offlineBeamSpot::RECO",  # source of the event beamspot (in the ALCARECO files)
    bsFromDB = "offlineBeamSpot::@currentProcess", # source of the DB beamspot (from Global Tag) NOTE: only if dbFromEvent is True!
    dbFromEvent = True,
    warningThr = 3, # significance threshold to emit a warning message
    errorThr = 5    # significance threshold to abort the job
)

###################################################################
# Path
###################################################################
process.p = cms.Path(process.seqTrackselRefit                       +
                     #process.offlineBeamSpot                       +
                     #process.TrackRefitter                         +
                     process.BeamSpotChecker                        +
                     process.offlinePrimaryVerticesFromRefittedTrks +
                     process.PrimaryVertexResolution)
