import FWCore.ParameterSet.Config as cms

puppiCentral = cms.VPSet(
                 cms.PSet(
                  algoId           = cms.untracked.int32(5),  #0 is default Puppi
                  useCharged       = cms.untracked.bool(True),
                  applyLowPUCorr   = cms.untracked.bool(True),
                  combOpt          = cms.untracked.int32(0),
                  cone             = cms.untracked.double(0.2),
                  rmsPtMin         = cms.untracked.double(0.5),
                  rmsScaleFactor   = cms.untracked.double(1.0)
                 )
                )

puppiForward = cms.VPSet(
                cms.PSet(
                 algoId         = cms.untracked.int32(0),  #0 is default Puppi
                 useCharged     = cms.untracked.bool(False),
                 applyLowPUCorr = cms.untracked.bool(True),
                 combOpt        = cms.untracked.int32(0),
                 cone           = cms.untracked.double(0.2),
                 rmsPtMin       = cms.untracked.double(0.5),
                 rmsScaleFactor = cms.untracked.double(1.0)
                 ),
                cms.PSet(
                 algoId         = cms.untracked.int32(1),  #0 is default Puppi
                 useCharged     = cms.untracked.bool(False),
                 applyLowPUCorr = cms.untracked.bool(True),
                 combOpt        = cms.untracked.int32(0),
                 cone           = cms.untracked.double(0.2),
                 rmsPtMin       = cms.untracked.double(0.5),
                 rmsScaleFactor = cms.untracked.double(1.0)
                 )
                )

puppi = cms.PSet(#"PuppiProducer",
                       PuppiName      = cms.untracked.string("Puppi"),
                       UseDeltaZCut   = cms.untracked.bool  (True),
                       DeltaZCut      = cms.untracked.double(0.3),
                       candName       = cms.untracked.string('particleFlow'),
                       vertexName     = cms.untracked.string('offlinePrimaryVertices'),
                       applyCHS       = cms.untracked.bool  (True),
                       useExp         = cms.untracked.bool  (False),
                       MinNeutralPt   = cms.untracked.double(0.05),
                       MinPuppiWeight = cms.untracked.double(0.01),
                       algos          = cms.VPSet( 
                        cms.PSet( 
                         etaMin = cms.untracked.double(-2.5),
                         etaMax = cms.untracked.double( 2.5),
                         ptMin  = cms.untracked.double(-1.),
                         puppiAlgos = puppiCentral
                        ),
                        cms.PSet( 
                         etaMin = cms.untracked.double(-10.),
                         etaMax = cms.untracked.double(-2.5),
                         ptMin  = cms.untracked.double(-1.),
                         puppiAlgos = puppiForward
                        ),
                        cms.PSet( 
                         etaMin = cms.untracked.double(2.5),
                         etaMax = cms.untracked.double(10.0),
                         ptMin  = cms.untracked.double(-1.),
                         puppiAlgos = puppiForward
                        )
                       )
                      )
