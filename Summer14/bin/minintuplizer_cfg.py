import FWCore.ParameterSet.Config as cms

process = cms.Process("MiniNtuplizer")

process.Options = cms.PSet(
    maxEvents       = cms.int32(100),
    jetR            = cms.double(0.8),
    doCMSSWJets     = cms.bool(False),
    puppiConfig     = cms.string("Puppi_cff.py"),
    L1FastJetJEC    = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS162_V5_L1FastJet_AK7PF.txt"),
    L2RelativeJEC   = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS162_V5_L2Relative_AK7PF.txt"),
    L3AbsoluteJEC   = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS162_V5_L3Absolute_AK7PF.txt"),
    L2L3ResidualJEC = cms.string(""), 
    JECUncertainty  = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS162_V5_Uncertainty_AK7PF.txt")
)
