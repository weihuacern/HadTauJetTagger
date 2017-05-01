import FWCore.ParameterSet.Config as cms

prodPFCands = cms.EDFilter(
  "prodPFCands",
  debug  = cms.bool(False),
  pfCandSrc = cms.InputTag("packedPFCandidates"),
  #genDecayLVec = cms.InputTag("prodGenInfo:genDecayLVec"),
)

