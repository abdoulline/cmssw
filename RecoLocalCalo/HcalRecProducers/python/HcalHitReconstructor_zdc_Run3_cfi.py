import FWCore.ParameterSet.Config as cms

zdcreco = cms.EDProducer(
    "ZdcHitReconstructor_Run3",
    correctionPhaseNS = cms.double(0.0),
    digiLabelQIE10ZDC = cms.InputTag("hcalDigis:ZDC"),
    Subdetector = cms.string('ZDC'),
    correctForPhaseContainment = cms.bool(False),
    correctForTimeslew = cms.bool(False),
    dropZSmarkedPassed = cms.bool(True),
    ignoreRPD = cms.bool(True),
    recoMethod = cms.int32(1),
        
    #Tags for calculating status flags
    # None of the flag algorithms have been implemented for zdc, so these booleans do nothing
    correctTiming = cms.bool(True),
    setNoiseFlags = cms.bool(True),
    setHSCPFlags  = cms.bool(True),
    setSaturationFlags = cms.bool(True),
    setTimingTrustFlags = cms.bool(False), # timing flags currently only implemented for HF
    
    saturationParameters=  cms.PSet(maxADCvalue=cms.int32(255))
    ) # zdcreco


