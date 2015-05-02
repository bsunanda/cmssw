import FWCore.ParameterSet.Config as cms

RecAnalyzerMinbias = cms.EDAnalyzer("RecAnalyzerMinbias",
                                    HistOutFile = cms.untracked.string('analysis_minbias.root'),
                                    hbheInputMB = cms.InputTag("hbherecoMB"),
                                    hfInputMB   = cms.InputTag("hfrecoMB"),
                                    Recalib     = cms.bool(False),
                                    IgnoreL1    = cms.untracked.bool(False),
                                    RunNZS      = cms.bool(True),
                                    ELowHB      = cms.double(4),
                                    EHighHB     = cms.double(100),
                                    ELowHE      = cms.double(4),
                                    EHighHE     = cms.double(150),
                                    ELowHF      = cms.double(10),
                                    EHighHF     = cms.double(150),
                                    )
