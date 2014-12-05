import FWCore.ParameterSet.Config as cms

process = cms.Process("Tree")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out0459FA4C-E16B-E411-A5F8-02163E0103BD.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out0475675C-DF6B-E411-B9CF-02163E00E8CB.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out0644761D-E46B-E411-A7B7-02163E00EB2B.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out0A689926-E06B-E411-9B97-001E67AC06E0.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out0C70FBD2-E06B-E411-9173-02163E00E886.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out0C888F70-EC6B-E411-8A38-02163E00FFB5.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out0E558B37-DD6B-E411-AAB0-02163E00F46B.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out12B258AE-E26B-E411-BB49-02163E00F355.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out1618E3ED-E56B-E411-AE29-02163E010340.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out1A122B72-E96B-E411-8E1B-0025904B201E.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out1CC92F5D-E16B-E411-9307-001E675053A5.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out2A186B2B-E46B-E411-9428-02163E00F377.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out2A2E9E6C-E16B-E411-8343-001E67ABF6FC.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out36193A3A-E06B-E411-AEF1-003048FEB8FE.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out40DC3A54-EE6B-E411-BF56-02163E0103BD.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out4CB71478-E16B-E411-875D-02163E00E5EE.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out4E538362-DB6B-E411-B254-02163E00EAF3.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out5275D24D-DF6B-E411-97E6-0025B3244016.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out64230499-E36B-E411-99B9-02163E00E8F7.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out74D83D88-E16B-E411-B8C1-02163E00E7B7.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out7C0BD880-E16B-E411-984D-02163E00B752.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out82AA88F0-EE6B-E411-800E-02163E00E968.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out96EF7753-DF6B-E411-999B-0025904B26A8.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/Out9C30CC4E-DF6B-E411-BD95-0025B320378A.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutA8D5BF81-E46B-E411-B2E9-02163E00F4C2.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutB00CDC41-DD6B-E411-8466-02163E00E5A6.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutB0D7CE09-E96B-E411-B990-02163E00E86F.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutB45CBD30-E16B-E411-8F5A-002590495076.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutB6619371-E16B-E411-8C2E-003048FEBA12.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutB8904E27-FE6B-E411-AA23-0025901D627E.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutBC5E9A9E-E66B-E411-ABBC-02163E00E9CE.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutC432C720-DE6B-E411-BB1F-02163E00E96E.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutD09AE332-E16B-E411-97E1-02163E00F3BF.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutD654F076-E46B-E411-A06A-0025904B1FC0.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutD6FADC12-E86B-E411-AD90-02163E00E8AE.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutD855698A-E16B-E411-8E4A-02163E00EB6E.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutD8EC5EA2-E26B-E411-9BC8-00215AEDFFA6.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutDACD381B-CB6B-E411-97BD-FA163EC657C1.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutDC4392F6-EB6B-E411-9692-02163E00BC02.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutDE636D72-E16B-E411-B7F4-02163E00E6D1.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutE4AC81C6-DE6B-E411-85C2-0025904B1FB8.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutEAD788B2-E26B-E411-8107-002590494DD2.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutF4B3F896-DF6B-E411-BC64-00215AEDFD8A.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutF6FBFF84-DC6B-E411-AF9B-0025B3203616.root',
'root://eoscms//eos/cms/store/group/phys_jetmet/pharris/Hcal2/config2correct.py/TT_PU25ns_PRE_LS172_V15-v1/OutFAD849E5-DD6B-E411-B588-02163E00E81B.root',
    )
)

process.tree = cms.EDAnalyzer('SimpleJetResponseTreeMaker',  
  pvSrc = cms.InputTag('offlinePrimaryVertices'),
  rhoSrc = cms.InputTag('fixedGridRhoAll'),#'ak4PFJets', 'rho'),
  caloRhoSrc = cms.InputTag('fixedGridRhoFastjetAllCalo'),#'ak4PFJets', 'rho'),
  caloJetSrc = cms.InputTag('ak4CaloJets'),
  recoJetSrc = cms.InputTag('ak4PFJets'),
  genJetSrc = cms.InputTag('ak4GenJets'),
  genJetMatchDeltaR  = cms.double(0.2),
  genJetPtThreshold  = cms.double(2.0),
  recoJetPtThreshold = cms.double(2.0),
  jecPayloads = cms.vstring([
    'CSA14_V4_MC_L1FastJet_AK4PF.txt',
    'CSA14_V4_MC_L2L3Residual_AK4PF.txt',
    'CSA14_V4_MC_L2Relative_AK4PF.txt',
    'CSA14_V4_MC_L3Absolute_AK4PF.txt',
    'CSA14_V4_MC_Uncertainty_AK4PF.txt'  
  ])
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("/tmp/jdolen/tt25_Hcal2_config2correct.root"),
      closeFileFast = cms.untracked.bool(True)
  )

process.p = cms.Path(process.tree)
