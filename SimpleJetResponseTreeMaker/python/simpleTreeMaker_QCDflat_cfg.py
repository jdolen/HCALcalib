import FWCore.ParameterSet.Config as cms

process = cms.Process("Tree")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_2_0_pre8/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v1/00000/00EBC831-0354-E411-B30E-0025905A60D0.root',
      'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_2_0_pre8/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v1/00000/06A04FBB-0554-E411-929B-003048FFD76E.root',
      'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_2_0_pre8/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v1/00000/14EFB9FC-0154-E411-825B-0025905964C0.root',
      'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_2_0_pre8/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v1/00000/1EC900FB-0054-E411-94AB-0025905B85EE.root',
      'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_2_0_pre8/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v1/00000/34D1538D-0254-E411-AC3A-0025905A48BC.root',
      'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_2_0_pre8/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v1/00000/3A04D180-0854-E411-B38D-0025905A60F8.root',
      'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_2_0_pre8/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v1/00000/40179BD3-FF53-E411-8BD0-0025905964B4.root',
      'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_2_0_pre8/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v1/00000/74B5AE1A-0054-E411-80DC-0025905A6056.root',
      'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_2_0_pre8/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v1/00000/9AECAEFC-0154-E411-A3DC-0025905A6118.root',
      'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_2_0_pre8/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v1/00000/9CE31B67-0454-E411-A0B6-0025905A6088.root',
      'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_2_0_pre8/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v1/00000/A48ECF8E-0254-E411-9DA2-00248C55CC62.root'
    )
)

process.tree = cms.EDAnalyzer('SimpleJetResponseTreeMaker',  
  pvSrc = cms.InputTag('offlinePrimaryVertices'),
  rhoSrc = cms.InputTag('fixedGridRhoAll'),#'ak4PFJets', 'rho'),
  recoJetSrc = cms.InputTag('ak4PFJets'),
  genJetSrc = cms.InputTag('ak5GenJets'),
  jecPayloads = cms.vstring([
    'CSA14_V4_MC_L1FastJet_AK4PF.txt',
    'CSA14_V4_MC_L2L3Residual_AK4PF.txt',
    'CSA14_V4_MC_L2Relative_AK4PF.txt',
    'CSA14_V4_MC_L3Absolute_AK4PF.txt',
    'CSA14_V4_MC_Uncertainty_AK4PF.txt'  
  ])
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("OutTreeQCDflat.root"),
      closeFileFast = cms.untracked.bool(True)
  )

process.p = cms.Path(process.tree)
