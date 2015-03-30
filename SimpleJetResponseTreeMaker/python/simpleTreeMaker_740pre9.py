import FWCore.ParameterSet.Config as cms

process = cms.Process("Tree")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/00847E71-97D1-E411-A46E-002618FDA259.root',
'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/08BCA6E5-9CD1-E411-BF2D-003048FFD75C.root',
'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/0C44E9EF-9BD1-E411-9192-0025905A6132.root',
'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/1CCF8EDA-9DD1-E411-999E-0025905A6136.root',
'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/1E44A3A6-AAD1-E411-95E6-002354EF3BDF.root',
'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/30467D71-97D1-E411-A8CC-0025905B8596.root',
'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/3C300513-9BD1-E411-BCDD-002618943964.root',
'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/76A010AA-9ED1-E411-B2A6-002618943961.root',
'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/80880103-9AD1-E411-83E3-0025905A60DE.root',
'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/8AF6D1A4-AAD1-E411-B430-002618943908.root',
'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/B6357DE1-98D1-E411-8CA2-0025905A60D6.root',
'root://cmsxrootd.fnal.gov///store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/FC3C5EAD-9ED1-E411-8537-0025905A60AA.root', 
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
      fileName = cms.string("/uscmst1b_scratch/lpc1/lpcphys/jdolen/JMAR/740pre9MassResolution/qcdflat.root"),
      closeFileFast = cms.untracked.bool(True)
  )

process.p = cms.Path(process.tree)
