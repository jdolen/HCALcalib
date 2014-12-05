// -*- C++ -*-
//
// Package:    HCALcalib/SimpleJetResponseTreeMaker
// Class:      SimpleJetResponseTreeMaker
// 
/**\class SimpleJetResponseTreeMaker SimpleJetResponseTreeMaker.cc HCALcalib/SimpleJetResponseTreeMaker/plugins/SimpleJetResponseTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  James Dolen
//         Created:  Fri, 24 Oct 2014 23:27:05 GMT
//
//

// system include files
#include <memory>

// core
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// DataFormats
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h" 
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
// JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// utilities
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

// TFile
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// root
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"


class SimpleJetResponseTreeMaker : public edm::EDAnalyzer {
   public:
      explicit SimpleJetResponseTreeMaker(const edm::ParameterSet&);
      ~SimpleJetResponseTreeMaker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::InputTag   pvSrc_;
      edm::InputTag   rhoSrc_;
      edm::InputTag   caloRhoSrc_;
      edm::InputTag   caloJetSrc_;
      edm::InputTag   recoJetSrc_;
      edm::InputTag   genJetSrc_;
      std::vector<std::string>  jecPayloads_; /// files for JEC payloads
      double genJetMatchDeltaR_   ;
      double genJetPtThreshold_   ;
      double recoJetPtThreshold_  ;

      boost::shared_ptr<FactorizedJetCorrector> jec_;
      boost::shared_ptr<JetCorrectionUncertainty> jecUnc_;

      // Tree filled once per event
      TTree *EventTree;
      Int_t     Event_Run                        ;
      Int_t     Event_Lumi                       ;
      Int_t     Event_Event                      ;
      Int_t     Event_Nvtx                       ;
      Int_t     Event_NvtxGood                   ;
      Float_t   Event_Rho                        ;
      Float_t   Event_RhoCalo                    ;
      Float_t   Event_RhoCaloCentral             ;
      Int_t     Event_NPU_BXm3                   ; 
      Int_t     Event_NPU_BXm2                   ; 
      Int_t     Event_NPU_BXm1                   ; 
      Int_t     Event_NPU_BX0                    ; 
      Int_t     Event_NPU_BX1                    ; 
      Int_t     Event_NPU_BX2                    ; 
      Float_t   Event_CorrSumJetHadE             ;
      Float_t   Event_SumJetPt                   ;
      Float_t   Event_SumJetPt_Pt10              ;
      Float_t   Event_SumJetPt_Pt20              ;
      Float_t   Event_SumJetPt_Pt10Eta2p4        ;
      Float_t   Event_SumJetPt_Barrel            ;
      Float_t   Event_SumJetPt_Endcap            ;

      Int_t     Event_Njets                      ;
      Int_t     Event_Njets_Pt10                 ;
      Int_t     Event_Njets_Pt20                 ;
      Int_t     Event_Njets_Pt10Eta2p4           ;
      Float_t   Event_CorrSumJetPt               ;
      Float_t   Event_CorrSumJetPt_Pt10          ;
      Float_t   Event_CorrSumJetPt_Pt20          ;
      Float_t   Event_CorrSumJetPt_Pt10Eta2p4    ;
      Int_t     Event_CorrNjets                  ;
      Int_t     Event_CorrNjets_Pt10             ;
      Int_t     Event_CorrNjets_Pt20             ;
      Int_t     Event_CorrNjets_Pt10Eta2p4       ;

      Float_t   Event_CaloSumJetHadE             ;
      Float_t   Event_CaloSumJetPt               ;
      Float_t   Event_CaloSumJetPt_Pt10          ;
      Float_t   Event_CaloSumJetPt_Pt20          ;
      Float_t   Event_CaloSumJetPt_Pt10Eta2p4    ;
      Float_t   Event_CaloSumJetPt_Barrel        ;
      Float_t   Event_CaloSumJetPt_Endcap        ;

      Int_t     Event_CaloNjets                  ;
      Int_t     Event_CaloNjets_Pt10             ;
      Int_t     Event_CaloNjets_Pt20             ;
      Int_t     Event_CaloNjets_Pt10Eta2p4       ;

      Float_t   Event_GenPart1Pt                 ;
      Float_t   Event_GenPart2Pt                 ;


      // Tree filled for each calojet
      TTree *CaloJetTree;

      Float_t CaloJet_MatchedGenJet_Mass             ;
      Float_t CaloJet_MatchedGenJet_Energy           ;
      Float_t CaloJet_MatchedGenJet_Pt               ;
      Float_t CaloJet_MatchedGenJet_Eta              ;
      Float_t CaloJet_MatchedGenJet_Phi              ;
      Float_t CaloJet_MatchedGenJet_DeltaR           ;
      Float_t CaloJet_MatchedGenJet_HadPt            ;                               
      Float_t CaloJet_MatchedGenJet_HadEnergy        ;   

      Float_t CaloJet_Mass                           ;
      Float_t CaloJet_Energy                         ;
      Float_t CaloJet_Pt                             ;
      Float_t CaloJet_Eta                            ;
      Float_t CaloJet_Phi                            ;

      Float_t CaloJet_CorrPt                         ;
      Float_t CaloJet_CorrPtRhoArea                  ;
      Float_t CaloJet_CorrPtRhocentralArea           ;
      Float_t CaloJet_Area                           ;
      Float_t CaloJet_RhoArea                        ;
      Float_t CaloJet_RhoCentralArea                 ;
      Float_t CaloJet_RhoAreaOfficial                ;
                                                
      Float_t CaloJet_HadEnergy                      ;                                
      Float_t CaloJet_HadPt                          ;                                
      Float_t CaloJet_HadPt_Over_GenJet_HadPt      ;                                             
      Float_t CaloJet_HadE_Over_GenJet_HadE        ;                                           

      Float_t CaloJet_Pt_Over_GenJet_Pt            ;
      Float_t CaloJet_Pt_Minus_GenJet_Pt           ;
      Float_t CaloJet_CorrPt_Over_GenJet_Pt        ;
      Float_t CaloJet_CorrPt_Minus_GenJet_Pt       ;

      Float_t CaloJet_EMEB                           ;
      Float_t CaloJet_EfracHad                       ;
      Float_t CaloJet_hadHB                          ;
      Float_t CaloJet_hadHO                          ;
      Float_t CaloJet_hadHE                          ;
      Float_t CaloJet_hadHF                          ;
      Float_t CaloJet_Nconst                         ;


      // Tree filled for each PF jet
      TTree *PFJetTree;

      Int_t   PFJet_Run                          ;
      Int_t   PFJet_Lumi                         ;
      Int_t   PFJet_Event                        ;  
      Int_t   PFJet_Nvtx                         ;
      Int_t   PFJet_NvtxGood                     ;
      Float_t PFJet_Rho                          ;

      Float_t PFJet_MatchedGenJet_Mass                      ;
      Float_t PFJet_MatchedGenJet_Energy                    ;
      Float_t PFJet_MatchedGenJet_Pt                        ;
      Float_t PFJet_MatchedGenJet_Eta                       ;
      Float_t PFJet_MatchedGenJet_Phi                       ;
      Float_t PFJet_MatchedGenJet_DeltaR                    ;
      Float_t PFJet_MatchedGenJet_HadPt                     ;                       
      Float_t PFJet_MatchedGenJet_HadEnergy                 ; 

      Float_t PFJet_Area                         ;
      Float_t PFJet_CorrFactor                   ;

      Float_t PFJet_Mass                         ;
      Float_t PFJet_Energy                       ;
      Float_t PFJet_Pt                           ;
      Float_t PFJet_Eta                          ;
      Float_t PFJet_Phi                          ;

      Float_t PFJet_CorrMass                     ;
      Float_t PFJet_CorrEnergy                   ;
      Float_t PFJet_CorrPt                       ;
      Float_t PFJet_CorrEta                      ;
      Float_t PFJet_CorrPhi                      ;

      Float_t PFJet_CorrPtRhoArea                ;
      Float_t PFJet_RhoArea                      ;
      Float_t PFJet_RhoAreaOfficial              ;

      Float_t PFJet_Pt_Over_GenJet_Pt          ;
      Float_t PFJet_Pt_Minus_GenJet_Pt         ;
      Float_t PFJet_CorrPt_Over_GenJet_Pt      ;
      Float_t PFJet_CorrPt_Minus_GenJet_Pt     ;
                        
      Float_t PFJet_HadEnergy                    ;                      
      Float_t PFJet_HadPt                        ;                      
      Float_t PFJet_HadPt_Over_GenJet_HadPt    ;                                     
      Float_t PFJet_HadE_Over_GenJet_HadE      ;                                     
                                        
      Float_t PFJet_Nconst                       ; 
      Float_t PFJet_chargedEmEnergy              ;
      Float_t PFJet_chargedEmEnergyFraction      ;
      Float_t PFJet_chargedHadronEnergy          ;
      Float_t PFJet_chargedHadronEnergyFraction  ;
      Float_t PFJet_chargedHadronMultiplicity    ;
      Float_t PFJet_chargedMuEnergy              ;
      Float_t PFJet_chargedMuEnergyFraction      ;
      Float_t PFJet_chargedMultiplicity          ;
      Float_t PFJet_electronEnergy               ;
      Float_t PFJet_electronEnergyFraction       ;
      Float_t PFJet_electronMultiplicity         ;
      Float_t PFJet_HFEMEnergy                   ;
      Float_t PFJet_HFEMEnergyFraction           ;
      Float_t PFJet_HFEMMultiplicity             ;
      Float_t PFJet_HFHadronEnergy               ;
      Float_t PFJet_HFHadronEnergyFraction       ;
      Float_t PFJet_HFHadronMultiplicity         ;
      Float_t PFJet_muonEnergy                   ;
      Float_t PFJet_muonEnergyFraction           ;
      Float_t PFJet_muonMultiplicity             ;
      Float_t PFJet_neutralEmEnergy              ;
      Float_t PFJet_neutralEmEnergyFraction      ;
      Float_t PFJet_neutralHadronEnergy          ;
      Float_t PFJet_neutralHadronEnergyFraction  ;
      Float_t PFJet_neutralHadronMultiplicity    ;
      Float_t PFJet_neutralMultiplicity          ;
      Float_t PFJet_photonEnergy                 ;
      Float_t PFJet_photonEnergyFraction         ;
      Float_t PFJet_photonMultiplicity           ;

};

SimpleJetResponseTreeMaker::SimpleJetResponseTreeMaker(const edm::ParameterSet& iConfig) :
pvSrc_              (iConfig.getParameter<edm::InputTag>("pvSrc") ),
rhoSrc_             (iConfig.getParameter<edm::InputTag>("rhoSrc") ),
caloRhoSrc_         (iConfig.getParameter<edm::InputTag>("caloRhoSrc") ),
caloJetSrc_         (iConfig.getParameter<edm::InputTag>("caloJetSrc") ),
recoJetSrc_         (iConfig.getParameter<edm::InputTag>("recoJetSrc") ),
genJetSrc_          (iConfig.getParameter<edm::InputTag>("genJetSrc") ),
jecPayloads_        (iConfig.getParameter<std::vector<std::string> >  ("jecPayloads") ),
genJetMatchDeltaR_  (iConfig.getParameter<double>( "genJetMatchDeltaR" ) ),      
genJetPtThreshold_  (iConfig.getParameter<double>( "genJetPtThreshold" ) ),                 
recoJetPtThreshold_ (iConfig.getParameter<double>( "recoJetPtThreshold" ) )                   
{

  std::vector<JetCorrectorParameters> vPar;
  for ( std::vector<std::string>::const_iterator ipayload = jecPayloads_.begin(),
    ipayloadEnd = jecPayloads_.end(); ipayload != ipayloadEnd - 1; ++ipayload ) {
    std::cout << "Adding payload " << *ipayload << std::endl;
    JetCorrectorParameters pars(*ipayload);
    vPar.push_back(pars);
  }

  jec_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  jecUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloads_.back()));


}


SimpleJetResponseTreeMaker::~SimpleJetResponseTreeMaker()
{
}

// ------------ method called for each event  ------------
void
SimpleJetResponseTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;

  int run   = iEvent.id().run();
  int event = iEvent.id().event();
  int lumi  = iEvent.id().luminosityBlock();

  // Count reconstructed verties
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(pvSrc_, vertices);
  int count_vertex      = 0;
  int count_good_vertex = 0;
  for(reco::VertexCollection::const_iterator v=vertices->begin();v!=vertices->end(); ++v)
  {
    count_vertex++;
    if ( v->ndof() > 4 && fabs(v->z()) < 24 && fabs(v->position().rho()) < 2 ) count_good_vertex++;
  }

  // Get actual number of MC pileup (in time and out of time)
  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
  std::vector<PileupSummaryInfo>::const_iterator PVI;

  float npIT_BX0  =-1.;
  float npIT_BX1  =-1.;
  float npIT_BX2  =-1.;
  float npIT_BXm1 =-1.;
  float npIT_BXm2 =-1.;
  float npIT_BXm3 =-1.;

  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
    int BX = PVI->getBunchCrossing();
    //double npT = PVI->getTrueNumInteractions();
    double npIT = PVI->getPU_NumInteractions();
    if(BX == -3) npIT_BXm3 = npIT;
    if(BX == -2) npIT_BXm2 = npIT;
    if(BX == -1) npIT_BXm1 = npIT;
    if(BX ==  0) npIT_BX0  = npIT;
    if(BX ==  1) npIT_BX1  = npIT;
    if(BX ==  2) npIT_BX2  = npIT;
  }


  // get PF rho
  edm::Handle<double> rhoH;
  iEvent.getByLabel(rhoSrc_, rhoH );
  double rhoVal = *rhoH;

  // get calo rho
  edm::Handle<double> rhoCaloH;
  iEvent.getByLabel(caloRhoSrc_, rhoCaloH );
  double rhoValCalo = *rhoCaloH;

  // get central calo rho
  edm::Handle<double> rhoCaloCentralH;
  iEvent.getByLabel("fixedGridRhoFastjetCentralCalo", rhoCaloCentralH );
  double rhoValCaloCentral = *rhoCaloCentralH;


  // Get gen particle info (for particle gun samples)
  reco::Candidate::LorentzVector gen_pi1;
  reco::Candidate::LorentzVector gen_pi2;

  Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  for(size_t i = 0; i < genParticles->size(); ++ i) 
  {
    const reco::GenParticle & p = (*genParticles)[i];
    if ( p.pdgId() ==  211 ) gen_pi1 = p.p4();
    if ( p.pdgId() == -211 ) gen_pi2 = p.p4();
  }

  // GenJets
  edm::Handle<reco::GenJetCollection> genJetH;
  iEvent.getByLabel(genJetSrc_, genJetH); 
 
  //-------------------------------------------------
  // Calo Jets
  //-------------------------------------------------
  edm::Handle<reco::CaloJetCollection> calojetH;
  iEvent.getByLabel(caloJetSrc_, calojetH); 

  double calo_sumJetHadE             = 0;

  double calo_sumJetPt             = 0;
  double calo_sumJetPt_Pt10        = 0;
  double calo_sumJetPt_Pt20        = 0;
  double calo_sumJetPt_Pt10_Eta2p4 = 0;

  double calo_sumJetPt_Barrel = 0;
  double calo_sumJetPt_Endcap = 0;

  int calo_count_jets             = 0;
  int calo_count_jets_Pt10        = 0;
  int calo_count_jets_Pt20        = 0;
  int calo_count_jets_Pt10_Eta2p4 = 0;
  
  for ( reco::CaloJetCollection::const_iterator jet = calojetH->begin(); jet != calojetH->end(); ++jet ) {
    // JEC by hand
    reco::Candidate::LorentzVector uncorrJet = jet->p4();
    jec_->setJetEta( uncorrJet.eta() );
    jec_->setJetPt ( uncorrJet.pt() );
    jec_->setJetE  ( uncorrJet.energy() );
    jec_->setJetA  ( jet->jetArea() );
    jec_->setRho   ( rhoValCalo );
    vector<float> vcor;
    vcor = jec_->getSubCorrections();
    jec_->setJetEta( uncorrJet.eta() );
    jec_->setJetPt ( uncorrJet.pt() );
    jec_->setJetE  ( uncorrJet.energy() );
    jec_->setJetA  ( jet->jetArea() );
    jec_->setRho   ( rhoValCalo );
    reco::Candidate::PolarLorentzVector corrJet (uncorrJet.pt(), uncorrJet.eta(), uncorrJet.phi(), uncorrJet.mass());
    corrJet *=  ( vcor[0] );  //use only L1 correction

    // Eventwise calojet variables
    calo_sumJetHadE += jet->hadEnergyInHB() +  jet->hadEnergyInHO() +jet->hadEnergyInHE() +jet->hadEnergyInHF();
    calo_sumJetPt += uncorrJet.pt() ;       
    if (uncorrJet.pt()>10) calo_sumJetPt_Pt10 += uncorrJet.pt() ;             
    if (uncorrJet.pt()>20) calo_sumJetPt_Pt20 += uncorrJet.pt() ;             
    if (uncorrJet.pt()>10 && fabs(jet->eta()) < 2.4 ) calo_sumJetPt_Pt10_Eta2p4 += uncorrJet.pt() ;  

    calo_count_jets ++;            
    if (uncorrJet.pt()>10) calo_count_jets_Pt10  ++;     
    if (uncorrJet.pt()>20) calo_count_jets_Pt20  ++;     
    if (uncorrJet.pt()>10 && fabs(jet->eta()) < 2.4 ) calo_count_jets_Pt10_Eta2p4 ++;

    if ( fabs(jet->eta()) < 1.3)  calo_sumJetPt_Barrel += uncorrJet.pt() ;  
    if ( fabs(jet->eta()) > 1.3 && fabs(jet->eta()) < 3.0)  calo_sumJetPt_Endcap += uncorrJet.pt() ;  

    // GenJetMatch
    if (corrJet.pt() < recoJetPtThreshold_) continue;

    double closest_genjet_dR = 9999;
    reco::GenJet theMatchingGenJet;
    for (GenJetCollection::const_iterator genjet=genJetH->begin(); genjet!=genJetH->end(); genjet++) {
      if ( genjet->pt() < genJetPtThreshold_ || fabs(genjet->eta())>3 ) continue; 
      double deltar = deltaR( uncorrJet.eta(), uncorrJet.phi(), genjet->eta(), genjet->phi() );
      if ( deltar > genJetMatchDeltaR_ ) continue;
      if ( deltar < closest_genjet_dR ){
        closest_genjet_dR =  deltar;
        theMatchingGenJet = (*genjet);
      }
    }
    
    // CaloJetTree Fill
    if (closest_genjet_dR < genJetMatchDeltaR_)
    {
      double genHadPt = theMatchingGenJet.hadEnergy() * theMatchingGenJet.pt() / theMatchingGenJet.energy() ;

      CaloJet_MatchedGenJet_HadPt     = genHadPt;     
      CaloJet_MatchedGenJet_HadEnergy = theMatchingGenJet.hadEnergy(); 

      CaloJet_MatchedGenJet_Mass    = theMatchingGenJet.mass();     
      CaloJet_MatchedGenJet_Energy  = theMatchingGenJet.energy();   
      CaloJet_MatchedGenJet_Pt      = theMatchingGenJet.pt();        
      CaloJet_MatchedGenJet_Eta     = theMatchingGenJet.eta();       
      CaloJet_MatchedGenJet_Phi     = theMatchingGenJet.phi();    
      CaloJet_MatchedGenJet_DeltaR  = deltaR( uncorrJet.eta(), uncorrJet.phi(), theMatchingGenJet.eta(), theMatchingGenJet.phi() );

      CaloJet_Eta        = uncorrJet.eta();        
      CaloJet_Mass       = uncorrJet.mass();     
      CaloJet_Energy     = uncorrJet.energy();   
      CaloJet_Pt         = uncorrJet.pt();         
    
      double rhoCaloArea        = jet->jetArea() * rhoValCalo ;
      double rhoCaloCentralArea = jet->jetArea() * rhoValCaloCentral ;

      CaloJet_CorrPt                = corrJet.pt() ;         
      CaloJet_CorrPtRhoArea         = uncorrJet.pt()-rhoCaloArea ;         
      CaloJet_CorrPtRhocentralArea  = uncorrJet.pt()-rhoCaloCentralArea ;         
      CaloJet_Area                  = jet->jetArea();
      CaloJet_RhoArea               = jet->jetArea() * rhoValCalo ;
      CaloJet_RhoCentralArea        = jet->jetArea() * rhoCaloCentralArea ;
      CaloJet_RhoAreaOfficial       = uncorrJet.pt()-corrJet.pt() ;

      CaloJet_EMEB         = jet->emEnergyInEB();
      CaloJet_EfracHad     = jet->energyFractionHadronic();
      CaloJet_hadHB        = jet->hadEnergyInHB();
      CaloJet_hadHO        = jet->hadEnergyInHO();
      CaloJet_hadHE        = jet->hadEnergyInHE();
      CaloJet_hadHF        = jet->hadEnergyInHF();
      CaloJet_Nconst       = jet->getCaloConstituents().size();

      double hadE  = jet->hadEnergyInHB() +  jet->hadEnergyInHO() +jet->hadEnergyInHE() +jet->hadEnergyInHF();
      double hadPt = hadE * uncorrJet.pt() / uncorrJet.energy() ;

      CaloJet_HadEnergy     = hadE;   
      CaloJet_HadPt         = hadPt;   

      if (genHadPt> 0) CaloJet_HadPt_Over_GenJet_HadPt = hadPt/genHadPt;
      if (theMatchingGenJet.hadEnergy() > 0) CaloJet_HadE_Over_GenJet_HadE = hadE/theMatchingGenJet.hadEnergy();


      if (theMatchingGenJet.pt() > 0) CaloJet_Pt_Over_GenJet_Pt = uncorrJet.pt()/theMatchingGenJet.pt();
      CaloJet_Pt_Minus_GenJet_Pt = uncorrJet.pt() - theMatchingGenJet.pt();
      CaloJetTree->Fill();

    }
  }

  //-------------------------------------------------
  // PF Jets
  //-------------------------------------------------
  edm::Handle<reco::PFJetCollection> pfjetH;
  iEvent.getByLabel(recoJetSrc_, pfjetH); 

  double sumJetHadE             = 0;

  double sumJetPt               = 0;
  double sumJetPt_Pt10          = 0;
  double sumJetPt_Pt20          = 0;
  double sumJetPt_Pt10_Eta2p4   = 0;
  double sumJetPt_Barrel        = 0;
  double sumJetPt_Endcap        = 0;
  int count_jets                = 0;
  int count_jets_Pt10           = 0;
  int count_jets_Pt20           = 0;
  int count_jets_Pt10_Eta2p4    = 0;

  double sumJetCorrPt             = 0;
  double sumJetCorrPt_Pt10        = 0;
  double sumJetCorrPt_Pt20        = 0;
  double sumJetCorrPt_Pt10_Eta2p4 = 0;

  int count_jets_corr             = 0;
  int count_jets_corr_Pt10        = 0;
  int count_jets_corr_Pt20        = 0;
  int count_jets_corr_Pt10_Eta2p4 = 0;


  for ( reco::PFJetCollection::const_iterator jet = pfjetH->begin(); jet != pfjetH->end(); ++jet ) {

    // Get uncorrected and corrected jet
    reco::Candidate::LorentzVector uncorrJet = jet->p4();
    jec_->setJetEta( uncorrJet.eta() );
    jec_->setJetPt ( uncorrJet.pt() );
    jec_->setJetE  ( uncorrJet.energy() );
    jec_->setJetA  ( jet->jetArea() );
    jec_->setRho   ( rhoVal );
    // weird behavior...have to setJet multiple times
    double corr = jec_->getCorrection();
    jec_->setJetEta( uncorrJet.eta() );
    jec_->setJetPt ( uncorrJet.pt() );
    jec_->setJetE  ( uncorrJet.energy() );
    jec_->setJetA  ( jet->jetArea() );
    jec_->setRho   ( rhoVal );
    vector<float> vcor;
    vcor = jec_->getSubCorrections();
    jec_->setJetEta( uncorrJet.eta() );
    jec_->setJetPt ( uncorrJet.pt() );
    jec_->setJetE  ( uncorrJet.energy() );
    jec_->setJetA  ( jet->jetArea() );
    jec_->setRho   ( rhoVal );
    reco::Candidate::PolarLorentzVector corrJet (uncorrJet.pt(), uncorrJet.eta(), uncorrJet.phi(), uncorrJet.mass());
    corrJet *=  ( vcor[0] );  // Use only L1 correction

    // Calculate eventwise jet variables
    sumJetHadE += jet->chargedHadronEnergy() + jet->neutralHadronEnergy();
    sumJetPt += uncorrJet.pt() ;       
    if (uncorrJet.pt()>10) sumJetPt_Pt10 += uncorrJet.pt() ;             
    if (uncorrJet.pt()>20) sumJetPt_Pt20 += uncorrJet.pt() ;             
    if (uncorrJet.pt()>10 && fabs(jet->eta()) < 2.4 ) sumJetPt_Pt10_Eta2p4 += uncorrJet.pt() ;  

    sumJetCorrPt += corrJet.pt() ;       
    if (corrJet.pt()>10) sumJetCorrPt_Pt10 += corrJet.pt() ;             
    if (corrJet.pt()>20) sumJetCorrPt_Pt20 += corrJet.pt() ;             
    if (corrJet.pt()>10 && fabs(jet->eta()) < 2.4 ) sumJetCorrPt_Pt10_Eta2p4 += corrJet.pt() ;  

    if ( fabs(jet->eta()) < 1.3)  sumJetPt_Barrel += uncorrJet.pt() ;  
    if ( fabs(jet->eta()) > 1.3 && fabs(jet->eta()) < 3.0)  sumJetPt_Endcap += uncorrJet.pt() ;  

    count_jets ++;            
    if (uncorrJet.pt()>10) count_jets_Pt10  ++;     
    if (uncorrJet.pt()>20) count_jets_Pt20  ++;     
    if (uncorrJet.pt()>10 && fabs(jet->eta()) < 2.4 ) count_jets_Pt10_Eta2p4 ++;

    if (corrJet.pt()>10) count_jets_corr_Pt10  ++;     
    if (corrJet.pt()>20) count_jets_corr_Pt20  ++;     
    if (corrJet.pt()>10 && fabs(jet->eta()) < 2.4 ) count_jets_corr_Pt10_Eta2p4 ++;

    // Gen Matching

    if (corrJet.pt() < recoJetPtThreshold_) continue;  

    double closest_genjet_dR = 9999;
    reco::GenJet theMatchingGenJet;
    for (GenJetCollection::const_iterator genjet=genJetH->begin(); genjet!=genJetH->end(); genjet++) {
      if ( genjet->pt() < genJetPtThreshold_ || fabs(genjet->eta())>3 ) continue; 
      double deltar = deltaR( uncorrJet.eta(), uncorrJet.phi(), genjet->eta(), genjet->phi() );
      if ( deltar > genJetMatchDeltaR_ ) continue;
      if ( deltar < closest_genjet_dR ){
        closest_genjet_dR =  deltar;
        theMatchingGenJet = (*genjet);
      }
    }

    // PFJetTree Fill for each matched jet
    if (closest_genjet_dR < genJetMatchDeltaR_)
    {
      PFJet_Run        = run;
      PFJet_Lumi       = lumi;
      PFJet_Event      = event;
      PFJet_Nvtx       = count_vertex;     
      PFJet_NvtxGood   = count_good_vertex;
      PFJet_Rho        = rhoVal;

      double genHadPt = theMatchingGenJet.hadEnergy() * theMatchingGenJet.pt()  / theMatchingGenJet.energy() ;

      PFJet_MatchedGenJet_HadPt     = genHadPt;     
      PFJet_MatchedGenJet_HadEnergy = theMatchingGenJet.hadEnergy();     
      PFJet_MatchedGenJet_Mass      = theMatchingGenJet.mass();     
      PFJet_MatchedGenJet_Energy    = theMatchingGenJet.energy();   
      PFJet_MatchedGenJet_Pt        = theMatchingGenJet.pt();        
      PFJet_MatchedGenJet_Eta       = theMatchingGenJet.eta();       
      PFJet_MatchedGenJet_Phi       = theMatchingGenJet.phi();    
      PFJet_MatchedGenJet_DeltaR    = deltaR( uncorrJet.eta(), uncorrJet.phi(), theMatchingGenJet.eta(), theMatchingGenJet.phi() );

      PFJet_Area       = jet->jetArea() ;      
      PFJet_CorrFactor = corr;

      double hadE = jet->chargedHadronEnergy() + jet->neutralHadronEnergy();
      double hadPt = hadE * uncorrJet.pt() / uncorrJet.E();
       
      PFJet_HadEnergy = hadE;
      PFJet_HadPt     = hadPt;

      PFJet_Mass       = uncorrJet.mass();     
      PFJet_Energy     = uncorrJet.energy();   
      PFJet_Pt         = uncorrJet.pt();         
      PFJet_Eta        = uncorrJet.eta();        
      PFJet_Phi        = uncorrJet.phi();
      PFJet_CorrMass   = corrJet.mass();     
      PFJet_CorrEnergy = corrJet.energy();   
      PFJet_CorrPt     = corrJet.pt();         
      PFJet_CorrEta    = corrJet.eta();        
      PFJet_CorrPhi    = corrJet.phi();  
     
      double rhoArea        = jet->jetArea() * rhoVal ;
      PFJet_CorrPtRhoArea         = uncorrJet.pt()-rhoArea ;         
      PFJet_RhoArea               = jet->jetArea() * rhoValCalo ;
      PFJet_RhoAreaOfficial       = uncorrJet.pt()-corrJet.pt() ;

      if (theMatchingGenJet.hadEnergy() > 0) PFJet_HadPt_Over_GenJet_HadPt = hadPt/theMatchingGenJet.hadEnergy() ;
      if (genHadPt > 0)                      PFJet_HadE_Over_GenJet_HadE   = hadE/genHadPt;

      if (theMatchingGenJet.pt() > 0) PFJet_Pt_Over_GenJet_Pt = uncorrJet.pt()/theMatchingGenJet.pt();
      PFJet_Pt_Minus_GenJet_Pt = uncorrJet.pt() - theMatchingGenJet.pt();
      if (theMatchingGenJet.pt() > 0) PFJet_CorrPt_Over_GenJet_Pt = corrJet.pt()/theMatchingGenJet.pt();
      PFJet_CorrPt_Minus_GenJet_Pt = corrJet.pt() - theMatchingGenJet.pt();

      PFJet_Nconst                      = jet->numberOfDaughters            ();  
      PFJet_chargedEmEnergy             = jet->chargedEmEnergy              ();  
      PFJet_chargedEmEnergyFraction     = jet->chargedEmEnergyFraction      ();  
      PFJet_chargedHadronEnergy         = jet->chargedHadronEnergy          ();  
      PFJet_chargedHadronEnergyFraction = jet->chargedHadronEnergyFraction  ();  
      PFJet_chargedHadronMultiplicity   = jet->chargedHadronMultiplicity    ();  
      PFJet_chargedMuEnergy             = jet->chargedMuEnergy              ();  
      PFJet_chargedMuEnergyFraction     = jet->chargedMuEnergyFraction      ();  
      PFJet_chargedMultiplicity         = jet->chargedMultiplicity          ();  
      PFJet_electronEnergy              = jet->electronEnergy               ();  
      PFJet_electronEnergyFraction      = jet->electronEnergyFraction       ();  
      PFJet_electronMultiplicity        = jet->electronMultiplicity         ();  
      PFJet_HFEMEnergy                  = jet->HFEMEnergy                   ();  
      PFJet_HFEMEnergyFraction          = jet->HFEMEnergyFraction           ();  
      PFJet_HFEMMultiplicity            = jet->HFEMMultiplicity             ();  
      PFJet_HFHadronEnergy              = jet->HFHadronEnergy               ();  
      PFJet_HFHadronEnergyFraction      = jet->HFHadronEnergyFraction       ();  
      PFJet_HFHadronMultiplicity        = jet->HFHadronMultiplicity         ();  
      PFJet_muonEnergy                  = jet->muonEnergy                   ();  
      PFJet_muonEnergyFraction          = jet->muonEnergyFraction           ();  
      PFJet_muonMultiplicity            = jet->muonMultiplicity             ();  
      PFJet_neutralEmEnergy             = jet->neutralEmEnergy              ();  
      PFJet_neutralEmEnergyFraction     = jet->neutralEmEnergyFraction      ();  
      PFJet_neutralHadronEnergy         = jet->neutralHadronEnergy          ();  
      PFJet_neutralHadronEnergyFraction = jet->neutralHadronEnergyFraction  ();  
      PFJet_neutralHadronMultiplicity   = jet->neutralHadronMultiplicity    ();  
      PFJet_neutralMultiplicity         = jet->neutralMultiplicity          ();  
      PFJet_photonEnergy                = jet->photonEnergy                 ();  
      PFJet_photonEnergyFraction        = jet->photonEnergyFraction         ();  
      PFJet_photonMultiplicity          = jet->photonMultiplicity           ();  

      PFJetTree->Fill();
    }
  }

  //-------------------------------------------------
  // Event tree
  //-------------------------------------------------
  Event_Run                        = run;
  Event_Lumi                       = lumi;
  Event_Event                      = event;
  Event_Nvtx                       = count_vertex;     
  Event_NvtxGood                   = count_good_vertex;
  Event_Rho                        = rhoVal;
  Event_RhoCalo                    = rhoValCalo;
  Event_RhoCaloCentral             = rhoValCaloCentral;

  Event_NPU_BXm3                   = npIT_BXm3 ;
  Event_NPU_BXm2                   = npIT_BXm2 ;
  Event_NPU_BXm1                   = npIT_BXm1 ;
  Event_NPU_BX0                    = npIT_BX0  ;
  Event_NPU_BX1                    = npIT_BX1  ;
  Event_NPU_BX2                    = npIT_BX2  ;

  Event_SumJetPt                   = sumJetPt;
  Event_SumJetPt_Pt10              = sumJetPt_Pt10;
  Event_SumJetPt_Pt20              = sumJetPt_Pt20;
  Event_SumJetPt_Pt10Eta2p4        = sumJetPt_Pt10_Eta2p4;
  Event_SumJetPt_Barrel            = sumJetPt_Barrel;
  Event_SumJetPt_Endcap            = sumJetPt_Endcap;
  Event_Njets                      = count_jets;
  Event_Njets_Pt10                 = count_jets_Pt10;
  Event_Njets_Pt10Eta2p4           = count_jets_Pt10_Eta2p4;

  Event_CorrSumJetHadE             = sumJetHadE;
  Event_CorrSumJetPt               = sumJetCorrPt;
  Event_CorrSumJetPt_Pt10          = sumJetCorrPt_Pt10;
  Event_CorrSumJetPt_Pt20          = sumJetCorrPt_Pt20;
  Event_CorrSumJetPt_Pt10Eta2p4    = sumJetCorrPt_Pt10_Eta2p4;
  Event_CorrNjets                  = count_jets_corr;
  Event_CorrNjets_Pt10             = count_jets_corr_Pt10;
  Event_CorrNjets_Pt10Eta2p4       = count_jets_corr_Pt10_Eta2p4;

  Event_CaloSumJetHadE             = calo_sumJetHadE;
  Event_CaloSumJetPt               = calo_sumJetPt;
  Event_CaloSumJetPt_Pt10          = calo_sumJetPt_Pt10;
  Event_CaloSumJetPt_Pt20          = calo_sumJetPt_Pt20;
  Event_CaloSumJetPt_Pt10Eta2p4    = calo_sumJetPt_Pt10_Eta2p4;
  Event_CaloNjets                  = calo_count_jets;
  Event_CaloNjets_Pt10             = calo_count_jets_Pt10;
  Event_CaloNjets_Pt10Eta2p4       = calo_count_jets_Pt10_Eta2p4;
  Event_CaloSumJetPt_Barrel        = calo_sumJetPt_Barrel;
  Event_CaloSumJetPt_Endcap        = calo_sumJetPt_Endcap;

  Event_GenPart1Pt                 = gen_pi1.pt();
  Event_GenPart2Pt                 = gen_pi2.pt();

  EventTree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
SimpleJetResponseTreeMaker::beginJob()
{

  edm::Service<TFileService> fs;
  TFileDirectory subDir3 = fs->mkdir("Trees");

  EventTree = new TTree("EventTree","EventTree");
  EventTree->Branch("Event_Run"                         ,  & Event_Run                       , "Event_Run/I"                        );
  EventTree->Branch("Event_Lumi"                        ,  & Event_Lumi                      , "Event_Lumi/I"                       );
  EventTree->Branch("Event_Event"                       ,  & Event_Event                     , "Event_Event/I"                      );
  EventTree->Branch("Event_Nvtx"                        ,  & Event_Nvtx                      , "Event_Nvtx/I"                       );
  EventTree->Branch("Event_NvtxGood"                    ,  & Event_NvtxGood                  , "Event_NvtxGood/I"                   );
  EventTree->Branch("Event_Rho"                         ,  & Event_Rho                       , "Event_Rho/F"                        );
  EventTree->Branch("Event_RhoCalo"                     ,  & Event_RhoCalo                   , "Event_RhoCalo/F"                    );
  EventTree->Branch("Event_RhoCaloCentral"              ,  & Event_RhoCaloCentral            , "Event_RhoCaloCentral/F"             );
      
  EventTree->Branch("Event_NPU_BXm3"                    ,  & Event_NPU_BXm3                  , "Event_NPU_BXm3/I"                        );
  EventTree->Branch("Event_NPU_BXm2"                    ,  & Event_NPU_BXm2                  , "Event_NPU_BXm2/I"                        );
  EventTree->Branch("Event_NPU_BXm1"                    ,  & Event_NPU_BXm1                  , "Event_NPU_BXm1/I"                        );
  EventTree->Branch("Event_NPU_BX0"                     ,  & Event_NPU_BX0                   , "Event_NPU_BX0/I"                         );
  EventTree->Branch("Event_NPU_BX1"                     ,  & Event_NPU_BX1                   , "Event_NPU_BX1/I"                         );
  EventTree->Branch("Event_NPU_BX2"                     ,  & Event_NPU_BX2                   , "Event_NPU_BX2/I"                         );
    
  EventTree->Branch("Event_CorrSumJetHadE"              ,  & Event_CorrSumJetHadE            , "Event_CorrSumJetHadE/F"             );
  EventTree->Branch("Event_SumJetPt"                    ,  & Event_SumJetPt                  , "Event_SumJetPt/F"                   );
  EventTree->Branch("Event_SumJetPt_Pt10"               ,  & Event_SumJetPt_Pt10             , "Event_SumJetPt_Pt10/F"              );
  EventTree->Branch("Event_SumJetPt_Pt20"               ,  & Event_SumJetPt_Pt20             , "Event_SumJetPt_Pt20/F"              );
  EventTree->Branch("Event_SumJetPt_Pt10Eta2p4"         ,  & Event_SumJetPt_Pt10Eta2p4       , "Event_SumJetPt_Pt10Eta2p4/F"        );
  EventTree->Branch("Event_SumJetPt_Barrel"             ,  & Event_SumJetPt_Barrel           , "Event_SumJetPt_Barrel/F"            );
  EventTree->Branch("Event_SumJetPt_Endcap"             ,  & Event_SumJetPt_Endcap           , "Event_SumJetPt_Endcap/F"            );
  EventTree->Branch("Event_Njets"                       ,  & Event_Njets                     , "Event_Njets/I"                      );
  EventTree->Branch("Event_Njets_Pt10"                  ,  & Event_Njets_Pt10                , "Event_Njets_Pt10/I"                 );
  EventTree->Branch("Event_Njets_Pt20"                  ,  & Event_Njets_Pt20                , "Event_Njets_Pt20/I"                 );
  EventTree->Branch("Event_Njets_Pt10Eta2p4"            ,  & Event_Njets_Pt10Eta2p4          , "Event_Njets_Pt10Eta2p4/I"           );

  EventTree->Branch("Event_CorrSumJetPt"                ,  & Event_CorrSumJetPt              , "Event_CorrSumJetPt/F"                   );
  EventTree->Branch("Event_CorrSumJetPt_Pt10"           ,  & Event_CorrSumJetPt_Pt10         , "Event_CorrSumJetPt_Pt10/F"              );
  EventTree->Branch("Event_CorrSumJetPt_Pt20"           ,  & Event_CorrSumJetPt_Pt20         , "Event_CorrSumJetPt_Pt20/F"              );
  EventTree->Branch("Event_CorrSumJetPt_Pt10Eta2p4"     ,  & Event_CorrSumJetPt_Pt10Eta2p4   , "Event_CorrSumJetPt_Pt10Eta2p4/F"        );
  EventTree->Branch("Event_CorrNjets"                   ,  & Event_CorrNjets                 , "Event_CorrNjets/I"                      );
  EventTree->Branch("Event_CorrNjets_Pt10"              ,  & Event_CorrNjets_Pt10            , "Event_CorrNjets_Pt10/I"                 );
  EventTree->Branch("Event_CorrNjets_Pt20"              ,  & Event_CorrNjets_Pt20            , "Event_CorrNjets_Pt20/I"                 );
  EventTree->Branch("Event_CorrNjets_Pt10Eta2p4"        ,  & Event_CorrNjets_Pt10Eta2p4      , "Event_CorrNjets_Pt10Eta2p4/I"           );

  EventTree->Branch("Event_CaloSumJetHadE"              ,  & Event_CaloSumJetHadE            , "Event_CaloSumJetHadE/F"                 );
  EventTree->Branch("Event_CaloSumJetPt"                ,  & Event_CaloSumJetPt              , "Event_CaloSumJetPt/F"                   );
  EventTree->Branch("Event_CaloSumJetPt_Pt10"           ,  & Event_CaloSumJetPt_Pt10         , "Event_CaloSumJetPt_Pt10/F"              );
  EventTree->Branch("Event_CaloSumJetPt_Pt20"           ,  & Event_CaloSumJetPt_Pt20         , "Event_CaloSumJetPt_Pt20/F"              );
  EventTree->Branch("Event_CaloSumJetPt_Pt10Eta2p4"     ,  & Event_CaloSumJetPt_Pt10Eta2p4   , "Event_CaloSumJetPt_Pt10Eta2p4/F"        );
  EventTree->Branch("Event_CaloSumJetPt_Barrel"         ,  & Event_CaloSumJetPt_Barrel       , "Event_CaloSumJetPt_Barrel/F"            );
  EventTree->Branch("Event_CaloSumJetPt_Endcap"         ,  & Event_CaloSumJetPt_Endcap       , "Event_CaloSumJetPt_Endcap/F"            );
  EventTree->Branch("Event_CaloNjets"                   ,  & Event_CaloNjets                 , "Event_CaloNjets/I"                      );
  EventTree->Branch("Event_CaloNjets_Pt10"              ,  & Event_CaloNjets_Pt10            , "Event_CaloNjets_Pt10/I"                 );
  EventTree->Branch("Event_CaloNjets_Pt20"              ,  & Event_CaloNjets_Pt20            , "Event_CaloNjets_Pt20/I"                 );
  EventTree->Branch("Event_CaloNjets_Pt10Eta2p4"        ,  & Event_CaloNjets_Pt10Eta2p4      , "Event_CaloNjets_Pt10Eta2p4/I"           );

  EventTree->Branch("Event_GenPart1Pt"                  ,  & Event_GenPart1Pt                , "Event_GenPart1Pt/F"                     );
  EventTree->Branch("Event_GenPart2Pt"                  ,  & Event_GenPart2Pt                , "Event_GenPart2Pt/F"                     );


  CaloJetTree = new TTree("CaloJetTree","CaloJetTree");

  CaloJetTree->Branch("CaloJet_MatchedGenJet_Mass"            , & CaloJet_MatchedGenJet_Mass            , "CaloJet_MatchedGenJet_Mass/F"                       );
  CaloJetTree->Branch("CaloJet_MatchedGenJet_Energy"          , & CaloJet_MatchedGenJet_Energy          , "CaloJet_MatchedGenJet_Energy/F"                     ); 
  CaloJetTree->Branch("CaloJet_MatchedGenJet_Pt"              , & CaloJet_MatchedGenJet_Pt              , "CaloJet_MatchedGenJet_Pt/F"                         ); 
  CaloJetTree->Branch("CaloJet_MatchedGenJet_Eta"             , & CaloJet_MatchedGenJet_Eta             , "CaloJet_MatchedGenJet_Eta/F"                        ); 
  CaloJetTree->Branch("CaloJet_MatchedGenJet_Phi"             , & CaloJet_MatchedGenJet_Phi             , "CaloJet_MatchedGenJet_Phi/F"                        ); 
  CaloJetTree->Branch("CaloJet_MatchedGenJet_DeltaR"          , & CaloJet_MatchedGenJet_DeltaR          , "CaloJet_MatchedGenJet_DeltaR/F"                     ); 
  CaloJetTree->Branch("CaloJet_MatchedGenJet_HadPt"           , & CaloJet_MatchedGenJet_HadPt           , "CaloJet_MatchedGenJet_HadPt/F"                      ); 
  CaloJetTree->Branch("CaloJet_MatchedGenJet_HadEnergy"       , & CaloJet_MatchedGenJet_HadEnergy       , "CaloJet_MatchedGenJet_HadEnergy/F"                  ); 

  CaloJetTree->Branch("CaloJet_Mass"                          , & CaloJet_Mass                          , "CaloJet_Mass/F"                          );
  CaloJetTree->Branch("CaloJet_Energy"                        , & CaloJet_Energy                        , "CaloJet_Energy/F"                        ); 
  CaloJetTree->Branch("CaloJet_Pt"                            , & CaloJet_Pt                            , "CaloJet_Pt/F"                            ); 
  CaloJetTree->Branch("CaloJet_Eta"                           , & CaloJet_Eta                           , "CaloJet_Eta/F"                           ); 
  CaloJetTree->Branch("CaloJet_Phi"                           , & CaloJet_Phi                           , "CaloJet_Phi/F"                           ); 

  CaloJetTree->Branch("CaloJet_CorrPt"                        , & CaloJet_CorrPt                        , "CaloJet_CorrPt/F"                        ); 
  CaloJetTree->Branch("CaloJet_CorrPtRhoArea"                 , & CaloJet_CorrPtRhoArea                 , "CaloJet_CorrPtRhoArea/F"                 ); 
  CaloJetTree->Branch("CaloJet_CorrPtRhocentralArea"          , & CaloJet_CorrPtRhocentralArea          , "CaloJet_CorrPtRhocentralArea/F"          ); 
  CaloJetTree->Branch("CaloJet_Area"                          , & CaloJet_Area                          , "CaloJet_Area/F"                          ); 
  CaloJetTree->Branch("CaloJet_RhoArea"                       , & CaloJet_RhoArea                       , "CaloJet_RhoArea/F"                       ); 
  CaloJetTree->Branch("CaloJet_RhoCentralArea"                , & CaloJet_RhoCentralArea                , "CaloJet_RhoCentralArea/F"                ); 
  CaloJetTree->Branch("CaloJet_RhoAreaOfficial"               , & CaloJet_RhoAreaOfficial               , "CaloJet_RhoAreaOfficial/F"               ); 
   
  CaloJetTree->Branch("CaloJet_Pt_Over_GenJet_Pt"             , & CaloJet_Pt_Over_GenJet_Pt             , "CaloJet_Pt_Over_GenJet_Pt/F"           ); 
  CaloJetTree->Branch("CaloJet_Pt_Minus_GenJet_Pt"            , & CaloJet_Pt_Minus_GenJet_Pt            , "CaloJet_Pt_Minus_GenJet_Pt/F"          ); 

  CaloJetTree->Branch("CaloJet_EMEB"                          , & CaloJet_EMEB                          , "CaloJet_EMEB/F"                          );
  CaloJetTree->Branch("CaloJet_EfracHad"                      , & CaloJet_EfracHad                      , "CaloJet_EfracHad/F"                      ); 
  CaloJetTree->Branch("CaloJet_hadHB"                         , & CaloJet_hadHB                         , "CaloJet_hadHB/F"                         ); 
  CaloJetTree->Branch("CaloJet_hadHO"                         , & CaloJet_hadHO                         , "CaloJet_hadHO/F"                         ); 
  CaloJetTree->Branch("CaloJet_hadHE"                         , & CaloJet_hadHE                         , "CaloJet_hadHE/F"                         ); 
  CaloJetTree->Branch("CaloJet_hadHF"                         , & CaloJet_hadHF                         , "CaloJet_hadHF/F"                         ); 
  CaloJetTree->Branch("CaloJet_Nconst"                        , & CaloJet_Nconst                        , "CaloJet_Nconst/F"                        ); 

  CaloJetTree->Branch("CaloJet_HadEnergy"                     , & CaloJet_HadEnergy                     , "CaloJet_HadEnergy/F"                     ); 
  CaloJetTree->Branch("CaloJet_HadPt"                         , & CaloJet_HadPt                         , "CaloJet_HadPt/F"                         ); 
  CaloJetTree->Branch("CaloJet_HadPt_Over_GenJet_HadPt"       , & CaloJet_HadPt_Over_GenJet_HadPt       , "CaloJet_HadPt_Over_GenJet_HadPt/F"     ); 
  CaloJetTree->Branch("CaloJet_HadE_Over_GenJet_HadE"         , & CaloJet_HadE_Over_GenJet_HadE         , "CaloJet_HadE_Over_GenJet_HadE/F"       ); 



  PFJetTree = new TTree("PFJetTree","PFJetTree");

  PFJetTree->Branch("PFJet_Run"                         , & PFJet_Run                         , "PFJet_Run/I"                         );
  PFJetTree->Branch("PFJet_Lumi"                        , & PFJet_Lumi                        , "PFJet_Lumi/I"                        );
  PFJetTree->Branch("PFJet_Event"                       , & PFJet_Event                       , "PFJet_Event_/I"                      );
  PFJetTree->Branch("PFJet_Nvtx"                        , & PFJet_Nvtx                        , "PFJet_Nvtx/I"                        );
  PFJetTree->Branch("PFJet_NvtxGood"                    , & PFJet_NvtxGood                    , "PFJet_NvtxGood/I"                    );
  PFJetTree->Branch("PFJet_Rho"                         , & PFJet_Rho                         , "PFJet_Rho/F"                         );
 
  PFJetTree->Branch("PFJet_MatchedGenJet_Mass"          , & PFJet_MatchedGenJet_Mass          , "PFJet_MatchedGenJet_Mass/F"                     );
  PFJetTree->Branch("PFJet_MatchedGenJet_Energy"        , & PFJet_MatchedGenJet_Energy        , "PFJet_MatchedGenJet_Energy/F"                   ); 
  PFJetTree->Branch("PFJet_MatchedGenJet_Pt"            , & PFJet_MatchedGenJet_Pt            , "PFJet_MatchedGenJet_Pt/F"                       ); 
  PFJetTree->Branch("PFJet_MatchedGenJet_Eta"           , & PFJet_MatchedGenJet_Eta           , "PFJet_MatchedGenJet_Eta/F"                      ); 
  PFJetTree->Branch("PFJet_MatchedGenJet_Phi"           , & PFJet_MatchedGenJet_Phi           , "PFJet_MatchedGenJet_Phi/F"                      ); 
  PFJetTree->Branch("PFJet_MatchedGenJet_DeltaR"        , & PFJet_MatchedGenJet_DeltaR        , "PFJet_MatchedGenJet_DeltaR/F"                   ); 
  PFJetTree->Branch("PFJet_MatchedGenJet_HadPt"         , & PFJet_MatchedGenJet_HadPt         , "PFJet_MatchedGenJet_HadPt/F"                    );
  PFJetTree->Branch("PFJet_MatchedGenJet_HadEnergy"     , & PFJet_MatchedGenJet_HadEnergy     , "PFJet_MatchedGenJet_HadEnergy/F"                );

  PFJetTree->Branch("PFJet_Area"                        , & PFJet_Area                        , "PFJet_Area/F"                        ); 
  PFJetTree->Branch("PFJet_CorrFactor"                  , & PFJet_CorrFactor                  , "PFJet_CorrFactor/F"                  ); 

  PFJetTree->Branch("PFJet_Mass"                        , & PFJet_Mass                        , "PFJet_Mass/F"                        );
  PFJetTree->Branch("PFJet_Energy"                      , & PFJet_Energy                      , "PFJet_Energy/F"                      ); 
  PFJetTree->Branch("PFJet_Pt"                          , & PFJet_Pt                          , "PFJet_Pt/F"                          ); 
  PFJetTree->Branch("PFJet_Eta"                         , & PFJet_Eta                         , "PFJet_Eta/F"                         ); 
  PFJetTree->Branch("PFJet_Phi"                         , & PFJet_Phi                         , "PFJet_Phi/F"                         ); 

  PFJetTree->Branch("PFJet_CorrMass"                    , & PFJet_CorrMass                    , "PFJet_CorrMass/F"                    ); 
  PFJetTree->Branch("PFJet_CorrEnergy"                  , & PFJet_CorrEnergy                  , "PFJet_CorrEnergy/F"                  ); 
  PFJetTree->Branch("PFJet_CorrPt"                      , & PFJet_CorrPt                      , "PFJet_CorrPt/F"                      ); 
  PFJetTree->Branch("PFJet_CorrEta"                     , & PFJet_CorrEta                     , "PFJet_CorrEta/F"                     ); 
  PFJetTree->Branch("PFJet_CorrPhi"                     , & PFJet_CorrPhi                     , "PFJet_CorrPhi/F"                     ); 

  PFJetTree->Branch("PFJet_CorrPtRhoArea"               , & PFJet_CorrPtRhoArea               , "PFJet_CorrPtRhoArea/F"               ); 
  PFJetTree->Branch("PFJet_RhoArea"                     , & PFJet_RhoArea                     , "PFJet_RhoArea/F"                     ); 
  PFJetTree->Branch("PFJet_RhoAreaOfficial"             , & PFJet_RhoAreaOfficial             , "PFJet_RhoAreaOfficial/F"             ); 

  PFJetTree->Branch("PFJet_Pt_Over_GenJet_Pt"           , & PFJet_Pt_Over_GenJet_Pt           , "PFJet_Pt_Over_GenJet_Pt/F"           ); 
  PFJetTree->Branch("PFJet_Pt_Minus_GenJet_Pt"          , & PFJet_Pt_Minus_GenJet_Pt          , "PFJet_Pt_Minus_GenJet_Pt/F"          ); 
  PFJetTree->Branch("PFJet_CorrPt_Over_GenJet_Pt"       , & PFJet_CorrPt_Over_GenJet_Pt       , "PFJet_CorrPt_Over_GenJet_Pt/F"       ); 
  PFJetTree->Branch("PFJet_CorrPt_Minus_GenJet_Pt"      , & PFJet_CorrPt_Minus_GenJet_Pt      , "PFJet_CorrPt_Minus_GenJet_Pt/F"      ); 

  PFJetTree->Branch("PFJet_Nconst"                      , & PFJet_Nconst                      , "PFJet_Nconst/F"                      );
  PFJetTree->Branch("PFJet_chargedEmEnergy"             , & PFJet_chargedEmEnergy             , "PFJet_chargedEmEnergy/F"             );
  PFJetTree->Branch("PFJet_chargedEmEnergyFraction"     , & PFJet_chargedEmEnergyFraction     , "PFJet_chargedEmEnergyFraction/F"     );
  PFJetTree->Branch("PFJet_chargedHadronEnergy"         , & PFJet_chargedHadronEnergy         , "PFJet_chargedHadronEnergy/F"         );
  PFJetTree->Branch("PFJet_chargedHadronEnergyFraction" , & PFJet_chargedHadronEnergyFraction , "PFJet_chargedHadronEnergyFraction/F" );
  PFJetTree->Branch("PFJet_chargedHadronMultiplicity"   , & PFJet_chargedHadronMultiplicity   , "PFJet_chargedHadronMultiplicity/F"   );
  PFJetTree->Branch("PFJet_chargedMuEnergy"             , & PFJet_chargedMuEnergy             , "PFJet_chargedMuEnergy/F"             );
  PFJetTree->Branch("PFJet_chargedMuEnergyFraction"     , & PFJet_chargedMuEnergyFraction     , "PFJet_chargedMuEnergyFraction/F"     );
  PFJetTree->Branch("PFJet_chargedMultiplicity"         , & PFJet_chargedMultiplicity         , "PFJet_chargedMultiplicity/F"         );
  PFJetTree->Branch("PFJet_electronEnergy"              , & PFJet_electronEnergy              , "PFJet_electronEnergy/F"              );
  PFJetTree->Branch("PFJet_electronEnergyFraction"      , & PFJet_electronEnergyFraction      , "PFJet_electronEnergyFraction/F"      );
  PFJetTree->Branch("PFJet_electronMultiplicity"        , & PFJet_electronMultiplicity        , "PFJet_electronMultiplicity/F"        );
  PFJetTree->Branch("PFJet_HFEMEnergy"                  , & PFJet_HFEMEnergy                  , "PFJet_HFEMEnergy/F"                  );
  PFJetTree->Branch("PFJet_HFEMEnergyFraction"          , & PFJet_HFEMEnergyFraction          , "PFJet_HFEMEnergyFraction/F"          );
  PFJetTree->Branch("PFJet_HFEMMultiplicity"            , & PFJet_HFEMMultiplicity            , "PFJet_HFEMMultiplicity/F"            );
  PFJetTree->Branch("PFJet_HFHadronEnergy"              , & PFJet_HFHadronEnergy              , "PFJet_HFHadronEnergy/F"              );
  PFJetTree->Branch("PFJet_HFHadronEnergyFraction"      , & PFJet_HFHadronEnergyFraction      , "PFJet_HFHadronEnergyFraction/F"      );
  PFJetTree->Branch("PFJet_HFHadronMultiplicity"        , & PFJet_HFHadronMultiplicity        , "PFJet_HFHadronMultiplicity/F"        );
  PFJetTree->Branch("PFJet_muonEnergy"                  , & PFJet_muonEnergy                  , "PFJet_muonEnergy/F"                  );
  PFJetTree->Branch("PFJet_muonEnergyFraction"          , & PFJet_muonEnergyFraction          , "PFJet_muonEnergyFraction/F"          );
  PFJetTree->Branch("PFJet_muonMultiplicity"            , & PFJet_muonMultiplicity            , "PFJet_muonMultiplicity/F"            );
  PFJetTree->Branch("PFJet_neutralEmEnergy"             , & PFJet_neutralEmEnergy             , "PFJet_neutralEmEnergy/F"             );
  PFJetTree->Branch("PFJet_neutralEmEnergyFraction"     , & PFJet_neutralEmEnergyFraction     , "PFJet_neutralEmEnergyFraction/F"     );
  PFJetTree->Branch("PFJet_neutralHadronEnergy"         , & PFJet_neutralHadronEnergy         , "PFJet_neutralHadronEnergy/F"         );
  PFJetTree->Branch("PFJet_neutralHadronEnergyFraction" , & PFJet_neutralHadronEnergyFraction , "PFJet_neutralHadronEnergyFraction/F" );
  PFJetTree->Branch("PFJet_neutralHadronMultiplicity"   , & PFJet_neutralHadronMultiplicity   , "PFJet_neutralHadronMultiplicity/F"   );
  PFJetTree->Branch("PFJet_neutralMultiplicity"         , & PFJet_neutralMultiplicity         , "PFJet_neutralMultiplicity/F"         );
  PFJetTree->Branch("PFJet_photonEnergy"                , & PFJet_photonEnergy                , "PFJet_photonEnergy/F"                );
  PFJetTree->Branch("PFJet_photonEnergyFraction"        , & PFJet_photonEnergyFraction        , "PFJet_photonEnergyFraction/F"        );
  PFJetTree->Branch("PFJet_photonMultiplicity"          , & PFJet_photonMultiplicity          , "PFJet_photonMultiplicity/F"          );

  PFJetTree->Branch("PFJet_HadEnergy"                   , & PFJet_HadEnergy                   , "PFJet_HadEnergy/F"                   );
  PFJetTree->Branch("PFJet_HadPt"                       , & PFJet_HadPt                       , "PFJet_HadPt/F"                       );
  PFJetTree->Branch("PFJet_HadPt_Over_GenJet_HadPt"     , & PFJet_HadPt_Over_GenJet_HadPt     , "PFJet_HadPt_Over_GenJet_HadPt/F"   );
  PFJetTree->Branch("PFJet_HadE_Over_GenJet_HadE"       , & PFJet_HadE_Over_GenJet_HadE       , "PFJet_HadE_Over_GenJet_HadE/F"     );

}

// ------------ method called once each job just after ending the event loop  ------------
void 
SimpleJetResponseTreeMaker::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimpleJetResponseTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleJetResponseTreeMaker);
