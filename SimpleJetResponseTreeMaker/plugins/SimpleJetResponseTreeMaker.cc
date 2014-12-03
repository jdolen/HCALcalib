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
      edm::InputTag   recoJetSrc_;
      edm::InputTag   genJetSrc_;
      std::vector<std::string>  jecPayloads_; /// files for JEC payloads

      boost::shared_ptr<FactorizedJetCorrector> jec_;
      boost::shared_ptr<JetCorrectionUncertainty> jecUnc_;


      TTree *EventTree;
      Int_t     EventRun                    ;
      Int_t     EventLumi                   ;
      Int_t     EventEvent                  ;
      Int_t     EventNvtx                   ;
      Int_t     EventNvtxGood               ;
      Float_t   EventRho                    ;
      Float_t   EventRhoCalo                    ;
      Float_t   EventRhoCaloCentral                    ;

      Int_t EventNPU_BXm3  ; 
      Int_t EventNPU_BXm2  ; 
      Int_t EventNPU_BXm1  ; 
      Int_t EventNPU_BX0   ; 
      Int_t EventNPU_BX1   ; 
      Int_t EventNPU_BX2   ; 



      Float_t   EventCorrSumJetHadE               ;
      Float_t   EventSumJetPt                   ;
      Float_t   EventSumJetPt_Pt10              ;
      Float_t   EventSumJetPt_Pt20              ;
      Float_t   EventSumJetPt_Pt10Eta2p4        ;
      Float_t   EventSumJetPt_Barrel            ;
      Float_t   EventSumJetPt_Endcap             ;

//EventTree->Branch("EventSumJetPt_Barrel"     ,  & EventSumJetPt_Barrel   , "EventSumJetPt_Barrel/F"        );
//  EventTree->Branch("EventSumJetPt_Endcap"     ,  & EventSumJetPt_Endcap   , "EventSumJetPt_Endcap/F"        );


      Int_t     EventNjets                      ;
      Int_t     EventNjets_Pt10                 ;
      Int_t     EventNjets_Pt20                 ;
      Int_t     EventNjets_Pt10Eta2p4           ;
      Float_t   EventCorrSumJetPt               ;
      Float_t   EventCorrSumJetPt_Pt10          ;
      Float_t   EventCorrSumJetPt_Pt20          ;
      Float_t   EventCorrSumJetPt_Pt10Eta2p4    ;
      Int_t     EventCorrNjets                  ;
      Int_t     EventCorrNjets_Pt10             ;
      Int_t     EventCorrNjets_Pt20             ;
      Int_t     EventCorrNjets_Pt10Eta2p4       ;

      Float_t   EventCaloSumJetHadE             ;
      Float_t   EventCaloSumJetPt               ;
      Float_t   EventCaloSumJetPt_Pt10          ;
      Float_t   EventCaloSumJetPt_Pt20          ;
      Float_t   EventCaloSumJetPt_Pt10Eta2p4    ;
      Float_t   EventCaloSumJetPt_Barrel            ;
      Float_t   EventCaloSumJetPt_Endcap             ;

      Int_t     EventCaloNjets                  ;
      Int_t     EventCaloNjets_Pt10             ;
      Int_t     EventCaloNjets_Pt20             ;
      Int_t     EventCaloNjets_Pt10Eta2p4       ;

      Float_t   EventGenPart1Pt             ;
      Float_t   EventGenPart2Pt             ;



      TTree *CaloJetTree;
      Float_t MyGenJet_Mass                        ;
      Float_t MyGenJet_Energy                      ;
      Float_t MyGenJet_Pt                          ;
      Float_t MyGenJet_Eta                         ;
      Float_t MyGenJet_Phi                         ;
      Float_t MyGenJet_DeltaR                      ;

      Float_t CaloJet_Mass                         ;
      Float_t CaloJet_Energy                       ;
      Float_t CaloJet_Pt                           ;
      Float_t CaloJet_Eta                          ;
      Float_t CaloJet_Phi                          ;

      Float_t CaloJet_CorrPt                    ;
      Float_t CaloJet_CorrPtRhoArea             ;
      Float_t CaloJet_CorrPtRhocentralArea      ;
      Float_t CaloJet_Area                      ;
      Float_t CaloJet_RhoArea                   ;
      Float_t CaloJet_RhoCentralArea            ;
      Float_t CaloJet_RhoAreaOfficial           ;
                                                
      Float_t MyGenJet_HadPt                    ;                               
      Float_t MyGenJet_HadEnergy                ;                                
      Float_t CaloJet_HadEnergy                 ;                                
      Float_t CaloJet_HadPt                     ;                                
      Float_t CaloJet_HadPt_Over_GenJet_HadPt   ;                                             
      Float_t CaloJet_HadE_Over_GenJet_HadE     ;                                           

      Float_t CaloJet_Pt_Over_GenJet_Pt            ;
      Float_t CaloJet_Pt_Minus_GenJet_Pt           ;
      Float_t CaloJet_CorrPt_Over_GenJet_Pt        ;
      Float_t CaloJet_CorrPt_Minus_GenJet_Pt       ;

      Float_t CaloJet_EMEB       ;
      Float_t CaloJet_EfracHad   ;
      Float_t CaloJet_hadHB      ;
      Float_t CaloJet_hadHO      ;
      Float_t CaloJet_hadHE      ;
      Float_t CaloJet_hadHF      ;
      Float_t CaloJet_Nconst     ;

      TTree *JetTree;
      Int_t   Run                                ;
      Int_t   Lumi                               ;
      Int_t   Event                              ;  
      Int_t   Nvtx                               ;
      Int_t   NvtxGood                           ;
      Float_t Rho                                ;

      Float_t GenJet_Mass                        ;
      Float_t GenJet_Energy                      ;
      Float_t GenJet_Pt                          ;
      Float_t GenJet_Eta                         ;
      Float_t GenJet_Phi                         ;
      Float_t GenJet_DeltaR                      ;

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


      Float_t PFJet_Pt_Over_GenJet_Pt            ;
      Float_t PFJet_Pt_Minus_GenJet_Pt           ;
      Float_t PFJet_CorrPt_Over_GenJet_Pt        ;
      Float_t PFJet_CorrPt_Minus_GenJet_Pt       ;
  
      Float_t  GenJet_HadPt                       ;                       
      Float_t  GenJet_HadEnergy                   ;                       
      Float_t  PFJet_HadEnergy                    ;                      
      Float_t  PFJet_HadPt                        ;                      
      Float_t  PFJet_HadPt_Over_GenJet_HadPt      ;                                     
      Float_t  PFJet_HadE_Over_GenJet_HadE        ;                                     
                                        
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



      TTree *JetTreeMatchThresholdPt;
      Float_t Match_GenJet_Mass                        ;
      Float_t Match_GenJet_Energy                      ;
      Float_t Match_GenJet_Pt                          ;
      Float_t Match_GenJet_Eta                         ;
      Float_t Match_GenJet_Phi                         ;
      Float_t Match_GenJet_DeltaR                      ;

      Float_t Match_PFJet_Area                         ;
      Float_t Match_PFJet_CorrFactor                   ;

      Float_t Match_PFJet_Mass                         ;
      Float_t Match_PFJet_Energy                       ;
      Float_t Match_PFJet_Pt                           ;
      Float_t Match_PFJet_Eta                          ;
      Float_t Match_PFJet_Phi                          ;

      Float_t Match_PFJet_CorrMass                     ;
      Float_t Match_PFJet_CorrEnergy                   ;
      Float_t Match_PFJet_CorrPt                       ;
      Float_t Match_PFJet_CorrEta                      ;
      Float_t Match_PFJet_CorrPhi                      ;

      Float_t Match_PFJet_Pt_Over_GenJet_Pt            ;
      Float_t Match_PFJet_Pt_Minus_GenJet_Pt           ;
      Float_t Match_PFJet_CorrPt_Over_GenJet_Pt        ;
      Float_t Match_PFJet_CorrPt_Minus_GenJet_Pt       ;

      Float_t Match_PFJet_Nconst                       ; 
      Float_t Match_PFJet_chargedEmEnergy              ;
      Float_t Match_PFJet_chargedEmEnergyFraction      ;
      Float_t Match_PFJet_chargedHadronEnergy          ;
      Float_t Match_PFJet_chargedHadronEnergyFraction  ;
      Float_t Match_PFJet_chargedHadronMultiplicity    ;
      Float_t Match_PFJet_chargedMuEnergy              ;
      Float_t Match_PFJet_chargedMuEnergyFraction      ;
      Float_t Match_PFJet_chargedMultiplicity          ;
      Float_t Match_PFJet_electronEnergy               ;
      Float_t Match_PFJet_electronEnergyFraction       ;
      Float_t Match_PFJet_electronMultiplicity         ;
      Float_t Match_PFJet_HFEMEnergy                   ;
      Float_t Match_PFJet_HFEMEnergyFraction           ;
      Float_t Match_PFJet_HFEMMultiplicity             ;
      Float_t Match_PFJet_HFHadronEnergy               ;
      Float_t Match_PFJet_HFHadronEnergyFraction       ;
      Float_t Match_PFJet_HFHadronMultiplicity         ;
      Float_t Match_PFJet_muonEnergy                   ;
      Float_t Match_PFJet_muonEnergyFraction           ;
      Float_t Match_PFJet_muonMultiplicity             ;
      Float_t Match_PFJet_neutralEmEnergy              ;
      Float_t Match_PFJet_neutralEmEnergyFraction      ;
      Float_t Match_PFJet_neutralHadronEnergy          ;
      Float_t Match_PFJet_neutralHadronEnergyFraction  ;
      Float_t Match_PFJet_neutralHadronMultiplicity    ;
      Float_t Match_PFJet_neutralMultiplicity          ;
      Float_t Match_PFJet_photonEnergy                 ;
      Float_t Match_PFJet_photonEnergyFraction         ;
      Float_t Match_PFJet_photonMultiplicity           ;


};

SimpleJetResponseTreeMaker::SimpleJetResponseTreeMaker(const edm::ParameterSet& iConfig) :
pvSrc_        (iConfig.getParameter<edm::InputTag>("pvSrc") ),
rhoSrc_       (iConfig.getParameter<edm::InputTag>("rhoSrc") ),
recoJetSrc_   (iConfig.getParameter<edm::InputTag>("recoJetSrc") ),
genJetSrc_    (iConfig.getParameter<edm::InputTag>("genJetSrc") ),
jecPayloads_  (iConfig.getParameter<std::vector<std::string> >  ("jecPayloads"))
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

  int run = iEvent.id().run();
  int event = iEvent.id().event();
  int lumi = iEvent.id().luminosityBlock();

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(pvSrc_, vertices);
  int count_vertex      = 0;
  int count_good_vertex = 0;
  for(reco::VertexCollection::const_iterator v=vertices->begin();v!=vertices->end(); ++v)
  {
    count_vertex++;
    if ( v->ndof() > 4 && fabs(v->z()) < 24 && fabs(v->position().rho()) < 2 ) count_good_vertex++;
  }
    // cout<<"count_vertex "<<count_vertex<<endl;
    // cout<<"count_good_vertex "<<count_good_vertex<<endl;

  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
  std::vector<PileupSummaryInfo>::const_iterator PVI;

  float npIT_BX0  =-1.;
  float npIT_BX1  =-1.;
  float npIT_BX2  =-1.;
  float npIT_BXm1 =-1.;
  float npIT_BXm2 =-1.;
  float npIT_BXm3 =-1.;
  //int countme =0;
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
    // if(BX == -1) {
    //   //cout<<countme<<"  BX "<<BX<<" npT "<<npT<<" npIT "<<npIT<<endl;
    // }
    // if(BX == 0) {
    //   //cout<<countme<<"  BX "<<BX<<" npT "<<npT<<" npIT "<<npIT<<endl;
    //   npIT_BX0 = npIT;
    // }
    //countme++;
  }



  edm::Handle<double> rhoH;
  iEvent.getByLabel(rhoSrc_, rhoH );
  double rhoVal = *rhoH;

  // double                                "fixedGridRhoAll"           ""                "RECO"    
  // double                                "fixedGridRhoFastjetAll"    ""                "RECO"    
  // double                                "fixedGridRhoFastjetAllCalo"   ""                "RECO"    
  // double                                "fixedGridRhoFastjetCentralCalo"   ""                "RECO"    
  // double                                "fixedGridRhoFastjetCentralChargedPileUp"   ""                "RECO"    
  // double                                "fixedGridRhoFastjetCentralNeutral"   ""                "RECO"    
  edm::Handle<double> rhoCaloH;
  iEvent.getByLabel("fixedGridRhoFastjetAllCalo", rhoCaloH );
  double rhoValCalo = *rhoCaloH;

  edm::Handle<double> rhoCaloCentralH;
  iEvent.getByLabel("fixedGridRhoFastjetCentralCalo", rhoCaloCentralH );
  double rhoValCaloCentral = *rhoCaloCentralH;



  reco::Candidate::LorentzVector gen_pi1;
  reco::Candidate::LorentzVector gen_pi2;

  Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  for(size_t i = 0; i < genParticles->size(); ++ i) 
  {
    const reco::GenParticle & p = (*genParticles)[i];
    //int id = p.pdgId();
    //int st = p.status();  
    //const Candidate * mom = p.mother();
    //int momId = mom->pdgId();
    //double px = p.px(), py = p.py(), pz = p.pz(), e = p.energy();
    //double pt = p.pt();//, eta = p.eta(), phi = p.phi(), mass = p.mass();
    //double vx = p.vx(), vy = p.vy(), vz = p.vz();
    //int charge = p.charge();
    //int n = p.numberOfDaughters();
    //cout<<"id = "<<id<<" status = "<<st<<" mass = "<<mass<<" eta = "<<eta<<" phi = "<<phi<<endl;      
    //cout<<"id = "<<id<<" status = "<<st<<" pt = "<<pt<<" ndaughters "<<n<<endl;      

    if ( p.pdgId() ==  211 ) gen_pi1 = p.p4();
    if ( p.pdgId() == -211 ) gen_pi2 = p.p4();

  }
//cout<<"gen_pi1 "<<gen_pi1.pt()<<endl;

  edm::Handle<reco::GenJetCollection> genJetH;
  iEvent.getByLabel(genJetSrc_, genJetH); 
// for (GenJetCollection::const_iterator genjet=genJetH->begin(); genjet!=genJetH->end(); genjet++) {
// cout<<"genjet "<<genjet->pt()<<endl;
// }
  // reco::GenJet genjet1;
  // reco::GenJet genjet2;
  // double closest_genjet_dR1 = 9999;
  // double closest_genjet_dR2 = 9999;

  // for (GenJetCollection::const_iterator genjet=genJetH->begin(); genjet!=genJetH->end(); genjet++) {
  //     //if ( genjet->pt() < 10 || fabs(genjet->eta())>5 ) continue; 
  //     double deltar1 = deltaR( gen_pi1.eta(), gen_pi1.phi(), genjet->eta(), genjet->phi() );
  //     double deltar2 = deltaR( gen_pi2.eta(), gen_pi2.phi(), genjet->eta(), genjet->phi() );
  //     if ( deltar1 < closest_genjet_dR1 && deltar1 < 0.2  ){
  //       closest_genjet_dR1 =  deltar1;
  //       genjet1 = (*genjet);
  //     }
  //     if ( deltar2 < closest_genjet_dR2 && deltar2 < 0.2  ){
  //       closest_genjet_dR2 =  deltar2;
  //       genjet2 = (*genjet);
  //     }
  //     cout<<"Genjet pt "<<genjet->pt()<<" deltar1 "<<deltar1<<" deltar2 "<<deltar2<<endl;
  // }
  // cout<<"  closest_genjet_dR1  "<<closest_genjet_dR1<<endl;
  // cout<<"  closest_genjet_dR2  "<<closest_genjet_dR2<<endl;
  // cout<<"  gen_pi1.pt() "<<gen_pi1.pt()<<" genjet1->pt() "<<genjet1.pt() <<endl;
  // cout<<"  gen_pi2.pt() "<<gen_pi2.pt()<<" genjet2->pt() "<<genjet2.pt() <<endl;

// vector<reco::CaloJet>                 "ak4CaloJets"               ""                "HLT"     
// vector<reco::CaloJet>                 "ak5CaloJets"               ""                "HLT"     
// vector<reco::CaloJet>                 "ak7CaloJets"               ""                "HLT"     

  edm::Handle<reco::CaloJetCollection> calojetH;
  iEvent.getByLabel("ak4CaloJets", calojetH); 

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
    reco::Candidate::LorentzVector uncorrJet = jet->p4();
    //cout<<"calo "<<jet->pt()<<endl;
    //reco::Candidate::LorentzVector uncorrJet = jet->p4();

    // JEC by hand
    jec_->setJetEta( uncorrJet.eta() );
    jec_->setJetPt ( uncorrJet.pt() );
    jec_->setJetE  ( uncorrJet.energy() );
    jec_->setJetA  ( jet->jetArea() );
    jec_->setRho   ( rhoValCalo );

    // rhoValCalo
    // rhoValCaloCentral

    // double corr = jec_->getCorrection();
    // cout<<"corr "<<corr<<endl;
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

   
    // cout<<"vcor.size() "<<vcor.size()<<endl;
    // cout<<"Correction applied to jet after L1Offset = "<<vcor[0]<<endl;

    reco::Candidate::PolarLorentzVector corrJet (uncorrJet.pt(), uncorrJet.eta(), uncorrJet.phi(), uncorrJet.mass());
    corrJet *=  ( vcor[0] );


   

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

    if (corrJet.pt() < 2) continue;

    double closest_genjet_dR = 9999;
    reco::GenJet theMatchingGenJet;
    for (GenJetCollection::const_iterator genjet=genJetH->begin(); genjet!=genJetH->end(); genjet++) {
      if ( genjet->pt() < 2 || fabs(genjet->eta())>3 ) continue; 
      double deltar = deltaR( uncorrJet.eta(), uncorrJet.phi(), genjet->eta(), genjet->phi() );
      if ( deltar > 0.2 ) continue;
      if ( deltar < closest_genjet_dR ){
        closest_genjet_dR =  deltar;
        theMatchingGenJet = (*genjet);
      }
    }

    if (closest_genjet_dR<0.2)
    {

      // CaloJetTree Fill
      double genHadPt = theMatchingGenJet.hadEnergy() * theMatchingGenJet.pt() / theMatchingGenJet.energy() ;

      MyGenJet_HadPt     = genHadPt;     
      MyGenJet_HadEnergy = theMatchingGenJet.hadEnergy(); 

      MyGenJet_Mass    = theMatchingGenJet.mass();     
      MyGenJet_Energy  = theMatchingGenJet.energy();   
      MyGenJet_Pt      = theMatchingGenJet.pt();        
      MyGenJet_Eta     = theMatchingGenJet.eta();       
      MyGenJet_Phi     = theMatchingGenJet.phi();    
      MyGenJet_DeltaR  = deltaR( uncorrJet.eta(), uncorrJet.phi(), theMatchingGenJet.eta(), theMatchingGenJet.phi() );

      CaloJet_Eta        = uncorrJet.eta();        
      CaloJet_Mass       = uncorrJet.mass();     
      CaloJet_Energy     = uncorrJet.energy();   
      CaloJet_Pt         = uncorrJet.pt();         
    
      double rhoCaloArea        = jet->jetArea() * rhoValCalo ;
      double rhoCaloCentralArea = jet->jetArea() * rhoValCaloCentral ;
      // cout<<"rhoCaloArea "<<rhoCaloArea<<endl;
      // cout<<"rhoCaloCentralArea "<<rhoCaloCentralArea<<endl;

      // cout<<"uncorrJet.pt() - corrJet.pt() "<<uncorrJet.pt()-corrJet.pt()<<endl;
      // cout<<"uncorrJet.pt() "<<uncorrJet.pt()<<endl;
      // cout<<"corrJet.pt() "<<corrJet.pt()<<endl;
      // cout<<"uncorrJet.pt()-rhoCaloArea "<<uncorrJet.pt()-rhoCaloArea<<endl;
      // cout<<"uncorrJet.pt()-rhoCaloCentralArea "<<uncorrJet.pt()-rhoCaloCentralArea<<endl;


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

  edm::Handle<reco::PFJetCollection> pfjetH;
  iEvent.getByLabel(recoJetSrc_, pfjetH); 

  double sumJetHadE             = 0;

  double sumJetPt             = 0;
  double sumJetPt_Pt10        = 0;
  double sumJetPt_Pt20        = 0;
  double sumJetPt_Pt10_Eta2p4 = 0;
  double sumJetPt_Barrel =0;
  double sumJetPt_Endcap =0;
  int count_jets             = 0;
  int count_jets_Pt10        = 0;
  int count_jets_Pt20        = 0;
  int count_jets_Pt10_Eta2p4 = 0;

  double sumJetCorrPt             = 0;
  double sumJetCorrPt_Pt10        = 0;
  double sumJetCorrPt_Pt20        = 0;
  double sumJetCorrPt_Pt10_Eta2p4 = 0;

  int count_jets_corr             = 0;
  int count_jets_corr_Pt10        = 0;
  int count_jets_corr_Pt20        = 0;
  int count_jets_corr_Pt10_Eta2p4 = 0;


  for ( reco::PFJetCollection::const_iterator jet = pfjetH->begin(); jet != pfjetH->end(); ++jet ) {
  //cout<<"pf "<<jet->pt()<<endl;

    reco::Candidate::LorentzVector uncorrJet = jet->p4();
    jec_->setJetEta( uncorrJet.eta() );
    jec_->setJetPt ( uncorrJet.pt() );
    jec_->setJetE  ( uncorrJet.energy() );
    jec_->setJetA  ( jet->jetArea() );
    jec_->setRho   ( rhoVal );
    // jec_->setNPV   ( npv );
    
    double corr = jec_->getCorrection();
//cout<<"corr "<<corr<<endl;
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

//	cout<<"vcor.size() "<<vcor.size()<<endl;
  //  cout<<"Correction applied to jet after L1Offset = "<<vcor[0]<<endl;
   // cout<<"Correction applied to jet after 1  = "<<vcor[1]<<endl;
   // cout<<"Correction applied to jet after 2  = "<<vcor[2]<<endl;


    reco::Candidate::PolarLorentzVector corrJet (uncorrJet.pt(), uncorrJet.eta(), uncorrJet.phi(), uncorrJet.mass());
    corrJet *=  ( vcor[0] );

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

   if (corrJet.pt() < 2) continue;

    double closest_genjet_dR = 9999;
    reco::GenJet theMatchingGenJet;
    for (GenJetCollection::const_iterator genjet=genJetH->begin(); genjet!=genJetH->end(); genjet++) {
      if ( genjet->pt() < 2 || fabs(genjet->eta())>5 ) continue; 
      //double deltar = deltaR( corrJet.eta(), corrJet.phi(), genjet->eta(), genjet->phi() );
      double deltar = deltaR( uncorrJet.eta(), uncorrJet.phi(), genjet->eta(), genjet->phi() );
      if ( deltar > 0.2 ) continue;
      if ( deltar < closest_genjet_dR ){
        closest_genjet_dR =  deltar;
        theMatchingGenJet = (*genjet);
      }
    }

    if (closest_genjet_dR<0.2)
    {

      // JetTree Fill

      Run        = run;
      Lumi       = lumi;
      Event      = event;
      Nvtx       = count_vertex;     
      NvtxGood   = count_good_vertex;
      Rho        = rhoVal;

      double genHadPt = theMatchingGenJet.hadEnergy() * theMatchingGenJet.pt()  / theMatchingGenJet.energy() ;

      GenJet_HadPt     = genHadPt;     
      GenJet_HadEnergy = theMatchingGenJet.hadEnergy();     
      GenJet_Mass    = theMatchingGenJet.mass();     
      GenJet_Energy  = theMatchingGenJet.energy();   
      GenJet_Pt      = theMatchingGenJet.pt();        
      GenJet_Eta     = theMatchingGenJet.eta();       
      GenJet_Phi     = theMatchingGenJet.phi();    
      GenJet_DeltaR  = deltaR( uncorrJet.eta(), uncorrJet.phi(), theMatchingGenJet.eta(), theMatchingGenJet.phi() );

      PFJet_Area       = jet->jetArea() ;      
      PFJet_CorrFactor = corr;

      double hadE = jet->chargedHadronEnergy() + jet->neutralHadronEnergy();
      double hadPt = hadE * uncorrJet.pt() / uncorrJet.E();
       
      PFJet_HadEnergy = hadE;
      PFJet_HadPt     = hadPt;



      //cout<<"hadPt "<<hadPt<<" Pt "<<uncorrJet.pt()<<" HadE "<<hadE<<" E "<<uncorrJet.energy()<<endl;
    
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
      // cout<<"rhoCaloArea "<<rhoCaloArea<<endl;
      // cout<<"rhoCaloCentralArea "<<rhoCaloCentralArea<<endl;

      // cout<<"uncorrJet.pt() - corrJet.pt() "<<uncorrJet.pt()-corrJet.pt()<<endl;
      // cout<<"uncorrJet.pt() "<<uncorrJet.pt()<<endl;
      // cout<<"corrJet.pt() "<<corrJet.pt()<<endl;
      // cout<<"uncorrJet.pt()-rhoCaloArea "<<uncorrJet.pt()-rhoCaloArea<<endl;
      // cout<<"uncorrJet.pt()-rhoCaloCentralArea "<<uncorrJet.pt()-rhoCaloCentralArea<<endl;

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

      JetTree->Fill();
    }
  }

  //EventTree fill
  EventRun                        = run;
  EventLumi                       = lumi;
  EventEvent                      = event;
  EventNvtx                       = count_vertex;     
  EventNvtxGood                   = count_good_vertex;
  EventRho                        = rhoVal;
  EventRhoCalo                    = rhoValCalo;
  EventRhoCaloCentral                    = rhoValCaloCentral;

  EventNPU_BXm3   = npIT_BXm3 ;
  EventNPU_BXm2   = npIT_BXm2 ;
  EventNPU_BXm1   = npIT_BXm1 ;
  EventNPU_BX0    = npIT_BX0  ;
  EventNPU_BX1    = npIT_BX1  ;
  EventNPU_BX2    = npIT_BX2  ;


  EventSumJetPt                   = sumJetPt;
  EventSumJetPt_Pt10              = sumJetPt_Pt10;
  EventSumJetPt_Pt20              = sumJetPt_Pt20;
  EventSumJetPt_Pt10Eta2p4        = sumJetPt_Pt10_Eta2p4;
  EventSumJetPt_Barrel        = sumJetPt_Barrel;
  EventSumJetPt_Endcap        = sumJetPt_Endcap;
  EventNjets                      = count_jets;
  EventNjets_Pt10                 = count_jets_Pt10;
  EventNjets_Pt10Eta2p4           = count_jets_Pt10_Eta2p4;

  EventCorrSumJetHadE                   = sumJetHadE;
  EventCorrSumJetPt                   = sumJetCorrPt;
  EventCorrSumJetPt_Pt10              = sumJetCorrPt_Pt10;
  EventCorrSumJetPt_Pt20              = sumJetCorrPt_Pt20;
  EventCorrSumJetPt_Pt10Eta2p4        = sumJetCorrPt_Pt10_Eta2p4;
  EventCorrNjets                      = count_jets_corr;
  EventCorrNjets_Pt10                 = count_jets_corr_Pt10;
  EventCorrNjets_Pt10Eta2p4           = count_jets_corr_Pt10_Eta2p4;

  EventCaloSumJetHadE                   = calo_sumJetHadE;
  EventCaloSumJetPt                   = calo_sumJetPt;
  EventCaloSumJetPt_Pt10              = calo_sumJetPt_Pt10;
  EventCaloSumJetPt_Pt20              = calo_sumJetPt_Pt20;
  EventCaloSumJetPt_Pt10Eta2p4        = calo_sumJetPt_Pt10_Eta2p4;
  EventCaloNjets                      = calo_count_jets;
  EventCaloNjets_Pt10                 = calo_count_jets_Pt10;
  EventCaloNjets_Pt10Eta2p4           = calo_count_jets_Pt10_Eta2p4;
  EventCaloSumJetPt_Barrel        = calo_sumJetPt_Barrel;
  EventCaloSumJetPt_Endcap        = calo_sumJetPt_Endcap;


  EventGenPart1Pt                 = gen_pi1.pt();
  EventGenPart2Pt                 = gen_pi2.pt();
  EventTree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
SimpleJetResponseTreeMaker::beginJob()
{

  edm::Service<TFileService> fs;
  TFileDirectory subDir3 = fs->mkdir("Trees");

  EventTree = new TTree("EventTree","EventTree");
  EventTree->Branch("EventRun"                     ,  & EventRun                   , "EventRun/I"                        );
  EventTree->Branch("EventLumi"                    ,  & EventLumi                  , "EventLumi/I"                       );
  EventTree->Branch("EventEvent"                   ,  & EventEvent                 , "EventEvent/I"                      );
  EventTree->Branch("EventNvtx"                    ,  & EventNvtx                  , "EventNvtx/I"                       );
  EventTree->Branch("EventNvtxGood"                ,  & EventNvtxGood              , "EventNvtxGood/I"                   );
  EventTree->Branch("EventRho"                     ,  & EventRho                   , "EventRho/F"                        );
  EventTree->Branch("EventRhoCalo"                     ,  & EventRhoCalo                   , "EventRhoCalo/F"                        );
  EventTree->Branch("EventRhoCaloCentral"                     ,  & EventRhoCaloCentral                   , "EventRhoCaloCentral/F"                        );
  
  EventTree->Branch("EventNPU_BXm3"                     ,  & EventNPU_BXm3                   , "EventNPU_BXm3/I"                        );
  EventTree->Branch("EventNPU_BXm2"                     ,  & EventNPU_BXm2                   , "EventNPU_BXm2/I"                        );
  EventTree->Branch("EventNPU_BXm1"                     ,  & EventNPU_BXm1                   , "EventNPU_BXm1/I"                        );
  EventTree->Branch("EventNPU_BX0"                      ,  & EventNPU_BX0                    , "EventNPU_BX0/I"                         );
  EventTree->Branch("EventNPU_BX1"                      ,  & EventNPU_BX1                    , "EventNPU_BX1/I"                         );
  EventTree->Branch("EventNPU_BX2"                      ,  & EventNPU_BX2                    , "EventNPU_BX2/I"                         );


 
 
 
 
 
 


  EventTree->Branch("EventCorrSumJetHadE"                ,  & EventCorrSumJetHadE              , "EventCorrSumJetHadE/F"                   );
  EventTree->Branch("EventSumJetPt"                ,  & EventSumJetPt              , "EventSumJetPt/F"                   );
  EventTree->Branch("EventSumJetPt_Pt10"           ,  & EventSumJetPt_Pt10         , "EventSumJetPt_Pt10/F"              );
  EventTree->Branch("EventSumJetPt_Pt20"           ,  & EventSumJetPt_Pt20         , "EventSumJetPt_Pt20/F"              );
  EventTree->Branch("EventSumJetPt_Pt10Eta2p4"     ,  & EventSumJetPt_Pt10Eta2p4   , "EventSumJetPt_Pt10Eta2p4/F"        );
  EventTree->Branch("EventSumJetPt_Barrel"     ,  & EventSumJetPt_Barrel   , "EventSumJetPt_Barrel/F"        );
  EventTree->Branch("EventSumJetPt_Endcap"     ,  & EventSumJetPt_Endcap   , "EventSumJetPt_Endcap/F"        );
  EventTree->Branch("EventNjets"                   ,  & EventNjets                 , "EventNjets/I"                      );
  EventTree->Branch("EventNjets_Pt10"              ,  & EventNjets_Pt10            , "EventNjets_Pt10/I"                 );
  EventTree->Branch("EventNjets_Pt20"              ,  & EventNjets_Pt20            , "EventNjets_Pt20/I"                 );
  EventTree->Branch("EventNjets_Pt10Eta2p4"        ,  & EventNjets_Pt10Eta2p4      , "EventNjets_Pt10Eta2p4/I"           );

  EventTree->Branch("EventCorrSumJetPt"                ,  & EventCorrSumJetPt              , "EventCorrSumJetPt/F"                   );
  EventTree->Branch("EventCorrSumJetPt_Pt10"           ,  & EventCorrSumJetPt_Pt10         , "EventCorrSumJetPt_Pt10/F"              );
  EventTree->Branch("EventCorrSumJetPt_Pt20"           ,  & EventCorrSumJetPt_Pt20         , "EventCorrSumJetPt_Pt20/F"              );
  EventTree->Branch("EventCorrSumJetPt_Pt10Eta2p4"     ,  & EventCorrSumJetPt_Pt10Eta2p4   , "EventCorrSumJetPt_Pt10Eta2p4/F"        );
  EventTree->Branch("EventCorrNjets"                   ,  & EventCorrNjets                 , "EventCorrNjets/I"                      );
  EventTree->Branch("EventCorrNjets_Pt10"              ,  & EventCorrNjets_Pt10            , "EventCorrNjets_Pt10/I"                 );
  EventTree->Branch("EventCorrNjets_Pt20"              ,  & EventCorrNjets_Pt20            , "EventCorrNjets_Pt20/I"                 );
  EventTree->Branch("EventCorrNjets_Pt10Eta2p4"        ,  & EventCorrNjets_Pt10Eta2p4      , "EventCorrNjets_Pt10Eta2p4/I"           );


  EventTree->Branch("EventCaloSumJetHadE"                ,  & EventCaloSumJetHadE              , "EventCaloSumJetHadE/F"                   );
  EventTree->Branch("EventCaloSumJetPt"                ,  & EventCaloSumJetPt              , "EventCaloSumJetPt/F"                   );
  EventTree->Branch("EventCaloSumJetPt_Pt10"           ,  & EventCaloSumJetPt_Pt10         , "EventCaloSumJetPt_Pt10/F"              );
  EventTree->Branch("EventCaloSumJetPt_Pt20"           ,  & EventCaloSumJetPt_Pt20         , "EventCaloSumJetPt_Pt20/F"              );
  EventTree->Branch("EventCaloSumJetPt_Pt10Eta2p4"     ,  & EventCaloSumJetPt_Pt10Eta2p4   , "EventCaloSumJetPt_Pt10Eta2p4/F"        );
  EventTree->Branch("EventCaloSumJetPt_Barrel"     ,  & EventCaloSumJetPt_Barrel   , "EventCaloSumJetPt_Barrel/F"        );
  EventTree->Branch("EventCaloSumJetPt_Endcap"     ,  & EventCaloSumJetPt_Endcap   , "EventCaloSumJetPt_Endcap/F"        );
  EventTree->Branch("EventCaloNjets"                   ,  & EventCaloNjets                 , "EventCaloNjets/I"                      );
  EventTree->Branch("EventCaloNjets_Pt10"              ,  & EventCaloNjets_Pt10            , "EventCaloNjets_Pt10/I"                 );
  EventTree->Branch("EventCaloNjets_Pt20"              ,  & EventCaloNjets_Pt20            , "EventCaloNjets_Pt20/I"                 );
  EventTree->Branch("EventCaloNjets_Pt10Eta2p4"        ,  & EventCaloNjets_Pt10Eta2p4      , "EventCaloNjets_Pt10Eta2p4/I"           );


  EventTree->Branch("EventGenPart1Pt"              ,  & EventGenPart1Pt            , "EventGenPart1Pt/F"                 );
  EventTree->Branch("EventGenPart2Pt"              ,  & EventGenPart2Pt            , "EventGenPart2Pt/F"                 );


  CaloJetTree = new TTree("CaloJetTree","CaloJetTree");

  CaloJetTree->Branch("MyGenJet_Mass"                       , & MyGenJet_Mass                       , "MyGenJet_Mass/F"                        );
  CaloJetTree->Branch("MyGenJet_Energy"                     , & MyGenJet_Energy                     , "MyGenJet_Energy/F"                      ); 
  CaloJetTree->Branch("MyGenJet_Pt"                         , & MyGenJet_Pt                         , "MyGenJet_Pt/F"                          ); 
  CaloJetTree->Branch("MyGenJet_Eta"                        , & MyGenJet_Eta                        , "MyGenJet_Eta/F"                         ); 
  CaloJetTree->Branch("MyGenJet_Phi"                        , & MyGenJet_Phi                        , "MyGenJet_Phi/F"                         ); 
  CaloJetTree->Branch("MyGenJet_DeltaR"                     , & MyGenJet_DeltaR                     , "MyGenJet_DeltaR/F"                      ); 

  CaloJetTree->Branch("CaloJet_Mass"                        , & CaloJet_Mass                        , "CaloJet_Mass/F"                        );
  CaloJetTree->Branch("CaloJet_Energy"                      , & CaloJet_Energy                      , "CaloJet_Energy/F"                      ); 
  CaloJetTree->Branch("CaloJet_Pt"                          , & CaloJet_Pt                          , "CaloJet_Pt/F"                          ); 
  CaloJetTree->Branch("CaloJet_Eta"                         , & CaloJet_Eta                         , "CaloJet_Eta/F"                         ); 
  CaloJetTree->Branch("CaloJet_Phi"                         , & CaloJet_Phi                         , "CaloJet_Phi/F"                         ); 

  CaloJetTree->Branch("CaloJet_CorrPt"                      , & CaloJet_CorrPt                      , "CaloJet_CorrPt/F"                       ); 
  CaloJetTree->Branch("CaloJet_CorrPtRhoArea"               , & CaloJet_CorrPtRhoArea               , "CaloJet_CorrPtRhoArea/F"                ); 
  CaloJetTree->Branch("CaloJet_CorrPtRhocentralArea"        , & CaloJet_CorrPtRhocentralArea        , "CaloJet_CorrPtRhocentralArea/F"         ); 
  CaloJetTree->Branch("CaloJet_Area"                        , & CaloJet_Area                        , "CaloJet_Area/F"                         ); 
  CaloJetTree->Branch("CaloJet_RhoArea"                     , & CaloJet_RhoArea                     , "CaloJet_RhoArea/F"                      ); 
  CaloJetTree->Branch("CaloJet_RhoCentralArea"              , & CaloJet_RhoCentralArea              , "CaloJet_RhoCentralArea/F"               ); 
  CaloJetTree->Branch("CaloJet_RhoAreaOfficial"             , & CaloJet_RhoAreaOfficial             , "CaloJet_RhoAreaOfficial/F"              ); 
   
   
   
   
   
   

  CaloJetTree->Branch("CaloJet_Pt_Over_GenJet_Pt"           , & CaloJet_Pt_Over_GenJet_Pt           , "CaloJet_Pt_Over_GenJet_Pt/F"           ); 
  CaloJetTree->Branch("CaloJet_Pt_Minus_GenJet_Pt"          , & CaloJet_Pt_Minus_GenJet_Pt          , "CaloJet_Pt_Minus_GenJet_Pt/F"          ); 

  CaloJetTree->Branch("CaloJet_EMEB"                          , & CaloJet_EMEB                      , "CaloJet_EMEB/F"                          );
  CaloJetTree->Branch("CaloJet_EfracHad"                      , & CaloJet_EfracHad                  , "CaloJet_EfracHad/F"                      ); 
  CaloJetTree->Branch("CaloJet_hadHB"                         , & CaloJet_hadHB                     , "CaloJet_hadHB/F"                         ); 
  CaloJetTree->Branch("CaloJet_hadHO"                         , & CaloJet_hadHO                     , "CaloJet_hadHO/F"                         ); 
  CaloJetTree->Branch("CaloJet_hadHE"                         , & CaloJet_hadHE                     , "CaloJet_hadHE/F"                         ); 
  CaloJetTree->Branch("CaloJet_hadHF"                         , & CaloJet_hadHF                     , "CaloJet_hadHF/F"                         ); 
  CaloJetTree->Branch("CaloJet_Nconst"                        , & CaloJet_Nconst                    , "CaloJet_Nconst/F"                        ); 


  CaloJetTree->Branch("MyGenJet_HadPt"                      , & MyGenJet_HadPt                                    , "MyGenJet_HadPt/F"                        ); 
  CaloJetTree->Branch("MyGenJet_HadEnergy"                  , & MyGenJet_HadEnergy                                , "MyGenJet_HadEnergy/F"                        ); 
  CaloJetTree->Branch("CaloJet_HadEnergy"                   , & CaloJet_HadEnergy                                 , "CaloJet_HadEnergy/F"                        ); 
  CaloJetTree->Branch("CaloJet_HadPt"                       , & CaloJet_HadPt                                     , "CaloJet_HadPt/F"                        ); 
  CaloJetTree->Branch("CaloJet_HadPt_Over_GenJet_HadPt"     , & CaloJet_HadPt_Over_GenJet_HadPt                   , "CaloJet_HadPt_Over_GenJet_HadPt/F"                        ); 
  CaloJetTree->Branch("CaloJet_HadE_Over_GenJet_HadE"       , & CaloJet_HadE_Over_GenJet_HadE                     , "CaloJet_HadE_Over_GenJet_HadE/F"                        ); 


  
  
  
  
  
  


  JetTree = new TTree("JetTree","JetTree");
  JetTree->Branch("Run"      ,  & Run,        "Run/I");
  JetTree->Branch("Lumi"     ,  & Lumi,       "Lumi/I");
  JetTree->Branch("Event"    ,  & Event,      "Event/I");
  JetTree->Branch("Nvtx"     ,  & Nvtx,       "Nvtx/I");
  JetTree->Branch("NvtxGood" ,  & NvtxGood,   "NvtxGood/I");
  JetTree->Branch("Rho"      ,  & Rho,        "Rho/F");
 
  JetTree->Branch("GenJet_Mass"                       , & GenJet_Mass                       , "GenJet_Mass/F"                        );
  JetTree->Branch("GenJet_Energy"                     , & GenJet_Energy                     , "GenJet_Energy/F"                      ); 
  JetTree->Branch("GenJet_Pt"                         , & GenJet_Pt                         , "GenJet_Pt/F"                          ); 
  JetTree->Branch("GenJet_Eta"                        , & GenJet_Eta                        , "GenJet_Eta/F"                         ); 
  JetTree->Branch("GenJet_Phi"                        , & GenJet_Phi                        , "GenJet_Phi/F"                         ); 
  JetTree->Branch("GenJet_DeltaR"                     , & GenJet_DeltaR                     , "GenJet_DeltaR/F"                      ); 

  JetTree->Branch("PFJet_Area"                        , & PFJet_Area                        , "PFJet_Area/F"                        ); 
  JetTree->Branch("PFJet_CorrFactor"                  , & PFJet_CorrFactor                  , "PFJet_CorrFactor/F"                  ); 

  JetTree->Branch("PFJet_Mass"                        , & PFJet_Mass                        , "PFJet_Mass/F"                        );
  JetTree->Branch("PFJet_Energy"                      , & PFJet_Energy                      , "PFJet_Energy/F"                      ); 
  JetTree->Branch("PFJet_Pt"                          , & PFJet_Pt                          , "PFJet_Pt/F"                          ); 
  JetTree->Branch("PFJet_Eta"                         , & PFJet_Eta                         , "PFJet_Eta/F"                         ); 
  JetTree->Branch("PFJet_Phi"                         , & PFJet_Phi                         , "PFJet_Phi/F"                         ); 

  JetTree->Branch("PFJet_CorrMass"                    , & PFJet_CorrMass                    , "PFJet_CorrMass/F"                    ); 
  JetTree->Branch("PFJet_CorrEnergy"                  , & PFJet_CorrEnergy                  , "PFJet_CorrEnergy/F"                  ); 
  JetTree->Branch("PFJet_CorrPt"                      , & PFJet_CorrPt                      , "PFJet_CorrPt/F"                      ); 
  JetTree->Branch("PFJet_CorrEta"                     , & PFJet_CorrEta                     , "PFJet_CorrEta/F"                     ); 
  JetTree->Branch("PFJet_CorrPhi"                     , & PFJet_CorrPhi                     , "PFJet_CorrPhi/F"                     ); 


  JetTree->Branch("PFJet_CorrPtRhoArea"                      , & PFJet_CorrPtRhoArea                      , "PFJet_CorrPtRhoArea/F"                      ); 
  JetTree->Branch("PFJet_RhoArea"                      , & PFJet_RhoArea                      , "PFJet_RhoArea/F"                      ); 
  JetTree->Branch("PFJet_RhoAreaOfficial"                      , & PFJet_RhoAreaOfficial                      , "PFJet_RhoAreaOfficial/F"                      ); 



  JetTree->Branch("PFJet_Pt_Over_GenJet_Pt"           , & PFJet_Pt_Over_GenJet_Pt           , "PFJet_Pt_Over_GenJet_Pt/F"           ); 
  JetTree->Branch("PFJet_Pt_Minus_GenJet_Pt"          , & PFJet_Pt_Minus_GenJet_Pt          , "PFJet_Pt_Minus_GenJet_Pt/F"          ); 
  JetTree->Branch("PFJet_CorrPt_Over_GenJet_Pt"       , & PFJet_CorrPt_Over_GenJet_Pt       , "PFJet_CorrPt_Over_GenJet_Pt/F"       ); 
  JetTree->Branch("PFJet_CorrPt_Minus_GenJet_Pt"      , & PFJet_CorrPt_Minus_GenJet_Pt      , "PFJet_CorrPt_Minus_GenJet_Pt/F"      ); 

  JetTree->Branch("PFJet_Nconst"                      , & PFJet_Nconst                      , "PFJet_Nconst/F"                      );
  JetTree->Branch("PFJet_chargedEmEnergy"             , & PFJet_chargedEmEnergy             , "PFJet_chargedEmEnergy/F"             );
  JetTree->Branch("PFJet_chargedEmEnergyFraction"     , & PFJet_chargedEmEnergyFraction     , "PFJet_chargedEmEnergyFraction/F"     );
  JetTree->Branch("PFJet_chargedHadronEnergy"         , & PFJet_chargedHadronEnergy         , "PFJet_chargedHadronEnergy/F"         );
  JetTree->Branch("PFJet_chargedHadronEnergyFraction" , & PFJet_chargedHadronEnergyFraction , "PFJet_chargedHadronEnergyFraction/F" );
  JetTree->Branch("PFJet_chargedHadronMultiplicity"   , & PFJet_chargedHadronMultiplicity   , "PFJet_chargedHadronMultiplicity/F"   );
  JetTree->Branch("PFJet_chargedMuEnergy"             , & PFJet_chargedMuEnergy             , "PFJet_chargedMuEnergy/F"             );
  JetTree->Branch("PFJet_chargedMuEnergyFraction"     , & PFJet_chargedMuEnergyFraction     , "PFJet_chargedMuEnergyFraction/F"     );
  JetTree->Branch("PFJet_chargedMultiplicity"         , & PFJet_chargedMultiplicity         , "PFJet_chargedMultiplicity/F"         );
  JetTree->Branch("PFJet_electronEnergy"              , & PFJet_electronEnergy              , "PFJet_electronEnergy/F"              );
  JetTree->Branch("PFJet_electronEnergyFraction"      , & PFJet_electronEnergyFraction      , "PFJet_electronEnergyFraction/F"      );
  JetTree->Branch("PFJet_electronMultiplicity"        , & PFJet_electronMultiplicity        , "PFJet_electronMultiplicity/F"        );
  JetTree->Branch("PFJet_HFEMEnergy"                  , & PFJet_HFEMEnergy                  , "PFJet_HFEMEnergy/F"                  );
  JetTree->Branch("PFJet_HFEMEnergyFraction"          , & PFJet_HFEMEnergyFraction          , "PFJet_HFEMEnergyFraction/F"          );
  JetTree->Branch("PFJet_HFEMMultiplicity"            , & PFJet_HFEMMultiplicity            , "PFJet_HFEMMultiplicity/F"            );
  JetTree->Branch("PFJet_HFHadronEnergy"              , & PFJet_HFHadronEnergy              , "PFJet_HFHadronEnergy/F"              );
  JetTree->Branch("PFJet_HFHadronEnergyFraction"      , & PFJet_HFHadronEnergyFraction      , "PFJet_HFHadronEnergyFraction/F"      );
  JetTree->Branch("PFJet_HFHadronMultiplicity"        , & PFJet_HFHadronMultiplicity        , "PFJet_HFHadronMultiplicity/F"        );
  JetTree->Branch("PFJet_muonEnergy"                  , & PFJet_muonEnergy                  , "PFJet_muonEnergy/F"                  );
  JetTree->Branch("PFJet_muonEnergyFraction"          , & PFJet_muonEnergyFraction          , "PFJet_muonEnergyFraction/F"          );
  JetTree->Branch("PFJet_muonMultiplicity"            , & PFJet_muonMultiplicity            , "PFJet_muonMultiplicity/F"            );
  JetTree->Branch("PFJet_neutralEmEnergy"             , & PFJet_neutralEmEnergy             , "PFJet_neutralEmEnergy/F"             );
  JetTree->Branch("PFJet_neutralEmEnergyFraction"     , & PFJet_neutralEmEnergyFraction     , "PFJet_neutralEmEnergyFraction/F"     );
  JetTree->Branch("PFJet_neutralHadronEnergy"         , & PFJet_neutralHadronEnergy         , "PFJet_neutralHadronEnergy/F"         );
  JetTree->Branch("PFJet_neutralHadronEnergyFraction" , & PFJet_neutralHadronEnergyFraction , "PFJet_neutralHadronEnergyFraction/F" );
  JetTree->Branch("PFJet_neutralHadronMultiplicity"   , & PFJet_neutralHadronMultiplicity   , "PFJet_neutralHadronMultiplicity/F"   );
  JetTree->Branch("PFJet_neutralMultiplicity"         , & PFJet_neutralMultiplicity         , "PFJet_neutralMultiplicity/F"         );
  JetTree->Branch("PFJet_photonEnergy"                , & PFJet_photonEnergy                , "PFJet_photonEnergy/F"                );
  JetTree->Branch("PFJet_photonEnergyFraction"        , & PFJet_photonEnergyFraction        , "PFJet_photonEnergyFraction/F"        );
  JetTree->Branch("PFJet_photonMultiplicity"          , & PFJet_photonMultiplicity          , "PFJet_photonMultiplicity/F"          );

  JetTree->Branch("GenJet_HadPt"          , & GenJet_HadPt          , "GenJet_HadPt/F"          );
  JetTree->Branch("GenJet_HadEnergy"          , & GenJet_HadEnergy          , "GenJet_HadEnergy/F"          );
  JetTree->Branch("PFJet_HadEnergy"          , & PFJet_HadEnergy          , "PFJet_HadEnergy/F"          );
  JetTree->Branch("PFJet_HadPt"          , & PFJet_HadPt          , "PFJet_HadPt/F"          );
  JetTree->Branch("PFJet_HadPt_Over_GenJet_HadPt"          , & PFJet_HadPt_Over_GenJet_HadPt          , "PFJet_HadPt_Over_GenJet_HadPt/F"          );
  JetTree->Branch("PFJet_HadE_Over_GenJet_HadE"          , & PFJet_HadE_Over_GenJet_HadE          , "PFJet_HadE_Over_GenJet_HadE/F"          );



  JetTreeMatchThresholdPt = new TTree("JetTreeMatchThresholdPt","JetTreeMatchThresholdPt");

  JetTreeMatchThresholdPt->Branch("Match_GenJet_Mass"                       , & Match_GenJet_Mass                       , "Match_GenJet_Mass/F"                        );
  JetTreeMatchThresholdPt->Branch("Match_GenJet_Energy"                     , & Match_GenJet_Energy                     , "Match_GenJet_Energy/F"                      ); 
  JetTreeMatchThresholdPt->Branch("Match_GenJet_Pt"                         , & Match_GenJet_Pt                         , "Match_GenJet_Pt/F"                          ); 
  JetTreeMatchThresholdPt->Branch("Match_GenJet_Eta"                        , & Match_GenJet_Eta                        , "Match_GenJet_Eta/F"                         ); 
  JetTreeMatchThresholdPt->Branch("Match_GenJet_Phi"                        , & Match_GenJet_Phi                        , "Match_GenJet_Phi/F"                         ); 
  JetTreeMatchThresholdPt->Branch("Match_GenJet_DeltaR"                     , & Match_GenJet_DeltaR                     , "Match_GenJet_DeltaR/F"                      ); 

  JetTreeMatchThresholdPt->Branch("Match_PFJet_Area"                        , & Match_PFJet_Area                        , "Match_PFJet_Area/F"                        ); 
  JetTreeMatchThresholdPt->Branch("Match_PFJet_CorrFactor"                  , & Match_PFJet_CorrFactor                  , "Match_PFJet_CorrFactor/F"                  ); 

  JetTreeMatchThresholdPt->Branch("Match_PFJet_Mass"                        , & Match_PFJet_Mass                        , "Match_PFJet_Mass/F"                        );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_Energy"                      , & Match_PFJet_Energy                      , "Match_PFJet_Energy/F"                      ); 
  JetTreeMatchThresholdPt->Branch("Match_PFJet_Pt"                          , & Match_PFJet_Pt                          , "Match_PFJet_Pt/F"                          ); 
  JetTreeMatchThresholdPt->Branch("Match_PFJet_Eta"                         , & Match_PFJet_Eta                         , "Match_PFJet_Eta/F"                         ); 
  JetTreeMatchThresholdPt->Branch("Match_PFJet_Phi"                         , & Match_PFJet_Phi                         , "Match_PFJet_Phi/F"                         ); 

  JetTreeMatchThresholdPt->Branch("Match_PFJet_CorrMass"                    , & Match_PFJet_CorrMass                    , "Match_PFJet_CorrMass/F"                    ); 
  JetTreeMatchThresholdPt->Branch("Match_PFJet_CorrEnergy"                  , & Match_PFJet_CorrEnergy                  , "Match_PFJet_CorrEnergy/F"                  ); 
  JetTreeMatchThresholdPt->Branch("Match_PFJet_CorrPt"                      , & Match_PFJet_CorrPt                      , "Match_PFJet_CorrPt/F"                      ); 
  JetTreeMatchThresholdPt->Branch("Match_PFJet_CorrEta"                     , & Match_PFJet_CorrEta                     , "Match_PFJet_CorrEta/F"                     ); 
  JetTreeMatchThresholdPt->Branch("Match_PFJet_CorrPhi"                     , & Match_PFJet_CorrPhi                     , "Match_PFJet_CorrPhi/F"                     ); 

  JetTreeMatchThresholdPt->Branch("Match_PFJet_Pt_Over_GenJet_Pt"           , & Match_PFJet_Pt_Over_GenJet_Pt           , "Match_PFJet_Pt_Over_GenJet_Pt/F"           ); 
  JetTreeMatchThresholdPt->Branch("Match_PFJet_Pt_Minus_GenJet_Pt"          , & Match_PFJet_Pt_Minus_GenJet_Pt          , "Match_PFJet_Pt_Minus_GenJet_Pt/F"          ); 
  JetTreeMatchThresholdPt->Branch("Match_PFJet_CorrPt_Over_GenJet_Pt"       , & Match_PFJet_CorrPt_Over_GenJet_Pt       , "Match_PFJet_CorrPt_Over_GenJet_Pt/F"       ); 
  JetTreeMatchThresholdPt->Branch("Match_PFJet_CorrPt_Minus_GenJet_Pt"      , & Match_PFJet_CorrPt_Minus_GenJet_Pt      , "Match_PFJet_CorrPt_Minus_GenJet_Pt/F"      ); 

  JetTreeMatchThresholdPt->Branch("Match_PFJet_Nconst"                      , & Match_PFJet_Nconst                      , "Match_PFJet_Nconst/F"                      );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_chargedEmEnergy"             , & Match_PFJet_chargedEmEnergy             , "Match_PFJet_chargedEmEnergy/F"             );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_chargedEmEnergyFraction"     , & Match_PFJet_chargedEmEnergyFraction     , "Match_PFJet_chargedEmEnergyFraction/F"     );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_chargedHadronEnergy"         , & Match_PFJet_chargedHadronEnergy         , "Match_PFJet_chargedHadronEnergy/F"         );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_chargedHadronEnergyFraction" , & Match_PFJet_chargedHadronEnergyFraction , "Match_PFJet_chargedHadronEnergyFraction/F" );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_chargedHadronMultiplicity"   , & Match_PFJet_chargedHadronMultiplicity   , "Match_PFJet_chargedHadronMultiplicity/F"   );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_chargedMuEnergy"             , & Match_PFJet_chargedMuEnergy             , "Match_PFJet_chargedMuEnergy/F"             );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_chargedMuEnergyFraction"     , & Match_PFJet_chargedMuEnergyFraction     , "Match_PFJet_chargedMuEnergyFraction/F"     );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_chargedMultiplicity"         , & Match_PFJet_chargedMultiplicity         , "Match_PFJet_chargedMultiplicity/F"         );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_electronEnergy"              , & Match_PFJet_electronEnergy              , "Match_PFJet_electronEnergy/F"              );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_electronEnergyFraction"      , & Match_PFJet_electronEnergyFraction      , "Match_PFJet_electronEnergyFraction/F"      );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_electronMultiplicity"        , & Match_PFJet_electronMultiplicity        , "Match_PFJet_electronMultiplicity/F"        );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_HFEMEnergy"                  , & Match_PFJet_HFEMEnergy                  , "Match_PFJet_HFEMEnergy/F"                  );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_HFEMEnergyFraction"          , & Match_PFJet_HFEMEnergyFraction          , "Match_PFJet_HFEMEnergyFraction/F"          );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_HFEMMultiplicity"            , & Match_PFJet_HFEMMultiplicity            , "Match_PFJet_HFEMMultiplicity/F"            );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_HFHadronEnergy"              , & Match_PFJet_HFHadronEnergy              , "Match_PFJet_HFHadronEnergy/F"              );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_HFHadronEnergyFraction"      , & Match_PFJet_HFHadronEnergyFraction      , "Match_PFJet_HFHadronEnergyFraction/F"      );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_HFHadronMultiplicity"        , & Match_PFJet_HFHadronMultiplicity        , "Match_PFJet_HFHadronMultiplicity/F"        );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_muonEnergy"                  , & Match_PFJet_muonEnergy                  , "Match_PFJet_muonEnergy/F"                  );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_muonEnergyFraction"          , & Match_PFJet_muonEnergyFraction          , "Match_PFJet_muonEnergyFraction/F"          );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_muonMultiplicity"            , & Match_PFJet_muonMultiplicity            , "Match_PFJet_muonMultiplicity/F"            );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_neutralEmEnergy"             , & Match_PFJet_neutralEmEnergy             , "Match_PFJet_neutralEmEnergy/F"             );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_neutralEmEnergyFraction"     , & Match_PFJet_neutralEmEnergyFraction     , "Match_PFJet_neutralEmEnergyFraction/F"     );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_neutralHadronEnergy"         , & Match_PFJet_neutralHadronEnergy         , "Match_PFJet_neutralHadronEnergy/F"         );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_neutralHadronEnergyFraction" , & Match_PFJet_neutralHadronEnergyFraction , "Match_PFJet_neutralHadronEnergyFraction/F" );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_neutralHadronMultiplicity"   , & Match_PFJet_neutralHadronMultiplicity   , "Match_PFJet_neutralHadronMultiplicity/F"   );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_neutralMultiplicity"         , & Match_PFJet_neutralMultiplicity         , "Match_PFJet_neutralMultiplicity/F"         );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_photonEnergy"                , & Match_PFJet_photonEnergy                , "Match_PFJet_photonEnergy/F"                );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_photonEnergyFraction"        , & Match_PFJet_photonEnergyFraction        , "Match_PFJet_photonEnergyFraction/F"        );
  JetTreeMatchThresholdPt->Branch("Match_PFJet_photonMultiplicity"          , & Match_PFJet_photonMultiplicity          , "Match_PFJet_photonMultiplicity/F"          );



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
