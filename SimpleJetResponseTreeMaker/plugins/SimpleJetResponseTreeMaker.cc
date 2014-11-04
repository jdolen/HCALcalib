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

      TTree *JetTree;
      Int_t Run         ;
      Int_t Lumi        ;
      Int_t Event       ;
      Int_t Nvtx        ;
      Int_t NvtxGood    ;
      Float_t Rho       ;
      Float_t Mass      ;
      Float_t Energy    ;
      Float_t Pt        ;
      Float_t Eta       ;
      Float_t Phi       ;
      Float_t GenMass   ;
      Float_t GenEnergy ;
      Float_t GenPt     ;
      Float_t GenEta    ;
      Float_t GenPhi    ;
      Float_t GenDeltaR ;
      Float_t RecoOverGenPt;
      Float_t RecoMinusGenPt;
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
  
  edm::Handle<double> rhoH;
  iEvent.getByLabel(rhoSrc_, rhoH );
  double rhoVal = *rhoH;

  edm::Handle<reco::GenJetCollection> genJetH;
  iEvent.getByLabel(genJetSrc_, genJetH); 

  edm::Handle<reco::PFJetCollection> pfjetH;
  iEvent.getByLabel(recoJetSrc_, pfjetH); 

  for ( reco::PFJetCollection::const_iterator jet = pfjetH->begin(); jet != pfjetH->end(); ++jet ) {

    reco::Candidate::LorentzVector uncorrJet = jet->p4();
    jec_->setJetEta( uncorrJet.eta() );
    jec_->setJetPt ( uncorrJet.pt() );
    jec_->setJetE  ( uncorrJet.energy() );
    jec_->setJetA  ( jet->jetArea() );
    jec_->setRho   ( rhoVal );
    // jec_->setNPV   ( npv );
    
    double corr = jec_->getCorrection();
    reco::Candidate::PolarLorentzVector corrJet (uncorrJet.pt(), uncorrJet.eta(), uncorrJet.phi(), uncorrJet.mass());
    corrJet *=  (corr );

    if (corrJet.pt() < 10) continue;

    double closest_genjet_dR = 9999;
    reco::GenJet theMatchingGenJet;
    for (GenJetCollection::const_iterator genjet=genJetH->begin(); genjet!=genJetH->end(); genjet++) {
      if ( genjet->pt() < 10 || fabs(genjet->eta())>5 ) continue; 
      double deltar = deltaR( corrJet.eta(), corrJet.phi(), genjet->eta(), genjet->phi() );
      if ( deltar > 0.2 ) continue;
      if ( deltar < closest_genjet_dR ){
        closest_genjet_dR =  deltar;
        theMatchingGenJet = (*genjet);
      }
    }

    if (closest_genjet_dR<9000)
    {
      Run        = run;
      Lumi       = lumi;
      Event      = event;
      Nvtx       = count_vertex;     
      NvtxGood   = count_good_vertex;
      Rho        = rhoVal;
      Mass       = corrJet.mass();     
      Energy     = corrJet.energy();   
      Pt         = corrJet.pt();         
      Eta        = corrJet.eta();        
      Phi        = corrJet.phi();        
      GenMass    = theMatchingGenJet.mass();     
      GenEnergy  = theMatchingGenJet.energy();   
      GenPt      = theMatchingGenJet.pt();        
      GenEta     = theMatchingGenJet.eta();       
      GenPhi     = theMatchingGenJet.phi();    
      GenDeltaR  = deltaR( corrJet.eta(), corrJet.phi(), theMatchingGenJet.eta(), theMatchingGenJet.phi() );
      if (theMatchingGenJet.pt() > 0) RecoOverGenPt = corrJet.pt()/theMatchingGenJet.pt();
      RecoMinusGenPt = corrJet.pt() - theMatchingGenJet.pt();

      PFJet_Nconst                = jet->numberOfDaughters            ();  
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


}

// ------------ method called once each job just before starting event loop  ------------
void 
SimpleJetResponseTreeMaker::beginJob()
{

  edm::Service<TFileService> fs;
  TFileDirectory subDir3 = fs->mkdir("Trees");
  JetTree = new TTree("JetTree","JetTree");
  JetTree->Branch("Run"      ,  & Run,        "Run/I");
  JetTree->Branch("Lumi"     ,  & Lumi,       "Lumi/I");
  JetTree->Branch("Event"    ,  & Event,      "Event/I");
  JetTree->Branch("Nvtx"     ,  & Nvtx,       "Nvtx/I");
  JetTree->Branch("NvtxGood" ,  & NvtxGood,   "NvtxGood/I");
  JetTree->Branch("Rho"      ,  & Rho,        "Rho/F");
  JetTree->Branch("Mass"     ,  & Mass,       "Mass/F");
  JetTree->Branch("Energy"   ,  & Energy,     "Energy/F"); 
  JetTree->Branch("Pt"       ,  & Pt,         "Pt/F"); 
  JetTree->Branch("Eta"      ,  & Eta,        "Eta/F"); 
  JetTree->Branch("Phi"      ,  & Phi,        "Phi/F"); 
  JetTree->Branch("GenMass"  ,  & GenMass,    "GenMass/F");
  JetTree->Branch("GenEnergy",  & GenEnergy,  "GenEnergy/F"); 
  JetTree->Branch("GenPt"    ,  & GenPt,      "GenPt/F"); 
  JetTree->Branch("GenEta"   ,  & GenEta,     "GenEta/F"); 
  JetTree->Branch("GenPhi"   ,  & GenPhi,     "GenPhi/F"); 
  JetTree->Branch("GenDeltaR"   ,  & GenDeltaR,     "GenDeltaR/F"); 
  JetTree->Branch("RecoOverGenPt"   ,  & RecoOverGenPt,     "RecoOverGenPt/F"); 
  JetTree->Branch("RecoMinusGenPt"   ,  & RecoMinusGenPt,     "RecoMinusGenPt/F"); 

  JetTree->Branch("PFJet_Nconst"                , & PFJet_Nconst                , "PFJet_Nconst/F"                );
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
