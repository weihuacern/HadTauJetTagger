#include <memory>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "TLorentzVector.h"

class prodPFCands : public edm::EDFilter 
{
 public:
  explicit prodPFCands(const edm::ParameterSet & iConfig);
  ~prodPFCands();
 private:
  virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);
  bool isData_;
  bool debug_;

  edm::InputTag pfCandSrc_;
  edm::Handle<pat::PackedCandidateCollection> pfCandHandle_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandTok_;

  //double GetTrackActivity(const edm::Handle<pat::PackedCandidateCollection> & other_pfcands, const pat::PackedCandidate* track);
};

prodPFCands::prodPFCands(const edm::ParameterSet & iConfig)
{
  //input variables from python file
  isData_ = true;
  debug_ = iConfig.getParameter<bool>("debug");

  pfCandSrc_ = iConfig.getParameter<edm::InputTag>("pfCandSrc");
  pfCandTok_ = consumes<pat::PackedCandidateCollection> (pfCandSrc_);

  //genDecayLVec_Src_ = iConfig.getParameter<edm::InputTag>("genDecayLVec");
  //GenDecayLVec_Tok_=consumes<std::vector<TLorentzVector> > (genDecayLVec_Src_);

  //output variables
  //pf candidate, 
  //the momentum 4-vector, accessible through the usual methods like pt(), eta(), phi(), mass(), energy(), px(), p4(), polarP4(). An additional method phiAtVtx() is provided, returning the phi of the candidate's track at the vertex; this is identical to phi() for the vast majority of the particles, but the two might differ for some of them if the calorimeters had contributed significantly in defining the 4-vector of the particle.
  //the particle charge and pdgId: 11, 13, 22 for ele/mu/gamma, 211 for charged hadrons, 130 for neutral hadrons, 1 and 2 for hadronic and em particles in HF.
  //quality flags:
  //isolation: nominal, mini
  produces<std::vector<TLorentzVector> >("pfCandsLVec");
  
  /*
  produces<std::vector<TLorentzVector> >("looseisoTrksLVec");
  produces<std::vector<double> >("looseisoTrkscharge");
  produces<std::vector<double> >("looseisoTrksdz");
  produces<std::vector<int> >("looseisoTrkspdgId");
  produces<std::vector<int> >("looseisoTrksidx");
  produces<std::vector<double> >("looseisoTrksiso");
  produces<std::vector<double> >("looseisoTrksmtw");
  produces<std::vector<double> >("looseisoTrkspfActivity");
  */
}

prodPFCands::~prodPFCands() 
{
}

bool prodPFCands::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  if( !iEvent.isRealData() ) isData_ = false;

  std::auto_ptr<std::vector<TLorentzVector> > pfCandsLVec(new std::vector<TLorentzVector>());

  /*
  std::auto_ptr<std::vector<TLorentzVector> > loose_isoTrksLVec(new std::vector<TLorentzVector>());
  std::auto_ptr<std::vector<double> > loose_isoTrks_charge(new std::vector<double>());
  std::auto_ptr<std::vector<double> > loose_isoTrks_dz(new std::vector<double>());
  std::auto_ptr<std::vector<int> > loose_isoTrks_pdgId(new std::vector<int>());
  if( !isData_ )
  {
  }
  */
  /*
  edm::Handle< std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(VtxTok_, vertices);
  //reco::Vertex::Point vtxpos = (vertices->size() > 0 ? (*vertices)[0].position() : reco::Vertex::Point());
  edm::Handle<edm::View<reco::MET> > met;
  iEvent.getByToken(MetTok_, met);
  */

  iEvent.getByToken(pfCandTok_, pfCandHandle_);
  if( pfCandHandle_.isValid() )//if pf candidate is validated
  {
    for(unsigned int ip=0; ip<pfCandHandle_->size(); ip++)
    {
      if( std::isnan((*pfCandHandle_)[ip].pt()) || std::isinf((*pfCandHandle_)[ip].pt()) ) continue;//get rid of negative and infinite pt in the pf candidate 

      TLorentzVector perLVec; perLVec.SetPtEtaPhiE( (*pfCandHandle_)[ip].pt(), (*pfCandHandle_)[ip].eta(), (*pfCandHandle_)[ip].phi(), (*pfCandHandle_)[ip].energy() );
      //int perCharge = pfCandHandle_->at(ip).charge(); 
      //if( perCharge ==0 ) continue;
      //double dz = (*pfCandHandle_)[ip].dz();
      //if( fabs(dz) > isotrk_dz_ ) continue;
      //double pfActivity = GetTrackActivity(pfCandHandle_, &(*pfCandHandle_)[ip]);

      pfCandsLVec->push_back(perLVec);
      //trksForIsoVeto_charge->push_back(perCharge);
      //trksForIsoVeto_dz->push_back(dz);
      //trksForIsoVeto_pdgId->push_back(pfCandHandle_->at(ip).pdgId());
      //trksForIsoVeto_idx->push_back(ip);
      //trksForIsoVeto_iso->push_back(perIso);
      //trksForIsoVeto_pfActivity->push_back(pfActivity);
    }//loop over all the pf candidates
  }//if pf candidate is validated

  if( debug_ )
  {
    //put debug infomation here!
  }
  
  iEvent.put(pfCandsLVec, "pfCandsLVec");

  /*
  iEvent.put(loose_isoTrksLVec, "looseisoTrksLVec");
  iEvent.put(loose_isoTrks_charge, "looseisoTrkscharge");
  iEvent.put(loose_isoTrks_dz, "looseisoTrksdz");
  iEvent.put(loose_isoTrks_pdgId, "looseisoTrkspdgId");
  iEvent.put(loose_isoTrks_idx, "looseisoTrksidx");
  iEvent.put(loose_isoTrks_iso, "looseisoTrksiso");
  iEvent.put(loose_isoTrks_mtw, "looseisoTrksmtw");
  iEvent.put(loose_isoTrks_pfActivity, "looseisoTrkspfActivity");
  */
  return true;
}
/*
double prodPFCands::GetTrackActivity(const edm::Handle<pat::PackedCandidateCollection> & other_pfcands, const pat::PackedCandidate* track) 
{
  if (track->pt()<5.) return -1.0;
  double trkiso(0.); 
  double r_iso = 0.3;
  for (const pat::PackedCandidate &other_pfc : *other_pfcands) 
  {
    if (other_pfc.charge()==0) continue;
    double dr = deltaR(other_pfc, *track);
    if (dr < r_iso || dr > 0.4) continue; // activity annulus
    float dz_other = other_pfc.dz();
    if( fabs(dz_other) > 0.1 ) continue;
    trkiso += other_pfc.pt();
  }
  double activity = trkiso/track->pt();
  return activity;
}
*/
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(prodPFCands);
