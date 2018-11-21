#ifndef EcalTrackExtra_H
#define EcalTrackExtra_H
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "FastSimulation/CaloGeometryTools/interface/CaloSegment.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "CLHEP/GenericFunctions/IncompleteGamma.hh"
#include "FastSimulation/Event/interface/FSimTrack.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include <iostream>
#include <string>
#include <map>

class CaloGeometryHelper;
class FSimEvent;

typedef math::XYZVector XYZVector;
typedef math::XYZVector XYZPoint;

class EcalTrackExtra : public edm::stream::EDAnalyzer<>
{
 public:
  explicit EcalTrackExtra(const edm::ParameterSet&);
  ~EcalTrackExtra();
  virtual void beginRun(edm::Run const&, edm::EventSetup const& );
  virtual void analyze(const edm::Event & iEvent,const edm::EventSetup & c);
  virtual void endRun();

 private:
  void doExtrapolation(const math::XYZVectorD& pos,const math::XYZTLorentzVectorD & mom,int charge,std::vector<CaloSegment>& segments  );
  float getEnergy(const CaloSegment& seg, const EcalRecHitCollection & barrel, const EcalRecHitCollection &endcap);
  int doEnergyDeposits(const std::vector<CaloSegment>& segments, const EcalRecHitCollection & barrel, const EcalRecHitCollection &endcap,std::vector<float> & deposits,std::vector<float>&);
  bool  showerParam(const FSimTrack & mySimTrack, double &a, double &b);
  double deposit(double t, double a, double b, double dt) ;
  const reco::Candidate * mcMatch(const edm::View<reco::Candidate> *,const reco::GsfTrackRef,float &);

 private:
  edm::InputTag inputTagGSFTracks_;
  edm::InputTag inputTagBarrelRecHits_;
  edm::InputTag inputTagEndcapRecHits_;
  edm::InputTag inputTagPFCandidates_;
  edm::InputTag inputTagTruth_;

  edm::EDGetTokenT<reco::GsfTrackCollection> gsfTracksToken_;
  edm::EDGetTokenT<EBRecHitCollection> ebRecHitsToken_;
  edm::EDGetTokenT<EERecHitCollection> eeRecHitsToken_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatesToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > mcToken_;

  Genfun::IncompleteGamma myIncompleteGamma;
  bool matchWithMC_;
  double deltarMatch_ ;
  double maxPt_;

  CaloGeometryHelper * myGeometry;
  FSimEvent * mySimEvent;
  DQMStore * dbe;
  TH2F *h0,*h2,*h4;
  TH1F *h10,*h12,*h14;
  TH2F* h42,*h43;
  TH2F* h22,*h23;
  TH2F* h100,*h102,*h104;
  TH2F* h142,*h143;
  TH2F* h122,*h123;
  TH2F* h200, *h202, *h203;
  TH1F* h30;
  TH2F *h32;
  TH2F* h502, *h503;
  TH2F* h602, *h603;
  TH1F * h700,* h702,*h704;
};
#endif
