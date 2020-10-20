#ifndef BSISOLATION_H
#define BSISOLATION_H

#include <string>
#include <vector>

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"


class BsIsolation
{
 public:
 
  BsIsolation( reco::Vertex, 
	       edm::Handle<reco::TrackCollection>,
	       edm::ESHandle<MagneticField> , 
	       edm::ESHandle<TransientTrackBuilder> ,
	       const AnalyticalImpactPointExtrapolator*,
	       reco::BeamSpot ,
	       const pat::Muon*, 
	       const pat::Muon*, 	       	       
	       uint , // trk minus index
	       uint ,  // trk plus index
	       const reco::TransientTrack, 
	       const reco::TransientTrack,
	       const reco::TransientTrack,
	       const reco::TransientTrack,
	       const RefCountedKinematicTree 
	       );
  ~BsIsolation() {};
 
  std::vector<double> mum_isopts, mup_isopts, trkm_isopts, trkp_isopts; 
  std::vector<double> mum_isodr,  mup_isodr,  trkm_isodr,  trkp_isodr; 
  std::vector<double> bmass_isopts,bmass_isodr;
  std::vector<double> docatrks;
  float summpt;
  float sumppt;
  float sumkmpt;
  float sumkppt;


 private:

  const ParticleMass muonMass =   0.10565837;
  const ParticleMass pionMass =   0.13957018;
  const ParticleMass kaonMass =   0.493677;

  float mumasserr = 3.5e-9;

};
#endif
