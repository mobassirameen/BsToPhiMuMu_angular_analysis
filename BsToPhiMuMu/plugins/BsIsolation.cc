#include "BsIsolation.h"
#include "TLorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"


#define TRKMAXR 110.0 // [cm]
#define TRKMAXZ 280.0 // [cm]


BsIsolation::BsIsolation(reco::Vertex bestVtx, 
                         edm::Handle<reco::TrackCollection> tracks, 
                         edm::ESHandle<MagneticField> bFieldHandle,
			 edm::ESHandle<TransientTrackBuilder> theTTBuilder_,
			 const AnalyticalImpactPointExtrapolator* impactPointExtrapolator_,
                         reco::BeamSpot beamSpot,
                         const pat::Muon* mum,
                         const pat::Muon* mup,
                         uint itrkm, // trk minus index
                         uint itrkp,  // trk plus index
                         const reco::TransientTrack refitMumTT, 
                         const reco::TransientTrack refitMupTT,
                         const reco::TransientTrack refitTrkmTT,
                         const reco::TransientTrack refitTrkpTT,
                         const RefCountedKinematicTree vertexFitTree			 
			 ){
  ClosestApproachInRPhi ClosestApp;
  float summpt_(0);
  float sumppt_(0);
  float sumkmpt_(0);
  float sumkppt_(0);

  vertexFitTree->movePointerToTheTop(); // Bs --> phi(KK) mu+ mu-                                                          
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex vertex = vertexFitTree->currentDecayVertex();

  double  bpx = b_KP->currentState().globalMomentum().x();
  double  bpy = b_KP->currentState().globalMomentum().y();
  double  bpz = b_KP->currentState().globalMomentum().z();
  double masstmp = b_KP->currentState().mass();
  TLorentzVector tmp_bmeson_lv;
  tmp_bmeson_lv.SetXYZM(bpx,bpy, bpz, masstmp);

  // for (std::vector<reco::Vertex>::const_iterator iVertex = vertices->begin(); iVertex != vertices->end(); iVertex++) { 
  //   bestVtx = *(iVertex); if (bestVtx.isValid() == true) break; 
  // }
  
  for (uint itrkiso =0 ;  itrkiso < tracks->size(); itrkiso++)
    {

      reco::TrackRef tkiso(tracks,itrkiso) ;                                                
      if ( itrkiso == itrkm  || itrkiso == itrkp)           continue;
      //       if (!tkiso->quality(reco::TrackBase::highPurity))     continue;
      if ( tkiso->pt() < 0.8 )                              continue;


      // check if the track is one of the two muons
      bool skip_this_track = false;              
      for (unsigned int i = 0; i < mum->numberOfSourceCandidatePtrs(); ++i) {
	const edm::Ptr<reco::Candidate> & source = mum->sourceCandidatePtr(i);
	if (! ( (mum->sourceCandidatePtr(i)).isNonnull() &&  (mum->sourceCandidatePtr(i)).isAvailable() ))   continue;
	const reco::Candidate & cand = *(source);
	if (cand.charge() == 0 || cand.bestTrack() == nullptr)      continue;
	try{ cand.bestTrack()->eta();}
	catch(...) { continue;}
	if ( deltaR(tkiso->eta(),tkiso->phi(),cand.bestTrack()->eta(),cand.bestTrack()->phi()) < 0.00001 ) {
	  skip_this_track = true;
	  break;
	}
      }
      if (skip_this_track) continue;
      for (unsigned int i = 0; i < mup->numberOfSourceCandidatePtrs(); ++i) {
	const edm::Ptr<reco::Candidate> & source = mup->sourceCandidatePtr(i);
	if (! ( (mup->sourceCandidatePtr(i)).isNonnull() &&  (mup->sourceCandidatePtr(i)).isAvailable() ))   continue;
	const reco::Candidate & cand = *(source);
	if (cand.charge() == 0 || cand.bestTrack() == nullptr)      continue;
	try{ cand.bestTrack()->eta();}
	catch(...) { continue;}
	if ( deltaR(tkiso->eta(),tkiso->phi(),cand.bestTrack()->eta(),cand.bestTrack()->phi()) < 0.00001 ) {
	  skip_this_track = true;
	  break;
	}
      }
      if (skip_this_track) continue;


      // requirement that the track is not associated to any PV
      // if (findPV(trk_index, recVtxColl) == 1 ) continue;
      
      const reco::TransientTrack TrackIsoTT((*tkiso), &(*bFieldHandle));
      /// Bs Meson Isolation 

      VertexDistance3D distance3D;      
      const reco::TransientTrack transTrack = theTTBuilder_->build(*tkiso);
      assert(impactPointExtrapolator_);
      //if (!vertex->vertexIsValid()) {continue;}
      const GlobalPoint BVP = GlobalPoint( vertex->position() );

      //std::cout<<"secondary vertex globalpoint "<< BVP.x()<<" \t y "<<BVP.y()<<" \t z "<<BVP.z()<<std::endl;
      //std::cout<<"Beamspot globalpoint "<< beamSpot.position().x()<<" \t y "<<beamSpot.position().y()<<" \t z "<<beamSpot.position().z()<<std::endl;
      //double rho = TrackIsoTT.initialFreeState().transverseCurvature();     

      if(transTrack.isValid() && vertex->vertexIsValid()){
	
	//std::cout<<" It is an valid transient track with rho "<<rho<<std::endl;
	//TrajectoryStateOnSurface  tsos = impactPointExtrapolator_->extrapolate(TrackIsoTT.initialFreeState(), BVP);
	TrajectoryStateOnSurface  tsos = impactPointExtrapolator_->extrapolate(transTrack.initialFreeState(), BVP);
	if (tsos.isValid()) {      	
	  
	  Measurement1D doca = distance3D.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), vertex->vertexState());
	  double svDoca = doca.value();
	  //std::cout<<"doca test "<<svDoca<<std::endl;
	  
	  docatrks.push_back(svDoca);
	// check to ensure the goodness of the track
	  if( svDoca<0.05){
	    bmass_isopts.push_back( tkiso->pt() );
	    bmass_isodr.push_back ( deltaR(tkiso ->eta(), tkiso->phi(), tmp_bmeson_lv.Eta(), tmp_bmeson_lv.Phi() ));
	  }
	}
      }
      if (! (TrackIsoTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),
                                                                  beamSpot.position().y(),
                                                                  beamSpot.position().z())).perigeeError().transverseImpactParameterError() >0) ) continue;
      
    
      // add new iso
      ClosestApp.calculate(TrackIsoTT.initialFreeState(), refitMumTT.initialFreeState());
      if (ClosestApp.status() != false)
	{
	  if ( ClosestApp.distance() < 0.1 ) {
	    mum_isopts.push_back( tkiso->pt() );	    
	    mum_isodr.push_back ( deltaR(tkiso ->eta(), tkiso->phi(), refitMumTT.track().eta(), refitMumTT.track().phi() ));   
	    summpt_ += tkiso->pt();
	  }                 
	}
      ClosestApp.calculate(TrackIsoTT.initialFreeState(), refitMupTT.initialFreeState());
      if (ClosestApp.status() != false)
	{
	  if ( ClosestApp.distance() < 0.1 ) {
	    mup_isopts.push_back( tkiso->pt() );
	    mup_isodr.push_back ( deltaR(tkiso ->eta(), tkiso->phi(), refitMupTT.track().eta(), refitMupTT.track().phi() ));   
	    sumppt_ += tkiso->pt();
	  }                 
	}
      ClosestApp.calculate(TrackIsoTT.initialFreeState(), refitTrkmTT.initialFreeState());
      if (ClosestApp.status() != false)
	{
	  if ( ClosestApp.distance() < 0.1 ) {
	    trkm_isopts.push_back( tkiso->pt() );
	    trkm_isodr.push_back ( deltaR(tkiso ->eta(), tkiso->phi(), refitTrkmTT.track().eta(), refitTrkmTT.track().phi() ));   
	    sumkppt_ += tkiso->pt();
	  }                 
	}
      ClosestApp.calculate(TrackIsoTT.initialFreeState(), refitTrkpTT.initialFreeState());
      if (ClosestApp.status() != false)
	{
	  if ( ClosestApp.distance() < 0.1 ) {
	    trkp_isopts.push_back( tkiso->pt() );
	    trkp_isodr.push_back ( deltaR(tkiso ->eta(), tkiso->phi(), refitTrkpTT.track().eta(), refitTrkpTT.track().phi() ));   	                   
	    sumkmpt_ += tkiso->pt();
	  }
	}  
    }           
  summpt = summpt_;
  sumppt = sumppt_;
  sumkmpt = sumkmpt_;
  sumkppt = sumkppt_;
             
  //std::cout<<"summpt "<<summpt<<" sumppt "<<sumppt<<" sumkppt "<<sumkppt<<" sumkmpt "<<sumkmpt<<std::endl;
}
