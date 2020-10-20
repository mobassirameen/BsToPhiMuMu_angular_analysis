// -*- C++ -*-
//
// Package:    BsToPhiMuMu
// Class:      BsToPhiMuMu
// 
/**\class BsToPhiMuMu BsToPhiMuMu.cc BpHaNA/BsToPhiMuMu/src/BsToPhiMuMu.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//=====================================================================
// original author:  Niladribihari Sahoo,42 3-024,+41227662373,        |
//         copyright @ N.Sahoo, NISER, Bhubaneswar                     |
//         created:  Sat Nov 28 07:33:44 CET 2015                      |
//         added saveGenInfo  (sat 16 jan 2016)                        |
//         added soft muon id info, trig matching info (sat 16 jan'16) | 
//         added truth-matching vars (fri 28 oct'16)                   |
//         removed unnecessary vars (fri 28 oct'16)                    |
//         added quality variables (Tue 28 Feb'17)                     |
//                                                                     |
// Edit by: Jhovanny Mejia <jhovanny.andres.mejia.guisao@cern.ch>      | 
// Edit date <2016-08-11>                                              | 
// Editted by : Chandiprasad Kar                                       |
// Rewritten for MINIAOD                                               |
// Added pileup weight, Impact parameter variable w.r.t vtx            | 
//=====================================================================
// $Id$
//
//


//-----------------------                                                                                                                      
// system include files                                                                                                                             
//-----------------------                                                                                                                     
#include <memory>

//----------------------                                                                                                                       
// user include files                                                                                                                                
//----------------------                                                                                                                         
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"

#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include "myanalyzer/BsToPhiMuMu/plugins/BsIsolation.cc"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

using namespace std;
using namespace reco;
using namespace edm;


const int MUONMINUS_PDG_ID = 13;
const int KAONPLUS_PDG_ID = 321;
const int PHI_PDG_ID = 333;      // phi(1020)
const int BS_PDG_ID = 531;
const int JPSI_PDG_ID = 443;
const int PSI2S_PDG_ID = 100443;

const double PI = 3.141592653589793;



//-----------------------
// class declaration
//-----------------------

class BsToPhiMuMu : public edm::EDAnalyzer {
   public:
      explicit BsToPhiMuMu(const edm::ParameterSet&);
      ~BsToPhiMuMu();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);



  bool buildBsToPhiMuMu(const edm::Event &);

  void calLS (double, double, double, double, double, double, double,
              double, double,  double, double, double, double, double,
              double, double, double, double, double*, double*);

  void calCosAlpha (double, double, double, double, double,
                    double, double, double, double, double,
                    double, double, double, double,
                    double, double, double, double,
                    double*, double*);

  void calCosAlpha2d (double, double, double, double, double,
                      double, double, double, double, double,
                      double, double, double, double,
                      double, double, double, double,
                      double*, double*);

  void calCtau(RefCountedKinematicTree, double &, double &);
  double calEta(double, double, double);
  double calPhi(double, double, double);
  double calEtaPhiDistance (double, double, double, double, double, double);
  void clearVariables();

  bool hasBeamSpot(const edm::Event&);

  bool calClosestApproachTracks(const reco::TransientTrack,
                                const reco::TransientTrack,
                                double&, double &, double &);

  bool hasGoodPhiVertex(const reco::TransientTrack, const reco::TransientTrack, 
			reco::TransientTrack &, reco::TransientTrack &, 
			double &, double &);

  bool hasGoodMuonDcaBs (const reco::TransientTrack, double &, double &);
  bool hasGoodTrackDcaBs (const reco::TransientTrack, double &, double &);
  bool hasGoodTrackDcaPoint (const reco::TransientTrack, const GlobalPoint,
                             double, double &, double &);

  bool hasGoodBsMass(RefCountedKinematicTree, double &);  

  bool hasGoodBsVertex(const reco::TransientTrack, const reco::TransientTrack,
		       const reco::TransientTrack, const reco::TransientTrack,
		       double &, double &, double &,  RefCountedKinematicTree &);
  bool calClosestApproachBs (const RefCountedKinematicTree ,  double & , double & );
  bool hasGoodMuMuVertex (const reco::TransientTrack, const reco::TransientTrack,
                          reco::TransientTrack &, reco::TransientTrack &,
                          double &, double &, double &, double &, double &,
                          double &, double &, double &);

  //bool IsTheSame(const reco::TrackRef& , const pat::Muon& );
  bool hasGoodTrack(const edm::Event&, const reco::TrackRef , double &);

  bool hasPrimaryVertex(const edm::Event &);

  //void hltReport(const edm::Event&);

  bool matchMuonTrack (const edm::Event&, const reco::TrackRef);

  bool matchPrimaryVertexTracks ();
  std::pair<double,double> pionImpactParameter(reco::TransientTrack , reco::Vertex    );

  void saveBsToPhiMuMu(const RefCountedKinematicTree);
  void saveBsVertex(RefCountedKinematicTree);
  void saveBsCosAlpha(RefCountedKinematicTree);
  void saveBsCosAlpha2d(RefCountedKinematicTree);
  void saveBsLsig(RefCountedKinematicTree);
  void saveBsCtau(RefCountedKinematicTree);
  
  void saveIPvtx(const edm::Event&, const reco::TransientTrack, const reco::TransientTrack,
		 const reco::TransientTrack , const reco::TransientTrack );

  void saveGenInfo(const edm::Event&);
  void savePhiVariables(RefCountedKinematicTree,
			   reco::VertexCompositeCandidate);

  void saveDimuVariables(double, double, double, double, double, double,
                         double, double, double, double, double, double,
                         double, double);
  //void saveMuonTriggerMatches(const pat::Muon, const pat::Muon);
  
  void savePUinMC(const edm::Event& iEvent);
  bool skipOscillations (const reco::GenParticle &, edm::Handle<reco::GenParticleCollection>);
  bool isAncestor(const reco::Candidate*, const reco::Candidate * );
  
  void saveTruthMatch(const edm::Event& iEvent);

      // ----------member data ---------------------------


  // --- begin input from python file ---                                                                                                           
  string OutputFileName_;
  bool BuildBsToPhiMuMu_;

  //----------------------                                                                                                                                 
  // particle properties                                                                                                                                   
  //----------------------                                                                                                                                    
  ParticleMass MuonMass_;
  float MuonMassErr_;
  ParticleMass KaonMass_;
  float KaonMassErr_;
  ParticleMass PhiMass_;
  float PhiMassErr_;
  double BsMass_;

  //----------                                                                                                                                                
  // labels                                                                                                                                               
  //----------                                                                                                                                               
 
  edm::EDGetTokenT<reco::GenParticleCollection>         prunedGenToken_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection>    packedGenToken_;

  edm::EDGetTokenT<edm::TriggerResults> TriggerResultsLabel_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales>            triggerPrescales_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

  std::vector<std::string> trigTable_;
  //std::vector<std::string> l1Table_;
 
  //edm::EDGetTokenT<GlobalAlgBlkBxCollection> l1results_;
  //edm::EDGetTokenT<GlobalExtBlkBxCollection> l1ext_;

  edm::EDGetTokenT<reco::BeamSpot> BeamSpotLabel_;
  edm::EDGetTokenT<reco::VertexCollection> VertexLabel_;
  //edm::EDGetTokenT<pat::MuonCollection> MuonLabel_;
  edm::EDGetTokenT<edm::View<pat::Muon> > MuonLabel_;
  edm::EDGetTokenT<reco::TrackCollection> TrackLabel_; 
  //edm::EDGetTokenT<edm::View<pat::PackedCandidate>> TrackLabel_;

  edm::EDGetTokenT<std::vector< PileupSummaryInfo>> puToken_;  
  const AnalyticalImpactPointExtrapolator* impactPointExtrapolator_;
  vector<string> TriggerNames_;
  vector<string> LastFilterNames_;

  //---------------                                                                                                                                        
  // gen particle                                                                                                                                     
  //---------------                                                                                                                                          
  bool   IsMonteCarlo_;
  bool   KeepGENOnly_;
  double TruthMatchMuonMaxR_;
  double TruthMatchKaonMaxR_;
  //double TruthMatchPhiMaxVtx_;

  //---------------------                                                                                                                       
  // pre-selection cuts                                                                                                                                   
  //---------------------                                                                                                                              
  double MuonMinPt_;
  double MuonMaxEta_;
  double MuonMaxDcaBs_;
  double TrkMinPt_;
  double TrkMinDcaSigBs_;
  double TrkMaxR_;
  double TrkMaxZ_;
  double MuMuMaxDca_;
  double MuMuMinVtxCl_;
  double MuMuMinPt_;
  double MuMuMinInvMass_;
  double MuMuMaxInvMass_;
  double MuMuMinLxySigmaBs_;
  double MuMuMinCosAlphaBs_;
  double PhiMinMass_;
  double PhiMaxMass_;
  double BsMinVtxCl_;
  double BsMinMass_;
  double BsMaxMass_;


  //--------------------                                                                                                                       
  // Across the event                                                                                                                                   
  //--------------------                                                                                                                                     
  map<string, string> mapTriggerToLastFilter_;
  reco::BeamSpot beamSpot_;
  edm::ESHandle<TransientTrackBuilder> theTTBuilder_;
  edm::ESHandle<MagneticField> bFieldHandle_;
  reco::Vertex primaryVertex_;

  //-----------------                                                                                                                                          
  // Root Variables                                                                                                                                          
  //-----------------                                                                                                                                  
  TFile* fout_;
  TTree* tree_;

  unsigned int run, event, lumiblock, nprivtx;
  vector<string> *triggernames;
  vector<int> *triggerprescales;

  //----------                                                                                                                                              
  // dimuon                                                                                                                                                 
  //----------                                                                                                                                             
  vector<double> *mumdcabs, *mumdcabserr, *mumpx, *mumpy, *mumpz;
  vector<double> *mupdcabs, *mupdcabserr, *muppx, *muppy, *muppz;
  vector<double> *mumutrkr, *mumutrkz , *mumudca;
  vector<double> *mumuvtxcl, *mumulsbs, *mumulsbserr;
  vector<double> *mumucosalphabs, *mumucosalphabserr;
  vector<double> *mumumass, *mumumasserr;

  //-----------------------                                                                                                                                 
  // soft muon variables                                                                                                                                     
  //-----------------------                                                                                                                                  
  vector<bool>   *mumisgoodmuon, *mupisgoodmuon ;
  vector<int>    *mumcharge, *mupcharge ;
  vector<bool>   *mumtrackermuon, *muptrackermuon ;
  vector<bool>   *mumloosemuon, *muploosemuon ;
  vector<int>    *mumnpixhits, *mupnpixhits, *mumnpixlayers, *mupnpixlayers;
  vector<int>    *mumntrkhits, *mupntrkhits, *mumntrklayers, *mupntrklayers;
  vector<double> *mumnormchi2, *mupnormchi2;
  vector<double> *mumCL, *mupCL;
  vector<int>    *mumtrkqual, *muptrkqual;  /* added track quality vars */
  vector<double> *mumdxyvtx, *mupdxyvtx, *mumdzvtx, *mupdzvtx;
  vector<string> *mumtriglastfilter, *muptriglastfilter;
  vector<double> *mumpt, *muppt, *mumeta, *mupeta;

  std::vector<std::vector<double> > *mumIsoPt, *mumIsodR;
  std::vector<std::vector<double> > *mupIsoPt, *mupIsodR;

  std::vector<std::vector<double> > *kptrkIsoPt, *kptrkIsodR;
  std::vector<std::vector<double> > *kmtrkIsoPt, *kmtrkIsodR;
  vector<double> *rawmumpt, *rawmuppt, *rawmumphi, *rawmupphi,*rawmumeta, *rawmupeta;
  vector<double> *rawTrkmpt, *rawTrkppt, *rawTrkmphi, *rawTrkpphi,*rawTrkmeta, *rawTrkpeta;
  std::vector<std::vector<double> > *docatrk;
  std::vector<std::vector<double> > *bmassIsoPt, *bmassIsodR;

  vector<int>    *tri_LMNTk, *tri_PsipTk, *tri_JpsiTk; 
  vector<double> *tri_dr0, *tri_dr1;
  vector<double> *tri_dpt0, *tri_dpt1;

  vector<double> *kptrkMinIP2D, *kptrkMinIP2DE;
  vector<double> *kmtrkMinIP2D, *kmtrkMinIP2DE;
  vector<double> *mupMinIP2D, *mupMinIP2DE;
  vector<double> *mumMinIP2D, *mumMinIP2DE;
  vector<double> *kptrkMinIP, *kptrkMinIPE;
  vector<double> *kmtrkMinIP, *kmtrkMinIPE;
  vector<double> *mupMinIP, *mupMinIPE;
  vector<double> *mumMinIP, *mumMinIPE;
  //--------------
  // kaon track   
  //--------------                                                                                                                                  
  //vector<int> *trkchg; // +1 for K+, -1 for K-                                                                                                 
  //vector<double> *trkpx, *trkpy, *trkpz, *trkpt;
  vector<double> *kptrkdcabs, *kptrkdcabserr;
  vector<double> *kmtrkdcabs, *kmtrkdcabserr;
  
  //--------------
  // phi(1020)
  //--------------
  vector<double> *kmtrknormchi2, *kptrknormchi2;
  vector<double> *kmtrkCL, *kptrkCL;
  vector<int>    *kmtrkqual, *kptrkqual;  /* added track quality vars */
 
  vector<int>    *kpchg;
  vector<double> *kppx, *kppy, *kppz ;
  vector<int>    *kmchg;
  vector<double> *kmpx, *kmpy, *kmpz ;
  //  vector<double> *phipx, *phipy, *phipz;
  //  vector<double>  *phivtxx, *phivtxy, *phivtxz, *phivtxcl, *philsbs, *philsbserr;

  vector<double> *phimass /*, *phimasserr, *phibarmass, *phibarmasserr */ ;

  //-----------------
  // Bs and Bsbar
  //-----------------
  int nb;
  vector<double> *bpx, *bpxerr, *bpy, *bpyerr, *bpz, *bpzerr ;
  vector<double> *bmass, *bmasserr;
  vector<double> *bvtxcl, *bvtxx, *bvtxxerr, *bvtxy, *bvtxyerr, *bvtxz, *bvtxzerr;
  vector<double> *bcosalphabs, *bcosalphabserr, *bcosalphabs2d, *bcosalphabs2derr, *blsbs, *blsbserr, *bctau, *bctauerr; 
  vector<double> *bdcabs, *bdcabserr;
  // vector<double> *bbarmass, *bbarmasserr;

  //----------
  // For MC   
  //----------                                                                                                                                               
  double genbpx, genbpy, genbpz, genbpid;
  double genphipx, genphipy, genphipz;
  double genphivtxx, genphivtxy, genphivtxz;

  //                                                                                                                                                        
  int genkpchg;
  double genkppx, genkppy, genkppz;
  int genkmchg;
  double genkmpx, genkmpy, genkmpz;

  double genmumpx, genmumpy, genmumpz;
  double genmuppx, genmuppy, genmuppz;
  
  string decname;

  std::vector<double>       *bunchXingMC, *numInteractionsMC, *trueNumInteractionsMC;

  vector<bool> *istruemum, *istruemup, *istruekp, *istruekm, /**istruephi,*/ *istruebs;

  double fpuw8;
  // ############################
  // # Pileup information in MC #
  // ############################
  
  // Comment:
  // - PileupSummaryInfo::getTrueNumInteractions() gives the distribution of the mean number of interactions per crossing.
  // Since this is the mean value of the poisson distribution from which the number of interactions in- and out-of-time are
  // generated, no additional information should be required for reweighting if these values are matched in data and Monte Carlo.
  // - PileupSummaryInfo::getPU_NumInteractions() gives the expected mean number of interactions per crossing for each LumiSection.
  // Therefore the pileup histogram will contain the distribution of the number of interactions one would actually observe given
  // a poisson of that mean. So, this distribution is what one would see if one counted the number of events seen in a given beam
  // crossing (by looking at the number of vertices in data, for example. This would be appropriate for pileup reweighting based
  // on in-time-only distributions.

  //-----------------------
  // variables to monitor  
  //-----------------------                                                                                                                         
  TDatime t_begin_ , t_now_ ;
  int n_processed_, n_selected_;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BsToPhiMuMu::BsToPhiMuMu(const edm::ParameterSet& iConfig):

  OutputFileName_(iConfig.getParameter<string>("OutputFileName")),
  BuildBsToPhiMuMu_(iConfig.getUntrackedParameter<bool>("BuildBsToPhiMuMu")),

  //************************
  // particle properties    
  //************************                                                                                                                           

  MuonMass_(iConfig.getUntrackedParameter<double>("MuonMass")),
  MuonMassErr_(iConfig.getUntrackedParameter<double>("MuonMassErr")),
  KaonMass_(iConfig.getUntrackedParameter<double>("KaonMass")),
  KaonMassErr_(iConfig.getUntrackedParameter<double>("KaonMassErr")),
  BsMass_(iConfig.getUntrackedParameter<double>("BsMass")),

  //***********
  // labels    
  //***********                                                                                                                                         
  prunedGenToken_(consumes<reco::GenParticleCollection >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packed"))),

  TriggerResultsLabel_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResultsLabel"))),  
  triggerPrescales_ (consumes<pat::PackedTriggerPrescales>            (iConfig.getParameter<edm::InputTag>("prescales"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  
  trigTable_( iConfig.getParameter<std::vector<std::string> >("TriggerNames")),   

  BeamSpotLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotLabel"))),
  VertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VertexLabel"))),
// MuonLabel_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("MuonLabel"))),
  MuonLabel_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("MuonLabel"))),
  TrackLabel_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("TrackLabel"))),
  //TrackLabel_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("TrackLabel"))),
  
  //pileup
  puToken_(consumes<std::vector< PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PuInfoTag"))), 
  impactPointExtrapolator_(0),
 //TriggerNames_(iConfig.getParameter< vector<string> >("TriggerNames")),
  LastFilterNames_(iConfig.getParameter< vector<string> >("LastFilterNames")),
  
  //***************
  // gen particle  
  //***************                                                                                                                                  
  IsMonteCarlo_(iConfig.getUntrackedParameter<bool>("IsMonteCarlo")),
  KeepGENOnly_(iConfig.getUntrackedParameter<bool>("KeepGENOnly")),
  
  TruthMatchMuonMaxR_(iConfig.getUntrackedParameter<double>("TruthMatchMuonMaxR")),
  TruthMatchKaonMaxR_(iConfig.getUntrackedParameter<double>("TruthMatchKaonMaxR")),
 
  //*********************
  // pre-selection cuts                                                                                                                                
  //*********************
  MuonMinPt_(iConfig.getUntrackedParameter<double>("MuonMinPt")),
  MuonMaxEta_(iConfig.getUntrackedParameter<double>("MuonMaxEta")),
  MuonMaxDcaBs_(iConfig.getUntrackedParameter<double>("MuonMaxDcaBs")),

  TrkMinPt_(iConfig.getUntrackedParameter<double>("TrkMinPt")),
  TrkMinDcaSigBs_(iConfig.getUntrackedParameter<double>("TrkMinDcaSigBs")),
  TrkMaxR_(iConfig.getUntrackedParameter<double>("TrkMaxR")),
  TrkMaxZ_(iConfig.getUntrackedParameter<double>("TrkMaxZ")),

  MuMuMaxDca_(iConfig.getUntrackedParameter<double>("MuMuMaxDca")),
  MuMuMinVtxCl_(iConfig.getUntrackedParameter<double>("MuMuMinVtxCl")),
  MuMuMinPt_(iConfig.getUntrackedParameter<double>("MuMuMinPt")),
  MuMuMinInvMass_(iConfig.getUntrackedParameter<double>("MuMuMinInvMass")),
  MuMuMaxInvMass_(iConfig.getUntrackedParameter<double>("MuMuMaxInvMass")),
  MuMuMinLxySigmaBs_(iConfig.getUntrackedParameter<double>("MuMuMinLxySigmaBs")),
  MuMuMinCosAlphaBs_(iConfig.getUntrackedParameter<double>("MuMuMinCosAlphaBs")),

  PhiMinMass_(iConfig.getUntrackedParameter<double>("PhiMinMass")),
  PhiMaxMass_(iConfig.getUntrackedParameter<double>("PhiMaxMass")),
  BsMinVtxCl_(iConfig.getUntrackedParameter<double>("BsMinVtxCl")),
  BsMinMass_(iConfig.getUntrackedParameter<double>("BsMinMass")),
  BsMaxMass_(iConfig.getUntrackedParameter<double>("BsMaxMass")),
  
  tree_(0),
  triggernames(0), triggerprescales(0),
  mumdcabs(0), mumdcabserr(0), mumpx(0), mumpy(0), mumpz(0),
  mupdcabs(0),  mupdcabserr(0), muppx(0),  muppy(0), muppz(0),
  mumutrkr(0), mumutrkz(0), mumudca(0),  mumuvtxcl(0),  mumulsbs(0),
  mumulsbserr(0), mumucosalphabs(0),  mumucosalphabserr(0),
  mumumass(0), mumumasserr(0),
  mumisgoodmuon(0), mupisgoodmuon(0),
  mumcharge(0), mupcharge(0), 
  mumtrackermuon(0), muptrackermuon(0), 
  mumloosemuon(0), muploosemuon(0), 

  mumnpixhits(0), mupnpixhits(0), mumnpixlayers(0), mupnpixlayers(0),
  mumntrkhits(0), mupntrkhits(0), mumntrklayers(0), mupntrklayers(0),
  mumnormchi2(0), mupnormchi2(0), 
  mumCL(0), mupCL(0), 
  mumtrkqual(0), muptrkqual(0),         /* added */
  mumdxyvtx(0), mupdxyvtx(0),
  mumdzvtx(0), mupdzvtx(0), mumtriglastfilter(0), muptriglastfilter(0),
  mumpt(0), muppt(0), mumeta(0), mupeta(0),

//trkchg(0), trkpx(0), trkpy(0), trkpz(0), trkpt(0),
  mumIsoPt(0), mumIsodR(0),
  mupIsoPt(0), mupIsodR(0),

  kptrkIsoPt(0), kptrkIsodR(0),
  kmtrkIsoPt(0), kmtrkIsodR(0),

  rawmumpt(0), rawmuppt(0), rawmumphi(0), rawmupphi(0),rawmumeta(0), rawmupeta(0),
  rawTrkmpt(0), rawTrkppt(0), rawTrkmphi(0), rawTrkpphi(0),rawTrkmeta(0), rawTrkpeta(0),
  docatrk(0),
  bmassIsoPt(0), bmassIsodR(0),
  tri_LMNTk(0), tri_PsipTk(0), tri_JpsiTk(0),
  tri_dr0(0), tri_dr1(0),
  tri_dpt0(0), tri_dpt1(0),

  kptrkMinIP2D(0), kptrkMinIP2DE(0),
  kmtrkMinIP2D(0), kmtrkMinIP2DE(0),
  mupMinIP2D(0), mupMinIP2DE(0),
  mumMinIP2D(0), mumMinIP2DE(0),

  kptrkMinIP(0), kptrkMinIPE(0),
  kmtrkMinIP(0), kmtrkMinIPE(0),
  mupMinIP(0), mupMinIPE(0),
  mumMinIP(0), mumMinIPE(0),

  kptrkdcabs(0), kptrkdcabserr(0),
  kmtrkdcabs(0), kmtrkdcabserr(0),
  kmtrknormchi2(0), kptrknormchi2(0),
  kmtrkCL(0), kptrkCL(0),
  kmtrkqual(0), kptrkqual(0),         /* added */
 
  kpchg(0),
  kppx(0), kppy(0), kppz(0),
  
  kmchg(0),
  kmpx(0), kmpy(0), kmpz(0),


  //phipx(0), phipy(0), phipz(0),
  //phivtxx(0), phivtxy(0), phivtxz(0),

  phimass(0), /* phimasserr(0), phibarmass(0), phibarmasserr(0), */

  nb(0), bpx(0), bpxerr(0), bpy(0), bpyerr(0), bpz(0), bpzerr(0), bmass(0), bmasserr(0),
  bvtxcl(0), bvtxx(0), bvtxxerr(0), bvtxy(0), bvtxyerr(0), bvtxz(0), bvtxzerr(0),
  bcosalphabs(0), bcosalphabserr(0), bcosalphabs2d(0), bcosalphabs2derr(0), blsbs(0), blsbserr(0), bctau(0), bctauerr(0),

  bdcabs(0), bdcabserr(0),
  //bbarmass(0), bbarmasserr(0),

  genbpx(0), genbpy(0), genbpz(0), genbpid(0),
  genphipx(0), genphipy(0), genphipz(0), genphivtxx(0), genphivtxy(0), genphivtxz(0),

  genkpchg(0),
  genkppx(0), genkppy(0), genkppz(0),

  genkmchg(0),
  genkmpx(0), genkmpy(0), genkmpz(0),

  genmumpx(0), genmumpy(0), genmumpz(0),
  genmuppx(0), genmuppy(0), genmuppz(0),

  decname(""),

  bunchXingMC(0), numInteractionsMC(0), trueNumInteractionsMC(0),
  istruemum(0), istruemup(0), istruekp(0), istruekm(0), /*istruephi(0),*/ istruebs(0),
  fpuw8(0)



{
   //now do what ever initialization is needed
  //assert(TriggerNames_.size() == LastFilterNames_.size());
  //for (size_t i = 0; i < TriggerNames_.size(); ++i)
  //  mapTriggerToLastFilter_[TriggerNames_[i]] = LastFilterNames_[i];
  
}


BsToPhiMuMu::~BsToPhiMuMu()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
BsToPhiMuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  clearVariables();

  run = iEvent.id().run() ;
  event = iEvent.id().event() ;
  lumiblock = iEvent.luminosityBlock();

  n_processed_ += 1;


  if (IsMonteCarlo_) saveGenInfo(iEvent);

  ////hltReport(iEvent);

  ///if ( KeepGENOnly_){
  if ( IsMonteCarlo_ && KeepGENOnly_){
    tree_->Fill();
    n_selected_ += 1;
  }else{
    //hltReport(iEvent);
    if ( hasBeamSpot(iEvent) ) {
      iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle_);      
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder_);      
      if ( bFieldHandle_.isValid() && hasPrimaryVertex(iEvent) ) {
        buildBsToPhiMuMu(iEvent) ;
        if (IsMonteCarlo_){ 
	  savePUinMC(iEvent);
	  saveTruthMatch(iEvent);
	  
	}
        n_selected_ += 1;
      }
    }

    if (IsMonteCarlo_ || nb > 0){ // Keep failed events for MC to calculate reconstruction efficiency.                                                               
      tree_->Fill();
    }
  }


  clearVariables();


}


// ------------ method called once each job just before starting event loop  ------------
void 
BsToPhiMuMu::beginJob()
{

  t_begin_.Set();
  printf("\n ---------- Begin Job ---------- \n");
  t_begin_.Print();

  n_processed_ = 0;
  n_selected_ = 0;

  fout_ = new TFile(OutputFileName_.c_str(), "RECREATE");
  fout_->cd();
  ///edm::Service<TFileService> fs;



  tree_ = new TTree ("tree", "BsToPhiMuMu");
  ///tree_ = fs->make<TTree>("tree","Bs->J/psi kaskey menos ntuple");

  tree_->Branch("run", &run, "run/i");
  tree_->Branch("event", &event, "event/i");
  tree_->Branch("lumiblock", &lumiblock, "lumiblock/i");
  tree_->Branch("nprivtx", &nprivtx, "nprivtx/i");
  tree_->Branch("triggernames", &triggernames);
  tree_->Branch("triggerprescales", &triggerprescales);
  tree_->Branch("mumdcabs", &mumdcabs);
  tree_->Branch("mumdcabserr", &mumdcabserr);
  tree_->Branch("mumpx", &mumpx);
  tree_->Branch("mumpy", &mumpy);
  tree_->Branch("mumpz", &mumpz);
  tree_->Branch("mupdcabs", &mupdcabs);
  tree_->Branch("mupdcabserr", &mupdcabserr);
  tree_->Branch("muppx", &muppx);
  tree_->Branch("muppy", &muppy);
  tree_->Branch("muppz", &muppz);
  tree_->Branch("mumutrkr", &mumutrkr);
  tree_->Branch("mumutrkz", &mumutrkz);
  tree_->Branch("mumudca", &mumudca);
  tree_->Branch("mumuvtxcl", &mumuvtxcl);
  tree_->Branch("mumulsbs", &mumulsbs);
  tree_->Branch("mumulsbserr", &mumulsbserr);
  tree_->Branch("mumucosalphabs", &mumucosalphabs);
  tree_->Branch("mumucosalphabserr", &mumucosalphabserr);
  tree_->Branch("mumumass", &mumumass);
  tree_->Branch("mumumasserr", &mumumasserr);
  tree_->Branch("mumisgoodmuon", &mumisgoodmuon);
  tree_->Branch("mupisgoodmuon", &mupisgoodmuon);

  tree_->Branch("mumtrackermuon", &mumtrackermuon);
  tree_->Branch("muptrackermuon", &muptrackermuon);
  tree_->Branch("mumloosemuon", &mumloosemuon);
  tree_->Branch("muploosemuon", &muploosemuon);
  tree_->Branch("mumcharge", &mumcharge);
  tree_->Branch("mupcharge", &mupcharge);

  tree_->Branch("mumnpixhits", &mumnpixhits);
  tree_->Branch("mupnpixhits", &mupnpixhits);
  tree_->Branch("mumnpixlayers", &mumnpixlayers);
  tree_->Branch("mupnpixlayers", &mupnpixlayers);
  tree_->Branch("mumntrkhits", &mumntrkhits);
  tree_->Branch("mupntrkhits", &mupntrkhits);
  tree_->Branch("mumntrklayers", &mumntrklayers);
  tree_->Branch("mupntrklayers", &mupntrklayers);
  tree_->Branch("mumnormchi2", &mumnormchi2);
  tree_->Branch("mupnormchi2", &mupnormchi2);
  tree_->Branch("mumCL", &mumCL);
  tree_->Branch("mupCL", &mupCL);
  tree_->Branch("mumtrkqual", &mumtrkqual);  /* added quality vars */
  tree_->Branch("muptrkqual", &muptrkqual);
  tree_->Branch("mumdxyvtx", &mumdxyvtx);
  tree_->Branch("mupdxyvtx", &mupdxyvtx);
  tree_->Branch("mumdzvtx", &mumdzvtx);
  tree_->Branch("mupdzvtx", &mupdzvtx);
  tree_->Branch("mumtriglastfilter", &mumtriglastfilter);
  tree_->Branch("muptriglastfilter", &muptriglastfilter);
  tree_->Branch("mumpt", &mumpt);
  tree_->Branch("muppt", &muppt);
  tree_->Branch("mumeta", &mumeta);
  tree_->Branch("mupeta", &mupeta);

  tree_->Branch("mupIsoPt", &mupIsoPt);
  tree_->Branch("mupIsodR", &mupIsodR);
  tree_->Branch("mumIsoPt", &mumIsoPt);
  tree_->Branch("mumIsodR", &mumIsodR);
  tree_->Branch("kptrkIsoPt", &kptrkIsoPt);
  tree_->Branch("kptrkIsodR", &kptrkIsodR);
  tree_->Branch("kmtrkIsoPt", &kmtrkIsoPt);
  tree_->Branch("kmtrkIsodR", &kmtrkIsodR);
  tree_->Branch("docatrk", &docatrk);
  tree_->Branch("bmassIsoPt", &bmassIsoPt);
  tree_->Branch("bmassIsodR", &bmassIsodR);

  tree_->Branch("rawmumpt",&rawmumpt);
  tree_->Branch("rawmuppt",&rawmuppt);
  tree_->Branch("rawmumeta",&rawmumeta);
  tree_->Branch("rawmupeta",&rawmupeta);
  tree_->Branch("rawmumphi",&rawmumphi);
  tree_->Branch("rawmupphi",&rawmupphi);
  
  tree_->Branch("rawTrkmpt",&rawTrkmpt);
  tree_->Branch("rawTrkppt",&rawTrkppt);
  tree_->Branch("rawTrkmeta",&rawTrkmeta);
  tree_->Branch("rawTrkpeta",&rawTrkpeta);
  tree_->Branch("rawTrkmphi",&rawTrkmphi);
  tree_->Branch("rawTrkpphi",&rawTrkpphi);

  
  tree_->Branch("tri_LMNTk", &tri_LMNTk);
  tree_->Branch("tri_PsipTk", &tri_PsipTk);
  tree_->Branch("tri_JpsiTk", &tri_JpsiTk);
  tree_->Branch("tri_dr0", &tri_dr0);
  tree_->Branch("tri_dr1", &tri_dr1);
  tree_->Branch("tri_dpt0", &tri_dpt0);
  tree_->Branch("tri_dpt1", &tri_dpt1);


  tree_->Branch("kptrkMinIP2D", &kptrkMinIP2D);
  tree_->Branch("kptrkMinIP2DE", &kptrkMinIP2DE);
  tree_->Branch("kmtrkMinIP2D", &kmtrkMinIP2D);
  tree_->Branch("kmtrkMinIP2DE", &kmtrkMinIP2DE);
  tree_->Branch("mupMinIP2D", &mupMinIP2D);
  tree_->Branch("mupMinIP2DE", &mupMinIP2DE);
  tree_->Branch("mumMinIP2D", &mumMinIP2D);
  tree_->Branch("mumMinIP2DE", &mumMinIP2DE);
  tree_->Branch("kptrkMinIP", &kptrkMinIP);
  tree_->Branch("kptrkMinIPE", &kptrkMinIPE);
  tree_->Branch("kmtrkMinIP", &kmtrkMinIP);
  tree_->Branch("kmtrkMinIPE", &kmtrkMinIPE);
  tree_->Branch("mupMinIP", &mupMinIP);
  tree_->Branch("mupMinIPE", &mupMinIPE);
  tree_->Branch("mumMinIP", &mumMinIP);
  tree_->Branch("mumMinIPE", &mumMinIPE);


  //tree_->Branch("trkchg", &trkchg);
  //tree_->Branch("trkpx", &trkpx);
  //tree_->Branch("trkpy", &trkpy);
  //tree_->Branch("trkpz", &trkpz);
  //tree_->Branch("trkpt", &trkpt);
  tree_->Branch("kptrkdcabs", &kptrkdcabs);
  tree_->Branch("kptrkdcabserr", &kptrkdcabserr);
  tree_->Branch("kmtrkdcabs", &kmtrkdcabs);
  tree_->Branch("kmtrkdcabserr", &kmtrkdcabserr);
  tree_->Branch("kpchg", &kpchg);
  tree_->Branch("kppx", &kppx);
  tree_->Branch("kppy", &kppy);
  tree_->Branch("kppz", &kppz);
  tree_->Branch("kmchg", &kmchg);
  tree_->Branch("kmpx", &kmpx);
  tree_->Branch("kmpy", &kmpy);
  tree_->Branch("kmpz", &kmpz);
  tree_->Branch("kmtrknormchi2", &kmtrknormchi2);
  tree_->Branch("kptrknormchi2", &kptrknormchi2);
  tree_->Branch("kmtrkCL", &kmtrkCL);
  tree_->Branch("kptrkCL", &kptrkCL);
  tree_->Branch("kmtrkqual", &kmtrkqual);  /* added quality vars */
  tree_->Branch("kptrkqual", &kptrkqual);
  //tree_->Branch("phipx", &phipx);
  //tree_->Branch("phipy", &phipy);
  //tree_->Branch("phipz", &phipz);
  //tree_->Branch("phivtxx", &phivtxx);
  //tree_->Branch("phivtxy", &phivtxy);
  //tree_->Branch("phivtxz", &phivtxz);
  tree_->Branch("phimass", &phimass);
  //tree_->Branch("phimasserr", &phimasserr);
  //tree_->Branch("phibarmass", &phibarmass);
  //tree_->Branch("phibarmasserr", &phibarmasserr);
  tree_->Branch("nb", &nb, "nb/I");
  tree_->Branch("bpx", &bpx);
  tree_->Branch("bpxerr", &bpxerr);
  tree_->Branch("bpy", &bpy);
  tree_->Branch("bpyerr", &bpyerr);
  tree_->Branch("bpz", &bpz);
  tree_->Branch("bpzerr", &bpzerr);
  tree_->Branch("bmass", &bmass);
  tree_->Branch("bmasserr", &bmasserr);
  tree_->Branch("bvtxcl", &bvtxcl);
  tree_->Branch("bvtxx", &bvtxx);
  tree_->Branch("bvtxxerr", &bvtxxerr);
  tree_->Branch("bvtxy", &bvtxy);
  tree_->Branch("bvtxyerr", &bvtxyerr);
  tree_->Branch("bvtxz", &bvtxz);
  tree_->Branch("bvtxzerr", &bvtxzerr);
  tree_->Branch("bcosalphabs", &bcosalphabs);
  tree_->Branch("bcosalphabserr", &bcosalphabserr);
  tree_->Branch("bcosalphabs2d", &bcosalphabs2d);
  tree_->Branch("bcosalphabs2derr", &bcosalphabs2derr);
  tree_->Branch("blsbs", &blsbs);
  tree_->Branch("blsbserr", &blsbserr);
  tree_->Branch("bctau", &bctau);
  tree_->Branch("bctauerr", &bctauerr);
  //tree_->Branch("bbarmass", &bbarmass);
  //tree_->Branch("bbarmasserr", &bbarmasserr);
  tree_->Branch("bdcabs", &bdcabs);
  tree_->Branch("bdcabserr", &bdcabserr);

  if (IsMonteCarlo_) {
    tree_->Branch("genbpx",      &genbpx     , "genbpx/D"    );
    tree_->Branch("genbpy",      &genbpy     , "genbpy/D"    );
    tree_->Branch("genbpz",      &genbpz     , "genbpz/D"    );
    tree_->Branch("genbpid",     &genbpid    , "genbpid/D"    );
    tree_->Branch("genphipx",    &genphipx   , "genphipx/D"  );
    tree_->Branch("genphipy",    &genphipy   , "genphipy/D"  );
    tree_->Branch("genphipz",    &genphipz   , "genphipz/D"  );
    tree_->Branch("genphivtxx",  &genphivtxx    , "genphivtxx/D"   );
    tree_->Branch("genphivtxy",  &genphivtxy    , "genphivtxy/D"   );
    tree_->Branch("genphivtxz",  &genphivtxz    , "genphivtxz/D"   );
    tree_->Branch("genkpchg",   &genkpchg   , "genkpchg/I"   );
    tree_->Branch("genkppx",    &genkppx    , "genkppx/D"   );
    tree_->Branch("genkppy",    &genkppy    , "genkppy/D"   );
    tree_->Branch("genkppz",    &genkppz    , "genkppz/D"   );
    tree_->Branch("genkmchg",   &genkmchg   , "genkmchg/I"   );
    tree_->Branch("genkmpx",    &genkmpx    , "genkmpx/D"   );
    tree_->Branch("genkmpy",    &genkmpy    , "genkmpy/D"   );
    tree_->Branch("genkmpz",    &genkmpz    , "genkmpz/D"   );
    tree_->Branch("genmumpx",   &genmumpx   , "genmumpx/D"  );
    tree_->Branch("genmumpy",   &genmumpy   , "genmumpy/D"  );
    tree_->Branch("genmumpz",   &genmumpz   , "genmumpz/D"  );
    tree_->Branch("genmuppx",   &genmuppx   , "genmuppx/D"  );
    tree_->Branch("genmuppy",   &genmuppy   , "genmuppy/D"  );
    tree_->Branch("genmuppz",   &genmuppz   , "genmuppz/D"  );

    tree_->Branch("decname",    &decname);
    tree_->Branch("istruemum",  &istruemum );
    tree_->Branch("istruemup",  &istruemup );
    tree_->Branch("istruekp",   &istruekp  );
    tree_->Branch("istruekm",   &istruekm  );
    //tree_->Branch("istruephi",  &istruephi );
    tree_->Branch("istruebs",   &istruebs  );

    tree_->Branch("bunchXingMC",&bunchXingMC);
    tree_->Branch("numInteractionsMC",&numInteractionsMC);
    tree_->Branch("trueNumInteractionsMC",&trueNumInteractionsMC);
    tree_->Branch("fpuw8",   &fpuw8  );

  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
BsToPhiMuMu::endJob() 
{


  fout_->cd();
  //tree_->Write();
  /////tree_->GetDirectory()->cd();
  tree_->Write();
 
  fout_->Close();

  t_now_.Set();
  printf(" \n ---------- End Job ---------- \n" ) ;
  t_now_.Print();
  printf(" processed: %i \n selected: %i \n \
 duration: %i sec \n rate: %g evts/sec\n",
	 n_processed_, n_selected_,
	 t_now_.Convert() - t_begin_.Convert(),
	 float(n_processed_)/(t_now_.Convert()-t_begin_.Convert()) );

}

// ------------ method called when starting to processes a run  ------------
void 
BsToPhiMuMu::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
BsToPhiMuMu::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
BsToPhiMuMu::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
BsToPhiMuMu::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BsToPhiMuMu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void 
BsToPhiMuMu::clearVariables(){

  run = 0;
  event = 0;
  lumiblock = 0;
  nprivtx = 0;
  triggernames->clear();
  triggerprescales->clear();
  mumdcabs->clear();  mumdcabserr->clear();  mumpx->clear();   mumpy->clear();  mumpz->clear();
  mupdcabs->clear();  mupdcabserr->clear();  muppx->clear();   muppy->clear();  muppz->clear();
  mumutrkr->clear(); mumutrkz->clear();
  mumudca->clear();  mumuvtxcl->clear();   mumulsbs->clear();  mumulsbserr->clear();
  mumucosalphabs->clear();  mumucosalphabserr->clear();
  mumumass->clear(); mumumasserr->clear();
  mumisgoodmuon->clear();  mupisgoodmuon->clear();
  mumcharge->clear(); mupcharge->clear() ;
  mumtrackermuon->clear(); muptrackermuon->clear() ;
  mumloosemuon->clear(); muploosemuon->clear() ;
  mumnpixhits->clear();  mupnpixhits->clear();  mumnpixlayers->clear();  mupnpixlayers->clear();
  mumntrkhits->clear();  mupntrkhits->clear();  mumntrklayers->clear();  mupntrklayers->clear();

  mumnormchi2->clear(); mupnormchi2->clear();
  mumCL->clear(); mupCL->clear();
  mumtrkqual->clear(); muptrkqual->clear();  
  mumdxyvtx->clear(); mupdxyvtx->clear();
  mumdzvtx->clear(); mupdzvtx->clear();
  mumtriglastfilter->clear(); muptriglastfilter->clear();
  mumpt->clear(); muppt->clear();
  mumeta->clear(); mupeta->clear();

  mumIsoPt->clear(); mumIsodR->clear();
  mupIsoPt->clear(); mupIsodR->clear();
  kptrkIsoPt->clear(); kptrkIsodR->clear();
  kmtrkIsoPt->clear(); kmtrkIsodR->clear();
  rawmumpt->clear(); rawmuppt->clear(); rawmumphi->clear(); rawmupphi->clear();rawmumeta->clear(); rawmupeta->clear();
  rawTrkmpt->clear(); rawTrkppt->clear(); rawTrkmphi->clear(); rawTrkpphi->clear();rawTrkmeta->clear(); rawTrkpeta->clear();
  docatrk->clear();
  bmassIsoPt->clear(); bmassIsodR->clear();
  tri_LMNTk->clear(); tri_PsipTk->clear(); tri_JpsiTk->clear(); 
  tri_dr0->clear(); tri_dr1->clear();
  tri_dpt0->clear(); tri_dpt1->clear();

  kptrkMinIP2D->clear(); kptrkMinIP2DE->clear();
  kmtrkMinIP2D->clear(); kmtrkMinIP2DE->clear();
  mupMinIP2D->clear(); mupMinIP2DE->clear();
  mumMinIP2D->clear(); mumMinIP2DE->clear();
  kptrkMinIP->clear(); kptrkMinIPE->clear();
  kmtrkMinIP->clear(); kmtrkMinIPE->clear();
  mupMinIP->clear(); mupMinIPE->clear();
  mumMinIP->clear(); mumMinIPE->clear();

  //trkchg->clear(); trkpx->clear(); trkpy->clear(); trkpz->clear(); trkpt->clear();
  kptrkdcabs->clear(); kptrkdcabserr->clear();
  kmtrkdcabs->clear(); kmtrkdcabserr->clear();

  kmtrknormchi2->clear(); kptrknormchi2->clear();
  kmtrkCL->clear(); kptrkCL->clear();
  kmtrkqual->clear(); kptrkqual->clear();    
  kpchg->clear();
  kppx->clear(); kppy->clear(); kppz->clear();
 
  kmchg->clear();
  kmpx->clear(); kmpy->clear(); kmpz->clear();

  //phipx->clear(); phipy->clear(); phipz->clear();
  //phivtxx->clear(); phivtxy->clear(); phivtxz->clear();

  phimass->clear(); 
  //phimasserr->clear();
  //phibarmass->clear(); phibarmasserr->clear();

  nb = 0;

  bpx->clear(); bpxerr->clear(); bpy->clear();  bpyerr->clear();
  bpz->clear(); bpzerr->clear();

  bmass->clear(); bmasserr->clear();
  bvtxcl->clear(); bvtxx->clear(); bvtxxerr->clear(); bvtxy->clear(); bvtxyerr->clear();
  bvtxz->clear(); bvtxzerr->clear(); bcosalphabs->clear(); bcosalphabserr->clear();
  bcosalphabs2d->clear(); bcosalphabs2derr->clear();
  blsbs->clear(); blsbserr->clear(); bctau->clear(); bctauerr->clear();
  bdcabs->clear(); bdcabserr->clear();
  //bbarmass->clear(); bbarmasserr->clear();

  if (IsMonteCarlo_) {

    genbpx = 0;  genbpy = 0;  genbpz = 0; genbpid = 0;
    genphipx = 0;  genphipy = 0;  genphipz = 0;
    genphivtxx = 0; genphivtxy = 0; genphivtxz = 0;

    genkpchg = 0;
    genkppx = 0;  genkppy = 0;  genkppz = 0;
    genkmchg = 0;
    genkmpx = 0;  genkmpy = 0;  genkmpz = 0;

    genmumpx = 0;  genmumpy = 0;  genmumpz = 0;
    genmuppx = 0;  genmuppy = 0;  genmuppz = 0;

    decname = "";
    bunchXingMC->clear();
    numInteractionsMC->clear();
    trueNumInteractionsMC->clear();
    istruemum->clear(); istruemup->clear(); istruekp->clear();
    istruekm->clear(); /*istruephi->clear();*/ istruebs->clear();

    fpuw8 = 0;
  }

}
bool
BsToPhiMuMu::hasBeamSpot(const edm::Event& iEvent)
{
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  //iEvent.getByLabel(BeamSpotLabel_, beamSpotHandle);
  iEvent.getByToken(BeamSpotLabel_, beamSpotHandle);

  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("myBeam") << "No beam spot available from EventSetup" ;
    return false;
  }

  beamSpot_ = *beamSpotHandle;
  return true;
}

bool
BsToPhiMuMu::hasPrimaryVertex(const edm::Event& iEvent)
{
  edm::Handle<reco::VertexCollection> recVtxs;
  //iEvent.getByLabel(VertexLabel_, recVtxs);
  iEvent.getByToken(VertexLabel_, recVtxs);
  nprivtx = recVtxs->size();

  for (std::vector<reco::Vertex>::const_iterator iVertex = recVtxs->begin();
       iVertex != recVtxs->end(); iVertex++) {
    primaryVertex_ = *(iVertex);
    if (primaryVertex_.isValid()) break;
  }

  if (!primaryVertex_.isValid()) return false;

  return true;
}

//-------------------------------------------------------
//  main function to build BsToPhiMuMu candidate
//-------------------------------------------------------

bool
BsToPhiMuMu::buildBsToPhiMuMu(const edm::Event& iEvent)
{

  // init variables
  // edm::Handle<pat::MuonCollection> patMuonHandle;
  //iEvent.getByToken(MuonLabel_, patMuonHandle);
  edm::Handle< edm::View<pat::Muon> > patMuonHandle;
  iEvent.getByToken(MuonLabel_,patMuonHandle);

  if( patMuonHandle->size() < 2 ) return false;
  //std::cout<<" Pat Muon size "<<patMuonHandle->size()<<std::endl;

  // edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  // iEvent.getByToken(TrackLabel_, thePATTrackHandle);
   
  AnalyticalImpactPointExtrapolator extrapolator(bFieldHandle_.product());      
  impactPointExtrapolator_ = &extrapolator;

  edm::Handle<reco::TrackCollection> thePATTrackHandle;
  iEvent.getByToken(TrackLabel_, thePATTrackHandle);
    
  bool passed;
  double DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr;
  double mumutrk_R, mumutrk_Z, DCAmumu;
  double trk_R, trk_Z, trk_DCA;
  reco::TransientTrack refitMupTT, refitMumTT;
  double mu_mu_vtx_cl, mu_mu_pt, mu_mu_mass, mu_mu_mass_err;
  double MuMuLSBS, MuMuLSBSErr;
  double MuMuCosAlphaBS, MuMuCosAlphaBSErr;
  double kaon_trk_pt;
  reco::TransientTrack refitKmTT, refitKpTT;
  double phi_mass /*phi_mass_err*/; 
  double phi_vtx_cl;
  double b_vtx_chisq, b_vtx_cl, b_mass;
  //double bbar_mass;
  double DCABsBS, DCABsBSErr;
  double DCAPhiTrkpBS, DCAPhiTrkpBSErr;
  double DCAPhiTrkmBS, DCAPhiTrkmBSErr;
  RefCountedKinematicTree vertexFitTree, barVertexFitTree;

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales>            triggerPrescales;
  edm::Handle<edm::TriggerResults> hltTriggerResults;

  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  iEvent.getByToken( TriggerResultsLabel_, hltTriggerResults ); 
 
  std::string cctk_0 = "hltJpsiTkAllConeTracksIter";
  std::string cctk_1 = "hltPsiPrimeTkAllConeTracksIter";
  std::string cctk_2 = "hltLowMassNonResonantTkAllConeTracksIter";
  std::stringstream myString;
  std::vector<float> obj_eta, obj_phi, obj_pt;
  std::vector<int> obj_charge;

  //bool foundOneTrig = false;
  const edm::TriggerNames &names = iEvent.triggerNames(*hltTriggerResults);
  for (unsigned int i = 0, n = hltTriggerResults->size(); i < n; ++i) {
    for (unsigned int it = 0; it < trigTable_.size(); it++){
      if (names.triggerName(i).find(trigTable_[it]) != std::string::npos && hltTriggerResults->accept(i))
	{	 
	  triggernames->push_back(names.triggerName(i));
	  triggerprescales->push_back(triggerPrescales->getPrescaleForIndex(i));
	  // cout<<"Trigger name "<<names.triggerName(i)<<" Prescale "<<triggerPrescales->getPrescaleForIndex(i)<<endl;
	  //foundOneTrig = true;
	}
    }
  }
  //if ( iEvent.isRealData() && !foundOneTrig) return;
     
  if (triggerObjects.isValid()) {
    //std::cout << "will try to match trigger object with track " << triggerObjects->size() << endl;    
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
        obj.unpackPathNames(names);
        obj.unpackFilterLabels(iEvent,*hltTriggerResults);
    	std::string cc1 = obj.collection();
	
    	for (unsigned int i = 0; i < trigTable_.size(); i++)
        {
    	  myString.clear(); myString.str(""); myString << trigTable_[i] << "*";
    	  if ( obj.hasPathName(myString.str().c_str(), true, true) )
            {

    	      std::vector<std::string> pathNamesAll  = obj.pathNames(false);
    	      std::vector<std::string> pathNamesLast = obj.pathNames(true);
	      
    	      //std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
	      
    	      if ( (cc1.find(cctk_0) != std::string::npos ) || (cc1.find(cctk_1) != std::string::npos) || ( cc1.find(cctk_2) != std::string::npos)) {
		
		obj_eta.push_back(obj.eta());
		obj_phi.push_back(obj.phi());
		obj_pt.push_back(obj.pt());
		obj_charge.push_back(obj.charge());
		
		// std::cout << myString.str().c_str() << std::endl;
		// std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
		// std::cout << "\t   Charge:   "<<obj.charge()<<std::endl;
		// std::cout << "\t   Collection: " << obj.collection() << std::endl;
		// std::cout << "\t   Type IDs:   ";
		// for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
		// std::cout << "\t   Filters:    ";
		// for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
		// std::cout << std::endl;
	      }
    	    }
    	}
    }
  }
    
  // ---------------------------------
  // loop 1: mu-
  // ---------------------------------
  //for(View<pat::Muon>::const_iterator iMuonM = patMuonHandle->begin(); iMuonM != patMuonHandle->end(); ++iMuonM) {
  for(edm::View<pat::Muon>::const_iterator iMuon1 = patMuonHandle->begin(); iMuon1 != patMuonHandle->end(); ++iMuon1) 
    {      
      for(edm::View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != patMuonHandle->end(); ++iMuon2) 
	{
	    
	  if(iMuon1==iMuon2) continue;
	    
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue; // <-------------------------------------
	    
	  TrackRef muTrackp;  
	  TrackRef muTrackm;  
	    
	  if(iMuon1->charge() == 1){muTrackp = iMuon1->track();}
	  if(iMuon1->charge() == -1){muTrackm = iMuon1->track();}
	    
	  if(iMuon2->charge() == 1) {muTrackp = iMuon2->track();}
	  if(iMuon2->charge() == -1){muTrackm = iMuon2->track();}
	    
	  if( muTrackp.isNull() || muTrackm.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  // for (const pat::Muon &iMuonM : *patMuonHandle) 
	  //   {  
	  //     reco::TrackRef muTrackm = iMuonM.innerTrack();
	  //     //reco::TrackRef muTrackm = iMuonM->track();
	  //     if ( muTrackm.isNull() ) continue;
	  
	  
	  if((muTrackm->pt() < MuonMinPt_) ||
	     (fabs(muTrackm->eta()) > MuonMaxEta_)) continue;
	  // check mu- DCA to beam spot
	  const reco::TransientTrack muTrackmTT(muTrackm, &(*bFieldHandle_));
	  passed = hasGoodMuonDcaBs(muTrackmTT, DCAmumBS, DCAmumBSErr) ;
	  if ( ! passed ) continue;

	  if((muTrackp->pt() < MuonMinPt_) ||
	     (fabs(muTrackp->eta()) > MuonMaxEta_)) continue;	
	  
	  // check mu+ DCA to beam spot
	  const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle_));
	  passed = hasGoodMuonDcaBs(muTrackpTT, DCAmupBS, DCAmupBSErr);
	  if ( ! passed ) continue;
	  
	  // check goodness of muons closest approach and the 3D-DCA	
	  if ( !calClosestApproachTracks(muTrackpTT, muTrackmTT,
					 mumutrk_R, mumutrk_Z, DCAmumu)) continue;
	  
	  if ( mumutrk_R > TrkMaxR_ ||
	       mumutrk_Z > TrkMaxZ_ ||
	       DCAmumu > MuMuMaxDca_ ) continue;
	  
	  // check dimuon vertex
	  passed = hasGoodMuMuVertex(muTrackpTT, muTrackmTT, refitMupTT, refitMumTT,
				     mu_mu_vtx_cl, mu_mu_pt,
				     mu_mu_mass, mu_mu_mass_err,
				     MuMuLSBS, MuMuLSBSErr,
				     MuMuCosAlphaBS, MuMuCosAlphaBSErr);
	  	  
	  if ( !passed) continue;
	  
	  // ---------------------------------
	  // loop 3: track-
	  // ---------------------------------
	  //for(View<pat::PackedCandidate>::const_iterator iTrackM = thePATTrackHandle->begin();
	  //	  iTrackM != thePATTrackHandle->end(); ++iTrackM )
	  //{
	  //for ( vector<pat::GenericParticle>::const_iterator iTrackM
	  //   = thePATTrackHandle->begin();
	  //iTrackM != thePATTrackHandle->end(); ++iTrackM ) {
	  
	  for (uint iTrack1 =0 ;  iTrack1 < thePATTrackHandle->size(); iTrack1++){
	    for (uint iTrack2 =iTrack1+1 ;  iTrack2 < thePATTrackHandle->size(); iTrack2++){
	    
	      
	      reco::TrackRef Track1(thePATTrackHandle,iTrack1) ;
	      reco::TrackRef Track2(thePATTrackHandle,iTrack2) ;
	      
	      //reco::TrackRef Trackm = iTrackM->track();
	      if( (Track1->charge() * Track2->charge()) ==1)continue;
	      if ( Track1.isNull() || Track2.isNull() ) continue;
	      //if (! Trackm->quality(reco::Track::highPurity)) continue;			 
	      TrackRef Trackm;
	      TrackRef Trackp;
	      if(Track1->charge() == 1){Trackp = Track1; }
	      if(Track1->charge() == -1){Trackm = Track1; }
	        
	      if(Track2->charge() == 1) {Trackp = Track2;}
	      if(Track2->charge() == -1){Trackm = Track2;}
	        

	      if ( Trackm->pt() < TrkMinPt_)  continue;
	      if (fabs(Trackm->eta()) > 2.5)  continue;
	      passed = hasGoodTrack(iEvent, Trackm, kaon_trk_pt);
	      if(!passed) continue;

	      // compute track DCA to beam spot
	      const reco::TransientTrack theTrackmTT((*Trackm), &(*bFieldHandle_));
	      if(!theTrackmTT.isValid()) continue;
	      passed = hasGoodTrackDcaBs(theTrackmTT, DCAPhiTrkpBS, DCAPhiTrkpBSErr);
	      
	      if (!passed) continue;
	      
	      // ---------------------------------
	      // loop 4: track+
	      // ---------------------------------
	      // for ( vector<pat::GenericParticle>::const_iterator iTrackP
	      // 	= thePATTrackHandle->begin();
	      //       iTrackP != thePATTrackHandle->end(); ++iTrackP ) {
	      
	      //reco::TrackRef Trackp = iTrackP.track();
	      //if (! Trackp->quality(reco::Track::highPurity)) continue;
	      if ( Trackp->pt() < TrkMinPt_)  continue;
	      if (fabs(Trackp->eta()) > 2.5)  continue;
	      passed = hasGoodTrack(iEvent, Trackp, kaon_trk_pt);
	      if (!passed) continue;
	      
	      // compute track DCA to beam spot
	      const reco::TransientTrack theTrackpTT((*Trackp), &(*bFieldHandle_));
	      if(!theTrackpTT.isValid()) continue;
	      passed = hasGoodTrackDcaBs(theTrackpTT, DCAPhiTrkmBS, DCAPhiTrkmBSErr);
	      if (!passed) continue;
	      
	      
	      // check goodness of two tracks closest approach and the 3D-DCA
	      if (! calClosestApproachTracks(theTrackpTT, theTrackmTT,
					     trk_R, trk_Z, trk_DCA)) continue ;
	      if ( trk_R > TrkMaxR_ || trk_Z > TrkMaxZ_ ) continue;
	      
	      
	      TLorentzVector pion14V,pion24V,mu14V, mu24V, phi4V,Bs4V; 
	      pion14V.SetXYZM(theTrackpTT.track().px(),theTrackpTT.track().py(),theTrackpTT.track().pz(),KaonMass_);
	      pion24V.SetXYZM(theTrackmTT.track().px(),theTrackmTT.track().py(),theTrackmTT.track().pz(),KaonMass_);
	      phi4V=pion14V+pion24V;
	      if(phi4V.M()<0.950 || phi4V.M()>1.09) continue;
	      
	      // check two tracks vertex for Bs
	      if ( ! hasGoodPhiVertex(theTrackmTT, theTrackpTT, refitKmTT, refitKpTT, phi_vtx_cl, phi_mass) ) continue;
	      if ( phi_mass < PhiMinMass_ || phi_mass > PhiMaxMass_ ) continue;
	      
	      mu14V.SetXYZM(iMuon1->track()->px(),iMuon1->track()->py(),iMuon1->track()->pz(),MuonMass_);
	      mu24V.SetXYZM(iMuon2->track()->px(),iMuon2->track()->py(),iMuon2->track()->pz(),MuonMass_);
	      Bs4V=pion14V+pion24V+mu14V+mu24V;
	      if(Bs4V.M()<4.0 || Bs4V.M()>6.7) continue;
	      //cout<<"It clculated the Bs CL "<<b_vtx_cl<<endl;
	      // fit Bs vertex  mu- mu+ K- K+
	      if ( ! hasGoodBsVertex(muTrackmTT, muTrackpTT, theTrackmTT, theTrackpTT,
				     b_vtx_chisq, b_vtx_cl, b_mass,
				     vertexFitTree) ) continue;
	      //cout<<"It passed good vertex fit "<<endl;
	      if ( (b_vtx_cl < BsMinVtxCl_) || (b_mass < BsMinMass_) || (b_mass > BsMaxMass_) ) continue;
	      
	      //cout<<"It clculated the Bs mass "<<Bs4V.M()<<" CL "<<b_vtx_cl<<endl;
	      if( !calClosestApproachBs(vertexFitTree, DCABsBS, DCABsBSErr)) continue;
	      
	      
	      // need to check with primaryVertex tracks?
	    
	      nb++;
	      
	      //cout<<"number of B: "<<nb<<" mass of the B is : "<<Bs4V.M()<<endl;	  
	      // save the tree variables
	      saveDimuVariables(DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr,
				mumutrk_R, mumutrk_Z, DCAmumu, mu_mu_vtx_cl,
				MuMuLSBS, MuMuLSBSErr,
				MuMuCosAlphaBS, MuMuCosAlphaBSErr,
				mu_mu_mass, mu_mu_mass_err);
	      
	      
	      bdcabs->push_back(DCABsBS);
	      bdcabserr->push_back(DCABsBSErr);
	      
	      mumisgoodmuon->push_back(muon::isGoodMuon(*iMuon1, muon::TMOneStationTight));
	      mupisgoodmuon->push_back(muon::isGoodMuon(*iMuon2, muon::TMOneStationTight));
	      
	      mumcharge->push_back(iMuon1->charge());
	      mupcharge->push_back(iMuon2->charge());
	      mumnpixhits->push_back(muTrackm->hitPattern().numberOfValidPixelHits());
	      mupnpixhits->push_back(muTrackp->hitPattern().numberOfValidPixelHits());
	      mumnpixlayers->push_back(muTrackm->hitPattern().pixelLayersWithMeasurement());
	      mupnpixlayers->push_back(muTrackp->hitPattern().pixelLayersWithMeasurement());
	      
	      mumntrkhits->push_back(muTrackm->hitPattern().numberOfValidTrackerHits());
	      mupntrkhits->push_back(muTrackp->hitPattern().numberOfValidTrackerHits());
	      mumntrklayers->push_back(muTrackm->hitPattern().trackerLayersWithMeasurement());
	      mupntrklayers->push_back(muTrackp->hitPattern().trackerLayersWithMeasurement());
	      
	      mumnormchi2->push_back(muTrackm->normalizedChi2());
	      mupnormchi2->push_back(muTrackp->normalizedChi2());
	      
	      mumtrkqual->push_back(muTrackm->quality(reco::TrackBase::highPurity));   
	      muptrkqual->push_back(muTrackp->quality(reco::TrackBase::highPurity));
	      
	      mumdxyvtx->push_back(muTrackm->dxy(primaryVertex_.position()));
	      mupdxyvtx->push_back(muTrackp->dxy(primaryVertex_.position()));
	      
	      mumdzvtx->push_back(muTrackm->dz(primaryVertex_.position()));
	      mupdzvtx->push_back(muTrackp->dz(primaryVertex_.position()));
	      
	      mumpt->push_back(muTrackm->pt());
	      muppt->push_back(muTrackp->pt());
	    
	      mumeta->push_back(muTrackm->eta());
	      mupeta->push_back(muTrackp->eta());
	      
	    
	      kptrkdcabs->push_back(DCAPhiTrkpBS);
	      kptrkdcabserr->push_back(DCAPhiTrkpBSErr);
	      kmtrkdcabs->push_back(DCAPhiTrkmBS);
	      kmtrkdcabserr->push_back(DCAPhiTrkmBSErr);
	      phimass->push_back(phi_mass);
	      //phimasserr->push_back(phi_mass_err);
	      bvtxcl->push_back(b_vtx_cl);
	      //bmass->push_back(b_mass);
	      kmtrkqual->push_back( (int)Trackm->quality(reco::Track::highPurity));
	      kmtrkCL->push_back(TMath::Prob(theTrackmTT.chi2(), static_cast<int>(rint(theTrackmTT.ndof()))));
	      kmtrknormchi2->push_back(Trackm->normalizedChi2());
	      kptrkqual->push_back( (int)Trackp->quality(reco::Track::highPurity));
	      kptrkCL->push_back(TMath::Prob(theTrackpTT.chi2(), static_cast<int>(rint(theTrackpTT.ndof()))));
	      kptrknormchi2->push_back(Trackp->normalizedChi2());
	    
	      mupCL->push_back(TMath::Prob(muTrackpTT.chi2(), static_cast<int>(rint(muTrackpTT.ndof()))));
	      mumCL->push_back(TMath::Prob(muTrackmTT.chi2(), static_cast<int>(rint(muTrackmTT.ndof()))));
	    
	      float dr0 = 99999.;
	      float dpt0 = 99999.;
	      for (uint ii=0 ; ii < obj_eta.size(); ii++) {
		float dp = Trackp->phi() - obj_phi[ii];
		float de = Trackp->eta() - obj_eta[ii];	      
		if (dp>float(M_PI)) dp-=float(2*M_PI);  
		float dr = std::sqrt(de*de + dp*dp);
		if (dr < dr0) dr0 = dr;
		float dpt = Trackp->pt() - obj_pt[ii];
		if (abs(dpt) < dpt0) dpt0 = abs(dpt);
		//std::cout << "\tTrigger object: eta " << obj_eta[ii] << " phi " << obj_phi[ii] << " pt "<< obj_pt[ii] <<" -> " << dr0 << " dr "<<dr<< std::endl;
		//std::cout << "\tTrack object: eta " << Trackp->eta() << " phi " << Trackp->phi() <<" pt "<< Trackp->pt() << " -> " << dr0 << " dpt "<<dpt<< " dpt0 "<<dpt0<< std::endl;
	      }
	      //std::cout<<std::endl;
	      float dr1 = 99999.;
	      float dpt1 = 99999.;
	      for (uint ii=0 ; ii < obj_eta.size(); ii++) {
		float dp = Trackm->phi() - obj_phi[ii];
		float de = Trackm->eta() - obj_eta[ii];	    
		if (dp>float(M_PI)) dp-=float(2*M_PI);  
		float dr = std::sqrt(de*de + dp*dp);
		if (dr < dr1) dr1 = dr;
		float dpt = Trackm->pt() - obj_pt[ii];
		if (abs(dpt) < dpt1) dpt1 =abs(dpt);
		
		//std::cout << "\nTrigger object: eta " << obj_eta[ii] << " phi " << obj_phi[ii] << " pt "<< obj_pt[ii] <<" -> " << dr1 << " dr "<<dr<< std::endl;
		//std::cout << "\t m Track object: eta " << Trackm->eta() << " phi " << Trackm->phi() <<" pt "<< Trackm->pt() << " -> " << dr1 <<  " dpt "<<dpt<< " dpt0 "<<dpt1<<std::endl;
	      }
	      // std::cout <<" Dr0 "<< 1+int(dr0*1000) <<"\t dr1 "<< 1+int(dr1*1000)<<std::endl;
	      // std::cout <<" Dr0 "<< dr0 <<"\t dr1 "<< dr1<<std::endl;
	      // std::cout <<" dpt0 "<< dpt0 <<"\t dpt1 "<< dpt1<<std::endl;
	      
	      const pat::Muon* muon1 = &(*iMuon1);
	      const pat::Muon* muon2 = &(*iMuon2);
	      
	      // // Works for 2017	  
	      // const pat::TriggerObjectStandAloneCollection muHLTMatches1_t1 = muon1->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4PsiPrime");
	      // const pat::TriggerObjectStandAloneCollection muHLTMatches2_t1 = muon2->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4PsiPrime");
	      // const pat::TriggerObjectStandAloneCollection muHLTMatches1_t2 = muon1->triggerObjectMatchesByFilter("hltJpsiTkVertexFilter");
	      // const pat::TriggerObjectStandAloneCollection muHLTMatches2_t2 = muon2->triggerObjectMatchesByFilter("hltJpsiTkVertexFilter");
	      
	      // //const pat::TriggerObjectStandAloneCollection muHLTMatches1_t2 = muon1->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi");
	      // //const pat::TriggerObjectStandAloneCollection muHLTMatches2_t2 = muon2->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi");
	      
	      // const pat::TriggerObjectStandAloneCollection muHLTMatches1_t3 = muon1->triggerObjectMatchesByFilter("hltLowMassNonResonantTkVertexFilter");
	      // const pat::TriggerObjectStandAloneCollection muHLTMatches2_t3 = muon2->triggerObjectMatchesByFilter("hltLowMassNonResonantTkVertexFilter");
	      // std::cout<<"Muonm t1 "<<muHLTMatches1_t1.size()<<" Muonp "<<  muHLTMatches2_t1.size()<<std::endl;
	      // std::cout<<"Muonm t2 "<<muHLTMatches1_t2.size()<<" Muonp "<<  muHLTMatches2_t2.size()<<std::endl;
	      // std::cout<<"Muonm t3 "<<muHLTMatches1_t3.size()<<" Muonp "<<  muHLTMatches2_t3.size()<<std::endl;
	      
	      int tri_PsipTk_tmp = 0, tri_JpsiTk_tmp = 0, tri_LMNTk_tmp = 0;
	      
	      
	      //works for 2018, 2016
	      if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v*")!=nullptr) tri_PsipTk_tmp = 1;
	      if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr) tri_JpsiTk_tmp = 1;
	      if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v*")!=nullptr) tri_LMNTk_tmp = 1;
	      
	      // std::cout<<"Muonpsiprime "<<tri_PsipTk_tmp <<std::endl;
	      // std::cout<<"Muonm jpsitrk "<<tri_JpsiTk_tmp<<std::endl;
	      // std::cout<<"Muonm lowmasstrk "<<tri_LMNTk_tmp <<std::endl;
	      
	      
	      //   if (muHLTMatches1_t1.size() > 0 && muHLTMatches2_t1.size() > 0) tri_PsipTk_tmp = 1;
	      //   if (muHLTMatches1_t2.size() > 0 && muHLTMatches2_t2.size() > 0) tri_JpsiTk_tmp = 1;
	      //   if (muHLTMatches1_t3.size() > 0 && muHLTMatches2_t3.size() > 0) tri_LMNTk_tmp = 1;
	      // //if (muHLTMatches1_t4.size() > 0 && muHLTMatches2_t4.size() > 0) cond "<<std::endl;}//+int(1000*dr0);
	      
	      
	      tri_LMNTk->push_back( tri_LMNTk_tmp );	       
	      tri_PsipTk->push_back( tri_PsipTk_tmp );	       
	      tri_JpsiTk->push_back( tri_JpsiTk_tmp );
	      
	      tri_dr0->push_back(dr0);
	      tri_dr1->push_back(dr1);
	      tri_dpt0->push_back(abs(dpt0));
	      tri_dpt1->push_back(abs(dpt1));
	      mumloosemuon->push_back(muon1->isLooseMuon());
	      muploosemuon->push_back(muon2->isLooseMuon());
	      mumtrackermuon->push_back(muon1->isTrackerMuon());
	      muptrackermuon->push_back(muon2->isTrackerMuon());
	      
	      saveBsToPhiMuMu(vertexFitTree);
	      saveBsVertex(vertexFitTree);
	      saveBsCosAlpha(vertexFitTree);
	      saveBsCosAlpha2d(vertexFitTree);
	      saveBsLsig(vertexFitTree);
	      saveBsCtau(vertexFitTree);
	      saveIPvtx(iEvent, muTrackmTT, muTrackpTT, theTrackmTT, theTrackpTT);
	      
	      BsIsolation BsIso = BsIsolation(primaryVertex_,//bestVtx,
					      thePATTrackHandle, 
					      bFieldHandle_,					  
					      theTTBuilder_,
					      impactPointExtrapolator_,
					      beamSpot_ , 
					      muon1, muon2,
					      iTrack1, iTrack2,
					      muTrackmTT, muTrackpTT,
					      theTrackmTT, theTrackpTT,
					      vertexFitTree
					      );
	      mumIsoPt       -> push_back(BsIso.mum_isopts);
	      mupIsoPt       -> push_back(BsIso.mup_isopts);
	      kmtrkIsoPt     -> push_back(BsIso.trkm_isopts);
	      kptrkIsoPt     -> push_back(BsIso.trkp_isopts);
	      mumIsodR       -> push_back(BsIso.mum_isodr);
	      mupIsodR       -> push_back(BsIso.mup_isodr);
	      kmtrkIsodR     -> push_back(BsIso.trkm_isodr);
	      kptrkIsodR     -> push_back(BsIso.trkp_isodr);
	      rawmumpt        -> push_back(muTrackm -> pt() );
	      rawmumphi       -> push_back(muTrackm -> phi());
	      rawmumeta       -> push_back(muTrackm -> eta());
	      rawmuppt        -> push_back(muTrackp -> pt() );
	      rawmupphi       -> push_back(muTrackp -> phi());
	      rawmupeta       -> push_back(muTrackp -> eta());
	      rawTrkmpt    -> push_back(Trackm -> pt() );
	      rawTrkmphi   -> push_back(Trackm -> phi());
	      rawTrkmeta   -> push_back(Trackm -> eta());
	      rawTrkppt    -> push_back(Trackp -> pt() );
	      rawTrkpphi   -> push_back(Trackp -> phi());
	      rawTrkpeta   -> push_back(Trackp -> eta());
	      
	      bmassIsoPt    -> push_back(BsIso.bmass_isopts);
	      bmassIsodR    -> push_back(BsIso.bmass_isodr);
	      docatrk       -> push_back(BsIso.docatrks);
	      //std::cout<<"track charge "<<Track1->charge() <<" tack2 "<< Track2->charge()<<std::endl;
	    } // close track+ loop
	  } // close track- loop*/
	} // close mu+ loop
    } // close mu- loop

  
  if ( nb > 0) {
    //edm::LogInfo("myBs") << "Found " << nb << " Bs -> phi(KK) mu+ mu- ";
    std::cout<<"-------Number of Bs ------------ "<<nb<<std::endl;
    return true;
  }
  
  return false;
  
  
}

void
BsToPhiMuMu::calLS (double Vx, double Vy, double Vz,
		     double Wx, double Wy, double Wz,
		     double VxErr2, double VyErr2, double VzErr2,
		     double VxyCov, double VxzCov, double VyzCov,
		     double WxErr2, double WyErr2, double WzErr2,
		     double WxyCov, double WxzCov, double WyzCov,
		     double* deltaD, double* deltaDErr)
{
  *deltaD = sqrt((Vx-Wx) * (Vx-Wx) + (Vy-Wy) * (Vy-Wy) + (Vz-Wz) * (Vz-Wz));
  if (*deltaD > 0.)
    *deltaDErr = sqrt((Vx-Wx) * (Vx-Wx) * VxErr2 +
		      (Vy-Wy) * (Vy-Wy) * VyErr2 +
		      (Vz-Wz) * (Vz-Wz) * VzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*VxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*VxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*VyzCov +
		      
		      (Vx-Wx) * (Vx-Wx) * WxErr2 +
		      (Vy-Wy) * (Vy-Wy) * WyErr2 +
		      (Vz-Wz) * (Vz-Wz) * WzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*WxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*WxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*WyzCov) / *deltaD;
  else *deltaDErr = 0.;
}


void
BsToPhiMuMu::calCosAlpha (double Vx, double Vy, double Vz,
			   double Wx, double Wy, double Wz,
			   double VxErr2, double VyErr2, double VzErr2,
			   double VxyCov, double VxzCov, double VyzCov,
			   double WxErr2, double WyErr2, double WzErr2,
			   double WxyCov, double WxzCov, double WyzCov,
			   double* cosAlpha, double* cosAlphaErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;

  if ((Vnorm > 0.) && (Wnorm > 0.)) {
    *cosAlpha = VdotW / (Vnorm * Wnorm);
    *cosAlphaErr = sqrt( (
			  (Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			  (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			  (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +

			  (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			  (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			  (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			 (Wnorm*Wnorm*Wnorm*Wnorm) +
			 
			 ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			  (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			  (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +
			  
			  (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			  (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			  (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			 (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
  }  else {
    *cosAlpha = 0.;
    *cosAlphaErr = 0.;
  }

}


void
BsToPhiMuMu::calCosAlpha2d (double Vx, double Vy, double Vz,
			     double Wx, double Wy, double Wz,
			     double VxErr2, double VyErr2, double VzErr2,
			     double VxyCov, double VxzCov, double VyzCov,
			     double WxErr2, double WyErr2, double WzErr2,
			     double WxyCov, double WxzCov, double WyzCov,
			     double* cosAlpha2d, double* cosAlpha2dErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;

  if ((Vnorm > 0.) && (Wnorm > 0.)) {
    *cosAlpha2d = VdotW / (Vnorm * Wnorm);
    *cosAlpha2dErr = sqrt( (
			    (Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			    (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			    (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +

			    (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			    (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			    (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			   (Wnorm*Wnorm*Wnorm*Wnorm) +

			   ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			    (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			    (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +

			    (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			    (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			    (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			   (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
  }  else {
    *cosAlpha2d = 0.;
    *cosAlpha2dErr = 0.;
  }

}


bool
BsToPhiMuMu::hasGoodMuonDcaBs (const reco::TransientTrack muTrackTT,
				double &muDcaBs, double &muDcaBsErr)
{
  TrajectoryStateClosestToPoint theDCAXBS =
    muTrackTT.trajectoryStateClosestToPoint(
					    GlobalPoint(beamSpot_.position().x(),
							beamSpot_.position().y(),beamSpot_.position().z()));

  if ( !theDCAXBS.isValid() )  return false;

  muDcaBs = theDCAXBS.perigeeParameters().transverseImpactParameter();
  muDcaBsErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  if ( fabs(muDcaBs) > MuonMaxDcaBs_ )   return false;
  return true;
}

bool
BsToPhiMuMu::hasGoodTrackDcaBs (const reco::TransientTrack TrackTT,
				 double &DcaBs, double &DcaBsErr)
{
  TrajectoryStateClosestToPoint theDCAXBS =
    TrackTT.trajectoryStateClosestToPoint(
					  GlobalPoint(beamSpot_.position().x(),
						      beamSpot_.position().y(),beamSpot_.position().z()));

  if ( !theDCAXBS.isValid() )  return false;

  DcaBs = theDCAXBS.perigeeParameters().transverseImpactParameter();
  DcaBsErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  //if ( fabs(DcaBs/DcaBsErr) < TrkMinDcaSigBs_ )   return false;
  if (!(DcaBsErr > 0)) return false;
  //cout<<" DcaBsErr "<<DcaBsErr<<endl;
  return true;
}


bool
BsToPhiMuMu::hasGoodTrackDcaPoint (const reco::TransientTrack track,
				    const GlobalPoint p,
				    double maxdca, double &dca, double &dcaerr)
{
  TrajectoryStateClosestToPoint theDCAX = track.trajectoryStateClosestToPoint(p);
  if ( !theDCAX.isValid() ) return false;

  dca = theDCAX.perigeeParameters().transverseImpactParameter();
  dcaerr = theDCAX.perigeeError().transverseImpactParameterError();
  if ( dca > maxdca ) return false;

  return true;
}


bool
BsToPhiMuMu::calClosestApproachTracks (const reco::TransientTrack trackpTT,
					const reco::TransientTrack trackmTT,
					double & trk_R,
					double & trk_Z,
					double & trk_DCA)
{
  ClosestApproachInRPhi ClosestApp;
  ClosestApp.calculate(trackpTT.initialFreeState(),
		       trackmTT.initialFreeState());
  if (! ClosestApp.status() )  return false ;

  GlobalPoint XingPoint = ClosestApp.crossingPoint();

  trk_R = sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y());
  trk_Z = fabs(XingPoint.z());

  // if ((sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) >
  //      TrkMaxR_) || (fabs(XingPoint.z()) > TrkMaxZ_))  return false;

  trk_DCA = ClosestApp.distance();
  // if (DCAmumu > MuMuMaxDca_) return false;

  return true;
}


bool
BsToPhiMuMu::hasGoodMuMuVertex ( const reco::TransientTrack muTrackpTT,
				 const reco::TransientTrack muTrackmTT,
				 reco::TransientTrack &refitMupTT,
				 reco::TransientTrack &refitMumTT,
				 double & mu_mu_vtx_cl, double & mu_mu_pt,
				 double & mu_mu_mass, double & mu_mu_mass_err,
				 double & MuMuLSBS, double & MuMuLSBSErr,
				 double & MuMuCosAlphaBS,
				 double & MuMuCosAlphaBSErr)
{
  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;

  vector<RefCountedKinematicParticle> muonParticles;
  double chi = 0.;
  double ndf = 0.;
  muonParticles.push_back(partFactory.particle(muTrackmTT,
					       MuonMass_,chi,ndf,MuonMassErr_));
  muonParticles.push_back(partFactory.particle(muTrackpTT,
					       MuonMass_,chi,ndf,MuonMassErr_));

  RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles);

  if ( !mumuVertexFitTree->isValid())  return false;

  mumuVertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();
  RefCountedKinematicVertex mumu_KV = mumuVertexFitTree->currentDecayVertex();

  if ( !mumu_KV->vertexIsValid()) return false;

  mu_mu_vtx_cl = TMath::Prob((double)mumu_KV->chiSquared(),
			     int(rint(mumu_KV->degreesOfFreedom())));

  if (mu_mu_vtx_cl < MuMuMinVtxCl_)  return false;

  // extract the re-fitted tracks
  mumuVertexFitTree->movePointerToTheTop();

  mumuVertexFitTree->movePointerToTheFirstChild();
  RefCountedKinematicParticle refitMum = mumuVertexFitTree->currentParticle();
  refitMumTT = refitMum->refittedTransientTrack();

  mumuVertexFitTree->movePointerToTheNextChild();
  RefCountedKinematicParticle refitMup = mumuVertexFitTree->currentParticle();
  refitMupTT = refitMup->refittedTransientTrack();

  TLorentzVector mymum, mymup, mydimu;

  mymum.SetXYZM(refitMumTT.track().momentum().x(),
		refitMumTT.track().momentum().y(),
		refitMumTT.track().momentum().z(), MuonMass_);

  mymup.SetXYZM(refitMupTT.track().momentum().x(),
		refitMupTT.track().momentum().y(),
		refitMupTT.track().momentum().z(), MuonMass_);

  mydimu = mymum + mymup;
  mu_mu_pt = mydimu.Perp();

  mu_mu_mass = mumu_KP->currentState().mass();
  mu_mu_mass_err = sqrt(mumu_KP->currentState().kinematicParametersError().
			matrix()(6,6));

  if ((mu_mu_pt < MuMuMinPt_) || 
      (mu_mu_mass < MuMuMinInvMass_) ||
      (mu_mu_mass > MuMuMaxInvMass_))  return false;
  //(mu_mu_mass < MuMuMinInvMassl1_) || (mu_mu_mass > MuMuMaxInvMassl2_) || (mu_mu_mass < MuMuMinInvMassl2_ && mu_mu_mass > MuMuMaxInvMassl1_)) return false;

  //cout<<"mumumass in goodmumucvertex "<<mu_mu_mass<<endl;
  // compute the distance between mumu vtx and beam spot
  calLS (mumu_KV->position().x(),mumu_KV->position().y(),0.0,
	 beamSpot_.position().x(),beamSpot_.position().y(),0.0,
	 mumu_KV->error().cxx(),mumu_KV->error().cyy(),0.0,
	 mumu_KV->error().matrix()(0,1),0.0,0.0,
	 beamSpot_.covariance()(0,0),beamSpot_.covariance()(1,1),0.0,
	 beamSpot_.covariance()(0,1),0.0,0.0,
	 &MuMuLSBS,&MuMuLSBSErr);

  if (MuMuLSBS/MuMuLSBSErr < MuMuMinLxySigmaBs_)  return false;

  calCosAlpha(mumu_KP->currentState().globalMomentum().x(),
	      mumu_KP->currentState().globalMomentum().y(),
	      0.0,
	      mumu_KV->position().x() - beamSpot_.position().x(),
	      mumu_KV->position().y() - beamSpot_.position().y(),
	      0.0,
	      mumu_KP->currentState().kinematicParametersError().matrix()(3,3),
	      mumu_KP->currentState().kinematicParametersError().matrix()(4,4),
	      0.0,
	      mumu_KP->currentState().kinematicParametersError().matrix()(3,4),
	      0.0,
	      0.0,
	      mumu_KV->error().cxx() + beamSpot_.covariance()(0,0),
	      mumu_KV->error().cyy() + beamSpot_.covariance()(1,1),
	      0.0,
	      mumu_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),
	      0.0,
	      0.0,
	      &MuMuCosAlphaBS,&MuMuCosAlphaBSErr);

  if (MuMuCosAlphaBS < MuMuMinCosAlphaBs_)  return false;

  return true;
}


bool
BsToPhiMuMu::matchMuonTrack (const edm::Event& iEvent,
			      const reco::TrackRef theTrackRef)
{
  if ( theTrackRef.isNull() ) return false;

  // edm::Handle<pat::MuonCollection> thepatMuonHandle;
  // iEvent.getByToken(MuonLabel_, thepatMuonHandle);
  edm::Handle< edm::View<pat::Muon> > patMuonHandle;
  iEvent.getByToken(MuonLabel_,patMuonHandle);  
  // reco::TrackRef muTrackRef;
  for (const pat::Muon &iMuontmp : *patMuonHandle) {    
    for (unsigned int i = 0; i < iMuontmp.numberOfSourceCandidatePtrs(); ++i) {      
      const edm::Ptr<reco::Candidate> & source = iMuontmp.sourceCandidatePtr(i);
      if (! (iMuontmp.sourceCandidatePtr(i)).isNonnull())   continue;
      if (! (iMuontmp.sourceCandidatePtr(i)).isAvailable()) continue;
      
      const reco::Candidate & cand = *(source);
      if (cand.charge() == 0) continue;
      if (cand.bestTrack() == nullptr) continue;
      try{ cand.bestTrack()->eta();}
      catch(...) { std::cout << "should continue: " << std::endl; continue;}
      
      if ( iMuontmp.charge() == -1 && theTrackRef->charge() == -1 &&
  	   deltaR(theTrackRef->eta(),theTrackRef->phi(),cand.bestTrack()->eta(),cand.bestTrack()->phi()) < 0.00001
  	   ){
	//std::cout<<"Eta "<<theTrackRef->eta()<<" phi "<<theTrackRef->phi()<<" muon eta "<<cand.bestTrack()->eta()<<" muon phi "<<cand.bestTrack()->phi()<<std::endl;
	//std::cout<<"it Matches with Negative Muon charge "<<std::endl;
  	return true;
      }
      else if ( iMuontmp.charge() == 1 && theTrackRef->charge() == 1 &&
  		deltaR(theTrackRef->eta(),theTrackRef->phi(),cand.bestTrack()->eta(),cand.bestTrack()->phi()) < 0.00001
  		){	
	//std::cout<<"it Matches with positive Muon charge "<<std::endl;
  	return true;
      }
    }
  }
  //std::cout<<"it does not Matches with Muon charge "<<std::endl;
  return false;
}
 

// bool BsToPhiMuMu::IsTheSame(const reco::TrackRef& tk, const pat::Muon& mu){
//   double DeltaEta = fabs(mu.eta()-tk.eta());
//   double DeltaP   = fabs(mu.p()-tk.p());
//   if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
//   return false;
// }

bool
BsToPhiMuMu::hasGoodTrack(const edm::Event& iEvent,
			  const reco::TrackRef iTrack,
			  double & kaon_trk_pt)
{
  //reco::TrackRef theTrackRef = iTrack.track();
  if ( iTrack.isNull() ) return false;

  // veto muon tracks
  if ( matchMuonTrack(iEvent, iTrack) ) return false;

  // veto pion tracks from Kshort
  //if ( matchKshortTrack(iEvent, theTrackRef) ) return false;

  // check the track kinematics

  return true;
}


bool
BsToPhiMuMu::hasGoodPhiVertex( const reco::TransientTrack kaonmTT,
			       const reco::TransientTrack kaonpTT,
			       reco::TransientTrack &refitKmTT,
			       reco::TransientTrack &refitKpTT,
			       double & phi_vtx_cl, double & phi_mass)
{
  KinematicParticleFactoryFromTransientTrack pFactory;

  float chi = 0.;
  float ndf = 0.;

  vector<RefCountedKinematicParticle> phiParticles;
  phiParticles.push_back(pFactory.particle(kaonmTT,KaonMass_,chi,ndf,KaonMassErr_));
  phiParticles.push_back(pFactory.particle(kaonpTT,KaonMass_,chi,ndf,KaonMassErr_));

  KinematicParticleVertexFitter fitter;
  RefCountedKinematicTree phiVertexFitTree = fitter.fit(phiParticles);
  if ( ! phiVertexFitTree->isValid() ) return false ;

  phiVertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle phi_KP = phiVertexFitTree->currentParticle();
  RefCountedKinematicVertex phi_KV   = phiVertexFitTree->currentDecayVertex();
  if ( !phi_KV->vertexIsValid() ) return false;

  phi_vtx_cl = TMath::Prob((double)phi_KV->chiSquared(),
			   int(rint(phi_KV->degreesOfFreedom())));

  // NS @ added 2016-10-28
  //phivtxx->push_back(phi_KV->position().x());
  //phivtxy->push_back(phi_KV->position().y());
  //phivtxz->push_back(phi_KV->position().z());

  // extract the re-fitted tracks
  phiVertexFitTree->movePointerToTheTop();

  phiVertexFitTree->movePointerToTheFirstChild();
  RefCountedKinematicParticle refitKm = phiVertexFitTree->currentParticle();
  refitKmTT = refitKm->refittedTransientTrack();

  phiVertexFitTree->movePointerToTheNextChild();
  RefCountedKinematicParticle refitKp = phiVertexFitTree->currentParticle();
  refitKpTT = refitKp->refittedTransientTrack();

  TLorentzVector mykm, mykp, myphi;

  mykm.SetXYZM(refitKmTT.track().momentum().x(),
	       refitKmTT.track().momentum().y(),
	       refitKmTT.track().momentum().z(), KaonMass_);

  mykp.SetXYZM(refitKpTT.track().momentum().x(),
	       refitKpTT.track().momentum().y(),
	       refitKpTT.track().momentum().z(), KaonMass_);

  myphi = mykm + mykp;
  //phi_pt = myphi.Perp();

  phi_mass = phi_KP->currentState().mass();
  //phi_mass_err = sqrt(phi_KP->currentState().kinematicParametersError().matrix()(6,6));

  //if ( (phi_pt < PhiMinPt_) || (phi_mass < PhiMinMass_) || (phi_mass > PhiMaxMass_) )  return false;

  return true;

}

bool
BsToPhiMuMu::hasGoodBsVertex(const reco::TransientTrack mu1TT,
			     const reco::TransientTrack mu2TT,
			     const reco::TransientTrack kaonmTT,
			     const reco::TransientTrack kaonpTT,
			     double & b_vtx_chisq, double & b_vtx_cl,
			     double & b_mass, 
			     RefCountedKinematicTree & vertexFitTree)
{

  KinematicParticleFactoryFromTransientTrack pFactory;
  float chi = 0.;
  float ndf = 0.;

  // Bs -> mu+ mu- phi (K+ K-)
  vector<RefCountedKinematicParticle> vFitMCParticles;
  vFitMCParticles.push_back(pFactory.particle(mu1TT,MuonMass_,
					      chi,ndf,MuonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(mu2TT,MuonMass_,
					      chi,ndf,MuonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(kaonmTT, KaonMass_, chi,
					      ndf, KaonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(kaonpTT, KaonMass_, chi,
					      ndf, KaonMassErr_));


  KinematicParticleVertexFitter fitter;
  vertexFitTree = fitter.fit(vFitMCParticles);
  if (!vertexFitTree->isValid()) return false;
  ////cout << "particles fitted to a single vertex found: " << boolalpha << vertexFitTree->isValid() << endl; 

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();

  if ( !b_KV->vertexIsValid()) return false;
  ////cout << "Bs decay vertex found: " << boolalpha << b_KV->vertexIsValid() << endl;

  b_vtx_cl = TMath::Prob((double)b_KV->chiSquared(),
			 int(rint(b_KV->degreesOfFreedom())));

  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  b_mass = b_KP->currentState().mass();

  //if ( (b_vtx_cl < BsMinVtxCl_) || (b_mass < BsMinMass_) || (b_mass > BsMaxMass_) ) return false;

  //printf("reco Bs cand vtxcl: %6.4f , mass: %6.4f \n", b_vtx_cl, b_mass); 
  
  return true;

}
bool
BsToPhiMuMu::calClosestApproachBs (const RefCountedKinematicTree vertexFitTree,  double & DCABsBS, double & DCABsBSErr){
  vertexFitTree->movePointerToTheTop(); // Bs --> phi(KK) mu+ mu-                                                                                                
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();

  TrajectoryStateClosestToPoint theDCAXBS =b_KP->refittedTransientTrack().trajectoryStateClosestToPoint(GlobalPoint(beamSpot_.position().x(),beamSpot_.position().y(),beamSpot_.position().z()));
  if (theDCAXBS.isValid() == false) return false ;

  DCABsBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
  DCABsBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();

  return true; 
}

void
BsToPhiMuMu::saveBsToPhiMuMu(const RefCountedKinematicTree vertexFitTree){

  vertexFitTree->movePointerToTheTop(); // Bs --> phi(KK) mu+ mu-                                                                                         
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();

  bpx->push_back(b_KP->currentState().globalMomentum().x());
  bpxerr->push_back( sqrt( b_KP->currentState().kinematicParametersError().matrix()(3,3) ) );
  bpy->push_back(b_KP->currentState().globalMomentum().y());
  bpyerr->push_back( sqrt( b_KP->currentState().kinematicParametersError().matrix()(4,4) ) );
  bpz->push_back(b_KP->currentState().globalMomentum().z());
  bpzerr->push_back( sqrt( b_KP->currentState().kinematicParametersError().matrix()(5,5) ) );
  bmass->push_back(b_KP->currentState().mass());
  bmasserr->push_back( sqrt( b_KP->currentState().kinematicParametersError().matrix()(6,6) ) );

  vertexFitTree->movePointerToTheFirstChild(); // mu1                                                                                                      
  RefCountedKinematicParticle mu1_KP = vertexFitTree->currentParticle();
  vertexFitTree->movePointerToTheNextChild();  // mu2                                                                                                         
  RefCountedKinematicParticle mu2_KP = vertexFitTree->currentParticle();

  RefCountedKinematicParticle mup_KP, mum_KP ;

  if ( mu1_KP->currentState().particleCharge() > 0 ) mup_KP = mu1_KP;
  if ( mu1_KP->currentState().particleCharge() < 0 ) mum_KP = mu1_KP;
  if ( mu2_KP->currentState().particleCharge() > 0 ) mup_KP = mu2_KP;
  if ( mu2_KP->currentState().particleCharge() < 0 ) mum_KP = mu2_KP;

  muppx->push_back(mup_KP->currentState().globalMomentum().x());
  muppy->push_back(mup_KP->currentState().globalMomentum().y());
  muppz->push_back(mup_KP->currentState().globalMomentum().z());

  mumpx->push_back(mum_KP->currentState().globalMomentum().x());
  mumpy->push_back(mum_KP->currentState().globalMomentum().y());
  mumpz->push_back(mum_KP->currentState().globalMomentum().z());

  // add the variables for K+ and K-  ??

  vertexFitTree->movePointerToTheNextChild();
  RefCountedKinematicParticle k1_KP = vertexFitTree->currentParticle();
  vertexFitTree->movePointerToTheNextChild();
  RefCountedKinematicParticle k2_KP = vertexFitTree->currentParticle();

  RefCountedKinematicParticle kp_KP, km_KP ;

  if ( k1_KP->currentState().particleCharge() > 0 ) kp_KP = k1_KP;
  if ( k1_KP->currentState().particleCharge() < 0 ) km_KP = k1_KP;
  if ( k2_KP->currentState().particleCharge() > 0 ) kp_KP = k2_KP;
  if ( k2_KP->currentState().particleCharge() < 0 ) km_KP = k2_KP;

  kpchg->push_back(kp_KP->currentState().particleCharge());
  kppx->push_back(kp_KP->currentState().globalMomentum().x());
  kppy->push_back(kp_KP->currentState().globalMomentum().y());
  kppz->push_back(kp_KP->currentState().globalMomentum().z());

  kmchg->push_back(km_KP->currentState().particleCharge());
  kmpx->push_back(km_KP->currentState().globalMomentum().x());
  kmpy->push_back(km_KP->currentState().globalMomentum().y());
  kmpz->push_back(km_KP->currentState().globalMomentum().z());


}

void
BsToPhiMuMu::saveBsVertex(RefCountedKinematicTree vertexFitTree){
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();
  bvtxx->push_back((*b_KV).position().x());
  bvtxxerr->push_back(sqrt( abs(b_KV->error().cxx()) ));
  bvtxy->push_back((*b_KV).position().y());
  bvtxyerr->push_back(sqrt( abs(b_KV->error().cyy()) ));
  bvtxz->push_back((*b_KV).position().z());
  bvtxzerr->push_back(sqrt( abs(b_KV->error().czz()) ));

}

void 
BsToPhiMuMu::saveBsCosAlpha(RefCountedKinematicTree vertexFitTree)
{
  // alpha is the angle in the transverse plane between the B0 momentum                                                                             
  // and the seperation between the B0 vertex and the beamspot                                                                                           

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex   b_KV = vertexFitTree->currentDecayVertex();

  double cosAlphaBS, cosAlphaBSErr;
  calCosAlpha(b_KP->currentState().globalMomentum().x(),
              b_KP->currentState().globalMomentum().y(),
              b_KP->currentState().globalMomentum().z(),
              b_KV->position().x() - beamSpot_.position().x(),
              b_KV->position().y() - beamSpot_.position().y(),
              b_KV->position().z() - beamSpot_.position().z(),
              b_KP->currentState().kinematicParametersError().matrix()(3,3),
              b_KP->currentState().kinematicParametersError().matrix()(4,4),
              b_KP->currentState().kinematicParametersError().matrix()(5,5),
              b_KP->currentState().kinematicParametersError().matrix()(3,4),
              b_KP->currentState().kinematicParametersError().matrix()(3,5),
              b_KP->currentState().kinematicParametersError().matrix()(4,5),
              b_KV->error().cxx() + beamSpot_.covariance()(0,0),
              b_KV->error().cyy() + beamSpot_.covariance()(1,1),
              b_KV->error().czz() + beamSpot_.covariance()(2,2),
              b_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),
              b_KV->error().matrix()(0,2) + beamSpot_.covariance()(0,2),
              b_KV->error().matrix()(1,2) + beamSpot_.covariance()(1,2),
              &cosAlphaBS,&cosAlphaBSErr);


  bcosalphabs->push_back(cosAlphaBS);
  bcosalphabserr->push_back(cosAlphaBSErr);

}


void
BsToPhiMuMu::saveBsCosAlpha2d(RefCountedKinematicTree vertexFitTree)
{

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex   b_KV = vertexFitTree->currentDecayVertex();

  double cosAlphaBS2d, cosAlphaBS2dErr;
  calCosAlpha2d(b_KP->currentState().globalMomentum().x(),
                b_KP->currentState().globalMomentum().y(),0.0,                   
                b_KV->position().x() - beamSpot_.position().x(),
                b_KV->position().y() - beamSpot_.position().y(),0.0,
                b_KP->currentState().kinematicParametersError().matrix()(3,3),
                b_KP->currentState().kinematicParametersError().matrix()(4,4),0.0,
                b_KP->currentState().kinematicParametersError().matrix()(3,4),0.0,0.0,
                b_KV->error().cxx() + beamSpot_.covariance()(0,0),
                b_KV->error().cyy() + beamSpot_.covariance()(1,1),0.0,
                b_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),0.0,0.0,
                &cosAlphaBS2d,&cosAlphaBS2dErr);


  bcosalphabs2d->push_back(cosAlphaBS2d);
  bcosalphabs2derr->push_back(cosAlphaBS2dErr);

}

void 
BsToPhiMuMu::saveBsLsig(RefCountedKinematicTree vertexFitTree)
{
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();
  double LSBS, LSBSErr;

  calLS (b_KV->position().x(), b_KV->position().y(), 0.0,
         beamSpot_.position().x(), beamSpot_.position().y(), 0.0,
         b_KV->error().cxx(), b_KV->error().cyy(), 0.0,
         b_KV->error().matrix()(0,1), 0.0, 0.0,
         beamSpot_.covariance()(0,0), beamSpot_.covariance()(1,1), 0.0,
         beamSpot_.covariance()(0,1), 0.0, 0.0,
         &LSBS,&LSBSErr);

  blsbs->push_back(LSBS);
  blsbserr->push_back(LSBSErr);

}

void
BsToPhiMuMu::calCtau(RefCountedKinematicTree vertexFitTree,
		     double &bctau, double &bctauerr)
{
  //calculate ctau = (mB*(Bvtx-Pvtx)*pB)/(|pB|**2)                                                                                                     

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex   b_KV = vertexFitTree->currentDecayVertex();

  double betagamma = (b_KP->currentState().globalMomentum().mag()/BsMass_);

  // calculate ctau error. Momentum error is negligible compared to                                                                                       
  // the vertex errors, so don't worry about it                                                                                                         

  GlobalPoint BVP = GlobalPoint( b_KV->position() );
  GlobalPoint PVP = GlobalPoint( primaryVertex_.position().x(),
                                 primaryVertex_.position().y(),
                                 primaryVertex_.position().z() );
  GlobalVector sep3D = BVP-PVP;
  GlobalVector pBV = b_KP->currentState().globalMomentum();
  bctau = (BsMass_* (sep3D.dot(pBV)))/(pBV.dot(pBV));

  GlobalError BVE = b_KV->error();
  GlobalError PVE = GlobalError( primaryVertex_.error() );
  VertexDistance3D theVertexDistance3D;
  Measurement1D TheMeasurement = theVertexDistance3D.distance( VertexState(BVP, BVE), VertexState(PVP, PVE) );
  double myError = TheMeasurement.error();

  //  ctau is defined by the portion of the flight distance along                                                                  
  //  the compoenent of the B momementum, so only consider the error                                                                                      
  //  of that component, too, which is accomplished by scaling by                                                                                        
  //  ((VB-VP)(dot)PB)/|VB-VP|*|PB|                                                                                                                       

  double scale = abs( (sep3D.dot(pBV))/(sep3D.mag()*pBV.mag()) );
  bctauerr =  (myError*scale)/betagamma;

}

double
BsToPhiMuMu::calEta (double Px, double Py, double Pz)
{
  double P = sqrt(Px*Px + Py*Py + Pz*Pz);
  return 0.5*log((P + Pz) / (P - Pz));
}

double
BsToPhiMuMu::calPhi (double Px, double Py, double Pz)
{
  double phi = atan(Py / Px);
  if (Px < 0 && Py < 0) phi = phi - PI;
  if (Px < 0 && Py > 0) phi = phi + PI;
  return phi;
}

double
BsToPhiMuMu::calEtaPhiDistance (double Px1, double Py1, double Pz1,
				double Px2, double Py2, double Pz2)
{
  double phi1 = calPhi (Px1,Py1,Pz1);
  double eta1 = calEta (Px1,Py1,Pz1);
  double phi2 = calPhi (Px2,Py2,Pz2);
  double eta2 = calEta (Px2,Py2,Pz2);
  return sqrt((eta1-eta2) * (eta1-eta2) + (phi1-phi2) * (phi1-phi2));
}

void 
BsToPhiMuMu::saveBsCtau(RefCountedKinematicTree vertexFitTree)
{
  double bctau_temp, bctauerr_temp;
  calCtau(vertexFitTree, bctau_temp, bctauerr_temp);
  bctau->push_back(bctau_temp);
  bctauerr->push_back(bctauerr_temp);
}
void
BsToPhiMuMu::saveIPvtx(const edm::Event& iEvent, const reco::TransientTrack mu1TT,
		    const reco::TransientTrack mu2TT,
		    const reco::TransientTrack kaonmTT,
		    const reco::TransientTrack kaonpTT)
{
  edm::Handle<reco::VertexCollection> recVtx;
  iEvent.getByToken(VertexLabel_, recVtx);
  if (recVtx->empty()) return; // skip the event if no PV found
  double mumMind0  = 100;
  double mupMind0  = 100;
  double TrkmMind0 = 100;
  double TrkpMind0 = 100;
  double mumMind0E, mupMind0E, TrkmMind0E, TrkpMind0E;
  double mumMinip  = 100;
  double mupMinip  = 100;
  double TrkmMinip = 100;
  double TrkpMinip = 100;
  double mumMinipE, mupMinipE, TrkmMinipE, TrkpMinipE;
  GlobalPoint vert;
  TrajectoryStateClosestToPoint traj;

  mumMind0  = mupMind0  = 100;
  TrkmMind0 = TrkpMind0 = 100;
  mumMinip  = mupMinip  = 100;
  TrkmMinip = TrkpMinip = 100;
  std::pair<double,double>  IPPair;
  
  for (std::vector<reco::Vertex>::const_iterator ipv = recVtx->begin(); ipv != recVtx->end(); ipv++) { 
    if (! ipv->isValid() ) continue; 
    vert = GlobalPoint(ipv->x(), ipv->y(), ipv->z());

    traj = mu1TT.trajectoryStateClosestToPoint(vert );
    if (fabs(traj.perigeeParameters().transverseImpactParameter()) < mumMind0){
      mumMind0  = fabs(traj.perigeeParameters().transverseImpactParameter());
      mumMind0E = traj.perigeeError().transverseImpactParameterError();
    }  
    
    traj = mu2TT.trajectoryStateClosestToPoint(vert );
    if (fabs(traj.perigeeParameters().transverseImpactParameter()) < mupMind0){
      mupMind0  = fabs(traj.perigeeParameters().transverseImpactParameter());
      mupMind0E = traj.perigeeError().transverseImpactParameterError();
    }  
    
    traj = kaonmTT.trajectoryStateClosestToPoint(vert );
    if (fabs(traj.perigeeParameters().transverseImpactParameter()) < TrkmMind0){
      TrkmMind0  = fabs(traj.perigeeParameters().transverseImpactParameter());
      TrkmMind0E = traj.perigeeError().transverseImpactParameterError();
    }  
    
    traj = kaonpTT.trajectoryStateClosestToPoint(vert );
    if (fabs(traj.perigeeParameters().transverseImpactParameter()) < TrkpMind0){
      TrkpMind0  = fabs(traj.perigeeParameters().transverseImpactParameter());
      TrkpMind0E = traj.perigeeError().transverseImpactParameterError();
    }  
    
    IPPair = pionImpactParameter(mu1TT,*ipv);
    if (IPPair.first < mumMinip){
      mumMinip  = IPPair.first;
      mumMinipE = IPPair.second;
    }
    IPPair = pionImpactParameter(mu2TT,*ipv);
    if (IPPair.first < mupMinip){
      mupMinip  = IPPair.first;
      mupMinipE = IPPair.second;
    }
    IPPair = pionImpactParameter(kaonmTT,*ipv);
    if (IPPair.first < TrkmMinip){
      TrkmMinip  = IPPair.first;
      TrkmMinipE = IPPair.second;
    }
    IPPair = pionImpactParameter(kaonpTT,*ipv);
    if (IPPair.first < TrkpMinip){
      TrkpMinip  = IPPair.first;
      TrkpMinipE = IPPair.second;
    }
  }
  mumMinIP2D -> push_back(mumMind0);
  mumMinIP2DE     -> push_back(mumMind0E);
  mupMinIP2D      -> push_back(mupMind0);
  mupMinIP2DE     -> push_back(mupMind0E);
  kmtrkMinIP2D  -> push_back(TrkmMind0);
  kmtrkMinIP2DE -> push_back(TrkmMind0E);
  kptrkMinIP2D  -> push_back(TrkpMind0);
  kptrkMinIP2DE -> push_back(TrkpMind0E);

  mumMinIP      -> push_back(mumMinip);
  mumMinIPE     -> push_back(mumMinipE);
  mupMinIP      -> push_back(mupMinip);
  mupMinIPE     -> push_back(mupMinipE);
  kmtrkMinIP  -> push_back(TrkmMinip);
  kmtrkMinIPE -> push_back(TrkmMinipE);
  kptrkMinIP  -> push_back(TrkpMinip);
  kptrkMinIPE -> push_back(TrkpMinipE);
 
}
std::pair<double,double> BsToPhiMuMu::pionImpactParameter(reco::TransientTrack piTT, reco::Vertex myVtx)
{
  std::pair<double,double> measure;
  std::pair<bool,Measurement1D>  piIP_pair = IPTools::absoluteImpactParameter3D(piTT, myVtx);
  if (piIP_pair.first)
    {
      measure.first  = piIP_pair.second.value();
      measure.second = piIP_pair.second.significance();
    }
  else 
    {
      //if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 3D" << std::endl;
      measure.first  = 0;
      measure.second = 0;
    } 
  return measure;
}
void
BsToPhiMuMu::saveGenInfo(const edm::Event& iEvent){

  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(prunedGenToken_,pruned);
  
  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packedGenToken_,packed);
  
  const reco::Candidate* genMum = NULL;
  const reco::Candidate* genMup = NULL;

  const reco::Candidate* genTrkm = NULL;
  const reco::Candidate* genTrkp = NULL;

  
  bool found_mum  = false;
  bool found_mup  = false;
  bool found_trkp = false;
  bool found_trkm = false;

  for (const reco::GenParticle &bMeson : *pruned) {
    if(abs( bMeson.pdgId()) == BS_PDG_ID){
      if (skipOscillations(bMeson, pruned)) continue;
      
      //   if (printMsg)
      //std::cout << "PdgID: " << bMeson.pdgId() << " pt " << bMeson.pt() << " eta: " << bMeson.eta() << " phi: " << bMeson.phi()  << "mother: " << bMeson.mother(0)->pdgId() << std::endl;

      genMum    = NULL;
      genMup    = NULL;
      genTrkm   = NULL;
      genTrkp   = NULL;
      
      found_mum   = false;
      found_mup   = false;
      found_trkp  = false;
      found_trkm  = false;
      
      for (const pat::PackedGenParticle &dau : *packed) {
	//get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
	const reco::Candidate * motherInPrunedCollection = dau.mother(0) ;
	if(motherInPrunedCollection != nullptr && isAncestor( &bMeson , motherInPrunedCollection))
	  {
	    //if (printMsg) 
	    //  std::cout << " After check  PdgID: " << dau.pdgId() << " pt " << dau.pt() << " eta: " << dau.eta() << " phi: " << dau.phi() << std::endl;
	    
	    if (dau.pdgId() == 13){
	      found_mum = true;
	      genMum = &dau;
	    }
	    else if (dau.pdgId() == -13) {
	      found_mup = true;
	      genMup = &dau;
	    }
	    else if (dau.pdgId() == 321 ) {
	      found_trkp = true;
	      genTrkp = &dau;
	    }
	    else if (dau.pdgId() == -321 ) {
	      found_trkm = true;
	      genTrkm = &dau;
	    }
            
	  }
      }
      

      if (found_mup && found_mum && found_trkp && found_trkm ){
	// save gen info                                                                                                                                     
	genbpx = bMeson.px();
	genbpy = bMeson.py();
	genbpz = bMeson.pz();
	genbpid = bMeson.pdgId();
	//std::cout<<"checking pdgid final "<<bMeson.pdgId()<<std::endl;
	// genphipx = phi.px();
	// genphipy = phi.py();
	// genphipz = phi.pz();
	
	// genphivtxx = phi.vx();
	// genphivtxy = phi.vy();
	// genphivtxz = phi.vz();
	
	genkpchg = genTrkp->charge();
	genkppx  = genTrkp->px();
	genkppy  = genTrkp->py();
	genkppz  = genTrkp->pz();
	
	genkmchg = genTrkm->charge();
	genkmpx  = genTrkm->px();
	genkmpy  = genTrkm->py();
	genkmpz  = genTrkm->pz();
	
	genmumpx = genMum->px();
	genmumpy = genMum->py();
	genmumpz = genMum->pz();
	
	genmuppx = genMup->px();
	genmuppy = genMup->py();
	genmuppz = genMup->pz();
      }
    }
  }
}





void
BsToPhiMuMu::saveDimuVariables(double DCAmumBS, double DCAmumBSErr,
			       double DCAmupBS, double DCAmupBSErr,
			       double mumutrk_R, double mumutrk_Z,
			       double DCAmumu,  double mu_mu_vtx_cl,
			       double MuMuLSBS, double MuMuLSBSErr,
			       double MuMuCosAlphaBS, double MuMuCosAlphaBSErr,
			       double mu_mu_mass, double mu_mu_mass_err)

{
  mumdcabs->push_back(DCAmumBS);
  mumdcabserr->push_back(DCAmumBSErr);

  mupdcabs->push_back(DCAmupBS);
  mupdcabserr->push_back(DCAmupBSErr);

  mumutrkr->push_back(mumutrk_R);
  mumutrkz->push_back(mumutrk_Z);
  mumudca->push_back(DCAmumu);
  mumuvtxcl->push_back(mu_mu_vtx_cl);
  mumulsbs->push_back(MuMuLSBS);
  mumulsbserr->push_back(MuMuLSBSErr);
  mumucosalphabs->push_back(MuMuCosAlphaBS);
  mumucosalphabserr->push_back(MuMuCosAlphaBSErr);

  mumumass->push_back(mu_mu_mass);
  //cout<<" mumumass in saving "<<mu_mu_mass<<endl;
  mumumasserr->push_back(mu_mu_mass_err);
}


void BsToPhiMuMu::savePUinMC(const edm::Event& iEvent){
  // #################################
  // # Save pileup information in MC #
  // #################################
  edm::Handle< std::vector<PileupSummaryInfo> > PupInfo;
  iEvent.getByToken(puToken_, PupInfo);
  edm::LumiReWeighting lumi_weights;
  //lumi_weights      = edm::LumiReWeighting("PileupMC_2018.root", "DataPileupHistogram2018_rereco.root", "input_Event/N_TrueInteractions", "pileup");
  lumi_weights      = edm::LumiReWeighting("PileupMC_2016.root", "DataPileupHistogram2016_rereco.root", "input_Event/N_TrueInteractions", "pileup");
  //lumi_weights      = edm::LumiReWeighting("PileupMC_2017.root", "DataPileupHistogram2017_rereco.root", "input_Event/N_TrueInteractions", "pileup");
  float tnpv = -1;       // True number of primary vertices
  float wpu = 1;         // Pile-up re-weight factor
  
  for (std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); PVI++)
    {
      bunchXingMC->push_back(PVI->getBunchCrossing());
      numInteractionsMC->push_back(PVI->getPU_NumInteractions());
      trueNumInteractionsMC->push_back(PVI->getTrueNumInteractions());
      int bx = PVI->getBunchCrossing();
      if (bx == 0)tnpv = PVI->getTrueNumInteractions();        
    }

  wpu = lumi_weights.weight(tnpv);
  fpuw8 = wpu ;
}

bool BsToPhiMuMu::skipOscillations (const reco::GenParticle &bMeson, edm::Handle<reco::GenParticleCollection> pruned)
{
  for (unsigned int i = 0; i < bMeson.numberOfDaughters(); i++){
    ////if (bMeson.daughter(i)->pdgId() == 511 || bMeson.daughter(i)->pdgId() == 531 || bMeson.daughter(i)->pdgId() == 5122)    return true; 
    if (bMeson.daughter(i)->pdgId() == 511  || bMeson.daughter(i)->pdgId() == 5122)    return true; 
    //std::cout << "oscillating to:     PdgID: " << bMeson.daughter(i)->pdgId() << " pt " << bMeson.daughter(i)->pt() << " eta: " << bMeson.daughter(i)->eta() << " phi: " << bMeson.daughter(i)->phi() << std::endl;
  }
  
  for (const reco::GenParticle &bMother : *pruned) {
    if ( fabs(bMother.pdgId()) == 531){
      const reco::Candidate * mother = bMother.mother(0) ;
      //std::cout<<"mother PDGID  "<<bMother.mother(0)->pdgId()<<std::endl;
      if(mother != nullptr && isAncestor( &bMeson , mother)) return true;
    }
  }
  return false;
}


bool BsToPhiMuMu::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
    //particle is already the ancestor
    if(ancestor == particle ) return true;

    //otherwise loop on mothers, if any and return true if the ancestor is found
    for(size_t i=0;i< particle->numberOfMothers();i++)
    {
        if(isAncestor(ancestor,particle->mother(i))) return true;
    }
    //if we did not return yet, then particle and ancestor are not relatives
    return false;
}


void
BsToPhiMuMu::saveTruthMatch(const edm::Event& iEvent){
  double deltaEtaPhi;

  for (vector<int>::size_type i = 0; i < bmass->size(); i++) {//{{{
   
    //-----------------------
    // truth match with mu-
    //-----------------------
    deltaEtaPhi = calEtaPhiDistance(genmumpx, genmumpy, genmumpz,
				    mumpx->at(i), mumpy->at(i), mumpz->at(i));
    if (deltaEtaPhi < TruthMatchMuonMaxR_) {
      istruemum->push_back(true);
    } else {
      istruemum->push_back(false);
    }

    //-----------------------
    // truth match with mu+
    //-----------------------
    deltaEtaPhi = calEtaPhiDistance(genmuppx, genmuppy, genmuppz,
				    muppx->at(i), muppy->at(i), muppz->at(i));

    if (deltaEtaPhi < TruthMatchMuonMaxR_) {
      istruemup->push_back(true);
    }
    else {
      istruemup->push_back(false);
    }

    //---------------------------------
    // truth match with kaon+ track   
    //---------------------------------                                                                                                                       
    deltaEtaPhi = calEtaPhiDistance(genkppx, genkppy, genkppz,
                                    kppx->at(i), kppy->at(i), kppz->at(i));
    if (deltaEtaPhi < TruthMatchKaonMaxR_){
      istruekp->push_back(true);
    } else {
      istruekp->push_back(false);
    }

    //---------------------------------                                                                                                                           
    // truth match with kaon- track                                                                                                                                 
    //---------------------------------                                                                                                                              
    deltaEtaPhi = calEtaPhiDistance(genkmpx, genkmpy, genkmpz,
                                    kmpx->at(i), kmpy->at(i), kmpz->at(i));
    if (deltaEtaPhi < TruthMatchKaonMaxR_){
      istruekm->push_back(true);
    } else {
      istruekm->push_back(false);
    }


    //---------------------------------------
    // truth match with Bs or Bs bar 
    //---------------------------------------                                                                                                
    if ( istruemum->back() && istruemup->back() && istruekm->back() && istruekp->back() ) {
      istruebs->push_back(true);
      //std::cout<<"Getting correct Bs "<<std::endl;
    } else {
      istruebs->push_back(false);
      // std::cout<<"Getting wrong Bs "<<std::endl;
    }



  }//}}}

}

//define this as a plug-in
DEFINE_FWK_MODULE(BsToPhiMuMu);
