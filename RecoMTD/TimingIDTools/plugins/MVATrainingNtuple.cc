#include <memory>
#include "TTree.h"
#include "TFile.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
// #include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToTPAssociatorBaseImpl.h"
#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociationMap.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/PrimaryVertexProducer/interface/HITrackFilterForPVFinding.h"

using reco::TrackCollection;

class MVATrainingNtuple : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  typedef math::XYZTLorentzVector LorentzVector;

  // auxiliary class holding simulated vertices (originally from Primary4DVertexValidation)
  struct simPrimaryVertex {
    simPrimaryVertex(double x1, double y1, double z1, double t1, int k1) : x(x1), y(y1), z(z1), t(t1), key(k1) {};
    double x, y, z, t;
    int key;
    int eventId;
    int bunchCrossing;
    TrackingVertexRef sim_vertex;
    int OriginalIndex = -1;
    bool is_LV;
  };

public:
  explicit MVATrainingNtuple(const edm::ParameterSet&);
  ~MVATrainingNtuple() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  const bool trkTPSelAll(const TrackingParticle&);
  const bool trkRecSel(const reco::TrackBase&);

  const edm::Ref<std::vector<TrackingParticle>>* getAnyMatchedTP(const reco::TrackBaseRef&);
  double timeFromTrueMass(double, double, double, double);

  std::vector<MVATrainingNtuple::simPrimaryVertex> getSimPVs(const edm::Handle<TrackingVertexCollection>&);

  edm::Service<TFileService> fs_;

  // ----------member data ---------------------------
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;
  TrackFilterForPVFindingBase* theTrackFilter;

  static constexpr double simUnit_ = 1e9;     //sim time in s while reco time in ns
  static constexpr double c_ = 2.99792458e1;  //c in cm/ns
  static constexpr double etacutGEN_ = 4.;    // |eta| < 4;
  static constexpr double etacutREC_ = 3.;    // |eta| < 3;
  static constexpr double etacutBTL_ = 1.5;   // |eta| < 1.5;
  static constexpr double pTcutBTL_ = 0.7;    // PT > 0.7 GeV
  static constexpr double pTcutETL_ = 0.2;    // PT > 0.2 GeV
  static constexpr double rBTL_ = 110.0;
  static constexpr double zETL_ = 290.0;
  std::string fileName_;

  bool saveNtupleforBDT_;
  bool saveNtupleforGNN_;

  const reco::RecoToSimCollection* r2s_;

  // GNN input variables
  std::vector<double> gnn_pt, gnn_eta, gnn_phi, gnn_z_pca, gnn_dz, gnn_t_Pi, gnn_t_K, gnn_t_P, gnn_t0safe, gnn_t0pid,
      gnn_sigma_t0safe, gnn_mtdTime, gnn_sigma_tmtd, gnn_mva_qual, gnn_btlMatchChi2, gnn_btlMatchTimeChi2,
      gnn_etlMatchChi2, gnn_etlMatchTimeChi2, gnn_pathLength, gnn_probPi, gnn_probK, gnn_probP, gnn_trk_chi2,
      gnn_trk_ndof, gnn_sigma_tof_Pi, gnn_sigma_tof_K, gnn_sigma_tof_P, gnn_sim_vertex_z, gnn_sim_vertex_t, gnn_tp_tEst,
      gnn_outermostHitPosition;
  std::vector<int> gnn_npixBarrel, gnn_npixEndcap, gnn_sim_vertex_evID, gnn_sim_vertex_BX, gnn_sim_vertex_index,
      gnn_tp_pdgId, gnn_trk_validhits;
  std::vector<bool> gnn_is_matched_tp, gnn_sim_vertex_isLV;

  // BDT input variables
  std::vector<double> Ttrack_pt, Ttrack_eta, Ttrack_phi, Ttrack_dz, Ttrack_dxy, Ttrack_chi2, Ttrack_BTLchi2,
      Ttrack_BTLtime_chi2, Ttrack_ETLchi2, Ttrack_ETLtime_chi2, Ttrack_t0, Ttrack_sigmat0, Ttrack_Tmtd,
      Ttrack_sigmaTmtd, Ttrack_length, Ttrack_MtdMVA, Ttrack_lHitPos;
  std::vector<int> Ttrack_ndof, Ttrack_nValidHits, Ttrack_npixBarrelValidHits, Ttrack_npixEndcapValidHits;
  std::vector<bool> Ttrack_Signal, Ttrack_Associated, Ttrack_TPDirClu, Ttrack_TPOtherClu, Ttrack_isBTL, Ttrack_isETL,
      Ttrack_withMTD;

  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchTimeChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> etlMatchChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> etlMatchTimeChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<int>> npixBarrelToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> npixEndcapToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> outermostHitPositionToken_;

  edm::EDGetTokenT<reco::TrackCollection> RecTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> RecMTDTrackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> RecVertexToken_;
  edm::EDGetTokenT<reco::BeamSpot> RecBeamSpotToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;
  edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoAssociationToken_;
  edm::EDGetTokenT<TrackingVertexCollection> trackingVertexCollectionToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> trackAssocToken_;
  edm::EDGetTokenT<reco::TPToSimCollectionMtd> tp2SimAssociationMapToken_;
  edm::EDGetTokenT<MtdRecoClusterToSimLayerClusterAssociationMap> r2sAssociationMapToken_;
  edm::EDGetTokenT<FTLClusterCollection> btlRecCluToken_;
  edm::EDGetTokenT<FTLClusterCollection> etlRecCluToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> pathLengthToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> momentumToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatimeToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmat0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmat0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmat0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> trackMVAQualToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofpiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofkToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofpToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probPToken_;
};

MVATrainingNtuple::MVATrainingNtuple(const edm::ParameterSet& iConfig)
    : theTTBToken(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      fileName_(iConfig.getUntrackedParameter<std::string>("fileName")) {
  RecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTracks"));
  RecMTDTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagT"));
  RecVertexToken_ = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("inputTagV"));
  tp2SimAssociationMapToken_ =
      consumes<reco::TPToSimCollectionMtd>(iConfig.getParameter<edm::InputTag>("tp2SimAssociationMapTag"));
  r2sAssociationMapToken_ = consumes<MtdRecoClusterToSimLayerClusterAssociationMap>(
      iConfig.getParameter<edm::InputTag>("r2sAssociationMapTag"));
  trackAssocToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackAssocSrc"));
  RecBeamSpotToken_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("offlineBS"));
  trackingParticleCollectionToken_ =
      consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  recoToSimAssociationToken_ =
      consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  simToRecoAssociationToken_ =
      consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  trackingVertexCollectionToken_ = consumes<TrackingVertexCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  btlRecCluToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("recCluTagBTL"));
  etlRecCluToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("recCluTagETL"));
  pathLengthToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pathLengthSrc"));
  momentumToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("momentumSrc"));
  sigmatimeToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmaSrc"));
  t0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"));
  sigmat0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"));
  t0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0PID"));
  sigmat0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0PID"));
  t0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0SafePID"));
  sigmat0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0SafePID"));
  trackMVAQualToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMVAQual"));
  tmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtd"));
  tofPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofPi"));
  tofKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofK"));
  tofPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofP"));
  probPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probPi"));
  probKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probK"));
  probPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probP"));
  sigmatofpiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofpiSrc"));
  sigmatofkToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofkSrc"));
  sigmatofpToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofpSrc"));
  btlMatchChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchChi2Src"));
  btlMatchTimeChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchTimeChi2Src"));
  etlMatchChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchChi2Src"));
  etlMatchTimeChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchTimeChi2Src"));
  npixBarrelToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("npixBarrelSrc"));
  npixEndcapToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("npixEndcapSrc"));
  outermostHitPositionToken_ =
      consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("outermostHitPositionSrc"));
  saveNtupleforBDT_ = iConfig.getParameter<bool>("ntupleforBDT");
  saveNtupleforGNN_ = iConfig.getParameter<bool>("ntupleforGNN");
  // select and configure the track selection
  std::string trackSelectionAlgorithm =
      iConfig.getParameter<edm::ParameterSet>("TkFilterParameters").getParameter<std::string>("algorithm");
  if (trackSelectionAlgorithm == "filter") {
    theTrackFilter = new TrackFilterForPVFinding(iConfig.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else if (trackSelectionAlgorithm == "filterWithThreshold") {
    theTrackFilter = new HITrackFilterForPVFinding(iConfig.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else {
    edm::LogWarning("MVATrainingNtuple: unknown track selection algorithm: " + trackSelectionAlgorithm);
  }
}

MVATrainingNtuple::~MVATrainingNtuple() {
  if (theTrackFilter)
    delete theTrackFilter;
}

const edm::Ref<std::vector<TrackingParticle>>* MVATrainingNtuple::getAnyMatchedTP(const reco::TrackBaseRef& recoTrack) {
  auto found = r2s_->find(recoTrack);

  // reco track not matched to any TP
  if (found == r2s_->end())
    return nullptr;

  //matched TP equal to any TP
  for (const auto& tp : found->val) {
    return &tp.first;
  }

  // reco track not matched to any TP from vertex
  return nullptr;
}

double MVATrainingNtuple::timeFromTrueMass(double mass, double pathlength, double momentum, double time) {
  if (time > 0 && pathlength > 0 && mass > 0) {
    double gammasq = 1. + momentum * momentum / (mass * mass);
    double v = c_ * std::sqrt(1. - 1. / gammasq);  // cm / ns
    double t_est = time - (pathlength / v);

    return t_est;
  } else {
    return -1;
  }
}

std::vector<MVATrainingNtuple::simPrimaryVertex> MVATrainingNtuple::getSimPVs(
    const edm::Handle<TrackingVertexCollection>& tVC) {
  std::vector<MVATrainingNtuple::simPrimaryVertex> simpv;

  int current_event = -1;
  int s = -1;
  for (TrackingVertexCollection::const_iterator v = tVC->begin(); v != tVC->end(); ++v) {
    // LV is the first vertex in each event, keep only at BX=0
    int eventId = v->eventId().event();
    int bunchCrossing = v->eventId().bunchCrossing();

    if (bunchCrossing != 0)
      continue;

    bool is_LV = true;
    if (eventId != current_event) {
      current_event = eventId;
    } else {
      is_LV = false;
    }
    s++;

    // could be a new vertex, check  all primaries found so far to avoid multiple entries
    int key = std::distance(tVC->begin(), v);
    simPrimaryVertex sv(v->position().x(), v->position().y(), v->position().z(), v->position().t(), key);
    sv.eventId = eventId;
    sv.bunchCrossing = bunchCrossing;
    sv.sim_vertex = TrackingVertexRef(tVC, key);
    sv.OriginalIndex = s;
    sv.is_LV = is_LV;

    simPrimaryVertex* vp = nullptr;  // will become non-NULL if a vertex is found and then point to it
    for (std::vector<simPrimaryVertex>::iterator v0 = simpv.begin(); v0 != simpv.end(); v0++) {
      if ((sv.eventId == v0->eventId) && (sv.bunchCrossing == v0->bunchCrossing) && (std::abs(sv.x - v0->x) < 1e-5) &&
          (std::abs(sv.y - v0->y) < 1e-5) && (std::abs(sv.z - v0->z) < 1e-5)) {
        vp = &(*v0);
        break;
      }
    }
    if (!vp) {
      // this is a new vertex, add it to the list of sim-vertices
      simpv.push_back(sv);
    }

  }  // End of for loop on tracking vertices

  // In case of no simulated vertices, break here
  if (simpv.empty())
    return simpv;

  return simpv;
}

// ------------ method called for each event  ------------
void MVATrainingNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using edm::Handle;
  using std::vector;
  using namespace reco;

  edm::Handle<reco::TrackCollection> tracksH;
  iEvent.getByToken(RecTrackToken_, tracksH);

  const auto& theB = &iSetup.getData(theTTBToken);
  std::vector<reco::TransientTrack> t_tks;

  edm::Handle<TrackingParticleCollection> TPCollectionH;
  iEvent.getByToken(trackingParticleCollectionToken_, TPCollectionH);
  if (!TPCollectionH.isValid())
    edm::LogWarning("MVATrainingNtuple") << "TPCollectionH is not valid";

  edm::Handle<reco::RecoToSimCollection> recoToSimH;
  iEvent.getByToken(recoToSimAssociationToken_, recoToSimH);
  if (recoToSimH.isValid())
    r2s_ = recoToSimH.product();
  else
    edm::LogWarning("MVATrainingNtuple") << "recoToSimH is not valid";

  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> BeamSpotH;
  iEvent.getByToken(RecBeamSpotToken_, BeamSpotH);
  if (!BeamSpotH.isValid())
    edm::LogWarning("MVATrainingNtuple") << "BeamSpotH is not valid";
  beamSpot = *BeamSpotH;

  edm::Handle<TrackingVertexCollection> TVCollectionH;
  iEvent.getByToken(trackingVertexCollectionToken_, TVCollectionH);
  if (!TVCollectionH.isValid())
    edm::LogWarning("MVATrainingNtuple") << "TVCollectionH is not valid";

  std::vector<simPrimaryVertex> simpv;
  simpv = getSimPVs(TVCollectionH);

  const auto& trackAssoc = iEvent.get(trackAssocToken_);

  const auto& tp2SimAssociationMap = iEvent.get(tp2SimAssociationMapToken_);
  const auto& r2sAssociationMap = iEvent.get(r2sAssociationMapToken_);

  // auto RecTrackHandle = makeValid(iEvent.getHandle(RecTrackToken_));
  auto RecTrackHandle = iEvent.getHandle(RecTrackToken_);

  const auto& btlRecCluHandle = iEvent.getHandle(btlRecCluToken_);
  const auto& etlRecCluHandle = iEvent.getHandle(etlRecCluToken_);

  const auto& pathLength = iEvent.get(pathLengthToken_);
  const auto& momentum = iEvent.get(momentumToken_);
  const auto& sigmatimemtd = iEvent.get(sigmatimeToken_);
  const auto& t0Src = iEvent.get(t0SrcToken_);
  const auto& sigmat0Src = iEvent.get(sigmat0SrcToken_);
  const auto& t0Pid = iEvent.get(t0PidToken_);
  const auto& sigmat0Pid = iEvent.get(sigmat0PidToken_);
  const auto& t0Safe = iEvent.get(t0SafePidToken_);
  const auto& sigmat0Safe = iEvent.get(sigmat0SafePidToken_);
  const auto& mtdQualMVA = iEvent.get(trackMVAQualToken_);
  const auto& tMtd = iEvent.get(tmtdToken_);
  const auto& tofPi = iEvent.get(tofPiToken_);
  const auto& tofK = iEvent.get(tofKToken_);
  const auto& tofP = iEvent.get(tofPToken_);
  const auto& probPi = iEvent.get(probPiToken_);
  const auto& probK = iEvent.get(probKToken_);
  const auto& probP = iEvent.get(probPToken_);
  const auto& sigmatofpi = iEvent.get(sigmatofpiToken_);
  const auto& sigmatofk = iEvent.get(sigmatofkToken_);
  const auto& sigmatofp = iEvent.get(sigmatofpToken_);
  const auto& btlMatchChi2 = iEvent.get(btlMatchChi2Token_);
  const auto& btlMatchTimeChi2 = iEvent.get(btlMatchTimeChi2Token_);
  const auto& etlMatchChi2 = iEvent.get(etlMatchChi2Token_);
  const auto& etlMatchTimeChi2 = iEvent.get(etlMatchTimeChi2Token_);
  const auto& npixBarrel = iEvent.get(npixBarrelToken_);
  const auto& npixEndcap = iEvent.get(npixEndcapToken_);
  const auto& outermostHitPosition = iEvent.get(outermostHitPositionToken_);

  // Fill TTree with input variables for GNN
  if (saveNtupleforGNN_) {
    std::string GNNtreeName = "GNNtree_" + std::to_string(iEvent.id().event());
    TTree* GNNtree = fs_->make<TTree>(GNNtreeName.c_str(), "Tree for GNN tracks");

    GNNtree->Branch("gnn_pt", &gnn_pt);
    GNNtree->Branch("gnn_eta", &gnn_eta);
    GNNtree->Branch("gnn_phi", &gnn_phi);
    GNNtree->Branch("gnn_z_pca", &gnn_z_pca);
    GNNtree->Branch("gnn_dz", &gnn_dz);
    GNNtree->Branch("gnn_t_Pi", &gnn_t_Pi);
    GNNtree->Branch("gnn_t_K", &gnn_t_K);
    GNNtree->Branch("gnn_t_P", &gnn_t_P);
    GNNtree->Branch("gnn_sigma_t0safe", &gnn_sigma_t0safe);
    GNNtree->Branch("gnn_sigma_tmtd", &gnn_sigma_tmtd);
    GNNtree->Branch("gnn_t0safe", &gnn_t0safe);
    GNNtree->Branch("gnn_t0pid", &gnn_t0pid);
    GNNtree->Branch("gnn_mva_qual", &gnn_mva_qual);
    GNNtree->Branch("gnn_btlMatchChi2", &gnn_btlMatchChi2);
    GNNtree->Branch("gnn_btlMatchTimeChi2", &gnn_btlMatchTimeChi2);
    GNNtree->Branch("gnn_etlMatchChi2", &gnn_etlMatchChi2);
    GNNtree->Branch("gnn_etlMatchTimeChi2", &gnn_etlMatchTimeChi2);
    GNNtree->Branch("gnn_pathLength", &gnn_pathLength);
    GNNtree->Branch("gnn_npixBarrel", &gnn_npixBarrel);
    GNNtree->Branch("gnn_npixEndcap", &gnn_npixEndcap);
    GNNtree->Branch("gnn_outermostHitPosition", &gnn_outermostHitPosition);
    GNNtree->Branch("gnn_mtdTime", &gnn_mtdTime);
    GNNtree->Branch("gnn_is_matched_tp", &gnn_is_matched_tp);
    GNNtree->Branch("gnn_tp_tEst", &gnn_tp_tEst);
    GNNtree->Branch("gnn_tp_pdgId", &gnn_tp_pdgId);
    GNNtree->Branch("gnn_probPi", &gnn_probPi);
    GNNtree->Branch("gnn_probK", &gnn_probK);
    GNNtree->Branch("gnn_probP", &gnn_probP);
    GNNtree->Branch("gnn_sigma_tof_Pi", &gnn_sigma_tof_Pi);
    GNNtree->Branch("gnn_sigma_tof_K", &gnn_sigma_tof_K);
    GNNtree->Branch("gnn_sigma_tof_P", &gnn_sigma_tof_P);
    GNNtree->Branch("gnn_trk_chi2", &gnn_trk_chi2);
    GNNtree->Branch("gnn_trk_ndof", &gnn_trk_ndof);
    GNNtree->Branch("gnn_trk_validhits", &gnn_trk_validhits);
    GNNtree->Branch("gnn_sim_vertex_evID", &gnn_sim_vertex_evID);
    GNNtree->Branch("gnn_sim_vertex_BX", &gnn_sim_vertex_BX);
    GNNtree->Branch("gnn_sim_vertex_index", &gnn_sim_vertex_index);
    GNNtree->Branch("gnn_sim_vertex_z", &gnn_sim_vertex_z);
    GNNtree->Branch("gnn_sim_vertex_t", &gnn_sim_vertex_t);
    GNNtree->Branch("gnn_sim_vertex_isLV", &gnn_sim_vertex_isLV);

    gnn_pt.clear();
    gnn_eta.clear();
    gnn_phi.clear();
    gnn_z_pca.clear();
    gnn_dz.clear();
    gnn_t_Pi.clear();
    gnn_t_K.clear();
    gnn_t_P.clear();
    gnn_sigma_t0safe.clear();
    gnn_sigma_tmtd.clear();
    gnn_t0safe.clear();
    gnn_t0pid.clear();
    gnn_mva_qual.clear();
    gnn_btlMatchChi2.clear();
    gnn_btlMatchTimeChi2.clear();
    gnn_etlMatchChi2.clear();
    gnn_etlMatchTimeChi2.clear();
    gnn_pathLength.clear();
    gnn_npixBarrel.clear();
    gnn_npixEndcap.clear();
    gnn_outermostHitPosition.clear();
    gnn_mtdTime.clear();
    gnn_is_matched_tp.clear();
    gnn_tp_tEst.clear();
    gnn_tp_pdgId.clear();
    gnn_probPi.clear();
    gnn_probK.clear();
    gnn_probP.clear();
    gnn_sigma_tof_Pi.clear();
    gnn_sigma_tof_K.clear();
    gnn_sigma_tof_P.clear();
    gnn_trk_chi2.clear();
    gnn_trk_ndof.clear();
    gnn_trk_validhits.clear();
    gnn_sim_vertex_evID.clear();
    gnn_sim_vertex_BX.clear();
    gnn_sim_vertex_index.clear();
    gnn_sim_vertex_z.clear();
    gnn_sim_vertex_t.clear();
    gnn_sim_vertex_isLV.clear();

    // build TransientTracks
    t_tks = (*theB).build(tracksH, beamSpot, t0Safe, sigmat0Safe);

    // track filter
    std::vector<reco::TransientTrack>&& seltks = theTrackFilter->select(t_tks);

    for (std::vector<reco::TransientTrack>::const_iterator itk = seltks.begin(); itk != seltks.end(); itk++) {
      reco::TrackBaseRef trackref = (*itk).trackBaseRef();

      gnn_pt.push_back((*itk).track().pt());
      gnn_eta.push_back((*itk).track().eta());
      gnn_phi.push_back((*itk).track().phi());
      gnn_z_pca.push_back((*itk).track().vz());
      gnn_dz.push_back((*itk).track().dzError());
      gnn_t_Pi.push_back(tMtd[trackref] - tofPi[trackref]);
      gnn_t_K.push_back(tMtd[trackref] - tofK[trackref]);
      gnn_t_P.push_back(tMtd[trackref] - tofP[trackref]);
      gnn_sigma_t0safe.push_back(sigmat0Safe[trackref]);
      gnn_sigma_tmtd.push_back(sigmatimemtd[trackref]);
      gnn_t0safe.push_back(t0Safe[trackref]);
      gnn_t0pid.push_back(t0Pid[trackref]);
      gnn_mva_qual.push_back(mtdQualMVA[trackref]);
      gnn_btlMatchChi2.push_back(btlMatchChi2[trackref]);
      gnn_btlMatchTimeChi2.push_back(btlMatchTimeChi2[trackref]);
      gnn_etlMatchChi2.push_back(etlMatchChi2[trackref]);
      gnn_etlMatchTimeChi2.push_back(etlMatchTimeChi2[trackref]);
      gnn_pathLength.push_back(pathLength[trackref]);
      gnn_npixBarrel.push_back(npixBarrel[trackref]);
      gnn_npixEndcap.push_back(npixEndcap[trackref]);
      gnn_outermostHitPosition.push_back(outermostHitPosition[trackref]);
      gnn_mtdTime.push_back(tMtd[trackref]);
      gnn_probPi.push_back(probPi[trackref]);
      gnn_probK.push_back(probK[trackref]);
      gnn_probP.push_back(probP[trackref]);
      gnn_sigma_tof_Pi.push_back(sigmatofpi[trackref]);
      gnn_sigma_tof_K.push_back(sigmatofk[trackref]);
      gnn_sigma_tof_P.push_back(sigmatofp[trackref]);
      gnn_trk_chi2.push_back((*itk).track().chi2());
      gnn_trk_ndof.push_back((*itk).track().ndof());
      gnn_trk_validhits.push_back((*itk).track().numberOfValidHits());

      auto anytp_info = getAnyMatchedTP(trackref);
      if (anytp_info != nullptr) {
        gnn_is_matched_tp.push_back(true);
        double anytp_mass = (*anytp_info)->mass();
        gnn_tp_tEst.push_back(timeFromTrueMass(anytp_mass, pathLength[trackref], momentum[trackref], tMtd[trackref]));
        gnn_tp_pdgId.push_back(std::abs((*anytp_info)->pdgId()));

        TrackingVertexRef parentVertexRef = (*anytp_info)->parentVertex();

        // Loop on TV Collection to retrive info on sim vertices
        bool vertex_match = false;
        for (const auto& vsim : simpv) {
          if (vsim.sim_vertex == parentVertexRef) {
            vertex_match = true;
            // Found the matching simPrimaryVertex
            gnn_sim_vertex_z.push_back(vsim.z);
            gnn_sim_vertex_t.push_back(vsim.t * simUnit_);
            gnn_sim_vertex_evID.push_back(vsim.eventId);
            gnn_sim_vertex_BX.push_back(vsim.bunchCrossing);
            gnn_sim_vertex_index.push_back(vsim.key);
            gnn_sim_vertex_isLV.push_back(vsim.is_LV);
          }
        }
        if (vertex_match == false) {
          gnn_sim_vertex_z.push_back(-999.);
          gnn_sim_vertex_t.push_back(-999.);
          gnn_sim_vertex_evID.push_back(-999);
          gnn_sim_vertex_BX.push_back(-999);
          gnn_sim_vertex_index.push_back(-999);
          gnn_sim_vertex_isLV.push_back(false);
        }

      } else {
        gnn_is_matched_tp.push_back(false);
        gnn_tp_tEst.push_back(-999.);
        gnn_tp_pdgId.push_back(-999);
        gnn_sim_vertex_z.push_back(-999.);
        gnn_sim_vertex_t.push_back(-999.);
        gnn_sim_vertex_evID.push_back(-999);
        gnn_sim_vertex_BX.push_back(-999);
        gnn_sim_vertex_index.push_back(-999);
        gnn_sim_vertex_isLV.push_back(false);
      }

    }  // loop on sel tracks

    GNNtree->Fill();

  }  // ntuple for GNN

  // Fill TTree with input variables for BDT
  if (saveNtupleforBDT_) {
    std::string BDTtreeName = "BDTtree_" + std::to_string(iEvent.id().event());
    TTree* BDTtree = fs_->make<TTree>(BDTtreeName.c_str(), "Tree for BDT tracks");

    BDTtree->Branch("Track_pt", &Ttrack_pt);
    BDTtree->Branch("Track_eta", &Ttrack_eta);
    BDTtree->Branch("Track_phi", &Ttrack_phi);
    BDTtree->Branch("Track_dz", &Ttrack_dz);
    BDTtree->Branch("Track_dxy", &Ttrack_dxy);
    BDTtree->Branch("Track_chi2", &Ttrack_chi2);
    BDTtree->Branch("Track_ndof", &Ttrack_ndof);
    BDTtree->Branch("Track_nValidHits", &Ttrack_nValidHits);
    BDTtree->Branch("Track_npixBarrelValidHits", &Ttrack_npixBarrelValidHits);
    BDTtree->Branch("Track_npixEndcapValidHits", &Ttrack_npixEndcapValidHits);
    BDTtree->Branch("Track_Signal", &Ttrack_Signal);
    BDTtree->Branch("Track_Associated", &Ttrack_Associated);
    BDTtree->Branch("Track_BTLchi2", &Ttrack_BTLchi2);
    BDTtree->Branch("Track_BTLtime_chi2", &Ttrack_BTLtime_chi2);
    BDTtree->Branch("Track_ETLchi2", &Ttrack_ETLchi2);
    BDTtree->Branch("Track_ETLtime_chi2", &Ttrack_ETLtime_chi2);
    BDTtree->Branch("Track_t0", &Ttrack_t0);
    BDTtree->Branch("Track_sigmat0", &Ttrack_sigmat0);
    BDTtree->Branch("Track_Tmtd", &Ttrack_Tmtd);
    BDTtree->Branch("Track_MtdMVA", &Ttrack_MtdMVA);
    BDTtree->Branch("Track_lHitPos", &Ttrack_lHitPos);
    BDTtree->Branch("Track_sigmaTmtd", &Ttrack_sigmaTmtd);
    BDTtree->Branch("Track_length", &Ttrack_length);
    BDTtree->Branch("Track_DirClu", &Ttrack_TPDirClu);
    BDTtree->Branch("Track_OtherClu", &Ttrack_TPOtherClu);
    BDTtree->Branch("Track_isBTL", &Ttrack_isBTL);
    BDTtree->Branch("Track_isETL", &Ttrack_isETL);
    BDTtree->Branch("Track_withMTD", &Ttrack_withMTD);

    Ttrack_pt.clear();
    Ttrack_eta.clear();
    Ttrack_phi.clear();
    Ttrack_dz.clear();
    Ttrack_dxy.clear();
    Ttrack_chi2.clear();
    Ttrack_ndof.clear();
    Ttrack_nValidHits.clear();
    Ttrack_npixBarrelValidHits.clear();
    Ttrack_npixEndcapValidHits.clear();
    Ttrack_Signal.clear();
    Ttrack_Associated.clear();
    Ttrack_BTLchi2.clear();
    Ttrack_BTLtime_chi2.clear();
    Ttrack_ETLchi2.clear();
    Ttrack_ETLtime_chi2.clear();
    Ttrack_t0.clear();
    Ttrack_sigmat0.clear();
    Ttrack_Tmtd.clear();
    Ttrack_MtdMVA.clear();
    Ttrack_lHitPos.clear();
    Ttrack_sigmaTmtd.clear();
    Ttrack_length.clear();
    Ttrack_TPDirClu.clear();
    Ttrack_TPOtherClu.clear();
    Ttrack_isBTL.clear();
    Ttrack_isETL.clear();
    Ttrack_withMTD.clear();

    unsigned int index = 0;

    // --- Loop over all RECO tracks ---
    for (const auto& trackGen : *RecTrackHandle) {
      const reco::TrackRef trackref(iEvent.getHandle(RecTrackToken_), index);
      index++;

      if (trackAssoc[trackref] == -1) {
        LogInfo("mtdTracks") << "Extended track not associated";
        continue;
      }

      const reco::TrackRef mtdTrackref = reco::TrackRef(iEvent.getHandle(RecMTDTrackToken_), trackAssoc[trackref]);
      const reco::Track& track = *mtdTrackref;

      bool withMTD = false;
      bool isBTL = false;
      bool isETL = false;
      bool twoETLdiscs = false;
      bool isTPDirClu = false;
      bool isTPOtherClu = false;

      bool isTPmtdDirectBTL = false, isTPmtdOtherBTL = false, isTPmtdDirectETL = false, isTPmtdOtherETL = false,
           isTPmtdCorrect = false;

      if (trkRecSel(trackGen)) {
        if (std::round(sigmatimemtd[trackref] - sigmat0Pid[trackref]) != 0) {
          LogWarning("mtdTracks")
              << "TimeError associated to refitted track is different from TimeError stored in tofPID "
                 "sigmat0 ValueMap: this should not happen";
        }

        std::vector<edm::Ref<edmNew::DetSetVector<FTLCluster>, FTLCluster>> recoClustersRefs;

        if (std::abs(trackGen.eta()) < etacutBTL_) {
          // --- all BTL tracks (with and without hit in MTD) ---

          for (const auto hit : track.recHits()) {
            if (hit->isValid() == false)
              continue;
            MTDDetId Hit = hit->geographicalId();
            if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 1)) {
              isBTL = true;
              const auto* mtdhit = static_cast<const MTDTrackingRecHit*>(hit);
              const auto& hitCluster = mtdhit->mtdCluster();
              if (hitCluster.size() != 0) {
                auto recoClusterRef = edmNew::makeRefTo(btlRecCluHandle, &hitCluster);
                recoClustersRefs.push_back(recoClusterRef);
              }
            }
          }
        }  //loop over (geometrical) BTL tracks
        else {
          // --- all ETL tracks (with and without hit in MTD) ---
          for (const auto hit : track.recHits()) {
            if (hit->isValid() == false)
              continue;
            MTDDetId Hit = hit->geographicalId();
            if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 2)) {
              isETL = true;

              const auto* mtdhit = static_cast<const MTDTrackingRecHit*>(hit);
              const auto& hitCluster = mtdhit->mtdCluster();
              if (hitCluster.size() != 0) {
                auto recoClusterRef = edmNew::makeRefTo(etlRecCluHandle, &hitCluster);
                recoClustersRefs.push_back(recoClusterRef);
              }
            }
          }
        }

        LogDebug("MVATrainingNtuple") << "Track p/pt = " << trackGen.p() << " " << trackGen.pt() << " eta "
                                      << trackGen.eta() << " BTL " << isBTL << " ETL " << isETL << " 2disks "
                                      << twoETLdiscs;

        // == TrackingParticle based matching
        const reco::TrackBaseRef trkrefb(trackref);
        auto tp_info = getAnyMatchedTP(trkrefb);
        if (tp_info != nullptr && trkTPSelAll(**tp_info)) {
          // ==  MC truth matching

          auto simClustersRefsIt = tp2SimAssociationMap.find(*tp_info);
          withMTD = (simClustersRefsIt != tp2SimAssociationMap.end());

          // If there is a mtdSimLayerCluster from the tracking particle
          if (withMTD) {
            // -- Get the refs to MtdSimLayerClusters associated to the TP
            std::vector<edm::Ref<MtdSimLayerClusterCollection>> simClustersRefs;
            for (const auto& ref : simClustersRefsIt->val) {
              simClustersRefs.push_back(ref);
            }
            // === ETL
            // -- Sort ETL sim clusters by time
            if (std::abs(trackGen.eta()) >= etacutBTL_ && !simClustersRefs.empty()) {
              std::sort(simClustersRefs.begin(), simClustersRefs.end(), [](const auto& a, const auto& b) {
                return a->simLCTime() < b->simLCTime();
              });
              // Check if TP has direct or other sim cluster for BTL
              for (const auto& simClusterRef : simClustersRefs) {
                if (simClusterRef->trackIdOffset() == 0) {
                  isTPmtdDirectETL = true;
                } else if (simClusterRef->trackIdOffset() != 0) {
                  isTPmtdOtherETL = true;
                }
              }
            }
            // === BTL
            // -- Sort BTL sim clusters by time
            if (std::abs(trackGen.eta()) < etacutBTL_ && !simClustersRefs.empty()) {
              std::sort(simClustersRefs.begin(), simClustersRefs.end(), [](const auto& a, const auto& b) {
                return a->simLCTime() < b->simLCTime();
              });
              // Check if TP has direct or other sim cluster for BTL
              for (const auto& simClusterRef : simClustersRefs) {
                if (simClusterRef->trackIdOffset() == 0) {
                  isTPmtdDirectBTL = true;
                } else if (simClusterRef->trackIdOffset() != 0) {
                  isTPmtdOtherBTL = true;
                }
              }
            }

            // ==  Check if the track-cluster association is correct: Track->RecoClus->SimClus == Track->TP->SimClus
            for (const auto& recClusterRef : recoClustersRefs) {
              if (recClusterRef.isNonnull()) {
                auto itp = r2sAssociationMap.equal_range(recClusterRef);
                if (itp.first != itp.second) {
                  auto& simClustersRefs_RecoMatch = (*itp.first).second;
                  for (const auto& simClusterRef_RecoMatch : simClustersRefs_RecoMatch) {
                    // Check if simClusterRef_RecoMatch  exists in SimClusters
                    auto simClusterIt =
                        std::find(simClustersRefs.begin(), simClustersRefs.end(), simClusterRef_RecoMatch);
                    // SimCluster found in SimClusters
                    if (simClusterIt != simClustersRefs.end()) {
                      isTPmtdCorrect = true;
                    }
                  }
                }
              }
            }  /// end loop over reco clusters associated to this track.
          }  // --- end "withMTD"
        }  // TP matching
        Ttrack_pt.push_back(trackGen.pt());
        Ttrack_phi.push_back(trackGen.phi());
        Ttrack_eta.push_back(trackGen.eta());
        Ttrack_dz.push_back(std::abs(trackGen.dz(beamSpot.position())));
        Ttrack_dxy.push_back(std::abs(trackGen.dxy(beamSpot.position())));
        Ttrack_chi2.push_back(trackGen.chi2());
        Ttrack_ndof.push_back(trackGen.ndof());
        Ttrack_nValidHits.push_back(trackGen.numberOfValidHits());

        Ttrack_npixBarrelValidHits.push_back(npixBarrel[trackref]);
        Ttrack_npixEndcapValidHits.push_back(npixEndcap[trackref]);
        Ttrack_BTLchi2.push_back(btlMatchChi2[trackref]);
        Ttrack_BTLtime_chi2.push_back(btlMatchTimeChi2[trackref]);
        Ttrack_ETLchi2.push_back(etlMatchChi2[trackref]);
        Ttrack_ETLtime_chi2.push_back(etlMatchTimeChi2[trackref]);

        Ttrack_t0.push_back(t0Src[trackref]);
        Ttrack_sigmat0.push_back(sigmat0Src[trackref]);
        Ttrack_Tmtd.push_back(tMtd[trackref]);
        Ttrack_sigmaTmtd.push_back(sigmatimemtd[trackref]);
        Ttrack_length.push_back(pathLength[trackref]);
        Ttrack_MtdMVA.push_back(mtdQualMVA[trackref]);
        Ttrack_lHitPos.push_back(outermostHitPosition[trackref]);

        Ttrack_isBTL.push_back(isBTL);
        Ttrack_isETL.push_back(isETL);
        Ttrack_withMTD.push_back(withMTD);

        if ((tp_info != nullptr) && (*tp_info)->eventId().bunchCrossing() == 0 &&
            (*tp_info)->eventId().event() == 0) {  // Signal vs PU seperation
          Ttrack_Signal.push_back(true);           // Signal track
        } else {
          Ttrack_Signal.push_back(false);  // PU track?
        }

        if (isTPmtdCorrect) {
          Ttrack_Associated.push_back(true);  // Associated track
        } else {
          Ttrack_Associated.push_back(false);  // Not associated track
        }
        if (isTPmtdDirectBTL || isTPmtdDirectETL) {
          isTPDirClu = true;  // Direct cluster
        }
        if (isTPmtdOtherBTL || isTPmtdOtherETL) {
          isTPOtherClu = true;  // Not direct cluster
        }

        Ttrack_TPDirClu.push_back(isTPDirClu);
        Ttrack_TPOtherClu.push_back(isTPOtherClu);
      }  // trkRecSel

    }  // RECO tracks loop

    BDTtree->Fill();

  }  // ntuple for BDT
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MVATrainingNtuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("inputTracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("inputTagT", edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>("inputTagV", edm::InputTag("offlinePrimaryVertices4D"));
  desc.add<edm::InputTag>("TPtoRecoTrackAssoc", edm::InputTag("trackingParticleRecoTrackAsssociation"));
  desc.add<edm::InputTag>("tp2SimAssociationMapTag", edm::InputTag("mtdSimLayerClusterToTPAssociation"));
  desc.add<edm::InputTag>("r2sAssociationMapTag", edm::InputTag("mtdRecoClusterToSimLayerClusterAssociation"));
  desc.add<edm::InputTag>("SimTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("offlineBS", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"))
      ->setComment("Association between General and MTD Extended tracks");
  desc.add<edm::InputTag>("recCluTagBTL", edm::InputTag("mtdClusters", "FTLBarrel"));
  desc.add<edm::InputTag>("recCluTagETL", edm::InputTag("mtdClusters", "FTLEndcap"));
  desc.add<edm::InputTag>("pathLengthSrc", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("momentumSrc", edm::InputTag("trackExtenderWithMTD:generalTrackp"));
  desc.add<edm::InputTag>("tmtd", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmaSrc", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("t0Src", edm::InputTag("trackExtenderWithMTD:generalTrackt0"));
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("trackExtenderWithMTD:generalTracksigmat0"));
  desc.add<edm::InputTag>("t0PID", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("sigmat0PID", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("t0SafePID", edm::InputTag("tofPID:t0safe"));
  desc.add<edm::InputTag>("sigmat0SafePID", edm::InputTag("tofPID:sigmat0safe"));
  desc.add<edm::InputTag>("trackMVAQual", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<edm::InputTag>("tofPi", edm::InputTag("trackExtenderWithMTD:generalTrackTofPi"));
  desc.add<edm::InputTag>("tofK", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"));
  desc.add<edm::InputTag>("tofP", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"));
  desc.add<edm::InputTag>("probPi", edm::InputTag("tofPID:probPi"));
  desc.add<edm::InputTag>("probK", edm::InputTag("tofPID:probK"));
  desc.add<edm::InputTag>("probP", edm::InputTag("tofPID:probP"));
  desc.add<edm::InputTag>("sigmatofpiSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"));
  desc.add<edm::InputTag>("sigmatofkSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"));
  desc.add<edm::InputTag>("sigmatofpSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"));
  desc.add<edm::InputTag>("btlMatchChi2Src", edm::InputTag("trackExtenderWithMTD", "btlMatchChi2"));
  desc.add<edm::InputTag>("btlMatchTimeChi2Src", edm::InputTag("trackExtenderWithMTD", "btlMatchTimeChi2"));
  desc.add<edm::InputTag>("etlMatchChi2Src", edm::InputTag("trackExtenderWithMTD", "etlMatchChi2"));
  desc.add<edm::InputTag>("etlMatchTimeChi2Src", edm::InputTag("trackExtenderWithMTD", "etlMatchTimeChi2"));
  desc.add<edm::InputTag>("npixBarrelSrc", edm::InputTag("trackExtenderWithMTD", "npixBarrel"));
  desc.add<edm::InputTag>("npixEndcapSrc", edm::InputTag("trackExtenderWithMTD", "npixEndcap"));
  desc.add<edm::InputTag>("outermostHitPositionSrc",
                          edm::InputTag("trackExtenderWithMTD", "generalTrackOutermostHitPosition"));
  desc.addUntracked<std::string>("fileName", "file.root");
  desc.add<bool>("ntupleforBDT", true);
  desc.add<bool>("ntupleforGNN", false);
  {
    edm::ParameterSetDescription psd0;
    HITrackFilterForPVFinding::fillPSetDescription(psd0);  // extension of TrackFilterForPVFinding
    desc.add<edm::ParameterSetDescription>("TkFilterParameters", psd0);
  }

  descriptions.add("mvaTrainingNtuple", desc);
}

const bool MVATrainingNtuple::trkTPSelAll(const TrackingParticle& tp) {
  bool match = false;

  auto x_pv = tp.parentVertex()->position().x();
  auto y_pv = tp.parentVertex()->position().y();
  auto z_pv = tp.parentVertex()->position().z();

  auto r_pv = std::sqrt(x_pv * x_pv + y_pv * y_pv);

  match = tp.charge() != 0 && std::abs(tp.eta()) < etacutGEN_ &&
          ((std::abs(tp.eta()) < etacutBTL_ && tp.pt() > pTcutBTL_) ||
           (std::abs(tp.eta()) > etacutBTL_ && tp.pt() > pTcutETL_)) &&
          r_pv < rBTL_ && std::abs(z_pv) < zETL_;
  return match;
}

const bool MVATrainingNtuple::trkRecSel(const reco::TrackBase& trk) {
  bool match = false;
  match = std::abs(trk.eta()) <= etacutREC_ && ((std::abs(trk.eta()) <= etacutBTL_ && trk.pt() > pTcutBTL_) ||
                                                (std::abs(trk.eta()) > etacutBTL_ && trk.pt() > pTcutETL_));
  return match;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MVATrainingNtuple);
