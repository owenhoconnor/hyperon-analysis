////////////////////////////////////////////////////////////////////////
// Class:       hyperon:AnalyzeEvents
// Plugin Type: analyzer (Unknown Unknown)
// File:        hyperon:AnalyzeEvents_module.cc
//
// Generated at Fri Oct  4 15:30:11 2024 by Jarek Nowak using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////



#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"


#include "art_root_io/TFileService.h"
// LArSoft Includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

// // Root Includes
#include <iostream>
#include <vector>
#include <TTree.h>
#include <TH1.h>
#include <string>
#include <exception>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>


// Anonymous namespace for helper functions and structs

namespace {

// Helper functions for MCTruth Saving

bool IsDecayProcess(std::string process)
{
    std::transform(
        process.begin(),
        process.end(),
        process.begin(),
        [](unsigned char c) {
            return std::tolower(c);
        }
    );

    return process.find("decay") != std::string::npos;
}

// Structs for Track/Shower Reco Info
struct TrackRecoInfo
{
    bool valid = false;
    bool uniqueTrack = false;
    int id = -1;

    float length = -9999.f;

    float startX = -9999.f;
    float startY = -9999.f;
    float startZ = -9999.f;

    float endX = -9999.f;
    float endY = -9999.f;
    float endZ = -9999.f;

    float startDirX = -9999.f;
    float startDirY = -9999.f;
    float startDirZ = -9999.f;

    float endDirX = -9999.f;
    float endDirY = -9999.f;
    float endDirZ = -9999.f;

    float theta = -9999.f;
    float phi   = -9999.f;
};

struct ShowerRecoInfo
{
    bool valid = false;
    bool uniqueShower = false;
    int id = -1;

    float length = -9999.f;

    float startX = -9999.f;
    float startY = -9999.f;
    float startZ = -9999.f;

    float dirX = -9999.f;
    float dirY = -9999.f;
    float dirZ = -9999.f;
};


TrackRecoInfo GetTrackRecoInfo(const std::vector<art::Ptr<recob::Track>>& tracks)
{
    TrackRecoInfo info;

    if (tracks.empty()) {
        return info;
    }
    if(tracks.size() == 1){
        info.uniqueTrack = true;
    }

    const art::Ptr<recob::Track>& track = tracks.front();

    info.valid = true;
    info.id = track->ID();
    info.length = track->Length();

    const auto& start = track->Vertex();
    const auto& end   = track->End();

    info.startX = start.X();
    info.startY = start.Y();
    info.startZ = start.Z();

    info.endX = end.X();
    info.endY = end.Y();
    info.endZ = end.Z();

    const auto& startDir = track->StartDirection();
    const auto& endDir   = track->EndDirection();

    info.startDirX = startDir.X();
    info.startDirY = startDir.Y();
    info.startDirZ = startDir.Z();

    info.endDirX = endDir.X();
    info.endDirY = endDir.Y();
    info.endDirZ = endDir.Z();

    info.theta = track->Theta();
    info.phi   = track->Phi();

    return info;
}

ShowerRecoInfo GetShowerRecoInfo(const std::vector<art::Ptr<recob::Shower>>& showers)
{
    ShowerRecoInfo info;

    if (showers.empty()) {
        return info;
    }
    if (showers.size() == 1){
        info.uniqueShower = true;
    }

    const art::Ptr<recob::Shower>& shower = showers.front();

    info.valid = true;
    info.id = shower->ID();
    info.length = shower->Length();

    const auto& start = shower->ShowerStart();

    info.startX = start.X();
    info.startY = start.Y();
    info.startZ = start.Z();

    const auto& dir = shower->Direction();

    info.dirX = dir.X();
    info.dirY = dir.Y();
    info.dirZ = dir.Z();

    return info;
}


float GetPFPTrackScore(
    const art::Ptr<recob::PFParticle>& pfp,
    const art::FindManyP<larpandoraobj::PFParticleMetadata>& metadataAssoc)
{
    const auto metadataVec = metadataAssoc.at(pfp.key());

    for (const auto& metadata : metadataVec)
    {
        const auto& properties = metadata->GetPropertiesMap();

        auto it = properties.find("TrackScore");

        if (it != properties.end()) {
            return it->second;
        }
    }

    return -9999.f;
}

// Struct for PFP truth matching

struct PFPTruthMatch
{
    bool valid = false;
    int trackID = -9999;
    int pdg = -9999;

    int nHits = 0;
    int nMatchedHits = 0;
    float purity = -1.f;
};

} // anonymous namespace


namespace hyperon {
  class AnalyzeEvents;
}

class hyperon::AnalyzeEvents : public art::EDAnalyzer {
public:
  explicit AnalyzeEvents(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeEvents(AnalyzeEvents const&) = delete;
  AnalyzeEvents(AnalyzeEvents&&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents const&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void endSubRun(art::SubRun const& sr) override;

private:

  // Declare member data here.

  int fVerbose;

  std::string fHitLabel;
  std::string fGenieGenModuleLabel;
  std::string fPFParticleLabel;
  std::string fSliceLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fClusterLabel;

  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

// Neutrino Interaction variables
   int intType, CCNC, neutrinoPDG, numProtons, numNeutrons, numPi, numPi0, numTrueHits;
   float W, X, Y, QSqr, Pt, Theta, neutrinoE, leptonP;
   float trueVertexX, trueVertexY, trueVertexZ;


   int lambda=0,mu=0, mubar=0,kaonp=0,kaonm=0,kaon0=0,proton=0,neutron=0,pip=0,pim=0,pi0=0;
   int sigmap=0, sigma0=0, sigmam=0, gamma=0, goodSigma=0, goodLambda=0;
   int nMultiEvents = 0;
   int nSingleEvents = 0;
   int nTotEvents = 0;
   int nMCParticles = 0;

  TTree *fTree;
  unsigned int fEventID;

  std::vector<int> fNSlices;
  size_t totalSlices = 0;
  size_t totalNeutrinos = 0;
  int nuID = -1;
  int nuSliceKey = -1;
  int fNPrimaryParticles;
  int fNPrimaryChildren;

  // MC Truth parameters
  int fRun;       // Run number
  int fSubRun; // Sub run number
 
  // MCTruth Level
  std::vector<int> trueOrigin;
  std::vector<float> trueW;
  std::vector<float> trueX;
  std::vector<float> trueY;
  std::vector<float> trueQSqr;
  std::vector<float> truePt;
  std::vector<float> trueTheta;
  std::vector<float> trueNuVtxX;
  std::vector<float> trueNuVtxY;
  std::vector<float> trueNuVtxZ;
  std::vector<int> trueNuPDG;
  std::vector<int> trueNuTrackID;
  std::vector<float> trueNuEnergy;
  std::vector<int> trueCCNC;
  std::vector<int> trueIntMode;
  std::vector<int> trueIntType;
  std::vector<int> trueTargetPDG;

  // MCParticle Level
  std::vector<int> truePDG;
  std::vector<int> trueTrackID; 
  std::vector<int> trueMotherPDG;
  std::vector<int> trueMotherTrackID;
  std::vector<int> trueGeneration;
  std::vector<int> trueIsPrimary;
  std::vector<int> trueIsDecayProduct;
  std::vector<int> trueMCTruthIndex;
  std::vector<int> trueNHits;
  std::vector<int> trueNHitsU;
  std::vector<int> trueNHitsV;
  std::vector<int> trueNHitsZ;
  std::vector<float> trueP;
  std::vector<float> truePx;
  std::vector<float> truePy;
  std::vector<float> truePz;
  std::vector<float> trueMass;
  std::vector<float> trueE;
  std::vector<float> trueStartX;
  std::vector<float> trueStartY;
  std::vector<float> trueStartZ;
  std::vector<float> trueStartT;
  std::vector<float> trueEndX;
  std::vector<float> trueEndY;
  std::vector<float> trueEndZ;
  std::vector<float> trueEndT;
  std::vector<float> trueEndPx;
  std::vector<float> trueEndPy;
  std::vector<float> trueEndPz;
  std::vector<float> trueEndE;
  std::vector<int> trueNTrajectoryPoints;
  std::vector<float> trueTrajectoryLength;
  std::vector<std::string> trueProcess;
  std::vector<std::string> trueEndProcess;
  std::vector<int> trueStatusCode;
  std::vector<int> trueNDaughters;
  std::vector<int> trueNStoredDecayDaughters;
  std::vector<int> trueNPrimaryParticles;
  std::vector<int> trueNSavedParticles;
  std::vector<int> trueParticleStartIndex;
  

  // Reconstructed parameters

  // Slice level

  std::vector<int> fSliceKey;
  std::vector<int> fSliceID;
  std::vector<float> fSliceNuScore;
  std::vector<int> fSliceTotalHits;
  std::vector<int> fSliceTrueNuHits;
  std::vector<int> fSliceTrueOrigin;
  int              fEventTotalTrueNuHits;
  std::vector<float> fSliceVtxX;
  std::vector<float> fSliceVtxY;
  std::vector<float> fSliceVtxZ;
  std::vector<float> fSliceOpt0Score;

  // PFP level

  std::vector<int> fPfpKey;
  std::vector<int> fPfpSelfID;
  std::vector<int> fPfpParentID;
  std::vector<int> fPfpRecoPDG;
  std::vector<int> fPfpSliceKey;
  std::vector<int> fPfpIsNuSlice;
  std::vector<int> fPfpIsPrimary;
  std::vector<int> fPfpNDaughters;
  std::vector<float> fPfpTrackScore;
  std::vector<int> fPfpHasTrackScore;

  std::vector<int> fPfpNTracks;
  std::vector<int> fPfpNShowers;
  std::vector<int> fPfpHasTrack;
  std::vector<int> fPfpHasUniqueTrack;
  std::vector<int> fPfpHasShower;
  std::vector<int> fPfpHasUniqueShower;

  // Reco track parameters
  std::vector<int> fPfpTrackID;
  std::vector<float> fPfpTrackLength;
  std::vector<float> fPfpTrackStartX;
  std::vector<float> fPfpTrackStartY;
  std::vector<float> fPfpTrackStartZ;
  std::vector<float> fPfpTrackEndX;
  std::vector<float> fPfpTrackEndY;
  std::vector<float> fPfpTrackEndZ;
  std::vector<float> fPfpTrackStartDirX;
  std::vector<float> fPfpTrackStartDirY;
  std::vector<float> fPfpTrackStartDirZ;
  std::vector<float> fPfpTrackEndDirX;
  std::vector<float> fPfpTrackEndDirY;
  std::vector<float> fPfpTrackEndDirZ;
  std::vector<float> fPfpTrackVertexDirX;
  std::vector<float> fPfpTrackVertexDirY;
  std::vector<float> fPfpTrackVertexDirZ;
  std::vector<float> fPfpTrackTheta;
  std::vector<float> fPfpTrackPhi;

  std::vector<int> fPfpTrackSliceID;
  std::vector<int> fPfpTrackTrueG4ID;
  std::vector<int> fPfpTrackIsPrimary;

  // Reco shower parameters
  std::vector<int> fPfpShowerID;
  std::vector<float> fPfpShowerLength;
  std::vector<float> fPfpShowerStartX;
  std::vector<float> fPfpShowerStartY;
  std::vector<float> fPfpShowerStartZ;
  std::vector<float> fPfpShowerDirX;
  std::vector<float> fPfpShowerDirY;
  std::vector<float> fPfpShowerDirZ;
 
  std::vector<int> fPfpShowerSliceID;
  std::vector<int> fPfpShowerTrueG4ID;
  std::vector<int> fPfpShowerIsPrimary;

  // Pfp vertex Params
  std::vector<int> fPfpNVertices;
  std::vector<int> fPfpHasVertex;
  std::vector<int> fPfpHasUniqueVertex;
  std::vector<float> fPfpVertexX;
  std::vector<float> fPfpVertexY;
  std::vector<float> fPfpVertexZ;

  // Truth matching parameters

  std::vector<int> fPfpTrueTrackID;
  std::vector<int> fPfpTruePDG;
  std::vector<int> fPfpNHits;
  std::vector<int> fPfpNMatchedHits;
  std::vector<float> fPfpTruthPurity;
  std::vector<int> fTrueIsReconstructed;
  std::vector<int> fTrueIsReconstructedInNuSlice;
  std::vector<int> fTrueNMatchedPfps;
  std::vector<int> fTrueBestRecoPfpIdx;
  std::vector<float> fTrueBestRecoTrackScore;
  std::vector<int> fTrueBestRecoHasTrack;
  std::vector<int> fTrueBestRecoHasShower;

  // other stuff

  std::vector<float> fnuScore;
  std::map<int, int> trueTrackPDGMap;

  std::vector<float> fNeutrinoNuScores;
  std::vector<float> fCosmicNuScores;
  
  float highestNuScore = -1;


  // Run and Subrun Information

  int fRun_sr;
  int fSubRun_sr;
  double fPOT;
  TTree *fSubRunTree;
  
};


hyperon::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset} 
  , fHitLabel(pset.get<std::string>("HitLabel"))
  , fGenieGenModuleLabel(pset.get<std::string>("GenieGenModuleLabel"))
  , fPFParticleLabel(pset.get<std::string>("PFParticleLabel"))
  , fSliceLabel(pset.get<std::string>("SliceLabel"))
  , fTrackLabel(pset.get<std::string>("TrackLabel"))
  //, fCalorimetryLabel(pset.get<std::string>("CalorimetryLabel"))
  , fShowerLabel(pset.get<std::string>("ShowerLabel"))
  , fClusterLabel(pset.get<std::string>("ClusterLabel"))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void hyperon::AnalyzeEvents::analyze(art::Event const& evt)
{
  // Implementation of required member function here.
    fEventID = evt.id().event(); 
    std::cout<<"Event# "<<evt.id().event()<<std::endl;
    fRun = evt.run();
    fSubRun = evt.subRun();

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    art::ServiceHandle<cheat::ParticleInventoryService> piService;
    art::ServiceHandle<cheat::BackTrackerService> btService;

 // Clear reco parameters

  fSliceKey.clear();
  fSliceID.clear();
  fSliceNuScore.clear();
  fSliceTotalHits.clear();
  fSliceTrueNuHits.clear();
  fSliceTrueOrigin.clear();
  fEventTotalTrueNuHits = 0;
  fSliceVtxX.clear();
  fSliceVtxY.clear();
  fSliceVtxZ.clear();
  fSliceOpt0Score.clear();

  fPfpKey.clear();
  fPfpSelfID.clear();
  fPfpParentID.clear();
  fPfpRecoPDG.clear();
  fPfpSliceKey.clear();
  fPfpIsNuSlice.clear();
  fPfpIsPrimary.clear();
   fPfpNDaughters.clear();
  fPfpTrackScore.clear();
  fPfpHasTrackScore.clear();
  fPfpNTracks.clear();
  fPfpNShowers.clear();
  fPfpHasTrack.clear();
  fPfpHasUniqueTrack.clear();
  fPfpHasShower.clear();
  fPfpHasUniqueShower.clear();
  fPfpNVertices.clear();
    fPfpHasVertex.clear();
    fPfpHasUniqueVertex.clear();
    fPfpVertexX.clear();
    fPfpVertexY.clear();
    fPfpVertexZ.clear();
  fPfpTrueTrackID.clear();
  fPfpTruePDG.clear();
  fPfpNHits.clear();
  fPfpNMatchedHits.clear();
  fPfpTruthPurity.clear();
  fTrueIsReconstructed.clear();
  fTrueIsReconstructedInNuSlice.clear();
  fTrueNMatchedPfps.clear(); 
  fTrueBestRecoPfpIdx.clear();
  fTrueBestRecoTrackScore.clear();
  fTrueBestRecoHasTrack.clear();
  fTrueBestRecoHasShower.clear();
 fPfpTrackID.clear();
 fPfpTrackLength.clear();
 fPfpTrackStartX.clear();
 fPfpTrackStartY.clear();
 fPfpTrackStartZ.clear();
 fPfpTrackEndX.clear();
 fPfpTrackEndY.clear();
 fPfpTrackEndZ.clear();
 fPfpTrackStartDirX.clear();
 fPfpTrackStartDirY.clear();
 fPfpTrackStartDirZ.clear();
 fPfpTrackEndDirX.clear();
 fPfpTrackEndDirY.clear();
 fPfpTrackEndDirZ.clear();
 fPfpTrackVertexDirX.clear();
 fPfpTrackVertexDirY.clear();
 fPfpTrackVertexDirZ.clear();
 fPfpTrackTheta.clear();
 fPfpTrackPhi.clear();
 fPfpShowerID.clear();
 fPfpShowerLength.clear();
 fPfpShowerStartX.clear();
 fPfpShowerStartY.clear();
 fPfpShowerStartZ.clear();
 fPfpShowerDirX.clear();
 fPfpShowerDirY.clear();
 fPfpShowerDirZ.clear();

 fnuScore.clear();
 fNeutrinoNuScores.clear();
 fCosmicNuScores.clear();
 highestNuScore = -1.f;
 nuSliceKey = -1;
 nuID = -1;


 fNPrimaryParticles = 0;



 // =====================================================
 // Reconstructed Slices Analysis
 // =====================================================

 // Get event slices
 
   art::ValidHandle<std::vector<recob::Slice>> sliceHandle = evt.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
   art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
   std::vector<art::Ptr<recob::Slice>> sliceVector;
   if (sliceHandle.isValid()){art::fill_ptr_vector(sliceVector, sliceHandle);}

// Get associations between slices and PFParticles, Get associations between PFParticles and Clusters, Clusters and Hits, and between Hits and MCParticles, // Get associations between PFParticle and Vertex
// Get associations between PFPParticles and PFPMetaData

   art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, evt, fSliceLabel);
   art::FindManyP<recob::Hit> sliceHitAssoc(sliceHandle, evt, fSliceLabel);
   art::FindManyP<recob::Vertex> pfpVertexAssoc(pfpHandle, evt, fPFParticleLabel);
   art::FindManyP<larpandoraobj::PFParticleMetadata> pfpMetadataAssoc(pfpHandle, evt, fPFParticleLabel);
   art::FindManyP<recob::Cluster> pfpClusterAssoc(pfpHandle, evt, fPFParticleLabel);
   art::FindManyP<recob::Hit> clusterHitAssoc(evt.getValidHandle<std::vector<recob::Cluster>>(fClusterLabel), evt, fClusterLabel);

   art::Handle<std::vector<recob::Hit>> globalHitHandle;
   evt.getByLabel(fHitLabel, globalHitHandle);
   art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> hitTruthAssns(globalHitHandle, evt, "gaushitTruthMatch");
   
   std::vector<art::Ptr<recob::Hit>> globalHits;
   art::fill_ptr_vector(globalHits, globalHitHandle);

   //std::unique_ptr<art::FindManyP<sbn::OpT0Finder>> opt0Assns;
   //try { opt0Assns = std::make_unique<art::FindManyP<sbn::OpT0Finder>>(sliceHandle, e, fOpT0Label); } catch(...) {}


      // art::FindManyP<simb::MCParticle> hitMCParticleAssoc(evt.getValidHandle<std::vector<recob::Hit>>(fHitLabel), evt, fHitLabel);

   // Define helper function to get hits from PFP
   auto getPFPHits = [&](const art::Ptr<recob::PFParticle>& pfp){
        // Get clusters associated with PFP and loop over them
        std::vector<art::Ptr<recob::Hit>> pfpHits;
        const std::vector<art::Ptr<recob::Cluster>> clusters = pfpClusterAssoc.at(pfp.key());

        for (const art::Ptr<recob::Cluster>& cluster : clusters) {

            // get hits associated with cluster
            const std::vector<art::Ptr<recob::Hit>> clusterHits = clusterHitAssoc.at(cluster.key());
            pfpHits.insert(pfpHits.end(), clusterHits.begin(), clusterHits.end());
        }

        return pfpHits;
    };

    // Define helper to define whether a hit is truth matched to a beam neutrino
    /*auto hitIsBeamNeutrino = [&](const art::Ptr<recob::Hit>& hit){
        if (!hitTruthAssns.isValid()) {
            return false;
        }

        const auto& particles = hitTruthAssns.at(hit.key());

        for (const auto& truePart : particles)
        {
            const art::Ptr<simb::MCTruth> hitMCTruth = piService->TrackIdToMCTruth_P(std::abs(truePart->TrackId()));

            if (hitMCTruth.isNonnull() && hitMCTruth->Origin() == simb::kBeamNeutrino){
                return true;
            }
        }

        return false;
    };
    
*/
    constexpr bool rollupUnsavedIDs = true;

    struct HitTruthInfo
    {
        bool valid = false;
        int trackID = -9999;
        int origin = -1;
        bool isBeamNeutrino = false;
        bool isCosmic = false;
    };

    auto getHitTruthInfo = [&](const art::Ptr<recob::Hit>& hit)
    {
        HitTruthInfo info;
        const TruthMatchUtils::G4ID g4ID = TruthMatchUtils::TrueParticleID(clockData, hit, rollupUnsavedIDs);

        if (!TruthMatchUtils::Valid(g4ID)) {
            return info;
        }

        const art::Ptr<simb::MCTruth> hitMCTruth = piService->TrackIdToMCTruth_P(g4ID);

        if (!hitMCTruth.isNonnull()) {
                return info;
        }

        info.valid = true;
        info.trackID = g4ID;
        info.origin = static_cast<int>(hitMCTruth->Origin());
        info.isBeamNeutrino = (hitMCTruth->Origin() == simb::kBeamNeutrino);
        info.isCosmic = (hitMCTruth->Origin() == simb::kCosmicRay);

        return info;
    };

    // Loop over global hits

    struct TrueHitCounts{
        int total = 0;
        int U = 0;
        int V = 0;
        int Z = 0;
    };
    std::unordered_map<int, TrueHitCounts> trueHitCounts;

    for (const art::Ptr<recob::Hit>& hit : globalHits){
        //bool isBeamNeutrino = hitIsBeamNeutrino(iHit);

        const HitTruthInfo truthInfo = getHitTruthInfo(hit);
        if(!truthInfo.valid){continue;}

        if(truthInfo.isBeamNeutrino){
            ++fEventTotalTrueNuHits;
        }

        TrueHitCounts& counts = trueHitCounts[truthInfo.trackID];
        ++counts.total;

        switch (hit->View()){
            case geo::kU:
                ++counts.U;
                break;
            
            case geo::kV:
                ++counts.V;
                break;
            
            case geo::kZ:
                ++counts.Z;
                break;
            
            default:
                std::cerr<<"Unexpected hit view: "<<static_cast<int>(hit->View())<<std::endl;
                break;
        }

        std::cout
        << "TPC " << hit->WireID().TPC
        << " plane " << hit->WireID().Plane
        << " view " << static_cast<int>(hit->View())
        << std::endl;

    }

// Filling our neutrino hierarchy variables by looping over slices in event:
// What are the properties of the slices, the nuScores, the truth origin, and what is the nuSlice?

   std::cout<<"Event "<<fEventID<<" has "<<sliceVector.size()<<" slices."<<std::endl;
   if (sliceVector.size() == 0){
	   std::cerr<<"No slices found in this event!"<<std::endl;
   }

   for (const art::Ptr<recob::Slice> &slice : sliceVector){

    float nuScore = -1;
	totalSlices++;

	if (slice.key() >= slicePFPAssoc.size()) {
		std::cerr<<"Error: Slice key "<<slice.key()<<" is out of bounds for slicePFPAssoc (size =  "<<slicePFPAssoc.size()<<")"<<std::endl;
	       continue; // skip this slice
	}	       

	std::vector<art::Ptr<recob::PFParticle>> slicePFPs(slicePFPAssoc.at(slice.key()));
	if (slicePFPs.empty()){
		std::cerr<<"Warning: No PFPParticles associated with slice key "<< slice.key()<<std::endl;
		continue; // skip this slice
	}

	std::cout<<"Slice key: "<< slice.key()<<", Number of PFPs: "<< slicePFPs.size() << std::endl;
	std::cout<<"nuSliceKey = "<<nuSliceKey<<std::endl;

    float sliceVtxX = -999.;
    float sliceVtxY = -999.;
    float sliceVtxZ = -999.;

    // Loop over PFPs in slice
	for (const art::Ptr<recob::PFParticle> &slicePFP : slicePFPs){

        // Only care about primary neutrinos
		const bool isPrimary = (slicePFP->IsPrimary());
		const bool isNeutrino = (std::abs(slicePFP->PdgCode()) == 14);

		if (!(isPrimary && isNeutrino)){
			std::cout<<"Not a primary neutrino, skipping PFP!"<<std::endl;
			continue;
		}

		std::cout<<"Is a primary neutrino, not skipping!"<<std::endl;

        // Get metadata for nuScore, IsClearCosmic properties
		if (slicePFP.key() >= pfpMetadataAssoc.size()){
			std::cerr<<"PFP key "<<slicePFP.key()<<" out of bounds for pfpMetadataAssoc, size = "<<pfpMetadataAssoc.size()<<std::endl;
			continue;
		}

		std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metadataVec = pfpMetadataAssoc.at(slicePFP.key());

		if (metadataVec.empty()){
			std::cerr<<"No meta data found for PFParticle with key :"<<slice.key()<<std::endl;
			continue;
		}

		bool isClearCosmic = false;
		bool foundNuScore = false;
		bool foundClearCosmic = false;

		for (const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata : metadataVec){
			const auto &propertiesMap = metadata->GetPropertiesMap();

			if (propertiesMap.find("NuScore") != propertiesMap.end()){
				nuScore = propertiesMap.at("NuScore");
				foundNuScore = true;
			}

			else {
				std::cerr<<"NuScore not found for PFParticle with key "<<slice.key()<<std::endl;
			}

			if (propertiesMap.find("IsClearCosmic") != propertiesMap.end()){
				isClearCosmic = static_cast<bool>(propertiesMap.at("IsClearCosmic"));
				foundClearCosmic = true;
			}
		}

		if (foundNuScore) {
			std::cout<<"Found nuScore for this PFP: "<<nuScore<<". Current highestNuScore: "<<highestNuScore<<std::endl;
			fnuScore.push_back(nuScore);

			if (isClearCosmic){
				fCosmicNuScores.push_back(nuScore);
				std::cout<<"Found clear cosmic with nuScore "<<nuScore<<std::endl;
			}

			else {
				fNeutrinoNuScores.push_back(nuScore);
			}
		}

		else {
			std::cerr<<"Skipping PFParticle due to missing nuScore. Key:  "<<slicePFP.key()<<std::endl;
		}

		if (!foundClearCosmic) {
			std::cerr<<"Warning: IsClearCosmic property missing for PFParticle with key: "<<slicePFP.key()<<std::endl;
		}

		std::cout << "Neutrino slice detected! PDG codes of particles in this slice: ";	

		for (const auto &pfp : slicePFPs){
		    std::cout << pfp->PdgCode()<<", ";
		}
		std::cout<<std::endl;
	
		// We have found our neutrino!!
		totalNeutrinos++;

		if (nuScore > highestNuScore) {
			highestNuScore = nuScore;
			nuSliceKey = slice.key();
			nuID = slicePFP->Self();
			fNPrimaryParticles = slicePFPs.size();
			std::cout<<"Highest nuScore overwritten! New highestNuScore: "<<highestNuScore<<std::endl;
			std::cout<<"new nuSliceKey = "<<nuSliceKey<<std::endl;
			std::cout<<"new nuID = "<<nuID<<std::endl;
		}

		else {
			std::cout<<"This PFP nuScore is not the highest. Highest nuScore remains: "<<highestNuScore<<std::endl;
		}

        // Get vertex of PFP
        auto vertices = pfpVertexAssoc.at(slicePFP.key());

		if (!vertices.empty()) {
			const recob::Vertex& vertex = *vertices.at(0);
			auto const& vertexPos = vertex.position();
			
			sliceVtxX = vertexPos.X();
			sliceVtxY = vertexPos.Y();
			sliceVtxZ = vertexPos.Z();
			break;
	    }

	} // end loop over slice PFPs

    fSliceKey.push_back(slice.key());
    fSliceID.push_back(slice.id());
    fSliceNuScore.push_back(nuScore);
    //fSliceOpt0Score.push_back(opt0Score)
    fSliceVtxX.push_back(sliceVtxX);
    fSliceVtxY.push_back(sliceVtxY);
    fSliceVtxZ.push_back(sliceVtxZ);

    // Now, get hits in slice and loop over these hits

    std::vector<art::Ptr<recob::Hit>> sliceHits(sliceHitAssoc.at(slice.key()));
    int totalSliceHits = sliceHits.size();
    int trueNuSliceHits = 0;
    std::map<int, int> originVotes;

    for(const art::Ptr<recob::Hit> &hit : sliceHits){

        const HitTruthInfo truthInfo = getHitTruthInfo(hit);

        if(!truthInfo.valid){continue;}
        ++originVotes[truthInfo.origin];

        if(truthInfo.isBeamNeutrino){++trueNuSliceHits;}

        /*if (hitTruthAssns.isValid()){
            auto const& particles = hitTruthAssns.at(hit.key());

            for (const auto& truePart : particles){
                const art::Ptr<simb::MCTruth> hitMCTruth = piService->TrackIdToMCTruth_P(std::abs(truePart->TrackId()));

                if(hitMCTruth.isNonnull()){
                    int hitOrigin = static_cast<int>(hitMCTruth->Origin());
                    originVotes[hitOrigin]++;

                    if (hitOrigin == simb::kBeamNeutrino){
                        hitBelongsToNu = true;
                    }
                }
            }
        }*/

    } // end of loop over slice hits

    int bestOrigin = 0;
    int maxVotes = -1;
    for (auto const& [originType, count] : originVotes){
        if (count > maxVotes){
            maxVotes = count;
            bestOrigin = originType;
        }
    }

    fSliceTotalHits.push_back(totalSliceHits);
    fSliceTrueNuHits.push_back(trueNuSliceHits);
    fSliceTrueOrigin.push_back(bestOrigin);

   } // end loop over slices

// Define vector of PFPs in nuSlice
   //std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs(slicePFPAssoc.at(nuSliceKey));
   //std::cout<<"!!! Now looping through nuSlicePFPs of size = "<<nuSlicePFPs.size()<<" !!!"<<std::endl;
   //std::cout<<"Equal to nPFParticles = "<<fNPrimaryParticles<<std::endl;


// Track and Shower reco diagnostics

   art::ValidHandle<std::vector<recob::Track>> trackHandle = evt.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
   if (!trackHandle.isValid()){
	   std::cerr<<"Warning: recob::Track collection not found in event. "<<std::endl;
   }
   art::ValidHandle<std::vector<recob::Shower>> showerHandle = evt.getValidHandle<std::vector<recob::Shower>>(fShowerLabel);
   art::FindManyP<recob::Track> pfpTrackAssoc(pfpHandle, evt, fTrackLabel);
   art::FindManyP<recob::Shower> pfpShowerAssoc(pfpHandle, evt, fShowerLabel);
   art::FindManyP<recob::PFParticle> trackToPFPAssoc(trackHandle, evt, fTrackLabel);
   art::FindManyP<recob::PFParticle> showerToPFPAssoc(showerHandle, evt, fShowerLabel);

    for (const art::Ptr<recob::Slice>& slice : sliceVector)
    {
        const bool isNuSlice = (nuSliceKey >= 0 && static_cast<int>(slice.key()) == nuSliceKey);
        const auto slicePFPs = slicePFPAssoc.at(slice.key());

        for (const art::Ptr<recob::PFParticle>& pfp : slicePFPs)
        {
            // ------------------------------------------------
            // Basic PFP information
            // ------------------------------------------------

            fPfpKey.push_back(pfp.key());
            fPfpSelfID.push_back(pfp->Self());
            fPfpParentID.push_back(pfp->Parent());
            fPfpRecoPDG.push_back(pfp->PdgCode());
            fPfpSliceKey.push_back(slice.key());
            fPfpIsNuSlice.push_back(isNuSlice);
            fPfpIsPrimary.push_back(pfp->IsPrimary());
            fPfpNDaughters.push_back(pfp->NumDaughters());

            // ------------------------------------------------
            // TrackScore
            // ------------------------------------------------

            const float trackScore = GetPFPTrackScore(pfp, pfpMetadataAssoc);

            fPfpTrackScore.push_back(trackScore);
            fPfpHasTrackScore.push_back(trackScore > -9990.f);

            // ------------------------------------------------
            // Reco associations
            // ------------------------------------------------

            const auto tracks = pfpTrackAssoc.at(pfp.key());
            const auto showers = pfpShowerAssoc.at(pfp.key());
            const auto vertices = pfpVertexAssoc.at(pfp.key());

            fPfpNTracks.push_back(tracks.size());
            fPfpNShowers.push_back(showers.size());
            fPfpNVertices.push_back(vertices.size());

            // ------------------------------------------------
            // Track interpretation
            // ------------------------------------------------

            const TrackRecoInfo trackInfo = GetTrackRecoInfo(tracks);

            fPfpHasTrack.push_back(trackInfo.valid);
            fPfpHasUniqueTrack.push_back(trackInfo.uniqueTrack);
            fPfpTrackID.push_back(trackInfo.id);
            fPfpTrackLength.push_back(trackInfo.length);

            fPfpTrackStartX.push_back(trackInfo.startX);
            fPfpTrackStartY.push_back(trackInfo.startY);
            fPfpTrackStartZ.push_back(trackInfo.startZ);
            fPfpTrackEndX.push_back(trackInfo.endX);
            fPfpTrackEndY.push_back(trackInfo.endY);
            fPfpTrackEndZ.push_back(trackInfo.endZ);
            fPfpTrackStartDirX.push_back(trackInfo.startDirX);
            fPfpTrackStartDirY.push_back(trackInfo.startDirY);
            fPfpTrackStartDirZ.push_back(trackInfo.startDirZ);
            fPfpTrackEndDirX.push_back(trackInfo.endDirX);
            fPfpTrackEndDirY.push_back(trackInfo.endDirY);
            fPfpTrackEndDirZ.push_back(trackInfo.endDirZ);
            fPfpTrackTheta.push_back(trackInfo.theta);
            fPfpTrackPhi.push_back(trackInfo.phi);

            // ------------------------------------------------
            // Shower interpretation
            // ------------------------------------------------

            const ShowerRecoInfo showerInfo = GetShowerRecoInfo(showers);

            fPfpHasShower.push_back(showerInfo.valid);
            fPfpHasUniqueShower.push_back(showerInfo.uniqueShower);
            fPfpShowerID.push_back(showerInfo.id);
            fPfpShowerLength.push_back(showerInfo.length);

            fPfpShowerStartX.push_back(showerInfo.startX);
            fPfpShowerStartY.push_back(showerInfo.startY);
            fPfpShowerStartZ.push_back(showerInfo.startZ);
            fPfpShowerDirX.push_back(showerInfo.dirX);
            fPfpShowerDirY.push_back(showerInfo.dirY);
            fPfpShowerDirZ.push_back(showerInfo.dirZ);

            // ------------------------------------------------
            // Vertex interpretation
            // ------------------------------------------------

            fPfpHasVertex.push_back(!vertices.empty());
            fPfpHasUniqueVertex.push_back(vertices.size() == 1);

            if (vertices.size() == 1){
                const auto& vertex = vertices.front();
                const auto& pos = vertex->position();

                float x = pos.X();
                float y = pos.Y();
                float z = pos.Z();
                fPfpVertexX.push_back(x);
                fPfpVertexY.push_back(y);
                fPfpVertexZ.push_back(z);
            }
            else{
                fPfpVertexX.push_back(-9999.f);
                fPfpVertexY.push_back(-9999.f);
                fPfpVertexZ.push_back(-9999.f);
            }
            
            // ------------------------------------------------
            // Truth match THIS PFP
            // ------------------------------------------------

            PFPTruthMatch match;

            // get vector of pfpHits from helper function
            const std::vector<art::Ptr<recob::Hit>> pfpHits = getPFPHits(pfp);

            // Always store one value per selected PFP so the vectors stay aligned.

            if (!pfpHits.empty()) {
                const TruthMatchUtils::G4ID g4ID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHits, rollupUnsavedIDs);

                if (TruthMatchUtils::Valid(g4ID)) {
                    const simb::MCParticle* trueParticle = particleInventory->TrackIdToParticle_P(g4ID);

                    if(trueParticle){
                        match.valid = true;
                        match.trackID = trueParticle->TrackId();
                        match.pdg = trueParticle->PdgCode();

                         // Hit-count purity consistent with TruthMatchUtils.

                        std::size_t nMatchedHits = 0;

                        for (const art::Ptr<recob::Hit>& hit : pfpHits) {

                            const TruthMatchUtils::G4ID hitG4ID = TruthMatchUtils::TrueParticleID(clockData, hit, rollupUnsavedIDs);

                            if (TruthMatchUtils::Valid(hitG4ID) && hitG4ID == g4ID) {
                            ++nMatchedHits;
                            }
                        }
                        
                        match.nMatchedHits = nMatchedHits;
                        // How many hits match out of the total number of hits associated with this PFP?
                        match.purity = static_cast<float>(nMatchedHits) / static_cast<float>(pfpHits.size());
                        
                    }
                }
            }

            match.nHits = static_cast<int>(pfpHits.size());

            fPfpTrueTrackID.push_back(match.trackID);
            fPfpTruePDG.push_back(match.pdg);
            fPfpNHits.push_back(match.nHits);
            fPfpNMatchedHits.push_back(match.nMatchedHits);
            fPfpTruthPurity.push_back(match.purity);

        } // end loop over PFPs in slice 
    } // end loop over slices


// Tracks reco -> truth matching
/* std::cout<<"========= Track Truth Matching =========="<<std::endl;
fTruePfpPDG.clear();
fTrueTrackPDG.clear();
fMuonTrackScores.clear();
fProtonTrackScores.clear();
fPionTrackScores.clear();
std::map<int, TVector3> trackTrueVertices;

art::Handle<art::Assns<simb::MCParticle, recob::Hit>> hitTruthAssn;
evt.getByLabel("gaushitTruthMatch", hitTruthAssn);
art::FindManyP<recob::Hit> trackHitAssoc(trackHandle, evt, fTrackLabel);

for (const art::Ptr<recob::PFParticle>& slicePFP : nuSlicePFPs) { // loop through PFPs in nuSlice
    art::FindManyP<recob::Track> pfpTrackAssoc(pfpHandle, evt, fTrackLabel);

    std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(slicePFP.key());

    // Get clusters from pfp to cluster association
    //std::vector<art::Ptr<recob::Cluster>> clusters = pfpClusterAssoc.at(slicePFP.key());
    //art::Ptr<recob::Cluster> cluster = clusters[0];

    // Get hits associated with cluster (hits from PFP)
    //std::vector<art::Ptr<recob::Hit>> clusterHits = clusterHitAssoc.at(cluster.ID());

    // Only care about neutrino children (same condition for saving track/shower reco params)
    if (slicePFP->Parent() != static_cast<long unsigned int>(nuID)){
		std::cout<<"skipping this slicePFP (reason: not a neutrino child)"<<std::endl;
		continue;
    }

    if (tracks.empty()) 
	continue;
   
    if (tracks.size() != 1) // ensure there is only 1 track associated with this PFP
        continue;

for (const auto& track : tracks) { // loop over tracks (is this really necessary???)
  
    //std::vector<art::Ptr<recob::PFParticle>> trackPFPs = trackToPFPAssoc.at(track.key());
    //art::Ptr<recob::PFParticle> trackPFP = trackPFPs.front();
    //std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> trackMetadataVec = pfpMetadataAssoc.at(trackPFP.key());
    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadataVec = pfpMetadataAssoc.at(slicePFP.key());
    float trackScore = -1.0;
    if (pfpMetadataVec.empty()) {
	std::cerr << "No metadata found for PFP with key: " << slicePFP.key() << std::endl;
       }

        for (const auto& metadata : pfpMetadataVec) {
            const auto& propertiesMap = metadata->GetPropertiesMap();
      		//  std::cout << "\nTrack PFParticle Metadata Properties (ID " << track->ID() << "):" << std::endl;	 
	 if (propertiesMap.find("TrackScore") != propertiesMap.end()) {
                        trackScore = propertiesMap.at("TrackScore");
			std::cout<<"trackScore = "<<trackScore<<std::endl;
	 
    	 } 
	 else {
        	std::cerr << "No metadata found for PFParticle with key: " << slicePFP.key() << std::endl;
   	 }

	}

    std::vector<art::Ptr<recob::Hit>> trackHits = trackHitAssoc.at(track.key());    

    art::FindManyP<simb::MCParticle> hitToMCParticleAssn(trackHits, evt, "gaushitTruthMatch");
    std::map<int, int> mcParticleHitCount;
    
    for (size_t i_hit = 0; i_hit < trackHits.size(); i_hit++) { // loop over hits
        std::vector<art::Ptr<simb::MCParticle>> mcParticles = hitToMCParticleAssn.at(i_hit);
        for (const auto& mcParticle : mcParticles) {
            mcParticleHitCount[mcParticle->TrackId()]++;
        }
    }
   
    int bestTrackID = -1;
    int maxHitCount = 0;
    for (const auto& [trackID, hitCount] : mcParticleHitCount) { // find track with highest hit count
        if (hitCount > maxHitCount) {
            maxHitCount = hitCount;
            bestTrackID = trackID;
        }
    } 
    
    if (bestTrackID != -1) {
	art::ServiceHandle<cheat::ParticleInventoryService> piService;
        const simb::MCParticle* mcParticle = piService->TrackIdToParticle_P(bestTrackID);
	if (mcParticle && trackScore > 0.5) { // only truth match track if we classify this PFP as a track
	trueTrackPDGMap[track->ID()] = mcParticle->PdgCode();	

	 int pdgCode = mcParticle->PdgCode();
            if (pdgCode == -13) {
                fMuonTrackScores.push_back(trackScore);
            } else if (pdgCode == 2212) {
                fProtonTrackScores.push_back(trackScore);
            } else if (pdgCode == -211) {
                fPionTrackScores.push_back(trackScore);
            }

		float matchQuality = static_cast<float>(maxHitCount) / trackHits.size();
		//std::cout<<"match quality = "<<matchQuality<<std::endl;

		//size_t trackIndex = std::find(fTrackIDs.begin(), fTrackIDs.end(), track->ID()) - fTrackIDs.begin();
                //if (trackIndex < fTrackIDs.size()) {
		fTrueTrackPDG.push_back(mcParticle->PdgCode());
		fTruePfpPDG.push_back(mcParticle->PdgCode());

                std::cout << "Best matched MC Particle for Track: "
                          << "PDG: " << mcParticle->PdgCode()
                          << ", Track ID: " << mcParticle->TrackId()
                          << ", Hit Count: " << maxHitCount
			  << ", trackHits.size(): "<<trackHits.size()
			  <<", Match Quality: "<< matchQuality
			  <<", trackScore: "<<trackScore<< std::endl;

		// fIsLongestTrack.push_back(track->ID() == fLongestTrackID);
               // fIsClosestTrack.push_back(track->ID() == fClosestTrackID);

*/		/* if (track->ID() == fLongestTrackID) {
                        fLongestTrackTruePDG = mcParticle->PdgCode();
                        longestTrackCount++;
                        longestTrackPDGCounts[mcParticle->PdgCode()]++;

		                  
                        if (mcParticle->PdgCode() == 13) { 
                            longestTrackMuonCount++;
                        }
                        
                        std::cout << "=========================" << std::endl;
                        std::cout << "Best matched MC Particle for Longest Track: "
                                  << "MCParticle Track ID: " << mcParticle->TrackId()
                                  << ", PDG: " << mcParticle->PdgCode()
                                  << ", Hit Count: " << maxHitCount
                                  << ", Match Quality: " << matchQuality
                                  << std::endl;
                    }
*/
		/*	 if (track->ID() == fClosestTrackID) {
                        fClosestTrackTruePDG = mcParticle->PdgCode();
                        closestTrackCount++;
                        closestTrackPDGCounts[mcParticle->PdgCode()]++;
                        
                        if (abs(mcParticle->PdgCode()) == 321) { // Kaon
                            closestTrackKaonCount++;
                        }
                        
                        std::cout << "=========================" << std::endl;
                        std::cout << "Best matched MC Particle for Closest Track: "
                                  << "MCParticle Track ID: " << mcParticle->TrackId()
                                  << ", PDG: " << mcParticle->PdgCode()
                                  << ", Hit Count: " << maxHitCount
                                  << ", Match Quality: " << matchQuality
                                  << std::endl;
         	           	}	
*//*
				//}
			}
		}
	} // end loop over tracks
} // end loop over slice PFPs
*/


/*
//Showers reco -> truth matching
std::cout<<"======== Shower Truth Matching ==========="<<std::endl;
fTrueShowerPDG.clear();

art::FindManyP<recob::Hit> showerHitAssoc(showerHandle, evt, fShowerLabel);
for (const art::Ptr<recob::PFParticle>& slicePFP : nuSlicePFPs) {
    art::FindManyP<recob::Shower> pfpShowerAssoc(pfpHandle, evt, fShowerLabel);

    std::vector<art::Ptr<recob::Shower>> showers = pfpShowerAssoc.at(slicePFP.key());

    // Only care about neutrino children (same as track logic)
    if (slicePFP->Parent() != static_cast<long unsigned int>(nuID)){
	std::cout<<"skipping this slicePFP (reason: not a neutrino child)"<<std::endl;
	continue;
    }

    if (showers.empty())
        continue;

    if (showers.size() != 1)
        continue;

    art::Ptr<recob::Shower> shower = showers[0];

    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadataVec = pfpMetadataAssoc.at(slicePFP.key());
    float trackScore = -1.0;
    if (pfpMetadataVec.empty()) {
	std::cerr << "No metadata found for PFP with key: " << slicePFP.key() << std::endl;
       }

        for (const auto& metadata : pfpMetadataVec) {
            const auto& propertiesMap = metadata->GetPropertiesMap();
      		//  std::cout << "\nTrack PFParticle Metadata Properties (ID " << track->ID() << "):" << std::endl;	 
	 if (propertiesMap.find("TrackScore") != propertiesMap.end()) {
                        trackScore = propertiesMap.at("TrackScore");
			std::cout<<"trackScore = "<<trackScore<<std::endl;
	 
    	 } 
	 else {
        	std::cerr << "No metadata found for PFParticle with key: " << slicePFP.key() << std::endl;
   	 }

	}
	
    std::vector<art::Ptr<recob::Hit>> showerHits = showerHitAssoc.at(shower.key());
    art::FindManyP<simb::MCParticle> showerHitToMCParticleAssn(showerHits, evt, "gaushitTruthMatch");
    
    std::map<int, int> mcParticleHitCount;
    
for (size_t i_hit = 0; i_hit < showerHits.size(); i_hit++) { // loop over shower hits
    std::vector<art::Ptr<simb::MCParticle>> mcParticles = showerHitToMCParticleAssn.at(i_hit); 

 for (const auto& mcParticle : mcParticles) {
        mcParticleHitCount[mcParticle->TrackId()]++;
    }
}

    int bestTrackID = -1;
    int maxHitCount = 0;
    for (const auto& [trackID, hitCount] : mcParticleHitCount) { // find shower with highest number of hits
        if (hitCount > maxHitCount) {
            maxHitCount = hitCount;
            bestTrackID = trackID;
        }
    }
    if (bestTrackID != -1) {
           art::ServiceHandle<cheat::ParticleInventoryService> piService;
    const simb::MCParticle* mcParticle = piService->TrackIdToParticle_P(bestTrackID);    

	if (mcParticle && trackScore < 0.5) { // only truth match shower if we classify this PFP as a shower

		float matchQuality = static_cast<float>(maxHitCount) / showerHits.size();
                //frecoTruePDG[shower->ID()] = mcParticle->PdgCode();
		fTrueShowerPDG.push_back(mcParticle->PdgCode());
		fTruePfpPDG.push_back(mcParticle->PdgCode());
                std::cout << "Best matched MC Particle for Shower: "
                          << "PDG: " << mcParticle->PdgCode()
                          << ", Track ID: " << mcParticle->TrackId()
                          << ", Hit Count: " << maxHitCount 
			  << ", showerHits.size(): "<<showerHits.size()
			  << ", Match Quality: "<<matchQuality << std::endl;
            }
        }
}*/

// ============================================================================
// MC TRUTH PARAMETERS
// ============================================================================

std::cout<< "------------ MC TRUTH Parameters ----------------"<<std::endl;

art::Handle<std::vector<simb::MCTruth> > mctruthListHandle;
std::vector<art::Ptr<simb::MCTruth> > mclist;
if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle)){
      art::fill_ptr_vector(mclist, mctruthListHandle);
}

art::FindManyP<simb::MCParticle> fmpart( mctruthListHandle, evt, "largeant" );

// corsika MC truth information

/*art::ValidHandle<std::vector<simb::MCTruth>> cosmicMCTruthListHandle = evt.getValidHandle<std::vector<simb::MCTruth>>("corsika");
std::vector<art::Ptr<simb::MCTruth>> cosmicMCTruthVector;
if (cosmicMCTruthListHandle.isValid()){
    art::fill_ptr_vector(cosmicMCTruthVector, cosmicMCTruthListHandle);
}

art::FindManyP<simb::MCParticle> cosmicGeantAssoc(cosmicMCTruthListHandle, evt, "largeant");*/

++nTotEvents;

// ---------------------------------------------------------------------------
// Number of MCTruth interactions in this art::Event
// ---------------------------------------------------------------------------

if (mclist.size() > 1) {
    ++nMultiEvents;
}

if (mclist.size() == 1) {
    ++nSingleEvents;
}

std::cout<< "mclist size = "<< mclist.size()<< std::endl;


// ===========================================================================
// Loop over neutrino interactions
// ===========================================================================

for (size_t i_truth = 0; i_truth < mclist.size(); ++i_truth)
{
    art::Ptr<simb::MCTruth> truth = mclist.at(i_truth);

    const auto& neutrino = truth->GetNeutrino();
    const auto& nu       = neutrino.Nu();
    const auto origin = truth->Origin();

    // -----------------------------------------------------------------------
    // Save interaction-level information
    // -----------------------------------------------------------------------

    std::cout
        << "\n========== MCTruth " << i_truth << " =========="<< std::endl;

    std::cout<< "MCTruth Origin: "<<origin<< std::endl;
    std::cout<< "Neutrino PDG: "<< nu.PdgCode()<< std::endl;
    std::cout<< "Interaction mode: "<< neutrino.Mode()<< std::endl;
    std::cout<< "Interaction type: "<< neutrino.InteractionType()<< std::endl;
    std::cout<< "CCNC: "<< neutrino.CCNC()<< std::endl;
    std::cout<< "Neutrino energy: "<< nu.E()<< std::endl;
    std::cout<< "Neutrino vertex: ("<< nu.Vx() << ", "<< nu.Vy() << ", "<< nu.Vz() << ")"<< std::endl;
    std::cout<< "Target: "<< neutrino.Target()<< std::endl;
    std::cout<< "Generator neutrino TrackID: "<< nu.TrackId()<< std::endl;

    trueOrigin.push_back(origin);
    trueW.push_back(neutrino.W());
    trueX.push_back(neutrino.X());
    trueY.push_back(neutrino.Y());
    trueQSqr.push_back(neutrino.QSqr());
    truePt.push_back(neutrino.Pt());
    trueTheta.push_back(neutrino.Theta());
    trueNuPDG.push_back(nu.PdgCode());
    trueNuEnergy.push_back(nu.E());
    trueCCNC.push_back(neutrino.CCNC());
    trueIntMode.push_back(neutrino.Mode());
    trueIntType.push_back(neutrino.InteractionType());
    trueTargetPDG.push_back(neutrino.Target());
    trueNuTrackID.push_back(nu.TrackId());
    trueNuVtxX.push_back(nu.Vx());
    trueNuVtxY.push_back(nu.Vy());
    trueNuVtxZ.push_back(nu.Vz());

    // -----------------------------------------------------------------------
    // All simulated particles associated with THIS MCTruth interaction.
    //
    // This collection can be huge because it contains descendants from
    // detector interactions, ionisation, scattering, etc.
    // -----------------------------------------------------------------------

    std::vector<art::Ptr<simb::MCParticle>> assocParticles = fmpart.at(i_truth);

    std::cout<< "Total associated Geant4 MCParticles = "<< assocParticles.size()<< std::endl;

    // --------------------------------------------
    // Debugging: print out process counts for assocParticles in this MCTruth
    // --------------------------------------------

    std::map<std::string, int> processCounts;

    for (const auto& particle : assocParticles)
    {
        if (!particle) {
            continue;
        }

        ++processCounts[
            particle->Process()
        ];
    }

    std::cout<< "\nProcesses in MCTruth "<< i_truth<< ":"<< std::endl;

    for (const auto& [process, count] : processCounts)
    {
        std::cout
            << "    "
            << process
            << " : "
            << count
            << std::endl;
    }

    // -----------------------------------------------------------------------
    // Build fast genealogy lookup tables
    //
    // particleByTrackID:
    //
    //     TrackID -> particle
    //
    // childrenByMotherID:
    //
    //     mother TrackID -> all stored children
    //
    // -----------------------------------------------------------------------

    std::unordered_map<
        int,
        art::Ptr<simb::MCParticle>
    > particleByTrackID;


    std::unordered_map<
        int,
        std::vector<art::Ptr<simb::MCParticle>>
    > childrenByMotherID;


    for (const auto& particle : assocParticles)
    {
        if (!particle) {
            continue;
        }

        particleByTrackID[
            particle->TrackId()
        ] = particle;

        childrenByMotherID[
            particle->Mother()
        ].push_back(particle);
    }


    // -----------------------------------------------------------------------
    // Save where this MCTruth starts in the flat particle vectors.
    //
    // This isn't essential because trueMCTruthIndex also identifies the
    // interaction, but it makes later indexing extremely convenient.
    // -----------------------------------------------------------------------

    const int firstParticleIndex = static_cast<int>(truePDG.size());
    trueParticleStartIndex.push_back(firstParticleIndex);

    // -----------------------------------------------------------------------
    // Find all primary particles for this interaction.
    //
    // LArSoft defines detector-simulation primaries with Process()=="primary".
    // -----------------------------------------------------------------------

    std::queue<
        std::pair<
            art::Ptr<simb::MCParticle>,
            int
        >
    > particlesToProcess;

    int nPrimariesThisTruth = 0;

    for (const auto& particle : assocParticles)
    {
        if (!particle) {
            continue;
        }

        if (particle->Process() != "primary") {
            continue;
        }

        // generation = 0 for an interaction primary
        particlesToProcess.push({
            particle,
            0
        });

        ++nPrimariesThisTruth;
    }

    std::cout<< "Primary particles = "<< nPrimariesThisTruth<< std::endl;

    // -----------------------------------------------------------------------
    // Avoid duplicate saving even if some bizarre genealogy occurs.
    // -----------------------------------------------------------------------

    std::unordered_set<int> savedTrackIDs;
    int nSavedThisTruth = 0;

    // =======================================================================
    // Traverse:
    //
    // primary
    //   -> decay child
    //       -> decay child
    //           -> ...
    //
    // No fixed depth.
    // =======================================================================

    while (!particlesToProcess.empty())
    {
        const auto current = particlesToProcess.front();

        particlesToProcess.pop();

        const art::Ptr<simb::MCParticle>& particle = current.first;

        const int generation = current.second;

        if (!particle) {
            continue;
        }

        const int trackID = particle->TrackId();

        int nHits = 0;
        int nHitsU = 0;
        int nHitsV = 0;
        int nHitsZ = 0;
        auto it = trueHitCounts.find(trackID);

        if (it != trueHitCounts.end()){
            nHits = it->second.total;
            nHitsU = it->second.U;
            nHitsV = it->second.V;
            nHitsZ = it->second.Z;
        }

        if (nHits != nHitsU + nHitsV + nHitsZ){
        std::cerr
            << "Hit-view mismatch for TrackID "
            << trackID
            << ": total = " << nHits
            << ", U = " << nHitsU
            << ", V = " << nHitsV
            << ", Z = " << nHitsZ
            << std::endl;
        }

        // Already handled?
        if (savedTrackIDs.count(trackID) != 0) {
            continue;
        }

        savedTrackIDs.insert(trackID);


        // -------------------------------------------------------------------
        // Mother information
        // -------------------------------------------------------------------

        const int motherTrackID = particle->Mother();
        int motherPDGCode = 0;

        auto motherIt = particleByTrackID.find(motherTrackID); 

        if (motherIt != particleByTrackID.end()) {
            motherPDGCode = motherIt->second->PdgCode();
        }


        // -------------------------------------------------------------------
        // Determine how many STORED decay children this particle has.
        //
        // NumberDaughters() is still worth saving separately because that is
        // the total daughter information stored on MCParticle.
        // -------------------------------------------------------------------

        int nStoredDecayDaughters = 0;

        auto childrenIt = childrenByMotherID.find(trackID);

        if (childrenIt != childrenByMotherID.end())
        {
            for (const auto& child : childrenIt->second)
            {
                if (!child) {
                    continue;
                }

                if (IsDecayProcess(child->Process())) {
                    ++nStoredDecayDaughters;
                }
            }
        }


        // ===================================================================
        // SAVE THIS PARTICLE
        // ===================================================================

        truePDG.push_back(particle->PdgCode());
        trueTrackID.push_back(trackID);
        trueMotherTrackID.push_back(motherTrackID);
        trueMotherPDG.push_back(motherPDGCode);
        trueMCTruthIndex.push_back(static_cast<int>(i_truth));
        trueGeneration.push_back(generation);
        trueIsPrimary.push_back(generation == 0);
        trueIsDecayProduct.push_back(generation > 0);
        trueNHits.push_back(nHits);
        trueNHitsU.push_back(nHitsU);
        trueNHitsV.push_back(nHitsV);
        trueNHitsZ.push_back(nHitsZ);

        // -------------------------------------------------------------------
        // Process information
        // -------------------------------------------------------------------

        trueProcess.push_back(particle->Process());

        trueEndProcess.push_back(particle->EndProcess());

        trueStatusCode.push_back(particle->StatusCode());

        // -------------------------------------------------------------------
        // Daughter information
        // -------------------------------------------------------------------

        trueNDaughters.push_back(particle->NumberDaughters());

        trueNStoredDecayDaughters.push_back(nStoredDecayDaughters);


        // -------------------------------------------------------------------
        // Particle mass
        // -------------------------------------------------------------------

        trueMass.push_back(particle->Mass());

        // -------------------------------------------------------------------
        // Trajectory / kinematic information
        // -------------------------------------------------------------------

        const unsigned int nTrajectoryPoints = particle->NumberTrajectoryPoints();

        trueNTrajectoryPoints.push_back(nTrajectoryPoints);


        if (nTrajectoryPoints > 0)
        {
            // Start position

            trueStartX.push_back(particle->Vx());
            trueStartY.push_back(particle->Vy());
            trueStartZ.push_back(particle->Vz());
            trueStartT.push_back(particle->T());

            // Initial momentum

            truePx.push_back(particle->Px());
            truePy.push_back(particle->Py());
            truePz.push_back(particle->Pz());
            trueP.push_back(particle->P());
            trueE.push_back(particle->E());

            // End position

            trueEndX.push_back(particle->EndX());
            trueEndY.push_back(particle->EndY());
            trueEndZ.push_back(particle->EndZ());
            trueEndT.push_back(particle->EndT());

            // End momentum

            trueEndPx.push_back(particle->EndPx());
            trueEndPy.push_back(particle->EndPy());
            trueEndPz.push_back(particle->EndPz());
            trueEndE.push_back(particle->EndE());


            // ---------------------------------------------------------------
            // Actual simulated trajectory length
            // ---------------------------------------------------------------

            double trajectoryLength = 0.0;

            for (unsigned int iPoint = 1; iPoint < nTrajectoryPoints; ++iPoint){
                const auto& previous = particle->Position(iPoint - 1);

                const auto& currentPoint = particle->Position(iPoint);

                const double dx = currentPoint.X() - previous.X();

                const double dy = currentPoint.Y() - previous.Y();

                const double dz = currentPoint.Z() - previous.Z();

                trajectoryLength +=
                    std::sqrt(
                        dx * dx +
                        dy * dy +
                        dz * dz
                    );
            }

            trueTrajectoryLength.push_back(trajectoryLength);
        }

        else
        {
            // ---------------------------------------------------------------
            // Keep ALL flat vectors aligned.
            // ---------------------------------------------------------------

            trueStartX.push_back(-9999.);
            trueStartY.push_back(-9999.);
            trueStartZ.push_back(-9999.);
            trueStartT.push_back(-9999.);

            truePx.push_back(-9999.);
            truePy.push_back(-9999.);
            truePz.push_back(-9999.);
            trueP.push_back(-9999.);
            trueE.push_back(-9999.);

            trueEndX.push_back(-9999.);
            trueEndY.push_back(-9999.);
            trueEndZ.push_back(-9999.);
            trueEndT.push_back(-9999.);

            trueEndPx.push_back(-9999.);
            trueEndPy.push_back(-9999.);
            trueEndPz.push_back(-9999.);
            trueEndE.push_back(-9999.);

            trueTrajectoryLength.push_back(-9999.);
        }

        ++nSavedThisTruth;

        // -------------------------------------------------------------------
        // DEBUG OUTPUT
        // -------------------------------------------------------------------

        std::cout
            << std::string(2 * generation, ' ')
            << "PDG = " << particle->PdgCode()
            << " TrackID = " << trackID
            << " Mother = " << motherTrackID
            << " MotherPDG = " << motherPDGCode
            << " generation = " << generation
            << " process = " << particle->Process()
            << " end process = " << particle->EndProcess()
            << " daughters = " << particle->NumberDaughters()
            << " retained decay daughters = "
            << nStoredDecayDaughters
            << std::endl;

        // ===================================================================
        // Look for children to retain.
        // ===================================================================

        if (childrenIt == childrenByMotherID.end()) {
            continue;
        }


        for (const auto& child : childrenIt->second)
        {
            if (!child) {
                continue;
            }


            if (!IsDecayProcess(child->Process()))
            {
                continue;
            }

            particlesToProcess.push({
                child,
                generation + 1
            });
        }
    }


    // -----------------------------------------------------------------------
    // Interaction-level particle counts
    // -----------------------------------------------------------------------

    trueNPrimaryParticles.push_back(nPrimariesThisTruth);

    trueNSavedParticles.push_back(nSavedThisTruth);


    std::cout
        << "Saved "
        << nSavedThisTruth
        << " particles for MCTruth "
        << i_truth
        << " ("
        << nPrimariesThisTruth
        << " primaries)"
        << std::endl;
}

// ===========================================================================
// Truth Matching: True MCParticles -> Reconstructed PFPs
// ===========================================================================
for (size_t iTrue = 0; iTrue < trueTrackID.size(); ++iTrue)
{
    int nMatchedPfps = 0;
    int bestPfpIndex = -1;
    int bestNMatchedHits = -1;

    bool reconstructedInDefaultNuSlice = false;

    for (size_t iPfp = 0; iPfp < fPfpTrueTrackID.size(); ++iPfp){
        if (fPfpTrueTrackID.at(iPfp) != trueTrackID.at(iTrue)) {
            continue;
        }

        ++nMatchedPfps;

        if (fPfpIsNuSlice.at(iPfp)) {
            reconstructedInDefaultNuSlice = true;
        }

        if (fPfpNMatchedHits.at(iPfp) > bestNMatchedHits)
        {
            bestNMatchedHits = fPfpNMatchedHits.at(iPfp);
            bestPfpIndex = static_cast<int>(iPfp);
        }
    }

    fTrueIsReconstructed.push_back(nMatchedPfps > 0);
    fTrueIsReconstructedInNuSlice.push_back(reconstructedInDefaultNuSlice);
    fTrueNMatchedPfps.push_back(nMatchedPfps);
    fTrueBestRecoPfpIdx.push_back(bestPfpIndex);

    if (bestPfpIndex >= 0)
    {
        fTrueBestRecoTrackScore.push_back(fPfpTrackScore.at(bestPfpIndex));
        fTrueBestRecoHasTrack.push_back(fPfpHasTrack.at(bestPfpIndex));
        fTrueBestRecoHasShower.push_back(fPfpHasShower.at(bestPfpIndex));
    }
    else
    {
        fTrueBestRecoTrackScore.push_back(-9999.f);
        fTrueBestRecoHasTrack.push_back(0);
        fTrueBestRecoHasShower.push_back(0);
    }
}


//std::cin.get();
//  int nGeniePrimaries = 0, nGEANTparticles = 0, nMCNeutrinos = 0;


/*
 const std::vector<art::Ptr<simb::MCTruth>> truthVec = particleInventory->MCTruthVector_Ps();

  std::cout << std::setprecision(1) << std::fixed;
  if (fVerbose) {
    for (auto const &truth : truthVec) {
      std::cout << "Truth: " << truth << std::endl;
      if (truth->NeutrinoSet()) {
        const simb::MCNeutrino neutrino = truth->GetNeutrino();
        std::cout << "Neutrino: " << neutrino << std::endl;

        const simb::MCParticle nu = neutrino.Nu();
        std::cout << "X: " << nu.Vx() << " Y: " << nu.Vy() << " Z " << nu.Vz() << std::endl;
      } // truth->NeutrinoSet
    } // fVerbose
  } // truth: truthVec
  std::cout << std::setprecision(2) << std::fixed;
*/
/// std::cout<<nuSliceKey<<std::endl;

// DIAGNOSTIC OUT
/*const size_t nPfp = fPfpKey.size();

std::cout
    << "\n========== PFP VECTOR SANITY ==========\n"
    << "pfpKey:              " << nPfp << "\n"
    << "pfpTrackScore:       " << fPfpTrackScore.size() << "\n"
    << "pfpNTracks:          " << fPfpNTracks.size() << "\n"
    << "pfpHasTrack:         " << fPfpHasTrack.size() << "\n"
    << "pfpHasUniqueTrack:   " << fPfpHasUniqueTrack.size() << "\n"
    << "pfpTrackLength:      " << fPfpTrackLength.size() << "\n"
    << "pfpNShowers:         " << fPfpNShowers.size() << "\n"
    << "pfpHasShower:        " << fPfpHasShower.size() << "\n"
    << "pfpHasUniqueShower:  " << fPfpHasUniqueShower.size() << "\n"
    << "pfpShowerLength:     " << fPfpShowerLength.size() << "\n"
    << "pfpNVertices:        " << fPfpNVertices.size() << "\n"
    << "pfpVertexX:          " << fPfpVertexX.size() << "\n"
    << "pfpTrueTrackID:      " << fPfpTrueTrackID.size() << "\n"
    << "pfpTruthPurity:      " << fPfpTruthPurity.size() << "\n"
    << "=======================================\n";*/

fTree->Fill();
 
std::cout<<"trueNuEnergy size = "<<trueNuEnergy.size()<<std::endl;
std::cout<<"trueIntMode size = "<<trueIntMode.size()<<std::endl;
std::cout<<"trueIntType size = "<<trueIntType.size()<<std::endl;
std::cout<<"trueCCNC size = "<<trueCCNC.size()<<std::endl;

    trueOrigin.clear();
    trueW.clear();
    trueX.clear();
    trueY.clear();
    trueQSqr.clear();
    truePt.clear();
    trueTheta.clear();
    trueNuPDG.clear();
    trueNuTrackID.clear();
    trueNuVtxX.clear();
    trueNuVtxY.clear();
    trueNuVtxZ.clear();
    trueNuEnergy.clear();
    trueIntMode.clear();
    trueIntType.clear();
    trueCCNC.clear();
    trueTargetPDG.clear();

    truePDG.clear();
    trueTrackID.clear();
    trueMotherPDG.clear();
    trueMotherTrackID.clear();
    trueGeneration.clear();
    trueIsPrimary.clear();
    trueIsDecayProduct.clear();
    trueMCTruthIndex.clear();
    trueNHits.clear();
    trueNHitsU.clear();
    trueNHitsV.clear();
    trueNHitsZ.clear();
    trueP.clear();
    trueMass.clear();
    trueStartX.clear();
    trueStartY.clear();
    trueStartZ.clear();
    trueStartT.clear();
    truePx.clear();
    truePy.clear();
    truePz.clear();
    trueE.clear();
    trueEndX.clear();
    trueEndY.clear();
    trueEndZ.clear();
    trueEndT.clear();
    trueEndPx.clear();
    trueEndPy.clear();
    trueEndPz.clear();
    trueEndE.clear();
    trueNTrajectoryPoints.clear();
    trueTrajectoryLength.clear();
    trueProcess.clear();
    trueEndProcess.clear();
    trueStatusCode.clear();
    trueNDaughters.clear();
    trueNStoredDecayDaughters.clear();
    trueNPrimaryParticles.clear();
    trueNSavedParticles.clear();
    trueParticleStartIndex.clear();
}


void hyperon::AnalyzeEvents::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree> ("tree", "Output TTree");

  //add branches here
  
  // MC truth parameters
  fTree->Branch("eventID", &fEventID);
  fTree ->Branch("run", &fRun, "run/I");
  fTree ->Branch("subrun", &fSubRun, "subrun/I");

  fTree->Branch("trueOrigin", &trueOrigin);
  fTree->Branch("trueW", &trueW);
  fTree->Branch("trueX", &trueX);
  fTree->Branch("trueY", &trueY);
  fTree->Branch("trueQSqr", &trueQSqr);
  fTree->Branch("truePt", &truePt);
  fTree->Branch("trueTheta", &trueTheta);
  fTree->Branch("trueNuPDG", &trueNuPDG);
  fTree->Branch("trueNuTrackID", &trueNuTrackID);
  fTree->Branch("trueNuVtxX", &trueNuVtxX);
  fTree->Branch("trueNuVtxY", &trueNuVtxY);
  fTree->Branch("trueNuVtxZ", &trueNuVtxZ);
  fTree->Branch("trueNuEnergy", &trueNuEnergy);
  fTree->Branch("trueCCNC", &trueCCNC);
  fTree->Branch("trueIntMode", &trueIntMode);
  fTree->Branch("trueIntType", &trueIntType);
  fTree->Branch("trueTargetPDG", &trueTargetPDG);

  fTree->Branch("truePDG", &truePDG);
  fTree->Branch("trueTrackID", &trueTrackID);
  fTree->Branch("trueMotherPDG", &trueMotherPDG);
  fTree->Branch("trueMotherTrackID", &trueMotherTrackID);
  fTree->Branch("trueGeneration", &trueGeneration);
  fTree->Branch("trueIsPrimary", &trueIsPrimary);
  fTree->Branch("trueIsDecayProduct", &trueIsDecayProduct);
  fTree->Branch("trueMCTruthIndex", &trueMCTruthIndex);
  fTree->Branch("trueNHits", &trueNHits);
  fTree->Branch("trueNHitsU", &trueNHitsU);
  fTree->Branch("trueNHitsV", &trueNHitsV);
  fTree->Branch("trueNHitsZ", &trueNHitsZ);
  fTree->Branch("trueP", &trueP);
  fTree->Branch("trueMass", &trueMass);
  fTree->Branch("trueStartX", &trueStartX);
  fTree->Branch("trueStartY", &trueStartY);
  fTree->Branch("trueStartZ", &trueStartZ);
  fTree->Branch("trueStartT", &trueStartT);
  fTree->Branch("truePx", &truePx);
  fTree->Branch("truePy", &truePy);
  fTree->Branch("truePz", &truePz);
  fTree->Branch("trueE", &trueE);
  fTree->Branch("trueEndX", &trueEndX);
  fTree->Branch("trueEndY", &trueEndY);
  fTree->Branch("trueEndZ", &trueEndZ);
  fTree->Branch("trueEndT", &trueEndT);
  fTree->Branch("trueEndPx", &trueEndPx);
  fTree->Branch("trueEndPy", &trueEndPy);
  fTree->Branch("trueEndPz", &trueEndPz);
  fTree->Branch("trueEndE", &trueEndE);
  fTree->Branch("trueNTrajectoryPoints", &trueNTrajectoryPoints);
  fTree->Branch("trueTrajectoryLength", &trueTrajectoryLength);
  fTree->Branch("trueProcess", &trueProcess);
  fTree->Branch("trueEndProcess", &trueEndProcess);
  fTree->Branch("trueStatusCode", &trueStatusCode);
  fTree->Branch("trueNDaughters", &trueNDaughters);
  fTree->Branch("trueNStoredDecayDaughters", &trueNStoredDecayDaughters);
  fTree->Branch("trueNPrimaryParticles", &trueNPrimaryParticles);
  fTree->Branch("trueNSavedParticles", &trueNSavedParticles);
  fTree->Branch("trueParticleStartIndex", &trueParticleStartIndex);

  // reco parameters
  fTree->Branch("sliceKey", &fSliceKey);
  fTree->Branch("sliceID", &fSliceID);
  fTree->Branch("sliceNuScore", &fSliceNuScore);
  fTree->Branch("sliceTotalHits", &fSliceTotalHits);
  fTree->Branch("sliceTrueNuHits", &fSliceTrueNuHits);
  fTree->Branch("sliceTrueOrigin", &fSliceTrueOrigin);
  fTree->Branch("eventTotalTrueNuHits", &fEventTotalTrueNuHits, "eventTotalTrueNuHits/I");
  fTree->Branch("sliceVtxX", &fSliceVtxX);
  fTree->Branch("sliceVtxY", &fSliceVtxY);
  fTree->Branch("sliceVtxZ", &fSliceVtxZ);

  fTree->Branch("pfpKey", &fPfpKey);
    fTree->Branch("pfpSelfID", &fPfpSelfID);
    fTree->Branch("pfpParentID", &fPfpParentID);
    fTree->Branch("pfpRecoPDG", &fPfpRecoPDG);
    fTree->Branch("pfpSliceKey", &fPfpSliceKey);
    fTree->Branch("pfpIsNuSlice", &fPfpIsNuSlice);
    fTree->Branch("pfpIsPrimary", &fPfpIsPrimary);
    fTree->Branch("pfpNDaughters", &fPfpNDaughters);
    fTree->Branch("pfpTrackScore", &fPfpTrackScore);
    fTree->Branch("pfpHasTrackScore", &fPfpHasTrackScore);
    fTree->Branch("pfpNTracks", &fPfpNTracks);
    fTree->Branch("pfpNShowers", &fPfpNShowers);
    fTree->Branch("pfpHasTrack", &fPfpHasTrack);
    fTree->Branch("pfpHasShower", &fPfpHasShower);
    fTree->Branch("pfpHasUniqueTrack", &fPfpHasUniqueTrack);
    fTree->Branch("pfpHasUniqueShower", &fPfpHasUniqueShower);
    fTree->Branch("pfpNVertices", &fPfpNVertices);
    fTree->Branch("pfpHasVertex", &fPfpHasVertex);
    fTree->Branch("pfpHasUniqueVertex", &fPfpHasUniqueVertex);
    fTree->Branch("pfpTrackID", &fPfpTrackID);
  fTree->Branch("pfpTrackLength", &fPfpTrackLength);
  fTree->Branch("pfpTrackStartX", &fPfpTrackStartX);
  fTree->Branch("pfpTrackStartY", &fPfpTrackStartY);
  fTree->Branch("pfpTrackStartZ", &fPfpTrackStartZ);
  fTree->Branch("pfpTrackEndX", &fPfpTrackEndX);
  fTree->Branch("pfpTrackEndY", &fPfpTrackEndY);
  fTree->Branch("pfpTrackEndZ", &fPfpTrackEndZ);
  fTree->Branch("pfpTrackStartDirX", &fPfpTrackStartDirX);
  fTree->Branch("pfpTrackStartDirY", &fPfpTrackStartDirY);
  fTree->Branch("pfpTrackStartDirZ", &fPfpTrackStartDirZ);
  fTree->Branch("pfpTrackEndDirX", &fPfpTrackEndDirX);
  fTree->Branch("pfpTrackEndDirY", &fPfpTrackEndDirY);
  fTree->Branch("pfpTrackEndDirZ", &fPfpTrackEndDirZ);
  fTree->Branch("pfpTrackVertexDirX", &fPfpTrackVertexDirX);
  fTree->Branch("pfpTrackVertexDirY", &fPfpTrackVertexDirY);
  fTree->Branch("pfpTrackVertexDirZ", &fPfpTrackVertexDirZ);
  fTree->Branch("pfpTrackTheta", &fPfpTrackTheta);
  fTree->Branch("pfpTrackPhi", &fPfpTrackPhi);
  fTree->Branch("pfpShowerID", &fPfpShowerID);
  fTree->Branch("pfpShowerLength", &fPfpShowerLength);
  fTree->Branch("pfpShowerStartX", &fPfpShowerStartX);
  fTree->Branch("pfpShowerStartY", &fPfpShowerStartY);
  fTree->Branch("pfpShowerStartZ", &fPfpShowerStartZ);
  fTree->Branch("pfpShowerDirX", &fPfpShowerDirX);
  fTree->Branch("pfpShowerDirY", &fPfpShowerDirY);
  fTree->Branch("pfpShowerDirZ", &fPfpShowerDirZ);
  fTree->Branch("pfpVertexX", &fPfpVertexX);
  fTree->Branch("pfpVertexY", &fPfpVertexY);
  fTree->Branch("pfpVertexZ", &fPfpVertexZ);
  fTree->Branch("pfpTrueTrackID", &fPfpTrueTrackID);
  fTree->Branch("pfpTruePDG", &fPfpTruePDG);
  fTree->Branch("pfpNHits", &fPfpNHits);
  fTree->Branch("pfpNMatchedHits", &fPfpNMatchedHits);
  fTree->Branch("pfpTruthPurity", &fPfpTruthPurity);
  fTree->Branch("trueIsReconstructed", &fTrueIsReconstructed);
  fTree->Branch("trueIsReconstructedInNuSlice", &fTrueIsReconstructedInNuSlice);
  fTree->Branch("trueNMatchedPfps", &fTrueNMatchedPfps);
  fTree->Branch("trueBestRecoPfpIdx", &fTrueBestRecoPfpIdx);
  fTree->Branch("trueBestRecoTrackScore", &fTrueBestRecoTrackScore);
  fTree->Branch("trueBestRecoHasTrack", &fTrueBestRecoHasTrack);
  fTree->Branch("trueBestRecoHasShower", &fTrueBestRecoHasShower);

  fTree->Branch("nuScores", &fnuScore);
  fTree->Branch("NeutrinoNuScores", &fNeutrinoNuScores);

  fSubRunTree = tfs->make<TTree> ("subRunTree", "SubRun Level Info TTree");
  fSubRunTree->Branch("run", &fRun_sr);
  fSubRunTree->Branch("subRun", &fSubRun_sr);
  fSubRunTree->Branch("pot", &fPOT);

  //Histograms
}

void hyperon::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.


}

void hyperon::AnalyzeEvents::endSubRun(art::SubRun const& sr) {
    fRun_sr = sr.run();
    fSubRun_sr = sr.subRun();
    fPOT = 0.0;

    art::Handle<sumdata::POTSummary> potHandle;
    if (sr.getByLabel("generator", potHandle)){
        fPOT = potHandle->totpot;
    }
    else {
        mf::LogVerbatim("hyperon::AnalyzeEvents") << "Warning: No POTSummary found in this SubRun!";
    }

    fSubRunTree->Fill();
}

DEFINE_ART_MODULE(hyperon::AnalyzeEvents)
