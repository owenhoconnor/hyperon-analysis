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
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbndcode/CosmicId/Algs/PandoraNuScoreCosmicIdAlg.h"

// // Root Includes
#include "TLorentzVector.h"
#include <iostream>
#include <vector>
#include <TTree.h>
#include <TH1.h>
#include <fstream>
#include <string>
#include <exception>

//Function to calculate distance between tracks
float ComputeTrackDistance(const art::Ptr<recob::Track>& track1, const art::Ptr<recob::Track>& track2)
{
    TVector3 start1 = track1->Vertex<TVector3>();
        TVector3 start2 = track2->Vertex<TVector3>();

            return (start1 - start2).Mag(); // Euclidean distance
            }
            

void DumpTrackInfo(const recob::Track& track) {
    std::cout << "Longest Track Reco DUMP INFO ****";
    std::cout << "Reco Track ID: " << track.ID() << "\n";
    std::cout << "Track Length: " << track.Length() << " cm\n";
    std::cout << "Start Point (Vertex): (" 
              << track.Vertex().X() << ", " 
              << track.Vertex().Y() << ", " 
              << track.Vertex().Z() << ")\n";
    std::cout << "End Point: (" 
              << track.End().X() << ", " 
              << track.End().Y() << ", " 
              << track.End().Z() << ")\n";
    std::cout << "Theta (Polar Angle): " << track.Theta() << " radians\n";
    std::cout << "Phi (Azimuthal Angle): " << track.Phi() << " radians\n";
    std::cout << "Zenith Angle: " << track.ZenithAngle() << " radians\n";
    std::cout << "Azimuth Angle: " << track.AzimuthAngle() << " radians\n";
}





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

private:
//sbnd::PandoraNuScoreCosmicIdAlg fCosmicIdAlg;
  // Declare member data here.

  int fVerbose;
 // std::vector<bool> fChildTrackIsLongest;
//  std::vector<float> fChildTrackLengths;
  std::string fHitLabel;
  std::string fGenieGenModuleLabel;
  std::string fSliceLabel;
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fCalorimetryLabel;
  int fNPrimaryParticles;
  int fNPrimaryChildren;
  std::vector<int> motherPDG;
  std::vector<int> daughterPDG;
  float vertexX, vertexY, vertexZ;
  std::vector<float> momentumX;
  std::vector<float> momentumY;
  std::vector<float> momentumZ;
//Calorimetry
  std::vector<std::vector<float>> fChildTrackdEdx;
  std::vector<std::vector<float>> fChildTrackResRange;
  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

// Neutrino Interaction variables
   int intType, CCNC, neutrinoPDG, numProtons, numNeutrons, numPi, numPi0, numTrueHits;
   float W, X, Y, QSqr, Pt, Theta, neutrinoE, leptonP;
   float trueVertexX, trueVertexY, trueVertexZ;
// Showers
   std::vector<int> fShowerPDG;
   std::string fShowerLabel; // Label for recob::Shower (pandoraShowerSBN)
   std::vector<float> fShowerEnergy; // Store shower energies
   std::vector<float> fShowerDirectionX; // Shower direction X
   std::vector<float> fShowerDirectionY; // Shower direction Y
   std::vector<float> fShowerDirectionZ; // Shower direction Z
   std::vector<double> fLength; 
 TTree *fTree;
  unsigned int fEventID;
 std::vector<float> trueP;
 std::vector<int> truePDG;
 std::string fOpT0FinderLabel;
 std::vector<double> ftime;
 std::vector<double> fscore;
 std::vector<float> fnuScore;
 std::vector<int> fNPfpSlices;
 std::vector<int> fNSlices;
 size_t totalSlices = 0;
 size_t totalNeutrinos = 0;
 int nuID = -1;
 int nuSliceKey = -1;
 std::vector<int> fnuSlicePDG;
 std::vector<std::vector<int>> fNeutrinoDaughterPDGs;
 std::vector<std::vector<int>> ftruthNeutrinoDaughterPDGs;
 std::string fClusterLabel;
 std::vector<int> fTrueTrackPDG;
 std::vector<int> fTrueShowerPDG;
// std::vector<bool> fChildTrackIsClosest;
// std::vector<float> fClosestTrackDistance;
// int longestID;
// float longestLength;
// int longestTrackPDG;
// int closestID;
// float minDistance;
// int closestTrackPDG;
// std::vector<float> fDistanceToRecoVertex;
 int fRun;       // Run number
 int fSubRun;    // Subrun number
// std::vector<double> fTrueTrackMomentum;
// std::vector<double> fTrueTrackEnergy;
// std::vector<double> fTrueTrackStartX;
// std::vector<double> fTrueTrackStartY;
// std::vector<double> fTrueTrackStartZ;
// std::vector<float> fLambdaDecayVertexX; // Decay vertex X position
// std::vector<float> fLambdaDecayVertexY; // Decay vertex Y position
// std::vector<float> fLambdaDecayVertexZ; // Decay vertex Z position
// std::vector<float> fProtonMomentum;    // Proton momentum
// std::vector<float> fPionMomentum;      // Negative pion momentum
// std::vector<float> fLambdaDecayDistance; // Distance from neutrino vertex to decay vertex
// std::vector<float> fLambdaInvariantMass; // Invariant mass of proton and pion
// int longestTrackMuonCount = 0;
// int longestTrackCount = 0;
// std::map<int, int> longestTrackPDGCounts;
 std::vector<int> fPFParticleIndex;
// std::vector<int> fTrackID;
 std::vector<std::string> fAlgorithmName;
 std::vector<int> fAssumedPdg;
// std::vector<float>fMuonScore;
// std::vector<float>fProtonScore;
// std::vector<float>fKaonScore;
// std::vector<float>fPionScore;
// std::vector<double> muonChi2, protonChi2, kaonChi2, pionChi2;

 std::map<int, int> trueTrackPDGMap;
 std::map<int, std::vector<std::pair<int, double>>> pidScoreMap;
 int fPlaneID;
 std::vector<int> fassumedPDGs;  
 std::vector<double> fchi2Scores; 
 std::vector<int> ftruePDGs;
 std::vector<int> frecoTruePDG; 

	 std::string fBadFilesOutput;  // Path to write problematic files to
  std::vector<std::string> fBadFiles;  // List of bad files encountered
//NEW VARIABLES
std::vector<int> fTrackIDs;
std::vector<float> fTrackLengths;
std::vector<int> fTrackPDGs;
std::vector<bool> fIsLongestTrack;
std::vector<bool> fIsClosestTrack;
std::vector<float> fDistanceToLongestTrack;
std::vector<float> fDistanceToRecoVertex;
std::vector<TVector3> fTrackStartPositions;
std::vector<TVector3> fTrackEndPositions;

std::vector<int> fTrueTrackPDGs;
std::vector<float> fTrueTrackEnergies;
std::vector<TVector3> fTrueTrackMomenta;
std::vector<TVector3> fTrueTrackStartPositions;
std::vector<float> fMatchQuality;

int fLongestTrackID;
float fLongestTrackLength;
int fLongestTrackPDG;
int fLongestTrackTruePDG;
float fLongestTrackDistanceToVertex;
int fClosestTrackID;
float fClosestTrackDistance;
int fClosestTrackPDG;
int fClosestTrackTruePDG;
TVector3 fRecoVertex;
bool fFoundRecoVertex;
bool foundLongestTrack;

std::map<int, int> longestTrackPDGCounts;
int longestTrackCount = 0;
int longestTrackMuonCount = 0;
std::map<int, int> closestTrackPDGCounts;
int closestTrackCount = 0;
int closestTrackKaonCount = 0;

std::vector<float> fMuonTrackScores;
std::vector<float> fKaonTrackScores;
std::vector<float> fProtonTrackScores;
std::vector<float> fPionTrackScores;
std::vector<float> fTrackScores;

std::vector<float> fAllPFPMuonScores;
std::vector<float> fAllPFPKaonScores;
std::vector<float> fAllPFPProtonScores;
std::vector<float> fAllPFPPionScores;
std::vector<float> fAllPFPElectronScores;

std::vector<float> fAllSliceT0s;
std::vector<int> fAllSliceTriggerTypes;

std::vector<double> fFlashMatchScores;
std::vector<double> fFlashPEs;

std::vector<float> fCosmicNuScores;  // nuScores for clear cosmics
std::vector<float> fNeutrinoNuScores; // nuScores for neutrinos
//
};


hyperon::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset}
    , fHitLabel(pset.get<std::string>("HitLabel"))
    , fGenieGenModuleLabel(pset.get<std::string>("GenieGenModuleLabel"))
    , fSliceLabel(pset.get<std::string>("SliceLabel")) // More initializers here.
    , fPFParticleLabel(pset.get<std::string>("PFParticleLabel"))
    , fTrackLabel(pset.get<std::string>("TrackLabel"))
    , fCalorimetryLabel(pset.get<std::string>("CalorimetryLabel"))
    , fShowerLabel(pset.get<std::string>("ShowerLabel"))
    , fOpT0FinderLabel(pset.get<std::string>("OpT0FinderLabel"))
    , fClusterLabel(pset.get<std::string>("ClusterLabel"))
    , fBadFilesOutput(pset.get<std::string>("BadFilesOutput", "bad_files.txt"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void hyperon::AnalyzeEvents::analyze(art::Event const& evt)
{
  // Implementation of required member function here.
 fEventID = evt.id().event(); 
 std::cout<<"Event# "<<evt.id().event()<<std::endl;

// NEW VARIABLES
fTrackIDs.clear();
    fTrackLengths.clear();
    fTrackPDGs.clear();
    fIsLongestTrack.clear();
    fIsClosestTrack.clear();
    fDistanceToLongestTrack.clear();
    fDistanceToRecoVertex.clear();
    fTrackStartPositions.clear();
    fTrackEndPositions.clear();
    fTrueTrackPDGs.clear();
    fTrueTrackEnergies.clear();
    fTrueTrackMomenta.clear();
    fTrueTrackStartPositions.clear();
    fMatchQuality.clear();
fTrackScores.clear();

fLongestTrackID = std::numeric_limits<int>::lowest();
    fLongestTrackLength = 0;
    fLongestTrackPDG = 0;
    fLongestTrackTruePDG = 0;
    fLongestTrackDistanceToVertex = -9999.0f;
    fClosestTrackID = -1;
    fClosestTrackDistance = std::numeric_limits<float>::max();
    fClosestTrackPDG = 0;
    fClosestTrackTruePDG = 0;
    fRecoVertex.SetXYZ(-9999., -9999., -9999.);
    fFoundRecoVertex = false;
	foundLongestTrack = false;
//





fRun = evt.run();
  fSubRun = evt.subRun();



   fnuScore.clear();
   fNPfpSlices.clear();
   fNSlices.clear();	
   fNPrimaryParticles = 0; 	
   fNPrimaryChildren = 0;
 //  fChildTrackLengths.clear();
 //  fChildTrackIsLongest.clear();
   fnuSlicePDG.clear();
	fAllSliceT0s.clear();
fAllSliceTriggerTypes.clear();
 fCosmicNuScores.clear();
    fNeutrinoNuScores.clear();
    // Get event slices
    art::ValidHandle<std::vector<recob::Slice>> sliceHandle = evt.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
    std::vector<art::Ptr<recob::Slice>> sliceVector;
    if (sliceHandle.isValid())
        art::fill_ptr_vector(sliceVector, sliceHandle);

    // Get associations between slices and PFParticles
    art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, evt, fSliceLabel);
	art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle =
  	evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
  
    // Get assocations between PFParticles and PFPMetaData
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfpMetadataAssoc(pfpHandle, evt, fPFParticleLabel);

    // Get MCTruth information
    art::ValidHandle<std::vector<simb::MCTruth>> mcTruthHandle = evt.getValidHandle<std::vector<simb::MCTruth>>("generator");
    std::vector<art::Ptr<simb::MCTruth>> mcTruthVector;
    if (mcTruthHandle.isValid())
        art::fill_ptr_vector(mcTruthVector, mcTruthHandle);
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
        art::fill_ptr_vector(mclist, mctruthListHandle);
    // simb::MCParticle to simb::MCTruth assocaition
    art::FindManyP< simb::MCParticle > fmpart( mctruthListHandle, evt, "largeant" );    

    // Get associations between PFParticles and Clusters, and between Hits and MCParticles
    art::FindManyP<recob::Cluster> pfpClusterAssoc(pfpHandle, evt, fPFParticleLabel);
    art::FindManyP<recob::Hit> clusterHitAssoc(evt.getValidHandle<std::vector<recob::Cluster>>(fClusterLabel), evt, fClusterLabel);
    art::FindManyP<simb::MCParticle> hitMCParticleAssoc(evt.getValidHandle<std::vector<recob::Hit>>(fHitLabel), evt, fHitLabel);

   // Get assocations between pfParticle and Vertex
   art::FindManyP<recob::Vertex> pfpVertexAssoc(pfpHandle, evt, fPFParticleLabel);

	std::cout << "Event " << fEventID << " has " << sliceVector.size() << " slices." << std::endl;
    fNSlices.push_back(sliceVector.size());

art::FindManyP<anab::T0> sliceToT0(sliceHandle, evt, "opt0finder");
art::FindManyP<sbn::SimpleFlashMatch> pfpToFlashMatch(pfpHandle, evt, "fmatch");

//	for (size_t i = 0; i < pfpHandle->size(); ++i) {
//    pfParticleMap[i] = art::Ptr<recob::PFParticle>(pfpHandle, i);
//}

    // Filling our neutrino hierarchy variables
    for (const art::Ptr<recob::Slice> &slice : sliceVector) {

	totalSlices++;

	if (slice.key() >= slicePFPAssoc.size()) {
    	std::cerr << "Error: Slice key " << slice.key() 
              << " is out of bounds for slicePFPAssoc (size = " 
              << slicePFPAssoc.size() << ")" << std::endl;
    	continue; // Skip this slice
	}

	  std::vector<art::Ptr<anab::T0>> t0vec = sliceToT0.at(slice.key());
	if (!t0vec.empty()) {
        const art::Ptr<anab::T0> t0 = t0vec.front();
        fAllSliceT0s.push_back(t0->Time());             //Store T0 time (in ns)
        fAllSliceTriggerTypes.push_back(t0->TriggerType()); //Store trigger type (cosmic, beam, etc.)
    } else {
        fAllSliceT0s.push_back(-999);   // Placeholder for missing T0
        fAllSliceTriggerTypes.push_back(-1);
    }



	std::vector<art::Ptr<recob::PFParticle>> slicePFPs(slicePFPAssoc.at(slice.key()));
       
//	auto pfpHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
//const auto& pfParticleVec = *pfpHandle;
	  if (slicePFPs.empty()) {
        std::cerr << "Warning: No PFParticles associated with slice key " 
                  << slice.key() << std::endl;
        	continue; // Skip this slice
    }
	
/*
	for (const art::Ptr<recob::PFParticle> &slicePFP : slicePFPs) {

    // Get the neutrino PFP ancestor of this slice
    recob::PFParticle pfpNu = fCosmicIdAlg.GetPFPNeutrino(slicePFP, pfParticleVec);

    // Get the NuScore of the neutrino PFP
    float nuScore = fCosmicIdAlg.GetPandoraNuScore(pfpNu, pfpMetadataAssoc);

  //   Determine if it's a clear cosmic
    bool isClearCosmic = fCosmicIdAlg.PandoraNuScoreCosmicId(pfpNu, pfParticleMap, evt);

    // Fill the overall NuScore vector
    fnuScore.push_back(nuScore);

//     Separate into cosmic or neutrino NuScore vectors
    if (isClearCosmic) {
        fCosmicNuScores.push_back(nuScore);
        std::cout << "Clear cosmic! NuScore: " << nuScore << std::endl;
    } else {
        fNeutrinoNuScores.push_back(nuScore);
    }
}

*/


 
	fNPfpSlices.push_back(slicePFPs.size());
	std::cout << "Slice key: " << slice.key() << ", Number of PFPs: " << slicePFPs.size() << std::endl;

// THROWS ERROR HERE
	for (const art::Ptr<recob::PFParticle> &slicePFP : slicePFPs) {

            const bool isPrimary = (slicePFP->IsPrimary());
    //        const bool isNeutrino ((std::abs(slicePFP->PdgCode()) == 12) || (std::abs(slicePFP->PdgCode()) == 14));

          if (!(isPrimary)) //&& isNeutrino))
                continue;



	  std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metadataVec = pfpMetadataAssoc.at(slicePFP.key());

          if (metadataVec.empty()) {
                   std::cerr << "No metadata found for PFParticle with key: " << slice.key() << std::endl;
                 continue;
                }

		float nuScore = -1.0;
		bool isClearCosmic = false;
		bool foundNuScore = false;
		bool foundClearCosmic = false;
       		
		for (const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata : metadataVec) {
            	const auto &propertiesMap = metadata->GetPropertiesMap();
            
		if (propertiesMap.find("NuScore") != propertiesMap.end()) {
                	nuScore = propertiesMap.at("NuScore");
        		 foundNuScore = true; 
	   	} else {
                	std::cerr << "NuScore not found for PFParticle with key: " << slice.key() << std::endl;
           	}


		if (propertiesMap.find("IsClearCosmic") != propertiesMap.end()) {
                isClearCosmic = static_cast<bool>(propertiesMap.at("IsClearCosmic"));
        	        foundClearCosmic = true;
    
		}

		}
       
		if (foundNuScore) {
    fnuScore.push_back(nuScore);

    if (isClearCosmic) {
        fCosmicNuScores.push_back(nuScore);
        std::cout << "Found clear cosmic with nuScore: " << nuScore << std::endl;
    } else {
        fNeutrinoNuScores.push_back(nuScore);
    }
} else {
    std::cerr << "Skipping PFParticle due to missing NuScore. Key: " << slicePFP.key() << std::endl;
}

if (!foundClearCosmic) {
    std::cerr << "Warning: IsClearCosmic property missing for PFParticle with key: " << slicePFP.key() << std::endl;
}
/*
	fnuScore.push_back(nuScore);  // Store all NuScores


	  if (isClearCosmic) {
            fCosmicNuScores.push_back(nuScore);
       		std::cout << "Found clear cosmic with nuScore: " << nuScore << std::endl;	
	 } else {
            fNeutrinoNuScores.push_back(nuScore);
        }
*/
	std::cout << "Neutrino slice detected! PDG codes of particles in this slice: ";
        for (const auto &particle : slicePFPs) {
            std::cout << particle->PdgCode() << " ";
	 }
        std::cout << std::endl;

		
            // We have found our neutrino!
            totalNeutrinos++;
            nuSliceKey = slice.key();
            nuID = slicePFP->Self();
            fNPrimaryParticles = slicePFPs.size();
            fNPrimaryChildren = slicePFP->NumDaughters();
	    fnuSlicePDG.push_back(slicePFP->PdgCode());
           
	}
	}




	std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs(slicePFPAssoc.at(nuSliceKey));

for (const art::Ptr<recob::PFParticle>& nuSlicePFP : nuSlicePFPs) {
        if (nuSlicePFP->IsPrimary() && std::abs(nuSlicePFP->PdgCode()) == 14) {
            auto vertices = pfpVertexAssoc.at(nuSlicePFP.key());
            
            if (!vertices.empty()) {
                const recob::Vertex& vertex = *vertices.at(0);
                auto const& vertexPos = vertex.position();
                fRecoVertex.SetXYZ(vertexPos.X(), vertexPos.Y(), vertexPos.Z());
                fFoundRecoVertex = true;
                break;
            }
        }
    }
    


	 if (!fFoundRecoVertex) {
        std::cerr << "Error: RecoVertex not found for the neutrino PFParticle!" << std::endl;
        fRecoVertex.SetXYZ(-9999., -9999., -9999.); // Use placeholder values
    }
    
    std::cout << "RecoVertex: (" << fRecoVertex.X() << ", " << fRecoVertex.Y() << ", " << fRecoVertex.Z() << ")" << std::endl;


	//Debug - should never return this because every slice will have a neutrino
//    if(nuID < 0){
//	std::cout << "No neutrino found in this event. Skipping to the next event." << std::endl;
//        return;

//	}
		//Sanity numbers check
	     //   std::cout << "Total number of slices: " << totalSlices << std::endl;
             //   std::cout << "Total number of neutrinos: " << totalNeutrinos << std::endl;



//Longest Track PID

art::FindManyP<recob::Track> pfpTrackAssoc(pfpHandle, evt, fTrackLabel);

/*
    longestID = std::numeric_limits<int>::lowest();
    longestLength = 0;
    longestTrackPDG = 0;
    closestID = std::numeric_limits<int>::lowest();
    minDistance = std::numeric_limits<float>::max();
    closestTrackPDG = 0;
    fChildTrackLengths.clear();
    fChildTrackIsLongest.clear();
    fChildTrackIsClosest.clear();
    fClosestTrackDistance.clear();
*/

art::Ptr<recob::Track> longestTrack;
TVector3 longestTrackStart; // Store the start point of the longest track
for (const art::Ptr<recob::PFParticle>& nuSlicePFP : nuSlicePFPs)
{
    // We are only interested in neutrino children particles
    if (nuSlicePFP->Parent() != static_cast<long unsigned int>(nuID))
        continue;

    // Get tracks associated with this PFParticle
    std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(nuSlicePFP.key());

    // There should only be 0 or 1 tracks associated with a PFP
    if (tracks.size() != 1)
        continue;

    // Get the track
    art::Ptr<recob::Track> track = tracks.at(0);
    float trackLength = track->Length();
 
	
	    auto const& trackStart = track->Vertex();
        TVector3 trackStartPos(trackStart.X(), trackStart.Y(), trackStart.Z());

	float distanceToVertex = -9999.0f;
        if (fFoundRecoVertex) {
            distanceToVertex = (trackStartPos - fRecoVertex).Mag();
        }

   // Check if this track is longer than the current longest
//	 longestID = track->ID();
//        longestLength = trackLength;
//        longestTrack = track;
//        longestTrackStart = track->Vertex<TVector3>();
//        longestTrackPDG = nuSlicePFP->PdgCode();       
 

	fTrackIDs.push_back(track->ID());
        fTrackLengths.push_back(trackLength);
        //fTrackPDGs.push_back(nuSlicePFP->PdgCode()); useless comes back as 13 for tracks and 11 for showers 
        fTrackStartPositions.push_back(track->Vertex<TVector3>());
        fTrackEndPositions.push_back(track->End<TVector3>());
	fDistanceToRecoVertex.push_back(distanceToVertex);
  
	//recob::DumpTracks trackDumper;
 //       std::cout << "\nNew longest track found:\n";
 //       DumpTrackInfo(*longestTrack);	

	if (trackLength > fLongestTrackLength) {
            fLongestTrackID = track->ID();
            fLongestTrackLength = trackLength;
            //fLongestTrackPDG = nuSlicePFP->PdgCode();
 	    fLongestTrackDistanceToVertex = distanceToVertex;
            longestTrack = track;
            longestTrackStart = track->Vertex<TVector3>();
		foundLongestTrack = true;


     }

}

//Finding closest track   
//int closestID = std::numeric_limits<int>::lowest();
//float minDistance = std::numeric_limits<float>::max();
//Now loop through the PFPs again to fill the track variables for the tree


if(foundLongestTrack && longestTrack){
for (const art::Ptr<recob::PFParticle> &nuSlicePFP : nuSlicePFPs)
{
    // We are only interested in neutrino children particles
    if (nuSlicePFP->Parent() != static_cast<long unsigned int>(nuID))
        continue;

    // Get tracks associated with this PFParticle
    std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(nuSlicePFP.key());

    // There should only be 0 or 1 tracks associated with a PFP
    if (tracks.size() != 1)
        continue;

    // Get the track
    art::Ptr<recob::Track> track = tracks.at(0);
    float distance = std::numeric_limits<float>::max();


 	  if (track->ID() == fLongestTrackID)
                continue;

    
        auto const& trackStart = track->Vertex();
        auto const& longestTrackStartVrt = longestTrack->Vertex();
        distance = (trackStart - longestTrackStartVrt).R();

	   
	if (distance < fClosestTrackDistance) {
                fClosestTrackDistance = distance;
                fClosestTrackID = track->ID();
                //fClosestTrackPDG = nuSlicePFP->PdgCode();
                
                std::cout << "\nNew closest track found:" << std::endl;
                std::cout << "Track ID: " << track->ID() 
                          << ", Distance: " << distance << " cm"
                          << std::endl;
            }

        }

}

	 if (fClosestTrackID != -1 && fClosestTrackDistance < 100) {
	size_t trackIndex = std::find(fTrackIDs.begin(), fTrackIDs.end(), fClosestTrackID) - fTrackIDs.begin();
        if (trackIndex < fTrackIDs.size()) {
            fDistanceToLongestTrack.push_back(fClosestTrackDistance);
        }
    }




/*

//Finding Distance to reconstruction vertex
fDistanceToRecoVertex.clear();
bool foundRecoVertex = false;

TVector3 recoVertex(-9999., -9999., -9999.);
for (const art::Ptr<recob::PFParticle> &nuSlicePFP : nuSlicePFPs)
{
    if (nuSlicePFP->IsPrimary() && std::abs(nuSlicePFP->PdgCode()) == 14)
    {
        auto vertices = pfpVertexAssoc.at(nuSlicePFP.key());
//        if (vertices.empty())
//	continue;
	
	if (!vertices.empty())
        {
            const recob::Vertex &vertex = *vertices.at(0);
            auto const& vertexPos = vertex.position();
            recoVertex.SetXYZ(vertexPos.X(), vertexPos.Y(), vertexPos.Z());
        foundRecoVertex = true;    
	break;
        }
    }
}

if (!foundRecoVertex)
{
    std::cerr << "Error: RecoVertex not found for the neutrino PFParticle!" << std::endl;
	   recoVertex.SetXYZ(-9999., -9999., -9999.); // Use placeholder values
//    fDistanceToRecoVertex.push_back(-9999.0f);
}


std::cout << "RecoVertex: (" << recoVertex.X() << ", " << recoVertex.Y() << ", " << recoVertex.Z() << ")" << std::endl;
float distanceToRecoVertex = -9999.0f;
if (longestTrack)
{
//    std::cerr << "Error: longestTrack is not valid!" << std::endl;
//	    longestLength = -1.0;
 //   longestTrackPDG = 0;
//    fDistanceToRecoVertex.push_back(-9999.0f);
 

auto const& longestTrackStartVrt = longestTrack->Vertex();
TVector3 longestTrackStartVrt_TV3(longestTrackStartVrt.X(), longestTrackStartVrt.Y(), longestTrackStartVrt.Z());
if(foundRecoVertex){
float distanceToRecoVertex = (longestTrackStartVrt_TV3 - recoVertex).Mag();
fDistanceToRecoVertex.push_back(distanceToRecoVertex);

}
}

std::cout << "Distance to RecoVertex: " << distanceToRecoVertex << " cm" << std::endl;
*/








// std::cout << "\nFinal longest track - ID: " << longestID 
//              << ", Length: " << longestLength 
//              << ", PDG: " << longestTrackPDG << std::endl;
    





//Tracks reco->Truth
frecoTruePDG.clear();
fTrueTrackPDG.clear();
fMuonTrackScores.clear();
fKaonTrackScores.clear();
fProtonTrackScores.clear();
fPionTrackScores.clear();
std::map<int, TVector3> trackTrueVertices;
art::ValidHandle<std::vector<recob::Track>> trackHandle = 
    evt.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
if (!trackHandle.isValid()) {
    std::cerr << "WARNING: recob::Track collection not found in event." << std::endl;
    return;
}

art::Handle<art::Assns<simb::MCParticle, recob::Hit>> hitTruthAssn;
evt.getByLabel("gaushitTruthMatch", hitTruthAssn);
art::FindManyP<recob::Hit> trackHitAssoc(trackHandle, evt, fTrackLabel);
 art::FindManyP<recob::PFParticle> trackToPFPAssoc(trackHandle, evt, fTrackLabel);

for (const art::Ptr<recob::PFParticle>& slicePFP : nuSlicePFPs) {
    art::FindManyP<recob::Track> pfpTrackAssoc(pfpHandle, evt, fTrackLabel);

    std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(slicePFP.key());
    if (tracks.empty()) 
	continue;
   
    if (tracks.size() != 1)
        continue;
//art::Ptr<recob::Track> track = tracks[0];

for (const auto& track : tracks) {
 
//    art::Ptr<recob::Track> track = tracks[0];
  
    std::vector<art::Ptr<recob::PFParticle>> trackPFPs = trackToPFPAssoc.at(track.key());
	art::Ptr<recob::PFParticle> trackPFP = trackPFPs.front();
	 std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> trackMetadataVec = pfpMetadataAssoc.at(trackPFP.key());
    float trackScore = -1.0;
    if (trackMetadataVec.empty()) {
	std::cerr << "No metadata found for track with ID: " << track->ID() << std::endl;
       
                }

        for (const auto& metadata : trackMetadataVec) {
            const auto& propertiesMap = metadata->GetPropertiesMap();
      		  std::cout << "\nTrack PFParticle Metadata Properties (ID " << track->ID() << "):" << std::endl;	 
	 for (const auto& [propertyName, propertyValue] : propertiesMap) {
            std::cout << "  " << propertyName << " = " << propertyValue << std::endl;
		}
	 if (propertiesMap.find("TrackScore") != propertiesMap.end()) {
                        trackScore = propertiesMap.at("TrackScore");
	 
    } else {
        std::cerr << "No metadata found for track PFParticle ID: " << track->ID() << std::endl;
    }
	          fTrackScores.push_back(trackScore);

	}



//}







  
    std::vector<art::Ptr<recob::Hit>> trackHits = trackHitAssoc.at(track.key());
    

    art::FindManyP<simb::MCParticle> hitToMCParticleAssn(trackHits, evt, "gaushitTruthMatch");

    std::map<int, int> mcParticleHitCount;
    
    for (size_t i_hit = 0; i_hit < trackHits.size(); i_hit++) {
        std::vector<art::Ptr<simb::MCParticle>> mcParticles = hitToMCParticleAssn.at(i_hit);
        for (const auto& mcParticle : mcParticles) {
            mcParticleHitCount[mcParticle->TrackId()]++;
        }
    }
   
    int bestTrackID = -1;
    int maxHitCount = 0;
    for (const auto& [trackID, hitCount] : mcParticleHitCount) {
        if (hitCount > maxHitCount) {
            maxHitCount = hitCount;
            bestTrackID = trackID;
        }
    }
    
    if (bestTrackID != -1) {
	art::ServiceHandle<cheat::ParticleInventoryService> piService;
        const simb::MCParticle* mcParticle = piService->TrackIdToParticle_P(bestTrackID);
	if (mcParticle) {
	trueTrackPDGMap[track->ID()] = mcParticle->PdgCode();	

	 int pdgCode = mcParticle->PdgCode();
            if (pdgCode == 13) {
                fMuonTrackScores.push_back(trackScore);
            } else if (pdgCode == 321) {
                fKaonTrackScores.push_back(trackScore);
            } else if (pdgCode == 2212) {
                fProtonTrackScores.push_back(trackScore);
            } else if (pdgCode == 211 || pdgCode == -211) {
                fPionTrackScores.push_back(trackScore);
            }



		float matchQuality = static_cast<float>(maxHitCount) / trackHits.size();

		size_t trackIndex = std::find(fTrackIDs.begin(), fTrackIDs.end(), track->ID()) - fTrackIDs.begin();
                if (trackIndex < fTrackIDs.size()) {
		fTrueTrackPDGs.push_back(mcParticle->PdgCode());

		 fIsLongestTrack.push_back(track->ID() == fLongestTrackID);
                    fIsClosestTrack.push_back(track->ID() == fClosestTrackID);

		 if (track->ID() == fLongestTrackID) {
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

			 if (track->ID() == fClosestTrackID) {
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

				}
			}
		}
	}
}

/*
//		 frecoTruePDG[track->ID()] = mcParticle->PdgCode();
		 TLorentzVector momentum = mcParticle->Momentum();
                    double px = momentum.Px();
                    double py = momentum.Py();
                    double pz = momentum.Pz();
                    double energy = momentum.E();
		 //   double p = momentum.P();
		 trackTrueVertices[track->ID()] = TVector3(mcParticle->Vx(), mcParticle->Vy(), mcParticle->Vz());
	       	 fTrueTrackPDG.push_back(mcParticle->PdgCode());
		 //fTrueTrackMomentum.push_back(p);        
            	 //fTrueTrackEnergy.push_back(energy);           
             	 fTrueTrackStartX.push_back(mcParticle->Vx());  
            	 fTrueTrackStartY.push_back(mcParticle->Vy());   
            	 fTrueTrackStartZ.push_back(mcParticle->Vz());   
			
		 std::cout<< "=========================="<<"\n" << std::endl;
		 std::cout<< "Best Matched MC Particle:"
			      << "Reco Track ID: " << track->ID()
                              << ", PDG: " << mcParticle->PdgCode()
                              << ", MCParticle Track ID: " << mcParticle->TrackId()
                              << ", Hit Count: " << maxHitCount
                              << ", Momentum: (" << px << ", " << py << ", " << pz << ")"
                              << ", Energy: " << energy << " GeV"
                              << std::endl;

		 if (isLongestTrack) {
                        std::cout <<"************************"<<"\n" << std::endl;
			std::cout << "Best matched MC Particle for Longest Track: "
                                  << "MCParticle Track ID: " << mcParticle->TrackId()
				  << ", PDG: " << mcParticle->PdgCode()
                                  << ", Hit Count: " << maxHitCount
                                  << ", Momentum: (" << px << ", " << py << ", " << pz << ")"
                                  << ", Energy: " << energy << " GeV"
                                  << std::endl;
			 longestTrackCount++;
			 longestTrackPDGCounts[mcParticle->PdgCode()]++;
			 if (mcParticle->PdgCode() == 13) {
               	 	 longestTrackMuonCount++; // Increment the counter
            		}
                    }
		} else {
        std::cerr << "ERROR: Could not find MCParticle for track ID " << bestTrackID << std::endl;
                }
*/           


  std::cout << "\nFinal longest track - ID: " << fLongestTrackID 
              << ", Length: " << fLongestTrackLength 
              << ", PDG: " << fLongestTrackPDG
              << ", True PDG: " << fLongestTrackTruePDG
              << ", Distance to Vertex: " << fLongestTrackDistanceToVertex << " cm" << std::endl;
              
    std::cout << "\nFinal closest track - ID: " << fClosestTrackID 
              << ", Distance: " << fClosestTrackDistance << " cm" 
              << ", PDG: " << fClosestTrackPDG
              << ", True PDG: " << fClosestTrackTruePDG << std::endl;

/*


std::cout << "Number of times the longest track is a muon (PDG ID 13): " << longestTrackMuonCount << "Number of longest Track : "  << longestTrackCount<< std::endl;
double fractionMuon = static_cast<double>(longestTrackMuonCount) / longestTrackCount;
std::cout << "Fraction of times the longest track is a muon: " << fractionMuon << std::endl;

std::cout << "\n==============================" << std::endl;
std::cout << "Accumulated counts of longest track PDG codes:" << std::endl;
for (const auto& [pdg, count] : longestTrackPDGCounts) {
    std::cout << "PDG Code: " << pdg << ", Count: " << count << std::endl;
}
std::cout << "==============================" << std::endl;

*/

//Showers reco->Truth
fTrueShowerPDG.clear();

art::ValidHandle<std::vector<recob::Shower>> showerHandle = 
	evt.getValidHandle<std::vector<recob::Shower>>(fShowerLabel);

art::FindManyP<recob::Hit> showerHitAssoc(showerHandle, evt, fShowerLabel);
for (const art::Ptr<recob::PFParticle>& slicePFP : nuSlicePFPs) {
    art::FindManyP<recob::Shower> pfpShowerAssoc(pfpHandle, evt, fShowerLabel);

    std::vector<art::Ptr<recob::Shower>> showers = pfpShowerAssoc.at(slicePFP.key());

    if (showers.empty())
        continue;

    if (showers.size() != 1)
        continue;


    art::Ptr<recob::Shower> shower = showers[0];

    std::vector<art::Ptr<recob::Hit>> showerHits = showerHitAssoc.at(shower.key());
	art::FindManyP<simb::MCParticle> showerHitToMCParticleAssn(showerHits, evt, "gaushitTruthMatch");
    
    std::map<int, int> mcParticleHitCount;
    
    
for (size_t i_hit = 0; i_hit < showerHits.size(); i_hit++) {
    std::vector<art::Ptr<simb::MCParticle>> mcParticles = showerHitToMCParticleAssn.at(i_hit); 

 for (const auto& mcParticle : mcParticles) {
        mcParticleHitCount[mcParticle->TrackId()]++;
    }
}
    int bestTrackID = -1;
    int maxHitCount = 0;
    for (const auto& [trackID, hitCount] : mcParticleHitCount) {
        if (hitCount > maxHitCount) {
            maxHitCount = hitCount;
            bestTrackID = trackID;
        }
    }
    if (bestTrackID != -1) {
           art::ServiceHandle<cheat::ParticleInventoryService> piService;
    const simb::MCParticle* mcParticle = piService->TrackIdToParticle_P(bestTrackID);    

	if (mcParticle) {
                //frecoTruePDG[shower->ID()] = mcParticle->PdgCode();
		fTrueShowerPDG.push_back(mcParticle->PdgCode());
                std::cout << "Best matched MC Particle for Shower: "
                          << "PDG: " << mcParticle->PdgCode()
                          << ", Track ID: " << mcParticle->TrackId()
                          << ", Hit Count: " << maxHitCount << std::endl;
            }
        }
}

/// TESTING
    fAllPFPMuonScores.clear();
    fAllPFPKaonScores.clear();
    fAllPFPProtonScores.clear();
    fAllPFPPionScores.clear();
    fAllPFPElectronScores.clear();
	fFlashMatchScores.clear();
fFlashPEs.clear();


	art::FindManyP<recob::SpacePoint> pfpToSP(pfpHandle, evt, fPFParticleLabel);
art::Handle<std::vector<recob::SpacePoint>> spHandle;
evt.getByLabel(fPFParticleLabel, spHandle);
art::FindManyP<recob::Hit> spToHit(spHandle, evt, fPFParticleLabel);

    for (size_t i = 0; i < pfpHandle->size(); ++i) {
        art::Ptr<recob::PFParticle> pfp(pfpHandle, i);
    


std::vector<art::Ptr<sbn::SimpleFlashMatch>> flashMatches = pfpToFlashMatch.at(pfp.key());
	std::cout << "PFParticle key: " << pfp.key() 
          << ", nFlashMatches: " << flashMatches.size() << std::endl;

	
        if (!flashMatches.empty()) {
        const sbn::SimpleFlashMatch& match = *flashMatches.front();  // Use first match

	std::cout << "  present: " << match.present 
              << ", score: " << match.score.total 
              << ", PE: " << match.light.pe << std::endl;

        if (match.present) {
            double totalScore = match.score.total;
            double flashPE = match.light.pe;

            fFlashMatchScores.push_back(totalScore);
            fFlashPEs.push_back(flashPE);
        } else {
       //      No valid match present
            fFlashMatchScores.push_back(-999);
            fFlashPEs.push_back(-999);
        }
    }



    if (pfp->IsPrimary() && abs(pfp->PdgCode()) == 14) continue;

        float trackScore = -1.0;
	std::vector<art::Ptr<recob::Hit>> pfpHits;
	std::vector<art::Ptr<recob::SpacePoint>> spVec = pfpToSP.at(pfp.key());
for (const auto& sp : spVec) {
    std::vector<art::Ptr<recob::Hit>> hits = spToHit.at(sp.key());
    pfpHits.insert(pfpHits.end(), hits.begin(), hits.end());
}

        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metadataVec = pfpMetadataAssoc.at(pfp.key());
        if (metadataVec.empty()) {
            std::cerr << "No metadata found for PFP with ID: " << pfp->Self() << std::endl;
            continue;
        }

        for (const auto& metadata : metadataVec) {
            const auto& propertiesMap = metadata->GetPropertiesMap();
            if (propertiesMap.find("TrackScore") != propertiesMap.end()) {
                trackScore = propertiesMap.at("TrackScore");
                break;
            }
        }

        if (trackScore < 0) continue;

        if (pfpHits.empty()) continue;

        art::FindManyP<simb::MCParticle> hitToMCParticleAssn(pfpHits, evt, "gaushitTruthMatch");
        std::map<int, int> mcParticleHitCount;
	for (size_t i_hit = 0; i_hit < pfpHits.size(); i_hit++) {
            std::vector<art::Ptr<simb::MCParticle>> mcParticles = hitToMCParticleAssn.at(i_hit);
            for (const auto& mcParticle : mcParticles) {
		     mcParticleHitCount[mcParticle->TrackId()]++;
            }
        }

        int bestTrackID = -1;
        int maxHitCount = 0;
        for (const auto& [trackID, hitCount] : mcParticleHitCount) {
            if (hitCount > maxHitCount) {
                maxHitCount = hitCount;
                bestTrackID = trackID;
            }
        }

        if (bestTrackID == -1) continue;

        art::ServiceHandle<cheat::ParticleInventoryService> piService;
        const simb::MCParticle* mcParticle = piService->TrackIdToParticle_P(bestTrackID);
        if (!mcParticle) continue;

	int correctedHitCount = 0;
for (size_t i_hit = 0; i_hit < pfpHits.size(); i_hit++) {
    std::vector<art::Ptr<simb::MCParticle>> mcParticles = hitToMCParticleAssn.at(i_hit);
    for (const auto& mc : mcParticles) {
        if (mc->TrackId() == bestTrackID) {
            correctedHitCount++;
            break; // Only count once per hit
        }
    }
}



        float matchQuality = static_cast<float>(correctedHitCount) / pfpHits.size();
        int pdgCode = mcParticle->PdgCode();
	fMatchQuality.push_back(matchQuality);

       // if (matchQuality > 0.5) {
            if (pdgCode == 13) {
                fAllPFPMuonScores.push_back(trackScore);
            } else if (pdgCode == 321) {
                fAllPFPKaonScores.push_back(trackScore);
            } else if (pdgCode == 2212) {
                fAllPFPProtonScores.push_back(trackScore);
            } else if (pdgCode == -211) {
                fAllPFPPionScores.push_back(trackScore);
            } else if (pdgCode == 11) {
                fAllPFPElectronScores.push_back(trackScore);
            }


            std::cout << "PFP ID: " << pfp->Self() 
                      << ", Track Score: " << trackScore 
                      << ", True PDG: " << pdgCode 
                      << ", Match Quality: " << matchQuality 
                      << std::endl;
       // }
    }



//



//Calorimetry
/*
art::FindManyP<anab::Calorimetry> trackCaloAssoc(trackHandle, evt, fCalorimetryLabel);
for(const auto& nuSlicePFP : nuSlicePFPs) 
{
	if (nuSlicePFP->Parent() != static_cast<long unsigned int>(nuID))
		continue;

	std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(nuSlicePFP.key());
	if (tracks.size() != 1) continue;

	art::Ptr<recob::Track> track = tracks.at(0);


//Get the Calorimetry object
std::vector<art::Ptr<anab::Calorimetry>> calos = trackCaloAssoc.at(track.key());

for(auto const& calo : calos)
	{
	const int plane = calo->PlaneID().Plane;
	
	//collection plane (2)
	if(plane !=2)
		continue;
	fChildTrackdEdx.push_back(calo->dEdx());
	fChildTrackResRange.push_back(calo->ResidualRange());
	}
}
*/
//Distance to reco vertex for Lambda daughters 
/*
fLambdaDecayVertexX.clear();
fLambdaDecayVertexY.clear();
fLambdaDecayVertexZ.clear();
fLambdaDecayDistance.clear();
fProtonMomentum.clear();
fPionMomentum.clear();
fLambdaInvariantMass.clear();

for (size_t i_truth = 0; i_truth < mclist.size(); i_truth++) {
   art::Ptr<simb::MCTruth> truth(mclist.at(i_truth));
   std::vector<art::Ptr<simb::MCParticle>> assocParticles(fmpart.at(i_truth));
   for (size_t i_mcpart = 0; i_mcpart < assocParticles.size(); i_mcpart++) {
       art::Ptr<simb::MCParticle> mcParticle(assocParticles.at(i_mcpart));
       if (mcParticle->PdgCode() == 3122) {
           for (int i_daughter = 0; i_daughter < mcParticle->NumberDaughters(); ++i_daughter) {
               int daughterTrackID = mcParticle->Daughter(i_daughter);
               for (size_t i = 0; i < assocParticles.size(); ++i) {
                   if (assocParticles.at(i)->TrackId() == daughterTrackID) {
                       int daughterPdgCode = assocParticles.at(i)->PdgCode();
                       if (daughterPdgCode == 2212 || daughterPdgCode == -211) {
                           TVector3 decayVertex(
                               assocParticles.at(i)->Vx(),
                               assocParticles.at(i)->Vy(),
                               assocParticles.at(i)->Vz()
                           );
                           TVector3 neutrinoVertex(vertexX, vertexY, vertexZ);
                           float distance = (decayVertex - neutrinoVertex).Mag();
                           
                           fLambdaDecayVertexX.push_back(decayVertex.X());
                           fLambdaDecayVertexY.push_back(decayVertex.Y());
                           fLambdaDecayVertexZ.push_back(decayVertex.Z());
                           fLambdaDecayDistance.push_back(distance);
                           
                           TLorentzVector daughterMomentum = assocParticles.at(i)->Momentum();
                           if (daughterPdgCode == 2212) {
                               fProtonMomentum.push_back(daughterMomentum.P());
                           } else if (daughterPdgCode == -211) {
                               fPionMomentum.push_back(daughterMomentum.P());
                           }
                       }
                   }
               }
           }
       }
   }
}
*/
//PID anab
art::FindManyP<anab::ParticleID> trackParticleIDAssoc(trackHandle, evt, "pandoraPid");

fPFParticleIndex.clear();
//fTrackID.clear();
fAlgorithmName.clear();
fAssumedPdg.clear();
//fMuonScore.clear();
//fProtonScore.clear();
//fKaonScore.clear();
//fPionScore.clear();



for (size_t i = 0; i < pfpHandle->size(); ++i)
{
    art::Ptr<recob::PFParticle> pfParticle(pfpHandle, i);
 

    auto tracks = pfpTrackAssoc.at(pfParticle.key());
    if (!tracks.empty())
    {
        for (const art::Ptr<recob::Track> &track : tracks)
        {
	     int trackID = track->ID();

              
                auto it = trueTrackPDGMap.find(trackID);
                if (it == trueTrackPDGMap.end()) continue;

            auto particleIDs = trackParticleIDAssoc.at(track.key());
            if (!particleIDs.empty())
            {
			 pidScoreMap[trackID].clear();

			std::cout << "[DEBUG] Processing track " << trackID << " with true PDG " 
                          << it->second << std::endl;
		
		           for (auto it_pid = particleIDs.rbegin(); it_pid != particleIDs.rend(); ++it_pid) {
                    const art::Ptr<anab::ParticleID>& particleID = *it_pid;
                    const int planeID = particleID->PlaneID().Plane;
		if(planeID !=2){
                continue;
		}
		

		std::map<int, double> tempScores;  // PDG -> Score
    		bool allScoresValid = true;


                    const std::vector<anab::sParticleIDAlgScores> &algScores = particleID->ParticleIDAlgScores();
                	 std::cout << "[DEBUG] Found " << algScores.size() 
                              << " algorithm scores for track " << trackID 
                              << " on plane " << planeID << std::endl;  


		 for (const anab::sParticleIDAlgScores &score : algScores)
                    {  	
			std::cout << "[DEBUG] Algorithm Name: " << score.fAlgName 
              << ", Assumed PDG: " << score.fAssumedPdg 
              << ", Score Value: " << score.fValue << std::endl;	

			if(score.fAlgName == "Chi2" && score.fValue != 0.0) {
            tempScores[score.fAssumedPdg] = score.fValue;
      		  }
		    }
			const std::set<int> requiredPDGs = {13, 2212, 321, 211};  // mu, proton, kaon, pion
    			allScoresValid = true;

			for (int pdg : requiredPDGs) {
        if(tempScores.find(pdg) == tempScores.end() || tempScores[pdg] == 0.0) {
            allScoresValid = false;
            break;
        }
    }


	if(allScoresValid) {
        pidScoreMap[trackID].clear();


		        for (const auto& [pdg, score] : tempScores) {
            pidScoreMap[trackID].emplace_back(pdg, score);
            std::cout << "[DEBUG] Storing PDG " << pdg 
                      << " with score " << score 
                      << " for track " << trackID << std::endl;


		}
		break;
		} else {

			  std::cout << "[WARNING] Incomplete/Invalid Chi2 set found for track " 
                  << trackID << " in ParticleID " << particleID.key()
                  << ". Skipping to earlier ParticleID." << std::endl;
    }
}
}
}                       
}     
}
		

	



size_t nTracks = pidScoreMap.size();
std::cout << "[DEBUG] Number of tracks with PID scores: " << nTracks << std::endl;



std::map<int, int> truePDGCounts;
for (const auto& entry : trueTrackPDGMap) {
    truePDGCounts[entry.second]++;
}

std::cout << "\nTrue PDG Distribution:\n";
std::cout << "----------------------\n";
for (const auto& pdgCount : truePDGCounts) {
    std::cout << "PDG " << pdgCount.first << ": " << pdgCount.second << " tracks\n";
}
std::cout << "----------------------\n";



size_t totalScores = 0;
fassumedPDGs.clear();
fchi2Scores.clear();
ftruePDGs.clear();

for (const auto& pair : pidScoreMap) {
    int trackID = pair.first;
    int truePDG = trueTrackPDGMap[trackID]; 
    size_t nScores = pair.second.size();
    totalScores += nScores;

    
    for (const auto& scorePair : pair.second) {
        ftruePDGs.push_back(truePDG);
        fassumedPDGs.push_back(scorePair.first);
        fchi2Scores.push_back(scorePair.second);
    }
}

std::cout << "[DEBUG] Total PID scores across all tracks: " << totalScores << std::endl;


std::cout << "Vector sizes: truePDGs=" << ftruePDGs.size()
          << ", assumedPDGs=" << fassumedPDGs.size()
          << ", chi2Scores=" << fchi2Scores.size() << std::endl;




/*
			 fPFParticleIndex.push_back(i);
                        fTrackID.push_back(track->ID());
                        fAlgorithmName.push_back(score.fAlgName);
          		 int assumedPdg = score.fAssumedPdg;         

	if (score.fAlgName == "Chi2") {
                            if (assumedPdg == 13) {  // Muon
                                fMuonScore.push_back(score.fValue);
                            } else if (assumedPdg == 2212) {  // Proton
                                fProtonScore.push_back(score.fValue);
                            } else if (assumedPdg == 211) {  // Pion
                                fPionScore.push_back(score.fValue);
                            } else if (assumedPdg == 321) {  // Kaon
                                fKaonScore.push_back(score.fValue);
                            }
                        }
	std::cout << "PFParticle " << i << ", Track ID: " << track->ID()
          << ", Algorithm: " << score.fAlgName
          << ", Assumed PDG: " << assumedPdg
          << ", Score: " << score.fValue << std::endl;
		    }
                }
            }
            else
            {
                std::cout << "No ParticleID associated with Track " << track->ID() << std::endl;
            }
        }
    }
    else
    {
        std::cout << "No Tracks associated with PFParticle " << i << 
*/



//####
trueP.clear();
truePDG.clear();
motherPDG.clear();
daughterPDG.clear();
momentumX.clear();
momentumY.clear();
momentumZ.clear();
vertexX = 0;
vertexY = 0;
vertexZ = 0; 



// * MC truth information
   for (size_t i_truth = 0; i_truth < mclist.size(); i_truth++)
   {
       art::Ptr<simb::MCTruth> truth(mclist.at(i_truth));
//       std::cout<<"Neutrino at index " << i_truth << " has pdg: " << truth->GetNeutrino().Nu().PdgCode() << std::endl;


       std::vector< art::Ptr<simb::MCParticle> > assocParticles(fmpart.at(i_truth));
       for (size_t i_mcpart = 0; i_mcpart < assocParticles.size(); i_mcpart++)
       {
           art::Ptr<simb::MCParticle> mcParticle(assocParticles.at(i_mcpart));
           if(mcParticle->Mother()!=10000000) continue;

//           std::cout<<"--Particle at index " << i_mcpart << " has pdg: " << mcParticle->PdgCode() 
//		    <<" has momentum: " << mcParticle->P()<<" MotherID: "<<mcParticle->Mother()
//		    <<" has status: "<< mcParticle->StatusCode()<<" number trajectory points: "<< mcParticle->NumberTrajectoryPoints()
//		    <<" TrackId: " << mcParticle->TrackId()<< std::endl;
	//FOR FINAL STATE PARTICLES
           trueP.push_back(mcParticle->P()); 

	   truePDG.push_back(mcParticle->PdgCode());

 	   // New Variables 21.10  /Fill new Variables 
	   motherPDG.push_back(mcParticle->Mother());
	   if (mcParticle->NumberTrajectoryPoints()>0) {
		auto const& startPoint = mcParticle->Position(0);
		vertexX = startPoint.X();
		vertexY = startPoint.Y();
		vertexZ = startPoint.Z();
		}
		else {
            // If no trajectory points, store placeholders for consistency
                vertexX = -9999.; 
                vertexY = -9999.;
                vertexZ = -9999.;
                      }	   
	   momentumX.push_back(mcParticle->Px());
	   momentumY.push_back(mcParticle->Py());
  	   momentumZ.push_back(mcParticle->Pz());
		
	// Daughter		
        for (int i_daughter = 0; i_daughter < mcParticle->NumberDaughters(); ++i_daughter) {
            int daughterTrackID = mcParticle->Daughter(i_daughter);
	for (size_t i = 0; i < assocParticles.size(); ++i) {
                if (assocParticles.at(i)->TrackId() == daughterTrackID) {
                    int daughterPdgCode = assocParticles.at(i)->PdgCode();
                    daughterPDG.push_back(daughterPdgCode);
//                    std::cout << "---- Daughter Particle: PDG Code: " << daughterPdgCode 
//                              << ", Track ID: " << daughterTrackID 
//                              << std::endl;
                    break;
                }
            }
        }

      

}


   }


//Daughter truth matching
ftruthNeutrinoDaughterPDGs.clear();

//for (const art::Ptr<recob::PFParticle> &nuSlicePFP : nuSlicePFPs) {
	std::vector<int> truthMatchedDaughterPDGs;
    
    for (size_t i_truth = 0; i_truth < mclist.size(); i_truth++) {
        art::Ptr<simb::MCTruth> truth(mclist.at(i_truth));
	// if (std::abs(truth->GetNeutrino().Nu().PdgCode()) == std::abs(nuSlicePFP->PdgCode())) {
            std::vector<art::Ptr<simb::MCParticle>> assocParticles(fmpart.at(i_truth));
                       int neutrinoTrackId = -1;
            for (const auto& mcParticle : assocParticles) {
                if (mcParticle->PdgCode() == truth->GetNeutrino().Nu().PdgCode()) {
                    neutrinoTrackId = mcParticle->TrackId();
                    break;
                }
            }
 	 if (neutrinoTrackId != -1) {
                for (const auto& mcParticle : assocParticles) {
		if (mcParticle->Mother() == neutrinoTrackId) {
                        truthMatchedDaughterPDGs.push_back(mcParticle->PdgCode());
                    }
                }
            }
       // }
    }
//	if (!truthMatchedDaughterPDGs.empty()) {
        ftruthNeutrinoDaughterPDGs.push_back(truthMatchedDaughterPDGs);
   // }
//}    




/* THIS WONT WORK
// OpT0Finder
 //Define associaiton
 art::ValidHandle<std::vector<sbn::OpT0Finder>> opt0FinderHandle = 
    evt.getValidHandle<std::vector<sbn::OpT0Finder>>(fOpT0FinderLabel);

 art::FindManyP<sbn::OpT0Finder> sliceOpT0Assoc(sliceHandle, evt, fOpT0FinderLabel);
 
 for (const art::Ptr<recob::PFParticle>& nuSlicePFP : nuSlicePFPs)
{
    // We are only interested in neutrino children particles
         if (nuSlicePFP->Parent() != static_cast<long unsigned int>(nuID))
                 continue;
                  
	 // store associaitons  
	 std::vector<art::Ptr<sbn::OpT0Finder>> opT0s = sliceOpT0Assoc.at(nuSliceKey);
	//Loop over opT0s we just found
	for (const art::Ptr<sbn::OpT0Finder>& opT0 : opT0s)
	{
	
	//Store time and score
	fscore.push_back(opT0->score);
        ftime.push_back(opT0->time);	
                }
}
*/

// Track this change
// clearly that didnt work
// bosh 
std::cout<<"***************************"<<std::endl;

// NOT MY CODE
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


 
 fTree->Fill();
}


void hyperon::AnalyzeEvents::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree> ("tree", "Output TTree");
  fTree ->Branch("nPFParticle", &fNPrimaryParticles);
  fTree ->Branch("nPrimaryChildren", &fNPrimaryChildren);
 // fTree ->Branch("childTrackLengths", &fChildTrackLengths);
 // fTree ->Branch("childTrackdEdx", &fChildTrackdEdx);
 // fTree ->Branch("childTrackResRange", &fChildTrackResRange);
//  fTree ->Branch("childTrackIsLongest", &fChildTrackIsLongest); 
  fTree ->Branch("eventID", &fEventID);
  fTree ->Branch("trueP", &trueP);
  fTree ->Branch("truePDG", &truePDG);
//  fTree ->Branch("NeutrinoDaughterPDGs", &fNeutrinoDaughterPDGs);
//  fTree ->Branch("truthNeutrinoDaughterPDGs", &ftruthNeutrinoDaughterPDGs);
  fTree ->Branch("motherPDG", &motherPDG);
  fTree ->Branch("daughterPDG", &daughterPDG);
  fTree ->Branch("vertexX", &vertexX);
  fTree ->Branch("vertexY", &vertexY);
  fTree ->Branch("vertexZ", &vertexZ);
  fTree ->Branch("momentumX", &momentumX);
  fTree ->Branch("momentumY", &momentumY);
  fTree ->Branch("momentumZ", &momentumZ); 
  fTree ->Branch("nuID", &nuID);
//  fTree ->Branch("Length", &fLength);
  //fTree ->Branch("ShowerEnergy", &fShowerEnergy);
  //fTree ->Branch("ShowerDirectionX", &fShowerDirectionX);
  //fTree ->Branch("ShowerDirectionY", &fShowerDirectionY);
  //fTree ->Branch("ShowerDirectionZ", &fShowerDirectionZ); 
//  fTree ->Branch("ShowerPDG", &fShowerPDG);
  fTree ->Branch("nuScore", &fnuScore);
//  fTree ->Branch("NPfpSlices", &fNPfpSlices);
//  fTree ->Branch("NSlices", &fNSlices);
//  fTree ->Branch("nuSlicePDG", &fnuSlicePDG);
  fTree ->Branch("TrueTrackPDG", &fTrueTrackPDG);
  fTree ->Branch("TrueShowerPDG", &fTrueShowerPDG);
//  fTree ->Branch("longestID", &longestID);
 // fTree ->Branch("longestLength", &longestLength);
//  fTree ->Branch("longestTrackPDG", &longestTrackPDG);
//  fTree ->Branch("closestID", &closestID);
 // fTree ->Branch("minDistance", &minDistance);
//  fTree ->Branch("closestTrackPDG", &closestTrackPDG);
//  fTree ->Branch("ChildTrackIsClosest", &fChildTrackIsClosest);
//  fTree ->Branch("ClosestTrackDistance", &fClosestTrackDistance);
 // fTree ->Branch("DistanceToRecoVertex", &fDistanceToRecoVertex);
  fTree ->Branch("run", &fRun, "run/I");
  fTree ->Branch("subrun", &fSubRun, "subrun/I");
 // fTree ->Branch("LambdaDecayVertexX", &fLambdaDecayVertexX);
 // fTree ->Branch("LambdaDecayVertexY", &fLambdaDecayVertexY);
 // fTree ->Branch("LambdaDecayVertexZ", &fLambdaDecayVertexZ);
 // fTree ->Branch("LambdaDecayDistance", &fLambdaDecayDistance);
 // fTree ->Branch("ProtonMomentum", &fProtonMomentum);
 // fTree ->Branch("PionMomentum", &fPionMomentum);
 // fTree ->Branch("LambdaInvariantMass", &fLambdaInvariantMass);
  //fTree ->Branch("TrueTrackMomentum", &fTrueTrackMomentum);
  //fTree ->Branch("TrueTrackEnergy", &fTrueTrackEnergy);
 // fTree ->Branch("TrueTrackStartX", &fTrueTrackStartX);
 // fTree ->Branch("TrueTrackStartY", &fTrueTrackStartY);
 // fTree ->Branch("TrueTrackStartZ", &fTrueTrackStartZ);
 // fTree ->Branch("PIDPFParticleIndex", &fPFParticleIndex);
 // fTree ->Branch("PIDTrackID", &fTrackID);
 // fTree ->Branch("AlgorithmName", &fAlgorithmName);
 // fTree ->Branch("PIDAssumedPdg", &fAssumedPdg);
 // fTree->Branch("MuonScore", &fMuonScore);
 // fTree->Branch("ProtonScore", &fProtonScore);
 // fTree->Branch("KaonScore", &fKaonScore);
 // fTree->Branch("PionScore", &fPionScore);
fTree->Branch("truePDGs",& ftruePDGs);
fTree->Branch("assumedPDGs", &fassumedPDGs);
fTree->Branch("chi2Scores", &fchi2Scores);
//NEW VARIABLES


fTree->Branch("TrackIDs", &fTrackIDs);
fTree->Branch("TrackLengths", &fTrackLengths);
fTree->Branch("TrackPDGs", &fTrackPDGs);
fTree->Branch("IsLongestTrack", &fIsLongestTrack);
fTree->Branch("IsClosestTrack", &fIsClosestTrack);
fTree->Branch("DistanceToLongestTrack", &fDistanceToLongestTrack);
fTree->Branch("DistanceToRecoVertex", &fDistanceToRecoVertex);


fTree->Branch("TrueTrackPDGs", &fTrueTrackPDGs);

fTree->Branch("MatchQuality", &fMatchQuality);


fTree->Branch("LongestTrackID", &fLongestTrackID);
fTree->Branch("LongestTrackLength", &fLongestTrackLength);

fTree->Branch("LongestTrackTruePDG", &fLongestTrackTruePDG);
fTree->Branch("LongestTrackDistanceToVertex", &fLongestTrackDistanceToVertex);
fTree->Branch("ClosestTrackID", &fClosestTrackID);
fTree->Branch("ClosestTrackDistance", &fClosestTrackDistance);

fTree->Branch("ClosestTrackTruePDG", &fClosestTrackTruePDG);
fTree->Branch("TrackScore", &fTrackScores);
fTree->Branch("MuonTrackScores", &fMuonTrackScores);
fTree->Branch("KaonTrackScores", &fKaonTrackScores);
fTree->Branch("ProtonTrackScores", &fProtonTrackScores);
fTree->Branch("PionTrackScores", &fPionTrackScores);

fTree->Branch("AllPFPMuonScores", &fAllPFPMuonScores);
fTree->Branch("AllPFPKaonScores", &fAllPFPKaonScores);
fTree->Branch("AllPFPProtonScores", &fAllPFPProtonScores);
fTree->Branch("AllPFPPionScores", &fAllPFPPionScores);
fTree->Branch("AllPFPElectronScores", &fAllPFPElectronScores);

fTree->Branch("SliceT0s", &fAllSliceT0s);
fTree->Branch("SliceTriggerTypes", &fAllSliceTriggerTypes);

    fTree->Branch("FlashMatchScores", &fFlashMatchScores);
    fTree->Branch("FlashPEs", &fFlashPEs);
 fTree->Branch("CosmicNuScores", &fCosmicNuScores);
    fTree->Branch("NeutrinoNuScores", &fNeutrinoNuScores);
}


void hyperon::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.

}



DEFINE_ART_MODULE(hyperon::AnalyzeEvents)
