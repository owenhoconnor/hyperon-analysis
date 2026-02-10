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
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"

// // Root Includes
#include <iostream>
#include <vector>
#include <TTree.h>
#include <TH1.h>
#include <string>
#include <exception>
#include <fstream>

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

  // Declare member data here.

  int fVerbose;

  //std::string fHitLabel;
  std::string fGenieGenModuleLabel;
  std::string fPFParticleLabel;
  std::string fSliceLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  //std::string fClusterLabel;

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
   float isoVertexX, isoVertexY, isoVertexZ;
   bool fFoundRecoVertex;

  TTree *fTree;
  unsigned int fEventID;
  std::vector<int> fNPfpSlices;
  std::vector<int> fNSlices;
  size_t totalSlices = 0;
  size_t totalNeutrinos = 0;
  int nuID = -1;
  int nuSliceKey = -1;
  int fNPrimaryParticles;
  int fNPrimaryChildren;
  // MC Truth parameters
  std::vector<float> trueP;
  std::vector<int> truePDG; 
  std::vector<int> daughterPDG;
  std::vector<int> motherPDG;
  std::vector<float> vertexX;
  std::vector<float> vertexY;
  std::vector<float> vertexZ;
  std::vector<int> vertexSize;
  std::vector<int> daughterSize;

  // reco parameters
  std::vector<int> fTrackIDs;
  std::vector<float> fTrackLengths;
  std::vector<int> fTrackPDGs;
  float fRecoVertexX;
  float fRecoVertexY;
  float fRecoVertexZ;
  std::vector<float> fDistanceToRecoVertex;
  std::vector<float> fTrackStartPositionX;
  std::vector<float> fTrackStartPositionY;
  std::vector<float> fTrackStartPositionZ;
  std::vector<float> fTrackEndPositionX;
  std::vector<float> fTrackEndPositionY;
  std::vector<float> fTrackEndPositionZ;
  std::vector<int> fShowerPDG;
  std::vector<float> fnuScore;
  std::vector<float> fNeutrinoNuScores;
  std::vector<float> fCosmicNuScores;
  std::vector<float> fTrackScores;
  int fTrackCount = 0;
  int fShowerCount = 0;
  
  float highestNuScore = -1;
  TVector3 fRecoVertex;

  bool isSignal = false; // define bools for signal and background
  bool isBkg = false;
};


hyperon::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset} 
  //, fHitLabel(pset.get<std::string>("HitLabel"))
  , fGenieGenModuleLabel(pset.get<std::string>("GenieGenModuleLabel"))
  , fPFParticleLabel(pset.get<std::string>("PFParticleLabel"))
  , fSliceLabel(pset.get<std::string>("SliceLabel"))
  , fTrackLabel(pset.get<std::string>("TrackLabel"))
  //, fCalorimetryLabel(pset.get<std::string>("CalorimetryLabel"))
  , fShowerLabel(pset.get<std::string>("ShowerLabel"))
  //, fClusterLabel(pset.get<std::string>("ClusterLabel"))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void hyperon::AnalyzeEvents::analyze(art::Event const& evt)
{
  // Implementation of required member function here.
 fEventID = evt.id().event(); 
 std::cout<<"Event# "<<evt.id().event()<<std::endl;

 // Clear reco parameters
 fTrackIDs.clear();
 fTrackLengths.clear();
 fTrackPDGs.clear();
 fDistanceToRecoVertex.clear();
 fTrackStartPositionX.clear();
 fTrackStartPositionY.clear();
 fTrackStartPositionZ.clear();
 fTrackEndPositionX.clear();
 fTrackEndPositionY.clear();
 fTrackEndPositionZ.clear();
 fShowerPDG.clear();
 fnuScore.clear();
 fNeutrinoNuScores.clear();
 fCosmicNuScores.clear();
 fTrackScores.clear();
 highestNuScore = 0.;
 fFoundRecoVertex = false;
 fRecoVertex.SetXYZ(-9999., -9999., -9999.);
 fRecoVertexX = -9999.;
 fRecoVertexY = -9999.;
 fRecoVertexZ = -9999.;

 fNPrimaryParticles = 0;
 fNPrimaryChildren = 0;
 fTrackCount = 0;
 fShowerCount = 0;


 // Get event slices
 
   art::ValidHandle<std::vector<recob::Slice>> sliceHandle = evt.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
   std::vector<art::Ptr<recob::Slice>> sliceVector;
   if (sliceHandle.isValid())
	   art::fill_ptr_vector(sliceVector, sliceHandle);

// Get associations between slices and PFParticles

   art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, evt, fSliceLabel);
   art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

// Get associations between PFPParticles and PFPMetaData

   art::FindManyP<larpandoraobj::PFParticleMetadata> pfpMetadataAssoc(pfpHandle, evt, fPFParticleLabel);

// * MC truth information
   art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
   std::vector<art::Ptr<simb::MCTruth> > mclist;
   if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);

   art::FindManyP< simb::MCParticle > fmpart( mctruthListHandle, evt, "largeant" );

// Get associations between PFParticles and Clusters, and between Hits and MCParticles

   art::FindManyP<recob::Cluster> pfpClusterAssoc(pfpHandle, evt, fPFParticleLabel);
  // art::FindManyP<recob::Hit> clusterHitAssoc(evt.getValidHandle<std::vector<recob::Cluster>>(fClusterLabel), evt, fClusterLabel);
  // art::FindManyP<simb::MCParticle> hitMCParticleAssoc(evt.getValidHandle<std::vector<recob::Hit>>(fHitLabel), evt, fHitLabel);

// Get associations between PFParticle and Vertex

   art::FindManyP<recob::Vertex> pfpVertexAssoc(pfpHandle, evt, fPFParticleLabel);
   std::cout<<"Event "<<fEventID<<" has "<<sliceVector.size()<<" slices."<<std::endl;

// Filling our neutrino hierarchy variables

   for (const art::Ptr<recob::Slice> &slice : sliceVector){

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


	fNPfpSlices.push_back(slicePFPs.size());
	std::cout<<"Slice key: "<< slice.key()<<", Number of PFPs: "<< slicePFPs.size() << std::endl;


	for (const art::Ptr<recob::PFParticle> &slicePFP : slicePFPs){

		const bool isPrimary = (slicePFP->IsPrimary());
		const bool isNeutrino = (std::abs(slicePFP->PdgCode()) == 14);

		if (!(isPrimary && isNeutrino)){
			//std::cout<<"Not a primary neutrino, skipping PFP!"<<std::endl;
			continue;
		}

		std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metadataVec = pfpMetadataAssoc.at(slicePFP.key());

		if (metadataVec.empty()){
			std::cerr<<"No meta data found for PFParticle with key :"<<slice.key()<<std::endl;
		}

		float nuScore = -1;
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

		for (const auto &particle : slicePFPs){
		std::cout << particle->PdgCode()<<", ";
		}
		std::cout<<std::endl;
	
		// We have found our neutrino!!
		totalNeutrinos++;

		if (nuScore > highestNuScore) {
			highestNuScore = nuScore;
			nuSliceKey = slice.key();
			nuID = slicePFP->Self();
			fNPrimaryParticles = slicePFPs.size();
			fNPrimaryChildren = slicePFP->NumDaughters();
			std::cout<<"Highest nuScore overwritten! New highestNuScore: "<<highestNuScore<<std::endl;
			std::cout<<"new nuSliceKey = "<<nuSliceKey<<std::endl;
			std::cout<<"new nuID = "<<nuID<<std::endl;
		}

		else {
			std::cout<<"This PFP nuScore is not the highest. Highest nuScore remains: "<<highestNuScore<<std::endl;
		}
	}	

   } 


// Define vector of PFPs in nuSlice
   std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs(slicePFPAssoc.at(nuSliceKey));
   std::cout<<"!!! Now looping through nuSlicePFPs of size = "<<nuSlicePFPs.size()<<" !!!"<<std::endl;

// try to find reco vertex by looping over nuSlicePFPs

   for (const art::Ptr<recob::PFParticle>& nuSlicePFP : nuSlicePFPs) {
	if (nuSlicePFP->IsPrimary() && std::abs(nuSlicePFP->PdgCode()) == 14) {
		auto vertices = pfpVertexAssoc.at(nuSlicePFP.key());

		if (!vertices.empty()) {
			const recob::Vertex& vertex = *vertices.at(0);
			auto const& vertexPos = vertex.position();
			fRecoVertex.SetXYZ(vertexPos.X(), vertexPos.Y(), vertexPos.Z());
			fRecoVertexX = fRecoVertex.X();
			fRecoVertexY = fRecoVertex.Y();
			fRecoVertexZ = fRecoVertex.Z();
			fFoundRecoVertex = true;
			break;
		}	
	}
   }


   if (!fFoundRecoVertex) {
	std::cerr<<"Error: Reco Vertex not found for the neutrino PFPParticle!"<<std::endl;
	fRecoVertex.SetXYZ(-9999., -9999., -9999.); // store placeholder values
   }

   std::cout<<"RecoVertex: (" << fRecoVertex.X() << ", "<< fRecoVertex.Y() << ", "<< fRecoVertex.Z() << ")"<<std::endl;

// Save Reco Parameters

   art::ValidHandle<std::vector<recob::Track>> trackHandle = evt.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
   if (!trackHandle.isValid()){
	   std::cerr<<"Warning: recob::Track collection not found in event. "<<std::endl;
   }
   art::ValidHandle<std::vector<recob::Shower>> showerHandle = evt.getValidHandle<std::vector<recob::Shower>>(fShowerLabel);
   art::FindManyP<recob::Track> pfpTrackAssoc(pfpHandle, evt, fTrackLabel);
   art::FindManyP<recob::Shower> pfpShowerAssoc(pfpHandle, evt, fShowerLabel);
   art::FindManyP<recob::PFParticle> trackToPFPAssoc(trackHandle, evt, fTrackLabel);
   art::FindManyP<recob::PFParticle> showerToPFPAssoc(showerHandle, evt, fShowerLabel);

   // Loop over nuSlice PFPs to count tracks and showers
   /*for (const art::Ptr<recob::PFParticle>& nuSlicePFP : nuSlicePFPs) {
	// Only interested in neutrino children	
	std::cout<<"nuSlicePFP Parent = "<<nuSlicePFP->Parent()<<" and nuID = "<<nuID<<std::endl;
	if (nuSlicePFP->Parent() != static_cast<long unsigned int>(nuID)){
		std::cout<<"skipping this slicePFP (reason: not a neutrino child)"<<std::endl;
		continue;
	}

	// Get pfp meta data
	std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadataVec = pfpMetadataAssoc.at(nuSlicePFP.key());

	float trackScore = -1.0;

	if (pfpMetadatVec.empty()){
		std::cerr<<"No metadata found for PFParticle with key "<<nuSlicePFP.key()<<std::endl;
	}

	for (const auto& metadata : pfpMetadataVec){
		const auto& propertiesMap = metadata->GetPropertiesMap();

		if (propertiesMap.find("TrackScore") != propertiesMap.end()){
			trackScore = propertiesMap.at("TrackScore");
			fTrackScores.push_back(trackScore);
			std::cout<<"Found track score for PFParticle: "<<trackScore<<std::endl;
		}

		else {
			std::cerr<<"No metadata found for PFParticle with key "<<nuSlicePFP.key()<<std::endl;
		}

	}

	if (trackScore > 0.5){
		fTrackCount++;
		std::cout<<"trackScore > 0.5, +1 to tracks, nTracks = "<<fTrackCount<<std::endl;
	}

	else{
		fShowerCount++;
		std::cout<<"trackScore < 0.5, +1 to showers, nShowers = "<<fShowerCount<<std::endl;

	}
   }*/

   for (const art::Ptr<recob::PFParticle>& nuSlicePFP : nuSlicePFPs) {
	// Only interested in neutrino children
	std::cout<<"nuSlicePFP Parent = "<<nuSlicePFP->Parent()<<" and nuID = "<<nuID<<std::endl;
	if (nuSlicePFP->Parent() != static_cast<long unsigned int>(nuID)){
		std::cout<<"skipping this slicePFP (reason: not a neutrino child)"<<std::endl;
		continue;
	}

	// Get tracks associated with this particle
	
	std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(nuSlicePFP.key());
	std::vector<art::Ptr<recob::Shower>> showers = pfpShowerAssoc.at(nuSlicePFP.key());

	// There should only be 0 or 1 tracks / showers associated with a PFP
	if (tracks.size() != 1 && showers.size() != 1){
		std::cout<<"Skipping this slicePFP (reason: tracks.size() != 1 and showers.size() != 1)"<<std::endl;
		continue;
	}

	// Get shower if there is one
	if (showers.size() == 1){
	    art::Ptr<recob::Shower> shower = showers.at(0); // call it track for easier code later
	    std::cout<<"Got shower: "<<shower<<std::endl;
	    float showerLength = shower->Length();
	    auto const& showerStart = shower->ShowerStart();
	    TVector3 showerStartPos(showerStart.X(), showerStart.Y(), showerStart.Z());

	    float distanceToVertex = -9999.0f;
	    if (fFoundRecoVertex) {
		distanceToVertex = (showerStartPos - fRecoVertex).Mag();
	    }

	    fTrackLengths.push_back(showerLength);
	    fTrackStartPositionX.push_back(showerStart.X());
	    fTrackStartPositionY.push_back(showerStart.Y());
	    fTrackStartPositionZ.push_back(showerStart.Z());
	    fDistanceToRecoVertex.push_back(distanceToVertex);
	}
	// Get PFPs associated with track
	//std::vector<art::Ptr<recob::PFParticle>> trackPFPs = trackToPFPAssoc.at(track.key());
	//art::Ptr<recob::PFParticle> trackPFP = trackPFPs.front();
	//std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> trackMetadataVec = pfpMetadataAssoc.at(trackPFP.key());
	
	// Get PFP Metadata
	std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadataVec = pfpMetadataAssoc.at(nuSlicePFP.key());
	float trackScore = -1.0;
	//if (trackMetadataVec.empty()){
	//	std::cerr<<"No meta found for track with ID = "<<track->ID()<<std::endl;
	//}
	
	if (pfpMetadataVec.empty()){
		std::cerr<<"No metadata found for PFParticle with key "<<nuSlicePFP.key()<<std::endl;
	}

	for (const auto& metadata : pfpMetadataVec) { // find trackScore in PFP metadata

		const auto& propertiesMap = metadata->GetPropertiesMap();
		//std::cout<<"PFParticle Metadata Properties (ID )"<< track->ID()<<"):\n"<<std::endl;
		//for (const auto& [propertyName, propertyValue] : propertiesMap){
		//	std::cout<<""<<propertyName<<" = "<<propertyValue<<std::endl;
		//}

		if (propertiesMap.find("TrackScore") != propertiesMap.end()){
			trackScore = propertiesMap.at("TrackScore");
			fTrackScores.push_back(trackScore);
			std::cout<<"Found track score for PFParticle: "<<trackScore<<std::endl;
		}

		else {
			std::cerr<<"No metadata found for PFParticle with key "<<nuSlicePFP.key()<<std::endl;
		}

	}

	// Classify PFP as track or shower with trackScore
	if (trackScore > 0.5){
		fTrackCount++;
		std::cout<<"trackScore > 0.5, +1 to tracks, nTracks = "<<fTrackCount<<std::endl;
	}

	else{
		fShowerCount++;
		std::cout<<"trackScore < 0.5, +1 to showers, nShowers = "<<fShowerCount<<std::endl;
	}

	// if there's a track (not a shower)
	if (tracks.size() == 1){
	    art::Ptr<recob::Track> track = tracks.at(0);
	    float trackLength = track->Length();
	    auto const& trackStart = track->Vertex();
	    auto const& trackEnd = track->End();
	    TVector3 trackStartPos(trackStart.X(), trackStart.Y(), trackStart.Z());
	    TVector3 trackEndPos(trackEnd.X(), trackEnd.Y(), trackEnd.Z());


	    float distanceToVertex = -9999.0f;
	    if (fFoundRecoVertex) {
		distanceToVertex = (trackStartPos - fRecoVertex).Mag();
	    }

	     std::cout<<"Pushing back reco parameters now"<<std::endl;
	    fTrackIDs.push_back(track->ID());
	    fTrackLengths.push_back(trackLength);
	    fTrackStartPositionX.push_back(trackStart.X());
	    fTrackStartPositionY.push_back(trackStart.Y());
	    fTrackStartPositionZ.push_back(trackStart.Z());
	    fTrackEndPositionX.push_back(trackEnd.X());
	    fTrackEndPositionY.push_back(trackEnd.Y());
	    fTrackEndPositionZ.push_back(trackEnd.Z());
	    fDistanceToRecoVertex.push_back(distanceToVertex);
	}

   }


// reco Track -> MC truth matching




//reco Shower -> MC truth matching

// MCTRUTH PARAMETERS

   std::cout<<"------------ MC TRUTH Parameters ----------------"<<std::endl;  
   nTotEvents++;
   vertexSize.push_back(0);

   //if (mclist.size() > 1){
//	nMultiEvents++;
  // }

  //if (mclist.size() == 1)
//	nSingleEvents++;

   for (size_t i_truth = 0; i_truth < mclist.size(); i_truth++)
   {
		
	   
       art::Ptr<simb::MCTruth> truth(mclist.at(i_truth));
       std::cout<<"Neutrino at index " << i_truth << " has pdg: " << truth->GetNeutrino().Nu().PdgCode() << std::endl;

 //      int lambda=0,mu=0, mubar=0,kaonp=0,kaonm=0,kaon0=0,proton=0,neutron=0,pip=0,pim=0,pi0=0;
 //      int sigmap=0, sigma0=0, sigmam=0;
       std::vector< art::Ptr<simb::MCParticle> > assocParticles(fmpart.at(i_truth));
       for (size_t i_mcpart = 0; i_mcpart < assocParticles.size(); i_mcpart++)
       {
           art::Ptr<simb::MCParticle> mcParticle(assocParticles.at(i_mcpart));
           if(mcParticle->Mother()!=10000000) continue;

           std::cout<<"--Particle at index " << i_mcpart <<"at mclist index "<<i_truth<< " has pdg: " << mcParticle->PdgCode() 
		    <<" has momentum: " << mcParticle->P()<<" Mother: "<<mcParticle->Mother()
		    <<" has status: "<< mcParticle->StatusCode()<<" number trajectory points: "<< mcParticle->NumberTrajectoryPoints()
		    <<" TrackId: " << mcParticle->TrackId()<< std::endl;
           trueP.push_back(mcParticle->P()); 
	   truePDG.push_back(mcParticle->PdgCode());
	   motherPDG.push_back(mcParticle->Mother());
	   nMCParticles++; // count num of saved mc particles in each interaction vertex
	  // Get interaction vertices
	  if (mcParticle->NumberTrajectoryPoints()>0) {
		 auto const& startPoint = mcParticle->Position(0);
		 vertexX.push_back(startPoint.X());
		 vertexY.push_back(startPoint.Y());
		 vertexZ.push_back(startPoint.Z());

		 isoVertexX = startPoint.X();
                 isoVertexY = startPoint.Y();
 		 isoVertexZ = startPoint.Z();
	  }
	  else{ // Store placeholders if no trajectory points for consistency
		  vertexX.push_back(-9999.);
		  vertexY.push_back(-9999.);
		  vertexZ.push_back(-9999.);

		  isoVertexX = -9999.;
		  isoVertexY = -9999.;
		  isoVertexZ = -9999.;
	  }

	  // Daughter information
	  daughterSize.push_back(mcParticle->NumberDaughters()); // save num of daughters for each particle for later indexing
	  std::cout<<"num daughters = "<<mcParticle->NumberDaughters()<<std::endl;
	  for (int i_daughter = 0; i_daughter < mcParticle->NumberDaughters(); ++i_daughter){
		  int daughterTrackID = mcParticle->Daughter(i_daughter);
	  for (size_t i = 0; i < assocParticles.size(); ++i){
		  if (assocParticles.at(i)->TrackId() == daughterTrackID){
			  int daughterPdgCode = assocParticles.at(i)->PdgCode();
			  daughterPDG.push_back(daughterPdgCode); // save all interaction particle daughter pdgs in one array
			  break;
		  }
	  }
	  }
           
	  if(mcParticle->PdgCode()==  13) mu++;   
	  if(mcParticle->PdgCode()== -13) mubar++;
	  if(mcParticle->PdgCode()== 211) pip++; 
	  if(mcParticle->PdgCode()==-211) pim++;
          if(mcParticle->PdgCode()== 111) pi0++;
	  if(mcParticle->PdgCode() == 22) gamma++;

	  
          if(mcParticle->PdgCode()== 321) kaonp++;
          if(mcParticle->PdgCode()==-321) kaonm++;
          if(mcParticle->PdgCode()== 311 || mcParticle->PdgCode() == 130 || mcParticle->PdgCode() == 310) kaon0++;
          if(mcParticle->PdgCode()== 2212){
		proton++;

	  for (int i_daughter = 0; i_daughter < mcParticle->NumberDaughters(); ++i_daughter){
		  int daughterTrackID = mcParticle->Daughter(i_daughter);
	  for (size_t i = 0; i < assocParticles.size(); ++i){
		  if (assocParticles.at(i)->TrackId() == daughterTrackID){
			  int daughterPdgCode = assocParticles.at(i)->PdgCode();
			  std::cout<<"PROTON DAUGHTER PDG = "<<daughterPdgCode<<std::endl;
			  break;
		 	 }
	  	    }
	  	}

	  }
          if(mcParticle->PdgCode()== 2112) neutron++;
          if(mcParticle->PdgCode()== 3122){
		  lambda++;
		  
	  for (int i_daughter = 0; i_daughter < mcParticle->NumberDaughters(); ++i_daughter){
		  int daughterTrackID = mcParticle->Daughter(i_daughter);
	  for (size_t i = 0; i < assocParticles.size(); ++i){
		  if (assocParticles.at(i)->TrackId() == daughterTrackID){
			  int daughterPdgCode = assocParticles.at(i)->PdgCode();
			  std::cout<<"LAMBDA DAUGHTER PDG = "<<daughterPdgCode<<std::endl;

			  if (daughterPdgCode == 2212){
				  std::cout<<"LAMBDA w/ PROTON DAUGHTER"<<std::endl;
				  goodLambda++;
			  }
			  break;
		 	 }
	  	    }
	  	}
	  }
          if(mcParticle->PdgCode()== 3222) sigmap++;
          if(mcParticle->PdgCode()== 3212){
		  sigma0++;

	  for (int i_daughter = 0; i_daughter < mcParticle->NumberDaughters(); ++i_daughter){
		  int daughterTrackID = mcParticle->Daughter(i_daughter);
	  for (size_t i = 0; i < assocParticles.size(); ++i){
		  if (assocParticles.at(i)->TrackId() == daughterTrackID){
			  int daughterPdgCode = assocParticles.at(i)->PdgCode();
			  std::cout<<"SIGMA0 DAUGHTER PDG = "<<daughterPdgCode<<std::endl;

			  if (daughterPdgCode == 3122){
				std::cout<<"SIGMA0 w/ LAMBDA DAUGHTER"<<std::endl;
				goodSigma++;
			  }
			  break;
		 	 }
	  	    }
	  	}
	  }
          if(mcParticle->PdgCode()== 3112) sigmam++;

	


      }
	vertexSize.push_back(nMCParticles); // save the number of particles at each interaction vertex
	nMCParticles = 0; // reset counter

   }
   
//std::cout<<"***************************"<<std::endl;
//std::cout<<vertexSize.size()<<std::endl;
//std::cout<<"isoVertex X is: "<<isoVertexX<<", isoVertexY is "<<isoVertexY<<", isoVertexZ is: "<<isoVertexZ<<std::endl;

//double multiRatio = double(nMultiEvents) / double(nTotEvents);
//double singleRatio = double(nSingleEvents) / double(nTotEvents);
//double sanityCheck = multiRatio + singleRatio;
//std::cout<<"num of multi int events = "<<nMultiEvents<<std::endl;
//std::cout<<"num of single int events = "<<nSingleEvents<<std::endl;
//std::cout<<"num of total events = "<<nTotEvents<<std::endl;
//std::cout<<"running ratio of multi int events = "<<multiRatio<<std::endl;
//std::cout<<"running ratio of single int events = "<<singleRatio<<std::endl;
//std::cout<<"sanity check 1 = "<<sanityCheck<<std::endl;

std::cout<<"num of good lambdas = "<<goodLambda<<std::endl;
std::cout<<"num of good sigmas = "<<goodSigma<<std::endl;
std::cout<<"****************************************************************"<<std::endl;


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
 //if isSignal{
   fTree->Fill();
 //}
 trueP.clear();
 truePDG.clear();
 daughterPDG.clear();
 motherPDG.clear();
 vertexX.clear();
 vertexY.clear();
 vertexZ.clear();
 vertexSize.clear();
 daughterSize.clear();
}


void hyperon::AnalyzeEvents::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree> ("tree", "Output TTree");

  //add branches here
  
  // MC truth parameters
  fTree->Branch("eventID", &fEventID);
  fTree->Branch("trueP", &trueP);
  fTree->Branch("truePDG", &truePDG);
  fTree->Branch("daughterPDG", &daughterPDG);
  fTree->Branch("motherPDG", &motherPDG);
  fTree->Branch("vertexX", &vertexX);
  fTree->Branch("vertexY", &vertexY);
  fTree->Branch("vertexZ", &vertexZ);
  fTree->Branch("vertexSize", &vertexSize);
  fTree->Branch("daughterSize", &daughterSize);
  fTree->Branch("isoVertexX", &isoVertexX);
  fTree->Branch("isoVertexY", &isoVertexY);
  fTree->Branch("isoVertexZ", &isoVertexZ);


  // reco parameters
  fTree->Branch("nPFParticles", &fNPrimaryParticles);
  fTree->Branch("nPrimaryChildren", &fNPrimaryChildren);
  fTree->Branch("trackCount", &fTrackCount);
  fTree->Branch("showerCount", &fShowerCount);
  fTree->Branch("TrackIDs", &fTrackIDs);
  fTree->Branch("TrackLengths", &fTrackLengths);
  fTree->Branch("RecoVertexX", &fRecoVertexX);
  fTree->Branch("RecoVertexY", &fRecoVertexY);
  fTree->Branch("RecoVertexZ", &fRecoVertexZ);
  fTree->Branch("DistanceToRecoVertex", &fDistanceToRecoVertex);
  fTree->Branch("nuScores", &fnuScore);
  fTree->Branch("trackScores", &fTrackScores);
  fTree->Branch("NeutrinoNuScores", &fNeutrinoNuScores);
  fTree->Branch("TrackStartPositionX", &fTrackStartPositionX);
  fTree->Branch("TrackStartPositionY", &fTrackStartPositionY);
  fTree->Branch("TrackStartPositionZ", &fTrackStartPositionZ);
  fTree->Branch("TrackEndPositionX", &fTrackEndPositionX);
  fTree->Branch("TrackEndPositionY", &fTrackEndPositionY);
  fTree->Branch("TrackEndPositionZ", &fTrackEndPositionZ);

  //Histograms
}

void hyperon::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.


}

DEFINE_ART_MODULE(hyperon::AnalyzeEvents)
