#define tmvaPrep_cxx
#include "tmvaPrep.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void tmvaPrep::Loop()
{
//   In a ROOT session, you can do:
//      root> .L tmvaPrep.C
//      root> tmvaPrep t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   // Create new files and TTrees 
 
   TFile *sigFile = new TFile("/data/ooconnor/sbnd/hyperons/preselection_output/tmvaSample_sig.root", "RECREATE");
   TTree *sigTree = fChain->CloneTree(0);
   sigTree->SetName("sigTree");
   sigTree->SetDirectory(sigFile);

   // Initialize TMVA variables
   float track1Length;
   float track2Length;
   float track3Length;
   float shower1Length;

   float track1StartPosX;
   float track1StartPosY;
   float track1StartPosZ;
   float track2StartPosX;
   float track2StartPosY;
   float track2StartPosZ;
   float track3StartPosX;
   float track3StartPosY;
   float track3StartPosZ;
   float shower1StartPosX;
   float shower1StartPosY;
   float shower1StartPosZ;

   float track1StartDirX;
   float track1StartDirY;
   float track1StartDirZ;
   float track2StartDirX;
   float track2StartDirY;
   float track2StartDirZ;
   float track3StartDirX;
   float track3StartDirY;
   float track3StartDirZ;
   float shower1DirX;
   float shower1DirY;
   float shower1DirZ;

   float track1DistRecoVtx;
   float track2DistRecoVtx;
   float track3DistRecoVtx;
   float shower1DistRecoVtx;

   float track1Track2Angle;
   float track1Track3Angle;
   float track2Track3Angle;
   float track1Shower1Angle;
   float track2Shower1Angle;
   float track3Shower1Angle;

   float track1Track2Dist;
   float track1Track3Dist;
   float track2Track3Dist;
   float track1Shower1Dist;
   float track2Shower1Dist;
   float track3Shower1Dist;

   // Tracks, shower lengths
   sigTree->Branch("track1Length", &track1Length);
   sigTree->Branch("track2Length", &track2Length);
   sigTree->Branch("track3Length", &track3Length);
   sigTree->Branch("shower1Length", &shower1Length);

   // Track, shower start positions
   sigTree->Branch("track1StartPosX", &track1StartPosX, "track1StartPosX/F");
   sigTree->Branch("track1StartPosY", &track1StartPosY, "track1StartPosY/F");
   sigTree->Branch("track1StartPosZ", &track1StartPosZ, "track1StartPosZ/F");

   sigTree->Branch("track2StartPosX", &track2StartPosX, "track2StartPosX/F");
   sigTree->Branch("track2StartPosY", &track2StartPosY, "track2StartPosY/F");
   sigTree->Branch("track2StartPosZ", &track2StartPosZ, "track2StartPosZ/F");

   sigTree->Branch("track3StartPosX", &track3StartPosX, "track3StartPosX/F");
   sigTree->Branch("track3StartPosY", &track3StartPosY, "track3StartPosY/F");
   sigTree->Branch("track3StartPosZ", &track3StartPosZ, "track3StartPosZ/F");

   sigTree->Branch("shower1StartPosX", &shower1StartPosX, "shower1StartPosX/F");
   sigTree->Branch("shower1StartPosY", &shower1StartPosY, "shower1StartPosY/F");
   sigTree->Branch("shower1StartPosZ", &shower1StartPosZ, "shower1StartPosZ/F");


// Track, shower direction vectors
   sigTree->Branch("track1StartDirX", &track1StartDirX, "track1StartDirX/F");
   sigTree->Branch("track1StartDirY", &track1StartDirY, "track1StartDirY/F");
   sigTree->Branch("track1StartDirZ", &track1StartDirZ, "track1StartDirZ/F");

   sigTree->Branch("track2StartDirX", &track2StartDirX, "track2StartDirX/F");
   sigTree->Branch("track2StartDirY", &track2StartDirY, "track2StartDirY/F");
   sigTree->Branch("track2StartDirZ", &track2StartDirZ, "track2StartDirZ/F");

   sigTree->Branch("track3StartDirX", &track3StartDirX, "track3StartDirX/F");
   sigTree->Branch("track3StartDirY", &track3StartDirY, "track3StartDirY/F");
   sigTree->Branch("track3StartDirZ", &track3StartDirZ, "track3StartDirZ/F");

   sigTree->Branch("shower1DirX", &shower1DirX);
   sigTree->Branch("shower1DirY", &shower1DirY);
   sigTree->Branch("shower1DirZ", &shower1DirZ);

// Distances to reconstructed vertex
   sigTree->Branch("track1DistRecoVtx", &track1DistRecoVtx, "track1DistRecoVtx/F");
   sigTree->Branch("track2DistRecoVtx", &track2DistRecoVtx, "track2DistRecoVtx/F");
   sigTree->Branch("track3DistRecoVtx", &track3DistRecoVtx, "track3DistRecoVtx/F");

   sigTree->Branch("shower1DistRecoVtx", &shower1DistRecoVtx, "shower1DistRecoVtx/F");

// Track-track, track-shower angles
   sigTree->Branch("track1Track2Angle", &track1Track2Angle);
   sigTree->Branch("track1Track3Angle", &track1Track3Angle);
   sigTree->Branch("track2Track3Angle", &track2Track3Angle);
   sigTree->Branch("track1Shower1Angle", &track1Shower1Angle);
   sigTree->Branch("track2Shower1Angle", &track2Shower1Angle);
   sigTree->Branch("track3Shower1Angle", &track3Shower1Angle);

// Track-track, track-shower distances
   sigTree->Branch("track1Track2Dist", &track1Track2Dist);
   sigTree->Branch("track1Track3Dist", &track1Track3Dist);
   sigTree->Branch("track2Track3Dist", &track2Track3Dist);
   sigTree->Branch("track1Shower1Dist", &track1Shower1Dist);
   sigTree->Branch("track2Shower1Dist", &track2Shower1Dist);
   sigTree->Branch("track3Shower1Dist", &track3Shower1Dist);

   // Clone tree structure for background TTree
   TFile *bkgFile = new TFile("/data/ooconnor/sbnd/hyperons/preselection_output/tmvaSample_bkg.root", "RECREATE");
   TTree *bkgTree = sigTree->CloneTree(0);
   bkgTree->SetName("bkgTree");
   bkgTree->SetTitle("Background Tree");
   bkgTree->SetDirectory(bkgFile);

   int nSig = 0;
   int nBkg = 0;
   int nBadSig = 0;
   int nBadBkg = 0;

   // define distance cut thresholds based on visual inspection of distributions

   int track1DistRecoVtxCut = 30;
   int track2DistRecoVtxCut = 30;
   int track3DistRecoVtxCut = 30;
   int shower1DistRecoVtxCut = 150;

   int track1Track2DistCut = 25;
   int track1Track3DistCut = 25;
   int track2Track3DistCut = 25;
   int track1Shower1DistCut = 150;
   int track2Shower1DistCut = 150;
   int track3Shower1DistCut = 150;

   // define length cut thresholds

   int track1LengthCut = 500;
   int track2LengthCut = 125;
   int track3LengthCut = 60;
   int shower1LengthCut = 80;

   // define counters for cuts

   int nTrk1DistRecoVtxCut = 0;
   int nTrk2DistRecoVtxCut = 0;
   int nTrk3DistRecoVtxCut = 0;
   int nShwr1DistRecoVtxCut = 0;
   int nTrk1Trk2DistCut = 0;
   int nTrk1Trk3DistCut = 0;
   int nTrk2Trk3DistCut = 0;
   int nTrk1Shwr1DistCut = 0;
   int nTrk2Shwr1DistCut = 0;
   int nTrk3Shwr1DistCut = 0;
   int nTrk1LengthCut = 0;
   int nTrk2LengthCut = 0;
   int nTrk3LengthCut = 0;
   int nShwr1LengthCut = 0;

   int nBadTrack = 0;
   int nBadShower = 0;

   int nBadTrkStartSig = 0;
   int nBadTrkStartBkg = 0;
   int nBadTrkDirSig = 0;
   int nBadTrkDirBkg = 0;
   int nBadShwrDirSig = 0;
   int nBadShwrDirBkg = 0;
   int nBadShwrStartSig = 0;
   int nBadShwrStartBkg = 0;
   int nBadTrkShwrSig = 0;
   int nBadTrkShwrBkg = 0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      std::cout<<"********* Event #"<<jentry<<" ***************"<<std::endl;
      std::cout<<"sampleType = "<<sampleType<<std::endl;
      std::cout<<"trackCount = "<<trackCount<<std::endl;
      std::cout<<"trackLengths size = "<<trackLengths->size()<<std::endl;
      std::cout<<"showerCount = "<<showerCount<<std::endl;
      std::cout<<"showerLengths size = "<<showerLengths->size()<<std::endl;

      // apply stricter cuts

      bool IsBadTrack = false;
      bool IsBadShower = false;

      // flag if any tracks or shower start pos is outside reco FV
      for (int i = 0; i < trackLengths->size(); i++){
	//      std::cout<<"track "<<i<<" start position = "<<trackStartPositionX->at(i)<<", "<<trackStartPositionY->at(i)<<", "<<trackStartPositionZ->at(i)<<std::endl;
	      if(trackStartPositionX->at(i) > 180 || std::abs(trackStartPositionY->at(i)) > 180
			      || trackStartPositionZ->at(i) < 0 || trackStartPositionZ->at(i) > 450){
		      if(sampleType == 0){nBadTrkStartSig++;}
		      if(sampleType == 2){nBadTrkStartBkg++;} 
		      IsBadTrack = true;
	      }
      }

      //std::cout<<"shower start position is "<<showerStartPositionX->at(0)<<", "<<showerStartPositionY->at(0)<<", "<<showerStartPositionZ<<std::endl; 
      if(std::abs(showerStartPositionX->at(0)) > 180 || std::abs(showerStartPositionY->at(0)) > 180 
			      || showerStartPositionZ->at(0) < 10 || showerStartPositionZ->at(0) > 450){
	      IsBadShower = true;
	      if(sampleType == 0){nBadShwrStartSig++;}
	      if(sampleType == 2){nBadShwrStartBkg++;}
      }

      // flag if any tracks or shower directions are outside cosine range
      for (int i = 0; i < trackLengths->size(); i++){
	  std::cout<<"track "<<i<<" direction is "<<trackStartDirX->at(i)<<", "<<trackStartDirY->at(i)<<", "<<trackStartDirZ->at(i)<<std::endl;
	  if(std::abs(trackStartDirX->at(i)) > 1 || std::abs(trackStartDirY->at(i)) > 1 || std::abs(trackStartDirZ->at(i)) > 1){
		  IsBadTrack = true;
		  if(sampleType == 0){nBadTrkDirSig++;}
		  if(sampleType == 2){nBadTrkDirBkg++;}
		  //std::cout<<"bad track due to direction out of cosine range"<<std::endl;
	  }
      }

      if(std::abs(showerDirX->at(0)) > 1 || std::abs(showerDirY->at(0)) > 1 || std::abs(showerDirZ->at(0)) > 1){
	      if(sampleType == 0){nBadShwrDirSig++;}
	      if(sampleType == 2){nBadShwrDirBkg++;} 
	      IsBadShower = true;
      }

      if(IsBadTrack || IsBadShower){
	      if(sampleType==0){nBadTrkShwrSig++; nBadSig++;}
	      if(sampleType==2){nBadTrkShwrBkg++; nBadBkg++;} 
	      continue;
      } // exclude events that have any tracks or showers outside reco fv

      // sort tracks in descending order
      std::vector<int> idx = {0, 1, 2};
      std::sort(idx.begin(), idx.end(), [&](int a, int b){return trackLengths->at(a) > trackLengths->at(b);});
      
      // assign index to tracks in descending order: LONGESTEST = track1, SHORTEST = track3
      int trk1Idx = idx[0];
      int trk2Idx = idx[1];
      int trk3Idx = idx[2];
      int shwr1Idx = 0; // only 1 shower

      std::cout<<"track1 index = "<<trk1Idx<<std::endl;
      std::cout<<"track2 index = "<<trk2Idx<<std::endl;
      std::cout<<"track3 index = "<<trk3Idx<<std::endl;

      // define recoVtx
      TVector3 recoVtx(RecoVertexX, RecoVertexY, RecoVertexZ);

      // define ordered lengths
      track1Length = trackLengths->at(trk1Idx);
      track2Length = trackLengths->at(trk2Idx);
      track3Length = trackLengths->at(trk3Idx);
      shower1Length = showerLengths->at(shwr1Idx);

      std::cout<<"track1 Length = "<<track1Length<<std::endl;
      std::cout<<"track2 Length = "<<track2Length<<std::endl;
      std::cout<<"track3 Length = "<<track3Length<<std::endl;

      TVector3 track1StartPos(trackStartPositionX->at(trk1Idx), trackStartPositionY->at(trk1Idx), trackStartPositionZ->at(trk1Idx));
      TVector3 track2StartPos(trackStartPositionX->at(trk2Idx), trackStartPositionY->at(trk2Idx), trackStartPositionZ->at(trk2Idx));
      TVector3 track3StartPos(trackStartPositionX->at(trk3Idx), trackStartPositionY->at(trk3Idx), trackStartPositionZ->at(trk3Idx));
      TVector3 shower1StartPos(showerStartPositionX->at(shwr1Idx), showerStartPositionY->at(shwr1Idx), showerStartPositionZ->at(shwr1Idx));

      track1StartPosX = track1StartPos.X();
      track1StartPosY = track1StartPos.Y();
      track1StartPosZ = track1StartPos.Z();

      track2StartPosX = track2StartPos.X();
      track2StartPosY = track2StartPos.Y();
      track2StartPosZ =  track2StartPos.Z();

      track3StartPosX = track3StartPos.X();
      track3StartPosY = track3StartPos.Y();
      track3StartPosZ = track3StartPos.Z();

      shower1StartPosX = shower1StartPos.X();
      shower1StartPosY = shower1StartPos.Y();
      shower1StartPosZ = shower1StartPos.Z();
      
      TVector3 track1StartDir(trackStartDirX->at(trk1Idx), trackStartDirY->at(trk1Idx), trackStartDirZ->at(trk1Idx));
      TVector3 track2StartDir(trackStartDirX->at(trk2Idx), trackStartDirY->at(trk2Idx), trackStartDirZ->at(trk2Idx));
      TVector3 track3StartDir(trackStartDirX->at(trk3Idx), trackStartDirY->at(trk3Idx), trackStartDirZ->at(trk3Idx));
      TVector3 shower1Dir(showerDirX->at(shwr1Idx), showerDirY->at(shwr1Idx), showerDirZ->at(shwr1Idx));

      track1StartDirX = track1StartDir.X();
      track1StartDirY = track1StartDir.Y();
      track1StartDirZ = track1StartDir.Z();
      track2StartDirX = track2StartDir.X();
      track2StartDirY = track2StartDir.Y();
      track2StartDirZ = track2StartDir.Z();
      track3StartDirX = track3StartDir.X();
      track3StartDirY = track3StartDir.Y();
      track3StartDirZ = track3StartDir.Z();
      shower1DirX = shower1Dir.X();
      shower1DirY = shower1Dir.Y();
      shower1DirZ = shower1Dir.Z();

      track1DistRecoVtx = (track1StartPos - recoVtx).Mag();
      track2DistRecoVtx = (track2StartPos - recoVtx).Mag();
      track3DistRecoVtx = (track3StartPos - recoVtx).Mag();
      shower1DistRecoVtx = (shower1StartPos - recoVtx).Mag();

      track1Track2Angle = track1StartDir.Angle(track2StartDir);
      track1Track3Angle = track1StartDir.Angle(track3StartDir);
      track2Track3Angle = track2StartDir.Angle(track3StartDir);
      track1Shower1Angle = track1StartDir.Angle(shower1Dir);
      track2Shower1Angle = track2StartDir.Angle(shower1Dir);
      track3Shower1Angle = track3StartDir.Angle(shower1Dir);

      track1Track2Dist = (track1StartPos - track2StartPos).Mag();
      track1Track3Dist = (track1StartPos - track3StartPos).Mag();
      track2Track3Dist = (track2StartPos - track3StartPos).Mag();
      track1Shower1Dist = (track1StartPos - shower1StartPos).Mag();
      track2Shower1Dist = (track2StartPos - shower1StartPos).Mag();
      track3Shower1Dist = (track3StartPos - shower1StartPos).Mag();

      // apply visual inspec distance cuts
      
      if (track1DistRecoVtx > track1DistRecoVtxCut){
	     if(sampleType==0){nBadSig++; nTrk1DistRecoVtxCut++;}
	     if(sampleType==2){nBadBkg++;}
	      continue;
      }
      if(track2DistRecoVtx > track2DistRecoVtxCut){
	      if(sampleType==0){nBadSig++; nTrk2DistRecoVtxCut++;}
	      if(sampleType==2){nBadBkg++;}
	      continue;}
      if(track3DistRecoVtx > track3DistRecoVtxCut){
	      if(sampleType==0){nBadSig++; nTrk3DistRecoVtxCut++;}
	      if(sampleType==2){nBadBkg++;}
	      continue;}
      if(shower1DistRecoVtx > shower1DistRecoVtxCut){
	      if(sampleType==0){nBadSig++; nShwr1DistRecoVtxCut++;}
	      if(sampleType==2){nBadBkg++;}
	      continue;}
      
      if(track1Track2Dist > track1Track2DistCut){
	      if(sampleType==0){nBadSig++; nTrk1Trk2DistCut++;}
	      if(sampleType==2){nBadBkg++;}
	      continue;}
      if(track1Track3Dist > track1Track3DistCut){
	      if(sampleType==0){nBadSig++; nTrk1Trk3DistCut++;}
	      if(sampleType==2){nBadBkg++;}
	      continue;}
      if(track2Track3Dist > track2Track3DistCut){
	      if(sampleType==0){nBadSig++; nTrk2Trk3DistCut++;}
	      if(sampleType==2){nBadBkg++;}
	      continue;}
      if(track1Shower1Dist > track1Shower1DistCut){
	      if(sampleType==0){nBadSig++; nTrk1Shwr1DistCut++;}
	      if(sampleType==2){nBadBkg++;}
	      continue;}
      if(track2Shower1Dist > track2Shower1DistCut){
	      if(sampleType==0){nBadSig++; nTrk2Shwr1DistCut++;}
	      if(sampleType==2){nBadBkg++;}
	      continue;}
      if(track3Shower1Dist > track3Shower1DistCut){
	     if(sampleType==0){nBadSig++; nTrk3Shwr1DistCut++;}
	     if(sampleType==2){nBadBkg++;}
	      continue;}

      // apply inspec length cuts

      if(track1Length > track1LengthCut){
	     if(sampleType==0){nBadSig++; nTrk1LengthCut++;}
	     if(sampleType==2){nBadBkg++;}
	      continue;}
      if(track2Length > track2LengthCut){
	     if(sampleType==0){nBadSig++; nTrk2LengthCut++;}
	     if(sampleType==2){nBadBkg++;}
	      continue;}
      if(track3Length > track3LengthCut){
	      if(sampleType==0){nBadSig++; nTrk3LengthCut++;}
	      if(sampleType==2){nBadBkg++;}
	      continue;}
      if(shower1Length > shower1LengthCut){
	      if(sampleType==0){nBadSig++; nShwr1LengthCut++;}
	      if(sampleType==2){nBadBkg++;}
	      continue;}

      // Fill signal and background trees:

      std::cout<<"Filling trees for event ID "<<eventID<<" with sampleType "<<sampleType<<std::endl;
      if(sampleType == 0){nSig++; sigTree->Fill();}
      if(sampleType == 2){nBkg++; bkgTree->Fill();}
      
   }

   int nBad = nBadSig + nBadBkg;
   int nBadTrkStart = nBadTrkStartSig + nBadTrkStartBkg;
   int nBadTrkDir = nBadTrkDirSig + nBadTrkDirBkg;
   int nBadShwrStart = nBadShwrStartSig + nBadShwrStartBkg;
   int nBadShwrDir = nBadShwrDirSig + nBadShwrDirBkg;
   int nBadLengths = nTrk1LengthCut + nTrk2LengthCut + nTrk3LengthCut + nShwr1LengthCut;
   int nBadDistRecoVtx = nTrk1DistRecoVtxCut + nTrk2DistRecoVtxCut + nTrk3DistRecoVtxCut + nShwr1DistRecoVtxCut;

   std::cout<<"Dropped "<<nBad<<" events ("<<nBadSig<<" signal, "<<nBadBkg<<" background)"<<std::endl;
   std::cout<<"Dropped "<<nBadTrkStart<<" events due to track start outside FV("<<nBadTrkStartSig<< "signal, "<<nBadTrkStartBkg<<" background)"<<std::endl;
   std::cout<<"Dropped "<<nBadShwrStart<<"events due to shower start outside FV ("<<nBadShwrStartSig<<" signal, "<<nBadShwrStartBkg<<" background)"<<std::endl;
   std::cout<<"Dropped "<<nBadTrkDir<<" events due to track dir outside cosine range("<<nBadTrkDirSig<<"signal, "<<nBadTrkDirBkg<<" background)"<<std::endl;
   std::cout<<"Dropped "<<nBadShwrDir<<" events due to shower dir outside cosine range("<<nBadShwrDirSig<<"signal, "<<nBadShwrDirBkg<<" background)"<<std::endl;
   std::cout<<"Dropped "<<nBadLengths<<" events because of track/shower lengths cut"<<std::endl;
   std::cout<<"Dropped "<<nBadDistRecoVtx<<" events because of distance to reco vertex cut"<<std::endl;
   std::cout<<"nTrk1DistRecoVtxCut =  "<<nTrk1DistRecoVtxCut<<std::endl;
   std::cout<<"nTrk2DistRecoVtxCut = "<<nTrk2DistRecoVtxCut<<std::endl;
   std::cout<<"nTrk3DistRecoVtxCut = "<<nTrk3DistRecoVtxCut<<std::endl;
   std::cout<<"nShwr1DistRecoVtxCut = "<<nShwr1DistRecoVtxCut<<std::endl;
   std::cout<<"nTrk1Trk2DistCut = "<<nTrk1Trk2DistCut<<std::endl;
   std::cout<<"nTrk1Trk3DistCut = "<<nTrk1Trk3DistCut<<std::endl;
   std::cout<<"nTrk2Trk3DistCut = "<<nTrk2Trk3DistCut<<std::endl;
   std::cout<<"nTrk1Shwr1DistCut = "<<nTrk1Shwr1DistCut<<std::endl;
   std::cout<<"nTrk2Shwr1DistCut = "<<nTrk2Shwr1DistCut<<std::endl;
   std::cout<<"nTrk3Shwr1DistCut = "<<nTrk3Shwr1DistCut<<std::endl;
   std::cout<<"nTrk1LengthCut = "<<nTrk1LengthCut<<std::endl;
   std::cout<<"nTrk2LengthCut = "<<nTrk2LengthCut<<std::endl;
   std::cout<<"nTrk3LengthCut = "<<nTrk3LengthCut<<std::endl;
   std::cout<<"nShwr1LengthCut = "<<nShwr1LengthCut<<std::endl;
   std::cout<<"Signal left after track/shower start outside FV = "<<std::endl;
   std::cout<<"Total signal left = "<<nSig<<" events"<<std::endl;
   std::cout<<"Total background left = "<<nBkg<<" events"<<std::endl; 

   sigFile->cd();
   sigTree->Write("sigTree");
   sigFile->Close();
   delete sigFile;

   bkgFile->cd();
   bkgTree->Write("bkgTree");
   bkgFile->Close();
   delete bkgFile;

}
