//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 12 15:43:48 2026 by ROOT version 6.40.02
// from TTree bkgTree/Background Tree
// found on file: /data/ooconnor/sbnd/hyperons/preselection_output/tmvaSample_bkg.root
//////////////////////////////////////////////////////////

#ifndef backgroundPlots_h
#define backgroundPlots_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class backgroundPlots {
public :
   TTree          *fChain;   ///<!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; ///<!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          eventID;
   Int_t           run;
   Int_t           subrun;
   vector<float>   *trueP;
   vector<int>     *truePDG;
   vector<int>     *daughterPDG;
   vector<int>     *motherPDG;
   vector<float>   *vertexX;
   vector<float>   *vertexY;
   vector<float>   *vertexZ;
   vector<int>     *vertexSize;
   vector<int>     *daughterSize;
   vector<float>    *trueNuEnergy;
   vector<int>     *trueCCNC;
   vector<int>     *trueIntMode;
   vector<int>     *trueIntType;
   Int_t           nPFParticles;
   Int_t           nPrimaryChildren;
   Int_t           trackCount;
   Int_t           showerCount;
   vector<int>     *TrackIDs;
   vector<float>   *trackLengths;
   Float_t         RecoVertexX;
   Float_t         RecoVertexY;
   Float_t         RecoVertexZ;
   vector<float>   *DistanceToRecoVertex;
   vector<float>   *nuScores;
   vector<float>   *trackScores;
   vector<float>   *muonTrackScores;
   vector<float>   *protonTrackScores;
   vector<float>   *pionTrackScores;
   vector<float>   *NeutrinoNuScores;
   vector<float>   *trackStartPositionX;
   vector<float>   *trackStartPositionY;
   vector<float>   *trackStartPositionZ;
   vector<float>   *trackEndPositionX;
   vector<float>   *trackEndPositionY;
   vector<float>   *trackEndPositionZ;
   vector<float>   *trackStartDirX;
   vector<float>   *trackStartDirY;
   vector<float>   *trackStartDirZ;
   vector<float>   *trackEndDirX;
   vector<float>   *trackEndDirY;
   vector<float>   *trackEndDirZ;
   vector<float>   *trackVertexDirX;
   vector<float>   *trackVertexDirY;
   vector<float>   *trackVertexDirZ;
   vector<float>   *trackTheta;
   vector<float>   *trackPhi;
   vector<float>   *showerLengths;
   vector<float>   *showerStartPositionX;
   vector<float>   *showerStartPositionY;
   vector<float>   *showerStartPositionZ;
   vector<float>   *showerDirX;
   vector<float>   *showerDirY;
   vector<float>   *showerDirZ;
   vector<int>     *pfpPDG;
   Int_t           sampleType;
   Int_t           chosenMCTruthIdx;
   Float_t         track1Length;
   Float_t         track2Length;
   Float_t         track3Length;
   Float_t         shower1Length;
   Float_t         track1StartPosX;
   Float_t         track1StartPosY;
   Float_t         track1StartPosZ;
   Float_t         track2StartPosX;
   Float_t         track2StartPosY;
   Float_t         track2StartPosZ;
   Float_t         track3StartPosX;
   Float_t         track3StartPosY;
   Float_t         track3StartPosZ;
   Float_t         shower1StartPosX;
   Float_t         shower1StartPosY;
   Float_t         shower1StartPosZ;
   Float_t         track1StartDirX;
   Float_t         track1StartDirY;
   Float_t         track1StartDirZ;
   Float_t         track2StartDirX;
   Float_t         track2StartDirY;
   Float_t         track2StartDirZ;
   Float_t         track3StartDirX;
   Float_t         track3StartDirY;
   Float_t         track3StartDirZ;
   Float_t         shower1DirX;
   Float_t         shower1DirY;
   Float_t         shower1DirZ;
   Float_t         track1DistRecoVtx;
   Float_t         track2DistRecoVtx;
   Float_t         track3DistRecoVtx;
   Float_t         shower1DistRecoVtx;
   Float_t         track1Track2Angle;
   Float_t         track1Track3Angle;
   Float_t         track2Track3Angle;
   Float_t         track1Shower1Angle;
   Float_t         track2Shower1Angle;
   Float_t         track3Shower1Angle;
   Float_t         track1Track2Dist;
   Float_t         track1Track3Dist;
   Float_t         track2Track3Dist;
   Float_t         track1Shower1Dist;
   Float_t         track2Shower1Dist;
   Float_t         track3Shower1Dist;

   // List of branches
   TBranch        *b_eventID;   ///<!
   TBranch        *b_run;   ///<!
   TBranch        *b_subrun;   ///<!
   TBranch        *b_trueP;   ///<!
   TBranch        *b_truePDG;   ///<!
   TBranch        *b_daughterPDG;   ///<!
   TBranch        *b_motherPDG;   ///<!
   TBranch        *b_vertexX;   ///<!
   TBranch        *b_vertexY;   ///<!
   TBranch        *b_vertexZ;   ///<!
   TBranch        *b_vertexSize;   ///<!
   TBranch        *b_daughterSize;   ///<!
   TBranch        *b_trueNuEnergy;   ///<!
   TBranch        *b_trueCCNC;   ///<!
   TBranch        *b_trueIntMode;   ///<!
   TBranch        *b_trueIntType;   ///<!
   TBranch        *b_nPFParticles;   ///<!
   TBranch        *b_nPrimaryChildren;   ///<!
   TBranch        *b_trackCount;   ///<!
   TBranch        *b_showerCount;   ///<!
   TBranch        *b_TrackIDs;   ///<!
   TBranch        *b_trackLengths;   ///<!
   TBranch        *b_RecoVertexX;   ///<!
   TBranch        *b_RecoVertexY;   ///<!
   TBranch        *b_RecoVertexZ;   ///<!
   TBranch        *b_DistanceToRecoVertex;   ///<!
   TBranch        *b_nuScores;   ///<!
   TBranch        *b_trackScores;   ///<!
   TBranch        *b_muonTrackScores;   ///<!
   TBranch        *b_protonTrackScores;   ///<!
   TBranch        *b_pionTrackScores;   ///<!
   TBranch        *b_NeutrinoNuScores;   ///<!
   TBranch        *b_trackStartPositionX;   ///<!
   TBranch        *b_trackStartPositionY;   ///<!
   TBranch        *b_trackStartPositionZ;   ///<!
   TBranch        *b_trackEndPositionX;   ///<!
   TBranch        *b_trackEndPositionY;   ///<!
   TBranch        *b_trackEndPositionZ;   ///<!
   TBranch        *b_trackStartDirX;   ///<!
   TBranch        *b_trackStartDirY;   ///<!
   TBranch        *b_trackStartDirZ;   ///<!
   TBranch        *b_trackEndDirX;   ///<!
   TBranch        *b_trackEndDirY;   ///<!
   TBranch        *b_trackEndDirZ;   ///<!
   TBranch        *b_trackVertexDirX;   ///<!
   TBranch        *b_trackVertexDirY;   ///<!
   TBranch        *b_trackVertexDirZ;   ///<!
   TBranch        *b_trackTheta;   ///<!
   TBranch        *b_trackPhi;   ///<!
   TBranch        *b_showerLengths;   ///<!
   TBranch        *b_showerStartPositionX;   ///<!
   TBranch        *b_showerStartPositionY;   ///<!
   TBranch        *b_showerStartPositionZ;   ///<!
   TBranch        *b_showerDirX;   ///<!
   TBranch        *b_showerDirY;   ///<!
   TBranch        *b_showerDirZ;   ///<!
   TBranch        *b_pfpPDG;   ///<!
   TBranch        *b_sampleType;   ///<!
   TBranch        *b_chosenMCTruthIdx;   ///<!
   TBranch        *b_track1Length;   ///<!
   TBranch        *b_track2Length;   ///<!
   TBranch        *b_track3Length;   ///<!
   TBranch        *b_shower1Length;   ///<!
   TBranch        *b_track1StartPosX;   ///<!
   TBranch        *b_track1StartPosY;   ///<!
   TBranch        *b_track1StartPosZ;   ///<!
   TBranch        *b_track2StartPosX;   ///<!
   TBranch        *b_track2StartPosY;   ///<!
   TBranch        *b_track2StartPosZ;   ///<!
   TBranch        *b_track3StartPosX;   ///<!
   TBranch        *b_track3StartPosY;   ///<!
   TBranch        *b_track3StartPosZ;   ///<!
   TBranch        *b_shower1StartPosX;   ///<!
   TBranch        *b_shower1StartPosY;   ///<!
   TBranch        *b_shower1StartPosZ;   ///<!
   TBranch        *b_track1StartDirX;   ///<!
   TBranch        *b_track1StartDirY;   ///<!
   TBranch        *b_track1StartDirZ;   ///<!
   TBranch        *b_track2StartDirX;   ///<!
   TBranch        *b_track2StartDirY;   ///<!
   TBranch        *b_track2StartDirZ;   ///<!
   TBranch        *b_track3StartDirX;   ///<!
   TBranch        *b_track3StartDirY;   ///<!
   TBranch        *b_track3StartDirZ;   ///<!
   TBranch        *b_shower1DirX;   ///<!
   TBranch        *b_shower1DirY;   ///<!
   TBranch        *b_shower1DirZ;   ///<!
   TBranch        *b_track1DistRecoVtx;   ///<!
   TBranch        *b_track2DistRecoVtx;   ///<!
   TBranch        *b_track3DistRecoVtx;   ///<!
   TBranch        *b_shower1DistRecoVtx;   ///<!
   TBranch        *b_track1Track2Angle;   ///<!
   TBranch        *b_track1Track3Angle;   ///<!
   TBranch        *b_track2Track3Angle;   ///<!
   TBranch        *b_track1Shower1Angle;   ///<!
   TBranch        *b_track2Shower1Angle;   ///<!
   TBranch        *b_track3Shower1Angle;   ///<!
   TBranch        *b_track1Track2Dist;   ///<!
   TBranch        *b_track1Track3Dist;   ///<!
   TBranch        *b_track2Track3Dist;   ///<!
   TBranch        *b_track1Shower1Dist;   ///<!
   TBranch        *b_track2Shower1Dist;   ///<!
   TBranch        *b_track3Shower1Dist;   ///<!

   backgroundPlots(TTree *tree=0);
   virtual ~backgroundPlots();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef backgroundPlots_cxx
backgroundPlots::backgroundPlots(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/ooconnor/sbnd/hyperons/preselection_output/tmvaSample_bkg.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/ooconnor/sbnd/hyperons/preselection_output/tmvaSample_bkg.root");
      }
      f->GetObject("bkgTree",tree);

   }
   Init(tree);
}

backgroundPlots::~backgroundPlots()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t backgroundPlots::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t backgroundPlots::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void backgroundPlots::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.

   // Set object pointer
   trueP = 0;
   truePDG = 0;
   daughterPDG = 0;
   motherPDG = 0;
   vertexX = 0;
   vertexY = 0;
   vertexZ = 0;
   vertexSize = 0;
   daughterSize = 0;
   trueNuEnergy = 0;
   trueCCNC = 0;
   trueIntMode = 0;
   trueIntType = 0;
   TrackIDs = 0;
   trackLengths = 0;
   DistanceToRecoVertex = 0;
   nuScores = 0;
   trackScores = 0;
   muonTrackScores = 0;
   protonTrackScores = 0;
   pionTrackScores = 0;
   NeutrinoNuScores = 0;
   trackStartPositionX = 0;
   trackStartPositionY = 0;
   trackStartPositionZ = 0;
   trackEndPositionX = 0;
   trackEndPositionY = 0;
   trackEndPositionZ = 0;
   trackStartDirX = 0;
   trackStartDirY = 0;
   trackStartDirZ = 0;
   trackEndDirX = 0;
   trackEndDirY = 0;
   trackEndDirZ = 0;
   trackVertexDirX = 0;
   trackVertexDirY = 0;
   trackVertexDirZ = 0;
   trackTheta = 0;
   trackPhi = 0;
   showerLengths = 0;
   showerStartPositionX = 0;
   showerStartPositionY = 0;
   showerStartPositionZ = 0;
   showerDirX = 0;
   showerDirY = 0;
   showerDirZ = 0;
   pfpPDG = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("trueP", &trueP, &b_trueP);
   fChain->SetBranchAddress("truePDG", &truePDG, &b_truePDG);
   fChain->SetBranchAddress("daughterPDG", &daughterPDG, &b_daughterPDG);
   fChain->SetBranchAddress("motherPDG", &motherPDG, &b_motherPDG);
   fChain->SetBranchAddress("vertexX", &vertexX, &b_vertexX);
   fChain->SetBranchAddress("vertexY", &vertexY, &b_vertexY);
   fChain->SetBranchAddress("vertexZ", &vertexZ, &b_vertexZ);
   fChain->SetBranchAddress("vertexSize", &vertexSize, &b_vertexSize);
   fChain->SetBranchAddress("daughterSize", &daughterSize, &b_daughterSize);
   fChain->SetBranchAddress("trueNuEnergy", &trueNuEnergy, &b_trueNuEnergy);
   fChain->SetBranchAddress("trueCCNC", &trueCCNC, &b_trueCCNC);
   fChain->SetBranchAddress("trueIntMode", &trueIntMode, &b_trueIntMode);
   fChain->SetBranchAddress("trueIntType", &trueIntType, &b_trueIntType);
   fChain->SetBranchAddress("nPFParticles", &nPFParticles, &b_nPFParticles);
   fChain->SetBranchAddress("nPrimaryChildren", &nPrimaryChildren, &b_nPrimaryChildren);
   fChain->SetBranchAddress("trackCount", &trackCount, &b_trackCount);
   fChain->SetBranchAddress("showerCount", &showerCount, &b_showerCount);
   fChain->SetBranchAddress("TrackIDs", &TrackIDs, &b_TrackIDs);
   fChain->SetBranchAddress("trackLengths", &trackLengths, &b_trackLengths);
   fChain->SetBranchAddress("RecoVertexX", &RecoVertexX, &b_RecoVertexX);
   fChain->SetBranchAddress("RecoVertexY", &RecoVertexY, &b_RecoVertexY);
   fChain->SetBranchAddress("RecoVertexZ", &RecoVertexZ, &b_RecoVertexZ);
   fChain->SetBranchAddress("DistanceToRecoVertex", &DistanceToRecoVertex, &b_DistanceToRecoVertex);
   fChain->SetBranchAddress("nuScores", &nuScores, &b_nuScores);
   fChain->SetBranchAddress("trackScores", &trackScores, &b_trackScores);
   fChain->SetBranchAddress("muonTrackScores", &muonTrackScores, &b_muonTrackScores);
   fChain->SetBranchAddress("protonTrackScores", &protonTrackScores, &b_protonTrackScores);
   fChain->SetBranchAddress("pionTrackScores", &pionTrackScores, &b_pionTrackScores);
   fChain->SetBranchAddress("NeutrinoNuScores", &NeutrinoNuScores, &b_NeutrinoNuScores);
   fChain->SetBranchAddress("trackStartPositionX", &trackStartPositionX, &b_trackStartPositionX);
   fChain->SetBranchAddress("trackStartPositionY", &trackStartPositionY, &b_trackStartPositionY);
   fChain->SetBranchAddress("trackStartPositionZ", &trackStartPositionZ, &b_trackStartPositionZ);
   fChain->SetBranchAddress("trackEndPositionX", &trackEndPositionX, &b_trackEndPositionX);
   fChain->SetBranchAddress("trackEndPositionY", &trackEndPositionY, &b_trackEndPositionY);
   fChain->SetBranchAddress("trackEndPositionZ", &trackEndPositionZ, &b_trackEndPositionZ);
   fChain->SetBranchAddress("trackStartDirX", &trackStartDirX, &b_trackStartDirX);
   fChain->SetBranchAddress("trackStartDirY", &trackStartDirY, &b_trackStartDirY);
   fChain->SetBranchAddress("trackStartDirZ", &trackStartDirZ, &b_trackStartDirZ);
   fChain->SetBranchAddress("trackEndDirX", &trackEndDirX, &b_trackEndDirX);
   fChain->SetBranchAddress("trackEndDirY", &trackEndDirY, &b_trackEndDirY);
   fChain->SetBranchAddress("trackEndDirZ", &trackEndDirZ, &b_trackEndDirZ);
   fChain->SetBranchAddress("trackVertexDirX", &trackVertexDirX, &b_trackVertexDirX);
   fChain->SetBranchAddress("trackVertexDirY", &trackVertexDirY, &b_trackVertexDirY);
   fChain->SetBranchAddress("trackVertexDirZ", &trackVertexDirZ, &b_trackVertexDirZ);
   fChain->SetBranchAddress("trackTheta", &trackTheta, &b_trackTheta);
   fChain->SetBranchAddress("trackPhi", &trackPhi, &b_trackPhi);
   fChain->SetBranchAddress("showerLengths", &showerLengths, &b_showerLengths);
   fChain->SetBranchAddress("showerStartPositionX", &showerStartPositionX, &b_showerStartPositionX);
   fChain->SetBranchAddress("showerStartPositionY", &showerStartPositionY, &b_showerStartPositionY);
   fChain->SetBranchAddress("showerStartPositionZ", &showerStartPositionZ, &b_showerStartPositionZ);
   fChain->SetBranchAddress("showerDirX", &showerDirX, &b_showerDirX);
   fChain->SetBranchAddress("showerDirY", &showerDirY, &b_showerDirY);
   fChain->SetBranchAddress("showerDirZ", &showerDirZ, &b_showerDirZ);
   fChain->SetBranchAddress("pfpPDG", &pfpPDG, &b_pfpPDG);
   fChain->SetBranchAddress("sampleType", &sampleType, &b_sampleType);
   fChain->SetBranchAddress("chosenMCTruthIdx", &chosenMCTruthIdx, &b_chosenMCTruthIdx);
   fChain->SetBranchAddress("track1Length", &track1Length, &b_track1Length);
   fChain->SetBranchAddress("track2Length", &track2Length, &b_track2Length);
   fChain->SetBranchAddress("track3Length", &track3Length, &b_track3Length);
   fChain->SetBranchAddress("shower1Length", &shower1Length, &b_shower1Length);
   fChain->SetBranchAddress("track1StartPosX", &track1StartPosX, &b_track1StartPosX);
   fChain->SetBranchAddress("track1StartPosY", &track1StartPosY, &b_track1StartPosY);
   fChain->SetBranchAddress("track1StartPosZ", &track1StartPosZ, &b_track1StartPosZ);
   fChain->SetBranchAddress("track2StartPosX", &track2StartPosX, &b_track2StartPosX);
   fChain->SetBranchAddress("track2StartPosY", &track2StartPosY, &b_track2StartPosY);
   fChain->SetBranchAddress("track2StartPosZ", &track2StartPosZ, &b_track2StartPosZ);
   fChain->SetBranchAddress("track3StartPosX", &track3StartPosX, &b_track3StartPosX);
   fChain->SetBranchAddress("track3StartPosY", &track3StartPosY, &b_track3StartPosY);
   fChain->SetBranchAddress("track3StartPosZ", &track3StartPosZ, &b_track3StartPosZ);
   fChain->SetBranchAddress("shower1StartPosX", &shower1StartPosX, &b_shower1StartPosX);
   fChain->SetBranchAddress("shower1StartPosY", &shower1StartPosY, &b_shower1StartPosY);
   fChain->SetBranchAddress("shower1StartPosZ", &shower1StartPosZ, &b_shower1StartPosZ);
   fChain->SetBranchAddress("track1StartDirX", &track1StartDirX, &b_track1StartDirX);
   fChain->SetBranchAddress("track1StartDirY", &track1StartDirY, &b_track1StartDirY);
   fChain->SetBranchAddress("track1StartDirZ", &track1StartDirZ, &b_track1StartDirZ);
   fChain->SetBranchAddress("track2StartDirX", &track2StartDirX, &b_track2StartDirX);
   fChain->SetBranchAddress("track2StartDirY", &track2StartDirY, &b_track2StartDirY);
   fChain->SetBranchAddress("track2StartDirZ", &track2StartDirZ, &b_track2StartDirZ);
   fChain->SetBranchAddress("track3StartDirX", &track3StartDirX, &b_track3StartDirX);
   fChain->SetBranchAddress("track3StartDirY", &track3StartDirY, &b_track3StartDirY);
   fChain->SetBranchAddress("track3StartDirZ", &track3StartDirZ, &b_track3StartDirZ);
   fChain->SetBranchAddress("shower1DirX", &shower1DirX, &b_shower1DirX);
   fChain->SetBranchAddress("shower1DirY", &shower1DirY, &b_shower1DirY);
   fChain->SetBranchAddress("shower1DirZ", &shower1DirZ, &b_shower1DirZ);
   fChain->SetBranchAddress("track1DistRecoVtx", &track1DistRecoVtx, &b_track1DistRecoVtx);
   fChain->SetBranchAddress("track2DistRecoVtx", &track2DistRecoVtx, &b_track2DistRecoVtx);
   fChain->SetBranchAddress("track3DistRecoVtx", &track3DistRecoVtx, &b_track3DistRecoVtx);
   fChain->SetBranchAddress("shower1DistRecoVtx", &shower1DistRecoVtx, &b_shower1DistRecoVtx);
   fChain->SetBranchAddress("track1Track2Angle", &track1Track2Angle, &b_track1Track2Angle);
   fChain->SetBranchAddress("track1Track3Angle", &track1Track3Angle, &b_track1Track3Angle);
   fChain->SetBranchAddress("track2Track3Angle", &track2Track3Angle, &b_track2Track3Angle);
   fChain->SetBranchAddress("track1Shower1Angle", &track1Shower1Angle, &b_track1Shower1Angle);
   fChain->SetBranchAddress("track2Shower1Angle", &track2Shower1Angle, &b_track2Shower1Angle);
   fChain->SetBranchAddress("track3Shower1Angle", &track3Shower1Angle, &b_track3Shower1Angle);
   fChain->SetBranchAddress("track1Track2Dist", &track1Track2Dist, &b_track1Track2Dist);
   fChain->SetBranchAddress("track1Track3Dist", &track1Track3Dist, &b_track1Track3Dist);
   fChain->SetBranchAddress("track2Track3Dist", &track2Track3Dist, &b_track2Track3Dist);
   fChain->SetBranchAddress("track1Shower1Dist", &track1Shower1Dist, &b_track1Shower1Dist);
   fChain->SetBranchAddress("track2Shower1Dist", &track2Shower1Dist, &b_track2Shower1Dist);
   fChain->SetBranchAddress("track3Shower1Dist", &track3Shower1Dist, &b_track3Shower1Dist);
   Notify();
}

bool backgroundPlots::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be for a new TTree in a TChain. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void backgroundPlots::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t backgroundPlots::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef backgroundPlots_cxx
