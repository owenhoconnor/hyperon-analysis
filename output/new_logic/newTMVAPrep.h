//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 12 19:41:40 2026 by ROOT version 6.40.02
// from TTree tree/Output TTree
// found on file: /data/ooconnor/sbnd/hyperons/preselection_output/new_logic/signalDef_output_bkg.root
//////////////////////////////////////////////////////////

#ifndef newTMVAPrep_h
#define newTMVAPrep_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class newTMVAPrep {
public :
   TTree          *fChain;   ///<!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; ///<!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          eventID;
   Int_t           run;
   Int_t           subrun;
   vector<int> *trueOrigin;
   vector<float> *trueW;
   vector<float> *trueX;
   vector<float> *trueY;
   vector<float> *trueQSqr;
   vector<float> *truePt;
   vector<float> *trueTheta;
   vector<int>     *trueNuPDG;
   vector<int>     *trueNuTrackID;
   vector<float>   *trueNuVtxX;
   vector<float>   *trueNuVtxY;
   vector<float>   *trueNuVtxZ;
   vector<float>     *trueNuEnergy;
   vector<int>     *trueCCNC;
   vector<int>     *trueIntMode;
   vector<int>     *trueIntType;
   vector<int>     *trueTargetPDG;
   vector<int>     *truePDG;
   vector<int>     *trueTrackID;
   vector<int>     *trueMotherPDG;
   vector<int>     *trueMotherTrackID;
   vector<int>     *trueGeneration;
   vector<int>     *trueIsPrimary;
   vector<int>     *trueIsDecayProduct;
   vector<int>     *trueMCTruthIndex;
   vector<float>   *trueP;
   vector<float>   *trueMass;
   vector<float>   *trueStartX;
   vector<float>   *trueStartY;
   vector<float>   *trueStartZ;
   vector<float>   *trueStartT;
   vector<float>   *truePx;
   vector<float>   *truePy;
   vector<float>   *truePz;
   vector<float>   *trueE;
   vector<float>   *trueEndX;
   vector<float>   *trueEndY;
   vector<float>   *trueEndZ;
   vector<float>   *trueEndT;
   vector<float>   *trueEndPx;
   vector<float>   *trueEndPy;
   vector<float>   *trueEndPz;
   vector<float>   *trueEndE;
   vector<int>     *trueNTrajectoryPoints;
   vector<float>   *trueTrajectoryLength;
   vector<string>  *trueProcess;
   vector<string>  *trueEndProcess;
   vector<int>     *trueStatusCode;
   vector<int>     *trueNDaughters;
   vector<int>     *trueNStoredDecayDaughters;
   vector<int>     *trueNPrimaryParticles;
   vector<int>     *trueNSavedParticles;
   vector<int>     *trueParticleStartIndex;
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
   Int_t           chosenTruthIdx;

   // List of branches
   TBranch        *b_eventID;   ///<!
   TBranch        *b_run;   ///<!
   TBranch        *b_subrun;   ///<!
   TBranch        *b_trueOrigin; ///<!
   TBranch        *b_trueW;      ///<!
   TBranch        *b_trueX;      ///<!
   TBranch        *b_trueY;      ///<!
   TBranch        *b_trueQSqr;   ///<!
   TBranch        *b_truePt;     ///<!
   TBranch        *b_trueTheta;  ///<!
   TBranch        *b_trueNuPDG;   ///<!
   TBranch        *b_trueNuTrackID;   ///<!
   TBranch        *b_trueNuVtxX;   ///<!
   TBranch        *b_trueNuVtxY;   ///<!
   TBranch        *b_trueNuVtxZ;   ///<!
   TBranch        *b_trueNuEnergy;   ///<!
   TBranch        *b_trueCCNC;   ///<!
   TBranch        *b_trueIntMode;   ///<!
   TBranch        *b_trueIntType;   ///<!
   TBranch        *b_trueTargetPDG;   ///<!
   TBranch        *b_truePDG;   ///<!
   TBranch        *b_trueTrackID;   ///<!
   TBranch        *b_trueMotherPDG;   ///<!
   TBranch        *b_trueMotherTrackID;   ///<!
   TBranch        *b_trueGeneration;   ///<!
   TBranch        *b_trueIsPrimary;   ///<!
   TBranch        *b_trueIsDecayProduct;   ///<!
   TBranch        *b_trueMCTruthIndex;   ///<!
   TBranch        *b_trueP;   ///<!
   TBranch        *b_trueMass;   ///<!
   TBranch        *b_trueStartX;   ///<!
   TBranch        *b_trueStartY;   ///<!
   TBranch        *b_trueStartZ;   ///<!
   TBranch        *b_trueStartT;   ///<!
   TBranch        *b_truePx;   ///<!
   TBranch        *b_truePy;   ///<!
   TBranch        *b_truePz;   ///<!
   TBranch        *b_trueE;   ///<!
   TBranch        *b_trueEndX;   ///<!
   TBranch        *b_trueEndY;   ///<!
   TBranch        *b_trueEndZ;   ///<!
   TBranch        *b_trueEndT;   ///<!
   TBranch        *b_trueEndPx;   ///<!
   TBranch        *b_trueEndPy;   ///<!
   TBranch        *b_trueEndPz;   ///<!
   TBranch        *b_trueEndE;   ///<!
   TBranch        *b_trueNTrajectoryPoints;   ///<!
   TBranch        *b_trueTrajectoryLength;   ///<!
   TBranch        *b_trueProcess;   ///<!
   TBranch        *b_trueEndProcess;   ///<!
   TBranch        *b_trueStatusCode;   ///<!
   TBranch        *b_trueNDaughters;   ///<!
   TBranch        *b_trueNStoredDecayDaughters;   ///<!
   TBranch        *b_trueNPrimaryParticles;   ///<!
   TBranch        *b_trueNSavedParticles;   ///<!
   TBranch        *b_trueParticleStartIndex;   ///<!
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
   TBranch        *b_chosenTruthIdx;   ///<!

   newTMVAPrep(TTree *tree=0);
   virtual ~newTMVAPrep();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef newTMVAPrep_cxx
newTMVAPrep::newTMVAPrep(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/ooconnor/sbnd/hyperons/preselection_output/new_logic/signalDef_output_bkg.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/ooconnor/sbnd/hyperons/preselection_output/new_logic/signalDef_output_bkg.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

newTMVAPrep::~newTMVAPrep()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t newTMVAPrep::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t newTMVAPrep::LoadTree(Long64_t entry)
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

void newTMVAPrep::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.

   // Set object pointer
   trueOrigin = 0;
   trueW = 0;
   trueX = 0;
   trueY = 0;
   trueQSqr = 0;
   truePt = 0;
   trueTheta = 0;
   trueNuPDG = 0;
   trueNuTrackID = 0;
   trueNuVtxX = 0;
   trueNuVtxY = 0;
   trueNuVtxZ = 0;
   trueNuEnergy = 0;
   trueCCNC = 0;
   trueIntMode = 0;
   trueIntType = 0;
   trueTargetPDG = 0;
   truePDG = 0;
   trueTrackID = 0;
   trueMotherPDG = 0;
   trueMotherTrackID = 0;
   trueGeneration = 0;
   trueIsPrimary = 0;
   trueIsDecayProduct = 0;
   trueMCTruthIndex = 0;
   trueP = 0;
   trueMass = 0;
   trueStartX = 0;
   trueStartY = 0;
   trueStartZ = 0;
   trueStartT = 0;
   truePx = 0;
   truePy = 0;
   truePz = 0;
   trueE = 0;
   trueEndX = 0;
   trueEndY = 0;
   trueEndZ = 0;
   trueEndT = 0;
   trueEndPx = 0;
   trueEndPy = 0;
   trueEndPz = 0;
   trueEndE = 0;
   trueNTrajectoryPoints = 0;
   trueTrajectoryLength = 0;
   trueProcess = 0;
   trueEndProcess = 0;
   trueStatusCode = 0;
   trueNDaughters = 0;
   trueNStoredDecayDaughters = 0;
   trueNPrimaryParticles = 0;
   trueNSavedParticles = 0;
   trueParticleStartIndex = 0;
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
   fChain->SetBranchAddress("trueNuPDG", &trueNuPDG, &b_trueNuPDG);
   fChain->SetBranchAddress("trueOrigin". &trueOrigin, &b_trueOrigin);
   fChain->SetBranchAddress("trueW", &trueW, &b_trueW);
   fChain->SetBranchAddress("trueX", &trueX, &b_trueX);
   fChain->SetBranchAddress("trueY", &trueY, &b_trueY);
   fChain->SetBranchAddress("trueQSqr", &trueQSqr, &b_trueQSqr);
   fChain->SetBranchAddress("truePt", &truePt, &b_truePt);
   fChain->SetBranchAddress("trueTheta", &trueTheta, &b_trueTheta);
   fChain->SetBranchAddress("trueNuTrackID", &trueNuTrackID, &b_trueNuTrackID);
   fChain->SetBranchAddress("trueNuVtxX", &trueNuVtxX, &b_trueNuVtxX);
   fChain->SetBranchAddress("trueNuVtxY", &trueNuVtxY, &b_trueNuVtxY);
   fChain->SetBranchAddress("trueNuVtxZ", &trueNuVtxZ, &b_trueNuVtxZ);
   fChain->SetBranchAddress("trueNuEnergy", &trueNuEnergy, &b_trueNuEnergy);
   fChain->SetBranchAddress("trueCCNC", &trueCCNC, &b_trueCCNC);
   fChain->SetBranchAddress("trueIntMode", &trueIntMode, &b_trueIntMode);
   fChain->SetBranchAddress("trueIntType", &trueIntType, &b_trueIntType);
   fChain->SetBranchAddress("trueTargetPDG", &trueTargetPDG, &b_trueTargetPDG);
   fChain->SetBranchAddress("truePDG", &truePDG, &b_truePDG);
   fChain->SetBranchAddress("trueTrackID", &trueTrackID, &b_trueTrackID);
   fChain->SetBranchAddress("trueMotherPDG", &trueMotherPDG, &b_trueMotherPDG);
   fChain->SetBranchAddress("trueMotherTrackID", &trueMotherTrackID, &b_trueMotherTrackID);
   fChain->SetBranchAddress("trueGeneration", &trueGeneration, &b_trueGeneration);
   fChain->SetBranchAddress("trueIsPrimary", &trueIsPrimary, &b_trueIsPrimary);
   fChain->SetBranchAddress("trueIsDecayProduct", &trueIsDecayProduct, &b_trueIsDecayProduct);
   fChain->SetBranchAddress("trueMCTruthIndex", &trueMCTruthIndex, &b_trueMCTruthIndex);
   fChain->SetBranchAddress("trueP", &trueP, &b_trueP);
   fChain->SetBranchAddress("trueMass", &trueMass, &b_trueMass);
   fChain->SetBranchAddress("trueStartX", &trueStartX, &b_trueStartX);
   fChain->SetBranchAddress("trueStartY", &trueStartY, &b_trueStartY);
   fChain->SetBranchAddress("trueStartZ", &trueStartZ, &b_trueStartZ);
   fChain->SetBranchAddress("trueStartT", &trueStartT, &b_trueStartT);
   fChain->SetBranchAddress("truePx", &truePx, &b_truePx);
   fChain->SetBranchAddress("truePy", &truePy, &b_truePy);
   fChain->SetBranchAddress("truePz", &truePz, &b_truePz);
   fChain->SetBranchAddress("trueE", &trueE, &b_trueE);
   fChain->SetBranchAddress("trueEndX", &trueEndX, &b_trueEndX);
   fChain->SetBranchAddress("trueEndY", &trueEndY, &b_trueEndY);
   fChain->SetBranchAddress("trueEndZ", &trueEndZ, &b_trueEndZ);
   fChain->SetBranchAddress("trueEndT", &trueEndT, &b_trueEndT);
   fChain->SetBranchAddress("trueEndPx", &trueEndPx, &b_trueEndPx);
   fChain->SetBranchAddress("trueEndPy", &trueEndPy, &b_trueEndPy);
   fChain->SetBranchAddress("trueEndPz", &trueEndPz, &b_trueEndPz);
   fChain->SetBranchAddress("trueEndE", &trueEndE, &b_trueEndE);
   fChain->SetBranchAddress("trueNTrajectoryPoints", &trueNTrajectoryPoints, &b_trueNTrajectoryPoints);
   fChain->SetBranchAddress("trueTrajectoryLength", &trueTrajectoryLength, &b_trueTrajectoryLength);
   fChain->SetBranchAddress("trueProcess", &trueProcess, &b_trueProcess);
   fChain->SetBranchAddress("trueEndProcess", &trueEndProcess, &b_trueEndProcess);
   fChain->SetBranchAddress("trueStatusCode", &trueStatusCode, &b_trueStatusCode);
   fChain->SetBranchAddress("trueNDaughters", &trueNDaughters, &b_trueNDaughters);
   fChain->SetBranchAddress("trueNStoredDecayDaughters", &trueNStoredDecayDaughters, &b_trueNStoredDecayDaughters);
   fChain->SetBranchAddress("trueNPrimaryParticles", &trueNPrimaryParticles, &b_trueNPrimaryParticles);
   fChain->SetBranchAddress("trueNSavedParticles", &trueNSavedParticles, &b_trueNSavedParticles);
   fChain->SetBranchAddress("trueParticleStartIndex", &trueParticleStartIndex, &b_trueParticleStartIndex);
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
   fChain->SetBranchAddress("chosenTruthIdx", &chosenTruthIdx, &b_chosenTruthIdx);
   Notify();
}

bool newTMVAPrep::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be for a new TTree in a TChain. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void newTMVAPrep::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t newTMVAPrep::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef newTMVAPrep_cxx
