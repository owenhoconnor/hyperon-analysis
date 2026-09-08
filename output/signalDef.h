//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep  8 11:13:28 2026 by ROOT version 6.40.02
// from TTree tree/Output TTree
// found on file: /data/ooconnor/sbnd/hyperons/analyzer_output/merged_anaOut_hyperons.root
//////////////////////////////////////////////////////////

#ifndef signalDef_h
#define signalDef_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class signalDef {
public :
   TTree          *fChain;   ///<!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; ///<!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          eventID;
   Int_t           run;
   Int_t           subrun;
   vector<int>     *trueOrigin;
   vector<float>   *trueW;
   vector<float>   *trueX;
   vector<float>   *trueY;
   vector<float>   *trueQSqr;
   vector<float>   *truePt;
   vector<float>   *trueTheta;
   vector<int>     *trueNuPDG;
   vector<int>     *trueNuTrackID;
   vector<float>   *trueNuVtxX;
   vector<float>   *trueNuVtxY;
   vector<float>   *trueNuVtxZ;
   vector<float>   *trueNuEnergy;
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
   vector<int>     *sliceKey;
   vector<int>     *sliceID;
   vector<float>   *sliceNuScore;
   vector<int>     *sliceTotalHits;
   vector<int>     *sliceTrueNuHits;
   vector<int>     *sliceTrueOrigin;
   Int_t           eventTotalTrueNuHits;
   vector<float>   *sliceVtxX;
   vector<float>   *sliceVtxY;
   vector<float>   *sliceVtxZ;
   vector<int>     *pfpKey;
   vector<int>     *pfpSelfID;
   vector<int>     *pfpParentID;
   vector<int>     *pfpRecoPDG;
   vector<int>     *pfpSliceKey;
   vector<int>     *pfpIsNuSlice;
   vector<int>     *pfpIsPrimary;
   vector<int>     *pfpNPrimaryChildren;
   vector<float>   *pfpTrackScore;
   vector<int>     *pfpHasTrackScore;
   vector<int>     *pfpNTracks;
   vector<int>     *pfpNShowers;
   vector<int>     *pfpHasTrack;
   vector<int>     *pfpHasShower;
   vector<int>     *pfpHasUniqueTrack;
   vector<int>     *pfpHasUniqueShower;
   vector<int>     *pfpNVertices;
   vector<int>     *pfpHasVertex;
   vector<int>     *pfpHasUniqueVertex;
   vector<int>     *pfpTrackID;
   vector<float>   *pfpTrackLength;
   vector<float>   *pfpTrackStartX;
   vector<float>   *pfpTrackStartY;
   vector<float>   *pfpTrackStartZ;
   vector<float>   *pfpTrackEndX;
   vector<float>   *pfpTrackEndY;
   vector<float>   *pfpTrackEndZ;
   vector<float>   *pfpTrackStartDirX;
   vector<float>   *pfpTrackStartDirY;
   vector<float>   *pfpTrackStartDirZ;
   vector<float>   *pfpTrackEndDirX;
   vector<float>   *pfpTrackEndDirY;
   vector<float>   *pfpTrackEndDirZ;
   vector<float>   *pfpTrackVertexDirX;
   vector<float>   *pfpTrackVertexDirY;
   vector<float>   *pfpTrackVertexDirZ;
   vector<float>   *pfpTrackTheta;
   vector<float>   *pfpTrackPhi;
   vector<int>     *pfpShowerID;
   vector<float>   *pfpShowerLength;
   vector<float>   *pfpShowerStartX;
   vector<float>   *pfpShowerStartY;
   vector<float>   *pfpShowerStartZ;
   vector<float>   *pfpShowerDirX;
   vector<float>   *pfpShowerDirY;
   vector<float>   *pfpShowerDirZ;
   vector<float>   *pfpVertexX;
   vector<float>   *pfpVertexY;
   vector<float>   *pfpVertexZ;
   vector<int>     *pfpTrueTrackID;
   vector<int>     *pfpTruePDG;
   vector<int>     *pfpNHits;
   vector<int>     *pfpNMatchedHits;
   vector<float>   *pfpTruthPurity;
   vector<int>     *trueIsReconstructed;
   vector<int>     *trueIsReconstructedInNuSlice;
   vector<int>     *trueNMatchedPfps;
   vector<int>     *trueBestRecoPfpIdx;
   vector<float>   *trueBestRecoTrackScore;
   vector<int>     *trueBestRecoHasTrack;
   vector<int>     *trueBestRecoHasShower;
   vector<float>   *nuScores;
   vector<float>   *NeutrinoNuScores;

   // List of branches
   TBranch        *b_eventID;   ///<!
   TBranch        *b_run;   ///<!
   TBranch        *b_subrun;   ///<!
   TBranch        *b_trueOrigin;   ///<!
   TBranch        *b_trueW;   ///<!
   TBranch        *b_trueX;   ///<!
   TBranch        *b_trueY;   ///<!
   TBranch        *b_trueQSqr;   ///<!
   TBranch        *b_truePt;   ///<!
   TBranch        *b_trueTheta;   ///<!
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
   TBranch        *b_sliceKey;   ///<!
   TBranch        *b_sliceID;   ///<!
   TBranch        *b_sliceNuScore;   ///<!
   TBranch        *b_sliceTotalHits;   ///<!
   TBranch        *b_sliceTrueNuHits;   ///<!
   TBranch        *b_sliceTrueOrigin;   ///<!
   TBranch        *b_eventTotalTrueNuHits;   ///<!
   TBranch        *b_sliceVtxX;   ///<!
   TBranch        *b_sliceVtxY;   ///<!
   TBranch        *b_sliceVtxZ;   ///<!
   TBranch        *b_pfpKey;   ///<!
   TBranch        *b_pfpSelfID;   ///<!
   TBranch        *b_pfpParentID;   ///<!
   TBranch        *b_pfpRecoPDG;   ///<!
   TBranch        *b_pfpSliceKey;   ///<!
   TBranch        *b_pfpIsNuSlice;   ///<!
   TBranch        *b_pfpIsPrimary;   ///<!
   TBranch        *b_pfpNPrimaryChildren;   ///<!
   TBranch        *b_pfpTrackScore;   ///<!
   TBranch        *b_pfpHasTrackScore;   ///<!
   TBranch        *b_pfpNTracks;   ///<!
   TBranch        *b_pfpNShowers;   ///<!
   TBranch        *b_pfpHasTrack;   ///<!
   TBranch        *b_pfpHasShower;   ///<!
   TBranch        *b_pfpHasUniqueTrack;   ///<!
   TBranch        *b_pfpHasUniqueShower;   ///<!
   TBranch        *b_pfpNVertices;   ///<!
   TBranch        *b_pfpHasVertex;   ///<!
   TBranch        *b_pfpHasUniqueVertex;   ///<!
   TBranch        *b_pfpTrackID;   ///<!
   TBranch        *b_pfpTrackLength;   ///<!
   TBranch        *b_pfpTrackStartX;   ///<!
   TBranch        *b_pfpTrackStartY;   ///<!
   TBranch        *b_pfpTrackStartZ;   ///<!
   TBranch        *b_pfpTrackEndX;   ///<!
   TBranch        *b_pfpTrackEndY;   ///<!
   TBranch        *b_pfpTrackEndZ;   ///<!
   TBranch        *b_pfpTrackStartDirX;   ///<!
   TBranch        *b_pfpTrackStartDirY;   ///<!
   TBranch        *b_pfpTrackStartDirZ;   ///<!
   TBranch        *b_pfpTrackEndDirX;   ///<!
   TBranch        *b_pfpTrackEndDirY;   ///<!
   TBranch        *b_pfpTrackEndDirZ;   ///<!
   TBranch        *b_pfpTrackVertexDirX;   ///<!
   TBranch        *b_pfpTrackVertexDirY;   ///<!
   TBranch        *b_pfpTrackVertexDirZ;   ///<!
   TBranch        *b_pfpTrackTheta;   ///<!
   TBranch        *b_pfpTrackPhi;   ///<!
   TBranch        *b_pfpShowerID;   ///<!
   TBranch        *b_pfpShowerLength;   ///<!
   TBranch        *b_pfpShowerStartX;   ///<!
   TBranch        *b_pfpShowerStartY;   ///<!
   TBranch        *b_pfpShowerStartZ;   ///<!
   TBranch        *b_pfpShowerDirX;   ///<!
   TBranch        *b_pfpShowerDirY;   ///<!
   TBranch        *b_pfpShowerDirZ;   ///<!
   TBranch        *b_pfpVertexX;   ///<!
   TBranch        *b_pfpVertexY;   ///<!
   TBranch        *b_pfpVertexZ;   ///<!
   TBranch        *b_pfpTrueTrackID;   ///<!
   TBranch        *b_pfpTruePDG;   ///<!
   TBranch        *b_pfpNHits;   ///<!
   TBranch        *b_pfpNMatchedHits;   ///<!
   TBranch        *b_pfpTruthPurity;   ///<!
   TBranch        *b_trueIsReconstructed;   ///<!
   TBranch        *b_trueIsReconstructedInNuSlice;   ///<!
   TBranch        *b_trueNMatchedPfps;   ///<!
   TBranch        *b_trueBestRecoPfpIdx;   ///<!
   TBranch        *b_trueBestRecoTrackScore;   ///<!
   TBranch        *b_trueBestRecoHasTrack;   ///<!
   TBranch        *b_trueBestRecoHasShower;   ///<!
   TBranch        *b_nuScores;   ///<!
   TBranch        *b_NeutrinoNuScores;   ///<!

   signalDef(TTree *tree=0);
   virtual ~signalDef();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef signalDef_cxx
signalDef::signalDef(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/ooconnor/sbnd/hyperons/analyzer_output/merged_anaOut_hyperons.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/ooconnor/sbnd/hyperons/analyzer_output/merged_anaOut_hyperons.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/data/ooconnor/sbnd/hyperons/analyzer_output/merged_anaOut_hyperons.root:/ana");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

signalDef::~signalDef()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t signalDef::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t signalDef::LoadTree(Long64_t entry)
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

void signalDef::Init(TTree *tree)
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
   sliceKey = 0;
   sliceID = 0;
   sliceNuScore = 0;
   sliceTotalHits = 0;
   sliceTrueNuHits = 0;
   sliceTrueOrigin = 0;
   sliceVtxX = 0;
   sliceVtxY = 0;
   sliceVtxZ = 0;
   pfpKey = 0;
   pfpSelfID = 0;
   pfpParentID = 0;
   pfpRecoPDG = 0;
   pfpSliceKey = 0;
   pfpIsNuSlice = 0;
   pfpIsPrimary = 0;
   pfpNPrimaryChildren = 0;
   pfpTrackScore = 0;
   pfpHasTrackScore = 0;
   pfpNTracks = 0;
   pfpNShowers = 0;
   pfpHasTrack = 0;
   pfpHasShower = 0;
   pfpHasUniqueTrack = 0;
   pfpHasUniqueShower = 0;
   pfpNVertices = 0;
   pfpHasVertex = 0;
   pfpHasUniqueVertex = 0;
   pfpTrackID = 0;
   pfpTrackLength = 0;
   pfpTrackStartX = 0;
   pfpTrackStartY = 0;
   pfpTrackStartZ = 0;
   pfpTrackEndX = 0;
   pfpTrackEndY = 0;
   pfpTrackEndZ = 0;
   pfpTrackStartDirX = 0;
   pfpTrackStartDirY = 0;
   pfpTrackStartDirZ = 0;
   pfpTrackEndDirX = 0;
   pfpTrackEndDirY = 0;
   pfpTrackEndDirZ = 0;
   pfpTrackVertexDirX = 0;
   pfpTrackVertexDirY = 0;
   pfpTrackVertexDirZ = 0;
   pfpTrackTheta = 0;
   pfpTrackPhi = 0;
   pfpShowerID = 0;
   pfpShowerLength = 0;
   pfpShowerStartX = 0;
   pfpShowerStartY = 0;
   pfpShowerStartZ = 0;
   pfpShowerDirX = 0;
   pfpShowerDirY = 0;
   pfpShowerDirZ = 0;
   pfpVertexX = 0;
   pfpVertexY = 0;
   pfpVertexZ = 0;
   pfpTrueTrackID = 0;
   pfpTruePDG = 0;
   pfpNHits = 0;
   pfpNMatchedHits = 0;
   pfpTruthPurity = 0;
   trueIsReconstructed = 0;
   trueIsReconstructedInNuSlice = 0;
   trueNMatchedPfps = 0;
   trueBestRecoPfpIdx = 0;
   trueBestRecoTrackScore = 0;
   trueBestRecoHasTrack = 0;
   trueBestRecoHasShower = 0;
   nuScores = 0;
   NeutrinoNuScores = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("trueOrigin", &trueOrigin, &b_trueOrigin);
   fChain->SetBranchAddress("trueW", &trueW, &b_trueW);
   fChain->SetBranchAddress("trueX", &trueX, &b_trueX);
   fChain->SetBranchAddress("trueY", &trueY, &b_trueY);
   fChain->SetBranchAddress("trueQSqr", &trueQSqr, &b_trueQSqr);
   fChain->SetBranchAddress("truePt", &truePt, &b_truePt);
   fChain->SetBranchAddress("trueTheta", &trueTheta, &b_trueTheta);
   fChain->SetBranchAddress("trueNuPDG", &trueNuPDG, &b_trueNuPDG);
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
   fChain->SetBranchAddress("sliceKey", &sliceKey, &b_sliceKey);
   fChain->SetBranchAddress("sliceID", &sliceID, &b_sliceID);
   fChain->SetBranchAddress("sliceNuScore", &sliceNuScore, &b_sliceNuScore);
   fChain->SetBranchAddress("sliceTotalHits", &sliceTotalHits, &b_sliceTotalHits);
   fChain->SetBranchAddress("sliceTrueNuHits", &sliceTrueNuHits, &b_sliceTrueNuHits);
   fChain->SetBranchAddress("sliceTrueOrigin", &sliceTrueOrigin, &b_sliceTrueOrigin);
   fChain->SetBranchAddress("eventTotalTrueNuHits", &eventTotalTrueNuHits, &b_eventTotalTrueNuHits);
   fChain->SetBranchAddress("sliceVtxX", &sliceVtxX, &b_sliceVtxX);
   fChain->SetBranchAddress("sliceVtxY", &sliceVtxY, &b_sliceVtxY);
   fChain->SetBranchAddress("sliceVtxZ", &sliceVtxZ, &b_sliceVtxZ);
   fChain->SetBranchAddress("pfpKey", &pfpKey, &b_pfpKey);
   fChain->SetBranchAddress("pfpSelfID", &pfpSelfID, &b_pfpSelfID);
   fChain->SetBranchAddress("pfpParentID", &pfpParentID, &b_pfpParentID);
   fChain->SetBranchAddress("pfpRecoPDG", &pfpRecoPDG, &b_pfpRecoPDG);
   fChain->SetBranchAddress("pfpSliceKey", &pfpSliceKey, &b_pfpSliceKey);
   fChain->SetBranchAddress("pfpIsNuSlice", &pfpIsNuSlice, &b_pfpIsNuSlice);
   fChain->SetBranchAddress("pfpIsPrimary", &pfpIsPrimary, &b_pfpIsPrimary);
   fChain->SetBranchAddress("pfpNPrimaryChildren", &pfpNPrimaryChildren, &b_pfpNPrimaryChildren);
   fChain->SetBranchAddress("pfpTrackScore", &pfpTrackScore, &b_pfpTrackScore);
   fChain->SetBranchAddress("pfpHasTrackScore", &pfpHasTrackScore, &b_pfpHasTrackScore);
   fChain->SetBranchAddress("pfpNTracks", &pfpNTracks, &b_pfpNTracks);
   fChain->SetBranchAddress("pfpNShowers", &pfpNShowers, &b_pfpNShowers);
   fChain->SetBranchAddress("pfpHasTrack", &pfpHasTrack, &b_pfpHasTrack);
   fChain->SetBranchAddress("pfpHasShower", &pfpHasShower, &b_pfpHasShower);
   fChain->SetBranchAddress("pfpHasUniqueTrack", &pfpHasUniqueTrack, &b_pfpHasUniqueTrack);
   fChain->SetBranchAddress("pfpHasUniqueShower", &pfpHasUniqueShower, &b_pfpHasUniqueShower);
   fChain->SetBranchAddress("pfpNVertices", &pfpNVertices, &b_pfpNVertices);
   fChain->SetBranchAddress("pfpHasVertex", &pfpHasVertex, &b_pfpHasVertex);
   fChain->SetBranchAddress("pfpHasUniqueVertex", &pfpHasUniqueVertex, &b_pfpHasUniqueVertex);
   fChain->SetBranchAddress("pfpTrackID", &pfpTrackID, &b_pfpTrackID);
   fChain->SetBranchAddress("pfpTrackLength", &pfpTrackLength, &b_pfpTrackLength);
   fChain->SetBranchAddress("pfpTrackStartX", &pfpTrackStartX, &b_pfpTrackStartX);
   fChain->SetBranchAddress("pfpTrackStartY", &pfpTrackStartY, &b_pfpTrackStartY);
   fChain->SetBranchAddress("pfpTrackStartZ", &pfpTrackStartZ, &b_pfpTrackStartZ);
   fChain->SetBranchAddress("pfpTrackEndX", &pfpTrackEndX, &b_pfpTrackEndX);
   fChain->SetBranchAddress("pfpTrackEndY", &pfpTrackEndY, &b_pfpTrackEndY);
   fChain->SetBranchAddress("pfpTrackEndZ", &pfpTrackEndZ, &b_pfpTrackEndZ);
   fChain->SetBranchAddress("pfpTrackStartDirX", &pfpTrackStartDirX, &b_pfpTrackStartDirX);
   fChain->SetBranchAddress("pfpTrackStartDirY", &pfpTrackStartDirY, &b_pfpTrackStartDirY);
   fChain->SetBranchAddress("pfpTrackStartDirZ", &pfpTrackStartDirZ, &b_pfpTrackStartDirZ);
   fChain->SetBranchAddress("pfpTrackEndDirX", &pfpTrackEndDirX, &b_pfpTrackEndDirX);
   fChain->SetBranchAddress("pfpTrackEndDirY", &pfpTrackEndDirY, &b_pfpTrackEndDirY);
   fChain->SetBranchAddress("pfpTrackEndDirZ", &pfpTrackEndDirZ, &b_pfpTrackEndDirZ);
   fChain->SetBranchAddress("pfpTrackVertexDirX", &pfpTrackVertexDirX, &b_pfpTrackVertexDirX);
   fChain->SetBranchAddress("pfpTrackVertexDirY", &pfpTrackVertexDirY, &b_pfpTrackVertexDirY);
   fChain->SetBranchAddress("pfpTrackVertexDirZ", &pfpTrackVertexDirZ, &b_pfpTrackVertexDirZ);
   fChain->SetBranchAddress("pfpTrackTheta", &pfpTrackTheta, &b_pfpTrackTheta);
   fChain->SetBranchAddress("pfpTrackPhi", &pfpTrackPhi, &b_pfpTrackPhi);
   fChain->SetBranchAddress("pfpShowerID", &pfpShowerID, &b_pfpShowerID);
   fChain->SetBranchAddress("pfpShowerLength", &pfpShowerLength, &b_pfpShowerLength);
   fChain->SetBranchAddress("pfpShowerStartX", &pfpShowerStartX, &b_pfpShowerStartX);
   fChain->SetBranchAddress("pfpShowerStartY", &pfpShowerStartY, &b_pfpShowerStartY);
   fChain->SetBranchAddress("pfpShowerStartZ", &pfpShowerStartZ, &b_pfpShowerStartZ);
   fChain->SetBranchAddress("pfpShowerDirX", &pfpShowerDirX, &b_pfpShowerDirX);
   fChain->SetBranchAddress("pfpShowerDirY", &pfpShowerDirY, &b_pfpShowerDirY);
   fChain->SetBranchAddress("pfpShowerDirZ", &pfpShowerDirZ, &b_pfpShowerDirZ);
   fChain->SetBranchAddress("pfpVertexX", &pfpVertexX, &b_pfpVertexX);
   fChain->SetBranchAddress("pfpVertexY", &pfpVertexY, &b_pfpVertexY);
   fChain->SetBranchAddress("pfpVertexZ", &pfpVertexZ, &b_pfpVertexZ);
   fChain->SetBranchAddress("pfpTrueTrackID", &pfpTrueTrackID, &b_pfpTrueTrackID);
   fChain->SetBranchAddress("pfpTruePDG", &pfpTruePDG, &b_pfpTruePDG);
   fChain->SetBranchAddress("pfpNHits", &pfpNHits, &b_pfpNHits);
   fChain->SetBranchAddress("pfpNMatchedHits", &pfpNMatchedHits, &b_pfpNMatchedHits);
   fChain->SetBranchAddress("pfpTruthPurity", &pfpTruthPurity, &b_pfpTruthPurity);
   fChain->SetBranchAddress("trueIsReconstructed", &trueIsReconstructed, &b_trueIsReconstructed);
   fChain->SetBranchAddress("trueIsReconstructedInNuSlice", &trueIsReconstructedInNuSlice, &b_trueIsReconstructedInNuSlice);
   fChain->SetBranchAddress("trueNMatchedPfps", &trueNMatchedPfps, &b_trueNMatchedPfps);
   fChain->SetBranchAddress("trueBestRecoPfpIdx", &trueBestRecoPfpIdx, &b_trueBestRecoPfpIdx);
   fChain->SetBranchAddress("trueBestRecoTrackScore", &trueBestRecoTrackScore, &b_trueBestRecoTrackScore);
   fChain->SetBranchAddress("trueBestRecoHasTrack", &trueBestRecoHasTrack, &b_trueBestRecoHasTrack);
   fChain->SetBranchAddress("trueBestRecoHasShower", &trueBestRecoHasShower, &b_trueBestRecoHasShower);
   fChain->SetBranchAddress("nuScores", &nuScores, &b_nuScores);
   fChain->SetBranchAddress("NeutrinoNuScores", &NeutrinoNuScores, &b_NeutrinoNuScores);
   Notify();
}

bool signalDef::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be for a new TTree in a TChain. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void signalDef::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t signalDef::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef signalDef_cxx
