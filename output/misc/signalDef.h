//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 13 15:25:29 2026 by ROOT version 6.38.00
// from TTree tree/Output TTree
// found on file: hyperonAnalysisSIG.root
//////////////////////////////////////////////////////////

#ifndef signalDef_h
#define signalDef_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class signalDef {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          eventID;
   vector<float>   *trueP;
   vector<int>     *truePDG;
   vector<int>     *daughterPDG;
   vector<int>     *motherPDG;
   vector<float>   *vertexX;
   vector<float>   *vertexY;
   vector<float>   *vertexZ;
   vector<int>     *vertexSize;
   vector<int>     *daughterSize;
   Float_t         isoVertexX;
   Float_t         isoVertexY;
   Float_t         isoVertexZ;
   Int_t           nPFParticles;
   Int_t           nPrimaryChildren;
   Int_t           trackCount;
   Int_t           showerCount;
   vector<int>     *TrackIDs;
   vector<float>   *TrackLengths;
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
   vector<float>   *TrackStartPositionX;
   vector<float>   *TrackStartPositionY;
   vector<float>   *TrackStartPositionZ;
   vector<float>   *TrackEndPositionX;
   vector<float>   *TrackEndPositionY;
   vector<float>   *TrackEndPositionZ;
   vector<int>     *pfpTrackPDG;
   vector<int>     *pfpShowerPDG;

   // List of branches
   TBranch        *b_eventID;   //!
   TBranch        *b_trueP;   //!
   TBranch        *b_truePDG;   //!
   TBranch        *b_daughterPDG;   //!
   TBranch        *b_motherPDG;   //!
   TBranch        *b_vertexX;   //!
   TBranch        *b_vertexY;   //!
   TBranch        *b_vertexZ;   //!
   TBranch        *b_vertexSize;   //!
   TBranch        *b_daughterSize;   //!
   TBranch        *b_isoVertexX;   //!
   TBranch        *b_isoVertexY;   //!
   TBranch        *b_isoVertexZ;   //!
   TBranch        *b_nPFParticles;   //!
   TBranch        *b_nPrimaryChildren;   //!
   TBranch        *b_trackCount;   //!
   TBranch        *b_showerCount;   //!
   TBranch        *b_TrackIDs;   //!
   TBranch        *b_TrackLengths;   //!
   TBranch        *b_RecoVertexX;   //!
   TBranch        *b_RecoVertexY;   //!
   TBranch        *b_RecoVertexZ;   //!
   TBranch        *b_DistanceToRecoVertex;   //!
   TBranch        *b_nuScores;   //!
   TBranch        *b_trackScores;   //!
   TBranch        *b_muonTrackScores;   //!
   TBranch        *b_protonTrackScores;   //!
   TBranch        *b_pionTrackScores;   //!
   TBranch        *b_NeutrinoNuScores;   //!
   TBranch        *b_TrackStartPositionX;   //!
   TBranch        *b_TrackStartPositionY;   //!
   TBranch        *b_TrackStartPositionZ;   //!
   TBranch        *b_TrackEndPositionX;   //!
   TBranch        *b_TrackEndPositionY;   //!
   TBranch        *b_TrackEndPositionZ;   //!
   TBranch        *b_pfpTrackPDG;   //!
   TBranch        *b_pfpShowerPDG;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("hyperonAnalysisSIG.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("hyperonAnalysisSIG.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("hyperonAnalysisSIG.root:/ana");
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
   trueP = 0;
   truePDG = 0;
   daughterPDG = 0;
   motherPDG = 0;
   vertexX = 0;
   vertexY = 0;
   vertexZ = 0;
   vertexSize = 0;
   daughterSize = 0;
   TrackIDs = 0;
   TrackLengths = 0;
   DistanceToRecoVertex = 0;
   nuScores = 0;
   trackScores = 0;
   muonTrackScores = 0;
   protonTrackScores = 0;
   pionTrackScores = 0;
   NeutrinoNuScores = 0;
   TrackStartPositionX = 0;
   TrackStartPositionY = 0;
   TrackStartPositionZ = 0;
   TrackEndPositionX = 0;
   TrackEndPositionY = 0;
   TrackEndPositionZ = 0;
   pfpTrackPDG = 0;
   pfpShowerPDG = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("trueP", &trueP, &b_trueP);
   fChain->SetBranchAddress("truePDG", &truePDG, &b_truePDG);
   fChain->SetBranchAddress("daughterPDG", &daughterPDG, &b_daughterPDG);
   fChain->SetBranchAddress("motherPDG", &motherPDG, &b_motherPDG);
   fChain->SetBranchAddress("vertexX", &vertexX, &b_vertexX);
   fChain->SetBranchAddress("vertexY", &vertexY, &b_vertexY);
   fChain->SetBranchAddress("vertexZ", &vertexZ, &b_vertexZ);
   fChain->SetBranchAddress("vertexSize", &vertexSize, &b_vertexSize);
   fChain->SetBranchAddress("daughterSize", &daughterSize, &b_daughterSize);
   fChain->SetBranchAddress("isoVertexX", &isoVertexX, &b_isoVertexX);
   fChain->SetBranchAddress("isoVertexY", &isoVertexY, &b_isoVertexY);
   fChain->SetBranchAddress("isoVertexZ", &isoVertexZ, &b_isoVertexZ);
   fChain->SetBranchAddress("nPFParticles", &nPFParticles, &b_nPFParticles);
   fChain->SetBranchAddress("nPrimaryChildren", &nPrimaryChildren, &b_nPrimaryChildren);
   fChain->SetBranchAddress("trackCount", &trackCount, &b_trackCount);
   fChain->SetBranchAddress("showerCount", &showerCount, &b_showerCount);
   fChain->SetBranchAddress("TrackIDs", &TrackIDs, &b_TrackIDs);
   fChain->SetBranchAddress("TrackLengths", &TrackLengths, &b_TrackLengths);
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
   fChain->SetBranchAddress("TrackStartPositionX", &TrackStartPositionX, &b_TrackStartPositionX);
   fChain->SetBranchAddress("TrackStartPositionY", &TrackStartPositionY, &b_TrackStartPositionY);
   fChain->SetBranchAddress("TrackStartPositionZ", &TrackStartPositionZ, &b_TrackStartPositionZ);
   fChain->SetBranchAddress("TrackEndPositionX", &TrackEndPositionX, &b_TrackEndPositionX);
   fChain->SetBranchAddress("TrackEndPositionY", &TrackEndPositionY, &b_TrackEndPositionY);
   fChain->SetBranchAddress("TrackEndPositionZ", &TrackEndPositionZ, &b_TrackEndPositionZ);
   fChain->SetBranchAddress("pfpTrackPDG", &pfpTrackPDG, &b_pfpTrackPDG);
   fChain->SetBranchAddress("pfpShowerPDG", &pfpShowerPDG, &b_pfpShowerPDG);
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
