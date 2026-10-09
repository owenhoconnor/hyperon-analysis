#define preselection_cxx
#include "preselection.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void preselection::Loop()
{
//   In a ROOT session, you can do:
//      root> .L preselection.C
//      root> preselection t
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

float efficiency = 0.0;
float purity = 0.0;

// preselection cut values
float recoVtxXHighCut = 180;
float recoVtxXLowCut = 180;
float recoVtxYHighCut = 180;
float recoVtxYLowCut = -180;
float recoVtxZHighCut = 450;
float recoVtxZLowCut = 10;
float nuScoreCut = 0.5;
float bcfmCut = 0.2;
float nTracksCut = 3;
float nShowersCut = 1;


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      std::cout<<"---------------------------------"
      std::cout<<"Processing event "<<jentry<<" out of "<<nentries<<std::endl;

      nuSliceNuScore = sliceNuScore->at(nuSliceIdx);
      nuSliceBCFMScore = sliceBaryFlashScore->at(nuSliceIdx);

      // Reco FV cut
      if (std::abs(RecoVertexX) < 180 && std::abs(RecoVertexY) < 180 && RecoVertexZ < 450 && RecoVertexZ > 10){
         isInRecoFV = true;
      }

      // nuScore cut
      if(nuSliceNuScore < nuScoreCut){
         continue;
      }

      // BCFM cut
      if(nuSliceBCFMScore < bfcmCut){
         continue;
      }

      // Track/Shower multiplicty cut (3 Track + 1 Shower)

      if(nTracks != nTracksCut){
         continue;
      }
      if(nShowers != nShowersCut){
         continue;
      }
   }

   // Efficiency * Purity optimization

   float efficiency; // selected signal / total true signal
   float purity; // selected signal / total selected

   float step = 1;
   for (int i = 0; i < ; ++i){
      float effPur = efficiency * purity;
   }

   std::cout<<"=============== Preselection Results ==================="<<std::endl;
   std::cout<<"========================================================"std::endl;
}
