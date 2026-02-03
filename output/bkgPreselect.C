#define bkgPreselect_cxx
#include "bkgPreselect.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void bkgPreselect::Loop()
{
//   In a ROOT session, you can do:
//      root> .L bkgPreselect.C
//      root> bkgPreselect t
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

   TFile *outFile = TFile::Open("preTreeB.root", "RECREATE");
   if (!outFile || outFile->IsZombie()){
	   std::cerr<<"Could not open file!"<<std::endl;
   }

   outFile->cd();
   TTree *bkgPreTree = fChain->CloneTree(0);
   bkgPreTree->SetName("preTreeB");
   bkgPreTree->SetDirectory(outFile);

   int minShowers = 1;
   int maxShowers = 2;
   int nBkg = 0;
   int nPreBkg = 0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      bool IsInRecoFV = false;
      bool HasTrackRange = false;
      bool HasShowerRange = false;
   
      std::cout<<"Event in TreeB with reco vertex: "<<RecoVertexX<<", "<<RecoVertexY<<", "<<RecoVertexZ<<"and trackCount = "<<trackCount
<<"and showerCount = "<<showerCount<<std::endl;	
	nBkg++;      

      // Apply preselection

      // Only include events inside FV
      if (std::abs(RecoVertexX) < 180 && std::abs(RecoVertexY) < 180 && RecoVertexZ > 10 && RecoVertexZ < 450){IsInRecoFV = true;}

      // More than 3 but less than 6 tracks
      if (trackCount > 2 && trackCount < 6){HasTrackRange = true;}

      // 1 or 2 showers
      if (showerCount == minShowers || showerCount == maxShowers){HasShowerRange = true;}

      // Fill preselection tree
      if (IsInRecoFV && HasTrackRange && HasShowerRange){
	      nPreBkg++;
	      bkgPreTree->Fill();
      }
      
   }

   std::cout<<"num of total background = "<<nBkg<<", num of bkg after preselection = "<<nPreBkg<<std::endl;
   outFile->cd();
   bkgPreTree->Write("preTreeB");
   outFile->Close();
   delete outFile;
}
