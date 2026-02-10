#define signalPreselect_cxx
#include "signalPreselect.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void signalPreselect::Loop()
{
//   In a ROOT session, you can do:
//      root> .L signalPreselect.C
//      root> signalPreselect t
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

   const int Ncuts = 3;
   TString CutName[Ncuts];
   CutName[0] = "AllEvents";
   CutName[1] = "NShowers";
   CutName[2] = "NTracks";


   const int Nvars = 3;
   TString VarName[nVars];
   VarName[0] = "trackCount";
   VarName[1] = "showerCount";
   VarName[2] = "RecoVertexX";
   VarName[3] = "RecoVertexY";
   VarName[4] = "RecoVertexZ";


   for (int s = 0; s < 3; s++){ // samples
	   for (int v = 0; v < Nvars; v++){ // variables
		   for (int c = 0; c < Ncuts; c++){ // cuts
			   Hist[s][v][c] = new TH1F("", "", 100, -1, -1)

   TFile *outFile = TFile::Open("preTreeS.root", "RECREATE");
   if (!outFile || outFile->IsZombie()){
	   std::cerr<<"Could not open file!"<<std::endl;
	   return;
   }

   outFile->cd();
   TTree *signalPreTree = fChain->CloneTree(0);
   signalPreTree->SetName("preTreeS");
   signalPreTree->SetDirectory(outFile);
  
   int minShowers = 1;
   int maxShowers = 2;

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

      // Apply preselection cuts
      
      // Only include signal events that are in reco FV
      if (std::abs(RecoVertexX) < 180 && std::abs(RecoVertexY) < 180 && RecoVertexZ > 10 && RecoVertexZ < 450){IsInRecoFV = true;}
      
      // At least 3 tracks, less than 6
      if (trackCount > 2 && trackCount < 6){HasTrackRange = true;}	

      // 1 or 2 showers
      if (showerCount == minShowers || showerCount == maxShowers){HasShowerRange = true;};

      // Fill new tree
      if (IsInRecoFV){
	   if(HasTrackRange){
//		  if(HasShowerRange){
			  signalPreTree->Fill();
		  }
//	   }
      }
   }

   outFile->cd();
   signalPreTree->Write("preTreeS");
   outFile->Close();
   delete outFile;
}
