#define testingSIG_cxx
#include "testingSIG.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void testingSIG::Loop()
{
//   In a ROOT session, you can do:
//      root> .L testingSIG.C
//      root> testingSIG t
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
//
// DEFINE 
//
//
   
   int nSigma = 0;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      // DEFINE PER EVENT VARIABLES
      
      bool IsAntiMuon = false;
      bool IsLambda = false;
      bool IsSigma0 = false;
      bool IsPhoton = false;
      bool IsKaonp = false;
      bool IsKaonm = false;
      bool IsKaon0 = false;
      bool IsProton = false;
      bool IsPionm = false;
   
      bool IsInFV = false; // true vertex FV for signal def
      bool IsSignal = false;

      int sumCounter = 0;
      for (int i_vert = 0; i_vert < vertexSize->size(); i_vert++){ // Loop over truth vertices in event (may be more than 1)
      sumCounter += vertexSize->at(i_vert); // Keep track of what particle index we are at over all vertices

      for (int i = (sumCounter - vertexSize->at(i_vert)); i < sumCounter; i++){ // Loop through particles from previous vertex end to current vertex end  

          std::cout<<"Particle index "<<i<<" in vertex index  "<<i_vert<<" has vertex "<<vertexX->at(i)<<", "<<vertexY->at(i)<<", "
		  <<vertexZ->at(i)<<", has PDG "<<truePDG->at(i)<<
		  ", mother "<<motherPDG->at(i)<<" trueP "<<trueP->at(i)<<std::endl;

	  if (truePDG->at(i) == 3212){

	   std::cout<<"SIGMA FOUND"<<std::endl;
	   std::cout<<"sigma has "<<daughterPDG->size()<<" daughters"<<std::endl;
	  /* for (int i_daugh = 0; i_daugh < daughterPDG->size(); i++){
		std::cout<<"sigma daugher particle has PDG "<<daughterPDG->at(i_daugh)<<std::endl;
		nSigma++;
	   }*/

	  }

      }
      }


      // SIGNAL DEFINITION
      
      if (IsAntiMuon && IsInFV){
	if (!IsKaonp && !IsKaonm && !IsKaon0){
	    if(IsSigma0 && IsPhoton){
		std::cout<<"check if sigma has lambda daughter... or otherway around"<<std::endl;
	    }
	}
      }

     
   }

   std::cout<<"num of sigma: "<<nSigma<<std::endl;
}
