#define signalDef_cxx
#include "signalDef.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void signalDef::Loop()
{
//   In a ROOT session, you can do:
//      root> .L signalDef.C
//      root> signalDef t
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
// fChain is the TTree 
//
//
//
	
   
   int nSigma = 0;
   int nLambda = 0;
   int nGoodLambda = 0;
   int nGoodSigma = 0;

   if (fChain == 0) return;

   TString Samples[3];
   Samples[0] = "Sigma";
   Samples[1] = "Hyperons";
   Samples[2] = "Background";
   
   int nEvents[3] = {0};
   nEvents[0] = 1; // Sigma
   nEvents[1] = 1; // Hyperons
   nEvents[2] = 1; // BKG
   

   double PoT[3] = {0};
   PoT[0] = 1; // sigma
   PoT[1] = 1; // hyperons
   PoT[2] = 1; // bkg

   double scale[3] = {0};
   scale[0] = PoT[2] / PoT[0];
   scale[1] = PoT[2] / PoT[1];
   scale[2] = PoT[2] / PoT[2];

   const int nCuts = 4;
   TString CutName[nCuts];
   CutName[0] = "AllEvents";
   CutName[1] = "FV";
   CutName[2] = "nShowers";
   CutName[3] = "nTracks";

   const int nVars = 4;
   TString VarName[nVars];
   VarName[0] = "nTracks";
   VarName[1] = "nShowers";

   TH1F *Hist[nVars][3][nCuts];

   for (int s = 0; s < 4; s++){ // samples
	  for (int v = 0; v < nVars; v++){ // variables
		 for (int c = 0; c < nCuts; c++){ // cuts
			Hist[v][s][c] = new TH1F(VarName[v] + Samples[s] + CutName[c], "", 100, -1, -1);
		 }
	  }
   }
   

   TFile *outFile = TFile::Open("TreeS.root", "RECREATE");

   if (!outFile || outFile->IsZombie()){
	   std::cerr<<"Could not open file!!"<<std::endl;
	   return;
   }

   outFile->cd();
   TTree *signalTree = fChain->CloneTree(0);
   signalTree->SetName("TreeS");
   signalTree->SetDirectory(outFile);

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
      bool IsGoodLambda = false;
      bool IsSigma0 = false;
      bool IsGoodSigma = false;
      bool IsPhoton = false;
      bool IsKaonp = false;
      bool IsKaonm = false;
      bool IsKaon0 = false;
      bool IsProton = false;
      bool IsPionm = false;
   
      bool IsInFV = false; // true vertex FV for signal def
      bool IsInRecoFV = false;
      bool IsInTrackRange = false;
      bool IsInShowerRange = false;
      bool IsSignal = false;
      bool Cuts[nCuts] = {false};
      int s = 0; // sample index (0 = sigma signal, 1 = other hyperons, 2=bkg)

      int sumCounter = 0;
      int daughterCounter = 0;

      std::cout<<"****************** new event ********************"<<std::endl;
      std::cout<<"truePDG size = "<<truePDG->size()<<", daughterSize size = "<<daughterSize->size()<<std::endl;
      std::cout<<"trackScores size = "<<trackScores->size()<<" track start pos size = "<<TrackStartPositionX->size()<<std::endl;
      std::cout<<"distance to reco vertex size = "<<DistanceToRecoVertex->size()<<std::endl;
      std::cout<<"track lengths size = "<<TrackLengths->size()<<std::endl;
      std::cout<<"trackCount = "<<trackCount<<" showerCount = "<<showerCount<<std::endl;

      for (int i_vert = 0; i_vert < vertexSize->size(); i_vert++){ // Loop over truth vertices in event (may be more than 1)
      sumCounter += vertexSize->at(i_vert); // Keep track of what particle index we are at over all vertices

      for (int i = (sumCounter - vertexSize->at(i_vert)); i < sumCounter; i++){ // Loop through particles from previous vertex end to current vertex end   
    
	      std::cout<<"------------------- new particle ------------------"<<std::endl;    
          std::cout<<"Particle index "<<i<<" in vertex index  "<<i_vert<<" has vertex "<<vertexX->at(i)<<", "<<vertexY->at(i)<<", "
		  <<vertexZ->at(i)<<", has PDG "<<truePDG->at(i)<<
		  ", mother "<<motherPDG->at(i)<<", num of daughters: "<<daughterSize->at(i)<<" trueP "<<trueP->at(i)<<std::endl;

	  for (int i_daug = daughterCounter; i_daug < (daughterCounter + daughterSize->at(i)); i_daug++){
		  std::cout<<"Daughter index "<<i_daug<<" has PDG "<<daughterPDG->at(i_daug)<<std::endl;
	  } 

	  if (std::abs(vertexX->at(i)) < 180 && std::abs(vertexY->at(i)) < 180 && vertexZ->at(i) > 10 && vertexZ->at(i) < 450){
		  IsInFV = true;
	  }

	  if (truePDG->at(i) == 310 || truePDG->at(i) == 130 || truePDG->at(i) == 311){
		  IsKaon0 = true;
	  }
	  if (truePDG->at(i) == 321){
		  IsKaonp = true;
	  }
	  if (truePDG->at(i) == -321){
		  IsKaonm = true;
	  }
	  if (truePDG->at(i) == 2212){
		  IsProton = true;
	  }
	  if (truePDG->at(i) == 3122){
		  IsLambda = true;
		  bool HasProtonDaught = false;
		  bool HasPionmDaught = false;

		  for (int i_daught = daughterCounter; i_daught < (daughterCounter + daughterSize->at(i)); i_daught++){
			  std::cout<<"Lambda daughter num "<<i_daught<<" has PDG "<<daughterPDG->at(i_daught)<<std::endl;
			  if (daughterPDG->at(i_daught) == 2212){
				  HasProtonDaught = true;
			  }
			  if (daughterPDG->at(i_daught) == -211){
				  HasPionmDaught = true;
			  }
		  }

		  if (HasProtonDaught && HasPionmDaught){
			  IsGoodLambda = true;
		  }
	  }

	  if (truePDG->at(i) == 3212){
		IsSigma0 = true;
		bool HasPhotonDaughter = false;
		bool HasLambdaDaughter = false;

		for (int i_daught = daughterCounter; i_daught < (daughterCounter + daughterSize->at(i)); i_daught++){
			std::cout<<"Sigma0 daughter has PDG "<<daughterPDG->at(i_daught)<<std::endl;

			if(daughterPDG->at(i_daught) == 22){
				HasPhotonDaughter = true;
			}
			if(daughterPDG->at(i_daught) == 3122){
				HasLambdaDaughter = true;
			}
		}

		if (HasPhotonDaughter && HasLambdaDaughter){
			IsGoodSigma = true;
		}
	   //std::cout<<"SIGMA FOUND"<<std::endl;
	   //std::cout<<"sigma has "<<daughterPDG->size()<<" daughters"<<std::endl;
	  /* for (int i_daugh = 0; i_daugh < daughterPDG->size(); i++){
		std::cout<<"sigma daugher particle has PDG "<<daughterPDG->at(i_daugh)<<std::endl;
		nSigma++;
	   }*/

	  }
	  if (truePDG->at(i) == -13){
		  IsAntiMuon = true;
	  }
	  if (truePDG->at(i) == -211){
		  IsPionm = true;
	  }
	  if (truePDG->at(i) == 22){
		  IsPhoton = true;
	  }

	  daughterCounter += daughterSize->at(i); // keep track of cumulative index of daughter PDGs
      } // end of loop over particles in vertex
  
      } // end of loop over vertices in event


      // SIGNAL DEFINITION
      
      if (IsInFV && IsAntiMuon && !IsKaonp && !IsKaonm && !IsKaon0){
	      if (IsGoodSigma){
			nGoodSigma++;
			std::cout<<"signal event, nSignal = "<<nGoodSigma<<std::endl;
			IsSignal = true;
			signalTree->Fill();
	    }
	    if (IsLambda){
		nLambda++;
		std::cout<<"lambda event, nLambda = "<<nLambda<<std::endl; // check daughters of lambda for pion- and proton
		if (IsGoodLambda){
			nGoodLambda++;
			std::cout<<"lambda w proton/pion daughters event, n goodlambda = "<<nGoodLambda<<std::endl;
		}
	    }
	}

      // Preselection

      if (std::abs(RecoVertexX) < 180 && std::abs(RecoVertexY) < 180 && RecoVertexZ > 10 && RecoVertexZ < 450){IsInRecoFV = true;}

      if (trackCount > 2 && trackCount < 6){IsInTrackRange = true;}
      if (showerCount == 1 || showerCount == 2){IsInShowerRange = true;}

      Cuts[0] = true;
      Cuts[1] = Cuts[0] && IsInRecoFV;
      Cuts[2] = Cuts[1] && IsInTrackRange;
      Cuts[3] = Cuts[2] && IsInShowerRange;

      for (int c = 0; c<nCuts; c++){

	if(Cuts[c]){

		std::cout<<"fill hists here"<<std::endl;
	}

      }

     
   } // end event loop

   outFile->cd();
   signalTree->Write("TreeS"); // Write signal tree and close file
   outFile->Close();
   delete outFile;

   std::cout<<"num of sigma: "<<nSigma<<std::endl;
   std::cout<<"num of lambda: "<<nLambda<<std::endl;
   std::cout<<"num of good sigma: "<<nGoodSigma<<std::endl;
   std::cout<<"num of good lambda: "<<nGoodLambda<<std::endl;
}
