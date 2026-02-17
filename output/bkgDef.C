#define bkgDef_cxx
#include "bkgDef.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void bkgDef::Loop()
{
//   In a ROOT session, you can do:
//      root> .L bkgDef.C
//      root> bkgDef t
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
   
   int nSigma = 0;
   int nLambda = 0;
   int nGoodLambda = 0;
   int nGoodSigma = 0;
   int bkgCount = 0;

   if (fChain == 0) return;

   TFile *outFile = TFile::Open("TreeB.root", "RECREATE");

   if (!outFile || outFile->IsZombie()){
	   std::cerr<<"Could not open file!!"<<std::endl;
	   return;
   }

   outFile->cd();
   TTree *bkgTree = fChain->CloneTree(0);
   bkgTree->SetName("TreeB");
   bkgTree->SetDirectory(outFile);

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
      bool IsSignal = false;
      int signalCount = 0;

      int sumCounter = 0;
      int daughterCounter = 0;

      std::cout<<"****************** new event ********************"<<std::endl;
      std::cout<<"truePDG size = "<<truePDG->size()<<", daughterSize size = "<<daughterSize->size()<<std::endl;

      for (int i_vert = 0; i_vert < vertexSize->size(); i_vert++){ // Loop over truth vertices in event (may be more than 1)
      sumCounter += vertexSize->at(i_vert); // Keep track of what particle index we are at over all vertices

      for (int i = (sumCounter - vertexSize->at(i_vert)); i < sumCounter; i++){ // Loop through particles from previous vertex end to current vertex end   
    
	      std::cout<<"------------------- new particle ------------------"<<std::endl;    
          std::cout<<"Particle index "<<i<<" in vertex index  "<<i_vert<<" has vertex "<<vertexX->at(i)<<", "<<vertexY->at(i)<<", "
		  <<vertexZ->at(i)<<", has PDG "<<truePDG->at(i)<<
		  ", mother "<<motherPDG->at(i)<<", num of daughters: "<<daughterSize->at(i)<<" trueP "<<trueP->at(i)<<std::endl;

	 // for (int i_daug = daughterCounter; i_daug < (daughterCounter + daughterSize->at(i)); i_daug++){
	//	  std::cout<<"Daughter index "<<i_daug<<" has PDG "<<daughterPDG->at(i_daug)<<std::endl;
	  //} 

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
      
      if (IsAntiMuon && IsInFV){
	if (!IsKaonp && !IsKaonm && !IsKaon0){
	    if(IsSigma0){
		nSigma++;
		if (IsGoodSigma){
			nGoodSigma++;
			std::cout<<"signal event, nSignal = "<<nGoodSigma<<std::endl;
			signalCount++;
			IsSignal = true;
		}
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
      }

      if (!IsSignal){
	      bkgTree->Fill();
	      bkgCount++;
      }
     
   }

   outFile->cd();
   bkgTree->Write("TreeB"); // Write signal tree and close file
   outFile->Close();
   delete outFile;

   std::cout<<"num of sigma: "<<nSigma<<std::endl;
   std::cout<<"num of lambda: "<<nLambda<<std::endl;
   std::cout<<"num of good sigma: "<<nGoodSigma<<std::endl;
   std::cout<<"num of good lambda: "<<nGoodLambda<<std::endl;
   std::cout<<"num of bkg events: "<<bkgCount<<std::endl;
}
