#define hyperonAnalysis_cxx
#include "hyperonAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void hyperonAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L test.C
//      root> test t
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


   auto *hLambda = new TH1D("hLambda", "trueP Lambda Baryons", 100, -1, -1); // Declaration of histograms
   auto *hSigma = new TH1D("hSigma", "trueP Sigma Baryons", 100, -1, -1);
   auto *hProton = new TH1D("hProton", "trueP Protons", 100, -1, -1);
   auto *hMuon = new TH1D("hMuon", "trueP Muons", 100, -1, -1);
   auto *hPion = new TH1D("hPion", "trueP Pions", 100, -1, -1);
   auto *hCCQELambda = new TH1D("hCCQELambda", "trueP CCQE Lambdas", 100, -1, -1);
   auto *hCCQESigma0 = new TH1D("hCCQESigma0", "trueP CCQE Sigma0s", 100, -1, -1);
   auto *hPhoton = new TH1D("hPhoton", "trueP Photons", 100, -1, -1);
   auto *hTrackCount = new TH1D("hTrackCount", "Track Count", 100, -1, -1);
   auto *hShowerCount = new TH1D("hShowerCount", "Shower Count", 100, -1, -1);
   auto *hRecoVertexX = new TH1D("hRecoVertexX", "Reco Vertex X", 100, -1, -1);
   auto *hRecoVertexY = new TH1D("hRecoVertexY", "RecoVertex Y", 100, -1, -1);
   auto *hRecoVertexZ = new TH1D("hRecoVertexZ", "Reco Vertex Z", 100, -1, -1);
   auto *hTrackScores = new TH1D("hTrackScores", "Track Scores", 100, -1, -1);
   auto *hPFParticles = new TH1D("hPFParticles", "Num PFParticles", 100, -1, -1);

   int NCCQELambda = 0; // Declaration of counters for num of events
   int NCCQESigma0 = 0;

   int nfiles = 7; // num of hyperon files
   	
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
    /*  std::cout<<"truePDG size is "<<truePDG->size()<<std::endl;
      std::cout<<"trueP size is "<<trueP->size()<<std::endl;
      std::cout<<"daughterPDG size is "<<daughterPDG->size()<<std::endl;
      std::cout<<"motherPDG size is "<<motherPDG->size()<<std::endl;
     */// std::cout<<"vertices are "<<vertexX<<" "<<vertexY<<" "<<vertexZ<<std::endl;

      bool IsInFV = false;

      bool IsLambda = false;
      bool IsSigma0 = false;
      bool IsAntiMuon = false;
      bool IsPionm = false;
      bool IsProton = false;
      bool IsKaon0 = false;
      bool IsKaonp = false;
      int i_lambda = 0;
      int i_sigma0 = 0;
      int sumCounter = 0;
      int sumCounter2 = 0;
      std::cout<<"event start **************************"<<std::endl;
      std::cout<<"truePDG size is: "<<truePDG->size()<<std::endl;
      if (TrackIDs->size() > 0 && TrackLengths->size() > 0 && DistanceToRecoVertex->size() > 0){
      std::cout<<"TrackIDs size is: "<<TrackIDs->size()<<std::endl;
      std::cout<<"TrackLengths size is: "<<TrackLengths->size()<<std::endl;
      std::cout<<"Reco Vertex is "<<RecoVertexX<<" ,"<<RecoVertexY<<" ,"<<RecoVertexZ<<std::endl;
      std::cout<<"DistanceToRecoVertex size is "<<DistanceToRecoVertex->size()<<std::endl;
      //std::cout<<"TrackIDs first element is: "<<TrackIDs->at(0)<<std::endl;
      //std::cout<<"TrackLengths first element is is: "<<TrackLengths->at(0)<<std::endl;
      //std::cout<<"DistanceToRecoVertex first element is: "<<DistanceToRecoVertex->at(0)<<std::endl;
      }
  
      std::cout<<"num of primary particles = "<<nPFParticles<<std::endl;
      std::cout<<"num of tracks = "<<trackCount<<std::endl;
      std::cout<<"num of showers = "<<showerCount<<std::endl;
      std::cout<<"num of track score entries = "<<trackScores->size()<<std::endl;

      if (std::abs(RecoVertexX) < 180 && std::abs(RecoVertexY) < 180 && RecoVertexZ > 10 & RecoVertexZ < 450){
      	 hRecoVertexX->Fill(RecoVertexX);
     	 hRecoVertexY->Fill(RecoVertexY);
     	 hRecoVertexZ->Fill(RecoVertexZ);

	 hPFParticles->Fill(nPFParticles);
	 hTrackCount->Fill(trackCount);
	 hShowerCount->Fill(showerCount);

      	 for (int i = 0; i < trackScores->size(); i++){
             if(trackScores->at(i) > 0)
	          hTrackScores->Fill(trackScores->at(i));
      	}
      }
      //if (std::abs(isoVertexX) < 180 && std::abs(isoVertexY) < 180 && isoVertexZ > 10 & isoVertexZ < 450){IsInFV = true;}     

      for (int i_vert = 0; i_vert < vertexSize->size(); i_vert++){ // Loop over interaction vertices (may be more than 1)
//	std::cout<<"i_vert is: "<<i_vert<<std::endl;
//      std::cout<<"vertexSize at (i_vert) is: "<<vertexSize->at(i_vert)<<std::endl;
        sumCounter += vertexSize->at(i_vert);
//	std::cout<<"cumulative sum of vertexSize is: "<<sumCounter<<std::endl;
//	std::cout<<"loop from i = " <<(sumCounter - vertexSize->at(i_vert))<<"til i < "<<sumCounter<<std::endl;
 	      
      for (int i = (sumCounter - vertexSize->at(i_vert)); i < sumCounter; i++){ // Loop over particles 

	     // std::cout<<"particle loop start ------------"<<std::endl;
	      //std::cout<<"vertexX at i = "<<i<<" is "<<vertexX->at(i)<<std::endl;
	      //std::cout<<"vertexY at i = "<<i<<" is "<<vertexY->at(i)<<std::endl;
	      //std::cout<<"vertexZ at i = "<<i<<" is "<<vertexZ->at(i)<<std::endl;
     // std::cout<<"TrackID =  "<<TrackIDs->at(i)<<std::endl;
     // std::cout<<"track length = "<<TrackLengths->at(i)<<std::endl;
     // std::cout<<"distance to reco vertex = "<<DistanceToRecoVertex->at(i)<<std::endl;
         
      if (std::abs(vertexX->at(i)) < 180 && std::abs(vertexY->at(i)) < 180 && vertexZ->at(i) > 10 & vertexZ->at(i) < 450){IsInFV = true;}
 
      if (truePDG->at(i) == 2212){ // Flag if particle is a proton
	IsProton=true;
	//std::cout<<"Proton with P"<<"\t"<<trueP->at(i) <<std::endl;
	//hProton->Fill(trueP->at(i));	
      }

      if (truePDG->at(i) == 3122){ // Flag if particle is Lambda
	IsLambda=true;
	//std::cout<<"Lambda with P"<<"\t"<<trueP->at(i)<<std::endl;
	//hLambda->Fill(trueP->at(i));
	i_lambda = i;
	}

      if (truePDG->at(i) == 3212){ // Flag is particle is Sigma
	IsSigma0=true;
	//std::cout<<"Sigma0 with P"<<"\t"<<trueP->at(i)<<std::endl;
	//hSigma->Fill(trueP->at(i));
	i_sigma0 = i;
      }

      if (truePDG->at(i) == -13){ // Flag if particle is antimuon
	 IsAntiMuon=true;
	 //std::cout<<"antimuon with P"<<"\t"<<trueP->at(i)<<std::endl;
         //hMuon->Fill(trueP->at(i));
      }

      if (truePDG->at(i) == -211){
	 IsPionm=true;
	 //std::cout<<"pi- with P"<<"\t"<<trueP->at(i)<<std::endl;
         //hPion->Fill(trueP->at(i));
      }

      if (truePDG->at(i) == 311){
	      IsKaon0 = true;
      }

      if (truePDG->at(i) == 321){
	      IsKaonp = true;
      }
    }
  }

   if (IsAntiMuon && IsInFV && !IsKaon0 && !IsKaonp){ //  no kaons to exclude associated production
	std::cout<<"-------------------------------"<<std::endl;   
	if (IsLambda){NCCQELambda++; hCCQELambda->Fill(trueP->at(i_lambda)); std::cout<<"CCQE Lambda event"<<std::endl;}
	if (IsSigma0){NCCQESigma0++; hCCQESigma0->Fill(trueP->at(i_sigma0)); std::cout<<"CCQE Sigma0 event"<<std::endl;}
	// loop over vertices and particles again  here to see what i have in my interactions, then exclude what i dont need
	for (int i_vert = 1; i_vert < vertexSize->size(); ++i_vert){
		sumCounter2 += vertexSize->at(i_vert);
	for (int i = (sumCounter2 - vertexSize->at(i_vert)); i < sumCounter2; ++i){
		std::cout<<"particle at index: "<<i<<"has pdg code "<<truePDG->at(i)<<std::endl;
		//std::cout<<"particle is from vertex "<<i_vert<<std::endl;
		std::cout<<"particle has vertex: "<<vertexX->at(i)<<", "<<vertexY->at(i)<<", "<<vertexZ->at(i)<<std::endl;


		/*if (truePDG->at(i) == 22){
			std::cout<<"photon with P = "<<trueP->at(i)<<std::endl;
	 		hPhoton->Fill(trueP->at(i));
		}

		if (truePDG->at(i) == 2212){
			std::cout<<"proton with P = "<<trueP->at(i)<<std::endl;
			hProton->Fill(trueP->at(i));
		}

		if (truePDG->at(i) == 2112){
			std::cout<<"neutron with P = "<<trueP->at(i)<<std::endl;
		}
		
		if (truePDG->at(i) == -13){
			std::cout<<"Antimuon with P = "<<trueP->at(i)<<std::endl;
		}

		if (truePDG->at(i) == 3122){
			std::cout<<"Lambda with P = "<<trueP->at(i)<<std::endl;
		}

		if (truePDG->at(i) == 3212){
			std::cout<<"Sigma with P = "<<trueP->at(i)<<std::endl;
		}

		if (truePDG->at(i) == 14){
			std::cout<<"Numu with P = "<<trueP->at(i)<<std::endl;
		}

		if (truePDG->at(i) == 13){
			std::cout<<"Muon with P = "<<trueP->at(i)<<std::endl;
		}
		
		if (truePDG->at(i) == 211) {
			std::cout<<"Pionp with P = "<<trueP->at(i)<<std::endl;
		}
		if (truePDG->at(i) == -211){
			std::cout<<"Pionm with P = "<<trueP->at(i)<<std::endl;	
		}*///
	    }
	}

	std::cin.get();
	// look at neutrino energies and photon distribution in sigma0
   }   

  }

  double scale = 29/double(nfiles); // 1e21/3.4151e19 for total POT / POT per file

  std::cout<<"NCCQE Lambda =  "<<NCCQELambda<<" scaled to 1e21 POT = "<<NCCQELambda*scale<<std::endl;
  std::cout<<"NCCQE Sigma0 = "<<NCCQESigma0<<" scaled to 1e21 POT = "<<NCCQESigma0*scale<<std::endl;
   
  auto *c1 = new TCanvas("c1", "CCQE hyperons trueP", 200, 10, 700, 900);
  
  
  hCCQELambda->Scale(scale);
  hCCQESigma0->Scale(scale);
  hCCQELambda->Draw("HIST");
  c1->Print("plots/hCCQELambda.png");
  hCCQESigma0->Draw("HIST");
  c1->Print("plots/hCCQESigma0.png");
  hProton->Draw("HIST");
  c1->Print("plots/hProton.png");
  hPFParticles->Draw("HIST");
  c1->Print("plots/hPFParticles.png");
  hTrackCount->Draw("HIST");
  c1->Print("plots/hTrackCount.png");
  hShowerCount->Draw("HIST");
  c1->Print("plots/hShowerCount.png");
  hRecoVertexX->Draw("HIST");
  c1->Print("plots/hRecoVertexX.png");
  hRecoVertexY->Draw("HIST");
  c1->Print("plots/hRecoVertexY.png");
  hRecoVertexZ->Draw("HIST");
  c1->Print("plots/hRecoVertexZ.png");
  hTrackScores->Draw("HIST");
  c1->Print("plots/hTrackScores.png");
 // hPhoton->Draw("HIST");
 // c1->Print("plots/hPhoton.png");
 /* hProton->Draw();
  c1->Print("plots/hProton.png");
  hLambda->Draw();
  c1->Print("plots/hLambda.png");
  hSigma->Draw();
  c1->Print("hSigma.png");
  hMuon->Draw();
  c1->Print("hMuon.png");
  hPion->Draw();
  c1->Print("hPion.png");
*/

}
