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
	
   if (fChain == 0) return;
   
   int nSigma = 0;
   int nLambda = 0;
   int nGoodLambda = 0;
   int nGoodSigma = 0;

   TString Samples[3];
   Samples[0] = "Sigma";
   Samples[1] = "Hyperons";
   Samples[2] = "Background";
   
   int nEvents[3] = {0};
   nEvents[0] = 141422; // number of events in all hyperon files
   nEvents[1] = nEvents[0]; // hyperons
   nEvents[2] = 209009 + nEvents[0]; // num of events in bkg files + num of events in hyp files
   

   double PoT[4] = {0};
   PoT[0] = 2.69e19; // pot per hyperon file
   PoT[1] = PoT[0]; //  (hyperons)
   PoT[2] = 8.875e15; // pot per bkg file
   PoT[3] = 1e21; // total pot

   double nHyperonFiles = 1776.;
   double nBkgFiles = 2236.;

   double scale[3] = {0};
   scale[0] = 0.035131; //PoT[3] / (nHyperonFiles * PoT[0]); // sigma
   scale[1] = 0.035131; //PoT[3] / (nHyperonFiles * PoT[1]); // hyperons 
   scale[2] = PoT[3] / (nBkgFiles * PoT[2]); // bkg
   std::cout<<"SIGMA and HYPERON SCALE = "<<scale[0]<<std::endl;
   std::cout<<"BACKGROUND SCALE = "<<scale[2]<<std::endl;

   const int nCuts = 4;
   TString CutName[nCuts];
   CutName[0] = "AllEvents";
   CutName[1] = "FV";
   CutName[2] = "nShowers";
   CutName[3] = "nTracks";

   const int nVars = 15;
   TString VarName[nVars];
   VarName[0] = "nTracks";
   VarName[1] = "nShowers";
   VarName[2] = "trackScores";

   VarName[3] = "muonTrackScore";
   VarName[4] = "muonTrackLength";
   VarName[5] = "muonTrackDistToVtx";

   VarName[6] = "pionTrackScore";
   VarName[7] = "pionTrackLength";
   VarName[8] = "pionTrackDistToVtx";

   VarName[9] = "protonTrackScore";
   VarName[10] = "protonTrackLength";
   VarName[11] = "protonTrackDistToVtx";
   
   VarName[12] = "gammaTrackScore";
   VarName[13] = "gammaTrackLength";
   VarName[14] = "gammaTrackDistToVtx";

   TH1F *Hist[nVars][3][nCuts];

   int MaxTracks = 6;
   int MaxShowers = 6;

   TH2F *hTracksShowersSig = new TH2F("hTracksShowersSig", ";num of Tracks;num of Showers",
		   MaxTracks + 1, -0.5, MaxTracks + 0.5,
		   MaxShowers + 1, -0.5, MaxShowers + 0.5);
   hTracksShowersSig->SetDirectory(nullptr);

   TH2F *hTracksShowersBkg = new TH2F("nTracksShowersBkg", ";num of Tracks;num of Showers",
		   MaxTracks + 1, -0.5, MaxTracks + 0.5,
		   MaxShowers + 1, -0.5, MaxShowers + 0.5);
   hTracksShowersBkg->SetDirectory(nullptr); 

   for (int s = 0; s < 3; s++){ // samples
	  for (int v = 0; v < nVars; v++){ // variables
		 for (int c = 0; c < nCuts; c++){ // cuts
			Hist[v][s][c] = new TH1F(VarName[v] + Samples[s] + CutName[c], "", 100, -1, -1);
		 }
	  }
   }
   

   TFile *sigFile = TFile::Open("TreeS.root", "RECREATE");

   if (!sigFile || sigFile->IsZombie()){
	   std::cerr<<"Could not open file!"<<std::endl;
	   return;
   }

   sigFile->cd();
   TTree *signalTree = fChain->CloneTree(0);
   signalTree->SetName("preTree");
   signalTree->SetDirectory(sigFile);

   TFile *bkgFile = TFile::Open("TreeB.root", "RECREATE");

   if (!bkgFile || bkgFile->IsZombie()){
	   std::cerr<<"Could not open file!"<<std::endl;
	   return;
   }

   bkgFile->cd();
   TTree *bkgTree = fChain->CloneTree(0);
   bkgTree->SetName("preTree");
   bkgTree->SetDirectory(bkgFile);

   int sampleType;
   signalTree->Branch("sampleType", &sampleType);
   bkgTree->Branch("sampleType", &sampleType);

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
      bool IsSignal = false;
   
      bool IsInFV = false; // true vertex FV for signal def
      bool IsInRecoFV = false;
      bool IsInTrackRange = false;
      bool IsInShowerRange = false;
      bool IsBadShower = false;
      bool SelMuonCandidate = false;
      bool SelPionCandidate = false;
      bool SelProtonCandidate = false;
      int muonIndex = -1;
      int pionIndex = -1;
      int protonIndex = -1;
      int shower1Index = -1;
      int shower2Index = -1;
      bool Cuts[nCuts] = {false};
      int s = -1; // sample index (0 = sigma signal, 1 = other hyperons, 2=bkg)

      int sumCounter = 0;
      int daughterCounter = 0;

      std::cout<<"****************** new event ********************"<<std::endl;
      std::cout<<"trackStartPosition size = "<<trackStartPositionX->size()<<std::endl;
      std::cout<<"distance to reco vertex size = "<<DistanceToRecoVertex->size()<<std::endl;
      std::cout<<"TrackLengths size = "<<trackLengths->size()<<std::endl;
      std::cout<<"trackScores size = "<<trackScores->size()<<std::endl;
      std::cout<<"trackCount = "<<trackCount<<" showerCount = "<<showerCount<<std::endl;
      std::cout<<"pfpTrackPDG size = "<<pfpTrackPDG->size()<<std::endl;
      std::cout<<"pfpShowerPDG size = "<<pfpShowerPDG->size()<<std::endl;
      std::cout<<"nPFParticles = "<<nPFParticles<<std::endl;
      std::cout<<"pfpPDG size = "<<pfpPDG->size()<<std::endl;

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
	    }
	   // if (IsLambda){
	//	nLambda++;
	//	std::cout<<"lambda event, nLambda = "<<nLambda<<std::endl; // check daughters of lambda for pion- and proton
	//	if (IsGoodLambda){
	//		nGoodLambda++;
	//		std::cout<<"lambda w proton/pion daughters event, n goodlambda = "<<nGoodLambda<<std::endl;
	//	}
	  //  }
	}

      if(IsSignal && jentry < nEvents[0] + 1){ // if signal and we're still within range of hyperon file events
	  s=0;
	  sampleType=0;
	  hTracksShowersSig->Fill(trackCount, showerCount);
      }
      if(!IsSignal && jentry < nEvents[0] + 1){ // if not signal but we're still within range of hyperon file events (other hyperons)
	      s=1;
	      sampleType=1;
      }
      if(!IsSignal && jentry > nEvents[0]){ // now within range of bkg events, exclude rare signal
          s=2;
	  sampleType=2;
	  hTracksShowersBkg->Fill(trackCount, showerCount);
      }	
      if(IsSignal && jentry > nEvents[0]){continue;} // if there's signal in bkg range, ignore

      // Preselection

      TVector3 track1Start;
      TVector3 track2Start;
      TVector3 track3Start;
      TVector3 shower1Start;

      /*if(muonIndex>-1){muonStart.SetXYZ(TrackStartPositionX->at(muonIndex), TrackStartPositionY->at(muonIndex), TrackStartPositionZ->at(muonIndex));}
      if(pionIndex>-1){pionStart.SetXYZ(TrackStartPositionX->at(pionIndex), TrackStartPositionY->at(pionIndex), TrackStartPositionZ->at(pionIndex));}
      if(protonIndex>-1)  protonStart.SetXYZ(trk_start_x->at(proton_index), trk_start_y->at(proton_index),trk_start_z->at(proton_index));
    if(shower1Index>-1) shower1Start.SetXYZ(trk_start_x->at(shower1_index), trk_start_y->at(shower1_index),trk_start_z->at(shower1_index));
       if(shower2Index>-1) shower2Start.SetXYZ(trk_start_x->at(shower2_index), trk_start_y->at(shower2_index),trk_start_z->at(shower2_index));
*/
      if (std::abs(RecoVertexX) < 180 && std::abs(RecoVertexY) < 180 && RecoVertexZ > 10 && RecoVertexZ < 450){IsInRecoFV = true;}

      if (trackCount > 2 && trackCount < 6){IsInTrackRange = true;} // need to fill TMVA trees for specific topology
      if (showerCount == 1 || showerCount == 2){IsInShowerRange = true;}

      // Select topology:

      /*for (int i = 0; i < showerLengths->size(); i++){
	      if(showerLengths->at(i) < 0){IsBadShower = true;}
	      if(std::abs(showerStartPositionX->at(i)) > 180){IsBadShower = true;}
	      if(std::abs(showerStartPositionY->at(i)) > 180){IsBadShower = true;}
	      if(showerStartPositionZ->at(i) < 0){IsBadShower = true;} 
      }*/
      if (IsInRecoFV && trackCount == 3 && showerCount == 1 && !IsBadShower){
	      if(s==0){signalTree->Fill();}
	      if(s==2){bkgTree->Fill();}
      }

      Cuts[0] = true;
      Cuts[1] = Cuts[0] && IsInRecoFV;
      Cuts[2] = Cuts[1] && IsInTrackRange;
      Cuts[3] = Cuts[2] && IsInShowerRange;

    /*  for (int c = 0; c<nCuts; c++){

	if(Cuts[c]){

		Hist[0][s][c]->Fill(trackCount, scale[s]);
		Hist[1][s][c]->Fill(showerCount, scale[s]);

		for (int i = 0; i < trackScores->size(); i++){

			// Fill non-particle specific per PFP variables (e.g trackScore)
			Hist[2][s][c]->Fill(trackScores->at(i), scale[s]);


			if(pfpPDG->size() == trackScores->size()) {
			// Fill particle specific per PFP variables (need pfpPDG to be same size as trackScores)
			if(pfpPDG->at(i) == -13 && TrackLengths->at(i) > 0){ // Antimuon
				Hist[3][s][c]->Fill(trackScores->at(i), 1);
				Hist[4][s][c]->Fill(TrackLengths->at(i), 1);
				Hist[5][s][c]->Fill(DistanceToRecoVertex->at(i), 1);
			}
			if(pfpPDG->at(i) == -211 && TrackLengths->at(i) > 0){ // Pion (minus)
				Hist[6][s][c]->Fill(trackScores->at(i), 1);
				Hist[7][s][c]->Fill(TrackLengths->at(i), 1);
				Hist[8][s][c]->Fill(DistanceToRecoVertex->at(i), 1);
			}
			if(pfpPDG->at(i) == 2212 && TrackLengths->at(i) > 0){ //Proton
				Hist[9][s][c]->Fill(trackScores->at(i), 1);
				Hist[10][s][c]->Fill(TrackLengths->at(i), 1);
				Hist[11][s][c]->Fill(DistanceToRecoVertex->at(i), 1);
			}
			if(pfpPDG->at(i) == 22 && TrackLengths->at(i) > 0){ // Photon
				Hist[12][s][c]->Fill(trackScores->at(i), 1);
				Hist[13][s][c]->Fill(TrackLengths->at(i), 1);
				Hist[14][s][c]->Fill(DistanceToRecoVertex->at(i), 1);
			}
		    }

		}
	   }

      } // end loop over cuts
*/
     
   } // end event loop

   TCanvas *c1 = new TCanvas("c1", "Track-Shower Topology", 2000, 2000);
   c1->SetRightMargin(0.15);
   c1->SetLeftMargin(0.15);
   c1->SetBottomMargin(0.15);

   hTracksShowersSig->GetXaxis()->SetTitle("Number of tracks");
   hTracksShowersSig->GetYaxis()->SetTitle("Number of showers");
   hTracksShowersSig->GetZaxis()->SetTitle("Events");

   hTracksShowersSig->GetXaxis()->CenterTitle();
   hTracksShowersSig->GetYaxis()->CenterTitle();

   hTracksShowersSig->GetXaxis()->SetTitleSize(0.04);
   hTracksShowersSig->GetYaxis()->SetTitleSize(0.04);
   hTracksShowersSig->GetZaxis()->SetTitleSize(0.04);

   hTracksShowersSig->GetXaxis()->SetLabelSize(0.04);
   hTracksShowersSig->GetYaxis()->SetLabelSize(0.04);
   hTracksShowersSig->GetZaxis()->SetLabelSize(0.04);

   hTracksShowersSig->GetXaxis()->SetNdivisions(MaxTracks + 1, 0, 0, kTRUE);
   hTracksShowersSig->GetYaxis()->SetNdivisions(MaxShowers + 1, 0, 0, kTRUE);

   gStyle->SetPalette(kViridis);
   hTracksShowersSig->Draw("COLZ TEXT");
   c1->Print("plots/nTracksShowersSig.png");

   TCanvas *c2 = new TCanvas("c2", "Background Track-Shower Topology", 2000, 2000);
   c2->SetRightMargin(0.15);
   c2->SetLeftMargin(0.15);
   c2->SetBottomMargin(0.15);

   hTracksShowersBkg->GetXaxis()->SetTitle("Number of tracks");
   hTracksShowersBkg->GetYaxis()->SetTitle("Number of showers");
   hTracksShowersBkg->GetZaxis()->SetTitle("Events");

   hTracksShowersBkg->GetXaxis()->CenterTitle();
   hTracksShowersBkg->GetYaxis()->CenterTitle();

   hTracksShowersBkg->GetXaxis()->SetTitleSize(0.04);
   hTracksShowersBkg->GetYaxis()->SetTitleSize(0.04);
   hTracksShowersBkg->GetZaxis()->SetTitleSize(0.04);

   hTracksShowersBkg->GetXaxis()->SetLabelSize(0.04);
   hTracksShowersBkg->GetYaxis()->SetLabelSize(0.04);
   hTracksShowersBkg->GetZaxis()->SetLabelSize(0.04);

   hTracksShowersBkg->GetXaxis()->SetNdivisions(MaxTracks + 1, 0, 0, kTRUE);
   hTracksShowersBkg->GetYaxis()->SetNdivisions(MaxShowers + 1, 0, 0, kTRUE);

   gStyle->SetPalette(kViridis);
   hTracksShowersBkg->Draw("COLZ TEXT");
   c2->Print("plots/nTracksShowersBkg.png");

   sigFile->cd();
   signalTree->Write("preTree"); // Write signal tree and close file
   sigFile->Close();
   delete sigFile;

   bkgFile->cd();
   bkgTree->Write("preTree"); // write bkg tree
   bkgFile->Close();
   delete bkgFile;

/*   TCanvas *c2 = new TCanvas("c2","",5000,2000); // Print results to a canvas

   for(int v = 0; v< nVars; v++){
    c2->Divide(nCuts,2);// 2 rows and 4 columns

    for (int c = 0; c < nCuts; c++){
      for(int s = 0; s < 2; s++){
	c2->cd(c+1+s*nCuts);
	Hist[v][s][c]->Draw();
      }
    }
    c2->Print("plots/"+VarName[v]+".png");
    c2->Clear();
  }
*/
  // TCanvas *c3 = new TCanvas("c3", "", 3000, 3000);
   //hTracksShowers->Draw("text");
   //c3->Print("plots/ntracks_v_nshowers.png");

   std::cout<<"num of good sigma: "<<nGoodSigma<<std::endl;
   std::cout<<"sigma scale = "<<scale[0]<<std::endl;
   std::cout<<"bkg scale = "<<scale[2]<<std::endl;
}
