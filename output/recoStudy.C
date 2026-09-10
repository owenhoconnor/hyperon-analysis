#define recoStudy_cxx
#include "recoStudy.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

TGraphAsymmErrors* MakeEfficiencyGraph(
    TEfficiency* eff,
    Color_t color,
    Style_t markerStyle)
{
    TGraphAsymmErrors* gr = eff->CreateGraph();

    for (int i = 0; i < gr->GetN(); ++i)
    {
        // Remove horizontal bin-width error bars
        gr->SetPointEXlow(i, 0.0);
        gr->SetPointEXhigh(i, 0.0);
    }

    gr->SetMarkerColor(color);
    gr->SetLineColor(color);

    gr->SetMarkerStyle(markerStyle);
    gr->SetMarkerSize(1.3);

    gr->SetLineWidth(2);

    gr->SetTitle(eff->GetTitle());

    return gr;
}

void recoStudy::Loop()
{
//   In a ROOT session, you can do:
//      root> .L recoStudy.C
//      root> recoStudy t
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

  // Create Histograms

   TEfficiency* effMuon = new TEfficiency("effMuon", "Muon Reconstructed Efficiency;True Momentum [GeV/c]; Efficiency", 30, 0.0, 3.0);
   TEfficiency* effPhoton = new TEfficiency("effPhoton", "Photon Reconstruction Efficiency;True Momentum [GeV/c]; Efficiency", 30, 0.0, 0.3);
   TEfficiency* effProton = new TEfficiency("effProton", "Proton Reconstruction Efficiency;True Momentum [GeV/c]; Efficiency", 30, 0.0, 1.5);
   TEfficiency* effPion = new TEfficiency("effPion", "Pion Reconstruction Efficiency;True Momentum [GeV/C]; Efficiency", 30, 0.0, 0.3);

   TH1F* hMuonTrackScore = new TH1F("hMuonTrackScore", "Muon Track Score", 100, 0, 1);
   TH1F* hPhotonTrackScore = new TH1F("hPhotonTrackScore", "Photon Track Score", 100, 0, 1);
   TH1F* hProtonTrackScore = new TH1F("hProtonTrackScore", "Proton Track Score", 100, 0, 1);
   TH1F* hPionTrackScore = new TH1F("hPionTrackScore", "Pion Track Score", 100, 0, 1);

   // Define counters

   int nMuons = 0;
   int nPhotons = 0;
   int nPions = 0;
   int nProtons = 0;

   int nPrimaryPFPs = 0;
   int nSecondaryPFPs = 0;

   int nPrimaryMuons = 0;
   int nSecondaryMuons = 0;
   int nPrimaryPhotons = 0;
   int nSecondaryPhotons = 0;
   int nPrimaryProtons = 0;
   int nSecondaryProtons = 0;
   int nPrimaryPions = 0;
   int nSecondaryPions = 0;

   int nNotRecoMuons = 0;
   int nNotRecoPhotons = 0;
   int nNotRecoProtons = 0;
   int nNotRecoPions = 0;

            int nMultiNuRootEvents = 0;


   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // find PDG for each signal particle, fill TEfficiency for each
      // isReconstructed->at(that PDG)

      int signalTrueMuonID = -999;
      int signalTruePhotonID = -999;
      int signalTrueProtonID = -999;
      int signalTruePionID = -999;

      std::cout<<"New event # "<<jentry<<std::endl;
      if (sampleType == 0){ // is Signal

         // loop over reco PFPs to determine nu PFP ID
         int nNuRoots = 0;
         int nuPfpSelfID = -1;

         for (size_t iPfp = 0; iPfp < pfpKey->size(); ++iPfp)
         {
            if (pfpIsNuSlice->at(iPfp) && pfpIsPrimary->at(iPfp))
            {
               ++nNuRoots;
               nuPfpSelfID = pfpSelfID->at(iPfp);
            }
         }

         if (nNuRoots != 1)
         {
            std::cout<< "Event " << jentry<< " has " << nNuRoots<< " primary PFPs in selected nu slice"<< std::endl;
            nMultiNuRootEvents++;
         }

         // loop over true particles

         for (int i = 0; i < truePDG->size(); ++i){
            std::cout<<"PDG = "<<truePDG->at(i)<<" at generation "<<trueGeneration->at(i)<<std::endl;

            const bool reconstructed = (trueIsReconstructed->at(i) == 1);
            const float momentum = trueP->at(i);
            
            // define reconstructability w/ per view hits
            int nGoodViews = 0;

            if (trueNHitsU->at(i) >= 5){++nGoodViews;}
            if (trueNHitsV->at(i) >= 5){++nGoodViews;}
            if (trueNHitsZ->at(i) >= 5){++nGoodViews;}

            const bool reconstructable = trueNHits->at(i) >= 15 && nGoodViews >= 2;

            // muons
            if (truePDG->at(i) == -13 && trueGeneration->at(i) == 0){
               effMuon->Fill(reconstructed, momentum);
               
               nMuons++;
               signalTrueMuonID = trueTrackID->at(i);

               if(!reconstructed){
                  nNotRecoMuons++;
               }
               else
               {
                  const int bestPfpIdx = trueBestRecoPfpIdx->at(i);
                  if(bestPfpIdx < 0 || bestPfpIdx > static_cast<int>(pfpParentID->size())){
                     std::cout<<"Error: muon says recoed, but valid best PFP index "<<bestPfpIdx<<" event = "<<jentry<<std::endl;
                  }
                  else
                  {
                     const int parentID = pfpParentID->at(bestPfpIdx);
                     const bool isPrimary = parentID == nuPfpSelfID && nuPfpSelfID >= 0;

                     if(isPrimary){
                        nPrimaryMuons++;
                     }
                     else{
                        nSecondaryMuons++;
                     }
                  }

               }
            }
            // photons
            if (truePDG->at(i) == 22 && trueGeneration->at(i) == 1 && trueMotherPDG->at(i) == 3212){
               effPhoton->Fill(reconstructed, momentum);
               
               nPhotons++;
               signalTruePhotonID = trueTrackID->at(i);

               if(!reconstructed){
                  nNotRecoPhotons++;
               }
               else
               {
                  const int bestPfpIdx = trueBestRecoPfpIdx->at(i);
                  if(bestPfpIdx < 0 || bestPfpIdx > static_cast<int>(pfpParentID->size())){
                     std::cout<<"Error: photon says recoed, but valid best PFP index "<<bestPfpIdx<<" event = "<<jentry<<std::endl;
                  }
                  else
                  {
                     const int parentID = pfpParentID->at(bestPfpIdx);
                     const bool isPrimary = parentID == nuPfpSelfID && nuPfpSelfID >= 0;

                     if(isPrimary){
                        nPrimaryPhotons++;
                     }
                     else{
                        nSecondaryPhotons++;
                     }
                  }

               }
            }
            // protons
            if (truePDG->at(i) == 2212 && trueGeneration->at(i) == 2 && trueMotherPDG->at(i) == 3122){
               effProton->Fill(reconstructed, momentum);
               
               nProtons++;
               if (signalTrueProtonID != -999) {
                  std::cout<< "WARNING: multiple signal proton candidates in event "<< jentry << std::endl;
               }
               signalTrueProtonID = trueTrackID->at(i);

               if(!reconstructed){
                  nNotRecoProtons++;
               }
               else
               {
                  const int bestPfpIdx = trueBestRecoPfpIdx->at(i);
                  if(bestPfpIdx < 0 || bestPfpIdx > static_cast<int>(pfpParentID->size())){
                     std::cout<<"Error: proton says recoed, but valid best PFP index "<<bestPfpIdx<<" event = "<<jentry<<std::endl;
                  }
                  else
                  {
                     const int parentID = pfpParentID->at(bestPfpIdx);
                     const bool isPrimary = parentID == nuPfpSelfID && nuPfpSelfID >= 0;

                     if(isPrimary){
                        nPrimaryProtons++;
                     }
                     else{
                        nSecondaryProtons++;
                     }
                  }
               }
            }
            // pions
            if (truePDG->at(i) == -211 && trueGeneration->at(i) == 2 && trueMotherPDG->at(i) == 3122){
               effPion->Fill(reconstructed, momentum);
               
               nPions++;
               signalTruePionID = trueTrackID->at(i);

               if(!reconstructed){
                  nNotRecoPions++;
               }
               else
               {
                  const int bestPfpIdx = trueBestRecoPfpIdx->at(i);
                  if(bestPfpIdx < 0 || bestPfpIdx > static_cast<int>(pfpParentID->size())){
                     std::cout<<"Error: pion says recoed, but valid best PFP index "<<bestPfpIdx<<" event = "<<jentry<<std::endl;
                  }
                  else
                  {
                     const int parentID = pfpParentID->at(bestPfpIdx);
                     const bool isPrimary = parentID == nuPfpSelfID && nuPfpSelfID >= 0;

                     if(isPrimary){
                        nPrimaryPions++;
                     }
                     else{
                        nSecondaryPions++;
                     }
                  }

               }
            }

         } // end loop over truth 

         // loop over reco PFPs
      
         std::cout<<"num of PFPs in event = "<<pfpKey->size()<<std::endl;
         for (int iPfp = 0; iPfp < pfpKey->size(); ++iPfp){
            std::cout<<"PFP truth matched PDG ="<<pfpTruePDG->at(iPfp)<<std::endl;

            const float trackScore = pfpTrackScore->at(iPfp);
            const bool inNuSlice = pfpIsNuSlice->at(iPfp);
            const int parentID = pfpParentID->at(iPfp);
            const bool isPrimary = nuPfpSelfID >= 0 && parentID == nuPfpSelfID;
            const int matchedTrackID = pfpTrueTrackID->at(iPfp);
            // if this pfp->Parent is primary, need to be able to look up isPrimary for parent -> make map
            //if(trackScore<0){continue;}

            if(isPrimary){nPrimaryPFPs++;}
            if(!isPrimary){nSecondaryPFPs++;}

            // muon
            /*if (pfpTruePDG->at(iPfp) == -13){
               hMuonTrackScore->Fill(trackScore);
            }
            // photon
            if (pfpTruePDG->at(iPfp) == 22){
               hPhotonTrackScore->Fill(trackScore);
            }
            // proton
            if (pfpTruePDG->at(iPfp) == 2212){
               hProtonTrackScore->Fill(trackScore);
            }
            // pion
            if (pfpTruePDG->at(iPfp) == -211){
               hPionTrackScore->Fill(trackScore);
            }*/

         } // end loop over PFPs

      } // end check if event is signal

   } // end loop over events

   float fracPrimaryMuon = static_cast<float>(nPrimaryMuons) / nMuons;//(static_cast<float>(nPrimaryMuons) + static_cast<float>(nSecondaryMuons));
   float fracPrimaryPhoton = static_cast<float>(nPrimaryPhotons) / nPhotons; //(static_cast<float>(nPrimaryPhotons) + static_cast<float>(nSecondaryPhotons));
   float fracPrimaryProton = static_cast<float>(nPrimaryProtons) / nProtons; //(static_cast<float>(nPrimaryProtons) + static_cast<float>(nSecondaryProtons));
   float fracPrimaryPion = static_cast<float>(nPrimaryPions) / nPions; //(static_cast<float>(nPrimaryPions) + static_cast<float>(nSecondaryPions));

   float fracSecondaryMuon = static_cast<float>(nSecondaryMuons) / nMuons; // (static_cast<float>(nPrimaryMuons) + static_cast<float>(nSecondaryMuons));
   float fracSecondaryPhoton = static_cast<float>(nSecondaryPhotons) / nPhotons; // (static_cast<float>(nPrimaryPhotons) + static_cast<float>(nSecondaryPhotons));
   float fracSecondaryProton = static_cast<float>(nSecondaryProtons) / nProtons; // (static_cast<float>(nPrimaryProtons) + static_cast<float>(nSecondaryProtons));
   float fracSecondaryPion = static_cast<float>(nSecondaryPions) / nPions; // (static_cast<float>(nPrimaryPions) + static_cast<float>(nSecondaryPions));

   float fracNotRecoMuon = static_cast<float>(nNotRecoMuons) / nMuons; // (static_cast<float>(nPrimaryMuons) + static_cast<float>(nSecondaryMuons));
   float fracNotRecoPhoton = static_cast<float>(nNotRecoPhotons) / nPhotons; // (static_cast<float>(nPrimaryPhotons) + static_cast<float>(nSecondaryPhotons));
   float fracNotRecoProton = static_cast<float>(nNotRecoProtons) / nProtons; // (static_cast<float>(nPrimaryProtons) + static_cast<float>(nSecondaryProtons));
   float fracNotRecoPion = static_cast<float>(nNotRecoPions) / nPions; // (static_cast<float>(nPrimaryPions) + static_cast<float>(nSecondaryPions));



   std::cout
   <<"num of muons = "<<nMuons
   <<"num of photons = "<<nPhotons
   <<"num of protons = "<<nProtons
   <<"num of pions = "<<nPions
   <<std::endl;

   std::cout<<"============================================="<<std::endl;
   std::cout<<"Signal Particle | Frac Primary | Frac Secondary | Not Reco"<<std::endl;
   std::cout<<"Muon            | "<<fracPrimaryMuon<<" | "<<fracSecondaryMuon<<" | "<<fracNotRecoMuon<<std::endl;
   std::cout<<"Photon          | "<<fracPrimaryPhoton<<" | "<<fracSecondaryPhoton<<"| "<<fracNotRecoPhoton<<std::endl;
   std::cout<<"Proton          | "<<fracPrimaryProton<<" | "<<fracSecondaryProton<<" | "<<fracNotRecoProton<<std::endl;
   std::cout<<"Pion            | "<<fracPrimaryPion<<" | "<<fracSecondaryPion<<" | "<<fracNotRecoPion<<std::endl;
   std::cout<<"============================================="<<std::endl;

   std::cout<<"num of primary muons = "<<nPrimaryMuons<<std::endl;
   std::cout<<"num of secondary muons = "<<nSecondaryMuons<<std::endl;

   std::cout<<"num of primary photons = "<<nPrimaryPhotons<<std::endl;
   std::cout<<"num of secondary photons = "<<nSecondaryPhotons<<std::endl;
   std::cout<<" num of primary protons = "<<nPrimaryProtons<<std::endl;
   std::cout<<"num of secondary protons = "<<nSecondaryProtons<<std::endl;
   std::cout<<"num of primray pions = "<<nPrimaryPions<<std::endl;
   std::cout<<"num of secondary pions = "<<nSecondaryPions<<std::endl;
   std::cout<<"num of primary PFPs = "<<nPrimaryPFPs<<std::endl;
   std::cout<<"num of secondary PFPs = "<<nSecondaryPFPs<<std::endl;

   // -------------------------------------------------------
   // Track Score by Signal Particle
   // -------------------------------------------------------

   THStack* hsTrackScores = new THStack("hsTrackScores", "Track Scores by Signal Particles");
   TLegend* legTrackScore = new TLegend(0.75, 0.3, 0.97, 0.91);
   legTrackScore->SetBorderSize(0);
   legTrackScore->SetFillStyle(0);
   legTrackScore->SetTextSize(0.025);

   hMuonTrackScore->SetFillColorAlpha(kRed, 0.75);
   hMuonTrackScore->SetLineColor(kBlack);
   hsTrackScores->Add(hMuonTrackScore);
   legTrackScore->AddEntry(hMuonTrackScore, "Muons", "f");

   hPhotonTrackScore->SetFillColorAlpha(kBlue, 0.75);
   hPhotonTrackScore->SetLineColor(kBlack);
   hsTrackScores->Add(hPhotonTrackScore);
   legTrackScore->AddEntry(hPhotonTrackScore, "Photons", "f");

   hProtonTrackScore->SetFillColorAlpha(kYellow, 0.75);
   hProtonTrackScore->SetLineColor(kBlack);
   hsTrackScores->Add(hProtonTrackScore);
   legTrackScore->AddEntry(hProtonTrackScore, "Protons", "f");

   hPionTrackScore->SetFillColorAlpha(kOrange, 0.75);
   hPionTrackScore->SetLineColor(kBlack);
   hsTrackScores->Add(hPionTrackScore);
   legTrackScore->AddEntry(hPionTrackScore, "Pions", "f");

   TCanvas* cTrackScores = new TCanvas("cTrackScores", "Track Scores by Signal Particle", 1800, 1200);
   hsTrackScores->Draw("HIST");
   hsTrackScores->GetXaxis()->SetTitle(" PFP Track Score");
   hsTrackScores->GetYaxis()->SetTitle("PFPs");

   legTrackScore->Draw();
   cTrackScores->Print("plots/trackScoresStacked.png");


   // -------------------------------------------------------
   // Reconstruction Efficiencies
   // -------------------------------------------------------

   TCanvas* cEffMuon = new TCanvas("cEffMuon", "", 1800, 1200);
  cEffMuon->SetFillColor(kWhite);

   TGraphAsymmErrors* grMuon =
      MakeEfficiencyGraph(
         effMuon,
         kBlue + 1,
         20
      );

   grMuon->Draw("APZ");

   grMuon->GetXaxis()->SetLimits(0.0, 3.0);
   grMuon->GetYaxis()->SetRangeUser(0.0, 1.05);

   cEffMuon->Print("plots/efficiencyMuon.png");

   TCanvas* cEffPhoton =
      new TCanvas("cEffPhoton", "", 1800, 1200);

   TGraphAsymmErrors* grPhoton =
      MakeEfficiencyGraph(
         effPhoton,
         kBlue + 1,
         20
      );

   grPhoton->Draw("APZ");

   grPhoton->GetXaxis()->SetLimits(0.0, 0.3);
   grPhoton->GetYaxis()->SetRangeUser(0.0, 1.05);

   cEffPhoton->Print("plots/efficiencyPhoton.png");

   TCanvas* cEffProton = new TCanvas("cEffProton", "", 1800, 1200);

   TGraphAsymmErrors* grProton =
      MakeEfficiencyGraph(
         effProton,
         kBlue + 1,
         20
      );

   grProton->Draw("APZ");

   grProton->GetXaxis()->SetLimits(0.0, 1.5);
   grProton->GetYaxis()->SetRangeUser(0.0, 1.05);
   cEffProton->Print("plots/efficiencyProton.png");

   TCanvas* cEffPion = new TCanvas("cEffPion", "", 1800, 1200);
   cEffPion->SetFillStyle(1001);
   cEffPion->SetFillColor(kWhite);
   
   TGraphAsymmErrors* grPion =
      MakeEfficiencyGraph(
         effPion,
         kBlue + 1,
         20
      );
   grPion->Draw("APZ");
   cEffPion->Print("plots/efficiencyPion.png");
}
