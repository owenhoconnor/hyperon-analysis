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


   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // find PDG for each signal particle, fill TEfficiency for each
      // isReconstructed->at(that PDG)

      std::cout<<"New event # "<<jentry<<std::endl;
      if (sampleType == 0){ // is Signal

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
            }
            // photons
            if (truePDG->at(i) == 22 && trueGeneration->at(i) == 1 && trueMotherPDG->at(i) == 3212){
               effPhoton->Fill(reconstructed, momentum);
               nPhotons++;
            }
            // protons
            if (truePDG->at(i) == 2212 && trueGeneration->at(i) == 2 && trueMotherPDG->at(i) == 3122){
               effProton->Fill(reconstructed, momentum);
               nProtons++;
            }
            // pions
            if (truePDG->at(i) == -211 && trueGeneration->at(i) == 2 && trueMotherPDG->at(i) == 3122){
               effPion->Fill(reconstructed, momentum);
               nPions++;
            }
         }

         // loop over reco PFPs
         std::cout<<"num of PFPs in event = "<<pfpKey->size()<<std::endl;
         for (int iPfp = 0; iPfp < pfpKey->size(); ++iPfp){
            std::cout<<"PFP truth matched PDG ="<<pfpTruePDG->at(iPfp)<<std::endl;

            const float trackScore = pfpTrackScore->at(iPfp);
            const bool inNuSlice = pfpIsNuSlice->at(iPfp);
            const bool isPrimary = pfpIsPrimary->at(iPfp);
            //if(trackScore<0){continue;}

            if(isPrimary){

               // muon
               if (pfpTruePDG->at(iPfp) == -13){
                  nPrimaryMuons++;
               }
               // photon
               if (pfpTruePDG->at(iPfp) == 22){
                  nPrimaryPhotons++;
               }
               // proton
               if (pfpTruePDG->at(iPfp) == 2212){
                  nPrimaryProtons++;
               }
               // pion
               if (pfpTruePDG->at(iPfp) == -211){
                  nPrimaryPhotons++;
               }
            }

            if(!isPrimary){
               // muon
               if (pfpTruePDG->at(iPfp) == -13){
                  nSecondaryMuons++;
               }
               // photon
               if (pfpTruePDG->at(iPfp) == 22){
                  nSecondaryPhotons++;
               }
               // proton
               if (pfpTruePDG->at(iPfp) == 2212){
                  nSecondaryProtons++;
               }
               // pion
               if (pfpTruePDG->at(iPfp) == -211){
                  nSecondaryPions++;
               }
            }

            // muon
            if (pfpTruePDG->at(iPfp) == -13){
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
            }

            const int matchedTrackID = pfpTrueTrackID->at(iPfp);
            for (int iTrue = 0; iTrue < trueTrackID->size(); ++iTrue){

               if(trueTrackID->at(iTrue) != matchedTrackID){
                  continue;
               }

               if (trueGeneration->at(iTrue) == 0){
                  // isPrimary

                  // muon
                  if (pfpTruePDG->at(iPfp) == -13){
                     //nPrimaryMuons++;
                  }
                  // photon
                  if (pfpTruePDG->at(iPfp) == 22){
                     //nPrimaryPhotons++;
                  }
                  // proton
                  if (pfpTruePDG->at(iPfp) == 2212){
                      //nPrimaryProtons++;
                  }
                  // pion
                  if (pfpTruePDG->at(iPfp) == -211){
                     //nPrimaryPions++;
                  }
               }
               else if (trueGeneration->at(iTrue) != 0){
                  // isSecondary
                  // muon
                  if (pfpTruePDG->at(iPfp) == -13){
                     //nSecondaryMuons++;
                  }
                  // photon
                  if (pfpTruePDG->at(iPfp) == 22){
                     //nSecondaryPhotons++;
                  }
                  // proton
                  if (pfpTruePDG->at(iPfp) == 2212){
                     // nSecondaryProtons++;
                  }
                  // pion
                  if (pfpTruePDG->at(iPfp) == -211){
                     //nSecondaryPions++;
                  }
               }
            }

         }
      }

   }

  /* float fracPrimaryMuon = nPrimaryMuons / (nPrimaryMuons + nSecondaryMuons);
   float fracPrimaryPhoton = nPrimaryPhotons / (nPrimaryPhotons + nSecondaryPhotons);
   float fracPrimaryProton = nPrimaryProtons / (nPrimaryProtons + nSecondaryProtons);
   float fracPrimaryPion = nPrimaryPions / (nPrimaryPions + nSecondaryPions);*/

   std::cout
   <<"num of muons = "<<nMuons
   <<"num of photons = "<<nPhotons
   <<"num of protons = "<<nProtons
   <<"num of pions = "<<nPions
   <<std::endl;

   /*std::cout<<"Fraction of Primary Muons = "<<fracPrimaryMuon
   <<" Fraction of Primary Photons = "<<fracPrimaryPhoton
   <<" Fraction of Primary Protons = "<<fracPrimaryProton
   <<" Fraction of Primary Pions = "<<fracPrimaryPion<<std::endl;*/

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
