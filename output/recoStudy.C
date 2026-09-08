#define recoStudy_cxx
#include "recoStudy.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void StyleEfficiency(TEfficiency* eff, Color_t color, Style_t markerStyle, double xmin, double xmax)
{
    eff->Draw("AP");
    gPad->Update();

    TGraphAsymmErrors* gr = eff->GetPaintedGraph();

    for (int i = 0; i < gr->GetN(); ++i)
    {
        gr->SetPointEXlow(i, 0.0);
        gr->SetPointEXhigh(i, 0.0);
    }

    gr->SetMarkerColor(color);
    gr->SetLineColor(color);
    gr->SetMarkerStyle(markerStyle);
    gr->SetMarkerSize(5);
    gr->SetLineWidth(2);

    gr->GetXaxis()->SetRangeUser(xmin, xmax);
    gr->GetYaxis()->SetRangeUser(0.0, 1.05);

    gPad->Modified();
    gPad->Update();
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

  /* std::array<TH1F*, 4> hRecoEfficiencies;

   for (int i = 0; i < 4; ++i){
      hRecoEfficiencies->at(i) = new TH1F(Form("hRecoeff_%d", i), "", 100, 0, 0);
      hRecoEfficiencies->at(i)->SetDirectory(nullptr);
   }*/

   TEfficiency* effMuon = new TEfficiency("effMuon", "Muon Reconstructed Efficiency;True Momentum [GeV/c]; Efficiency", 20, 0.0, 3.0);
   TEfficiency* effPhoton = new TEfficiency("effPhoton", "Photon Reconstruction Efficiency;True Momentum [GeV/c]; Efficiency", 20, 0.0, 0.5);
   TEfficiency* effProton = new TEfficiency("effProton", "Proton Reconstruction Efficiency;True Momentum [GeV/c]; Efficiency", 20, 0.0, 1.5);
   TEfficiency* effPion = new TEfficiency("effPion", "Pion Reconstruction Efficiency;True Momentum [GeV/C]; Efficiency", 20, 0.0, 0.5);

   int nMuons = 0;
   int nPhotons = 0;
   int nPions = 0;
   int nProtons = 0;


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

         for (int i = 0; i < truePDG->size(); ++i){
            std::cout<<"PDG = "<<truePDG->at(i)<<" at generation "<<trueGeneration->at(i)<<std::endl;

            const bool reconstructed = (trueIsReconstructed->at(i) == 1);
            const float momentum = trueP->at(i);
            
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
      }

   }

   std::cout
   <<"num of muons = "<<nMuons
   <<"num of photons = "<<nPhotons
   <<"num of protons = "<<nProtons
   <<"num of pions = "<<nPions
   <<std::endl;

   TCanvas* cEffMuon = new TCanvas("cEffMuon", "", 1800, 1200);
   cEffMuon->SetFillStyle(1001);
   cEffMuon->SetFillColor(kWhite);
   effMuon->Draw("AP");
   StyleEfficiency(effMuon, kBlue, 20, 0, 3.0);
   cEffMuon->Print("plots/efficiencyMuon.png");

   TCanvas* cEffPhoton = new TCanvas("cEffPhoton", "", 1800, 1200);
   cEffPhoton->SetFillStyle(1001);
   cEffPhoton->SetFillColor(kWhite);
   effPhoton->Draw("AP");
   StyleEfficiency(effPhoton, kBlue, 20, 0, 1.0);
   cEffPhoton->Print("plots/efficiencyPhoton.png");

   TCanvas* cEffProton = new TCanvas("cEffProton", "", 1800, 1200);
   cEffProton->SetFillStyle(1001);
   cEffProton->SetFillColor(kWhite);
   effProton->Draw("AP");
   StyleEfficiency(effProton, kBlue, 20, 0, 2.0);
   cEffProton->Print("plots/efficiencyProton.png");

   TCanvas* cEffPion = new TCanvas("cEffPion", "", 1800, 1200);
   cEffPion->SetFillStyle(1001);
   cEffPion->SetFillColor(kWhite);
   StyleEfficiency(effPion, kBlue, 20, 0, 1.0);
   cEffPion->Print("plots/efficiencyPion.png");
}
