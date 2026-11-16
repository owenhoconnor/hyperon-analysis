#define recoStudy_cxx
#include "recoStudy.h"

#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <TAxis.h>

#include <array>
#include <iostream>
#include <string>
#include <vector>

namespace {

// Truth-level reconstructability: at least five hits in EACH wire plane.
bool IsReconstructable(int nHits, int nU, int nV, int nZ)
{
    (void)nHits; // Per-plane requirement already implies >= 15 total hits.
    return nU >= 5 && nV >= 5 && nZ >= 5;
}

struct ParticleStudy {
    const char* name;
    const char* filename;
    double maxMomentum;
    TEfficiency* all;
    TEfficiency* reconstructable;
    long long nTrue = 0;
    long long nReconstructable = 0;
    long long nMatched = 0;
    long long nMatchedReconstructable = 0;
    long long nInvalidBestPFP = 0;
    long long nPrimary = 0;
    long long nSecondary = 0;
};

TGraphAsymmErrors* MakeEfficiencyGraph(TEfficiency* eff, Color_t colour,
                                       Style_t marker)
{
    TGraphAsymmErrors* graph = eff->CreateGraph();
    graph->SetMarkerColor(colour);
    graph->SetLineColor(colour);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(1.2);
    graph->SetLineWidth(2);
    for (int i = 0; i < graph->GetN(); ++i) {
        graph->SetPointEXlow(i, 0.0);
        graph->SetPointEXhigh(i, 0.0);
    }
    return graph;
}

void PlotEfficiencies(const ParticleStudy& study)
{
    TCanvas* canvas = new TCanvas(Form("cEff%s", study.name), study.name, 1400, 1000);
    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.05);
    canvas->SetBottomMargin(0.12);
    canvas->SetTopMargin(0.08);

    TGraphAsymmErrors* grAll = MakeEfficiencyGraph(study.all, kBlue + 1, 20);
    TGraphAsymmErrors* grReco = MakeEfficiencyGraph(study.reconstructable, kOrange + 7, 22);

    grAll->SetTitle(Form("%s Reconstruction Study;True Momentum [GeV/c];Fraction of True Particles", study.name));
    grAll->Draw("APZ");
    grAll->GetXaxis()->SetLimits(0.0, study.maxMomentum);
    grAll->GetYaxis()->SetRangeUser(0.0, 1.05);
    grAll->GetXaxis()->SetTitleSize(0.045);
    grAll->GetYaxis()->SetTitleSize(0.045);
    grAll->GetXaxis()->SetLabelSize(0.04);
    grAll->GetYaxis()->SetLabelSize(0.04);
    grAll->GetYaxis()->SetTitleOffset(1.15);

    grReco->Draw("PZ SAME");

    TLegend* legend = new TLegend(0.51, 0.15, 0.91, 0.29);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
    legend->AddEntry(grAll, "Reconstructed", "lp");
    legend->AddEntry(grReco, "Reconstructable (5 hits/plane)", "lp");
    legend->Draw();

    canvas->Modified();
    canvas->Update();
    canvas->Print(Form("plots/%s", study.filename));
    delete canvas;
}

double Fraction(long long numerator, long long denominator)
{
    return denominator ? static_cast<double>(numerator) / denominator : 0.0;
}

} // namespace

void recoStudy::Loop()
{
    if (!fChain) return;

    gSystem->mkdir("plots", true);
    gStyle->SetOptStat(0);

    std::array<ParticleStudy, 4> studies = {{
        {"Muon",   "efficiencyMuon.png",   3.0, nullptr, nullptr},
        {"Photon", "efficiencyPhoton.png", 0.3, nullptr, nullptr},
        {"Proton", "efficiencyProton.png", 1.5, nullptr, nullptr},
        {"Pion",   "efficiencyPion.png",   0.3, nullptr, nullptr}
    }};

    for (auto& study : studies) {
        study.all = new TEfficiency(Form("eff%sAll", study.name), "",
                                    30, 0.0, study.maxMomentum);
        study.reconstructable = new TEfficiency(Form("eff%sReconstructable", study.name), "",
                                                30, 0.0, study.maxMomentum);
    }

    long long nSignalEvents = 0;
    long long nUnusualNuRoots = 0;
    long long nBadTruthVectors = 0;
    const Long64_t nEntries = fChain->GetEntriesFast();

    for (Long64_t jentry = 0; jentry < nEntries; ++jentry) {
        if (LoadTree(jentry) < 0) break;
        fChain->GetEntry(jentry);
        if (sampleType != 0) continue;
        ++nSignalEvents;

        // Identify the selected neutrino PFP root, as in the original script.
        int nNuRoots = 0;
        int nuPfpSelfID = -1;
        for (size_t i = 0; i < pfpKey->size(); ++i) {
            if (pfpIsNuSlice->at(i) && pfpIsPrimary->at(i)) {
                ++nNuRoots;
                nuPfpSelfID = pfpSelfID->at(i);
            }
        }
        if (nNuRoots != 1) ++nUnusualNuRoots;

        const size_t nTruth = truePDG->size();
        if (trueGeneration->size() != nTruth || trueMotherPDG->size() != nTruth ||
            trueP->size() != nTruth || trueIsReconstructed->size() != nTruth ||
            trueNHits->size() != nTruth || trueNHitsU->size() != nTruth ||
            trueNHitsV->size() != nTruth || trueNHitsZ->size() != nTruth ||
            trueBestRecoPfpIdx->size() != nTruth) {
            ++nBadTruthVectors;
            continue;
        }

        for (size_t i = 0; i < nTruth; ++i) {
            const int pdg = truePDG->at(i);
            const int generation = trueGeneration->at(i);
            const int mother = trueMotherPDG->at(i);

            // Original signal-particle definitions; do not loosen these.
            int category = -1;
            if (pdg == -13 && generation == 0) category = 0;
            else if (pdg == 22 && generation == 1 && mother == 3212) category = 1;
            else if (pdg == 2212 && generation == 2 && mother == 3122) category = 2;
            else if (pdg == -211 && generation == 2 && mother == 3122) category = 3;
            if (category < 0) continue;

            ParticleStudy& study = studies.at(category);
            const bool matched = (trueIsReconstructed->at(i) == 1);
            const bool reconstructable = IsReconstructable(
                trueNHits->at(i), trueNHitsU->at(i),
                trueNHitsV->at(i), trueNHitsZ->at(i));
            const double momentum = trueP->at(i);

            ++study.nTrue;
            if (matched) ++study.nMatched;
            study.all->Fill(matched, momentum);

            if (reconstructable) {
                ++study.nReconstructable;
                if (matched) ++study.nMatchedReconstructable;
            }
            // Fraction of ALL true particles satisfying the hit criterion.
            study.reconstructable->Fill(reconstructable, momentum);

            if (matched) {
                const int bestPfpIdx = trueBestRecoPfpIdx->at(i);
                if (bestPfpIdx < 0 ||
                    bestPfpIdx >= static_cast<int>(pfpParentID->size())) {
                    ++study.nInvalidBestPFP;
                } else {
                    const bool primary = (nuPfpSelfID >= 0 &&
                        pfpParentID->at(bestPfpIdx) == nuPfpSelfID);
                    if (primary) ++study.nPrimary;
                    else ++study.nSecondary;
                }
            }
        }
    }

    std::cout << "\n========== Reconstruction study ==========\n"
              << "Signal events: " << nSignalEvents << "\n"
              << "Events with != 1 selected nu root: " << nUnusualNuRoots << "\n"
              << "Events with inconsistent truth-vector sizes: " << nBadTruthVectors << "\n"
              << "Efficiencies below are integrated over all momenta;\n"
              << "plots use only the configured momentum ranges.\n\n";

    for (const auto& study : studies) {
        std::cout << study.name
                  << ": true=" << study.nTrue
                  << ", reconstructable=" << study.nReconstructable
                  << ", matched=" << study.nMatched
                  << ", matched & reconstructable=" << study.nMatchedReconstructable
                  << ", all efficiency=" << Fraction(study.nMatched, study.nTrue)
                  << ", reconstructable fraction="
                  << Fraction(study.nReconstructable, study.nTrue)
                  << ", matched given reconstructable="
                  << Fraction(study.nMatchedReconstructable, study.nReconstructable)
                  << ", primary=" << study.nPrimary
                  << ", secondary=" << study.nSecondary
                  << ", invalid best PFP=" << study.nInvalidBestPFP
                  << '\n';
    }

    for (const auto& study : studies) PlotEfficiencies(study);
    for (auto& study : studies) {
        delete study.all;
        delete study.reconstructable;
    }
    std::cout << "Plots saved in plots/efficiency{Muon,Photon,Proton,Pion}.png\n";
}
