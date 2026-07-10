// plotEffVsMomentum.C

#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"

void plotEffVsMomentum(
    const char* filename = "validTree_BDT.root",
    const char* treename = "validTree_BDT"
)
{
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: could not open " << filename << std::endl;
        return;
    }

    TTree *t = (TTree*)f->Get(treename);
    if (!t) {
        std::cerr << "ERROR: could not find tree " << treename << std::endl;
        f->ls();
        return;
    }

    std::vector<int>   *TruePDG = nullptr;
    std::vector<float> *TrueP   = nullptr;

    Bool_t IsSignal = false;
    Bool_t PassBDT  = false;

    t->SetBranchAddress("TruePDG", &TruePDG);
    t->SetBranchAddress("TrueP",   &TrueP);
    t->SetBranchAddress("IsSignal", &IsSignal);
    t->SetBranchAddress("PassBDT",  &PassBDT);

    gStyle->SetOptStat(0);

    // Adjust ranges/binning after looking at your truth distributions.
    const int nBins = 20;
    const double pMin = 0.0;
    const double pMax = 2.0; // GeV/c

    std::map<std::string, int> pdgs = {
        {"photon", 22},
        {"Lambda", 3122},
        {"proton", 2212},
        {"pion", 211}
    };

    std::map<std::string, TEfficiency*> effs;

    for (const auto &kv : pdgs) {
        const std::string &name = kv.first;

        effs[name] = new TEfficiency(
            Form("eff_%s", name.c_str()),
            Form("BDT efficiency vs true %s momentum;True p [GeV/c];Efficiency", name.c_str()),
            nBins, pMin, pMax
        );
    }

    Long64_t nEntries = t->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);

        // We want signal efficiency, so denominator is true signal events only.
        if (!IsSignal) continue;

        if (!TruePDG || !TrueP) continue;

        if (TruePDG->size() != TrueP->size()) {
            std::cerr << "WARNING: TruePDG and TrueP size mismatch in entry "
                      << i << ": " << TruePDG->size()
                      << " vs " << TrueP->size() << std::endl;
            continue;
        }

        bool foundPhoton = false;
        bool foundLambda = false;
        bool foundProton = false;
        bool foundPion   = false;

        float pPhoton = -999.f;
        float pLambda = -999.f;
        float pProton = -999.f;
        float pPion   = -999.f;

        for (size_t j = 0; j < TruePDG->size(); ++j) {
            int pdg = TruePDG->at(j);
            float p = TrueP->at(j);

            if (pdg == 22 && !foundPhoton) {
                pPhoton = p;
                foundPhoton = true;
            }

            if (pdg == 3122 && !foundLambda) {
                pLambda = p;
                foundLambda = true;
            }

            if (pdg == 2212 && !foundProton) {
                pProton = p;
                foundProton = true;
            }

            // Use abs if you want pi+ and pi- grouped together.
            if (std::abs(pdg) == 211 && !foundPion) {
                pPion = p;
                foundPion = true;
            }
        }

        if (foundPhoton) effs["photon"]->Fill(PassBDT, pPhoton);
        if (foundLambda) effs["Lambda"]->Fill(PassBDT, pLambda);
        if (foundProton) effs["proton"]->Fill(PassBDT, pProton);
        if (foundPion)   effs["pion"]->Fill(PassBDT, pPion);
    }

    TFile *out = TFile::Open("efficiency_vs_momentum.root", "RECREATE");

    for (auto &kv : effs) {
        const std::string &name = kv.first;
        TEfficiency *eff = kv.second;

        TCanvas *c = new TCanvas(
            Form("c_%s", name.c_str()),
            Form("Efficiency vs momentum: %s", name.c_str()),
            900, 700
        );

        eff->SetLineWidth(2);
        eff->SetMarkerStyle(20);
        eff->SetMarkerSize(1.0);

        eff->Draw("AP");
        c->Update();

        eff->Write();
        c->Write();
        c->SaveAs(Form("eff_vs_p_%s.png", name.c_str()));
    }

    out->Close();
    f->Close();

    std::cout << "Wrote efficiency_vs_momentum.root and PNG plots." << std::endl;
}
