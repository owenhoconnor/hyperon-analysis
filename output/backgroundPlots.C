#define backgroundPlots_cxx
#include "backgroundPlots.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>

// -----------------------------------------------------------------------------
// Define structures for particle classification
// -----------------------------------------------------------------------------

struct InteractionParticles {
    std::vector<int> pdg;
    std::vector<float> p;
};

struct ParticleSummary {
    int nProtons = 0;
    int nNeutrons = 0;
    int nPions = 0;
    int nPhotons = 0;
    int nStrange = 0;

    float leadingProtonP = -1.0;
    float leadingPionP = -1.0;
    float leadingPhotonP = -1.0;
};

// -----------------------------------------------------------------------------
// Interaction mode definitions
// -----------------------------------------------------------------------------

constexpr int kNModes = 14;

const std::array<std::string, kNModes> modeLabels = {
    "QE",
    "RES",
    "DIS",
    "Coh",
    "Coh Elastic",
    "e Scatter",
    "IMD Annih",
    "Inverse Beta",
    "Glashow Res",
    "AMNuGamma",
    "MEC",
    "Diffractive",
    "EM",
    "Weak Mix"
};

// -----------------------------------------------------------------------------
// Topology definitions
// -----------------------------------------------------------------------------

enum Topology {
    kCC0Pi0P = 0,
    kCC0Pi1P,
    kCC0Pi2PlusP,
    kCC1Pi0,
    kCC1PiPlus,
    kCC1PiMinus,
    kCCMultiPi,
    kCCStrange,

    kNC0Pi0P,
    kNC0Pi1P,
    kNC0Pi2PlusP,
    kNC1Pi0,
    kNC1PiPlus,
    kNC1PiMinus,
    kNCMultiPi,
    kNCStrange,

    kOther,

    kNTopologies
};

const std::array<std::string, kNTopologies> topologyLabels = {
    "CC 0#pi, 0p",
    "CC 0#pi, 1p",
    "CC 0#pi, #geq2p",
    "CC 1#pi^{0}",
    "CC 1#pi^{+}",
    "CC 1#pi^{-}",
    "CC multi-#pi",
    "CC strange",

    "NC 0#pi, 0p",
    "NC 0#pi, 1p",
    "NC 0#pi, #geq2p",
    "NC 1#pi^{0}",
    "NC 1#pi^{+}",
    "NC 1#pi^{-}",
    "NC multi-#pi",
    "NC strange",

    "Other"
};

// -----------------------------------------------------------------------------
// Identify strange-particle content
// -----------------------------------------------------------------------------

bool IsStrangeHadron(const int pdg){
    const int absPDG = std::abs(pdg);

    switch (absPDG) {
        // Kaons
        case 130:   // K_L
        case 310:   // K_S
        case 311:   // K0
        case 321:   // K+

        // Hyperons
        case 3122:  // Lambda
        case 3112:  // Sigma-
        case 3212:  // Sigma0
        case 3222:  // Sigma+
        case 3312:  // Xi-
        case 3322:  // Xi0
        case 3334:  // Omega-
            return true;

        default:
            return false;
    }
}

// -----------------------------------------------------------------------------
// Classify one interaction
// -----------------------------------------------------------------------------

int ClassifyTopology(
    const int ccnc,
    const std::vector<int>& truePDGs
)
{
    int nProtons = 0;
    int nNeutrons = 0;

    int nPiPlus = 0;
    int nPiMinus = 0;
    int nPiZero = 0;

    bool hasStrangeHadron = false;

    for (const int pdg : truePDGs) {

        switch (pdg) {
            case 2212:
                ++nProtons;
                break;

            case 2112:
                ++nNeutrons;
                break;

            case 211:
                ++nPiPlus;
                break;

            case -211:
                ++nPiMinus;
                break;

            case 111:
                ++nPiZero;
                break;

            default:
                break;
        }

        if (IsStrangeHadron(pdg)) {
            hasStrangeHadron = true;
        }
    }

    const int nChargedPions = nPiPlus + nPiMinus;
    const int nPions = nChargedPions + nPiZero;

    // In simb::MCNeutrino:
    //     CCNC == 0 -> charged current
    //     CCNC == 1 -> neutral current

    if (ccnc == 0) {

        // Give strange interactions their own category.
        if (hasStrangeHadron) {
            return kCCStrange;
        }

        if (nPions == 0) {
            if (nProtons == 0) {
                return kCC0Pi0P;
            }

            if (nProtons == 1) {
                return kCC0Pi1P;
            }

            return kCC0Pi2PlusP;
        }

        if (nPiZero == 1 && nChargedPions == 0) {
            return kCC1Pi0;
        }

        if (nPiPlus == 1 && nPiMinus == 0 && nPiZero == 0) {
            return kCC1PiPlus;
        }

        if (nPiMinus == 1 && nPiPlus == 0 && nPiZero == 0) {
            return kCC1PiMinus;
        }

        return kCCMultiPi;
    }

    if (ccnc == 1) {

        if (hasStrangeHadron) {
            return kNCStrange;
        }

        if (nPions == 0) {
            if (nProtons == 0) {
                return kNC0Pi0P;
            }

            if (nProtons == 1) {
                return kNC0Pi1P;
            }

            return kNC0Pi2PlusP;
        }

        if (nPiZero == 1 && nChargedPions == 0) {
            return kNC1Pi0;
        }

        if (nPiPlus == 1 && nPiMinus == 0 && nPiZero == 0) {
            return kNC1PiPlus;
        }

        if (nPiMinus == 1 && nPiPlus == 0 && nPiZero == 0) {
            return kNC1PiMinus;
        }

        return kNCMultiPi;
    }

    return kOther;
}

// Summarise particle content of an interaction

ParticleSummary SummariseParticles(
    const std::vector<int>& pdgs,
    const std::vector<float>& momenta
)
{
    ParticleSummary summary;

    const size_t n = std::min(pdgs.size(), momenta.size());

    for (size_t i = 0; i < n; ++i) {

        const int pdg = pdgs.at(i);
        const float p = momenta.at(i);

        if (pdg == 2212) {
            ++summary.nProtons;

            if (p > summary.leadingProtonP) {
                summary.leadingProtonP = p;
            }
        }

        if (pdg == 2112) {
            ++summary.nNeutrons;
        }

        if (
            pdg == 211 ||
            pdg == -211 ||
            pdg == 111
        ) {
            ++summary.nPions;

            if (p > summary.leadingPionP) {
                summary.leadingPionP = p;
            }
        }

        if (pdg == 22) {
            ++summary.nPhotons;

            if (p > summary.leadingPhotonP) {
                summary.leadingPhotonP = p;
            }
        }

        if (IsStrangeHadron(pdg)) {
            ++summary.nStrange;
        }
    }

    return summary;
}

// Get Old Interaction Particles

InteractionParticles GetOldInteractionParticles(
    const int truthIdx,
    const std::vector<int>& truePDGs,
    const std::vector<float>& truePs,
    const std::vector<int>& vertexSizes
)
{
    InteractionParticles particles;

    // Because vertexSize has an initial dummy zero.
    if (
        truthIdx < 0 ||
        static_cast<size_t>(truthIdx + 1) >= vertexSizes.size()
    ) {
        return particles;
    }

    size_t startIndex = 0;

    // Sum particle counts from all preceding MCTruths.
    for (int i = 0; i < truthIdx; ++i) {
        startIndex += vertexSizes.at(i + 1);
    }

    const size_t nParticles = vertexSizes.at(truthIdx + 1);

    const size_t endIndex =
        std::min(
            {
                startIndex + nParticles,
                truePDGs.size(),
                truePs.size()
            }
        );

    for ( size_t i = startIndex; i < endIndex; ++i){
        particles.pdg.push_back(truePDGs.at(i));
        particles.p.push_back(truePs.at(i));
    }

    return particles;
}

void backgroundPlots::Loop()
{
//   In a ROOT session, you can do:
//      root> .L backgroundPlots.C
//      root> backgroundPlots t
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

   TH1F *hTrueNuEnergy = new TH1F("hTrueNuEnergy", "", 30, 0.0, 3.0);
   TH1F *hTrueIntMode = new TH1F("hTrueIntMode", "", kNModes, -0.5, 13.5);
   TH1F *hTrueCCNC = new TH1F("hTrueCCNC", "", 2, -0.5, 1.5);

   // Mode vs neutrino energy

   std::array<TH1F*, kNModes> hNuEnergyByMode;

   for (int mode = 0; mode < kNModes; ++mode) {
      hNuEnergyByMode[mode] = new TH1F(Form("hNuEnergyMode_%d", mode), "", 100, 0.0, 5.0);
      hNuEnergyByMode[mode]->SetDirectory(nullptr);
   }

   // Topology vs neutrino energy

   std::array<TH1F*, kNTopologies> hNuEnergyByTopology;

   for (int topology = 0; topology < kNTopologies; ++topology) {
    hNuEnergyByTopology[topology] = new TH1F(Form("hNuEnergyTopology_%d", topology), "", 50, 0.0, 5.0);
    hNuEnergyByTopology[topology]->SetDirectory(nullptr);
   }

   // Mode vs neutrino energy 2D Matrix
   TH2F *hModeVsNuE = new TH2F("hModeVsNuE","", 30, 0.0, 3.0, kNModes, -0.5, 13.5);

   // Topology vs neutrino energy 2D Matrix
   TH2F *hTopologyVsNuE = new TH2F("hTopologyVsNuE", "", 30, 0.0, 3.0, kNTopologies, -0.5,kNTopologies - 0.5);

   // Interaction mode vs topology
   TH2I *hModeVsTopology = new TH2I("hModeVsTopology", "", kNTopologies, -0.5, kNTopologies - 0.5, kNModes, -0.5, 13.5);

   // Hadronic multiplicities
   TH1I *hNProtons = new TH1I("hNProtons","",6, -0.5, 5.5);

   TH1I *hNPions = new TH1I("hNPions", "", 6, -0.5, 5.5);

   TH1I *hNPhotons = new TH1I("hNPhotons", "", 6, -0.5, 5.5);

   // Leading-particle momenta
   TH1F *hLeadingProtonP = new TH1F("hLeadingProtonP", "", 30, 0.0, 2.0);

   TH1F *hLeadingPionP = new TH1F("hLeadingPionP", "", 30, 0.0, 2.0);

   TH1F *hLeadingPhotonP = new TH1F("hLeadingPhotonP", "", 30, 0.0, 1.0);

   for (int mode = 0; mode < kNModes; ++mode) {

    hTrueIntMode
        ->GetXaxis()
        ->SetBinLabel(
            mode + 1,
            modeLabels.at(mode).c_str()
        );

    hModeVsNuE
        ->GetYaxis()
        ->SetBinLabel(
            mode + 1,
            modeLabels.at(mode).c_str()
        );

    hModeVsTopology
        ->GetYaxis()
        ->SetBinLabel(
            mode + 1,
            modeLabels.at(mode).c_str()
        );
}

for (
    int topology = 0;
    topology < kNTopologies;
    ++topology
) {

    hTopologyVsNuE
        ->GetYaxis()
        ->SetBinLabel(
            topology + 1,
            topologyLabels.at(topology).c_str()
        );

    hModeVsTopology
        ->GetXaxis()
        ->SetBinLabel(
            topology + 1,
            topologyLabels.at(topology).c_str()
        );
}

hTrueCCNC->GetXaxis()->SetBinLabel(
    1,
    "CC"
);

hTrueCCNC->GetXaxis()->SetBinLabel(
    2,
    "NC"
);

   std::array<int, kNTopologies> topologyCounts{};

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // ---------------------------------------------------
      // Classify topologies and count them
      // ---------------------------------------------------

      const int truthIdx = chosenMCTruthIdx - 1;

      if (
         truthIdx < 0 ||
         static_cast<size_t>(truthIdx) >= trueIntMode->size() ||
         static_cast<size_t>(truthIdx) >= trueCCNC->size() ||
         static_cast<size_t>(truthIdx) >= trueNuEnergy->size()
      ) {
         continue;
        }

      const int intMode = trueIntMode->at(truthIdx);
      const int intType = trueIntType->at(truthIdx);
      const int ccnc = trueCCNC->at(truthIdx);
      const float nuEnergy = trueNuEnergy->at(truthIdx);

   // Extract ONLY the old primary particles belonging to the selected MCTruth.

      const InteractionParticles particles =
         GetOldInteractionParticles(
            truthIdx,
            *truePDG,
            *trueP,
            *vertexSize
         );

      const int topology =
         ClassifyTopology(
            ccnc,
            particles.pdg
         );

      const ParticleSummary summary =
         SummariseParticles(
            particles.pdg,
            particles.p
         );

      hTrueNuEnergy->Fill(nuEnergy);
      hTrueIntMode->Fill(intMode);
      hTrueCCNC->Fill(ccnc);

      if (intMode >= 0 && intMode < kNModes){

      hNuEnergyByMode[intMode]->Fill(nuEnergy);
      hNuEnergyByTopology[topology]->Fill(nuEnergy);

      hModeVsNuE->Fill(
         nuEnergy,
         intMode
      );

      hModeVsTopology->Fill(
         topology,
         intMode
      );
      }

      hTopologyVsNuE->Fill(nuEnergy, topology);

      hNProtons->Fill(summary.nProtons);

      hNPions->Fill(summary.nPions);

      hNPhotons->Fill(summary.nPhotons);

      if (summary.leadingProtonP >= 0.0) {
            hLeadingProtonP->Fill(summary.leadingProtonP);
         }

         if (summary.leadingPionP >= 0.0) {
            hLeadingPionP->Fill(summary.leadingPionP);
         }

         if (summary.leadingPhotonP >= 0.0) {
            hLeadingPhotonP->Fill(summary.leadingPhotonP);
         }

         ++topologyCounts.at(topology);

   } // End event loop

   // --------------------------------------------------------------------------
   // Prepare Histograms
   // --------------------------------------------------------------------------

   THStack *hsModeVsNuE = new THStack("hsModeVsNuE", "Surviving Background by Interaction Mode;True E_{#nu} [GeV];Interactions");
   THStack *hsTopologyVsNuE = new THStack("hsTopologyVsNuE", "Surviving Background by Interaction Topology;True E_{#nu} [GeV];Interactions");

   const std::array<int, kNModes> modeColours = {
    kAzure + 1,
    kOrange + 7,
    kGreen + 2,
    kRed - 4,
    kViolet - 4,
    kCyan + 1,
    kPink + 1,
    kYellow + 1,
    kTeal + 2,
    kMagenta - 2,
    kBlue - 7,
    kOrange - 3,
    kSpring + 5,
    kGray + 1
   };

   const std::array<int, kNTopologies> topologyColours = {
    kAzure + 1,
    kAzure + 6,
    kBlue - 7,

    kOrange + 7,
    kOrange + 1,
    kRed - 4,
    kPink + 1,
    kViolet - 4,

    kGreen + 2,
    kGreen - 6,
    kSpring + 5,

    kCyan + 1,
    kTeal + 2,
    kYellow + 1,
    kMagenta - 2,
    kViolet + 1,

    kGray + 1
   };

   TLegend *legMode = new TLegend(0.72, 0.55, 0.89, 0.89);
   legMode->SetBorderSize(0);
   legMode->SetFillStyle(0);

for (int mode = 0; mode < kNModes; ++mode) {

    // Don't add empty modes to the plot/legend.
    if (hNuEnergyByMode[mode]->GetEntries() == 0) {
        continue;
    }

    hNuEnergyByMode[mode]->SetFillColorAlpha(
        modeColours[mode],
        0.75
    );

    hNuEnergyByMode[mode]->SetLineColor(
        kBlack
    );

    hNuEnergyByMode[mode]->SetLineWidth(1);

    hsModeVsNuE->Add(
        hNuEnergyByMode[mode]
    );

    legMode->AddEntry(
        hNuEnergyByMode[mode],
        modeLabels[mode].c_str(),
        "f"
    );
}

TCanvas *cModeVsNuE = new TCanvas(
    "cModeVsNuE",
    "Interaction Mode vs Neutrino Energy",
    1800,
    1200
);

cModeVsNuE->SetRightMargin(0.22);

gStyle->SetOptStat(0);

hsModeVsNuE->Draw("HIST");

hsModeVsNuE->GetXaxis()->SetTitle(
    "True E_{#nu} [GeV]"
);

hsModeVsNuE->GetYaxis()->SetTitle(
    "Interactions"
);

hsModeVsNuE->GetXaxis()->CenterTitle();
hsModeVsNuE->GetYaxis()->CenterTitle();

legMode->Draw();

cModeVsNuE->Print(
    "plots/modeVsNuEnergyBkg.png"
);

cModeVsNuE->SaveAs(
    "plots/modeVsNuEnergyBkg.C"
);

TLegend *legTopology = new TLegend(
    0.73,
    0.40,
    0.93,
    0.90
);

legTopology->SetBorderSize(0);
legTopology->SetFillStyle(0);
legTopology->SetTextSize(0.025);

for (
    int topology = 0;
    topology < kNTopologies;
    ++topology
) {

    if (
        hNuEnergyByTopology[topology]
            ->GetEntries() == 0
    ) {
        continue;
    }

    hNuEnergyByTopology[topology]
        ->SetFillColorAlpha(
            topologyColours[topology],
            0.75
        );

    hNuEnergyByTopology[topology]
        ->SetLineColor(
            kBlack
        );

    hsTopologyVsNuE->Add(
        hNuEnergyByTopology[topology]
    );

    legTopology->AddEntry(
        hNuEnergyByTopology[topology],
        topologyLabels[topology].c_str(),
        "f"
    );
}

TCanvas *cTopologyVsNuE = new TCanvas(
    "cTopologyVsNuE",
    "Topology vs Neutrino Energy",
    1800,
    1200
);

cTopologyVsNuE->SetRightMargin(0.27);

gStyle->SetOptStat(0);

hsTopologyVsNuE->Draw("HIST");

hsTopologyVsNuE->GetXaxis()->SetTitle(
    "True E_{#nu} [GeV]"
);

hsTopologyVsNuE->GetYaxis()->SetTitle(
    "Interactions"
);

hsTopologyVsNuE->GetXaxis()->CenterTitle();
hsTopologyVsNuE->GetYaxis()->CenterTitle();

legTopology->Draw();

cTopologyVsNuE->Print(
    "plots/topologyVsNuEnergyBkg.png"
);

cTopologyVsNuE->SaveAs(
    "plots/topologyVsNuEnergyBkg.C"
);

   std::vector<int> nonEmptyBins;

   for (int bin = 1;
      bin <= hTrueIntMode->GetNbinsX();
      ++bin) {

      if (hTrueIntMode->GetBinContent(bin) != 0) {
         nonEmptyBins.push_back(bin);
      }

   } 

   TH1F* hTrueIntModeCompact = new TH1F(
      "hTrueIntModeCompact",
      "",
      static_cast<int>(nonEmptyBins.size()),
      0.0,
      static_cast<double>(nonEmptyBins.size())
   );

   for (std::size_t i = 0; i < nonEmptyBins.size(); ++i) {

      const int oldBin = nonEmptyBins[i];
      const int newBin = static_cast<int>(i) + 1;

      hTrueIntModeCompact->SetBinContent(
         newBin,
         hTrueIntMode->GetBinContent(oldBin)
      );

      hTrueIntModeCompact->SetBinError(
         newBin,
         hTrueIntMode->GetBinError(oldBin)
      );

      const char* oldLabel =
         hTrueIntMode->GetXaxis()->GetBinLabel(oldBin);

      hTrueIntModeCompact->GetXaxis()->SetBinLabel(
         newBin,
         oldLabel
      );
}

   /*TCanvas *c1 = new TCanvas("c1", "True Interaction Type", 2000, 2000);
   hTrueIntType->Draw("HIST");
   c1->Print("plots/trueBeamIntType.png");
   c1->SaveAs("plots/trueBeamIntType.C");*/

   TCanvas *c2 = new TCanvas("c2", "True Interaction Mode", 2000, 2000);

   //hTrueIntMode->GetXAxis()->SetTitle("Interaction Mode");
   hTrueIntMode->GetXaxis()->SetBinLabel(1, "QE");
   hTrueIntMode->GetXaxis()->SetBinLabel(2, "RES");
   hTrueIntMode->GetXaxis()->SetBinLabel(3, "DIS");
   hTrueIntMode->GetXaxis()->SetBinLabel(4, "Coh");
   hTrueIntMode->GetXaxis()->SetBinLabel(5, "Coh Elastic");
   hTrueIntMode->GetXaxis()->SetBinLabel(6, "e Scatter");
   hTrueIntMode->GetXaxis()->SetBinLabel(7, "IMD Annih");
   hTrueIntMode->GetXaxis()->SetBinLabel(8, "Inverse Beta");
   hTrueIntMode->GetXaxis()->SetBinLabel(9, "Glashow Res");
   hTrueIntMode->GetXaxis()->SetBinLabel(10, "AMNuGamma");
   hTrueIntMode->SetFillColorAlpha(kAzure+1, 0.7);
   hTrueIntMode->SetLineColor(kBlack);
   hTrueIntMode->SetLineWidth(2);
   hTrueIntMode->Draw("HIST");
   c2->Print("plots/trueBeamIntMode.png");
   c2->SaveAs("plots/trueBeamIntMode.C");

   TCanvas *c3 = new TCanvas("c3", "True CCNC", 2000, 2000);

  // hTrueCCNC->GetXAxis()->SetTitle("CCNC");
   hTrueCCNC->GetXaxis()->SetBinLabel(1, "CC");
   hTrueCCNC->GetXaxis()->SetBinLabel(2, "NC");

   hTrueCCNC->GetXaxis()->SetTitle("Interaction Current");
   hTrueCCNC->GetYaxis()->SetTitle("Interactions");

   hTrueCCNC->GetXaxis()->CenterTitle();
   hTrueCCNC->GetYaxis()->CenterTitle();

   hTrueCCNC->SetMinimum(0);

   gStyle->SetOptStat(0);
   hTrueCCNC->SetFillColorAlpha(kAzure + 1, 0.7);
   hTrueCCNC->SetLineColor(kBlack);
   hTrueCCNC->SetLineWidth(2);

   hTrueCCNC->Draw("HIST TEXT0");
   c3->Print("plots/trueBeamCCNC.png");
   c3->SaveAs("plots/trueBeamCCNC.C");

   TCanvas *c4 = new TCanvas("c4", "True Interaction Mode Compact", 2000, 2000);
   hTrueIntModeCompact->LabelsOption("v", "X");
   hTrueIntModeCompact->Draw("HIST TEXT0");
   c4->Print("plots/trueIntModeCompact.png");
   c4->SaveAs("plots/trueIntModeCompact.C");

   std::vector<int> populatedTopologies;

    for (int topology = 0;
        topology < kNTopologies;
        ++topology) {

        if (topologyCounts.at(topology) > 0) {
            populatedTopologies.push_back(topology);
        }
    }   

    TH1I* hTopology = new TH1I(
        "hTopology",
        "",
        static_cast<int>(populatedTopologies.size()),
        0.0,
        static_cast<double>(populatedTopologies.size())
    );

    for (std::size_t i = 0;
        i < populatedTopologies.size();
        ++i) {

        const int topology = populatedTopologies.at(i);
        const int bin = static_cast<int>(i) + 1;

        hTopology->SetBinContent(
            bin,
            topologyCounts.at(topology)
        );

        hTopology->GetXaxis()->SetBinLabel(
            bin,
            topologyLabels.at(topology).c_str()
        );
    }

        gStyle->SetOptStat(0);

        TCanvas* cTopology = new TCanvas(
            "cTopology",
            "Final-state topologies",
            1800,
            1100
        );

        cTopology->SetLeftMargin(0.12);
        cTopology->SetRightMargin(0.05);
        cTopology->SetBottomMargin(0.30);

        hTopology->SetTitle(
            "Background Final State Topology Distribution After Presel"
        );

        hTopology->GetXaxis()->SetTitle("True final state topology");
        hTopology->GetYaxis()->SetTitle("Interactions");

        hTopology->GetXaxis()->CenterTitle();
        hTopology->GetYaxis()->CenterTitle();

        hTopology->GetXaxis()->SetTitleOffset(3.5);
        hTopology->GetYaxis()->SetTitleOffset(1.3);

        hTopology->GetXaxis()->SetLabelSize(0.035);
        hTopology->GetYaxis()->SetLabelSize(0.035);

        hTopology->SetMinimum(0.0);
        hTopology->SetMaximum(
            1.20 * hTopology->GetMaximum()
        );

        // Vertical labels are usually clearest for long topology names.
        hTopology->LabelsOption("v", "X");

        hTopology->SetFillColorAlpha(kAzure +1, 0.7);
        hTopology->SetLineColor(kBlack);
        hTopology->SetLineWidth(2);
        hTopology->Draw("HIST TEXT0");

        cTopology->Modified();
        cTopology->Update();

        cTopology->Print("plots/finalStateTopologiesBkg.png");
        cTopology->SaveAs("plots/finalStateTopologiesBkg.C");
   

   TCanvas *cModeVsTopology =
    new TCanvas(
        "cModeVsTopology",
        "Mode vs topology",
        1900,
        1300
    );

cModeVsTopology->SetLeftMargin(0.15);
cModeVsTopology->SetRightMargin(0.15);
cModeVsTopology->SetBottomMargin(0.30);

hModeVsTopology->SetTitle(
    "Surviving Background: Interaction Mode vs Final-State Topology"
);

hModeVsTopology
    ->GetXaxis()
    ->SetTitle(
        "Final-State Topology"
    );

hModeVsTopology
    ->GetYaxis()
    ->SetTitle(
        "Interaction Mode"
    );

hModeVsTopology
    ->GetXaxis()
    ->LabelsOption(
        "v"
    );

hModeVsTopology->Draw(
    "COLZ TEXT"
);

cModeVsTopology->Print(
    "plots/modeVsTopologyBkg.png"
);

TCanvas *cTopologyVsNuE =
    new TCanvas(
        "cTopologyVsNuE",
        "Topology vs neutrino energy",
        1800,
        1300
    );

cTopologyVsNuE->SetLeftMargin(0.18);
cTopologyVsNuE->SetRightMargin(0.15);

hTopologyVsNuE->SetTitle(
    "Surviving Background: Final-State Topology vs True Neutrino Energy"
);

hTopologyVsNuE
    ->GetXaxis()
    ->SetTitle(
        "True E_{#nu} [GeV]"
    );

hTopologyVsNuE
    ->GetYaxis()
    ->SetTitle(
        "Final-State Topology"
    );

hTopologyVsNuE->Draw(
    "COLZ"
);

cTopologyVsNuE->Print(
    "plots/topologyVsNuEnergyBkg.png"
);
}
