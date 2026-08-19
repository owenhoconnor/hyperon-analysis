#define newBackgroundPlots_cxx

#include "newBackgroundPlots.h"

#include <TH1F.h>
#include <TH1I.h>
#include <TH2I.h>
#include <THStack.h>

#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>


// ============================================================================
// INTERACTION MODE DEFINITIONS
// ============================================================================

// MCNeutrino::Mode() can run from -1 (unknown) to 13.
//
// Array index:
//   0 -> mode -1
//   1 -> mode  0
//   ...
//  14 -> mode 13

constexpr int kMinMode = -1;
constexpr int kMaxMode = 13;
constexpr int kNModeCategories = kMaxMode - kMinMode + 1;

const std::array<
    std::string,
    kNModeCategories
> modeLabels = {

    "Unknown",       // -1
    "QE",            //  0
    "RES",           //  1
    "DIS",           //  2
    "Coh",           //  3
    "Coh Elastic",   //  4
    "e Scatter",     //  5
    "IMD Annih",     //  6
    "Inverse Beta",  //  7
    "Glashow Res",   //  8
    "AMNuGamma",     //  9
    "MEC",           // 10
    "Diffractive",   // 11
    "EM",            // 12
    "Weak Mix"       // 13
};


int ModeToArrayIndex(const int mode)
{
    if (
        mode < kMinMode ||
        mode > kMaxMode
    ) {
        return -1;
    }

    return mode - kMinMode;
}


// ============================================================================
// TOPOLOGY DEFINITIONS
// ============================================================================

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


const std::array<
    std::string,
    kNTopologies
> topologyLabels = {

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


// ============================================================================
// STRANGE-PARTICLE IDENTIFICATION
// ============================================================================

bool IsStrangeHadron(const int pdg)
{
    const int absPDG = std::abs(pdg);

    switch (absPDG) {

        // Kaons
        case 130:
        case 310:
        case 311:
        case 321:

        // Hyperons
        case 3122:
        case 3112:
        case 3212:
        case 3222:
        case 3312:
        case 3322:
        case 3334:

            return true;

        default:

            return false;
    }
}


// ============================================================================
// CLASSIFY FINAL-STATE TOPOLOGY
// ============================================================================

int ClassifyTopology(
    const int ccnc,
    const std::vector<int>& truePDGs
)
{
    int nProtons = 0;
    int nPiPlus = 0;
    int nPiMinus = 0;
    int nPiZero = 0;

    bool hasStrangeHadron = false;


    for (const int pdg : truePDGs)
    {
        switch (pdg) {

            case 2212:

                ++nProtons;
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

    // ------------------------------------------------------------------------
    // CC
    // ------------------------------------------------------------------------

    if (ccnc == 0)
    {
        if (hasStrangeHadron) {
            return kCCStrange;
        }

        if (nPions == 0)
        {
            if (nProtons == 0) {
                return kCC0Pi0P;
            }

            if (nProtons == 1) {
                return kCC0Pi1P;
            }

            return kCC0Pi2PlusP;
        }

        if (
            nPiZero == 1 &&
            nChargedPions == 0
        ) {
            return kCC1Pi0;
        }


        if (
            nPiPlus == 1 &&
            nPiMinus == 0 &&
            nPiZero == 0
        ) {
            return kCC1PiPlus;
        }

        if (
            nPiMinus == 1 &&
            nPiPlus == 0 &&
            nPiZero == 0
        ) {
            return kCC1PiMinus;
        }

        return kCCMultiPi;
    }


    // ------------------------------------------------------------------------
    // NC
    // ------------------------------------------------------------------------

    if (ccnc == 1)
    {
        if (hasStrangeHadron) {
            return kNCStrange;
        }

        if (nPions == 0)
        {
            if (nProtons == 0) {
                return kNC0Pi0P;
            }

            if (nProtons == 1) {
                return kNC0Pi1P;
            }

            return kNC0Pi2PlusP;
        }


        if (
            nPiZero == 1 &&
            nChargedPions == 0
        ) {
            return kNC1Pi0;
        }

        if (
            nPiPlus == 1 &&
            nPiMinus == 0 &&
            nPiZero == 0
        ) {
            return kNC1PiPlus;
        }

        if (
            nPiMinus == 1 &&
            nPiPlus == 0 &&
            nPiZero == 0
        ) {
            return kNC1PiMinus;
        }


        return kNCMultiPi;
    }

    return kOther;
}

// ============================================================================
// PRIMARY PARTICLES BELONGING TO ONE MCTRUTH
// ============================================================================

struct InteractionParticles {
    std::vector<int> pdg;
    std::vector<float> p;
};


InteractionParticles GetPrimaryInteractionParticles(
    const int truthIdx,
    const std::vector<int>& pdgs,
    const std::vector<float>& momenta,

    const std::vector<int>& mcTruthIndices,
    const std::vector<int>& generations
)
{
    InteractionParticles particles;

    const size_t n =
        std::min(
            {
                pdgs.size(),
                momenta.size(),
                mcTruthIndices.size(),
                generations.size()
            }
        );

    for (size_t i = 0; i < n; ++i)
    {
        // Particle must belong to our selected neutrino interaction.

        if (mcTruthIndices.at(i) != truthIdx){
            continue;
        }

        // generation == 0 means interaction primary.
        // Decay products NOT included when defining final-state topology.

        if (generations.at(i) != 0) {
            continue;
        }

        particles.pdg.push_back(
            pdgs.at(i)
        );

        particles.p.push_back(
            momenta.at(i)
        );
    }

    return particles;
}


// ============================================================================
// SIMPLE PARTICLE CONTENT SUMMARY
// ============================================================================

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


ParticleSummary SummariseParticles(
    const std::vector<int>& pdgs,
    const std::vector<float>& momenta
)
{
    ParticleSummary summary;

    const size_t n = std::min(pdgs.size(), momenta.size());

    for (size_t i = 0; i < n; ++i)
    {
        const int pdg = pdgs.at(i);
        const float p = momenta.at(i);

        // --------------------------------------------------------------------
        // Protons
        // --------------------------------------------------------------------

        if (pdg == 2212)
        {
            ++summary.nProtons;

            if (p > summary.leadingProtonP) {
                summary.leadingProtonP = p;
            }
        }

        // --------------------------------------------------------------------
        // Neutrons
        // --------------------------------------------------------------------

        if (pdg == 2112) {
            ++summary.nNeutrons;
        }

        // --------------------------------------------------------------------
        // Pions
        // --------------------------------------------------------------------

        if (pdg == 211 || pdg == -211 || pdg == 111)
        {
            ++summary.nPions;

            if (
                p > summary.leadingPionP
            ) {
                summary.leadingPionP = p;
            }
        }

        // --------------------------------------------------------------------
        // Photons
        // --------------------------------------------------------------------

        if (pdg == 22)
        {
            ++summary.nPhotons;

            if (
                p > summary.leadingPhotonP
            ) {
                summary.leadingPhotonP = p;
            }
        }

        // --------------------------------------------------------------------
        // Strange particles
        // --------------------------------------------------------------------

        if (IsStrangeHadron(pdg)) {
            ++summary.nStrange;
        }
    }

    return summary;
}


// ============================================================================
// MAIN LOOP
// ============================================================================

void newBackgroundPlots::Loop()
{
    if (fChain == 0) {return;}

    gSystem->mkdir("plots", true);
    gStyle->SetOptStat(0);

    // ========================================================================
    // HISTOGRAMS
    // ========================================================================

    // ------------------------------------------------------------------------
    // CC / NC
    // ------------------------------------------------------------------------

    TH1F* hTrueCCNC = new TH1F("hTrueCCNC", "", 2, -0.5, 1.5);
    hTrueCCNC->SetDirectory(nullptr);
    hTrueCCNC->GetXaxis()->SetBinLabel(1, "CC");
    hTrueCCNC->GetXaxis()->SetBinLabel(2, "NC");

    // ------------------------------------------------------------------------
    // Interaction mode
    // ------------------------------------------------------------------------

    TH1F* hTrueIntMode = new TH1F("hTrueIntMode", "", kNModeCategories, -1.5, 13.5);
    hTrueIntMode->SetDirectory(nullptr);

    for (int i = 0; i < kNModeCategories; ++i){
        hTrueIntMode->GetXaxis()->SetBinLabel(i + 1, modeLabels.at(i).c_str());
    }

    // ------------------------------------------------------------------------
    // MCTruth orgin
    // ------------------------------------------------------------------------

    std::map<int, int> originCounts;

    TH1F* hTrueOrigin = new TH1F("hTrueOrigin", "", 3, -0.5, 2.5);
    hTrueOrigin->SetDirectory(nullptr)
    hTrueOrigin->GetXaxis()->SetBinLabel(1, "Beam Neutrino")
    hTrueOrigin->GetXAxis()->SetBinLabel(2, "Cosmic")

    // ------------------------------------------------------------------------
    // Overall true neutrino energy
    // ------------------------------------------------------------------------

    TH1F* hTrueNuEnergy = new TH1F("hTrueNuEnergy", "", 30, 0.0, 3.0);
    hTrueNuEnergy->SetDirectory(nullptr);

    // ------------------------------------------------------------------------
    // Interaction type
    //
    // InteractionType codes are sparse, so collect counts first and make
    // a compact categorical histogram after the event loop.
    // ------------------------------------------------------------------------

    std::map<int, int> interactionTypeCounts;

    // ------------------------------------------------------------------------
    // Final-state topology counts
    // ------------------------------------------------------------------------

    std::array<
        int,
        kNTopologies
    > topologyCounts{};

    // ------------------------------------------------------------------------
    // Nu E distributions separated by interaction mode
    // ------------------------------------------------------------------------

    std::array<
        TH1F*,
        kNModeCategories
    > hNuEnergyByMode;

    for (int modeIndex = 0; modeIndex < kNModeCategories; ++modeIndex){
        hNuEnergyByMode.at(modeIndex) = new TH1F(Form("hNuEnergyMode_%d", modeIndex), "", 30, 0.0, 3.0);
        hNuEnergyByMode.at(modeIndex)->SetDirectory(nullptr);
    }

    // ------------------------------------------------------------------------
    // E_nu distributions separated by topology
    // ------------------------------------------------------------------------

    std::array<
        TH1F*,
        kNTopologies
    > hNuEnergyByTopology;

    for (int topology = 0; topology < kNTopologies; ++topology) {
        hNuEnergyByTopology.at(topology) = new TH1F(Form("hNuEnergyTopology_%d", topology), "", 30, 0.0, 3.0);
        hNuEnergyByTopology.at(topology)->SetDirectory(nullptr);
    }

    // ------------------------------------------------------------------------
    // Mode vs topology
    // ------------------------------------------------------------------------

    TH2I* hModeVsTopology = new TH2I("hModeVsTopology", "", kNTopologies, -0.5, kNTopologies - 0.5, kNModeCategories, -1.5, 13.5);
    hModeVsTopology->SetDirectory(nullptr);

    for (int topology = 0; topology < kNTopologies; ++topology) {
        hModeVsTopology->GetXaxis()->SetBinLabel(topology + 1, topologyLabels.at(topology).c_str());
    }

    for (int modeIndex = 0; modeIndex < kNModeCategories; ++modeIndex) {
        hModeVsTopology->GetYaxis()->SetBinLabel(modeIndex + 1, modeLabels.at(modeIndex).c_str());
    }

    // ------------------------------------------------------------------------
    // Final-state multiplicities
    // ------------------------------------------------------------------------

    TH1I* hNProtons = new TH1I("hNProtons", "", 7, -0.5, 6.5);
    TH1I* hNNeutrons = new TH1I("hNNeutrons", "", 7, -0.5, 6.5);
    TH1I* hNPions = new TH1I("hNPions", "", 7, -0.5, 6.5);
    TH1I* hNPhotons =new TH1I("hNPhotons","",7,-0.5,6.5);
    TH1I* hNStrange =new TH1I("hNStrange", "", 5, -0.5, 4.5);

    // ------------------------------------------------------------------------
    // Leading-particle momenta
    // ------------------------------------------------------------------------

    TH1F* hLeadingProtonP = new TH1F("hLeadingProtonP", "", 30, 0.0, 2.0);
    TH1F* hLeadingPionP =new TH1F("hLeadingPionP","", 30, 0.0, 2.0);
    TH1F* hLeadingPhotonP = new TH1F("hLeadingPhotonP", "", 30, 0.0, 1.0);

    // ========================================================================
    // EVENT LOOP
    // ========================================================================

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0;
    Long64_t nb = 0;

    int nBadTruthIndices = 0;

    for (Long64_t jentry = 0;jentry < nentries;++jentry)
    {
        Long64_t ientry = LoadTree(jentry);

        if (ientry < 0) {
            break;
        }

        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        // --------------------------------------------------------------------
        // Sanity checks
        // --------------------------------------------------------------------

        if (
            chosenTruthIdx < 0 ||

            static_cast<size_t>(chosenTruthIdx)
                >= trueIntMode->size() ||

            static_cast<size_t>(chosenTruthIdx)
                >= trueIntType->size() ||

            static_cast<size_t>(chosenTruthIdx)
                >= trueCCNC->size() ||

            static_cast<size_t>(chosenTruthIdx)
                >= trueNuEnergy->size()
        )
        {
            ++nBadTruthIndices;
            std::cout<< "WARNING: invalid chosen MCTruth index "<< chosenTruthIdx<< " in tree entry "<< jentry<< std::endl;
            continue;
        }

        // --------------------------------------------------------------------
        // Interaction-level quantities
        // --------------------------------------------------------------------

        const int intMode = trueIntMode->at(chosenTruthIdx);
        const int intType = trueIntType->at(chosenTruthIdx);
        const int ccnc = trueCCNC->at(chosenTruthIdx);
        const int origin = trueOrigin->at(chosenTruthIdx);
        const float nuEnergy = trueNuEnergy->at(chosenTruthIdx);

        // ====================================================================
        // SELECT PRIMARY PARTICLES BELONGING TO THIS MCTRUTH
        // ====================================================================

        const InteractionParticles particles =
        GetPrimaryInteractionParticles(
                chosenTruthIdx,
                *truePDG,
                *trueP,
                *trueMCTruthIndex,
                *trueGeneration
            );

        // ====================================================================
        // CLASSIFY EVENT
        // ====================================================================

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


        // ====================================================================
        // FILL HISTOGRAMS
        // ====================================================================

        hTrueCCNC->Fill(ccnc);
        hTrueIntMode->Fill(intMode);
        hTrueOrigin->Fill(origin);
        hTrueNuEnergy->Fill(nuEnergy);

        ++originCounts[origin];
        ++interactionTypeCounts[intType];
        ++topologyCounts.at(topology);

        // --------------------------------------------------------------------
        // Mode vs E_nu stack
        // --------------------------------------------------------------------

        const int modeIndex = ModeToArrayIndex(intMode);

        if (modeIndex >= 0)
        {
            hNuEnergyByMode.at(modeIndex)->Fill(nuEnergy);
            hModeVsTopology->Fill(topology, intMode);
        }

        // --------------------------------------------------------------------
        // Topology vs E_nu stack
        // --------------------------------------------------------------------

        hNuEnergyByTopology.at(topology)->Fill(nuEnergy);

        // --------------------------------------------------------------------
        // Multiplicity
        // --------------------------------------------------------------------

        hNProtons->Fill(summary.nProtons);
        hNNeutrons->Fill(summary.nNeutrons);
        hNPions->Fill(summary.nPions);
        hNPhotons->Fill(summary.nPhotons);
        hNStrange->Fill(summary.nStrange);

        // --------------------------------------------------------------------
        // Leading momenta
        // --------------------------------------------------------------------

        if (
            summary.leadingProtonP >= 0.0
        ) {
            hLeadingProtonP->Fill(
                summary.leadingProtonP
            );
        }


        if (
            summary.leadingPionP >= 0.0
        ) {
            hLeadingPionP->Fill(
                summary.leadingPionP
            );
        }


        if (
            summary.leadingPhotonP >= 0.0
        ) {
            hLeadingPhotonP->Fill(
                summary.leadingPhotonP
            );
        }
    }

    std::cout<< "Finished loop over "<< nentries<< " entries."<< std::endl;
    std::cout<< "Entries skipped because of bad chosenMCchosenTruthIdx: "<< nBadTruthIndices<< std::endl;
    std::cout<<"# of Events with Beam Neutrino Origin = "<<originCounts[0]<<std::endl;

    // ========================================================================
    // COMMON HISTOGRAM STYLE
    // ========================================================================

    const int fillColour = kAzure + 1;

    // ========================================================================
    // 1. CC / NC
    // ========================================================================

    TCanvas* cCCNC = new TCanvas("cCCNC", "True CCNC", 1600, 1200);

    hTrueCCNC->GetXaxis()->SetTitle("Interaction Current");
    hTrueCCNC->GetYaxis()->SetTitle("Interactions");
    hTrueCCNC->SetMinimum(0.0);
    hTrueCCNC->SetFillColorAlpha(fillColour, 0.7);
    hTrueCCNC->SetLineColor(kBlack);
    hTrueCCNC->SetLineWidth(2);
    hTrueCCNC->Draw("HIST TEXT0");
    cCCNC->Print("plots/new_logic/trueBeamCCNC_new.png");

    // ========================================================================
    // 2. INTERACTION MODE DISTRIBUTION
    // ========================================================================

    TCanvas* cMode = new TCanvas("cMode", "True Interaction Mode", 1800, 1200);

    cMode->SetBottomMargin(0.20);

    hTrueIntMode->GetXaxis()->SetTitle("Interaction Mode");
    hTrueIntMode->GetYaxis()->SetTitle("Interactions");
    hTrueIntMode->SetFillColorAlpha(fillColour,0.7);
    hTrueIntMode->SetLineColor(kBlack);
    hTrueIntMode->SetLineWidth(2);
    hTrueIntMode->LabelsOption("v","X");
    hTrueIntMode->Draw("HIST TEXT0");
    cMode->Print("plots/new_logic/trueBeamIntMode_new.png");

    // ========================================================================
    // 3. INTERACTION TYPE DISTRIBUTION
    // ========================================================================

    TH1I* hTrueIntType = new TH1I("hTrueIntType", "", static_cast<int>(interactionTypeCounts.size()), 0.0, static_cast<double>(interactionTypeCounts.size()));
    int intTypeBin = 1;

    for (const auto& entry : interactionTypeCounts)
    {
        hTrueIntType->SetBinContent(intTypeBin, entry.second);
        hTrueIntType->GetXaxis()->SetBinLabel(intTypeBin, std::to_string(entry.first).c_str());
        ++intTypeBin;
    }

    TCanvas* cIntType = new TCanvas("cIntType", "Interaction Type", 1800, 1200);
    cIntType->SetBottomMargin(0.20);

    hTrueIntType->GetXaxis()->SetTitle("Interaction Type Code");
    hTrueIntType->GetYaxis()->SetTitle("Interactions");
    hTrueIntType->SetFillColorAlpha(fillColour, 0.7);
    hTrueIntType->SetLineColor(kBlack);
    hTrueIntType->LabelsOption("v","X");
    hTrueIntType->Draw("HIST TEXT0");
    cIntType->Print("plots/new_logic/trueBeamIntType_new.png");

    // ========================================================================
    // 4. MCTRUTH ORIGIN
    // ========================================================================

    TCanvas* cNuOrigin = new TCanvas("cNuOrigin", "True Neutrino Origin", 1600, 1200);
    hTrueOrigin->GetXAxis()->SetTitle("True Neutrino Origin");
    hTrueOrigin->GetYAxis()->SetTitle("Interactions");
    hTrueOrigin->SetFillColorAlpha(fillColour, 0.7);
    hTrueOrigin->SetLineColor(kBlack);
    hTrueOrigin->LabelsOption("v", "X");
    hTrueOrigin->Draw("HIST TEXT0");
    cNuOrigin->Print("plots/trueNuOrigin.pdf")

    // ========================================================================
    // 5. OVERALL NEUTRINO ENERGY
    // ========================================================================

    TCanvas* cNuEnergy = new TCanvas("cNuEnergy", "True Neutrino Energy", 1600, 1200);

    hTrueNuEnergy->GetXaxis()->SetTitle("True E_{#nu} [GeV]");
    hTrueNuEnergy->GetYaxis()->SetTitle("Interactions");
    hTrueNuEnergy->SetFillColorAlpha(fillColour,0.7);
    hTrueNuEnergy->SetLineColor(kBlack);
    hTrueNuEnergy->Draw("HIST");
    cNuEnergy->Print("plots/new_logic/trueNuEnergy_new.png");

    // ========================================================================
    // 6. FINAL-STATE TOPOLOGY DISTRIBUTION
    // ========================================================================

    std::vector<int> populatedTopologies;

    for (int topology = 0; topology < kNTopologies; ++topology) {

        if (topologyCounts.at(topology) > 0) {
            populatedTopologies.push_back(topology);
        }
    }

    TH1I* hTopology = new TH1I("hTopology","", static_cast<int>(populatedTopologies.size()), 0.0, static_cast<double>(populatedTopologies.size()));

    for (size_t i = 0; i < populatedTopologies.size(); ++i)
    {
        const int topology = populatedTopologies.at(i);
        const int bin = static_cast<int>(i) + 1;
        hTopology->SetBinContent(bin, topologyCounts.at(topology));
        hTopology->GetXaxis()->SetBinLabel(bin,topologyLabels.at(topology).c_str());
    }

    TCanvas* cTopology =new TCanvas("cTopology","Final-state topologies", 1800, 1100);

    cTopology->SetLeftMargin( 0.12);
    cTopology->SetRightMargin(0.05);
    cTopology->SetBottomMargin(0.30);

    hTopology->SetTitle("Background Final-State Topology Distribution After Preselection");
    hTopology->GetXaxis() ->SetTitle("True Final-State Topology");
    hTopology->GetYaxis()->SetTitle("Interactions");
    hTopology->GetXaxis()->SetTitleOffset(3.5);
    hTopology->SetMinimum(0.0);
    hTopology->SetMaximum( 1.20 *hTopology->GetMaximum());
    hTopology->LabelsOption("v","X" );
    hTopology->SetFillColorAlpha(fillColour,0.7);
    hTopology->SetLineColor(kBlack);
    hTopology->SetLineWidth(2);
    hTopology->Draw("HIST TEXT0");
    cTopology->Print("plots/new_logic/finalStateTopologiesBkg_new.png");

    // ========================================================================
    // COLOURS FOR STACKS
    // ========================================================================

    const std::array<
        int,
        kNModeCategories
    > modeColours = {
        kGray + 1,
        kAzure + 1,
        kOrange + 7,
        kGreen + 2,
        kRed - 4,
        kViolet - 4,
        kCyan + 1,
        kPink + 1,
        kYellow + 1,
        kTeal + 2,
        kBlue - 7,
        kOrange - 3,
        kSpring + 5,
        kMagenta - 2,
        kViolet + 1
    };

    const std::array<
        int,
        kNTopologies
    > topologyColours = {

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

    // ========================================================================
    // STACKED INTERACTION MODE VS NEUTRINO ENERGY
    // ========================================================================

    THStack* hsModeVsNuE = new THStack("hsModeVsNuE","Surviving Background by Interaction Mode;True E_{#nu} [GeV];Interactions");
    TLegend* legMode =new TLegend(0.76, 0.45, 0.95, 0.90);

    legMode->SetBorderSize(0);
    legMode->SetFillStyle(0);

    for (int modeIndex = 0;modeIndex < kNModeCategories;++modeIndex)
    {
        if (hNuEnergyByMode.at(modeIndex)->GetEntries() == 0) {
            continue;
        }

        hNuEnergyByMode.at(modeIndex)->SetFillColorAlpha(modeColours.at(modeIndex), 0.75);
        hNuEnergyByMode.at(modeIndex)->SetLineColor( kBlack );
        hsModeVsNuE->Add(hNuEnergyByMode.at(modeIndex));
        legMode->AddEntry(hNuEnergyByMode.at(modeIndex), modeLabels.at(modeIndex).c_str(), "f");

    }

    TCanvas* cModeVsNuE = new TCanvas("cModeVsNuE", "Interaction Mode vs Neutrino Energy", 1800, 1200);
    cModeVsNuE->SetRightMargin(0.25);

    hsModeVsNuE->Draw("HIST");
    hsModeVsNuE->GetXaxis()->SetTitle("True E_{#nu} [GeV]" );
    hsModeVsNuE->GetYaxis()->SetTitle("Interactions");
    legMode->Draw();
    cModeVsNuE->Print("plots/new_logic/modeVsNuEnergyBkg_new.png");

    // ========================================================================
    // STACKED TOPOLOGY VS NEUTRINO ENERGY
    // ========================================================================

    THStack* hsTopologyVsNuE = new THStack("hsTopologyVsNuE", "Surviving Background by Final-State Topology;True E_{#nu} [GeV];Interactions");
    TLegend* legTopology = new TLegend(0.74, 0.32, 0.97, 0.91);

    legTopology->SetBorderSize(0);
    legTopology->SetFillStyle(0);
    legTopology->SetTextSize(0.025);

    for (int topology = 0;topology < kNTopologies;++topology)
    {
        if (hNuEnergyByTopology.at(topology)->GetEntries() == 0
        ) {
            continue;
        }

        hNuEnergyByTopology.at(topology)->SetFillColorAlpha(topologyColours.at(topology), 0.75);
        hNuEnergyByTopology.at(topology)->SetLineColor(kBlack);
        hsTopologyVsNuE->Add(hNuEnergyByTopology.at(topology));
        legTopology->AddEntry(hNuEnergyByTopology.at(topology), topologyLabels.at(topology).c_str(), "f");
    }

    TCanvas* cTopologyVsNuE =new TCanvas("cTopologyVsNuE", "Topology vs Neutrino Energy", 1800, 1200);
    cTopologyVsNuE->SetRightMargin(0.30);

    hsTopologyVsNuE->Draw("HIST");
    hsTopologyVsNuE->GetXaxis()->SetTitle("True E_{#nu} [GeV]" );
    hsTopologyVsNuE->GetYaxis()->SetTitle("Interactions");

    legTopology->Draw();
    cTopologyVsNuE->Print("plots/new_logic/topologyVsNuEnergyBkg_new.png");

    // ========================================================================
    // MODE VS TOPOLOGY MATRIX
    // ========================================================================

    TCanvas* cModeVsTopology = new TCanvas("cModeVsTopology","Mode vs Topology", 2000, 1400);

    cModeVsTopology->SetRightMargin(0.15);
    cModeVsTopology->SetBottomMargin(0.30);
    cModeVsTopology->SetLeftMargin(0.15);

    hModeVsTopology->GetXaxis()->LabelsOption("v");
    hModeVsTopology->GetXaxis()->SetTitle("Final-State Topology");
    hModeVsTopology->GetYaxis()->SetTitle("Interaction Mode");
    hModeVsTopology->Draw("COLZ TEXT");

    cModeVsTopology->Print("plots/new_logic/modeVsTopologyBkg_new.png");


    // ========================================================================
    // FINAL-STATE MULTIPLICITY PLOTS
    // ========================================================================

    std::array<
        TH1*,
        5
    > multiplicityHists = {

        hNProtons,
        hNNeutrons,
        hNPions,
        hNPhotons,
        hNStrange
    };

    const std::array<
        std::string,
        5
    > multiplicityNames = {
        "Proton",
        "Neutron",
        "Pion",
        "Photon",
        "Strange Hadron"
    };


    const std::array<
        std::string,
        5
    > multiplicityFiles = {
        "plots/new_logic/nPrimaryProtons_new.png",
        "plots/new_logic/nPrimaryNeutrons_new.png",
        "plots/new_logic/nPrimaryPions_new.png",
        "plots/new_logic/nPrimaryPhotons_new.png",
        "plots/new_logic/nPrimaryStrange_new.png"
    };


    for (size_t i = 0;i < multiplicityHists.size();++i)
    {
        TCanvas* c = new TCanvas(Form("cMultiplicity_%zu", i),"", 1400, 1000);

        multiplicityHists.at(i)->GetXaxis()->SetTitle(Form("Number of Primary %ss",multiplicityNames.at(i).c_str()));
        multiplicityHists.at(i)->GetYaxis()->SetTitle("Interactions");
        multiplicityHists.at(i)->SetFillColorAlpha(fillColour, 0.7);
        multiplicityHists.at(i)->SetLineColor(kBlack);
        multiplicityHists.at(i)->Draw("HIST TEXT0");
        c->Print(multiplicityFiles.at(i).c_str());
    }


    // ========================================================================
    // LEADING PARTICLE MOMENTA
    // ========================================================================

    std::array<
        TH1F*,
        3
    > momentumHists = {
        hLeadingProtonP,
        hLeadingPionP,
        hLeadingPhotonP
    };


    const std::array<
        std::string,
        3
    > momentumTitles = {

        "Leading Primary Proton Momentum",
        "Leading Primary Pion Momentum",
        "Leading Primary Photon Momentum"
    };


    const std::array<
        std::string,
        3
    > momentumFiles = {
        "plots/new_logic/leadingProtonP_new.png",
        "plots/new_logic/leadingPionP_new.png",
        "plots/new_logic/leadingPhotonP_new.png"
    };


    for (size_t i = 0;i < momentumHists.size();++i)
    {
        TCanvas* c = new TCanvas(Form("cMomentum_%zu", i), "", 1400, 1000);

        momentumHists.at(i)->SetTitle(momentumTitles.at(i).c_str());
        momentumHists.at(i)->GetXaxis()->SetTitle("True Momentum [GeV/c]");
        momentumHists.at(i)->GetYaxis()->SetTitle("Particles");
        momentumHists.at(i)->SetFillColorAlpha(fillColour, 0.7);
        momentumHists.at(i)->SetLineColor(kBlack);
        momentumHists.at(i)->Draw("HIST");

        c->Print(momentumFiles.at(i).c_str());
    }
}