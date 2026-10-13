#define preselection_cxx
#include "preselection.h"
#include <cmath>
#include <iostream>

PreselectionResults preselection::Loop(const PreselectionCuts& cuts)
{
    PreselectionResults results;

    if (!fChain) return results;

    const Long64_t nentries = fChain->GetEntries();

    for (Long64_t jentry = 0; jentry < nentries; ++jentry) {

        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;

        fChain->GetEntry(jentry);

        // Count all true signal events BEFORE applying cuts.
        if (sampleType == 0) {
            results.nTrueSignal++;
        }

        // Validate the selected neutrino slice.
        if (!sliceNuScore || !sliceBaryFlashScore)
            continue;

        if (nuSliceIdx < 0 ||
            static_cast<size_t>(nuSliceIdx) >= sliceNuScore->size() ||
            static_cast<size_t>(nuSliceIdx) >= sliceBaryFlashScore->size())
            continue;

        const float recoVtxX = sliceVtxX->at(nuSliceIdx);
        const float recoVtxY = sliceVtxY->at(nuSliceIdx);
        const float recoVtxZ = sliceVtxZ->at(nuSliceIdx);
        const double nuSliceNuScore = sliceNuScore->at(nuSliceIdx);
        const double nuSliceBCFMScore = sliceBaryFlashScore->at(nuSliceIdx);

        // Reconstructed fiducial volume
        if (std::abs(recoVtxX) >= cuts.recoVtxXMax || std::abs(recoVtxY) >= cuts.recoVtxYMax ||
            recoVtxZ <= cuts.recoVtxZMin ||
            recoVtxZ >= cuts.recoVtxZMax)
            continue;

        // Neutrino score
        if (nuSliceNuScore < cuts.nuScore.value)
            continue;

        // Barycentre flash matching score
        if (nuSliceBCFMScore < cuts.bcfm.value)
            continue;

        // Track / shower multiplicity
        if (nTracks != cuts.nTracks)
            continue;

        if (nShowers != cuts.nShowers)
            continue;

        // Event survives all preselection cuts
        results.nSelected++;

        if (sampleType == 0) {
            results.nSelectedSignal++;
        }
    }

    return results;
}