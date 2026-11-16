#define preselection_cxx
#include "preselection.h"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>

PreselectionResults preselection::Loop(const PreselectionCuts& cuts)
{
    PreselectionResults results;
    if (!fChain) return results;

    const Long64_t nentries = fChain->GetEntries();

    for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
        if (LoadTree(jentry) < 0) break;
        if (GetEntry(jentry) <= 0) {
            ++results.unreadableEntries;
            continue;
        }

        const bool knownCategory = sampleType >= 0 && sampleType < kNSamples;
        const double weight = static_cast<double>(eventWeight);
        const bool validWeight = std::isfinite(weight) && weight >= 0.0;

        if (!knownCategory) ++results.invalidSampleType;
        if (!validWeight) ++results.invalidWeight;
        if (!knownCategory) continue; // No safe sample category for bookkeeping.

        const int category = sampleType;
        if (!validWeight) {
            // A raw invalid-weight count; sumW and sumW2 cannot be meaningful.
            ++results.invalidEvents[kInvalidWeight][category].raw;
            continue;
        }

        results.cutflow[kInitial][category].fill(weight);

        // Check every applicable slice-validation failure, not just the first.
        const bool noSlice = nuSliceIdx < 0;
        const bool missingVector = !sliceNuScore || !sliceBaryFlashScore || !sliceBaryFlashTime ||
                                   !sliceVtxX || !sliceVtxY || !sliceVtxZ;
        bool outOfRange = false;
        bool nonFinite = false;

        if (!noSlice && !missingVector) {
            const std::size_t idx = static_cast<std::size_t>(nuSliceIdx);
            outOfRange = idx >= sliceNuScore->size() ||
                         idx >= sliceBaryFlashScore->size() ||
                         idx >= sliceVtxX->size() ||
                         idx >= sliceVtxY->size() ||
                         idx >= sliceVtxZ->size();
            if (!outOfRange) {
                nonFinite = !std::isfinite(static_cast<double>(sliceVtxX->at(idx))) ||
                            !std::isfinite(static_cast<double>(sliceVtxY->at(idx))) ||
                            !std::isfinite(static_cast<double>(sliceVtxZ->at(idx))) ||
                            !std::isfinite(static_cast<double>(sliceNuScore->at(idx))) ||
                            !std::isfinite(static_cast<double>(sliceBaryFlashScore->at(idx)));
            }
        }

        if (noSlice) results.invalidEvents[kNoNuSlice][category].fill(weight);
        if (missingVector) results.invalidEvents[kMissingSliceVector][category].fill(weight);
        if (outOfRange) results.invalidEvents[kSliceIndexOutOfRange][category].fill(weight);
        if (nonFinite) results.invalidEvents[kNonFiniteRecoValue][category].fill(weight);

        if (noSlice || missingVector || outOfRange || nonFinite) {
            ++results.invalidSlice;
            continue; // Physics scores/vertex cannot safely be evaluated.
        }

        ++results.validSlice;
        const std::size_t idx = static_cast<std::size_t>(nuSliceIdx);
        const double x = sliceVtxX->at(idx);
        const double y = sliceVtxY->at(idx);
        const double z = sliceVtxZ->at(idx);
        const double nuScore = sliceNuScore->at(idx);
        const double bcfm = sliceBaryFlashScore->at(idx);
        const double flashTime = sliceBaryFlashTime->at(idx);

        // Keep the previous RecoFV geometry (dent not applied).
        const bool pass[kNPhysicsCuts] = {
            std::abs(x) < cuts.recoVtxXMax &&
            std::abs(y) < cuts.recoVtxYMax &&
            z > cuts.recoVtxZMin && z < cuts.recoVtxZMax,
            nuScore > cuts.nuScore.value,
            bcfm > cuts.bcfm.value,
            flashTime >= 0.18 && flashTime <= 1.73,
            nTracks == cuts.nTracks,
            nShowers == cuts.nShowers
        };

        unsigned int pattern = 0;
        bool cumulativePass = true;
        for (int cut = 0; cut < kNPhysicsCuts; ++cut) {
            const int stage = cut + 1;
            if (pass[cut]) {
                results.independentPass[stage][category].fill(weight);
                pattern |= (1u << cut);
            } else {
                results.independentFail[stage][category].fill(weight);
            }
            cumulativePass = cumulativePass && pass[cut];
            if (cumulativePass) results.cutflow[stage][category].fill(weight);
        }
        results.cutPatterns[pattern][category].fill(weight);
    }

    const char* stageNames[kNStages] = {
        "Initial", "RecoFV", "nuScore", "BCFM", "FlashTime", "Tracks", "Showers"
    };
    const char* sampleNames[kNSamples] = {
        "Signal", "Other Hyperons", "Beam", "Dirt", "Cosmic"
    };

    std::cout << std::setprecision(10)
              << "\n=== CUMULATIVE CUT FLOW (expected events / 1e21 POT) ===\n";
    for (int stage = 0; stage < kNStages; ++stage) {
        std::cout << "\n" << stageNames[stage] << ":\n";
        for (int sample = 0; sample < kNSamples; ++sample) {
            const auto& c = results.cutflow[stage][sample];
            std::cout << "  " << std::setw(15) << sampleNames[sample]
                      << " raw=" << c.raw << " sumW=" << c.sumW
                      << " sumW2=" << c.sumW2 << '\n';
        }
        const auto s = static_cast<SelectionStage>(stage);
        std::cout << "  efficiency=" << results.efficiency(s)
                  << " purity=" << results.purity(s) << '\n';
    }
    std::cout << "\nInvalid sampleType=" << results.invalidSampleType
              << " invalid weight=" << results.invalidWeight
              << " invalid slice (unique)=" << results.invalidSlice
              << " unreadable entries=" << results.unreadableEntries
              << " evaluable slices=" << results.validSlice << '\n';
    return results;
}
