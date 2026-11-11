#include "preselection.h"
#include <TChain.h>
#include <TSystem.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
const char* sampleNames[kNSamples] = {"signal", "hyperon", "beam", "dirt", "cosmic"};
const char* stageNames[kNStages] = {"initial", "reco_fv", "nu_score", "bcfm", "tracks", "showers"};
const char* reasonNames[kNInvalidReasons] = {
    "no_nu_slice", "missing_slice_vector", "slice_index_out_of_range",
    "nonfinite_reco_value", "invalid_weight"
};

void writeCount(std::ostream& out, const WeightedCount& c) {
    out << c.raw << ',' << c.sumW << ',' << c.sumW2;
}

std::ofstream openCSV(const std::string& filename) {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot write " + filename);
    out << std::setprecision(17);
    return out;
}
}

void runPreselection()
{
    TChain chain("tree");
    chain.Add("/data/ooconnor/sbnd/hyperons/preselection_output/eventClassification_output.root");
    if (chain.GetEntries() <= 0) {
        std::cerr << "No input entries found.\n";
        return;
    }

    PreselectionCuts cuts;
    preselection selector(&chain);
    const PreselectionResults r = selector.Loop(cuts);

    // Change this path if the ROOT scripts are run outside the presentations repo.
    const std::string outputDir = "/data/ooconnor/sbnd/hyperons/preselection_output/csv";
    if (gSystem->mkdir(outputDir.c_str(), true) != 0 &&
        gSystem->AccessPathName(outputDir.c_str())) {
        throw std::runtime_error("Could not create output directory: " + outputDir);
    }

    // One row per stage/sample, with exact cuts and their parameters in metadata.
    {
        auto out = openCSV(outputDir + "/preselection_cumulative.csv");
        out << "stage,sample,raw,sumW,sumW2,efficiency,purity\n";
        for (int s = 0; s < kNStages; ++s) {
            for (int c = 0; c < kNSamples; ++c) {
                out << stageNames[s] << ',' << sampleNames[c] << ',';
                writeCount(out, r.cutflow[s][c]);
                out << ',' << r.efficiency(static_cast<SelectionStage>(s))
                    << ',' << r.purity(static_cast<SelectionStage>(s)) << '\n';
            }
        }
    }
    {
        auto out = openCSV(outputDir + "/preselection_independent.csv");
        out << "cut,sample,pass_raw,pass_sumW,pass_sumW2,fail_raw,fail_sumW,fail_sumW2\n";
        for (int s = 1; s < kNStages; ++s) {
            for (int c = 0; c < kNSamples; ++c) {
                out << stageNames[s] << ',' << sampleNames[c] << ',';
                writeCount(out, r.independentPass[s][c]);
                out << ',';
                writeCount(out, r.independentFail[s][c]);
                out << '\n';
            }
        }
    }
    {
        auto out = openCSV(outputDir + "/preselection_patterns.csv");
        out << "pattern,pass_reco_fv,pass_nu_score,pass_bcfm,pass_tracks,pass_showers,sample,raw,sumW,sumW2\n";
        for (int p = 0; p < kNCutPatterns; ++p) {
            for (int c = 0; c < kNSamples; ++c) {
                out << p;
                for (int bit = 0; bit < kNPhysicsCuts; ++bit)
                    out << ',' << ((p >> bit) & 1);
                out << ',' << sampleNames[c] << ',';
                writeCount(out, r.cutPatterns[p][c]);
                out << '\n';
            }
        }
    }
    {
        auto out = openCSV(outputDir + "/preselection_invalid.csv");
        out << "reason,sample,raw,sumW,sumW2\n";
        for (int reason = 0; reason < kNInvalidReasons; ++reason) {
            for (int c = 0; c < kNSamples; ++c) {
                out << reasonNames[reason] << ',' << sampleNames[c] << ',';
                writeCount(out, r.invalidEvents[reason][c]);
                out << '\n';
            }
        }
        // These are totals rather than per-category entries; weight is undefined.
        out << "invalid_sample_type,all," << r.invalidSampleType << ",,\n";
        out << "invalid_weight_all,all," << r.invalidWeight << ",,\n";
        out << "invalid_slice_unique,all," << r.invalidSlice << ",,\n";
        out << "unreadable_entries,all," << r.unreadableEntries << ",,\n";
        out << "valid_slice,all," << r.validSlice << ",,\n";
    }
    {
        auto out = openCSV(outputDir + "/preselection_metadata.csv");
        out << "parameter,value\n"
            << "target_POT,1e21\n"
            << "reco_vtx_x_abs_max," << cuts.recoVtxXMax << '\n'
            << "reco_vtx_y_abs_max," << cuts.recoVtxYMax << '\n'
            << "reco_vtx_z_min," << cuts.recoVtxZMin << '\n'
            << "reco_vtx_z_max," << cuts.recoVtxZMax << '\n'
            << "nu_score_min_exclusive," << cuts.nuScore.value << '\n'
            << "bcfm_min_exclusive," << cuts.bcfm.value << '\n'
            << "n_tracks_exact," << cuts.nTracks << '\n'
            << "n_showers_exact," << cuts.nShowers << '\n'
            << "dent_applied,false\n"
            << "flash_time_applied,false\n"
            << "weight_definition,eventWeight_POT_scale_only\n"
            << "sample_type_0,signal\n"
            << "sample_type_1,hyperon\n"
            << "sample_type_2,beam\n"
            << "sample_type_3,dirt\n"
            << "sample_type_4,cosmic\n"
            << "pattern_bit_0,reco_fv\n"
            << "pattern_bit_1,nu_score\n"
            << "pattern_bit_2,bcfm\n"
            << "pattern_bit_3,tracks\n"
            << "pattern_bit_4,showers\n";
    }
    // Compatibility with the existing presentations Python script: wide-format
    // CSV with one row per cumulative stage and expected weighted yields.
    {
        auto out = openCSV(outputDir + "/preselection.csv");
        out << "cut,signal,hyperon,beam,cosmic,dirt\n";
        for (int s = 0; s < kNStages; ++s) {
            out << stageNames[s] << ','
                << r.cutflow[s][kSignal].sumW << ','
                << r.cutflow[s][kHyperon].sumW << ','
                << r.cutflow[s][kBeam].sumW << ','
                << r.cutflow[s][kDirt].sumW << ','
                << r.cutflow[s][kCosmic].sumW << '\n';
        }
    }

    std::cout << "\nFinal efficiency=" << r.efficiency()
              << " purity=" << r.purity()
              << " eff*pur=" << r.effPur()
              << "\nWrote CSV outputs to " << outputDir << "/\n";
}
