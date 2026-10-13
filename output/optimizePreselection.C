#include "preselection.h"
#include <TChain.h>
#include <iostream>
#include <iomanip>
#include <cmath>

std::vector<Cut PreselectionCuts::*> parameters = {
    &PreselectionCuts::nuScore,
    &PreselectionCuts::bcfm,
    &PreselectionCuts::flashTime
};


void optimizePreselection()
{
    // Input data
    TChain chain("tree"); // Replace with your tree path
    chain.Add("/data/ooconnor/sbnd/hyperons/preselection_output/eventClassification_output.root");

    // create instance of preselection class
    preselection selector(&chain);

    // Starting cut values
    PreselectionCuts cuts;

    const int nSteps = 100;
    const int maxIters = 10;
    double highestEffPur = -1.0;

    // Optimize NuScore and BCFM sequentially
    for (const auto& parameter : parameters) {

        Cut& currentCut = cuts.&parameter;

        std::cout<<"------------------------------------------"<<std::endl;
        std::cout<<"Optimizing Cut: "<<currentCut.name <<std::endl;
        std::cout<<"Cut Min Value: "<<currentCut.minValue<<std::endl;
        std::cout<<"Cut Max Value: "<<currentCut.maxValue<<std::endl;

        double bestCutValue = currentCut.value
        double bestEffPur = -1.0;

        for (int iStep = 0; iStep <= currentCut.nSteps; ++iStep) {

            const double stepFrac = static_cast<double>(iStep) / currentCut.nSteps;
            const double stepCutValue = currentCut.minValue + stepFrac * (currentCut.maxValue - currentCut.minValue);

            std::cout<<"Step number = "<<iStep<<std::endl;
            std::cout<<"Step fraction = "<<stepFrac<<std::endl;
            std::cout<<"Step cut value = "<<stepCutValue<<std::endl;

            // Make a copy of the current configuration
            PreselectionCuts trialCuts = cuts;

            (trialCuts.*parameter).value = stepCutValue;

            // Run the entire preselection
            PreselectionResults results = selector.Loop(trialCuts);
            const double effPur = results.effPur();

            if (effPur > bestEffPur) {
                bestEffPur = effPur;
                bestCutValue = stepFrac;
                std::cout<<"New highest eff*pur, new bestCutValue"<<std::endl;
            }

            std::cout<<"Efficiency * Purity for this step = "<<effPur<<std::endl;
            std::cout<<"Best efficiency*purity = "<<bestEffPur<<std::endl;
            std::cout<<"Best cut value = "<<bestCutValue<<std::endl;
        }

        // Update the optimised cut
        (cuts.*parameter).value = bestCutValue; 

        std::cout <<(cuts.*parameter).name<<" best cut = " << bestCutValue<< ", Eff*Pur = " << bestEffPur<< std::endl;
    }

    PreselectionResults finalResults = selector.Loop(cuts);
    highestEffPur = finalResults.effPur();

    std::cout << "\n========= OPTIMAL CUTS =========\n";

    std::cout<< "NuScore: " << cuts.nuScore.value <<std::endl;;
    std::cout<< "BCFM: " << cuts.bcfm.value <<std::endl;
    std::cout<<"FlashTime: "<<std::endl;

    std::cout << "Selected signal: "<< finalResults.nSelectedSignal << "\n";
    std::cout << "Total signal: "<< finalResults.nTrueSignal << "\n";
    std::cout << "Total selected: "<< finalResults.nSelected << "\n";
    std::cout << "Efficiency: "<< finalResults.efficiency() << "\n";
    std::cout << "Purity: "<< finalResults.purity() << "\n";
    std::cout << "Eff * Pur: "<< finalResults.effPur() << "\n";
}