#include "preselection.h"
#include <TChain.h>
#include <iostream>
#include <iomanip>
#include <cmath>

void runPreselection()
{
    // Input data
    TChain chain("tree");
    chain.Add("/data/ooconnor/sbnd/hyperons/preselection_output/eventClassification_output.root");

    // create instance of preselection class
    preselection selector(&chain);

    PreselectionCuts cuts; // Starting cut values (use defaults in preselection.h)
    PreselectionResults finalResults = selector.Loop(cuts);

    std::cout << "\n========= CUTS =========\n";

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