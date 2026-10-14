#define preselection_cxx
#include "preselection.h"
#include <cmath>
#include <iostream>

PreselectionResults preselection::Loop(const PreselectionCuts& cuts)
{
    PreselectionResults results;

    if (!fChain) return results;

    float nPassNuScoreSig = 0.0;
    float nPassNuScoreHyp = 0.0;
    float nPassNuScoreBkg = 0.0;
    float nPassNuScoreDirt = 0.0;
    float nPassNuScoreCosmic = 0.0;
    float nPassBCFMSig = 0.0;
    float nPassBCFMHyp = 0.0;
    float nPassBCFMBkg = 0.0;
    float nPassBCFMDirt = 0.0;
    float nPassBCFMCosmic = 0.0;
    float nPassRecoFVSig = 0.0;
    float nPassRecoFVHyp = 0.0;
    float nPassRecoFVBkg = 0.0;
    float nPassRecoFVDirt = 0.0;
    float nPassRecoFVCosmic = 0.0;
    float nPassTrackMultiplicitySig = 0.0;
    float nPassTrackMultiplicityHyp = 0.0;
    float nPassTrackMultiplicityBkg = 0.0;
    float nPassTrackMultiplicityDirt = 0.0;
    float nPassTrackMultiplicityCosmic = 0.0;
    float nPassShowerMultiplicitySig = 0.0;
    float nPassShowerMultiplicityHyp = 0.0;
    float nPassShowerMultiplicityBkg = 0.0;
    float nPassShowerMultiplicityDirt = 0.0;
    float nPassShowerMultiplicityCosmic = 0.0;

    const Long64_t nentries = fChain->GetEntries();

    for (Long64_t jentry = 0; jentry < nentries; ++jentry) {

        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;

        fChain->GetEntry(jentry);

        bool passNuScore = false;
        bool passBCFM = false;
        bool passRecoFV = false;
        bool passTrackMultiplicity = false;
        bool passShowerMultiplicity = false;

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
        {
            passRecoFV = true;
        }

        // Neutrino score
        if (nuSliceNuScore > cuts.nuScore.value)
        {
            passNuScore = true;
            //nPassNuScore += eventWeight
        }

        // Barycentre flash matching score
        if (nuSliceBCFMScore > cuts.bcfm.value)
        {
                passBCFM = true;
                //nPassBCFM += eventWeight
        }

        // Track / shower multiplicity
        if (nTracks == cuts.nTracks)
        {
            passTrackMultiplicity = true;
        }

        if (nShowers == cuts.nShowers)
        {
            passShowerMultiplicity = true;
        }

        if(passNuScore)
        {
            if(sampleType==0){nPassNuScoreSig += eventWeight;}
            if(sampleType==4){nPassNuScoreHyp += eventWeight;}
            if(sampleType==1){nPassNuScoreBkg += eventWeight;}
            if(sampleType==2){nPassNuScoreDirt += eventWeight;}
            if(sampleType==3){nPassNuScoreCosmic += eventWeight;}
            
            if(passBCFM)
            {
                if(sampleType==0){nPassBCFMSig += eventWeight;}
                if(sampleType==4){nPassBCFMHyp += eventWeight;}
                if(sampleType==1){nPassBCFMBkg += eventWeight;}
                if(sampleType==2){nPassBCFMDirt += eventWeight;}
                if(sampleType==3){nPassBCFMCosmic += eventWeight;}

                if(passRecoFV)
                {
                    if(sampleType==0){nPassRecoFVSig += eventWeight;}
                    if(sampleType==4){nPassRecoFVHyp += eventWeight;}
                    if(sampleType==1){nPassRecoFVBkg += eventWeight;}
                    if(sampleType==2){nPassRecoFVDirt += eventWeight;}
                    if(sampleType==3){nPassRecoFVCosmic += eventWeight;}

                    if(passTrackMultiplicity)
                    {
                        if(sampleType==0){nPassTrackMultiplicitySig += eventWeight;}
                        if(sampleType==4){nPassTrackMultiplicityHyp += eventWeight;}
                        if(sampleType==1){nPassTrackMultiplicityBkg += eventWeight;}
                        if(sampleType==2){nPassTrackMultiplicityDirt += eventWeight;}
                        if(sampleType==3){nPassTrackMultiplicityCosmic += eventWeight;}

                        if(passShowerMultiplicity)
                        {
                            if(sampleType==0){nPassShowerMultiplicitySig += eventWeight;}
                            if(sampleType==4){nPassShowerMultiplicityHyp += eventWeight;}
                            if(sampleType==1){nPassShowerMultiplicityBkg += eventWeight;}
                            if(sampleType==2){nPassShowerMultiplicityDirt += eventWeight;}
                            if(sampleType==3){nPassShowerMultiplicityCosmic += eventWeight;}
                            //nPassShowerMultiplicty += eventWeight;
                        }
                    }
                }
            }

        }
        // Event survives all preselection cuts
        if(passNuScore && passBCFM && passRecoFV & passTrackMultiplicity && passShowerMultiplicity)
        {
            results.nSelected++;

            if (sampleType == 0)
            {    
                results.nSelectedSignal++;
            }
        }
    }

    std::cout<<"------------------- Results ----------------------------"<<std::endl;
    std::cout<<"# Signal passing NuScore Cut = "<<nPassNuScoreSig<<std::endl;
    std::cout<<"# Hyperons passing NuScore Cut = "<<nPassNuScoreHyp<<std::endl;
    std::cout<<"# Bkg passing NuScore cut = "<<nPassNuScoreBkg<<std::endl;
    std::cout<<"# Dirt passing NuScore cut = "<<nPassNuScoreDirt<<std::endl;
    std::cout<<"# Cosmic passing NuScore cut = "<<nPassNuScoreCosmic<<std::endl;

    return results;
}