#ifndef recoStudy_h
#define recoStudy_h

// Minimal MakeClass-style reader for recoStudy.C.
// Pass the TTree/TChain explicitly; no stale hard-coded input file.
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>
#include <vector>

class recoStudy {
public:
    TTree* fChain = nullptr;
    Int_t fCurrent = -1;

    Int_t sampleType = -1;
    std::vector<int>* truePDG = nullptr;
    std::vector<int>* trueGeneration = nullptr;
    std::vector<int>* trueMotherPDG = nullptr;
    std::vector<float>* trueP = nullptr;
    std::vector<int>* trueIsReconstructed = nullptr;
    std::vector<int>* trueNHits = nullptr;
    std::vector<int>* trueNHitsU = nullptr;
    std::vector<int>* trueNHitsV = nullptr;
    std::vector<int>* trueNHitsZ = nullptr;
    std::vector<int>* trueBestRecoPfpIdx = nullptr;
    std::vector<int>* pfpKey = nullptr;
    std::vector<int>* pfpIsNuSlice = nullptr;
    std::vector<int>* pfpIsPrimary = nullptr;
    std::vector<int>* pfpSelfID = nullptr;
    std::vector<int>* pfpParentID = nullptr;

    explicit recoStudy(TTree* tree = nullptr) { Init(tree); }
    virtual ~recoStudy() = default; // Input tree/file owned by caller.

    virtual Int_t GetEntry(Long64_t entry) {
        return fChain ? fChain->GetEntry(entry) : 0;
    }
    virtual Long64_t LoadTree(Long64_t entry) {
        if (!fChain) return -5;
        Long64_t localEntry = fChain->LoadTree(entry);
        if (localEntry >= 0 && fChain->GetTreeNumber() != fCurrent) {
            fCurrent = fChain->GetTreeNumber();
            Notify();
        }
        return localEntry;
    }
    virtual void Init(TTree* tree) {
        if (!tree) {
            std::cerr << "recoStudy: supply a TTree/TChain, e.g. recoStudy study(tree);\n";
            return;
        }
        fChain = tree;
        fCurrent = -1;
        fChain->SetMakeClass(1);

        // Fail safely when run on a filtered or differently structured tree.
        const char* required[] = {
            "sampleType", "truePDG", "trueGeneration", "trueMotherPDG",
            "trueP", "trueIsReconstructed", "trueNHits", "trueNHitsU",
            "trueNHitsV", "trueNHitsZ", "trueBestRecoPfpIdx",
            "pfpKey", "pfpIsNuSlice", "pfpIsPrimary", "pfpSelfID", "pfpParentID"
        };
        bool missing = false;
        for (const char* name : required) {
            if (!fChain->GetBranch(name)) {
                std::cerr << "recoStudy: missing required branch: " << name << '\n';
                missing = true;
            }
        }
        if (missing) {
            std::cerr << "recoStudy: incompatible input tree; analysis disabled.\n";
            fChain = nullptr;
            return;
        }

        fChain->SetBranchAddress("sampleType", &sampleType);
        fChain->SetBranchAddress("truePDG", &truePDG);
        fChain->SetBranchAddress("trueGeneration", &trueGeneration);
        fChain->SetBranchAddress("trueMotherPDG", &trueMotherPDG);
        fChain->SetBranchAddress("trueP", &trueP);
        fChain->SetBranchAddress("trueIsReconstructed", &trueIsReconstructed);
        fChain->SetBranchAddress("trueNHits", &trueNHits);
        fChain->SetBranchAddress("trueNHitsU", &trueNHitsU);
        fChain->SetBranchAddress("trueNHitsV", &trueNHitsV);
        fChain->SetBranchAddress("trueNHitsZ", &trueNHitsZ);
        fChain->SetBranchAddress("trueBestRecoPfpIdx", &trueBestRecoPfpIdx);
        fChain->SetBranchAddress("pfpKey", &pfpKey);
        fChain->SetBranchAddress("pfpIsNuSlice", &pfpIsNuSlice);
        fChain->SetBranchAddress("pfpIsPrimary", &pfpIsPrimary);
        fChain->SetBranchAddress("pfpSelfID", &pfpSelfID);
        fChain->SetBranchAddress("pfpParentID", &pfpParentID);
        Notify();
    }
    virtual bool Notify() { return true; }
    virtual void Show(Long64_t entry = -1) {
        if (fChain) fChain->Show(entry);
    }
    virtual Int_t Cut(Long64_t) { return 1; }
    virtual void Loop();
};

#endif // recoStudy_h
