root -b -l -x << EOF
.L eventClassification.C
TChain *chain = new TChain("ana/tree")
chain->Add("/data/ooconnor/sbnd/hyperons/analyzer_output/merged_anaOut_hyperons.root") // tree number = 0
chain->Add("/data/ooconnor/sbnd/hyperons/analyzer_output/merged_anaOut_2026.root") // tree number = 1
eventClassification t(chain);
t.Loop();
.q
EOF
