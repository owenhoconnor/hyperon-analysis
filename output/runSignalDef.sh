root -b -l -x << EOF
.L signalDef.C
TChain *chain = new TChain("ana/tree")

chain->Add("/data/ooconnor/sbnd/hyperons/analyzer_output/merged_anaOut_2026.root")
signalDef t(chain);
t.Loop();
.q
EOF
