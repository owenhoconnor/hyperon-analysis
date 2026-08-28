root -b -l -x << EOF
.L newSignalDef.C
TChain *chain = new TChain("ana/tree")

chain->Add("/data/ooconnor/sbnd/hyperons/analyzer_output/merged_anaOut_2026.root")
newSignalDef t(chain);
t.Loop();
.q
EOF
