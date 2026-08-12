root -b -l -x << EOF
.L signalDef.C
TChain *chain = new TChain("ana/tree")

chain->Add("/data/ooconnor/sbnd/hyperons/analyzer_output/analyser_output_2026_all.root")
signalDef t(chain);
t.Loop();
.q
EOF
