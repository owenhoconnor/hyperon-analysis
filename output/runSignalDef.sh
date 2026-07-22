root -b -l -x << EOF
.L signalDef.C
TChain *chain = new TChain("ana/tree")
chain->Add("mergedSig.root")
chain->Add("/data/ooconnor/sbnd/hyperons/hyperon_analyzer_output/2026_prod/analyzer_output_prod2026.root")
signalDef t(chain);
t.Loop();
.q
EOF
