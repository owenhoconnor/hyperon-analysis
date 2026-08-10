root -b -l -x << EOF
.L signalDef.C
TChain *chain = new TChain("ana/tree")

chain->Add("/data/ooconnor/sbnd/hyperons/analyzer_output/combined_anaOut_prod2026_new.root")
signalDef t(chain);
t.Loop();
.q
EOF
