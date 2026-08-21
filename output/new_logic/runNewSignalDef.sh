root -b -l -x << EOF
.L newSignalDef.C
TChain *chain = new TChain("ana/tree")
chain->Add("/data/ooconnor/sbnd/hyperons/analyzer_output/analyzer_output_hyperons_combined.root")
chain->Add("/data/ooconnor/sbnd/hyperons/analyzer_output/analyzer_output_2026_combined.root")
newSignalDef t(chain);
t.Loop();
.q
EOF
