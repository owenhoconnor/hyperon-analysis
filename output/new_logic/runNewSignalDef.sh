root -b -l -x << EOF
.L newSignalDef.C
TChain *chain = new TChain("ana/tree")
chain->Add("/data/ooconnor/sbnd/hyperons/analyzer_output/new_logic/analyzer_output_new_2026_all.root")
newSignalDef t(chain);
t.Loop();
.q
EOF
