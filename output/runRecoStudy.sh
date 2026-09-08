root -b -l -x << EOF
.L recoStudy.C
TChain *chain = new TChain("tree")
chain->Add("/data/ooconnor/sbnd/hyperons/preselection_output/signalDef_output_sig.root")
recoStudy t(chain);
t.Loop();
.q
EOF
