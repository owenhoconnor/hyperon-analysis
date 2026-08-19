root -l -b -x << EOF
.L tmvaPrep.C
TChain *chain = new TChain("tree")
chain->Add("/data/ooconnor/sbnd/hyperons/preselection_output/signalDef_output_sig.root")
chain->Add("/data/ooconnor/sbnd/hyperons/preselection_output/signalDef_output_bkg_test.root")
tmvaPrep t(chain)
t.Loop()
.q
EOF
