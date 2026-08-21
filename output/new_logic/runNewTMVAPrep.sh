root -l -b -x << EOF
.L newTMVAPrep.C
TChain *chain = new TChain("tree")
chain->Add("/data/ooconnor/sbnd/hyperons/preselection_output/signalDef_output_sig.root")
chain->Add("/data/ooconnor/sbnd/hyperons/preselection_output/signalDef_output_bkg.root")
newTMVAPrep t(chain)
t.Loop()
.q
EOF
