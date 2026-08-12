root -l -b -x << EOF
.L newBackgroundPlots.C
TChain *chain = new TChain("bkgTree")
chain->Add("/data/ooconnor/sbnd/hyperons/preselection_output/new_logic/tmvaSample_bkg.root")
newBackgroundPlots a(chain)
a.Loop()
.q
EOF