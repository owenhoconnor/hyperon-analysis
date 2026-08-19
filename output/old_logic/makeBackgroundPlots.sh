root -l -b -x << EOF
.L backgroundPlots.C
TChain *chain = new TChain("bkgTree")
chain->Add("/data/ooconnor/sbnd/hyperons/preselection_output/tmvaSample_bkg.root")
backgroundPlots a(chain)
a.Loop()
.q
EOF