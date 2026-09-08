root -l -b -x << EOF
.L backgroundPlots.C
TChain *chain = new TChain("tree")
chain->AddFile("/data/ooconnor/sbnd/hyperons/preselection_output/tmvaSample_sig.root", TChain::kBigNumber, "sigTree")
chain->AddFile("/data/ooconnor/sbnd/hyperons/preselection_output/tmvaSample_bkg.root", TChain::kBigNumber, "bkgTree")
chain->AddFile("/data/ooconnor/sbnd/hyperons/preselection_output/tmvaSample_cosmic.root", TChain::kBigNumber, "cosmicTree")
backgroundPlots a(chain)
a.Loop()
.q
EOF
