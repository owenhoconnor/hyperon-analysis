root -l -b -x << EOF
.L backgroundPlots.C
TChain *chain = new TChain("tree")
chain->AddFile("/data/ooconnor/sbnd/hyperons/preselection_output/eventClassification_output.root", TChain::kBigNumber, "tree")
backgroundPlots a(chain)
a.Loop()
.q
EOF
