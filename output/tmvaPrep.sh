root -l -b -x << EOF
.L varPrep.C
TChain *chain = new TChain("preTree")
chain->Add("TreeS.root")
chain->Add("TreeB.root")
varPrep t(chain)
t.Loop()
.q
EOF
