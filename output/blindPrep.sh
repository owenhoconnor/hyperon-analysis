root -l -b -x << EOF
.L blindPre.C
TChain *chain = new TChain("unlabTree")
chain->Add("unlabTree.root")
blindPre t(chain)
t.Loop()
.q
EOF
