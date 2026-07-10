root -l -b -x << EOF
.L plotEfficiencyP.C
TChain *chain = new TChain("unlabTree")
chain->Add("unlabTree.root")
blindPre t(chain)
t.Loop()
.q
EOF
