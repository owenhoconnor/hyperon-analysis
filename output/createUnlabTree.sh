root -b -l -x << EOF
.L signalDef.C
TChain *chain = new TChain("ana/tree")
chain->Add("hyperons.root")
signalDef t(chain);
t.Loop();
.q
EOF
