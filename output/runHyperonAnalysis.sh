root -b -l -x << EOF
.L signalDef.C
TChain *chain = new TChain("ana/tree")
chain->Add("mergedSig.root")
chain->Add("mergedBkg.root")
signalDef t(chain);
t.Loop();
.q
EOF
