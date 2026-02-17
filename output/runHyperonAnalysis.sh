root -b -l -x << EOF
.L signalDef.C
signalDef t;
t.Loop();
.q
EOF
