#!/usr/bin/env bash
set -euo pipefail
root -l -b <<'ROOT'
.L preselection.C+
.L runPreselection.C+
runPreselection()
.q
ROOT
