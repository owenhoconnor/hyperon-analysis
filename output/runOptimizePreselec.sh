#!/usr/bin/env bash
set -euo pipefail

echo "==================================="
echo " Running SBND Cut Optimisation"
echo "==================================="

# Run ROOT in batch mode
root -l -b <<'EOF'
.L preselection.C+
.L optimizePreselection.C+
optimizePreselection()
.q
EOF

echo "==================================="
echo " Optimisation finished"
echo "==================================="

