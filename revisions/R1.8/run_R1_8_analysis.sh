#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$repo_root"

python reviewer_response/02_analysis/pathway_reanalysis/run_pathway_reanalysis.py
python -m pytest reviewer_response/02_analysis/pathway_reanalysis/test_pathway_reanalysis.py -q

