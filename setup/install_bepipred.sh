#!/usr/bin/env bash
# Install BepiPred-3.0 (linear B-cell epitope prediction) into tools/bepipred3.
#
#   bash setup/install_bepipred.sh
#
# Source: https://github.com/UberClifford/BepiPred-3.0
#
# Despite what this project's notes said for a long time, BepiPred-3.0 does NOT
# require a licence form or registration. DTU Health Tech publish the code and
# the trained weights on GitHub, free for academic and other non-commercial use.
# Only for-profit use needs a separate licence (contact morni@dtu.dk). So it
# installs unattended, like IAPred.
#
# Two things are downloaded, and only the first happens here:
#
#   1. the repository, about 140 MB, because the trained ensemble weights
#      (bp3/BP3Models) are committed to it;
#   2. the ESM-2 language model, about 2.5 GB, fetched automatically by fair-esm
#      the first time a prediction runs - not by this script.
#
# Python dependencies are NOT installed here either: they belong to the rule's
# conda environment, workflow/envs/bepipred.yaml, which Snakemake builds on the
# first run. This script only places the source tree.

set -euo pipefail

REPO_URL="https://github.com/UberClifford/BepiPred-3.0"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="${PROJECT_ROOT}/tools/bepipred3"

echo "Installing BepiPred-3.0"
echo "  from: ${REPO_URL}"
echo "  into: ${TARGET}"
echo "  size: ~140 MB (the repository carries the trained model weights)"
echo

if ! command -v git >/dev/null 2>&1; then
    echo "ERROR: git is not installed." >&2
    echo "Inside the dev container it is already present; on the host, install git first." >&2
    exit 1
fi

if [ -d "${TARGET}/.git" ]; then
    echo "Already present - updating."
    git -C "${TARGET}" pull --ff-only
else
    mkdir -p "$(dirname "${TARGET}")"
    # --depth 1: the history is not needed and roughly doubles the transfer.
    git clone --depth 1 "${REPO_URL}" "${TARGET}"
fi

echo
echo "Checking the layout the adapter expects:"
missing=0
for required in "bepipred3_CLI.py" "bp3/bepipred3.py" "bp3/BP3Models"; do
    if [ -e "${TARGET}/${required}" ]; then
        echo "  ok       ${required}"
    else
        echo "  MISSING  ${required}"
        missing=1
    fi
done

# The weights are what make this repository large; a partial clone that silently
# dropped them would only fail much later, in the middle of a prediction run.
folds=$(find "${TARGET}/bp3/BP3Models" -name '*Fold*' 2>/dev/null | wc -l | tr -d ' ')
echo "  ${folds} model fold file(s) found"
if [ "${folds}" -eq 0 ]; then
    missing=1
fi

if [ "${missing}" -ne 0 ]; then
    echo
    echo "ERROR: the clone is not the layout the adapter expects." >&2
    echo "Check ${TARGET} by hand before running stage 3." >&2
    exit 1
fi

cat <<'EOF'

------------------------------------------------------------------------------
Installed.
------------------------------------------------------------------------------
The adapter was written against this repository's own source (bp3/bepipred3.py),
not its README, so the output columns and the command line should be right:

  raw_output.csv is  Accession,Residue,BepiPred-3.0 score,BepiPred-3.0 linear
  epitope score  - one row per residue, in sequence order, with NO position
  column. The pipeline reads the raw score, because that is the column
  BepiPred's own -t threshold is compared against.

Still to happen on the first run, unattended:

  * conda builds workflow/envs/bepipred.yaml (torch, fair-esm)
  * fair-esm downloads the ESM-2 weights, about 2.5 GB, into the torch hub
    cache. In the dev container TORCH_HOME points at a named volume, so this
    happens once and survives a container rebuild.

Then:  python setup/check_tools.py
------------------------------------------------------------------------------
EOF
