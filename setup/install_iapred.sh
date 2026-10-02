#!/usr/bin/env bash
# Install IAPred (antigenicity prediction) into tools/iapred.
#
#   bash setup/install_iapred.sh
#
# Unlike the other predictors in this pipeline, IAPred is on GitHub and can be
# fetched unattended:  https://github.com/sebamiles/IAPred
#
# After installing, run  python setup/check_tools.py  to confirm it resolves, and
# read the note it prints about the adapter's unverified command line.

set -euo pipefail

REPO_URL="https://github.com/sebamiles/IAPred"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="${PROJECT_ROOT}/tools/iapred"

echo "Installing IAPred"
echo "  from: ${REPO_URL}"
echo "  into: ${TARGET}"
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
    git clone --depth 1 "${REPO_URL}" "${TARGET}"
fi

echo
echo "Contents of ${TARGET}:"
ls -la "${TARGET}"

# The repository may declare its own dependencies. Install them into whichever
# environment is active, rather than assuming a particular one.
for requirements in "${TARGET}/requirements.txt" "${TARGET}/requirements-dev.txt"; do
    if [ -f "${requirements}" ]; then
        echo
        echo "Installing dependencies from $(basename "${requirements}")"
        python -m pip install -r "${requirements}"
    fi
done

if [ -f "${TARGET}/setup.py" ] || [ -f "${TARGET}/pyproject.toml" ]; then
    echo
    echo "The repository is packaged; installing it into the active environment."
    python -m pip install -e "${TARGET}"
fi

cat <<'EOF'

------------------------------------------------------------------------------
Installed.
------------------------------------------------------------------------------
The adapter has been confirmed against a real run of this tool. Two things it
found are worth knowing, because the repository's own README disagrees with both:

  * the score column is Intrinsic_Antigenicity_Score. The README documents an
    older name, IAscore, which no longer appears in the output.
  * IApred truncates deflines to 20 characters, which would collide ours. The
    adapter feeds it surrogate IDs and maps the results back, so this is
    handled - but it means the headers in the raw output are not the real ones.

The adapter tries several argument styles automatically and matches output
columns by name. If a future version breaks it, you have two options, neither
of which touches anything downstream:

  1. Correct _IAPRED_ENTRYPOINTS / _IAPRED_ARG_STYLES / _IAPRED_*_COLUMNS at the
     top of src/vaxpipe/adapters/screening.py

  2. Bypass the adapter entirely - in config.yaml set

         antigenicity:
           tool: "external"
           command: "python tools/iapred/IApred.py {input} {output}"

Also confirm the score threshold. config.yaml uses 0.3 for both the stage 2
screen and the stage 6 filter, which is IApred's own "High antigenicity"
boundary. Note its scores run roughly -3 to 3, NOT 0 to 1, so negative values are
normal. That value decides which candidates survive.

Then:  python setup/check_tools.py
------------------------------------------------------------------------------
EOF
