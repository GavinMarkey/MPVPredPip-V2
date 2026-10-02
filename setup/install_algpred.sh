#!/usr/bin/env bash
# Install AlgPred 2.0 (allergenicity prediction) into tools/algpred2.
#
#   bash setup/install_algpred.sh
#
# AlgPred 2.0 is freely available from the Raghava group (GPL-3.0):
#   https://github.com/raghavagps/algpred2
#   https://webs.iiitd.edu.in/raghava/algpred2/stand.html
#
# tools/README.md used to describe this as a manual download. It is not - and
# that same claim was already wrong for BepiPred and for all five IEDB tools.
#
# WHY NOT `pip install algpred2`
# ------------------------------
# A pip wheel exists and is the easiest route, but it is only enough for model 1
# (amino-acid-composition Random Forest). This pipeline runs model 2, the hybrid
# of RF + BLAST + MERCI, because that is the model config.yaml's threshold is
# characterised against. The hybrid needs two things that live in the repository
# rather than in the wheel:
#
#   Database/   BLAST database and IgE motif files
#   progs/      MERCI, a Perl program
#
# So the whole package is cloned. BLAST and Perl come from
# workflow/envs/screening.yaml. Switch screening.allergenicity.model to 1 if you
# would rather avoid all of that - the adapter handles either.

set -euo pipefail

REPO_URL="https://github.com/raghavagps/algpred2"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="${PROJECT_ROOT}/tools/algpred2"

echo "Installing AlgPred 2.0"
echo "  from: ${REPO_URL}"
echo "  into: ${TARGET}"
echo

if ! command -v git >/dev/null 2>&1; then
    echo "ERROR: git is not installed." >&2
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

# --- what the hybrid model needs --------------------------------------------
missing=0
for required in algpred2.py Database progs; do
    if [ ! -e "${TARGET}/${required}" ]; then
        echo "WARNING: ${required} is missing from the clone." >&2
        missing=1
    fi
done

# --- rf_model, which git clone alone does NOT get you ------------------------
# .gitattributes tracks rf_model with Git LFS, so a plain clone leaves either
# nothing or a ~130-byte pointer stub. Without the real file the Random Forest
# cannot load and AlgPred dies inside joblib with an error that never mentions
# LFS. Fetch it properly here rather than discovering it at stage 6.
MODEL="${TARGET}/rf_model"

is_real_model() {
    [ -f "${MODEL}" ] || return 1
    # A pointer stub starts with this line; the model is a joblib pickle.
    ! head -c 23 "${MODEL}" 2>/dev/null | grep -q "version https://git-lfs"
}

if is_real_model; then
    echo
    echo "rf_model present ($(wc -c < "${MODEL}") bytes)."
else
    echo
    echo "rf_model is missing or is an unfetched LFS pointer - retrieving it."

    if command -v git-lfs >/dev/null 2>&1 || git lfs version >/dev/null 2>&1; then
        echo "  trying: git lfs pull"
        git -C "${TARGET}" lfs install --local >/dev/null 2>&1 || true
        git -C "${TARGET}" lfs pull || true
    else
        echo "  git-lfs is not installed, skipping to the zip mirror."
    fi

    if ! is_real_model; then
        ZIP_URL="https://webs.iiitd.edu.in/raghava/algpred2/algpred2_new.zip"
        echo "  trying: ${ZIP_URL}"
        TMP="$(mktemp -d)"
        trap 'rm -rf "${TMP}"' EXIT
        if command -v curl >/dev/null 2>&1; then
            curl -fsSL "${ZIP_URL}" -o "${TMP}/algpred2.zip" || true
        elif command -v wget >/dev/null 2>&1; then
            wget -q "${ZIP_URL}" -O "${TMP}/algpred2.zip" || true
        fi
        if [ -s "${TMP}/algpred2.zip" ]; then
            # The zip's internal layout varies between mirrors, so find the
            # model wherever it sits rather than assuming a prefix.
            unzip -q -o "${TMP}/algpred2.zip" -d "${TMP}/x" || true
            FOUND="$(find "${TMP}/x" -name 'rf_model' -type f | head -1)"
            if [ -n "${FOUND}" ]; then
                cp "${FOUND}" "${MODEL}"
                echo "  recovered rf_model from the zip."
            fi
        fi
    fi

    if is_real_model; then
        echo "  rf_model is now in place ($(wc -c < "${MODEL}") bytes)."
    else
        echo "ERROR: could not obtain rf_model." >&2
        echo "       Download this zip by hand and copy rf_model into ${TARGET}:" >&2
        echo "       https://webs.iiitd.edu.in/raghava/algpred2/algpred2_new.zip" >&2
        missing=1
    fi
fi

if [ -f "${TARGET}/envfile" ]; then
    echo
    echo "envfile ships pointing at the author's own machine (/home/neelam/...)."
    echo "Leave it alone: the adapter writes a correct one per run, resolving"
    echo "blastp from PATH inside the conda environment. Nothing to edit here."
fi

cat <<'EOF'

------------------------------------------------------------------------------
Installed.
------------------------------------------------------------------------------
Verified CLI (from the upstream README):

    algpred2.py -i INPUT [-o OUTPUT] [-t THRESHOLD] [-m {1,2}] [-d {1,2}]

The adapter passes -m 2 -d 2 explicitly. Both matter:

  -m 2  hybrid model (RF + BLAST + MERCI). AlgPred defaults to 1.

  -d 2  report EVERY peptide. AlgPred defaults to 1, which writes out only the
        peptides it calls allergenic - under that default the peptides that
        pass the screen would be missing from the output entirely, and "nothing
        was allergenic" would look identical to "the tool never ran". The
        adapter cross-checks the output against the submitted FASTA and raises
        if any peptide is absent, so this cannot fail quietly.

Output columns were read off the source (algpred2.py, hybrid()):

    Subject,ML Score,MERCI Score,BLAST Score,Hybrid Score,Prediction

The adapter takes "Hybrid Score", not "ML Score", and trusts the tool's own
Allergen / Non-Allergen verdict rather than recomputing one.

Two things the adapter handles so you do not have to:

  * envfile is rewritten per run against this machine. The shipped copy points
    at /home/neelam/... and would fail; leave it as it is.
  * each run happens in its own temporary directory. algpred2.py reads and
    writes about ten files by bare relative name, so screening two proteins at
    once from one directory would corrupt both.

Then:  python setup/check_tools.py
------------------------------------------------------------------------------
EOF

if [ "${missing}" -ne 0 ]; then
    echo "Finished with warnings - see the ERROR/WARNING lines above." >&2
    exit 1
fi
