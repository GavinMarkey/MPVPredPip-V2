#!/usr/bin/env bash
# Install the five IEDB standalone tools used by stages 4, 5 and 8.
#
#   bash setup/install_iedb.sh              # all five
#   bash setup/install_iedb.sh mhc_i        # just one, by target directory name
#
# tools/README.md claimed for a long time that these were "gated behind a
# registration or licence click, so they cannot be scripted". That is wrong -
# the same mistake it made about BepiPred. All five are plain HTTP downloads
# from downloads.iedb.org, each with an MD5SUM file beside it. They are free for
# non-commercial use.
#
# About 1.9 GB in total, nearly all of it the two binding predictors, which
# bundle a separate trained model set per method.
#
# The configure step differs per tool and is taken from each one's own README:
# mhc_i has `configure`, mhc_ii has `configure.py`, population_coverage has both
# but its README names `configure`. immunogenicity and cluster have none. The
# cluster tool is downloaded for reference only - stage 5 no longer uses it, see
# src/vaxpipe/cluster.py - so it can be skipped with an explicit tool list.

set -euo pipefail

BASE="https://downloads.iedb.org/tools"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS="${PROJECT_ROOT}/tools"

# target_dir | remote subdirectory | archive | configure command ("" = none)
# The target directory names are the ones config.yaml points at under `tools:`.
TOOLS_LIST=(
    "mhc_i|mhci|IEDB_MHC_I-3.1.7.tar.gz|./configure"
    "mhc_ii|mhcii|IEDB_MHC_II-3.1.12.tar.gz|./configure.py"
    "immunogenicity|immunogenicity|IEDB_Immunogenicity-3.0.tar.gz|"
    "cluster|cluster|IEDB_Cluster-1.0.tar.gz|"
    # Its README.txt does say to run configure, despite the tool being small.
    "population_coverage|population|IEDB_Population_Coverage-3.0.2.tar.gz|./configure"
)

WANTED=("$@")

want() {
    [ ${#WANTED[@]} -eq 0 ] && return 0
    for w in "${WANTED[@]}"; do [ "$w" = "$1" ] && return 0; done
    return 1
}

install_one() {
    local target="$1" remote="$2" archive="$3" configure="$4"
    local dest="${TOOLS}/${target}"
    local url="${BASE}/${remote}/LATEST/${archive}"

    echo
    echo "=============================================================="
    echo "${target}"
    echo "  from: ${url}"
    echo "  into: ${dest}"

    if [ -d "${dest}" ] && [ -n "$(ls -A "${dest}" 2>/dev/null)" ]; then
        echo "  Already present - skipping. Delete the directory to reinstall."
        return 0
    fi

    local staging
    staging="$(mktemp -d)"
    # Clean up the staging directory whatever happens, so a failed download
    # cannot leave a half-unpacked tool behind that would then be "already
    # present" on the next run.
    trap 'rm -rf "${staging}"' RETURN

    echo "  Downloading."
    curl -fL --progress-bar -o "${staging}/${archive}" "${url}"

    # The MD5SUM file covers every archive in that directory, so pick our line
    # out of it. These are large binary downloads over a long transfer; a silent
    # truncation would surface much later as an unreadable model file.
    echo "  Verifying checksum."
    if curl -fsL "${BASE}/${remote}/LATEST/MD5SUM" -o "${staging}/MD5SUM" 2>/dev/null; then
        local expected actual
        expected="$(grep -F "${archive}" "${staging}/MD5SUM" | awk '{print $1}' | head -1)"
        actual="$(md5sum "${staging}/${archive}" | awk '{print $1}')"
        if [ -z "${expected}" ]; then
            echo "  WARNING: no entry for ${archive} in MD5SUM; cannot verify."
        elif [ "${expected}" != "${actual}" ]; then
            echo "  ERROR: checksum mismatch." >&2
            echo "    expected ${expected}" >&2
            echo "    got      ${actual}" >&2
            return 1
        else
            echo "  Checksum ok."
        fi
    else
        echo "  WARNING: MD5SUM not retrievable; cannot verify."
    fi

    echo "  Unpacking."
    tar -xzf "${staging}/${archive}" -C "${staging}"

    # Each archive contains a single top-level directory, but the name is not
    # always the one config.yaml expects, so find it rather than assume it.
    local unpacked
    unpacked="$(find "${staging}" -mindepth 1 -maxdepth 1 -type d ! -name 'LATEST' | head -1)"
    if [ -z "${unpacked}" ]; then
        echo "  ERROR: the archive did not contain a directory." >&2
        return 1
    fi
    echo "  Archive contains: $(basename "${unpacked}")"

    mkdir -p "${TOOLS}"
    mv "${unpacked}" "${dest}"

    if [ -n "${configure}" ]; then
        local script_name="${configure#./}"
        if [ -f "${dest}/${script_name}" ]; then
            echo "  Running ${configure} (sets absolute paths to the trained models)."
            # Must run from inside the tool directory: configure writes the paths
            # it finds relative to its own working directory.
            #
            # A .py configure script is run through `python` rather than executed
            # directly, because mhc_ii's carries the shebang #!/usr/bin/python3,
            # which does not exist in a conda environment - executing it fails
            # with "bad interpreter" and the tool is then silently left
            # unconfigured. Whatever python is on PATH here is the rule's own
            # interpreter, which is the one that will run the tool later anyway.
            case "${script_name}" in
                *.py) ( cd "${dest}" && python "${script_name}" ) ;;
                *)    ( cd "${dest}" && chmod +x "${script_name}" 2>/dev/null || true
                        cd "${dest}" && "./${script_name}" ) ;;
            esac
            echo "  Configured."
        else
            echo "  NOTE: expected ${configure} but it is not there; skipping."
            echo "        Check ${dest}/README before running stage 4."
        fi
    else
        # If a future release adds one, say so rather than silently skipping it.
        for candidate in configure configure.py; do
            if [ -f "${dest}/${candidate}" ]; then
                echo "  NOTE: ${dest}/${candidate} exists but this script did not"
                echo "        expect one here. Read the README and run it if needed."
            fi
        done
    fi

    echo "  Done."
}

echo "Installing the IEDB standalone tools (~1.9 GB in total)"

for entry in "${TOOLS_LIST[@]}"; do
    IFS='|' read -r target remote archive configure <<< "${entry}"
    want "${target}" || continue
    install_one "${target}" "${remote}" "${archive}" "${configure}"
done

cat <<'EOF'

------------------------------------------------------------------------------
Installed. Expect a shakedown.
------------------------------------------------------------------------------
None of the IEDB adapters has ever been executed, and the scripts that drive
them were converted from Snakemake's script: directive to argparse without ever
being run. BepiPred needed four corrections on its first real run even after its
adapter was written from the tool's own source. Assume the same here.

The things most likely to be wrong, in order:

  1. the class I command line. The adapter assumes
        ./src/predict_binding.py <method> <alleles> <lengths> <fasta>
     with alleles and lengths as PARALLEL comma-separated lists of equal
     length. That matches this version's README, but has not been run.
  2. the output column names, resolved by name in
     src/vaxpipe/adapters/iedb.py. This is what cost the most time with IAPred:
     its README documented a column the tool no longer emitted.
  3. the class II tool uses configure.py where class I uses configure - already
     handled above, but a sign the two are not as symmetrical as they look.

Run one protein through stage 4 first and read the raw output before trusting a
full run.

Then:  python setup/check_tools.py
------------------------------------------------------------------------------
EOF
