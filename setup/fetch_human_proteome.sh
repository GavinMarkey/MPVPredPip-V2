#!/usr/bin/env bash
# Fetch the human reference proteome for the stage 6 autoimmunity screen.
#
#   bash setup/fetch_human_proteome.sh
#
# Stage 6 flags any consensus peptide sharing an exact substring of
# screening.autoimmunity.min_exact_match_length residues (default 9) with a human
# protein, on the reasoning that such a peptide risks provoking an autoimmune
# response rather than a useful one. That check needs the human proteome.
#
# This is UniProt's own reference-proteome file, so the name config.yaml expects
# is the name upstream uses - UP000005640 is the human reference proteome and
# 9606 is the NCBI taxon ID. One canonical protein per gene, about 20,600
# sequences.
#
# Not included: UP000005640_9606_additional.fasta.gz, a further 41 MB of isoforms
# and unreviewed entries. Screening against the canonical set is the convention.
# Adding the isoforms would make the screen stricter, and can only remove
# candidates, never add them - if you want that, download it alongside and
# concatenate the two.

set -euo pipefail

BASE_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640"
ARCHIVE="UP000005640_9606.fasta.gz"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET_DIR="${PROJECT_ROOT}/resources/human"
TARGET="${TARGET_DIR}/UP000005640_9606.fasta"

echo "Fetching the human reference proteome"
echo "  from: ${BASE_URL}/${ARCHIVE}"
echo "  into: ${TARGET}"
echo "  size: ~7.7 MB compressed, ~75 MB on disk"
echo

if [ -s "${TARGET}" ]; then
    echo "Already present: ${TARGET}"
    echo "Delete it to re-download."
    exit 0
fi

mkdir -p "${TARGET_DIR}"

# Download beside the target and move into place only once complete, so an
# interrupted transfer cannot leave a truncated proteome that would silently
# under-report autoimmune matches - a screen that quietly passes everything is
# worse than one that fails.
TMP="${TARGET}.partial.gz"
trap 'rm -f "${TMP}"' EXIT

if command -v curl >/dev/null 2>&1; then
    curl -fL --progress-bar -o "${TMP}" "${BASE_URL}/${ARCHIVE}"
elif command -v wget >/dev/null 2>&1; then
    wget -O "${TMP}" "${BASE_URL}/${ARCHIVE}"
else
    echo "ERROR: neither curl nor wget is available." >&2
    exit 1
fi

echo "Decompressing."
gunzip -c "${TMP}" > "${TARGET}"

count=$(grep -c '^>' "${TARGET}" || true)
echo
echo "  ${count} sequences in $(basename "${TARGET}")"

# The human reference proteome is ~20,600 canonical entries. An order-of-
# magnitude miss means the download or the decompression went wrong, and the
# autoimmunity screen would then pass peptides it should have flagged.
if [ "${count}" -lt 15000 ]; then
    echo >&2
    echo "ERROR: expected roughly 20,600 sequences, got ${count}." >&2
    echo "The file looks truncated. Delete it and re-run." >&2
    rm -f "${TARGET}"
    exit 1
fi

echo
echo "Done. Next:  python setup/check_tools.py"
