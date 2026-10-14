#!/usr/bin/env python3
"""Fetch NCBI protein reference sets for conservancy analysis.

    python setup/fetch_ncbi_conservancy.py

Reads the stage 1 manifest to discover every (protein, strain) pair, then
queries NCBI's protein database for each one.  The search term is the manifest's
``protein_name`` and ``strain`` columns verbatim, filtered to sequences within
±``length_margin`` residues of the input sequence length.  Downloaded records
are filtered (drop sequences containing ``X``, deduplicate by sequence) and
written as a clean FASTA plus a provenance JSON.

This is a **setup script**, not a DAG rule.  NCBI contents change over time;
fetching inside the DAG would make two runs of the same pipeline silently
incomparable.  Running the fetch once and committing the provenance file makes
the target set explicit and reproducible.

Output layout (under ``resources/conservancy/``)::

    resources/conservancy/Spike/Sudan_ebolavirus.fasta
    resources/conservancy/Spike/Sudan_ebolavirus.provenance.json

Rate limits: 3 requests/second anonymous, 10/second with an NCBI API key.
NCBI asks automated clients to send ``tool`` and ``email``.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Project root & paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESOURCE_DIR = PROJECT_ROOT / "resources" / "conservancy"
MANIFEST_PATH = PROJECT_ROOT / "1 - Sequence_retrieval" / "manifest.tsv"


# ---------------------------------------------------------------------------
# Config loading (no pyyaml dependency — only stdlib is available on the host)
# ---------------------------------------------------------------------------

def _load_config() -> dict:
    """Load config.yaml with PyYAML if available, otherwise fall back to a
    minimal parser for the fields we need.
    """
    config_path = PROJECT_ROOT / "config.yaml"
    try:
        import yaml  # type: ignore[import-untyped]

        with open(config_path, encoding="utf-8") as fh:
            return yaml.safe_load(fh)
    except ImportError:
        pass

    return _parse_reference_sets_block(config_path)


def _parse_reference_sets_block(path: Path) -> dict:
    """Rough YAML reader — enough for the reference_sets control keys."""
    text = path.read_text(encoding="utf-8")
    enabled = False
    email = ""
    api_key = ""
    length_margin = 5

    in_ref = False
    indent_ref = 0
    for line in text.splitlines():
        stripped = line.lstrip()
        indent = len(line) - len(stripped)

        if stripped.startswith("reference_sets:"):
            in_ref = True
            indent_ref = indent
            continue
        if in_ref and indent <= indent_ref and stripped and not stripped.startswith("#"):
            in_ref = False
            continue
        if not in_ref:
            continue

        if stripped.startswith("enabled:"):
            enabled = stripped.split(":", 1)[1].strip().lower() in ("true", "yes")
        elif stripped.startswith("email:"):
            email = stripped.split(":", 1)[1].strip().strip('"').strip("'")
        elif stripped.startswith("api_key:"):
            api_key = stripped.split(":", 1)[1].strip().strip('"').strip("'")
        elif stripped.startswith("length_margin:"):
            try:
                length_margin = int(stripped.split(":", 1)[1].strip())
            except ValueError:
                pass

    return {
        "conservation": {
            "reference_sets": {
                "enabled": enabled,
                "email": email,
                "api_key": api_key,
                "length_margin": length_margin,
            }
        }
    }


# ---------------------------------------------------------------------------
# Manifest reading
# ---------------------------------------------------------------------------

def _read_manifest() -> list[dict[str, str]]:
    """Read the stage 1 manifest TSV."""
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(
            f"Manifest not found: {MANIFEST_PATH}\n"
            "Run  snakemake --use-conda --cores 1 "
            f"'{MANIFEST_PATH.relative_to(PROJECT_ROOT)}'\n"
            "to generate it (stage 1 only), then re-run this script."
        )
    with open(MANIFEST_PATH, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _discover_pairs(
    manifest: list[dict[str, str]],
) -> list[dict[str, str | int]]:
    """Derive one fetch entry per unique (protein_label, strain) from the manifest.

    Each entry carries the protein_name (for the search term), the strain, and
    the sequence length (for the ±margin window).  When two manifest rows share
    the same (protein_label, strain) — which shouldn't happen, but is handled —
    the first one wins.
    """
    seen: set[tuple[str, str]] = set()
    pairs: list[dict[str, str | int]] = []
    for row in manifest:
        key = (row["protein_label"], row["strain"])
        if key in seen:
            continue
        seen.add(key)
        pairs.append({
            "protein_label": row["protein_label"],
            "protein_name": row["protein_name"],
            "strain": row["strain"],
            "length": int(row["length"]),
        })
    return pairs


# ---------------------------------------------------------------------------
# NCBI E-utilities helpers
# ---------------------------------------------------------------------------

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL_NAME = "mpvpredpip"
BATCH_SIZE = 500           # efetch retmax
MAX_RETRIES = 3


def _request_json(url: str) -> dict:
    """GET *url*, parse JSON, respect rate limits."""
    for attempt in range(MAX_RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                return json.loads(resp.read())
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt == MAX_RETRIES - 1:
                raise
            print(f"  retry {attempt + 1}/{MAX_RETRIES}: {exc}")
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def _request_text(url: str) -> str:
    """GET *url*, return text."""
    for attempt in range(MAX_RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                return resp.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt == MAX_RETRIES - 1:
                raise
            print(f"  retry {attempt + 1}/{MAX_RETRIES}: {exc}")
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def _rate_pause(api_key: str) -> None:
    """Sleep to stay within NCBI rate limits: 3/s anonymous, 10/s with key."""
    time.sleep(0.12 if api_key else 0.35)


def esearch(term: str, email: str, api_key: str) -> dict:
    """Search the protein database and return the search result with history."""
    params = {
        "db": "protein",
        "term": term,
        "retmode": "json",
        "retmax": 0,
        "usehistory": "y",
        "tool": TOOL_NAME,
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    url = f"{EUTILS}/esearch.fcgi?{urllib.parse.urlencode(params)}"
    data = _request_json(url)
    result = data.get("esearchresult", {})

    if "ERROR" in result:
        raise RuntimeError(f"NCBI esearch error: {result['ERROR']}")

    return {
        "count": int(result.get("count", 0)),
        "webenv": result.get("webenv", ""),
        "query_key": result.get("querykey", ""),
        "query_translation": result.get("querytranslation", ""),
        "id_list": result.get("idlist", []),
    }


def efetch_fasta(
    webenv: str,
    query_key: str,
    count: int,
    email: str,
    api_key: str,
) -> str:
    """Fetch all records as FASTA using the history server."""
    chunks: list[str] = []
    for start in range(0, count, BATCH_SIZE):
        params = {
            "db": "protein",
            "query_key": query_key,
            "WebEnv": webenv,
            "rettype": "fasta",
            "retmode": "text",
            "retstart": start,
            "retmax": BATCH_SIZE,
            "tool": TOOL_NAME,
        }
        if email:
            params["email"] = email
        if api_key:
            params["api_key"] = api_key

        url = f"{EUTILS}/efetch.fcgi?{urllib.parse.urlencode(params)}"
        text = _request_text(url)
        chunks.append(text)
        fetched = min(start + BATCH_SIZE, count)
        print(f"    fetched {fetched}/{count}")
        _rate_pause(api_key)

    return "\n".join(chunks)


# ---------------------------------------------------------------------------
# FASTA filtering — mirrors the user's manual fasta_filter.py
# ---------------------------------------------------------------------------

def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Parse a FASTA string into ``(header, sequence)`` pairs."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks).upper()))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line)

    if header is not None:
        records.append((header, "".join(chunks).upper()))
    return records


def filter_fasta(
    records: list[tuple[str, str]],
) -> tuple[list[tuple[str, str]], dict[str, int]]:
    """Apply the user's manual filtering rules.

    1. Drop sequences containing ``X`` (ambiguous residue — partial or
       uncertain translations that would depress conservancy scores).
    2. Drop exact duplicate sequences (first occurrence wins).

    Returns the filtered list and a stats dict.
    """
    stats: dict[str, int] = {"input": len(records), "dropped_x": 0, "dropped_dup": 0}

    # Pass 1: drop sequences containing X
    no_x: list[tuple[str, str]] = []
    for header, seq in records:
        if "X" in seq:
            stats["dropped_x"] += 1
        else:
            no_x.append((header, seq))

    # Pass 2: deduplicate by sequence (first occurrence wins)
    seen_md5: set[str] = set()
    unique: list[tuple[str, str]] = []
    for header, seq in no_x:
        md5 = hashlib.md5(seq.encode()).hexdigest()
        if md5 in seen_md5:
            stats["dropped_dup"] += 1
        else:
            seen_md5.add(md5)
            unique.append((header, seq))

    stats["output"] = len(unique)
    return unique, stats


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 80) -> None:
    """Write records as FASTA, wrapped at *width* columns."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def fetch_one(
    protein: str,
    strain: str,
    protein_name: str,
    seq_length: int,
    length_margin: int,
    email: str,
    api_key: str,
) -> Path:
    """Fetch, filter and write one reference set.  Returns the FASTA path."""
    min_len = seq_length - length_margin
    max_len = seq_length + length_margin
    base_term = f"{protein_name} {strain}"
    full_term = f"{base_term} AND {min_len}:{max_len}[SLEN]"

    # Slugify the strain for the filename (spaces -> underscores)
    strain_slug = re.sub(r"[^A-Za-z0-9]+", "_", strain).strip("_")

    out_dir = RESOURCE_DIR / protein
    fasta_path = out_dir / f"{strain_slug}.fasta"
    prov_path = out_dir / f"{strain_slug}.provenance.json"

    if fasta_path.exists():
        print(f"  already present: {fasta_path.relative_to(PROJECT_ROOT)}")
        print("  delete it to re-fetch.")
        return fasta_path

    print(f"\n  searching: {full_term}")
    result = esearch(full_term, email, api_key)
    count = result["count"]
    print(f"  NCBI count: {count}")
    print(f"  query translation: {result['query_translation']}")

    if count == 0:
        raise RuntimeError(
            f"No records found for {protein}/{strain}.\n"
            f"  term: {full_term}\n"
            "  Check that the protein name and strain produce results on NCBI."
        )

    _rate_pause(api_key)

    # Fetch all records
    print(f"  fetching {count} records...")
    raw_fasta = efetch_fasta(
        result["webenv"], result["query_key"], count, email, api_key,
    )

    # Parse and extract accessions before filtering
    raw_records = parse_fasta(raw_fasta)
    if len(raw_records) != count:
        print(
            f"  warning: expected {count} records from NCBI, parsed {len(raw_records)}"
        )

    raw_accessions = [h.split()[0] for h, _ in raw_records]

    # Filter
    filtered, stats = filter_fasta(raw_records)
    print(
        f"  filtering: {stats['input']} raw -> "
        f"{stats['dropped_x']} with X, {stats['dropped_dup']} duplicates -> "
        f"{stats['output']} kept"
    )

    if not filtered:
        raise RuntimeError(
            f"All {count} records were filtered out for {protein}/{strain}.\n"
            "  This likely means the search term is returning wrong organism records."
        )

    # Write
    write_fasta(fasta_path, filtered)
    print(f"  wrote: {fasta_path.relative_to(PROJECT_ROOT)}")

    # Provenance
    provenance = {
        "protein": protein,
        "strain": strain,
        "protein_name": protein_name,
        "input_length": seq_length,
        "length_margin": length_margin,
        "term": full_term,
        "base_term": base_term,
        "min_length": min_len,
        "max_length": max_len,
        "query_translation": result["query_translation"],
        "ncbi_count": count,
        "post_filter_count": len(filtered),
        "dropped_x": stats["dropped_x"],
        "dropped_duplicate": stats["dropped_dup"],
        "accessions_raw": raw_accessions,
        "accessions_kept": [h.split()[0] for h, _ in filtered],
        "fetch_date": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "tool": TOOL_NAME,
        "email": email or "(not set)",
    }
    with open(prov_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(provenance, fh, indent=2)
        fh.write("\n")
    print(f"  wrote: {prov_path.relative_to(PROJECT_ROOT)}")

    return fasta_path


def main() -> int:
    config = _load_config()
    ref_cfg = config.get("conservation", {}).get("reference_sets", {})

    if not ref_cfg.get("enabled"):
        print("conservation.reference_sets.enabled is false — nothing to fetch.")
        print("Set it to true in config.yaml to enable NCBI reference sets.")
        return 0

    email = ref_cfg.get("email", "")
    api_key = ref_cfg.get("api_key", "")
    length_margin = int(ref_cfg.get("length_margin", 5))

    if not email:
        print(
            "WARNING: no email set in conservation.reference_sets.email.\n"
            "  NCBI asks automated clients to identify themselves.\n"
            "  The fetch will proceed, but please add your email.\n"
        )

    # Discover all (protein, strain) pairs from the manifest
    manifest = _read_manifest()
    pairs = _discover_pairs(manifest)

    print("Fetching NCBI conservancy reference sets")
    print(f"  source: {MANIFEST_PATH.relative_to(PROJECT_ROOT)}")
    print(f"  output: {RESOURCE_DIR.relative_to(PROJECT_ROOT)}/")
    print(f"  length margin: ±{length_margin} residues")
    if api_key:
        print("  rate limit: 10 req/s (API key set)")
    else:
        print("  rate limit: 3 req/s (no API key)")
    print(f"  pairs to fetch: {len(pairs)}")
    print()

    errors: list[str] = []
    for pair in pairs:
        protein = pair["protein_label"]
        strain = pair["strain"]
        print(f"--- {protein} / {strain} (input length {pair['length']}) ---")
        try:
            fetch_one(
                protein=protein,
                strain=strain,
                protein_name=str(pair["protein_name"]),
                seq_length=int(pair["length"]),
                length_margin=length_margin,
                email=email,
                api_key=api_key,
            )
        except Exception as exc:
            msg = f"{protein}/{strain}: {exc}"
            errors.append(msg)
            print(f"  ERROR: {exc}")
        _rate_pause(api_key)

    print()
    if errors:
        print(f"{len(errors)} error(s):")
        for e in errors:
            print(f"  - {e}")
        return 1

    print("Done.  Next:  snakemake --use-conda --cores 8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
