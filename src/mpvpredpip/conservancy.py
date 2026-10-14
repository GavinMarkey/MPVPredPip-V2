"""Stage 7 - epitope conservancy across the supplied strains.

The IEDB conservancy analysis tool is web-only. Its algorithm is simple and fully
specified, so it is reimplemented here to keep the pipeline unattended, while
``mode: both`` in ``config.yaml`` also writes upload-ready files so the web tool
can be used as a cross-check.

Algorithm, as defined by IEDB:

  For an epitope of length *L* and one target protein, slide a window of length
  *L* along the target and compute the fraction of identical positions at each
  offset. The identity of the epitope against that protein is the **maximum** over
  all offsets. The *degree of conservancy* is the fraction of target proteins
  whose identity reaches the chosen threshold; at a threshold of 1.0 this is the
  fraction of strains containing the epitope exactly.

Note the difference from the alignment-based view: conservancy is computed on
ungapped sliding windows, independently of any MSA, so it is unaffected by how
the alignment in stage 7 was built.

When ``conservation.reference_sets.enabled`` is true, each strain is analysed
against the corresponding NCBI reference set in ``resources/conservancy/``
(fetched by ``setup/fetch_ncbi_conservancy.py``), producing one row per
(consensus, strain).  When disabled, the fallback is the user's own input FASTAs
from the manifest — which usually means 1 sequence per strain and 100% for
everything.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

from . import fasta as fasta_mod
from . import schemas, tables


def identity(epitope: str, target: str) -> float:
    """Maximum fractional identity of *epitope* against any window of *target*."""
    if not epitope:
        return 0.0
    length = len(epitope)
    if len(target) < length:
        # Compare against the whole short target, scoring over the epitope length
        # so that a truncated record cannot inflate identity.
        best = sum(1 for a, b in zip(epitope, target) if a == b)
        return best / length

    best = 0
    for offset in range(len(target) - length + 1):
        window = target[offset : offset + length]
        matches = 0
        for a, b in zip(epitope, window):
            if a == b:
                matches += 1
        if matches > best:
            best = matches
            if best == length:
                break
    return best / length


def analyse(
    peptides: Sequence[tuple[str, str]],
    targets: Sequence[tuple[str, str]],
    threshold: float = 1.0,
) -> list[dict[str, object]]:
    """Compute conservancy of each ``(peptide_id, sequence)`` across *targets*.

    *targets* are ``(label, protein_sequence)`` pairs.  When called with a
    per-strain reference set the labels are accession IDs; when called with the
    manifest fallback they are strain names (legacy behaviour).
    """
    rows: list[dict[str, object]] = []
    for peptide_id, peptide in peptides:
        identities: list[tuple[str, float]] = [
            (label, identity(peptide, sequence)) for label, sequence in targets
        ]
        matched = [label for label, value in identities if value >= threshold]
        values = [value for _, value in identities] or [0.0]
        fraction = len(matched) / len(targets) if targets else 0.0

        rows.append(
            {
                "consensus_id": peptide_id,
                "sequence": peptide,
                "n_target_sequences": len(targets),
                "n_matched": len(matched),
                "conservancy_fraction": fraction,
                "conservancy_percent": round(fraction * 100.0, 2),
                "min_identity": round(min(values) * 100.0, 2),
                "max_identity": round(max(values) * 100.0, 2),
                "matched_accessions": "|".join(sorted(matched)),
            }
        )
    return rows


def deduplicate_records(records: Iterable[fasta_mod.Record]) -> list[fasta_mod.Record]:
    """Collapse identical sequences so conservancy is not inflated by duplicates.

    This is the filtering step called for in the pipeline specification: if the
    same protein sequence was submitted twice under different accessions, counting
    it twice would overstate how conserved an epitope is.
    """
    seen: dict[str, fasta_mod.Record] = {}
    for record in records:
        seen.setdefault(record.md5, record)
    return list(seen.values())


def _load_manifest_targets(
    manifest_path: str,
    protein_label: str,
) -> dict[str, list[fasta_mod.Record]]:
    """Load target sequences from the manifest, grouped by strain.

    Returns ``{strain: [records]}`` so the caller can iterate strains the same
    way it would iterate reference FASTAs.  This is the fallback when NCBI
    reference sets are not configured.
    """
    by_strain: dict[str, list[fasta_mod.Record]] = {}
    for row in tables.read_table(manifest_path):
        if row["protein_label"] != protein_label:
            continue
        for record in fasta_mod.read_fasta(row["path"]):
            record.strain = row["strain"]
            by_strain.setdefault(row["strain"], []).append(record)
    return by_strain


def run(
    consensus_fasta: str,
    manifest_path: str,
    protein_label: str,
    out_table: str,
    out_targets_fasta: str,
    out_epitopes_txt: str,
    threshold: float = 1.0,
    deduplicate: bool = True,
    reference_fastas: dict[str, str] | None = None,
) -> list[dict[str, object]]:
    """Full stage-7 conservancy run for one protein.

    When *reference_fastas* is ``None`` (the default), target sequences come
    from the manifest — the user's own inputs.  When it is a dict mapping
    ``{strain_name: fasta_path}``, each strain is analysed against its own
    NCBI reference set, producing one row per ``(consensus, strain)``.

    Also writes the two files the IEDB web tool expects, so a manual cross-check
    needs no reformatting: the non-redundant target FASTA and a plain list of
    epitope sequences.
    """
    peptides = [(record.id, record.sequence) for record in fasta_mod.read_fasta(consensus_fasta)]

    # Write the epitope list (same regardless of reference mode)
    with open(out_epitopes_txt, "w", encoding="utf-8", newline="\n") as handle:
        for _, sequence in peptides:
            handle.write(sequence + "\n")

    # --- Build per-strain target sets ----------------------------------------
    if reference_fastas:
        # NCBI reference sets — one FASTA per strain, already filtered by the
        # fetch script.  Deduplication is still applied because the fetch
        # script's dedupe is by raw sequence, and a record that differs only
        # in header could appear across merged batches.
        strain_targets: dict[str, list[fasta_mod.Record]] = {}
        for strain, fasta_path in sorted(reference_fastas.items()):
            records = fasta_mod.read_fasta(fasta_path)
            for r in records:
                r.strain = strain
            if deduplicate:
                before = len(records)
                records = deduplicate_records(records)
                if before != len(records):
                    print(
                        f"[mpvpredpip] conservancy: collapsed {before} -> "
                        f"{len(records)} non-redundant {protein_label}/{strain} "
                        f"reference sequences",
                    )
            strain_targets[strain] = records
    else:
        # Fallback: the user's own manifest entries
        by_strain = _load_manifest_targets(manifest_path, protein_label)
        if deduplicate:
            for strain, records in by_strain.items():
                before = len(records)
                by_strain[strain] = deduplicate_records(records)
                after = len(by_strain[strain])
                if before != after:
                    print(
                        f"[mpvpredpip] conservancy: collapsed {before} -> "
                        f"{after} non-redundant {protein_label}/{strain} sequences",
                    )
        strain_targets = by_strain

    # --- Write the combined targets FASTA for IEDB cross-check ---------------
    all_targets: list[fasta_mod.Record] = []
    for strain in sorted(strain_targets):
        all_targets.extend(strain_targets[strain])

    fasta_mod.write_fasta(
        out_targets_fasta,
        [(f"{record.strain}|{record.accession}", record.sequence) for record in all_targets],
    )

    # --- Analyse per strain --------------------------------------------------
    all_rows: list[dict[str, object]] = []
    for strain in sorted(strain_targets):
        targets = strain_targets[strain]
        target_pairs = [(record.accession or record.id, record.sequence) for record in targets]
        rows = analyse(peptides, target_pairs, threshold=threshold)
        for row in rows:
            row["protein_label"] = protein_label
            row["strain"] = strain
        all_rows.extend(rows)

    tables.write_table(out_table, schemas.CONSERVANCY, all_rows)
    return all_rows
