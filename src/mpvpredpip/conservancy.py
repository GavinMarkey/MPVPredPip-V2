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

    *targets* are ``(strain_label, protein_sequence)`` pairs.
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
                "matched_strains": "|".join(sorted(matched)),
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


def run(
    consensus_fasta: str,
    manifest_path: str,
    protein_label: str,
    out_table: str,
    out_targets_fasta: str,
    out_epitopes_txt: str,
    threshold: float = 1.0,
    deduplicate: bool = True,
) -> list[dict[str, object]]:
    """Full stage-7 conservancy run for one protein.

    Also writes the two files the IEDB web tool expects, so a manual cross-check
    needs no reformatting: the non-redundant target FASTA and a plain list of
    epitope sequences.
    """
    peptides = [(record.id, record.sequence) for record in fasta_mod.read_fasta(consensus_fasta)]

    targets: list[fasta_mod.Record] = []
    for row in tables.read_table(manifest_path):
        if row["protein_label"] != protein_label:
            continue
        for record in fasta_mod.read_fasta(row["path"]):
            record.strain = row["strain"]
            targets.append(record)

    if deduplicate:
        before = len(targets)
        targets = deduplicate_records(targets)
        if before != len(targets):
            print(
                f"[mpvpredpip] conservancy: collapsed {before} -> {len(targets)} "
                f"non-redundant {protein_label} sequences",
            )

    fasta_mod.write_fasta(
        out_targets_fasta,
        [(f"{record.strain}|{record.accession}", record.sequence) for record in targets],
    )
    with open(out_epitopes_txt, "w", encoding="utf-8", newline="\n") as handle:
        for _, sequence in peptides:
            handle.write(sequence + "\n")

    rows = analyse(
        peptides,
        [(record.strain, record.sequence) for record in targets],
        threshold=threshold,
    )
    for row in rows:
        row["protein_label"] = protein_label

    tables.write_table(out_table, schemas.CONSERVANCY, rows)
    return rows
