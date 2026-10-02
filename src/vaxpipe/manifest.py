"""Stage 1 - build the sequence manifest.

Scans the input directory, parses protein and strain from each record, validates
the sequences and writes ``manifest.tsv`` plus one normalised per-protein FASTA
that every later stage consumes. Doing the parsing once, here, is what lets the
rest of the pipeline treat "protein" and "strain" as first-class wildcards.
"""

from __future__ import annotations

import os
import re
import sys
from collections import defaultdict

from . import fasta, inputs, schemas, tables


def collect(input_dir: str, extensions: list[str], aliases: dict[str, str]) -> list[fasta.Record]:
    """Read every FASTA under *input_dir* and annotate each record."""
    records: list[fasta.Record] = []
    if not os.path.isdir(input_dir):
        raise SystemExit(
            f"Input directory does not exist: {input_dir}\n\n"
            + inputs.short_rules(input_dir, extensions)
        )

    for path in inputs.find_inputs(input_dir, extensions):
        for record in fasta.read_fasta(path):
            record.protein_label = fasta.normalise_protein_label(record.name, aliases)
            records.append(record)

    if not records:
        raise SystemExit(
            f"No FASTA files with extensions {sorted(e.lower() for e in extensions)} "
            f"found in {input_dir!r}.\n\n" + inputs.short_rules(input_dir, extensions)
        )
    return records


def build(
    input_dir: str,
    out_manifest: str,
    out_dir: str,
    extensions: list[str],
    aliases: dict[str, str],
    min_length: int = 50,
) -> list[fasta.Record]:
    """Validate the inputs and write the manifest plus per-protein FASTAs."""
    records = collect(input_dir, extensions, aliases)

    # Shared with the preflight report (setup/init_input.py), so what that
    # check calls acceptable is exactly what this stage accepts.
    problems = inputs.find_problems(records, min_length=min_length)

    if problems:
        raise SystemExit(
            "Stage 1 input validation failed:\n" + "\n".join(problems) +
            "\n\nFix the inputs above and re-run. Nothing downstream has been touched."
        )

    tables.write_table(
        out_manifest,
        schemas.MANIFEST,
        [
            {
                "accession": record.accession,
                "protein_label": record.protein_label,
                "protein_name": record.name,
                "strain": record.strain,
                "length": len(record),
                "md5": record.md5,
                "path": record.source_path.replace("\\", "/"),
            }
            for record in records
        ],
    )

    # One FASTA per protein, holding every strain, with stable deflines of the
    # form ">Protein|Strain|Accession" that later stages can split on.
    by_protein: dict[str, list[fasta.Record]] = defaultdict(list)
    for record in records:
        by_protein[record.protein_label].append(record)

    for protein, group in by_protein.items():
        fasta.write_fasta(
            os.path.join(out_dir, f"{protein}.fasta"),
            [(make_header(r), r.sequence) for r in sorted(group, key=lambda r: r.strain)],
        )

    # Per protein *and* strain, because the predictors run one sequence at a time
    # and their outputs must stay attributable to a single strain.
    for record in records:
        fasta.write_fasta(
            os.path.join(out_dir, "by_strain", f"{record.protein_label}__{slug(record.strain)}.fasta"),
            [(make_header(record), record.sequence)],
        )

    print(
        f"Stage 1: {len(records)} sequences across {len(by_protein)} proteins "
        f"({', '.join(sorted(by_protein))})",
        file=sys.stderr,
    )
    return records


def _no_spaces(value: str) -> str:
    """Collapse internal whitespace to underscores.

    A record's identifier is the first whitespace-delimited token of its defline
    (:attr:`fasta.Record.id`). Strain names routinely contain a space - "Zaire
    ebolavirus", "Bundibugyo virus" - so without this the defline is severed at
    that space: the identifier becomes "Nucleoprotein|Zaire" and the accession is
    dropped from every downstream table.
    """
    return re.sub(r"\s+", "_", value.strip())


def make_header(record: fasta.Record) -> str:
    """Stable defline: ``Protein|Strain|Accession``.

    Each field has its whitespace collapsed to underscores so the whole defline
    survives as a single token; see :func:`_no_spaces`. The manifest keeps the
    strain's original spelling, so nothing is lost - this affects identifiers
    only.
    """
    return "|".join(
        _no_spaces(part)
        for part in (record.protein_label, record.strain, record.accession)
    )


def parse_header(header: str) -> tuple[str, str, str]:
    """Inverse of :func:`make_header`. Tolerates extra description text.

    Strain comes back underscored, as it is written; ``manifest.tsv`` holds the
    original spelling for display.
    """
    first = header.split()[0] if header else ""
    parts = first.split("|")
    while len(parts) < 3:
        parts.append("")
    return parts[0], parts[1], parts[2]


def slug(text: str) -> str:
    """Filesystem-safe token. Used for strain names, which contain spaces."""
    keep = [c if c.isalnum() else "_" for c in text.strip()]
    out = "".join(keep)
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_") or "unknown"


def proteins(manifest_path: str) -> list[str]:
    """Distinct protein labels in the manifest, sorted."""
    return sorted({row["protein_label"] for row in tables.read_table(manifest_path)})


def strains(manifest_path: str, protein: str | None = None) -> list[str]:
    rows = tables.read_table(manifest_path)
    if protein is not None:
        rows = [row for row in rows if row["protein_label"] == protein]
    return sorted({row["strain"] for row in rows})


def pairs(manifest_path: str) -> list[tuple[str, str]]:
    """``(protein_label, strain)`` pairs, which are the wildcards of stages 2-4."""
    return sorted(
        {(row["protein_label"], row["strain"]) for row in tables.read_table(manifest_path)}
    )
