"""Stage 4 - merge immunogenic class I binders with class II binders."""

# No `from __future__ import annotations` - see 01_manifest.py.

import _common

_common.bootstrap()

from mpvpredpip import fasta as fasta_mod  # noqa: E402
from mpvpredpip import schemas, tables  # noqa: E402

smk = snakemake  # noqa: F821

with _common.tee_log(_common.log_path(smk)):
    rows = tables.read_table(str(smk.input.mhc_i)) + tables.read_table(str(smk.input.mhc_ii))
    tables.write_table(str(smk.output.table), schemas.EPITOPES, rows)

    # One record per unique peptide: the same peptide binds many alleles, and
    # writing it once per allele would distort the clustering input.
    #
    # The strain is aggregated over every row for that peptide, not taken from
    # whichever row came first. An identical 9-mer predicted in all three species
    # is one clustering input but three species' worth of evidence, and naming
    # only the first would understate the conservation that makes it a good
    # candidate. It goes after a space so the identifier remains the first
    # whitespace token that read_fasta uses as Record.id.
    strains_by_sequence: dict[str, set[str]] = {}
    first_id: dict[str, str] = {}
    order: list[str] = []
    for row in rows:
        sequence = str(row["sequence"])
        if sequence not in strains_by_sequence:
            strains_by_sequence[sequence] = set()
            first_id[sequence] = str(row["epitope_id"])
            order.append(sequence)
        strain = str(row.get("strain", "") or "").strip()
        if strain:
            strains_by_sequence[sequence].add(strain)

    unique_records = [
        (f"{first_id[s]} {'|'.join(sorted(strains_by_sequence[s]))}", s) for s in order
    ]
    fasta_mod.write_fasta(str(smk.output.fasta), unique_records, width=0)

    by_tool: dict[str, list[dict]] = {}
    for row in rows:
        by_tool.setdefault(row["source_tool"], []).append(row)

    lines = [f"Total T-cell epitope rows (peptide x allele): **{len(rows)}**", ""]
    for tool, items in sorted(by_tool.items()):
        peptides = {item["sequence"] for item in items}
        alleles = {item["allele"] for item in items if item["allele"]}
        label = "MHC class I" if tool == "mhc_i" else "MHC class II"
        lines.append(
            f"- **{label}**: {len(peptides)} unique peptides across {len(alleles)} alleles "
            f"({len(items)} peptide-allele pairs)"
        )

    if "mhc_i" not in by_tool or "mhc_ii" not in by_tool:
        lines += [
            "",
            "> One of the two MHC classes produced nothing. Stage 5 requires at least one",
            "> epitope of each class per cluster, so no candidate can pass. Check the",
            "> relevant log before continuing.",
        ]

    _common.write_summary(str(smk.output.summary), "Stage 4 - T-cell epitopes", lines)
    print(f"Merged {len(rows)} T-cell epitope rows; {len(order)} unique peptides.")
