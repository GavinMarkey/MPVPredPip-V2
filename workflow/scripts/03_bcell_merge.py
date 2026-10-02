"""Stage 3 - merge BepiPred and EpiDope epitopes into one per-protein table."""

# No `from __future__ import annotations` - see 01_manifest.py.

import _common

_common.bootstrap()

from mpvpredpip import fasta as fasta_mod  # noqa: E402
from mpvpredpip import schemas, tables  # noqa: E402

smk = snakemake  # noqa: F821

with _common.tee_log(_common.log_path(smk)):
    rows = tables.read_table(str(smk.input.bepipred)) + tables.read_table(str(smk.input.epidope))
    tables.write_table(str(smk.output.table), schemas.EPITOPES, rows)

    # The strain goes in the defline after a space, so the identifier stays the
    # first whitespace token and every join keyed on it still works - read_fasta
    # takes Record.id from that token. Without it this file says which protein
    # and which predictor an epitope came from but not which species, which is
    # the one thing a three-species design most needs to see.
    fasta_mod.write_fasta(
        str(smk.output.fasta),
        [(f"{row['epitope_id']} {row['strain']}", row["sequence"]) for row in rows],
        width=0,
    )

    by_tool: dict[str, list[dict]] = {}
    for row in rows:
        by_tool.setdefault(row["source_tool"], []).append(row)

    # An epitope supported by both predictors is the strongest B-cell evidence
    # available here, and stage 5 requires at least one from each, so it is worth
    # surfacing the overlap now rather than discovering it is empty later.
    sequences_by_tool = {tool: {r["sequence"] for r in items} for tool, items in by_tool.items()}
    shared = set.intersection(*sequences_by_tool.values()) if len(sequences_by_tool) > 1 else set()

    lines = [f"Total B-cell epitopes: **{len(rows)}**", ""]
    for tool, items in sorted(by_tool.items()):
        lengths = [int(item["length"]) for item in items]
        lines.append(
            f"- **{tool}**: {len(items)} epitopes, "
            f"{len(sequences_by_tool[tool])} unique peptides, "
            f"length {min(lengths) if lengths else 0}-{max(lengths) if lengths else 0}"
        )
    lines += [
        "",
        f"Peptides predicted identically by both tools: **{len(shared)}**",
    ]
    if len(by_tool) < 2:
        lines += [
            "",
            "> Only one B-cell predictor produced epitopes. Stage 5 requires one epitope",
            "> from BepiPred **and** one from EpiDope by default, so no cluster can pass.",
            "> Check the other predictor's log, or relax",
            "> `clustering.require_both_bcell_tools` in `config.yaml`.",
        ]

    _common.write_summary(str(smk.output.summary), "Stage 3 - B-cell epitopes", lines)
    print(f"Merged {len(rows)} B-cell epitopes from {len(by_tool)} predictor(s).")
