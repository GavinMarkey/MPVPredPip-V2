"""Stage 5 - assemble and deduplicate the clustering input.

Pools the B-cell and T-cell epitopes for one protein. Identical peptides from the
*same* predictor collapse into one clustering entry; identical peptides from
*different* predictors are deliberately kept apart, because the whole point of
the stage 5 rule is to see which predictors support a given cluster.

No ``from __future__ import annotations`` - see 01_manifest.py.
"""

import _common

_common.bootstrap()

from mpvpredpip import epitopes as epitope_mod  # noqa: E402
from mpvpredpip import fasta as fasta_mod  # noqa: E402
from mpvpredpip import schemas, tables  # noqa: E402

smk = snakemake  # noqa: F821

with _common.tee_log(_common.log_path(smk)):
    rows = tables.read_table(str(smk.input.bcell)) + tables.read_table(str(smk.input.tcell))
    if not rows:
        raise SystemExit(
            "No epitopes reached stage 5. Check the stage 3 and stage 4 summaries - "
            "the thresholds in config.yaml are the usual cause."
        )

    tables.write_table(str(smk.output.combined), schemas.EPITOPES, rows)

    unique = epitope_mod.deduplicate(rows)
    tables.write_table(str(smk.output.table), schemas.EPITOPES_NR, unique)

    # deduplicate() has already aggregated the species into "strains", so put it
    # in the defline. The identifier stays the first whitespace token, which is
    # what read_fasta uses as Record.id and what stage 5 clusters and joins on.
    fasta_mod.write_fasta(
        str(smk.output.fasta),
        [(f"{row['epitope_id']} {row['strains']}", row["sequence"]) for row in unique],
        width=0,
    )

    counts: dict[str, int] = {}
    for row in unique:
        counts[str(row["source_tool"])] = counts.get(str(row["source_tool"]), 0) + 1

    print(f"Clustering input: {len(unique)} unique peptides from {len(rows)} epitope rows.")
    for tool in sorted(counts):
        print(f"  {tool}: {counts[tool]}")

    missing = [tool for tool in ("bepipred", "epidope", "mhc_i", "mhc_ii") if tool not in counts]
    if missing:
        print(
            f"WARNING: no peptides from {', '.join(missing)}. With the default "
            "composition rule, no cluster can pass stage 5."
        )
