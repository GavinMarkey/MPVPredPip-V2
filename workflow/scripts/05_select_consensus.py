"""Stage 5 - apply the B/T composition rule and emit consensus peptides."""

# No `from __future__ import annotations` - see 01_manifest.py.

import _common

_common.bootstrap()

from mpvpredpip import cluster as cluster_mod  # noqa: E402
from mpvpredpip import tables  # noqa: E402

smk = snakemake  # noqa: F821

with _common.tee_log(_common.log_path(smk)):
    protein = str(smk.wildcards.protein)
    require = dict(smk.params.require)
    require_both = bool(smk.params.require_both)

    total, passing = cluster_mod.select(
        cluster_path=str(smk.input.raw),
        epitope_table=str(smk.input.epitopes),
        protein_label=protein,
        out_clusters=str(smk.output.clusters),
        out_members=str(smk.output.members),
        out_fasta=str(smk.output.fasta),
        require=require,
        require_both_bcell_tools=require_both,
        reject_ambiguous=bool(smk.params.reject_ambiguous),
    )

    rows = tables.read_table(str(smk.output.clusters))
    kept = [row for row in rows if tables.as_bool(row["passes"])]

    rule_text = (
        f"at least {require.get('min_mhc_i', 1)} MHC-I, {require.get('min_mhc_ii', 1)} MHC-II"
        + (
            f", {require.get('min_bepipred', 1)} BepiPred and {require.get('min_epidope', 1)} EpiDope"
            if require_both
            else f", and {require.get('min_bepipred', 1) + require.get('min_epidope', 1)} B-cell epitopes from either tool"
        )
    )

    lines = [
        f"Clusters found: **{total}**  ",
        f"Clusters passing the composition rule: **{len(kept)}**",
        "",
        f"Rule applied: {rule_text}.",
        "",
        "| Cluster | Consensus | Len | BepiPred | EpiDope | MHC-I | MHC-II | Species | Disagree | Passes |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | :---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['cluster_id']} | `{row['consensus_sequence']}` | {row['consensus_length']} | "
            f"{row['n_bepipred']} | {row['n_epidope']} | {row['n_mhc_i']} | {row['n_mhc_ii']} | "
            f"{row['strains'].replace('|', ', ') or '-'} | {row['n_disagreeing']} | "
            f"{'yes' if tables.as_bool(row['passes']) else 'no'} |"
        )

    rejected = [row for row in rows if not tables.as_bool(row["passes"])]
    if rejected:
        lines += ["", "### Why clusters were rejected", ""]
        for row in rejected:
            lines.append(f"- **{row['cluster_id']}**: {row['reason']}")

    _common.write_summary(
        str(smk.output.summary), f"Stage 5 - {protein} cluster selection", lines
    )

    print(f"{protein}: {len(kept)}/{total} clusters passed the composition rule.")
    if not kept:
        print(
            "WARNING: no clusters passed. The most common causes are a clustering "
            "threshold that is too high to merge B- and T-cell epitopes into shared "
            "clusters, or one predictor contributing nothing. See the summary for the "
            "per-cluster reason."
        )
