"""Stage 5 - sweep the identity thresholds, apply the B/T rule, pool the result.

Clustering runs once per threshold in ``clustering.threshold`` (``"90:50"`` ->
90, 80, 70, 60, 50). Each is a complete, independent clustering of the same
epitopes: a tighter threshold merges less, so its consensus is closer to a real
peptide, while a looser one is likelier to gather all four epitope types into
one cluster and so to pass the composition rule at all.

Everything that passes, at any threshold, is pooled and screened. Peptides are
deduplicated by sequence first, keeping the tightest threshold that produced
them, because a peptide passing at 90% usually passes again at every looser
threshold below it.

Clusters whose consensus contains an X are discarded outright rather than
listed as failures - see select_at_threshold in mpvpredpip.cluster.
"""

# No `from __future__ import annotations` - see 01_manifest.py.

import _common

_common.bootstrap()

from mpvpredpip import cluster as cluster_mod  # noqa: E402
from mpvpredpip import fasta as fasta_mod  # noqa: E402
from mpvpredpip import schemas, tables  # noqa: E402

smk = snakemake  # noqa: F821

#: Rejection reasons carry counts ("only 0 MHC-II epitope(s), needs 1"), so the
#: raw strings would never group. The leading phrase is what distinguishes them.
_REASON_WIDTH = 60

with _common.tee_log(_common.log_path(smk)):
    protein = str(smk.wildcards.protein)
    require = dict(smk.params.require)
    require_both = bool(smk.params.require_both)
    reject_ambiguous = bool(smk.params.reject_ambiguous)
    thresholds = [int(t) for t in smk.params.thresholds]

    raw_paths = list(smk.input.raw)

    # expand() preserves the order of the list it is given, so these line up
    # with `thresholds`. Checked rather than assumed: a mismatch would file one
    # threshold's clusters under another's, which nothing downstream could
    # detect - every pooled row carries the threshold it was labelled with, so
    # a mislabelling would look like a legitimate result.
    if len(raw_paths) != len(thresholds):
        raise SystemExit(
            f"Stage 5 was given {len(thresholds)} thresholds but "
            f"{len(raw_paths)} cluster files. These come from the same expand() "
            "over CLUSTER_THRESHOLDS and must match."
        )

    epitope_rows = tables.read_table(str(smk.input.epitopes))

    results = []
    for threshold, raw_path in zip(thresholds, raw_paths):
        result = cluster_mod.select_at_threshold(
            cluster_path=str(raw_path),
            epitope_rows=epitope_rows,
            protein_label=protein,
            threshold=threshold,
            require=require,
            require_both_bcell_tools=require_both,
            reject_ambiguous=reject_ambiguous,
        )
        results.append(result)
        print(
            f"{protein} @ {threshold}%: {result.n_found} clusters, "
            f"{result.n_ambiguous} dropped for an ambiguous consensus, "
            f"{result.n_passing} passed."
        )

    # The two tables stage 5 produces. Every row carries the threshold it came
    # from, which is what the per-threshold files used to provide by splitting
    # the same rows across ten files.
    tables.write_table(
        str(smk.output.clusters),
        schemas.CLUSTERS,
        [row for result in results for row in result.cluster_rows],
    )
    tables.write_table(
        str(smk.output.members),
        schemas.CLUSTER_MEMBERS,
        [row for result in results for row in result.member_rows],
    )

    entries, duplicates = cluster_mod.pool_passing(results)
    fasta_mod.write_fasta(str(smk.output.fasta), entries, width=0)

    # --- summary -------------------------------------------------------------
    rule_text = (
        f"at least {require.get('min_mhc_i', 1)} MHC-I, {require.get('min_mhc_ii', 1)} MHC-II"
        + (
            f", {require.get('min_bepipred', 1)} BepiPred and {require.get('min_epidope', 1)} EpiDope"
            if require_both
            else f", and {require.get('min_bepipred', 1) + require.get('min_epidope', 1)} B-cell epitopes from either tool"
        )
        + (
            f", spanning at least {require['min_strains']} strain(s)"
            if require.get("min_strains", 0)
            else ""
        )
    )

    total_found = sum(result.n_found for result in results)
    total_ambiguous = sum(result.n_ambiguous for result in results)
    total_passing = sum(result.n_passing for result in results)

    lines = [
        f"Thresholds swept: **{', '.join(str(t) + '%' for t in thresholds)}**  ",
        f"Clusters found across the sweep: **{total_found}**  ",
        f"Discarded for an ambiguous (X) consensus: **{total_ambiguous}**  ",
        f"Passing the composition rule: **{total_passing}**  ",
        f"Unique peptides sent to stage 6: **{len(entries)}** "
        f"({duplicates} duplicate(s) collapsed)",
        "",
        f"Rule applied: {rule_text}.",
        "",
        "Each threshold is an independent clustering of the same epitopes, so a",
        "cluster numbered 1.1 at 90% is unrelated to 1.1 at 80%. Everything that",
        "passes at any threshold is screened; where one peptide passed at several,",
        "the tightest is the one kept.",
        "",
        "A consensus containing an X is not a synthesisable peptide and nothing",
        "downstream can score it, so those clusters are discarded here rather than",
        "listed below as failures.",
        "",
    ]

    passing_rows = [
        row
        for result in results
        for row in result.cluster_rows
        if tables.as_bool(row["passes"])
    ]
    kept_ids = {consensus_id for consensus_id, _ in entries}

    lines += ["## Passing clusters", ""]
    if passing_rows:
        lines += [
            "`Kept` marks the row whose peptide went forward; a `dup` passed at a",
            "tighter threshold too and was collapsed into it.",
            "",
            "| Threshold | Cluster | Consensus | Len | BepiPred | EpiDope | MHC-I | MHC-II | Species | Dropped | Kept |",
            "| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | :---: |",
        ]
        for row in passing_rows:
            lines.append(
                f"| {row['threshold']}% | {row['cluster_id']} | `{row['consensus_sequence']}` | "
                f"{row['consensus_length']} | {row['n_bepipred']} | {row['n_epidope']} | "
                f"{row['n_mhc_i']} | {row['n_mhc_ii']} | "
                f"{row['strains'].replace('|', ', ') or '-'} | {row['n_members_dropped']} | "
                f"{'yes' if row['consensus_id'] in kept_ids else 'dup'} |"
            )
    else:
        lines.append("None. The per-threshold breakdown below shows what stopped them.")
    lines.append("")

    lines += [
        "## Every threshold",
        "",
        "| Threshold | Found | Ambiguous | Examined | Passed |",
        "| ---: | ---: | ---: | ---: | ---: |",
    ]
    for result in results:
        lines.append(
            f"| {result.threshold}% | {result.n_found} | {result.n_ambiguous} | "
            f"{result.n_examined} | {result.n_passing} |"
        )
    lines.append("")

    lines += [
        "### Why clusters were rejected",
        "",
        "Counted per threshold. A cluster failing several criteria is counted",
        "under each, so these total more than the rejected count above.",
        "",
    ]
    for result in results:
        rejected = result.n_examined - result.n_passing
        lines.append(f"**{result.threshold}%** - {rejected} of {result.n_examined} rejected")
        if result.reasons:
            for reason, count in result.reasons.most_common():
                text = reason if len(reason) <= _REASON_WIDTH else reason[:_REASON_WIDTH] + "..."
                lines.append(f"- {text}: {count}")
        else:
            lines.append("- (none)")
        lines.append("")

    lines += [
        "Full per-cluster detail, including the clusters summarised above, is in",
        "`clusters.tsv`, and each cluster's constituent epitopes are in",
        "`cluster_members.tsv`. Both carry every threshold in the sweep, tagged",
        "by a `threshold` column - filter on it to read one at a time, e.g.",
        "`awk -F'\\t' 'NR==1 || $2==90' clusters.tsv`.",
    ]

    _common.write_summary(
        str(smk.output.summary), f"Stage 5 - {protein} cluster selection", lines
    )

    print(
        f"{protein}: {total_passing} passing cluster(s) across {len(thresholds)} "
        f"threshold(s) -> {len(entries)} unique peptide(s) for stage 6."
    )
    if not entries:
        print(
            "WARNING: no clusters passed at any threshold. The most common causes "
            "are one predictor contributing nothing, or a composition rule no "
            "cluster can satisfy. The per-threshold breakdown in the summary names "
            "what each rejection was for."
        )
