"""Stage 9 - assemble the final candidate list and the run report.

Joins the per-protein tables into one ranked list and records how many candidates
were lost at each stage. The attrition table is the fastest way to tell a real
biological result from a mis-set threshold: a stage that removes everything is
almost always a configuration problem.

No ``from __future__ import annotations`` - see 01_manifest.py.
"""

import _common

_common.bootstrap()

from mpvpredpip import fasta as fasta_mod  # noqa: E402
from mpvpredpip import schemas, tables  # noqa: E402

smk = snakemake  # noqa: F821


def load_all(paths) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        rows.extend(tables.read_table(str(path)))
    return rows


with _common.tee_log(_common.log_path(smk)):
    proteins = list(smk.params.proteins)

    manifest = tables.read_table(str(smk.input.manifest))
    bcell = load_all(smk.input.bcell)
    tcell = load_all(smk.input.tcell)
    clusters = load_all(smk.input.clusters)
    members = load_all(smk.input.members)
    screening = load_all(smk.input.screening)
    # Conservancy now produces one row per (consensus, strain).  Group by
    # consensus and report the minimum across strains — the weakest strain
    # bounds a multi-strain claim (100% in Zaire and 40% in Sudan is a 40%
    # candidate, not a 100% one).  Per-strain values are stored alongside.
    conservancy_rows = load_all(smk.input.conservancy)
    conservancy_by_id: dict[str, list[dict]] = {}
    for row in conservancy_rows:
        conservancy_by_id.setdefault(row["consensus_id"], []).append(row)

    # World coverage is each candidate's OWN figure. Stage 8 submits every
    # peptide to the IEDB tool separately for exactly this reason: the set
    # figure describes what the peptides achieve together and is identical for
    # all of them, so ranking on it compared nothing. The set figure is still
    # in the stage 8 table under consensus_id "ALL", and is reported below.
    #
    # "combined" is the class I + class II calculation - see
    # population_coverage.mhc_class in config.yaml. Rows for class I and class
    # II alone are in the stage 8 table too.
    coverage_rows = load_all(smk.input.coverage)
    world_coverage: dict[str, str] = {}       # consensus_id -> its own coverage
    set_world_coverage: dict[str, str] = {}   # protein_label -> the "ALL" row
    alleles_by_consensus: dict[str, str] = {}
    for row in coverage_rows:
        world_combined = (
            row["population"].strip().lower() == "world"
            and str(row["mhc_class"]).strip().lower() == "combined"
        )
        if row["consensus_id"] == "ALL":
            if world_combined:
                set_world_coverage[row["protein_label"]] = row["coverage_percent"]
        else:
            alleles_by_consensus[row["consensus_id"]] = row["alleles"]
            if world_combined:
                world_coverage[row["consensus_id"]] = row["coverage_percent"]

    cluster_by_id = {row["consensus_id"]: row for row in clusters}

    # Split each consensus peptide's alleles back into class I and class II using
    # the epitope class recorded in stage 5 membership.
    class_by_epitope = {row["epitope_id"]: row["epitope_class"] for row in (bcell + tcell)}
    mhc_i_by_consensus: dict[str, set[str]] = {}
    mhc_ii_by_consensus: dict[str, set[str]] = {}
    allele_by_epitope: dict[str, set[str]] = {}
    for row in bcell + tcell:
        if row.get("allele"):
            allele_by_epitope.setdefault(row["epitope_id"], set()).add(row["allele"])
    for row in members:
        alleles = allele_by_epitope.get(row["epitope_id"], set())
        if not alleles:
            continue
        target = (
            mhc_i_by_consensus
            if class_by_epitope.get(row["epitope_id"]) == "MHC-I"
            else mhc_ii_by_consensus
        )
        target.setdefault(row["consensus_id"], set()).update(alleles)

    passing = [row for row in screening if tables.as_bool(row["passes"])]

    final_rows: list[dict] = []
    for row in passing:
        consensus_id = row["consensus_id"]
        cluster = cluster_by_id.get(consensus_id, {})

        # Aggregate conservancy: minimum across strains, with per-strain detail.
        cons_rows = conservancy_by_id.get(consensus_id, [])
        if cons_rows:
            min_cons = min(tables.as_float(r.get("conservancy_percent", ""), 0.0) for r in cons_rows)
            per_strain = "|".join(
                f"{r.get('strain', '?')}:{r.get('conservancy_percent', '?')}%"
                for r in sorted(cons_rows, key=lambda r: str(r.get("strain", "")))
            )
        else:
            min_cons = ""
            per_strain = ""

        final_rows.append(
            {
                "consensus_id": consensus_id,
                "protein_label": row["protein_label"],
                "cluster_id": row.get("cluster_id", ""),
                "sequence": row["sequence"],
                "length": row["length"],
                # Prefer stage 6's copy, falling back to the cluster row: both
                # carry it, and a candidate must never reach the final table
                # without saying which species it covers.
                "strains": row.get("strains", "") or cluster.get("strains", ""),
                "n_strains": row.get("n_strains", "") or cluster.get("n_strains", ""),
                "antigenicity_score": row.get("antigenicity_score", ""),
                "allergen_score": row.get("allergen_score", ""),
                "allergen": tables.as_bool(row.get("allergen", False)),
                "allergen_identity_percent": row.get("allergen_identity_percent", ""),
                "allergen_homology_risk": tables.as_bool(
                    row.get("allergen_homology_risk", False)
                ),
                "human_match_count": row.get("human_match_count", 0),
                "conservancy_percent": min_cons,
                "conservancy_per_strain": per_strain,
                "world_coverage_percent": world_coverage.get(consensus_id, ""),
                "n_bepipred": cluster.get("n_bepipred", ""),
                "n_epidope": cluster.get("n_epidope", ""),
                "n_mhc_i": cluster.get("n_mhc_i", ""),
                "n_mhc_ii": cluster.get("n_mhc_ii", ""),
                "mhc_i_alleles": "|".join(sorted(mhc_i_by_consensus.get(consensus_id, set()))),
                "mhc_ii_alleles": "|".join(sorted(mhc_ii_by_consensus.get(consensus_id, set()))),
            }
        )

    rank_by = str(smk.params.rank_by)
    final_rows.sort(key=lambda r: tables.as_float(r.get(rank_by, ""), float("-inf")), reverse=True)
    for index, row in enumerate(final_rows, start=1):
        row["rank"] = index

    tables.write_table(str(smk.output.table), schemas.FINAL, final_rows)
    fasta_mod.write_fasta(
        str(smk.output.fasta),
        [(row["consensus_id"], row["sequence"]) for row in final_rows],
        width=0,
    )

    # --- attrition -----------------------------------------------------------
    attrition: list[dict] = []
    for protein in proteins:
        n_strains = len({r["strain"] for r in manifest if r["protein_label"] == protein})
        n_bcell = len({r["sequence"] for r in bcell if r["protein_label"] == protein})
        n_tcell = len({r["sequence"] for r in tcell if r["protein_label"] == protein})
        protein_clusters = [r for r in clusters if r["protein_label"] == protein]
        n_pass_cluster = sum(1 for r in protein_clusters if tables.as_bool(r["passes"]))
        protein_screened = [r for r in screening if r["protein_label"] == protein]
        n_pass_screen = sum(1 for r in protein_screened if tables.as_bool(r["passes"]))

        attrition.append(
            {
                "protein": protein,
                "strains": n_strains,
                "bcell_peptides": n_bcell,
                "tcell_peptides": n_tcell,
                "clusters": len(protein_clusters),
                "clusters_passing": n_pass_cluster,
                "screened": len(protein_screened),
                "final_candidates": n_pass_screen,
            }
        )

    tables.write_table(
        str(smk.output.attrition),
        ["protein", "strains", "bcell_peptides", "tcell_peptides", "clusters",
         "clusters_passing", "screened", "final_candidates"],
        attrition,
    )

    # --- report --------------------------------------------------------------
    stage_dirs = dict(smk.params.stage_dirs)
    # Said from the table rather than from the config: allergenicity is advisory
    # by default, so a candidate can reach this list having been CALLED an
    # allergen. The old wording asserted every candidate was "non-allergenic",
    # which silently became false the moment that screen stopped dropping
    # peptides - the most load-bearing sentence in the report, wrong, with
    # nothing to flag it. Counting the rows cannot drift from what is in them.
    flagged = [row for row in final_rows if row["allergen"]]
    lines = [
        f"**{len(final_rows)}** consensus epitopes passed every stage of the pipeline.",
        "",
        "Each of them contains, within one cluster, at least one MHC class I epitope, one",
        "MHC class II epitope and B-cell epitopes from both BepiPred and EpiDope; is",
        "predicted antigenic; and has no exact match to the human proteome.",
        "",
    ]
    if flagged:
        lines.extend([
            f"**{len(flagged)} of the {len(final_rows)} were predicted allergenic** "
            "by AlgPred and are on this list anyway, because allergenicity is "
            "advisory rather than a hard filter at this peptide length - see "
            "`screening.allergenicity.drop_on_allergen` in config.yaml and the "
            "caveats in each protein's stage 6 summary. The `allergen` and "
            "`allergen_score` columns below carry the call, and these candidates "
            "need an allergenicity assessment that does not rest on amino-acid "
            "composition alone before being taken forward: "
            + ", ".join(f"`{row['consensus_id']}`" for row in flagged) + ".",
            "",
        ])
    lines.extend([
        "## Attrition by stage",
        "",
        "| Protein | Strains | B-cell peptides | T-cell peptides | Clusters | Passed stage 5 | Passed stage 6 |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for row in attrition:
        lines.append(
            f"| {row['protein']} | {row['strains']} | {row['bcell_peptides']} | "
            f"{row['tcell_peptides']} | {row['clusters']} | {row['clusters_passing']} | "
            f"{row['final_candidates']} |"
        )

    lines += ["", "## Final candidates", ""]
    if final_rows:
        lines += [
            "`World coverage` is each candidate's own, class I and II combined.",
            "",
            "| Rank | ID | Protein | Sequence | Len | Antigenicity | Min conservancy | Per strain | World coverage |",
            "| ---: | --- | --- | --- | ---: | ---: | ---: | --- | ---: |",
        ]
        for row in final_rows:
            cons_pct = row["conservancy_percent"]
            cons_str = f"{cons_pct}%" if cons_pct != "" else "-"
            per_strain = _common.md_cell(str(row.get("conservancy_per_strain", ""))) or "-"
            lines.append(
                f"| {row['rank']} | `{_common.md_cell(str(row['consensus_id']))}` "
                f"| {row['protein_label']} | "
                f"`{row['sequence']}` | {row['length']} | {row['antigenicity_score']} | "
                f"{cons_str} | {per_strain} | {row['world_coverage_percent']}% |"
            )
        if set_world_coverage:
            lines += [
                "",
                "Taken together, per protein (the `ALL` rows of the stage 8 table) - "
                "the fraction of the World population responding to at least one of "
                "that protein's candidates:",
                "",
            ]
            for protein_label in sorted(set_world_coverage):
                lines.append(f"- **{protein_label}**: {set_world_coverage[protein_label]}%")
    else:
        lines += [
            "No candidates survived every filter. Work back through the per-stage",
            "summaries to find where the list emptied:",
            "",
            f"- `{stage_dirs['S3']}/<protein>/summary.md` - B-cell epitope counts per tool",
            f"- `{stage_dirs['S4']}/<protein>/summary.md` - T-cell binders per MHC class",
            f"- `{stage_dirs['S5']}/<protein>/summary.md` - per-cluster rejection reasons",
            f"- `{stage_dirs['S6']}/<protein>/summary.md` - which screen removed each peptide",
            "",
            "A stage that removes everything usually means a threshold in `config.yaml`",
            "needs adjusting, not that there are no candidates.",
        ]

    lines += [
        "",
        "## Where the outputs are",
        "",
        f"- `{stage_dirs['S1']}/manifest.tsv` - inputs, proteins and strains",
        f"- `{stage_dirs['S2']}/whole_protein_antigenicity.tsv` - whole-protein screen",
        f"- `{stage_dirs['S3']}/<protein>/bcell_epitopes.tsv` - BepiPred + EpiDope",
        f"- `{stage_dirs['S4']}/<protein>/tcell_epitopes.tsv` - MHC-I (immunogenic) + MHC-II",
        f"- `{stage_dirs['S5']}/<protein>/clusters.tsv` - clusters and why each passed or failed",
        f"- `{stage_dirs['S6']}/<protein>/screening.tsv` - antigenicity, allergenicity, autoimmunity",
        f"- `{stage_dirs['S7']}/<protein>/conservancy.tsv` - conservancy across strains",
        f"- `{stage_dirs['S8']}/<protein>/population_coverage.tsv` - IEDB population coverage",
        f"- `{stage_dirs['S9']}/final_candidates.tsv` - this list, ranked by {rank_by}",
        "",
        "Each stage also keeps the untouched native output of the tool that produced it,",
        "under that stage's `raw/` directory.",
    ]

    _common.write_summary(str(smk.output.report), "Vaccine candidate pipeline - run report", lines)

    print(f"Final candidates: {len(final_rows)}")
    for row in final_rows:
        print(f"  {row['rank']}. {row['consensus_id']}  {row['sequence']}")
