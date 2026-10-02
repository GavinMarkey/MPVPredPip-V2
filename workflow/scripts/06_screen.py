"""Stage 6 - antigenicity, allergenicity and autoimmunity screening.

Hard filters applied to the consensus peptides that passed stage 5. A peptide is
carried into stage 9 only if it is predicted antigenic, is not predicted
allergenic, and does not share an exact window with the human proteome.

Every peptide keeps a row in the output table whether it passed or not, with the
reason recorded, so the candidate list can always be traced back.

No ``from __future__ import annotations`` - see 01_manifest.py.
"""

import os

import _common

_common.bootstrap()

from mpvpredpip import autoimmunity  # noqa: E402
from mpvpredpip import fasta as fasta_mod  # noqa: E402
from mpvpredpip import schemas, tables  # noqa: E402
from mpvpredpip.adapters import screening  # noqa: E402

smk = snakemake  # noqa: F821

with _common.tee_log(_common.log_path(smk)):
    protein = str(smk.wildcards.protein)
    raw_dir = str(smk.params.raw_dir)
    os.makedirs(raw_dir, exist_ok=True)

    antigenicity_cfg = dict(smk.params.antigenicity)
    allergenicity_cfg = dict(smk.params.allergenicity)
    autoimmunity_cfg = dict(smk.params.autoimmunity)

    records = fasta_mod.read_fasta(str(smk.input.fasta))
    cluster_rows = {row["consensus_id"]: row for row in tables.read_table(str(smk.input.clusters))}

    if not records:
        print(f"{protein}: no consensus peptides reached stage 6.")
        tables.write_table(str(smk.output.table), schemas.SCREENING, [])
        fasta_mod.write_fasta(str(smk.output.passed), [], width=0)
        _common.write_summary(
            str(smk.output.summary),
            f"Stage 6 - {protein} screening",
            ["No consensus peptides passed stage 5, so there was nothing to screen."],
        )
        raise SystemExit(0)

    # --- antigenicity --------------------------------------------------------
    antigenicity: dict[str, tuple[float, bool]] = {}
    if antigenicity_cfg.get("enabled", True):
        rows = screening.predict_antigenicity(
            fasta_path=str(smk.input.fasta),
            out_path=os.path.join(raw_dir, "antigenicity.tsv"),
            threshold=float(antigenicity_cfg["threshold"]),
            tool=str(smk.params.antigenicity_tool),
            command=str(smk.params.antigenicity_command or ""),
            tool_dir=str(smk.params.antigenicity_dir),
        )
        tables.write_table(os.path.join(raw_dir, "antigenicity.tsv"), schemas.ANTIGENICITY, rows)
        # A blank score means the predictor was reached but declined to score
        # this peptide (IApred refuses sequences under 20 aa, for instance). It
        # is carried as None so the reason can be reported as itself rather than
        # as a score below the threshold.
        antigenicity = {
            str(row["sequence_id"]): (
                float(row["antigenicity_score"]) if row["antigenicity_score"] != "" else None,
                bool(row["antigenic"]),
                str(row.get("note", "")),
            )
            for row in rows
        }

    # --- allergenicity -------------------------------------------------------
    allergenicity: dict[str, tuple[float, bool]] = {}
    inert_components: list[str] = []
    if allergenicity_cfg.get("enabled", True):
        # AlgPred already ran, in its own rule under its own conda environment
        # (see workflow/envs/algpred.yaml). Only its CSV is read here.
        allergenicity = screening.parse_allergenicity(
            out_path=str(smk.input.allergenicity),
            fasta_path=str(smk.input.fasta),
            threshold=float(allergenicity_cfg["threshold"]),
        )
        # Which components of the hybrid actually contributed. Reported because
        # a degenerate hybrid is graded against a threshold characterised for
        # the full one - see config.yaml.
        inert_components = screening.inert_hybrid_components(
            str(smk.input.allergenicity)
        )
        if inert_components:
            print(
                "NOTE allergenicity: " + " and ".join(inert_components) +
                " were zero for every peptide, so the hybrid score equals the "
                "ML score and each call rests on composition alone."
            )

    # --- autoimmunity --------------------------------------------------------
    window = int(autoimmunity_cfg.get("min_exact_match_length", 9))
    human_hits: dict[str, list[str]] = {}
    unscreenable: dict[str, list[str]] = {}
    if autoimmunity_cfg.get("enabled", True):
        proteome = autoimmunity.HumanProteome(str(smk.input.proteome))
        print(
            f"Human proteome: {proteome.n_proteins} proteins, "
            f"{proteome.n_residues} residues, {window}-mer windows."
        )
        screened = autoimmunity.find_matches(
            peptides={record.id: record.sequence for record in records},
            proteome=proteome,
            length=window,
        )
        for consensus_id, result in screened.items():
            if result.matches:
                human_hits[consensus_id] = sorted(
                    {match.describe() for match in result.matches}
                )
            if result.skipped_windows:
                unscreenable[consensus_id] = result.skipped_windows

        # Windows holding an ambiguous residue cannot be searched for literally.
        # Say so loudly: a peptide only partly screened must not be read off the
        # table as one that came back clean.
        for consensus_id, skipped in sorted(unscreenable.items()):
            print(
                f"WARNING {consensus_id}: {len(skipped)} window(s) skipped - "
                f"ambiguous residues, so this peptide is only partly screened."
            )

    # --- combine -------------------------------------------------------------
    out_rows: list[dict] = []
    passed: list[tuple[str, str]] = []

    for record in records:
        reasons: list[str] = []

        score, is_antigenic, note = antigenicity.get(record.id, (None, True, ""))
        if antigenicity_cfg.get("enabled", True) and not is_antigenic:
            if score is None:
                reasons.append(note or "antigenicity could not be scored")
            else:
                reasons.append(
                    f"antigenicity {float(score):.4f} below threshold {antigenicity_cfg['threshold']}"
                )

        allergen_score, is_allergen = allergenicity.get(record.id, ("", False))
        if allergenicity_cfg.get("enabled", True) and is_allergen:
            reasons.append("predicted allergen")

        hits = human_hits.get(record.id, [])
        autoimmune = bool(hits)
        if autoimmune and autoimmunity_cfg.get("drop_on_match", True):
            reasons.append(f"{len(hits)} exact {window}-mer match(es) to the human proteome")

        # Recorded, but deliberately not a failure: an unsearchable window is
        # missing evidence, not adverse evidence. It rides in the reason text so
        # it cannot be missed when reading the table.
        skipped = unscreenable.get(record.id, [])
        if skipped:
            reasons.append(
                f"NOTE {len(skipped)} window(s) unscreened for autoimmunity "
                "(ambiguous residue)"
            )

        cluster = cluster_rows.get(record.id, {})
        # A peptide failing only on the note above has no real failure, so it
        # still passes - test the substantive reasons, not the list's length.
        passes = not [reason for reason in reasons if not reason.startswith("NOTE ")]
        out_rows.append(
            {
                "consensus_id": record.id,
                "protein_label": protein,
                "cluster_id": cluster.get("cluster_id", ""),
                "sequence": record.sequence,
                "length": len(record.sequence),
                "strains": cluster.get("strains", ""),
                "n_strains": cluster.get("n_strains", ""),
                "antigenicity_score": "" if score is None else score,
                "antigenic": is_antigenic,
                "allergen_score": allergen_score,
                "allergen": is_allergen,
                "human_match_count": len(hits),
                "human_match_peptides": ";".join(sorted(hits)[:20]),
                "autoimmune_risk": autoimmune,
                "windows_skipped": len(skipped),
                "passes": passes,
                "reason": "; ".join(reasons),
            }
        )
        if passes:
            passed.append((record.id, record.sequence))

    tables.write_table(str(smk.output.table), schemas.SCREENING, out_rows)
    fasta_mod.write_fasta(str(smk.output.passed), passed, width=0)

    lines = [
        f"Consensus peptides screened: **{len(records)}**  ",
        f"Passed every enabled screen: **{len(passed)}**",
        "",
        "| Consensus | Len | Species | Antigenicity | Allergen | Human matches | Passes | Reason |",
        "| --- | ---: | --- | ---: | :---: | ---: | :---: | --- |",
    ]
    for row in out_rows:
        species = str(row["strains"]).replace("|", ", ") or "-"
        lines.append(
            f"| `{row['consensus_id']}` | {row['length']} | {species} | "
            f"{row['antigenicity_score']} | "
            f"{'yes' if row['allergen'] else 'no'} | {row['human_match_count']} | "
            f"{'yes' if row['passes'] else 'no'} | {row['reason']} |"
        )

    # Caveats that change how the table above should be read. Kept out of the
    # per-peptide `reason` column because they apply to the whole run, not to
    # any one candidate.
    caveats: list[str] = []
    if inert_components:
        caveats.append(
            f"**Allergenicity rests on composition alone.** "
            f"{' and '.join(inert_components)} were zero for every peptide, so "
            "AlgPred's hybrid score equals its ML score throughout. This is "
            "expected at 20-30 aa - BLAST finds no allergen homologs and the "
            "IgE motifs do not hit - but the 0.3 threshold is characterised "
            "for the full hybrid, so treat an allergen call here as weak "
            "evidence rather than a settled verdict."
        )
    if unscreenable:
        caveats.append(
            f"**{len(unscreenable)} peptide(s) only partly screened for "
            "autoimmunity.** Windows holding an ambiguous residue cannot be "
            "searched for literally; see `windows_skipped`. A zero human-match "
            "count for those peptides means incomplete, not clean."
        )
    if caveats:
        lines.extend(["", "## Caveats", ""])
        lines.extend(f"- {caveat}" for caveat in caveats)

    _common.write_summary(str(smk.output.summary), f"Stage 6 - {protein} screening", lines)
    print(f"{protein}: {len(passed)}/{len(records)} consensus peptides passed screening.")
