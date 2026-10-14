"""Stage 6 - antigenicity, allergenicity and autoimmunity screening.

Screens applied to the consensus peptides that passed stage 5. A peptide is
carried into stage 9 only if it is predicted antigenic and does not share an
exact window with the human proteome.

Allergenicity is the exception: it is recorded but advisory by default, because
at 20-30 aa AlgPred's hybrid degenerates to its composition model alone. See
screening.allergenicity.drop_on_allergen in config.yaml, which restores the
hard filter.

Every peptide keeps a row in the output table whether it passed or not, with the
reason recorded, so the candidate list can always be traced back.

No ``from __future__ import annotations`` - see 01_manifest.py.
"""

import os

import _common

_common.bootstrap()

from mpvpredpip import allergen_homology  # noqa: E402
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
    homology_cfg = dict(smk.params.allergen_homology or {})
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

    # --- allergen homology ---------------------------------------------------
    # Sequence evidence to set beside AlgPred's composition score. Two criteria
    # on two different inputs - exact matching on the constituent epitopes,
    # which are real predicted sequences, and percent identity on the
    # consensus, which is synthesised. See mpvpredpip.allergen_homology.
    homology: dict = {}
    exact_length = int(homology_cfg.get("exact_match_length", 8))
    identity_threshold = float(homology_cfg.get("identity_threshold", 30))
    if homology_cfg.get("enabled", False):
        allergens = allergen_homology.AllergenSet(str(smk.input.allergens))
        member_rows = tables.read_table(str(smk.input.members))
        wanted = {record.id for record in records}
        consensus_seq = {record.id: record.sequence for record in records}
        members: dict[str, dict[str, str]] = {}
        for row in member_rows:
            consensus_id = row.get("consensus_id", "")
            if consensus_id not in wanted:
                continue
            # The row carrying the consensus itself is not a constituent
            # epitope; '-' is how stage 5 writes an absent field.
            if row.get("member_index", "") == "Consensus":
                continue
            sequence = row.get("sequence", "")
            epitope_id = row.get("epitope_id", "")
            if not sequence or sequence == "-" or not epitope_id or epitope_id == "-":
                continue
            members.setdefault(consensus_id, {})[epitope_id] = sequence

        print(
            f"Allergen set: {allergens.n_sequences} sequences, "
            f"{allergens.n_residues} residues; {exact_length}-mer exact match on "
            f"{sum(len(v) for v in members.values())} constituent epitopes, "
            f"{identity_threshold:g}% identity on {len(consensus_seq)} consensus "
            "sequence(s)."
        )
        homology = allergen_homology.screen(
            consensus=consensus_seq,
            members=members,
            allergens=allergens,
            exact_match_length=exact_length,
            identity_threshold=identity_threshold,
        )
        # Said once per run, because it governs how the identity column should
        # be read and is easy to lose sight of per-peptide.
        flagged_ident = [
            key for key, result in homology.items()
            if allergen_homology.identity_flagged(result, identity_threshold)
        ]
        if flagged_ident and not homology_cfg.get("drop_on_identity", False):
            print(
                f"NOTE allergen homology: {len(flagged_ident)} of {len(homology)} "
                f"consensus sequence(s) reached {identity_threshold:g}% identity, "
                "which is at chance for this peptide length - recorded, not "
                "filtered. See screening.allergen_homology in config.yaml."
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
            # Advisory by default - see screening.allergenicity.drop_on_allergen
            # in config.yaml. The NOTE prefix is what keeps the call in the
            # reason column without failing the peptide, exactly as it does for
            # an unscreenable autoimmunity window below.
            if allergenicity_cfg.get("drop_on_allergen", True):
                reasons.append("predicted allergen")
            else:
                reasons.append(
                    "NOTE predicted allergen (advisory, not filtered)"
                )

        # Allergen homology. The two criteria are kept apart in the reason text
        # as well as in the columns: an exact run of residues shared with a
        # known allergen is real evidence, while a percent identity at this
        # length is at chance, and collapsing them into one "homology" verdict
        # would let the weaker one borrow the stronger one's authority.
        homology_result = homology.get(record.id)
        exact_epitopes = 0
        exact_matches = ""
        identity_percent = ""
        identity_match = ""
        homology_risk = False
        if homology_result is not None:
            exact_epitopes = homology_result.n_exact_epitopes
            exact_matches = homology_result.describe_exact()
            identity_percent = round(homology_result.identity.percent, 1)
            identity_match = homology_result.identity.allergen_id
            identity_hit = allergen_homology.identity_flagged(
                homology_result, identity_threshold
            )
            homology_risk = homology_result.exact_flagged or identity_hit

            if homology_result.exact_flagged:
                detail = (
                    f"{exact_epitopes} constituent epitope(s) share an exact "
                    f"{exact_length}-mer with a known allergen"
                )
                if homology_cfg.get("drop_on_exact_match", True):
                    reasons.append(detail)
                else:
                    reasons.append(f"NOTE {detail} (advisory, not filtered)")

            if identity_hit:
                detail = (
                    f"consensus {identity_percent}% identical to "
                    f"{identity_match or 'a known allergen'} "
                    f"(threshold {identity_threshold:g}%)"
                )
                if homology_cfg.get("drop_on_identity", False):
                    reasons.append(detail)
                else:
                    reasons.append(f"NOTE {detail} (at chance for this length)")

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
                "allergen_exact_epitopes": exact_epitopes,
                "allergen_exact_matches": exact_matches,
                "allergen_identity_percent": identity_percent,
                "allergen_identity_match": identity_match,
                "allergen_homology_risk": homology_risk,
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
        "| Consensus | Len | Species | Antigenicity | Allergen | Exact | Identity | Human matches | Passes | Reason |",
        "| --- | ---: | --- | ---: | :---: | ---: | ---: | ---: | :---: | --- |",
    ]
    for row in out_rows:
        species = str(row["strains"]).replace("|", ", ") or "-"
        # "-" rather than blank when the homology screen is off, so an absent
        # screen reads differently from one that found nothing.
        exact = "-" if not homology else str(row["allergen_exact_epitopes"])
        identity = "-" if not homology else f"{row['allergen_identity_percent']}%"
        lines.append(
            f"| `{row['consensus_id']}` | {row['length']} | {species} | "
            f"{row['antigenicity_score']} | "
            f"{'yes' if row['allergen'] else 'no'} | {exact} | {identity} | "
            f"{row['human_match_count']} | "
            f"{'yes' if row['passes'] else 'no'} | "
            # Escaped, because the reason now carries allergen accessions and
            # those are pipe-delimited - `sp|P01574|IFNB_HUMAN` would open two
            # extra columns and stop the table rendering as a table. Backticks
            # do not protect it; markdown splits the row before parsing spans.
            f"{_common.md_cell(row['reason'])} |"
        )

    # Caveats that change how the table above should be read. Kept out of the
    # per-peptide `reason` column because they apply to the whole run, not to
    # any one candidate.
    caveats: list[str] = []
    if homology:
        ident_hits = [
            row for row in out_rows
            if row["allergen_identity_percent"] != ""
            and float(row["allergen_identity_percent"]) >= identity_threshold
        ]
        exact_hits = [row for row in out_rows if row["allergen_exact_epitopes"]]
        caveats.append(
            f"**Allergen homology was screened two ways.** The "
            f"{exact_length}-mer exact-match criterion ran on the constituent "
            f"epitopes and flagged {len(exact_hits)} of {len(out_rows)} "
            "candidate(s); this is the FAO/WHO short-sequence criterion and is "
            "the stronger evidence. The "
            f"{identity_threshold:g}% identity criterion ran on the consensus "
            f"sequence and flagged {len(ident_hits)} of {len(out_rows)}. "
            "Those two numbers are not comparable: at 20-30 aa a "
            f"{identity_threshold:g}% identity cut flags composition-matched "
            "shuffled controls about as often as real peptides, so it is "
            "recorded rather than acted on. See "
            "`screening.allergen_homology` in config.yaml for the measurement."
        )
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
    if (
        allergenicity_cfg.get("enabled", True)
        and not allergenicity_cfg.get("drop_on_allergen", True)
    ):
        flagged = [row for row in out_rows if row["allergen"]]
        # Said as "no longer removes them on that basis", not "were not
        # filtered out": a peptide called allergen can still fail antigenicity
        # or autoimmunity, and on this run most of them do. The blunter wording
        # reads as though all of them passed.
        also_failed = [row for row in flagged if not row["passes"]]
        detail = ""
        if also_failed:
            detail = (
                f" {len(also_failed)} of those still fail another screen, so "
                "this did not rescue every one."
            )
        caveats.append(
            f"**Allergenicity is advisory on this run.** "
            f"{len(flagged)} of {len(out_rows)} peptide(s) were called allergen; "
            "that call is recorded but no longer removes them on its own, and "
            "the score stays in the `allergen` and `allergen_score` columns."
            + detail +
            " Set `screening.allergenicity.drop_on_allergen: true` in "
            "config.yaml to make it a hard filter again."
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
