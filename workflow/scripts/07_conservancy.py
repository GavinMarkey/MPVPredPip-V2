"""Stage 7 - conservancy of the screened consensus peptides across strains."""

# No `from __future__ import annotations` - see 01_manifest.py.

import _common

_common.bootstrap()

from mpvpredpip import conservancy as conservancy_mod  # noqa: E402

smk = snakemake  # noqa: F821

with _common.tee_log(_common.log_path(smk)):
    protein = str(smk.wildcards.protein)
    threshold = float(smk.params.threshold)
    mode = str(smk.params.mode)
    reference_fastas = smk.params.reference_fastas  # dict or None

    rows = conservancy_mod.run(
        consensus_fasta=str(smk.input.consensus),
        manifest_path=str(smk.input.manifest),
        protein_label=protein,
        out_table=str(smk.output.table),
        out_targets_fasta=str(smk.output.targets),
        out_epitopes_txt=str(smk.output.epitopes),
        threshold=threshold,
        deduplicate=bool(smk.params.deduplicate),
        reference_fastas=reference_fastas,
    )

    # --- summary ---
    using_refs = reference_fastas is not None
    strains_seen = sorted({str(r["strain"]) for r in rows})

    lines = [
        f"Peptides analysed: **{len({r['consensus_id'] for r in rows})}**  ",
        f"Identity threshold: **{threshold:.0%}**  ",
        f"Mode: `{mode}`  ",
        f"Target source: **{'NCBI reference sets' if using_refs else 'manifest (input FASTAs)'}**",
        "",
    ]

    if using_refs:
        lines += [
            f"Strains: {', '.join(strains_seen)}",
            "",
        ]

    # Per-strain tables when we have multiple strains
    if len(strains_seen) > 1:
        for strain in strains_seen:
            strain_rows = [r for r in rows if r["strain"] == strain]
            n_targets = strain_rows[0]["n_target_sequences"] if strain_rows else 0
            lines += [
                f"### {strain} ({n_targets} target sequences)",
                "",
                "| Consensus | Len | Conservancy | Matched | Min id | Max id |",
                "| --- | ---: | ---: | ---: | ---: | ---: |",
            ]
            for row in sorted(strain_rows, key=lambda r: -float(r["conservancy_percent"])):
                lines.append(
                    f"| `{_common.md_cell(str(row['consensus_id']))}` "
                    f"| {len(str(row['sequence']))} | "
                    f"{row['conservancy_percent']}% | {row['n_matched']}/{row['n_target_sequences']} | "
                    f"{row['min_identity']}% | {row['max_identity']}% |"
                )
            lines.append("")
    else:
        # Single strain or manifest fallback — flat table
        lines += [
            "| Consensus | Len | Conservancy | Matched accessions | Min id | Max id |",
            "| --- | ---: | ---: | --- | ---: | ---: |",
        ]
        for row in sorted(rows, key=lambda r: -float(r["conservancy_percent"])):
            lines.append(
                f"| `{_common.md_cell(str(row['consensus_id']))}` "
                f"| {len(str(row['sequence']))} | "
                f"{row['conservancy_percent']}% | {row['matched_accessions'] or '-'} | "
                f"{row['min_identity']}% | {row['max_identity']}% |"
            )

    lines += [
        "",
        "Conservancy is computed on ungapped sliding windows, independently of the MSA,",
        "following the IEDB definition: the fraction of non-redundant target sequences in",
        "which the peptide is found at or above the identity threshold.",
    ]

    if using_refs:
        lines += [
            "",
            "Target sequences were fetched from NCBI by ``setup/fetch_ncbi_conservancy.py``.",
            "See the ``.provenance.json`` files in ``resources/conservancy/`` for exact",
            "counts, search terms and accession lists.",
        ]

    if mode in {"export", "both"}:
        lines += [
            "",
            "### Cross-checking against the IEDB web tool",
            "",
            "Upload-ready files have been written to `iedb_upload/`:",
            "",
            f"- `non_redundant_sequences.fasta` - the deduplicated {protein} sequences",
            "- `epitope_list.txt` - one peptide per line",
            "",
            "Paste these into the IEDB conservancy analysis web tool to confirm the",
            "numbers above. The local calculation is what the pipeline uses so that runs",
            "stay unattended and reproducible.",
        ]

    _common.write_summary(str(smk.output.summary), f"Stage 7 - {protein} conservancy", lines)
    n_peptides = len({r["consensus_id"] for r in rows})
    n_strains = len(strains_seen)
    print(f"{protein}: conservancy computed for {n_peptides} peptides across {n_strains} strain(s).")
