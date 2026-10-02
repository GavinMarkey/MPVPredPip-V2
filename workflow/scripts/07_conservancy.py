"""Stage 7 - conservancy of the screened consensus peptides across strains."""

# No `from __future__ import annotations` - see 01_manifest.py.

import _common

_common.bootstrap()

from vaxpipe import conservancy as conservancy_mod  # noqa: E402

smk = snakemake  # noqa: F821

with _common.tee_log(_common.log_path(smk)):
    protein = str(smk.wildcards.protein)
    threshold = float(smk.params.threshold)
    mode = str(smk.params.mode)

    rows = conservancy_mod.run(
        consensus_fasta=str(smk.input.consensus),
        manifest_path=str(smk.input.manifest),
        protein_label=protein,
        out_table=str(smk.output.table),
        out_targets_fasta=str(smk.output.targets),
        out_epitopes_txt=str(smk.output.epitopes),
        threshold=threshold,
        deduplicate=bool(smk.params.deduplicate),
    )

    lines = [
        f"Peptides analysed: **{len(rows)}**  ",
        f"Identity threshold: **{threshold:.0%}**  ",
        f"Mode: `{mode}`",
        "",
        "| Consensus | Len | Conservancy | Matched strains | Min id | Max id |",
        "| --- | ---: | ---: | --- | ---: | ---: |",
    ]
    for row in sorted(rows, key=lambda r: -float(r["conservancy_percent"])):
        lines.append(
            f"| `{row['consensus_id']}` | {len(str(row['sequence']))} | "
            f"{row['conservancy_percent']}% | {row['matched_strains'] or '-'} | "
            f"{row['min_identity']}% | {row['max_identity']}% |"
        )

    lines += [
        "",
        "Conservancy is computed on ungapped sliding windows, independently of the MSA,",
        "following the IEDB definition: the fraction of non-redundant target sequences in",
        "which the peptide is found at or above the identity threshold.",
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
    print(f"{protein}: conservancy computed for {len(rows)} peptides.")
