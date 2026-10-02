"""Stage 2 - whole-protein antigenicity screen."""

# No `from __future__ import annotations` - see 01_manifest.py.

import os

import _common

_common.bootstrap()

from mpvpredpip import schemas, tables  # noqa: E402
from mpvpredpip.adapters import screening  # noqa: E402

smk = snakemake  # noqa: F821

with _common.tee_log(_common.log_path(smk)):
    threshold = float(smk.params.threshold)
    raw_dir = str(smk.params.raw_dir)
    os.makedirs(raw_dir, exist_ok=True)

    rows: list[dict] = []
    for fasta_path in smk.input.proteins:
        protein = os.path.splitext(os.path.basename(str(fasta_path)))[0]
        rows.extend(
            screening.predict_antigenicity(
                fasta_path=str(fasta_path),
                out_path=os.path.join(raw_dir, f"{protein}_antigenicity.tsv"),
                threshold=threshold,
                tool=str(smk.params.tool),
                command=str(smk.params.command or ""),
                tool_dir=str(smk.params.tool_dir),
            )
        )

    tables.write_table(str(smk.output.table), schemas.ANTIGENICITY, rows)

    passed = [row for row in rows if row["antigenic"]]
    failed = [row for row in rows if not row["antigenic"]]

    lines = [
        f"Predictor: `{rows[0]['tool'] if rows else 'n/a'}`  ",
        f"Threshold: **{threshold}**  ",
        f"Sequences scored: **{len(rows)}** - {len(passed)} at or above threshold, "
        f"{len(failed)} below.",
        "",
        "| Sequence | Length | Score | Antigenic |",
        "| --- | ---: | ---: | :---: |",
    ]
    for row in sorted(rows, key=lambda r: -float(r["antigenicity_score"])):
        mark = "yes" if row["antigenic"] else "no"
        lines.append(
            f"| `{row['sequence_id']}` | {row['length']} | "
            f"{float(row['antigenicity_score']):.4f} | {mark} |"
        )

    if failed:
        lines += [
            "",
            "Sequences below the threshold are reported but **not** dropped: this stage is",
            "informational, and antigenicity is applied as a hard filter to the final",
            "consensus peptides in stage 6. Set `antigenicity.drop_failing_proteins: true`",
            "in `config.yaml` to stop here instead.",
        ]

    if smk.params.drop_failing and failed:
        raise SystemExit(
            "antigenicity.drop_failing_proteins is set and these sequences are below "
            f"the threshold of {threshold}:\n"
            + "\n".join(f"  {row['sequence_id']}: {row['antigenicity_score']}" for row in failed)
        )

    _common.write_summary(str(smk.output.summary), "Stage 2 - whole-protein antigenicity", lines)
    print(f"Scored {len(rows)} sequences; {len(passed)} at or above {threshold}.")
