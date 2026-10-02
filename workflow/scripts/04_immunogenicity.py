"""Stage 4 - class I immunogenicity screen.

Every MHC class I binder is scored with the IEDB immunogenicity predictor and
only those above ``tcell.immunogenicity.min_score`` continue, as the pipeline
specification requires. Peptides that are dropped are written to their own table
rather than discarded silently, so the loss at this step is auditable.

No ``from __future__ import annotations`` - see 01_manifest.py. This script runs
under iedb.yaml, which pins **Python 3.8**, so annotations here must not use
builtin generics or ``X | Y`` at runtime.
"""

import argparse

import _common

_common.bootstrap()

from mpvpredpip import schemas, tables  # noqa: E402
from mpvpredpip.adapters import iedb  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--binders", required=True)
parser.add_argument("--table", required=True, help="surviving binders to write")
parser.add_argument("--dropped", required=True, help="binders cut by this filter")
parser.add_argument("--raw", required=True, help="the predictor's own output")
parser.add_argument("--tool-dir", required=True)
parser.add_argument("--min-score", type=float, required=True)
# Snakemake renders booleans as "True"/"False", which bool() would both accept.
parser.add_argument("--enabled", required=True, choices=["True", "False"])
parser.add_argument("--log", default=None)
args = parser.parse_args()

with _common.tee_log(args.log):
    binders = tables.read_table(args.binders)
    min_score = args.min_score
    enabled = args.enabled == "True"

    if not enabled:
        print("Immunogenicity screening disabled; passing all class I binders through.")
        tables.write_table(args.table, schemas.EPITOPES, binders)
        tables.write_table(args.dropped, schemas.EPITOPES, [])
        with open(args.raw, "w", encoding="utf-8") as handle:
            handle.write("# immunogenicity screening disabled in config.yaml\n")
        raise SystemExit(0)

    if not binders:
        print("No class I binders to screen.")
        tables.write_table(args.table, schemas.EPITOPES, [])
        tables.write_table(args.dropped, schemas.EPITOPES, [])
        with open(args.raw, "w", encoding="utf-8") as handle:
            handle.write("# no class I binders reached this stage\n")
        raise SystemExit(0)

    peptides = sorted({row["sequence"] for row in binders})
    scores = iedb.predict_immunogenicity(
        peptides=peptides,
        out_path=args.raw,
        tool_dir=args.tool_dir,
    )

    kept = []     # list[dict]
    dropped = []  # list[dict]
    for row in binders:
        score = scores.get(row["sequence"])
        enriched = dict(row)
        # Keep the immunogenicity score alongside the binding score so the reason
        # a peptide survived or was cut stays visible in the table itself.
        enriched["percentile"] = row.get("percentile", "")
        enriched["score"] = row.get("score", "")
        if score is not None and score > min_score:
            kept.append(enriched)
        else:
            dropped.append(enriched)

    tables.write_table(args.table, schemas.EPITOPES, kept)
    tables.write_table(args.dropped, schemas.EPITOPES, dropped)

    unique_kept = {row["sequence"] for row in kept}
    print(
        f"Immunogenicity: kept {len(kept)}/{len(binders)} binder rows "
        f"({len(unique_kept)} unique peptides) with score > {min_score}; "
        f"dropped {len(dropped)}."
    )
    if not kept:
        print(
            "WARNING: every class I binder was dropped by the immunogenicity filter. "
            f"No cluster can satisfy the MHC-I requirement in stage 5. Raw scores: {args.raw}"
        )
