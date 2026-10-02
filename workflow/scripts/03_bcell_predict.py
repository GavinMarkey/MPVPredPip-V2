"""Stage 3 - run one B-cell predictor and convert per-residue scores to epitopes.

Shared by the BepiPred and EpiDope rules: which predictor runs is decided by
``params.source_tool``, and the interval-building rules are identical for both so
their epitopes stay comparable when stage 5 clusters them together.

*** Invoked through ``shell:``, not ``script:``. ***

EpiDope's environment pins Python 3.8, and Snakemake's ``script:`` directive
imports Snakemake itself inside the rule's environment to unpickle the job
object. Snakemake 9 is written in Python 3.10+ syntax, so that import fails with
``TypeError: unsupported operand type(s) for |``. Taking parameters on the
command line instead means this script never imports Snakemake and runs happily
on 3.8. BepiPred's rule calls it the same way, so both predictors stay on one
code path.

For the same reason there is no ``from __future__ import annotations`` here, and
annotations must not use builtin generics or ``X | Y`` at runtime.
"""

import argparse
import os

import _common

_common.bootstrap()

from vaxpipe import epitopes as epitope_mod  # noqa: E402
from vaxpipe import fasta as fasta_mod  # noqa: E402
from vaxpipe import manifest as manifest_mod  # noqa: E402
from vaxpipe import schemas, tables  # noqa: E402
from vaxpipe.adapters import bcell  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fasta", required=True)
parser.add_argument("--residues", required=True, help="per-residue score table to write")
parser.add_argument("--epitopes", required=True, help="epitope table to write")
parser.add_argument("--raw-dir", required=True, help="where the predictor's own output is kept")
parser.add_argument("--tool-dir", default="", help="install directory, or executable for EpiDope")
parser.add_argument("--source-tool", required=True, choices=["bepipred", "epidope"])
parser.add_argument("--threshold", type=float, required=True)
parser.add_argument("--length", type=int, required=True, help="epitope length in residues")
parser.add_argument("--step", type=int, required=True, help="offset between overlapping epitopes")
parser.add_argument(
    "--min-region",
    type=int,
    required=True,
    help="shortest above-threshold region worth padding out to --length",
)
parser.add_argument("--merge-gap", type=int, required=True)
parser.add_argument("--threads", type=int, default=1)
parser.add_argument(
    "--esm-dir",
    default="",
    help="BepiPred only: where to cache ESM-2 encodings. Kept outside --raw-dir so "
         "that clearing the raw output does not discard them.",
)
parser.add_argument(
    "--skip-predict",
    action="store_true",
    help="parse an existing --raw-dir instead of running the predictor. EpiDope's "
         "environment pins Python 3.6 and cannot import this package, so its rule "
         "runs the tool by plain shell command and parses in a second rule.",
)
parser.add_argument("--log", default=None)
args = parser.parse_args()

with _common.tee_log(args.log):
    source_tool = args.source_tool
    fasta_path = args.fasta
    raw_dir = args.raw_dir
    os.makedirs(raw_dir, exist_ok=True)

    records = {record.id: record for record in fasta_mod.read_fasta(fasta_path)}
    if not records:
        raise SystemExit(f"No sequences in {fasta_path}")

    if args.skip_predict:
        # The predictor already ran in its own environment; only parse.
        if source_tool == "bepipred":
            scores = bcell.parse_bepipred(raw_dir)
        else:
            scores = bcell.parse_epidope(raw_dir)
    elif source_tool == "bepipred":
        scores = bcell.predict_bepipred(
            fasta_path=fasta_path,
            out_dir=raw_dir,
            tool_dir=args.tool_dir,
            threshold=args.threshold,
            esm_dir=args.esm_dir,
        )
    else:
        scores = bcell.predict_epidope(
            fasta_path=fasta_path,
            out_dir=raw_dir,
            executable=args.tool_dir or "epidope",
            threads=args.threads,
        )

    bcell.save_per_residue(args.residues, scores)

    # Predictors vary in how much of the defline they keep. Match on the exact ID
    # first, then on a unique prefix, so a truncated identifier is still resolved
    # rather than silently dropping a whole protein's predictions.
    def resolve(sequence_id):
        # -> str | None. Unannotated: this runs under Python 3.8 (epidope.yaml).
        if sequence_id in records:
            return sequence_id
        matches = [key for key in records if key.startswith(sequence_id) or sequence_id.startswith(key)]
        return matches[0] if len(matches) == 1 else None

    found = []        # list[epitope_mod.Epitope]
    unresolved = []   # list[str]

    for sequence_id, vector in scores.items():
        key = resolve(sequence_id)
        if key is None:
            unresolved.append(sequence_id)
            continue

        record = records[key]
        if len(vector) != len(record.sequence):
            raise SystemExit(
                f"{source_tool} returned {len(vector)} residue scores for {key!r} but the "
                f"sequence is {len(record.sequence)} residues. Refusing to guess an "
                "alignment - check the raw output in " + raw_dir
            )

        protein, strain, accession = manifest_mod.parse_header(record.header)
        for start, end, score in epitope_mod.residues_to_epitopes(
            sequence=record.sequence,
            scores=vector,
            threshold=args.threshold,
            length=args.length,
            step=args.step,
            merge_gap=args.merge_gap,
            min_region=args.min_region,
        ):
            found.append(
                epitope_mod.Epitope(
                    protein_label=protein,
                    strain=strain,
                    accession=accession,
                    source_tool=source_tool,
                    start=start,
                    end=end,
                    sequence=record.sequence[start - 1 : end],
                    score=score,
                )
            )

    if unresolved:
        raise SystemExit(
            f"{source_tool} reported sequence IDs that do not match the input FASTA:\n"
            + "\n".join(f"  {name}" for name in unresolved[:10])
            + f"\n\nInput IDs: {sorted(records)[:10]}\nRaw output kept in {raw_dir}"
        )

    found = epitope_mod.assign_ids(found)
    tables.write_table(args.epitopes, schemas.EPITOPES, [e.as_row() for e in found])

    print(
        f"{source_tool}: {len(found)} epitopes across {len(records)} sequence(s) "
        f"at threshold {args.threshold}"
    )
    if not found:
        print(
            f"WARNING: {source_tool} produced no epitopes. The threshold "
            f"({args.threshold}) may be too strict for this protein - check "
            f"{args.residues} for the score distribution."
        )
