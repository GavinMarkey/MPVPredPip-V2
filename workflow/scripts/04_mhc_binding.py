"""Stage 4 - MHC class I or class II binding prediction.

Shared by both binding rules; ``params.mhc_class`` selects which. Reproduces the
Tepitool selection locally: predict against an explicit allele panel, then keep
peptides at or below the configured percentile-rank cut-off.

No ``from __future__ import annotations`` - see 01_manifest.py. This script runs
under iedb.yaml, which pins **Python 3.8**, so annotations here must not use
builtin generics or ``X | Y`` at runtime.
"""

import argparse

import _common

_common.bootstrap()

from mpvpredpip import epitopes as epitope_mod  # noqa: E402
from mpvpredpip import fasta as fasta_mod  # noqa: E402
from mpvpredpip import manifest as manifest_mod  # noqa: E402
from mpvpredpip import schemas, tables  # noqa: E402
from mpvpredpip.adapters import iedb  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fasta", required=True)
parser.add_argument("--table", required=True, help="epitope table to write")
parser.add_argument("--raw", required=True, help="the predictor's own output")
parser.add_argument("--tool-dir", required=True)
parser.add_argument("--allele-file", required=True)
parser.add_argument("--lengths", required=True, help="comma-separated peptide lengths")
parser.add_argument("--cutoff", type=float, required=True)
parser.add_argument("--method", required=True)
parser.add_argument("--mhc-class", required=True)
parser.add_argument("--source-tool", required=True)
parser.add_argument("--log", default=None)
args = parser.parse_args()

with _common.tee_log(args.log):
    fasta_path = args.fasta
    records = fasta_mod.read_fasta(fasta_path)
    if not records:
        raise SystemExit(f"No sequences in {fasta_path}")

    # The IEDB predictors report seq_num as the 1-based index of the sequence in
    # the submitted FASTA, which is the only link back to strain and accession.
    id_by_seqnum = {index + 1: record.id for index, record in enumerate(records)}
    by_id = {record.id: record for record in records}

    alleles = iedb.read_alleles(args.allele_file)
    lengths = [int(length) for length in args.lengths.split(",") if length.strip()]
    cutoff = args.cutoff
    mhc_class = args.mhc_class

    print(
        f"MHC class {mhc_class}: {len(alleles)} alleles x lengths {lengths} "
        f"over {len(records)} sequence(s), keeping percentile <= {cutoff}"
    )

    iedb.predict_binding(
        fasta_path=fasta_path,
        out_path=args.raw,
        tool_dir=args.tool_dir,
        alleles=alleles,
        lengths=lengths,
        method=args.method,
        mhc_class=mhc_class,
    )

    hits = iedb.parse_binding(
        path=args.raw,
        percentile_cutoff=cutoff,
        id_by_seqnum=id_by_seqnum,
    )

    found = []  # list[epitope_mod.Epitope]
    for hit in hits:
        record = by_id.get(str(hit["sequence_id"]))
        if record is None:
            raise SystemExit(
                f"Binding output references sequence number {hit['seq_num']} which is not "
                f"in the input FASTA ({len(records)} sequences). Raw output: {args.raw}"
            )
        protein, strain, accession = manifest_mod.parse_header(record.header)
        found.append(
            epitope_mod.Epitope(
                protein_label=protein,
                strain=strain,
                accession=accession,
                source_tool=args.source_tool,
                start=int(hit["start"]),
                end=int(hit["end"]),
                sequence=str(hit["peptide"]),
                score=float(hit["score"]),
                allele=str(hit["allele"]),
                percentile=str(hit["percentile"]),
            )
        )

    found = epitope_mod.assign_ids(found)
    tables.write_table(args.table, schemas.EPITOPES, [e.as_row() for e in found])

    unique = {epitope.sequence for epitope in found}
    print(f"Kept {len(found)} binder rows ({len(unique)} unique peptides) at percentile <= {cutoff}.")
    if not found:
        print(
            f"WARNING: no MHC class {mhc_class} binders passed. Check the raw output at "
            f"{args.raw} and consider relaxing the percentile cut-off."
        )
