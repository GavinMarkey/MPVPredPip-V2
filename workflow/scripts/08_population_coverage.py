"""Stage 8 - population coverage of the alleles behind each consensus peptide.

Alleles are taken from the T-cell epitopes that are members of each consensus
peptide's cluster, traced through the stage 5 membership table. That is
deliberate: the coverage figure then describes the alleles that actually drove
the cluster's selection, rather than a fresh prediction against the consensus
sequence that would not correspond to the evidence the candidate rests on.

No ``from __future__ import annotations`` - see 01_manifest.py. This script runs
under iedb.yaml, which pins **Python 3.8**, so annotations here must not use
builtin generics or ``X | Y`` at runtime.
"""

import argparse
import os

import _common

_common.bootstrap()

from mpvpredpip import schemas, tables  # noqa: E402
from mpvpredpip.adapters import iedb  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--protein", required=True)
parser.add_argument("--screening", required=True)
parser.add_argument("--members", required=True)
parser.add_argument("--epitopes", required=True)
parser.add_argument("--table", required=True, help="coverage table to write")
parser.add_argument("--input-file", required=True, help="the input built for the IEDB tool")
parser.add_argument("--raw", required=True, help="the tool's own output")
parser.add_argument("--summary", required=True)
parser.add_argument("--tool-dir", required=True)
parser.add_argument("--populations", required=True, help="comma-separated population names")
parser.add_argument("--mhc-class", required=True)
parser.add_argument("--log", default=None)
args = parser.parse_args()

with _common.tee_log(args.log):
    protein = args.protein
    mhc_class = args.mhc_class

    screened = tables.read_table(args.screening)
    members = tables.read_table(args.members)
    epitopes = tables.read_table(args.epitopes)

    passing = [row for row in screened if tables.as_bool(row["passes"])]
    if not passing:
        print(f"{protein}: no peptides passed stage 6; nothing to compute coverage for.")
        tables.write_table(args.table, schemas.POPULATION_COVERAGE, [])
        for path in (args.input_file, args.raw):
            os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("# no peptides reached stage 8\n")
        _common.write_summary(
            args.summary,
            f"Stage 8 - {protein} population coverage",
            ["No peptides passed stage 6, so population coverage was not computed."],
        )
        raise SystemExit(0)

    # epitope_id -> alleles, restricted to T-cell epitopes
    alleles_by_epitope = {}  # dict[str, set[str]]
    class_by_epitope = {}    # dict[str, str]
    for row in epitopes:
        if row["epitope_class"] not in {"MHC-I", "MHC-II"}:
            continue
        if row["allele"]:
            alleles_by_epitope.setdefault(row["epitope_id"], set()).add(row["allele"])
        class_by_epitope[row["epitope_id"]] = row["epitope_class"]

    # consensus_id -> alleles contributed by its cluster members
    per_consensus = {}        # dict[str, set[str]]
    per_consensus_class = {}  # dict[str, dict[str, set[str]]]
    for row in members:
        consensus_id = row["consensus_id"]
        epitope_id = row["epitope_id"]
        alleles = alleles_by_epitope.get(epitope_id)
        if not alleles:
            continue
        per_consensus.setdefault(consensus_id, set()).update(alleles)
        bucket = per_consensus_class.setdefault(consensus_id, {"MHC-I": set(), "MHC-II": set()})
        bucket.setdefault(class_by_epitope.get(epitope_id, "MHC-I"), set()).update(alleles)

    entries = []  # list[tuple[str, list[str]]]
    for row in passing:
        alleles = sorted(per_consensus.get(row["consensus_id"], set()))
        if alleles:
            entries.append((str(row["sequence"]), alleles))
        else:
            print(
                f"WARNING: {row['consensus_id']} has no T-cell alleles traced from its "
                "cluster members and is excluded from the coverage calculation."
            )

    if not entries:
        raise SystemExit(
            f"{protein}: none of the passing peptides could be linked back to MHC alleles. "
            "Check that stage 5 cluster membership resolved correctly."
        )

    iedb.write_population_input(args.input_file, entries)
    iedb.run_population_coverage(
        input_path=args.input_file,
        out_path=args.raw,
        tool_dir=args.tool_dir,
        populations=[p for p in args.populations.split(",") if p.strip()],
        mhc_class=mhc_class,
    )

    coverage_rows = iedb.parse_population_coverage(args.raw)

    # The tool reports coverage over the whole submitted allele set, so the result
    # is per protein rather than per peptide. The per-peptide allele lists are
    # carried through so the final table can show what each candidate contributed.
    all_alleles = sorted({allele for _, alleles in entries for allele in alleles})
    out_rows = [
        {
            "consensus_id": "ALL",
            "protein_label": protein,
            "population": row["population"],
            "mhc_class": mhc_class,
            "coverage_percent": row["coverage_percent"],
            "average_hit": row.get("average_hit", ""),
            "pc90": row.get("pc90", ""),
            "n_alleles": len(all_alleles),
            "alleles": "|".join(all_alleles),
        }
        for row in coverage_rows
    ]

    for row in passing:
        alleles = sorted(per_consensus.get(row["consensus_id"], set()))
        buckets = per_consensus_class.get(row["consensus_id"], {})
        out_rows.append(
            {
                "consensus_id": row["consensus_id"],
                "protein_label": protein,
                "population": "",
                "mhc_class": mhc_class,
                "coverage_percent": "",
                "average_hit": "",
                "pc90": "",
                "n_alleles": len(alleles),
                "alleles": "|".join(alleles),
            }
        )

    tables.write_table(args.table, schemas.POPULATION_COVERAGE, out_rows)

    lines = [
        f"Peptides included: **{len(entries)}**  ",
        f"Distinct alleles: **{len(all_alleles)}**  ",
        f"MHC class: `{mhc_class}`",
        "",
        "| Population | Coverage | Average hit | PC90 |",
        "| --- | ---: | ---: | ---: |",
    ]
    for row in coverage_rows:
        lines.append(
            f"| {row['population']} | {row['coverage_percent']}% | "
            f"{row.get('average_hit', '-')} | {row.get('pc90', '-')} |"
        )

    _common.write_summary(
        args.summary, f"Stage 8 - {protein} population coverage", lines
    )
    print(f"{protein}: coverage computed over {len(all_alleles)} alleles from {len(entries)} peptides.")
