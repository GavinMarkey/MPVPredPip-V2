"""Stage 8 - population coverage of the alleles behind each consensus peptide.

Alleles are taken from the T-cell epitopes that are members of each consensus
peptide's cluster, traced through the stage 5 membership table. That is
deliberate: the coverage figure then describes the alleles that actually drove
the cluster's selection, rather than a fresh prediction against the consensus
sequence that would not correspond to the evidence the candidate rests on.

------------------------------------------------------------------------------
One tool run per peptide, not one for the whole set
------------------------------------------------------------------------------
Given a file of several epitopes the IEDB tool answers a *set* question: what
fraction of the population responds to **at least one** epitope in the file.
That is the right figure for a finished multi-epitope construct, but it is not
a property of any single candidate, and stage 9 ranks candidates. Submitting
all peptides at once therefore produced one number that was then repeated
against every candidate as though it were theirs.

So each peptide is submitted on its own, and the whole set is submitted once
more under the id ``ALL``. Both are kept: per-peptide rows let candidates be
compared, the ``ALL`` row describes what they achieve together. Every run
covers all configured populations and classes, so this costs one tool
invocation per peptide plus one.

A peptide with no alleles of a given class produces no row for that class -
the tool omits the block rather than reporting zero. A missing row means "no
alleles of that class", which is not the same as no coverage.

No ``from __future__ import annotations`` - see 01_manifest.py. This script runs
under iedb.yaml, which pins **Python 3.8**, so annotations here must not use
builtin generics or ``X | Y`` at runtime.
"""

import argparse
import os
import tempfile

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

    def _norm_allele(a):
        # type: (str) -> str
        """Strip the beta chain from MHC-II paired alleles.

        The IEDB population coverage tool expects only the alpha chain,
        e.g. ``HLA-DPA1*02:01`` not ``HLA-DPA1*02:01/DPB1*14:01``.
        """
        return a.split("/", 1)[0] if "/" in a else a

    # Index alleles by *sequence* rather than epitope_id.  The membership
    # table carries only the .1 epitope_id per unique peptide (by design —
    # stage 4 deduplicates to avoid inflating cluster sizes), but the same
    # peptide can bind multiple alleles (.1, .2, .3, …).  Keying by
    # sequence ensures every allele variant is captured.
    sequence_by_epitope = {}   # type: dict[str, str]
    alleles_by_sequence = {}   # type: dict[str, set[str]]
    for row in epitopes:
        if row["epitope_class"] not in ("MHC-I", "MHC-II"):
            continue
        sequence_by_epitope[row["epitope_id"]] = row["sequence"]
        if row["allele"]:
            alleles_by_sequence.setdefault(row["sequence"], set()).add(
                _norm_allele(row["allele"])
            )

    # consensus_id -> alleles contributed by its cluster members. The tool
    # decides each allele's class itself, so they are not split up here.
    per_consensus = {}  # type: dict[str, set[str]]
    for row in members:
        consensus_id = row["consensus_id"]
        seq = sequence_by_epitope.get(row["epitope_id"])
        if not seq:
            continue
        alleles = alleles_by_sequence.get(seq)
        if not alleles:
            continue
        per_consensus.setdefault(consensus_id, set()).update(alleles)

    entries = []  # list[tuple[str, str, list[str]]] - consensus_id, sequence, alleles
    for row in passing:
        alleles = sorted(per_consensus.get(row["consensus_id"], set()))
        if alleles:
            entries.append((str(row["consensus_id"]), str(row["sequence"]), alleles))
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

    populations = [p for p in args.populations.split(",") if p.strip()]
    all_alleles = sorted({allele for _, _, alleles in entries for allele in alleles})

    # The submitted set, kept as the record of what went in and as an
    # upload-ready file for the IEDB web tool.
    iedb.write_population_input(
        args.input_file, [(sequence, alleles) for _, sequence, alleles in entries]
    )

    raw_chunks = []  # list[str] - the tool's own reports, kept verbatim

    def _coverage_for(label, submission, scratch):
        # type: (str, list, str) -> list
        """Run the tool over *submission* and return its parsed coverage rows."""
        submission_path = os.path.join(scratch, "input.tsv")
        report_path = os.path.join(scratch, "report.txt")
        iedb.write_population_input(submission_path, submission)
        iedb.run_population_coverage(
            input_path=submission_path,
            out_path=report_path,
            tool_dir=args.tool_dir,
            populations=populations,
            mhc_class=mhc_class,
        )
        with open(report_path, encoding="utf-8", errors="replace") as handle:
            raw_chunks.append("{0}\n{1}\n{0}\n{2}".format("=" * 78, label, handle.read()))
        return iedb.parse_population_coverage(report_path)

    out_rows = []
    with tempfile.TemporaryDirectory(prefix="popcov_") as scratch_root:
        for index, (consensus_id, sequence, alleles) in enumerate(entries):
            scratch = os.path.join(scratch_root, str(index))
            os.makedirs(scratch)
            label = (
                f"PER-EPITOPE {consensus_id}\nsequence: {sequence}\n"
                f"alleles ({len(alleles)}): {','.join(alleles)}"
            )
            for row in _coverage_for(label, [(sequence, alleles)], scratch):
                out_rows.append(
                    {
                        "consensus_id": consensus_id,
                        "protein_label": protein,
                        "population": row["population"],
                        "mhc_class": row["mhc_class"],
                        "coverage_percent": row["coverage_percent"],
                        "average_hit": row.get("average_hit", ""),
                        "pc90": row.get("pc90", ""),
                        "n_alleles": len(alleles),
                        "alleles": "|".join(alleles),
                    }
                )

        scratch = os.path.join(scratch_root, "all")
        os.makedirs(scratch)
        label = f"ALL EPITOPES TOGETHER ({len(entries)} peptides, {len(all_alleles)} alleles)"
        set_rows = _coverage_for(
            label, [(sequence, alleles) for _, sequence, alleles in entries], scratch
        )
        for row in set_rows:
            out_rows.append(
                {
                    "consensus_id": "ALL",
                    "protein_label": protein,
                    "population": row["population"],
                    "mhc_class": row["mhc_class"],
                    "coverage_percent": row["coverage_percent"],
                    "average_hit": row.get("average_hit", ""),
                    "pc90": row.get("pc90", ""),
                    "n_alleles": len(all_alleles),
                    "alleles": "|".join(all_alleles),
                }
            )

    with open(args.raw, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n\n".join(raw_chunks))

    tables.write_table(args.table, schemas.POPULATION_COVERAGE, out_rows)

    classes = sorted({str(row["mhc_class"]) for row in out_rows})
    lines = [
        f"Peptides analysed individually: **{len(entries)}**  ",
        f"Distinct alleles across all peptides: **{len(all_alleles)}**  ",
        f"MHC classes reported: `{'`, `'.join(classes)}`",
        "",
        "Each peptide is submitted to the IEDB tool on its own, so each row below",
        "is that peptide's own coverage. `ALL` is the whole set submitted together:",
        "the fraction responding to at least one peptide, which is the figure for a",
        "construct carrying all of them.",
        "",
        "| Peptide | Population | Class | Coverage | Average hit | PC90 |",
        "| --- | --- | --- | ---: | ---: | ---: |",
    ]
    for row in out_rows:
        lines.append(
            f"| {row['consensus_id']} | {row['population']} | {row['mhc_class']} | "
            f"{row['coverage_percent']}% | {row['average_hit']} | {row['pc90']} |"
        )

    # The same figures again, peptides down and regions across, one table per
    # class. The long form above is the complete record - it also carries
    # average_hit and pc90 - but it runs to peptides x regions x classes rows,
    # which is already 54 for two peptides. Comparing candidates across regions
    # is the question this stage exists to answer, and that is a read down a
    # column, not a scan through a list.
    #
    # Columns come from the rows rather than from the configured population
    # list so that a name the tool spells differently still lands in a column
    # instead of silently leaving every cell blank.
    column_order = []
    for row in out_rows:
        if row["population"] not in column_order:
            column_order.append(row["population"])
    row_order = [consensus_id for consensus_id, _, _ in entries] + ["ALL"]

    lines += ["", "## Coverage by region", ""]
    for mhc in classes:
        cells = {
            (row["consensus_id"], row["population"]): row["coverage_percent"]
            for row in out_rows
            if row["mhc_class"] == mhc
        }
        lines += [
            f"**MHC class {mhc}**",
            "",
            "| Peptide | " + " | ".join(column_order) + " |",
            "| --- |" + " ---: |" * len(column_order),
        ]
        for consensus_id in row_order:
            # A missing cell means the tool emitted no block for this class,
            # i.e. the peptide has no alleles of it - not zero coverage.
            values = [
                f"{cells[(consensus_id, population)]}%"
                if (consensus_id, population) in cells
                else "-"
                for population in column_order
            ]
            label = "**ALL together**" if consensus_id == "ALL" else consensus_id
            lines.append(f"| {label} | " + " | ".join(values) + " |")
        lines.append("")

    lines += [
        "`-` means the tool reported no block for that class, which happens when "
        "the peptide has no alleles of it. It does not mean zero coverage.",
    ]

    _common.write_summary(
        args.summary, f"Stage 8 - {protein} population coverage", lines
    )
    print(
        f"{protein}: coverage computed for {len(entries)} peptides individually "
        f"plus the combined set, over {len(populations)} populations."
    )
