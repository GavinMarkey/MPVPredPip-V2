"""Adapters for the IEDB standalone tools.

Covers stages 4, 5 and 8: MHC class I and class II binding prediction, class I
immunogenicity, epitope cluster analysis and population coverage.

------------------------------------------------------------------------------
On Tepitool
------------------------------------------------------------------------------
Tepitool is a guided *web* front end; it has no standalone distribution. What it
actually does is call the IEDB binding predictors with a curated allele panel and
a percentile cut-off. This module reproduces that locally: the allele panels live
in ``workflow/config/alleles/`` and the cut-offs in ``config.yaml``
(``tcell.mhc_i.percentile_cutoff`` and ``tcell.mhc_ii.percentile_cutoff``), so
the selection is explicit and versioned instead of being a hidden default on a
web form. Adjust the panels to match whatever Tepitool preset the study needs.

------------------------------------------------------------------------------
ADAPTER CONTRACT
------------------------------------------------------------------------------
The IEDB standalones have changed their CLI and column names across releases.
Every parser below resolves columns **by name** and every invocation keeps the
tool's raw stdout on disk. If a call signature does not match your installed
version, fix it in the ``_command_*`` helper for that tool - the parsing and all
downstream logic are independent of it.

Native outputs, tab-separated with a header. The binding ones are verified
against the installed tools - class I from its source, class II from the
expected values in its own configure.py test suite:

  MHC-I          allele, seq_num, start, end, length, peptide, core, icore,
                 netmhcpan_el    Score, rank
  MHC-II         allele, seq_num, start, end, length, core_peptide, peptide,
                 netmhciipan_el  ic50, percentile_rank, adjusted_rank
  immunogenicity peptide, length, score
  pop. coverage  population, coverage, average_hit, pc90

Three things there are easy to get wrong:

* the percentile column is ``rank`` for class I but ``percentile_rank`` for
  class II. Both resolve, but only because ``rank`` is in the candidate list;
* both have a ``peptide`` and class II also has a ``core_peptide``. Exact-match
  runs before prefix and substring in base.find_column, which is what keeps the
  binding peptide from being read out of the core column - do not "simplify"
  that ordering;
* for the EL (eluted-ligand) methods the ``ic50`` column is **0.00 on every
  row** - see the expected values in mhc_ii/configure.py. The eluted-ligand
  score has no IC50 to report, so the column is retained but empty of meaning.
  Filtering and ranking must use the percentile, which is what
  ``percentile_cutoff`` in config.yaml does. A score of 0 from this column is
  not a weak binder; it is no information at all.

  cluster        *** THIS ONE IS WRONG - see below. ***
                 Cluster.Sub-Cluster Number, Peptide Number, Alignment, Position,
                 Description, Peptide

The cluster entry above describes the IEDB **web** tool. The standalone
(IEDB_Cluster-1.0) emits something else entirely, verified by reading the
installed source:

    cluster_number, num_epitopes_in_cluster, epitope_num, epitope_name, epitope_seq

There is no consensus row, no alignment and no position - the whole package is
two files, and the words "consensus", "alignment" and "cluster-break" do not
appear anywhere in it. It only groups peptides. Its CLI is also not the one
``run_cluster`` builds: the input is positional, there is no -m/-i/-f, and
``--threshold`` is on a **0-100 scale** (default 80), not the 0-1 fraction in
config.yaml. It is Python 2 source and will not import under Python 3.

This matters beyond the adapter: stage 5 exists to produce consensus peptides,
and stages 6-9 screen and rank them. That consensus has to come from somewhere
the standalone does not provide.
"""

from __future__ import annotations

import glob
import os
from collections.abc import Sequence

from . import base

# --- column candidates -------------------------------------------------------
_ALLELE_COLUMNS = ["allele", "mhc", "Allele"]
_PEPTIDE_COLUMNS = ["peptide", "sequence", "Peptide", "peptide sequence"]
_START_COLUMNS = ["start", "seq_start", "peptide start", "Peptide start"]
_END_COLUMNS = ["end", "seq_end", "peptide end", "Peptide end"]
_SEQNUM_COLUMNS = ["seq_num", "seq #", "sequence number", "Seq #"]
_PERCENTILE_COLUMNS = [
    "percentile_rank", "percentile rank", "rank", "percentile",
    "el_rank", "ann_percentile", "consensus_percentile_rank",
]
_SCORE_COLUMNS = ["score", "el_score", "ic50", "affinity", "Score"]


def _entry(tool_dir: str, names: Sequence[str], what: str) -> str:
    """Locate a tool's entry point. Always returns an **absolute** path.

    Absolute is not a nicety. Every IEDB tool is invoked with ``cwd=tool_dir``,
    because each one resolves its own bundled data relative to where it is run
    from. A path like ``tools/mhc_i/src/predict_binding.py`` is relative to the
    project root, so handing it to a process started in ``tools/mhc_i`` makes it
    resolve against that directory a second time and the file is not found.

    This exact mistake has now been made three times in this codebase - IAPred,
    BepiPred and here - because the symptom is a "No such file or directory" for
    a path that plainly exists. Resolving it here means no caller has to
    remember.
    """
    for name in names:
        direct = os.path.join(tool_dir, name)
        if os.path.exists(direct):
            return os.path.abspath(direct)
        matches = sorted(glob.glob(os.path.join(tool_dir, "**", name), recursive=True))
        if matches:
            return os.path.abspath(matches[0])
    raise base.ToolMissing(
        f"{what} not found under {tool_dir}.\n  Looked for: {list(names)}\n\n"
        "IEDB standalone tools are downloaded per-tool; see tools/README.md."
    )


def read_alleles(path: str) -> list[str]:
    """Read an allele panel: one allele per line, ``#`` comments and blanks ignored."""
    base.require_path(path, "Allele panel")
    alleles: list[str] = []
    with open(path, encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.split("#", 1)[0].strip()
            if line:
                alleles.append(line)
    if not alleles:
        raise base.ToolError(f"Allele panel {path} is empty.")
    return alleles


# ---------------------------------------------------------------------------
# MHC binding
# ---------------------------------------------------------------------------

def predict_binding(
    fasta_path: str,
    out_path: str,
    tool_dir: str,
    alleles: Sequence[str],
    lengths: Sequence[int],
    method: str,
    mhc_class: str,
    python_executable: str = "python",
) -> str:
    """Run an IEDB binding predictor; returns the path to the raw output.

    The IEDB CLI takes allele and length as *parallel comma-separated lists*, so
    every (allele, length) combination has to be expanded explicitly - passing
    one allele and three lengths silently predicts only the first pairing.
    """
    base.require_path(tool_dir, f"IEDB MHC class {mhc_class} standalone")

    if mhc_class.upper() in {"I", "1"}:
        script = _entry(tool_dir, ["src/predict_binding.py", "predict_binding.py"],
                        "IEDB MHC-I predict_binding.py")
    else:
        script = _entry(tool_dir, ["mhc_II_binding.py", "src/mhc_II_binding.py",
                                   "predict_binding.py"],
                        "IEDB MHC-II binding script")

    allele_field = ",".join(a for a in alleles for _ in lengths)
    length_field = ",".join(str(length) for _ in alleles for length in lengths)

    command = [python_executable, script, method, allele_field]
    # The class I CLI takes lengths as a separate positional argument; class II
    # infers length from the method for most predictors but accepts it when given.
    if mhc_class.upper() in {"I", "1"}:
        command += [length_field, os.path.abspath(fasta_path)]
    else:
        command += [os.path.abspath(fasta_path), length_field]

    base.run_command(command, stdout_path=out_path, cwd=tool_dir)
    _check_binding_output(out_path, mhc_class, command)
    return out_path


def _check_binding_output(path: str, mhc_class: str, command: Sequence[str]) -> None:
    """Fail loudly when a predictor exits 0 without producing a table.

    The IEDB binding tools catch their own import errors, print the message to
    **stdout**, and exit 0. Stdout is where this adapter puts the results, so a
    broken install writes a one-line file that parses as a table with no data
    rows - which is indistinguishable, downstream, from "no peptide passed the
    cut-off".

    That is not hypothetical. With ``python=3.8`` in envs/iedb.yaml, class I
    failed on ``from importlib.resources import files`` and the entire stage 4
    run completed green, reporting 0 class I binders and merging 4231 class II
    epitopes with no class I support at all. Nothing in the run said the
    predictor had never executed; the only clue was that it finished in two
    seconds.

    So the header is checked rather than trusted. An empty result with a valid
    header is left alone - that really is a legitimate "nothing passed".
    """
    with open(path, encoding="utf-8", errors="replace") as handle:
        first = handle.readline().strip()
        second = handle.readline().strip()

    columns = {part.strip().lower() for part in first.split("\t")}
    if {"allele", "peptide"} <= columns:
        return

    printable = " ".join(str(part) for part in command)
    raise base.ToolError(
        f"The MHC class {mhc_class} predictor exited successfully but did not write a "
        f"results table to {path}.\n\n"
        f"  first line:  {first[:300] or '(the file is empty)'}\n"
        f"  second line: {second[:300] or '(there is no second line)'}\n\n"
        "A binding table must have at least 'allele' and 'peptide' columns. These "
        "tools report their own startup failures on stdout and still exit 0, so what "
        "is above is most likely the real error - an import failure, a missing "
        "bundled binary, or an allele the method does not support.\n\n"
        f"  command: {printable}\n\n"
        "This is raised rather than returning zero binders because the two are "
        "impossible to tell apart downstream, and zero binders would quietly empty "
        "stage 5."
    )


def parse_binding(
    path: str,
    percentile_cutoff: float,
    id_by_seqnum: dict[int, str] | None = None,
) -> list[dict[str, object]]:
    """Parse binding output and keep rows at or below *percentile_cutoff*."""
    rows = base.read_delimited(path)
    if not rows:
        return []

    header = list(rows[0].keys())
    c_allele = base.find_column(header, _ALLELE_COLUMNS, "allele", path)
    c_peptide = base.find_column(header, _PEPTIDE_COLUMNS, "peptide", path)
    c_start = base.find_column(header, _START_COLUMNS, "start", path, required=False)
    c_end = base.find_column(header, _END_COLUMNS, "end", path, required=False)
    c_seqnum = base.find_column(header, _SEQNUM_COLUMNS, "sequence number", path, required=False)
    c_pct = base.find_column(header, _PERCENTILE_COLUMNS, "percentile rank", path)
    c_score = base.find_column(header, _SCORE_COLUMNS, "score", path, required=False)

    kept: list[dict[str, object]] = []
    for row in rows:
        percentile = base.to_float(row[c_pct], "percentile rank")
        if percentile > percentile_cutoff:
            continue
        peptide = str(row[c_peptide]).strip().upper()
        if not peptide:
            continue

        start = int(base.to_float(row[c_start], "start")) if c_start and row.get(c_start) else 0
        end = (
            int(base.to_float(row[c_end], "end"))
            if c_end and row.get(c_end)
            else (start + len(peptide) - 1 if start else 0)
        )
        seq_num = int(base.to_float(row[c_seqnum], "seq_num")) if c_seqnum and row.get(c_seqnum) else 1

        kept.append(
            {
                "allele": str(row[c_allele]).strip(),
                "peptide": peptide,
                "start": start,
                "end": end,
                "seq_num": seq_num,
                "sequence_id": (id_by_seqnum or {}).get(seq_num, ""),
                "percentile": percentile,
                "score": base.to_float(row[c_score], "score") if c_score and row.get(c_score) else percentile,
            }
        )
    return kept


# ---------------------------------------------------------------------------
# class I immunogenicity
# ---------------------------------------------------------------------------

def predict_immunogenicity(
    peptides: Sequence[str],
    out_path: str,
    tool_dir: str,
    python_executable: str = "python",
) -> dict[str, float]:
    """Score peptides with the IEDB class I immunogenicity predictor."""
    base.require_path(tool_dir, "IEDB immunogenicity standalone")
    script = _entry(
        tool_dir,
        ["predict_immunogenicity.py", "src/predict_immunogenicity.py", "immunogenicity.py"],
        "IEDB immunogenicity script",
    )

    unique = sorted({peptide.strip().upper() for peptide in peptides if peptide.strip()})
    if not unique:
        return {}

    peptide_file = os.path.join(os.path.dirname(os.path.abspath(out_path)), "immunogenicity_input.txt")
    os.makedirs(os.path.dirname(peptide_file), exist_ok=True)
    with open(peptide_file, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(unique) + "\n")

    base.run_command([python_executable, script, peptide_file], stdout_path=out_path, cwd=tool_dir)

    # Two things about this tool's output that the others do not do:
    #
    #   * it prints its settings first - "masking: default", "masked variables:
    #     ..." and a blank line - so the table does not start at line 1;
    #   * the table is comma-separated, even though this adapter stores it with a
    #     .tsv extension. That extension is ours, not the tool's; the delimiter
    #     is sniffed, so it does not matter, but do not "fix" the content to
    #     match the name.
    #
    # header_contains pins the real header row and fails loudly if there is no
    # table at all, rather than reading the preamble as a header.
    rows = base.read_delimited(out_path, header_contains="peptide")
    if not rows:
        raise base.ToolError(f"Immunogenicity output {out_path} is empty.")
    header = list(rows[0].keys())
    c_peptide = base.find_column(header, _PEPTIDE_COLUMNS, "peptide", out_path)
    c_score = base.find_column(header, ["score", "immunogenicity", "Score"], "score", out_path)

    scores: dict[str, float] = {}
    for row in rows:
        peptide = str(row[c_peptide]).strip().upper()
        if peptide:
            scores[peptide] = base.to_float(row[c_score], "immunogenicity score")

    missing = [peptide for peptide in unique if peptide not in scores]
    if missing:
        raise base.ToolError(
            f"The immunogenicity predictor returned no score for {len(missing)} of "
            f"{len(unique)} peptides (first few: {missing[:5]}).\n"
            "Scoring every class I binder is required before the positive-score filter "
            "can be applied, so this is treated as an error rather than dropping them."
        )
    return scores


# ---------------------------------------------------------------------------
# epitope clustering - deliberately absent
# ---------------------------------------------------------------------------
#
# There was a ``run_cluster`` here. It has been removed rather than corrected,
# because there is nothing left for it to do: stage 5 clusters in-process with
# vaxpipe.cluster and never invokes the standalone. Leaving it would have left a
# function whose command line matched no released version of the tool, which is
# worse than no function at all - see this module's docstring for what the
# standalone actually is, and src/vaxpipe/cluster.py for what replaced it.


# ---------------------------------------------------------------------------
# population coverage
# ---------------------------------------------------------------------------

def run_population_coverage(
    input_path: str,
    out_path: str,
    tool_dir: str,
    populations: Sequence[str],
    mhc_class: str = "I,II",
    python_executable: str = "python",
) -> str:
    """Run the IEDB population coverage tool.

    *input_path* must be a tab-separated file of ``peptide<TAB>allele[,allele...]``
    - the format the tool expects. :func:`write_population_input` builds it.
    """
    base.require_path(tool_dir, "IEDB population coverage standalone")
    script = _entry(
        tool_dir,
        ["calculate_population_coverage.py", "src/calculate_population_coverage.py"],
        "IEDB population coverage script",
    )
    command = [
        python_executable, script,
        "-p", ",".join(populations),
        "-c", mhc_class,
        "-f", os.path.abspath(input_path),
    ]
    base.run_command(command, stdout_path=out_path, cwd=tool_dir)
    return out_path


def write_population_input(path: str, entries: Sequence[tuple[str, Sequence[str]]]) -> str:
    """Write the ``peptide<TAB>comma-separated alleles`` file the tool expects."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for peptide, alleles in entries:
            unique = sorted({allele for allele in alleles if allele.strip()})
            if unique:
                handle.write(f"{peptide}\t{','.join(unique)}\n")
    return path


def parse_population_coverage(path: str) -> list[dict[str, object]]:
    """Parse population-coverage output into ``population``/``coverage`` rows.

    The tool's output is a short human-readable report rather than a clean table,
    so this is a line scanner: it looks for the section header and then reads
    whitespace- or tab-separated rows until the block ends.
    """
    if not os.path.exists(path):
        raise base.ToolError(f"Population coverage output not found: {path}")

    rows: list[dict[str, object]] = []
    with open(path, encoding="utf-8-sig", errors="replace") as handle:
        lines = [line.rstrip("\n") for line in handle]

    in_table = False
    for line in lines:
        stripped = line.strip()
        if not stripped:
            in_table = False
            continue

        lowered = stripped.lower()
        if lowered.startswith("population") and "coverage" in lowered:
            in_table = True
            continue
        if not in_table:
            continue

        parts = [part.strip() for part in (stripped.split("\t") if "\t" in stripped else stripped.split("  ")) if part.strip()]
        if len(parts) < 2:
            continue

        coverage = parts[1].replace("%", "").strip()
        try:
            coverage_value = float(coverage)
        except ValueError:
            continue

        rows.append(
            {
                "population": parts[0],
                "coverage_percent": coverage_value,
                "average_hit": base.to_float(parts[2], "average hit") if len(parts) > 2 and _numeric(parts[2]) else "",
                "pc90": base.to_float(parts[3], "pc90") if len(parts) > 3 and _numeric(parts[3]) else "",
            }
        )

    if not rows:
        raise base.ToolError(
            f"Could not parse any population rows from {path}.\n"
            "The raw report has been kept - check its layout and adjust "
            "parse_population_coverage() in src/vaxpipe/adapters/iedb.py."
        )
    return rows


def _numeric(text: str) -> bool:
    try:
        float(text)
        return True
    except ValueError:
        return False
