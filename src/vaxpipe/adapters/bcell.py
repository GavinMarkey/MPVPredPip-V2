"""Stage 3 adapters - linear B-cell epitope prediction (BepiPred and EpiDope).

Both predictors score **every residue** of a protein. Neither emits epitopes
directly in a form this pipeline can use, so each adapter normalises its native
output into a per-residue score vector and hands that to
:func:`vaxpipe.epitopes.residues_to_epitopes`, which applies the thresholds and
length rules from ``config.yaml``. Keeping the interval logic in one shared place
means BepiPred and EpiDope epitopes are built by identical rules and stay
directly comparable in stage 5.

------------------------------------------------------------------------------
ADAPTER CONTRACT - read this before changing either parser
------------------------------------------------------------------------------
``predict_*`` must return ``{sequence_id: [score_per_residue, ...]}`` where the
list length equals the ungapped protein length and index 0 is residue 1.

Both parsers resolve columns **by name**, accept CSV or TSV, and search a list of
candidate output filenames. They were written against the output layouts
described below. If your installed version differs, the raw output is preserved
under the stage's ``raw/`` directory - adjust the candidate lists and column
names at the top of each function, not the interval logic.

  BepiPred-3.0  VERIFIED against the upstream source (bp3/bepipred3.py,
                ``create_csvfile``). Writes ``raw_output.csv``::

                    Accession,Residue,BepiPred-3.0 score,BepiPred-3.0 linear epitope score
                    Spike|Zaire_ebolavirus|NP_066246.1,M,0.0871, 0.1043

                There is **no position column** - one row per residue, in
                sequence order, grouped by accession. ``Residue`` is the amino
                acid letter, not an index. The parser handles this: no candidate
                resolves to a position column, so it falls back to file order.
  EpiDope       VERIFIED by running it. A *sectioned* tab-separated file, not a
                flat table - see :func:`parse_epidope`.
"""

from __future__ import annotations

import glob
import os
from collections import defaultdict

from . import base

# ---------------------------------------------------------------------------
# BepiPred
# ---------------------------------------------------------------------------

_BEPIPRED_OUTPUT_CANDIDATES = [
    "raw_output.csv",
    "Bepipred_raw_output.csv",
    "bepipred3_raw_output.csv",
    "*raw_output*.csv",
    "*residue*.csv",
    "*.csv",
]

_BEPIPRED_ID_COLUMNS = ["Accession", "Entry", "ID", "Protein", "Header", "seq_id", "Sequence ID"]

#: BepiPred emits no position column. These are kept for a future version that
#: adds one; none of them matches the real header, which is deliberate - the
#: "Residue" column holds the amino acid letter, and resolving it as a position
#: would feed "M" to float(). Verified safe: none of these is a prefix or
#: substring of "residue" under base._norm.
_BEPIPRED_POS_COLUMNS = ["Residue Position", "Position", "Pos", "Index", "Residue Index", "AA Position"]

#: "BepiPred-3.0 score" is the raw per-residue ensemble probability;
#: "BepiPred-3.0 linear epitope score" is that smoothed by a rolling mean
#: (default window 9). We take the RAW score, because BepiPred's own epitope call
#: thresholds the raw score too - bp3_pred_variable_threshold compares
#: ``avg_prob`` against ``var_threshold``, not the rolling mean. So our 0.1512 in
#: config.yaml means the same thing BepiPred's -t means, and our epitopes agree
#: with the Bcell_epitope_preds.fasta kept alongside them. The exact-match rule
#: in base.find_column picks the raw column even though both start with the same
#: words; do not reorder this list.
_BEPIPRED_SCORE_COLUMNS = [
    "BepiPred-3.0 score",
    "BepiPred-3.0 linear epitope score",
    "score",
    "Score",
    "Prediction",
    "epitope score",
]


def predict_bepipred(
    fasta_path: str,
    out_dir: str,
    tool_dir: str,
    threshold: float | None = None,
    esm_dir: str = "",
    python_executable: str = "python",
    extra_args: list[str] | None = None,
) -> dict[str, list[float]]:
    """Run BepiPred and return per-residue scores keyed by sequence ID.

    *threshold* is passed through as ``-t``. It does not affect ``raw_output.csv``,
    which is what this pipeline reads - it decides the epitope call in BepiPred's
    own ``Bcell_epitope_preds.fasta``. Passing it keeps the preserved raw output
    consistent with the analysis built from it, for the same reason the EpiDope
    rule passes its slice length and shift.

    *esm_dir* is where the ESM-2 encodings are cached. BepiPred defaults to a
    directory beside its own source, which would write into ``tools/``; the
    workflow points it somewhere the encodings survive a re-run instead, since
    producing them is the slow part.
    """
    base.require_path(
        tool_dir,
        "BepiPred installation",
        "BepiPred-3.0 is free for academic use and installs unattended:\n"
        "    bash setup/install_bepipred.sh",
    )
    os.makedirs(out_dir, exist_ok=True)

    entry = _find_entrypoint(
        tool_dir,
        ["bepipred3_CLI.py", "BepiPred3_CLI.py", "bepipred3.py", "CLI.py"],
        "BepiPred command-line entry point",
    )

    command = [
        python_executable,
        # Absolute, because the command runs with cwd=tool_dir so that the local
        # bp3 package - which carries the trained model weights in
        # bp3/BP3Models - is the one imported. A relative path would be resolved
        # against tool_dir a second time and not be found.
        os.path.abspath(entry),
        "-i", os.path.abspath(fasta_path),
        "-o", os.path.abspath(out_dir),
    ]
    if esm_dir:
        # expanduser: the configured default is under ~, because the cache has to
        # be outside the workspace (see config.yaml) and the run user does not
        # own /opt.
        esm_dir = os.path.abspath(os.path.expanduser(esm_dir))
        os.makedirs(esm_dir, exist_ok=True)
        command += ["-esm_dir", esm_dir]
    if threshold is not None:
        command += ["-t", str(threshold)]
    # -pred is required by BepiPred's own argument parser. "vt_pred" thresholds
    # the averaged ensemble probability, which is the column we read.
    command += extra_args or ["-pred", "vt_pred"]
    base.run_command(command, cwd=tool_dir)

    return parse_bepipred(out_dir)


def parse_bepipred(out_dir: str) -> dict[str, list[float]]:
    """Read BepiPred's own output from *out_dir*, without running it."""
    return _parse_per_residue(
        out_dir,
        _BEPIPRED_OUTPUT_CANDIDATES,
        _BEPIPRED_ID_COLUMNS,
        _BEPIPRED_POS_COLUMNS,
        _BEPIPRED_SCORE_COLUMNS,
        tool="BepiPred",
    )


# ---------------------------------------------------------------------------
# EpiDope
# ---------------------------------------------------------------------------

#: ``epidope_scores.csv`` is the per-residue file and must come first. The others
#: EpiDope writes - ``predicted_epitopes.csv`` and ``predicted_epitopes_sliced.faa``
#: - are epitope-level, and reading one of those as per-residue scores would
#: silently produce nonsense. The loose patterns are kept only as a fallback for
#: a version that renames things.
_EPIDOPE_OUTPUT_CANDIDATES = [
    "epidope_scores.csv",
    "*scores*.csv",
    "*score*.tsv",
]

_EPIDOPE_ID_COLUMNS = ["#Header", "Header", "protein", "Protein", "ID", "seq_id", "Accession"]
_EPIDOPE_POS_COLUMNS = ["position", "Position", "pos", "index", "residue_position"]
_EPIDOPE_SCORE_COLUMNS = ["score", "Score", "epidope_score", "value", "prediction"]


def predict_epidope(
    fasta_path: str,
    out_dir: str,
    executable: str = "epidope",
    threads: int = 4,
    extra_args: list[str] | None = None,
) -> dict[str, list[float]]:
    """Run EpiDope and return per-residue scores keyed by sequence ID."""
    os.makedirs(out_dir, exist_ok=True)
    command = [
        executable,
        "-i", os.path.abspath(fasta_path),
        "-o", os.path.abspath(out_dir),
    ]
    if threads:
        command += ["-p", str(threads)]
    command += extra_args or []
    base.run_command(command)

    return parse_epidope(out_dir)


def parse_epidope(out_dir: str) -> dict[str, list[float]]:
    """Read EpiDope's own output from *out_dir*, without running it.

    Kept separate from :func:`predict_epidope` because EpiDope's environment pins
    **Python 3.6**, which cannot even parse this package (``from __future__ import
    annotations`` requires 3.7). The workflow therefore runs the tool in its own
    environment with a bare shell command and parses the result here, under the
    orchestrator's environment.

    ``epidope_scores.csv`` is *sectioned*, not a flat table, and is tab-separated
    despite the extension::

        position/header	aminoacid	score
        >Spike|Bundibugyo_virus|AYI50316.1
        1	M	0.5627997517585754
        2	V	0.5722602605819702
        >Spike|Sudan_ebolavirus|WEY06878.1
        1	M	0.5301142334938049

    A ``>`` line opens a block; the rows beneath it belong to that sequence. The
    deflines survive intact here, pipes included - note that the per-sequence
    files EpiDope writes alongside it have the pipes stripped from their names.
    """
    path = _locate_output(out_dir, _EPIDOPE_OUTPUT_CANDIDATES, "EpiDope")

    ordered: dict[str, list[tuple[int, float]]] = {}
    current: str | None = None
    sectioned = False

    with open(path, encoding="utf-8") as handle:
        for number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue

            if line.startswith(">"):
                sectioned = True
                current = line[1:].strip().split()[0]
                ordered.setdefault(current, [])
                continue

            fields = line.split("\t") if "\t" in line else line.split(",")
            if len(fields) < 3:
                continue

            position_text = fields[0].strip()
            if not position_text.lstrip("-").isdigit():
                continue  # the column-header line

            if current is None:
                raise base.ToolError(
                    f"{path} line {number}: residue scores before any '>' header. "
                    "The file layout is not the one this parser expects."
                )
            ordered[current].append(
                (int(position_text), base.to_float(fields[2], "EpiDope score"))
            )

    if not sectioned:
        # A version that emits a flat table with an identifier column instead.
        return _parse_per_residue(
            out_dir,
            _EPIDOPE_OUTPUT_CANDIDATES,
            _EPIDOPE_ID_COLUMNS,
            _EPIDOPE_POS_COLUMNS,
            _EPIDOPE_SCORE_COLUMNS,
            tool="EpiDope",
        )

    if not ordered:
        raise base.ToolError(f"No per-residue scores found in {path}.")

    scores: dict[str, list[float]] = {}
    for sequence_id, pairs in ordered.items():
        pairs.sort(key=lambda pair: pair[0])
        scores[sequence_id] = [score for _, score in pairs]
    return scores


# ---------------------------------------------------------------------------
# shared parsing
# ---------------------------------------------------------------------------

def _find_entrypoint(tool_dir: str, names: list[str], what: str) -> str:
    """Locate a tool's entry point. Always returns an **absolute** path.

    These tools run with ``cwd=tool_dir`` so they can find their own bundled
    data, which means a project-relative path would be resolved against that
    directory a second time. See the same note on ``_entry`` in iedb.py - this
    has caught three adapters now.
    """
    for name in names:
        candidate = os.path.join(tool_dir, name)
        if os.path.exists(candidate):
            return os.path.abspath(candidate)
        matches = sorted(glob.glob(os.path.join(tool_dir, "**", name), recursive=True))
        if matches:
            return os.path.abspath(matches[0])
    raise base.ToolMissing(
        f"{what} not found under {tool_dir}.\n  Looked for: {names}\n\n"
        "See tools/README.md for the expected layout."
    )


def _locate_output(out_dir: str, candidates: list[str], tool: str) -> str:
    for pattern in candidates:
        matches = sorted(glob.glob(os.path.join(out_dir, "**", pattern), recursive=True))
        # Ignore anything this pipeline wrote itself.
        matches = [m for m in matches if os.path.basename(m) != "normalised.tsv"]
        if matches:
            return matches[0]
    listing = "\n".join(
        f"    {os.path.relpath(p, out_dir)}"
        for p in sorted(glob.glob(os.path.join(out_dir, "**", "*"), recursive=True))
        if os.path.isfile(p)
    ) or "    (the tool produced no files)"
    raise base.ToolError(
        f"Could not find the {tool} per-residue score file in {out_dir}.\n"
        f"  Looked for: {candidates}\n  Files present:\n{listing}\n\n"
        f"Update the candidate list in src/vaxpipe/adapters/bcell.py for your {tool} version."
    )


def _parse_per_residue(
    out_dir: str,
    candidates: list[str],
    id_columns: list[str],
    pos_columns: list[str],
    score_columns: list[str],
    tool: str,
) -> dict[str, list[float]]:
    path = _locate_output(out_dir, candidates, tool)
    rows = base.read_delimited(path)
    if not rows:
        raise base.ToolError(f"{tool} output {path} contains no data rows.")

    header = list(rows[0].keys())
    id_column = base.find_column(header, id_columns, f"{tool} sequence-ID", path)
    score_column = base.find_column(header, score_columns, f"{tool} score", path)
    pos_column = base.find_column(header, pos_columns, f"{tool} position", path, required=False)

    ordered: dict[str, list[tuple[int, float]]] = defaultdict(list)
    for index, row in enumerate(rows):
        sequence_id = str(row[id_column]).strip().lstrip(">").split()[0]
        score = base.to_float(row[score_column], f"{tool} score")
        position = base.to_float(row[pos_column], "position") if pos_column else float(index)
        ordered[sequence_id].append((int(position), score))

    scores: dict[str, list[float]] = {}
    for sequence_id, pairs in ordered.items():
        # Sorting by the reported position protects against a tool that emits
        # residues out of order; without a position column the file order is the
        # only ordering available and is preserved.
        pairs.sort(key=lambda pair: pair[0])
        scores[sequence_id] = [score for _, score in pairs]
    return scores


def load_precomputed(path: str) -> dict[str, list[float]]:
    """Read a per-residue score table this pipeline wrote earlier.

    Lets stage 3 be re-run from cached predictions, and gives a documented escape
    hatch for plugging in any other per-residue B-cell predictor: emit a TSV with
    ``sequence_id``, ``position`` and ``score`` and point the adapter at it.
    """
    rows = base.read_delimited(path)
    ordered: dict[str, list[tuple[int, float]]] = defaultdict(list)
    for row in rows:
        ordered[str(row["sequence_id"]).strip()].append(
            (int(float(row["position"])), float(row["score"]))
        )
    return {
        sequence_id: [score for _, score in sorted(pairs)]
        for sequence_id, pairs in ordered.items()
    }


def save_per_residue(path: str, scores: dict[str, list[float]]) -> None:
    """Persist per-residue scores so stage 3 is reproducible without re-running."""
    from .. import tables

    tables.write_table(
        path,
        ["sequence_id", "position", "score"],
        (
            {"sequence_id": sequence_id, "position": index + 1, "score": score}
            for sequence_id, vector in sorted(scores.items())
            for index, score in enumerate(vector)
        ),
    )
