"""Stage 2 and 6 adapters - antigenicity, allergenicity and autoimmunity.

------------------------------------------------------------------------------
Antigenicity predictor backends
------------------------------------------------------------------------------
``antigenicity.tool`` in ``config.yaml`` selects one of three backends:

  ``iapred``    IAPred from https://github.com/sebamiles/IAPred, installed into
                ``tools/iapred`` by ``setup/install_iapred.sh``. This is the
                default and the one the pipeline specification calls for.

                The invocation in :func:`_predict_iapred` was originally inferred
                without access to the repository, and has since been corrected
                against a real run: the entry point is ``IApred.py``, the CLI is
                positional (``IApred.py in.fasta [out.csv]``) and the score column
                is ``Intrinsic_Antigenicity_Score``. Note that IApred's own README
                is stale on that last point - it documents ``IAscore``.

                It is still resolved leniently - any ``*.py`` entry point is
                discovered, several argument styles are attempted, and columns are
                matched by name - and the raw output is always preserved. If a
                future release changes the CLI, edit ``_IAPRED_*`` below, or fall
                back to ``external``, which needs no code change at all.

                *** IApred will not score sequences under 20 aa. *** It writes
                "Sequence too short" into the score column instead. Stage 6
                screens cluster consensus peptides, which are typically 20-50 aa,
                so most are scored - but short ones do occur and are recorded as
                an explicit screening failure, never as a low score. See
                :data:`_IAPRED_SENTINELS`.

  ``external``  run any command you configure (``antigenicity.command``),
                formatted with ``{input}`` and ``{output}``. Use this for IAPred
                if its real CLI turns out not to match, or for any other
                predictor. It must read a FASTA and write a delimited file with
                an identifier column and a score column; nothing else is assumed.

  ``kolaskar``  a built-in, dependency-free implementation of the
                Kolaskar & Tongaonkar (1990) semi-empirical antigenic-propensity
                method, so the workflow runs end to end before anything is
                installed. It is **not** equivalent to IAPred or VaxiJen and must
                not be reported as such.
"""

from __future__ import annotations

import os
from collections.abc import Sequence

from .. import fasta as fasta_mod
from . import base

# ---------------------------------------------------------------------------
# antigenicity - IAPred
# ---------------------------------------------------------------------------
# https://github.com/sebamiles/IAPred
#
# The entry point, argument order and column names below were checked against
# the repository's own README after installing it (tools/iapred/README.md). They
# have NOT been confirmed by an actual run - do that before trusting stage 2 or
# stage 6 output.
#
# Everything here is still used leniently - the entry point falls back to any
# *.py in the install, the flag styles are tried in order until one succeeds, and
# output columns are matched by name. Correct them here if the real CLI differs;
# nothing outside this module depends on them.
#
# Scores run roughly -3 to 3, NOT 0 to 1: negative values are legitimate and mean
# low antigenicity. The README's category boundaries are High > 0.3,
# Moderate -0.3 to 0.3, Low < -0.3.

#: Candidate entry-point filenames, most specific first. The real name is
#: ``IApred.py`` - lowercase "p" - which matters on a case-sensitive filesystem.
_IAPRED_ENTRYPOINTS = [
    "IApred.py", "iapred.py", "IAPred.py", "run_iapred.py", "main.py", "predict.py", "*.py",
]

#: Argument styles attempted in order. ``{input}``/``{output}`` are substituted.
#: Positional comes first: it is what the repository's README documents
#: (``python IApred.py input.fasta [output.csv]``). The flag styles are kept
#: behind it as fallbacks in case a future release grows a different CLI.
_IAPRED_ARG_STYLES = [
    ["{input}", "{output}"],
    ["-i", "{input}", "-o", "{output}"],
    ["--input", "{input}", "--output", "{output}"],
    ["-f", "{input}", "-o", "{output}"],
]

_IAPRED_ID_COLUMNS = ["id", "seq_id", "sequence_id", "name", "accession", "header", "protein", "query"]

#: IApred's CSV is ``Header, Sequence_Length, Intrinsic_Antigenicity_Score,
#: Antigenicity_Category`` - established by RUNNING it, not from its README,
#: which documents a stale ``IAscore`` column the code no longer emits.
#:
#: The real name comes first so it wins on the exact-match pass. "antigenicity"
#: is deliberately LAST: it is a prefix of ``Antigenicity_Category``, which holds
#: the text label (High/Moderate/Low), and resolving to that column instead of
#: the numeric one makes the score unparseable.
_IAPRED_SCORE_COLUMNS = [
    "intrinsic_antigenicity_score", "iascore", "ia_score", "score",
    "iapred", "iapred_score", "prediction", "value",
    "antigenicity",
]

#: IApred writes these strings INTO the score column in place of a number
#: (IApred.py lines 113, 118, 150): for sequences under 20 aa, for sequences with
#: no valid residues, and for anything that raised inside the tool.
#:
#: They are not scores and must not be parsed as one - but they must not abort
#: the stage either. Stage 6 screens cluster consensus peptides, which are
#: usually 20-50 aa, so most are scored; a single short one would otherwise take
#: down a whole protein's run. They become ``None`` here, which the caller
#: records as an explicit screening failure rather than a low score.
_IAPRED_SENTINELS = {"sequence too short", "invalid sequence", "error", "n/a"}

#: Recorded against a peptide IApred declined to score. The raw output is always
#: preserved next to the table, so the exact sentinel is recoverable.
IAPRED_UNSCORABLE = "IApred returned no score (sequence too short, invalid, or an internal error)"


def _predict_iapred(
    fasta_path: str,
    raw_path: str,
    tool_dir: str,
    python_executable: str = "python",
) -> dict[str, float | None]:
    """Run IAPred and return ``{sequence_id: score}``.

    Each argument style in :data:`_IAPRED_ARG_STYLES` is tried until one both
    succeeds and leaves a parseable file. If every style fails, the error lists
    what was attempted along with each tool's own stderr, so the real signature
    can be read off the failure rather than guessed at again.

    IApred truncates every defline to 20 characters (``IApred.py`` line 105), so
    the pipeline's own identifiers do not survive the round trip -
    ``Nucleoprotein|Bundibugyo`` is 24. The sequences are therefore written to a
    temporary FASTA under short surrogate names, and the scores mapped back here.
    A ``.ids.tsv`` sidecar records the mapping so the preserved raw output stays
    readable.
    """
    import glob as _glob

    base.require_path(
        tool_dir,
        "IAPred installation",
        "Install it with:  bash setup/install_iapred.sh\n"
        "(clones https://github.com/sebamiles/IAPred into tools/iapred)",
    )

    entry = None
    for pattern in _IAPRED_ENTRYPOINTS:
        matches = sorted(_glob.glob(os.path.join(tool_dir, pattern)))
        matches += sorted(_glob.glob(os.path.join(tool_dir, "**", pattern), recursive=True))
        # Ignore setup and test scripts, which would otherwise match "*.py".
        matches = [
            m for m in matches
            if not os.path.basename(m).startswith(("setup", "test_", "conftest", "_"))
        ]
        if matches:
            entry = matches[0]
            break

    if entry is None:
        raise base.ToolMissing(
            f"No IAPred entry script found under {tool_dir}.\n"
            f"  Looked for: {_IAPRED_ENTRYPOINTS}\n\n"
            "Set antigenicity.tool to 'external' and give the real command in "
            "antigenicity.command instead."
        )

    records = fasta_mod.read_fasta(fasta_path)
    if not records:
        raise base.ToolError(f"No sequences to score in {fasta_path}.")

    stem = os.path.splitext(os.path.abspath(raw_path))[0]
    surrogate_path = stem + ".input.fasta"
    surrogates = {f"s{index}": record for index, record in enumerate(records)}

    os.makedirs(os.path.dirname(stem) or ".", exist_ok=True)
    fasta_mod.write_fasta(
        surrogate_path, [(key, record.sequence) for key, record in surrogates.items()]
    )
    with open(stem + ".ids.tsv", "w", encoding="utf-8") as handle:
        handle.write("surrogate_id\tsequence_id\n")
        for key, record in surrogates.items():
            handle.write(f"{key}\t{record.id}\n")

    attempts: list[str] = []
    for style in _IAPRED_ARG_STYLES:
        args = [
            part.format(input=surrogate_path, output=os.path.abspath(raw_path))
            for part in style
        ]
        command = [python_executable, os.path.abspath(entry)] + args
        try:
            result = base.run_command(command, cwd=tool_dir, check=False)
        except base.ToolError as exc:
            attempts.append(f"  {' '.join(args)}\n    {exc}")
            continue

        if result.returncode != 0:
            attempts.append(
                f"  {' '.join(args)}\n    exit {result.returncode}: "
                f"{(result.stderr or '').strip()[:300]}"
            )
            continue

        # Some scripts write to stdout rather than to an -o path.
        if not os.path.exists(raw_path) and result.stdout:
            os.makedirs(os.path.dirname(os.path.abspath(raw_path)) or ".", exist_ok=True)
            with open(raw_path, "w", encoding="utf-8") as handle:
                handle.write(result.stdout)

        try:
            scored = _parse_iapred(raw_path)
        except base.ToolError as exc:
            attempts.append(f"  {' '.join(args)}\n    ran, but output unusable: {exc}")
            continue

        unknown = sorted(set(scored) - set(surrogates))
        if unknown:
            raise base.ToolError(
                f"IAPred returned identifiers this run did not submit: {unknown[:5]}.\n"
                f"Expected the surrogate names in {stem}.ids.tsv."
            )
        return {surrogates[key].id: value for key, value in scored.items()}

    raise base.ToolError(
        f"Could not drive IAPred at {entry}. Attempts:\n" + "\n".join(attempts) +
        "\n\nIf a release has changed the CLI, correct _IAPRED_ARG_STYLES in "
        "src/mpvpredpip/adapters/screening.py, or set antigenicity.tool: 'external' "
        "and put the real command in antigenicity.command."
    )


def _parse_iapred(raw_path: str) -> dict[str, float | None]:
    """Parse IAPred output, resolving the identifier and score columns by name.

    A value of ``None`` means IApred returned one of its :data:`_IAPRED_SENTINELS`
    for that sequence - it was *present but not scored*. That is deliberately
    distinct from an identifier being absent altogether, which indicates a
    plumbing fault (a truncated or mismatched ID) and is still an error.
    """
    rows = base.read_delimited(raw_path)
    if not rows:
        raise base.ToolError(f"no rows in {raw_path}")

    header = list(rows[0].keys())
    c_id = base.find_column(header, _IAPRED_ID_COLUMNS, "IAPred sequence-ID", raw_path)
    c_score = base.find_column(header, _IAPRED_SCORE_COLUMNS, "IAPred score", raw_path)

    scores: dict[str, float | None] = {}
    for row in rows:
        identifier = str(row.get(c_id, "")).strip()
        if not identifier:
            continue
        identifier = identifier.lstrip(">").split()[0]

        raw = str(row.get(c_score, "")).strip()
        if raw.lower() in _IAPRED_SENTINELS:
            scores[identifier] = None
        else:
            scores[identifier] = base.to_float(raw, "IAPred score")
    return scores


# ---------------------------------------------------------------------------
# antigenicity - built-in fallback
# ---------------------------------------------------------------------------

#: Kolaskar & Tongaonkar (1990) antigenic propensity values, one per residue.
#: Cross-check these against the original paper before publishing results.
KOLASKAR_TONGAONKAR = {
    "A": 1.064, "C": 1.412, "D": 0.866, "E": 0.851, "F": 1.091,
    "G": 0.874, "H": 1.105, "I": 1.152, "K": 0.930, "L": 1.250,
    "M": 0.826, "N": 0.776, "P": 1.064, "Q": 1.015, "R": 0.873,
    "S": 1.012, "T": 0.909, "V": 1.383, "W": 0.893, "Y": 1.161,
}


def antigenicity_kolaskar(sequence: str, window: int = 7) -> float:
    """Mean antigenic propensity of *sequence*.

    The published method smooths propensity over a sliding window and calls
    residues antigenic where the smoothed value exceeds the sequence average. For
    a whole-peptide score the pipeline needs one number, so this returns the mean
    propensity across the sequence; the window is applied first so that short
    peptides are treated the same way as long proteins.
    """
    residues = [character for character in sequence.upper() if character in KOLASKAR_TONGAONKAR]
    if not residues:
        return 0.0
    values = [KOLASKAR_TONGAONKAR[residue] for residue in residues]

    if window > 1 and len(values) >= window:
        half = window // 2
        smoothed = []
        for index in range(len(values)):
            lo = max(0, index - half)
            hi = min(len(values), index + half + 1)
            smoothed.append(sum(values[lo:hi]) / (hi - lo))
        values = smoothed

    return sum(values) / len(values)


def predict_antigenicity(
    fasta_path: str,
    out_path: str,
    threshold: float,
    tool: str = "iapred",
    command: str = "",
    tool_dir: str = "tools/iapred",
    id_columns: Sequence[str] = ("id", "sequence_id", "name", "accession", "header", "protein"),
    score_columns: Sequence[str] = ("score", "antigenicity", "prediction", "vaxijen", "value"),
) -> list[dict[str, object]]:
    """Score every record in *fasta_path*; returns rows in the ANTIGENICITY schema."""
    records = fasta_mod.read_fasta(fasta_path)

    if tool == "iapred":
        scores = _predict_iapred(
            fasta_path=fasta_path,
            raw_path=os.path.splitext(out_path)[0] + ".raw.txt",
            tool_dir=tool_dir,
        )
        backend = "iapred"

    elif tool == "kolaskar":
        scores = {record.id: antigenicity_kolaskar(record.sequence) for record in records}
        backend = "kolaskar-tongaonkar-1990"

    elif tool == "external":
        if not command:
            raise base.ToolError(
                "antigenicity.tool is 'external' but antigenicity.command is empty in "
                "config.yaml.\nSet it to the command line for your antigenicity "
                "predictor, using {input} and {output} placeholders, e.g.\n"
                '  command: "python tools/iapred/iapred.py -i {input} -o {output}"'
            )
        raw_path = os.path.splitext(out_path)[0] + ".raw.txt"
        formatted = command.format(input=os.path.abspath(fasta_path), output=os.path.abspath(raw_path))
        base.run_command(formatted.split())

        rows = base.read_delimited(raw_path)
        if not rows:
            raise base.ToolError(f"Antigenicity tool produced no rows in {raw_path}.")
        header = list(rows[0].keys())
        c_id = base.find_column(header, list(id_columns), "antigenicity sequence-ID", raw_path)
        c_score = base.find_column(header, list(score_columns), "antigenicity score", raw_path)
        scores = {
            str(row[c_id]).strip().lstrip(">").split()[0]: base.to_float(row[c_score], "antigenicity")
            for row in rows
        }
        backend = "external"

    else:
        raise base.ToolError(
            f"Unknown antigenicity backend {tool!r}. Use 'iapred', 'external' or 'kolaskar'."
        )

    out_rows: list[dict[str, object]] = []
    for record in records:
        if record.id not in scores:
            raise base.ToolError(
                f"The antigenicity backend returned no score for {record.id!r}.\n"
                f"Scored: {sorted(scores)[:10]}"
            )
        score = scores[record.id]
        # `None` means the backend was reached and declined to score this
        # sequence (see _IAPRED_SENTINELS). It fails the screen, but with its own
        # reason - reporting it as "below threshold" would invent a number.
        scored = score is not None
        out_rows.append(
            {
                "sequence_id": record.id,
                "sequence": record.sequence,
                "length": len(record.sequence),
                "antigenicity_score": score if scored else "",
                "antigenic": bool(scored and score >= threshold),
                "threshold": threshold,
                "tool": backend,
                "note": "" if scored else IAPRED_UNSCORABLE,
            }
        )
    return out_rows


# ---------------------------------------------------------------------------
# allergenicity
# ---------------------------------------------------------------------------

#: AlgPred 2.0's hybrid CSV, VERIFIED against the upstream source
#: (``algpred2.py``, ``hybrid()``, which writes the file)::
#:
#:     Subject,ML Score,MERCI Score,BLAST Score,Hybrid Score,Prediction
#:     Spike_Consensus_20.1,0.41,0.0,0.0,0.41,Non-Allergen
#:
#: The identifier column is ``Subject`` - not any of the names this adapter
#: originally guessed at. The guessed list would have raised on the real file.
_ALGPRED_ID_COLUMNS = ["subject", "seq_id", "sequence_id", "id", "name", "header", "query"]

#: ``Hybrid Score`` first, deliberately. With ``-m 2`` the file carries the
#: component scores *and* the combined one; the machine-learning column is not
#: the value the screen is configured against, and silently grading on it would
#: apply the hybrid threshold to a different model's scale. Same failure mode as
#: ``Antigenicity_Category`` in IApred above - the wrong column parses fine.
_ALGPRED_SCORE_COLUMNS = [
    "hybrid score", "hybrid_score", "combined score",
    "ml score", "ml_score", "score", "prediction score",
]

#: Values are ``Allergen`` / ``Non-Allergen``. Note AlgPred calls allergenic on
#: ``score > threshold`` strictly; this adapter prefers its verdict over
#: recomputing one, so that boundary stays the tool's to define.
_ALGPRED_CALL_COLUMNS = ["prediction", "call", "result", "type"]

# Running AlgPred - entry-point discovery, the rf_model/LFS check, the private
# working directory, the regenerated envfile and the source patch - lives in
# mpvpredpip.algpred, not here. It executes under workflow/envs/algpred.yaml
# (Python 3.8, scikit-learn 0.22), which this module's environment is not.
# Only the parsing below runs in the screening environment.


#: Component columns of AlgPred's hybrid score. The hybrid is their sum, so a
#: component that is zero on every row contributed nothing to any verdict.
_ALGPRED_COMPONENTS = {
    "MERCI Score": ["merci score", "merci_score"],
    "BLAST Score": ["blast score", "blast_score"],
}


def inert_hybrid_components(out_path: str) -> list[str]:
    """Component columns that were zero for every peptide AlgPred scored.

    With ``-m 2`` the hybrid score is RF + MERCI + BLAST. At 20-30 aa the
    latter two routinely contribute nothing - BLAST finds no allergen homologs
    and the IgE motifs do not hit - so the "hybrid" degenerates to the
    composition model while still being reported, and graded, as a hybrid.

    That is not an error, but it means the call rests on one line of evidence
    instead of three, so stage 6 says which components were inert rather than
    leaving the reader to infer it. An empty list means every component
    contributed somewhere.
    """
    rows = base.read_delimited(out_path)
    if not rows:
        return []
    header = list(rows[0].keys())

    inert: list[str] = []
    for label, candidates in _ALGPRED_COMPONENTS.items():
        column = base.find_column(header, candidates, f"AlgPred {label}", out_path,
                                  required=False)
        if not column:
            continue
        values = [base.to_float(row.get(column, 0) or 0, label) for row in rows]
        if values and all(value == 0 for value in values):
            inert.append(label)
    return inert


def parse_allergenicity(
    out_path: str,
    fasta_path: str,
    threshold: float = 0.3,
) -> dict[str, tuple[float, bool]]:
    """Read AlgPred 2.0's CSV; returns ``{sequence_id: (score, is_allergen)}``.

    Parsing only. AlgPred itself runs in its own Snakemake rule under
    ``workflow/envs/algpred.yaml``, because its ``rf_model`` was pickled with
    scikit-learn 0.22 and will not unpickle under the 1.5.2 that IApred
    requires. :mod:`mpvpredpip.algpred` holds the invocation and the reasoning;
    the split keeps that ancient environment free of everything but the tool,
    and keeps output handling here with the rest of stage 6.

    The coverage check at the end is the load-bearing part. AlgPred's ``-d``
    default (1) writes out *only* the peptides it calls allergenic, so a
    missing peptide would otherwise read as one that passed. A peptide absent
    from the file has an **unknown** verdict, not a clean one - conflating the
    two is exactly what made the MHC-I failure invisible in stage 4.
    """
    rows = base.read_delimited(out_path)
    if not rows:
        raise base.ToolError(
            f"AlgPred produced no rows in {out_path}.\n"
            "With -d 2 it should write one row per input sequence even when "
            "none are allergenic, so an empty file means it did not run."
        )
    header = list(rows[0].keys())
    c_id = base.find_column(header, _ALGPRED_ID_COLUMNS, "AlgPred sequence-ID", out_path)
    c_score = base.find_column(header, _ALGPRED_SCORE_COLUMNS, "AlgPred score", out_path)
    c_call = base.find_column(header, _ALGPRED_CALL_COLUMNS, "AlgPred call", out_path,
                              required=False)

    results: dict[str, tuple[float, bool]] = {}
    for row in rows:
        sequence_id = str(row[c_id]).strip().lstrip(">").split()[0]
        score = base.to_float(row[c_score], "AlgPred score")
        if c_call and row.get(c_call):
            # "Allergen" vs "Non-Allergen": test for the negation first, since
            # "allergen" is a substring of both.
            call = str(row[c_call]).lower()
            is_allergen = "non" not in call and "allergen" in call
        else:
            is_allergen = score >= threshold
        results[sequence_id] = (score, is_allergen)

    # Coverage check. A peptide AlgPred did not report is not a peptide that
    # passed - it is one whose verdict is unknown, and the two must never be
    # conflated. This is the guard against -d reverting to its default.
    submitted = [record.id for record in fasta_mod.read_fasta(fasta_path)]
    missing = [identifier for identifier in submitted if identifier not in results]
    if missing:
        raise base.ToolError(
            f"AlgPred scored {len(results)} of {len(submitted)} submitted peptides; "
            f"{len(missing)} are absent from {out_path}.\n"
            f"  Missing: {missing[:5]}\n"
            f"  Header was: {header}\n\n"
            "AlgPred's -d default (1) reports only allergenic peptides. The run "
            "passes -d 2; if a release has renamed that flag, correct it in "
            "mpvpredpip.algpred rather than treating the gap as a pass."
        )
    return results


# ---------------------------------------------------------------------------
# autoimmunity
# ---------------------------------------------------------------------------
# Not here. The human-proteome screen is done in-process by mpvpredpip.autoimmunity
# rather than by the PIR Peptide Match standalone - that module's docstring
# records why, and carries the L/I-equivalence caveat that goes with the choice.
#
# Nothing in this file drives a Java tool any more, which is why openjdk has
# left workflow/envs/screening.yaml.
