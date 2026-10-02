"""Shared machinery for external-tool adapters.

Every third-party predictor in this pipeline is wrapped by an adapter with two
halves:

``run(...)``    invoke the tool and leave its *untouched* native output on disk
``parse(...)``  convert that native output into the canonical schema

Native output is always kept. When a predictor is upgraded and changes its
columns, only ``parse`` needs editing, and the raw file is still there to check
against. Adapters never delete or rewrite what the tool produced.

Parsers here locate columns **by name**, not by position, and fail loudly with
the actual header they saw. Third-party immunoinformatics tools reorder and
rename columns between releases; a positional parser silently mis-assigns scores
when that happens, which is the worst possible failure mode for this pipeline.
"""

from __future__ import annotations

import os
import shlex
import subprocess
import sys
from collections.abc import Iterable, Mapping, Sequence


class ToolError(RuntimeError):
    """Raised when an external tool is missing, fails, or emits an unusable format."""


class ToolMissing(ToolError):
    """Raised when the tool is not installed. Carries install guidance."""


def run_command(
    command: Sequence[str],
    stdout_path: str | None = None,
    cwd: str | None = None,
    env: Mapping[str, str] | None = None,
    check: bool = True,
    timeout: int | None = None,
) -> subprocess.CompletedProcess:
    """Run *command*, optionally teeing stdout to a file.

    Raises :class:`ToolError` with the tool's own stderr on failure, because a
    truncated traceback from Snakemake is useless when the real message is three
    frames down inside a predictor's own script.
    """
    printable = " ".join(shlex.quote(part) for part in command)
    print(f"[vaxpipe] $ {printable}", file=sys.stderr, flush=True)

    merged = dict(os.environ)
    if env:
        merged.update(env)

    if stdout_path:
        os.makedirs(os.path.dirname(os.path.abspath(stdout_path)) or ".", exist_ok=True)

    try:
        with open(stdout_path, "w", encoding="utf-8") if stdout_path else _Null() as sink:
            result = subprocess.run(
                list(command),
                stdout=sink if stdout_path else subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=cwd,
                env=merged,
                text=True,
                timeout=timeout,
            )
    except FileNotFoundError as exc:
        raise ToolMissing(
            f"Executable not found while running:\n  {printable}\n\n"
            f"{exc}\n\nRun 'python setup/check_tools.py' and see tools/README.md."
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise ToolError(f"Timed out after {timeout}s:\n  {printable}") from exc

    if check and result.returncode != 0:
        raise ToolError(
            f"Command failed with exit code {result.returncode}:\n  {printable}\n\n"
            f"--- stderr ---\n{(result.stderr or '').strip()[:4000]}"
        )
    return result


class _Null:
    """Context manager standing in for a file handle when stdout is captured."""

    def __enter__(self):
        return None

    def __exit__(self, *exc):
        return False


def require_path(path: str, what: str, hint: str = "") -> str:
    """Assert that a tool directory or file exists, with actionable guidance."""
    if os.path.exists(path):
        return path
    message = f"{what} not found at: {path}"
    if hint:
        message += f"\n\n{hint}"
    message += "\n\nRun 'python setup/check_tools.py' and see tools/README.md."
    raise ToolMissing(message)


# ---------------------------------------------------------------------------
# tolerant delimited-file reading
# ---------------------------------------------------------------------------

def sniff_delimiter(path: str) -> str:
    with open(path, encoding="utf-8-sig", errors="replace") as handle:
        sample = handle.read(8192)
    counts = {"\t": sample.count("\t"), ",": sample.count(","), ";": sample.count(";")}
    return max(counts, key=counts.get) if max(counts.values()) else "\t"


def read_delimited(
    path: str,
    delimiter: str | None = None,
    header_contains: str = "",
) -> list[dict[str, str]]:
    """Read a delimited file into dicts, auto-detecting the delimiter.

    *header_contains* names a column that the header row must contain. When
    given, leading lines are skipped until one whose fields include it, compared
    case-insensitively. Without it the first line is the header, which is wrong
    for any tool that prints settings before its table - the IEDB immunogenicity
    predictor emits::

        masking: default
        masked variables: [1, 2, 9]

        peptide,length,score
        AAAAAAAAA,9,0.12345

    and its first line was being read as a one-column header called
    ``masking: default``.
    """
    import csv

    if not os.path.exists(path):
        raise ToolError(f"Expected tool output does not exist: {path}")
    delimiter = delimiter or sniff_delimiter(path)

    with open(path, encoding="utf-8-sig", errors="replace", newline="") as handle:
        lines = handle.readlines()

    start = 0
    if header_contains:
        wanted = header_contains.strip().lower()
        for index, line in enumerate(lines):
            fields = [f.strip().lower() for f in line.rstrip("\r\n").split(delimiter)]
            if wanted in fields:
                start = index
                break
        else:
            # Do not fall back to line 0. A tool that printed a message instead
            # of a table would otherwise parse as an empty result, which reads
            # downstream as "nothing qualified" rather than "this never ran".
            preview = "".join(f"    {line.rstrip()}\n" for line in lines[:5]) or "    (empty)\n"
            raise ToolError(
                f"No header row containing {header_contains!r} was found in {path}.\n"
                f"  First lines were:\n{preview}\n"
                "Either the tool's output format has changed, or it failed and printed "
                "a message where the table should be."
            )

    reader = csv.DictReader(lines[start:], delimiter=delimiter)
    return [row for row in reader if _has_content(row)]


def _has_content(row: dict) -> bool:
    """True when a parsed row holds anything but whitespace.

    Values are not always strings. When a line has more fields than the header -
    which happens with ragged tool output, and with any file whose header this
    reader has mis-identified - csv.DictReader collects the surplus into a
    **list** under a ``None`` key. Calling .strip() on that raises
    AttributeError, turning a recoverable parse into a crash with a traceback
    that points nowhere near the real problem.
    """
    for value in row.values():
        if isinstance(value, (list, tuple)):
            if any((item or "").strip() for item in value):
                return True
        elif (value or "").strip():
            return True
    return False


def find_column(
    header: Iterable[str],
    candidates: Sequence[str],
    what: str,
    path: str = "",
    required: bool = True,
) -> str | None:
    """Resolve one of *candidates* against *header*, ignoring case and separators.

    Matching is exact-normalised first, then prefix, then substring. Returns the
    real header string so the caller can index rows with it.
    """
    columns = [column for column in header if column is not None]
    norm = {_norm(column): column for column in columns}

    for candidate in candidates:
        key = _norm(candidate)
        if key in norm:
            return norm[key]
    for candidate in candidates:
        key = _norm(candidate)
        for normalised, original in norm.items():
            if normalised.startswith(key):
                return original
    for candidate in candidates:
        key = _norm(candidate)
        for normalised, original in norm.items():
            if key in normalised:
                return original

    if not required:
        return None
    raise ToolError(
        f"Could not find the {what} column"
        + (f" in {path}" if path else "")
        + f".\n  Tried: {list(candidates)}\n  Header was: {columns}\n\n"
        "The predictor's output format has probably changed. Update the matching "
        "adapter in src/vaxpipe/adapters/ - the raw output has been preserved."
    )


def _norm(text: str) -> str:
    return "".join(character for character in str(text).lower() if character.isalnum())


def to_float(value: str, what: str = "value") -> float:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ToolError(f"Could not read {what} as a number: {value!r}") from exc
