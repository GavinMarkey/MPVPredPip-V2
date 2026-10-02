"""Tab-separated table I/O against the schemas in :mod:`vaxpipe.schemas`.

Dependency-free on purpose: the legacy predictor environments cannot all carry
pandas, and every stage boundary in this pipeline is a plain TSV so that any
intermediate file can be opened in a text editor or Excel without tooling.
"""

from __future__ import annotations

import csv
import os
from collections.abc import Iterable, Sequence
from typing import Any

#: Written into empty cells so that a TSV never has ambiguous blank columns.
NA = ""


def write_table(path: str, columns: Sequence[str], rows: Iterable[dict[str, Any]]) -> int:
    """Write *rows* to *path* with exactly *columns*, in order.

    Keys absent from a row become empty strings; keys not in *columns* are
    dropped. Returns the number of data rows written.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    count = 0
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([_fmt(row.get(column, NA)) for column in columns])
            count += 1
    return count


def read_table(path: str) -> list[dict[str, str]]:
    """Read a TSV written by :func:`write_table` into a list of dicts."""
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _fmt(value: Any) -> str:
    if value is None:
        return NA
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, float):
        # Keep tables diffable: fixed precision, no scientific notation, and no
        # trailing zeros that would churn between runs.
        text = f"{value:.6f}".rstrip("0").rstrip(".")
        return text or "0"
    return str(value)


def as_bool(value: Any) -> bool:
    """Parse the strings this module writes (and common variants) back to bool."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "t", "yes", "y", "1"}


def as_float(value: Any, default: float = float("nan")) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def as_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default
