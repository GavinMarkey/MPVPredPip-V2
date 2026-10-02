"""Boilerplate shared by every Snakemake script in this workflow.

Snakemake runs scripts with the working directory set to the project root, so
``src`` is located relative to that. A fallback walks up from this file in case
the workflow is invoked with an explicit ``--directory``.
"""

from __future__ import annotations

import os
import sys
from contextlib import contextmanager


def bootstrap() -> str:
    """Put ``src`` on ``sys.path`` and return the project root."""
    candidates = [os.getcwd()]

    here = os.path.dirname(os.path.abspath(__file__))
    # workflow/scripts -> workflow -> project root
    candidates.append(os.path.dirname(os.path.dirname(here)))

    for root in candidates:
        src = os.path.join(root, "src")
        if os.path.isdir(os.path.join(src, "vaxpipe")):
            if src not in sys.path:
                sys.path.insert(0, src)
            return root

    raise SystemExit(
        "Could not locate the 'src/vaxpipe' package.\n"
        f"Looked under: {candidates}\n"
        "Run snakemake from the project root, or pass --directory <project root>."
    )


@contextmanager
def tee_log(log_path: str | None):
    """Mirror stdout and stderr into *log_path* while still showing them live."""
    if not log_path:
        yield
        return

    os.makedirs(os.path.dirname(os.path.abspath(log_path)) or ".", exist_ok=True)
    handle = open(log_path, "w", encoding="utf-8")

    class _Tee:
        def __init__(self, *streams):
            self._streams = streams

        def write(self, text):
            for stream in self._streams:
                stream.write(text)
                stream.flush()

        def flush(self):
            for stream in self._streams:
                stream.flush()

    original_out, original_err = sys.stdout, sys.stderr
    sys.stdout = _Tee(original_out, handle)
    sys.stderr = _Tee(original_err, handle)
    try:
        yield
    finally:
        sys.stdout, sys.stderr = original_out, original_err
        handle.close()


def log_path(smk) -> str | None:
    try:
        return str(smk.log[0]) if smk.log else None
    except (IndexError, AttributeError):
        return None


def write_summary(path: str, title: str, lines: list[str]) -> None:
    """Write a short per-stage markdown summary next to the stage's tables."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"# {title}\n\n")
        for line in lines:
            handle.write(line.rstrip() + "\n")
