"""Minimal FASTA reading and writing plus NCBI defline parsing.

Deliberately dependency-free so that every conda environment in the workflow can
import it, including the ones pinned to old interpreters for legacy predictors.
"""

from __future__ import annotations

import hashlib
import os
import re
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass

#: ``>ACCESSION some protein description [organism or strain]``
_DEFLINE = re.compile(r"^(?P<accession>\S+)\s*(?P<name>.*?)\s*(?:\[(?P<strain>[^\]]+)\])?\s*$")

#: ``ACCESSION protein name - strain name`` - the filename convention for the
#: files in "0 - Input". Only a FALLBACK, used when the header lacks a protein
#: name or a [strain]; note the whitespace required on both sides of the dash.
#: The notice users see (vaxpipe.inputs) is tested against this pattern.
_FILENAME = re.compile(r"^(?P<accession>\S+)\s+(?P<name>.+?)\s+-\s+(?P<strain>.+)$")

_VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")


@dataclass
class Record:
    """A single FASTA record."""

    header: str
    sequence: str
    accession: str = ""
    name: str = ""
    strain: str = ""
    source_path: str = ""
    protein_label: str = ""

    @property
    def id(self) -> str:
        return self.header.split()[0] if self.header else ""

    @property
    def md5(self) -> str:
        return hashlib.md5(self.sequence.upper().encode()).hexdigest()

    def __len__(self) -> int:
        return len(self.sequence)


def read_fasta(path: str) -> list[Record]:
    """Read every record in *path*.

    Sequence lines are concatenated, whitespace stripped and upper-cased. ``*``
    and ``-`` are removed so that downstream coordinates always refer to the
    ungapped protein.
    """
    records: list[Record] = []
    header: str | None = None
    chunks: list[str] = []

    with open(path, encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(_build(header, chunks, path))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)

    if header is not None:
        records.append(_build(header, chunks, path))
    return records


def _build(header: str, chunks: Sequence[str], path: str) -> Record:
    sequence = "".join(chunks).upper().replace("-", "").replace("*", "")
    record = Record(header=header, sequence=sequence, source_path=path)
    _annotate(record, path)
    return record


def _annotate(record: Record, path: str) -> None:
    """Fill in accession/name/strain from the defline, falling back to the filename."""
    match = _DEFLINE.match(record.header)
    if match:
        record.accession = match.group("accession") or ""
        record.name = (match.group("name") or "").strip()
        record.strain = (match.group("strain") or "").strip()

    if not record.strain or not record.name:
        stem = os.path.splitext(os.path.basename(path))[0]
        fallback = _FILENAME.match(stem)
        if fallback:
            record.accession = record.accession or fallback.group("accession")
            record.name = record.name or fallback.group("name").strip()
            record.strain = record.strain or fallback.group("strain").strip()

    if not record.name:
        record.name = "unknown protein"
    if not record.strain:
        record.strain = "unknown strain"


def write_fasta(path: str, records: Iterable[tuple[str, str]], width: int = 60) -> int:
    """Write ``(header, sequence)`` pairs to *path*. Returns the number written."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    count = 0
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            if width and width > 0:
                for i in range(0, len(sequence), width):
                    handle.write(sequence[i : i + width] + "\n")
            else:
                handle.write(sequence + "\n")
            count += 1
    return count


def iter_fasta(path: str) -> Iterator[Record]:
    yield from read_fasta(path)


def validate(record: Record, min_length: int = 1) -> list[str]:
    """Return a list of human-readable problems with *record*; empty means valid."""
    problems: list[str] = []
    if len(record.sequence) < min_length:
        problems.append(
            f"sequence is {len(record.sequence)} residues, below the minimum of {min_length}"
        )
    if not record.sequence:
        problems.append("sequence is empty")
        return problems

    bad = sorted(set(record.sequence) - _VALID_AA)
    if bad:
        problems.append(f"contains non-amino-acid characters: {''.join(bad)}")
    if set(record.sequence) <= set("ACGTUN"):
        problems.append(
            "looks like a nucleotide sequence - this pipeline requires protein FASTA"
        )
    return problems


def normalise_protein_label(name: str, aliases: dict[str, str]) -> str:
    """Map a free-text protein name onto a short label used in epitope IDs.

    Matching is case-insensitive and falls back to the longest alias that occurs
    as a substring of *name*, so "spike glycoprotein precursor" still resolves to
    "Spike". Unmatched names are CamelCased so the pipeline never silently drops
    a protein.
    """
    lowered = name.lower().strip()
    if lowered in aliases:
        return aliases[lowered]

    best: tuple[int, str] | None = None
    for alias, label in aliases.items():
        if alias.lower() in lowered:
            if best is None or len(alias) > best[0]:
                best = (len(alias), label)
    if best:
        return best[1]

    cleaned = re.sub(r"[^A-Za-z0-9 ]+", " ", name).split()
    return "".join(word.capitalize() for word in cleaned) or "Unknown"
