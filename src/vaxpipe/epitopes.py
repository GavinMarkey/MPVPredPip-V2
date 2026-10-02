"""Epitope construction, identifier assignment and deduplication.

Two jobs live here:

1. Turning the *per-residue* scores that BepiPred and EpiDope emit into discrete
   epitope intervals (:func:`residues_to_epitopes`).
2. Assigning the stable identifiers that stage 5 relies on to decide whether a
   cluster contains the required mix of epitope types.

Identifier format::

    {ProteinLabel}_{B_cell|T_cell}_{BepiPred|EpiDope|MHCI|MHCII}_{n}.{m}

``n`` indexes distinct peptide sequences within a (protein, tool) pair; ``m``
indexes repeated occurrences of that same peptide (a different start position, or
the same peptide found in another strain). Keeping the ``_B_cell`` / ``_T_cell``
infix means the identifiers stay readable in the IEDB cluster output, exactly as
in the earlier SARS-CoV-2 analysis, while the added tool token is what lets the
"one epitope from BepiPred *and* one from EpiDope" rule be enforced mechanically
rather than by eye.
"""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from . import schemas

EPITOPE_ID = re.compile(
    r"^(?P<protein>.+?)_(?P<cls>B_cell|T_cell)_(?P<tool>BepiPred|EpiDope|MHCI|MHCII)_"
    r"(?P<n>\d+)\.(?P<m>\d+)$"
)

#: Maps the tool token inside an identifier back to the canonical source_tool.
_TOOL_TOKEN = {
    "BepiPred": "bepipred",
    "EpiDope": "epidope",
    "MHCI": "mhc_i",
    "MHCII": "mhc_ii",
}


@dataclass
class Epitope:
    protein_label: str
    strain: str
    accession: str
    source_tool: str
    start: int           # 1-based inclusive
    end: int             # 1-based inclusive
    sequence: str
    score: float
    allele: str = ""
    percentile: str = ""
    epitope_id: str = ""

    @property
    def epitope_class(self) -> str:
        return schemas.CLASS_TAGS[self.source_tool][2]

    @property
    def length(self) -> int:
        return len(self.sequence)

    def as_row(self) -> dict[str, object]:
        return {
            "epitope_id": self.epitope_id,
            "protein_label": self.protein_label,
            "strain": self.strain,
            "accession": self.accession,
            "source_tool": self.source_tool,
            "epitope_class": self.epitope_class,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "sequence": self.sequence,
            "score": self.score,
            "allele": self.allele,
            "percentile": self.percentile,
        }


# ---------------------------------------------------------------------------
# per-residue scores -> epitope intervals
# ---------------------------------------------------------------------------

def residues_to_epitopes(
    sequence: str,
    scores: Sequence[float],
    threshold: float,
    length: int = 14,
    step: int = 2,
    merge_gap: int = 1,
    min_region: int = 6,
) -> list[tuple[int, int, float]]:
    """Convert per-residue scores into ``(start, end, mean_score)`` epitopes.

    Residues scoring at or above *threshold* form above-threshold runs. Runs of
    sub-threshold residues no longer than *merge_gap* that sit *between* two such
    runs are absorbed, which stops a single dip from splitting one real epitope
    in two.

    Each run is then **tiled** with overlapping peptides of exactly *length*
    residues, advancing *step* residues at a time. Linear B-cell epitopes are
    typically 12-16 residues, so a long above-threshold region is not one long
    epitope - it is a stretch that contains several, and which offsets are the
    real ones is not knowable from a per-residue score. Emitting every offset
    hands that decision to stage 5, where clustering collapses the overlaps, and
    to the downstream screens.

    A 39-residue run at ``length=14, step=2`` yields 13 peptides::

        PHDWTKNITDKIDQ   (1-14)
        DWTKNITDKIDQII   (3-16)
        TKNITDKIDQIIHD   (5-18)
        ...

    Where the stride does not land exactly on the end of a run, a final window
    anchored at that end is added, so the tail of a run is never lost.

    A run **shorter** than *length* is padded out to *length* with the residues
    flanking it, centred on the run, and clamped at the ends of the sequence. The
    flanks are real sequence that merely scores below the threshold, and a
    10-residue predicted core is still evidence - discarding it would mean
    discarding every epitope for a protein whose regions all fall short, and with
    it the protein itself at stage 5. Runs shorter than *min_region* are dropped
    as too small to be credible, and nothing is emitted at all if the whole
    sequence is shorter than *length*.

    Coordinates returned are 1-based and inclusive.
    """
    if len(sequence) != len(scores):
        raise ValueError(
            f"sequence length {len(sequence)} does not match score count {len(scores)}"
        )
    if length <= 0:
        raise ValueError(f"length must be positive, got {length}")
    if step <= 0:
        raise ValueError(f"step must be positive, got {step}")

    mask = [score >= threshold for score in scores]
    mask = _absorb_gaps(mask, merge_gap)

    out: list[tuple[int, int, float]] = []
    for start, end in _runs(mask):
        run_length = end - start + 1

        if run_length < length:
            if run_length < min_region or len(sequence) < length:
                continue
            windows = [_pad(start, end, length, len(sequence))]
        else:
            windows = _tile(start, end, length, step)

        for lo, hi in windows:
            window = scores[lo : hi + 1]
            mean = sum(window) / len(window) if window else 0.0
            out.append((lo + 1, hi + 1, mean))
    return out


def _pad(start: int, end: int, length: int, sequence_length: int) -> tuple[int, int]:
    """Centre a window of *length* on the run ``[start, end]``.

    Clamped to the sequence, so a run near either terminus produces a window that
    still has the full length but is no longer centred on it.
    """
    short_by = length - (end - start + 1)
    lo = start - short_by // 2
    lo = max(0, min(lo, sequence_length - length))
    return (lo, lo + length - 1)


def _absorb_gaps(mask: list[bool], merge_gap: int) -> list[bool]:
    """Fill interior False runs of length <= *merge_gap*."""
    if merge_gap <= 0:
        return list(mask)
    mask = list(mask)
    n = len(mask)
    i = 0
    while i < n:
        if mask[i]:
            i += 1
            continue
        j = i
        while j < n and not mask[j]:
            j += 1
        # Only interior gaps qualify: a leading or trailing run is not "between"
        # two epitopes and must not extend the prediction past its evidence.
        if i > 0 and j < n and (j - i) <= merge_gap:
            for k in range(i, j):
                mask[k] = True
        i = j
    return mask


def _runs(mask: Sequence[bool]) -> list[tuple[int, int]]:
    """Maximal runs of True as 0-based inclusive ``(start, end)`` pairs."""
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for i, value in enumerate(mask):
        if value and start is None:
            start = i
        elif not value and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(mask) - 1))
    return runs


def _tile(start: int, end: int, length: int, step: int) -> list[tuple[int, int]]:
    """Tile ``[start, end]`` with windows of *length*, advancing by *step*.

    Every window is exactly *length* residues. A run shorter than that yields
    nothing. If the stride does not land on *end*, one further window ending at
    *end* is appended, so no tail is dropped - it may overlap its predecessor by
    more than the usual step.
    """
    run_length = end - start + 1
    if run_length < length:
        return []

    out = [(lo, lo + length - 1) for lo in range(start, end - length + 2, step)]
    last_end = out[-1][1]
    if last_end < end:
        out.append((end - length + 1, end))
    return out


# ---------------------------------------------------------------------------
# identifiers
# ---------------------------------------------------------------------------

def assign_ids(epitopes: Iterable[Epitope]) -> list[Epitope]:
    """Assign ``{Protein}_{class}_{Tool}_{n}.{m}`` identifiers in place.

    Ordering is deterministic - by protein, tool, sequence, strain then start -
    so that re-running the pipeline on unchanged inputs reproduces identical
    identifiers and the outputs stay diffable.
    """
    items = sorted(
        epitopes,
        key=lambda e: (e.protein_label, e.source_tool, e.sequence, e.strain, e.start),
    )
    counters: dict[tuple[str, str], dict[str, int]] = defaultdict(dict)
    occurrences: dict[tuple[str, str, str], int] = defaultdict(int)

    for epitope in items:
        key = (epitope.protein_label, epitope.source_tool)
        index_map = counters[key]
        if epitope.sequence not in index_map:
            index_map[epitope.sequence] = len(index_map) + 1
        n = index_map[epitope.sequence]

        occ_key = (epitope.protein_label, epitope.source_tool, epitope.sequence)
        occurrences[occ_key] += 1
        m = occurrences[occ_key]

        cls_tag, tool_tag, _ = schemas.CLASS_TAGS[epitope.source_tool]
        epitope.epitope_id = f"{epitope.protein_label}_{cls_tag}_{tool_tag}_{n}.{m}"
    return items


def parse_id(epitope_id: str) -> tuple[str, str, str] | None:
    """``(protein_label, source_tool, epitope_class)`` or ``None`` if unparseable."""
    match = EPITOPE_ID.match(epitope_id.strip())
    if not match:
        return None
    tool = _TOOL_TOKEN[match.group("tool")]
    return match.group("protein"), tool, schemas.CLASS_TAGS[tool][2]


# ---------------------------------------------------------------------------
# deduplication
# ---------------------------------------------------------------------------

def deduplicate(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    """Collapse identical peptides from the same tool into one clustering input.

    The earlier pipeline did this with a script that concatenated headers with
    hyphens; here the collapsed strains and alleles are kept in dedicated columns
    so the provenance of a peptide survives into the final tables instead of
    being buried in a defline.

    Peptides are keyed on (protein, tool, sequence): the *same* peptide predicted
    by two different tools stays as two rows, because stage 5 must be able to see
    that both BepiPred and EpiDope support it.
    """
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        key = (
            str(row["protein_label"]),
            str(row["source_tool"]),
            str(row["sequence"]),
        )
        grouped[key].append(row)

    out: list[dict[str, object]] = []
    for (protein, tool, sequence), members in sorted(grouped.items()):
        # The representative is the first identifier in sorted order, which is
        # the ".1" occurrence assigned by assign_ids.
        members = sorted(members, key=lambda r: str(r["epitope_id"]))
        scores = [float(r["score"]) for r in members if str(r.get("score", "")).strip()]
        strains = sorted({str(r["strain"]) for r in members if str(r.get("strain", "")).strip()})
        alleles = sorted({str(r["allele"]) for r in members if str(r.get("allele", "")).strip()})

        out.append(
            {
                "epitope_id": members[0]["epitope_id"],
                "protein_label": protein,
                "source_tool": tool,
                "epitope_class": schemas.CLASS_TAGS[tool][2],
                "sequence": sequence,
                "length": len(sequence),
                "best_score": max(scores) if scores else "",
                "strains": "|".join(strains),
                "alleles": "|".join(alleles),
                "n_occurrences": len(members),
            }
        )
    return out
