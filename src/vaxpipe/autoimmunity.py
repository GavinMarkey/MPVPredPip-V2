"""Stage 6 autoimmunity screen - exact matching against the human proteome.

A candidate peptide that shares a run of residues with a human protein risks
provoking a self-reactive response, so stage 6 rejects any consensus peptide
containing an exact window of ``min_exact_match_length`` residues that also
occurs in the human reference proteome.

------------------------------------------------------------------------------
Why this is not the PIR Peptide Match standalone
------------------------------------------------------------------------------
The pipeline specification names PIR Peptide Match, and stage 6 was originally
written to drive it: a rule built a Lucene index over the proteome, and the
screen queried that index through a ``.jar``.

That tool is free (MIT, github.com/udel-cbcb/PeptideMatch) but awkward to obtain
in the form this pipeline needs. It publishes no releases and ships no built
jar; ``peptidematch_cmd`` is an Eclipse project with no ``pom.xml``, so using it
means compiling Java against Lucene by hand. The maintained path has moved to an
Elasticsearch service needing Docker, Java 17, Maven and a multi-gigabyte heap.

Set against what the screen actually computes - does any k-mer of a ~20-30 aa
peptide occur verbatim in a 20,000-protein FASTA - that is a great deal of
machinery. A dozen candidates yield a few hundred windows; the proteome is about
11 MB of residues. :func:`find_matches` answers it in roughly a second.

So the check is reimplemented here, exactly as the IEDB cluster tool is
reimplemented in :mod:`vaxpipe.cluster` and stage 7's conservancy is
reimplemented from the IEDB definition - and for the same reason: the dependency
cost of the original outweighs what it contributes. The index rule and the
``java`` environment are gone with it.

**Matching is exact.** PIR Peptide Match offers an L/I equivalence mode that
treats isobaric leucine and isoleucine as the same residue; it is off by default
there and is not implemented here. Enabling it would make the screen stricter
(more peptides flagged), so its absence is not a silent relaxation - but it does
mean a candidate differing from a human window only at an L/I position is not
flagged.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass

from . import fasta as fasta_mod

#: Joins proteins into one searchable blob. Not a residue, so no window can
#: straddle two proteins and produce a match that exists in neither.
SEPARATOR = "\n"

#: The 20 standard residues. A window containing anything else - most often the
#: ``X`` that :mod:`vaxpipe.cluster` writes where a cluster's members disagree -
#: is not a sequence that can be searched for literally: ``X`` in a candidate
#: means "unresolved", while ``X`` in the proteome means "unsequenced", and
#: matching one against the other asserts a correspondence neither file claims.
#: Such windows are skipped and counted, never silently treated as clean.
STANDARD_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class Match:
    """One exact occurrence of *window* inside a human protein."""

    window: str
    protein_id: str
    offset: int  # 0-based, within that protein

    def describe(self) -> str:
        return f"{self.window}@{self.protein_id}:{self.offset + 1}"


def windows(sequence: str, length: int) -> list[str]:
    """All overlapping windows of *length* residues in *sequence*.

    A candidate is a concern if *any part* of it occurs in the human proteome,
    so the comparison slides across the peptide rather than testing it whole. A
    peptide shorter than *length* is returned unchanged: it cannot contain a
    window of that size, and dropping it would exempt the shortest candidates
    from the screen entirely.
    """
    sequence = sequence.upper()
    if length <= 0 or len(sequence) < length:
        return [sequence] if sequence else []
    return [sequence[index : index + length] for index in range(len(sequence) - length + 1)]


class HumanProteome:
    """The proteome held as one blob, with offsets mapped back to protein IDs.

    Concatenating is what makes this fast enough to do in-process: one
    :meth:`str.find` per query window scans all 11 MB at C speed, instead of
    looping over 20,000 separate sequences in Python for every window.
    """

    def __init__(self, path: str) -> None:
        records = fasta_mod.read_fasta(path)
        if not records:
            raise ValueError(
                f"No sequences in the human proteome FASTA at {path!r}.\n"
                "Fetch it with:  bash setup/fetch_human_proteome.sh"
            )

        self.path = path
        self.n_proteins = len(records)
        self._ids: list[str] = []
        self._starts: list[int] = []

        pieces: list[str] = []
        cursor = 0
        for record in records:
            sequence = record.sequence.upper()
            self._ids.append(record.id)
            self._starts.append(cursor)
            pieces.append(sequence)
            cursor += len(sequence) + len(SEPARATOR)
        self._blob = SEPARATOR.join(pieces)
        self.n_residues = sum(len(piece) for piece in pieces)

    def find(self, window: str) -> list[Match]:
        """Every exact occurrence of *window*, as :class:`Match` records."""
        if not window:
            return []
        found: list[Match] = []
        cursor = self._blob.find(window)
        while cursor != -1:
            # bisect_right - 1 gives the protein whose start is the greatest one
            # at or before this offset, i.e. the protein the hit falls inside.
            index = bisect.bisect_right(self._starts, cursor) - 1
            found.append(
                Match(window=window, protein_id=self._ids[index], offset=cursor - self._starts[index])
            )
            cursor = self._blob.find(window, cursor + 1)
        return found


@dataclass
class ScreenResult:
    """Outcome of screening one peptide."""

    matches: list[Match]
    skipped_windows: list[str]

    @property
    def flagged(self) -> bool:
        return bool(self.matches)

    def describe(self) -> str:
        return ";".join(sorted({match.describe() for match in self.matches}))


def find_matches(
    peptides: dict[str, str],
    proteome: HumanProteome,
    length: int,
) -> dict[str, ScreenResult]:
    """Screen ``{peptide_id: sequence}`` against *proteome*.

    Returns one :class:`ScreenResult` per input peptide, including those with no
    matches - the caller needs a verdict for every candidate, and an absent key
    would be indistinguishable from a peptide that was never screened.

    Windows are deduplicated across peptides before searching, since overlapping
    candidates from the same cluster frequently share them.
    """
    results = {key: ScreenResult(matches=[], skipped_windows=[]) for key in peptides}

    owners: dict[str, set[str]] = {}
    for key, sequence in peptides.items():
        for window in windows(sequence, length):
            if not set(window) <= STANDARD_RESIDUES:
                results[key].skipped_windows.append(window)
                continue
            owners.setdefault(window, set()).add(key)

    for window, keys in owners.items():
        hits = proteome.find(window)
        if not hits:
            continue
        for key in keys:
            results[key].matches.extend(hits)

    return results
