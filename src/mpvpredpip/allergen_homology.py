"""Allergen homology screening - sequence evidence, not composition.

AlgPred's hybrid score is the pipeline's primary allergenicity call, but at the
20-30 aa this pipeline produces it degenerates: BLAST finds no allergen homologs
and the IgE motifs do not hit, so the "hybrid" is its amino-acid-composition
model alone and nothing else. See ``inert_hybrid_components`` in
:mod:`mpvpredpip.adapters.screening`, and the stage 6 caveats it feeds.

This module adds the missing kind of evidence - does the candidate actually
*look like* a known allergen - by two separate criteria against a FASTA of known
allergens. Two, because one number cannot serve both peptide lengths this
pipeline handles.

------------------------------------------------------------------------------
Criterion 1: contiguous exact match, on the constituent epitopes
------------------------------------------------------------------------------
A run of identical residues shared with a known allergen, searched for
literally. This is the criterion the FAO/WHO and EFSA allergenicity guidelines
use for short sequences, where a percentage is meaningless (below), and it is
length-honest: 8 identical residues in a row means the same thing whether the
query is 8 residues or 800.

Applied to the *constituent epitopes* rather than the consensus, because the
consensus is a synthesised sequence. :mod:`mpvpredpip.cluster` builds it from a
cluster's aligned members, so it need not occur in any real protein - which is
exactly why stage 7's conservancy comes out at 0% for candidates whose members
are conserved. Searching the consensus for literal matches therefore asks
whether a sequence that may exist nowhere shares a window with an allergen. The
constituent epitopes are real predicted sequences, so they are what the exact
criterion can speak about.

------------------------------------------------------------------------------
Criterion 2: percent identity, on the consensus sequence
------------------------------------------------------------------------------
Best ungapped identity over the whole allergen set, against a configurable
threshold. ``30`` is the usual figure quoted for protein homology, and it is
sound for full-length proteins and domains - but that relationship is steeply
length-dependent, and at 20-30 aa the threshold is at chance.

Measured on this pipeline's own output, against the 2426-sequence allergen set
(586,536 residues), with each candidate's *shuffled* self as a null - shuffling
preserves length and composition while destroying homology::

                            real (n=5)   shuffled (n=10)
    mean identity %         5.9          5.9
    95th pct identity %     13.5         13.5
    max identity %          35.9         34.5
    longest exact match     5.0 aa       5.0 aa

    max identity >= 30%:    real 5/5     shuffled 8/10

A 30% cut flags composition-matched noise about as readily as the real thing,
so at this length it is measuring length and amino-acid composition, not
homology. That is why ``drop_on_identity`` defaults to false: the screen is
computed and recorded, and the number is there to be read, but it does not
remove a candidate on its own. Raising the threshold or lengthening the query
is what makes it mean something - the FAO/WHO percentage criterion pairs 35%
with an 80-residue window for this reason, and that window is longer than
anything stage 5 emits.

Both criteria reuse :class:`mpvpredpip.autoimmunity.SequenceBlob`, which is the
same index the human-proteome screen runs on.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import autoimmunity


class AllergenSet(autoimmunity.SequenceBlob):
    """A FASTA of known allergens, indexed for both criteria."""

    def __init__(self, path: str) -> None:
        super().__init__(
            path,
            empty_message=(
                f"No sequences in the allergen FASTA at {path!r}.\n"
                "screening.allergen_homology.allergens in config.yaml must point at "
                "a FASTA of known allergens, e.g. workflow/allergens_dataset.fasta."
            ),
        )


@dataclass
class IdentityHit:
    """Best ungapped identity found for one query over the whole set."""

    percent: float = 0.0
    allergen_id: str = ""

    def describe(self) -> str:
        if not self.allergen_id:
            return ""
        return f"{self.percent:.1f}%@{self.allergen_id}"


@dataclass
class HomologyResult:
    """Both criteria's outcome for one consensus candidate."""

    #: Exact-match hits, keyed by the constituent epitope that carried them.
    exact: dict[str, list[autoimmunity.Match]] = field(default_factory=dict)
    #: Windows skipped because they held a non-standard residue.
    skipped_windows: list[str] = field(default_factory=list)
    #: Best percent identity for the consensus sequence.
    identity: IdentityHit = field(default_factory=IdentityHit)

    @property
    def exact_flagged(self) -> bool:
        return any(self.exact.values())

    @property
    def n_exact_epitopes(self) -> int:
        """How many constituent epitopes carried a match."""
        return sum(1 for hits in self.exact.values() if hits)

    def describe_exact(self, limit: int = 8) -> str:
        seen: set[str] = set()
        for hits in self.exact.values():
            for hit in hits:
                seen.add(hit.describe())
        return ";".join(sorted(seen)[:limit])


def max_identity_matches(query: str, target: str) -> int:
    """Most identical residues over any full-length ungapped placement of *query*.

    An offset histogram rather than a comparison at every offset: for each query
    residue, every position of that residue in *target* votes for the one
    alignment offset that would pair them. The best offset's tally is the match
    count. That is O(len(target) x len(query) / 20) on average against
    O(len(target) x len(query)) for the naive scan, which is what makes an
    exhaustive sweep of a 586,000-residue set affordable in-process.

    Only placements where *query* fits entirely inside *target* are considered.
    A partial overhang would be scored against the full query length and so
    report an identity the alignment does not support.
    """
    span, length = len(target), len(query)
    if length == 0 or span < length:
        return 0

    positions: dict[str, list[int]] = {}
    for index, residue in enumerate(target):
        positions.setdefault(residue, []).append(index)

    limit = span - length
    tally: dict[int, int] = {}
    for index, residue in enumerate(query):
        for found in positions.get(residue, ()):
            offset = found - index
            if 0 <= offset <= limit:
                tally[offset] = tally.get(offset, 0) + 1
    return max(tally.values()) if tally else 0


def best_identity(queries: dict[str, str], allergens: AllergenSet) -> dict[str, IdentityHit]:
    """Best percent identity for each query over every allergen.

    Loops allergens outermost and queries innermost on purpose: the residue
    index for a given allergen is built once and reused by every query, instead
    of being rebuilt per (allergen, query) pair.
    """
    best = {key: IdentityHit() for key in queries}
    prepared = {key: value.upper() for key, value in queries.items() if value}

    for allergen_id, target in allergens.records:
        for key, query in prepared.items():
            if len(target) < len(query):
                continue
            matches = max_identity_matches(query, target)
            if not matches:
                continue
            percent = 100.0 * matches / len(query)
            if percent > best[key].percent:
                best[key] = IdentityHit(percent=percent, allergen_id=allergen_id)
    return best


def screen(
    consensus: dict[str, str],
    members: dict[str, dict[str, str]],
    allergens: AllergenSet,
    exact_match_length: int,
    identity_threshold: float,
) -> dict[str, HomologyResult]:
    """Run both criteria.

    *consensus* is ``{consensus_id: sequence}``; *members* is
    ``{consensus_id: {epitope_id: sequence}}``, the constituent epitopes. Every
    key of *consensus* appears in the result whether it was flagged or not - the
    caller needs a verdict for every candidate, and an absent key could not be
    told apart from one that was never screened.
    """
    results = {key: HomologyResult() for key in consensus}

    # --- criterion 1, on the constituent epitopes ---------------------------
    for consensus_id in consensus:
        epitopes = members.get(consensus_id, {})
        if not epitopes:
            continue
        screened = autoimmunity.find_matches(epitopes, allergens, exact_match_length)
        for epitope_id, outcome in screened.items():
            results[consensus_id].exact[epitope_id] = outcome.matches
            results[consensus_id].skipped_windows.extend(outcome.skipped_windows)

    # --- criterion 2, on the consensus sequence ----------------------------
    for consensus_id, hit in best_identity(consensus, allergens).items():
        results[consensus_id].identity = hit

    return results


def identity_flagged(result: HomologyResult, identity_threshold: float) -> bool:
    """Whether the consensus met or exceeded the identity threshold.

    Stated as >= rather than >, so a threshold of 30 flags a candidate at
    exactly 30.0%. The boundary belongs to the flagged side: the threshold is
    the point at which the screen is meant to have something to say.
    """
    return result.identity.percent >= identity_threshold
