"""Stage 5 - cluster epitopes and select B/T consensus peptides.

Clustering is done by :func:`cluster_peptides` in this module, not by the IEDB
cluster standalone. That is a deliberate change, and the reason matters:

  The pipeline was written against the IEDB **web** cluster tool, which returns a
  consensus sequence and a gapped alignment per member. The **standalone**
  (IEDB_Cluster-1.0) is a different and much smaller program - two files - that
  only groups peptides. It emits
  ``cluster_number, num_epitopes_in_cluster, epitope_num, epitope_name, epitope_seq``
  and the words "consensus" and "alignment" appear nowhere in it. It is also
  Python 2 source and will not import under Python 3.

  Stage 5 exists to produce consensus peptides and stages 6-9 screen and rank
  them, so that consensus has to come from somewhere. It is computed here, from
  the tool's own identity rule, which :func:`percent_identity` reimplements
  faithfully from the standalone's source. This project already reimplements the
  IEDB conservancy algorithm for the same reason, so the precedent is not new.

:func:`parse` is kept for reading a file produced by the IEDB web tool, whose
output looks like this - one block per sub-cluster, a ``Consensus`` row carrying
the merged sequence, then one row per member with its gapped alignment::

    Cluster.Sub-Cluster Number,Peptide Number,Alignment,Position,Description,Peptide
    1.1,Consensus,KKDAPYIVGDVVQEGVLTAVVIPTKKA,-,-,-
    1.1,1,KKDAPYIVGDVVQEGV-----------,1,Spike_B_cell_BepiPred_4.2,KKDAPYIVGDVVQEGV
    1.1,2,----------VVQEGVLTAVVIPTKKA,11,Spike_T_cell_MHCI_23.1,VVQEGVLTAVVIPTKKA

Selection rule, from the pipeline specification: a cluster is carried forward
only if it contains at least one MHC class I epitope, at least one MHC class II
epitope, and at least one B-cell epitope from BepiPred *and* one from EpiDope.
Membership is resolved by joining each ``Description`` against the epitope table
from stages 3 and 4, so the rule is enforced on recorded provenance rather than
on string matching against the identifier. Parsing the identifier is kept only as
a fallback for rows the join cannot resolve.
"""

from __future__ import annotations

import csv
import os
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

from . import epitopes as epitope_mod

CONSENSUS_MARKERS = {"consensus", "consensus sequence"}

#: Written into a consensus column whose members disagree, or that nothing
#: covers. config.yaml's ``reject_ambiguous_consensus`` discards any consensus
#: containing it, reproducing the previous pipeline's "'X' not in alignment
#: column" check.
AMBIGUOUS = "X"


@dataclass
class Member:
    cluster_id: str
    member_index: str
    alignment: str
    position: str
    description: str
    sequence: str
    source_tool: str = ""
    epitope_class: str = ""
    #: Pipe-delimited species this peptide was predicted in, e.g.
    #: "Sudan_ebolavirus|Zaire_ebolavirus". Plural because deduplication collapses
    #: a peptide found in several strains into one clustering input, and which
    #: species a consensus actually covers is the point of a multi-species design.
    strains: str = ""


@dataclass
class Cluster:
    cluster_id: str
    consensus: str
    members: list[Member] = field(default_factory=list)

    def counts(self) -> Counter:
        return Counter(member.source_tool for member in self.members if member.source_tool)


def percent_identity(first: str, second: str) -> float:
    """Percent identity between two peptides, as the IEDB cluster tool defines it.

    Reimplemented from IEDB_Cluster-1.0, ``predict_cluster.py`` - functions
    ``homologous_sequence`` and ``compute_homology``. Three properties of that
    definition are easy to get wrong by assuming a conventional identity measure,
    and all three are deliberate here:

    * the shorter peptide is slid along the longer one with **no gaps**, and the
      best offset wins. There is no alignment and no substitution matrix;
    * the match count is divided by the length of the **longer** peptide, not by
      the shorter one that was slid. A short peptide sitting inside a long one
      therefore scores *low*, not high;
    * except that an exact substring short-circuits to 99.9 - never 100 - so two
      peptides of different lengths are never reported as identical.

    The asymmetry in the second and third points is the tool's, not a mistake in
    transcription: without the substring rule a 9-mer inside a 15-mer would score
    60% and usually fall below the threshold.
    """
    short, long_ = (first, second) if len(first) <= len(second) else (second, first)
    if not short or not long_:
        return 0.0
    if short in long_:
        return 99.9
    best = 0
    for offset in range(len(long_) - len(short) + 1):
        matches = sum(
            1 for index, residue in enumerate(short) if residue == long_[offset + index]
        )
        if matches > best:
            best = matches
    return 100.0 * best / len(long_)


def _components(peptides: Sequence[str], threshold: float) -> list[list[int]]:
    """Group peptide indices into connected components of the homology graph.

    The original builds, for each peptide, the list of its homologues, then
    repeatedly merges any two lists sharing a member until none do. Every such
    list contains its own seed, so two lists share a member exactly when their
    seeds are joined by a path of at-or-above-threshold pairs - that is, the
    result is the connected components of the homology graph. Union-find
    computes them directly; the original re-scans every pair on every merge
    round.
    """
    parent = list(range(len(peptides)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left, right = find(left), find(right)
        if left != right:
            # Always point at the lower index so components stay keyed by their
            # earliest member, which keeps the output order input order.
            parent[max(left, right)] = min(left, right)

    for i in range(len(peptides)):
        for j in range(i + 1, len(peptides)):
            if percent_identity(peptides[i], peptides[j]) >= threshold:
                union(i, j)

    groups: dict[int, list[int]] = {}
    for index in range(len(peptides)):
        groups.setdefault(find(index), []).append(index)
    return [groups[key] for key in sorted(groups)]


def _best_offset(anchor: str, peptide: str) -> int:
    """Where *peptide* sits relative to *anchor*, by best ungapped overlap.

    May be negative, meaning the peptide starts before the anchor does. Ties are
    broken toward the larger overlap and then the smaller shift, so the result
    depends only on the two sequences and not on iteration order.
    """
    best_key: tuple[int, int, int] | None = None
    best_offset = 0
    for offset in range(-len(peptide) + 1, len(anchor)):
        score = overlap = 0
        for index, residue in enumerate(peptide):
            position = offset + index
            if 0 <= position < len(anchor):
                overlap += 1
                if anchor[position] == residue:
                    score += 1
        if overlap == 0:
            continue
        key = (score, overlap, -abs(offset))
        if best_key is None or key > best_key:
            best_key, best_offset = key, offset
    return best_offset


def build_consensus(peptides: Sequence[str]) -> tuple[str, list[str], list[int]]:
    """Merge a cluster's peptides into one consensus.

    Returns ``(consensus, alignments, positions)``, the latter two in the order
    the peptides were given. Each alignment is the peptide padded with ``-`` to
    the consensus width, and each position is its 1-based start in it - the same
    shape the IEDB web tool reports, so the rest of stage 5 is unchanged.

    Every peptide is placed against the longest one, and the consensus spans
    from the leftmost start to the rightmost end. Each column is decided by
    **majority vote**: the most common residue wins, and only a genuine tie - or
    a column nothing covers - becomes :data:`AMBIGUOUS`.

    It used to require unanimity, and that was wrong for this project. The
    workflow this reimplements analysed a single species, where demanding no
    ambiguity cost nothing. Here three Ebola species are compared on purpose, so
    unanimity puts an X at every position where any two of them differ - which is
    the whole reason a cross-species consensus is interesting. It rejected the
    only Spike cluster that satisfied every other requirement, whose consensus
    was ``DXXLPXQXXXDNWWTGWRQWXPAGIG``.

    A majority call is a real, synthesisable peptide, which an X-laden string is
    not - nothing downstream can score or screen an X. How well that peptide
    actually covers each strain is a question of degree, and stage 7's
    conservancy analysis measures it directly; it is not something a hard gate
    here can express. Use :func:`disagreeing_positions` to see how much voting a
    given consensus required.
    """
    if not peptides:
        return "", [], []

    anchor = max(peptides, key=len)
    offsets = [_best_offset(anchor, peptide) for peptide in peptides]
    start = min(offsets)
    width = max(
        offset + len(peptide) for offset, peptide in zip(offsets, peptides)
    ) - start

    columns: list[list[str]] = [[] for _ in range(width)]
    for offset, peptide in zip(offsets, peptides):
        for index, residue in enumerate(peptide):
            columns[offset - start + index].append(residue)

    residues = []
    for column in columns:
        if not column:
            # Nothing covers this position. Only reachable if a member aligned
            # far enough off the anchor to leave a hole, which the span
            # calculation above otherwise prevents.
            residues.append(AMBIGUOUS)
            continue
        ranked = Counter(column).most_common()
        # A tie is the one case with no defensible answer: two residues are
        # equally supported and picking either would invent evidence.
        if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
            residues.append(AMBIGUOUS)
        else:
            residues.append(ranked[0][0])
    consensus = "".join(residues)

    alignments = []
    positions = []
    for offset, peptide in zip(offsets, peptides):
        lead = offset - start
        alignments.append("-" * lead + peptide + "-" * (width - lead - len(peptide)))
        positions.append(lead + 1)
    return consensus, alignments, positions


def disagreeing_positions(cluster: Cluster) -> int:
    """How many consensus positions the members did not all agree on.

    With majority voting the consensus no longer advertises this - a position
    where two strains differ now shows the winning residue rather than an X - so
    it is reported separately. A consensus needing many votes still covers its
    members less well than one that needed none, and that difference is invisible
    in the sequence alone.

    Reconstructed from the members' stored alignments, which are already padded
    to the consensus width, so this costs nothing to keep.
    """
    count = 0
    for index in range(len(cluster.consensus)):
        seen = {
            member.alignment[index]
            for member in cluster.members
            if index < len(member.alignment)
        }
        seen.discard("-")
        if len(seen) > 1:
            count += 1
    return count


def member_matches_consensus(cluster: Cluster, member: Member) -> bool:
    """Does this member's sequence appear verbatim in the consensus it sits under?

    The consensus is a per-column **majority vote**, so a residue only two of
    three members agree on is written into it and the third member is overruled.
    The result is a sequence that need not exist in any strain, carrying members
    that do not actually contain it: in the Ebola spike run, cluster 32.1's
    consensus took S at position 10 where many members have A, and 21 of its 34
    members did not match the consensus they were filed under.

    That matters because every downstream stage treats the consensus as the
    candidate and its members as the evidence for it. A member that differs from
    the consensus is evidence for a *different* peptide.

    An unreadable position counts as a mismatch rather than raising: the member
    cannot be shown to support the consensus, which is the question being asked.
    """
    try:
        start = int(member.position) - 1
    except (TypeError, ValueError):
        return False
    if start < 0:
        return False
    return cluster.consensus[start:start + len(member.sequence)] == member.sequence


def prune_to_consensus(cluster: Cluster) -> list[Member]:
    """Drop members that are not an exact match to their consensus; return them.

    Mutates ``cluster.members``. Everything downstream - the composition counts,
    the strain set, the disagreement count, pass/fail, and the membership table -
    is then computed from the members that genuinely support the consensus.

    The consensus itself is deliberately NOT rebuilt from the survivors. Doing so
    would change it, which could orphan further members, and the loop has no
    natural fixed point. The consensus stays as the clustering produced it; what
    changes is that the evidence recorded against it is now evidence for it.
    """
    kept: list[Member] = []
    dropped: list[Member] = []
    for member in cluster.members:
        (kept if member_matches_consensus(cluster, member) else dropped).append(member)
    cluster.members = kept
    return dropped


def cluster_peptides(
    entries: Sequence[tuple[str, str]],
    threshold: float,
    unique: bool = True,
) -> list[Cluster]:
    """Cluster ``(identifier, peptide)`` pairs and build each cluster's consensus.

    *threshold* is a **percentage**, 0-100, matching the IEDB tool's own scale
    where the default is 80. It is not a fraction; passing 0.6 would cluster
    every peptide in the run into a single blob, so the caller is expected to
    have rejected that - see ``workflow/scripts/05_run_cluster.py``.

    *unique* mirrors the tool's ``--unique=on`` default: identical peptides are
    collapsed to their first occurrence before clustering. Two tools predicting
    the same peptide is a real occurrence, not a duplicate record, so the member
    that survives keeps only one identifier - which is why stage 5's composition
    rule is evaluated on the epitope table rather than on cluster membership
    alone.
    """
    if unique:
        collapsed: dict[str, tuple[str, str]] = {}
        for identifier, peptide in entries:
            collapsed.setdefault(peptide, (identifier, peptide))
        entries = list(collapsed.values())

    peptides = [peptide for _, peptide in entries]
    clusters: list[Cluster] = []

    for number, indices in enumerate(_components(peptides, threshold), start=1):
        group = [entries[index] for index in indices]
        consensus, alignments, positions = build_consensus([p for _, p in group])
        # "N.1": the second component is the sub-cluster, which the IEDB web tool
        # uses when it splits a cluster further. Nothing here sub-clusters, so it
        # is always 1 - kept so identifiers stay the shape the rest of the
        # pipeline and its tests already expect.
        cluster_id = f"{number}.1"
        cluster = Cluster(cluster_id=cluster_id, consensus=consensus)
        # Members are emitted in consensus order - ascending start position -
        # rather than in whatever order the clustering happened to visit them, so
        # the table reads as an alignment block: every member sits under the slice
        # of the consensus it covers. member_index is assigned after the sort, so
        # it numbers the rows as they appear.
        #
        # Python's sort is stable, so peptides starting at the same position keep
        # the order they were given in - which is the order of the stage 3/4
        # tables, and the only ordering those rows ever had.
        ordered = sorted(zip(group, alignments, positions), key=lambda item: item[2])
        for member_index, ((identifier, peptide), alignment, position) in enumerate(
            ordered, start=1
        ):
            cluster.members.append(
                Member(
                    cluster_id=cluster_id,
                    member_index=str(member_index),
                    alignment=alignment,
                    position=str(position),
                    description=identifier,
                    sequence=peptide,
                )
            )
        clusters.append(cluster)
    return clusters


CLUSTER_HEADER = [
    "Cluster.Sub-Cluster Number",
    "Peptide Number",
    "Alignment",
    "Position",
    "Description",
    "Peptide",
]


def write_clusters(path: str, clusters: Sequence[Cluster]) -> None:
    """Write clusters in the IEDB web tool's own layout.

    Deliberately the format :func:`parse` reads, rather than anything neater.
    Stage 5 stays two rules with a file between them, as every other stage
    boundary in this pipeline is a file: only the producer has changed. It also
    means the raw output can still be compared against a run of the IEDB web
    tool on the same input, which is the one check that would catch this
    module's clustering drifting away from the tool it reimplements.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(CLUSTER_HEADER)
        for cluster in clusters:
            writer.writerow([cluster.cluster_id, "Consensus", cluster.consensus, "-", "-", "-"])
            for member in cluster.members:
                writer.writerow([
                    member.cluster_id,
                    member.member_index,
                    member.alignment,
                    member.position,
                    member.description,
                    member.sequence,
                ])


def parse(path: str) -> list[Cluster]:
    """Read an IEDB cluster output file. Accepts comma- or tab-delimited input."""
    if not os.path.exists(path):
        raise SystemExit(f"Cluster output not found: {path}")

    with open(path, encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "\t" if sample.count("\t") > sample.count(",") else ","
        rows = [row for row in csv.reader(handle, delimiter=delimiter) if row]

    clusters: list[Cluster] = []
    current: Cluster | None = None

    for row in rows:
        cells = [cell.strip() for cell in row]
        if len(cells) < 3:
            continue
        # Skip the header wherever it appears.
        if cells[0].lower().startswith("cluster") and "number" in cells[1].lower():
            continue

        cluster_id, marker = cells[0], cells[1]

        if marker.lower() in CONSENSUS_MARKERS:
            current = Cluster(cluster_id=cluster_id, consensus=cells[2])
            clusters.append(current)
            continue

        if current is None or cluster_id != current.cluster_id:
            # A member row without a preceding consensus row: the IEDB tool emits
            # singletons this way. Treat the peptide as its own consensus.
            current = Cluster(cluster_id=cluster_id, consensus=cells[5] if len(cells) > 5 else cells[2])
            clusters.append(current)

        current.members.append(
            Member(
                cluster_id=cluster_id,
                member_index=marker,
                alignment=cells[2] if len(cells) > 2 else "",
                position=cells[3] if len(cells) > 3 else "",
                description=cells[4] if len(cells) > 4 else "",
                sequence=cells[5] if len(cells) > 5 else "",
            )
        )
    return clusters


def annotate(clusters: Iterable[Cluster], epitope_rows: Sequence[dict[str, str]]) -> None:
    """Attach ``source_tool``/``epitope_class``/``strains`` to every member, in place."""
    lookup = {str(row["epitope_id"]): row for row in epitope_rows}
    # Secondary index so a member can still be resolved if the description was
    # truncated or reformatted but the peptide survived intact.
    by_sequence: dict[str, dict[str, str]] = {}
    for row in epitope_rows:
        by_sequence.setdefault(str(row["sequence"]), row)

    # Strains are collected by peptide, across every row, not taken from the one
    # row the identifier resolves to. Deduplication keeps a single representative
    # identifier for a peptide predicted in several species, so that row names
    # only the first of them - reading strain from it would silently attribute a
    # peptide shared by all three species to whichever sorted first.
    strains_by_sequence: dict[str, set[str]] = {}
    for row in epitope_rows:
        sequence = str(row.get("sequence", ""))
        if not sequence:
            continue
        found = strains_by_sequence.setdefault(sequence, set())
        # EPITOPES has "strain"; the deduplicated EPITOPES_NR has pipe-joined
        # "strains". Accept either, so this works whichever table is passed.
        for key in ("strain", "strains"):
            value = str(row.get(key, "") or "").strip()
            if value:
                found.update(part for part in value.split("|") if part)

    for cluster in clusters:
        for member in cluster.members:
            member.strains = "|".join(sorted(strains_by_sequence.get(member.sequence, ())))

            row = lookup.get(member.description)
            if row is None and member.sequence:
                row = by_sequence.get(member.sequence)
            if row is not None:
                member.source_tool = str(row["source_tool"])
                member.epitope_class = str(row["epitope_class"])
                continue

            parsed = epitope_mod.parse_id(member.description)
            if parsed:
                _, member.source_tool, member.epitope_class = parsed


def cluster_strains(cluster: Cluster) -> list[str]:
    """Every species represented among a cluster's members, sorted.

    This is what makes a consensus peptide a *multi-species* candidate or not.
    A cluster whose members all come from one species yields a consensus for that
    species alone, which is a materially different result from one spanning all
    three - and before this was recorded, the two were indistinguishable in every
    stage 5 output.
    """
    found: set[str] = set()
    for member in cluster.members:
        found.update(part for part in member.strains.split("|") if part)
    return sorted(found)


def evaluate(
    cluster: Cluster,
    require: dict[str, int],
    require_both_bcell_tools: bool = True,
    reject_ambiguous: bool = True,
) -> tuple[bool, str, Counter]:
    """Apply the selection rule. Returns ``(passes, reason, counts)``."""
    counts = cluster.counts()
    reasons: list[str] = []

    if reject_ambiguous and "X" in cluster.consensus.upper():
        reasons.append("consensus contains an ambiguous residue (X)")

    n_mhc_i = counts.get("mhc_i", 0)
    n_mhc_ii = counts.get("mhc_ii", 0)
    n_bepipred = counts.get("bepipred", 0)
    n_epidope = counts.get("epidope", 0)

    if n_mhc_i < require.get("min_mhc_i", 1):
        reasons.append(f"only {n_mhc_i} MHC-I epitope(s), needs {require.get('min_mhc_i', 1)}")
    if n_mhc_ii < require.get("min_mhc_ii", 1):
        reasons.append(f"only {n_mhc_ii} MHC-II epitope(s), needs {require.get('min_mhc_ii', 1)}")

    if require_both_bcell_tools:
        if n_bepipred < require.get("min_bepipred", 1):
            reasons.append(
                f"only {n_bepipred} BepiPred epitope(s), needs {require.get('min_bepipred', 1)}"
            )
        if n_epidope < require.get("min_epidope", 1):
            reasons.append(
                f"only {n_epidope} EpiDope epitope(s), needs {require.get('min_epidope', 1)}"
            )
    else:
        needed = require.get("min_bepipred", 1) + require.get("min_epidope", 1)
        if n_bepipred + n_epidope < needed:
            reasons.append(
                f"only {n_bepipred + n_epidope} B-cell epitope(s) across both tools, needs {needed}"
            )

    # Strain coverage. A cluster whose members all come from one species yields a
    # consensus for that species alone; requiring every strain is what makes the
    # result a genuinely multi-species candidate. Default 0 keeps this off, so a
    # config written before this option behaves exactly as it did.
    min_strains = require.get("min_strains", 0)
    if min_strains:
        n_strains = len(cluster_strains(cluster))
        if n_strains < min_strains:
            reasons.append(
                f"only {n_strains} strain(s) represented, needs {min_strains}"
            )

    unresolved = sum(1 for member in cluster.members if not member.source_tool)
    if unresolved:
        reasons.append(
            f"{unresolved} member(s) could not be traced back to a predictor - "
            "check that the clustering input was built from the stage 3/4 tables"
        )

    return (not reasons), "; ".join(reasons), counts


#: How far the sweep drops between thresholds, in identity percentage points.
#: Fixed rather than configurable: "90:50" is meant to read as a range, and a
#: third number would make the common case harder to write for a step nobody
#: has asked to vary. The minimum is always included even when the step does
#: not land on it, so "90:55" ends at 55 rather than overshooting to 50.
SWEEP_STEP = 10


def parse_threshold_sweep(value: object) -> list[int]:
    """Read ``clustering.threshold`` into the list of thresholds to run.

    ``"90:50"`` sweeps 90, 80, 70, 60, 50 - tightest first. A single number is
    still accepted and means one threshold, so a config written before the
    sweep existed behaves exactly as it did.

    Descending order is load-bearing downstream: when the same consensus
    sequence passes at several thresholds the pooled output keeps the first
    one, which is then the tightest clustering that produced it.
    """
    text = str(value).strip()
    if ":" not in text:
        return [_valid_threshold(text)]

    parts = text.split(":")
    if len(parts) != 2:
        raise ValueError(
            f"clustering.threshold {value!r} is not a range. Use 'MAX:MIN', "
            "for example '90:50', or a single number."
        )
    high, low = (_valid_threshold(part) for part in parts)
    if high < low:
        raise ValueError(
            f"clustering.threshold {value!r} runs backwards: it is written "
            f"MAX:MIN and sweeps from tight to loose, so {high} must not be "
            f"below {low}."
        )

    sweep = list(range(high, low, -SWEEP_STEP))
    if not sweep or sweep[-1] != low:
        sweep.append(low)
    return sweep


def _valid_threshold(text: str) -> int:
    """One threshold on the IEDB standalone's 0-100 scale."""
    try:
        value = float(text)
    except (TypeError, ValueError):
        raise ValueError(
            f"clustering.threshold {text!r} is not a number."
        ) from None
    # The same guard 05_run_cluster.py applies: the scale is a percentage, and
    # a value at or below 1 is a fraction left over from the old config, which
    # would cluster every epitope in the run together.
    if not 1 < value <= 100:
        raise ValueError(
            f"clustering.threshold {text!r} is out of range. It is a PERCENTAGE "
            "on the IEDB standalone's 0-100 scale, so give 60 rather than 0.6."
        )
    return int(value)


def _position_key(member: Member) -> tuple[int, int]:
    """Sort key putting members in ascending consensus position.

    Anything unparseable sorts last rather than raising: ``parse`` accepts a
    cluster file this pipeline did not necessarily write, and an odd value in
    one column is not a reason to lose the whole stage.
    """
    try:
        return (0, int(member.position))
    except (TypeError, ValueError):
        return (1, 0)


@dataclass
class ThresholdResult:
    """What one identity threshold produced, before anything is written."""

    threshold: int
    cluster_rows: list[dict[str, object]] = field(default_factory=list)
    member_rows: list[dict[str, object]] = field(default_factory=list)
    #: ``(consensus_id, sequence)`` for the clusters that passed.
    passing: list[tuple[str, str]] = field(default_factory=list)
    #: Clusters the tool produced, before the ambiguous ones were discarded.
    n_found: int = 0
    #: Discarded for an X in the consensus, so absent from every row above.
    n_ambiguous: int = 0
    #: Rejection reason -> how many clusters gave it. One cluster can fail on
    #: several criteria and is counted under each.
    reasons: Counter = field(default_factory=Counter)

    @property
    def n_examined(self) -> int:
        return self.n_found - self.n_ambiguous

    @property
    def n_passing(self) -> int:
        return len(self.passing)


def select_at_threshold(
    cluster_path: str,
    epitope_rows: Sequence[dict[str, object]],
    protein_label: str,
    threshold: int,
    require: dict[str, int],
    require_both_bcell_tools: bool = True,
    reject_ambiguous: bool = True,
) -> ThresholdResult:
    """Apply the selection rule to one threshold's clustering.

    Returns the rows rather than writing them: the caller runs this once per
    threshold in the sweep and pools the results, and pooling needs the rows in
    memory anyway.

    *epitope_rows* is passed in already read, because the sweep would otherwise
    re-read the same table once per threshold.

    A cluster whose consensus contains an X is **discarded outright** when
    *reject_ambiguous* is set - it contributes no cluster row, no member rows
    and no summary line, only to :attr:`ThresholdResult.n_ambiguous`. It used
    to be kept and marked as failing, which put a quarter of all clusters into
    the output purely to say they were unusable. Discarding them before the
    members are pruned also skips the most expensive step for the clusters that
    can never be used.
    """
    result = ThresholdResult(threshold=threshold)

    clusters = parse(cluster_path)
    annotate(clusters, epitope_rows)
    result.n_found = len(clusters)

    if reject_ambiguous:
        usable = [c for c in clusters if AMBIGUOUS not in c.consensus.upper()]
        result.n_ambiguous = len(clusters) - len(usable)
        clusters = usable

    # A cluster's consensus is a majority vote and can differ from the members
    # that produced it; those members are not evidence for it. Pruning here -
    # rather than at write time - means the composition counts, strain
    # coverage, pass/fail, the pooled FASTA and the summary are all derived
    # from members that match their consensus exactly.
    dropped_counts = {
        cluster.cluster_id: len(prune_to_consensus(cluster)) for cluster in clusters
    }

    for cluster in clusters:
        passes, reason, counts = evaluate(
            cluster, require, require_both_bcell_tools, reject_ambiguous
        )
        # The threshold is part of the identity, not decoration: cluster ids
        # restart at 1.1 for every threshold, so without it the sweep would
        # produce several different peptides all called
        # "Nucleoprotein_Consensus_1.1" and stage 6 onwards would join them
        # together.
        consensus_id = f"{protein_label}_Consensus_{threshold}_{cluster.cluster_id}"
        strains = cluster_strains(cluster)

        if not passes:
            for part in reason.split("; "):
                if part:
                    result.reasons[part] += 1

        result.cluster_rows.append(
            {
                "protein_label": protein_label,
                "threshold": threshold,
                "cluster_id": cluster.cluster_id,
                "consensus_id": consensus_id,
                "consensus_sequence": cluster.consensus,
                "consensus_length": len(cluster.consensus),
                "n_members": len(cluster.members),
                "n_members_dropped": dropped_counts.get(cluster.cluster_id, 0),
                "n_bepipred": counts.get("bepipred", 0),
                "n_epidope": counts.get("epidope", 0),
                "n_mhc_i": counts.get("mhc_i", 0),
                "n_mhc_ii": counts.get("mhc_ii", 0),
                "strains": "|".join(strains),
                "n_strains": len(strains),
                "n_disagreeing": disagreeing_positions(cluster),
                "passes": passes,
                "reason": reason,
            }
        )

        # The consensus itself, ahead of its members, so the file reads as an
        # alignment block: the consensus on top and each member padded to the same
        # width beneath it. Stages 8 and 9 join these rows on epitope_id to pull
        # alleles through; "-" matches no epitope, so both skip this row rather
        # than counting it as a member.
        result.member_rows.append(
            {
                "protein_label": protein_label,
                "threshold": threshold,
                "cluster_id": cluster.cluster_id,
                "consensus_id": consensus_id,
                "member_index": "Consensus",
                "epitope_id": "-",
                "source_tool": "-",
                "epitope_class": "-",
                "strains": "|".join(strains),
                "position": "-",
                "alignment": cluster.consensus,
                "sequence": cluster.consensus,
            }
        )

        # Ascending consensus position, so each block reads as an alignment under
        # its consensus. cluster_peptides already orders them this way, but this
        # must not rely on that: parse() accepts a cluster file this pipeline did
        # not write - the IEDB standalone's own output - which carries whatever
        # order that file had. member_index is the row's ordinal in this order,
        # so the numbering always matches the rows as printed.
        for member_index, member in enumerate(
            sorted(cluster.members, key=_position_key), start=1
        ):
            result.member_rows.append(
                {
                    "protein_label": protein_label,
                    "threshold": threshold,
                    "cluster_id": cluster.cluster_id,
                    "consensus_id": consensus_id,
                    "member_index": str(member_index),
                    "epitope_id": member.description,
                    "source_tool": member.source_tool,
                    "epitope_class": member.epitope_class,
                    "strains": member.strains,
                    "position": member.position,
                    "alignment": member.alignment,
                    "sequence": member.sequence,
                }
            )

        if passes:
            result.passing.append((consensus_id, cluster.consensus))

    return result


def pool_passing(
    results: Sequence[ThresholdResult],
) -> tuple[list[tuple[str, str]], int]:
    """Collapse the sweep's passing peptides to one entry per sequence.

    Returns ``(entries, n_duplicates)``. A peptide that survives at 90% is
    usually still there at 80% and below, so without this the same sequence
    would be screened several times over and appear repeatedly among the final
    candidates.

    *results* is expected tightest-first, and the first occurrence wins, so the
    survivor is attributed to the tightest clustering that produced it.
    """
    seen: dict[str, str] = {}
    duplicates = 0
    for result in results:
        for consensus_id, sequence in result.passing:
            if sequence in seen:
                duplicates += 1
                continue
            seen[sequence] = consensus_id
    return [(consensus_id, sequence) for sequence, consensus_id in seen.items()], duplicates
