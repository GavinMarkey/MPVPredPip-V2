"""Tests for the pipeline logic that does not depend on any external tool.

Run with:

    python -m pytest tests/ -v

These cover the parts that decide which candidates survive - interval building,
identifier assignment, the stage 5 composition rule and conservancy - so a change
to a threshold or a parser can be checked without installing a single predictor.
"""

from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from mpvpredpip import algpred, autoimmunity, inputs  # noqa: E402
from mpvpredpip import cluster as cluster_mod  # noqa: E402
from mpvpredpip import conservancy, epitopes, fasta, tables  # noqa: E402
from mpvpredpip import manifest as manifest_mod  # noqa: E402
from mpvpredpip.adapters import bcell, iedb, screening  # noqa: E402

# ---------------------------------------------------------------------------
# per-residue scores -> epitope intervals
# ---------------------------------------------------------------------------

class TestResiduesToEpitopes:
    """Epitopes are fixed-length peptides tiled across each above-threshold run.

    Most tests use small `length`/`step` values for legibility; the production
    defaults are 14 and 2 (config.yaml), exercised by
    `test_real_world_tiling_matches_the_documented_example`.
    """

    def test_single_run_gives_one_peptide_when_it_fits_exactly(self):
        sequence = "A" * 10
        scores = [0.0, 0.0, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.0, 0.0]
        result = epitopes.residues_to_epitopes(
            sequence, scores, threshold=0.5, length=6, step=2
        )
        assert result == [(3, 8, pytest.approx(0.9))]

    def test_coordinates_are_one_based_inclusive(self):
        sequence = "ABCDE"
        scores = [0.9, 0.9, 0.9, 0.0, 0.0]
        (start, end, _), = epitopes.residues_to_epitopes(
            sequence, scores, threshold=0.5, length=3, step=2
        )
        assert (start, end) == (1, 3)
        assert sequence[start - 1 : end] == "ABC"

    def test_short_region_is_padded_out_with_flanking_sequence(self):
        # The Ebola nucleoprotein case: every above-threshold region is shorter
        # than the epitope length. Dropping them would remove the protein from
        # the analysis entirely, so the region is centred in a full-length window.
        sequence = "ABCDEFGHIJKLMNOP"
        scores = [0.0] * 16
        for index in range(6, 12):  # a 6-residue region at positions 7-12
            scores[index] = 0.9
        result = epitopes.residues_to_epitopes(
            sequence, scores, threshold=0.5, length=10, step=2, min_region=4
        )
        (start, end, _), = result
        assert end - start + 1 == 10
        assert start <= 7 and end >= 12  # the region is inside the window

    def test_padding_is_clamped_at_the_sequence_start(self):
        sequence = "ABCDEFGHIJ"
        scores = [0.9, 0.9, 0.9] + [0.0] * 7  # region at the very beginning
        (start, end, _), = epitopes.residues_to_epitopes(
            sequence, scores, threshold=0.5, length=6, step=2, min_region=3
        )
        assert (start, end) == (1, 6)

    def test_padding_is_clamped_at_the_sequence_end(self):
        sequence = "ABCDEFGHIJ"
        scores = [0.0] * 7 + [0.9, 0.9, 0.9]  # region at the very end
        (start, end, _), = epitopes.residues_to_epitopes(
            sequence, scores, threshold=0.5, length=6, step=2, min_region=3
        )
        assert (start, end) == (5, 10)

    def test_regions_below_min_region_are_dropped(self):
        scores = [0.9, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert (
            epitopes.residues_to_epitopes(
                "A" * 8, scores, threshold=0.5, length=5, step=2, min_region=4
            )
            == []
        )

    def test_nothing_is_emitted_when_the_sequence_is_shorter_than_the_epitope(self):
        scores = [0.9] * 5
        assert (
            epitopes.residues_to_epitopes(
                "ABCDE", scores, threshold=0.5, length=14, step=2, min_region=3
            )
            == []
        )

    def test_long_run_is_tiled_with_overlapping_peptides(self):
        # 10 above-threshold residues, length 6, step 2 -> offsets 1, 3, 5.
        scores = [0.9] * 10
        result = epitopes.residues_to_epitopes(
            "A" * 10, scores, threshold=0.5, length=6, step=2
        )
        assert [(s, e) for s, e, _ in result] == [(1, 6), (3, 8), (5, 10)]
        assert all(e - s + 1 == 6 for s, e, _ in result)

    def test_real_world_tiling_matches_the_documented_example(self):
        # The Bundibugyo spike epitope from the first EpiDope run, at the
        # production settings.
        sequence = "PHDWTKNITDKIDQIIHDFIDKPLPDQTDNDNWWTGWRQ"
        scores = [0.9] * len(sequence)
        result = epitopes.residues_to_epitopes(
            sequence, scores, threshold=0.5, length=14, step=2
        )
        peptides = [sequence[s - 1 : e] for s, e, _ in result]
        assert peptides[:3] == [
            "PHDWTKNITDKIDQ",
            "DWTKNITDKIDQII",
            "TKNITDKIDQIIHD",
        ]
        assert all(len(peptide) == 14 for peptide in peptides)

    def test_a_tail_left_by_the_stride_is_still_covered(self):
        # 11 residues, length 6, step 2 -> offsets 1, 3, 5 reach residue 10, so a
        # final window anchored at the end is added for residue 11.
        scores = [0.9] * 11
        result = epitopes.residues_to_epitopes(
            "A" * 11, scores, threshold=0.5, length=6, step=2
        )
        assert [(s, e) for s, e, _ in result] == [(1, 6), (3, 8), (5, 10), (6, 11)]
        assert result[-1][1] == 11

    def test_interior_gap_is_absorbed(self):
        # One sub-threshold residue between two runs should not split the region.
        scores = [0.9, 0.9, 0.9, 0.1, 0.9, 0.9, 0.9]
        result = epitopes.residues_to_epitopes(
            "A" * 7, scores, threshold=0.5, length=7, step=2, merge_gap=1
        )
        assert result == [(1, 7, pytest.approx(sum(scores) / 7))]

    def test_gap_wider_than_merge_gap_splits(self):
        scores = [0.9, 0.9, 0.9, 0.1, 0.1, 0.9, 0.9, 0.9]
        result = epitopes.residues_to_epitopes(
            "A" * 8, scores, threshold=0.5, length=3, step=2, merge_gap=1
        )
        assert [(s, e) for s, e, _ in result] == [(1, 3), (6, 8)]

    def test_leading_gap_is_not_absorbed(self):
        # A gap at the edge is not "between" two regions; absorbing it would
        # extend the prediction past its evidence.
        scores = [0.1, 0.9, 0.9, 0.9]
        result = epitopes.residues_to_epitopes(
            "ABCD", scores, threshold=0.5, length=3, step=2, merge_gap=1
        )
        assert [(s, e) for s, e, _ in result] == [(2, 4)]

    def test_every_epitope_has_exactly_the_configured_length(self):
        # The regression this replaced: a 42-residue run used to become one
        # 42-residue epitope. Length is now fixed, whatever the run length.
        scores = [0.9] * 42
        result = epitopes.residues_to_epitopes(
            "A" * 42, scores, threshold=0.5, length=14, step=2
        )
        assert all(e - s + 1 == 14 for s, e, _ in result)
        assert len(result) == 15  # offsets 1,3,...,29

    def test_step_of_one_gives_every_offset(self):
        scores = [0.9] * 8
        result = epitopes.residues_to_epitopes(
            "A" * 8, scores, threshold=0.5, length=6, step=1
        )
        assert [(s, e) for s, e, _ in result] == [(1, 6), (2, 7), (3, 8)]

    def test_invalid_length_or_step_is_rejected(self):
        scores = [0.9] * 8
        with pytest.raises(ValueError):
            epitopes.residues_to_epitopes("A" * 8, scores, threshold=0.5, length=0, step=2)
        with pytest.raises(ValueError):
            epitopes.residues_to_epitopes("A" * 8, scores, threshold=0.5, length=6, step=0)

    def test_threshold_is_inclusive(self):
        scores = [0.5, 0.5, 0.5]
        assert epitopes.residues_to_epitopes("ABC", scores, threshold=0.5, length=3, step=2)

    def test_length_mismatch_is_an_error(self):
        with pytest.raises(ValueError, match="does not match"):
            epitopes.residues_to_epitopes("ABC", [0.9, 0.9], threshold=0.5)


# ---------------------------------------------------------------------------
# identifiers
# ---------------------------------------------------------------------------

class TestIdentifiers:
    def _make(self, tool, sequence, strain="Zaire", start=1):
        return epitopes.Epitope(
            protein_label="Spike", strain=strain, accession="X1",
            source_tool=tool, start=start, end=start + len(sequence) - 1,
            sequence=sequence, score=1.0,
        )

    def test_format(self):
        assigned = epitopes.assign_ids([self._make("bepipred", "ACDEFGHI")])
        assert assigned[0].epitope_id == "Spike_B_cell_BepiPred_1.1"

    def test_class_infix_matches_tool(self):
        assigned = epitopes.assign_ids([
            self._make("bepipred", "AAAAAAAA"),
            self._make("epidope", "CCCCCCCC"),
            self._make("mhc_i", "DDDDDDDDD"),
            self._make("mhc_ii", "EEEEEEEEEEEEEEE"),
        ])
        ids = {e.source_tool: e.epitope_id for e in assigned}
        assert "_B_cell_BepiPred_" in ids["bepipred"]
        assert "_B_cell_EpiDope_" in ids["epidope"]
        assert "_T_cell_MHCI_" in ids["mhc_i"]
        assert "_T_cell_MHCII_" in ids["mhc_ii"]

    def test_repeat_peptide_increments_occurrence(self):
        assigned = epitopes.assign_ids([
            self._make("bepipred", "ACDEFGHI", strain="Zaire", start=1),
            self._make("bepipred", "ACDEFGHI", strain="Sudan", start=40),
        ])
        assert sorted(e.epitope_id for e in assigned) == [
            "Spike_B_cell_BepiPred_1.1",
            "Spike_B_cell_BepiPred_1.2",
        ]

    def test_distinct_peptides_get_distinct_indices(self):
        assigned = epitopes.assign_ids([
            self._make("bepipred", "AAAAAAAA"),
            self._make("bepipred", "CCCCCCCC"),
        ])
        assert len({e.epitope_id for e in assigned}) == 2

    def test_assignment_is_deterministic(self):
        def build():
            return epitopes.assign_ids([
                self._make("bepipred", "CCCCCCCC", strain="Sudan"),
                self._make("bepipred", "AAAAAAAA", strain="Zaire"),
            ])
        assert [e.epitope_id for e in build()] == [e.epitope_id for e in build()]

    @pytest.mark.parametrize("tool,expected", [
        ("Spike_B_cell_BepiPred_1.1", ("Spike", "bepipred", "B")),
        ("Spike_B_cell_EpiDope_2.1", ("Spike", "epidope", "B")),
        ("Nucleoprotein_T_cell_MHCI_3.2", ("Nucleoprotein", "mhc_i", "MHC-I")),
        ("Nucleoprotein_T_cell_MHCII_4.1", ("Nucleoprotein", "mhc_ii", "MHC-II")),
    ])
    def test_round_trip_parse(self, tool, expected):
        assert epitopes.parse_id(tool) == expected

    def test_parse_rejects_junk(self):
        assert epitopes.parse_id("not an epitope id") is None


# ---------------------------------------------------------------------------
# deduplication
# ---------------------------------------------------------------------------

class TestDeduplicate:
    def _row(self, tool, sequence, strain, epitope_id, allele=""):
        return {
            "epitope_id": epitope_id, "protein_label": "Spike", "strain": strain,
            "source_tool": tool, "sequence": sequence, "score": "1.0", "allele": allele,
        }

    def test_same_tool_same_peptide_collapses(self):
        result = epitopes.deduplicate([
            self._row("bepipred", "ACDEFGHI", "Zaire", "Spike_B_cell_BepiPred_1.1"),
            self._row("bepipred", "ACDEFGHI", "Sudan", "Spike_B_cell_BepiPred_1.2"),
        ])
        assert len(result) == 1
        assert result[0]["n_occurrences"] == 2
        assert result[0]["strains"] == "Sudan|Zaire"

    def test_different_tools_stay_separate(self):
        # Critical: stage 5 must be able to see that both B-cell tools support
        # the same peptide, so these must not collapse into one entry.
        result = epitopes.deduplicate([
            self._row("bepipred", "ACDEFGHI", "Zaire", "Spike_B_cell_BepiPred_1.1"),
            self._row("epidope", "ACDEFGHI", "Zaire", "Spike_B_cell_EpiDope_1.1"),
        ])
        assert len(result) == 2
        assert {row["source_tool"] for row in result} == {"bepipred", "epidope"}

    def test_alleles_are_aggregated(self):
        result = epitopes.deduplicate([
            self._row("mhc_i", "ACDEFGHIK", "Zaire", "Spike_T_cell_MHCI_1.1", "HLA-A*01:01"),
            self._row("mhc_i", "ACDEFGHIK", "Zaire", "Spike_T_cell_MHCI_1.2", "HLA-A*02:01"),
        ])
        assert result[0]["alleles"] == "HLA-A*01:01|HLA-A*02:01"


# ---------------------------------------------------------------------------
# stage 5 selection rule
# ---------------------------------------------------------------------------

CLUSTER_CSV = """Cluster.Sub-Cluster Number,Peptide Number,Alignment,Position,Description,Peptide
1.1,Consensus,ACDEFGHIKLMNPQRST,-,-,-
1.1,1,ACDEFGHIK--------,1,Spike_B_cell_BepiPred_1.1,ACDEFGHIK
1.1,2,--DEFGHIKLMN-----,3,Spike_B_cell_EpiDope_1.1,DEFGHIKLMN
1.1,3,----FGHIKLMNP----,5,Spike_T_cell_MHCI_1.1,FGHIKLMNP
1.1,4,---EFGHIKLMNPQRST,4,Spike_T_cell_MHCII_1.1,EFGHIKLMNPQRST
2.1,Consensus,WWWWWWWWWWWW,-,-,-
2.1,1,WWWWWWWWW---,1,Spike_B_cell_BepiPred_2.1,WWWWWWWWW
2.1,2,---WWWWWWWWW,4,Spike_T_cell_MHCI_2.1,WWWWWWWWW
3.1,Consensus,YYYYXYYYYYYY,-,-,-
3.1,1,YYYYXYYYY---,1,Spike_B_cell_BepiPred_3.1,YYYYXYYYY
3.1,2,--YYXYYYYYYY,3,Spike_B_cell_EpiDope_3.1,YYXYYYYYYY
3.1,3,---YXYYYYYY-,4,Spike_T_cell_MHCI_3.1,YXYYYYYY
3.1,4,----XYYYYYYY,5,Spike_T_cell_MHCII_3.1,XYYYYYYY
"""

EPITOPE_ROWS = [
    {"epitope_id": "Spike_B_cell_BepiPred_1.1", "source_tool": "bepipred",
     "epitope_class": "B", "sequence": "ACDEFGHIK"},
    {"epitope_id": "Spike_B_cell_EpiDope_1.1", "source_tool": "epidope",
     "epitope_class": "B", "sequence": "DEFGHIKLMN"},
    {"epitope_id": "Spike_T_cell_MHCI_1.1", "source_tool": "mhc_i",
     "epitope_class": "MHC-I", "sequence": "FGHIKLMNP"},
    {"epitope_id": "Spike_T_cell_MHCII_1.1", "source_tool": "mhc_ii",
     "epitope_class": "MHC-II", "sequence": "EFGHIKLMNPQRST"},
    {"epitope_id": "Spike_B_cell_BepiPred_2.1", "source_tool": "bepipred",
     "epitope_class": "B", "sequence": "WWWWWWWWW"},
    {"epitope_id": "Spike_T_cell_MHCI_2.1", "source_tool": "mhc_i",
     "epitope_class": "MHC-I", "sequence": "WWWWWWWWW"},
    {"epitope_id": "Spike_B_cell_BepiPred_3.1", "source_tool": "bepipred",
     "epitope_class": "B", "sequence": "YYYYXYYYY"},
    {"epitope_id": "Spike_B_cell_EpiDope_3.1", "source_tool": "epidope",
     "epitope_class": "B", "sequence": "YYXYYYYYYY"},
    {"epitope_id": "Spike_T_cell_MHCI_3.1", "source_tool": "mhc_i",
     "epitope_class": "MHC-I", "sequence": "YXYYYYYY"},
    {"epitope_id": "Spike_T_cell_MHCII_3.1", "source_tool": "mhc_ii",
     "epitope_class": "MHC-II", "sequence": "XYYYYYYY"},
]

REQUIRE = {"min_mhc_i": 1, "min_mhc_ii": 1, "min_bepipred": 1, "min_epidope": 1}


@pytest.fixture
def parsed(tmp_path):
    path = tmp_path / "cluster.csv"
    path.write_text(CLUSTER_CSV, encoding="utf-8")
    clusters = cluster_mod.parse(str(path))
    cluster_mod.annotate(clusters, EPITOPE_ROWS)
    return {c.cluster_id: c for c in clusters}


class TestClusterSelection:
    def test_parses_all_clusters(self, parsed):
        assert set(parsed) == {"1.1", "2.1", "3.1"}
        assert parsed["1.1"].consensus == "ACDEFGHIKLMNPQRST"
        assert len(parsed["1.1"].members) == 4

    def test_complete_cluster_passes(self, parsed):
        passes, reason, counts = cluster_mod.evaluate(parsed["1.1"], REQUIRE)
        assert passes, reason
        assert counts == {"bepipred": 1, "epidope": 1, "mhc_i": 1, "mhc_ii": 1}

    def test_missing_epidope_and_mhc_ii_fails(self, parsed):
        passes, reason, _ = cluster_mod.evaluate(parsed["2.1"], REQUIRE)
        assert not passes
        assert "MHC-II" in reason and "EpiDope" in reason

    def test_ambiguous_consensus_is_rejected(self, parsed):
        # Cluster 3.1 is complete but its consensus contains X.
        passes, reason, _ = cluster_mod.evaluate(parsed["3.1"], REQUIRE, reject_ambiguous=True)
        assert not passes
        assert "ambiguous" in reason

    def test_ambiguity_check_can_be_disabled(self, parsed):
        passes, _, _ = cluster_mod.evaluate(parsed["3.1"], REQUIRE, reject_ambiguous=False)
        assert passes

    def test_relaxed_bcell_rule_accepts_single_tool(self, parsed):
        # With require_both False, 2.1 still fails on MHC-II but its B-cell
        # requirement is now satisfied by BepiPred alone.
        _, reason, _ = cluster_mod.evaluate(
            parsed["2.1"], REQUIRE, require_both_bcell_tools=False
        )
        assert "EpiDope" not in reason
        assert "MHC-II" in reason

    def test_unresolved_member_fails_loudly(self, tmp_path):
        # A member that cannot be traced to a predictor must not be silently
        # ignored - that would let a cluster pass on evidence that was never checked.
        path = tmp_path / "c.csv"
        path.write_text(CLUSTER_CSV, encoding="utf-8")
        clusters = cluster_mod.parse(str(path))
        cluster_mod.annotate(clusters, [])  # nothing resolves
        by_id = {c.cluster_id: c for c in clusters}
        # Identifier parsing is the fallback, so these still resolve; blank the
        # descriptions to simulate a truly untraceable member.
        for member in by_id["1.1"].members:
            member.description = "???"
            member.sequence = "???"
            member.source_tool = ""
        passes, reason, _ = cluster_mod.evaluate(by_id["1.1"], REQUIRE)
        assert not passes
        assert "could not be traced" in reason

    def test_identifier_parsing_is_the_fallback(self, tmp_path):
        path = tmp_path / "c.csv"
        path.write_text(CLUSTER_CSV, encoding="utf-8")
        clusters = cluster_mod.parse(str(path))
        cluster_mod.annotate(clusters, [])  # no join table at all
        by_id = {c.cluster_id: c for c in clusters}
        passes, reason, _ = cluster_mod.evaluate(by_id["1.1"], REQUIRE)
        assert passes, reason


# ---------------------------------------------------------------------------
# conservancy
# ---------------------------------------------------------------------------

class TestConservancy:
    def test_exact_match_is_full_identity(self):
        assert conservancy.identity("ACDEF", "XXACDEFXX") == 1.0

    def test_single_mismatch(self):
        assert conservancy.identity("ACDEF", "XXACDEGXX") == pytest.approx(4 / 5)

    def test_best_window_is_taken(self):
        # The second occurrence is exact and must win over the first.
        assert conservancy.identity("ACDEF", "ACDEGXXXACDEF") == 1.0

    def test_target_shorter_than_epitope(self):
        assert conservancy.identity("ACDEFGH", "ACD") == pytest.approx(3 / 7)

    def test_conservancy_fraction(self):
        rows = conservancy.analyse(
            [("ep1", "ACDEF")],
            [("Zaire", "XXACDEFXX"), ("Sudan", "XXACDEFXX"), ("Bundibugyo", "XXXXXXX")],
            threshold=1.0,
        )
        assert rows[0]["n_matched"] == 2
        assert rows[0]["conservancy_percent"] == pytest.approx(66.67, abs=0.01)
        assert rows[0]["matched_strains"] == "Sudan|Zaire"

    def test_deduplicate_collapses_identical_sequences(self):
        records = [
            fasta.Record(header="a", sequence="ACDEF"),
            fasta.Record(header="b", sequence="ACDEF"),
            fasta.Record(header="c", sequence="GHIKL"),
        ]
        assert len(conservancy.deduplicate_records(records)) == 2


# ---------------------------------------------------------------------------
# FASTA handling
# ---------------------------------------------------------------------------

class TestFasta:
    def test_ncbi_defline_is_parsed(self, tmp_path):
        path = tmp_path / "NP_066246.1 spike glycoprotein - Zaire ebolavirus.fasta"
        path.write_text(">NP_066246.1 spike glycoprotein [Zaire ebolavirus]\nMGVTG\nILQLP\n")
        record, = fasta.read_fasta(str(path))
        assert record.accession == "NP_066246.1"
        assert record.name == "spike glycoprotein"
        assert record.strain == "Zaire ebolavirus"
        assert record.sequence == "MGVTGILQLP"

    def test_filename_is_the_fallback(self, tmp_path):
        path = tmp_path / "AAB81001.1 nucleoprotein - Zaire ebolavirus.fasta"
        path.write_text(">AAB81001.1\nMDSRPQ\n")
        record, = fasta.read_fasta(str(path))
        assert record.name == "nucleoprotein"
        assert record.strain == "Zaire ebolavirus"

    def test_gaps_and_stops_are_stripped(self, tmp_path):
        path = tmp_path / "x.fasta"
        path.write_text(">x test [strain]\nAC-DE*F\n")
        record, = fasta.read_fasta(str(path))
        assert record.sequence == "ACDEF"

    def test_nucleotide_input_is_caught(self):
        record = fasta.Record(header="x", sequence="ACGTACGTACGT")
        assert any("nucleotide" in problem for problem in fasta.validate(record))

    def test_short_sequence_is_caught(self):
        record = fasta.Record(header="x", sequence="ACDEF")
        assert any("below the minimum" in problem for problem in fasta.validate(record, 50))

    @pytest.mark.parametrize("name,expected", [
        ("spike glycoprotein", "Spike"),
        ("Spike Glycoprotein", "Spike"),
        ("spike glycoprotein precursor", "Spike"),
        ("nucleoprotein", "Nucleoprotein"),
        ("some novel protein", "SomeNovelProtein"),
    ])
    def test_protein_label_normalisation(self, name, expected):
        aliases = {"spike glycoprotein": "Spike", "nucleoprotein": "Nucleoprotein"}
        assert fasta.normalise_protein_label(name, aliases) == expected

    def test_longest_alias_wins(self):
        # "glycoprotein" also matches, but the more specific alias must win.
        aliases = {"glycoprotein": "GP", "spike glycoprotein": "Spike"}
        assert fasta.normalise_protein_label("spike glycoprotein", aliases) == "Spike"


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------

class TestTables:
    def test_round_trip(self, tmp_path):
        path = tmp_path / "t.tsv"
        rows = [{"a": 1, "b": "x", "c": True}, {"a": 2, "b": "y", "c": False}]
        tables.write_table(str(path), ["a", "b", "c"], rows)
        back = tables.read_table(str(path))
        assert [row["a"] for row in back] == ["1", "2"]
        assert [tables.as_bool(row["c"]) for row in back] == [True, False]

    def test_missing_keys_become_empty(self, tmp_path):
        path = tmp_path / "t.tsv"
        tables.write_table(str(path), ["a", "b"], [{"a": 1}])
        assert tables.read_table(str(path))[0]["b"] == ""

    def test_floats_are_stable(self, tmp_path):
        path = tmp_path / "t.tsv"
        tables.write_table(str(path), ["v"], [{"v": 0.5000000}, {"v": 1.0}])
        assert [row["v"] for row in tables.read_table(str(path))] == ["0.5", "1"]


# ---------------------------------------------------------------------------
# screening helpers
# ---------------------------------------------------------------------------

class TestScreening:
    def test_sliding_windows(self):
        # Moved to mpvpredpip.autoimmunity along with the rest of the human-proteome
        # screen; the behaviour is unchanged.
        assert autoimmunity.windows("ACDEFGH", 9) == ["ACDEFGH"]
        assert autoimmunity.windows("ACDEFGHIKL", 9) == ["ACDEFGHIK", "CDEFGHIKL"]

    def test_antigenicity_is_deterministic(self):
        assert screening.antigenicity_kolaskar("ACDEFGHIK") == screening.antigenicity_kolaskar("ACDEFGHIK")

    def test_antigenicity_ignores_unknown_residues(self):
        assert screening.antigenicity_kolaskar("XXXX") == 0.0

    def test_antigenicity_responds_to_composition(self):
        # Cysteine and valine carry the highest propensity values; a run of them
        # must score above a run of the lowest (asparagine).
        assert screening.antigenicity_kolaskar("CVCVCVCVC") > screening.antigenicity_kolaskar("NNNNNNNNN")


# ---------------------------------------------------------------------------
# autoimmunity - the in-process human proteome screen
# ---------------------------------------------------------------------------
# This replaces the PIR Peptide Match standalone, so it is the only thing
# standing between a self-matching peptide and the final candidate list. The
# cases that matter are the ones where a bug looks like a clean result: a window
# spanning two proteins, an offset attributed to the wrong protein, or an
# ambiguous window quietly counted as "no match".

class TestAutoimmunity:
    def _proteome(self, tmp_path, records, name="proteome.fasta"):
        path = tmp_path / name
        fasta.write_fasta(str(path), records, width=0)
        return autoimmunity.HumanProteome(str(path))

    def test_finds_an_exact_window(self, tmp_path):
        proteome = self._proteome(tmp_path, [("HUMAN1", "MKWVTFISLLLLFSSAYS")])
        results = autoimmunity.find_matches({"c1": "FISLLLLFS"}, proteome, 9)
        assert results["c1"].flagged
        assert results["c1"].matches[0].protein_id == "HUMAN1"

    def test_reports_a_clean_peptide_rather_than_omitting_it(self, tmp_path):
        """Every submitted peptide gets a result, matched or not.

        An absent key would be indistinguishable from a peptide that was never
        screened - the same conflation that made the MHC-I failure invisible.
        """
        proteome = self._proteome(tmp_path, [("HUMAN1", "MKWVTFISLLLLFSSAYS")])
        results = autoimmunity.find_matches({"c1": "WWWWWWWWW"}, proteome, 9)
        assert "c1" in results
        assert not results["c1"].flagged

    def test_window_cannot_span_two_proteins(self, tmp_path):
        """The separator is why this is safe.

        "AAAAA" ends protein 1 and "CCCCC" starts protein 2. Their concatenation
        exists in neither, and must not be reported as a human match.
        """
        proteome = self._proteome(tmp_path, [("P1", "MKWVTAAAAA"), ("P2", "CCCCCMKWVT")])
        results = autoimmunity.find_matches({"c1": "AAAAACCCCC"}, proteome, 10)
        assert not results["c1"].flagged

    def test_offset_maps_to_the_right_protein(self, tmp_path):
        """A hit inside the second protein must not be attributed to the first."""
        proteome = self._proteome(
            tmp_path, [("FIRST", "MKWVTFISLL"), ("SECOND", "QQQQQYNDYLE")]
        )
        results = autoimmunity.find_matches({"c1": "QYNDYLE"}, proteome, 7)
        match = results["c1"].matches[0]
        assert match.protein_id == "SECOND"
        assert match.offset == 4                      # 0-based inside SECOND
        assert match.describe().endswith(":5")        # 1-based when reported

    def test_finds_every_occurrence(self, tmp_path):
        proteome = self._proteome(tmp_path, [("P1", "ACDEFACDEF"), ("P2", "WWACDEFWW")])
        results = autoimmunity.find_matches({"c1": "ACDEF"}, proteome, 5)
        assert {m.protein_id for m in results["c1"].matches} == {"P1", "P2"}
        assert len(results["c1"].matches) == 3

    def test_ambiguous_window_is_skipped_not_cleared(self, tmp_path):
        """An X window is unscreenable, which is not the same as unmatched.

        Stage 5 can emit X where a cluster's members disagree. Searching for it
        literally would assert that an unresolved residue corresponds to an
        unsequenced one, so the window is set aside - and counted, so the
        partial coverage is visible in the output.
        """
        proteome = self._proteome(tmp_path, [("HUMAN1", "MKWVTFISLLLLFSSAYS")])
        results = autoimmunity.find_matches({"c1": "FISXLLLFS"}, proteome, 9)
        assert not results["c1"].flagged
        assert results["c1"].skipped_windows == ["FISXLLLFS"]

    def test_peptide_shorter_than_the_window_is_still_screened(self, tmp_path):
        """Short peptides must not fall through the screen untested."""
        proteome = self._proteome(tmp_path, [("HUMAN1", "MKWVTFISLLLLFSSAYS")])
        results = autoimmunity.find_matches({"c1": "MKWVT"}, proteome, 9)
        assert results["c1"].flagged

    def test_matching_is_case_insensitive(self, tmp_path):
        proteome = self._proteome(tmp_path, [("HUMAN1", "mkwvtfisllllfssays")])
        results = autoimmunity.find_matches({"c1": "fisllllfs"}, proteome, 9)
        assert results["c1"].flagged

    def test_empty_proteome_is_an_error(self, tmp_path):
        path = tmp_path / "empty.fasta"
        path.write_text("", encoding="utf-8")
        with pytest.raises(ValueError, match="No sequences"):
            autoimmunity.HumanProteome(str(path))


# ---------------------------------------------------------------------------
# AlgPred 2.0 - the parts of its installation that silently do not work
# ---------------------------------------------------------------------------
# The tool cannot be driven in tests, but three things about it fail in ways
# that do not name their own cause, and all three are decided before it runs:
# an LFS-stub model, an envfile pointing at the author's machine, and a parser
# that takes exactly four positionally-ordered lines.

class TestAlgPredSetup:
    def _install(self, tmp_path, model=b"\x80\x04fake joblib pickle"):
        root = tmp_path / "algpred2"
        (root / "progs").mkdir(parents=True)
        (root / "Database").mkdir(parents=True)
        (root / "progs" / "MERCI_motif_locator.pl").write_text("#!perl\n", encoding="utf-8")
        (root / "Database" / "pos_ige_motifs.txt").write_text("motif\n", encoding="utf-8")
        (root / "Database" / "data.phr").write_bytes(b"\x00")
        if model is not None:
            (root / "rf_model").write_bytes(model)
        return str(root)

    def test_missing_model_is_named_as_lfs(self, tmp_path):
        root = self._install(tmp_path, model=None)
        with pytest.raises(algpred.AlgPredError, match="LFS"):
            algpred.check_model(root)

    def test_lfs_pointer_is_not_mistaken_for_the_model(self, tmp_path):
        """The stub is a valid file of the right name - and useless.

        Without this check joblib raises an unpickling error that says nothing
        about why, which is a long way from the real cause.
        """
        pointer = b"version https://git-lfs.github.com/spec/v1\noid sha256:abc\nsize 1\n"
        root = self._install(tmp_path, model=pointer)
        with pytest.raises(algpred.AlgPredError, match="pointer"):
            algpred.check_model(root)

    def test_real_model_passes(self, tmp_path):
        algpred.check_model(self._install(tmp_path))

    def test_envfile_matches_the_parser_contract(self, tmp_path, monkeypatch):
        """algpred2.py takes exactly 4 non-comment lines, split on first colon.

        It reads them positionally - BLAST, database, MERCI, motifs - so the
        order is part of the format, not a convention.
        """
        monkeypatch.setattr(algpred.shutil, "which", lambda name: "/usr/bin/blastp")
        root = self._install(tmp_path)
        workdir = tmp_path / "work"
        workdir.mkdir()
        algpred.write_envfile(str(workdir), root)

        lines = (workdir / "envfile").read_text(encoding="utf-8").splitlines()
        payload = [line for line in lines if "#" not in line]
        assert len(payload) == 4                       # the parser rejects any other count

        paths = [line.split(":")[1] for line in payload]   # exactly how AlgPred parses it
        assert paths[0] == "/usr/bin/blastp"
        assert paths[1].endswith(os.path.join("Database", "data"))
        assert paths[2].endswith(os.path.join("progs", "MERCI_motif_locator.pl"))
        assert paths[3].endswith(os.path.join("Database", "pos_ige_motifs.txt"))

    def test_envfile_replaces_a_symlink_to_the_shipped_one(self, tmp_path, monkeypatch):
        """The run directory is built from symlinks; writing through one would
        corrupt the installed copy."""
        monkeypatch.setattr(algpred.shutil, "which", lambda name: "/usr/bin/blastp")
        root = self._install(tmp_path)
        shipped = os.path.join(root, "envfile")
        with open(shipped, "w", encoding="utf-8") as handle:
            handle.write("BLAST:/home/neelam/blastp\n")

        workdir = tmp_path / "work"
        workdir.mkdir()
        algpred.link_or_copy(shipped, str(workdir / "envfile"))
        algpred.write_envfile(str(workdir), root)

        assert "neelam" not in (workdir / "envfile").read_text(encoding="utf-8")
        assert "neelam" in open(shipped, encoding="utf-8").read()

    def test_missing_blastp_is_reported_as_such(self, tmp_path, monkeypatch):
        monkeypatch.setattr(algpred.shutil, "which", lambda name: None)
        root = self._install(tmp_path)
        workdir = tmp_path / "work"
        workdir.mkdir()
        with pytest.raises(algpred.AlgPredError, match="blastp"):
            algpred.write_envfile(str(workdir), root)

    def test_sklearn_externals_import_is_patched(self, tmp_path):
        """AlgPred will not even start under modern scikit-learn without this.

        ``sklearn.externals.joblib`` was removed in 0.23; the environment pins
        1.5.2 for IApred's sake, so the import has to be rewritten.
        """
        source = tmp_path / "algpred2.py"
        source.write_text(
            "import numpy as np\nfrom sklearn.externals import joblib\n"
            "clf = joblib.load('rf_model')\n"
            "data_test = np.loadtxt(file_name, delimiter=',')\n",
            encoding="utf-8",
        )
        destination = tmp_path / "work.py"
        applied = algpred.stage_script(str(source), str(destination))

        patched = destination.read_text(encoding="utf-8")
        assert "sklearn.externals" not in patched
        assert "import joblib" in patched
        assert "joblib.load('rf_model')" in patched     # the call site still works
        assert "joblib-import" in applied

    def test_single_sequence_bug_is_patched(self, tmp_path):
        """AlgPred cannot score one peptide without this.

        ``np.loadtxt`` yields a 1-D array for a single-row file, so
        ``predict_proba`` receives ``(20,)`` rather than ``(1, 20)`` and
        raises. Stage 5 regularly produces exactly one consensus peptide for a
        protein - it did for the nucleoprotein - so this is a normal input,
        not an edge case.
        """
        source = tmp_path / "algpred2.py"
        source.write_text(
            "import numpy as np\nimport joblib\n"
            "data_test = np.loadtxt(file_name, delimiter=',')\n",
            encoding="utf-8",
        )
        destination = tmp_path / "work.py"
        applied = algpred.stage_script(str(source), str(destination))

        patched = destination.read_text(encoding="utf-8")
        assert "np.atleast_2d(np.loadtxt(file_name, delimiter=','))" in patched
        assert "single-sequence" in applied

    def test_reworded_loadtxt_is_caught(self, tmp_path):
        """The failure that actually matters: the patch stops matching.

        If upstream renames the variable or changes the spacing, the
        replacement silently does nothing and AlgPred goes back to being
        unable to score a single sequence. The check is written as an outcome
        - the fixed form must be present - precisely so that this is caught
        here rather than as a ValueError from inside scikit-learn.
        """
        source = tmp_path / "algpred2.py"
        source.write_text(
            "import numpy as np\nimport joblib\n"
            "data_test = np.loadtxt(infile, delimiter=',')\n",   # renamed arg
            encoding="utf-8",
        )
        with pytest.raises(algpred.AlgPredError, match="two dimensions"):
            algpred.stage_script(str(source), str(tmp_path / "out.py"))

    def test_already_wrapped_source_is_accepted(self, tmp_path):
        """An upstream fix in the same shape must not be reported as a problem.

        ``atleast_2d`` is idempotent, so double-wrapping would be harmless
        anyway - but the check has to pass, or a fixed upstream would look
        broken.
        """
        source = tmp_path / "algpred2.py"
        source.write_text(
            "import numpy as np\nimport joblib\n"
            "data_test = np.atleast_2d(np.loadtxt(file_name, delimiter=','))\n",
            encoding="utf-8",
        )
        algpred.stage_script(str(source), str(tmp_path / "out.py"))

    def test_installed_script_is_never_modified(self, tmp_path):
        """The run directory is symlinks; patching must not write through one."""
        source = tmp_path / "algpred2.py"
        original = (
            "from sklearn.externals import joblib\n"
            "data_test = np.loadtxt(file_name, delimiter=',')\n"
        )
        source.write_text(original, encoding="utf-8")

        workdir = tmp_path / "work"
        workdir.mkdir()
        destination = workdir / "algpred2.py"
        algpred.link_or_copy(str(source), str(destination))
        algpred.stage_script(str(source), str(destination))

        assert source.read_text(encoding="utf-8") == original
        assert "sklearn.externals" not in destination.read_text(encoding="utf-8")

    def test_unknown_sklearn_externals_use_is_not_silently_passed(self, tmp_path):
        """A form the patch list does not cover must fail loudly, not quietly."""
        source = tmp_path / "algpred2.py"
        source.write_text("import sklearn.externals.joblib as jl\n", encoding="utf-8")
        with pytest.raises(algpred.AlgPredError, match="sklearn.externals"):
            algpred.stage_script(str(source), str(tmp_path / "out.py"))


# ---------------------------------------------------------------------------
# IAPred output parsing
# ---------------------------------------------------------------------------
# The adapter's *invocation* of IAPred cannot be tested without the tool
# installed. Its *parsing* can be, and that is the half that would silently
# corrupt results: a parser reading the wrong column produces plausible numbers.
#
# The first test below uses IApred's REAL output header, taken from its README.
# The rest cover layouts it might plausibly emit if the format changes.

class TestIAPredParsing:
    def _write(self, tmp_path, text, name="iapred.out"):
        path = tmp_path / name
        path.write_text(text, encoding="utf-8")
        return str(path)

    def test_real_iapred_csv_format(self, tmp_path):
        """IApred's actual output header, captured from a real run.

        This is NOT what its README documents - the README says ``IAscore``, the
        code writes ``Intrinsic_Antigenicity_Score``. The score must come from
        that column and NOT from ``Antigenicity_Category``, which holds a text
        label; an earlier candidate list matched the category column by prefix
        and made every score unparseable.

        The negative value is deliberate: IApred scores run roughly -3 to 3, so
        anything assuming a 0-1 probability is wrong about this tool.
        """
        path = self._write(
            tmp_path,
            "Header,Sequence_Length,Intrinsic_Antigenicity_Score,Antigenicity_Category\n"
            "Spike|Zaire|NP_1,676,0.69,High\n"
            "NP|Zaire|AAB_1,739,-0.42,Low\n",
        )
        assert screening._parse_iapred(path) == {
            "Spike|Zaire|NP_1": pytest.approx(0.69),
            "NP|Zaire|AAB_1": pytest.approx(-0.42),
        }

    def test_stale_readme_column_name_still_works(self, tmp_path):
        """``IAscore`` stays a fallback, in case the code is realigned to the
        README or an older IApred build is installed."""
        path = self._write(tmp_path, "Header,IAscore\nA,1.25\n")
        assert screening._parse_iapred(path) == {"A": pytest.approx(1.25)}

    def test_category_column_alone_is_not_used_as_score(self, tmp_path):
        """With no numeric column present, fail rather than parse the label."""
        path = self._write(tmp_path, "Header,Antigenicity_Category\nA,High\n")
        with pytest.raises(screening.base.ToolError):
            screening._parse_iapred(path)

    def test_sentinels_become_none_not_an_error(self, tmp_path):
        """IApred writes text into the score column for sequences it won't score.

        A peptide under 20 aa gets "Sequence too short". One of those among fifty
        consensus peptides must not abort the stage - it becomes None, which the
        caller reports as its own screening failure.
        """
        path = self._write(
            tmp_path,
            "Header,Sequence_Length,Intrinsic_Antigenicity_Score,Antigenicity_Category\n"
            "Good_1,676,0.69,High\n"
            "Short_1,14,Sequence too short,N/A\n"
            "Empty_1,0,Invalid sequence,N/A\n"
            "Broke_1,120,Error,N/A\n",
        )
        assert screening._parse_iapred(path) == {
            "Good_1": pytest.approx(0.69),
            "Short_1": None,
            "Empty_1": None,
            "Broke_1": None,
        }

    def test_a_real_non_numeric_value_still_fails(self, tmp_path):
        """Only the known sentinels are tolerated; anything else is a bug."""
        path = self._write(
            tmp_path,
            "Header,Intrinsic_Antigenicity_Score\nA,not-a-number\n",
        )
        with pytest.raises(screening.base.ToolError):
            screening._parse_iapred(path)

    def test_tab_separated(self, tmp_path):
        path = self._write(tmp_path, "id\tscore\nSpike|Zaire|NP_1\t0.71\nNP|Zaire|AAB_1\t0.32\n")
        assert screening._parse_iapred(path) == {
            "Spike|Zaire|NP_1": pytest.approx(0.71),
            "NP|Zaire|AAB_1": pytest.approx(0.32),
        }

    def test_comma_separated(self, tmp_path):
        path = self._write(tmp_path, "sequence_id,iapred_score\nA,0.9\nB,0.1\n")
        assert screening._parse_iapred(path) == {"A": pytest.approx(0.9), "B": pytest.approx(0.1)}

    def test_alternative_column_names(self, tmp_path):
        path = self._write(tmp_path, "Protein\tAntigenicity\nX\t0.55\n")
        assert screening._parse_iapred(path) == {"X": pytest.approx(0.55)}

    def test_extra_columns_are_ignored(self, tmp_path):
        path = self._write(tmp_path, "id\tlength\tscore\tcall\nA\t120\t0.8\tantigenic\n")
        assert screening._parse_iapred(path) == {"A": pytest.approx(0.8)}

    def test_leading_gt_is_stripped(self, tmp_path):
        path = self._write(tmp_path, "id\tscore\n>A desc here\t0.4\n")
        assert screening._parse_iapred(path) == {"A": pytest.approx(0.4)}

    def test_blank_rows_are_skipped(self, tmp_path):
        path = self._write(tmp_path, "id\tscore\nA\t0.4\n\t\n")
        assert screening._parse_iapred(path) == {"A": pytest.approx(0.4)}

    def test_unrecognised_columns_fail_loudly(self, tmp_path):
        # Must raise, not guess. A wrong column silently produces valid-looking
        # scores, which is the worst outcome for this pipeline.
        path = self._write(tmp_path, "alpha\tbeta\n1\t2\n")
        with pytest.raises(screening.base.ToolError) as exc:
            screening._parse_iapred(path)
        assert "alpha" in str(exc.value)  # the real header is reported back

    def test_empty_file_fails(self, tmp_path):
        path = self._write(tmp_path, "id\tscore\n")
        with pytest.raises(screening.base.ToolError):
            screening._parse_iapred(path)


class TestReadDelimitedHeaderSkipping:
    """Tools that print settings before their table must still parse.

    The IEDB immunogenicity predictor prints "masking: default", "masked
    variables: ..." and a blank line before its (comma-separated) table. Reading
    line 1 as the header gave a single column literally named
    "masking: default", and every column lookup failed.
    """

    def _write(self, tmp_path, text, name="raw.tsv"):
        path = tmp_path / name
        path.write_text(text, encoding="utf-8")
        return str(path)

    REAL_OUTPUT = (
        "masking: default\n"
        "masked variables: [1, 2, 9]\n"
        "\n"
        "peptide,length,score\n"
        "AAAAAAAAA,9,0.12345\n"
        "CCCCCCCCC,9,-0.54321\n"
    )

    def test_preamble_is_skipped(self, tmp_path):
        path = self._write(tmp_path, self.REAL_OUTPUT)
        rows = bcell.base.read_delimited(path, header_contains="peptide")
        assert [row["peptide"] for row in rows] == ["AAAAAAAAA", "CCCCCCCCC"]
        assert rows[0]["score"] == "0.12345"

    def test_negative_scores_survive(self, tmp_path):
        # Immunogenicity scores are routinely negative; min_score is 0.0, so a
        # sign lost here would silently pass every peptide.
        path = self._write(tmp_path, self.REAL_OUTPUT)
        rows = bcell.base.read_delimited(path, header_contains="peptide")
        assert float(rows[1]["score"]) < 0

    def test_without_the_hint_the_preamble_wins(self, tmp_path):
        # Documents why the parameter exists: this is the old behaviour, and it
        # is exactly how stage 4 failed. It must not raise either - a line with
        # more fields than the header puts a list in the row, which used to
        # crash the whitespace filter with an AttributeError.
        path = self._write(tmp_path, self.REAL_OUTPUT)
        rows = bcell.base.read_delimited(path)
        assert all("peptide" not in row for row in rows)

    def test_a_table_with_no_preamble_still_works(self, tmp_path):
        path = self._write(tmp_path, "peptide,length,score\nAAAAAAAAA,9,0.1\n")
        rows = bcell.base.read_delimited(path, header_contains="peptide")
        assert rows[0]["peptide"] == "AAAAAAAAA"

    def test_a_missing_header_is_an_error_not_an_empty_result(self, tmp_path):
        # The important half: a tool that printed a message instead of a table
        # must not read as "nothing qualified".
        path = self._write(tmp_path, "some error the tool printed\nand another line\n")
        with pytest.raises(bcell.base.ToolError) as exc:
            bcell.base.read_delimited(path, header_contains="peptide")
        assert "some error the tool printed" in str(exc.value)


class TestBindingOutputValidation:
    """A predictor that fails must not look like a predictor that found nothing.

    The IEDB binding tools print startup failures to stdout and exit 0, and this
    adapter sends stdout to the results file. Without this check, a broken
    install produced a green stage 4 run with zero class I binders and 4231
    class II epitopes merged with no class I support whatsoever.
    """

    def _write(self, tmp_path, text):
        path = tmp_path / "raw_binding.tsv"
        path.write_text(text, encoding="utf-8")
        return str(path)

    HEADER = "allele\tseq_num\tstart\tend\tlength\tpeptide\tcore\ticore\tScore\trank"

    def test_a_valid_header_with_no_rows_is_accepted(self):
        # A real "nothing passed the cut-off". Must NOT raise - that is a
        # biological result, and the caller reports it as a warning.
        import tempfile

        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "raw.tsv")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(self.HEADER + "\n")
            iedb._check_binding_output(path, "I", ["python", "predict_binding.py"])

    def test_a_valid_table_is_accepted(self, tmp_path):
        path = self._write(
            tmp_path,
            self.HEADER + "\nHLA-A*02:01\t1\t1\t9\t9\tMGVTGILQL\tMGVTGILQL\tMGVTGILQL\t0.5\t0.8\n",
        )
        iedb._check_binding_output(path, "I", ["python", "predict_binding.py"])

    def test_an_error_message_on_stdout_is_rejected(self, tmp_path):
        # The exact failure seen: python 3.8 could not satisfy the bundled
        # netmhcpan-4.1-executable, and this one line was the whole file.
        path = self._write(
            tmp_path,
            "cannot import name 'files' from 'importlib.resources' "
            "(/opt/conda/lib/python3.8/importlib/resources.py)\n",
        )
        with pytest.raises(iedb.base.ToolError) as exc:
            iedb._check_binding_output(path, "I", ["python", "predict_binding.py"])
        # The tool's own message must survive into the error, or the next person
        # debugging this learns nothing from it.
        assert "importlib.resources" in str(exc.value)

    def test_an_empty_file_is_rejected(self, tmp_path):
        path = self._write(tmp_path, "")
        with pytest.raises(iedb.base.ToolError):
            iedb._check_binding_output(path, "II", ["python", "mhc_II_binding.py"])

    def test_the_class_ii_header_is_accepted(self, tmp_path):
        # Different columns from class I - core_peptide, ic50, percentile_rank -
        # so the check must not be written against class I's shape alone.
        path = self._write(
            tmp_path,
            "allele\tseq_num\tstart\tend\tlength\tcore_peptide\tpeptide\tic50\t"
            "percentile_rank\tadjusted_rank\n",
        )
        iedb._check_binding_output(path, "II", ["python", "mhc_II_binding.py"])


class TestEntryPointsAreAbsolute:
    """Every adapter's entry-point finder must return an absolute path.

    This is a regression test for a bug made three separate times - IAPred,
    BepiPred and the IEDB tools. All of them run with ``cwd=tool_dir`` so the
    tool can find its own bundled data, so a project-relative path gets resolved
    against that directory twice and the process dies with "No such file or
    directory" for a path that is plainly there. It costs a full pipeline run to
    notice, so it is pinned here instead.
    """

    def _tool(self, tmp_path, *relative):
        for part in relative:
            path = tmp_path / part
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("# entry point\n", encoding="utf-8")
        return str(tmp_path)

    def test_iedb_entry_direct_hit(self, tmp_path):
        tool_dir = self._tool(tmp_path, "predict_binding.py")
        found = iedb._entry(tool_dir, ["predict_binding.py"], "test tool")
        assert os.path.isabs(found)
        assert found.endswith("predict_binding.py")

    def test_iedb_entry_found_by_glob(self, tmp_path):
        # The nested case is the one that actually bit: the name resolves through
        # the recursive glob rather than the direct join.
        tool_dir = self._tool(tmp_path, "src/predict_binding.py")
        found = iedb._entry(tool_dir, ["src/predict_binding.py"], "test tool")
        assert os.path.isabs(found)

    def test_iedb_entry_absolute_even_for_a_relative_tool_dir(self, tmp_path, monkeypatch):
        # config.yaml gives tool paths relative to the project root, so this is
        # how the adapters are actually called.
        self._tool(tmp_path, "tools/mhc_i/src/predict_binding.py")
        monkeypatch.chdir(tmp_path)
        found = iedb._entry("tools/mhc_i", ["src/predict_binding.py"], "test tool")
        assert os.path.isabs(found)
        assert os.path.exists(found)

    def test_bcell_entry_is_absolute(self, tmp_path, monkeypatch):
        self._tool(tmp_path, "tools/bepipred3/bepipred3_CLI.py")
        monkeypatch.chdir(tmp_path)
        found = bcell._find_entrypoint(
            "tools/bepipred3", ["bepipred3_CLI.py"], "test tool"
        )
        assert os.path.isabs(found)
        assert os.path.exists(found)

    def test_missing_entry_point_is_actionable(self, tmp_path):
        with pytest.raises(bcell.base.ToolMissing):
            iedb._entry(str(tmp_path), ["nothing_here.py"], "test tool")


class TestPercentIdentity:
    """The IEDB cluster tool's identity rule, reimplemented in mpvpredpip.cluster.

    Transcribed from IEDB_Cluster-1.0's own source. These tests pin the three
    properties that a conventional identity measure would get wrong, so that a
    later "tidy-up" of the formula fails loudly instead of silently reshaping
    every cluster in the run.
    """

    def test_exact_substring_scores_99_9_not_100(self):
        # Never 100: the tool reserves that, so a short peptide inside a long one
        # is not reported as identical to it.
        assert cluster_mod.percent_identity("TKNITDKID", "PHDWTKNITDKIDQII") == 99.9

    def test_identical_peptides_also_score_99_9(self):
        # A consequence of the substring rule, not an oversight - a string is a
        # substring of itself, so the short-circuit fires first.
        assert cluster_mod.percent_identity("PHDWTKNITDKIDQ", "PHDWTKNITDKIDQ") == 99.9

    def test_match_count_is_divided_by_the_longer_peptide(self):
        # 9 of 9 residues match, but the denominator is the longer sequence's 18,
        # so this is 50% and not 100%. Only the substring short-circuit rescues
        # such a pair, and it does not apply here because the match is gapped.
        short, long_ = "AAAAAAAAA", "AAAAAAAAACCCCCCCCC"
        assert cluster_mod.percent_identity(short, long_) == 99.9  # substring wins
        # Same lengths, but now not a substring:
        assert cluster_mod.percent_identity("AAAAACAAAA", "AAAAAGAAAACCCCCCCCC") == pytest.approx(
            100.0 * 9 / 19
        )

    def test_ungapped_best_offset_is_taken(self):
        assert cluster_mod.percent_identity("CCAAAAC", "GGGGAAAAGGG") == pytest.approx(
            100.0 * 4 / 11
        )

    def test_single_mismatch(self):
        assert cluster_mod.percent_identity(
            "PHDWTKNITDKIDQ", "PHDWTKNITDKIDA"
        ) == pytest.approx(100.0 * 13 / 14)

    def test_unrelated_peptides_score_zero(self):
        assert cluster_mod.percent_identity("AAAAAAAAAAAAAA", "CCCCCCCCCCCCCC") == 0.0

    def test_symmetric(self):
        first, second = "PHDWTKNITDKIDQ", "WTKNITDKIDQIIHD"
        assert cluster_mod.percent_identity(first, second) == cluster_mod.percent_identity(
            second, first
        )

    def test_empty_is_not_an_error(self):
        assert cluster_mod.percent_identity("", "ABC") == 0.0


class TestClusterPeptides:
    def test_identical_peptides_collapse_when_unique(self):
        clusters = cluster_mod.cluster_peptides(
            [("a", "PEPTIDEPEPTIDE"), ("b", "PEPTIDEPEPTIDE")], threshold=60
        )
        assert len(clusters) == 1
        assert len(clusters[0].members) == 1

    def test_unique_off_keeps_both(self):
        clusters = cluster_mod.cluster_peptides(
            [("a", "PEPTIDEPEPTIDE"), ("b", "PEPTIDEPEPTIDE")], threshold=60, unique=False
        )
        assert len(clusters) == 1
        assert len(clusters[0].members) == 2

    def test_unrelated_peptides_form_separate_clusters(self):
        clusters = cluster_mod.cluster_peptides(
            [("a", "AAAAAAAAAAAAAA"), ("b", "CCCCCCCCCCCCCC")], threshold=60
        )
        assert len(clusters) == 2
        assert [c.cluster_id for c in clusters] == ["1.1", "2.1"]

    def test_clustering_is_transitive(self):
        # a~b and b~c, but a and c are not similar enough on their own. The tool
        # merges any lists that share a member, so all three land together.
        entries = [
            ("a", "AAAAAAAAAAAACC"),
            ("b", "AAAAAAACCCCCCC"),
            ("c", "AACCCCCCCCCCCC"),
        ]
        clusters = cluster_mod.cluster_peptides(entries, threshold=50)
        assert len(clusters) == 1
        assert len(clusters[0].members) == 3

    def test_threshold_is_inclusive(self):
        # 13/14 = 92.857...; clustering at exactly that value must still join them.
        entries = [("a", "PHDWTKNITDKIDQ"), ("b", "PHDWTKNITDKIDA")]
        assert len(cluster_mod.cluster_peptides(entries, threshold=100.0 * 13 / 14)) == 1

    def test_member_order_follows_input(self):
        entries = [("first", "PHDWTKNITDKIDQ"), ("second", "PHDWTKNITDKIDA")]
        clusters = cluster_mod.cluster_peptides(entries, threshold=60)
        assert [m.description for m in clusters[0].members] == ["first", "second"]


class TestBuildConsensus:
    def test_agreeing_members_give_a_clean_consensus(self):
        consensus, alignments, positions = cluster_mod.build_consensus(
            ["PHDWTKNITDKIDQ", "PHDWTKNITDKIDQ"]
        )
        assert consensus == "PHDWTKNITDKIDQ"
        assert cluster_mod.AMBIGUOUS not in consensus
        assert positions == [1, 1]
        assert alignments == ["PHDWTKNITDKIDQ", "PHDWTKNITDKIDQ"]

    def test_a_two_way_tie_is_ambiguous(self):
        # Two members, two different residues: equally supported, so there is no
        # defensible winner and X is correct.
        consensus, _, _ = cluster_mod.build_consensus(
            ["PHDWTKNITDKIDQ", "PHDWTKNITDKIDA"]
        )
        assert consensus == "PHDWTKNITDKID" + cluster_mod.AMBIGUOUS

    def test_the_majority_residue_wins(self):
        # The change that unblocked stage 5. Under the old unanimity rule this
        # was X, which is not a peptide anyone can synthesise or screen; three
        # Ebola species differing at one position is the normal case here, not an
        # error.
        consensus, _, _ = cluster_mod.build_consensus(
            ["PHDWTKNITDKIDQ", "PHDWTKNITDKIDQ", "PHDWTKNITDKIDA"]
        )
        assert consensus == "PHDWTKNITDKIDQ"
        assert cluster_mod.AMBIGUOUS not in consensus

    def test_a_three_way_split_with_a_clear_winner(self):
        consensus, _, _ = cluster_mod.build_consensus(
            ["AAAAAAAAAAAAAQ", "AAAAAAAAAAAAAQ", "AAAAAAAAAAAAAC", "AAAAAAAAAAAAAD"]
        )
        assert consensus.endswith("Q")

    def test_a_three_way_tie_is_ambiguous(self):
        consensus, _, _ = cluster_mod.build_consensus(
            ["AAAAAAAAAAAAAQ", "AAAAAAAAAAAAAC", "AAAAAAAAAAAAAD"]
        )
        assert consensus.endswith(cluster_mod.AMBIGUOUS)

    def test_a_position_only_one_member_covers_is_not_a_tie(self):
        # The overhang of a longer member is supported by one peptide and
        # contradicted by none, so it is not ambiguous.
        consensus, _, _ = cluster_mod.build_consensus(["AAAAAAAAAAAAAA", "AAAAAAAAAAAAAAWWW"])
        assert cluster_mod.AMBIGUOUS not in consensus
        assert consensus.endswith("WWW")

    def test_overlapping_members_extend_the_consensus(self):
        # The case the IEDB web tool's own example shows: two peptides that agree
        # where they overlap merge into a span longer than either.
        consensus, alignments, positions = cluster_mod.build_consensus(
            ["KKDAPYIVGDVVQEGV", "VVQEGVLTAVVIPTKKA"]
        )
        assert consensus == "KKDAPYIVGDVVQEGVLTAVVIPTKKA"
        assert positions == [1, 11]
        assert alignments[0] == "KKDAPYIVGDVVQEGV" + "-" * 11
        assert alignments[1] == "-" * 10 + "VVQEGVLTAVVIPTKKA"

    def test_every_alignment_is_the_consensus_width(self):
        consensus, alignments, _ = cluster_mod.build_consensus(
            ["KKDAPYIVGDVVQEGV", "VVQEGVLTAVVIPTKKA", "DAPYIVGD"]
        )
        assert all(len(a) == len(consensus) for a in alignments)

    def test_single_peptide_is_its_own_consensus(self):
        consensus, alignments, positions = cluster_mod.build_consensus(["PHDWTKNITDKIDQ"])
        assert consensus == "PHDWTKNITDKIDQ"
        assert positions == [1]
        assert alignments == ["PHDWTKNITDKIDQ"]

    def test_empty_input(self):
        assert cluster_mod.build_consensus([]) == ("", [], [])


class TestStrainProvenance:
    """Which species a consensus covers must survive into stage 5.

    The pipeline compares three Ebola species on purpose, so a consensus built
    from one species is a different result from one spanning all three. Stage 5
    had no strain column at all, which made those two outcomes identical in every
    output it produced.
    """

    EPITOPE_ROWS = [
        {"epitope_id": "Spike_B_cell_BepiPred_1.1", "sequence": "PHDWTKNITDKIDQ",
         "source_tool": "bepipred", "epitope_class": "B", "strain": "Zaire_ebolavirus"},
        # Same peptide, another species. Deduplication keeps only the first
        # identifier, so reading strain from the resolved row alone would lose
        # this one.
        {"epitope_id": "Spike_B_cell_BepiPred_1.2", "sequence": "PHDWTKNITDKIDQ",
         "source_tool": "bepipred", "epitope_class": "B", "strain": "Sudan_ebolavirus"},
        {"epitope_id": "Spike_T_cell_MHCI_4.1", "sequence": "TKNITDKID",
         "source_tool": "mhc_i", "epitope_class": "MHC-I", "strain": "Zaire_ebolavirus"},
    ]

    def _annotated(self):
        clusters = cluster_mod.cluster_peptides(
            [("Spike_B_cell_BepiPred_1.1", "PHDWTKNITDKIDQ"),
             ("Spike_T_cell_MHCI_4.1", "TKNITDKID")],
            threshold=60,
        )
        cluster_mod.annotate(clusters, self.EPITOPE_ROWS)
        return clusters

    def test_strains_are_aggregated_over_every_row_for_a_peptide(self):
        clusters = self._annotated()
        member = next(
            m for c in clusters for m in c.members
            if m.sequence == "PHDWTKNITDKIDQ"
        )
        assert member.strains == "Sudan_ebolavirus|Zaire_ebolavirus"

    def test_cluster_strains_unions_its_members(self):
        clusters = self._annotated()
        # The 9-mer is an exact substring of the 14-mer, so they cluster.
        combined = clusters[0]
        assert cluster_mod.cluster_strains(combined) == [
            "Sudan_ebolavirus",
            "Zaire_ebolavirus",
        ]

    def test_a_single_species_cluster_reports_only_that_species(self):
        clusters = cluster_mod.cluster_peptides(
            [("Spike_T_cell_MHCI_4.1", "TKNITDKID")], threshold=60
        )
        cluster_mod.annotate(clusters, self.EPITOPE_ROWS)
        assert cluster_mod.cluster_strains(clusters[0]) == ["Zaire_ebolavirus"]

    def test_the_deduplicated_plural_column_is_also_accepted(self):
        # annotate is given all_epitopes.tsv ("strain") today, but must not break
        # if handed the deduplicated table ("strains", pipe-joined).
        clusters = cluster_mod.cluster_peptides(
            [("Spike_B_cell_BepiPred_1.1", "PHDWTKNITDKIDQ")], threshold=60
        )
        cluster_mod.annotate(clusters, [{
            "epitope_id": "Spike_B_cell_BepiPred_1.1",
            "sequence": "PHDWTKNITDKIDQ",
            "source_tool": "bepipred",
            "epitope_class": "B",
            "strains": "Bundibugyo_virus|Zaire_ebolavirus",
        }])
        assert clusters[0].members[0].strains == "Bundibugyo_virus|Zaire_ebolavirus"

    def test_a_defline_with_a_strain_still_yields_the_bare_identifier(self, tmp_path):
        # The whole reason the strain goes after a space: Record.id must stay the
        # identifier that every downstream join is keyed on.
        path = tmp_path / "epitopes.fasta"
        path.write_text(
            ">Spike_B_cell_BepiPred_1.1 Sudan_ebolavirus|Zaire_ebolavirus\n"
            "PHDWTKNITDKIDQ\n",
            encoding="utf-8",
        )
        record = fasta.read_fasta(str(path))[0]
        assert record.id == "Spike_B_cell_BepiPred_1.1"
        assert record.sequence == "PHDWTKNITDKIDQ"


class TestDisagreeingPositions:
    """Majority voting hides disagreement, so it has to be reported separately.

    A consensus that needed eight votes covers its members far less faithfully
    than one that needed none, and after voting both look like ordinary peptides.
    """

    def _cluster(self, peptides):
        entries = [(f"id{index}", peptide) for index, peptide in enumerate(peptides)]
        # threshold 0 puts everything in one cluster regardless of similarity,
        # which is what these cases need.
        return cluster_mod.cluster_peptides(entries, threshold=0, unique=False)[0]

    def test_identical_members_disagree_nowhere(self):
        assert cluster_mod.disagreeing_positions(
            self._cluster(["PHDWTKNITDKIDQ", "PHDWTKNITDKIDQ"])
        ) == 0

    def test_one_differing_position_is_counted(self):
        assert cluster_mod.disagreeing_positions(
            self._cluster(["PHDWTKNITDKIDQ", "PHDWTKNITDKIDA"])
        ) == 1

    def test_a_majority_call_still_counts_as_disagreement(self):
        # The point of the function: the consensus says "Q" and looks clean, but
        # the members did not agree here.
        cluster = self._cluster(
            ["PHDWTKNITDKIDQ", "PHDWTKNITDKIDQ", "PHDWTKNITDKIDA"]
        )
        assert cluster_mod.AMBIGUOUS not in cluster.consensus
        assert cluster_mod.disagreeing_positions(cluster) == 1

    def test_a_singleton_disagrees_nowhere(self):
        assert cluster_mod.disagreeing_positions(self._cluster(["PHDWTKNITDKIDQ"])) == 0


class TestClusterRoundTrip:
    """write_clusters must emit exactly what parse reads.

    Stage 5 is two rules with a file between them, so a mismatch here would not
    show up until select_consensus produced empty output.
    """

    def test_round_trip(self, tmp_path):
        entries = [
            ("Spike_B_cell_BepiPred_1.1", "KKDAPYIVGDVVQEGV"),
            ("Spike_T_cell_MHCI_23.1", "VVQEGVLTAVVIPTKKA"),
            ("Spike_B_cell_EpiDope_9.1", "WWWWWWWWWWWWWW"),
        ]
        original = cluster_mod.cluster_peptides(entries, threshold=60)
        path = str(tmp_path / "cluster_raw.csv")
        cluster_mod.write_clusters(path, original)
        reparsed = cluster_mod.parse(path)

        assert len(reparsed) == len(original)
        for before, after in zip(original, reparsed):
            assert after.cluster_id == before.cluster_id
            assert after.consensus == before.consensus
            assert [m.description for m in after.members] == [
                m.description for m in before.members
            ]
            assert [m.sequence for m in after.members] == [
                m.sequence for m in before.members
            ]
            assert [m.alignment for m in after.members] == [
                m.alignment for m in before.members
            ]

    def test_written_file_has_the_web_tool_header(self, tmp_path):
        path = str(tmp_path / "cluster_raw.csv")
        cluster_mod.write_clusters(
            path, cluster_mod.cluster_peptides([("a", "PEPTIDEPEPTIDE")], threshold=60)
        )
        with open(path, encoding="utf-8") as handle:
            assert handle.readline().strip().split(",")[0] == "Cluster.Sub-Cluster Number"


class TestBepiPredParsing:
    """The real ``raw_output.csv`` layout, taken from bp3/bepipred3.py.

    Written from the tool's source rather than its README, because with IAPred
    the README named a score column that no longer existed and cost a whole
    debugging cycle. The exact line written there is::

        f"{acc},{seq[i].upper()},{avg_prob[i]}, {avg_prob_rolling_mean[i]}"

    note the space before the fourth field, which is reproduced below.
    """

    HEADER = "Accession,Residue,BepiPred-3.0 score,BepiPred-3.0 linear epitope score"

    def _write(self, tmp_path, text):
        (tmp_path / "raw_output.csv").write_text(text, encoding="utf-8")
        return str(tmp_path)

    def test_real_bepipred_csv_format(self, tmp_path):
        out_dir = self._write(
            tmp_path,
            self.HEADER + "\n"
            "Spike|Zaire_ebolavirus|NP_066246.1,M,0.0871, 0.1043\n"
            "Spike|Zaire_ebolavirus|NP_066246.1,G,0.1533, 0.1210\n"
            "Spike|Zaire_ebolavirus|NP_066246.1,V,0.2011, 0.1502\n"
            "Nucleoprotein|Sudan_ebolavirus|WEY07207.1,M,0.3000, 0.2900\n",
        )
        assert bcell.parse_bepipred(out_dir) == {
            "Spike|Zaire_ebolavirus|NP_066246.1": [
                pytest.approx(0.0871),
                pytest.approx(0.1533),
                pytest.approx(0.2011),
            ],
            "Nucleoprotein|Sudan_ebolavirus|WEY07207.1": [pytest.approx(0.3000)],
        }

    def test_the_raw_score_is_taken_not_the_rolling_mean(self, tmp_path):
        # The two columns differ only by smoothing, so picking the wrong one
        # produces entirely plausible numbers and would never be noticed. The raw
        # score is the one BepiPred's own -t is compared against.
        out_dir = self._write(
            tmp_path, self.HEADER + "\nA,M,0.9000, 0.1000\n"
        )
        assert bcell.parse_bepipred(out_dir) == {"A": [pytest.approx(0.9)]}

    def test_the_residue_letter_is_not_read_as_a_position(self, tmp_path):
        # "Residue" holds an amino acid, not an index. If a position candidate
        # ever matched it, float("M") would raise - so this asserts the absence
        # of a match, which is what lets the parser fall back to file order.
        assert (
            bcell.base.find_column(
                self.HEADER.split(","),
                bcell._BEPIPRED_POS_COLUMNS,
                "position",
                required=False,
            )
            is None
        )

    def test_file_order_is_the_residue_order(self, tmp_path):
        # There is no position column, so nothing can reorder these. A parser
        # that sorted on the ID or the score would scramble the sequence.
        descending = "".join(
            f"A,M,{score}, {score}\n" for score in ("0.9", "0.5", "0.7", "0.1")
        )
        out_dir = self._write(tmp_path, self.HEADER + "\n" + descending)
        assert bcell.parse_bepipred(out_dir) == {
            "A": [pytest.approx(v) for v in (0.9, 0.5, 0.7, 0.1)]
        }

    def test_pipes_in_the_defline_survive(self, tmp_path):
        # The whole defline is the accession - BepiPred does not split at
        # whitespace or truncate, unlike IAPred - so stage 3 can attribute each
        # epitope to a strain without a surrogate-ID mapping.
        out_dir = self._write(
            tmp_path, self.HEADER + "\nSpike|Bundibugyo_virus|AYI50316.1,M,0.5, 0.5\n"
        )
        assert list(bcell.parse_bepipred(out_dir)) == ["Spike|Bundibugyo_virus|AYI50316.1"]

    def test_missing_output_is_actionable(self, tmp_path):
        with pytest.raises(bcell.base.ToolError) as exc:
            bcell.parse_bepipred(str(tmp_path))
        assert "raw_output.csv" in str(exc.value)


class TestDeflineIdentifiers:
    """Identifiers must survive as a single whitespace-delimited token.

    ``Record.id`` is ``header.split()[0]``, and every strain in this project's
    inputs contains a space ("Zaire ebolavirus"), so an unsanitised defline is cut
    at that space and the accession is lost from every downstream table.
    """

    def _record(self, protein, strain, accession):
        record = fasta.Record(header="", sequence="ACDEF")
        record.protein_label = protein
        record.strain = strain
        record.accession = accession
        return record

    def test_spaces_in_strain_do_not_split_the_identifier(self):
        header = manifest_mod.make_header(
            self._record("Nucleoprotein", "Zaire ebolavirus", "AAB81001.1")
        )
        assert header == "Nucleoprotein|Zaire_ebolavirus|AAB81001.1"
        assert len(header.split()) == 1  # one token: nothing is lost to a split

    def test_header_round_trips_with_the_accession_intact(self):
        record = self._record("Spike", "Bundibugyo virus", "AYI50316.1")
        protein, strain, accession = manifest_mod.parse_header(
            manifest_mod.make_header(record)
        )
        assert (protein, strain, accession) == ("Spike", "Bundibugyo_virus", "AYI50316.1")

    def test_record_id_keeps_the_whole_header(self):
        record = fasta.Record(
            header=manifest_mod.make_header(
                self._record("Nucleoprotein", "Bundibugyo virus", "AYI50379.1")
            ),
            sequence="ACDEF",
        )
        assert record.id == "Nucleoprotein|Bundibugyo_virus|AYI50379.1"


class TestIAPredSurrogateIdentifiers:
    """IApred truncates deflines to 20 characters, so it is fed surrogates."""

    def test_scores_are_mapped_back_to_the_real_identifiers(self, tmp_path, monkeypatch):
        fasta_path = tmp_path / "in.fasta"
        fasta_path.write_text(
            ">Nucleoprotein|Bundibugyo_virus|AYI50379.1\nACDEFGHIKLMNPQRSTVWY\n"
            ">Nucleoprotein|Zaire_ebolavirus|AAB81001.1\nACDEFGHIKLMNPQRSTVWYACDEF\n"
        )
        raw_path = tmp_path / "raw" / "iapred.csv"
        tool_dir = tmp_path / "tools" / "iapred"
        tool_dir.mkdir(parents=True)
        (tool_dir / "IApred.py").write_text("# stub\n")

        def fake_run(command, cwd=None, check=False):
            # IApred truncates to 20 chars; these surrogates are immune to that.
            os.makedirs(os.path.dirname(str(raw_path)), exist_ok=True)
            raw_path.write_text(
                "Header,Sequence_Length,Intrinsic_Antigenicity_Score,Antigenicity_Category\n"
                "s0,20,0.69,High\n"
                "s1,25,-0.42,Low\n"
            )
            return type("R", (), {"returncode": 0, "stdout": "", "stderr": ""})()

        monkeypatch.setattr(screening.base, "run_command", fake_run)

        scores = screening._predict_iapred(
            fasta_path=str(fasta_path), raw_path=str(raw_path), tool_dir=str(tool_dir)
        )
        assert scores == {
            "Nucleoprotein|Bundibugyo_virus|AYI50379.1": pytest.approx(0.69),
            "Nucleoprotein|Zaire_ebolavirus|AAB81001.1": pytest.approx(-0.42),
        }

    def test_the_mapping_sidecar_is_written(self, tmp_path, monkeypatch):
        fasta_path = tmp_path / "in.fasta"
        fasta_path.write_text(">Spike|Zaire_ebolavirus|NP_066246.1\nACDEFGHIKLMNPQRSTVWY\n")
        raw_path = tmp_path / "iapred.csv"
        tool_dir = tmp_path / "tools"
        tool_dir.mkdir()
        (tool_dir / "IApred.py").write_text("# stub\n")

        def fake_run(command, cwd=None, check=False):
            raw_path.write_text("Header,Intrinsic_Antigenicity_Score\ns0,0.5\n")
            return type("R", (), {"returncode": 0, "stdout": "", "stderr": ""})()

        monkeypatch.setattr(screening.base, "run_command", fake_run)
        screening._predict_iapred(
            fasta_path=str(fasta_path), raw_path=str(raw_path), tool_dir=str(tool_dir)
        )

        sidecar = (tmp_path / "iapred.ids.tsv").read_text()
        assert "s0\tSpike|Zaire_ebolavirus|NP_066246.1" in sidecar


class TestAntigenicityBackendSelection:
    def test_unscorable_peptide_is_recorded_not_scored(self, tmp_path, monkeypatch):
        """A peptide IApred declines to score fails the screen with its own note.

        It must not abort the stage, and must not be reported as "below
        threshold" - there is no number to compare.
        """
        path = tmp_path / "x.fasta"
        path.write_text(
            ">Good_1 spike [Zaire]\nACDEFGHIKLMNPQRSTVWYACDEF\n"
            ">Short_1 spike [Zaire]\nACDEFGHIK\n"
        )
        monkeypatch.setattr(
            screening, "_predict_iapred", lambda **kwargs: {"Good_1": 0.69, "Short_1": None}
        )

        good, short = screening.predict_antigenicity(
            str(path), str(tmp_path / "out.tsv"), threshold=0.3, tool="iapred"
        )

        assert good["antigenicity_score"] == pytest.approx(0.69)
        assert good["antigenic"] is True
        assert good["note"] == ""

        assert short["antigenicity_score"] == ""
        assert short["antigenic"] is False
        assert short["note"] == screening.IAPRED_UNSCORABLE

    def test_unknown_backend_is_rejected(self, tmp_path):
        path = tmp_path / "x.fasta"
        path.write_text(">x test [strain]\nACDEF\n")
        with pytest.raises(screening.base.ToolError, match="Unknown antigenicity backend"):
            screening.predict_antigenicity(
                str(path), str(tmp_path / "out.tsv"), threshold=0.4, tool="nonsense"
            )

    def test_external_backend_requires_a_command(self, tmp_path):
        path = tmp_path / "x.fasta"
        path.write_text(">x test [strain]\nACDEF\n")
        with pytest.raises(screening.base.ToolError, match="antigenicity.command is empty"):
            screening.predict_antigenicity(
                str(path), str(tmp_path / "out.tsv"), threshold=0.4, tool="external", command=""
            )

    def test_missing_iapred_install_is_actionable(self, tmp_path):
        path = tmp_path / "x.fasta"
        path.write_text(">x test [strain]\nACDEF\n")
        with pytest.raises(screening.base.ToolMissing, match="install_iapred"):
            screening.predict_antigenicity(
                str(path), str(tmp_path / "out.tsv"), threshold=0.4,
                tool="iapred", tool_dir=str(tmp_path / "does-not-exist"),
            )


# ---------------------------------------------------------------------------
# the input folder, and the rules it documents
# ---------------------------------------------------------------------------
# The naming rule used to be written out by hand in five places, and the copies
# had drifted from the parser: the documented filename convention is only a
# fallback, and not one of the project's own six input files satisfied it. These
# tests push the notice's *own examples* through the real parser, so the text
# users read cannot describe something the code does not do.

_PROTEIN = "ACDEFGHIKLMNPQRSTVWY" * 3      # 60 aa, valid and clearly not nucleotide
_EXTENSIONS = [".fasta", ".fa", ".faa"]
_ALIASES = {
    "spike glycoprotein": "Spike",
    "glycoprotein": "Spike",
    "nucleoprotein": "Nucleoprotein",
}


class TestInputFolder:
    def _write(self, directory, name, text):
        path = directory / name
        path.write_text(text, encoding="utf-8")
        return path

    def _only_record(self, path):
        (record,) = fasta.read_fasta(str(path))
        return record

    # --- the documented examples parse as documented --------------------------

    def test_documented_header_parses_as_documented(self, tmp_path):
        path = self._write(tmp_path, "any name.fasta", inputs.EXAMPLE_DEFLINE + "\n" + _PROTEIN + "\n")
        record = self._only_record(path)
        assert record.accession == "NP_066246.1"
        assert record.name == "spike glycoprotein"
        assert record.strain == "Zaire ebolavirus"

    def test_documented_filename_is_the_fallback_for_a_bare_header(self, tmp_path):
        path = self._write(tmp_path, inputs.EXAMPLE_FILENAME, ">NP_066246.1\n" + _PROTEIN + "\n")
        record = self._only_record(path)
        assert (record.accession, record.name, record.strain) == (
            "NP_066246.1", "spike glycoprotein", "Zaire ebolavirus",
        )

    def test_a_complete_header_makes_the_filename_irrelevant(self, tmp_path):
        """Why the real input files work despite not matching the filename rule."""
        path = self._write(
            tmp_path, "AAB81001.1 nucleoprotein -Zaire ebolavirus.fasta",   # no space after '-'
            ">AAB81001.1 nucleoprotein [Zaire ebolavirus]\n" + _PROTEIN + "\n",
        )
        record = self._only_record(path)
        assert record.name == "nucleoprotein"
        assert record.strain == "Zaire ebolavirus"

    def test_filename_fallback_needs_whitespace_on_both_sides_of_the_dash(self, tmp_path):
        """The notice says "both sides". This is the proof, and the way it bites.

        "-Zaire" has no space after the dash, so a bare header gets nothing from
        the filename and both fields fall to the sentinels.
        """
        path = self._write(
            tmp_path, "AAB81001.1 nucleoprotein -Zaire ebolavirus.fasta",
            ">AAB81001.1\n" + _PROTEIN + "\n",
        )
        record = self._only_record(path)
        assert record.name == inputs.UNKNOWN_PROTEIN
        assert record.strain == inputs.UNKNOWN_STRAIN

    # --- creating the folder --------------------------------------------------

    def test_creates_the_folder_and_its_readme(self, tmp_path):
        target = tmp_path / "0 - Input"
        assert inputs.ensure_input_dir(str(target), _EXTENSIONS, _ALIASES) is True
        readme = (target / inputs.README_NAME).read_text(encoding="utf-8")
        assert inputs.EXAMPLE_DEFLINE in readme
        assert inputs.EXAMPLE_FILENAME in readme

    def test_running_again_changes_nothing(self, tmp_path):
        """Stage 1 declares this folder as an input, and Snakemake watches a
        directory's modification time. A repeat run that touched the folder
        would send the whole pipeline back to the start."""
        target = tmp_path / "0 - Input"
        inputs.ensure_input_dir(str(target), _EXTENSIONS, _ALIASES)
        readme = target / inputs.README_NAME
        before = (target.stat().st_mtime_ns, readme.stat().st_mtime_ns)

        assert inputs.ensure_input_dir(str(target), _EXTENSIONS, _ALIASES) is False

        assert (target.stat().st_mtime_ns, readme.stat().st_mtime_ns) == before

    def test_never_touches_sequence_files(self, tmp_path):
        target = tmp_path / "0 - Input"
        target.mkdir()
        original = inputs.EXAMPLE_DEFLINE + "\n" + _PROTEIN + "\n"
        sequence = self._write(target, inputs.EXAMPLE_FILENAME, original)

        inputs.ensure_input_dir(str(target), _EXTENSIONS, _ALIASES)

        assert sequence.read_text(encoding="utf-8") == original

    def test_readme_is_not_mistaken_for_a_sequence(self, tmp_path):
        target = tmp_path / "0 - Input"
        inputs.ensure_input_dir(str(target), _EXTENSIONS, _ALIASES)
        assert inputs.find_inputs(str(target), _EXTENSIONS) == []

    def test_subfolders_are_ignored(self, tmp_path):
        target = tmp_path / "0 - Input"
        (target / "old").mkdir(parents=True)
        self._write(target / "old", "x.fasta", ">a b [c]\n" + _PROTEIN + "\n")
        self._write(target, "keep.fasta", ">a b [c]\n" + _PROTEIN + "\n")
        assert [os.path.basename(p) for p in inputs.find_inputs(str(target), _EXTENSIONS)] == [
            "keep.fasta"
        ]

    def test_readme_follows_the_configured_aliases(self, tmp_path):
        """It is generated from config.yaml, so it must track it."""
        target = tmp_path / "0 - Input"
        inputs.ensure_input_dir(str(target), _EXTENSIONS, _ALIASES)
        assert '"glycoprotein"  ->  Spike' in (target / inputs.README_NAME).read_text(encoding="utf-8")

        inputs.ensure_input_dir(str(target), _EXTENSIONS, {**_ALIASES, "matrix protein": "VP40"})
        assert '"matrix protein"  ->  VP40' in (target / inputs.README_NAME).read_text(encoding="utf-8")

    # --- the preflight report -------------------------------------------------

    def test_describe_reports_how_a_file_will_be_read(self, tmp_path):
        self._write(tmp_path, inputs.EXAMPLE_FILENAME, inputs.EXAMPLE_DEFLINE + "\n" + _PROTEIN + "\n")
        lines, ok = inputs.describe(str(tmp_path), _EXTENSIONS, _ALIASES)
        text = "\n".join(lines)
        assert ok is True
        assert "protein Spike | strain Zaire ebolavirus | accession NP_066246.1" in text

    def test_describe_warns_about_an_unknown_strain_without_failing(self, tmp_path):
        """Stage 1 accepts it, so the report must too - but it is worth saying."""
        self._write(tmp_path, "mystery.fasta", ">AAB81001.1\n" + _PROTEIN + "\n")
        lines, ok = inputs.describe(str(tmp_path), _EXTENSIONS, _ALIASES)
        text = "\n".join(lines)
        assert ok is True
        assert "WARNINGS" in text
        assert inputs.UNKNOWN_STRAIN in text

    def test_describe_fails_on_a_duplicate_pair(self, tmp_path):
        self._write(tmp_path, "a.fasta", ">A1 spike glycoprotein [Zaire ebolavirus]\n" + _PROTEIN + "\n")
        self._write(tmp_path, "b.fasta", ">B2 spike glycoprotein [Zaire ebolavirus]\n" + _PROTEIN + "\n")
        lines, ok = inputs.describe(str(tmp_path), _EXTENSIONS, _ALIASES)
        assert ok is False
        assert "duplicate protein/strain pair" in "\n".join(lines)

    def test_describe_fails_on_a_nucleotide_sequence(self, tmp_path):
        self._write(tmp_path, "a.fasta", ">A1 spike glycoprotein [Zaire ebolavirus]\n" + "ACGT" * 20 + "\n")
        lines, ok = inputs.describe(str(tmp_path), _EXTENSIONS, _ALIASES)
        assert ok is False
        assert "nucleotide" in "\n".join(lines)

    def test_describe_treats_an_empty_folder_as_fine(self, tmp_path):
        lines, ok = inputs.describe(str(tmp_path), _EXTENSIONS, _ALIASES)
        assert ok is True
        assert "No sequence files" in lines[0]

    def test_describe_and_stage_1_refuse_the_same_inputs(self, tmp_path):
        """One definition of "acceptable", used by both - not two that can drift."""
        folder = tmp_path / "in"
        folder.mkdir()
        self._write(folder, "a.fasta", ">A1 spike glycoprotein [Zaire ebolavirus]\n" + _PROTEIN + "\n")
        self._write(folder, "b.fasta", ">B2 spike glycoprotein [Zaire ebolavirus]\n" + _PROTEIN + "\n")

        _, ok = inputs.describe(str(folder), _EXTENSIONS, _ALIASES)
        with pytest.raises(SystemExit, match="duplicate protein/strain pair"):
            manifest_mod.build(
                str(folder), str(tmp_path / "manifest.tsv"), str(tmp_path / "norm"),
                _EXTENSIONS, _ALIASES,
            )
        assert ok is False

    # --- the errors users actually meet ---------------------------------------

    def test_an_empty_folder_error_carries_the_rules(self, tmp_path):
        folder = tmp_path / "0 - Input"
        folder.mkdir()
        with pytest.raises(SystemExit) as error:
            manifest_mod.collect(str(folder), _EXTENSIONS, _ALIASES)
        assert inputs.EXAMPLE_DEFLINE in str(error.value)
        assert "init_input.py" in str(error.value)

    def test_a_missing_folder_error_carries_the_rules(self, tmp_path):
        with pytest.raises(SystemExit) as error:
            manifest_mod.collect(str(tmp_path / "nope"), _EXTENSIONS, _ALIASES)
        assert inputs.EXAMPLE_DEFLINE in str(error.value)

    def test_kolaskar_backend_scores_every_record(self, tmp_path):
        path = tmp_path / "x.fasta"
        path.write_text(">a p [s]\nACDEFGHIK\n>b p [s]\nNNNNNNNNN\n")
        rows = screening.predict_antigenicity(
            str(path), str(tmp_path / "out.tsv"), threshold=1.0, tool="kolaskar"
        )
        assert [row["sequence_id"] for row in rows] == ["a", "b"]
        assert all(row["tool"] == "kolaskar-tongaonkar-1990" for row in rows)
        # Threshold of 1.0 sits inside the propensity range, so the two records
        # must fall on opposite sides of it.
        assert rows[0]["antigenic"] is True
        assert rows[1]["antigenic"] is False
