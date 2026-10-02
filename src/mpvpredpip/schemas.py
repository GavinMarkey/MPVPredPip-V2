"""Canonical table schemas passed between pipeline stages.

Every stage reads and writes tab-separated files with these exact columns. Tool
specific formats are converted at the edges (see ``mpvpredpip.adapters``) so that
no downstream stage ever has to know which predictor produced a row.
"""

from __future__ import annotations

# --- stage 1 -----------------------------------------------------------------
MANIFEST = [
    "accession",       # NCBI accession, e.g. NP_066246.1
    "protein_label",   # short normalised label used in epitope IDs, e.g. Spike
    "protein_name",    # full name from the defline, e.g. spike glycoprotein
    "strain",          # e.g. Zaire ebolavirus
    "length",          # residues
    "md5",             # checksum of the ungapped uppercase sequence
    "path",            # source FASTA, relative to the project root
]

# --- stages 3 and 4 ----------------------------------------------------------
# One row per predicted epitope, from any predictor.
EPITOPES = [
    "epitope_id",      # e.g. Spike_B_cell_BepiPred_12.1
    "protein_label",
    "strain",
    "accession",
    "source_tool",     # bepipred | epidope | mhc_i | mhc_ii
    "epitope_class",   # B | MHC-I | MHC-II
    "start",           # 1-based, inclusive, on the source protein
    "end",             # 1-based, inclusive
    "length",
    "sequence",
    "score",           # tool-native score (see source_tool for meaning)
    "allele",          # MHC allele, empty for B-cell epitopes
    "percentile",      # IEDB percentile rank, empty for B-cell epitopes
]

# Deduplicated epitopes, ready for clustering. ``strains`` and ``alleles`` are
# pipe-delimited aggregates of every row that collapsed into this peptide.
EPITOPES_NR = [
    "epitope_id",
    "protein_label",
    "source_tool",
    "epitope_class",
    "sequence",
    "length",
    "best_score",
    "strains",
    "alleles",
    "n_occurrences",
]

# --- stage 5 -----------------------------------------------------------------
CLUSTERS = [
    "protein_label",
    "cluster_id",          # e.g. 2.2  (cluster.sub-cluster)
    "consensus_id",        # e.g. Spike_Consensus_2.2
    "consensus_sequence",
    "consensus_length",
    "n_members",
    "n_bepipred",
    "n_epidope",
    "n_mhc_i",
    "n_mhc_ii",
    # Which species this consensus is built from. A candidate supported by one
    # species is a different result from one spanning all three, and this is a
    # three-species design - without these two columns the outputs could not tell
    # them apart.
    "strains",             # pipe-delimited, e.g. Sudan_ebolavirus|Zaire_ebolavirus
    "n_strains",
    # How many consensus positions the members did not all agree on. The
    # consensus is a per-column majority vote, so disagreement no longer shows up
    # as an X in the sequence; 0 means every member matched the consensus exactly.
    "n_disagreeing",
    "passes",              # True | False
    "reason",              # why it failed, empty when it passes
]

CLUSTER_MEMBERS = [
    "protein_label",
    "cluster_id",
    "consensus_id",
    "member_index",
    "epitope_id",
    "source_tool",
    "epitope_class",
    "strains",             # pipe-delimited species this peptide was predicted in
    "position",            # offset of the member within the consensus, 1-based
    "alignment",           # gapped alignment string
    "sequence",
]

# --- stages 2 and 6 ----------------------------------------------------------
ANTIGENICITY = [
    "sequence_id",
    "sequence",
    "length",
    "antigenicity_score",  # empty when the predictor declined to score it
    "antigenic",           # True | False against the configured threshold
    "threshold",
    "tool",
    "note",                # why there is no score, empty when there is one
]

SCREENING = [
    "consensus_id",
    "protein_label",
    "cluster_id",
    "sequence",
    "length",
    # Carried through from clusters.tsv. Stage 5 resolves which species each
    # consensus covers, and that is the single most important property of a
    # candidate here - without it the final table cannot distinguish a
    # pan-species peptide from a single-strain one.
    "strains",
    "n_strains",
    "antigenicity_score",
    "antigenic",
    "allergen_score",
    "allergen",            # True means predicted allergen -> fails
    "human_match_count",
    "human_match_peptides",
    "autoimmune_risk",     # True means it matches the human proteome -> fails
    # Windows that could not be searched because they hold an ambiguous
    # residue. Non-zero means the autoimmunity verdict is incomplete, NOT that
    # the peptide came back clean - see mpvpredpip.autoimmunity.
    "windows_skipped",
    "passes",
    "reason",
]

# --- stage 7 -----------------------------------------------------------------
CONSERVANCY = [
    "consensus_id",
    "protein_label",
    "sequence",
    "n_target_sequences",
    "n_matched",
    "conservancy_fraction",
    "conservancy_percent",
    "min_identity",
    "max_identity",
    "matched_strains",
]

# --- stage 8 -----------------------------------------------------------------
POPULATION_COVERAGE = [
    "consensus_id",
    "protein_label",
    "population",
    "mhc_class",
    "coverage_percent",
    "average_hit",
    "pc90",
    "n_alleles",
    "alleles",
]

# --- stage 9 -----------------------------------------------------------------
FINAL = [
    "rank",
    "consensus_id",
    "protein_label",
    "cluster_id",
    "sequence",
    "length",
    # Which species the candidate's cluster actually spanned, carried from
    # stage 5 through stage 6. For a pan-ebolavirus design this is the headline
    # property of a candidate: a peptide covering all three species is a
    # different proposition from one seen in a single strain, and without these
    # the final table cannot tell them apart.
    "strains",
    "n_strains",
    "antigenicity_score",
    "allergen_score",
    "human_match_count",
    "conservancy_percent",
    "world_coverage_percent",
    "n_bepipred",
    "n_epidope",
    "n_mhc_i",
    "n_mhc_ii",
    "mhc_i_alleles",
    "mhc_ii_alleles",
]

#: Epitope class tags used when building epitope IDs.
CLASS_TAGS = {
    "bepipred": ("B_cell", "BepiPred", "B"),
    "epidope": ("B_cell", "EpiDope", "B"),
    "mhc_i": ("T_cell", "MHCI", "MHC-I"),
    "mhc_ii": ("T_cell", "MHCII", "MHC-II"),
}
