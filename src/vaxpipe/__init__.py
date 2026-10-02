"""vaxpipe - linear peptide vaccine-candidate prediction pipeline.

Stage-by-stage helper library backing the Snakemake workflow in ``workflow/``.
Each module maps onto one numbered stage folder:

    manifest     stage 1  sequence retrieval and validation
    adapters     stages 2-4, 6, 8  external predictor wrappers
    epitopes     stages 3-4  per-residue scores to epitope intervals, IDs, dedup
    cluster      stage 5  IEDB cluster parsing and B/T consensus selection
    conservancy  stage 7  conservancy across strains
    report       stages 1-9  human-readable summaries

The design rule throughout: external tools are wrapped at the edge and their
native output is kept untouched; every stage boundary is a plain TSV following
:mod:`vaxpipe.schemas`.
"""

__version__ = "1.0.0"

__all__ = [
    "adapters",
    "cluster",
    "conservancy",
    "epitopes",
    "fasta",
    "manifest",
    "report",
    "schemas",
    "tables",
]
