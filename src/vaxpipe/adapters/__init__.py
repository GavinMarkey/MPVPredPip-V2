"""External predictor adapters.

Each module wraps third-party tools that this pipeline does not control:

    bcell      BepiPred, EpiDope                       (stage 3)
    iedb       MHC-I/II binding, immunogenicity,
               epitope clustering, population coverage (stages 4, 5, 8)
    screening  antigenicity, AlgPred, PIR Peptide Match (stages 2, 6)

All of them follow the contract documented in :mod:`vaxpipe.adapters.base`:
``run`` keeps native output untouched on disk, ``parse`` resolves columns by name
and fails loudly with the header it actually saw.
"""

from . import base, bcell, iedb, screening

__all__ = ["base", "bcell", "iedb", "screening"]
