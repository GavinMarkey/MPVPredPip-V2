"""Stage 1 - build the sequence manifest and normalised FASTA files.

No ``from __future__ import annotations`` here, or in any script run through
Snakemake's ``script:`` directive: Snakemake prepends its own preamble to inject
the ``snakemake`` object, so a future import is no longer the first statement and
Python rejects the file. Annotations in these scripts must therefore be valid at
runtime under the oldest interpreter the rule's environment provides.
"""

import _common

_common.bootstrap()

from vaxpipe import manifest as manifest_mod  # noqa: E402

smk = snakemake  # noqa: F821  (injected by Snakemake)

with _common.tee_log(_common.log_path(smk)):
    records = manifest_mod.build(
        input_dir=str(smk.input.directory),
        out_manifest=str(smk.output.manifest),
        out_dir=str(smk.params.out_dir),
        extensions=list(smk.params.extensions),
        aliases=dict(smk.params.aliases),
        min_length=int(smk.params.min_length),
    )

    by_protein: dict[str, list[str]] = {}
    for record in records:
        by_protein.setdefault(record.protein_label, []).append(record.strain)

    print(f"Wrote manifest with {len(records)} sequences.")
    for protein, strains in sorted(by_protein.items()):
        print(f"  {protein}: {len(strains)} strain(s) - {', '.join(sorted(strains))}")
