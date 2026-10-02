# Stage 1 - sequence retrieval, validation and normalisation.
#
# The user supplies NCBI protein FASTA files in "0 - Input" (sequences.input_dir
# in config.yaml) - a folder this pipeline only ever reads, which is why the
# outputs below go to the separate, wholly generated stage 1 folder. This rule
# validates the inputs (protein sequence, not nucleotide; no duplicate
# protein/strain pairs; above the minimum length) and emits the manifest plus
# normalised per-protein and per-protein-per-strain FASTA files with stable
# deflines of the form "Protein|Strain|Accession". Everything downstream keys
# off those.
#
# The naming rules are in vaxpipe.inputs, which is also the source of the
# messages users see - run `python setup/init_input.py` to read or check them.


rule sequence_manifest:
    input:
        directory=INPUT_DIR,
    output:
        manifest=f"{S1}/manifest.tsv",
        proteins=expand(f"{NORMALISED}/{{protein}}.fasta", protein=PROTEINS),
    params:
        out_dir=NORMALISED,
        extensions=config["sequences"]["extensions"],
        aliases=config["sequences"].get("protein_aliases", {}),
        min_length=config["sequences"].get("min_length", 50),
    log:
        f"{S1}/logs/manifest.log",
    conda:
        "../envs/core.yaml"
    script:
        "../scripts/01_manifest.py"
