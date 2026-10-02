# Stage 4 - T-cell epitope prediction (local Tepitool equivalent).
#
# Tepitool has no standalone release; it is a guided web front end over the IEDB
# binding predictors. Here the same job is done locally by calling the IEDB
# standalones with an explicit allele panel (workflow/config/alleles/) and an
# explicit percentile cut-off (config.yaml), so the selection criteria are
# versioned with the analysis instead of being hidden defaults on a web form.
#
# Class I binders are additionally screened with the IEDB immunogenicity
# predictor; only peptides scoring above tcell.immunogenicity.min_score are
# carried into stage 5, as required by the pipeline specification.
#
# These rules use `shell:` rather than `script:`: the IEDB environment pins
# Python 3.8, which cannot import Snakemake 9 (its source uses 3.10+ syntax), and
# `script:` requires that import. See 03_bcell_predict.py. Paths are quoted
# because the stage folders contain spaces, and list parameters are passed
# comma-separated.


rule mhc_i_binding:
    input:
        fasta=f"{NORMALISED}/{{protein}}.fasta",
    output:
        raw=f"{S4}/{{protein}}/MHC-I/raw_binding.tsv",
        table=f"{S4}/{{protein}}/MHC-I/binders.tsv",
    params:
        tool_dir=config["tools"]["iedb_mhc_i"],
        allele_file=config["tcell"]["mhc_i"]["allele_file"],
        method=config["tcell"]["mhc_i"]["method"],
        lengths=",".join(str(n) for n in config["tcell"]["mhc_i"]["peptide_lengths"]),
        cutoff=config["tcell"]["mhc_i"]["percentile_cutoff"],
        mhc_class="I",
        source_tool="mhc_i",
    threads: config["resources"]["threads_per_job"]
    log:
        f"{S4}/{{protein}}/logs/mhc_i.log",
    conda:
        "../envs/iedb.yaml"
    shell:
        "python workflow/scripts/04_mhc_binding.py"
        " --fasta '{input.fasta}'"
        " --table '{output.table}'"
        " --raw '{output.raw}'"
        " --tool-dir '{params.tool_dir}'"
        " --allele-file '{params.allele_file}'"
        " --lengths '{params.lengths}'"
        " --cutoff {params.cutoff}"
        " --method '{params.method}'"
        " --mhc-class '{params.mhc_class}'"
        " --source-tool '{params.source_tool}'"
        " --log '{log}'"


rule mhc_ii_binding:
    input:
        fasta=f"{NORMALISED}/{{protein}}.fasta",
    output:
        raw=f"{S4}/{{protein}}/MHC-II/raw_binding.tsv",
        table=f"{S4}/{{protein}}/MHC-II/binders.tsv",
    params:
        tool_dir=config["tools"]["iedb_mhc_ii"],
        allele_file=config["tcell"]["mhc_ii"]["allele_file"],
        method=config["tcell"]["mhc_ii"]["method"],
        lengths=",".join(str(n) for n in config["tcell"]["mhc_ii"]["peptide_lengths"]),
        cutoff=config["tcell"]["mhc_ii"]["percentile_cutoff"],
        mhc_class="II",
        source_tool="mhc_ii",
    threads: config["resources"]["threads_per_job"]
    log:
        f"{S4}/{{protein}}/logs/mhc_ii.log",
    conda:
        "../envs/iedb.yaml"
    shell:
        "python workflow/scripts/04_mhc_binding.py"
        " --fasta '{input.fasta}'"
        " --table '{output.table}'"
        " --raw '{output.raw}'"
        " --tool-dir '{params.tool_dir}'"
        " --allele-file '{params.allele_file}'"
        " --lengths '{params.lengths}'"
        " --cutoff {params.cutoff}"
        " --method '{params.method}'"
        " --mhc-class '{params.mhc_class}'"
        " --source-tool '{params.source_tool}'"
        " --log '{log}'"


rule mhc_i_immunogenicity:
    input:
        binders=f"{S4}/{{protein}}/MHC-I/binders.tsv",
    output:
        raw=f"{S4}/{{protein}}/Immunogenicity/raw_immunogenicity.tsv",
        table=f"{S4}/{{protein}}/Immunogenicity/immunogenic_binders.tsv",
        dropped=f"{S4}/{{protein}}/Immunogenicity/non_immunogenic_dropped.tsv",
    params:
        tool_dir=config["tools"]["iedb_immunogenicity"],
        enabled=config["tcell"]["immunogenicity"]["enabled"],
        min_score=config["tcell"]["immunogenicity"]["min_score"],
    log:
        f"{S4}/{{protein}}/logs/immunogenicity.log",
    conda:
        "../envs/iedb.yaml"
    shell:
        "python workflow/scripts/04_immunogenicity.py"
        " --binders '{input.binders}'"
        " --table '{output.table}'"
        " --dropped '{output.dropped}'"
        " --raw '{output.raw}'"
        " --tool-dir '{params.tool_dir}'"
        " --min-score {params.min_score}"
        " --enabled '{params.enabled}'"
        " --log '{log}'"


rule tcell_merge:
    input:
        mhc_i=f"{S4}/{{protein}}/Immunogenicity/immunogenic_binders.tsv",
        mhc_ii=f"{S4}/{{protein}}/MHC-II/binders.tsv",
    output:
        table=f"{S4}/{{protein}}/tcell_epitopes.tsv",
        fasta=f"{S4}/{{protein}}/tcell_epitopes.fasta",
        summary=f"{S4}/{{protein}}/summary.md",
    log:
        f"{S4}/{{protein}}/logs/merge.log",
    conda:
        "../envs/core.yaml"
    script:
        "../scripts/04_tcell_merge.py"
