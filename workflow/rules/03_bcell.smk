# Stage 3 - linear B-cell epitope prediction.
#
# BepiPred and EpiDope each run once per protein (the per-protein FASTA already
# holds every strain, so one invocation covers them all and keeps per-strain
# attribution through the "Protein|Strain|Accession" deflines).
#
# The two predictors need mutually incompatible dependency stacks, so each rule
# declares its own conda environment. That is the main reason this workflow uses
# Snakemake rather than a shell script: the environments are resolved per rule,
# not per pipeline.
#
# Both tools emit per-residue scores. Converting those to epitope intervals is
# done by the shared code in mpvpredpip.epitopes so that BepiPred and EpiDope
# epitopes are built by identical rules and remain comparable in stage 5.
#
# Both predictor rules use `shell:` rather than `script:`. Snakemake's script
# directive imports Snakemake inside the rule's environment, and EpiDope's
# environment pins Python 3.8 - too old for Snakemake 9's own syntax. Passing
# parameters on the command line avoids that import entirely. Every path is
# quoted because the stage folders contain spaces.


rule bepipred:
    input:
        fasta=f"{NORMALISED}/{{protein}}.fasta",
    output:
        residues=f"{S3}/{{protein}}/bepipred/per_residue_scores.tsv",
        epitopes=f"{S3}/{{protein}}/bepipred/epitopes.tsv",
    params:
        raw_dir=f"{S3}/{{protein}}/bepipred/raw",
        # Deliberately outside the project directory - see the comment on
        # bcell.bepipred.esm_cache_dir in config.yaml. BepiPred names each
        # encoding after the defline, and "|" cannot be written to a Windows
        # bind mount.
        esm_dir=f"{config['bcell']['bepipred']['esm_cache_dir']}/{{protein}}",
        tool_dir=config["tools"]["bepipred"],
        threshold=config["bcell"]["bepipred"]["threshold"],
        length=config["bcell"]["bepipred"]["epitope_length"],
        step=config["bcell"]["bepipred"]["epitope_step"],
        min_region=config["bcell"]["bepipred"]["min_region_length"],
        merge_gap=config["bcell"]["bepipred"]["merge_gap"],
        source_tool="bepipred",
    threads: config["resources"]["threads_per_job"]
    log:
        f"{S3}/{{protein}}/logs/bepipred.log",
    conda:
        "../envs/bepipred.yaml"
    shell:
        # Start from an empty raw directory for the same reason the EpiDope rule
        # does: a file left by a partial earlier run must not be parsed as if it
        # belonged to this one. The ESM cache is deliberately not under it.
        "rm -rf '{params.raw_dir}' && mkdir -p '{params.raw_dir}' &&"
        " python workflow/scripts/03_bcell_predict.py"
        " --fasta '{input.fasta}'"
        " --residues '{output.residues}'"
        " --epitopes '{output.epitopes}'"
        " --raw-dir '{params.raw_dir}'"
        " --esm-dir '{params.esm_dir}'"
        " --tool-dir '{params.tool_dir}'"
        " --source-tool '{params.source_tool}'"
        " --threshold {params.threshold}"
        " --length {params.length}"
        " --step {params.step}"
        " --min-region {params.min_region}"
        " --merge-gap {params.merge_gap}"
        " --threads {threads}"
        " --log '{log}'"


# EpiDope is split across two rules. Its environment (envs/epidope.yml, the
# upstream project's own export) pins Python 3.6, which cannot even parse this
# codebase - `from __future__ import annotations` requires 3.7. So the tool runs
# there as a bare shell command, and its output is parsed in the next rule under
# the orchestrator's environment. Nothing of ours enters the 3.6 environment.
rule epidope_run:
    input:
        fasta=f"{NORMALISED}/{{protein}}.fasta",
    output:
        marker=f"{S3}/{{protein}}/epidope/raw/.done",
    params:
        raw_dir=f"{S3}/{{protein}}/epidope/raw",
        executable=config["tools"]["epidope"],
        threshold=config["bcell"]["epidope"]["threshold"],
        slice_length=config["bcell"]["epidope"]["slice_length"],
        slice_shift=config["bcell"]["epidope"]["slice_shift"],
    threads: config["resources"]["threads_per_job"]
    log:
        f"{S3}/{{protein}}/logs/epidope.log",
    conda:
        "../envs/epidope.yml"
    shell:
        # Start from an empty directory. EpiDope opens its outputs with 'w' and
        # fails outright if it cannot truncate one left by an earlier run, and a
        # stale file surviving a partial run could otherwise be parsed as if it
        # belonged to this one. The directory is this rule's own output.
        "rm -rf '{params.raw_dir}' && mkdir -p '{params.raw_dir}' &&"
        " {params.executable} -i '{input.fasta}' -o '{params.raw_dir}'"
        " -t {params.threshold}"
        " -l {params.slice_length} -s {params.slice_shift}"
        " -p {threads} > '{log}' 2>&1 &&"
        " touch '{output.marker}'"


rule epidope_epitopes:
    input:
        fasta=f"{NORMALISED}/{{protein}}.fasta",
        marker=f"{S3}/{{protein}}/epidope/raw/.done",
    output:
        residues=f"{S3}/{{protein}}/epidope/per_residue_scores.tsv",
        epitopes=f"{S3}/{{protein}}/epidope/epitopes.tsv",
    params:
        raw_dir=f"{S3}/{{protein}}/epidope/raw",
        threshold=config["bcell"]["epidope"]["threshold"],
        length=config["bcell"]["epidope"]["epitope_length"],
        step=config["bcell"]["epidope"]["epitope_step"],
        min_region=config["bcell"]["epidope"]["min_region_length"],
        merge_gap=config["bcell"]["epidope"]["merge_gap"],
        source_tool="epidope",
    log:
        f"{S3}/{{protein}}/logs/epidope_epitopes.log",
    conda:
        "../envs/core.yaml"
    shell:
        "python workflow/scripts/03_bcell_predict.py"
        " --skip-predict"
        " --fasta '{input.fasta}'"
        " --residues '{output.residues}'"
        " --epitopes '{output.epitopes}'"
        " --raw-dir '{params.raw_dir}'"
        " --source-tool '{params.source_tool}'"
        " --threshold {params.threshold}"
        " --length {params.length}"
        " --step {params.step}"
        " --min-region {params.min_region}"
        " --merge-gap {params.merge_gap}"
        " --log '{log}'"


rule bcell_merge:
    input:
        bepipred=f"{S3}/{{protein}}/bepipred/epitopes.tsv",
        epidope=f"{S3}/{{protein}}/epidope/epitopes.tsv",
    output:
        table=f"{S3}/{{protein}}/bcell_epitopes.tsv",
        fasta=f"{S3}/{{protein}}/bcell_epitopes.fasta",
        summary=f"{S3}/{{protein}}/summary.md",
    log:
        f"{S3}/{{protein}}/logs/merge.log",
    conda:
        "../envs/core.yaml"
    script:
        "../scripts/03_bcell_merge.py"
