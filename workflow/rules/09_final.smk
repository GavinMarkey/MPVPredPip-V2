# Stage 9 - final candidate assembly.
#
# Joins every per-protein table into one ranked list of consensus epitopes that
# passed every enabled screen, and writes a report recording how many candidates
# were lost at each stage and why. The attrition table is the quickest way to see
# whether a threshold is doing something unintended - a stage that drops
# everything usually means a mis-set cut-off rather than a real biological result.


rule final_candidates:
    input:
        manifest=f"{S1}/manifest.tsv",
        antigenicity=f"{S2}/whole_protein_antigenicity.tsv",
        bcell=expand(f"{S3}/{{protein}}/bcell_epitopes.tsv", protein=PROTEINS),
        tcell=expand(f"{S4}/{{protein}}/tcell_epitopes.tsv", protein=PROTEINS),
        clusters=expand(f"{S5}/{{protein}}/clusters.tsv", protein=PROTEINS),
        members=expand(f"{S5}/{{protein}}/cluster_members.tsv", protein=PROTEINS),
        screening=expand(f"{S6}/{{protein}}/screening.tsv", protein=PROTEINS),
        conservancy=expand(f"{S7}/{{protein}}/conservancy.tsv", protein=PROTEINS),
        coverage=expand(f"{S8}/{{protein}}/population_coverage.tsv", protein=PROTEINS),
    output:
        table=f"{S9}/final_candidates.tsv",
        fasta=f"{S9}/final_candidates.fasta",
        report=f"{S9}/PIPELINE_REPORT.md",
        attrition=f"{S9}/attrition.tsv",
    params:
        proteins=PROTEINS,
        rank_by=config["final"]["rank_by"],
        require_all=config["final"]["require_all_screens"],
        stage_dirs={
            "S1": S1, "S2": S2, "S3": S3, "S4": S4, "S5": S5,
            "S6": S6, "S7": S7, "S8": S8, "S9": S9,
        },
    log:
        f"{S9}/logs/final.log",
    conda:
        "../envs/core.yaml"
    script:
        "../scripts/09_final.py"
