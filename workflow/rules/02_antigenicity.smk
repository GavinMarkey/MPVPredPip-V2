# Stage 2 - whole-protein antigenicity.
#
# Informational by default: every input protein is scored and reported, but
# nothing is dropped unless antigenicity.drop_failing_proteins is set. The same
# predictor is re-applied as a hard filter to the final consensus peptides in
# stage 6, which is where antigenicity actually gates the candidate list.


rule whole_protein_antigenicity:
    input:
        manifest=f"{S1}/manifest.tsv",
        proteins=expand(f"{NORMALISED}/{{protein}}.fasta", protein=PROTEINS),
    output:
        table=f"{S2}/whole_protein_antigenicity.tsv",
        summary=f"{S2}/summary.md",
    params:
        tool=config["antigenicity"]["tool"],
        command=config["antigenicity"].get("command", ""),
        tool_dir=config["tools"]["antigenicity"],
        threshold=config["antigenicity"]["threshold"],
        drop_failing=config["antigenicity"].get("drop_failing_proteins", False),
        raw_dir=f"{S2}/raw",
    log:
        f"{S2}/logs/antigenicity.log",
    conda:
        # Not core.yaml: this rule runs IApred, whose dependencies live in the
        # screening environment alongside the stage 6 antigenicity filter that
        # uses the same predictor.
        "../envs/screening.yaml"
    script:
        "../scripts/02_antigenicity.py"
