# Stage 5 - epitope clustering and B/T consensus selection.
#
# B-cell (BepiPred + EpiDope) and T-cell (immunogenic MHC-I + MHC-II) epitopes
# for one protein are pooled, deduplicated and clustered with the IEDB epitope
# cluster analysis tool at clustering.threshold.
#
# A cluster is carried forward only if it contains at least one MHC-I epitope,
# one MHC-II epitope and one B-cell epitope from BepiPred *and* one from EpiDope
# (relax with clustering.require_both_bcell_tools: false). Membership is resolved
# by joining the cluster output against the stage 3/4 epitope tables, so the rule
# is enforced on recorded provenance rather than on parsing identifiers by eye.


rule clustering_input:
    input:
        bcell=f"{S3}/{{protein}}/bcell_epitopes.tsv",
        tcell=f"{S4}/{{protein}}/tcell_epitopes.tsv",
    output:
        fasta=f"{S5}/{{protein}}/clustering_input.fasta",
        table=f"{S5}/{{protein}}/clustering_input.tsv",
        combined=f"{S5}/{{protein}}/all_epitopes.tsv",
    log:
        f"{S5}/{{protein}}/logs/clustering_input.log",
    conda:
        "../envs/core.yaml"
    script:
        "../scripts/05_clustering_input.py"


# Clustering is done in-process by mpvpredpip.cluster, not by the IEDB cluster
# standalone. The standalone only groups peptides - it emits no consensus and no
# alignment, and stage 5 exists to produce a consensus. Its identity rule is
# reimplemented faithfully; see src/mpvpredpip/cluster.py. This rule therefore needs
# no tool directory and runs in the orchestrator's own environment.
rule run_cluster:
    input:
        fasta=f"{S5}/{{protein}}/clustering_input.fasta",
    output:
        raw=f"{S5}/{{protein}}/cluster_raw.csv",
    params:
        # A PERCENTAGE, 0-100 - the IEDB tool's scale, not a fraction. The script
        # refuses a value <= 1 rather than silently clustering everything.
        threshold=config["clustering"]["threshold"],
        unique=config["clustering"]["unique"],
    log:
        f"{S5}/{{protein}}/logs/cluster.log",
    conda:
        "../envs/core.yaml"
    shell:
        "python workflow/scripts/05_run_cluster.py"
        " --fasta '{input.fasta}'"
        " --raw '{output.raw}'"
        " --threshold {params.threshold}"
        " --unique {params.unique}"
        " --log '{log}'"


rule select_consensus:
    input:
        raw=f"{S5}/{{protein}}/cluster_raw.csv",
        epitopes=f"{S5}/{{protein}}/all_epitopes.tsv",
    output:
        clusters=f"{S5}/{{protein}}/clusters.tsv",
        members=f"{S5}/{{protein}}/cluster_members.tsv",
        fasta=f"{S5}/{{protein}}/passing_consensus.fasta",
        summary=f"{S5}/{{protein}}/summary.md",
    params:
        require=config["clustering"]["require"],
        require_both=config["clustering"]["require_both_bcell_tools"],
        reject_ambiguous=config["clustering"]["reject_ambiguous_consensus"],
    log:
        f"{S5}/{{protein}}/logs/select.log",
    conda:
        "../envs/core.yaml"
    script:
        "../scripts/05_select_consensus.py"
