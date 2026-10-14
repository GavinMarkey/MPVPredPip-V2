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
#
# One job per threshold in the sweep. The thresholds come from
# clustering.threshold ("90:50" -> 90, 80, 70, 60, 50) and are known when the
# DAG is built, so this is a plain expand() rather than a checkpoint.
rule run_cluster:
    input:
        fasta=f"{S5}/{{protein}}/clustering_input.fasta",
    output:
        # temp(), so Snakemake removes it once select_consensus below has
        # consumed it. It is the channel between these two rules rather than a
        # result: nothing outside stage 5 reads it, and everything it carries
        # ends up in clusters.tsv and cluster_members.tsv. Declaring it as an
        # output is still what makes the dependency explicit - dropping it
        # entirely would mean the two rules passing data by convention.
        raw=temp(f"{S5}/{{protein}}/cluster_raw_{{threshold}}.csv"),
    params:
        # A PERCENTAGE, 0-100 - the IEDB tool's scale, not a fraction. The script
        # refuses a value <= 1 rather than silently clustering everything.
        unique=config["clustering"]["unique"],
    log:
        f"{S5}/{{protein}}/logs/cluster_{{threshold}}.log",
    conda:
        "../envs/core.yaml"
    shell:
        "python workflow/scripts/05_run_cluster.py"
        " --fasta '{input.fasta}'"
        " --raw '{output.raw}'"
        " --threshold {wildcards.threshold}"
        " --unique {params.unique}"
        " --log '{log}'"


# Pools the sweep into clusters.tsv and cluster_members.tsv, each carrying every
# threshold with a `threshold` column. Stages 6, 8 and 9 read those two.
#
# They used to be written alongside per-threshold clusters_<t>.tsv and
# cluster_members_<t>.tsv, for reading one threshold at a time. Those were
# dropped: they were row-for-row identical to the pooled tables, which is what a
# `threshold` column already gives you, so a five-threshold sweep wrote ten
# files duplicating these two. Filter on `threshold` instead - e.g.
#     awk -F'\t' 'NR==1 || $2==90' clusters.tsv
#
# passing_consensus.fasta is deduplicated by sequence, because a peptide that
# passes at 90% usually passes again at every looser threshold and would
# otherwise be screened five times over.
rule select_consensus:
    input:
        raw=expand(
            f"{S5}/{{{{protein}}}}/cluster_raw_{{threshold}}.csv",
            threshold=CLUSTER_THRESHOLDS,
        ),
        epitopes=f"{S5}/{{protein}}/all_epitopes.tsv",
    output:
        clusters=f"{S5}/{{protein}}/clusters.tsv",
        members=f"{S5}/{{protein}}/cluster_members.tsv",
        fasta=f"{S5}/{{protein}}/passing_consensus.fasta",
        summary=f"{S5}/{{protein}}/summary.md",
    params:
        thresholds=CLUSTER_THRESHOLDS,
        require=config["clustering"]["require"],
        require_both=config["clustering"]["require_both_bcell_tools"],
        reject_ambiguous=config["clustering"]["reject_ambiguous_consensus"],
    log:
        f"{S5}/{{protein}}/logs/select.log",
    conda:
        "../envs/core.yaml"
    script:
        "../scripts/05_select_consensus.py"
