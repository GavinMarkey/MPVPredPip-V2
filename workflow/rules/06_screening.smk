# Stage 6 - antigenicity, allergenicity and autoimmunity screening.
#
# Applied to the consensus peptides that passed stage 5. Unlike stage 2, these
# are hard filters: a peptide must be predicted antigenic, must not be predicted
# allergenic, and must not share an exact window of
# screening.autoimmunity.min_exact_match_length residues with the human proteome.
#
# The autoimmunity check is done in-process by mpvpredpip.autoimmunity rather than
# by the PIR Peptide Match standalone. That tool is free but publishes no built
# jar, and the work it would do here - exact k-mer lookup over a 20,000-protein
# FASTA for about a dozen peptides - takes roughly a second in Python. See the
# module docstring for the full reasoning and the L/I-equivalence caveat.
#
# This is why there is no index-building rule: there is no index. The proteome
# is read straight in as an input, so Snakemake reruns the screen if it changes.


# AlgPred runs in its own rule and its own environment. Its rf_model was
# pickled with scikit-learn 0.22 and will not unpickle under the 1.5.2 that
# IApred requires, so the two cannot share one environment - see
# workflow/envs/algpred.yaml. `shell:` rather than `script:` because Snakemake
# cannot be imported under that environment's Python 3.8, exactly as for
# EpiDope in stage 3.
rule algpred_run:
    input:
        fasta=f"{S5}/{{protein}}/passing_consensus.fasta",
    output:
        raw=f"{S6}/{{protein}}/raw/allergenicity.csv",
    params:
        tool_dir=config["tools"]["algpred"],
        threshold=config["screening"]["allergenicity"]["threshold"],
        model=config["screening"]["allergenicity"].get("model", 2),
    log:
        f"{S6}/{{protein}}/logs/allergenicity.log",
    conda:
        "../envs/algpred.yaml"
    shell:
        "python workflow/scripts/06_allergenicity.py"
        " --fasta '{input.fasta}'"
        " --out '{output.raw}'"
        " --tool-dir '{params.tool_dir}'"
        " --threshold {params.threshold}"
        " --model {params.model}"
        " --log '{log}'"


def _allergenicity_input(wildcards):
    """AlgPred's output, or nothing when the screen is switched off.

    Returned as a list so that disabling the screen removes the dependency
    entirely rather than demanding a file no rule will be asked to build.
    """
    if not config["screening"]["allergenicity"].get("enabled", True):
        return []
    return f"{S6}/{wildcards.protein}/raw/allergenicity.csv"


rule screen_consensus:
    input:
        fasta=f"{S5}/{{protein}}/passing_consensus.fasta",
        clusters=f"{S5}/{{protein}}/clusters.tsv",
        proteome=config["screening"]["autoimmunity"]["human_proteome"],
        allergenicity=_allergenicity_input,
    output:
        table=f"{S6}/{{protein}}/screening.tsv",
        passed=f"{S6}/{{protein}}/passed.fasta",
        summary=f"{S6}/{{protein}}/summary.md",
    params:
        antigenicity=config["screening"]["antigenicity"],
        antigenicity_tool=config["antigenicity"]["tool"],
        antigenicity_command=config["antigenicity"].get("command", ""),
        antigenicity_dir=config["tools"]["antigenicity"],
        allergenicity=config["screening"]["allergenicity"],
        autoimmunity=config["screening"]["autoimmunity"],
        raw_dir=f"{S6}/{{protein}}/raw",
    log:
        f"{S6}/{{protein}}/logs/screening.log",
    conda:
        "../envs/screening.yaml"
    script:
        "../scripts/06_screen.py"
