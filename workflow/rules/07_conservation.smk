# Stage 7 - conservation analysis.
#
# Two independent views of conservation:
#
#   conservancy  the IEDB conservancy calculation, reimplemented locally so the
#                pipeline runs unattended. Redundant input sequences are collapsed
#                first, so identical strains submitted twice cannot inflate the
#                percentages. With conservation.conservancy.mode set to "both",
#                upload-ready files for the IEDB web tool are written as well, to
#                allow a manual cross-check without reformatting anything.
#
#   MSA          a multiple alignment of the non-redundant strain sequences,
#                giving the alignment-based context for the conservancy numbers
#                and the input ConSurf needs.
#
# ConSurf is off by default: its standalone is licence-gated and pulls in a large
# dependency stack. Enable it with conservation.consurf.enabled once installed.


rule conservancy:
    input:
        consensus=f"{S6}/{{protein}}/passed.fasta",
        manifest=f"{S1}/manifest.tsv",
    output:
        table=f"{S7}/{{protein}}/conservancy.tsv",
        targets=f"{S7}/{{protein}}/iedb_upload/non_redundant_sequences.fasta",
        epitopes=f"{S7}/{{protein}}/iedb_upload/epitope_list.txt",
        summary=f"{S7}/{{protein}}/summary.md",
    params:
        threshold=config["conservation"]["conservancy"]["identity_threshold"],
        deduplicate=config["conservation"]["deduplicate_inputs"],
        mode=config["conservation"]["conservancy"]["mode"],
    log:
        f"{S7}/{{protein}}/logs/conservancy.log",
    conda:
        "../envs/core.yaml"
    script:
        "../scripts/07_conservancy.py"


rule msa:
    input:
        fasta=f"{S7}/{{protein}}/iedb_upload/non_redundant_sequences.fasta",
    output:
        alignment=f"{S7}/{{protein}}/MSA/{{protein}}_aligned.fasta",
    params:
        args=config["conservation"]["msa"]["args"],
    threads: config["resources"]["threads_per_job"]
    log:
        f"{S7}/{{protein}}/logs/msa.log",
    conda:
        "../envs/conservation.yaml"
    shell:
        "mafft {params.args} --thread {threads} '{input.fasta}' > '{output.alignment}' 2> '{log}'"
