# Stage 8 - population coverage.
#
# For each consensus peptide that survived stage 6, the MHC class I and class II
# alleles are collected from the T-cell epitopes that fall inside that consensus
# (traced through the stage 5 cluster membership table), and the IEDB population
# coverage tool is run over them.
#
# Taking the alleles from cluster membership rather than re-predicting against
# the consensus peptide matters: the coverage figure then describes the alleles
# that actually drove the cluster's selection, which is what the candidate is
# built from.


rule population_coverage:
    input:
        screening=f"{S6}/{{protein}}/screening.tsv",
        members=f"{S5}/{{protein}}/cluster_members.tsv",
        epitopes=f"{S5}/{{protein}}/all_epitopes.tsv",
    output:
        table=f"{S8}/{{protein}}/population_coverage.tsv",
        input_file=f"{S8}/{{protein}}/coverage_input.tsv",
        raw=f"{S8}/{{protein}}/raw_population_coverage.txt",
        summary=f"{S8}/{{protein}}/summary.md",
    params:
        tool_dir=config["tools"]["iedb_population_coverage"],
        populations=",".join(config["population_coverage"]["populations"]),
        mhc_class=config["population_coverage"]["mhc_class"],
    log:
        f"{S8}/{{protein}}/logs/population_coverage.log",
    conda:
        "../envs/iedb.yaml"
    # shell:, not script: - the IEDB environment is Python 3.8 and cannot import
    # Snakemake 9. See 03_bcell_predict.py.
    shell:
        "python workflow/scripts/08_population_coverage.py"
        " --protein '{wildcards.protein}'"
        " --screening '{input.screening}'"
        " --members '{input.members}'"
        " --epitopes '{input.epitopes}'"
        " --table '{output.table}'"
        " --input-file '{output.input_file}'"
        " --raw '{output.raw}'"
        " --summary '{output.summary}'"
        " --tool-dir '{params.tool_dir}'"
        " --populations '{params.populations}'"
        " --mhc-class '{params.mhc_class}'"
        " --log '{log}'"
