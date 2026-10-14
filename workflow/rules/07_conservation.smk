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
#                When conservation.reference_sets.enabled is true, each strain is
#                analysed against its own NCBI reference set (fetched beforehand
#                by ``python setup/fetch_ncbi_conservancy.py``), producing one
#                row per (consensus, strain).  This is the correct mode: a
#                peptide's conservancy is measured against natural variation, not
#                the user's 3 input FASTAs.
#
#   MSA          a multiple alignment of the non-redundant strain sequences,
#                giving the alignment-based context for the conservancy numbers
#                and the input ConSurf needs.
#
# ConSurf is off by default: its standalone is licence-gated and pulls in a large
# dependency stack. Enable it with conservation.consurf.enabled once installed.

def _strains_for_protein(protein: str) -> list[tuple[str, str]]:
    """Return ``[(strain, slug)]`` for *protein*, from the parse-time input scan.

    STRAIN_PAIRS is built in the Snakefile by scanning the input files, which are
    static. That matters here: this function feeds both the ``input`` and the
    ``params`` of ``conservancy``, and Snakemake evaluates those at different
    points in a run.

    It used to read the stage 1 manifest and return ``[]`` when that file was
    absent. On a full run the manifest does not exist while the DAG is being
    built but does by the time ``params`` is evaluated, so ``input`` saw no
    reference FASTAs while ``params`` saw all of them. Stage 7 was therefore
    scheduled without them among its inputs - nothing told Snakemake to build
    them - and then died inside the script on the first one it opened:

        FileNotFoundError: 'resources/conservancy/Spike/Bundibugyo_virus.fasta'

    Reading static input files instead gives the same answer at both moments, so
    the two can no longer disagree.
    """
    return [(strain, slug) for label, strain, slug in STRAIN_PAIRS if label == protein]


def _reference_fastas_for(wildcards):
    """Return the list of reference FASTA paths this protein needs.

    When reference_sets is enabled, the FASTAs become rule inputs so a
    missing fetch fails loudly at DAG time rather than silently falling
    back to the 3 manifest sequences.
    """
    ref_cfg = config.get("conservation", {}).get("reference_sets", {})
    if not ref_cfg.get("enabled"):
        return []
    pairs = _strains_for_protein(wildcards.protein)
    return [
        f"resources/conservancy/{wildcards.protein}/{slug}.fasta"
        for _, slug in pairs
    ]


def _reference_fastas_dict(wildcards):
    """Return a ``{strain: path}`` dict for the script, or None."""
    ref_cfg = config.get("conservation", {}).get("reference_sets", {})
    if not ref_cfg.get("enabled"):
        return None
    pairs = _strains_for_protein(wildcards.protein)
    if not pairs:
        return None
    return {
        strain: f"resources/conservancy/{wildcards.protein}/{slug}.fasta"
        for strain, slug in pairs
    }


# Provisions the NCBI reference sets that `conservancy` below reads, so a
# checkout without them builds them rather than failing on an input that nothing
# produces. That failure is not hypothetical: the FASTAs are declared inputs, no
# rule made them, and deleting them to refresh the sets broke stage 7.
#
# This does not reintroduce the drift problem config.yaml warns about.
# fetch_ncbi_conservancy.py never overwrites a FASTA that is already there - it
# prints "already present ... delete it to re-fetch" and moves on - so this rule
# provisions on first use and is a no-op afterwards. Refreshing a set stays a
# deliberate delete-then-refetch, which is what keeps two runs comparable.
#
# One job rather than one per protein: the script reads the whole manifest and
# fetches every pair in a single pass, so a per-protein rule would repeat the
# same fetch and race on the same output paths.
#
# Only the FASTAs are outputs. The script writes a .provenance.json only for the
# pairs it actually fetched, so listing those here would leave the rule short an
# output on any run where some sets were already present.
#
# Defined only when there is something to fetch - a rule with an empty `output`
# is an error, and REFERENCE_SET_PATHS is empty when reference sets are off.
if REFERENCE_SET_PATHS:

    rule fetch_conservancy_reference_sets:
        input:
            manifest=f"{S1}/manifest.tsv",
        output:
            fastas=list(REFERENCE_SET_PATHS.values()),
        log:
            f"{S7}/logs/fetch_reference_sets.log",
        conda:
            "../envs/core.yaml"
        shell:
            "python setup/fetch_ncbi_conservancy.py > '{log}' 2>&1"


rule conservancy:
    input:
        consensus=f"{S6}/{{protein}}/passed.fasta",
        manifest=f"{S1}/manifest.tsv",
        reference_fastas=_reference_fastas_for,
    output:
        table=f"{S7}/{{protein}}/conservancy.tsv",
        targets=f"{S7}/{{protein}}/iedb_upload/non_redundant_sequences.fasta",
        epitopes=f"{S7}/{{protein}}/iedb_upload/epitope_list.txt",
        summary=f"{S7}/{{protein}}/summary.md",
    params:
        threshold=config["conservation"]["conservancy"]["identity_threshold"],
        deduplicate=config["conservation"]["deduplicate_inputs"],
        mode=config["conservation"]["conservancy"]["mode"],
        reference_fastas=_reference_fastas_dict,
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
