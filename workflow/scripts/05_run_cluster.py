"""Stage 5 - cluster the pooled epitopes and build each cluster's consensus.

This no longer runs the IEDB cluster standalone. That tool only groups peptides -
it emits no consensus sequence and no alignment, is Python 2 source, and its CLI
is not the one this script used to build. The clustering and the consensus are
computed by vaxpipe.cluster instead, from the identity rule reimplemented from
the standalone's own source. See the module docstring of src/vaxpipe/cluster.py.

The output file is unchanged: the same layout the IEDB *web* tool produces, which
05_select_consensus.py reads. Only its producer has changed.

Invoked through ``shell:``, not ``script:`` - see 03_bcell_predict.py.
"""

import argparse

import _common

_common.bootstrap()

from vaxpipe import cluster as cluster_mod  # noqa: E402
from vaxpipe import fasta as fasta_mod  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fasta", required=True)
parser.add_argument("--raw", required=True, help="cluster output to write")
parser.add_argument(
    "--threshold",
    type=float,
    required=True,
    help="percent sequence identity, 0-100 (the IEDB tool's own scale; its default is 80)",
)
parser.add_argument(
    "--unique",
    choices=["True", "False"],
    default="True",
    help="collapse identical peptides before clustering, as the tool's --unique=on does",
)
parser.add_argument("--log", default=None)
args = parser.parse_args()

with _common.tee_log(args.log):
    # The threshold used to be documented as a 0-1 fraction, because the pipeline
    # was written against the IEDB web tool. The scale here is the standalone's:
    # 0-100. Passing the old 0.6 would put every peptide in the run into one
    # cluster and quietly produce a meaningless consensus, so refuse it outright
    # rather than let it through.
    if args.threshold <= 1.0:
        raise SystemExit(
            f"clustering.threshold is {args.threshold}, which looks like a fraction.\n"
            "It is a PERCENTAGE, 0-100 - the IEDB cluster tool's own scale, where the\n"
            "default is 80. A value at or below 1 would cluster every epitope together.\n"
            "Set clustering.threshold in config.yaml to e.g. 60, not 0.6."
        )

    records = fasta_mod.read_fasta(args.fasta)
    if not records:
        raise SystemExit(f"No sequences in {args.fasta}")

    entries = [(record.id, record.sequence) for record in records]
    clusters = cluster_mod.cluster_peptides(
        entries,
        threshold=args.threshold,
        unique=args.unique == "True",
    )
    cluster_mod.write_clusters(args.raw, clusters)

    singletons = sum(1 for cluster in clusters if len(cluster.members) == 1)
    ambiguous = sum(1 for cluster in clusters if cluster_mod.AMBIGUOUS in cluster.consensus)
    voted = [cluster_mod.disagreeing_positions(cluster) for cluster in clusters]
    needed_a_vote = sum(1 for count in voted if count)

    print(
        f"Clustered {len(entries)} peptide(s) into {len(clusters)} cluster(s) "
        f"at {args.threshold:.0f}% identity."
    )
    print(f"  {singletons} singleton(s)")
    # Each column is a majority vote, so a difference between strains no longer
    # shows up in the consensus itself. Report it here, because a consensus that
    # needed many votes represents its members less faithfully than one that
    # needed none - stage 7's conservancy analysis is what quantifies that.
    print(f"  {needed_a_vote} cluster(s) whose members disagree somewhere")
    if voted:
        print(f"  {max(voted)} disagreeing position(s) in the worst cluster")
    # X now means only that a column was a genuine tie, which is rare.
    print(f"  {ambiguous} cluster(s) left with an ambiguous residue after voting")
    print(f"Cluster output written to {args.raw}")
