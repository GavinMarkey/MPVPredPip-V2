"""Stage 6 - run AlgPred 2.0 and keep its raw CSV.

Driven by `shell:`, not `script:`, and run under workflow/envs/algpred.yaml,
which pins Python 3.8 and scikit-learn 0.22 so that AlgPred's rf_model will
unpickle. Snakemake 9 cannot be imported at that Python version, which is why
this takes its arguments on the command line - the same arrangement stage 3
uses for EpiDope.

Only mpvpredpip.algpred is imported here, and that module is standard-library-only
for this reason. Parsing the CSV happens afterwards, in the ordinary screening
environment.
"""

import argparse
import os
import sys


def main():
    parser = argparse.ArgumentParser(description="Run AlgPred 2.0 for one protein.")
    parser.add_argument("--fasta", required=True, help="consensus peptides to screen")
    parser.add_argument("--out", required=True, help="where to write AlgPred's CSV")
    parser.add_argument("--tool-dir", required=True)
    parser.add_argument("--threshold", type=float, default=0.3)
    parser.add_argument("--model", type=int, choices=[1, 2], default=2)
    parser.add_argument("--log", default="")
    args = parser.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(here)), "src"))
    from mpvpredpip import algpred  # noqa: E402

    handle = None
    if args.log:
        log_dir = os.path.dirname(os.path.abspath(args.log))
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)
        handle = open(args.log, "w")
        sys.stdout = handle
        sys.stderr = handle

    try:
        # An empty candidate list is a legitimate outcome of stage 5, not a
        # failure. AlgPred would reject the empty file, so short-circuit and
        # leave a header-only CSV for the parser to find.
        if os.path.getsize(args.fasta) == 0:
            out_dir = os.path.dirname(os.path.abspath(args.out))
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            with open(args.out, "w") as out:
                out.write("Subject,ML Score,MERCI Score,BLAST Score,Hybrid Score,Prediction\n")
            print("No consensus peptides to screen; wrote an empty AlgPred table.")
            return 0

        algpred.run(
            fasta_path=args.fasta,
            out_path=args.out,
            tool_dir=args.tool_dir,
            threshold=args.threshold,
            model=args.model,
        )
        print("AlgPred wrote {0}".format(args.out))
        return 0
    except algpred.AlgPredError as error:
        print("ERROR: {0}".format(error), file=sys.stderr)
        return 1
    finally:
        if handle is not None:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            handle.close()


if __name__ == "__main__":
    raise SystemExit(main())
