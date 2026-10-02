#!/usr/bin/env python3
"""Preflight check: report which external tools resolve before a run starts.

Run this after installing the predictors and before the first pipeline run:

    python setup/check_tools.py

It never executes a predictor - it only checks that the paths in ``config.yaml``
exist and that a plausible entry point is present inside each, then prints what
it found. That turns "stage 4 failed after forty minutes" into a ten-second
answer, which matters because several of these tools are slow and the failures
otherwise surface deep inside a Snakemake job.

Exit status is 0 if every *enabled* tool resolved, 1 otherwise.
"""

from __future__ import annotations

import glob
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from vaxpipe import inputs  # noqa: E402

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required: conda install -c conda-forge pyyaml")


GREEN, RED, YELLOW, DIM, RESET = "\033[32m", "\033[31m", "\033[33m", "\033[2m", "\033[0m"
if os.name == "nt" and not os.environ.get("WT_SESSION"):
    GREEN = RED = YELLOW = DIM = RESET = ""

OK, MISSING, SKIPPED = f"{GREEN}ok{RESET}", f"{RED}MISSING{RESET}", f"{YELLOW}skipped{RESET}"
EMPTY = f"{YELLOW}empty{RESET}"


def entry_points(root: str, names: list[str]) -> list[str]:
    """Find any of *names* anywhere under *root*, honouring their order.

    The first pattern that matches anything wins, and only its matches are
    returned. Callers pass names most-specific-first - ``["IApred.py", "*.py"]``
    - so that a precise name beats the catch-all.

    This used to merge every pattern's matches and ``sorted()`` the lot, which
    silently discarded that priority: the alphabetical winner for IAPred was
    ``IApred-training/10fold_CV.py``, a training script, because ``-`` sorts
    before ``.``. The tool then reported a working install while naming an
    entry point that would not run - the exact class of false green this
    script exists to prevent.
    """
    for name in names:
        found = set(glob.glob(os.path.join(root, name)))
        found |= set(glob.glob(os.path.join(root, "**", name), recursive=True))
        if found:
            return sorted(found)
    return []


def check(label: str, path: str, names: list[str], enabled: bool = True, note: str = "") -> bool:
    if not enabled:
        print(f"  {SKIPPED:<18} {label:<26} {DIM}disabled in config.yaml{RESET}")
        return True

    if not os.path.exists(path):
        print(f"  {MISSING:<18} {label:<26} {path}")
        if note:
            print(f"  {'':<18} {DIM}{note}{RESET}")
        return False

    hits = entry_points(path, names) if names else [path]
    if not hits:
        print(f"  {MISSING:<18} {label:<26} {path}")
        print(f"  {'':<18} {DIM}directory exists but no entry point matching {names}{RESET}")
        return False

    print(f"  {OK:<18} {label:<26} {os.path.relpath(hits[0])}")
    return True


def check_algpred_model(tool_dir: str, enabled: bool, model: int) -> bool:
    """Separate check: AlgPred's entry point existing does not mean it can run.

    ``rf_model`` is Git LFS-tracked, so a plain clone leaves it absent or as a
    pointer stub, and the tool then fails inside joblib. Reporting AlgPred as
    "ok" on the strength of algpred2.py alone is exactly the kind of false
    green this script exists to avoid.
    """
    if not enabled:
        print(f"  {SKIPPED:<18} {'AlgPred rf_model':<26} {DIM}disabled in config.yaml{RESET}")
        return True

    path = os.path.join(tool_dir, "rf_model")
    advice = f"Fetch it with:  git -C {tool_dir} lfs pull   (or re-run setup/install_algpred.sh)"

    if not os.path.exists(path):
        print(f"  {MISSING:<18} {'AlgPred rf_model':<26} {path}")
        print(f"  {'':<18} {DIM}LFS-tracked; a plain clone does not fetch it. {advice}{RESET}")
        return False

    with open(path, "rb") as handle:
        if handle.read(23) == b"version https://git-lfs":
            print(f"  {MISSING:<18} {'AlgPred rf_model':<26} {path}")
            print(f"  {'':<18} {DIM}this is an LFS pointer, not the model. {advice}{RESET}")
            return False

    size = os.path.getsize(path)
    print(f"  {OK:<18} {'AlgPred rf_model':<26} {os.path.relpath(path)} ({size} bytes)")

    if model == 2:
        for label, relative in (
            ("MERCI", os.path.join("progs", "MERCI_motif_locator.pl")),
            ("BLAST database", os.path.join("Database", "data.phr")),
            ("IgE motifs", os.path.join("Database", "pos_ige_motifs.txt")),
        ):
            if not os.path.exists(os.path.join(tool_dir, relative)):
                print(f"  {MISSING:<18} {('AlgPred ' + label):<26} {relative}")
                return False
        blastp = shutil.which("blastp")
        if blastp:
            print(f"  {OK:<18} {'blastp (hybrid model)':<26} {blastp}")
        else:
            # Not a failure here: blastp comes from workflow/envs/screening.yaml,
            # which Snakemake builds at run time and which is not active now.
            print(f"  {OK:<18} {'blastp (hybrid model)':<26} "
                  f"{DIM}not on PATH here; provided by envs/screening.yaml{RESET}")
    return True


def main() -> int:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)

    config_path = os.path.join(root, "config.yaml")
    if not os.path.exists(config_path):
        sys.exit(f"config.yaml not found at {config_path}")
    with open(config_path, encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    tools = config["tools"]
    screening = config["screening"]
    results: list[bool] = []

    print("\nExternal tool check\n" + "=" * 70)

    print("\nStage 2 / 6 - antigenicity")
    backend = config["antigenicity"]["tool"]
    if backend == "iapred":
        # IApred.py first: a bare "*.py" matches whichever file sorts first,
        # which is a training script, and reporting that as the entry point is
        # misleading when someone is checking their install.
        ok = check("IAPred", tools["antigenicity"], ["IApred.py", "*.py"],
                   note="Install with: bash setup/install_iapred.sh")
        results.append(ok)
        if ok:
            print(
                f"  {'':<18} {DIM}Threshold is {config['antigenicity']['threshold']}; IApred's "
                f"own 'High' boundary is 0.3 on a{RESET}"
            )
            print(
                f"  {'':<18} {DIM}roughly -3 to 3 scale, so negative scores are normal.{RESET}"
            )
    elif backend == "kolaskar":
        print(f"  {OK:<18} {'built-in (kolaskar)':<26} {DIM}no external tool needed{RESET}")
        print(
            f"  {'':<18} {YELLOW}NOTE{RESET} this is a propensity-scale stand-in, not "
            "IAPred or VaxiJen."
        )
        print(f"  {'':<18} {DIM}See tools/README.md section 1 before reporting results.{RESET}")
    else:
        command = config["antigenicity"].get("command", "")
        if command:
            print(f"  {OK:<18} {'external command':<26} {DIM}{command}{RESET}")
        else:
            print(f"  {MISSING:<18} {'external command':<26} antigenicity.command is empty")
            results.append(False)

    print("\nStage 3 - B-cell prediction")
    results.append(check("BepiPred", tools["bepipred"],
                         ["bepipred3_CLI.py", "BepiPred3_CLI.py", "*.py"],
                         note="Install with: bash setup/install_bepipred.sh"))
    # epidope.yml, not .yaml - the .yaml beside it is a superseded stub. The .yml
    # is the upstream project's own conda env export.
    print(f"  {DIM}EpiDope is installed by workflow/envs/epidope.yml at run time.{RESET}")

    print("\nStage 4 - T-cell prediction")
    iedb_note = "Install with: bash setup/install_iedb.sh"
    results.append(check("IEDB MHC-I", tools["iedb_mhc_i"], ["predict_binding.py"],
                         note=iedb_note))
    results.append(check("IEDB MHC-II", tools["iedb_mhc_ii"],
                         ["mhc_II_binding.py", "predict_binding.py"],
                         note=iedb_note))
    results.append(check("IEDB immunogenicity", tools["iedb_immunogenicity"],
                         ["predict_immunogenicity.py"],
                         enabled=config["tcell"]["immunogenicity"]["enabled"],
                         note=iedb_note))

    print("\nStage 5 - clustering")
    # Not a check: nothing here is required any more. The IEDB cluster standalone
    # only groups peptides - it emits no consensus, which is what stage 5 exists
    # to produce - and it is Python 2 source. vaxpipe.cluster reimplements its
    # identity rule and builds the consensus, so stage 5 needs no external tool.
    print(f"  {OK:<18} {'clustering':<26} {DIM}in-process (vaxpipe.cluster){RESET}")
    print(f"  {'':<18} {DIM}The IEDB cluster standalone is not used; see tools/README.md.{RESET}")

    print("\nStage 6 - screening")
    results.append(check("AlgPred 2.0", tools["algpred"], ["algpred*.py"],
                         enabled=screening["allergenicity"]["enabled"],
                         note="Install with: bash setup/install_algpred.sh"))
    results.append(check_algpred_model(
        tools["algpred"],
        enabled=screening["allergenicity"]["enabled"],
        model=int(screening["allergenicity"].get("model", 2)),
    ))
    # Not a check, for the same reason as clustering above: the PIR Peptide
    # Match standalone publishes no built jar, and the exact k-mer lookup it
    # would perform is done in-process by vaxpipe.autoimmunity.
    print(f"  {OK:<18} {'autoimmunity':<26} {DIM}in-process (vaxpipe.autoimmunity){RESET}")
    print(f"  {'':<18} {DIM}PIR Peptide Match is not used; see tools/README.md.{RESET}")
    results.append(check("Human proteome", screening["autoimmunity"]["human_proteome"], [],
                         enabled=screening["autoimmunity"]["enabled"],
                         note="Fetch with: bash setup/fetch_human_proteome.sh"))

    print("\nStage 7 - conservation")
    results.append(check("ConSurf", tools["consurf"], [],
                         enabled=config["conservation"]["consurf"]["enabled"]))

    print("\nStage 8 - population coverage")
    results.append(check("IEDB population coverage", tools["iedb_population_coverage"],
                         ["calculate_population_coverage.py"],
                         note=iedb_note))

    print("\nAllele panels")
    for label, key in (("MHC class I", "mhc_i"), ("MHC class II", "mhc_ii")):
        path = config["tcell"][key]["allele_file"]
        if os.path.exists(path):
            with open(path, encoding="utf-8-sig") as handle:
                count = sum(1 for line in handle if line.split("#")[0].strip())
            print(f"  {OK:<18} {label:<26} {count} alleles")
        else:
            print(f"  {MISSING:<18} {label:<26} {path}")
            results.append(False)

    print("\nInputs")
    input_dir = config["sequences"]["input_dir"]
    extensions = list(config["sequences"]["extensions"])
    if os.path.isdir(input_dir):
        files = inputs.find_inputs(input_dir, extensions)
        if files:
            print(f"  {OK:<18} {'sequence directory':<26} {len(files)} FASTA file(s)")
            for path in files[:10]:
                print(f"  {'':<18} {DIM}{os.path.basename(path)}{RESET}")
        else:
            # Not a failure - an empty folder is the normal state of a fresh
            # clone - but it is the one thing the pipeline cannot do without.
            print(f"  {EMPTY:<18} {'sequence directory':<26} {input_dir}")
            print(f"  {'':<18} {DIM}Add your sequences - naming rules are in "
                  f"{input_dir}/{inputs.README_NAME}{RESET}")
    else:
        print(f"  {MISSING:<18} {'sequence directory':<26} {input_dir}")
        print(f"  {'':<18} {DIM}Create it, and see the naming rules, with: "
              f"python setup/init_input.py{RESET}")
        results.append(False)

    failures = results.count(False)
    print("\n" + "=" * 70)
    if failures:
        print(
            f"{RED}{failures} tool(s) unresolved.{RESET} The pipeline will fail at the "
            "corresponding stage.\nSee tools/README.md for how to install each one."
        )
        return 1

    print(f"{GREEN}All enabled tools resolved.{RESET}")
    print(
        f"{DIM}Note: this checks that entry points exist, not that their CLI matches "
        f"what the adapters expect.\nRun a single protein first and check the raw output "
        f"in each stage folder.{RESET}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
