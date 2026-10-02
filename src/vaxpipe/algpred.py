"""Driving AlgPred 2.0 - everything that must happen before it will start.

This module is imported inside ``workflow/envs/algpred.yaml``, which pins
**Python 3.8 and scikit-learn 0.22**. It therefore uses the standard library
only, imports nothing else from :mod:`vaxpipe`, and must stay 3.8-compatible.
Parsing AlgPred's output is deliberately *not* here - that happens afterwards
in the ordinary screening environment, in
:func:`vaxpipe.adapters.screening.parse_allergenicity`.

------------------------------------------------------------------------------
Why AlgPred needs an environment of its own
------------------------------------------------------------------------------
``rf_model`` was pickled with scikit-learn 0.22 or older. Loading it under a
modern release fails in ``sklearn.tree._tree``::

    ValueError: node array from the pickle has an incompatible dtype
    - expected: [... 'missing_go_to_left'] itemsize 64
    - got:      [... 'weighted_n_node_samples']

``missing_go_to_left`` was added to the tree node dtype in scikit-learn 1.3, so
the model predates it. The pin cannot simply be loosened in
``envs/screening.yaml``, because IApred's models require 1.5.2 and the two would
then be mutually exclusive - a collision that file's own comment predicted. So
AlgPred runs in a separate rule under a separate environment, exactly as EpiDope
does in stage 3.

Switching ``screening.allergenicity.model`` to 1 is not an escape: model 1 loads
the same ``rf_model``.

------------------------------------------------------------------------------
Three further things the tool's state forces on the caller
------------------------------------------------------------------------------
* ``rf_model`` is **Git LFS-tracked**, so a plain clone leaves it absent or as a
  pointer stub; joblib then fails with an error that never mentions LFS.
* ``envfile`` ships pointing at the author's own machine, and AlgPred reads it
  from the current directory.
* ``algpred2.py`` reads and writes about ten files by bare relative name, so two
  proteins screened at once from one directory would overwrite each other's
  intermediates - producing wrong scores rather than a crash.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile

#: Files ``algpred2.py`` reads or writes by bare relative name, i.e. resolved
#: against the current directory. The reason each run needs a private one.
SCRATCH = (
    "seq.aac", "seq.pred", "seq.out", "Sequence_1", "RES_1_6_6.out",
    "merci.txt", "merci_output.csv", "merci_hybrid.csv", "blast_hybrid.csv",
    "final_output",
)

#: Candidate entry-point names, most specific first.
ENTRYPOINTS = ("algpred2.py", "algpred_2.py", "algpred.py", "python_scripts/algpred2.py")

#: First line of a Git LFS pointer file.
LFS_POINTER_PREFIX = b"version https://git-lfs"

#: Source fixes applied to a *copy* of the entry script, as
#: ``(name, old, new)``. The installation is never modified.
#:
#: ``joblib-import``
#:     ``sklearn.externals.joblib`` exists in the pinned 0.22 but emits a
#:     deprecation warning, and rewriting it to the standalone package keeps
#:     the script working if the pin is ever raised. ``joblib`` is a declared
#:     dependency of the environment either way.
#:
#: ``single-sequence``
#:     An upstream bug. ``np.loadtxt`` returns a 1-D array when the file holds
#:     a single row and a 2-D one otherwise, so scoring exactly one peptide
#:     reaches ``predict_proba`` with shape ``(20,)`` instead of ``(1, 20)``::
#:
#:         ValueError: Expected 2D array, got 1D array instead
#:
#:     AlgPred therefore cannot score a single sequence, in either model. That
#:     is not a rare edge case here - stage 5 routinely yields one consensus
#:     peptide for a protein, which is exactly what happened to the
#:     nucleoprotein. ``atleast_2d`` is a no-op for every other input.
PATCHES = (
    (
        "joblib-import",
        "from sklearn.externals import joblib",
        "import joblib",
    ),
    (
        "single-sequence",
        "data_test = np.loadtxt(file_name, delimiter=',')",
        "data_test = np.atleast_2d(np.loadtxt(file_name, delimiter=','))",
    ),
)


class AlgPredError(RuntimeError):
    """Raised for anything that stops AlgPred running, with how to fix it."""


def find_entrypoint(tool_dir):
    """Absolute path to AlgPred's entry script inside *tool_dir*."""
    for name in ENTRYPOINTS:
        candidate = os.path.join(tool_dir, name)
        if os.path.exists(candidate):
            return candidate
    raise AlgPredError(
        "No AlgPred entry script under {0}.\n"
        "Install it with:  bash setup/install_algpred.sh".format(tool_dir)
    )


def check_model(tool_dir):
    """Fail early and legibly if ``rf_model`` is missing or an unfetched stub."""
    path = os.path.join(tool_dir, "rf_model")
    advice = (
        "\nFetch it with either:\n"
        "    git -C {0} lfs pull\n"
        "  or download the zip and copy the file in:\n"
        "    https://webs.iiitd.edu.in/raghava/algpred2/algpred2_new.zip"
    ).format(tool_dir)

    if not os.path.exists(path):
        raise AlgPredError(
            "AlgPred's rf_model is missing from {0}.\n"
            "It is tracked with Git LFS, so a plain clone does not fetch it."
            .format(tool_dir) + advice
        )
    with open(path, "rb") as handle:
        if handle.read(len(LFS_POINTER_PREFIX)) == LFS_POINTER_PREFIX:
            raise AlgPredError(
                "AlgPred's rf_model at {0} is a Git LFS pointer, not the model.\n"
                "The clone recorded the file but never downloaded its contents."
                .format(path) + advice
            )


def link_or_copy(source, destination):
    """Symlink *source* into the run directory, copying if links are unavailable."""
    try:
        os.symlink(source, destination)
    except (OSError, NotImplementedError):
        if os.path.isdir(source):
            shutil.copytree(source, destination)
        else:
            shutil.copy2(source, destination)


def write_envfile(workdir, tool_dir):
    """Write an ``envfile`` pointing at this machine rather than the author's.

    The parser is positional and splits each line on its first colon
    (``algpred2.py`` lines 292-306), taking exactly four non-comment lines in
    this order. Do not reorder them, and keep colons out of the paths.

    ``blastp`` comes from PATH because it is provided by the conda environment,
    whose location is decided at environment-build time and is not knowable
    when the tool is installed.
    """
    blastp = shutil.which("blastp")
    if not blastp:
        raise AlgPredError(
            "blastp is not on PATH, and AlgPred's hybrid model (-m 2) needs it.\n"
            "It is declared in workflow/envs/algpred.yaml; if you are running "
            "outside that environment, install BLAST or set "
            "screening.allergenicity.model: 1 in config.yaml."
        )
    merci = os.path.join(tool_dir, "progs", "MERCI_motif_locator.pl")
    database = os.path.join(tool_dir, "Database", "data")
    motifs = os.path.join(tool_dir, "Database", "pos_ige_motifs.txt")

    for label, path in (("MERCI", merci), ("IgE motif file", motifs)):
        if not os.path.exists(path):
            raise AlgPredError("AlgPred's {0} is missing: {1}".format(label, path))

    target = os.path.join(workdir, "envfile")
    if os.path.islink(target) or os.path.exists(target):
        os.remove(target)          # replace the symlink to the shipped one
    with open(target, "w") as handle:
        handle.write("#Generated per run by vaxpipe - see vaxpipe.algpred\n")
        handle.write("BLAST:{0}\n".format(blastp))
        handle.write("BLAST database:{0}\n".format(database))
        handle.write("MERCI:{0}\n".format(merci))
        handle.write("MERCI motif file:{0}\n".format(motifs))


def stage_script(source, destination):
    """Write a compatibility-patched copy of *source* at *destination*.

    Returns the names of the patches that applied. A patch matching nothing is
    not fatal on its own - upstream may have fixed the problem - but each
    known defect is re-checked afterwards against the patched text, so a
    rewrite that silently stopped matching is caught here rather than as a
    crash inside AlgPred.
    """
    with open(source) as handle:
        text = handle.read()

    applied = []
    for name, old, new in PATCHES:
        if old in text:
            text = text.replace(old, new)
            applied.append(name)

    if "sklearn.externals" in text:
        raise AlgPredError(
            "{0} still imports from sklearn.externals after patching.\n"
            "Add the required rewrite to PATCHES in src/vaxpipe/algpred.py."
            .format(source)
        )
    # Stated as an outcome, not as "did a replacement match". Searching the
    # patched text for the bare call cannot work - the fix contains it - and
    # counting occurrences cannot work either, because every occurrence is
    # replaced, so a bare one never survives. Neither catches the failure that
    # matters: upstream rewording the line so the patch quietly stops applying.
    # Requiring the fixed form to be present does catch it.
    if "np.atleast_2d(np.loadtxt(" not in text:
        raise AlgPredError(
            "{0} does not force its feature table to two dimensions after "
            "patching.\nAlgPred then cannot score a single sequence, which "
            "stage 5 regularly produces - the nucleoprotein hit this with one "
            "consensus peptide while the spike's eight were fine.\n"
            "Most likely the line was reworded upstream and the "
            "'single-sequence' entry in PATCHES no longer matches; update it "
            "in src/vaxpipe/algpred.py. If upstream has instead fixed this "
            "another way (loadtxt's ndmin=2, say), relax this check."
            .format(source)
        )

    # The run directory is built from symlinks, so writing naively would edit
    # the installed tool through the link.
    parent = os.path.dirname(destination)
    if os.path.islink(parent):
        raise AlgPredError(
            "Refusing to write {0}: its parent {1} is a symlink into the "
            "installation.".format(destination, parent)
        )
    if os.path.islink(destination) or os.path.exists(destination):
        os.remove(destination)
    with open(destination, "w") as handle:
        handle.write(text)
    return applied


def run(fasta_path, out_path, tool_dir, threshold=0.3, model=2, python_executable=None):
    """Run AlgPred 2.0, writing its CSV to *out_path*.

    ``-d 2`` is **load-bearing** and is never configurable: AlgPred's default,
    ``-d 1``, writes out only the peptides it calls allergenic. Stage 6 needs a
    verdict for every candidate, so under the default the peptides that pass
    would simply be absent, and "nothing was allergenic" would look identical
    to "the tool never ran". The caller cross-checks the output against the
    submitted FASTA for exactly this reason.
    """
    python_executable = python_executable or sys.executable
    tool_dir = os.path.abspath(tool_dir)
    if not os.path.isdir(tool_dir):
        raise AlgPredError(
            "AlgPred is not installed at {0}.\n"
            "Install it with:  bash setup/install_algpred.sh".format(tool_dir)
        )

    script = find_entrypoint(tool_dir)
    check_model(tool_dir)

    out_path = os.path.abspath(out_path)
    parent = os.path.dirname(out_path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    workdir = tempfile.mkdtemp(prefix="algpred-")
    try:
        for entry in os.listdir(tool_dir):
            if entry in SCRATCH or entry == ".git":
                continue
            link_or_copy(os.path.join(tool_dir, entry), os.path.join(workdir, entry))

        staged = os.path.join(workdir, os.path.relpath(script, tool_dir))
        applied = stage_script(script, staged)
        print("[algpred] source patches applied: " + (", ".join(applied) or "none"))

        if int(model) == 2:
            write_envfile(workdir, tool_dir)

        command = [
            python_executable, staged,
            "-i", os.path.abspath(fasta_path),
            "-o", out_path,
            "-t", str(threshold),
            "-m", str(model),
            "-d", "2",      # every peptide, not only the allergenic ones
        ]
        print("[algpred] " + " ".join(command))
        result = subprocess.run(
            command, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        if result.stdout:
            print(result.stdout)
        if result.returncode != 0:
            raise AlgPredError(
                "AlgPred exited {0}.\n--- stderr ---\n{1}".format(
                    result.returncode, (result.stderr or "").strip()
                )
            )
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    return out_path
