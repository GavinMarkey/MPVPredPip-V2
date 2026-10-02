"""Stage 0 - the folder users put their sequences in, and the rules for them.

Everything the pipeline tells a user about their input files comes from this one
module: the ``README.txt`` written into the folder, the notice printed when the
container is set up, the error shown when the folder is empty, and the per-file
preflight report. They share a source on purpose.

The rule used to be written out separately in the README, in two error messages,
in a config comment and in a code comment, and the copies had drifted from what
the parser does. The documented filename convention is only a *fallback*, and
not one of the project's own six input files satisfied it - they worked because
NCBI headers already carry the strain.

What is actually enforced, read off :mod:`mpvpredpip.fasta` and pinned by tests
(``TestInputFolder``) so that this text cannot drift again:

* Protein, strain and accession come from the **header line**,
  ``>ACCESSION protein name [strain]``.
* The **filename** is consulted only when the header lacks a protein name or a
  ``[strain]``, and then only if it matches ``ACCESSION protein name - strain``
  with whitespace on *both* sides of the dash.
* Each ``(protein, strain)`` pair may appear once.
* Sequences must be protein, non-empty, and at least ``min_length`` residues.
"""

from __future__ import annotations

import os
from collections import defaultdict
from collections.abc import Iterable, Sequence

from . import fasta

#: Shown in the notice, and parsed through the real parser by the tests - so the
#: documentation cannot describe something the code does not do.
EXAMPLE_DEFLINE = ">NP_066246.1 spike glycoprotein [Zaire ebolavirus]"
EXAMPLE_FILENAME = "NP_066246.1 spike glycoprotein - Zaire ebolavirus.fasta"

#: Name of the notice written into the input folder. Its extension is not an
#: accepted sequence extension, which is what keeps it from being read as one.
README_NAME = "README.txt"

#: What :func:`mpvpredpip.fasta._annotate` falls back to when it finds nothing.
UNKNOWN_PROTEIN = "unknown protein"
UNKNOWN_STRAIN = "unknown strain"


def find_inputs(input_dir: str, extensions: Iterable[str]) -> list[str]:
    """Paths of the sequence files directly inside *input_dir*, sorted.

    Only top-level files with an accepted extension count. Subfolders are
    ignored, and so is ``README.txt``.
    """
    if not os.path.isdir(input_dir):
        return []
    wanted = {extension.lower() for extension in extensions}
    found: list[str] = []
    for entry in sorted(os.listdir(input_dir)):
        path = os.path.join(input_dir, entry)
        if os.path.isfile(path) and os.path.splitext(entry)[1].lower() in wanted:
            found.append(path)
    return found


def find_problems(records: Sequence[fasta.Record], min_length: int = 50) -> list[str]:
    """Everything stage 1 would reject about *records*, one line per problem.

    Shared by the stage 1 rule and the preflight report so the two cannot
    disagree about what is acceptable.
    """
    problems: list[str] = []
    seen: dict[tuple[str, str], str] = {}
    for record in records:
        for problem in fasta.validate(record, min_length=min_length):
            problems.append(f"  {record.source_path} [{record.id}]: {problem}")

        key = (record.protein_label, record.strain)
        if key in seen:
            problems.append(
                f"  duplicate protein/strain pair {key}: already supplied by {seen[key]}"
                f", now also {record.source_path}"
            )
        else:
            seen[key] = record.source_path
    return problems


def short_rules(input_dir: str, extensions: Iterable[str]) -> str:
    """The rule in a few lines, for error messages."""
    return (
        f"Put protein FASTA files directly in {input_dir!r} "
        f"({' '.join(extensions)}), one sequence per file.\n"
        f"  header : >ACCESSION protein name [strain]\n"
        f"           e.g. {EXAMPLE_DEFLINE}\n"
        f"  file   : <ACCESSION> <protein name> - <strain>.fasta\n"
        f"           e.g. {EXAMPLE_FILENAME}\n"
        f"Full rules and a check of your files:  python setup/init_input.py"
    )


def notice(
    input_dir: str,
    extensions: Sequence[str],
    aliases: dict[str, str],
    min_length: int = 50,
) -> str:
    """The full naming rules. Printed at setup and written to ``README.txt``."""
    alias_lines = "\n".join(
        f'       "{alias}"  ->  {label}' for alias, label in aliases.items()
    ) or "       (none configured)"
    accepted = "  ".join(extensions)

    return f"""\
INPUT FILES
===========

Put your protein sequences directly in this folder:  {input_dir}/
Subfolders are ignored.  Accepted extensions:  {accepted}

1. WHAT TO PUT HERE
   One protein sequence per file, downloaded from NCBI Protein in FASTA format.
   Amino acids only - nucleotide FASTA is rejected. Minimum {min_length} residues.

2. THE HEADER LINE IS WHAT THE PIPELINE READS
       >ACCESSION protein name [strain]
   for example
       {EXAMPLE_DEFLINE}

   accession  =  the first word
   protein    =  the text before the [ ... ]
   strain     =  the text inside the [ ... ]

   NCBI normally puts the SPECIES (its organism name) in the brackets, so by
   default each species counts as one "strain". Files downloaded from NCBI
   already look like this, so you rarely need to edit them.

3. NAME THE FILE
       <ACCESSION> <protein name> - <strain>.fasta
   for example
       {EXAMPLE_FILENAME}

   The filename is only a FALLBACK. It is used when a header has no [strain] or
   no protein name, and it needs a space on BOTH sides of the dash. When the
   header is complete the filename can be anything - but a consistent name keeps
   this folder readable and keeps the fallback working.

4. THINGS THAT STOP THE RUN AT STAGE 1
   * Every (protein, strain) pair must be unique. Two Zaire spike sequences
     with identical [bracket text] are rejected as duplicates. To supply
     several isolates of one species, make the bracketed text differ:
         [Zaire ebolavirus Makona]      [Zaire ebolavirus Mayinga]
   * A sequence that is nucleotide, empty, or shorter than {min_length} residues.

5. HOW A FILE IS ASSIGNED TO A PROTEIN
   The protein name is matched against sequences.protein_aliases in config.yaml.
   Currently configured:
{alias_lines}

   Matching is by substring, longest alias first, so a name that merely
   CONTAINS an alias gets that label. The Ebola GP gene also encodes a "secreted
   glycoprotein"; that name contains "glycoprotein" and would be labelled
   Spike. A name matching no alias becomes a protein of its own (CamelCased),
   so add an alias for every protein you analyse.

CHECK YOUR FILES BEFORE RUNNING
   python setup/init_input.py

THEN RUN
   snakemake --use-conda --cores 8
"""


def readme_text(
    input_dir: str,
    extensions: Sequence[str],
    aliases: dict[str, str],
    min_length: int = 50,
) -> str:
    """The notice with a header saying where it came from."""
    return (
        "(Written by setup/init_input.py from config.yaml. It is rewritten if the\n"
        " configuration changes, so put your own notes somewhere else.)\n\n"
        + notice(input_dir, extensions, aliases, min_length)
    )


def ensure_input_dir(
    input_dir: str,
    extensions: Sequence[str],
    aliases: dict[str, str],
    min_length: int = 50,
) -> bool:
    """Create *input_dir* and its ``README.txt``. Returns True if it was created.

    Idempotent, and deliberately quiet about it: the README is written only when
    it is missing or out of date. That matters because stage 1 declares this
    directory as an input, and Snakemake compares a directory's modification
    time, which moves whenever an entry is created or removed. Rewriting an
    unchanged README in place does not move it; creating one does, so running
    this repeatedly must not.

    Never touches sequence files.
    """
    created = not os.path.isdir(input_dir)
    os.makedirs(input_dir, exist_ok=True)

    path = os.path.join(input_dir, README_NAME)
    text = readme_text(input_dir, extensions, aliases, min_length)
    current = None
    if os.path.exists(path):
        with open(path, encoding="utf-8") as handle:
            current = handle.read()
    if current != text:
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
    return created


def describe(
    input_dir: str,
    extensions: Sequence[str],
    aliases: dict[str, str],
    min_length: int = 50,
) -> tuple[list[str], bool]:
    """Report how each input file will be read, and whether stage 1 will accept it.

    Returns ``(lines, ok)``. ``ok`` is False only when stage 1 would reject the
    inputs; a warning (an unrecognised strain, say) does not make it False,
    because stage 1 accepts those - it is just worth knowing about.
    """
    paths = find_inputs(input_dir, extensions)
    if not paths:
        return [f"No sequence files in {input_dir!r} yet."], True

    lines: list[str] = []
    warnings: list[str] = []
    unreadable: list[str] = []
    records: list[fasta.Record] = []

    for path in paths:
        name = os.path.basename(path)
        try:
            batch = fasta.read_fasta(path)
        except (OSError, UnicodeDecodeError) as error:
            unreadable.append(f"  {name}: cannot be read as text ({error})")
            continue
        if not batch:
            unreadable.append(f"  {name}: contains no FASTA record")
            continue

        for record in batch:
            record.protein_label = fasta.normalise_protein_label(record.name, aliases)
            records.append(record)
            lines.append(name)
            lines.append(
                f"      protein {record.protein_label} | strain {record.strain} | "
                f"accession {record.accession} | {len(record)} aa"
            )
            if record.strain == UNKNOWN_STRAIN:
                warnings.append(
                    f"  {name}: no [strain] in the header, and the filename does not "
                    "match '<ACCESSION> <protein> - <strain>'. It will be treated as "
                    f"'{UNKNOWN_STRAIN}'."
                )
            if record.name == UNKNOWN_PROTEIN:
                warnings.append(
                    f"  {name}: no protein name in the header, and the filename does "
                    "not match '<ACCESSION> <protein> - <strain>'."
                )

    by_protein: dict[str, set[str]] = defaultdict(set)
    for record in records:
        by_protein[record.protein_label].add(record.strain)
    if by_protein:
        lines.append("")
        lines.append(
            f"{len(paths)} file(s), {len(by_protein)} protein(s): "
            + ", ".join(
                f"{label} ({len(strains)} strain{'s' if len(strains) != 1 else ''})"
                for label, strains in sorted(by_protein.items())
            )
        )

    problems = unreadable + find_problems(records, min_length)
    if warnings:
        lines.append("")
        lines.append("WARNINGS (accepted, but check they are what you intended):")
        lines.extend(warnings)
    if problems:
        lines.append("")
        lines.append("PROBLEMS - stage 1 would reject these inputs:")
        lines.extend(problems)

    return lines, not problems
