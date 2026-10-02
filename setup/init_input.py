#!/usr/bin/env python
"""Create the input folder, explain how to fill it, and check what is in it.

    python setup/init_input.py

Run automatically when the dev container is created, and safe to run again at any
time: it creates only what is missing and never touches your sequence files.

It exits 0 when the folder is empty or its files are acceptable, and 1 only when
files are present that stage 1 would reject - so it doubles as a quick preflight
before a long run.

The rules themselves live in ``src/vaxpipe/inputs.py`` and are the same text the
pipeline prints in its own errors, so there is one place to correct them.
"""

from __future__ import annotations

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

try:
    import yaml
except ImportError:  # pragma: no cover - the dev container always has it
    sys.exit("PyYAML is required: conda install -c conda-forge pyyaml")

from vaxpipe import inputs  # noqa: E402

GREEN, RED, YELLOW, DIM, RESET = "\033[32m", "\033[31m", "\033[33m", "\033[2m", "\033[0m"


def main() -> int:
    os.chdir(ROOT)
    with open("config.yaml", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    sequences = config["sequences"]
    input_dir = sequences["input_dir"]
    extensions = list(sequences["extensions"])
    aliases = dict(sequences.get("protein_aliases", {}))
    min_length = int(sequences.get("min_length", 50))

    created = inputs.ensure_input_dir(input_dir, extensions, aliases, min_length)

    bar = "=" * 72
    print(bar)
    if created:
        print(f"{GREEN}Created {input_dir!r}{RESET} - this is where your sequences go.")
    else:
        print(f"{input_dir!r} already exists.")
    print(f"{DIM}The same rules are saved in {input_dir}/{inputs.README_NAME}{RESET}")
    print(bar)
    print()
    print(inputs.notice(input_dir, extensions, aliases, min_length))

    print(bar)
    print("Files currently in the folder")
    print(bar)
    lines, ok = inputs.describe(input_dir, extensions, aliases, min_length)
    for line in lines:
        print(line)

    if not ok:
        print(f"\n{RED}Fix the problems above before running the pipeline.{RESET}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
