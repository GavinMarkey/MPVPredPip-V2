# External tools

The pipeline orchestration, all parsing and all filtering logic are self
contained. The **predictors are not** — most are licence-gated, and none can be
redistributed here. This file records what each stage needs, and, importantly,
**how confident I am in the invocation details**.

> ### Verification status — read this first
>
> I was unable to reach the internet while building this pipeline, so the CLI
> flags, output filenames and column headers encoded in `src/vaxpipe/adapters/`
> are **written from prior knowledge and have not been checked against the
> current releases**. They are deliberately isolated: every adapter resolves
> columns *by name*, keeps the tool's native output untouched, and fails with the
> header it actually saw rather than mis-assigning a score.
>
> Treat the adapter layer as a first draft to be confirmed against your installs.
> Everything downstream of it — schemas, clustering rule, screening, conservancy,
> reporting — is independent of these details and does not need revisiting.
>
> Run `python setup/check_tools.py` after installing to see what resolves.

---

## Summary

| Stage | Tool | Obtainable unattended? | Confidence in adapter |
| --- | --- | --- | --- |
| 2, 6 | IAPred | **Yes** — `bash setup/install_iapred.sh` | **Confirmed by a real run** |
| 3 | BepiPred | **Yes** — `bash setup/install_bepipred.sh` | Written against its source |
| 3 | EpiDope | **Yes** — conda, via the `flomock` channel | **Confirmed by a real run** |
| 4 | IEDB MHC-I binding | **Yes** — direct download, 326 MB | Never run |
| 4 | IEDB MHC-II binding | **Yes** — direct download, 1.59 GB | Never run |
| 4 | IEDB immunogenicity | **Yes** — direct download, 4 KB | Never run |
| 5 | IEDB cluster | Not used — reimplemented in `vaxpipe.cluster` | n/a |
| 8 | IEDB population coverage | **Yes** — direct download, 0.8 MB | Never run |
| 6 | AlgPred 2.0 | **Yes** — `bash setup/install_algpred.sh` | Low, verify columns |
| 6 | PIR Peptide Match | Not used — reimplemented in `vaxpipe.autoimmunity` | n/a |
| 7 | ConSurf | Manual, optional | Off by default |

> **On "manual download".** Earlier revisions of this file listed BepiPred, all
> five IEDB tools and AlgPred 2.0 as manual or licence-gated. Every one of those
> claims turned out to be false when actually checked — each is freely
> downloadable without registration. Treat any remaining "manual" label here as
> unverified rather than established; ConSurf is the only one still carrying it.

---

## 1. IAPred — `tools/iapred/`

Source: <https://github.com/sebamiles/IAPred>. Unlike everything else here it is
on GitHub and installs unattended:

```bash
bash setup/install_iapred.sh
```

That clones the repository, installs any `requirements.txt` it declares, and
pip-installs it if it is packaged.

### The adapter has been checked against the README, but not yet run

The IAPred adapter was originally written **without access to the repository**.
Once installed, its assumptions in `_predict_iapred()` were checked against
`tools/iapred/README.md` and corrected:

| | Originally assumed | Actually |
| --- | --- | --- |
| Entry point | `iapred.py` / `IAPred.py` | `IApred.py` — lowercase `p`, which matters on Linux |
| Arguments | `-i {input} -o {output}` first | positional: `python IApred.py input.fasta [output.csv]` |
| Score column | `score` / `antigenicity` / … | `IAscore` |
| Score range | — | roughly **-3 to 3**, not 0–1; negatives are normal |

The score-column correction was the load-bearing one. `Antigenicity_Category`
holds a text label (`High`/`Moderate`/`Low`), and the old candidate list matched
it by prefix ahead of the numeric column, which would have made every score
unparseable. `iascore` now resolves first.

Two further behaviours were found by running it, and are handled in the adapter:

* **Deflines are truncated to 20 characters** (`IApred.py` line 105). This
  pipeline's identifiers are longer than that — `Nucleoprotein|Bundibugyo_virus|AYI50379.1`
  — so sequences are submitted under short surrogate names (`s0`, `s1`, …) and
  the scores mapped back. A `.ids.tsv` sidecar beside the preserved raw output
  records the mapping.
* **Sequences under 20 aa are not scored.** IApred writes `Sequence too short`
  into the score column instead of a number. Those peptides fail the stage 6
  screen with that as their recorded reason, rather than aborting the run or
  being reported as scoring below the threshold.

**Still to do: run it.** These corrections come from reading the README, not from
a successful execution. The adapter remains lenient — it discovers any `*.py`
entry point, tries four argument styles in turn, accepts output on stdout as well
as to a file, and matches columns by name — and on failure the error lists every
style it tried and what each one said, so the real signature can be read straight
off the failure.

Two ways to fix it, neither touching anything downstream:

1. Correct `_IAPRED_ENTRYPOINTS`, `_IAPRED_ARG_STYLES` and `_IAPRED_*_COLUMNS` at
   the top of `src/vaxpipe/adapters/screening.py`.
2. Bypass the adapter — set `antigenicity.tool: "external"` and give the real
   command in `antigenicity.command`, using `{input}` and `{output}`:

   ```yaml
   antigenicity:
     tool: "external"
     command: "python tools/iapred/IApred.py {input} {output}"
   ```

### The threshold

`config.yaml` now uses `0.3` for both the stage 2 screen and the stage 6 hard
filter. That is IApred's own "High antigenicity" boundary as documented in its
README:

| Category | Score |
| --- | --- |
| High | > 0.3 |
| Moderate | -0.3 to 0.3 |
| Low | < -0.3 |

It replaces a `0.4` placeholder. The two are close, but that was coincidence: the
placeholder was chosen as if the scale were a 0–1 probability, whereas IApred
scores run roughly **-3 to 3**. Anything that assumed non-negative scores was
wrong about this tool.

This value decides which candidates survive stage 6, so confirm it against the
paper before publishing. If you change `antigenicity.tool`, re-set it — the
`kolaskar` fallback scores around 1.0 and would pass everything at 0.3.

### Fallback backend

`antigenicity.tool: "kolaskar"` selects a built-in, dependency-free implementation
of the Kolaskar & Tongaonkar (1990) antigenic-propensity method, so the workflow
can be exercised end to end before IAPred is installed or if its CLI is still
being sorted out.

> It is **not** equivalent to IAPred and must not be reported as such. The
> propensity values in `screening.py` should also be checked against the original
> paper before being relied on.

## 2. BepiPred — `tools/bepipred3/`

Source: <https://github.com/UberClifford/BepiPred-3.0>. Installs unattended:

```bash
bash setup/install_bepipred.sh
```

**This file previously said BepiPred needed an academic licence and could not be
fetched unattended. That was wrong** and it blocked stage 3 for longer than it
should have. DTU Health Tech publish the code *and the trained weights* on
GitHub. They are free for academic and other non-commercial use with no form to
sign; only for-profit use needs a separate licence (contact morni@dtu.dk). The
repository is ~140 MB because the ensemble weights are committed to it.

Two downloads, not one. The ESM-2 language model (~2.5 GB) is fetched separately
by `fair-esm` the first time a prediction runs, into the torch hub cache.

The adapter was written against the tool's own source rather than its README,
after the IAPred experience below. Verified there:

- the entry point is `bepipred3_CLI.py`, and `-pred` is **required**;
- `raw_output.csv` is `Accession,Residue,BepiPred-3.0 score,BepiPred-3.0 linear
  epitope score` — one row per residue in sequence order, with **no position
  column**, and `Residue` holding the amino acid letter rather than an index;
- the pipeline reads the raw score, not the rolling-mean "linear epitope score"
  beside it, because the raw score is what BepiPred's own `-t` is compared
  against. So `bcell.bepipred.threshold` in `config.yaml` means the same thing
  `-t` means, and our epitopes agree with the `Bcell_epitope_preds.fasta` kept
  alongside them.

One caveat with a Windows host: BepiPred names each cached ESM-2 encoding after
the sequence defline, and ours contain `|`, which Windows will not accept in a
filename. The cache therefore has to live outside the bind-mounted project
directory — `bcell.bepipred.esm_cache_dir` in `config.yaml` points at a named
volume in the dev container. Every actual output file has a fixed name, so this
affects nothing else.

## 3. EpiDope

Nothing to install by hand. Snakemake builds `workflow/envs/epidope.yml` — the
upstream project's own `conda env export`, taken from
<https://github.com/rnajena/EpiDope> — on the first run.

It **is** conda-installable, but from the **`flomock`** channel, not bioconda.
This file and `workflow/envs/epidope.yaml` both claimed otherwise for a while;
`conda search -c bioconda -c conda-forge epidope` finding nothing is what that
claim rested on, and the conclusion drawn from it was wrong. `epidope.yaml` is
now a comment-only stub and can be deleted; the `.yml` is the real one.

That environment pins **Python 3.6**, which cannot parse this codebase. So the
EpiDope rule is split in two: the tool runs there as a bare shell command, and
its output is parsed by a second rule under the orchestrator's environment.
Nothing of ours is ever imported into the 3.6 environment.

The adapter expects a per-residue score CSV in the output directory; its
candidate filenames and column names in `src/vaxpipe/adapters/bcell.py` are
**unverified**, like IAPred's were.

## 4. IEDB standalone tools

`mhc_i`, `mhc_ii`, `immunogenicity`, `cluster`, `population_coverage`.

**This file previously said these were "gated behind a registration or licence
click, so they cannot be scripted". That is wrong** — the same mistake it made
about BepiPred. All five are plain HTTP downloads from `downloads.iedb.org`,
each with an `MD5SUM` beside it. Verified by request:

| Tool | Path under `https://downloads.iedb.org/tools/` | Size |
| --- | --- | --- |
| MHC-I | `mhci/LATEST/IEDB_MHC_I-3.1.7.tar.gz` | 326 MB |
| MHC-II | `mhcii/LATEST/IEDB_MHC_II-3.1.12.tar.gz` | 1.59 GB |
| Immunogenicity | `immunogenicity/LATEST/IEDB_Immunogenicity-3.0.tar.gz` | 4 KB |
| Cluster | `cluster/LATEST/IEDB_Cluster-1.0.tar.gz` | 5 KB |
| Population coverage | `population/LATEST/IEDB_Population_Coverage-3.0.2.tar.gz` | 0.8 MB |

Free for non-commercial use. The binding predictors are large because they
bundle every method; the other three are thin scripts. Install them with:

```bash
bash setup/install_iedb.sh
```

It verifies each download against the published `MD5SUM` and runs each tool's
configure step, which is not the same command in any two of them: `./configure`
for class I, `configure.py` for class II, `./configure` for population coverage
(which ships both), and none for immunogenicity.

A `.py` configure script has to be run as `python configure.py`, **not**
executed directly: class II's carries the shebang `#!/usr/bin/python3`, which
does not exist in a conda environment. Executed directly it fails with `bad
interpreter` and leaves the tool unconfigured, which then shows up much later as
a missing-model error in the middle of a stage 4 run.

**The cluster tool is not used.** Stage 5 clusters in-process — see
`src/vaxpipe/cluster.py` for why, and the docstring of
`src/vaxpipe/adapters/iedb.py` for what the standalone actually emits. It is
still downloaded by default, because it is 5 KB and it is the reference this
project's reimplementation is checked against. Skip it with
`bash setup/install_iedb.sh mhc_i mhc_ii immunogenicity population_coverage`.

Verified against the installed source rather than the READMEs, after IAPred:

- class I is `predict_binding.py <method> <alleles> <lengths> <fasta>` and class
  II is `mhc_II_binding.py <method> <allele> <fasta> <length>` — the input and
  length are the other way round between them;
- `netmhcpan_el` and `netmhciipan_el`, the methods in `config.yaml`, both exist
  in these versions (NetMHCpan 4.1 / NetMHCIIpan 4.1);
- class I's output columns are
  `allele, seq_num, start, end, length, peptide, core, icore, Score, rank` —
  note the percentile column is `rank`, not `percentile_rank`, for this method;
- immunogenicity takes one peptide per line, which is what the adapter writes;
- population coverage is `-p <populations> -c <class> -f <input>`.

Unpack each into the path named in `config.yaml` under `tools:`. Verify the CLI
signatures in `src/vaxpipe/adapters/iedb.py` against your versions — in
particular that the class I predictor takes
`<method> <alleles> <lengths> <fasta>` and that allele and length are parallel
comma-separated lists.

### Tepitool

Tepitool has no standalone release — it is a guided web front end over these same
predictors. The pipeline reproduces it locally with an explicit allele panel
(`workflow/config/alleles/`) and an explicit percentile cut-off
(`config.yaml`). **Check both against the Tepitool preset you intend to
reproduce**; the panels here are a reasonable "most frequent alleles" set, but
the exact membership of a given preset should be confirmed.

## 5. AlgPred 2.0 — `tools/algpred2/`

Raghava group standalone, GPL-3.0, freely available:
<https://github.com/raghavagps/algpred2>. Install with
`bash setup/install_algpred.sh`.

Verified CLI:

```
algpred2.py -i INPUT [-o OUTPUT] [-t THRESHOLD] [-m {1,2}] [-d {1,2}]
```

The adapter overrides two defaults, both deliberately:

* **`-d 2`** — AlgPred defaults to `-d 1`, which writes out *only* the peptides
  it calls allergenic. Under that default the peptides that pass the screen are
  absent from the file, and "nothing was allergenic" is indistinguishable from
  "the tool never ran" — the same silent-failure shape as the MHC-I bug in
  stage 4. `predict_allergenicity()` additionally cross-checks the output
  against the submitted FASTA and raises if any peptide is missing.
* **`-m 2`** — the hybrid model (RF + BLAST + MERCI), which is what
  `screening.allergenicity.threshold` is characterised against. AlgPred defaults
  to model 1 (composition RF only).

Model 2 is why the install clones the repository rather than using
`pip install algpred2`: the hybrid needs the bundled `Database/` and `progs/`
(MERCI, in Perl), which the wheel does not provide. Set
`screening.allergenicity.model: 1` to avoid all of it.

Output columns, read off `algpred2.py`'s `hybrid()` rather than guessed:

```
Subject,ML Score,MERCI Score,BLAST Score,Hybrid Score,Prediction
```

The identifier column is `Subject`, and the score taken is `Hybrid Score`, not
`ML Score`. The `Prediction` verdict (`Allergen` / `Non-Allergen`) is trusted
rather than recomputed, so the boundary stays the tool's to define.

### Four things the clone does not give you

These are handled in `predict_allergenicity()`; they are recorded here because
each one fails in a way that does not name its own cause.

1. **`rf_model` is Git LFS-tracked.** A plain `git clone` leaves it absent or
   as a ~130-byte pointer stub, and AlgPred then dies inside joblib with an
   error that never mentions LFS. `setup/install_algpred.sh` runs `git lfs pull`
   and falls back to the zip mirror; `check_tools.py` verifies the real file is
   there, since the entry point existing says nothing about whether it can run.
2. **`envfile` points at the author's machine** (`/home/neelam/...`). AlgPred
   reads it from the current directory, so the adapter writes a correct one per
   run, resolving `blastp` from PATH — which is necessary because BLAST comes
   from the conda environment, whose path is not known at install time.
3. **It writes about ten fixed-name files into the current directory**
   (`seq.aac`, `Sequence_1`, `RES_1_6_6.out`, `final_output`, …). Stage 6 fans
   out over proteins, so two concurrent jobs sharing a directory would overwrite
   each other's intermediates — producing wrong scores rather than a crash.
   Each call therefore runs in its own temporary directory of symlinks.
4. **`rf_model` cannot be unpickled by a modern scikit-learn.** Loading it
   under 1.5.2 fails in `sklearn.tree._tree`:

   ```
   ValueError: node array from the pickle has an incompatible dtype
   - expected: [... 'missing_go_to_left'] itemsize 64
   - got:      [... 'weighted_n_node_samples']
   ```

   `missing_go_to_left` entered the tree node dtype in scikit-learn 1.3, so the
   model predates it — consistent with `algpred2.py` importing
   `sklearn.externals.joblib`, which was removed in 0.23. IApred needs 1.5.2
   for *its* models, so one environment cannot serve both.

   AlgPred therefore has its own rule and its own `envs/algpred.yaml`
   (Python 3.8, scikit-learn 0.22), with the invocation in
   `src/vaxpipe/algpred.py` — standard-library-only so it can be imported
   there. Only the CSV parsing stays in `screening.yaml`. This is the same
   split stage 3 uses for EpiDope, and the collision was predicted by a comment
   in `screening.yaml` before it happened.

   Switching to `model: 1` is **not** an escape: model 1 loads the same
   `rf_model`.
5. **It cannot score a single sequence.** `prediction()` does
   `np.loadtxt(file_name, delimiter=',')`, which returns a 1-D array when the
   file holds one row and a 2-D array otherwise, so `predict_proba` gets
   `(20,)` instead of `(1, 20)`:

   ```
   ValueError: Expected 2D array, got 1D array instead
   ```

   This is an upstream bug affecting both models, and it is not an edge case
   here — stage 5 routinely yields one consensus peptide for a protein, which
   is exactly what happened to the nucleoprotein while the spike glycoprotein
   (8 peptides) ran fine. The `single-sequence` entry in `PATCHES` wraps the
   call in `np.atleast_2d`, a no-op for every other input.

## 6. Autoimmunity — no tool

PIR Peptide Match is **not used**. The screen is reimplemented in
`src/vaxpipe/autoimmunity.py`, exactly as the IEDB cluster tool is reimplemented
in `vaxpipe.cluster`.

The tool is free (MIT, <https://github.com/udel-cbcb/PeptideMatch>) but awkward
to obtain in the form this pipeline needs: it publishes no releases and ships no
built jar, `peptidematch_cmd` is an Eclipse project with no `pom.xml`, and the
maintained path has moved to an Elasticsearch service wanting Docker, Java 17,
Maven and a multi-gigabyte heap. What the screen computes — does any k-mer of a
~25 aa peptide occur verbatim in a 20,000-protein FASTA — takes about a second
in Python. The index rule and the JRE went with it.

You still need the **human reference proteome** at the path in
`screening.autoimmunity.human_proteome` (default
`resources/human/UP000005640_9606.fasta`): `bash setup/fetch_human_proteome.sh`.

Two properties of the reimplementation worth knowing:

* It queries every overlapping 9-mer (configurable) of each consensus peptide,
  not the whole peptide, since a candidate is a concern if *any* part of it
  occurs in the human proteome.
* Matching is **exact**. PIR's optional L/I-equivalence mode (isobaric leucine
  and isoleucine treated as identical) is off there by default and is not
  implemented here. Enabling it would flag *more* peptides, so its absence is
  not a silent relaxation — but a candidate differing from a human window only
  at an L/I position will not be flagged.
* A window containing an ambiguous residue cannot be searched for literally, so
  it is skipped and counted in `windows_skipped`. A non-zero count means that
  peptide's autoimmunity verdict is **incomplete, not clean**.

## 7. ConSurf — `tools/consurf/`

Optional and **off by default** (`conservation.consurf.enabled: false`). The
standalone is licence-gated and pulls a large dependency stack; its helpers
(MAFFT, CD-HIT, Rate4Site, RAxML) are already in
`workflow/envs/conservation.yaml`, so enabling it later needs no image rebuild.

Stage 7's conservancy calculation does **not** depend on ConSurf — it is
reimplemented locally from the IEDB definition so the pipeline runs unattended,
and it also writes upload-ready files for the IEDB web tool as a cross-check.

---

## Layout

```
tools/
├── bepipred3/          BepiPred (manual)
├── mhc_i/              IEDB MHC class I binding
├── mhc_ii/             IEDB MHC class II binding
├── immunogenicity/     IEDB class I immunogenicity
├── cluster/            IEDB epitope cluster analysis
├── population_coverage/ IEDB population coverage
├── algpred2/           AlgPred 2.0 (setup/install_algpred.sh)
├── iapred/             your antigenicity predictor (see section 1)
└── consurf/            ConSurf standalone (optional)
```

Paths are all configurable — `tools:` in `config.yaml` — so an existing
installation elsewhere can be pointed at instead of copied.
