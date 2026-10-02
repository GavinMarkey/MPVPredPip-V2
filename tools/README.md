# External Tools

The pipeline orchestration, parsing, and filtering logic are self-contained. The **predictors are not** — they must be installed separately. This file documents what each stage needs and how to install it.

---

## Summary

| Stage | Tool | Installation | Status |
| --- | --- | --- | --- |
| 2, 6 | IAPred | `bash setup/install_iapred.sh` | Automatic |
| 3 | BepiPred 3.0 | `bash setup/install_bepipred.sh` | Automatic |
| 3 | EpiDope | Automatic (conda) | Automatic |
| 4 | IEDB MHC-I binding | `bash setup/install_iedb.sh` | Automatic |
| 4 | IEDB MHC-II binding | `bash setup/install_iedb.sh` | Automatic |
| 4 | IEDB immunogenicity | `bash setup/install_iedb.sh` | Automatic |
| 5 | IEDB cluster | Reimplemented in `mpvpredpip.cluster` | Built-in |
| 8 | IEDB population coverage | `bash setup/install_iedb.sh` | Automatic |
| 6 | AlgPred 2.0 | `bash setup/install_algpred.sh` | Automatic |
| 6 | Autoimmunity screening | Reimplemented in `mpvpredpip.autoimmunity` | Built-in |
| 7 | ConSurf | Optional (off by default) | Manual |

All tools are freely available for non-commercial use. Run `python setup/check_tools.py` after installing to validate your setup.

---

## 1. IAPred

**Source:** <https://github.com/sebamiles/IAPred>

**Installation:**
```bash
bash setup/install_iapred.sh
```

### Configuration

The adapter was validated against the tool's README and tested in a real run:

| Parameter | Value | Notes |
| --- | --- | --- |
| Entry point | `IApred.py` | Lowercase `p` (matters on Linux) |
| Arguments | `python IApred.py input.fasta [output.csv]` | Positional, not flags |
| Score column | `IAscore` | Ranges from -3 to +3 (not 0–1) |
| Score range | -3 to +3 | Negative scores are normal |

### Known limitations

- **Deflines truncated to 20 characters:** The pipeline uses longer identifiers (e.g., `Nucleoprotein\|Bundibugyo_virus\|AYI50379.1`), so sequences are submitted under surrogate IDs (`s0`, `s1`, …) and mapped back. A `.ids.tsv` sidecar records the mapping.
- **Sequences under 20 aa are not scored:** IApred returns `Sequence too short` in place of a numeric score. These peptides fail stage 6 with that reason recorded.

### Threshold

The default threshold in `config.yaml` is `-0.3` (IApred's Moderate/Low boundary):

| Category | Score |
| --- | --- |
| High | > 0.3 |
| Moderate | -0.3 to 0.3 |
| Low | < -0.3 |

Verify this threshold against your biological hypothesis before reporting results. If you change `antigenicity.tool`, re-tune the threshold — the fallback `kolaskar` method scores differently.

### Fallback: Kolaskar method

If IAPred is unavailable, set `antigenicity.tool: "kolaskar"` to use a built-in, dependency-free implementation of Kolaskar & Tongaonkar (1990). This is useful for testing the workflow before full setup.

**Caveat:** This fallback is not equivalent to IAPred and must not be reported as such.

---

## 2. BepiPred 3.0

**Source:** <https://github.com/UberClifford/BepiPred-3.0>  
**License:** Free for academic and non-commercial use

**Installation:**
```bash
bash setup/install_bepipred.sh
```

This downloads BepiPred (~140 MB) and the ESM-2 language model (~2.5 GB, cached on first run).

### Configuration

| Parameter | Value | Notes |
| --- | --- | --- |
| Entry point | `bepipred3_CLI.py` | `-pred` flag is required |
| Output format | CSV: `Accession,Residue,BepiPred-3.0 score,...` | One row per residue |
| Score column | `BepiPred-3.0 score` (raw, not rolling mean) | Compared against `-t` threshold |

### Windows note

BepiPred caches ESM-2 encodings by sequence defline. The pipeline's IDs contain `|` (pipe), which Windows filesystems reject. The cache is stored in a named Docker volume instead (`bcell.bepipred.esm_cache_dir` in `config.yaml`).

---

## 3. EpiDope

**Installation:** Automatic via conda (no manual setup needed)

EpiDope installs from the `flomock` conda channel on first pipeline run. It runs in an isolated Python 3.6 environment (required by the tool) and its output is parsed in the main environment, so no downstream code depends on Python 3.6.

---

## 4. IEDB Standalone Tools

**Source:** <https://downloads.iedb.org/tools/>  
**License:** Free for non-commercial use

**Installation:**
```bash
bash setup/install_iedb.sh
```

This downloads and configures five tools:

| Tool | Size | CLI signature (key parameters) |
| --- | --- | --- |
| MHC-I binding | 326 MB | `predict_binding.py <method> <alleles> <lengths> <fasta>` |
| MHC-II binding | 1.59 GB | `mhc_II_binding.py <method> <allele> <fasta> <length>` |
| Immunogenicity | 4 KB | One peptide per line |
| Cluster | 5 KB | (Reimplemented; included for reference) |
| Population coverage | 0.8 MB | `-p <populations> -c <class> -f <input>` |

**Important:** The MHC-II configure script is `python configure.py`, not `./configure` (shebang points to `python3`, which doesn't exist in conda environments).

### Allele panels

The pipeline uses `workflow/config/alleles/` for MHC alleles. These are "most frequent" panels; verify they match your intended analysis (e.g., if reproducing a specific Tepitool preset).

### Skip the cluster tool

To skip downloading the cluster tool (5 KB, only used for reference):
```bash
bash setup/install_iedb.sh mhc_i mhc_ii immunogenicity population_coverage
```

---

## 5. AlgPred 2.0

**Source:** <https://github.com/raghavagps/algpred2>  
**License:** GPL-3.0

**Installation:**
```bash
bash setup/install_algpred.sh
```

### Configuration

| Parameter | Value | Why |
| --- | --- | --- |
| `-d 2` | All peptides (allergen + non-allergen) | Default `-d 1` omits non-allergen peptides, making failures silent |
| `-m 2` | Hybrid model (RF + BLAST + MERCI) | Recommended; model 1 (RF-only) available via config |

### Output format

```
Subject,ML Score,MERCI Score,BLAST Score,Hybrid Score,Prediction
```

The pipeline uses `Hybrid Score` (not `ML Score`) and trusts the `Prediction` verdict.

### Known issues

1. **`rf_model` is Git LFS-tracked:** A plain clone leaves it as a stub. The installer runs `git lfs pull` and falls back to a zip mirror.
2. **`envfile` points to author's machine:** The installer writes a correct one per run, resolving `blastp` from PATH.
3. **Fixed-name scratch files:** AlgPred writes files like `seq.aac`, `final_output` to the current directory. Concurrent jobs would collide. Each invocation runs in its own temporary directory.
4. **scikit-learn version mismatch:** The model was pickled with scikit-learn 0.22 and cannot load under 1.5.2 (required by IAPred). AlgPred therefore has its own conda environment (`workflow/envs/algpred.yaml`, Python 3.8, scikit-learn 0.22).
5. **Single-sequence scoring fails:** `np.loadtxt` returns 1-D for one row, 2-D for multiple. The installer patches this with `np.atleast_2d`.

---

## 6. Autoimmunity Screening

**No tool to install** — reimplemented locally in src/mpvpredpip/autoimmunity.py`.

The screen checks whether any 9-mer (configurable) of a consensus peptide occurs in the human proteome.

**You do need the human reference proteome:**
```bash
bash setup/fetch_human_proteome.sh
```

This downloads UniProt's human reference proteome (~75 MB) to `resources/human/UP000005640_9606.fasta`.

### Design choices

- Searches all overlapping 9-mers (any match = concern), not just the whole peptide.
- Matching is exact. L/I-equivalence (isobaric matching) is **not** implemented. Candidates differing from human windows only at L/I will not be flagged (a known limitation, not a silent relaxation).
- Windows with ambiguous residues (`X`, `N`, etc.) cannot be matched literally and are skipped. A non-zero `windows_skipped` count means that peptide's autoimmunity verdict is incomplete.

---

## 7. ConSurf

**Optional and off by default** (`conservation.consurf.enabled: false`)

The ConSurf standalone is licence-gated. Stage 7's conservancy calculation does **not** depend on it — it is reimplemented locally from the IEDB definition so the pipeline runs unattended. ConSurf outputs are available for cross-checking via the IEDB web tool.

---

## Troubleshooting

**"Tool not found" error:**
```bash
python setup/check_tools.py
```
This lists what is installed and what needs downloading.

**Tool produces unexpected output:**
Check the relevant section above for expected column names and CLI signatures. The raw output is always preserved in the stage's `raw/` folder for inspection.

**AlgPred fails with "bad interpreter":**
Ensure MHC-II was configured as `python configure.py`, not `./configure`.

---

## All tools in `config.yaml`

Paths are configurable under `tools:` in `config.yaml`, so existing installations elsewhere can be pointed at instead of copied.
