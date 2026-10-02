# MPVPredPip v2.0 — Multivalent Peptide-based Vaccine Prediction Pipeline

An automated, reproducible pipeline for predicting multivalent peptide-based vaccine candidates. Everything runs offline in a containerized Linux environment. Each stage writes its output into a numbered folder, and the output of each stage is the documented input of the next — so any stage can be inspected, re-run, or replaced independently.

Originally developed for **Ebola** virus strains, now generalized for any pathogen protein sequences.

**Code author:** [Gavin Markey](https://github.com/GavinMarkey)  
**Maintainer:** [Shukla Lab](https://shuklalab.github.io/)

---

## System Requirements

| Requirement | Windows | Mac | Linux |
| --- | --- | --- | --- |
| Docker | Docker Desktop | Docker Desktop | Docker Engine |
| Disk space | 30 GB (tools + environments + outputs) | 30 GB | 30 GB |
| RAM | 8 GB minimum (16 GB recommended) | 8 GB minimum | 8 GB minimum |
| Internet | ✅ (first run downloads ~4 GB) | ✅ | ✅ |
| **Optional: VSCode** | Dev Containers extension | Dev Containers extension | Dev Containers extension |

**First-run timing:** 30–60 minutes (includes downloading pre-trained models for BepiPred and IEDB tools)

---

## Quick Start

### Option 1: VSCode + Dev Containers (Easiest — Recommended)

1. **Install prerequisites:**
   - [Docker Desktop](https://www.docker.com/products/docker-desktop) (Windows/Mac) or [Docker Engine](https://docs.docker.com/engine/install/) (Linux)
   - [VSCode](https://code.visualstudio.com/)
   - VSCode extension: [Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

2. **Clone and open:**
   ```bash
   git clone https://github.com/YourName/MPVPredPip-V2.git
   cd MPVPredPip-V2
   code .
   ```

3. **Start the container:**
   - VSCode detects `.devcontainer/devcontainer.json` and prompts: **"Reopen in Container"**
   - Click the button (or use Command Palette: `Remote-Containers: Reopen in Container`)
   - First build takes ~5–10 minutes

4. **Setup happens automatically:**
   - The container installs the pipeline, validates tools, and creates `0 - Input/` folder
   - Watch the terminal for `✅` next to each tool
   - Any missing tools are reported with download instructions

5. **Add your sequences:**
   - Open `0 - Input/` folder in VSCode
   - Add your NCBI Protein FASTA files
   - Run `python setup/init_input.py` to validate filenames and preview how they'll be read

6. **Run the pipeline:**
   ```bash
   snakemake --use-conda --cores 8
   ```
   - First run downloads deep-learning models (~2 GB), subsequent runs use cache
   - Progress shown in terminal; intermediate results in numbered folders

7. **Check results:**
   - Final candidates: `9 - Final_Consensus_Epitopes/final_candidates.tsv`
   - Pipeline summary: `9 - Final_Consensus_Epitopes/PIPELINE_REPORT.md`
   - Attrition report (how many epitopes passed each stage): `9 - Final_Consensus_Epitopes/attrition.tsv`

---

### Option 2: Command Line (CLI) – No VSCode Needed

For users comfortable with Docker Compose commands:

1. **Install Docker** (see Option 1 above)

2. **Clone:**
   ```bash
   git clone https://github.com/YourName/MPVPredPip-V2.git
   cd MPVPredPip-V2
   ```

3. **Start the container:**
   ```bash
   docker-compose -f .devcontainer/docker-compose.yml up -d
   ```

4. **Run setup:**
   ```bash
   docker-compose -f .devcontainer/docker-compose.yml exec mpvpredpip \
     bash -c "python setup/init_input.py && python setup/check_tools.py"
   ```

5. **Add your sequences to `0 - Input/`**

6. **Run the pipeline:**
   ```bash
   docker-compose -f .devcontainer/docker-compose.yml exec mpvpredpip \
     snakemake --use-conda --cores 8
   ```

7. **Stop when done:**
   ```bash
   docker-compose -f .devcontainer/docker-compose.yml down
   ```

---

## Preparing Your Input Sequences

**The header line is what the pipeline reads, not the filename:**

```
>ACCESSION protein_name [strain]
>NP_066246.1 spike glycoprotein [Zaire ebolavirus]
```

- **Accession** = first word (NCBI accession number)
- **Protein name** = text before the `[...]` (matched against `config.yaml` aliases)
- **Strain** = text inside `[...]` (usually the organism/species name)

**Filename fallback** (optional, only used if header is incomplete):
```
NP_066246.1 spike glycoprotein - Zaire ebolavirus.fasta
```
Must have a space on **both sides** of the dash.

**Rules:**
- Each `(protein, strain)` pair must be unique
- Sequences must be protein (not nucleotide), non-empty, and ≥50 residues
- Protein aliases in `config.yaml` are matched by substring (longest first)

**Verify your files before running:**
```bash
python setup/init_input.py
```
This prints how each file will be read and validates the setup. Same rules are saved to `0 - Input/README.txt`.

---

## Understanding the Pipeline

### Stages

| # | Name | Input | Output | Time |
| --- | --- | --- | --- | --- |
| 0 | Input folder | Your FASTA files | – | – |
| 1 | Sequence retrieval | `0 - Input/*.fasta` | Validated manifest + normalised FASTA | <1 min |
| 2 | Whole-protein antigenicity | Sequences | Antigenicity scores (informational) | <1 min |
| 3 | B-cell epitopes | Sequences | Predicted epitopes from BepiPred + EpiDope | 10–20 min |
| 4 | T-cell epitopes | Sequences | MHC-I / MHC-II binding predictions | 15–30 min |
| 5 | Consensus clustering | B & T epitopes | Peptide clusters with both B & T epitopes | <5 min |
| 6 | Antigenicity / allergenicity / autoimmunity | Consensus peptides | Filtered candidates (screen for safety/efficacy) | 5–10 min |
| 7 | Conservation | Filtered peptides | Cross-strain conservancy analysis | 5–10 min |
| 8 | Population coverage | Conserved peptides | HLA coverage across global populations | 10–20 min |
| 9 | Final candidates | All screens | `final_candidates.tsv`, `PIPELINE_REPORT.md` | <1 min |

### Output Files

Each stage folder contains:
- **`summary.md`** — Human-readable overview of results
- **`<output>.tsv`** — Structured results (can open in Excel)
- **`raw/`** — Untouched predictor output (for debugging/verification)

### Selection Rule

A consensus peptide passes all stages only if it contains:
- **≥1 MHC-I epitope** AND
- **≥1 MHC-II epitope** AND
- **≥1 B-cell epitope** (from BOTH BepiPred AND EpiDope)

Override in `config.yaml`:
```yaml
clustering:
  require_both_bcell_tools: false  # Accept B-cell from either tool
```

---

## Configuration

All parameters live in `config.yaml` — **nothing is hard-coded**. Edit these before running:

```yaml
project:
  name: "MPVPredPip"              # Name for your analysis run

sequences:
  input_dir: "0 - Input"          # Where your FASTA files go
  min_length: 50                  # Minimum sequence length (residues)
  protein_aliases:                # Map NCBI names to short labels
    "spike glycoprotein": "Spike"
    "nucleoprotein": "Nucleoprotein"

screening:
  antigenicity:
    threshold: -0.3               # IApred Moderate/Low boundary (peptides)
  allergenicity:
    model: 2                       # AlgPred hybrid model (ML + BLAST + MERCI)
    threshold: 0.3                # Allergenicity threshold

tcell:
  mhc_i:
    percentile_cutoff: 1.0        # Top 1% of binders
  mhc_ii:
    percentile_cutoff: 10.0       # Top 10% of binders
```

See `config.yaml` for the complete list and detailed comments.

---

## Troubleshooting

### "No candidates at the end"
Read `9 - Final_Consensus_Epitopes/attrition.tsv` to see how many epitopes survived each stage per protein. A stage that drops everything is usually a threshold, not a biological signal. Then check that stage's `summary.md` for details.

### "Permission denied" on a file (Windows only)
The container owns the file, but Windows locks it if you have it open in Excel or another program. **Close the file** and re-run.

### "Tool not found" or "Download failed"
Run `python setup/check_tools.py` to see what's installed and what needs downloading. Some tools require manual download (e.g., AlgPred's model); check `tools/README.md` for instructions.

### "ModuleNotFoundError: No module named 'X'"
The container's conda environments may be incomplete. Rebuild:
```bash
docker-compose -f .devcontainer/docker-compose.yml down
docker-compose -f .devcontainer/docker-compose.yml up -d
```

### Re-running one stage
Delete that stage's output folder and re-run `snakemake --use-conda --cores 8`:
```bash
rm -rf "6 - Antigenicity_Allergenicity_AA_Checks"
snakemake --use-conda --cores 8
```
Snakemake redoes only what depends on the deleted output. Use `--rerun-incomplete` to recover from an interrupted run.

---

## What's New in v2.0

| Feature | Previous (SARS-CoV-2) | Now (MPVPredPip v2.0) |
| --- | --- | --- |
| **B-cell epitopes** | Manual Web form parsing | BepiPred + EpiDope locally |
| **T-cell epitopes** | Tepitool web form | IEDB standalone with versioned alleles |
| **Clustering** | IEDB web tool + hand-edited CSVs | Automated IEDB cluster + composition rule |
| **Provenance** | Epitope IDs in deflines | Typed TSV tables (tool, strain, allele, coords) |
| **Repetition** | Scripts copied per-protein | One workflow; discovers proteins from inputs |
| **Reproducibility** | Manual steps, no threshold record | Every parameter in `config.yaml`, all raw outputs saved |
| **Extensibility** | Hard-coded tools | Modular adapters; swap predictors without touching downstream |

Scientific logic is unchanged — same stages, same selection rule. The difference is full automation and transparency.

---

## File Structure

```
config.yaml                      Every tunable parameter
workflow/
  Snakefile                      Stage definitions
  rules/*.smk                    One rule file per stage
  scripts/*.py                   Thin wrappers calling mpvpredpip
  envs/*.yaml                    Per-rule conda environments
  config/alleles/                MHC allele panels
src/mpvpredpip/                     Core library
  adapters/                      Tool-specific wrappers
  schemas.py                     Data structures
  clustering.py, conservation.py etc.
tests/
  test_*.py                      Pytest suite
setup/
  init_input.py                  Creates 0 - Input/, validates filenames
  check_tools.py                 Preflight check: what tools are installed
tools/
  bepipred3/                     B-cell epitope predictor
  iapred/                        Antigenicity predictor
  iedb/                          T-cell + clustering tools
  algpred/                       Allergenicity predictor
resources/
  human/                         Human proteome (for autoimmunity screen)
  mhc_alleles/                   HLA allele panels

0 - Input/                       **YOUR sequences** (you create/populate)
1 - Sequence_retrieval/          Generated: manifest + normalised FASTA
2 - Antigenicity_prediction/     Generated: whole-protein scores
3 - B-cell_prediction/           Generated: BepiPred + EpiDope results
4 - T-cell_epitopes/             Generated: MHC-I/II predictions
5 - B_T-cell_consensus/          Generated: clusters + consensus peptides
6 - Antigenicity_Allergenicity_AA_Checks/   Generated: safety screens
7 - Conservation_Analysis/       Generated: cross-strain conservation
8 - Population_coverage/         Generated: HLA coverage estimates
9 - Final_Consensus_Epitopes/    Generated: final candidates + report
```

**Design principle:** External tools are wrapped at the edge, raw output is preserved, and every stage boundary is a plain TSV. Swap a predictor by editing one adapter; nothing downstream changes.

---

## Known Limitations & Caveats

- **First-run downloads:** BepiPred downloads ESM-2 models (~2 GB) on the first B-cell prediction run. Subsequent runs use the cache.
- **Tool adapters unverified:** Predictor invocations (src/mpvpredpip/adapters/`) are based on tool source code and documentation, not exhaustive testing. See `tools/README.md` for which tools have been confirmed on real data.
- **Antigenicity threshold:** The default (`-0.3`) is IApred's Moderate/Low boundary, suitable for peptides. Verify it matches your biological hypothesis before reporting results.
- **AlgPred model:** The hybrid model uses BLAST and sequence motifs, which have limited power at 20–30 aa length. The pipeline documents this in stage 6's output.
- **No wet-lab validation:** This is a **computational screen**, not a guarantee of immunogenicity or safety. Validation requires lab testing.
- **Stages 10+:** Docking, molecular dynamics, and manuscript generation are outside this workflow. Stage 9's FASTA is the handoff point.

---

## For Developers

Run the test suite before submitting changes:
```bash
python -m pytest tests/ -v
```

Code style: [Ruff](https://github.com/astral-sh/ruff) (auto-formatted on commit).

---

## Citation

If you use this pipeline, cite:
- **Markey et al.** (pending publication)
- **BepiPred 3.0** — Stricher et al., *Nucleic Acids Research* (2024)
- **EpiDope** — Kringelum et al., *Immunity* (2012)
- **IEDB** — Vita et al., *Nucleic Acids Research* (2018)

---

## License

Creative Commons Zero (CC0) — Public Domain. Use freely; no attribution required (but we'd appreciate it).

---

## Questions or Issues?

- **Setup/installation:** Check [Troubleshooting](#troubleshooting) above
- **Pipeline logic:** See `config.yaml` comments and stage `summary.md` files
- **Tool-specific questions:** Check `tools/README.md`
- **File naming/input validation:** Run `python setup/init_input.py`

---

**Last updated:** October 2026  
**Current version:** 2.0  
**Status:** Production-ready for vaccine candidate screening
