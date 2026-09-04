# HomoLogic 2
[![Python](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/1174208088.svg)](https://doi.org/10.5281/zenodo.20756173)

**HomoLogic 2** (`homologic2.py`) is a **single-file command-line driver** for the HomoLogic workflow.

It runs the same **integrative enzyme homology analysis** — sequence similarity, structural comparison, binding-site analysis and cavity characterization — but **installs its own environment, pins its own tools, and can be driven entirely from the shell** rather than from a hand-written Python driver script.

<img width="4757" height="3217" alt="image" src="https://github.com/user-attachments/assets/3d708e70-0325-432d-9c5a-e07b45a3396a" />

The eleven analysis stages are **behaviourally identical** to `homologic.py` — stages 1, 4 and 6 reproduce an existing run's outputs byte-for-byte. What is new sits around them:

- `--install` / `--check` — builds and verifies the whole environment
- **pinned executables** — the `boltz` that runs is the one in *this* environment, and it is logged
- **every Boltz entity type** — peptide, DNA, RNA, CCD and SMILES ligands, not SMILES alone
- **MSA reuse** — hand Boltz a precomputed MSA instead of calling the MSA server
- **ensemble generation** — many seeds × many diffusion samples per homolog
- **safe re-runs** — nothing is deleted without `--overwrite`; `--resume` skips finished work

---

# Overview

| file | what it is |
|----|----|
| `homologic.py` | the original `HomoLogic` class. No CLI — needs a hand-written driver script per run. |
| `homologic2.py` | the same pipeline as a self-contained CLI, plus the features above. Also importable: `from homologic2 import HomoLogic`. |

`homologic2.py` is a superset. Existing driver scripts that `from homologic import HomoLogic` keep working unchanged.

---

# Installation

```
python homologic2.py --install        # build the `homologic` conda environment
conda activate homologic
python homologic2.py --check          # verify it
```

| Flag | Description |
|----|----|
| `--install` | Creates the environment from a curated conda + pip spec. Resolves on any linux-64 machine. |
| `--install --pinned` | Reproduces the reference environment exactly from `homologic.yml`. |
| `--install --force` | Removes and rebuilds an existing environment. Without it, an existing environment is never touched. |
| `--install --dry-run-install` | Prints the exact commands and runs nothing. |
| `--install --env-name NAME` | Builds under a different environment name. |
| `--check` | Imports every required module, locates every binary, reports CUDA. Non-zero exit on any failure. |
| `--install-help` | Prints the equivalent manual commands. |

`--install`, `--check` and `--help` run under a **bare interpreter**, before the environment exists — the scientific imports are deferred until a pipeline stage actually runs.

**torch is installed before Boltz**, from the CUDA 12.8 wheel index. Letting Boltz resolve torch itself can pull a CPU-only or pre-cu128 wheel, which cannot drive an sm_120 (RTX 50-series) card.

Boltz downloads ~7 GB of weights into `~/.boltz` on first use.

---

# Inputs

Unchanged from HomoLogic: a **query** (reference enzyme) and a **haystack** (homolog candidates).

## Query (reference enzyme)

```
reference.fasta
reference.pdb
reference_ligands.pdb
```

| File | Description |
|----|----|
| `reference.fasta` | Reference enzyme amino acid sequence |
| `reference.pdb` | Reference protein structure |
| `reference_ligands.pdb` | Coordinates of bound ligand(s) |

## Haystack (homolog candidate sequences)

```
homologs.fasta
```

Putatively homologous sequences from similarity searches (BLAST, MMseqs2), structure search (Foldseek), pLM-BLAST, enzyme mining tools, orthology detection or curated libraries.

---

# Command-line interface

`-h` and `--help` deliberately differ:

| Command | Output |
|----|----|
| `homologic2.py -h` | Short — the usage line only. |
| `homologic2.py --help` | Full manual — the module docstring, then every option with its default. |
| `homologic2.py <stage> --help` | Every option for one stage. |

---

# Pipeline Workflow

Eleven stages. Each is a subcommand and each is one method on `HomoLogic`, so the pipeline can be driven from the shell or from Python.

| # | Subcommand | Method | Description |
|----|----|----|----|
| 1 | `seq-homology` | `calculate_sequence_homology()` | MMseqs2 search of reference vs. homologs, with a scikit-bio global-alignment fallback for sequences MMseqs2 reports no usable identity for |
| 2 | `boltz-input` | `generate_boltz_input()` | One Boltz input per homolog: chain A is the homolog, plus any co-folded ligand chains |
| 3 | `boltz-model` | `perform_boltz_modeling()` | Runs `boltz predict`, then converts each CIF to a protein-only PDB |
| 4 | `struct-homology` | `calculate_structure_homology()` | TM-align of every model against the reference → TM-score, RMSD |
| 5 | `superpose-structures` | `superpose_structures()` | PyMOL Cα superposition onto the reference, for a consistent coordinate system |
| 6 | `ref-bindingsite` | `extract_reference_binding_site()` | Reference residues within `--distance-cutoff` (default 6 Å) of the reference ligand |
| 7 | `homolog-bindingsites` | `extract_homolog_binding_sites()` | The same cut applied to each superposed homolog |
| 8 | `superpose-bindingsites` | `superpose_binding_sites()` | Binding-site-local superposition → binding-site RMSD, aligned Cα count |
| 9 | `analyze-bindingsites` | `analyze_binding_sites()` | Binding-site **metasequences** (residues within `--tolerated-misalignment`, default 1 Å), BLOSUM62 similarity, pLDDT statistics |
| 10 | `ref-cavity` | `analyze_reference_cavity_properties()` | pyKVFinder cavity detection on the reference binding site |
| 11 | `cavity-properties` | `analyze_cavity_properties()` | pyKVFinder per homolog + Open3D ICP cavity-geometry comparison |

Outputs follow the same numbered scheme as HomoLogic:

```
01_sequence_homology_results/    07_homolog_bindingsites/
02_boltz_input/                  08_homolog_bindingsites_superposed/
03_boltz_results/                09_bindingsite_metasequences/
04_structure_models/             10_bindingsite_similarity_results/
05_structure_homology_results/   11_detected_cavities/
06_superposed_structures/        12_cavity_analysis_results/
```

Run the whole pipeline, a subrange, or one stage:

```
python homologic2.py run-all -rf reference.fasta -rp reference.pdb \
    -rl reference_ligands.pdb -hf homologs.fasta --peptide KAFVQWLIAG --ccd CA

python homologic2.py run-all --dry-run                            # list stages, run nothing
python homologic2.py run-all -sa boltz-model -so struct-homology   # a subrange
python homologic2.py seq-homology -rf reference.fasta -hf homologs.fasta
```

---

# Co-folding entities

Every entity type Boltz accepts. Flags are **repeatable and order-preserving** — each adds one chain after the homolog's chain A, in the order given.

| Flag | Emits | Use for |
|----|----|----|
| `--smiles SMI` | `>X\|smiles` | small molecules |
| `--ccd CODE` | `>X\|ccd` | ions and standard cofactors — `--ccd CA` for calcium |
| `--peptide SEQ` | `>X\|protein\|empty` | peptide ligands, as a real protein chain |
| `--peptide-fasta FILE` | one chain per record | many or varying peptides |
| `--dna SEQ` / `--rna SEQ` | `>X\|dna` / `>X\|rna` | nucleic acids |
| `--entity TYPE:VALUE` | any of the five | escape hatch |
| `--copies N` | repeats the preceding entity | stoichiometry |

A peptide ligand is emitted as a **real protein chain** with MSA lookup disabled (`|empty`) — it no longer has to be hand-encoded as one long SMILES string, and ions are given by CCD code rather than as e.g. `[Ca+2]`.

SMILES (via RDKit), CCD codes (against `~/.boltz/mols`) and sequence alphabets are **validated before any file is written**, so a typo costs seconds instead of surfacing as failed GPU runs.

For constraints, templates, modifications or affinity — none of which the FASTA format can express — use YAML:

```
python homologic2.py boltz-input -hf homologs.fasta \
    --format yaml --yaml-template pocket.yaml
```

The template's target chain (default `A`) is replaced by each homolog in turn; everything else in it is carried through untouched.

---

# MSA reuse

Hand Boltz a precomputed MSA instead of letting it call the MSA server.

```
# ONE alignment containing the homologs being run
python homologic2.py boltz-input -hf homologs.fasta --use-custom-msa all.afa

# or a directory: <homolog_id>.<ext> per homolog
python homologic2.py boltz-input -hf homologs.fasta --use-custom-msa msas/

# or reuse whatever a previous run already computed (read-only)
python homologic2.py boltz-input -hf homologs.fasta --reuse-msa-from ../old_run/
```

Accepts **a3m, Boltz csv, aligned FASTA, Stockholm and Clustal**, all converted to Boltz csv. Boltz itself reads only `.a3m` and `.csv`.

## How one alignment serves every homolog

Boltz indexes an MSA row by the query chain's **residue index**, so the query row must be ungapped and exactly the chain's length. Reordering a shared alignment is therefore not sufficient: for each homolog the alignment is additionally **projected onto that homolog's ungapped columns**, with other sequences' residues in its gap columns re-emitted as a3m lowercase insertions.

Each homolog thus becomes the query of its own correctly-registered MSA, derived from the one shared file.

## Validation

Validation runs **before any stage executes**, and only when a Boltz stage is actually in the plan — `seq-homology` alone never requires an MSA. It refuses to start if:

- a homolog is absent from the alignment
- a non-query-centric alignment is ragged
- `--peptide-msa` is combined with a custom MSA (Boltz cannot mix custom and auto-generated MSAs in one input)

This matters because Boltz checks only that the MSA file **exists** and then consumes columns positionally — a mismatched alignment would otherwise produce a confident, meaningless prediction with no warning.

## Depth matters

An alignment of only the homolog set is **shallow**. Running 8 homologs gives each one a 6–8 sequence MSA, where the MSA server returns thousands, and Boltz confidence falls off steeply with depth.

Measured on the same two receptors, same peptide and Ca²⁺ ligand:

| MSA | `confidence_score` |
|----|----|
| server MSA | 0.82 – 0.91 |
| 6-sequence shared alignment | 0.40 – 0.61 |

The mechanism is not at fault — Boltz loads the custom MSA correctly (`processed/msa/*.npz` holds 6 sequences × 396 residues). The MSA is simply too thin to inform the prediction.

`--use-custom-msa` therefore **warns** when a converted MSA has fewer than 32 sequences. Use a deep alignment, or `--reuse-msa-from` / `--allow-msa-bootstrap` to reuse genuine server MSAs. A shallow alignment remains the right choice when a narrow, family-specific MSA is what is wanted — the warning exists so that it is a decision rather than an accident.

---

# Ensemble generation

`boltz-sweep` generates several seeds × many diffusion samples per homolog.

```
python homologic2.py boltz-sweep -i 02_boltz_input -o 03_sweep \
    --sweep-seeds 5 --diffusion-samples 100 --dry-run
```

This stage is **terminal**. It is not part of `run-all` and feeds none of stages 4–11, each of which assumes exactly one model per homolog. It exists to produce models for later analysis.

Two hard guards, both to avoid overloading the public MSA server:

1. **At most 5 homologs** per sweep.
2. **Every homolog needs a precomputed MSA.** `--allow-msa-bootstrap` relaxes this to Boltz building each homolog's MSA **once**, on its first seed, reused by the remaining seeds — one server call per homolog, never one per seed.

A sweep also reports what it would generate and refuses to start without `--yes`.

Outputs:

```
<out>/<homolog>/seed<N>/boltz_results_<homolog>/
<out>/<homolog>/seed<N>/<homolog>.fasta
<out>/ensemble/<homolog>/<homolog>_seed<N>_model_<M>.cif
<out>/manifest.jsonl
<out>/logs/<homolog>.seed<N>.boltz.log
```

Each seed's input file is kept deliberately — it records which MSA that seed used, the one thing not reconstructable from the output afterwards.

`manifest.jsonl` holds one row per model: homolog, seed, model index, paths, and the confidence metrics **Boltz itself wrote** (`confidence_score`, `ptm`, `iptm`, `complex_plddt`, `chains_ptm`). Nothing is recomputed. It is rebuilt from disk at the end of every sweep, so a resumed run still yields a complete index.

```
# best model per homolog
jq -s 'group_by(.homolog)[] | max_by(.confidence_score)' 03_sweep/manifest.jsonl
```

---

# Re-runs and safety

| Flag | Description |
|----|----|
| `--overwrite` | **Required** before a stage will delete an existing output folder, so re-running `boltz-model` cannot silently discard finished predictions |
| `--resume` | Reuse existing folders and skip completed work. For a sweep, a seed counts as done only when **every** diffusion sample has both its `.cif` and its confidence JSON — a seed interrupted partway is re-run, not left short |
| `--log FILE` | Run log as `timestamp LEVEL message` (default `homologic_run.log`); pass `""` to disable |

Boltz's stdout and stderr go to a per-homolog log file, and the tail is echoed on failure — an exit code alone is not enough to debug from.

---

# External executables

`boltz`, `mmseqs` and `TMalign` are resolved in this order:

1. an explicit `--boltz-bin` / `--mmseqs-bin` / `--tmalign-bin`
2. **this interpreter's own environment**
3. `$PATH`

Step 2 guarantees the Boltz that runs is the one in the same environment as the script, rather than whichever install sits first on `$PATH`. The resolved absolute path and version are logged at stage start, so any run stays traceable.

`--boltz-arg` passes anything else through to `boltz predict`. Use the `=` form so argparse does not read the value as one of HomoLogic 2's own flags:

```
--boltz-arg=--output_format=pdb
```

Commonly changed options are exposed directly: `--devices`, `--seed`, `--diffusion-samples`, `--recycling-steps`, `--use-potentials`, `--no-msa-server`.

---

# Example Workflow

A runnable, offline example (synthetic sequences, no GPU required) is in [`examples/`](examples/):

```
./examples/run_example.sh
```

For full worked examples on real data, see
https://github.com/Buller-Lab/Ssal-KRED_orthologs

## Command line

```bash
# once
python homologic2.py --install && conda activate homologic

# the whole pipeline, co-folding a peptide ligand and a calcium ion,
# reusing one alignment of the homolog set as the MSA
python homologic2.py run-all \
    -rf reference.fasta -rp reference.pdb -rl reference_ligands.pdb \
    -hf homologs.fasta \
    --peptide KAFVQWLIAG --ccd CA \
    --use-custom-msa homologs_aligned.afa \
    --overwrite

# an ensemble for a handful of homologs
python homologic2.py boltz-sweep -i 02_boltz_input -o 03_sweep \
    --sweep-seeds 5 --diffusion-samples 100 --yes
```

## Python

The class is unchanged, so the original driver style still works:

```python
from homologic2 import HomoLogic

hl = HomoLogic(reference_fasta="reference.fasta",
               reference_pdb="reference.pdb",
               reference_ligands_pdb="reference_ligands.pdb",
               homologs_fasta="homologs.fasta")

hl.calculate_sequence_homology(output_folder='01_sequence_homology_results')
hl.generate_boltz_input(
    entities=[{'type': 'protein', 'value': 'KAFVQWLIAG', 'msa': 'empty', 'copies': 1},
              {'type': 'ccd', 'value': 'CA', 'msa': None, 'copies': 1}],
    output_folder='02_boltz_input')
hl.perform_boltz_modeling(input_folder='02_boltz_input',
                          boltz_results_folder='03_boltz_results',
                          protein_only_folder='04_structure_models')
# ... stages 4-11 exactly as in HomoLogic
```

---

# Dependencies

All of these are installed and verified by `--install` / `--check`.

Core tools:

* MMseqs2
* Boltz-2
* TM-align
* PyMOL
* pyKVFinder

Python libraries:

* Biopython
* scikit-bio
* MDTraj
* Open3D
* RDKit
* NumPy
* Pandas
* PyTorch (CUDA 12.8)

For exact versions see `homologic.yml`, or run `python homologic2.py --check`.

---

# Acknowledgements

HomoLogic 2 is a command-line packaging of **HomoLogic**, whose pipeline,
method design and original implementation are by the Buller Lab — see
`homologic.py` and the project history. The analysis stages are unchanged;
this adds installation, tooling and reproducibility around them.

---

# Citation

If you use HomoLogic in your research, please cite:

```
Stockinger et al. Iterative ortholog mining enables data-driven discovery of stereoselective ketoreductases
```

(Manuscript in preparation)

---
