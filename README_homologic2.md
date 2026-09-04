# homologic2.py

Self-contained command-line driver for the HomoLogic homolog-structure
analysis pipeline: **sequence homology → Boltz-2 co-folding → structural
homology → binding-site extraction/superposition → binding-site
similarity/pLDDT → cavity/ICP analysis**.

One file. It installs its own environment, pins its own tools, and runs
either a single stage or the whole pipeline.

---

## Relationship to the other scripts here

| file | what it is |
|---|---|
| `homologic.py` | the original `HomoLogic` class, no CLI. Needs a hand-written `run_homologic.py` driver. Kept because existing run directories import it. |
| `homologic_gs.py` | the same class plus an argparse CLI. Kept for reference. |
| **`homologic2.py`** | **this.** Same pipeline, plus `--install`, pinned executables, full Boltz entity support, MSA reuse, an ensemble sweep, and safe re-runs. |

The `HomoLogic` class in `homologic2.py` is behaviourally identical to the
one in the other two for all eleven pipeline stages — stages 1, 4 and 6 were
checked to reproduce an existing run's outputs byte-for-byte. It is also
importable (`from homologic2 import HomoLogic`) without pulling in the heavy
scientific stack.

---

## Install

```bash
python homologic2.py --install            # build the `homologic` conda env
conda activate homologic
python homologic2.py --check              # verify it
```

| flag | effect |
|---|---|
| `--install` | Creates the env from a curated conda + pip spec. Resolves on any linux-64 machine. |
| `--install --pinned` | Reproduces the reference env exactly from `homologic.yml` (explicit build strings; this machine's arch/CUDA only). |
| `--install --force` | Remove and rebuild an existing env. Without it, an existing env is never touched. |
| `--install --dry-run-install` | Print the exact commands, run nothing. |
| `--install --env-name NAME` | Build somewhere other than `homologic`. |
| `--check` | Import every required module, locate every binary, report CUDA. Non-zero exit on any failure. |
| `--install-help` | Print the equivalent manual commands. |

`--install`, `--check` and `--help` all run under a bare interpreter, before
the environment exists — the heavy imports are deferred to
`_load_science_stack()`, which only runs for an actual pipeline stage.

**torch is installed before boltz**, from the CUDA 12.8 wheel index. Letting
boltz pull torch itself can land a CPU-only or pre-cu128 wheel, which cannot
drive an sm_120 (RTX 50-series) card.

Boltz downloads ~7 GB of weights into `~/.boltz` on first use.

---

## Help

`-h` and `--help` deliberately differ:

* **`-h`** — short. Usage line only.
* **`--help`** — full manual: the module docstring, then every option with its
  default. Works per stage too (`homologic2.py boltz-input --help`).

---

## The pipeline

Eleven stages, each a subcommand and each one method on `HomoLogic`. Run
individually, or all together with `run-all` (which uses a fixed, internally
consistent `01_…12_` folder chain).

| # | subcommand | method | what it does |
|---|---|---|---|
| 1 | `seq-homology` | `calculate_sequence_homology()` | MMseqs2 `easy-search` of the reference against the homologs, with a scikit-bio global-alignment fallback for anything MMseqs2 reports no usable identity for. → `sequence_homology.csv` |
| 2 | `boltz-input` | `generate_boltz_input()` | Splits the homologs FASTA into one Boltz input per homolog, chain A = the homolog, plus any co-folded ligand chains. |
| 3 | `boltz-model` | `perform_boltz_modeling()` | Runs `boltz predict` on each input, then converts each CIF to a protein-only PDB via mdtraj. |
| 4 | `struct-homology` | `calculate_structure_homology()` | TM-aligns every predicted model against the reference, appending `rmsd`/`tm_score` to the stage-1 CSV. |
| 5 | `superpose-structures` | `superpose_structures()` | PyMOL CA superposition of every model onto the reference, so all later geometry shares one frame. |
| 6 | `ref-bindingsite` | `extract_reference_binding_site()` | Reference residues with any atom within `--distance-cutoff` of the reference ligand. → `<reference>_bindingsite.pdb` |
| 7 | `homolog-bindingsites` | `extract_homolog_binding_sites()` | The same cut applied to each superposed homolog, using the reference ligand as the fixed frame. |
| 8 | `superpose-bindingsites` | `superpose_binding_sites()` | A second, binding-site-local PyMOL superposition, with quality gates (`--min-residues-total`, `--min-aligned-atoms`); too few aligned atoms → RMSD reported as NaN rather than trusted. |
| 9 | `analyze-bindingsites` | `analyze_binding_sites()` | Builds each homolog's binding-site "metasequence" (residues spatially matched to the reference within `--tolerated-misalignment`), scores its identity, and summarises pLDDT from Boltz's `.npz`. |
| 10 | `ref-cavity` | `analyze_reference_cavity_properties()` | pyKVFinder cavity detection on the reference binding site. |
| 11 | `cavity-properties` | `analyze_cavity_properties()` | pyKVFinder per homolog plus Open3D ICP shape comparison against the reference cavity. |

```bash
# whole pipeline
python homologic2.py run-all -rf reference.fasta -rp reference.pdb \
    -rl reference_ligands.pdb -hf homologs.fasta \
    --peptide KAFVQWLIAG --ccd CA

python homologic2.py run-all --dry-run          # list stages, run nothing
python homologic2.py run-all -sa boltz-model -so struct-homology   # a subrange
```

---

## Co-folding entities (stage 2)

Every entity type Boltz accepts. Flags are **repeatable and
order-preserving** — each adds one chain after the homolog's chain A, in the
order typed.

| flag | emits | use for |
|---|---|---|
| `--smiles SMI` | `>X\|smiles` | small molecules |
| `--ccd CODE` | `>X\|ccd` | ions and standard cofactors — `--ccd CA` for calcium |
| `--peptide SEQ` | `>X\|protein\|empty` | peptide ligands, as a real protein chain |
| `--peptide-fasta FILE` | one chain per record | many/varying peptides |
| `--dna SEQ` / `--rna SEQ` | `>X\|dna` / `>X\|rna` | nucleic acids |
| `--entity TYPE:VALUE` | any of the five | escape hatch |
| `--copies N` | repeats the preceding entity | stoichiometry |

`|empty` on a peptide means "no MSA lookup", which is what a short peptide
wants and is also what keeps a custom MSA on chain A legal (Boltz refuses to
mix custom and auto-generated MSAs in one input).

Everything is **validated before any file is written**: SMILES through RDKit,
CCD codes against `~/.boltz/mols`, sequences against their alphabet. A typo
costs seconds instead of surfacing as N failed GPU runs.

For constraints, templates, modifications or affinity — none of which the
FASTA format can express — use YAML:

```bash
python homologic2.py boltz-input -hf homologs.fasta \
    --format yaml --yaml-template pocket.yaml
```

Your template's target chain (default `A`) is replaced by each homolog in
turn; everything else in it is carried through untouched.

---

## MSA reuse

Hand Boltz a precomputed MSA instead of letting it call the MSA server.

```bash
# ONE alignment containing the homologs being run
python homologic2.py boltz-input -hf homologs.fasta --use-custom-msa all.afa

# or a directory: <homolog_id>.<ext> per homolog
python homologic2.py boltz-input -hf homologs.fasta --use-custom-msa msas/

# or pick up whatever a previous run already computed (read-only)
python homologic2.py boltz-input -hf homologs.fasta --reuse-msa-from ../old_run/
```

Accepts **a3m, Boltz csv, aligned FASTA, Stockholm, Clustal** — all converted
to Boltz csv. (Boltz itself reads only `.a3m` and `.csv`.)

**How one alignment serves every homolog.** Boltz indexes an MSA row by the
query chain's *residue index*, so the query row must be ungapped and exactly
the chain's length. Reordering the alignment is therefore not enough: for each
homolog the alignment is also **projected onto that homolog's ungapped
columns**, with other sequences' residues in its gap columns re-emitted as a3m
lowercase insertions. Each homolog thus becomes the query of its own
correctly-registered MSA, derived from the shared file.

**Validation runs before any stage executes**, and only when a Boltz stage is
actually in the plan (`seq-homology` alone never demands an MSA). It errors
if a homolog is absent from the alignment, if a non-query-centric alignment is
ragged, or if `--peptide-msa` is combined with a custom MSA. This matters:
Boltz checks only that the file *exists* and then consumes columns
positionally, so a mismatched MSA would produce a confident, meaningless
prediction with no warning.

Converted MSAs are cached beside the stage's output folder, as
`<parent>/msa/<homolog_id>.csv`, and referenced by absolute path.

### Depth matters — read this before using a shared alignment

An alignment of just the homolog set is **shallow**. Running 8 homologs gives
each one a 6–8 sequence MSA, where the MSA server returns thousands. Boltz
confidence falls off steeply with depth.

Measured on this machine, same two receptors, same peptide + Ca²⁺:

| MSA | `confidence_score` |
|---|---|
| server MSA (earlier sweeps) | 0.82 – 0.91 |
| 6-sequence shared alignment | 0.40 – 0.61 |

The mechanism is not at fault — Boltz loads the custom MSA correctly (verified:
`processed/msa/*.npz` holds 6 sequences × 396 residues). The MSA is simply too
thin to inform the prediction.

`--use-custom-msa` therefore warns when any converted MSA has fewer than 32
sequences. Use it with a **deep** alignment (e.g. one built from a real
database search over each homolog), or use `--reuse-msa-from` /
`--allow-msa-bootstrap` to reuse genuine server MSAs. A shallow alignment is
still the right choice when you deliberately want a narrow, family-specific
MSA — the warning exists so that is a decision, not an accident.

---

## `boltz-sweep` — ensembles

Generates several seeds × many diffusion samples per homolog.

```bash
python homologic2.py boltz-sweep -i 02_boltz_input -o 03_sweep \
    --sweep-seeds 5 --diffusion-samples 100 --dry-run
```

**This stage is terminal.** It is not part of `run-all` and feeds none of
stages 4–11, each of which assumes exactly one model per homolog. It exists to
produce models for later analysis, nothing more.

Two hard guards, both about not overloading the public MSA server:

1. **At most 5 homologs** per sweep.
2. **Every homolog needs a precomputed MSA.** `--allow-msa-bootstrap` relaxes
   this to Boltz building each homolog's MSA *once*, on its first seed, reused
   by the remaining seeds — one server call per homolog, never one per seed.

Plus a cost gate: the sweep reports what it would generate and refuses to
start without `--yes`. `--dry-run` previews without running.

Output:

```
<out>/<homolog>/seed<N>/boltz_results_<homolog>/…   # raw Boltz output
<out>/<homolog>/seed<N>/<homolog>.fasta             # that seed's exact input
<out>/ensemble/<homolog>/<homolog>_seed<N>_model_<M>.cif
<out>/manifest.jsonl
<out>/logs/<homolog>.seed<N>.boltz.log
```

Each seed's input file is kept deliberately: it records which MSA that seed
used, the one thing not reconstructable from the output afterwards.

`manifest.jsonl` is one row per model — homolog, seed, model index, paths, and
the confidence metrics **Boltz itself wrote** (`confidence_score`, `ptm`,
`iptm`, `complex_plddt`, `chains_ptm`, …). Nothing is recomputed; deriving
per-chain pLDDT/PAE from the token-level `.npz` arrays is separate analysis
and is not done here. It is rebuilt from disk at the end of every sweep, so a
resumed run still yields a complete index.

```bash
# best model per homolog
jq -s 'group_by(.homolog)[] | max_by(.confidence_score)' 03_sweep/manifest.jsonl
```

---

## Re-runs and safety

| flag | effect |
|---|---|
| `--overwrite` | **Required** before a stage will delete an existing output folder. Without it, re-running `boltz-model` cannot throw away finished predictions. |
| `--resume` | Reuse existing folders and skip work already done. For a sweep, a seed counts as done only when *every* diffusion sample has both its `.cif` and its confidence JSON — a seed interrupted partway is re-run, not silently left short. |
| `--log FILE` | Run log as `timestamp LEVEL message` (default `homologic_run.log`). Pass `""` to disable. |

Boltz's own stdout/stderr goes to a per-homolog log file, and the tail is
echoed on failure — an exit code alone is not enough to debug from.

---

## Executables

`boltz`, `mmseqs` and `TMalign` are resolved in this order:

1. an explicit `--boltz-bin` / `--mmseqs-bin` / `--tmalign-bin`,
2. **this interpreter's own environment**,
3. `$PATH`.

Step 2 is the important one: it guarantees the `boltz` that runs is the one in
the same env as the script, rather than whichever install happens to sit first
on `$PATH`. The resolved absolute path and version are logged at stage start,
so any run is traceable afterwards.

`--boltz-arg` passes anything through to `boltz predict`. Use the `=` form so
argparse does not read the value as one of ours:

```bash
--boltz-arg=--output_format=pdb
```

---

## Gotchas

* **Stage 2 needs stage 1's homologs FASTA**, not its CSV.
* **`--use-custom-msa` + `--peptide-msa` is rejected** — Boltz cannot mix
  custom and auto-generated MSAs in one input.
* **A sweep of 5 homologs × 5 seeds × 100 samples is 2500 structures** across
  25 `boltz predict` calls. Check `--dry-run` first.
* **`homologic2.py` is importable; `homologic2.0.py` was not** — `2.0` is not a
  valid Python module name. That is why this file is named as it is.
