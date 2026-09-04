#!/usr/bin/env python3
"""
Script written for HomoLogic v2
by gargs on 04.09.26
for bugs contact: sriram.garg@unibe.ch

homologic2.py
================
Self-contained command-line driver for the HomoLogic homolog-structure-
analysis pipeline (sequence homology -> Boltz-2 co-folding -> structural
homology -> binding-site extraction/superposition -> binding-site
similarity/pLDDT analysis -> cavity/ICP analysis).

Supersedes the homologic.py / homologic_gs.py pair kept alongside it:

  homologic.py     the HomoLogic class only, no CLI -- needs a hand-written
                   run_homologic.py driver. Kept because existing run
                   directories import it.
  homologic_gs.py  the same class plus an argparse CLI. Kept for reference.

What this version adds over homologic_gs.py:

  --install     builds the `homologic` conda environment and every tool the
                pipeline shells out to, so the script bootstraps itself.
  --check       verifies an existing environment (imports + binaries).
  pinned tools  boltz/mmseqs/TMalign are resolved to THIS interpreter's own
                environment and logged, instead of taking whatever `$PATH`
                happens to offer first.
  full Boltz    every entity type Boltz accepts -- protein (peptide ligands),
                dna, rna, ccd, smiles -- plus a YAML passthrough for
                constraints/templates/affinity. The old code could only emit
                up to four SMILES chains.
  safe re-runs  --overwrite is now required to wipe an existing output
                folder, and --resume skips homologs already predicted.

USAGE
-----
    python homologic2.py --install            # once, builds the env
    conda activate homologic
    python homologic2.py <stage> [options]
    python homologic2.py <stage> --help       # per-stage options

Pipeline stages, in the order they are normally run:
    seq-homology          calculate_sequence_homology()
    boltz-input            generate_boltz_input()
    boltz-model              perform_boltz_modeling()
    struct-homology            calculate_structure_homology()
    superpose-structures         superpose_structures()
    ref-bindingsite                extract_reference_binding_site()
    homolog-bindingsites             extract_homolog_binding_sites()
    superpose-bindingsites              superpose_binding_sites()
    analyze-bindingsites                   analyze_binding_sites()
    ref-cavity                               analyze_reference_cavity_properties()
    cavity-properties                          analyze_cavity_properties()

    run-all   runs all eleven stages above in sequence (or a --start-at /
              --stop-after subrange) using a fixed, internally consistent
              numbered-folder scheme (01_... through 12_...).

EXTERNAL TOOLS REQUIRED
-----------------------
    mmseqs   -> seq-homology            TMalign -> struct-homology
    boltz    -> boltz-model             PyMOL   -> the two superpose stages

All of these are installed by --install and verified by --check.
"""

# --- Standard library ONLY at module scope ----------------------------
# The heavy scientific stack is imported lazily by _load_science_stack(),
# because --install / --install-help / --check must run under a bare
# interpreter BEFORE the environment they build exists.
import argparse
import csv
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# Names filled in by _load_science_stack(); declared here so the pipeline
# method bodies below can reference them as plain globals, unchanged.
pd = np = md = o3d = pyKVFinder = None
PDBParser = NeighborSearch = PDBIO = is_aa = None
substitution_matrices = None
global_pairwise_align_protein = None
Protein = io = None

SCIENCE_MODULES = [
    # (import name, human label, how --install provides it)
    ('pandas', 'pandas', 'conda'),
    ('numpy', 'numpy', 'conda'),
    ('mdtraj', 'mdtraj', 'conda'),
    ('Bio', 'biopython', 'pip'),
    ('skbio', 'scikit-bio', 'conda'),
    ('pyKVFinder', 'pyKVFinder', 'pip'),
    ('open3d', 'open3d', 'pip'),
    ('pymol', 'pymol-open-source', 'conda'),
    ('boltz', 'boltz', 'pip'),
    ('torch', 'torch', 'pip'),
    ('rdkit', 'rdkit', 'pip'),
]

REQUIRED_BINARIES = ['boltz', 'mmseqs', 'TMalign']

ENV_NAME_DEFAULT = 'homologic'


def _load_science_stack():
    """
    Import the heavy dependencies into module globals.

    Called once, only by commands that actually run a pipeline stage --
    never by --install/--check/--help, which have to work in an interpreter
    where none of this is installed yet. Fails with one actionable message
    naming the missing module rather than a bare traceback.
    """
    global pd, np, md, o3d, pyKVFinder
    global PDBParser, NeighborSearch, PDBIO, is_aa, substitution_matrices
    global global_pairwise_align_protein, Protein, io

    missing = []
    try:
        import pandas as _pd; pd = _pd
        import numpy as _np; np = _np
    except ImportError as e:
        missing.append(str(e))
    try:
        import mdtraj as _md; md = _md
    except ImportError as e:
        missing.append(str(e))
    try:
        from Bio.PDB import PDBParser as _P, NeighborSearch as _N, PDBIO as _I, is_aa as _a
        from Bio.Align import substitution_matrices as _s
        PDBParser, NeighborSearch, PDBIO, is_aa, substitution_matrices = _P, _N, _I, _a, _s
    except ImportError as e:
        missing.append(str(e))
    try:
        import pyKVFinder as _k; pyKVFinder = _k
    except ImportError as e:
        missing.append(str(e))
    try:
        import open3d as _o; o3d = _o
    except ImportError as e:
        missing.append(str(e))
    try:
        from skbio.alignment import global_pairwise_align_protein as _g
        from skbio import Protein as _Pr, io as _io
        global_pairwise_align_protein, Protein, io = _g, _Pr, _io
    except ImportError as e:
        missing.append(str(e))

    if missing:
        sys.stderr.write("ERROR: the analysis environment is incomplete:\n")
        for m in missing:
            sys.stderr.write(f"    {m}\n")
        sys.stderr.write(
            "\nBuild or repair it with:\n"
            f"    python {os.path.basename(__file__)} --install\n"
            f"    python {os.path.basename(__file__)} --check\n")
        sys.exit(1)


def hms_string(sec_elapsed):
    """Computes a human-readable hms string."""
    h = int(sec_elapsed / (60 * 60))
    m = int((sec_elapsed % (60 * 60)) / 60)
    s = sec_elapsed % 60.
    return "{}:{:>02}:{:>05.2f}".format(h, m, s)


# ======================================================================
# Executable resolution
# ======================================================================

def resolve_exe(name, override=None, required=True):
    """
    Resolve an external binary to an absolute path, deterministically.

    Order:
      1. an explicit --<tool>-bin override,
      2. THIS interpreter's own environment (sys.executable's directory) --
         the important one: it guarantees `boltz` comes from the same env as
         the running script instead of whichever install happens to sit
         first on $PATH,
      3. $PATH,
      4. error (or None when required=False).
    """
    if override:
        p = Path(override).expanduser()
        if p.is_file() and os.access(p, os.X_OK):
            return str(p.resolve())
        sys.exit(f"ERROR: --{name.lower()}-bin '{override}' is not an executable file.")

    sibling = Path(sys.executable).resolve().parent / name
    if sibling.is_file() and os.access(sibling, os.X_OK):
        return str(sibling)

    found = shutil.which(name)
    if found:
        return str(Path(found).resolve())

    if required:
        sys.exit(
            f"ERROR: required executable '{name}' not found in this environment "
            f"({Path(sys.executable).parent}) or on $PATH.\n"
            f"       Build the environment with:  python {os.path.basename(__file__)} --install")
    return None


def exe_version(path):
    """Best-effort one-line version string for an external binary, for the log."""
    for flag in ('--version', '-version', '-h'):
        try:
            r = subprocess.run([path, flag], capture_output=True, text=True, timeout=60)
            out = (r.stdout or r.stderr).strip().splitlines()
            if out:
                return out[0].strip()
        except (OSError, subprocess.SubprocessError):
            continue
    return 'unknown version'


# ======================================================================
# Environment installation
# ======================================================================

CONDA_CHANNELS = ['conda-forge', 'bioconda']

# Curated spec: what the pipeline actually needs, without the explicit build
# strings of a full env dump, so it resolves on a machine that is not this
# one. Use --pinned for byte-identical reproduction of the reference env.
CONDA_SPEC = [
    'python=3.11',
    'mmseqs2',
    'tmalign',
    'pymol-open-source',
    'mdtraj',
    'openmm',
    'scikit-bio',
    'h5py',
    'pytables',
    'pandas',
    'numpy=1.26',
    'matplotlib-base',
    'seaborn',
    'openpyxl',
    'natsort',
    'pip',
]

# torch is installed first, from the CUDA 12.8 wheel index: boltz would
# otherwise pull whatever torch pip defaults to, and a CPU-only or pre-cu128
# wheel cannot drive an sm_120 card (RTX 50-series). Everything else is
# resolved normally afterwards.
PIP_TORCH = ['torch==2.10.0', '--index-url', 'https://download.pytorch.org/whl/cu128']
PIP_SPEC = [
    'boltz==2.2.1',
    'pyKVFinder',
    'open3d',
    'biopython==1.84',
    'rdkit',
]

PINNED_YML = 'homologic.yml'


def find_conda():
    """
    Locate a conda-compatible package manager, preferring the fast solvers.
    Returns (executable_path, kind) where kind is 'mamba' or 'conda'.
    """
    for name in ('mamba', 'micromamba'):
        found = shutil.which(name)
        if found:
            return found, 'mamba'
    found = shutil.which('conda') or os.environ.get('CONDA_EXE')
    if found and Path(found).is_file():
        return str(found), 'conda'
    sys.exit(
        "ERROR: no conda/mamba found on $PATH.\n"
        "       Install Miniforge or Miniconda first:\n"
        "       https://github.com/conda-forge/miniforge#install")


def conda_env_exists(conda, env_name):
    """True if an env of this name is already registered with conda."""
    try:
        r = subprocess.run([conda, 'env', 'list', '--json'],
                           capture_output=True, text=True, check=True)
        envs = json.loads(r.stdout).get('envs', [])
        return any(Path(e).name == env_name for e in envs)
    except (subprocess.SubprocessError, ValueError, OSError):
        return False


def env_python(conda, env_name):
    """Absolute path to the target env's interpreter."""
    try:
        r = subprocess.run([conda, 'env', 'list', '--json'],
                           capture_output=True, text=True, check=True)
        for e in json.loads(r.stdout).get('envs', []):
            if Path(e).name == env_name:
                return str(Path(e) / 'bin' / 'python')
    except (subprocess.SubprocessError, ValueError, OSError):
        pass
    return None


def run_step(cmd, dry_run=False, cwd=None):
    """Echo a command, then run it (unless --dry-run). Aborts on failure."""
    printable = ' '.join(str(c) for c in cmd)
    print(f"    $ {printable}")
    if dry_run:
        return
    r = subprocess.run([str(c) for c in cmd], cwd=cwd)
    if r.returncode != 0:
        sys.exit(f"ERROR: command failed (exit {r.returncode}): {printable}")


def do_install(args):
    """
    Build the analysis environment.

    Default path installs the curated CONDA_SPEC/PIP_SPEC above. --pinned
    instead replays homologic.yml (a full explicit-build dump sitting next to
    this script), which reproduces the reference environment exactly but only
    works on linux-64 with a matching CUDA stack.
    """
    conda, kind = find_conda()
    env_name = args.env_name
    print(f"HomoLogic v2 install")
    print(f"  package manager : {conda}  ({kind})")
    print(f"  target env      : {env_name}")
    print(f"  mode            : {'pinned (homologic.yml)' if args.pinned else 'curated spec'}")
    print()

    if conda_env_exists(conda, env_name):
        if not args.force:
            sys.exit(
                f"ERROR: conda env '{env_name}' already exists.\n"
                f"       Verify it with :  python {os.path.basename(__file__)} --check\n"
                f"       Rebuild it with:  --install --force\n"
                f"       Or build elsewhere: --install --env-name {env_name}_new")
        print(f"--force given: removing existing env '{env_name}'")
        run_step([conda, 'env', 'remove', '-n', env_name, '-y'], args.dry_run)

    if args.pinned:
        yml = Path(__file__).resolve().parent / PINNED_YML
        if not yml.is_file():
            sys.exit(f"ERROR: --pinned needs '{yml}', which is missing.\n"
                     f"       Drop the pinned environment dump there, or install "
                     f"without --pinned.")
        print(f"[1/2] Creating env from {yml}")
        run_step([conda, 'env', 'create', '-n', env_name, '-f', str(yml)], args.dry_run)
    else:
        print(f"[1/3] Creating env with conda packages")
        cmd = [conda, 'create', '-n', env_name, '-y']
        for ch in CONDA_CHANNELS:
            cmd += ['-c', ch]
        if kind == 'conda':
            cmd += ['--solver', 'libmamba']
        cmd += CONDA_SPEC
        run_step(cmd, args.dry_run)

        py = env_python(conda, env_name) if not args.dry_run else f'<{env_name}>/bin/python'
        if not args.dry_run and not py:
            sys.exit(f"ERROR: env '{env_name}' was created but its python could not be located.")

        print(f"\n[2/3] Installing torch (CUDA 12.8 wheels)")
        run_step([py, '-m', 'pip', 'install'] + PIP_TORCH, args.dry_run)

        print(f"\n[3/3] Installing remaining pip packages")
        run_step([py, '-m', 'pip', 'install'] + PIP_SPEC, args.dry_run)

    if args.dry_run:
        print("\n--dry-run: nothing was executed.")
        return 0

    print("\nInstall finished. Verifying:\n")
    rc = do_check(argparse.Namespace(env_name=env_name, check=True))
    if rc == 0:
        print(f"\nReady. Activate it with:\n    conda activate {env_name}")
    return rc


def do_check(args):
    """
    Verify an environment: import every required module and locate every
    required binary, inside the TARGET env's interpreter (not this one).
    Prints a PASS/FAIL table; returns non-zero if anything failed.
    """
    env_name = getattr(args, 'env_name', None) or ENV_NAME_DEFAULT

    # Check the running interpreter if it already is the target env,
    # otherwise reach into the named env.
    if Path(sys.prefix).name == env_name:
        py = sys.executable
    else:
        conda, _ = find_conda()
        py = env_python(conda, env_name)
        if not py or not Path(py).is_file():
            sys.exit(f"ERROR: conda env '{env_name}' not found.\n"
                     f"       Build it with:  python {os.path.basename(__file__)} --install")

    print(f"Checking env '{env_name}'")
    print(f"  interpreter: {py}\n")

    probe = r'''
import importlib, json, os, shutil, sys
from pathlib import Path
mods = json.loads(sys.argv[1])
bins = json.loads(sys.argv[2])
out = {"modules": [], "binaries": [], "cuda": None}
for name, label, src in mods:
    try:
        m = importlib.import_module(name)
        v = getattr(m, "__version__", "") or ""
        out["modules"].append([label, True, str(v)])
    except Exception as e:
        out["modules"].append([label, False, type(e).__name__ + ": " + str(e)[:60]])
here = Path(sys.executable).parent
for b in bins:
    p = here / b
    if not (p.is_file() and os.access(p, os.X_OK)):
        w = shutil.which(b)
        p = Path(w) if w else None
    out["binaries"].append([b, bool(p), str(p) if p else "not found"])
try:
    import torch
    out["cuda"] = [torch.cuda.is_available(), torch.version.cuda,
                   torch.cuda.get_device_name(0) if torch.cuda.is_available() else ""]
except Exception:
    pass
print(json.dumps(out))
'''
    r = subprocess.run([py, '-c', probe, json.dumps(SCIENCE_MODULES),
                        json.dumps(REQUIRED_BINARIES)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr.strip()[-2000:], file=sys.stderr)
        sys.exit(f"ERROR: could not probe env '{env_name}'.")

    try:
        res = json.loads(r.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        print(r.stdout[-2000:], file=sys.stderr)
        sys.exit("ERROR: unreadable probe output.")

    failed = 0
    print("  python modules")
    for label, ok, info in res['modules']:
        mark = 'PASS' if ok else 'FAIL'
        failed += 0 if ok else 1
        print(f"    [{mark}] {label:<22} {info}")
    print("\n  external binaries")
    for name, ok, info in res['binaries']:
        mark = 'PASS' if ok else 'FAIL'
        failed += 0 if ok else 1
        print(f"    [{mark}] {name:<22} {info}")

    if res.get('cuda'):
        avail, cuda_ver, dev = res['cuda']
        print("\n  gpu")
        if avail:
            print(f"    [PASS] torch CUDA {cuda_ver}      {dev}")
        else:
            print(f"    [WARN] torch reports no CUDA device (cuda={cuda_ver}); "
                  f"Boltz will run on CPU and be very slow")

    if failed:
        print(f"\n{failed} check(s) FAILED. Repair with:  "
              f"python {os.path.basename(__file__)} --install --force --env-name {env_name}")
        return 1
    print(f"\nAll checks passed.")
    return 0


INSTALL_HELP = f"""
HomoLogic v2 -- manual installation
====================================
`--install` runs all of this for you; these are the same commands, for when
you would rather drive it yourself.

    # 1. Create the environment
    conda create -n {ENV_NAME_DEFAULT} -y -c conda-forge -c bioconda \\
        {' '.join(CONDA_SPEC)}

    # 2. Activate it
    conda activate {ENV_NAME_DEFAULT}

    # 3. torch first, from the CUDA 12.8 wheel index. Order matters: let
    #    boltz pull torch itself and you may land on a CPU-only or pre-cu128
    #    wheel, which cannot drive an sm_120 (RTX 50-series) card.
    pip install {' '.join(PIP_TORCH)}

    # 4. The rest
    pip install {' '.join(PIP_SPEC)}

    # 5. Verify
    python {os.path.basename(__file__)} --check

To reproduce the reference environment exactly instead (linux-64 + matching
CUDA only), use the pinned dump next to this script:

    python {os.path.basename(__file__)} --install --pinned

Boltz downloads its weights (~7 GB) into ~/.boltz on first use, and
--use_msa_server needs outbound network access.
"""


# ======================================================================
# Boltz co-folding entities
# ======================================================================
#
# Boltz's FASTA format is  >CHAIN_ID|ENTITY_TYPE[|MSA_ID]  with
# ENTITY_TYPE in {protein, dna, rna, ccd, smiles}, and MSA_ID permitted on
# proteins only (boltz/data/parse/fasta.py). The MSA_ID "empty" is a
# sentinel meaning "do not look up an MSA for this chain"
# (boltz/data/parse/schema.py) -- which is what a short peptide ligand
# wants, since an MSA-server lookup on a 10-mer is both pointless and slow.
#
# Its YAML format instead uses protein/dna/rna/ligand entities (ligand
# taking smiles XOR ccd) and additionally carries constraints, templates,
# modifications and affinity. Anything not expressible as a flag below is
# reachable through --yaml-template.

ENTITY_TYPES = ('protein', 'dna', 'rna', 'ccd', 'smiles')

PROTEIN_ALPHABET = set('ACDEFGHIKLMNPQRSTVWYXBZUO')
DNA_ALPHABET = set('ACGTN')
RNA_ALPHABET = set('ACGUN')


def chain_ids():
    """
    Yield Boltz chain IDs: A..Z, then AA, AB, ... .

    homologic_gs.py used chr(66 + j) capped at four ligands, which both
    limited the number of co-folded entities and would have emitted
    non-alphabetic characters past chain E.
    """
    import itertools
    import string
    alphabet = string.ascii_uppercase
    for c in alphabet:
        yield c
    for n in itertools.count(2):
        for combo in itertools.product(alphabet, repeat=n):
            yield ''.join(combo)


def _abbrev(value, width=28):
    """Shorten a long SMILES/sequence for a log line."""
    value = str(value)
    return value if len(value) <= width else value[:width - 3] + '...'


def tail_of(path, n=20):
    """Last n lines of a file, as one string, for error reporting."""
    try:
        with open(path, errors='replace') as fh:
            lines = fh.readlines()
        return ''.join(lines[-n:]).rstrip() or '(log empty)'
    except OSError:
        return '(log unreadable)'


def boltz_cache_dir():
    """Where Boltz keeps its CCD component pickles."""
    return Path(os.environ.get('BOLTZ_CACHE', Path.home() / '.boltz'))


def validate_entities(entities):
    """
    Cheap pre-flight validation of the co-folding entities.

    Returns a list of human-readable problems (empty when all good). Run at
    input-generation time so a typo costs seconds rather than surfacing as
    N failed GPU predictions later.
    """
    problems = []
    mols = boltz_cache_dir() / 'mols'
    for ent in entities:
        etype, value = ent['type'], str(ent['value'])
        if etype not in ENTITY_TYPES:
            problems.append(f"unknown entity type '{etype}' "
                            f"(expected one of: {', '.join(ENTITY_TYPES)})")
            continue
        if not value:
            problems.append(f"empty value for {etype} entity")
            continue

        if etype == 'smiles':
            try:
                from rdkit import Chem
                from rdkit import RDLogger
                RDLogger.DisableLog('rdApp.*')
                if Chem.MolFromSmiles(value) is None:
                    problems.append(f"SMILES not parseable by RDKit: {_abbrev(value, 60)}")
            except ImportError:
                pass  # rdkit absent: skip rather than block
        elif etype == 'ccd':
            code = value.strip().upper()
            if not re.fullmatch(r'[A-Z0-9]{1,5}', code):
                problems.append(f"CCD code '{value}' is not 1-5 alphanumeric characters")
            elif mols.is_dir() and not (mols / f'{code}.pkl').is_file():
                problems.append(
                    f"CCD code '{code}' is not in the Boltz component cache ({mols}). "
                    f"Check the code, or pass it as --smiles instead.")
        else:
            alphabet = {'protein': PROTEIN_ALPHABET,
                        'dna': DNA_ALPHABET,
                        'rna': RNA_ALPHABET}[etype]
            bad = sorted(set(value.upper()) - alphabet)
            if bad:
                problems.append(f"{etype} sequence contains invalid residue(s) "
                                f"{''.join(bad)}: {_abbrev(value, 60)}")
            if ent.get('msa') and etype != 'protein':
                problems.append(f"MSA is only allowed on protein chains, not {etype}")
    return problems


def legacy_smiles_entities(smiles_code, owner=None):
    """
    Translate homologic.py's old `smiles_code=` argument into entities, so
    existing run_homologic.py drivers keep working unchanged.
    """
    if smiles_code in (False, None):
        return []
    if smiles_code == 'SMILES':
        values = list(getattr(owner, 'default_smiles', []) or [])
    elif isinstance(smiles_code, str):
        values = [smiles_code.strip()]
    elif isinstance(smiles_code, (list, tuple)):
        values = [str(s).strip() for s in smiles_code if str(s).strip()]
    else:
        values = [str(smiles_code).strip()]
    return [{'type': 'smiles', 'value': v, 'msa': None, 'copies': 1} for v in values]


def read_fasta_sequences(path):
    """[(id, sequence)] from a FASTA file."""
    records, cur_id, cur = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if cur_id:
                    records.append((cur_id, ''.join(cur)))
                cur_id, cur = line[1:].split()[0] if len(line) > 1 else 'seq', []
            elif line:
                cur.append(line)
    if cur_id:
        records.append((cur_id, ''.join(cur)))
    return records


def build_entities(args):
    """
    Assemble the ordered co-folding entity list from the CLI.

    Flags are collected in the order they were typed (see EntityAction), so
    chain assignment after the homolog's chain A is predictable:
        --peptide KAFVQWLIAG --ccd CA
    gives chain B = the peptide, chain C = the calcium ion.
    """
    entities = []
    peptide_msa = getattr(args, 'peptide_msa', False)
    for kind, value in getattr(args, 'entities', None) or []:
        if kind == 'copies':
            if not entities:
                sys.exit("ERROR: --copies must follow an entity flag.")
            entities[-1]['copies'] = int(value)
            continue
        if kind == 'peptide-fasta':
            for _sid, seq in read_fasta_sequences(value):
                entities.append({'type': 'protein', 'value': seq,
                                 'msa': None if peptide_msa else 'empty', 'copies': 1})
            continue
        if kind == 'entity':
            if ':' not in value:
                sys.exit(f"ERROR: --entity expects TYPE:VALUE, got '{value}'")
            etype, val = value.split(':', 1)
            etype = etype.strip().lower()
            if etype not in ENTITY_TYPES:
                sys.exit(f"ERROR: --entity type '{etype}' is not one of: "
                         f"{', '.join(ENTITY_TYPES)}")
            msa = 'empty' if (etype == 'protein' and not peptide_msa) else None
            entities.append({'type': etype, 'value': val.strip(), 'msa': msa, 'copies': 1})
            continue

        etype = {'peptide': 'protein'}.get(kind, kind)
        msa = 'empty' if (etype == 'protein' and not peptide_msa) else None
        entities.append({'type': etype, 'value': value.strip(), 'msa': msa, 'copies': 1})
    return entities


class EntityAction(argparse.Action):
    """
    Append (flag_kind, value) to a single shared `entities` list so that
    entities from different flags keep the order they were typed in.
    """
    def __init__(self, option_strings, dest, kind=None, **kwargs):
        self.kind = kind
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        items = list(getattr(namespace, self.dest, None) or [])
        items.append((self.kind, values))
        setattr(namespace, self.dest, items)


def add_entity_arguments(parser):
    """The co-folding entity flags, shared by boltz-input and run-all."""
    g = parser.add_argument_group(
        'Boltz co-folding entities',
        'Repeatable and order-preserving: each flag adds one chain after the '
        'homolog (chain A), in the order given. Every entity type Boltz '
        'accepts is available.')
    g.add_argument('-sm', '--smiles', dest='entities', action=EntityAction, kind='smiles',
                   metavar='SMILES', help='Small-molecule ligand as SMILES.')
    g.add_argument('--ccd', dest='entities', action=EntityAction, kind='ccd',
                   metavar='CODE',
                   help='Ligand/ion by CCD component code, e.g. --ccd CA for calcium. '
                        'Preferred over SMILES for standard ions and cofactors.')
    g.add_argument('--peptide', dest='entities', action=EntityAction, kind='peptide',
                   metavar='SEQ',
                   help='Peptide ligand as a one-letter sequence; emitted as a real '
                        'protein chain with MSA lookup disabled.')
    g.add_argument('--peptide-fasta', dest='entities', action=EntityAction,
                   kind='peptide-fasta', metavar='FILE',
                   help='Add one peptide-ligand chain per record in this FASTA.')
    g.add_argument('--dna', dest='entities', action=EntityAction, kind='dna', metavar='SEQ',
                   help='DNA chain.')
    g.add_argument('--rna', dest='entities', action=EntityAction, kind='rna', metavar='SEQ',
                   help='RNA chain.')
    g.add_argument('--entity', dest='entities', action=EntityAction, kind='entity',
                   metavar='TYPE:VALUE',
                   help='Escape hatch: any Boltz entity as TYPE:VALUE, where TYPE is '
                        'one of protein, dna, rna, ccd, smiles.')
    g.add_argument('--copies', dest='entities', action=EntityAction, kind='copies',
                   metavar='N', help='Number of copies of the entity flag that precedes it.')
    g.add_argument('--peptide-msa', action='store_true',
                   help='Do run MSA lookup for peptide/protein ligand chains '
                        '(default: skipped, which is what a short peptide wants).')
    g.add_argument('--format', choices=('fasta', 'yaml'), default='fasta',
                   help='Boltz input format to emit.')
    g.add_argument('--yaml-template', metavar='FILE',
                   help='With --format yaml: a Boltz YAML whose target chain is replaced '
                        'by each homolog. The only way to use constraints, templates, '
                        'modifications or affinity prediction.')
    g.add_argument('--yaml-target-chain', default='A', metavar='ID',
                   help='Chain in --yaml-template to substitute the homolog into.')
    g.add_argument('--no-validate', action='store_true',
                   help='Skip pre-flight validation of SMILES/CCD/sequences.')

    m = parser.add_argument_group(
        'MSA reuse',
        'Hand Boltz a precomputed MSA instead of letting it call the MSA server. '
        'Validated before any stage runs: Boltz itself only checks that the file '
        'exists, and consumes MSA columns positionally, so a mismatched alignment '
        'would give a confident but meaningless prediction.')
    m.add_argument('--use-custom-msa', '--use_custom_msa', dest='use_custom_msa',
                   metavar='PATH',
                   help='Either ONE alignment containing the homologs being run -- each '
                        'homolog\'s own row becomes its query and the rest become its '
                        'MSA -- or a DIRECTORY with one alignment per homolog named '
                        '<homolog_id>.<ext>. Accepts a3m, Boltz csv, aligned FASTA, '
                        'Stockholm or Clustal; all are converted to Boltz csv.')
    m.add_argument('--reuse-msa-from', dest='reuse_msa_from', metavar='RUN_DIR',
                   help='A previous run directory to take already-computed MSAs from '
                        '(read-only). Matched per homolog against Boltz\'s own '
                        'boltz_results_<id>/msa/ layout.')
    return parser


# -- YAML emission ------------------------------------------------------

def _yaml():
    try:
        import yaml
        return yaml
    except ImportError:
        sys.exit("ERROR: PyYAML is required for --format yaml (it ships with boltz).")


def load_yaml_template(path):
    """Load a Boltz YAML template, or synthesise a minimal one if none given."""
    if not path:
        return {'version': 1, 'sequences': [{'protein': {'id': 'A', 'sequence': ''}}]}
    if not os.path.isfile(path):
        sys.exit(f"ERROR: --yaml-template '{path}' not found.")
    doc = _yaml().safe_load(open(path))
    if not isinstance(doc, dict) or 'sequences' not in doc:
        sys.exit(f"ERROR: '{path}' is not a Boltz YAML (no top-level 'sequences' key).")
    return doc


def render_yaml_for_homolog(template, sequence, target_chain, entities, msa_path=None):
    """
    Deep-copy the template, substitute `sequence` into the protein entity
    whose id is `target_chain`, and append any extra entities on fresh
    chain IDs that do not collide with the template's own.
    """
    import copy
    doc = copy.deepcopy(template)
    seqs = doc.setdefault('sequences', [])

    used = set()
    for item in seqs:
        for body in item.values():
            ids = body.get('id')
            used.update(ids if isinstance(ids, list) else [ids])

    replaced = False
    for item in seqs:
        body = item.get('protein')
        if not body:
            continue
        ids = body.get('id')
        ids = ids if isinstance(ids, list) else [ids]
        if target_chain in ids:
            body['sequence'] = sequence
            if msa_path:
                body['msa'] = str(msa_path)
            replaced = True
            break
    if not replaced:
        sys.exit(f"ERROR: --yaml-template has no protein chain with id "
                 f"'{target_chain}' to substitute the homolog into.")

    if entities:
        ids = (c for c in chain_ids() if c not in used)
        for ent in entities:
            for _ in range(max(1, int(ent.get('copies', 1)))):
                cid = next(ids)
                if ent['type'] in ('ccd', 'smiles'):
                    seqs.append({'ligand': {'id': cid, ent['type']: ent['value']}})
                else:
                    body = {'id': cid, 'sequence': ent['value']}
                    if ent['type'] == 'protein' and ent.get('msa'):
                        body['msa'] = ent['msa']
                    seqs.append({ent['type']: body})
    return doc


def write_yaml(doc, path):
    with open(path, 'w') as fh:
        _yaml().safe_dump(doc, fh, sort_keys=False, default_flow_style=False, width=10 ** 6)



# ======================================================================
# MSA reuse
# ======================================================================
#
# Boltz accepts a precomputed MSA in exactly two formats -- ".a3m" and
# ".csv" (boltz/main.py: "MSA file ... not supported, only a3m or csv") --
# and reads it from the 3rd '|'-delimited field of a protein's FASTA header,
# or a protein entity's `msa:` key in YAML. Its ".csv" flavour is two
# columns, `key,sequence`, whose sequence column carries a3m insertion
# semantics (lowercase = insertion relative to the query, recorded as a
# deletion; boltz/data/parse/csv.py).
#
# Everything else an alignment tool emits -- aligned FASTA from MAFFT,
# Stockholm from hmmer/jackhmmer, Clustal -- has to be converted, which is
# what this section does. It always emits Boltz CSV.
#
# Two things Boltz itself does NOT check, and this section therefore must:
#
#   1. That the MSA's query row actually corresponds to the chain it is
#      attached to. Boltz validates only that the file exists and has a
#      known suffix; MSA columns are then consumed positionally. A
#      mismatched MSA yields a confident, wrong prediction with no warning.
#   2. That a custom MSA is not mixed with an auto-generated one in the
#      same input, which Boltz rejects late with "Cannot mix custom and
#      auto-generated MSAs in the same input!" (data/parse/schema.py).
#
# Both are enforced before any stage runs -- see preflight_msa().

MSA_SUFFIXES = ('.a3m', '.csv', '.sto', '.stk', '.stockholm', '.clustal',
                '.aln', '.afa', '.fasta', '.fa', '.fas', '.mfa', '.faa')

# Below this many sequences a custom MSA is worth flagging. The MSA server
# typically returns thousands; an alignment of just the homolog set being run
# may hold a handful, and Boltz confidence drops steeply with depth. Measured:
# the same two receptors scored ~0.82-0.91 on server MSAs and ~0.40-0.61 on a
# 6-sequence alignment. Not an error -- shallow is sometimes exactly what is
# wanted -- but never silent.
MSA_SHALLOW_WARN = 32

GAP_CHARS = '-.'


def _strip_gaps(seq):
    """The ungapped, upper-case residues of an aligned row."""
    return ''.join(c for c in seq if c not in GAP_CHARS).upper()


def sniff_msa_format(path):
    """
    Identify an alignment file's format from its first non-blank lines,
    falling back to the suffix. Returns one of
    'csv' | 'a3m' | 'aligned-fasta' | 'stockholm' | 'clustal'.
    """
    suffix = Path(path).suffix.lower()
    try:
        with open(path, errors='replace') as fh:
            head = [ln.rstrip('\n') for _, ln in zip(range(8), fh)]
    except OSError as e:
        sys.exit(f"ERROR: cannot read MSA '{path}': {e}")

    nonblank = [ln for ln in head if ln.strip()]
    if not nonblank:
        sys.exit(f"ERROR: MSA '{path}' is empty.")
    first = nonblank[0]

    if first.startswith('# STOCKHOLM'):
        return 'stockholm'
    if first.upper().startswith('CLUSTAL'):
        return 'clustal'
    if first.startswith('>'):
        # a3m and aligned FASTA are both '>'-headed. a3m is distinguished by
        # lowercase insertion characters, which aligned FASTA never carries.
        body = ''.join(ln for ln in nonblank if not ln.startswith('>'))
        return 'a3m' if any(c.islower() for c in body if c.isalpha()) else 'aligned-fasta'
    if suffix == '.csv' or ',' in first:
        return 'csv'
    if suffix in ('.sto', '.stk', '.stockholm'):
        return 'stockholm'
    if suffix in ('.clustal',):
        return 'clustal'
    sys.exit(f"ERROR: cannot identify the format of MSA '{path}'.\n"
             f"       Supported: a3m, Boltz csv, aligned FASTA, Stockholm, Clustal.")


def read_msa_rows(path):
    """
    Read any supported alignment into [(label, aligned_sequence)], in file
    order. a3m rows are kept verbatim (lowercase insertions preserved,
    which is what Boltz's parser expects); every other format is a
    fixed-width alignment.
    """
    fmt = sniff_msa_format(path)

    if fmt == 'csv':
        rows = []
        with open(path, newline='') as fh:
            reader = csv.DictReader(fh)
            cols = tuple(sorted(reader.fieldnames or []))
            if cols != ('key', 'sequence'):
                sys.exit(f"ERROR: '{path}' is not a Boltz MSA csv: columns are "
                         f"{list(reader.fieldnames or [])}, expected ['key', 'sequence'].")
            for r in reader:
                rows.append((str(r.get('key') or ''), (r.get('sequence') or '').strip()))
        return rows, fmt

    if fmt in ('a3m', 'aligned-fasta'):
        rows, label, buf = [], None, []
        with open(path, errors='replace') as fh:
            for line in fh:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if label is not None:
                        rows.append((label, ''.join(buf)))
                    label, buf = line[1:].strip(), []
                elif line.strip():
                    buf.append(line.strip())
        if label is not None:
            rows.append((label, ''.join(buf)))
        return rows, fmt

    if fmt == 'stockholm':
        seqs = {}
        order = []
        with open(path, errors='replace') as fh:
            for line in fh:
                line = line.rstrip('\n')
                if not line.strip() or line.startswith('#') or line.startswith('//'):
                    continue
                parts = line.split(None, 1)
                if len(parts) != 2:
                    continue
                name, chunk = parts[0], parts[1].replace(' ', '')
                if name not in seqs:
                    seqs[name] = []
                    order.append(name)
                seqs[name].append(chunk)
        return [(n, ''.join(seqs[n])) for n in order], fmt

    # clustal
    seqs, order = {}, []
    with open(path, errors='replace') as fh:
        for i, line in enumerate(fh):
            if i == 0 or not line.strip():
                continue
            if line[0].isspace():          # conservation line
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            name, chunk = parts[0], parts[1]
            if name not in seqs:
                seqs[name] = []
                order.append(name)
            seqs[name].append(chunk)
    return [(n, ''.join(seqs[n])) for n in order], fmt


def project_alignment_to_query(rows, query_idx):
    """
    Re-centre a multiple alignment on one of its own rows, producing the
    query-centric (a3m) form Boltz requires.

    Boltz indexes an MSA row by the query chain's residue index
    (`msa_residues[res_start + res_idx]` in data/feature/featurizer.py), so
    the query row must be ungapped and exactly as long as the chain. Simply
    reordering a MAFFT-style alignment is therefore not enough: every column
    in which the query has a gap is an insertion relative to the query and
    must be removed from the columns, re-emitted as a3m lowercase on the
    rows that have a residue there.

    This is what lets ONE alignment of the whole homolog set serve every
    homolog in the run: each homolog is projected onto its own columns in
    turn.
    """
    q = rows[query_idx][1]
    out = []
    for label, seq in rows:
        chars, pending = [], []
        for col, qc in enumerate(q):
            c = seq[col] if col < len(seq) else '-'
            if qc in GAP_CHARS:
                # Column the query does not occupy: an insertion in `seq`.
                if c not in GAP_CHARS:
                    pending.append(c.lower())
            else:
                if pending:
                    chars.append(''.join(pending))
                    pending = []
                chars.append('-' if c in GAP_CHARS else c.upper())
        out.append((label, ''.join(chars)))
    return out


def locate_query_row(rows, query_seq, query_id):
    """
    Find the alignment row belonging to this homolog: by exact ungapped
    sequence first (authoritative), then by header id. Returns an index, or
    None.
    """
    target = _strip_gaps(query_seq)
    for i, (_lab, seq) in enumerate(rows):
        if _strip_gaps(seq) == target:
            return i
    for i, (lab, _seq) in enumerate(rows):
        parts = lab.split()
        if parts and parts[0] == query_id:
            return i
    return None


_ALIGNMENT_CACHE = {}


def load_alignment(path):
    """Read and cache an alignment, so one shared file is parsed once."""
    key = str(Path(path).resolve())
    if key not in _ALIGNMENT_CACHE:
        _ALIGNMENT_CACHE[key] = read_msa_rows(path)
    return _ALIGNMENT_CACHE[key]


def convert_msa_to_boltz_csv(src, dest, query_seq, query_id):
    """
    Write this homolog's Boltz MSA csv, derived from `src`.

    `src` may be an alignment of just this homolog, or -- the common case --
    one alignment of the whole homolog set, from which this homolog's row is
    taken as the query and the rest become its MSA.

    Returns (dest_path, n_rows, query_row_index). Raises MSAError when the
    homolog is not present in the alignment or the alignment is ragged.
    """
    rows, fmt = load_alignment(src)
    if not rows:
        raise MSAError(f"'{src}' contains no sequences")

    target = _strip_gaps(query_seq)
    qi = locate_query_row(rows, query_seq, query_id)
    if qi is None:
        first_lab, first_seq = rows[0]
        raise MSAError(
            f"not found in '{Path(src).name}' ({len(rows)} sequences).\n"
            f"         Looked for a row whose ungapped sequence is this homolog's "
            f"({len(target)} aa)\n"
            f"         or whose header is '{query_id}'.\n"
            f"         query : {_abbrev(target, 50)}\n"
            f"         row 1 : '{_abbrev(first_lab, 24)}' "
            f"({len(_strip_gaps(first_seq))} aa) {_abbrev(_strip_gaps(first_seq), 40)}\n"
            f"         Boltz aligns MSA columns to the query by position and does not "
            f"check this itself,\n"
            f"         so the wrong row here gives a confident but meaningless "
            f"prediction.")

    # Ragged rows are legitimate in the query-centric formats (a3m, and
    # Boltz csv which uses the same insertion semantics), where lowercase
    # insertions genuinely vary the width. Anything else must be a real
    # fixed-width alignment, and is re-centred on the query below.
    if fmt not in ('a3m', 'csv'):
        widths = sorted({len(seq) for _lab, seq in rows})
        if len(widths) > 1:
            raise MSAError(
                f"'{Path(src).name}' is not a valid alignment: rows have "
                f"{len(widths)} different lengths ({', '.join(map(str, widths[:4]))}"
                f"{' ...' if len(widths) > 4 else ''}). Every column must line up.")
        projected = project_alignment_to_query(rows, qi)
    else:
        # Already query-centric; only the row order needs fixing.
        projected = rows

    ordered = [projected[qi]] + [r for i, r in enumerate(projected) if i != qi]

    # The query row must now be exactly the chain sequence, or Boltz will
    # read every row at the wrong offset.
    got = _strip_gaps(ordered[0][1])
    if got != target:
        raise MSAError(
            f"internal check failed for '{query_id}': after re-centring, the query "
            f"row is {len(got)} aa but the homolog is {len(target)} aa. "
            f"The alignment may be inconsistent.")

    Path(dest).parent.mkdir(parents=True, exist_ok=True)
    with open(dest, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['key', 'sequence'])
        for _lab, seq in ordered:
            # `key` is Boltz's taxonomy/pairing field; blank means unpaired,
            # which is what a single-chain MSA wants.
            w.writerow(['', seq])
    return str(Path(dest).resolve()), len(ordered), qi


class MSAError(Exception):
    """A custom MSA that Boltz would accept but that is wrong for this run."""


def find_msa_in_run(run_dir, stem):
    """
    Locate a homolog's MSA inside a previous run, in Boltz's own output
    layout (boltz_results_<stem>/msa/<stem>_<n>.csv). Read-only: a prior run
    is never written to. Also accepts this script's own converted-MSA
    folders. Returns an absolute path, or None.
    """
    run = Path(run_dir)
    patterns = [
        f'**/boltz_results_{stem}/msa/{stem}_*.csv',   # Boltz's own output
        f'**/msa/{stem}.csv',                          # our converted cache
        f'**/{stem}.a3m',
        f'**/{stem}.csv',
    ]
    for pat in patterns:
        hits = sorted(run.glob(pat))
        if hits:
            return str(hits[0].resolve())
    return None


def resolve_msa_sources(records, custom, reuse_from):
    """
    Map homolog id -> source alignment path.

    --use-custom-msa accepts either:
      * ONE alignment containing the homologs being run. Each homolog's own
        row becomes the query and the rest become its MSA (see
        convert_msa_to_boltz_csv / project_alignment_to_query). This is the
        normal case: the homolog set usually already has an alignment.
      * a DIRECTORY holding one alignment per homolog, named <homolog_id>.<ext>,
        for when each homolog has its own separately-built MSA.

    --reuse-msa-from points at a previous run directory and picks up whatever
    MSAs it already computed, per homolog.

    Returns (mapping, problems).
    """
    mapping, problems = {}, []
    ids = [pid for pid, _seq in records]

    if custom:
        p = Path(custom).expanduser()
        if not p.exists():
            problems.append(f"--use-custom-msa path '{custom}' does not exist")
        elif p.is_dir():
            for pid in ids:
                hits = [q for suf in MSA_SUFFIXES for q in sorted(p.glob(f'{pid}{suf}'))]
                if hits:
                    mapping[pid] = str(hits[0].resolve())
                else:
                    problems.append(
                        f"no MSA for homolog '{pid}' in {p}/ "
                        f"(expected {pid} with one of: "
                        f"{', '.join(MSA_SUFFIXES[:5])} ...)")
        else:
            # One shared alignment: every homolog is projected out of it.
            for pid in ids:
                mapping[pid] = str(p.resolve())

    if reuse_from:
        r = Path(reuse_from).expanduser()
        if not r.is_dir():
            problems.append(f"--reuse-msa-from '{reuse_from}' is not a directory")
        else:
            for pid in ids:
                if pid in mapping:
                    continue          # an explicit --use-custom-msa wins
                found = find_msa_in_run(r, pid)
                if found:
                    mapping[pid] = found
    return mapping, problems


def input_has_custom_msa(path):
    """
    True if a Boltz input file already points its protein chain at a
    precomputed MSA. FASTA carries it as the 3rd '|' field; YAML as an
    entity's `msa:` key. The literal "empty" is not a custom MSA -- it means
    single-sequence mode.
    """
    text = Path(path).read_text()
    if Path(path).suffix.lower() in ('.yaml', '.yml'):
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith('msa:'):
                val = stripped.split(':', 1)[1].strip()
                if val and val not in ('empty', '0', 'null', '~'):
                    return True
        return False
    for line in text.splitlines():
        if line.startswith('>'):
            parts = line[1:].split('|')
            if len(parts) >= 3 and parts[1].strip().lower() == 'protein':
                if parts[2].strip() and parts[2].strip() != 'empty':
                    return True
    return False


def with_msa_path(content, msa_path, is_yaml):
    """
    Point the first protein chain of a Boltz input at `msa_path`, leaving
    every other chain untouched. Used to hand seeds 2..N the MSA that seed 1
    generated, so a sweep makes one MSA-server call per homolog instead of
    one per seed.
    """
    if is_yaml:
        out, done = [], False
        for line in content.splitlines():
            out.append(line)
            if not done and line.strip().startswith('sequence:'):
                indent = ' ' * (len(line) - len(line.lstrip()))
                out.append(f"{indent}msa: {msa_path}")
                done = True
        return '\n'.join(out) + '\n'

    out, done = [], False
    for line in content.splitlines():
        if not done and line.startswith('>'):
            parts = line[1:].split('|')
            if len(parts) >= 2 and parts[1].strip().lower() == 'protein':
                out.append(f">{parts[0]}|{parts[1]}|{msa_path}")
                done = True
                continue
        out.append(line)
    return '\n'.join(out) + '\n'


def read_boltz_confidence(pred_dir, stem, idx):
    """
    Read Boltz's own confidence JSON for one predicted model.

    These are the scalars Boltz already wrote (confidence_score, ptm, iptm,
    complex_plddt, ...) plus its per-chain dicts. Nothing is recomputed here
    -- deriving per-chain pLDDT/PAE from the token-level .npz arrays is a
    separate analysis and deliberately not part of this script.
    """
    path = Path(pred_dir) / f'confidence_{stem}_model_{idx}.json'
    if not path.is_file():
        return {}
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def build_sweep_manifest(output_folder, ensemble_folder, stems, seeds):
    """
    Index every model a sweep produced, as one JSON object per line.

    Rebuilt from what is actually on disk at the end of each sweep, so a
    resumed run still yields a complete manifest rather than only the models
    written in that invocation.

    One row per model: which homolog, which seed, which diffusion sample,
    where the file is, and Boltz's confidence metrics for it.
    """
    rows = []
    for stem in stems:
        for seed in seeds:
            pred = (Path(output_folder) / stem / f'seed{seed}' /
                    f'boltz_results_{stem}' / 'predictions' / stem)
            if not pred.is_dir():
                continue
            for cif in sorted(pred.glob(f'{stem}_model_*.cif')):
                m = re.search(r'_model_(\d+)\.cif$', cif.name)
                if not m:
                    continue
                idx = int(m.group(1))
                row = {
                    'homolog': stem,
                    'seed': seed,
                    'model': idx,
                    'ensemble_cif': str((Path(ensemble_folder) / stem /
                                         f'{stem}_seed{seed}_model_{idx}.cif')),
                    'source_cif': str(cif),
                }
                row.update(read_boltz_confidence(pred, stem, idx))
                rows.append(row)

    manifest = Path(output_folder) / 'manifest.jsonl'
    with open(manifest, 'w') as fh:
        for row in rows:
            fh.write(json.dumps(row) + '\n')
    return manifest, len(rows)


def seed_is_complete(seed_out_dir, stem, n_samples):
    """
    True only when every diffusion sample of this seed has both its structure
    and its confidence file. Checking for model_0 alone would treat a seed
    interrupted partway through its samples as finished, and --resume would
    then leave that seed permanently short. Mirrors seed_is_complete() in
    bin/boltz-2_multi_seed_run/run_boltz_sweep.py.
    """
    pred = Path(seed_out_dir) / f'boltz_results_{stem}' / 'predictions' / stem
    if not pred.is_dir():
        return False
    for idx in range(n_samples):
        if not (pred / f'{stem}_model_{idx}.cif').exists():
            return False
        if not (pred / f'confidence_{stem}_model_{idx}.json').exists():
            return False
    return True


def find_generated_msa(seed_out_dir, stem):
    """
    Locate the MSA Boltz wrote during a --use_msa_server run, so the
    remaining seeds can reuse it. Boltz's layout is
    boltz_results_<stem>/msa/<stem>_<n>.csv.
    """
    msa_dir = Path(seed_out_dir) / f'boltz_results_{stem}' / 'msa'
    hits = sorted(msa_dir.glob(f'{stem}_*.csv')) if msa_dir.is_dir() else []
    return str(hits[0].resolve()) if hits else None


def preflight_msa(records, args, converted_dir, log=print):
    """
    Resolve, convert and validate every custom MSA BEFORE any stage runs,
    so a mismatched or unreadable alignment costs seconds instead of
    surfacing as a silently-wrong prediction after hours of GPU.

    Returns {homolog_id: absolute path to a Boltz csv}. Exits on any problem,
    reporting all of them at once.
    """
    custom = getattr(args, 'use_custom_msa', None)
    reuse = getattr(args, 'reuse_msa_from', None)
    if not custom and not reuse:
        return {}

    # Boltz rejects an input that mixes a custom MSA with an auto-generated
    # one, so catch the flag combination that would cause it.
    if getattr(args, 'peptide_msa', False):
        sys.exit("ERROR: --peptide-msa cannot be combined with a custom MSA.\n"
                 "       Boltz refuses to mix custom and auto-generated MSAs in one "
                 "input\n       (\"Cannot mix custom and auto-generated MSAs in the "
                 "same input!\").\n"
                 "       Drop --peptide-msa: peptide ligand chains default to an "
                 "explicit empty MSA,\n       which is both valid here and what a "
                 "short peptide wants.")

    sources, problems = resolve_msa_sources(records, custom, reuse)

    if custom:
        # --reuse-msa-from is best-effort (a partial prior run is normal and
        # the rest can be fetched); --use-custom-msa is a promise, so every
        # homolog must be covered.
        missing = [pid for pid, _ in records if pid not in sources]
        if missing and not problems:
            problems.append(f"no MSA found for {len(missing)} homolog(s): "
                            f"{', '.join(missing[:5])}"
                            f"{' ...' if len(missing) > 5 else ''}")

    resolved, seq_by_id, depths = {}, dict(records), []
    for pid, src in sorted(sources.items()):
        dest = Path(converted_dir) / f'{pid}.csv'
        try:
            path, n_rows, qrow = convert_msa_to_boltz_csv(src, dest, seq_by_id[pid], pid)
            resolved[pid] = path
            note = f", query taken from row {qrow + 1}" if qrow else ""
            log(f"MSA {pid}: {n_rows} sequences from {Path(src).name}{note}")
            depths.append((pid, n_rows))
        except MSAError as e:
            problems.append(f"{pid}: {e}")

    if problems:
        sys.stderr.write("ERROR: custom MSA validation failed:\n")
        for p in problems:
            sys.stderr.write(f"  - {p}\n")
        sys.exit(1)

    shallow = [(pid, n) for pid, n in depths if n < MSA_SHALLOW_WARN]
    if shallow:
        worst = min(n for _pid, n in shallow)
        log(f"WARNING: {len(shallow)}/{len(depths)} custom MSAs are shallow "
            f"(fewest: {worst} sequences, warning below {MSA_SHALLOW_WARN}). "
            f"The MSA server typically returns thousands; Boltz confidence drops "
            f"steeply with depth, so these predictions will be markedly worse than "
            f"server-MSA ones. Use a deeper alignment if that is not intended.")

    return resolved


# ======================================================================
# The HomoLogic pipeline
# ======================================================================

class HomoLogic:
    def __init__(self, reference_fasta="reference.fasta", reference_pdb="reference.pdb",
                 reference_ligands_pdb="reference_ligands.pdb", homologs_fasta="homologs.fasta",
                 boltz_bin=None, mmseqs_bin=None, tmalign_bin=None,
                 log_file=None, overwrite=False, resume=False):
        # These four paths default to plain filenames in the current working
        # directory. They are NOT auto-discovered anywhere in this class, so
        # in practice they must be explicitly supplied (via the CLI flags
        # --reference-fasta / --reference-pdb / --reference-ligands-pdb /
        # --homologs-fasta added below) unless the caller happens to be
        # running from a directory that already contains files with these
        # exact names.
        self.reference_fasta = reference_fasta            # reference protein sequence (FASTA)
        self.reference_pdb = reference_pdb                # reference protein structure (PDB)
        self.reference_ligands_pdb = reference_ligands_pdb  # reference-bound ligand/cofactor coords (PDB)
        self.homologs_fasta = homologs_fasta               # multi-sequence FASTA of candidate homologs

        # External binaries, resolved once to absolute paths (see
        # resolve_exe: this interpreter's own env wins over $PATH) so a run
        # can never silently pick up a different Boltz/MMseqs install than
        # the one this script was installed alongside. Resolved lazily --
        # a stage that does not shell out to a tool must not require it.
        self._bin_override = {'boltz': boltz_bin, 'mmseqs': mmseqs_bin, 'TMalign': tmalign_bin}
        self._bin_cache = {}

        # Re-run safety. safe_create_output_folder() deletes directories, so
        # it refuses to touch an existing one unless overwrite=True.
        self.overwrite = overwrite
        self.resume = resume

        # Run log. Console output is unchanged; every line is additionally
        # appended to log_file as "timestamp LEVEL message".
        self.log_file = log_file
        if self.log_file:
            os.makedirs(os.path.dirname(os.path.abspath(self.log_file)) or '.', exist_ok=True)
            with open(self.log_file, 'a') as fh:
                fh.write(f"{datetime.now():%Y-%m-%d %H:%M:%S} INFO  "
                         f"=== homologic2 run started: {' '.join(sys.argv)} ===\n")

    def tool(self, name):
        """Absolute path to an external binary, resolved and cached on first use."""
        if name not in self._bin_cache:
            path = resolve_exe(name, self._bin_override.get(name))
            self._bin_cache[name] = path
            self.print_progress(f"Using {name}: {path} ({exe_version(path)})")
        return self._bin_cache[name]

    # -- Lightweight timestamped logging. Console formatting is byte-for-byte
    #    what homologic_gs.py printed; the file sink is the addition. -------
    def _log(self, level, message):
        if not self.log_file:
            return
        try:
            with open(self.log_file, 'a') as fh:
                fh.write(f"{datetime.now():%Y-%m-%d %H:%M:%S} {level:<5} {message}\n")
        except OSError:
            pass  # never let logging kill a running pipeline

    def print_progress(self, message, start_time=None):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if start_time:
            elapsed = (datetime.now() - start_time).total_seconds()
            print(f"[{timestamp}] {message} ({elapsed:.1f}s)")
        else:
            print(f"[{timestamp}] {message}")
        self._log('INFO', message)

    def print_success(self, message, start_time=None):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if start_time:
            elapsed = (datetime.now() - start_time).total_seconds()
            print(f"[{timestamp}] SUCCESS: {message} ({elapsed:.1f}s)")
        else:
            print(f"[{timestamp}] SUCCESS: {message}")
        self._log('OK', message)

    def print_error(self, message):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] ERROR: {message}", file=sys.stderr)
        self._log('ERROR', message)

    def print_warning(self, message):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] WARNING: {message}")
        self._log('WARN', message)

    def safe_create_output_folder(self, folder_path):
        """
        Prepare a stage's output folder.

        homologic_gs.py deleted any existing folder of the same name without
        asking, which meant re-running `boltz-model` silently threw away
        hours of finished GPU work. Deleting now requires an explicit
        --overwrite; with --resume the existing folder is kept and reused.
        """
        if os.path.exists(folder_path):
            if self.resume and os.path.isdir(folder_path):
                self.print_progress(f"Output folder reused (--resume): {folder_path}")
                return
            if not self.overwrite:
                raise SystemExit(
                    f"ERROR: '{folder_path}' already exists and would be deleted.\n"
                    f"       Re-run with --overwrite to replace it, --resume to continue\n"
                    f"       into it, or point the stage at a different folder.")
            if os.path.isdir(folder_path):
                self.print_warning(f"Output folder '{folder_path}' exists → overwriting")
                shutil.rmtree(folder_path)
            else:
                self.print_warning(f"File '{folder_path}' exists → removing")
                os.remove(folder_path)
        os.makedirs(folder_path, exist_ok=True)
        self.print_progress(f"Output folder prepared: {folder_path}")
    @staticmethod
    def compute_identity_gaps_mismatches_openings(ref_aln_str, tgt_aln_str):
        # Given two already-aligned strings (same length, '-' = gap), walk
        # them column-by-column to derive percent identity, mismatch count,
        # and the number of distinct gap "openings" (a run of consecutive
        # gap columns on either strand counts once, not once per column).
        matches = 0
        aln_len = 0
        for a, b in zip(ref_aln_str, tgt_aln_str):
            if a != '-' and b != '-':
                aln_len += 1
                if a == b:
                    matches += 1
        identity = matches / aln_len if aln_len > 0 else 0.0
        mismatches = aln_len - matches
        opens = 0
        in_gap = False
        for a, b in zip(ref_aln_str, tgt_aln_str):
            is_gap = (a == '-' or b == '-')
            if is_gap and not in_gap:
                opens += 1
                in_gap = True
            elif not is_gap:
                in_gap = False
        return identity, mismatches, opens, aln_len

    def align_and_compute_identity_stats(self, ref_seq, tgt_seq, matrix_name='BLOSUM62',
                                         gap_open=11, gap_extend=1, return_alignment=False):
        # scikit-bio global (Needleman-Wunsch) protein alignment; used both
        # as the MMseqs2 fallback in calculate_sequence_homology() and for
        # binding-site metasequence comparison in analyze_binding_sites().
        alignment, _, _ = global_pairwise_align_protein(
            ref_seq, tgt_seq,
            gap_open_penalty=gap_open,
            gap_extend_penalty=gap_extend,
            substitution_matrix=substitution_matrices.load(matrix_name),
            penalize_terminal_gaps=False
        )
        if len(alignment) != 2:
            return None
        ref_aln = str(alignment[0])
        tgt_aln = str(alignment[1])
        identity, mismatches, gap_openings, aln_len = self.compute_identity_gaps_mismatches_openings(ref_aln, tgt_aln)
        result = {
            'identity': identity,
            'identity_pct_str': f"{identity:.3f}",
            'mismatches': mismatches,
            'gap_openings': gap_openings,
            'alignment_length': aln_len,
            'full_alignment_length': len(ref_aln)
        }
        if return_alignment:
            result['ref_aln'] = ref_aln
            result['tgt_aln'] = tgt_aln
        return result

    def calculate_sequence_homology(
        self,
        output_folder='01_sequence_homology_results',
        substitution_matrix_name='BLOSUM62',
        seqid_beyond_mmseqs=True
    ):
        """
        Stage 1: pairwise sequence identity of every homolog against the
        reference, primarily via MMseqs2 (self.reference_fasta vs.
        self.homologs_fasta -- both REQUIRED to already exist on disk),
        with a scikit-bio global-alignment fallback for any homolog MMseqs2
        does not report a usable percent-identity for. Requires the
        `mmseqs` binary on PATH.
        """
        start_time = datetime.now()
        self.print_progress("Starting sequence homology calculation", start_time)

        self.safe_create_output_folder(output_folder)
        tmp_dir = os.path.join(output_folder, 'tmp_mmseqs')
        self.safe_create_output_folder(tmp_dir)

        result_m8 = os.path.join(output_folder, 'result.m8')
        output_csv = os.path.join(output_folder, 'sequence_homology.csv')

        # Run MMseqs2 easy-search: reference (query) vs. homologs (target).
        cmd = [self.tool('mmseqs'), 'easy-search', '--max-seqs', '100000',
               self.reference_fasta, self.homologs_fasta, result_m8, tmp_dir]
        subprocess.run(cmd, check=True)

        homolog_seqs = {
            seq.metadata['id'].strip(): seq
            for seq in io.read(self.homologs_fasta, format='fasta', constructor=Protein)
        }
        self.print_progress(f"Loaded {len(homolog_seqs)} sequences")

        ref_seq_record = next(io.read(self.reference_fasta, format='fasta', constructor=Protein))
        ref_id = ref_seq_record.metadata['id'].strip()
        ref_seq = ref_seq_record

        headers = [
            'reference_id', 'homolog_id', 'sequence_identity', 'alignment_length',
            'mismatches', 'gap_openings', 'query_start', 'query_end',
            'target_start', 'target_end', 'e_value', 'bit_score'
        ]
        rows = []
        good_hits = set()       # MMseqs2 reported a sane 0-1 percent identity for these
        fallback_targets = set()  # MMseqs2 hit exists but pident looked malformed/out-of-range

        with open(result_m8, 'r') as f:
            for line in f:
                row = line.strip().split('\t')
                if len(row) != 12:
                    continue
                target_id = row[1].strip()
                try:
                    pident = float(row[2])
                    if 0 < pident <= 1:
                        good_hits.add(target_id)
                        rows.append(row)
                        continue
                except ValueError:
                    pass

                fallback_targets.add(target_id)
                rows.append(row)

        all_targets = set(homolog_seqs)
        unmatched = all_targets - good_hits - fallback_targets  # no MMseqs2 hit at all

        self.print_progress(
            f"MMseqs2: {len(good_hits)} good | {len(fallback_targets)} fallback | {len(unmatched)} unmatched"
        )

        # Re-derive identity/gap stats via global alignment for any hit whose
        # MMseqs2 pident column looked malformed, keeping MMseqs2's e-value
        # and bit-score columns (indices 10, 11) since those are still valid.
        if fallback_targets:
            self.print_progress(f"Global fallback for {len(fallback_targets)} sequences")
            for target_id in sorted(fallback_targets):
                if target_id not in homolog_seqs:
                    continue
                tgt_seq = homolog_seqs[target_id]
                stats = self.align_and_compute_identity_stats(ref_seq, tgt_seq, substitution_matrix_name)
                if stats:
                    for i, r in enumerate(rows):
                        if r[1].strip() == target_id:
                            rows[i] = [
                                ref_id, target_id,
                                stats['identity_pct_str'],
                                str(stats['full_alignment_length']),
                                str(stats['mismatches']),
                                str(stats['gap_openings']),
                                '1', str(len(ref_seq)),
                                '1', str(len(tgt_seq)),
                                r[10], r[11]
                            ]
                            break

        # Sequences MMseqs2 didn't report at all: optionally still compute
        # sequence identity via slower global alignment (seqid_beyond_mmseqs
        # controls this trade-off between completeness and runtime).
        if unmatched and seqid_beyond_mmseqs:
            self.print_progress(f"Global alignment for {len(unmatched)} unmatched sequences")
            for target_id in sorted(unmatched):
                if target_id not in homolog_seqs:
                    rows.append([ref_id, target_id] + ['n.d.'] * 10)
                    continue

                tgt_seq = homolog_seqs[target_id]
                stats = self.align_and_compute_identity_stats(ref_seq, tgt_seq, substitution_matrix_name)

                if stats:
                    rows.append([
                        ref_id, target_id,
                        stats['identity_pct_str'],
                        str(stats['full_alignment_length']),
                        str(stats['mismatches']),
                        str(stats['gap_openings']),
                        '1', str(len(ref_seq)),
                        '1', str(len(tgt_seq)),
                        'n.d.', 'n.d.'
                    ])
                else:
                    rows.append([ref_id, target_id] + ['n.d.'] * 10)

        elif unmatched and not seqid_beyond_mmseqs:
            self.print_progress(
                f"Skipping global alignment for {len(unmatched)} unmatched sequences "
                f"(seqid_beyond_mmseqs=False)"
            )
            for target_id in sorted(unmatched):
                rows.append([ref_id, target_id] + ['n.d.'] * 10)

        with open(output_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(rows)

        self.print_success(f"Results saved → {output_csv} ({len(rows)} entries)", start_time)

        try:
            os.remove(result_m8)
            shutil.rmtree(tmp_dir)
        except OSError:
            pass

        return output_csv

    def generate_boltz_input(self, entities=None, output_folder='02_boltz_input',
                             fmt='fasta', yaml_template=None, yaml_target_chain='A',
                             validate=True, smiles_code=None, msa_map=None):
        """
        Stage 2: split self.homologs_fasta (REQUIRED to exist) into one Boltz
        input file per homolog, with the homolog as chain A and every
        co-folded entity appended as its own chain.

        `entities` is an ordered list of dicts as produced by
        build_entities(): {'type': 'protein'|'dna'|'rna'|'ccd'|'smiles',
        'value': str, 'msa': str|None, 'copies': int}. Every entity type
        Boltz accepts is supported -- homologic_gs.py could only emit up to
        four SMILES chains, which is why peptide ligands previously had to be
        hand-encoded as one enormous SMILES string and metal ions as things
        like "[Ca+2]" instead of the CCD code.

        fmt='fasta' (default) writes Boltz's FASTA format:
            >CHAIN_ID|ENTITY_TYPE[|MSA_ID]
        fmt='yaml' clones `yaml_template` per homolog, substituting the
        homolog sequence into chain `yaml_target_chain`. YAML is the only
        format that can carry constraints/templates/affinity/modifications,
        so the template is the escape hatch for anything the flags above do
        not model.

        `msa_map` maps homolog id -> a precomputed Boltz MSA csv, produced
        and validated by preflight_msa(). When a homolog has one, its chain A
        carries that path (FASTA header field 3, or the YAML `msa:` key) and
        Boltz reuses it instead of calling the MSA server. Ligand chains keep
        their explicit empty MSA, which is what makes the combination legal:
        Boltz refuses to mix a custom MSA with an auto-generated one.

        `smiles_code` is accepted for backwards compatibility with
        homologic.py-era driver scripts and is folded into `entities`.
        """
        start_time = datetime.now()
        self.print_progress(f"Generating Boltz input: {self.homologs_fasta}", start_time)

        if not os.path.exists(self.homologs_fasta):
            self.print_error(f"Missing: {self.homologs_fasta}")
            return

        self.safe_create_output_folder(output_folder)

        records = []
        current_id = None
        current_seq = []
        # Minimal hand-rolled FASTA parser (deliberately not using skbio here,
        # since sequences are re-emitted verbatim rather than validated).
        with open(self.homologs_fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        records.append((current_id, ''.join(current_seq)))
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
        if current_id:
            records.append((current_id, ''.join(current_seq)))

        self.print_progress(f"Parsed {len(records)} sequences")

        # Legacy kwarg: fold any smiles_code= into the entity list.
        entities = list(entities or [])
        if smiles_code:
            entities = legacy_smiles_entities(smiles_code, self) + entities

        if entities:
            desc = ', '.join(f"{e['type']}:{_abbrev(e['value'])}"
                             + (f" x{e['copies']}" if e.get('copies', 1) > 1 else '')
                             for e in entities)
            self.print_progress(f"Co-folding entities: {desc}")

        # Validate before any GPU time is spent: a bad SMILES or an unknown
        # CCD code otherwise only surfaces after N failed Boltz runs.
        if validate and entities:
            problems = validate_entities(entities)
            for p in problems:
                self.print_error(f"Invalid entity: {p}")
            if problems:
                raise SystemExit("ERROR: fix the entities above before running boltz-model.")
            self.print_progress(f"Validated {len(entities)} co-folding entities")

        if msa_map:
            n = sum(1 for pid, _ in records if pid in msa_map)
            self.print_progress(f"Reusing precomputed MSAs for {n}/{len(records)} homologs")

        if fmt == 'yaml':
            template = load_yaml_template(yaml_template)
            written = 0
            for i, (prot_id, seq) in enumerate(records, 1):
                path = os.path.join(output_folder, f"{prot_id}.yaml")
                doc = render_yaml_for_homolog(template, seq, yaml_target_chain, entities,
                                              msa_path=(msa_map or {}).get(prot_id))
                write_yaml(doc, path)
                written += 1
                if i % 20 == 0 or i == len(records):
                    self.print_progress(f"Generated {i}/{len(records)} YAML files")
        else:
            written = 0
            for i, (prot_id, seq) in enumerate(records, 1):
                path = os.path.join(output_folder, f"{prot_id}.fasta")
                ids = chain_ids()
                msa_path = (msa_map or {}).get(prot_id)
                header = f">{next(ids)}|protein"
                if msa_path:
                    header += f"|{msa_path}"
                with open(path, 'w') as f:
                    f.write(f"{header}\n{seq}\n")
                    for ent in entities:
                        for _ in range(max(1, int(ent.get('copies', 1)))):
                            cid = next(ids)
                            header = f">{cid}|{ent['type']}"
                            if ent.get('msa'):
                                header += f"|{ent['msa']}"
                            f.write(f"{header}\n{ent['value']}\n")
                written += 1
                if i % 20 == 0 or i == len(records):
                    self.print_progress(f"Generated {i}/{len(records)} FASTA files")

        if fmt == 'yaml':
            n_chains = (len(template.get('sequences', []))
                        + sum(max(1, int(e.get('copies', 1))) for e in entities))
        else:
            n_chains = 1 + sum(max(1, int(e.get('copies', 1))) for e in entities)
        self.print_success(
            f"Generated {written} {fmt.upper()} files, {n_chains} chain(s) each", start_time)

    def perform_boltz_modeling(self, input_folder="02_boltz_input", boltz_results_folder='03_boltz_results',
                               protein_only_folder='04_structure_models', boltz_args=None,
                               use_msa_server=True, log_folder=None):
        """
        Stage 3: run `boltz predict` on every file produced by
        generate_boltz_input(), then convert each resulting CIF model to a
        protein-only PDB via mdtraj. Does not use any self.* reference
        attributes -- only the folder paths below.

        Changes from homologic_gs.py, all of which cost real runs there:
          - the Boltz binary is resolved to this environment's own (see
            self.tool) and logged, instead of whatever `$PATH` offered;
          - no os.chdir(): Boltz is given cwd= instead, so a crash mid-loop
            can no longer leave the process in the wrong directory and
            break every relative path after it;
          - the command is a list, not a shell string, so a homolog ID
            containing a shell metacharacter cannot misfire;
          - Boltz's stdout/stderr go to a per-homolog log file and the tail
            is echoed on failure, instead of being captured and discarded
            (which left only an exit code to debug from);
          - with --resume, homologs that already have a model are skipped
            rather than recomputed.
        """
        start_time = datetime.now()
        self.print_progress(f"Starting Boltz modeling: {input_folder}", start_time)

        if not os.path.isdir(input_folder):
            self.print_error(f"Missing input folder: {input_folder}")
            return

        input_files = sorted(f for f in os.listdir(input_folder)
                             if f.endswith(('.fasta', '.fa', '.yaml', '.yml')))
        if not input_files:
            self.print_progress("No Boltz input files → skipping")
            return

        self.print_progress(f"Found {len(input_files)} input files")
        boltz = self.tool('boltz')
        self.safe_create_output_folder(boltz_results_folder)
        self.safe_create_output_folder(protein_only_folder)

        log_folder = log_folder or os.path.join(boltz_results_folder, 'logs')
        os.makedirs(log_folder, exist_ok=True)

        successful = skipped = 0
        for i, fn in enumerate(input_files, 1):
            stem = os.path.splitext(fn)[0]
            if i % 10 == 0 or i == len(input_files):
                self.print_progress(f"Processing {i}/{len(input_files)}")

            # --resume: a finished prediction leaves a model CIF behind.
            if self.resume and glob.glob(os.path.join(
                    boltz_results_folder, '**', f'{stem}_model_0.cif'), recursive=True):
                skipped += 1
                continue

            src = os.path.join(input_folder, fn)
            dst = os.path.join(boltz_results_folder, fn)
            shutil.copy2(src, dst)

            cmd = [boltz, 'predict', fn]
            if use_msa_server:
                cmd.append('--use_msa_server')
            cmd += list(boltz_args or [])

            log_path = os.path.join(os.path.abspath(log_folder), f'{stem}.boltz.log')
            try:
                with open(log_path, 'w') as lf:
                    lf.write(f"# {' '.join(cmd)}\n# cwd={boltz_results_folder}\n\n")
                    lf.flush()
                    subprocess.run(cmd, cwd=boltz_results_folder, check=True,
                                   stdout=lf, stderr=subprocess.STDOUT, text=True)
                successful += 1
            except subprocess.CalledProcessError as e:
                self.print_error(f"Boltz failed {fn} (code {e.returncode}) → {log_path}")
                self.print_error(tail_of(log_path, 20))
            finally:
                try:
                    os.remove(dst)
                except OSError:
                    pass

        attempted = len(input_files) - skipped
        msg = f"Boltz: {successful}/{attempted} attempted successful"
        if skipped:
            msg += f", {skipped} already done and skipped (--resume)"
        self.print_progress(msg)

        # Prefer Boltz's primary model-0 output; fall back to any CIF found
        # if that naming pattern isn't present. dict.fromkeys() dedupes
        # while preserving order.
        cif_files = glob.glob(os.path.join(boltz_results_folder, "**/*_model_0.cif"), recursive=True) or \
                    glob.glob(os.path.join(boltz_results_folder, "**/*.cif"), recursive=True)
        cif_files = list(dict.fromkeys(cif_files))

        if not cif_files:
            self.print_error("No CIF files found after modeling")
            return

        self.print_progress(f"Found {len(cif_files)} CIF files")
        processed = 0
        for i, cif in enumerate(cif_files, 1):
            if i % 10 == 0 or i == len(cif_files):
                self.print_progress(f"Converting {i}/{len(cif_files)}")
            pid = os.path.basename(cif).split("_model_0")[0]
            out_pdb = os.path.join(protein_only_folder, f"{pid}.pdb")
            traj = md.load(cif)                  # CIF -> PDB (round-trip through disk)
            traj.save_pdb(out_pdb)
            traj = md.load(out_pdb)
            sel = traj.topology.select("protein")  # strip any co-folded ligand/hetero atoms
            protein_traj = traj.atom_slice(sel)
            protein_traj.save_pdb(out_pdb)
            processed += 1

        self.print_success(f"Generated {processed}/{len(cif_files)} protein-only PDBs", start_time)

    # -- Sweep ------------------------------------------------------------
    #
    # Deliberately terminal: it produces an ensemble of Boltz models per
    # homolog and stops there. It is NOT part of run-all and feeds none of
    # stages 4-11, which each assume exactly one model per homolog.

    SWEEP_MAX_HOMOLOGS = 5

    def perform_boltz_sweep(self, input_folder='02_boltz_input',
                            output_folder='03_boltz_sweep', seeds=(1, 2, 3, 4, 5),
                            diffusion_samples=100, boltz_args=None,
                            allow_msa_bootstrap=False, assume_yes=False,
                            ensemble_folder=None, dry_run=False):
        """
        Run `boltz predict` over several random seeds per homolog, each
        producing `diffusion_samples` structures, and collect everything into
        an ensemble folder.

        Two hard conditions, both about not melting the public MSA server:

          1. At most SWEEP_MAX_HOMOLOGS (5) homologs per sweep.
          2. Every homolog must have a precomputed MSA -- supplied via
             --use-custom-msa / --reuse-msa-from at the boltz-input stage, so
             it is baked into the input file. --allow-msa-bootstrap relaxes
             this to Boltz generating each homolog's MSA once, on its first
             seed, which the remaining seeds then reuse: that is one server
             call per homolog (<=5), never one per seed.

        Layout, following bin/boltz-2_multi_seed_run/run_boltz_sweep.py:
            <output_folder>/<stem>/seed<N>/boltz_results_<stem>/...
            <output_folder>/ensemble/<stem>/<stem>_seed<N>_<model>.cif
            <output_folder>/manifest.jsonl   -- one row per model, with the
                                                confidence metrics Boltz wrote
        The per-seed subdirectory keeps the input filename identical across
        seeds, which matters because Boltz derives its output folder name
        from that filename.
        """
        start_time = datetime.now()
        self.print_progress(f"Starting Boltz sweep: {input_folder}", start_time)

        if not os.path.isdir(input_folder):
            self.print_error(f"Missing input folder: {input_folder}")
            return

        input_files = sorted(f for f in os.listdir(input_folder)
                             if f.endswith(('.fasta', '.fa', '.yaml', '.yml')))
        if not input_files:
            self.print_progress("No Boltz input files → skipping")
            return

        seeds = list(seeds)

        # -- Hard condition 1: homolog count -------------------------------
        if len(input_files) > self.SWEEP_MAX_HOMOLOGS:
            raise SystemExit(
                f"ERROR: a sweep runs at most {self.SWEEP_MAX_HOMOLOGS} homologs; "
                f"'{input_folder}' has {len(input_files)}.\n"
                f"       {len(input_files)} homologs x {len(seeds)} seeds x "
                f"{diffusion_samples} samples = "
                f"{len(input_files) * len(seeds) * diffusion_samples} structures, and "
                f"one MSA-server call per homolog.\n"
                f"       Generate inputs for a smaller homolog set (boltz-input with a "
                f"trimmed --homologs-fasta) and sweep that.")

        # -- Hard condition 2: precomputed MSA -----------------------------
        without = [f for f in input_files
                   if not input_has_custom_msa(os.path.join(input_folder, f))]
        if without and not allow_msa_bootstrap:
            raise SystemExit(
                f"ERROR: a sweep needs a precomputed MSA for every homolog; "
                f"{len(without)} of {len(input_files)} inputs have none:\n"
                f"       {', '.join(without[:5])}"
                f"{' ...' if len(without) > 5 else ''}\n\n"
                f"       Without one, Boltz calls the MSA server once per seed per "
                f"homolog\n"
                f"       ({len(input_files) * len(seeds)} calls here instead of "
                f"{len(input_files)}).\n\n"
                f"       Either regenerate the inputs with a precomputed MSA:\n"
                f"           homologic2.py boltz-input ... --use-custom-msa <alignment>\n"
                f"           homologic2.py boltz-input ... --reuse-msa-from <old run>\n"
                f"       or let Boltz build each homolog's MSA once, on its first seed, "
                f"and reuse it\n"
                f"       for the rest (one call per homolog):\n"
                f"           --allow-msa-bootstrap")

        boltz = self.tool('boltz')
        total = len(input_files) * len(seeds) * diffusion_samples
        self.print_progress(
            f"Sweep: {len(input_files)} homologs x {len(seeds)} seeds "
            f"(={', '.join(map(str, seeds))}) x {diffusion_samples} diffusion samples "
            f"= {total} structures")
        if without:
            self.print_warning(
                f"MSA bootstrap: {len(without)} homolog(s) will each make ONE "
                f"MSA-server call on their first seed; later seeds reuse it")

        if dry_run:
            self.print_progress("Dry run -- what would happen:")
            for fn in input_files:
                stem = os.path.splitext(fn)[0]
                has = input_has_custom_msa(os.path.join(input_folder, fn))
                src = "precomputed MSA" if has else "MSA bootstrapped on seed " + str(seeds[0])
                print(f"    {stem}: seeds {', '.join(map(str, seeds))} "
                      f"x {diffusion_samples} samples  [{src}]")
            print(f"    -> {len(input_files) * len(seeds)} boltz predict calls, "
                  f"{total} structures")
            print(f"    -> ensemble in {ensemble_folder or os.path.join(output_folder, 'ensemble')}")
            print(f"    -> manifest at {os.path.join(output_folder, 'manifest.jsonl')}")
            return

        # This is a long GPU job -- confirm before starting one.
        if not assume_yes:
            raise SystemExit(
                f"ERROR: this sweep would generate {total} structures across "
                f"{len(input_files) * len(seeds)} `boltz predict` calls.\n"
                f"       That is hours of GPU time and a large amount of disk.\n"
                f"       Re-run with --yes to confirm, or use --dry-run to preview.")

        self.safe_create_output_folder(output_folder)
        ensemble_folder = ensemble_folder or os.path.join(output_folder, 'ensemble')
        os.makedirs(ensemble_folder, exist_ok=True)
        log_folder = os.path.join(output_folder, 'logs')
        os.makedirs(log_folder, exist_ok=True)

        successful = failed = skipped = 0
        collected = 0

        for fn in input_files:
            stem, suffix = os.path.splitext(fn)
            is_yaml = suffix.lower() in ('.yaml', '.yml')
            base_content = Path(os.path.join(input_folder, fn)).read_text()
            reuse_msa = None            # set after a bootstrap seed completes

            for i, seed in enumerate(seeds):
                seed_dir = os.path.join(output_folder, stem, f'seed{seed}')
                if self.resume and seed_is_complete(seed_dir, stem, diffusion_samples):
                    skipped += 1
                    if reuse_msa is None:
                        reuse_msa = find_generated_msa(seed_dir, stem)
                    continue

                os.makedirs(seed_dir, exist_ok=True)
                content = base_content
                if reuse_msa:
                    content = with_msa_path(base_content, reuse_msa, is_yaml)
                seed_input = os.path.join(seed_dir, fn)
                Path(seed_input).write_text(content)

                needs_server = not input_has_custom_msa(seed_input)
                cmd = [boltz, 'predict', fn,
                       '--seed', str(seed),
                       '--diffusion_samples', str(diffusion_samples)]
                if needs_server:
                    cmd.append('--use_msa_server')
                cmd += list(boltz_args or [])

                self.print_progress(f"{stem} seed {seed} ({i + 1}/{len(seeds)})")
                log_path = os.path.join(os.path.abspath(log_folder),
                                        f'{stem}.seed{seed}.boltz.log')
                try:
                    with open(log_path, 'w') as lf:
                        lf.write(f"# {' '.join(cmd)}\n# cwd={seed_dir}\n\n")
                        lf.flush()
                        subprocess.run(cmd, cwd=seed_dir, check=True,
                                       stdout=lf, stderr=subprocess.STDOUT, text=True)
                    successful += 1
                except subprocess.CalledProcessError as e:
                    failed += 1
                    self.print_error(f"{stem} seed {seed} failed (code {e.returncode}) "
                                     f"→ {log_path}")
                    self.print_error(tail_of(log_path, 20))
                    continue
                # The seed's exact input is deliberately kept: it records which
                # MSA that seed actually used, which is the one thing you cannot
                # reconstruct from the output afterwards.

                # Bootstrap: hand this homolog's freshly-built MSA to its
                # remaining seeds, so the server is called once, not per seed.
                if reuse_msa is None and needs_server:
                    reuse_msa = find_generated_msa(seed_dir, stem)
                    if reuse_msa:
                        self.print_progress(f"{stem}: MSA cached for the remaining "
                                            f"seeds → {reuse_msa}")
                    else:
                        self.print_warning(
                            f"{stem}: no MSA found after seed {seed}; the remaining "
                            f"seeds will each call the MSA server")

            # Collect this homolog's models into its ensemble folder.
            dest = os.path.join(ensemble_folder, stem)
            os.makedirs(dest, exist_ok=True)
            for seed in seeds:
                seed_dir = os.path.join(output_folder, stem, f'seed{seed}')
                for cif in sorted(glob.glob(os.path.join(seed_dir, '**', '*.cif'),
                                            recursive=True)):
                    model = os.path.basename(cif).replace(f'{stem}_', '')
                    out = os.path.join(dest, f'{stem}_seed{seed}_{model}')
                    if not os.path.exists(out):
                        shutil.copy2(cif, out)
                        collected += 1

        stems = [os.path.splitext(f)[0] for f in input_files]
        manifest, n_rows = build_sweep_manifest(output_folder, ensemble_folder,
                                                stems, seeds)
        self.print_progress(f"Manifest: {n_rows} models indexed → {manifest}")

        msg = (f"Sweep: {successful} succeeded, {failed} failed"
               + (f", {skipped} skipped (--resume)" if skipped else "")
               + f"; {collected} models collected in {ensemble_folder}")
        self.print_success(msg, start_time)

    def calculate_structure_homology(self, input_csv_file='01_sequence_homology_results/sequence_homology.csv',
                                     input_folder='04_structure_models', output_folder='05_structure_homology_results'):
        """
        Stage 4: TM-align every predicted structure in input_folder against
        self.reference_pdb (REQUIRED), appending rmsd/tm_score columns onto
        the Stage 1 CSV. Requires the `TMalign` binary on PATH.
        The default input_csv_file matches calculate_sequence_homology()'s
        default output folder; homologic_gs.py's two defaults disagreed on
        the '01_' prefix, so using both defaults silently failed.
        """
        start_time = datetime.now()
        self.print_progress("Starting structural homology (TM-align)", start_time)
        self.safe_create_output_folder(output_folder)

        if not os.path.exists(self.reference_pdb):
            self.print_error(f"Reference PDB missing: {self.reference_pdb}")
            return pd.DataFrame()

        df = pd.read_csv(input_csv_file)
        df['rmsd'] = np.nan
        df['tm_score'] = np.nan

        pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
        pdb_dict = {Path(f).stem: os.path.join(input_folder, f) for f in pdb_files}

        for idx, row in df.iterrows():
            tid = row.get('homolog_id', row.get('reference_id', 'n.d.'))
            if tid == 'n.d.':
                continue
            acc = Path(str(tid)).stem
            pdb_path = pdb_dict.get(acc)
            if not pdb_path or not os.path.exists(pdb_path):
                continue
            try:
                res = subprocess.run([self.tool('TMalign'), self.reference_pdb, pdb_path],
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
                rmsd = tm = None
                # TMalign prints RMSD and TM-score on separate lines of its
                # stdout report; scrape both with regex rather than parsing
                # its full fixed-width output format.
                for line in res.stdout.splitlines():
                    if 'RMSD=' in line:
                        m = re.search(r'RMSD[=:]\s*([\d.]+)', line)
                        if m:
                            rmsd = float(m.group(1))
                    if 'TM-score=' in line:
                        m = re.search(r'TM-score[=:]\s*([\d.]+)', line) or re.search(r'TM-score=\s*([\d.]+)', line)
                        if m:
                            tm = float(m.group(1))
                if rmsd is not None and tm is not None:
                    df.at[idx, 'rmsd'] = rmsd
                    df.at[idx, 'tm_score'] = tm
            except subprocess.CalledProcessError:
                continue

        out_csv = os.path.join(output_folder, 'structure_homology.csv')
        df.to_csv(out_csv, index=False)
        self.print_success(f"Results saved → {out_csv}", start_time)
        return df

    def superpose_structures(self, input_folder='04_structure_models', output_folder='06_superposed_structures'):
        """
        Stage 5: CA-atom superposition (PyMOL `align`) of every predicted
        structure in input_folder onto self.reference_pdb (REQUIRED).
        Requires a PyMOL installation importable as `from pymol import cmd`
        in this Python environment (imported lazily here, not at module
        top-level).
        """
        start_time = datetime.now()
        self.print_progress("Starting superposition to reference", start_time)
        self.safe_create_output_folder(output_folder)

        pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
        self.print_progress(f"Found {len(pdb_files)} PDB files")

        from pymol import cmd  # imported here so stages that don't need PyMOL can still run without it installed
        successful = 0

        for i, fn in enumerate(pdb_files, 1):
            in_path = os.path.join(input_folder, fn)
            out_path = os.path.join(output_folder, f"{os.path.splitext(fn)[0]}_superposed.pdb")
            try:
                cmd.reinitialize()
                cmd.load(self.reference_pdb, "reference")
                cmd.load(in_path, "target")
                cmd.align("target and name CA", "reference and name CA")
                cmd.save(out_path, "target")
                cmd.delete("all")
                successful += 1
                if i % 10 == 0 or i == len(pdb_files):
                    self.print_progress(f"Superposed {i}/{len(pdb_files)}")
            except Exception as e:
                self.print_warning(f"Superposition failed {fn}: {type(e).__name__}")

        self.print_success(f"Superposed {successful}/{len(pdb_files)} structures", start_time)

    class SelectBindingSiteResidues:
        # Bio.PDB.PDBIO Select subclass: writes out only the residues that
        # were found within distance_cutoff of a reference ligand atom.
        def __init__(self, binding_residues):
            self.binding_residues = binding_residues
        def accept_residue(self, residue):
            return residue in self.binding_residues
        def accept_chain(self, chain):
            return True
        def accept_model(self, model):
            return True
        def accept_atom(self, atom):
            return True

    def extract_reference_binding_site(self, distance_cutoff=6.0):
        """
        Stage 6: identify every reference protein residue with at least one
        atom within distance_cutoff (Angstrom) of any atom in
        self.reference_ligands_pdb, and write those residues out as
        "<reference_pdb-stem>_bindingsite.pdb". Both self.reference_pdb and
        self.reference_ligands_pdb are REQUIRED. This is not written into a
        --output-folder like most other stages; the output path is derived
        directly from self.reference_pdb's location.
        """
        start_time = datetime.now()
        self.print_progress(f"Extracting reference binding site ({distance_cutoff}Å)", start_time)
        parser = PDBParser(QUIET=True)
        ligand_struct = parser.get_structure("ligand", self.reference_ligands_pdb)
        ligand_atoms = list(ligand_struct.get_atoms())
        out_path = f"{os.path.splitext(self.reference_pdb)[0]}_bindingsite.pdb"
        prot_struct = parser.get_structure("protein", self.reference_pdb)
        prot_atoms = [a for a in prot_struct.get_atoms() if is_aa(a.parent)]  # exclude waters/hetero groups
        ns = NeighborSearch(prot_atoms)
        binding_res = set()
        for la in ligand_atoms:
            for a in ns.search(la.coord, distance_cutoff):
                binding_res.add(a.get_parent())  # a.get_parent() == the residue owning atom `a`
        io = PDBIO()
        io.set_structure(prot_struct)
        io.save(out_path, self.SelectBindingSiteResidues(binding_res))
        self.print_success(f"Reference binding site saved → {out_path}", start_time)

    def extract_homolog_binding_sites(self, distance_cutoff=6.0, input_folder='06_superposed_structures',
                                      output_folder='07_homolog_bindingsites'):
        """
        Stage 7: same idea as extract_reference_binding_site(), but applied
        to every already-superposed homolog structure in input_folder,
        using the reference ligand coordinates (self.reference_ligands_pdb,
        REQUIRED) as the fixed reference frame (this only works because the
        homolog structures were superposed onto the reference in Stage 5).
        """
        start_time = datetime.now()
        self.print_progress(f"Extracting homolog binding sites ({distance_cutoff}Å)", start_time)
        self.safe_create_output_folder(output_folder)
        parser = PDBParser(QUIET=True)
        ligand_struct = parser.get_structure("ligand", self.reference_ligands_pdb)
        ligand_atoms = list(ligand_struct.get_atoms())
        pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
        self.print_progress(f"Found {len(pdb_files)} superposed structures")
        successful = 0
        for i, fn in enumerate(pdb_files, 1):
            in_path = os.path.join(input_folder, fn)
            out_path = os.path.join(output_folder, f"{os.path.splitext(fn)[0]}_bindingsite.pdb")
            try:
                prot_struct = parser.get_structure("protein", in_path)
                prot_atoms = [a for a in prot_struct.get_atoms() if is_aa(a.parent)]
                ns = NeighborSearch(prot_atoms)
                binding_res = set()
                for la in ligand_atoms:
                    for a in ns.search(la.coord, distance_cutoff):
                        binding_res.add(a.get_parent())
                if not binding_res:
                    continue
                io = PDBIO()
                io.set_structure(prot_struct)
                io.save(out_path, self.SelectBindingSiteResidues(binding_res))
                successful += 1
                if i % 10 == 0 or i == len(pdb_files):
                    self.print_progress(f"Processed {i}/{len(pdb_files)}")
            except Exception as e:
                self.print_warning(f"Binding-site extraction failed for "
                                   f"{os.path.basename(pdb_path)}: {type(e).__name__}: {e}")
                continue
        if successful < len(pdb_files):
            self.print_warning(f"{len(pdb_files) - successful} of {len(pdb_files)} "
                               f"structures yielded no binding site")
        self.print_success(f"Extracted {successful}/{len(pdb_files)} binding sites", start_time)

    def superpose_binding_sites(self, csv_file='05_structure_homology_results/structure_homology.csv',
                                input_folder='07_homolog_bindingsites', output_folder='08_homolog_bindingsites_superposed',
                                min_residues_total=5, min_aligned_atoms=5):
        """
        Stage 8: refine the Stage 7 binding-site extraction with a
        second, binding-site-local PyMOL CA superposition (`cmd.super`)
        against the reference binding site (derived from self.reference_pdb,
        REQUIRED). min_residues_total/min_aligned_atoms are quality gates:
        too few residues -> skip; too few atoms actually aligned -> RMSD is
        reported as NaN rather than trusted. Requires PyMOL.
        """
        start_time = datetime.now()
        self.print_progress("Starting refined binding site superposition", start_time)
        self.safe_create_output_folder(output_folder)

        df = pd.read_csv(csv_file)
        for col in ['rmsd_bindingsite', 'n_bindingsite_residues', 'n_aligned_residues', 'aligned_fraction']:
            if col not in df.columns:
                df[col] = np.nan

        ref_bs_pdb = f"{os.path.splitext(self.reference_pdb)[0]}_bindingsite.pdb"
        ref_stem = Path(self.reference_pdb).stem
        if not os.path.exists(ref_bs_pdb):
            self.print_error(f"Reference binding site missing: {ref_bs_pdb}")
            return

        from pymol import cmd
        processed = valid_rmsd = 0

        for idx, row in df.iterrows():
            tid = row.get('homolog_id', row.get('reference_id', 'n.d.'))
            if tid == 'n.d.':
                continue
            acc = Path(str(tid)).stem
            in_pdb = os.path.join(input_folder, f"{acc}_superposed_bindingsite.pdb")
            if not os.path.exists(in_pdb):
                continue
            processed += 1
            if acc == ref_stem:
                df.at[idx, 'rmsd_bindingsite'] = 0.0  # reference-vs-itself trivially has zero RMSD
                continue
            try:
                cmd.reinitialize()
                cmd.load(ref_bs_pdb, "ref_bs")
                cmd.load(in_pdb, "tgt_bs")
                cmd.select("ref_ca", "name CA and ref_bs")
                cmd.select("tgt_ca", "name CA and tgt_bs")
                n_tgt = cmd.count_atoms("tgt_ca")
                df.at[idx, 'n_bindingsite_residues'] = n_tgt
                if n_tgt < min_residues_total:
                    cmd.delete("all")
                    continue
                # cmd.super() is a sequence-independent structural superposition
                # (more tolerant of local sequence divergence than cmd.align).
                rmsd_res = cmd.super("tgt_ca", "ref_ca")
                rmsd = rmsd_res[0] if isinstance(rmsd_res, tuple) and len(rmsd_res) >= 2 else np.nan
                n_aln = rmsd_res[1] if isinstance(rmsd_res, tuple) and len(rmsd_res) >= 2 else 0
                df.at[idx, 'n_aligned_residues'] = n_aln
                df.at[idx, 'aligned_fraction'] = n_aln / n_tgt if n_tgt > 0 else np.nan
                if n_aln < min_aligned_atoms:
                    rmsd = np.nan
                out_fn = f"{acc}_superposed_bindingsite_refined.pdb"
                cmd.save(os.path.join(output_folder, out_fn), "tgt_bs")
                df.at[idx, 'rmsd_bindingsite'] = rmsd
                if not np.isnan(rmsd):
                    valid_rmsd += 1
                cmd.delete("all")
            except Exception as e:
                self.print_warning(f"Binding-site superposition failed for {acc}: "
                                   f"{type(e).__name__}: {e}")
                df.at[idx, 'rmsd_bindingsite'] = np.nan
                df.at[idx, 'n_aligned_residues'] = np.nan
                df.at[idx, 'aligned_fraction'] = np.nan

        df = df.fillna('n.d.')
        out_csv = os.path.join(output_folder, 'bindingsite_geometry_homology.csv')
        df.to_csv(out_csv, index=False)
        self.print_success(f"Refinement done → {processed} processed, {valid_rmsd} valid RMSDs", start_time)

    def analyze_binding_sites(self,
                              csv_file='08_homolog_bindingsites_superposed/bindingsite_geometry_homology.csv',
                              input_folder='08_homolog_bindingsites_superposed',
                              bindingsite_metasequence_folder='09_bindingsite_metasequences',
                              bindingsite_similarity_results='10_bindingsite_similarity_results',
                              boltz_results_base='03_boltz_results',
                              tolerated_misalignment=1,
                              substitution_matrix_name='BLOSUM62'):
        """
        Stage 9: for every homolog, derive a "metasequence" (the one-letter
        residues making up its refined binding site, spatially matched
        residue-for-residue to the reference binding site within
        tolerated_misalignment Angstrom -- mdtraj's native nanometres are
        converted, see extract_metasequence) and compute:
          - a binding-site sequence-identity score, normalised against the
            reference-vs-itself score,
          - full-model and binding-site-local pLDDT summary statistics,
            read from the Boltz .npz outputs under boltz_results_base
            (REQUIRES that perform_boltz_modeling() was already run there).
        self.reference_pdb (REQUIRED) is used only to locate the reference
        binding-site PDB written by extract_reference_binding_site().
        """
        start_time = datetime.now()
        self.print_progress("Starting binding site metasequence & similarity analysis", start_time)

        ref_bs_pdb = f"{os.path.splitext(self.reference_pdb)[0]}_bindingsite.pdb"
        if not os.path.exists(ref_bs_pdb):
            self.print_error(f"Reference binding site missing: {ref_bs_pdb}")
            return

        aa_mapping = {
            'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
            'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'
        }

        def one_letter(res_dict):
            # 3-letter residue name -> 1-letter code; unrecognised (e.g.
            # non-standard) residues silently drop out as ''.
            return ''.join(aa_mapping.get(str(r)[:3], '') for r in res_dict.values())

        def extract_metasequence(ref_pdb, tgt_pdb, max_dist):
            # For each reference binding-site CA atom, find its nearest
            # target CA atom (nearest-neighbour matching, not a real
            # structural alignment); target residues whose nearest reference
            # CA is farther than max_dist are dropped as "non-equivalent".
            ref = md.load(ref_pdb)
            tgt = md.load(tgt_pdb)
            tgt_res = {i:res for i,res in enumerate(tgt.topology.residues)}
            orig_len = len(tgt_res)
            ref_ca_idx = [a.index for a in ref.topology.atoms if a.name == 'CA']
            tgt_ca_idx = [a.index for a in tgt.topology.atoms if a.name == 'CA']
            # mdtraj stores coordinates in NANOMETRES; every other distance in
            # this pipeline (Bio.PDB binding-site cutoff, TM-align RMSD, PyMOL
            # RMSD, pyKVFinder probe/volume, Open3D ICP) is in ANGSTROMS, and
            # max_dist is documented as Angstrom. Convert here so the units
            # match: without this, tolerated_misalignment=1.0 silently means
            # 1 nm = 10 A, ten times looser than intended.
            ref_xyz = ref.xyz[0][ref_ca_idx] * 10.0
            tgt_xyz = tgt.xyz[0][tgt_ca_idx] * 10.0
            dists = np.linalg.norm(ref_xyz[:,None] - tgt_xyz[None,:], axis=-1)  # pairwise, Angstrom
            closest = np.argmin(dists, axis=1)
            min_d = dists[np.arange(len(ref_ca_idx)), closest]
            # Take the target residue that is actually CLOSEST to each
            # reference CA -- closest[i], not i. Indexing tgt_ca_idx by the
            # reference's own position matched residues by rank rather than
            # proximity, and raised IndexError whenever the target site had
            # fewer CA atoms than the reference.
            matched_tgt = {tgt.topology.atom(tgt_ca_idx[closest[i]]).residue
                           for i, d in enumerate(min_d) if d < max_dist}
            keep = [a.index for a in tgt.topology.atoms if a.residue in matched_tgt]
            tgt = tgt.atom_slice(keep)
            ref_res = {i:res for i,res in enumerate(ref.topology.residues)}
            tgt_res = {i:res for i,res in enumerate(tgt.topology.residues)}
            return one_letter(ref_res), one_letter(tgt_res), orig_len - len(tgt_res)

        def write_fasta(seq, path, sid):
            with open(path, 'w') as f:
                f.write(f">{sid}\n{seq}\n")

        def plddt_stats(values):
            # Summary statistics over a set of per-residue pLDDT scores.
            if not values.size:
                return {k:np.nan for k in ['n_res_plddt','median_plddt','avg_plddt','std_plddt','iqr_plddt','p10_plddt']}
            return {
                'n_res_plddt': len(values),
                'median_plddt': np.median(values),
                'avg_plddt': np.mean(values),
                'std_plddt': np.std(values),
                'iqr_plddt': np.percentile(values,75) - np.percentile(values,25),
                'p10_plddt': np.percentile(values,10),
            }

        def get_plddt(target_id, mode='full'):
            # Reads Boltz's per-residue pLDDT array from the fixed output
            # path Boltz itself creates under boltz_results_base (see
            # perform_boltz_modeling). mode='full' summarises the whole
            # model; mode='bindingsite' summarises only binding-site
            # residues (matched back to array indices via PDB residue
            # numbering, 1-based -> 0-based).
            pid = os.path.splitext(str(target_id))[0]
            plddt_path = os.path.join(boltz_results_base, f'boltz_results_{pid}', 'predictions', pid, f'plddt_{pid}_model_0.npz')
            if not os.path.exists(plddt_path):
                return None
            plddt = np.load(plddt_path)['plddt']
            if mode == 'full':
                return plddt_stats(plddt)
            bs_pdb = os.path.join(input_folder, f"{pid}_superposed_bindingsite_refined.pdb")
            if not os.path.exists(bs_pdb):
                return None
            struct = PDBParser(QUIET=True).get_structure("bs", bs_pdb)
            bs_plddt = []
            for model in struct:
                for chain in model:
                    for res in chain:
                        if not is_aa(res): continue
                        idx = res.id[1] - 1
                        if 0 <= idx < len(plddt):
                            bs_plddt.append(plddt[idx])
            return plddt_stats(np.array(bs_plddt)) if bs_plddt else None

        self.safe_create_output_folder(bindingsite_metasequence_folder)
        self.safe_create_output_folder(bindingsite_similarity_results)

        df = pd.read_csv(csv_file)
        plddt_cols = ['n_res_plddt','median_plddt','avg_plddt','std_plddt','iqr_plddt','p10_plddt']
        for c in plddt_cols:
            df[c] = np.nan
            df[f'{c}_bs'] = np.nan
        df['bindingsite_similarity_score'] = np.nan
        df['non_equivalent_residues'] = np.nan

        processed = 0
        failed = []
        for idx, row in df.iterrows():
            tid = row.get('homolog_id', row.get('reference_id', 'n.d.'))
            if tid == 'n.d.':
                continue
            tgt_pdb = os.path.join(input_folder, f"{tid}_superposed_bindingsite_refined.pdb")
            if not os.path.exists(tgt_pdb):
                continue
            try:
                ref_seq, tgt_seq, non_eq = extract_metasequence(ref_bs_pdb, tgt_pdb, tolerated_misalignment)
                write_fasta(tgt_seq, os.path.join(bindingsite_metasequence_folder, f"{tid}.fasta"), tid)
                stats = self.align_and_compute_identity_stats(Protein(ref_seq), Protein(tgt_seq), substitution_matrix_name, gap_open=10, gap_extend=1)
                if stats:
                    max_stats = self.align_and_compute_identity_stats(Protein(ref_seq), Protein(ref_seq), substitution_matrix_name, gap_open=10, gap_extend=1)
                    max_sc = max_stats['identity'] if max_stats else 1.0
                    df.at[idx, 'bindingsite_similarity_score'] = stats['identity'] / max_sc if max_sc else np.nan
                df.at[idx, 'non_equivalent_residues'] = non_eq
                full_s = get_plddt(tid, 'full')
                if full_s:
                    for k,v in full_s.items():
                        df.at[idx, k] = round(v,3) if pd.notna(v) else np.nan
                bs_s = get_plddt(tid, 'bindingsite')
                if bs_s:
                    for k,v in bs_s.items():
                        df.at[idx, f'{k}_bs'] = round(v,3) if pd.notna(v) else np.nan
                processed += 1
                if processed % 10 == 0:
                    self.print_progress(f"Analyzed {processed} binding sites")
            except Exception as e:
                # Previously `except Exception: continue`, which dropped
                # homologs with no message at all -- a crash here left the
                # stage reporting success on a silently halved dataset.
                failed.append(tid)
                self.print_warning(f"Binding-site analysis failed for {tid}: "
                                   f"{type(e).__name__}: {e}")
                continue

        out_csv = os.path.join(bindingsite_similarity_results, 'bindingsite_homology.csv')
        df.to_csv(out_csv, index=False)
        if failed:
            self.print_warning(
                f"{len(failed)} of {processed + len(failed)} binding sites failed and "
                f"have no score: {', '.join(failed[:5])}"
                f"{' ...' if len(failed) > 5 else ''}")
        self.print_success(f"Processed {processed} binding sites → {out_csv}", start_time)
        return df

    def analyze_reference_cavity_properties(self, input='reference_bindingsite.pdb', output_results_folder='reference_cavity_analysis',
                                            probe_out=4.0, volume_cutoff=100.0):
        """
        Stage 10: pyKVFinder cavity detection on the reference binding site.
        Unlike most other stages, `input` here is a plain literal default
        ('reference_bindingsite.pdb') -- NOT derived from self.reference_pdb
        -- so it must be pointed explicitly at the file that
        extract_reference_binding_site() wrote out (self.reference_pdb's
        stem + '_bindingsite.pdb') unless that file happens to already be
        named 'reference_bindingsite.pdb' in the current directory.
        """
        start_time = datetime.now()
        self.print_progress("Starting reference cavity analysis", start_time)
        self.safe_create_output_folder(output_results_folder)
        input_path = os.path.join(input)
        results = pyKVFinder.run_workflow(
            input_path, probe_out=probe_out, volume_cutoff=volume_cutoff,
            include_depth=True, include_hydropathy=True, ignore_backbone=True
        )
        results.export_all(
            fn=os.path.join(output_results_folder, "reference_results.toml"),
            output=os.path.join(output_results_folder, "reference_cavities.pdb"),
            include_frequencies_pdf=True,
            pdf=os.path.join(output_results_folder, "reference_barplots.pdf")
        )
        self.print_success("Reference cavity analysis completed", start_time)

    def analyze_cavity_properties(self,
                                  input_folder='08_homolog_bindingsites_superposed',
                                  output_cavities_folder='11_detected_cavities',
                                  output_results_folder='12_cavity_analysis_results',
                                  probe_out=4.0, volume_cutoff=100.0,
                                  csv_input='10_bindingsite_similarity_results/bindingsite_homology.csv',
                                  reference_pdb=None, reference_cavity=None):
        """
        Stage 11: pyKVFinder cavity detection on every homolog binding site
        in input_folder, plus an Open3D ICP (iterative closest point) fit
        of each detected cavity point cloud onto the reference cavity, to
        get a shape-similarity RMSD.

        NOTE ON PARAMETER NAMING: the `reference_pdb` keyword argument here
        is LOCAL to this method and is NOT the same thing as
        self.reference_pdb (the constructor attribute). If left as None
        (the default) it falls back to a path derived from
        self.reference_pdb (REQUIRED in that case); if given explicitly, it
        instead OVERRIDES that derivation and must itself already be a path
        to a reference *binding-site* PDB (not the whole reference
        structure). Similarly, `reference_cavity` defaults to the file
        analyze_reference_cavity_properties() (Stage 10) writes out.
        """
        start_time = datetime.now()
        self.print_progress("Starting cavity detection + ICP", start_time)
        self.safe_create_output_folder(output_cavities_folder)
        self.safe_create_output_folder(output_results_folder)

        ref_pdb = reference_pdb or f"{os.path.splitext(self.reference_pdb)[0]}_bindingsite.pdb"
        ref_cav = reference_cavity or os.path.join("reference_cavity_analysis", "reference_cavities.pdb")

        def read_points(pdb_file):
            # Parse raw XYZ coordinates straight from PDB fixed-width
            # columns (faster than a full structure parser for this purpose).
            coords = []
            with open(pdb_file) as f:
                for line in f:
                    if line.startswith(("ATOM","HETATM")):
                        x,y,z = map(float, [line[30:38],line[38:46],line[46:54]])
                        coords.append([x,y,z])
            return np.array(coords)

        def run_icp(tgt_pdb, ref_pdb):
            # Rigid-body point-cloud alignment of a detected cavity onto the
            # reference cavity; inlier_rmse is the shape-fit quality metric.
            ref_pts = read_points(ref_pdb)
            tgt_pts = read_points(tgt_pdb)
            if len(ref_pts) < 3 or len(tgt_pts) < 3:
                return np.nan, 0
            pcd_ref = o3d.geometry.PointCloud()
            pcd_ref.points = o3d.utility.Vector3dVector(ref_pts)
            pcd_tgt = o3d.geometry.PointCloud()
            pcd_tgt.points = o3d.utility.Vector3dVector(tgt_pts)
            res = o3d.pipelines.registration.registration_icp(
                pcd_tgt, pcd_ref, max_correspondence_distance=2.0,
                init=np.eye(4),
                estimation_method=o3d.pipelines.registration.TransformationEstimationPointToPoint()
            )
            return round(float(res.inlier_rmse),3), len(res.correspondence_set)

        df = pd.read_csv(csv_input) if os.path.exists(csv_input) else pd.DataFrame()

        new_cols = [
            'cavity_volume','cavity_area','cavity_avg_depth','cavity_avg_hydropathy',
            'cavity_frequency_Alipathic_apolar','cavity_frequency_Aromatic',
            'cavity_frequency_Polar_uncharged','cavity_frequency_Negatively_charged',
            'cavity_frequency_Positively_charged','cavity_frequency_Non-standard',
            'cavity_rmsd','cavity_n_points'
        ]
        for c in new_cols:
            if c not in df.columns:
                df[c] = np.nan

        pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
        if not pdb_files:
            self.print_error(f"No PDB files in {input_folder}")
            df.fillna('n.d.').to_csv(os.path.join(output_results_folder, 'bindingsite_cavity_homology.csv'), index=False)
            return df

        self.print_progress(f"Found {len(pdb_files)} PDB files")
        processed = detected = 0

        for i, fn in enumerate(pdb_files, 1):
            tid = Path(fn).stem.split('_')[0]
            in_p = os.path.join(input_folder, fn)
            cav_pdb = os.path.join(output_cavities_folder, f"{tid}_cavities.pdb")
            if i % 5 == 0 or i == len(pdb_files):
                self.print_progress(f"Analyzing {i}/{len(pdb_files)}")
            try:
                res = pyKVFinder.run_workflow(
                    in_p, probe_out=probe_out, volume_cutoff=volume_cutoff,
                    include_depth=True, include_hydropathy=True, ignore_backbone=True
                )
                res.export_all(
                    fn=os.path.join(output_cavities_folder, f"{tid}_results.toml"),
                    output=cav_pdb,
                    include_frequencies_pdf=True,
                    pdf=os.path.join(output_cavities_folder, f"{tid}_barplots.pdf")
                )
                vol_dict = res.volume
                if not vol_dict:
                    continue
                detected += 1
                lbl = max(vol_dict, key=vol_dict.get)  # pyKVFinder can report multiple cavities; keep the largest
                rmsd, ncorr = run_icp(cav_pdb, ref_cav)
                mask = df['homolog_id'].isin([tid, f"{tid}_superposed_bindingsite_refined"])
                if not mask.any():
                    continue
                ridx = df.index[mask].tolist()[0]
                df.loc[ridx, 'cavity_volume'] = vol_dict[lbl]
                df.loc[ridx, 'cavity_area'] = res.area.get(lbl, np.nan)
                df.loc[ridx, 'cavity_avg_depth'] = res.avg_depth.get(lbl, np.nan)
                df.loc[ridx, 'cavity_avg_hydropathy'] = res.avg_hydropathy.get(lbl, np.nan)
                freq = res.frequencies.get(lbl, {}).get("CLASS", {})
                df.loc[ridx, 'cavity_frequency_Alipathic_apolar'] = freq.get("R1", np.nan)
                df.loc[ridx, 'cavity_frequency_Aromatic'] = freq.get("R2", np.nan)
                df.loc[ridx, 'cavity_frequency_Polar_uncharged'] = freq.get("R3", np.nan)
                df.loc[ridx, 'cavity_frequency_Negatively_charged'] = freq.get("R4", np.nan)
                df.loc[ridx, 'cavity_frequency_Positively_charged'] = freq.get("R5", np.nan)
                df.loc[ridx, 'cavity_frequency_Non-standard'] = freq.get("RX", np.nan)
                df.loc[ridx, 'cavity_rmsd'] = rmsd
                df.loc[ridx, 'cavity_n_points'] = ncorr
                processed += 1
            except Exception as e:
                self.print_warning(f"Cavity analysis failed for {acc}: "
                                   f"{type(e).__name__}: {e}")
                continue

        out_csv = os.path.join(output_results_folder, 'bindingsite_cavity_homology.csv')
        df.fillna('n.d.').to_csv(out_csv, index=False)
        self.print_success(f"Cavity+ICP done → {processed}/{len(pdb_files)} processed | {detected} cavities", start_time)
        return df

# =============================================================================
# Command-line interface
# =============================================================================
# Everything above this line is the original HomoLogic class -- unchanged in
# behaviour (only comments/docstrings and a few whitespace-only PEP8 tweaks
# were added). Everything below wires it up to argparse, one subcommand per
# pipeline stage, so each stage can be run independently from the shell.
# Every CLI default below is copied verbatim from the corresponding Python
# method's own default argument, so `homologic_gs.py <stage>` with no extra
# flags behaves identically to calling that method with no arguments.


# ======================================================================
# Command-line interface
# ======================================================================

def add_boltz_arguments(parser, skip=()):
    """Boltz prediction options, shared by boltz-model and run-all."""
    g = parser.add_argument_group(
        'Boltz prediction',
        'The options changed most often are named below; --boltz-arg reaches '
        'every other `boltz predict` flag.')
    g.add_argument('--devices', type=int, help='Number of GPUs for boltz predict.')
    g.add_argument('--seed', type=int, help='Boltz random seed.')
    if 'diffusion-samples' not in skip:
        g.add_argument('--diffusion-samples', type=int,
                       help='Boltz diffusion samples per input.')
    g.add_argument('--recycling-steps', type=int, help='Boltz recycling steps.')
    g.add_argument('--use-potentials', action='store_true',
                   help='Enable Boltz steering potentials.')
    g.add_argument('--no-msa-server', action='store_true',
                   help='Do not pass --use_msa_server (use when MSAs are supplied '
                        'in the input files, or when offline).')
    g.add_argument('--boltz-arg', action='append', default=[], metavar='ARG',
                   help='Repeatable raw passthrough to `boltz predict`. Use the "=" form '
                        'so argparse does not mistake the value for a flag of ours: '
                        '--boltz-arg=--output_format=pdb')
    return parser


def collect_boltz_args(args, skip_seed=False, skip_diffusion_samples=False):
    """
    Turn the flags above into an argv fragment for `boltz predict`.

    The sweep passes its own --seed and --diffusion_samples per call, so it
    skips those here rather than emitting each twice.
    """
    out = []
    pairs = [('--devices', 'devices'), ('--seed', 'seed'),
             ('--diffusion_samples', 'diffusion_samples'),
             ('--recycling_steps', 'recycling_steps')]
    if skip_seed:
        pairs = [p for p in pairs if p[1] != 'seed']
    if skip_diffusion_samples:
        pairs = [p for p in pairs if p[1] != 'diffusion_samples']
    for flag, attr in pairs:
        val = getattr(args, attr, None)
        if val is not None:
            out += [flag, str(val)]
    if getattr(args, 'use_potentials', False):
        out.append('--use_potentials')
    out += list(getattr(args, 'boltz_arg', None) or [])
    return out


# The 11 pipeline stages, in order. 'boltz-sweep' is deliberately absent: it
# is an ensemble generator, not a pipeline step, and stages 4-11 each assume
# exactly one model per homolog.
STAGE_ORDER = [
    'seq-homology', 'boltz-input', 'boltz-model', 'struct-homology',
    'superpose-structures', 'ref-bindingsite', 'homolog-bindingsites',
    'superpose-bindingsites', 'analyze-bindingsites', 'ref-cavity',
    'cavity-properties',
]


class ShortHelpAction(argparse.Action):
    """-h : the terse argparse usage summary."""

    def __init__(self, option_strings, dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings, dest=dest, default=default,
                         nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        parser.print_usage()
        if parser.prog.count(' ') == 0:          # top level
            print(f"\nRun '{parser.prog} --help' for the full manual, "
                  f"or '{parser.prog} <stage> -h' for one stage.")
        else:
            print(f"\nRun '{parser.prog} --help' for what each option does.")
        parser.exit()


class LongHelpAction(argparse.Action):
    """--help : the full manual -- this module's docstring, then every option."""

    def __init__(self, option_strings, dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings, dest=dest, default=default,
                         nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        if parser.prog.count(' ') == 0:      # top level only
            print(__doc__.strip())
            print()
        parser.print_help()
        parser.exit()


def add_help_arguments(parser):
    """Attach the split -h / --help pair to a parser built with add_help=False."""
    parser.add_argument('-h', action=ShortHelpAction,
                        help='Short usage summary.')
    parser.add_argument('--help', action=LongHelpAction,
                        help='Full manual: what every stage does, plus all options.')
    return parser


def build_arg_parser():
    """Assemble the top-level parser and its 11 pipeline-stage subcommands."""
    parser = argparse.ArgumentParser(
        prog='homologic2.py',
        description='HomoLogic v2: homolog structural-analysis pipeline, with a built-in installer.',
        epilog='First run:  python homologic2.py --install  &&  conda activate homologic\n'
               'Per-stage options:  python homologic2.py <stage> --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
    )
    add_help_arguments(parser)

    # -- Environment management. These run under a bare interpreter, before
    #    the environment they build exists, so they must not require any
    #    subcommand and must not touch the scientific stack. ---------------
    env = parser.add_argument_group('environment')
    env.add_argument('--install', action='store_true',
                     help='Build the conda environment and every tool the pipeline needs.')
    env.add_argument('--install-help', action='store_true',
                     help='Print the equivalent manual install commands and exit.')
    env.add_argument('--check', action='store_true',
                     help='Verify an environment: imports + external binaries + GPU.')
    env.add_argument('--env-name', default=ENV_NAME_DEFAULT,
                     help='Conda environment to build or check.')
    env.add_argument('--pinned', action='store_true',
                     help='With --install: reproduce the reference env exactly from '
                          'homologic.yml (linux-64 + matching CUDA only) instead of '
                          'resolving the curated spec.')
    env.add_argument('--force', action='store_true',
                     help='With --install: remove and rebuild an existing environment.')
    env.add_argument('--dry-run-install', dest='dry_run', action='store_true',
                     help='With --install: print the commands without running them.')

    # -- Cross-stage runtime options ------------------------------------------
    def add_runtime_arguments(target, suppress=False):
        default = argparse.SUPPRESS if suppress else None
        g = target.add_argument_group('runtime')
        g.add_argument('--boltz-bin', default=default, help='Override the boltz executable.')
        g.add_argument('--mmseqs-bin', default=default, help='Override the mmseqs executable.')
        g.add_argument('--tmalign-bin', default=default, help='Override the TMalign executable.')
        g.add_argument('--log', dest='log_file',
                       default=argparse.SUPPRESS if suppress else 'homologic_run.log',
                       help='Run log ("timestamp LEVEL message"). Pass "" to disable.')
        g.add_argument('--overwrite', action='store_true',
                       default=argparse.SUPPRESS if suppress else False,
                       help='Allow stages to delete an existing output folder. Required '
                            'before anything is wiped: re-running boltz-model without it '
                            'will not throw away finished predictions.')
        g.add_argument('--resume', action='store_true',
                       default=argparse.SUPPRESS if suppress else False,
                       help='Reuse existing output folders and skip homologs that already '
                            'have a Boltz model.')
        return target

    add_runtime_arguments(parser)

    # Same flags on every subcommand, with defaults suppressed so that giving
    # them before the subcommand still works and is not overwritten here.
    p_runtime = argparse.ArgumentParser(add_help=False)
    add_runtime_arguments(p_runtime, suppress=True)

    subparsers = parser.add_subparsers(dest='command', help='Pipeline stage to run')

    def add_stage(name, **kw):
        """Create a stage subparser carrying the same split -h / --help."""
        kw.setdefault('add_help', False)
        sub = subparsers.add_parser(name, **kw)
        add_help_arguments(sub)
        return sub

    # -- Shared/reusable argument groups, included only on the subcommands
    #    that actually consume them (parents=[...] defines each flag once
    #    and reuses it verbatim wherever it's relevant) ----------------------
    p_reference_fasta = argparse.ArgumentParser(add_help=False)
    p_reference_fasta.add_argument('-rf', '--reference-fasta', default='reference.fasta',
                                    help='Reference protein sequence (FASTA). REQUIRED input.')

    p_reference_pdb = argparse.ArgumentParser(add_help=False)
    p_reference_pdb.add_argument('-rp', '--reference-pdb', default='reference.pdb',
                                  help='Reference protein structure (PDB). REQUIRED input.')

    p_reference_ligands = argparse.ArgumentParser(add_help=False)
    p_reference_ligands.add_argument('-rl', '--reference-ligands-pdb', default='reference_ligands.pdb',
                                      help='Reference-bound ligand/cofactor coordinates (PDB). REQUIRED input.')

    p_homologs_fasta = argparse.ArgumentParser(add_help=False)
    p_homologs_fasta.add_argument('-hf', '--homologs-fasta', default='homologs.fasta',
                                   help='Multi-sequence FASTA of candidate homologs. REQUIRED input.')

    # -- Stage 1: sequence homology (MMseqs2 + global-alignment fallback) ----
    sp = add_stage(
        'seq-homology', parents=[p_runtime, p_reference_fasta, p_homologs_fasta],
        help='Stage 1: MMseqs2 + pairwise sequence identity vs. reference.')
    sp.add_argument('-o', '--output-folder', default='01_sequence_homology_results',
                     help='Folder for result.m8 / sequence_homology.csv.')
    sp.add_argument('-m', '--substitution-matrix', default='BLOSUM62',
                     help='Substitution matrix for the global-alignment fallback/rescue path.')
    sp.add_argument('-sb', '--skip-seqid-beyond-mmseqs', action='store_true',
                     help='Do NOT run global alignment for sequences MMseqs2 reported no hit for '
                          '(equivalent to seqid_beyond_mmseqs=False).')

    # -- Stage 2: Boltz input generation --------------------------------------
    sp = add_stage(
        'boltz-input', parents=[p_runtime, p_homologs_fasta],
        help='Stage 2: split homologs FASTA into per-protein Boltz input files.')
    add_entity_arguments(sp)
    sp.add_argument('-o', '--output-folder', default='02_boltz_input',
                     help='Folder to write one FASTA per homolog protein.')

    # -- Stage 3: Boltz structure prediction ----------------------------------
    sp = add_stage(
        'boltz-model', parents=[p_runtime],
        help='Stage 3: run `boltz predict` on each Stage-2 FASTA, convert CIF -> protein-only PDB.')
    sp.add_argument('-i', '--input-folder', default='02_boltz_input',
                     help='Folder of per-protein FASTA files from Stage 2.')
    sp.add_argument('-b', '--boltz-results-folder', default='03_boltz_results',
                     help='Folder Boltz writes its raw prediction output into.')
    add_boltz_arguments(sp)
    sp.add_argument('-p', '--protein-only-folder', default='04_structure_models',
                     help='Folder for the final protein-only (ligand-stripped) PDB models.')

    # -- Stage 4: TM-align structural homology --------------------------------
    sp = add_stage(
        'struct-homology', parents=[p_runtime, p_reference_pdb],
        help='Stage 4: TM-align each predicted model against the reference structure.')
    sp.add_argument('-c', '--input-csv-file', default='01_sequence_homology_results/sequence_homology.csv',
                     help="Stage-1 sequence_homology.csv to append rmsd/tm_score columns onto. NOTE: does "
                          "not match Stage 1's own default output folder name -- pass explicitly if you "
                          "used both stages' defaults.")
    sp.add_argument('-i', '--input-folder', default='04_structure_models',
                     help='Folder of predicted protein-only PDB models from Stage 3.')
    sp.add_argument('-o', '--output-folder', default='05_structure_homology_results',
                     help='Folder for structure_homology.csv.')

    # -- Stage 5: superpose full structures onto reference --------------------
    sp = add_stage(
        'superpose-structures', parents=[p_runtime, p_reference_pdb],
        help='Stage 5: PyMOL CA superposition of each model onto the reference structure.')
    sp.add_argument('-i', '--input-folder', default='04_structure_models',
                     help='Folder of predicted protein-only PDB models from Stage 3.')
    sp.add_argument('-o', '--output-folder', default='06_superposed_structures',
                     help='Folder for the superposed PDB models.')

    # -- Stage 6: extract reference binding site ------------------------------
    sp = add_stage(
        'ref-bindingsite', parents=[p_runtime, p_reference_pdb, p_reference_ligands],
        help='Stage 6: extract reference residues near the reference ligand.')
    sp.add_argument('-d', '--distance-cutoff', type=float, default=6.0,
                     help='Distance (Angstrom) from a ligand atom for a residue to count as binding-site.')

    # -- Stage 7: extract homolog binding sites -------------------------------
    sp = add_stage(
        'homolog-bindingsites', parents=[p_runtime, p_reference_ligands],
        help='Stage 7: extract homolog residues near the (superposed) reference ligand.')
    sp.add_argument('-d', '--distance-cutoff', type=float, default=6.0,
                     help='Distance (Angstrom) from a ligand atom for a residue to count as binding-site.')
    sp.add_argument('-i', '--input-folder', default='06_superposed_structures',
                     help='Folder of superposed homolog PDB models from Stage 5.')
    sp.add_argument('-o', '--output-folder', default='07_homolog_bindingsites',
                     help='Folder for the extracted homolog binding-site PDBs.')

    # -- Stage 8: refine binding-site superposition ---------------------------
    sp = add_stage(
        'superpose-bindingsites', parents=[p_runtime, p_reference_pdb],
        help='Stage 8: binding-site-local PyMOL superposition + RMSD quality gates.')
    sp.add_argument('-c', '--csv-file', default='05_structure_homology_results/structure_homology.csv',
                     help='Stage-4 structure_homology.csv to append binding-site RMSD columns onto.')
    sp.add_argument('-i', '--input-folder', default='07_homolog_bindingsites',
                     help='Folder of extracted homolog binding sites from Stage 7.')
    sp.add_argument('-o', '--output-folder', default='08_homolog_bindingsites_superposed',
                     help='Folder for the refined, re-superposed binding-site PDBs.')
    sp.add_argument('-mr', '--min-residues-total', type=int, default=5,
                     help='Minimum binding-site CA count required to attempt superposition.')
    sp.add_argument('-ma', '--min-aligned-atoms', type=int, default=5,
                     help='Minimum CA atoms actually aligned for the RMSD to be trusted (else NaN).')

    # -- Stage 9: binding-site metasequence, similarity, pLDDT ----------------
    sp = add_stage(
        'analyze-bindingsites', parents=[p_runtime, p_reference_pdb],
        help='Stage 9: binding-site metasequence identity + pLDDT summary stats.')
    sp.add_argument('-c', '--csv-file',
                     default='08_homolog_bindingsites_superposed/bindingsite_geometry_homology.csv',
                     help='Stage-8 CSV to append similarity-score/pLDDT columns onto.')
    sp.add_argument('-i', '--input-folder', default='08_homolog_bindingsites_superposed',
                     help='Folder of refined binding-site PDBs from Stage 8.')
    sp.add_argument('-ms', '--metasequence-folder', default='09_bindingsite_metasequences',
                     help='Folder for per-homolog binding-site metasequence FASTA files.')
    sp.add_argument('-sr', '--similarity-results-folder', default='10_bindingsite_similarity_results',
                     help='Folder for bindingsite_homology.csv.')
    sp.add_argument('-b', '--boltz-results-base', default='03_boltz_results',
                     help='Folder of raw Boltz output from Stage 3 (source of pLDDT .npz files).')
    sp.add_argument('-t', '--tolerated-misalignment', type=float, default=1,
                     help='Max CA-CA distance (Angstrom) for a target residue to be matched to a reference one.')
    sp.add_argument('-m', '--substitution-matrix', default='BLOSUM62',
                     help='Substitution matrix for the metasequence pairwise alignment.')

    # -- Stage 10: reference cavity properties --------------------------------
    sp = add_stage(
        'ref-cavity', parents=[p_runtime],
        help='Stage 10: pyKVFinder cavity detection on the reference binding site.')
    sp.add_argument('-i', '--input', default='reference_bindingsite.pdb',
                     help='Reference binding-site PDB (from Stage 6 -- NOT auto-derived from '
                          '--reference-pdb; point this explicitly at "<reference-pdb-stem>_bindingsite.pdb" '
                          'if that differs).')
    sp.add_argument('-o', '--output-folder', default='reference_cavity_analysis',
                     help='Folder for reference cavity results (results.toml, cavities.pdb, barplots.pdf).')
    sp.add_argument('-p', '--probe-out', type=float, default=4.0,
                     help='pyKVFinder outer probe radius (Angstrom).')
    sp.add_argument('-v', '--volume-cutoff', type=float, default=100.0,
                     help='Minimum cavity volume (cubic Angstrom) to be reported.')

    # -- Stage 11: homolog cavity properties + ICP ----------------------------
    sp = add_stage(
        'cavity-properties', parents=[p_runtime, p_reference_pdb],
        help='Stage 11: pyKVFinder cavity detection + ICP shape comparison for each homolog.')
    sp.add_argument('-i', '--input-folder', default='08_homolog_bindingsites_superposed',
                     help='Folder of refined homolog binding-site PDBs from Stage 8.')
    sp.add_argument('-oc', '--output-cavities-folder', default='11_detected_cavities',
                     help='Folder for per-homolog detected-cavity PDBs/results.')
    sp.add_argument('-o', '--output-folder', default='12_cavity_analysis_results',
                     help='Folder for the final bindingsite_cavity_homology.csv.')
    sp.add_argument('-p', '--probe-out', type=float, default=4.0,
                     help='pyKVFinder outer probe radius (Angstrom).')
    sp.add_argument('-v', '--volume-cutoff', type=float, default=100.0,
                     help='Minimum cavity volume (cubic Angstrom) to be reported.')
    sp.add_argument('-c', '--csv-input', default='10_bindingsite_similarity_results/bindingsite_homology.csv',
                     help='Stage-9 CSV to append cavity/ICP columns onto.')
    sp.add_argument('-ov', '--override-reference-bindingsite', default=None,
                     help='Explicit path to a reference binding-site PDB, OVERRIDING the path normally '
                          'derived from --reference-pdb. Leave unset to use that derivation.')
    sp.add_argument('-rc', '--reference-cavity', default=None,
                     help='Explicit path to the reference cavity PDB. Defaults to '
                          '"reference_cavity_analysis/reference_cavities.pdb" (Stage 10 output) if unset.')

    # -- boltz-sweep: an ensemble generator, deliberately outside the pipeline -
    sp = add_stage(
        'boltz-sweep', parents=[p_runtime],
        help='Generate an ensemble: several seeds x many diffusion samples per homolog.')
    sp.add_argument('-i', '--input-folder', default='02_boltz_input',
                    help='Folder of Boltz inputs from boltz-input.')
    sp.add_argument('-o', '--output-folder', default='03_boltz_sweep',
                    help='Root for per-seed Boltz output and the collected ensemble.')
    sp.add_argument('--sweep-seeds', type=int, default=5, metavar='N',
                    help='Number of random seeds per homolog (seeds 1..N).')
    sp.add_argument('--seed-list', metavar='A,B,C',
                    help='Explicit comma-separated seeds, overriding --sweep-seeds.')
    sp.add_argument('--diffusion-samples', type=int, default=100,
                    help='Boltz diffusion samples per seed.')
    sp.add_argument('--allow-msa-bootstrap', action='store_true',
                    help='Permit homologs without a precomputed MSA: Boltz builds each '
                         'one ONCE, on its first seed, and the remaining seeds reuse it '
                         '(one MSA-server call per homolog, never one per seed).')
    sp.add_argument('--ensemble-folder', default=None,
                    help='Where to collect every model. Default: <output-folder>/ensemble.')
    sp.add_argument('-y', '--yes', action='store_true',
                    help='Confirm the GPU cost and actually run. Without it the sweep '
                         'reports what it would generate and stops.')
    sp.add_argument('-n', '--dry-run', action='store_true',
                    help='Report the plan and the guard checks, then exit.')
    add_boltz_arguments(sp, skip=('diffusion-samples',))

    # -- run-all: the full 11-stage pipeline in one call ----------------------
    # Unlike the per-stage subcommands above (whose defaults mirror each
    # HomoLogic method's own signature verbatim, warts and all), run-all uses
    # a single, internally CONSISTENT numbered-folder scheme -- the same one
    # used in the original sequential driver script this class was written
    # for -- so each stage's output lands exactly where the next stage looks
    # for it. Per-stage folder names are not individually overridable here;
    # use the per-stage subcommands instead if you need custom folder names.
    sp = add_stage(
        'run-all', parents=[p_runtime, p_reference_fasta, p_reference_pdb, p_reference_ligands, p_homologs_fasta],
        help='Run all 11 stages (or a --start-at/--stop-after subrange) in sequence.')
    add_entity_arguments(sp)
    add_boltz_arguments(sp)
    sp.add_argument('-d', '--distance-cutoff', type=float, default=6.0,
                     help='Binding-site distance cutoff (Angstrom), used for BOTH the reference and '
                          'homolog binding-site extraction stages.')
    sp.add_argument('-m', '--substitution-matrix', default='BLOSUM62',
                     help='Substitution matrix, used for BOTH sequence homology and binding-site '
                          'metasequence alignment.')
    sp.add_argument('-sb', '--skip-seqid-beyond-mmseqs', action='store_true',
                     help='Do NOT run global alignment for sequences MMseqs2 reported no hit for.')
    sp.add_argument('-t', '--tolerated-misalignment', type=float, default=1,
                     help='Max CA-CA distance (Angstrom) for binding-site metasequence residue matching.')
    sp.add_argument('-mr', '--min-residues-total', type=int, default=5,
                     help='Minimum binding-site CA count required to attempt superposition.')
    sp.add_argument('-ma', '--min-aligned-atoms', type=int, default=5,
                     help='Minimum CA atoms actually aligned for a binding-site RMSD to be trusted.')
    sp.add_argument('-p', '--probe-out', type=float, default=4.0,
                     help='pyKVFinder outer probe radius (Angstrom), used for BOTH cavity stages.')
    sp.add_argument('-v', '--volume-cutoff', type=float, default=100.0,
                     help='Minimum cavity volume (cubic Angstrom), used for BOTH cavity stages.')
    sp.add_argument('-sa', '--start-at', choices=STAGE_ORDER, default=STAGE_ORDER[0],
                     help='First stage to run (skip everything before it, e.g. to resume a failed run).')
    sp.add_argument('-so', '--stop-after', choices=STAGE_ORDER, default=STAGE_ORDER[-1],
                     help='Last stage to run (stop before running anything after it).')
    sp.add_argument('-n', '--dry-run', action='store_true',
                     help='Print which stages would run, in order, and exit without running anything. '
                          'Recommended before a first full run, since safe_create_output_folder() wipes '
                          'any existing folder of the same name.')

    return parser


def _run_pipeline(hl, args):
    """
    Execute STAGE_ORDER[start_at .. stop_after] against `hl`, wiring each
    stage's output folder into the next stage's input folder using the
    fixed numbered-folder scheme documented in HOMOLOGIC_GS_README.md.
    Called only by the run-all subcommand.
    """
    start_idx = STAGE_ORDER.index(args.start_at)
    stop_idx = STAGE_ORDER.index(args.stop_after)
    if start_idx > stop_idx:
        print(f"ERROR: --start-at '{args.start_at}' comes after --stop-after "
              f"'{args.stop_after}' in the pipeline order.", file=sys.stderr)
        sys.exit(1)
    stages = STAGE_ORDER[start_idx:stop_idx + 1]

    if args.dry_run:
        hl.print_progress("Dry run -- stages that would execute, in order:")
        for s in stages:
            print(f"    {s}")
        return

    # Fixed folder chain (see the run-all help text above for why these are
    # not individually configurable here).
    seq_homology_folder = '01_sequence_homology_results'
    boltz_input_folder = '02_boltz_input'
    boltz_results_folder = '03_boltz_results'
    structure_models_folder = '04_structure_models'
    structure_homology_folder = '05_structure_homology_results'
    superposed_structures_folder = '06_superposed_structures'
    homolog_bindingsites_folder = '07_homolog_bindingsites'
    bindingsites_superposed_folder = '08_homolog_bindingsites_superposed'
    metasequence_folder = '09_bindingsite_metasequences'
    similarity_results_folder = '10_bindingsite_similarity_results'
    reference_cavity_folder = 'reference_cavity_analysis'
    detected_cavities_folder = '11_detected_cavities'
    cavity_analysis_results_folder = '12_cavity_analysis_results'
    # Matches exactly what extract_reference_binding_site() itself derives
    # from hl.reference_pdb, so this stays correct for any --reference-pdb.
    reference_bindingsite_pdb = f"{os.path.splitext(hl.reference_pdb)[0]}_bindingsite.pdb"

    entities = build_entities(args)
    msa_map = ({} if not any(st in BOLTZ_STAGES for st in stages)
               else msa_preflight_for(args, hl, args.homologs_fasta, boltz_input_folder))

    if 'seq-homology' in stages:
        hl.calculate_sequence_homology(
            output_folder=seq_homology_folder,
            substitution_matrix_name=args.substitution_matrix,
            seqid_beyond_mmseqs=not args.skip_seqid_beyond_mmseqs,
        )

    if 'boltz-input' in stages:
        hl.generate_boltz_input(
            entities=entities,
            msa_map=msa_map,
            fmt=args.format,
            yaml_template=args.yaml_template,
            yaml_target_chain=args.yaml_target_chain,
            validate=not args.no_validate,
            output_folder=boltz_input_folder,
        )

    if 'boltz-model' in stages:
        hl.perform_boltz_modeling(
            input_folder=boltz_input_folder,
            boltz_results_folder=boltz_results_folder,
            boltz_args=collect_boltz_args(args),
            use_msa_server=not args.no_msa_server,
            protein_only_folder=structure_models_folder,
        )

    if 'struct-homology' in stages:
        hl.calculate_structure_homology(
            input_csv_file=os.path.join(seq_homology_folder, 'sequence_homology.csv'),
            input_folder=structure_models_folder,
            output_folder=structure_homology_folder,
        )

    if 'superpose-structures' in stages:
        hl.superpose_structures(
            input_folder=structure_models_folder,
            output_folder=superposed_structures_folder,
        )

    if 'ref-bindingsite' in stages:
        hl.extract_reference_binding_site(
            distance_cutoff=args.distance_cutoff,
        )

    if 'homolog-bindingsites' in stages:
        hl.extract_homolog_binding_sites(
            distance_cutoff=args.distance_cutoff,
            input_folder=superposed_structures_folder,
            output_folder=homolog_bindingsites_folder,
        )

    if 'superpose-bindingsites' in stages:
        hl.superpose_binding_sites(
            csv_file=os.path.join(structure_homology_folder, 'structure_homology.csv'),
            input_folder=homolog_bindingsites_folder,
            output_folder=bindingsites_superposed_folder,
            min_residues_total=args.min_residues_total,
            min_aligned_atoms=args.min_aligned_atoms,
        )

    if 'analyze-bindingsites' in stages:
        hl.analyze_binding_sites(
            csv_file=os.path.join(bindingsites_superposed_folder, 'bindingsite_geometry_homology.csv'),
            input_folder=bindingsites_superposed_folder,
            bindingsite_metasequence_folder=metasequence_folder,
            bindingsite_similarity_results=similarity_results_folder,
            boltz_results_base=boltz_results_folder,
            tolerated_misalignment=args.tolerated_misalignment,
            substitution_matrix_name=args.substitution_matrix,
        )

    if 'ref-cavity' in stages:
        hl.analyze_reference_cavity_properties(
            input=reference_bindingsite_pdb,
            output_results_folder=reference_cavity_folder,
            probe_out=args.probe_out,
            volume_cutoff=args.volume_cutoff,
        )

    if 'cavity-properties' in stages:
        hl.analyze_cavity_properties(
            input_folder=bindingsites_superposed_folder,
            output_cavities_folder=detected_cavities_folder,
            output_results_folder=cavity_analysis_results_folder,
            probe_out=args.probe_out,
            volume_cutoff=args.volume_cutoff,
            csv_input=os.path.join(similarity_results_folder, 'bindingsite_homology.csv'),
        )

    hl.print_success(f"run-all complete -- stages executed: {', '.join(stages)}")


# Stages that consume a Boltz input file. A custom MSA is only meaningful --
# and only validated -- when one of these is going to run: `homologic2.py
# seq-homology` must stay usable with no MSA in sight.
BOLTZ_STAGES = ('boltz-input', 'boltz-model', 'boltz-sweep')


def msa_preflight_for(args, hl, homologs_fasta, output_root):
    """
    Resolve, convert and validate custom MSAs before any stage executes.
    Returns {homolog_id: boltz csv path}, or {} when no MSA flag was given.
    """
    if not (getattr(args, 'use_custom_msa', None) or getattr(args, 'reuse_msa_from', None)):
        return {}
    if not os.path.exists(homologs_fasta):
        sys.exit(f"ERROR: --use-custom-msa/--reuse-msa-from needs the homologs FASTA, "
                 f"but '{homologs_fasta}' does not exist.")
    records = read_fasta_sequences(homologs_fasta)
    # Converted MSAs live BESIDE the stage's output folder, not inside it:
    # the preflight runs before generate_boltz_input() prepares that folder,
    # and safe_create_output_folder() would otherwise either trip its own
    # overwrite guard on them or delete them.
    converted = os.path.join(os.path.dirname(os.path.abspath(output_root)), 'msa')
    return preflight_msa(records, args, converted, log=hl.print_progress)


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    # -- Environment management first: these must work in an interpreter
    #    that has none of the scientific stack installed. ------------------
    if args.install_help:
        print(INSTALL_HELP)
        return 0
    if args.install:
        return do_install(args)
    if args.check:
        return do_check(args)

    if not args.command:
        parser.print_help()
        return 1

    # Every remaining path runs a pipeline stage, so the heavy dependencies
    # are needed from here on.
    _load_science_stack()

    run_start = time.time()

    # Only the reference/homolog paths relevant to the chosen stage were
    # defined on its subparser above; anything else falls back to
    # HomoLogic's own constructor defaults via getattr().
    hl = HomoLogic(
        reference_fasta=getattr(args, 'reference_fasta', 'reference.fasta'),
        reference_pdb=getattr(args, 'reference_pdb', 'reference.pdb'),
        reference_ligands_pdb=getattr(args, 'reference_ligands_pdb', 'reference_ligands.pdb'),
        homologs_fasta=getattr(args, 'homologs_fasta', 'homologs.fasta'),
        boltz_bin=args.boltz_bin,
        mmseqs_bin=args.mmseqs_bin,
        tmalign_bin=args.tmalign_bin,
        log_file=args.log_file or None,
        overwrite=args.overwrite,
        resume=args.resume,
    )

    if args.command == 'seq-homology':
        hl.calculate_sequence_homology(
            output_folder=args.output_folder,
            substitution_matrix_name=args.substitution_matrix,
            seqid_beyond_mmseqs=not args.skip_seqid_beyond_mmseqs,
        )

    elif args.command == 'boltz-input':
        msa_map = msa_preflight_for(args, hl, args.homologs_fasta, args.output_folder)
        hl.generate_boltz_input(
            entities=build_entities(args),
            msa_map=msa_map,
            output_folder=args.output_folder,
            fmt=args.format,
            yaml_template=args.yaml_template,
            yaml_target_chain=args.yaml_target_chain,
            validate=not args.no_validate,
        )

    elif args.command == 'boltz-model':
        hl.perform_boltz_modeling(
            input_folder=args.input_folder,
            boltz_results_folder=args.boltz_results_folder,
            protein_only_folder=args.protein_only_folder,
            boltz_args=collect_boltz_args(args),
            use_msa_server=not args.no_msa_server,
        )

    elif args.command == 'boltz-sweep':
        if args.seed_list:
            try:
                seeds = [int(x) for x in args.seed_list.split(',') if x.strip()]
            except ValueError:
                sys.exit(f"ERROR: --seed-list '{args.seed_list}' is not a comma-separated "
                         f"list of integers.")
            if not seeds:
                sys.exit("ERROR: --seed-list is empty.")
        else:
            if args.sweep_seeds < 1:
                sys.exit("ERROR: --sweep-seeds must be at least 1.")
            seeds = list(range(1, args.sweep_seeds + 1))
        hl.perform_boltz_sweep(
            input_folder=args.input_folder,
            output_folder=args.output_folder,
            seeds=seeds,
            diffusion_samples=args.diffusion_samples,
            boltz_args=collect_boltz_args(args, skip_seed=True,
                                          skip_diffusion_samples=True),
            allow_msa_bootstrap=args.allow_msa_bootstrap,
            assume_yes=args.yes,
            ensemble_folder=args.ensemble_folder,
            dry_run=args.dry_run,
        )

    elif args.command == 'struct-homology':
        hl.calculate_structure_homology(
            input_csv_file=args.input_csv_file,
            input_folder=args.input_folder,
            output_folder=args.output_folder,
        )

    elif args.command == 'superpose-structures':
        hl.superpose_structures(
            input_folder=args.input_folder,
            output_folder=args.output_folder,
        )

    elif args.command == 'ref-bindingsite':
        hl.extract_reference_binding_site(
            distance_cutoff=args.distance_cutoff,
        )

    elif args.command == 'homolog-bindingsites':
        hl.extract_homolog_binding_sites(
            distance_cutoff=args.distance_cutoff,
            input_folder=args.input_folder,
            output_folder=args.output_folder,
        )

    elif args.command == 'superpose-bindingsites':
        hl.superpose_binding_sites(
            csv_file=args.csv_file,
            input_folder=args.input_folder,
            output_folder=args.output_folder,
            min_residues_total=args.min_residues_total,
            min_aligned_atoms=args.min_aligned_atoms,
        )

    elif args.command == 'analyze-bindingsites':
        hl.analyze_binding_sites(
            csv_file=args.csv_file,
            input_folder=args.input_folder,
            bindingsite_metasequence_folder=args.metasequence_folder,
            bindingsite_similarity_results=args.similarity_results_folder,
            boltz_results_base=args.boltz_results_base,
            tolerated_misalignment=args.tolerated_misalignment,
            substitution_matrix_name=args.substitution_matrix,
        )

    elif args.command == 'ref-cavity':
        hl.analyze_reference_cavity_properties(
            input=args.input,
            output_results_folder=args.output_folder,
            probe_out=args.probe_out,
            volume_cutoff=args.volume_cutoff,
        )

    elif args.command == 'cavity-properties':
        hl.analyze_cavity_properties(
            input_folder=args.input_folder,
            output_cavities_folder=args.output_cavities_folder,
            output_results_folder=args.output_folder,
            probe_out=args.probe_out,
            volume_cutoff=args.volume_cutoff,
            csv_input=args.csv_input,
            reference_pdb=args.override_reference_bindingsite,
            reference_cavity=args.reference_cavity,
        )

    elif args.command == 'run-all':
        _run_pipeline(hl, args)

    hl.print_success(f"{args.command} finished in {hms_string(time.time() - run_start)}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
