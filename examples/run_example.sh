#!/usr/bin/env bash
# Offline walkthrough of homologic2.py -- no GPU and no network needed.
# Everything here uses the synthetic sequences in this folder.
set -euo pipefail
cd "$(dirname "$0")"
H=../homologic2.py

echo "### 1. Is the environment sane?"
python "$H" --check || echo "  (run: python $H --install)"

echo
echo "### 2. Boltz inputs, protein only"
python "$H" boltz-input -hf homologs.faa -o out/01_plain --overwrite

echo
echo "### 3. Co-folding: a peptide ligand and a calcium ion"
#     --peptide becomes a real protein chain with MSA lookup disabled;
#     --ccd CA is preferred over the SMILES '[Ca+2]'.
python "$H" boltz-input -hf homologs.faa -o out/02_ligands \
    --peptide KAFVQWLIAG --ccd CA --overwrite

echo
echo "### 4. Reuse one alignment as every homolog's MSA"
#     Each homolog's own row becomes its query; the alignment is projected
#     onto that homolog's ungapped columns. Expect a shallow-MSA warning:
#     4 sequences is far below what an MSA server would return.
python "$H" boltz-input -hf homologs.faa -o out/03_msa \
    --peptide KAFVQWLIAG --ccd CA --use-custom-msa homologs_aligned.afa --overwrite

echo
echo "### 5. YAML input, for constraints/affinity"
python "$H" boltz-input -hf homologs.faa -o out/04_yaml \
    --format yaml --yaml-template pocket_template.yaml --overwrite

echo
echo "### 6. What a sweep would do (nothing is run without --yes)"
python "$H" boltz-sweep -i out/03_msa -o out/05_sweep \
    --sweep-seeds 3 --diffusion-samples 25 --dry-run

echo
echo "Done. Generated inputs are under ./out/ ; converted MSAs under ./out/msa/"
