# HomoLogic
[![Python](https://img.shields.io/badge/python-3.14%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**HomoLogic** is a modular Python workflow for **integrative enzyme homology analysis**. 

It combines **sequence similarity, structural comparison, binding-site analysis, and cavity characterization** into a unified computational pipeline designed for **enzyme discovery and comparative enzyme analysis**.

![Uploading image.png…]()


The pipeline generates a comprehensive set of **sequence-, structure-, and binding-site–derived descriptors** that can be used for:

- enzyme homology evaluation
- active-site comparison
- enzyme library analysis
- machine-learning driven enzyme discovery

---

# Overview

HomoLogic evaluates enzyme homologs by comparing a **reference enzyme ("query")** against a **pool of homologous sequences ("haystack")**.

The pipeline integrates:

- sequence homology
- structure prediction
- structural alignment
- binding-site extraction
- binding-site sequence similarity ("metasequences")
- cavity geometry and physicochemical properties

The result is a **unified descriptor table** that quantitatively characterizes enzyme homologs.

---

# Inputs

The workflow requires two input categories.

## Query (reference enzyme)

The reference enzyme defines the structural and functional template.

Required files:

```

reference.fasta
reference.pdb
reference_ligands.pdb

```

Optional:

- ligand SMILES strings for structure co-folding

| File | Description |
|----|----|
| `reference.fasta` | Reference enzyme amino acid sequence |
| `reference.pdb` | Reference protein structure |
| `reference_ligands.pdb` | Coordinates of bound ligand(s) |
| SMILES strings | Ligand definitions for structure prediction |

---

## Haystack (homolog candidate sequences)

```

homologs.fasta

```

This FASTA file contains **putatively homologous sequences** belonging to the same enzyme family or structural fold.

These sequences may originate from:

- sequence similarity searches (BLAST, MMseqs2)
- structure similarity searches (Foldseek)
- protein language model searches (pLM-BLAST)
- enzyme mining tools (BioCatNet, EnzymeMiner)
- orthology detection tools (OrthoFinder, SHOOT)
- curated enzyme libraries

---

# Pipeline Workflow

The HomoLogic workflow consists of the following steps.

---

## 1. Sequence homology analysis

Sequence similarity between the reference enzyme and all homologs is calculated using **MMseqs2**.

Outputs:

```

01_sequence_homology_results/sequence_homology.csv

```

Contains:

- sequence identity
- alignment length
- mismatches
- E-value
- bit score

---

## 2. Structure prediction

Structures for all homolog sequences are generated using **Boltz-2**.

The pipeline can optionally **co-fold ligands** to maintain the spatial context of the active site.

Input FASTA files for Boltz are generated automatically.

Outputs:

```

03_boltz_results/
04_structure_models/

```

Protein-only PDB files are extracted for downstream analysis.

### Alternative structure sources

If Boltz modeling is not desired, structures may also be obtained from:

- Protein Data Bank (PDB)
- AlphaFold DB
- homology modeling
- ESMFold

Requirement:

```

protein-only PDB files
filename must match FASTA header

```

Example:

```

> example
> example.pdb

```

---

## 3. Structural homology analysis

Structural similarity between reference and homolog models is quantified using **TM-align**.

Metrics computed:

- TM-score
- RMSD

Structures are then **superposed onto the reference structure using PyMOL** to ensure a consistent coordinate system.

Outputs:

```

05_structure_homology_results/structure_homology.csv
06_superposed_structures/

```

---

## 4. Binding-site extraction

Binding sites are defined as residues located within a **6 Å distance** from the reference ligand coordinates.

Binding sites are extracted from:

- the reference structure
- all homolog structures

Outputs:

```

reference_bindingsite.pdb
07_homolog_bindingsites/

```

---

## 5. Binding-site structural comparison

Homolog binding sites are superposed onto the reference binding site.

Metrics computed:

- binding-site RMSD
- number of aligned Cα atoms
- aligned residue fraction

Outputs:

```

08_homolog_bindingsites_superposed/
bindingsite_geometry_homology.csv

```

---

## 6. Binding-site metasequence analysis

To quantify sequence variation within structurally equivalent active-site regions, HomoLogic constructs **binding-site metasequences**.

These sequences represent spatially corresponding residues between reference and homolog binding sites.

Residues are considered equivalent if their **Cα atoms are within 1 Å distance** after superposition.

Metasequences are globally aligned using **scikit-bio** with substitution matrices from **Biopython** (default: **BLOSUM62**).

Outputs:

```

09_bindingsite_metasequences/
10_bindingsite_similarity_results/bindingsite_homology.csv

```

Metrics:

- normalized binding-site similarity score
- number of non-equivalent residues
- pLDDT statistics (full protein and binding site)

---

## 7. Cavity detection and analysis

Binding-site cavities are detected using **pyKVFinder**.

Default parameters:

```

probe radius = 4 Å
minimum cavity volume = 100 Å³

```

Descriptors calculated:

- cavity volume
- surface area
- average depth
- hydropathy
- residue class frequencies

Outputs:

```

reference_cavity_analysis/
11_detected_cavities/
12_cavity_analysis_results/

```

---

## 8. Cavity geometry comparison

Cavity geometries are compared using **point-cloud alignment** via **Open3D**.

The **Iterative Closest Point (ICP)** algorithm is used to align cavities to the reference pocket.

Metrics:

- cavity RMSD
- number of aligned cavity points

---

# Final Output

All descriptors are combined into the final table:

```

bindingsite_cavity_homology.csv

````

This dataset contains:

- sequence similarity descriptors
- structural similarity descriptors
- binding-site similarity scores
- structural confidence metrics
- cavity geometry descriptors
- cavity physicochemical properties

This unified feature set enables **quantitative enzyme comparison and machine-learning based enzyme discovery**.

---

# Example Workflow 
(please see detailed examples at https://github.com/Buller-Lab/Ssal-KRED_orthologs)

```python
from homologic import HomoLogic

# Define reference enzyme sequence, structure, and ligands, as well as homologous sequences (search space)
hl = HomoLogic(reference_fasta="reference.fasta", reference_pdb="reference.pdb", reference_ligands_pdb="reference_ligands.pdb", homologs_fasta="../Ssal-KRED_homologs_SwissProt/SwissProt_homologs/hits.fasta")

# Calculate sequence homology between reference and homologs via MMseqs2
hl.calculate_sequence_homology(output_folder='01_sequence_homology_results')

# Generate input files for Boltz-2 cofolding
hl.generate_boltz_input(smiles_code=["C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O"], output_folder='02_boltz_input')

# Perform Boltz-2 cofolding and generate protein only structure models
hl.perform_boltz_modeling(input_folder = "02_boltz_input", boltz_results_folder='03_boltz_results', protein_only_folder='04_structure_models')

# Calculate structure homology between reference and homologs via TM-align
hl.calculate_structure_homology(input_csv_file='01_sequence_homology_results/sequence_homology.csv', input_folder='04_structure_models', output_folder='05_structure_homology_results')

# Superpose all structure models onto the reference structure via PyMol for consistent spatial orientation
hl.superpose_structures(input_folder='04_structure_models', output_folder='06_superposed_structures')

# Based on given distance from reference ligand, extract binding sites from reference structure
hl.extract_reference_binding_site(distance_cutoff=6.0)

# Based on given distance from reference ligand, extract binding sites from superposed homolog structures
hl.extract_homolog_binding_sites(distance_cutoff=6.0, input_folder='06_superposed_structures', output_folder='07_homolog_bindingsites')

# Superpose all homolog binding sites onto the reference binding site via PyMol
hl.superpose_binding_sites(csv_file='05_structure_homology_results/structure_homology.csv', input_folder='07_homolog_bindingsites', output_folder='08_homolog_bindingsites_superposed')

# Analyze local "metasequences" of binding sites 
hl.analyze_binding_sites(csv_file='08_homolog_bindingsites_superposed/bindingsite_geometry_homology.csv', input_folder='08_homolog_bindingsites_superposed', bindingsite_metasequence_folder='09_bindingsite_metasequences', bindingsite_similarity_results='10_bindingsite_similarity_results', tolerated_misalignment=1.0)

# Analyze cavity properties of reference binding site
hl.analyze_reference_cavity_properties(input='reference_bindingsite.pdb', output_results_folder='reference_cavity_analysis', probe_out=4.0, volume_cutoff=100.0)

# And finally, detect cavities analyze their properties
hl.analyze_cavity_properties(input_folder='08_homolog_bindingsites_superposed', output_cavities_folder='11_detected_cavities', output_results_folder='12_cavity_analysis_results', probe_out=4.0, volume_cutoff=100.0, csv_input='10_bindingsite_similarity_results/bindingsite_homology.csv')
````

---

# Dependencies

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
* NumPy
* Pandas

For software versions, please anaconda .yml file
---

# Citation

If you use HomoLogic in your research, please cite:

```
Stockinger et al. Iterative ortholog mining enables data-driven discovery of stereoselective ketoreductases
```

(Manuscript in preparation)

---


