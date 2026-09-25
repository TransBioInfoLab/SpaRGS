# sGSVA: spatial gene set variation analysis

This repository contains code and reproducible analysis workflows for the sGSVA manuscript, 
including simulations, real-data analyses, figure generation, and supplementary tables.

## Project overview

sGSVA is a reference-free spatial pathway scoring framework for estimating
spatially coherent pathway activity across spot-based and cell-resolution spatial transcriptomics platforms. The pipeline integrates:

- Gaussian kernel smoothing on spatial neighborhoods
- correlation-based pathway directionality assignment
- signed and rank-based pathway scoring
- spatial significance testing using Moran's I / Geary's C
- permutation-based null models

The repository reproduces analyses for:

- breast cancer Visium
- human DLPFC Visium
- Xenium FFPE breast cancer
- simulation benchmarks and robustness

## Repository structure

- `scripts/` — dataset- and analysis-specific workflows
- `R/` — reusable helper functions for scoring, permutation, and evaluation
- `data/` — raw, processed, and example data
- `results/` — generated figures, machine-readable tables, summary CSVs, and session information
- `docs/` — manifest and supporting documentation
- `manuscript/` — manuscript text, formatted tables, and compiled supplementary files

## Requirements

This project was developed in R 4.4.1.

Major R package dependencies include:

- Seurat
- GSVA
- msigdbr
- org.Hs.eg.db
- spdep
- sf
- FNN
- hexbin
- Matrix
- matrixStats
- ggplot2
- patchwork
- dplyr
- tidyr
- readr
- spatialLIBD

Additional package-specific system dependencies may be required for some spatial packages, depending on the operating system.

### Install R dependencies

The code was tested using R 4.4.1.

```r
install.packages(c(
  "Matrix",
  "matrixStats",
  "spdep",
  "sf",
  "FNN",
  "msigdbr",
  "Seurat",
  "hexbin",
  "ggplot2",
  "patchwork",
  "dplyr",
  "tidyr",
  "readr"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "GSVA",
  "org.Hs.eg.db",
  "spatialLIBD"
))
```

## Session information and random seeds

Session information is saved in `results/session_info.txt`. A curated list of packages used in the sGSVA workflow is saved in `results/package_versions_sgsva.csv`. 
All stochastic analyses used fixed random seeds. Unless otherwise stated, manuscript analyses used seed = 1.

## Data access

### Public datasets

| Dataset | Source | Expected local path |
|---|---|---|
| Breast Visium | Zenodo accession 4739739 | `data/breast_visium/` |
| DLPFC Visium | `spatialLIBD` package | downloaded programmatically |
| Xenium FFPE breast cancer | 10x Genomics public Xenium FFPE human breast cancer dataset with 5k Human Pan-Tissue and Pathways panel plus 100 Custom Genes| `data/xenium/` |
Exact dataset accessions and vendor download links are provided in Supplementary Table S16.

### Breast Visium
Breast Visium input files are expected in `data/breast_visium/`.  
These include:
- filtered count matrices
- spatial image files and scale factors
- metadata and tissue annotations
- GSE176078 reference files

### Xenium FFPE breast cancer
Xenium input files are expected in `data/xenium/`.  
These include:
- cell-by-gene transcript count matrix
- cell centroid coordinates
- cell metadata and panel information
- Xenium gene-group marker table

### DLPFC
The DLPFC dataset is retrieved through the `spatialLIBD` package within the analysis scripts. No manual download into `data/` is required.

### Notes
Public raw datasets are expected to be stored locally under the data/ directory in the folder structure described above. 
Because these files can be large, they are not intended to be version-controlled within the repository. 
Users should obtain them from their original public sources and place them in the expected folders before running the full workflow.


## Key outputs

- Main and supplementary figures: `results/figures/`
- Machine-readable table sources: `results/tables/`
- Analysis-specific results: `results/simulation/`, `results/robustness/`, `results/breast_visium/`, `results/DLPFC/`, `results/xenium/`
- Session information: `results/session_info.txt`
- Script-to-output manifest: `docs/script_to_output_manifest.xlsx`


## How to run

### Full manuscript pipeline
To run the full analysis pipeline:

```bash
bash run_all.sh
```

### Minimal reproducibility example

A lightweight fixed-input example is provided to verify that the core sGSVA
implementation can run end-to-end from a clean installation. No external spatial
transcriptomics dataset download is required.

The example uses a small frozen input derived from the manuscript simulation
framework. The committed file `data/example/minimal_example_input.rds` 
contains the count matrix, spatial coordinates, representative KEGG Legacy gene
sets, simulated truth fields, and pathway-design information used for the test.
Using a fixed input avoids cross-platform differences introduced by regenerating
random simulations on different systems.

Run the example from the root of the repository:

```bash
Rscript scripts/run_minimal_example.R
```

### Large figure-source files

Some manuscript figures require large intermediate analysis objects, such as Seurat objects, full pathway-score matrices, or permutation result objects. 
These files are not included in the lightweight GitHub repository. If they are missing, run the corresponding full analysis script first. 





